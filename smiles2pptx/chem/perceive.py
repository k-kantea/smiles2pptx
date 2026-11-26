"""Perception (the second stage): working out what the SMILES leaves implicit.

perceive() takes the molecule from the reader and returns a copy with:

- the rings: the smallest set of smallest rings (Molecule.rings)
- the Kekule structure: each aromatic bond becomes a single or a double bond (Bond.order)
- the hydrogen count of every atom (Atom.h_count)
- the aromatic rings, by the 4n+2 rule on the Kekule structure (Molecule.aromatic_rings)
- which "@" and "@@" are on real stereocentres (Atom.stereocentre)
- the E/Z configuration of double bonds marked with / and \\ (Bond.cis_trans)

An aromatic system that has no Kekule structure (c1cccc1) raises PerceiveError.
"""
from __future__ import annotations

import copy
from collections import deque
from typing import Dict, List, Optional, Sequence, Set, Tuple

from .elements import NORMAL_VALENCE, normal_valences
from .molecule import IMPLICIT_H, Atom, Molecule

ORDER = {"single": 1, "double": 2, "triple": 3, "quadruple": 4}
# Rings with more atoms are not tested for aromaticity (the inner ring of a porphyrin).
LARGEST_AROMATIC_RING = 7
# Elements that can be part of an aromatic ring.
AROMATIC_ELEMENTS = {"B", "C", "N", "O", "Si", "P", "S", "As", "Se", "Te"}
# An atom with no double bond donates its lone pair when it has this many bonds (counting hydrogens).
LONE_PAIR_BONDS = {"N": 3, "P": 3, "As": 3, "O": 2, "S": 2, "Se": 2, "Te": 2}
# A ring carbon with a double bond to one of these outside the rings contributes no electrons (C=O).
PULLS_ELECTRONS = {"N", "O", "S", "Se", "Te"}
# Atoms whose lone pair makes them a stereocentre with three neighbours (sulfoxides, phosphines).
LONE_PAIR_STEREO = {"S", "Se", "P", "As"}
# Double bonds in smaller rings are not drawn with E/Z (as in RDKit).
SMALLEST_RING_WITH_EZ = 8


class PerceiveError(ValueError):
    """Error raised for a molecule whose aromatic atoms cannot be given a Kekule structure."""

    def __init__(self, molecule: Molecule, atoms: Sequence[int], reason: str):
        self.atoms = list(atoms)
        places = ", ".join(f"{molecule.atoms[i].element.lower() if molecule.atoms[i].aromatic else molecule.atoms[i].element}"
                           f" at position {molecule.atoms[i].position + 1}" for i in self.atoms)
        super().__init__(f"{reason}: {places}\n  {molecule.smiles}")


def perceive(molecule: Molecule) -> Molecule:
    """Return a copy of the molecule with the rings, the Kekule structure, the hydrogens
    and the real stereo marks filled in."""
    mol = copy.deepcopy(molecule)
    adjacency = _adjacency(mol)
    mol.rings = smallest_rings(len(mol.atoms), [(b.begin, b.end) for b in mol.bonds], len(mol.fragments))
    _kekulize(mol, adjacency)
    for i, atom in enumerate(mol.atoms):
        atom.h_count = _hydrogens(mol, i, adjacency)
    mol.aromatic_rings = aromatic_rings(mol)
    ranks = _symmetry_classes(mol, adjacency)
    _stereocentres(mol, adjacency, ranks)
    _cis_trans(mol, adjacency, ranks)
    return mol


def _adjacency(mol: Molecule) -> List[List[Tuple[int, int]]]:
    """Return a list of (neighbour, bond index) for each atom, in the order of the bonds."""
    adjacency: List[List[Tuple[int, int]]] = [[] for _ in mol.atoms]
    for k, b in enumerate(mol.bonds):
        adjacency[b.begin].append((b.end, k))
        adjacency[b.end].append((b.begin, k))
    return adjacency


# --- rings ------------------------------------------------------------------------

def smallest_rings(n_atoms: int, edges: Sequence[Tuple[int, int]], n_parts: int) -> List[List[int]]:
    """Return the smallest set of smallest rings (a minimum cycle basis).

    Horton's method: every ring in the set is a shortest path from an atom v to an atom x,
    a bond x-y and a shortest path from y back to v. The candidates are sorted by size and
    kept when they are independent of the rings already kept (as sets of bonds, added modulo 2).
    """
    wanted = len(edges) - n_atoms + n_parts
    if wanted <= 0:
        return []
    neighbours: List[List[Tuple[int, int]]] = [[] for _ in range(n_atoms)]
    for k, (a, b) in enumerate(edges):
        neighbours[a].append((b, k))
        neighbours[b].append((a, k))
    candidates: Dict[int, List[int]] = {}
    for v in range(n_atoms):
        if not neighbours[v]:
            continue
        parent: Dict[int, Tuple[int, int]] = {v: (-1, -1)}   # atom -> (previous atom, bond)
        queue = deque([v])
        while queue:
            x = queue.popleft()
            for y, k in neighbours[x]:
                if y not in parent:
                    parent[y] = (x, k)
                    queue.append(y)

        def path(x):
            atoms, bonds = [x], 0
            while x != v:
                x, k = parent[x]
                atoms.append(x)
                bonds |= 1 << k
            return atoms[::-1], bonds

        for k, (x, y) in enumerate(edges):
            if x not in parent or parent[x][1] == k or parent[y][1] == k:
                continue
            to_x, bonds_x = path(x)
            to_y, bonds_y = path(y)
            if set(to_x) & set(to_y) != {v}:
                continue   # the two paths meet before v: not a simple ring
            bonds = bonds_x | bonds_y | (1 << k)
            if bonds not in candidates:
                candidates[bonds] = to_x + to_y[:0:-1]
    rings: List[List[int]] = []
    basis: Dict[int, int] = {}   # highest bond -> ring (as a bit set of bonds); each highest bond differs
    for bonds, atoms in sorted(candidates.items(), key=lambda item: (len(item[1]), item[1])):
        reduced = bonds
        while reduced:
            top = reduced.bit_length() - 1
            if top not in basis:
                break
            reduced ^= basis[top]
        if reduced:
            basis[reduced.bit_length() - 1] = reduced
            rings.append(atoms)
            if len(rings) == wanted:
                break
    return rings


# --- Kekule structure ---------------------------------------------------------------

def _kekulize(mol: Molecule, adjacency) -> None:
    """Set the order of every bond; aromatic bonds become 1 or 2.

    An aromatic atom takes a double bond when its normal valence has room for one more
    bond than it has (counting aromatic bonds as single): c in benzene, n in pyridine,
    [n+] in pyridinium, but not [nH], o, or n with three bonds. The double bonds are then
    a perfect matching of those atoms over the aromatic bonds.
    """
    for b in mol.bonds:
        b.order = ORDER.get(b.kind, 1)
    needs = [_takes_double_bond(mol, i, adjacency) for i in range(len(mol.atoms))]
    nodes = [i for i, need in enumerate(needs) if need]
    if not nodes:
        return
    index = {atom: k for k, atom in enumerate(nodes)}
    graph: List[List[int]] = [[] for _ in nodes]
    bond_of: Dict[Tuple[int, int], int] = {}
    for k, b in enumerate(mol.bonds):
        if b.kind == "aromatic" and needs[b.begin] and needs[b.end]:
            i, j = index[b.begin], index[b.end]
            graph[i].append(j)
            graph[j].append(i)
            bond_of[(min(i, j), max(i, j))] = k
    match = _maximum_matching(graph)
    unmatched = [nodes[i] for i, m in enumerate(match) if m == -1]
    if unmatched:
        raise PerceiveError(mol, unmatched, "cannot kekulize: no double bond can be placed on these aromatic atoms")
    for i, j in enumerate(match):
        if i < j:
            mol.bonds[bond_of[(i, j)]].order = 2


def _takes_double_bond(mol: Molecule, i: int, adjacency) -> bool:
    atom = mol.atoms[i]
    if not atom.aromatic:
        return False
    used = sum(1 if mol.bonds[k].kind == "aromatic" else ORDER.get(mol.bonds[k].kind, 1) for _, k in adjacency[i])
    used += atom.hydrogens or 0
    valence = _lowest_at_least(normal_valences(atom.element, atom.charge), used)
    return valence is not None and valence > used


def _lowest_at_least(valences: Sequence[int], used: int) -> Optional[int]:
    return next((v for v in sorted(valences) if v >= used), None)


def _maximum_matching(graph: List[List[int]]) -> List[int]:
    """Return a maximum matching in a general graph (Edmonds' blossom method): the partner
    of each node, or -1."""
    n = len(graph)
    match = [-1] * n
    # Start from a greedy matching; the search below only has to fix what it leaves unmatched.
    for v in range(n):
        if match[v] == -1:
            for w in graph[v]:
                if match[w] == -1:
                    match[v], match[w] = w, v
                    break

    def augmenting_path(root):
        used = [False] * n
        parent = [-1] * n
        base = list(range(n))
        used[root] = True
        queue = [root]

        def common_base(a, b):
            seen = [False] * n
            while True:
                a = base[a]
                seen[a] = True
                if match[a] == -1:
                    break
                a = parent[match[a]]
            while True:
                b = base[b]
                if seen[b]:
                    return b
                b = parent[match[b]]

        def mark(v, b, child, blossom):
            while base[v] != b:
                blossom[base[v]] = blossom[base[match[v]]] = True
                parent[v] = child
                child = match[v]
                v = parent[match[v]]

        head = 0
        while head < len(queue):
            v = queue[head]
            head += 1
            for w in graph[v]:
                if base[v] == base[w] or match[v] == w:
                    continue
                if w == root or (match[w] != -1 and parent[match[w]] != -1):
                    b = common_base(v, w)
                    blossom = [False] * n
                    mark(v, b, w, blossom)
                    mark(w, b, v, blossom)
                    for i in range(n):
                        if blossom[base[i]]:
                            base[i] = b
                            if not used[i]:
                                used[i] = True
                                queue.append(i)
                elif parent[w] == -1:
                    parent[w] = v
                    if match[w] == -1:
                        return w, parent
                    used[match[w]] = True
                    queue.append(match[w])
        return -1, parent

    for root in range(n):
        if match[root] == -1:
            v, parent = augmenting_path(root)
            while v != -1:
                p = parent[v]
                after = match[p]
                match[v], match[p] = p, v
                v = after
    return match


# --- hydrogens ----------------------------------------------------------------------

def _hydrogens(mol: Molecule, i: int, adjacency) -> int:
    """Return the hydrogen count of atom i: as written for a bracket atom; otherwise what is
    left by the lowest normal valence not less than the sum of its bond orders (in the Kekule
    structure)."""
    atom = mol.atoms[i]
    if atom.bracket:
        return atom.hydrogens or 0
    used = sum(mol.bonds[k].order for _, k in adjacency[i])
    valence = _lowest_at_least(NORMAL_VALENCE.get(atom.element, ()), used)
    return 0 if valence is None else valence - used


# --- aromatic rings -------------------------------------------------------------------

def aromatic_rings(mol: Molecule) -> List[int]:
    """Return the rings (indices into mol.rings) that are aromatic by Hückel's 4n+2 rule.

    The rule is applied to the Kekule structure, so a ring written C1=CC=CC=C1 counts the
    same as c1ccccc1. Each ring atom contributes electrons:

    - 1 if it has one double bond and that bond lies in the ring or in an aromatic ring
      fused to it (the shared atoms of naphthalene)
    - 2 if it has no double bond and a lone pair: N with three bonds (pyrrole), O and S
      with two (furan), a carbon or nitrogen anion
    - 0 if it has no double bond and an empty orbital (a carbocation, boron), or if it is a
      carbon with a double bond outside the rings to N, O, S, Se or Te (the C=O of pyridone)

    Any other atom (CH2, an atom with two double bonds or a triple bond, a carbon with a
    C=C leaving the ring) rules the ring out. A ring is aromatic when its electron count is
    2, 6, 10, ... A ring that fails on its own is also tested together with a ring it shares
    a bond with (azulene; indolizine, where both rings share the lone pair of the nitrogen).
    Rings of more than LARGEST_AROMATIC_RING atoms are skipped.

    A ring is tested again whenever a ring fused to it turns out to be aromatic, until
    nothing changes. So no atom of an aromatic ring has a C=C into a ring that is not."""
    adjacency = _adjacency(mol)
    in_ring: Set[frozenset] = set()
    for ring in mol.rings:
        in_ring.update(_ring_pairs(ring))
    small = [r for r, ring in enumerate(mol.rings) if len(ring) <= LARGEST_AROMATIC_RING]
    pairs = [(r, s) for r in small for s in small
             if r < s and _ring_pairs(mol.rings[r]) & _ring_pairs(mol.rings[s])]
    found: Set[int] = set()
    aromatic_pairs: Set[frozenset] = set()

    def electrons(i: int, atoms: Set[int]) -> Optional[int]:
        atom = mol.atoms[i]
        if atom.element not in AROMATIC_ELEMENTS:
            return None
        multiple = [(j, k) for j, k in adjacency[i] if mol.bonds[k].order != 1]
        bonds = len(adjacency[i]) + (atom.h_count or 0)
        if not multiple:
            if atom.element == "C":
                return {-1: 2, 1: 0}.get(atom.charge) if bonds == 3 else None
            if atom.element == "B":
                return 0 if atom.charge == 0 and bonds == 3 else None
            if bonds - atom.charge == LONE_PAIR_BONDS.get(atom.element) and atom.charge in (0, -1):
                return 2
            return None
        if len(multiple) > 1 or mol.bonds[multiple[0][1]].order != 2:
            return None
        j = multiple[0][0]
        pair = frozenset((i, j))
        if pair in in_ring:
            return 1 if j in atoms or pair in aromatic_pairs else None
        return 0 if atom.element == "C" and mol.atoms[j].element in PULLS_ELECTRONS else None

    def huckel(rings: Sequence[int]) -> bool:
        atoms = {a for r in rings for a in mol.rings[r]}
        total = 0
        for i in atoms:
            n = electrons(i, atoms)
            if n is None:
                return False
            total += n
        return total % 4 == 2

    changed = True
    while changed:
        changed = False
        for unit in [(r,) for r in small] + pairs:
            if all(r in found for r in unit) or not huckel(unit):
                continue
            found.update(unit)
            for r in unit:
                aromatic_pairs.update(_ring_pairs(mol.rings[r]))
            changed = True
    return sorted(found)


def _ring_pairs(ring: Sequence[int]) -> Set[frozenset]:
    """Return the bonds of a ring, each as the set of its two atoms."""
    return {frozenset((a, b)) for a, b in zip(ring, list(ring[1:]) + list(ring[:1]))}


# --- equivalent atoms -----------------------------------------------------------------

def _symmetry_classes(mol: Molecule, adjacency, stereo: Optional[List[int]] = None) -> List[int]:
    """Return a class number for each atom; atoms with the same number are equivalent in the
    molecule. The classes start from the element, isotope, charge, hydrogens and number of
    bonds, and are refined by the bonds to the neighbours and their classes until they stop
    splitting. `stereo` adds a label to each atom (0 for none)."""
    def bond_label(k):
        b = mol.bonds[k]
        return 5 if b.kind == "aromatic" else b.order   # 5: not one of the orders 1 to 4

    labels = [(a.atomic_number, a.isotope or 0, a.charge, a.h_count, len(adjacency[i]),
               stereo[i] if stereo else 0) for i, a in enumerate(mol.atoms)]
    ranks = _number(labels)
    while True:
        labels = [(ranks[i], tuple(sorted((bond_label(k), ranks[j]) for j, k in adjacency[i])))
                  for i in range(len(mol.atoms))]
        refined = _number(labels)
        if len(set(refined)) == len(set(ranks)):
            return refined
        ranks = refined


def _number(labels) -> List[int]:
    order = {label: n for n, label in enumerate(sorted(set(labels)))}
    return [order[label] for label in labels]


# --- stereocentres -------------------------------------------------------------------------

def _stereocentres(mol: Molecule, adjacency, ranks: List[int]) -> None:
    """Keep "@" and "@@" only where they mean something.

    An atom is a stereocentre when its four neighbours (counting the hydrogens in brackets,
    or three neighbours and a lone pair on S, Se, P, As) are all different. Two equivalent
    neighbours may still differ in the stereocentres they contain (the middle carbon of
    2,3,4-trihydroxyglutaric acid), so the classes are recomputed with the sense of each
    stereocentre as a label, until nothing changes. A ring atom whose two ring neighbours
    are equivalent (cis-1,4-dimethylcyclohexane) is kept when one of its rings also contains
    another such marked atom or a stereocentre (also across two rings: trans-decalin).
    """
    marked = [i for i, a in enumerate(mol.atoms) if a.chirality in ("@", "@@") and _can_be_centre(mol, i, adjacency)]
    centres: set = set()
    classes = ranks
    for _ in range(len(mol.atoms) + 1):
        found = {i for i in marked if _all_different(mol, i, classes)}
        if found == centres:
            break
        centres = found
        labels = [0] * len(mol.atoms)
        for i in centres:
            labels[i] = 1 + (_sense(mol, i, classes) == "@")
        classes = _symmetry_classes(mol, adjacency, labels)
    ring_like = [i for i in marked if i not in centres and _ring_pair(mol, i, adjacency, classes)]
    kept = set(centres)
    for i in ring_like:
        others = (set(ring_like) | centres) - {i}
        if any(i in ring and others & set(ring) for ring in mol.rings):
            kept.add(i)
    for i in kept:
        mol.atoms[i].stereocentre = True


def _can_be_centre(mol: Molecule, i: int, adjacency) -> bool:
    atom = mol.atoms[i]
    heavy = len(adjacency[i])
    if atom.h_count > 1:
        return False
    if heavy + atom.h_count == 4:
        return True
    return heavy + atom.h_count == 3 and atom.element in LONE_PAIR_STEREO


def _neighbour_classes(mol: Molecule, i: int, classes: List[int]) -> List[int]:
    """Return the classes of the neighbours in Atom.neighbor_order; a hydrogen in brackets
    is -1, the same as a plain [H] atom."""
    out = []
    for j in mol.atoms[i].neighbor_order:
        if j == IMPLICIT_H or _plain_hydrogen(mol.atoms[j]):
            out.append(-1)
        else:
            out.append(classes[j])
    return out


def _plain_hydrogen(atom: Atom) -> bool:
    return atom.element == "H" and not atom.isotope and atom.charge == 0 and not atom.h_count


def _all_different(mol: Molecule, i: int, classes: List[int]) -> bool:
    c = _neighbour_classes(mol, i, classes)
    return len(set(c)) == len(c)


def _sense(mol: Molecule, i: int, classes: List[int]) -> str:
    """Return the mark of the same centre with its neighbours in the order of their classes,
    so that two equivalent centres with the same sense get the same mark."""
    atom = mol.atoms[i]
    keys = _neighbour_classes(mol, i, classes)
    inversions = sum(1 for a in range(len(keys)) for b in range(a + 1, len(keys)) if keys[a] > keys[b])
    if inversions % 2 == 0:
        return atom.chirality
    return "@@" if atom.chirality == "@" else "@"


def _ring_pair(mol: Molecule, i: int, adjacency, classes: List[int]) -> bool:
    """Return True if the only equivalent neighbours of atom i are two of its ring neighbours
    (in one ring, as in 1,4-dimethylcyclohexane, or one in each of two fused rings, as in decalin)."""
    c = _neighbour_classes(mol, i, classes)
    alike = [j for j, cj in zip(mol.atoms[i].neighbor_order, c) if c.count(cj) > 1]
    if len(alike) != 2 or IMPLICIT_H in alike:
        return False
    return all(any(i in ring and j in ring for ring in mol.rings) for j in alike)


# --- E/Z ---------------------------------------------------------------------------------

def _cis_trans(mol: Molecule, adjacency, ranks: List[int]) -> None:
    """Read E/Z from / and \\ on the bonds next to a double bond.

    "/" means the bond goes up from its begin atom to its end atom, "\\" down. On each
    side of the double bond, the lowest-numbered neighbour is the reference; the double
    bond is cis when both references are up or both are down. A double bond has no E/Z
    when a side has two equivalent neighbours or no mark, or when the bond lies in a ring
    of fewer than 8 atoms.
    """
    for k, bond in enumerate(mol.bonds):
        if bond.kind != "double":
            continue
        sides = []
        for atom, partner in ((bond.begin, bond.end), (bond.end, bond.begin)):
            side = _side(mol, atom, partner, k, adjacency, ranks)
            if side is None:
                break
            sides.append(side)
        else:
            if _smallest_ring_through(mol, k, adjacency) < SMALLEST_RING_WITH_EZ:
                continue
            (ref_a, up_a), (ref_b, up_b) = sides
            bond.cis_trans = (ref_a, ref_b, "cis" if up_a == up_b else "trans")


def _side(mol: Molecule, atom: int, partner: int, double: int, adjacency, ranks) -> Optional[Tuple[int, bool]]:
    """Return (the reference neighbour, whether it is up) for one end of a double bond, or None."""
    others = [(j, k) for j, k in adjacency[atom] if k != double]
    if not others or len(others) + mol.atoms[atom].h_count > 2:
        return None
    if len(others) == 2 and ranks[others[0][0]] == ranks[others[1][0]]:
        return None
    ups = {}
    for j, k in others:
        b = mol.bonds[k]
        if b.symbol in ("/", "\\"):
            goes_up = b.symbol == "/"
            ups[j] = goes_up if b.begin == atom else not goes_up
    if not ups or (len(ups) == 2 and len(set(ups.values())) == 1):
        return None   # no mark, or two marks that put both neighbours on the same side
    ref = min(j for j, _ in others)
    if ref in ups:
        return ref, ups[ref]
    (other, up), = ups.items()
    return ref, not up


def _smallest_ring_through(mol: Molecule, k: int, adjacency) -> float:
    """Return the size of the smallest ring through bond k (infinity if it is in no ring)."""
    start, goal = mol.bonds[k].begin, mol.bonds[k].end
    dist = {start: 0}
    queue = deque([start])
    while queue:
        x = queue.popleft()
        for y, kk in adjacency[x]:
            if kk == k or y in dist:
                continue
            dist[y] = dist[x] + 1
            if y == goal:
                return dist[y] + 1
            queue.append(y)
    return float("inf")
