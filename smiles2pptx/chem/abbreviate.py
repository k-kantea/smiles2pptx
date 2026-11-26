"""Substituent abbreviations: a group at the end of a bond collapsed into one labelled atom.

`abbreviate(mol, groups)` takes a perceived molecule and the groups to look for. Each group
is (a SMILES starting with "*", which stands for the atom the group is bonded to; the label
when the bond comes from the left; the label when it comes from the right), as in the table
in style.py:

    ("*C(=O)O", "CO2H", "HO2C")

It returns a copy of the molecule in which each group found is one atom (the one that was
bonded to "*") with both labels in Atom.label. The other atoms of the group are removed and
the atoms are renumbered. Layout and drawing then treat it like any other labelled atom.
A molecule without such a group is returned unchanged.

A group matches where the molecule has the same atoms (element, charge, isotope,
hydrogens), the same bonds (a bond in an aromatic ring matches one, however the ring is
written) and no bond leaving the group other than the one to "*". So CO2H does not match
an ester or a carboxylate. Other rules:

- a group bonded to an atom that already stands for a group is left as it is: in biphenyl
  one ring becomes Ph, not both
- of two groups that share atoms, the larger wins (OAc over acetyl)
- a group with a stereocentre, or with a double bond that has E/Z, is left as it is
"""
from __future__ import annotations

import copy
from functools import lru_cache
from typing import Dict, List, Sequence, Set, Tuple

from .molecule import IMPLICIT_H, Atom, Molecule
from .perceive import perceive
from .smiles import read_smiles

Group = Tuple[str, str, str]   # SMILES with "*", label for a bond from the left, label for one from the right


def abbreviate(mol: Molecule, groups: Sequence[Group]) -> Molecule:
    """Return the molecule with each of `groups` collapsed into one labelled atom (see the module)."""
    if not groups:
        return mol
    adjacency = _adjacency(mol)
    aromatic = mol.aromatic_bonds()
    ez = {a for b in mol.bonds if b.cis_trans is not None for a in (b.begin, b.end)}
    found = []
    for g, (smiles, left, right) in enumerate(groups):
        pattern = _pattern(smiles)
        for anchor, atoms in _places(mol, adjacency, aromatic, pattern):
            if not ez & set(atoms):
                found.append((-len(atoms), g, atoms[0], anchor, atoms, (left, right)))
    found.sort(key=lambda f: f[:4])
    labels: Dict[int, Tuple[str, str]] = {}
    removed: Set[int] = set()
    for _, _, attach, anchor, atoms, texts in found:
        taken = removed | set(labels)
        if taken & set(atoms) or anchor in taken:
            continue
        labels[attach] = texts
        removed.update(atoms[1:])
    if not labels:
        return mol
    return _without(mol, removed, labels)


@lru_cache(maxsize=None)
def _pattern(smiles: str) -> Molecule:
    """Return the group as a perceived molecule; atom 0 is "*"."""
    mol = perceive(read_smiles(smiles))
    stars = [i for i, a in enumerate(mol.atoms) if a.element == "*"]
    if stars != [0] or len(mol.atoms) < 2 or len(mol.neighbors(0)) != 1 or len(mol.fragments) != 1:
        raise ValueError(f"invalid abbreviation {smiles!r}: expected a single group as SMILES starting with '*' "
                         "(the attachment point) and no other '*'")
    return mol


def _adjacency(mol: Molecule) -> List[List[Tuple[int, int]]]:
    """Return a list of (neighbour, bond index) for each atom."""
    adjacency: List[List[Tuple[int, int]]] = [[] for _ in mol.atoms]
    for k, b in enumerate(mol.bonds):
        adjacency[b.begin].append((b.end, k))
        adjacency[b.end].append((b.begin, k))
    return adjacency


def _places(mol: Molecule, adjacency, aromatic: Set[int], pattern: Molecule):
    """Yield (the atom the group is bonded to, the atoms of the group with the bonded one
    first) for each match of the group. Matching a ring in both directions counts once."""
    p_adjacency = _adjacency(pattern)
    p_aromatic = pattern.aromatic_bonds()
    attach, star_bond = p_adjacency[0][0]
    order = [attach]   # atoms of the pattern, each after an atom it is bonded to
    for p in order:
        order += [q for q, _ in p_adjacency[p] if q != 0 and q not in order]

    def same_atom(p: int, i: int) -> bool:
        a, b = pattern.atoms[p], mol.atoms[i]
        return (a.element == b.element and a.charge == b.charge and a.isotope == b.isotope
                and a.h_count == b.h_count and len(p_adjacency[p]) == len(adjacency[i])
                and not b.stereocentre and b.label is None)

    def same_bond(pk: int, k: int) -> bool:
        if k < 0 or (pk in p_aromatic) != (k in aromatic):
            return False
        return pk in p_aromatic or pattern.bonds[pk].order == mol.bonds[k].order

    def bond_index(i: int, j: int) -> int:
        return next((k for n, k in adjacency[i] if n == j), -1)

    def extend(mapping: Dict[int, int], n: int) -> bool:
        if n == len(order):
            return True
        p = order[n]
        linked = [(q, pk) for q, pk in p_adjacency[p] if q in mapping]
        first = mapping[linked[0][0]]
        used = set(mapping.values())
        for i, _ in adjacency[first]:
            if i in used or not same_atom(p, i):
                continue
            if all(same_bond(pk, bond_index(i, mapping[q])) for q, pk in linked):
                mapping[p] = i
                if extend(mapping, n + 1):
                    return True
                del mapping[p]
        return False

    for i in range(len(mol.atoms)):
        if not same_atom(attach, i):
            continue
        for anchor, k in adjacency[i]:
            if not same_bond(star_bond, k):
                continue
            mapping = {0: anchor, attach: i}
            if extend(mapping, 1):
                yield anchor, [mapping[p] for p in order]


def _without(mol: Molecule, removed: Set[int], labels: Dict[int, Tuple[str, str]]) -> Molecule:
    """Return a copy without the `removed` atoms, renumbered, with the labels in `labels`."""
    new = {old: n for n, old in enumerate(i for i in range(len(mol.atoms)) if i not in removed)}
    atoms: List[Atom] = []
    for old in new:
        atom = copy.copy(mol.atoms[old])
        atom.neighbor_order = [j if j == IMPLICIT_H else new[j] for j in atom.neighbor_order
                               if j == IMPLICIT_H or j in new]
        if old in labels:
            atom.label = labels[old]
            atom.h_count = 0
        atoms.append(atom)
    bonds = []
    for b in mol.bonds:
        if b.begin in new and b.end in new:
            bond = copy.copy(b)
            bond.begin, bond.end = new[b.begin], new[b.end]
            if b.cis_trans is not None:
                bond.cis_trans = (new[b.cis_trans[0]], new[b.cis_trans[1]], b.cis_trans[2])
            bonds.append(bond)
    fragments = [[new[i] for i in part if i in new] for part in mol.fragments]
    kept = [r for r, ring in enumerate(mol.rings) if all(i in new for i in ring)]
    rings = [[new[i] for i in mol.rings[r]] for r in kept]
    aromatic_rings = [kept.index(r) for r in mol.aromatic_rings if r in kept]
    return Molecule(mol.smiles, atoms, bonds, fragments, rings, aromatic_rings)
