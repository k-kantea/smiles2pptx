"""Overlap removal after the first placement, by rotating parts of the molecule.

The first placement (place.py) works outwards from one atom and never looks back, so a long
chain can curl back onto itself and two branches can meet. This step counts overlaps the
way check.py does (atoms as points, and the bounding boxes of atom labels), then repeatedly
makes the move that lowers the count the most, until none are left or no move helps:

- flip: rotate one side of a single bond that is not in a ring by 180 degrees about the
  bond. In the plane this reflects that side in the line of the bond. Bond lengths, angles,
  rings and the E/Z of double bonds are kept
- swap: exchange two branches of an atom with three or more neighbours, by rotating each
  about the atom into the other's position. Not done at an atom of an E/Z double bond,
  where it would turn E into Z
- bend: rotate a branch about its atom into the middle of another empty angle there. This
  changes a bond angle (a substituent of a bridged ring that points into the rings moves
  outside), so it is chosen only when it clears more overlaps than any flip or swap. Not
  done at an atom of an E/Z double bond, or where two bonds are in a line

Among moves that clear the same number of overlaps, the one that leaves non-bonded atoms
less crowded (closer than a bond length) is chosen. Only moves that can separate an
overlap are tried: those that move some of its atoms but not all of them.
"""
from __future__ import annotations

import math
from collections import deque
from typing import Dict, FrozenSet, Iterable, List, Optional, Sequence, Set, Tuple

from ..check import CLASH, ON_BOND, parts_meet, segment_distance, segment_in_parts, segments_cross
from ..chem.molecule import Molecule
from ..depict.labels import atom_label, label_parts
from ..fonts import get_font
from ..style import Style
from .geometry import TAU, Point, add, direction, gaps, rotate, sub

MAX_MOVES = 100
CLOSE = 1.0   # non-bonded atoms closer than this (bond lengths) add to the tie-break


def untangle(mol: Molecule, points: Dict[int, Point], style: Style) -> Dict[int, Point]:
    """Return the points of one fragment, with overlaps cleared wherever a move can clear them."""
    scorer = _Scorer(mol, sorted(points), style)
    pos = dict(points)
    for _ in range(MAX_MOVES):
        hard, _, items = scorer.cost(pos, None, collect=True)
        if hard == 0:
            break
        best = None
        for moved, new, bends in scorer.moves(pos, items):
            before = scorer.cost(pos, moved)
            trial = dict(pos)
            trial.update(new)
            after = scorer.cost(trial, moved)
            # Fewer overlaps first, then a move that keeps the bond angles, then less crowding.
            gain = (after[0] - before[0], bends, after[1] - before[1])
            if best is None or gain < best[0]:
                best = (gain, new)
        if best is None:
            break
        fewer, bends, closer = best[0]
        if not (fewer < 0 or (fewer == 0 and not bends and closer < -1e-6)):
            break
        pos.update(best[1])
    return pos


class _Scorer:
    def __init__(self, mol: Molecule, atoms: List[int], style: Style):
        self.mol = mol
        self.atoms = atoms
        here = set(atoms)
        self.bonds: List[Tuple[int, int]] = [(b.begin, b.end) for b in mol.bonds if b.begin in here]
        self.orders = [b.order for b in mol.bonds if b.begin in here]
        self.bonded: Set[FrozenSet[int]] = {frozenset(b) for b in self.bonds}
        self.adjacency: Dict[int, List[int]] = {a: [] for a in atoms}
        for i, j in self.bonds:
            self.adjacency[i].append(j)
            self.adjacency[j].append(i)
        self.font = get_font(style.font)
        self.size = style.font_size / style.bond_length
        self.script = style.script_scale
        dummy = {a: (0.0, 0.0) for a in atoms}
        self.labelled = [a for a in atoms if atom_label(mol, a, dummy) is not None]
        self.ez_atoms = {a for b in mol.bonds if b.cis_trans is not None and b.begin in here
                         for a in (b.begin, b.end)}
        self.straight = set()
        for a in atoms:
            orders = [self.orders[k] for k, bond in enumerate(self.bonds) if a in bond]
            if len(orders) == 2 and (3 in orders or orders.count(2) == 2):
                self.straight.add(a)
        self.sides = self._acyclic_sides()

    # --- counting overlaps ------------------------------------------------------------

    def boxes(self, pos: Dict[int, Point], which: Iterable[int]):
        """Return {atom: [boxes]}: its label, plus its H if that is written above or below the symbol."""
        out = {}
        for a in which:
            label = atom_label(self.mol, a, pos)
            if label is not None:
                out[a] = label_parts(label, pos[a], self.font, self.size, self.script)
        return out

    def cost(self, pos: Dict[int, Point], moved: Optional[Set[int]], collect: bool = False):
        """Return (overlaps, closeness, items), counting only what involves an atom in `moved`
        (everything if it is None). With `collect`, `items` lists the atoms of each overlap;
        otherwise it is empty."""
        atoms, bonds = self.atoms, self.bonds
        every = moved is None
        M = set(atoms) if every else moved
        hard, soft, items = 0, 0.0, []
        # Atom against atom.
        for a in M:
            pa = pos[a]
            for b in atoms:
                if b == a or (b in M and b < a) or frozenset((a, b)) in self.bonded:
                    continue
                d = math.dist(pa, pos[b])
                if d < CLOSE:
                    soft += (CLOSE - d) ** 2
                    if d < CLASH:
                        hard += 1
                        if collect:
                            items.append((a, b))
        # Atom on a bond, and bond crossing bond.
        touched = [k for k, (i, j) in enumerate(bonds) if i in M or j in M]
        touched_set = set(touched)
        for k in touched:
            i, j = bonds[k]
            for a in atoms:
                if a not in (i, j) and self._on_bond(pos, a, i, j):
                    hard += 1
                    if collect:
                        items.append((a, i, j))
        for a in M:
            for k, (i, j) in enumerate(bonds):
                if k not in touched_set and a not in (i, j) and self._on_bond(pos, a, i, j):
                    hard += 1
                    if collect:
                        items.append((a, i, j))
        for k1 in touched:
            i1, j1 = bonds[k1]
            for k2, (i2, j2) in enumerate(bonds):
                if k2 == k1 or (k2 in touched_set and k2 < k1) or len({i1, j1, i2, j2}) < 4:
                    continue
                if segments_cross(pos[i1], pos[j1], pos[i2], pos[j2]):
                    hard += 1
                    if collect:
                        items.append((i1, j1, i2, j2))
        # Labels. A label can change when its atom or a neighbour has moved.
        if self.labelled:
            changed = M | {n for a in M for n in self.adjacency[a]}
            boxes = self.boxes(pos, self.labelled)
            near = [a for a in self.labelled if a in changed]
            near_set = set(near)
            for a in near:
                for i, j in bonds:
                    if a not in (i, j) and segment_in_parts(pos[i], pos[j], boxes[a]):
                        hard += 1
                        if collect:
                            items.append((a, i, j))
            for k in touched:
                i, j = bonds[k]
                for a in self.labelled:
                    if a not in near_set and a not in (i, j) and segment_in_parts(pos[i], pos[j], boxes[a]):
                        hard += 1
                        if collect:
                            items.append((a, i, j))
            for a in near:
                for b in self.labelled:
                    if b == a or (b in near_set and b < a):
                        continue
                    if parts_meet(boxes[a], boxes[b]):
                        hard += 1
                        if collect:
                            items.append((a, b))
        return hard, soft, items

    @staticmethod
    def _on_bond(pos, a, i, j) -> bool:
        p = pos[a]
        return (segment_distance(p, pos[i], pos[j]) < ON_BOND
                and math.dist(p, pos[i]) >= CLASH and math.dist(p, pos[j]) >= CLASH)

    # --- moves ------------------------------------------------------------------------

    def _acyclic_sides(self) -> Dict[Tuple[int, int], FrozenSet[int]]:
        """Return {(u, v): the atoms on v's side} for each bond u-v that is not in a ring."""
        sides = {}
        for i, j in self.bonds:
            side_j = self._reach(j, i)
            if i in side_j:
                continue
            side_i = frozenset(a for a in self.atoms if a not in side_j)
            sides[(i, j)] = side_j
            sides[(j, i)] = side_i
        return sides

    def _reach(self, start: int, without: int) -> FrozenSet[int]:
        """Return the atoms reachable from `start` without crossing the bond start-without."""
        seen = {start}
        queue = deque([start])
        while queue:
            a = queue.popleft()
            for b in self.adjacency[a]:
                if b in seen or (a == start and b == without):
                    continue
                seen.add(b)
                queue.append(b)
        return frozenset(seen)

    def moves(self, pos: Dict[int, Point], items: Sequence[Tuple[int, ...]]):
        """Yield (moved atoms, their new points, bends) for each move that can separate an
        overlap; `bends` is 1 for a bend and 0 otherwise."""
        involved = [set(it) for it in items]

        def parts(moved: Set[int]) -> bool:
            return any(0 < len(s & moved) < len(s) for s in involved)

        for k, (u, v) in enumerate(self.bonds):
            if self.orders[k] != 1 or (u, v) not in self.sides:
                continue
            # Reflect the smaller side; its atom on the bond stays where it is.
            far, near = self.sides[(u, v)], self.sides[(v, u)]
            side, pivot = (far, v) if len(far) <= len(near) else (near, u)
            moved = set(side) - {pivot}
            if not moved or not parts(moved):
                continue
            a, b = pos[u], pos[v]
            new = {x: _mirror(pos[x], a, b) for x in moved}
            if all(math.dist(new[x], pos[x]) < 1e-9 for x in moved):
                continue
            yield moved, new, 0
        for a in self.atoms:
            nbrs = self.adjacency[a]
            if len(nbrs) < 3 or a in self.ez_atoms:
                continue
            free = [b for b in nbrs if (a, b) in self.sides]
            for x in range(len(free)):
                for y in range(x + 1, len(free)):
                    b, c = free[x], free[y]
                    side_b, side_c = self.sides[(a, b)], self.sides[(a, c)]
                    moved = set(side_b) | set(side_c)
                    if not parts(moved):
                        continue
                    turn = direction(pos[a], pos[c]) - direction(pos[a], pos[b])
                    new = {}
                    for atom in side_b:
                        new[atom] = add(pos[a], rotate(sub(pos[atom], pos[a]), turn))
                    for atom in side_c:
                        new[atom] = add(pos[a], rotate(sub(pos[atom], pos[a]), -turn))
                    yield moved, new, 0
        for a in self.atoms:
            if a in self.ez_atoms or a in self.straight:
                continue
            for b in self.adjacency[a]:
                if (a, b) not in self.sides:
                    continue
                moved = set(self.sides[(a, b)])
                if not parts(moved):
                    continue
                now = direction(pos[a], pos[b])
                others = [direction(pos[a], pos[c]) for c in self.adjacency[a] if c != b]
                for start, size in gaps(others):
                    turn = math.remainder(start + size / 2 - now, TAU)
                    if abs(turn) < 1e-6:
                        continue
                    new = {atom: add(pos[a], rotate(sub(pos[atom], pos[a]), turn)) for atom in moved}
                    yield moved, new, 1


def _mirror(p: Point, a: Point, b: Point) -> Point:
    """Return p reflected in the line through a and b."""
    dx, dy = b[0] - a[0], b[1] - a[1]
    length2 = dx * dx + dy * dy
    t = ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / length2
    fx, fy = a[0] + t * dx, a[1] + t * dy
    return (2 * fx - p[0], 2 * fy - p[1])
