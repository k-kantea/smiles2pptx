"""Placement of the atoms of a perceived molecule in the plane.

The method follows SmilesDrawer (Probst and Reymond, 2018), reimplemented from its
description. Ring systems are placed first, as regular polygons that share bonds (fused)
or atoms (spiro) (rings.py). From there, placement spreads outwards one bond at a time:

- a chain zigzags, with the bond angle of the style
- the bonds from a branching atom share the space around it evenly; the longest
  branch continues the zigzag
- an atom with a triple bond or two double bonds keeps its two bonds in a straight line
- a substituent of a ring points away from the ring
- a ring system reached along the way is rotated so that the bond to it points away from it
- the / and \\ of a double bond decide which side its substituents go on

Each fragment of the SMILES (split by ".") is placed on its own, cleared of overlaps
(untangle.py) and rotated in steps of 30 degrees to make it as flat as possible, and the
fragments are set side by side, left to right, in the order they are written.
"""
from __future__ import annotations

import math
from collections import deque
from typing import Dict, List, Optional, Sequence, Tuple

from ..chem.molecule import Molecule
from ..style import Style
from .geometry import TAU, Point, add, at, direction, rotate, side, sub
from .rings import outward_gap, place_ring_system, ring_systems
from .lengthen import lengthen
from .untangle import untangle

# Maximum number of times a fragment is rotated flat and cleared of overlaps again.
FLAT_ROUNDS = 3

# A slot for a new atom: the direction of the bond to it, and the way the zigzag
# turns at that atom (1 or -1; None if nothing sets it).
Slot = Tuple[float, Optional[int]]


def layout(mol: Molecule, style: Optional[Style] = None) -> List[Point]:
    """Return one point per atom, in the order of mol.atoms (bond length 1, y pointing up).

    The molecule must come from perceive() (chem/perceive.py), which finds its rings."""
    style = style or Style()
    if any(atom.h_count is None for atom in mol.atoms):
        raise ValueError("call perceive() on the molecule before layout()")
    placer = _Placer(mol, math.radians(style.bond_angle))
    parts = []
    for part in mol.fragments:
        points = untangle(mol, placer.place_part(part), style)
        # Labels are boxes that do not rotate with the molecule: clear overlaps again after
        # rotating, and rotate flat again if that changed anything.
        for _ in range(FLAT_ROUNDS):
            flat = _flattest(points)
            points = untangle(mol, flat, style)
            if points == flat:
                break
        parts.append(lengthen(mol, points, style))
    return _side_by_side(len(mol.atoms), parts, style.fragment_gap)


class _Placer:
    def __init__(self, mol: Molecule, bond_angle: float):
        self.mol = mol
        self.bond_angle = bond_angle
        self.adjacency: List[List[Tuple[int, int]]] = [[] for _ in mol.atoms]
        for k, b in enumerate(mol.bonds):
            self.adjacency[b.begin].append((b.end, k))
            self.adjacency[b.end].append((b.begin, k))
        self.systems = ring_systems(mol.rings)
        self.system_of: Dict[int, int] = {}
        for s, ring_ids in enumerate(self.systems):
            for i in ring_ids:
                for atom in mol.rings[i]:
                    self.system_of[atom] = s
        self.pos: Dict[int, Point] = {}
        self.turn: Dict[int, Optional[int]] = {}

    def neighbours(self, atom: int) -> List[int]:
        return [b for b, _ in self.adjacency[atom]]

    # --- one fragment of the SMILES --------------------------------------------------

    def place_part(self, part: Sequence[int]) -> Dict[int, Point]:
        self.pos, self.turn = {}, {}
        systems_here = sorted({self.system_of[a] for a in part if a in self.system_of})
        queue: deque = deque()
        if systems_here:
            # Start from the largest ring system.
            start = max(systems_here, key=lambda s: (len(self._system_atoms(s)), -s))
            local = place_ring_system(self.mol.rings, self.systems[start], self._outside(start))
            self.pos.update(local)
            queue.extend(sorted(local))
        else:
            # Start from one end of the longest chain.
            first = self._chain_start(part)
            self.pos[first] = (0.0, 0.0)
            queue.append(first)
        while queue:
            queue.extend(self._expand(queue.popleft()))
        return {a: self.pos[a] for a in part}

    def _system_atoms(self, s: int) -> set:
        return {a for i in self.systems[s] for a in self.mol.rings[i]}

    def _outside(self, s: int) -> List[Tuple[int, int]]:
        """Return the bonds from atoms of ring system s to atoms outside it."""
        inside = self._system_atoms(s)
        return [(a, b) for a in sorted(inside) for b in self.neighbours(a) if b not in inside]

    def _chain_start(self, part: Sequence[int]) -> int:
        def farthest(source: int) -> int:
            dist = {source: 0}
            queue = deque([source])
            while queue:
                a = queue.popleft()
                for b in self.neighbours(a):
                    if b not in dist:
                        dist[b] = dist[a] + 1
                        queue.append(b)
            return max(dist, key=lambda a: (dist[a], -a))

        end = farthest(part[0])
        return min(end, farthest(end))

    # --- placing the neighbours of one atom --------------------------------------------

    def _expand(self, atom: int) -> List[int]:
        """Place the neighbours of `atom` that are not yet placed, and return the atoms placed."""
        children = [b for b in self.neighbours(atom) if b not in self.pos]
        if not children:
            return []
        placed = [b for b in self.neighbours(atom) if b in self.pos]
        slots, prefer = self._slots(atom, placed, len(children))
        order = sorted(children, key=lambda c: -self._branch_size(c, atom))
        sides = self._ez_sides(atom, children)
        free = list(prefer)
        chosen: Dict[int, int] = {}
        for child in sorted(order, key=lambda c: c not in sides):
            pick = free[0]
            if child in sides:
                partner, wanted = sides[child]
                for i in free:
                    if side(self.pos[partner], self.pos[atom], at(self.pos[atom], slots[i][0])) == wanted:
                        pick = i
                        break
            free.remove(pick)
            chosen[child] = pick
        new: List[int] = []
        for child in children:
            angle, turn = slots[chosen[child]]
            self.pos[child] = at(self.pos[atom], angle)
            self.turn[child] = turn
            new.append(child)
        for child in children:
            if child in self.system_of:
                new.extend(self._attach_system(child, atom))
        return new

    def _slots(self, atom: int, placed: List[int], k: int) -> Tuple[List[Slot], List[int]]:
        """Return the slots for k new bonds from `atom`, and the order in which to fill them."""
        here = self.pos[atom]
        if not placed:
            # First atom of a chain: the first bond slopes down so that the zigzag runs flat.
            if k == 1:
                return [(-(math.pi - self.bond_angle) / 2, -1)], [0]
            return [(i * TAU / k, None) for i in range(k)], list(range(k))
        if len(placed) == 1 and atom not in self.system_of:
            back = direction(here, self.pos[placed[0]])
            if k == 1 and self._straight(atom):
                return [(back + math.pi, self.turn.get(atom))], [0]
            if k <= 2:
                spread = self.bond_angle if k == 1 else TAU / 3
                slots = [(back + s * spread, -s) for s in (1, -1)]
                first = self._turn_here(atom, back, spread)
                return slots, ([0, 1] if first == 1 else [1, 0])
            slots = [(back + (i + 1) * TAU / (k + 1), None) for i in range(k)]
            straight = sorted(range(k), key=lambda i: abs(math.remainder(slots[i][0] - back - math.pi, TAU)))
            return slots, straight
        start, size = outward_gap(atom, self.pos, placed, [self.mol.rings[i] for i in self._rings_of(atom)])
        return [(start + size * (i + 1) / (k + 1), None) for i in range(k)], list(range(k))

    def _turn_here(self, atom: int, back: float, spread: float) -> int:
        """Return the way the zigzag turns at `atom`: as set when it was placed, or else
        towards the less crowded side."""
        turn = self.turn.get(atom)
        if turn is not None:
            return turn
        here = self.pos[atom]
        crowd = {s: self._crowding(at(here, back + s * spread)) for s in (1, -1)}
        return 1 if crowd[1] <= crowd[-1] + 1e-12 else -1

    def _crowding(self, point: Point) -> float:
        return sum(1 / max((p[0] - point[0]) ** 2 + (p[1] - point[1]) ** 2, 1e-9)
                   for p in self.pos.values())

    def _straight(self, atom: int) -> bool:
        """Return True if the two bonds of `atom` must be in a line: a triple bond, or two
        double bonds (allenes, CO2)."""
        orders = [self.mol.bonds[k].order for _, k in self.adjacency[atom]]
        return len(orders) == 2 and (3 in orders or orders.count(2) == 2)

    def _rings_of(self, atom: int) -> List[int]:
        return [i for i, ring in enumerate(self.mol.rings) if atom in ring]

    def _branch_size(self, child: int, atom: int) -> int:
        """Return the number of atoms reachable from `child` without going back through `atom`."""
        seen = {atom, child}
        queue = deque([child])
        while queue:
            a = queue.popleft()
            for b in self.neighbours(a):
                if b not in seen:
                    seen.add(b)
                    queue.append(b)
        return len(seen) - 1

    def _ez_sides(self, atom: int, children: List[int]) -> Dict[int, Tuple[int, int]]:
        """Return {child: (partner, side)} for the E/Z double bonds from `atom` to a placed
        partner: `side` is the side of the line partner -> atom that the reference neighbour of
        `atom` must go on, so that it is cis or trans to the reference neighbour of the partner."""
        wanted: Dict[int, Tuple[int, int]] = {}
        for partner, k in self.adjacency[atom]:
            bond = self.mol.bonds[k]
            if bond.cis_trans is None or partner not in self.pos:
                continue
            ref_begin, ref_end, relation = bond.cis_trans
            mine, theirs = (ref_begin, ref_end) if bond.begin == atom else (ref_end, ref_begin)
            if mine not in children or theirs not in self.pos:
                continue
            their_side = side(self.pos[partner], self.pos[atom], self.pos[theirs])
            wanted[mine] = (partner, their_side if relation == "cis" else -their_side)
        return wanted

    def _attach_system(self, first: int, parent: int) -> List[int]:
        """Place the ring system of `first` (already placed, bonded to `parent`), rotated so the
        bond to `parent` points away from the rings, and return the other atoms of the system."""
        ring_ids = self.systems[self.system_of[first]]
        local = place_ring_system(self.mol.rings, ring_ids, self._outside(self.system_of[first]))
        inside_ring = [b for b in self.neighbours(first) if b in local]
        outside = [b for b in self.neighbours(first) if b not in local]
        start, size = outward_gap(first, local, inside_ring, [self.mol.rings[i] for i in ring_ids])
        local_angle = start + size / (len(outside) + 1)
        turn = direction(self.pos[first], self.pos[parent]) - local_angle
        new = []
        for atom in sorted(local):
            if atom != first:
                self.pos[atom] = add(self.pos[first], rotate(sub(local[atom], local[first]), turn))
                self.turn[atom] = None
                new.append(atom)
        return new


def _flattest(points: Dict[int, Point]) -> Dict[int, Point]:
    """Return the fragment rotated by the multiple of 30 degrees that gives it the smallest
    height (on a tie, the first such rotation counting up from 0 degrees)."""
    best, best_height = points, None
    for i in range(12):
        turned = {a: rotate(p, i * math.pi / 6) for a, p in points.items()}
        ys = [p[1] for p in turned.values()]
        height = max(ys) - min(ys)
        if best_height is None or height < best_height - 1e-6:
            best, best_height = turned, height
    return best


def _side_by_side(n: int, parts: List[Dict[int, Point]], gap: float) -> List[Point]:
    """Return the fragments placed left to right, `gap` apart and centred on one line, with
    the whole centred on (0, 0)."""
    coords: List[Point] = [(0.0, 0.0)] * n
    x = 0.0
    for points in parts:
        xs = [p[0] for p in points.values()]
        ys = [p[1] for p in points.values()]
        dx, dy = x - min(xs), -(min(ys) + max(ys)) / 2
        for a, p in points.items():
            coords[a] = (p[0] + dx, p[1] + dy)
        x += max(xs) - min(xs) + gap
    if n:
        xs = [p[0] for p in coords]
        ys = [p[1] for p in coords]
        cx, cy = (min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2
        coords = [(p[0] - cx, p[1] - cy) for p in coords]
    return coords
