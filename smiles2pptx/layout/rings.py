"""Ring systems: rings that share atoms, placed together in their own coordinates.

A ring system is placed one ring at a time, with every bond of length 1. The first ring
is a regular polygon. Each following ring is the one that shares the most atoms with the
rings already placed:

- a ring that shares one atom (spiro) is a regular polygon on the far side of that atom
- a ring that shares two or more atoms has its other atoms placed on a circular arc
  between the two ends of the shared part, on the side away from the atoms already placed.
  When it shares one bond (fused rings), the ring comes out as a regular polygon.

When the arcs put atoms of the system on top of each other (bridged rings and cages), the
system is placed with springs instead (springs.py), and the placement with the fewest
overlaps, counting the first atom of each substituent, is kept.

The system is then rotated to fit where it joins the rest of the molecule (place.py).
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Set, Tuple

from .geometry import (TAU, Point, add, at, centroid, direction, distance, dot, inside,
                       largest_gap, mul, sub)
from .springs import PATIENCE, count_overlaps, spring_batches, worst_bond


def ring_systems(rings: Sequence[Sequence[int]]) -> List[List[int]]:
    """Return the rings of each ring system (rings joined by shared atoms) as indices into
    `rings`, with the systems in the order of their first ring."""
    parent = list(range(len(rings)))

    def root(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    owner: Dict[int, int] = {}
    for i, ring in enumerate(rings):
        for atom in ring:
            if atom in owner:
                parent[root(i)] = root(owner[atom])
            else:
                owner[atom] = i
    systems: Dict[int, List[int]] = {}
    for i in range(len(rings)):
        systems.setdefault(root(i), []).append(i)
    return sorted(systems.values(), key=lambda s: s[0])


def polygon_radius(n: int) -> float:
    """Return the distance from the centre to a vertex of a regular n-gon with sides of length 1."""
    return 1 / (2 * math.sin(math.pi / n))


def outward_gap(atom: int, pos: Dict[int, Point], neighbours: Sequence[int],
                rings: Sequence[Sequence[int]]):
    """Return the largest empty angle (start, size) around a placed atom, between the
    directions to its placed neighbours, skipping angles that open into a ring of the atom."""
    here = pos[atom]
    polygons = [[pos[a] for a in ring] for ring in rings
                if atom in ring and all(a in pos for a in ring)]

    def outside(start: float, size: float) -> bool:
        probe = at(here, start + size / 2, 0.5)
        return not any(inside(probe, polygon) for polygon in polygons)

    return largest_gap([direction(here, pos[b]) for b in neighbours if b in pos], outside)


def place_ring_system(rings: Sequence[Sequence[int]], system: Sequence[int],
                      outside: Sequence[Tuple[int, int]] = ()) -> Dict[int, Point]:
    """Return coordinates for every atom of a ring system, in the system's own frame.
    `outside` holds the bonds from atoms of the system to atoms outside it; the springs
    leave room for them."""
    pos: Dict[int, Point] = {}
    placed: List[int] = []
    left = list(system)
    while left:
        ring_id = _next_ring(rings, left, pos)
        left.remove(ring_id)
        _place_ring(rings, ring_id, placed, pos)
        placed.append(ring_id)
    bonds = ring_bonds(rings, system)
    if count_overlaps(pos, bonds) == 0:
        return pos
    # Bridged rings or a cage: the arcs put atoms on top of each other. Try springs, and
    # keep the placement with the fewest overlaps once the substituents are added.
    ring_list = [rings[i] for i in system]
    best = (_score(pos, bonds, outside, ring_list), pos)
    stale = 0
    for batch in spring_batches(pos, bonds, outside):
        stale += 1
        for candidate in batch:
            score = _score(candidate, bonds, outside, ring_list)
            if score[0] < best[0][0] or (score[0] == best[0][0] and score[1] < best[0][1] - 1e-9):
                best = (score, candidate)
                stale = 0
        if best[0][0] == 0 or stale >= PATIENCE:
            break
    return best[1]


def _score(pos: Dict[int, Point], bonds, outside, ring_list) -> Tuple[int, float]:
    """Return (overlaps, largest bond-length error) for a placement of a ring system, with the
    first atom of each substituent where place.py will put it: spread over the largest empty
    angle that does not open into a ring."""
    every = dict(pos)
    neighbours: Dict[int, List[int]] = {}
    for i, j in bonds:
        neighbours.setdefault(i, []).append(j)
        neighbours.setdefault(j, []).append(i)
    leaving: Dict[int, List[int]] = {}
    for a, b in outside:
        leaving.setdefault(a, []).append(b)
    for a, subs in leaving.items():
        start, size = outward_gap(a, pos, neighbours[a], ring_list)
        for k, b in enumerate(subs):
            every[b] = at(pos[a], start + size * (k + 1) / (len(subs) + 1))
    return count_overlaps(every, list(bonds) + list(outside)), worst_bond(pos, bonds)


def ring_bonds(rings: Sequence[Sequence[int]], system: Sequence[int]) -> List[Tuple[int, int]]:
    """Return the bonds of the rings of a system, each listed once."""
    seen = set()
    bonds = []
    for i in system:
        ring = rings[i]
        for k, a in enumerate(ring):
            b = ring[(k + 1) % len(ring)]
            key = frozenset((a, b))
            if key not in seen:
                seen.add(key)
                bonds.append((a, b))
    return bonds


def _next_ring(rings, left: List[int], pos: Dict[int, Point]) -> int:
    """Return the next ring to place: first the ring fused to the most other rings, then the
    ring that shares the most atoms with those placed (two or more shared before one)."""
    if not pos:
        def fused(i: int) -> int:
            return sum(1 for j in left if j != i and len(set(rings[i]) & set(rings[j])) >= 2)
        return max(left, key=lambda i: (fused(i), -left.index(i)))
    shared = {i: sum(1 for a in rings[i] if a in pos) for i in left}
    return max(left, key=lambda i: (shared[i] >= 2, shared[i], -left.index(i)))


def _place_ring(rings, ring_id: int, placed: List[int], pos: Dict[int, Point]) -> None:
    ring = list(rings[ring_id])
    n = len(ring)
    known = [i for i, a in enumerate(ring) if a in pos]
    if not known:
        # A regular polygon with its first atom at the top.
        r = polygon_radius(n)
        for k, atom in enumerate(ring):
            pos[atom] = at((0.0, 0.0), math.pi / 2 + k * TAU / n, r)
        return
    if len(known) == n:
        return
    if len(known) == 1:
        _place_spiro(rings, ring, known[0], placed, pos)
        return
    runs = _runs(ring, known)
    for u_i, run, v_i in runs:
        u, v = ring[u_i], ring[v_i]
        away_from = _reference(rings, ring, u_i, v_i, runs, placed, pos)
        _arc(pos, pos[u], pos[v], [ring[i] for i in run], away_from)


def _place_spiro(rings, ring: List[int], j: int, placed: List[int], pos: Dict[int, Point]) -> None:
    """Place a spiro ring as a regular polygon through the shared atom, on the far side of the
    rings that contain it."""
    n = len(ring)
    s = ring[j]
    neighbours = set()
    for i in placed:
        other = rings[i]
        if s in other:
            k = other.index(s)
            neighbours.update((other[k - 1], other[(k + 1) % len(other)]))
    start, size = outward_gap(s, pos, sorted(neighbours), [rings[i] for i in placed])
    toward = start + size / 2
    r = polygon_radius(n)
    centre = at(pos[s], toward, r)
    first = toward + math.pi
    for k in range(1, n):
        pos[ring[(j + k) % n]] = at(centre, first + k * TAU / n, r)


def _runs(ring: List[int], known: List[int]):
    """Return each run of atoms not yet placed as (index before it, its indices, index after
    it), going round the ring."""
    n = len(ring)
    known_set = set(known)
    runs = []
    for i in known:
        nxt = (i + 1) % n
        if nxt in known_set:
            continue
        run = []
        while nxt not in known_set:
            run.append(nxt)
            nxt = (nxt + 1) % n
        runs.append((i, run, nxt))
    return runs


def _reference(rings, ring, u_i, v_i, runs, placed, pos) -> Point:
    """Return a point on the side where the new atoms must not go: the centroid of the shared
    atoms between the two ends or, if the ring shares one bond, of the ring on the other side
    of that bond."""
    n = len(ring)
    if len(runs) == 1:
        between = []
        i = (v_i + 1) % n
        while i != u_i:
            between.append(ring[i])
            i = (i + 1) % n
        if between:
            return centroid(pos[a] for a in between)
        u, v = ring[u_i], ring[v_i]
        beside = [rings[i] for i in placed if _bonded_in(rings[i], u, v)]
        if beside:
            return centroid(pos[a] for a in beside[0])
    return centroid(pos.values())


def _bonded_in(ring: Sequence[int], u: int, v: int) -> bool:
    if u not in ring or v not in ring:
        return False
    k = ring.index(u)
    return v in (ring[k - 1], ring[(k + 1) % len(ring)])


def _arc(pos: Dict[int, Point], u: Point, v: Point, run: List[int], away_from: Point) -> None:
    """Put the atoms of `run` (in order from u) on a circular arc from u to v, with every bond
    of length 1, on the side of the line u-v away from `away_from`."""
    m = len(run) + 1
    d = distance(u, v)
    mid = mul(add(u, v), 0.5)
    chord = sub(v, u)
    normal = (-chord[1], chord[0])
    size = math.hypot(*normal)
    normal = (0.0, 1.0) if size < 1e-12 else (normal[0] / size, normal[1] / size)
    if dot(sub(away_from, mid), normal) > 0:
        normal = (-normal[0], -normal[1])
    if d >= m - 1e-9:
        # Too far apart for bonds of length 1: a straight line.
        for k, atom in enumerate(run, 1):
            pos[atom] = add(u, mul(chord, k / m))
        return
    phi = _arc_angle(d, m)
    r = 1 / (2 * math.sin(phi / (2 * m)))
    centre = sub(mid, mul(normal, r * math.cos(phi / 2)))
    start = direction(centre, u)
    # Go round in the direction that passes on the side of `normal`.
    sense = max((1, -1), key=lambda s: dot(sub(at(centre, start + s * phi / 2, r), mid), normal))
    for k, atom in enumerate(run, 1):
        pos[atom] = at(centre, start + sense * k * phi / m, r)


def _arc_angle(d: float, m: int) -> float:
    """Return the angle phi of an arc made of m chords of length 1 whose ends are d apart:
    sin(phi / 2) / sin(phi / 2m) = d, solved by bisection. The left-hand side falls from m to
    0 as phi goes from 0 to 2 pi."""
    lo, hi = 0.0, TAU
    for _ in range(80):
        phi = (lo + hi) / 2
        if math.sin(phi / 2) / math.sin(phi / (2 * m)) > d:
            lo = phi
        else:
            hi = phi
    return (lo + hi) / 2
