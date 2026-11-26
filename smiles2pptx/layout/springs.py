"""Spring layout for bridged and cage ring systems (norbornane, adamantane, cubane).

The arcs of rings.py cannot lay out a ring system flat when its rings share more than one
bond's worth of atoms: its atoms end up on top of each other. Here every pair of atoms in
the system is joined by a spring whose length grows with the number of bonds between them,
and the atoms are moved until the springs are as relaxed as possible (stress majorization:
Gansner, Koren and North, 2004).

The result depends on where the atoms start, so several starts are tried: the arc
placement, the atoms on a circle, and then fixed pseudo-random points. rings.py keeps the
result with the fewest overlaps.
"""
from __future__ import annotations

import math
from collections import deque
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

from ..check import point_overlaps
from .geometry import Point

ROUNDS = 500
STILL = 1e-4                                   # stop when no atom moved more than this in a round
SPRING_FACTORS = (math.sqrt(3) / 2, 0.75)      # spring length per bond between two atoms
RANDOM_STARTS = 8                              # tried when the first starts leave overlaps
PATIENCE = 4                                   # stop after this many batches that found nothing better

Bonds = Sequence[Tuple[int, int]]


def spring_batches(start: Dict[int, Point], bonds: Bonds, outside: Bonds = ()) -> Iterator[List[Dict[int, Point]]]:
    """Yield spring placements of the atoms of a ring system, one batch per start: the arc
    placement, a circle, then fixed pseudo-random points. Each start is relaxed with the
    system alone and with the first atom of each substituent (`outside` holds the bonds to
    them), which pushes the atoms that carry substituents outwards."""
    system = sorted(start)
    with_outside = dict(start)
    middle = (sum(p[0] for p in start.values()) / len(start), sum(p[1] for p in start.values()) / len(start))
    for a, b in outside:
        dx, dy = start[a][0] - middle[0], start[a][1] - middle[1]
        size = math.hypot(dx, dy) or 1.0
        with_outside[b] = (start[a][0] + dx / size, start[a][1] + dy / size)
    variants = [(list(bonds), sorted(start))]
    if outside:
        variants.append((list(bonds) + list(outside), sorted(with_outside)))

    def run(begin_of) -> List[Dict[int, Point]]:
        out = []
        for every, atoms in variants:
            begin = begin_of(atoms)
            for factor in SPRING_FACTORS:
                relaxed = relax(begin, every, factor)
                out.append({a: relaxed[a] for a in system})
        return out

    yield run(lambda atoms: {a: with_outside[a] for a in atoms})
    yield run(lambda atoms: {a: at_circle(k, len(atoms)) for k, a in enumerate(atoms)})
    numbers = _pseudo_random(len(system))
    for _ in range(RANDOM_STARTS):
        values = [next(numbers) for _ in range(2 * len(with_outside))]
        yield run(lambda atoms: {a: (values[2 * k] * len(atoms) / 3, values[2 * k + 1] * len(atoms) / 3)
                                 for k, a in enumerate(atoms)})


def at_circle(k: int, n: int) -> Point:
    r = max(n, 3) / 6
    return (r * math.cos(2 * math.pi * k / n), r * math.sin(2 * math.pi * k / n))


def count_overlaps(pos: Dict[int, Point], bonds: Bonds) -> int:
    """Count the overlaps among the given atoms, treating atoms as points (check.py)."""
    atoms = sorted(pos)
    index = {a: k for k, a in enumerate(atoms)}
    return point_overlaps([pos[a] for a in atoms], [(index[i], index[j]) for i, j in bonds]).points


def worst_bond(pos: Dict[int, Point], bonds: Bonds) -> float:
    return max((abs(math.dist(pos[i], pos[j]) - 1) for i, j in bonds), default=0.0)


def _pseudo_random(seed: int) -> Iterator[float]:
    """Yield numbers in [0, 1), the same sequence every time for the same seed."""
    x = seed
    while True:
        x = (1103515245 * x + 12345) % 2 ** 31
        yield x / 2 ** 31


def spring_length(bonds_apart: int, factor: float) -> float:
    """Return the spring length: 1 for a bond, and `factor` per bond for atoms further apart."""
    return 1.0 if bonds_apart == 1 else factor * bonds_apart


def relax(start: Dict[int, Point], bonds: Bonds, factor: float) -> Dict[int, Point]:
    """Return the atoms of `start` moved so that the springs between them are as relaxed as
    possible, scaled so that the median bond length is 1."""
    atoms = sorted(start)
    index = {a: k for k, a in enumerate(atoms)}
    n = len(atoms)
    adjacency: List[List[int]] = [[] for _ in atoms]
    for i, j in bonds:
        adjacency[index[i]].append(index[j])
        adjacency[index[j]].append(index[i])
    apart = [_bonds_apart(adjacency, k) for k in range(n)]
    x = [list(start[a]) for a in atoms]
    # Nudge every atom by a tiny fixed offset, so that atoms that start at the same point are
    # set apart, the same way every time.
    for k in range(n):
        x[k][0] += 1e-3 * math.cos(k + 1)
        x[k][1] += 1e-3 * math.sin(k + 1)
    springs = [[(j, 1.0 / apart[i][j] ** 2, spring_length(apart[i][j], factor))
                for j in range(n) if j != i and apart[i][j] is not None] for i in range(n)]
    for _ in range(ROUNDS):
        moved = 0.0
        for i in range(n):
            sx = sy = sw = 0.0
            xi = x[i]
            for j, w, length in springs[i]:
                xj = x[j]
                dx, dy = xi[0] - xj[0], xi[1] - xj[1]
                dist = math.hypot(dx, dy) or 1e-9
                sx += w * (xj[0] + length * dx / dist)
                sy += w * (xj[1] + length * dy / dist)
                sw += w
            if sw == 0:
                continue
            nx, ny = sx / sw, sy / sw
            moved = max(moved, abs(nx - xi[0]), abs(ny - xi[1]))
            x[i] = [nx, ny]
        if moved < STILL:
            break
    lengths = sorted(math.dist(x[index[i]], x[index[j]]) for i, j in bonds)
    unit = lengths[len(lengths) // 2] if lengths else 1.0
    unit = unit if unit > 1e-9 else 1.0
    return {a: (x[index[a]][0] / unit, x[index[a]][1] / unit) for a in atoms}


def _bonds_apart(adjacency: List[List[int]], source: int) -> List[Optional[int]]:
    dist: List[Optional[int]] = [None] * len(adjacency)
    dist[source] = 0
    queue = deque([source])
    while queue:
        a = queue.popleft()
        for b in adjacency[a]:
            if dist[b] is None:
                dist[b] = dist[a] + 1
                queue.append(b)
    return dist
