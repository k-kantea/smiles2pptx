"""Lengthening of bonds that would look too short between the labels of their atoms.

The lines of a bond stop short of its atom labels by the label margin, so little of a bond
between two labelled atoms is left visible: Br-Br, the bonds of FeBr3 and H2SO4, and P-Ph
with the phenyls written as Ph. A bond that is not in a ring is lengthened, by moving the
atoms on one side of it along the bond, until at least style.shortest_line of it is visible.
The side with fewer atoms moves. If that adds overlaps (counted the way untangle.py counts
them), the other side is tried, and if both do, the bond is left as it is. A bond in a ring
is always left as it is, since the ring would lose its shape.

This is the last step of the layout, after the molecule is rotated flat: labels do not rotate
with the molecule, so how much of a bond is visible depends on its direction.
"""
from __future__ import annotations

import math
from typing import Dict, FrozenSet, Optional, Tuple

from ..chem.molecule import Molecule
from ..depict.bonds import bond_lines
from ..depict.labels import atom_label, label_parts
from ..fonts import get_font
from ..style import Style
from .geometry import Point
from .untangle import _Scorer

# Extra length (points) added to the shortest visible length, so that rounding does not
# leave a bond a hair short of it.
OVER = 1e-6


def lengthen(mol: Molecule, points: Dict[int, Point], style: Style) -> Dict[int, Point]:
    """Return the points of one fragment, with bonds that are not in a ring made long enough
    to be visible (style.shortest_line) wherever that adds no overlaps."""
    if style.shortest_line <= 0:
        return points
    scorer = _Scorer(mol, sorted(points), style)
    pos = dict(points)
    for k, bond in enumerate(mol.bonds):
        i, j = bond.begin, bond.end
        if i not in pos or (i, j) not in scorer.sides:
            continue
        missing = _missing(mol, k, pos, style)
        if missing <= 0:
            continue
        before = scorer.cost(pos, None)[0]
        sides = sorted(((i, j), (j, i)), key=lambda e: (len(scorer.sides[e]), e != (i, j)))
        for near, far in sides:
            trial = _pushed(pos, near, far, scorer.sides[(near, far)], missing, style)
            if trial is not None and scorer.cost(trial, None)[0] <= before:
                pos = trial
                break
    return pos


def _pushed(pos: Dict[int, Point], near: int, far: int, side: FrozenSet[int], missing: float,
            style: Style) -> Optional[Dict[int, Point]]:
    """Return `pos` with the atoms of `side` (the side of `far`) moved `missing` points away
    from `near` along the bond, or None if the two atoms coincide.

    The lines of a bond are clipped where they meet the label boxes, and the boxes move with
    their atoms, so moving one end along the bond lengthens the visible part by the same
    amount: one step is enough."""
    (x0, y0), (x1, y1) = pos[near], pos[far]
    length = math.dist((x0, y0), (x1, y1))
    if length == 0:
        return None
    step = missing / style.bond_length
    ux, uy = (x1 - x0) / length, (y1 - y0) / length
    trial = dict(pos)
    for a in side:
        trial[a] = (trial[a][0] + ux * step, trial[a][1] + uy * step)
    return trial


def _missing(mol: Molecule, k: int, pos: Dict[int, Point], style: Style) -> float:
    """Return how much longer (in points) the visible part of bond k must be to reach
    style.shortest_line, or 0 if it is long enough. When none of the bond is visible, the
    labels of its atoms may overlap along the bond, and that overlap counts as missing too."""
    target = style.shortest_line * style.bond_length + OVER
    bond = mol.bonds[k]
    ends = (bond.begin, bond.end)
    font = get_font(style.font)
    L = style.bond_length
    P: Dict[int, Point] = {a: (pos[a][0] * L, -pos[a][1] * L) for a in _around(mol, ends, pos)}
    boxes = {}
    m = style.label_margin
    for a in ends:
        label = atom_label(mol, a, pos)
        if label is not None:
            boxes[a] = [(x0 - m, y0 - m, x1 + m, y1 + m) for x0, y0, x1, y1
                        in label_parts(label, P[a], font, style.font_size, style.script_scale, y_down=True)]
    if not boxes:
        return 0.0
    seen = max((math.dist(line.p0, line.p1) for line in bond_lines(mol, k, P, boxes, style)), default=0.0)
    if seen > 0:
        return max(0.0, target - seen)
    return max(0.0, target - _gap(P[ends[0]], P[ends[1]], boxes.get(ends[0]), boxes.get(ends[1])))


def _around(mol: Molecule, ends: Tuple[int, int], pos: Dict[int, Point]):
    """Return the atoms whose positions bond_lines may read: the two ends and their neighbours."""
    out = set(ends)
    for a in ends:
        out.update(n for n in mol.neighbors(a) if n in pos)
    return out


def _gap(p: Point, q: Point, parts_p, parts_q) -> float:
    """Return the distance between the labels along the line from p to q, negative if they
    overlap along it. `parts_p` and `parts_q` are the label boxes of p and q (or None)."""
    length = math.dist(p, q)
    t_out = _leave(p, q, parts_p) if parts_p else 0.0
    t_in = 1.0 - _leave(q, p, parts_q) if parts_q else 1.0
    return (t_in - t_out) * length


def _leave(p: Point, q: Point, parts) -> float:
    """Return where the line from p (inside the first box) towards q leaves the label, as a
    fraction of p-q (may exceed 1): after leaving the first box, it continues through each
    further box (an H above or below the symbol) that it is inside by then."""
    t = _leave_box(p, q, parts[0])
    for box in parts[1:]:
        t_in, t_out = _through(p, q, box)
        if t_in <= t + 1e-9 < t_out:
            t = t_out
    return t


def _leave_box(p: Point, q: Point, box) -> float:
    """Return where the line from p (inside the box) towards q leaves the box, as a fraction
    of p-q (may exceed 1)."""
    x0, y0, x1, y1 = box
    dx, dy = q[0] - p[0], q[1] - p[1]
    t = math.inf
    if dx > 0:
        t = min(t, (x1 - p[0]) / dx)
    elif dx < 0:
        t = min(t, (x0 - p[0]) / dx)
    if dy > 0:
        t = min(t, (y1 - p[1]) / dy)
    elif dy < 0:
        t = min(t, (y0 - p[1]) / dy)
    return max(0.0, t)


def _through(p: Point, q: Point, box) -> Tuple[float, float]:
    """Return where the line from p towards q enters and leaves the box, as fractions of p-q
    ((inf, -inf) if it misses the box)."""
    x0, y0, x1, y1 = box
    dx, dy = q[0] - p[0], q[1] - p[1]
    t_in, t_out = -math.inf, math.inf
    for d, lo, hi, s in ((dx, x0, x1, p[0]), (dy, y0, y1, p[1])):
        if d == 0:
            if not lo <= s <= hi:
                return math.inf, -math.inf
            continue
        a, b = (lo - s) / d, (hi - s) / d
        t_in, t_out = max(t_in, min(a, b)), min(t_out, max(a, b))
    return (t_in, t_out) if t_in < t_out else (math.inf, -math.inf)
