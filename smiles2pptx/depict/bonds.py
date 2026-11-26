"""Bond shapes, in points with y pointing down (the drawing stage).

- a single bond is one line; a triple bond is three parallel lines
- a double bond in a ring has its second line inside the ring, shortened at both ends.
  A double bond to an atom with no other bonds (C=O, =CH2) is two lines centred on the
  bond. Any other double bond has its second line on the side with more of the
  neighbouring bonds
- a wedge is a long, thin filled shape, narrow at the stereocentre; a hashed wedge is a
  row of short lines across the bond that get longer away from the stereocentre
- every line stops short of the labels of its atoms, by the label margin of the style
- when aromatic rings are drawn with circles, each bond of such a ring is a single line
  and the ring has a circle inside
- with coloured bonds, each line, wedge and hashed wedge is split at the midpoint of its
  bond, and each half takes the colour of the atom at its end
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Tuple

from ..chem.molecule import Molecule
from ..figure import Circle, Hash, Line, Wedge
from ..style import Style

Point = Tuple[float, float]
Box = Tuple[float, float, float, float]  # x0, y0, x1, y1
Segment = Tuple[Point, Point]
BLACK = "000000"
Colours = Tuple[str, str]   # the colours of the two atoms of a bond


def clip_out(p0: Point, p1: Point, box: Box) -> Optional[Segment]:
    """Return the segment without the part inside the box, when that part is at one end.
    Return None when the whole segment is inside. A segment that passes through the box is
    kept whole (that is an overlap, which check.py counts)."""
    x0, y0, x1, y1 = box
    dx, dy = p1[0] - p0[0], p1[1] - p0[1]
    t_in, t_out = 0.0, 1.0
    for p, q in ((-dx, p0[0] - x0), (dx, x1 - p0[0]), (-dy, p0[1] - y0), (dy, y1 - p0[1])):
        if abs(p) < 1e-12:
            if q < 0:
                return p0, p1
            continue
        t = q / p
        if p < 0:
            t_in = max(t_in, t)
        else:
            t_out = min(t_out, t)
    if t_in >= t_out:
        return p0, p1
    if t_in <= 1e-9 and t_out >= 1 - 1e-9:
        return None
    if t_in <= 1e-9:
        return _at(p0, p1, t_out), p1
    if t_out >= 1 - 1e-9:
        return p0, _at(p0, p1, t_in)
    return p0, p1


def clip_ends(seg: Optional[Segment], ends: Sequence[int], boxes: Dict[int, Sequence[Box]]) -> Optional[Segment]:
    """Return `seg` without the parts inside the labels of the atoms `ends`. `boxes` are
    {atom: [boxes]}: the symbol first, then the H written above or below it (labels.label_parts),
    so that a line that runs from the symbol into the H also stops at the H."""
    for k in ends:
        for box in boxes.get(k, ()):
            if seg is not None:
                seg = clip_out(seg[0], seg[1], box)
    return seg


def bond_lines(mol: Molecule, k: int, P: Sequence[Point], boxes: Dict[int, Sequence[Box]], style: Style,
               single: bool = False, colours: Colours = (BLACK, BLACK)) -> List[Line]:
    """Return the lines of bond k (not a wedge). `P` are the atom positions in points; `boxes`
    are the label boxes ({atom: [boxes]}, as in clip_ends), already widened by the label margin.
    `single` draws one line whatever the bond order (a bond in a ring drawn with a circle).
    `colours` are those of the begin and end atoms."""
    bond = mol.bonds[k]
    i, j = bond.begin, bond.end
    p, q = P[i], P[j]
    length = math.dist(p, q)
    if length == 0:
        return []
    nx, ny = -(q[1] - p[1]) / length, (q[0] - p[0]) / length
    gap = style.double_bond_gap * style.bond_length
    order = 1 if single else bond.order or 1
    segs: List[Segment] = []
    if order == 2:
        side = double_bond_side(mol, k, P)
        if side == 0:
            segs += [(_off(p, nx, ny, gap / 2), _off(q, nx, ny, gap / 2)),
                     (_off(p, nx, ny, -gap / 2), _off(q, nx, ny, -gap / 2))]
        else:
            trim = style.inner_line_trim * style.bond_length
            ux, uy = (q[0] - p[0]) / length, (q[1] - p[1]) / length
            a, b = _off(p, nx, ny, side * gap), _off(q, nx, ny, side * gap)
            segs += [(p, q), ((a[0] + ux * trim, a[1] + uy * trim), (b[0] - ux * trim, b[1] - uy * trim))]
    elif order >= 3:
        segs += [(p, q), (_off(p, nx, ny, gap), _off(q, nx, ny, gap)),
                 (_off(p, nx, ny, -gap), _off(q, nx, ny, -gap))]
    else:
        segs.append((p, q))
    lines = []
    for seg in segs:
        seg = clip_ends(seg, (i, j), boxes)
        if seg is not None:
            lines += [Line(a, b, style.line_width, colour=colour, name=f"Bond {k + 1}")
                      for (a, b), colour in halves(seg, p, q, colours)]
    return lines


def halves(seg: Segment, p: Point, q: Point, colours: Colours) -> List[Tuple[Segment, str]]:
    """Split `seg`, which runs along the bond from p to q, into parts coloured like the nearer
    atom: cut where it crosses the midpoint of the bond, or kept whole when the colours match."""
    if colours[0] == colours[1]:
        return [(seg, colours[0])]
    t0, t1 = _along(seg[0], p, q), _along(seg[1], p, q)
    if min(t0, t1) >= 0.5 - 1e-9 or max(t0, t1) <= 0.5 + 1e-9:
        return [(seg, colours[0] if t0 + t1 < 1 else colours[1])]
    middle = _at(seg[0], seg[1], (0.5 - t0) / (t1 - t0))
    first, second = (colours[0], colours[1]) if t0 < t1 else (colours[1], colours[0])
    return [((seg[0], middle), first), ((middle, seg[1]), second)]


def ring_circle(ring: Sequence[int], P: Sequence[Point], boxes: Dict[int, Sequence[Box]], style: Style,
                number: int) -> Circle:
    """Return the circle inside an aromatic ring. It is centred on the mean position of the
    ring atoms, and its radius is a fraction (circle_radius) of the mean distance from there
    to the bond midpoints. Where a label reaches into the ring (the + of N+), the circle
    shrinks to stay clear of `boxes` (the labels widened by the label margin, {atom: [boxes]}),
    but not below circle_min_radius."""
    n = len(ring)
    centre = (sum(P[a][0] for a in ring) / n, sum(P[a][1] for a in ring) / n)
    to_bonds = sum(math.dist(centre, ((P[a][0] + P[b][0]) / 2, (P[a][1] + P[b][1]) / 2))
                   for a, b in zip(ring, list(ring[1:]) + list(ring[:1]))) / n
    radius = style.circle_radius * to_bonds
    for box in (box for parts in boxes.values() for box in parts):
        nearest = (min(max(centre[0], box[0]), box[2]), min(max(centre[1], box[1]), box[3]))
        radius = min(radius, math.dist(centre, nearest))
    radius = max(radius, style.circle_min_radius * to_bonds)
    return Circle(centre, radius, style.line_width, name=f"Ring {number}")


def double_bond_side(mol: Molecule, k: int, P: Sequence[Point]) -> int:
    """Return the side for the second line of double bond k: 1 or -1 (along the left normal of
    the bond from begin to end, or against it), or 0 for two lines centred on the bond."""
    bond = mol.bonds[k]
    i, j = bond.begin, bond.end
    p, q = P[i], P[j]
    length = math.dist(p, q)
    nx, ny = -(q[1] - p[1]) / length, (q[0] - p[0]) / length
    mid = ((p[0] + q[0]) / 2, (p[1] + q[1]) / 2)
    rings = [r for r in mol.rings if _ring_has_bond(r, i, j)]
    if rings:
        ring = min(rings, key=lambda r: (abs(len(r) - 6), len(r)))
        cx = sum(P[a][0] for a in ring) / len(ring)
        cy = sum(P[a][1] for a in ring) / len(ring)
        return 1 if (cx - mid[0]) * nx + (cy - mid[1]) * ny > 0 else -1
    ni, nj = mol.neighbors(i), mol.neighbors(j)
    if len(ni) == 1 or len(nj) == 1:
        return 0
    s = 0.0
    for a, others in ((i, ni), (j, nj)):
        for n in others:
            if n in (i, j):
                continue
            s += (P[n][0] - mid[0]) * nx + (P[n][1] - mid[1]) * ny
    return 1 if s >= 0 else -1


def wedge_shapes(number: int, centre: Point, far: Point, kind: str, seg: Segment, style: Style,
                 colours: Colours = (BLACK, BLACK)) -> list:
    """Return the wedge or hashed wedge along `seg` as one shape, or as two when the atoms at
    its ends (`colours`: at the centre, at the far end) differ in colour."""
    if colours[0] == colours[1]:
        return [wedge_shape(number, centre, far, kind, seg, style, colours[0])]
    if kind == "wedge":
        return [wedge_shape(number, centre, far, kind, part, style, colour)
                for part, colour in halves(seg, centre, far, colours)]
    whole = wedge_shape(number, centre, far, kind, seg, style)
    near = [s for s in whole.segments if _along(_mid(s), centre, far) <= 0.5]
    rest = [s for s in whole.segments if _along(_mid(s), centre, far) > 0.5]
    return [Hash(part, whole.width, colour=colour, name=whole.name)
            for part, colour in ((near, colours[0]), (rest, colours[1])) if part]


def wedge_shape(number: int, centre: Point, far: Point, kind: str, seg: Segment, style: Style,
                colour: str = BLACK):
    """Return a wedge or hashed wedge along `seg` (the bond from `centre` to `far`, already
    clipped), named with `number` (the bond number, counting from 1).

    The width grows from the centre to the far end as it would along the whole bond, so a
    clipped wedge keeps the angle of an unclipped one."""
    length = math.dist(centre, far)
    nx, ny = -(far[1] - centre[1]) / length, (far[0] - centre[0]) / length
    wide = style.wedge_width
    s0, s1 = seg
    t0, t1 = math.dist(centre, s0) / length, math.dist(centre, s1) / length
    if kind == "wedge":
        w0, w1 = max(wide * t0, style.line_width) / 2, wide * t1 / 2
        return Wedge([_off(s0, nx, ny, w0), _off(s1, nx, ny, w1), _off(s1, nx, ny, -w1), _off(s0, nx, ny, -w0)],
                     colour=colour, name=f"Wedge {number}")
    step = style.hash_spacing
    n = max(2, int(math.dist(s0, s1) / step) + 1)
    lines = []
    for m in range(n):
        f = m / (n - 1)
        c = (s0[0] + (s1[0] - s0[0]) * f, s0[1] + (s1[1] - s0[1]) * f)
        half = max(wide * (t0 + (t1 - t0) * f), style.line_width * 2) / 2
        lines.append((_off(c, nx, ny, half), _off(c, nx, ny, -half)))
    return Hash(lines, style.line_width, colour=colour, name=f"Hash {number}")


def _ring_has_bond(ring: Sequence[int], i: int, j: int) -> bool:
    n = len(ring)
    return any({ring[m], ring[(m + 1) % n]} == {i, j} for m in range(n))


def _off(p: Point, nx: float, ny: float, d: float) -> Point:
    return (p[0] + nx * d, p[1] + ny * d)


def _at(p0: Point, p1: Point, t: float) -> Point:
    return (p0[0] + (p1[0] - p0[0]) * t, p0[1] + (p1[1] - p0[1]) * t)


def _mid(seg: Segment) -> Point:
    return _at(seg[0], seg[1], 0.5)


def _along(x: Point, p: Point, q: Point) -> float:
    """Return the position of x along the line from p (0) to q (1)."""
    dx, dy = q[0] - p[0], q[1] - p[1]
    return ((x[0] - p[0]) * dx + (x[1] - p[1]) * dy) / (dx * dx + dy * dy)
