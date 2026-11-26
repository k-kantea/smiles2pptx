"""Overlap counts for a molecule drawing, computed without PowerPoint.

Two counts, both of which should be 0:

- atoms as points: two non-bonded atoms closer than half a bond length (`clash`), an
  atom closer than 0.3 bond lengths to a bond it is not part of (`on_bond`), and two bonds
  that cross (`crossing`). The coordinates are first divided by the median bond length.
- with label bounding boxes: a bond that passes through the label of an atom other than
  its own two (`line_label`), and two labels whose boxes overlap (`label_label`).

`overlaps(mol, points, style)` returns both counts for a molecule and its layout.
"""
from __future__ import annotations

import math
import statistics
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple, Union

from .chem.molecule import Molecule
from .fonts import get_font
from .style import Style

Point = Tuple[float, float]
Box = Tuple[float, float, float, float]  # x0, y0, x1, y1

CLASH = 0.5      # two non-bonded atoms clash when closer than this (bond lengths)
ON_BOND = 0.3    # an atom closer than this to a bond it is not part of (bond lengths)


@dataclass
class Overlaps:
    clash: List[Tuple[int, int]] = field(default_factory=list)        # atom, atom
    on_bond: List[Tuple[int, int]] = field(default_factory=list)      # atom, bond
    crossing: List[Tuple[int, int]] = field(default_factory=list)     # bond, bond
    line_label: List[Tuple[int, int]] = field(default_factory=list)   # labelled atom, bond
    label_label: List[Tuple[int, int]] = field(default_factory=list)  # labelled atom, labelled atom

    @property
    def points(self) -> int:
        """Count overlaps, treating atoms as points."""
        return len(self.clash) + len(self.on_bond) + len(self.crossing)

    @property
    def labels(self) -> int:
        """Count overlaps that involve label bounding boxes."""
        return len(self.line_label) + len(self.label_label)


def point_overlaps(points: Sequence[Point], bonds: Sequence[Tuple[int, int]]) -> Overlaps:
    """Find overlaps, treating atoms as points. `bonds` are pairs of atom indices."""
    lengths = [math.dist(points[i], points[j]) for i, j in bonds]
    unit = statistics.median(lengths) if lengths else 1.0
    if unit <= 0:
        unit = 1.0
    P = [(x / unit, y / unit) for x, y in points]
    bonded = {frozenset(b) for b in bonds}
    result = Overlaps()
    n = len(P)
    for i in range(n):
        for j in range(i + 1, n):
            if frozenset((i, j)) not in bonded and math.dist(P[i], P[j]) < CLASH:
                result.clash.append((i, j))
    for k, (i, j) in enumerate(bonds):
        for a in range(n):
            if a in (i, j):
                continue
            if (segment_distance(P[a], P[i], P[j]) < ON_BOND
                    and math.dist(P[a], P[i]) >= CLASH and math.dist(P[a], P[j]) >= CLASH):
                result.on_bond.append((a, k))
    for k1 in range(len(bonds)):
        for k2 in range(k1 + 1, len(bonds)):
            s1, s2 = bonds[k1], bonds[k2]
            if set(s1) & set(s2):
                continue
            if segments_cross(P[s1[0]], P[s1[1]], P[s2[0]], P[s2[1]]):
                result.crossing.append((k1, k2))
    return result


def label_overlaps(points: Sequence[Point], bonds: Sequence[Tuple[int, int]],
                   boxes: Dict[int, Union[Box, Sequence[Box]]], result: Optional[Overlaps] = None) -> Overlaps:
    """Find overlaps that involve labels. `boxes` are {atom: box}, or {atom: [boxes]} for a
    label in parts (an H above or below its symbol), in the units of `points`. Each atom
    counts once, however many parts its label has."""
    result = result or Overlaps()
    parts = {a: as_parts(b) for a, b in boxes.items()}
    for k, (i, j) in enumerate(bonds):
        for a, box in parts.items():
            if a not in (i, j) and segment_in_parts(points[i], points[j], box):
                result.line_label.append((a, k))
    ids = sorted(parts)
    for x in range(len(ids)):
        for y in range(x + 1, len(ids)):
            if parts_meet(parts[ids[x]], parts[ids[y]]):
                result.label_label.append((ids[x], ids[y]))
    return result


def overlaps(mol: Molecule, points: Sequence[Point], style: Optional[Style] = None) -> Overlaps:
    """Return both counts for a perceived molecule and its layout (bond length 1)."""
    from .depict.labels import label_boxes

    style = style or Style()
    bonds = [(b.begin, b.end) for b in mol.bonds]
    result = point_overlaps(points, bonds)
    size = style.font_size / style.bond_length
    boxes = label_boxes(mol, points, get_font(style.font), size, style.script_scale)
    return label_overlaps(points, bonds, boxes, result)


# --- geometry --------------------------------------------------------------------------

def segment_distance(p: Point, a: Point, b: Point) -> float:
    dx, dy = b[0] - a[0], b[1] - a[1]
    length2 = dx * dx + dy * dy
    t = 0.0 if length2 == 0 else max(0.0, min(1.0, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / length2))
    return math.hypot(p[0] - (a[0] + t * dx), p[1] - (a[1] + t * dy))


def _cross(o: Point, a: Point, b: Point) -> float:
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def segments_cross(p1: Point, p2: Point, p3: Point, p4: Point) -> bool:
    """Return True if the segments cross at a point inside both (touching does not count)."""
    d1, d2 = _cross(p3, p4, p1), _cross(p3, p4, p2)
    d3, d4 = _cross(p1, p2, p3), _cross(p1, p2, p4)
    return d1 * d2 < -1e-9 and d3 * d4 < -1e-9


def segment_in_box(p0: Point, p1: Point, box: Box) -> bool:
    """Return True if a part of the segment with positive length lies inside the box (Liang-Barsky)."""
    x0, y0, x1, y1 = box
    dx, dy = p1[0] - p0[0], p1[1] - p0[1]
    t_in, t_out = 0.0, 1.0
    for p, q in ((-dx, p0[0] - x0), (dx, x1 - p0[0]), (-dy, p0[1] - y0), (dy, y1 - p0[1])):
        if abs(p) < 1e-12:
            if q < 0:
                return False
            continue
        t = q / p
        if p < 0:
            t_in = max(t_in, t)
        else:
            t_out = min(t_out, t)
    return t_out - t_in > 1e-6


def boxes_meet(a: Box, b: Box) -> bool:
    return a[0] < b[2] and b[0] < a[2] and a[1] < b[3] and b[1] < a[3]


def as_parts(boxes: Union[Box, Sequence[Box]]) -> List[Box]:
    """Return the boxes of one label, given either as a single box or as a list of parts."""
    return [boxes] if len(boxes) == 4 and not isinstance(boxes[0], (tuple, list)) else list(boxes)


def segment_in_parts(p0: Point, p1: Point, parts: Sequence[Box]) -> bool:
    return any(segment_in_box(p0, p1, box) for box in parts)


def parts_meet(a: Sequence[Box], b: Sequence[Box]) -> bool:
    return any(boxes_meet(x, y) for x in a for y in b)


def circle_meets_box(centre: Point, radius: float, box: Box) -> bool:
    """Return True if the disc (the circle and its interior) reaches into the box."""
    x = min(max(centre[0], box[0]), box[2])
    y = min(max(centre[1], box[1]), box[3])
    return math.dist(centre, (x, y)) < radius


# --- on a drawing ----------------------------------------------------------------------------

def drawing_overlaps(shapes, font_name: str = "Arial") -> Tuple[int, int]:
    """Return the (line_label, label_label) counts for drawn shapes (smiles2pptx.figure), in points.

    A line, an arrow, or the centre line of a wedge or hashed wedge that passes through the
    capital letters of a label counts once for each label. Bond lines already stop short of
    the labels of their own atoms, so every such crossing is an overlap. An aromatic ring
    circle that reaches the capital letters of a label counts the same way. Two labels whose
    capital letters overlap count once."""
    from .figure import Arrow, Circle, Hash, Line, Text, Wedge, walk

    font = get_font(font_name)
    segments, circles, boxes = [], [], []
    for s in walk(shapes):
        if isinstance(s, (Line, Arrow)):
            segments.append((s.p0, s.p1))
        elif isinstance(s, Wedge):
            p = s.points
            segments.append((_mid(p[0], p[3]), _mid(p[1], p[2])))
        elif isinstance(s, Hash):
            segments.append((_mid(*s.segments[0]), _mid(*s.segments[-1])))
        elif isinstance(s, Circle):
            circles.append((s.centre, s.radius))
        elif isinstance(s, Text):
            baseline = s.y + font.ascent_height(s.size)
            boxes.append((s.x, baseline - font.capital_height(s.size), s.x + s.width, baseline))
    line_label = sum(1 for a, b in segments for box in boxes if segment_in_box(a, b, box))
    line_label += sum(1 for centre, radius in circles for box in boxes if circle_meets_box(centre, radius, box))
    label_label = sum(1 for k, a in enumerate(boxes) for b in boxes[k + 1:] if boxes_meet(a, b))
    return line_label, label_label


def seen_lengths(group, bonds: int) -> List[float]:
    """Return the visible length of each bond in one molecule drawing (a group made by
    draw_molecule), in points: bonds 1 to `bonds` first, then the wedges to drawn hydrogens.

    Bond lines stop short of the labels of their atoms, so a bond between two labelled atoms
    may show only as a short line, or not at all. For a bond drawn as lines, the length is
    that of the longest line (a double bond has a shorter inner line; a line split at the
    middle into two colours counts as a whole). For a wedge or a hashed wedge, it is the
    length of the centre line. A bond with nothing drawn has length 0. A wedge to a drawn
    hydrogen with nothing drawn has no number of its own, so it is not in the list."""
    return _seen(group, bonds)[1]


def short_lines(group, bonds: int, minimum: float) -> List[int]:
    """Return the bonds that show less than `minimum` points (seen_lengths). Bonds are numbered
    from 1 as in the shape names, followed by the wedges to drawn hydrogens."""
    numbers, seen = _seen(group, bonds)
    return [n for n, length in zip(numbers, seen) if length < minimum]


def _seen(group, bonds: int) -> Tuple[List[int], List[float]]:
    from .figure import Hash, Line, Wedge, walk

    lines: Dict[int, List[Tuple[Point, Point]]] = {}
    dashes: Dict[int, List[Point]] = {}
    for s in walk([group]):
        parts = s.name.split()
        if len(parts) != 2 or not parts[1].isdigit():
            continue
        number = int(parts[1])
        if isinstance(s, Line) and parts[0] == "Bond":
            lines.setdefault(number, []).append((s.p0, s.p1))
        elif isinstance(s, Wedge) and parts[0] == "Wedge":
            p = s.points
            lines.setdefault(number, []).append((_mid(p[0], p[3]), _mid(p[1], p[2])))
        elif isinstance(s, Hash) and parts[0] == "Hash":
            # A two-colour hashed wedge is two shapes: measure the span of all their dashes.
            dashes.setdefault(number, []).extend(_mid(*d) for d in s.segments)
    for number, points in dashes.items():
        lines[number] = [max(((a, b) for a in points for b in points), key=lambda ab: math.dist(*ab))]
    numbers = list(range(1, bonds + 1)) + sorted(n for n in lines if n > bonds)
    return numbers, [_longest_line(lines.get(n, [])) for n in numbers]


def _longest_line(segments: List[Tuple[Point, Point]]) -> float:
    """Return the length of the longest line formed by the segments. Segments on the same
    straight line (the two colours of one line, or of one wedge) are added together."""
    totals: List[Tuple[Point, Point, float]] = []   # a point on the line, its direction, total length
    for a, b in segments:
        length = math.dist(a, b)
        if length == 0:
            continue
        u = ((b[0] - a[0]) / length, (b[1] - a[1]) / length)
        for k, (p, v, total) in enumerate(totals):
            parallel = abs(u[0] * v[1] - u[1] * v[0]) < 1e-6
            if parallel and abs((a[0] - p[0]) * v[1] - (a[1] - p[1]) * v[0]) < 1e-3:
                totals[k] = (p, v, total + length)
                break
        else:
            totals.append((a, u, length))
    return max((t for _, _, t in totals), default=0.0)


def outside(shapes, width: float, height: float) -> int:
    """Count the shapes (not groups) that extend beyond a slide of this size, in points."""
    from .figure import walk

    count = 0
    for s in walk(shapes):
        x0, y0, x1, y1 = s.bbox()
        if x0 < 0 or y0 < 0 or x1 > width or y1 > height:
            count += 1
    return count


def _mid(a: Point, b: Point) -> Point:
    return ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)
