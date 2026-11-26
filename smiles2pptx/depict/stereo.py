"""Wedges and hashed wedges for stereocentres.

Each stereocentre (Atom.stereocentre) gets one mark, with its narrow end on the centre: a
wedge (the bond points towards the viewer) or a hashed wedge (it points away). `candidates`
chooses where the mark goes: a substituent first, and at a ring fusion a hydrogen drawn for
the mark (`hydrogen_place`). Whether the mark is a wedge or a hashed wedge then follows from
the "@" or "@@" in the SMILES and the layout (`mark_kind`); only one of the two is correct.
The label of a centre whose hydrogen is drawn this way does not repeat that hydrogen
(`hydrogen_drawn`: N+, not NH+).

`drawn_chirality` reads a drawing back as a reader would and returns "@" or "@@" for a
centre, with its neighbours in SMILES order. The marks are checked against it.

Coordinates are those of the layout: the bond length is 1 and y points up.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

from ..check import boxes_meet, segment_in_box, segments_cross
from ..chem.molecule import IMPLICIT_H, Molecule
from ..figure import Run
from ..fonts import get_font
from ..style import Style
from .labels import Label, atom_label, label_box, label_boxes

Point = Tuple[float, float]
# Stands for the lone pair in the neighbour list of a centre with three neighbours (sulfoxides).
LONE_PAIR = -2
# A drawing whose signed volume is smaller than this cannot be read (bonds in a line).
FLAT = 1e-6
# Three bonds lie within half a turn when the widest gap between them is wider than half a
# turn minus this (radians, 10 degrees).
HALF_TURN_SLACK = math.radians(10)


@dataclass(frozen=True)
class Mark:
    """A wedge or a hashed wedge from `centre` to the atom `to`, or to a hydrogen drawn for
    the mark at `hydrogen` (then `to` is None)."""
    centre: int
    kind: str                         # "wedge" or "hash"
    to: Optional[int] = None
    hydrogen: Optional[Point] = None

    @property
    def z(self) -> int:
        return 1 if self.kind == "wedge" else -1


Target = Union[int, Point]   # a neighbour, or the position of a hydrogen to draw


def stereo_marks(mol: Molecule, points: Sequence[Point]) -> List[Mark]:
    """Return one mark for each stereocentre, on the first candidate that reads back correctly.
    A bond carries at most one mark, and a centre with no such candidate gets none. `used`
    holds the bonds with a mark (as frozensets) and the marks to drawn hydrogens."""
    marks: List[Mark] = []
    used = set()
    for c in range(len(mol.atoms)):
        if not mol.atoms[c].stereocentre:
            continue
        for target in candidates(mol, points, c, used):
            to, hydrogen = (target, None) if isinstance(target, int) else (None, tuple(target))
            kind = mark_kind(mol, points, c, to, hydrogen)
            if kind is not None:
                marks.append(Mark(c, kind, to, hydrogen))
                if to is not None:
                    used.add(frozenset((c, to)))
                else:
                    used.add(marks[-1])
                break
    return marks


def hydrogen_drawn(marks: Iterable[Mark]) -> set:
    """Return the centres whose hydrogen is drawn as a label of its own, at the end of their mark."""
    return {m.centre for m in marks if m.to is None and m.hydrogen is not None}


def candidates(mol: Molecule, points: Sequence[Point], centre: int, used: set) -> Iterable[Target]:
    """Yield the places where a mark of `centre` may go, best first, following textbook
    drawings of steroids (IUPAC 2006): the wedge goes on a substituent, or at a ring fusion on
    a hydrogen drawn for it, so that it lies outside the rings and its wide end is not on
    another stereocentre.

    1. a single bond in no ring, to an atom that is not a stereocentre
    2. the hydrogen in the bracket atom of the centre ([C@H]), drawn where `hydrogen_place`
       finds room (outside the rings when possible)
    3. a single bond in no ring, to another stereocentre
    4. a single bond in a ring (to an atom that is not a stereocentre first)

    2 comes before 3 because a wedge whose wide end is on another stereocentre can also be read
    as a mark of that centre (C17 of a steroid, whose only substituent is C20, gets its H).
    Within 1, 3 and 4, the bond to the atom with fewer neighbours comes first (a CH3 or OH at
    the end of a branch), then the lower atom index. A bond that already carries a mark is
    skipped."""
    free = []
    for j in mol.neighbors(centre):
        bond = mol.bond_between(centre, j)
        if bond.order == 1 and frozenset((centre, j)) not in used:
            free.append(j)
    ring = ring_bonds(mol)
    stereo = {j for j in free if mol.atoms[j].stereocentre}
    key = lambda j: (j in stereo, len(mol.neighbors(j)), j)   # noqa: E731
    chain = sorted((j for j in free if frozenset((centre, j)) not in ring), key=key)
    yield from (j for j in chain if j not in stereo)
    hydrogen = hydrogen_place(mol, points, centre, used)
    if hydrogen is not None:
        yield hydrogen
    yield from (j for j in chain if j in stereo)
    yield from sorted((j for j in free if frozenset((centre, j)) in ring), key=key)


def hydrogen_place(mol: Molecule, points: Sequence[Point], centre: int, used: set) -> Optional[Point]:
    """Return where to draw the hydrogen in the bracket atom of `centre` for its mark, or None.

    Only a centre with three neighbours and one hydrogen in its bracket atom gets one. The
    hydrogen goes in the middle of a gap between two of its bonds. A position is kept only
    when its label and wedge meet no line, label or other drawn hydrogen, and at least
    Style.shortest_line of the wedge shows between the label of the centre as drawn (N+) and
    the H. Outside the rings, distances from a labelled centre are measured from the edge of
    its label (`_Scene.reach`), so that the wedge shows as long as one from a carbon; inside a
    ring they are not, so that the H stays near its atom. Labels count with all their
    hydrogens, including the label of the centre (NH+), although it is drawn without this one
    (N+): a gap on the side where the label would have its H is not used (`_Scene.beside_h`),
    because the H would then line up with N+ and the + would sit on the wedge (N+···H, which
    reads as a hydrogen bond). An H that the label of the centre would have above or below its
    symbol does not count: the mark draws that H, and the + is not beside it. Positions are
    tried in this order, and the first gap with a clear position wins (widest gap first within
    each):

    1. a gap outside the rings, at OUTSIDE_LENGTHS (longest first, so that the wedge is as
       long as one to an OH). At a ring fusion, this is the gap facing away from the shared bond.
    2. a gap inside a ring, at INSIDE_LENGTHS. Where two rings bend round (C9 and C8 of a
       steroid on a hexagonal grid), the outside positions of two neighbouring atoms coincide,
       or fall on a methyl; textbooks then draw the hydrogen short, inside a ring.

    Return None when no position is clear; the mark then goes on a bond."""
    atom = mol.atoms[centre]
    neighbours = mol.neighbors(centre)
    if atom.h_count != 1 or IMPLICIT_H not in atom.neighbor_order or len(neighbours) != 3:
        return None
    c = points[centre]
    angles = sorted(math.atan2(points[j][1] - c[1], points[j][0] - c[0]) for j in neighbours)
    scene = _Scene(mol, points, used)
    best = None
    for k in range(3):
        width = (angles[(k + 1) % 3] - angles[k]) % (2 * math.pi)
        direction = angles[k] + width / 2
        ux, uy = math.cos(direction), math.sin(direction)
        probe = (c[0] + RING_PROBE * ux, c[1] + RING_PROBE * uy)
        inside = any(_inside(probe, [points[a] for a in r]) for r in mol.rings)
        if scene.beside_h(centre, (ux, uy)):
            continue
        reach = 0.0 if inside else scene.reach(centre, (ux, uy))
        for length in INSIDE_LENGTHS if inside else OUTSIDE_LENGTHS:
            h = (c[0] + (reach + length) * ux, c[1] + (reach + length) * uy)
            if scene.clear(centre, h) and scene.seen(centre, h) >= scene.shortest:
                score = (inside, -round(math.degrees(width)))
                if best is None or score < best[0]:
                    best = (score, h)
                break
    return None if best is None else best[1]


# Distances from a centre to a hydrogen drawn for its mark (bond lengths), in the order tried:
# outside the rings, and inside a ring (shorter, so that the H sits near its atom and not in
# the middle of the ring, where it would read as a label of the ring).
OUTSIDE_LENGTHS = (1.0, 0.8, 0.6)
INSIDE_LENGTHS = (0.8, 0.7, 0.6)
# A gap between two bonds is inside a ring when the point this far along its middle line is
# (bond lengths).
RING_PROBE = 0.5


class _Scene:
    """The lines and labels of a molecule in layout units (bond length 1), for testing whether
    a hydrogen drawn at a given position meets any of them. The sizes are those of the default
    style: they all grow with the scale, so in layout units only the font changes them, and
    only slightly. Labels include all their hydrogens (see hydrogen_place)."""

    def __init__(self, mol: Molecule, points: Sequence[Point], used: set):
        style = Style()
        self.font = get_font(style.font)
        self.size = style.font_size / style.bond_length
        self.scale = style.script_scale
        self.margin = style.label_margin / style.bond_length
        self.gap = style.double_bond_gap
        self.points = points
        self.bonds = [(b.begin, b.end, (b.order or 1) > 1) for b in mol.bonds]
        # {atom: [boxes]}: the label, and its H when that is written above or below the symbol.
        self.boxes = {a: [_grow(b, self.margin) for b in parts]
                      for a, parts in label_boxes(mol, points, self.font, self.size, self.scale).items()}
        self.drawn = [(m.centre, m.hydrogen) for m in used if isinstance(m, Mark) and m.hydrogen is not None]
        self.mol = mol
        self.shortest = style.shortest_line

    def centre_box(self, centre: int):
        """Return the label box of `centre` as drawn with its H at the end of the mark (N+, not
        NH+), including the margin, or None for an unlabelled centre."""
        label = atom_label(self.mol, centre, self.points, 1)
        if label is None:
            return None
        return _grow(label_box(label, self.points[centre], self.font, self.size, self.scale), self.margin)

    def reach(self, centre: int, u: Point) -> float:
        """Return how far the drawn label of `centre` extends along the unit vector u (0 if none)."""
        box = self.centre_box(centre)
        return 0.0 if box is None else _exit(self.points[centre], u, box)

    def beside_h(self, centre: int, u: Point) -> bool:
        """Return True if the label of `centre` with all its hydrogens (NH+) extends farther along
        u than the label as drawn (N+), that is, if u points to its H on the left or right."""
        parts = self.boxes.get(centre)
        if not parts:
            return False
        return _exit(self.points[centre], u, parts[0]) > self.reach(centre, u) + 1e-9

    def seen(self, centre: int, h: Point) -> float:
        """Return how much of the wedge from `centre` to a hydrogen at h shows between their labels."""
        c = self.points[centre]
        d = math.dist(c, h)
        u = ((h[0] - c[0]) / d, (h[1] - c[1]) / d)
        return d - self.reach(centre, u) - _exit(h, (-u[0], -u[1]), self.h_box(h))

    def h_box(self, h: Point):
        box = label_box(Label([Run("H")], 0), h, self.font, self.size, self.scale)
        return _grow(box, self.margin)

    def clear(self, centre: int, h: Point) -> bool:
        c, box = self.points[centre], self.h_box(h)
        # The second line of a double bond lies up to the gap away from the bond.
        wide = _grow(box, self.gap)
        for i, j, double in self.bonds:
            p, q = self.points[i], self.points[j]
            if segment_in_box(p, q, wide if double else box):
                return False
            if centre not in (i, j) and segments_cross(c, h, p, q):
                return False
        for a, parts in self.boxes.items():
            # An H above or below the symbol of the centre is the one its mark would draw.
            for b in (parts[:1] if a == centre else parts):
                if boxes_meet(box, b) or (a != centre and segment_in_box(c, h, b)):
                    return False
        for o, oh in self.drawn:
            ob = self.h_box(oh)
            if (boxes_meet(box, ob) or segment_in_box(c, h, ob) or segment_in_box(self.points[o], oh, box)
                    or segments_cross(c, h, self.points[o], oh)):
                return False
        return True


def _exit(p: Point, u: Point, box) -> float:
    """Return the distance from p, a point in the box, to the edge of the box along the unit vector u."""
    t = math.inf
    for k, (low, high) in enumerate(((box[0], box[2]), (box[1], box[3]))):
        if u[k] > 1e-12:
            t = min(t, (high - p[k]) / u[k])
        elif u[k] < -1e-12:
            t = min(t, (low - p[k]) / u[k])
    return max(0.0, t)


def _grow(box, d: float):
    return (box[0] - d, box[1] - d, box[2] + d, box[3] + d)


def _inside(p: Point, polygon: Sequence[Point]) -> bool:
    """Return True if p is inside the polygon (even-odd rule)."""
    inside = False
    n = len(polygon)
    for k in range(n):
        (x0, y0), (x1, y1) = polygon[k], polygon[(k + 1) % n]
        if (y0 > p[1]) != (y1 > p[1]) and p[0] < x0 + (p[1] - y0) * (x1 - x0) / (y1 - y0):
            inside = not inside
    return inside


def ring_bonds(mol: Molecule) -> set:
    """Return the ring bonds, as frozensets of their two atoms."""
    out = set()
    for r in mol.rings:
        out.update(frozenset((r[k], r[(k + 1) % len(r)])) for k in range(len(r)))
    return out


def mark_kind(mol: Molecule, points: Sequence[Point], centre: int, to: Optional[int],
              hydrogen: Optional[Point] = None) -> Optional[str]:
    """Return "wedge" or "hash", whichever mark from `centre` to `to` (or to a drawn hydrogen)
    reads back as the "@" or "@@" in the SMILES, or None if neither does (bonds in a line)."""
    want = mol.atoms[centre].chirality
    for kind in ("wedge", "hash"):
        mark = Mark(centre, kind, to, hydrogen)
        if drawn_chirality(mol, points, centre, [mark]) == want:
            return kind
    return None


def drawn_chirality(mol: Molecule, points: Sequence[Point], centre: int,
                    marks: Iterable[Mark]) -> Optional[str]:
    """Return "@" or "@@" as the drawing shows the centre, for its neighbours in SMILES order
    (Atom.neighbor_order). Only marks with their narrow end on the centre count. Return None
    when the drawing does not tell (no mark, bonds in a line, or the mark on the middle of
    three bonds within half a turn)."""
    mine = [m for m in marks if m.centre == centre]
    if not mine:
        return None
    if _middle_of_a_half_turn(mol, points, centre, mine):
        return None
    vectors = _vectors(mol, points, centre, mine)
    if vectors is None:
        return None
    v0, v1, v2, v3 = vectors
    volume = _det([_sub(v1, v0), _sub(v2, v0), _sub(v3, v0)])
    if abs(volume) < FLAT:
        return None
    # Looking from the first neighbour, the other three go anticlockwise for "@".
    return "@" if volume < 0 else "@@"


def _middle_of_a_half_turn(mol, points, centre, marks) -> bool:
    """Return True for three drawn bonds that lie within about half a turn, with the only mark
    on the middle one and the hydrogen or lone pair not drawn. Such a drawing can be read
    either way (RDKit reads some of them the other way), so it is not used."""
    if len(marks) != 1 or marks[0].to is None:
        return False
    neighbours = mol.neighbors(centre)
    if len(neighbours) != 3:
        return False
    c = points[centre]
    angles = sorted((math.atan2(points[j][1] - c[1], points[j][0] - c[0]), j) for j in neighbours)
    for k in range(3):
        a, b = angles[k][0], angles[(k + 1) % 3][0]
        gap = (b - a) % (2 * math.pi)
        if gap > math.pi - HALF_TURN_SLACK:
            middle = angles[(k + 2) % 3][1]
            return marks[0].to == middle
    return False


def neighbour_order(mol: Molecule, centre: int) -> List[int]:
    """Return Atom.neighbor_order, with LONE_PAIR inserted for a centre with three neighbours.
    The lone pair comes second, even when nothing is written before the centre, as RDKit reads
    it ([S@@](=O)(C)c1ccccc1 is C[S@](=O)c1ccccc1). The hydrogen in the bracket atom of a
    centre written first comes first instead ([C@@H](F)(Cl)Br is F[C@H](Cl)Br)."""
    order = list(mol.atoms[centre].neighbor_order)
    if len(order) == 3 and IMPLICIT_H not in order:
        order.insert(1, LONE_PAIR)
    return order


def _vectors(mol, points, centre, marks) -> Optional[List[Tuple[float, float, float]]]:
    order = neighbour_order(mol, centre)
    if len(order) != 4:
        return None
    c = points[centre]
    by_atom = {m.to: m for m in marks if m.to is not None}
    drawn_h = [m for m in marks if m.to is None and m.hydrogen is not None]
    vectors: List[Optional[Tuple[float, float, float]]] = []
    for j in order:
        if j >= 0:
            m = by_atom.get(j)
            vectors.append(_unit3(points[j], c, m.z if m else 0))
        elif j == IMPLICIT_H and drawn_h:
            vectors.append(_unit3(drawn_h[0].hydrogen, c, drawn_h[0].z))
        else:
            vectors.append(None)
    missing = [k for k, v in enumerate(vectors) if v is None]
    if len(missing) > 1:
        return None
    if missing:
        # A hydrogen or lone pair that is not drawn points away from the others, out of the
        # plane on the side opposite the marks.
        sx = sum(v[0] for v in vectors if v)
        sy = sum(v[1] for v in vectors if v)
        sz = sum(v[2] for v in vectors if v)
        z = -1.0 if sz > 0 else 1.0 if sz < 0 else 0.0
        vectors[missing[0]] = (-0.5 * sx, -0.5 * sy, z)
    return vectors


def _unit3(p: Point, c: Point, z: float) -> Tuple[float, float, float]:
    dx, dy = p[0] - c[0], p[1] - c[1]
    d = math.hypot(dx, dy) or 1.0
    return (dx / d, dy / d, float(z))


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _det(m) -> float:
    (a, b, c), (d, e, f), (g, h, i) = m
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def marks_by_bond(mol: Molecule, marks: Iterable[Mark]) -> Dict[int, Mark]:
    """Return {bond index: mark} for the marks that lie on bonds of the molecule."""
    out = {}
    for m in marks:
        if m.to is None:
            continue
        for k, b in enumerate(mol.bonds):
            if {b.begin, b.end} == {m.centre, m.to}:
                out[k] = m
                break
    return out
