"""Molecule drawing (the fourth stage): molecule + layout + style -> a group of shapes.

`draw_molecule(mol, points, style)` turns the layout (bond length 1, y up) into lines,
wedges, hashed wedges and text boxes in points (y down, as on a slide), in one group.
Atom labels are in labels.py, bonds in bonds.py and stereocentre wedges in stereo.py.

The style selects between drawing options: aromatic rings with alternating double bonds or
with a circle (Style.aromatic), and everything in black, atom labels in element colours,
or labels and bond halves in element colours (Style.colour). Groups drawn as text
(Style.abbreviate) are removed from the molecule before layout, by chem/abbreviate.py.
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Tuple

from ..chem.molecule import Molecule
from ..figure import Group, Run, Text
from ..fonts import Font, get_font
from ..style import Style
from .bonds import BLACK, bond_lines, clip_ends, ring_circle, wedge_shapes
from .labels import Label, atom_label, hydrogen_label, hydrogen_point, label_box
from .stereo import Mark, hydrogen_drawn, marks_by_bond, stereo_marks

Point = Tuple[float, float]
Box = Tuple[float, float, float, float]


def draw_molecule(mol: Molecule, points: Sequence[Point], style: Optional[Style] = None,
                  name: str = "Molecule", marks: Optional[List[Mark]] = None) -> Group:
    """Draw a perceived molecule from its layout, as one group with its top left corner at
    (0, 0). `marks` are the wedges; by default stereo_marks() chooses them."""
    style = style or Style()
    font = get_font(style.font)
    L = style.bond_length
    P = [(x * L, -y * L) for x, y in points]
    if marks is None:
        marks = stereo_marks(mol, points)
    # The colour of each atom: for its label, and for its half of a bond when bonds are coloured.
    colours = [BLACK if a.label else style.element_colour(a.element) for a in mol.atoms]
    hydrogen_colour = style.element_colour("H")
    coloured_bonds = style.colour == "atoms_and_bonds"

    def ends(i: int, j: int) -> Tuple[str, str]:
        return (colours[i], colours[j]) if coloured_bonds else (BLACK, BLACK)

    texts: List[Text] = []
    boxes: Dict[int, List[Box]] = {}   # {atom: [its label box, and the box of an H above or below]}
    bare: Dict[int, List[Box]] = {}    # the same boxes without the label margin
    # A hydrogen drawn at the end of a wedge is not repeated in the label of its atom.
    own_hydrogen = hydrogen_drawn(marks)
    for i, atom in enumerate(mol.atoms):
        label = atom_label(mol, i, points, 1 if i in own_hydrogen else 0)
        if label is None:
            continue
        said = atom.label[0] if atom.label else atom.element
        text, box = _text(label, P[i], font, style, f"Atom {said} {i + 1}", colours[i])
        texts.append(text)
        bare[i] = [box]
        boxes[i] = [_widen(box, style.label_margin)]
        h = hydrogen_point(label, P[i], style.font_size, y_down=True)
        if h is not None:
            # An H above or below the symbol is a text box of its own, in the colour of its atom,
            # as when it is written beside the symbol.
            text, box = _text(hydrogen_label(label), h, font, style, f"Atom {said} {i + 1} H", colours[i])
            texts.append(text)
            bare[i].append(box)
            boxes[i].append(_widen(box, style.label_margin))
    shapes = []
    on_bond = marks_by_bond(mol, marks)
    circled = mol.aromatic_bonds() if style.aromatic == "circle" else set()
    for k, bond in enumerate(mol.bonds):
        mark = on_bond.get(k)
        if mark is None:
            near = close_margin(mol, k, P, boxes, bare, style, single=k in circled)
            shapes += bond_lines(mol, k, P, near, style, single=k in circled,
                                 colours=ends(bond.begin, bond.end))
            continue
        centre, far = P[mark.centre], P[mark.to]
        seg = clip_ends((centre, far), (mark.centre, mark.to), boxes)
        if seg is not None:
            shapes += wedge_shapes(k + 1, centre, far, mark.kind, seg, style, ends(mark.centre, mark.to))
    if style.aromatic == "circle":
        shapes += [ring_circle(mol.rings[r], P, boxes, style, number + 1)
                   for number, r in enumerate(mol.aromatic_rings)]
    n = len(mol.atoms)
    for m in (m for m in marks if m.to is None and m.hydrogen is not None):
        n += 1
        h = (m.hydrogen[0] * L, -m.hydrogen[1] * L)
        label = Label([Run("H")], 0)
        text, box = _text(label, h, font, style, f"Atom H {n}", hydrogen_colour)
        texts.append(text)
        seg = clip_ends((P[m.centre], h), (m.centre, -1), {**boxes, -1: [_widen(box, style.label_margin)]})
        if seg is not None:
            both = (colours[m.centre], hydrogen_colour) if coloured_bonds else (BLACK, BLACK)
            shapes += wedge_shapes(len(mol.bonds) + n - len(mol.atoms), P[m.centre], h, m.kind, seg, style, both)
    group = Group(shapes + texts, name=name)
    x0, y0, _, _ = group.bbox()
    return group.moved(-x0, -y0)


# The smallest label margin close_margin gives the ends of a bond, as a fraction of Style.label_margin.
LEAST_MARGIN = 0.25


def close_margin(mol: Molecule, k: int, P: Sequence[Point], boxes: Dict[int, List[Box]],
                 bare: Dict[int, List[Box]], style: Style, single: bool = False) -> Dict[int, List[Box]]:
    """Return the boxes to clip bond k with. A bond between two labels that would show shorter
    than style.shortest_line gets a narrower margin at both ends, just enough to show that
    length, but no less than LEAST_MARGIN times the label margin. These are ring bonds (the
    N-N of a benzotriazole or a tetrazole, the O-O of a peroxide bridge): when such a bond is
    in no ring, the layout lengthens it instead (layout/lengthen.py). Every other bond keeps
    `boxes`."""
    bond = mol.bonds[k]
    i, j = bond.begin, bond.end
    if style.shortest_line <= 0 or i not in bare or j not in bare:
        return boxes
    target = style.shortest_line * style.bond_length

    def with_margin(m: float) -> Dict[int, List[Box]]:
        return {**boxes, i: [_widen(b, m) for b in bare[i]], j: [_widen(b, m) for b in bare[j]]}

    def seen(m: float) -> float:
        lines = bond_lines(mol, k, P, with_margin(m), style, single=single)
        return max((math.dist(line.p0, line.p1) for line in lines), default=0.0)

    high = style.label_margin
    if seen(high) >= target:
        return boxes
    low = LEAST_MARGIN * style.label_margin
    if seen(low) < target:
        return with_margin(low)
    for _ in range(30):
        middle = (low + high) / 2
        low, high = (middle, high) if seen(middle) >= target else (low, middle)
    return with_margin(low)


def _text(label: Label, point: Point, font: Font, style: Style, name: str, colour: str) -> Tuple[Text, Box]:
    size = style.font_size
    box = label_box(label, point, font, size, style.script_scale)
    baseline = point[1] + font.capital_height(size) / 2
    top = baseline - font.ascent_height(size)
    text = Text(box[0], top, box[2] - box[0], font.line_height(size), list(label.runs), size,
                font=style.font, colour=colour, name=name)
    return text, box


def _widen(box: Box, d: float) -> Box:
    return (box[0] - d, box[1] - d, box[2] + d, box[3] + d)


__all__ = ["draw_molecule", "Mark", "stereo_marks"]
