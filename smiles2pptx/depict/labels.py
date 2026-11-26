"""Atom labels: which atoms get one, what they say, and their bounding boxes.

Rules:

- carbon has no label, except a carbon with a charge or an isotope, a carbon with two
  double bonds (CO2), and a molecule of one atom
- the hydrogens of a labelled atom are written in its label (OH, NH2). When its only bond
  comes from the right, the hydrogens go on the left (HO)
- the hydrogens of an atom with two or more bonds go to the right of, left of, below or
  above the element symbol, whichever of the four is farthest from the nearest bond (ties
  are broken in the order right, left, below, above). An NH on the left side of a ring is
  HN, with its hydrogen outside the ring; an NH at the bottom corner of a ring is N with
  the H below it. Hydrogens above or below are a text box of their own, directly over or
  under the first letter of the symbol (style.HYDROGEN_LINE), and the charge stays on the
  right of the symbol (N+)
- a hydrogen drawn as a label of its own, at the wide end of a stereocentre wedge
  (stereo.py), is left out of the label of its atom: an NH+ at a ring fusion is labelled N+,
  and its one H is drawn at the wedge
- a molecule of one atom has its hydrogens in front for O, S, Se, Te and the halogens (H2O,
  HCl) and after the symbol for the others (NH3, CH4)
- a charge is a superscript on the right (+, 2+, the minus sign); an isotope is a
  superscript in front (13C)
- an atom that stands for a group (Atom.label, from chem/abbreviate.py) shows the text of
  the group: the one for a bond from the left (CO2H) or, when its bond comes from the
  right, the one for a bond from the right (HO2C)

The first letter of the element symbol is centred on the atom, so a label extends towards
the side its hydrogens are on. The text of a group has its first letter on the atom, or its
last letter when the bond comes from the right. The box is as wide as the text and as tall
as a capital letter. A label with its hydrogens above or below has two boxes (label_parts):
its bond lines stop at both, and overlaps are counted against both.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

from ..chem.molecule import Molecule
from ..figure import Run
from ..fonts import Font
from ..style import HYDROGEN_LINE

MINUS = "−"
# Elements whose hydrogens come first in a molecule of one atom (H2O, HCl).
H_IN_FRONT = {"O", "S", "Se", "Te", "F", "Cl", "Br", "I"}
# A neighbour this far to the right of the atom (in bond lengths) counts as on the right.
RIGHT = 1e-3
# The positions for the hydrogens of an atom with two or more bonds, in tie-break order, with
# their directions in the layout (radians, y up).
H_PLACES = (("right", 0.0), ("left", math.pi), ("below", -math.pi / 2), ("above", math.pi / 2))
# Two positions tie when their angles to the nearest bond differ by less than this.
SAME_ANGLE = math.radians(1)

Box = Tuple[float, float, float, float]  # x0, y0, x1, y1


@dataclass(frozen=True)
class Label:
    runs: List[Run]
    symbol_run: int   # the run that holds the element symbol; its first letter sits on the atom
    at_end: bool = False   # the last letter of that run sits on the atom instead (HO2C)
    # The hydrogens written above or below the symbol, as a separate label (H, H2), and their
    # side: 1 above the symbol and -1 below it, as seen on the page.
    h_runs: Optional[List[Run]] = None
    h_side: int = 0


def formula_runs(text: str) -> List[Run]:
    """Split text into runs, with each digit that follows a Latin letter (or such a digit)
    as a subscript: "H2SO4" -> H, 2, SO, 4."""
    runs: List[Run] = []
    prev, prev_kind = "", "normal"
    for ch in text:
        sub = ch.isdigit() and ((prev.isascii() and prev.isalpha()) or (prev.isdigit() and prev_kind == "sub"))
        kind = "sub" if sub else "normal"
        if runs and runs[-1].kind == kind:
            runs[-1] = Run(runs[-1].text + ch, kind)
        else:
            runs.append(Run(ch, kind))
        prev, prev_kind = ch, kind
    return runs


def atom_label(mol: Molecule, atom: int, points: Sequence[Tuple[float, float]],
               drawn_hydrogens: int = 0) -> Optional[Label]:
    """Return the label of one atom of a perceived molecule, or None if it has no label.
    `points` are the atom coordinates (any unit; only the directions to the neighbours
    matter). `drawn_hydrogens` of its hydrogens are drawn as labels of their own (at the end
    of a wedge) and are left out of this label."""
    a = mol.atoms[atom]
    neighbours = mol.neighbors(atom)
    from_right = len(neighbours) == 1 and points[neighbours[0]][0] - points[atom][0] > RIGHT
    if a.label is not None:
        runs = formula_runs(a.label[1] if from_right else a.label[0])
        normal = [n for n, run in enumerate(runs) if run.kind == "normal"]
        return Label(runs, normal[-1] if from_right else normal[0], at_end=from_right)
    if a.element == "C" and not _carbon_has_text(mol, atom, neighbours):
        return None
    h = max(0, (a.h_count or 0) - drawn_hydrogens)
    hydrogens = [Run("H")] + ([Run(str(h), "sub")] if h > 1 else []) if h else []
    front = [Run(str(a.isotope), "sup")] if a.isotope else []
    charge = []
    if a.charge:
        charge = [Run((str(abs(a.charge)) if abs(a.charge) > 1 else "") + ("+" if a.charge > 0 else MINUS), "sup")]
    symbol = Run(a.element)
    if len(neighbours) > 1:
        place = hydrogen_side(points, atom, neighbours)
    else:
        place = "left" if (from_right if neighbours else a.element in H_IN_FRONT) else "right"
    if hydrogens and place in ("above", "below"):
        return Label(front + [symbol] + charge, len(front), h_runs=hydrogens, h_side=1 if place == "above" else -1)
    if hydrogens and place == "left":
        runs = hydrogens + front + [symbol] + charge
        return Label(runs, len(hydrogens) + len(front))
    runs = front + [symbol] + hydrogens + charge
    return Label(runs, len(front))


def hydrogen_side(points, atom: int, neighbours: Sequence[int]) -> str:
    """Return where the hydrogens of an atom with several bonds go: "right", "left", "below"
    or "above", whichever is farthest from the nearest bond (H_PLACES). `points` are in the
    layout (y up). Return "right" when all the neighbours are at the position of the atom."""
    x, y = points[atom]
    bonds = [math.atan2(points[n][1] - y, points[n][0] - x) for n in neighbours
             if points[n][0] != x or points[n][1] != y]
    if not bonds:
        return "right"
    best, farthest = "right", -1.0
    for place, angle in H_PLACES:
        clear = min(abs(math.remainder(angle - b, 2 * math.pi)) for b in bonds)
        if clear > farthest + SAME_ANGLE:
            best, farthest = place, clear
    return best


def _carbon_has_text(mol: Molecule, atom: int, neighbours: List[int]) -> bool:
    a = mol.atoms[atom]
    if a.charge or a.isotope or not neighbours:
        return True
    doubles = sum(1 for b in mol.bonds if atom in (b.begin, b.end) and b.order == 2)
    return doubles >= 2


def label_box(label: Label, point: Tuple[float, float], font: Font, size: float,
              script_scale: float) -> Box:
    """Return the bounding box of the label, with the first letter of the symbol centred on
    `point` (the last letter if `label.at_end` is set). `size` is the font size in the units
    of `point`."""
    width = font.runs_width(label.runs, size, script_scale)
    half_height = font.capital_height(size) / 2
    symbol = label.runs[label.symbol_run].text
    if label.at_end:
        upto = font.runs_width(label.runs[:label.symbol_run + 1], size, script_scale)
        left = point[0] - upto + font.char_width(symbol[-1]) * size / 2
    else:
        before = font.runs_width(label.runs[:label.symbol_run], size, script_scale)
        left = point[0] - before - font.char_width(symbol[0]) * size / 2
    return (left, point[1] - half_height, left + width, point[1] + half_height)


def hydrogen_label(label: Label) -> Optional[Label]:
    """Return the hydrogens written above or below the symbol as a label of their own (its
    first letter, the H, sits on hydrogen_point), or None if they are not written there."""
    return Label(list(label.h_runs), 0) if label.h_runs else None


def hydrogen_point(label: Label, point: Tuple[float, float], size: float,
                   y_down: bool) -> Optional[Tuple[float, float]]:
    """Return the centre of the H written above or below the symbol of a label at `point`,
    or None. `y_down` means that y points down in the coordinates of `point` (points on a
    slide) rather than up (the layout)."""
    if not label.h_runs:
        return None
    up = label.h_side * HYDROGEN_LINE * size
    return (point[0], point[1] - up if y_down else point[1] + up)


def label_parts(label: Label, point: Tuple[float, float], font: Font, size: float, script_scale: float,
                y_down: bool = False) -> List[Box]:
    """Return the bounding boxes of the label: the one from label_box, then the box of the
    hydrogens when they are written above or below the symbol. `y_down` is as in
    hydrogen_point."""
    boxes = [label_box(label, point, font, size, script_scale)]
    h = hydrogen_point(label, point, size, y_down)
    if h is not None:
        boxes.append(label_box(hydrogen_label(label), h, font, size, script_scale))
    return boxes


def label_boxes(mol: Molecule, points: Sequence[Tuple[float, float]], font: Font, size: float,
                script_scale: float) -> dict:
    """Return {atom: [boxes]} (label_parts) for every labelled atom, with all its hydrogens in
    the label. `points` are in the layout (y up)."""
    boxes = {}
    for i in range(len(mol.atoms)):
        label = atom_label(mol, i, points)
        if label is not None:
            boxes[i] = label_parts(label, points[i], font, size, script_scale)
    return boxes
