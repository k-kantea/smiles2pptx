"""Tests for the overlap counts (smiles2pptx/check.py) and atom labels (depict/labels.py).

Each count is checked on small hand-made drawings, with one case that has the overlap
and one that just misses it. The label tests cover which atoms get a label, what it says
and where its hydrogens go.
"""
import math

import pytest

from smiles2pptx.check import (CLASH, ON_BOND, boxes_meet, label_overlaps, overlaps, point_overlaps,
                               segment_in_box, segments_cross)
from smiles2pptx.chem.perceive import perceive
from smiles2pptx.chem.smiles import read_smiles
from smiles2pptx.depict.labels import MINUS, atom_label, label_box, label_boxes, label_parts
from smiles2pptx.fonts import get_font
from smiles2pptx.layout import layout
from smiles2pptx.style import HYDROGEN_LINE, Style

ARIAL = get_font("Arial")


# --- atoms as points --------------------------------------------------------------------

def test_two_atoms_closer_than_half_a_bond_clash():
    bonds = [(0, 1), (2, 3)]
    near = [(0, 0), (1, 0), (0, CLASH - 0.01), (0, CLASH + 0.99)]
    apart = [(0, 0), (1, 0), (0, CLASH + 0.01), (0, CLASH + 1.01)]
    assert point_overlaps(near, bonds).clash == [(0, 2)]
    assert point_overlaps(apart, bonds).clash == []


def test_bonded_atoms_never_clash():
    assert point_overlaps([(0, 0), (1, 0)], [(0, 1)]).points == 0


def test_an_atom_close_to_the_middle_of_a_bond_is_on_the_bond():
    bonds = [(0, 1)]
    on = [(0, 0), (1, 0), (0.5, ON_BOND - 0.01)]
    off = [(0, 0), (1, 0), (0.5, ON_BOND + 0.01)]
    assert point_overlaps(on, bonds).on_bond == [(2, 0)]
    assert point_overlaps(off, bonds).on_bond == []


def test_an_atom_close_to_the_end_of_a_bond_clashes_and_is_not_also_on_the_bond():
    r = point_overlaps([(0, 0), (1, 0), (1.1, 0.1)], [(0, 1)])
    assert r.clash == [(1, 2)] and r.on_bond == []


def test_two_bonds_that_cross_are_counted_and_bonds_that_meet_at_an_atom_are_not():
    x = point_overlaps([(0, 0), (1, 1), (0, 1), (1, 0)], [(0, 1), (2, 3)])
    assert x.crossing == [(0, 1)]
    v = point_overlaps([(0, 0), (1, 0), (1, 1)], [(0, 1), (1, 2)])
    assert v.crossing == []
    assert not segments_cross((0, 0), (1, 0), (2, 1), (2, -1))   # would meet only if longer


def test_the_counts_do_not_depend_on_the_size_of_the_drawing():
    points = [(0, 0), (1, 1), (0, 1), (1, 0), (0.3, 0.55), (2, 2)]
    bonds = [(0, 1), (2, 3), (4, 5)]
    small = point_overlaps(points, bonds)
    big = point_overlaps([(20 * x, 20 * y) for x, y in points], bonds)
    assert (small.clash, small.on_bond, small.crossing) == (big.clash, big.on_bond, big.crossing)
    assert small.points > 0


# --- label bounding boxes -----------------------------------------------------------------

def test_a_bond_through_the_text_of_an_atom_it_is_not_on_is_counted():
    points = [(0, 0), (2, 0), (1, 0)]
    boxes = {2: (0.8, -0.2, 1.2, 0.2)}
    assert label_overlaps(points, [(0, 1)], boxes).line_label == [(2, 0)]
    assert label_overlaps(points, [(0, 2)], boxes).line_label == []          # the bond ends on that atom
    assert label_overlaps(points, [(0, 1)], {2: (0.8, 0.1, 1.2, 0.5)}).line_label == []


def test_two_texts_whose_boxes_meet_are_counted():
    boxes = {0: (0, 0, 1, 1), 1: (0.9, 0.5, 2, 1.5), 2: (1.1, 0, 2, 0.6)}
    assert label_overlaps([(0, 0)] * 3, [], boxes).label_label == [(0, 1), (1, 2)]
    assert boxes_meet((0, 0, 1, 1), (0.5, 0.5, 2, 2)) and not boxes_meet((0, 0, 1, 1), (1, 0, 2, 1))


def test_a_segment_that_passes_by_a_box_is_not_inside_it():
    assert not segment_in_box((0, 1.01), (2, 1.01), (0.5, 0, 1.5, 1))
    assert segment_in_box((0, 0.99), (2, 0.99), (0.5, 0, 1.5, 1))
    assert not segment_in_box((0, 0), (0.4, 0.5), (0.5, 0, 1.5, 1))   # stops short of it


# --- atom labels --------------------------------------------------------------------------

def label_text(smiles, atom, points=None):
    """Return the runs of an atom's label. Hydrogens written above or below the symbol come last,
    as ("above", runs) or ("below", runs)."""
    mol = perceive(read_smiles(smiles))
    label = atom_label(mol, atom, points or layout(mol))
    if label is None:
        return None
    runs = [(run.text, run.kind) for run in label.runs]
    if label.h_runs:
        runs.append(("above" if label.h_side == 1 else "below", [(run.text, run.kind) for run in label.h_runs]))
    return runs


def test_carbon_has_no_text_but_other_atoms_do():
    assert label_text("CCO", 0) is None and label_text("CCO", 1) is None
    assert label_text("CCN", 2) in ([("N", "normal"), ("H", "normal"), ("2", "sub")],
                                    [("H", "normal"), ("2", "sub"), ("N", "normal")])


def test_the_hydrogens_go_on_the_side_away_from_the_bond():
    assert label_text("CO", 1, [(0, 0), (1, 0)]) == [("O", "normal"), ("H", "normal")]
    assert label_text("CO", 1, [(1, 0), (0, 0)]) == [("H", "normal"), ("O", "normal")]
    # Two bonds in a straight line take right and left, so the hydrogen goes below.
    assert label_text("CNC", 1, [(1, 0), (0, 0), (-1, 0)]) == [("N", "normal"), ("below", [("H", "normal")])]


def test_the_hydrogens_of_an_atom_with_several_bonds_go_away_from_them():
    nh, hn = [("N", "normal"), ("H", "normal")], [("H", "normal"), ("N", "normal")]
    above, below = ("above", [("H", "normal")]), ("below", [("H", "normal")])
    n = ("N", "normal")
    # Both bonds to the right (the left corner of a ring): HN.
    assert label_text("CNC", 1, [(0.5, 0.87), (0, 0), (0.5, -0.87)]) == hn
    # One bond straight down and one to the upper right (the upper left corner of a ring): HN.
    assert label_text("CNC", 1, [(0, -1), (0, 0), (0.87, 0.5)]) == hn
    # Both bonds to the left: NH.
    assert label_text("CNC", 1, [(-0.5, 0.87), (0, 0), (-0.5, -0.87)]) == nh
    # One bond to each side, both pointing down (a chain, or the top corner of a ring): H above.
    assert label_text("CNC", 1, [(-0.87, -0.5), (0, 0), (0.87, -0.5)]) == [n, above]
    # Both pointing up (the bottom corner of a ring): H below.
    assert label_text("CNC", 1, [(-0.87, 0.5), (0, 0), (0.87, 0.5)]) == [n, below]
    # The bottom right corner of a ring (one bond to the left, one to the upper right): below,
    # which is 90 degrees from the nearest bond, whereas right is only 60 degrees from it.
    assert label_text("CNC", 1, [(-1, 0), (0, 0), (0.5, 0.87)]) == [n, below]
    # The size of the drawing does not matter.
    assert label_text("CNC", 1, [(50, 87), (0, 0), (50, -87)]) == hn
    assert label_text("CNC", 1, [(0.005, 0.0087), (0, 0), (0.005, -0.0087)]) == hn
    assert label_text("CNC", 1, [(-87, -50), (0, 0), (87, -50)]) == [n, above]
    # Three bonds and a charge: the hydrogen still goes away from the bonds, and the charge
    # stays to the right of the symbol.
    plus = ("+", "sup")
    assert label_text("C[NH+](C)C", 1, [(0, 1), (0, 0), (0.87, 0.5), (0.87, -0.5)]) == hn + [plus]
    assert label_text("C[NH+](C)C", 1, [(0, 1), (0, 0), (0.87, -0.5), (-0.87, -0.5)]) == [n, plus, below]
    assert label_text("C[NH+](C)C", 1, [(0, -1), (0, 0), (0.87, 0.5), (-0.87, 0.5)]) == [n, plus, above]
    assert label_text("C[NH+](C)C", 1, [(-1, 0), (0, 0), (0.5, 0.87), (0.5, -0.87)]) == nh + [plus]
    assert label_text("C[NH+](C)C", 1, [(1, 0), (0, 0), (-0.5, 0.87), (-0.5, -0.87)]) == hn + [plus]
    assert label_text("C[SiH](C)C", 1, [(0, 1), (0, 0), (0.87, -0.5), (-0.87, -0.5)]) == [("Si", "normal"), below]
    # Two hydrogens above or below are written together: H2.
    assert label_text("C[NH2+]C", 1, [(-0.87, -0.5), (0, 0), (0.87, -0.5)]) == [
        n, plus, ("above", [("H", "normal"), ("2", "sub")])]
    # Atoms at the same position (before layout) do not break it.
    assert label_text("CNC", 1, [(0, 0), (0, 0), (0, 0)]) == nh


def test_when_two_places_are_as_far_from_the_bonds_right_left_below_above_is_the_order():
    n, h = ("N", "normal"), [("H", "normal")]
    # Bonds up and down: right and left are both 90 degrees away; right first.
    assert label_text("CNC", 1, [(0, 1), (0, 0), (0, -1)]) == [n, ("H", "normal")]
    # Bonds right and left: below and above; below first.
    assert label_text("CNC", 1, [(1, 0), (0, 0), (-1, 0)]) == [n, ("below", h)]
    # Bonds to the upper right and lower left: all four are 45 degrees away; right first.
    assert label_text("CNC", 1, [(1, 1), (0, 0), (-1, -1)]) == [n, ("H", "normal")]
    # Bonds to the upper left and upper right of a square: left and right 45 degrees, below 135.
    assert label_text("CNC", 1, [(-1, 1), (0, 0), (1, 1)]) == [n, ("below", h)]
    # Half a degree counts as a tie, so the order decides; two degrees do not.
    tilt = math.radians(0.5)
    assert label_text("CNC", 1, [(math.cos(tilt), math.sin(tilt)), (0, 0), (-1, 0)]) == [n, ("below", h)]
    tilt = math.radians(2)
    assert label_text("CNC", 1, [(math.cos(tilt), math.sin(tilt)), (0, 0), (-1, 0)]) == [n, ("below", h)]
    assert label_text("CNC", 1, [(math.cos(tilt), -math.sin(tilt)), (0, 0), (-1, 0)]) == [n, ("above", h)]


@pytest.mark.parametrize("smiles,runs", [
    ("O", [("H", "normal"), ("2", "sub"), ("O", "normal")]),
    ("Cl", [("H", "normal"), ("Cl", "normal")]),
    ("S", [("H", "normal"), ("2", "sub"), ("S", "normal")]),
    ("N", [("N", "normal"), ("H", "normal"), ("3", "sub")]),
    ("C", [("C", "normal"), ("H", "normal"), ("4", "sub")]),
    ("[13CH4]", [("13", "sup"), ("C", "normal"), ("H", "normal"), ("4", "sub")]),
    ("[Na+]", [("Na", "normal"), ("+", "sup")]),
    ("[O-2]", [("O", "normal"), ("2" + MINUS, "sup")]),
])
def test_a_molecule_of_one_atom(smiles, runs):
    assert label_text(smiles, 0) == runs


def test_carbon_with_a_charge_an_isotope_or_two_double_bonds_has_text():
    assert label_text("O=C=O", 1) == [("C", "normal")]
    left = [(-0.5, 0.87), (0, 0), (-0.5, -0.87)]   # both bonds to the left: the H on the right
    assert label_text("C[CH-]C", 1, left) == [("C", "normal"), ("H", "normal"), (MINUS, "sup")]
    # In the layout the C is a lower corner of the zigzag, so the H goes below it.
    assert label_text("C[CH-]C", 1) == [("C", "normal"), (MINUS, "sup"), ("below", [("H", "normal")])]
    assert label_text("C[13CH2]C", 1)[0] == ("13", "sup")
    assert label_text("C[N+](C)(C)C", 1) == [("N", "normal"), ("+", "sup")]


def test_a_hydrogen_drawn_as_a_text_of_its_own_is_left_out_of_the_text_of_its_atom():
    def text(smiles, atom, drawn, points=((-0.87, -0.5), (0, 0), (0.87, -0.5), (0, 1))):
        mol = perceive(read_smiles(smiles))
        label = atom_label(mol, atom, list(points), drawn)
        if label is None:
            return None
        runs = [(run.text, run.kind) for run in label.runs]
        return runs + ([("below" if label.h_side == -1 else "above", "".join(r.text for r in label.h_runs))]
                       if label.h_runs else [])

    n, h, plus = ("N", "normal"), ("H", "normal"), ("+", "sup")
    # Three bonds, one pointing up: the hydrogen goes below the N+, or is left out when drawn separately.
    assert text("C[NH+](C)C", 1, 0) == [n, plus, ("below", "H")]
    assert text("C[NH+](C)C", 1, 1) == [n, plus]
    assert text("C[SiH](C)C", 1, 1) == [("Si", "normal")]
    assert text("C[BH-](C)C", 1, 1) == [("B", "normal"), (MINUS, "sup")]
    # Beside the symbol: the same.
    left = ((-0.87, -0.5), (0, 0), (-0.87, 0.5), (0, 1))   # all three bonds on the left: NH+
    assert text("C[NH+](C)C", 1, 0, left) == [n, h, plus]
    assert text("C[NH+](C)C", 1, 1, left) == [n, plus]
    # One of two drawn: the other stays in the label. More than the atom has: none is left.
    assert text("C[NH2+]C", 1, 1) == [n, plus, ("above", "H")]
    assert text("C[NH2+]C", 1, 1, left) == [n, h, plus]
    assert text("C[NH+](C)C", 1, 2) == [n, plus]
    # A carbon has no label, with or without its hydrogen.
    assert text("CC(C)C", 1, 0) is None and text("CC(C)C", 1, 1) is None


def test_the_first_letter_of_the_symbol_is_centred_on_the_atom():
    mol = perceive(read_smiles("CO"))
    o = ARIAL.char_width("O") * 1.0
    h = ARIAL.char_width("H") * 1.0
    right = label_box(atom_label(mol, 1, [(0, 0), (1, 0)]), (1, 0), ARIAL, 1.0, 0.66)
    left = label_box(atom_label(mol, 1, [(1, 0), (0, 0)]), (0, 0), ARIAL, 1.0, 0.66)
    assert right[0] == pytest.approx(1 - o / 2) and right[2] == pytest.approx(1 + o / 2 + h)
    assert left[0] == pytest.approx(-o / 2 - h) and left[2] == pytest.approx(o / 2)
    assert right[3] - right[1] == pytest.approx(ARIAL.capital_height(1.0))


def test_an_h_above_or_below_is_a_box_of_its_own_right_over_or_under_the_first_letter():
    mol = perceive(read_smiles("CNC"))
    label = atom_label(mol, 1, [(-0.87, 0.5), (0, 0), (0.87, 0.5)])   # both bonds above: the H below
    assert label.h_side == -1
    n, h = label_parts(label, (0, 0), ARIAL, 1.0, 0.66)                 # y up, as in the layout
    _, h_down = label_parts(label, (0, 0), ARIAL, 1.0, 0.66, y_down=True)  # y down, as on a slide
    assert n == label_box(label, (0, 0), ARIAL, 1.0, 0.66)
    assert (h[0] + h[2]) / 2 == pytest.approx(0)
    assert (h[1] + h[3]) / 2 == pytest.approx(-HYDROGEN_LINE) and h[3] - h[1] == pytest.approx(ARIAL.capital_height(1.0))
    assert h_down == pytest.approx((h[0], -h[3], h[2], -h[1]))
    assert not boxes_meet(n, h)
    # Beside the symbol, or no hydrogens: one box.
    assert len(label_parts(atom_label(mol, 1, [(-0.5, 0.87), (0, 0), (-0.5, -0.87)]), (0, 0), ARIAL, 1.0, 0.66)) == 1
    boxes = label_boxes(perceive(read_smiles("CN(C)C")), [(-0.87, 0.5), (0, 0), (0.87, 0.5), (0, -1)], ARIAL, 1.0, 0.66)
    assert [len(parts) for parts in boxes.values()] == [1]


def test_a_text_in_parts_counts_once_for_a_bond_or_a_text():
    # Atom 2 has its H below it: (0.8, -0.6, 1.2, -0.3).
    parts = {2: [(0.8, -0.2, 1.2, 0.2), (0.8, -0.6, 1.2, -0.3)]}
    points = [(0, -0.45), (2, -0.45), (1, 0)]
    assert label_overlaps(points, [(0, 1)], parts).line_label == [(2, 0)]         # through the H only
    assert label_overlaps([(1, -1), (1, 1), (1, 0)], [(0, 1)], parts).line_label == [(2, 0)]   # through both
    assert label_overlaps([(0, -0.25), (2, -0.25), (1, 0)], [(0, 1)], parts).line_label == []  # between them
    # A label that meets the H of another label.
    both = {**parts, 3: [(1.1, -0.7, 1.5, -0.5)]}
    assert label_overlaps([(0, 0)] * 4, [], both).label_label == [(2, 3)]


def test_overlaps_of_a_molecule_use_the_text_size_of_the_style():
    mol = perceive(read_smiles("OCCO"))
    # The two oxygens are close together: their labels meet, but as points the atoms do not clash.
    points = [(0.0, 0.0), (0.0, 1.0), (0.6, 1.3), (0.6, 0.3)]
    r = overlaps(mol, points)
    assert r.label_label == [(0, 3)] and r.clash == []
    small = overlaps(mol, points, Style(scale=0.5)).label_label
    assert small == [(0, 3)]   # the label size and the bond length scale together
    boxes = label_boxes(mol, points, ARIAL, Style().font_size / Style().bond_length, 0.66)
    assert set(boxes) == {0, 3}
