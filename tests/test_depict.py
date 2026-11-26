"""Tests for drawing a molecule (smiles2pptx/depict/): labels, bonds, double bonds, wedges.

- wedges: for every stereocentre of 23 molecules and every bond from it (and a hydrogen drawn
  in the widest gap), the wedge or hashed wedge is the one RDKit reads as the SMILES specifies
  (tests/data/stereo_reference.json, made by tools/make_stereo_reference.py)
- every stereocentre of the 15 molecules and the other stereo molecules gets one mark, no bond
  carries more than one, and the drawing reads back as the SMILES specifies
- mark placement: a substituent first (the end of a branch first), then a hydrogen drawn at a
  ring fusion (away from the shared bond, or inside a ring when that position is taken), then a
  bond to another stereocentre, then a ring bond; drawn hydrogens meet no line or label
- a hydrogen drawn for a mark is not repeated in the label of its atom (N+ at a ring fusion)
- the H of an NH goes where it is farthest from the bonds: beside the N, or above or below it
  as a separate label, where the lines stop
- the 15 molecules: no line through a label, no labels that meet, the group starts at (0, 0)
- a double bond in a ring has its second line inside the ring; C=O is centred; a chain
  double bond has its second line on the side of its neighbours
- lines stop short of the labels by the label margin; labels sit with their capitals on the atom
- charges, shape names, wedge and hashed wedge shapes
"""
import json
import math
from pathlib import Path

import pytest

from smiles2pptx.check import drawing_overlaps, outside, segment_in_box
from smiles2pptx.chem.perceive import perceive
from smiles2pptx.chem.smiles import read_smiles
from smiles2pptx.depict import draw_molecule
from smiles2pptx.depict.bonds import clip_ends, clip_out, double_bond_side
from smiles2pptx.depict.labels import Label, label_box
from smiles2pptx.depict.stereo import (INSIDE_LENGTHS, OUTSIDE_LENGTHS, Mark, _inside, _middle_of_a_half_turn,
                                       _Scene, drawn_chirality, hydrogen_place, mark_kind, ring_bonds,
                                       stereo_marks)
from smiles2pptx.figure import Hash, Line, Run, Text, Wedge, walk
from smiles2pptx.fonts import get_font
from smiles2pptx.layout import layout
from smiles2pptx.style import HYDROGEN_LINE, Style

DATA = Path(__file__).resolve().parent / "data"
MOLECULES = json.loads((DATA / "test_molecules.json").read_text(encoding="utf-8"))
TARGETS = {m["id"]: m["smiles"] for m in MOLECULES if m["target"]}
STEREO = json.loads((DATA / "stereo_reference.json").read_text(encoding="utf-8"))
STYLE = Style()
FONT = get_font(STYLE.font)


def drawn(smiles, style=STYLE):
    mol = perceive(read_smiles(smiles))
    points = layout(mol, style)
    return mol, points, draw_molecule(mol, points, style)


def texts_of(group):
    return [s for s in walk([group]) if isinstance(s, Text)]


def lines_named(group, name):
    return [s for s in walk([group]) if s.name == name]


# --- wedges ------------------------------------------------------------------------------

@pytest.mark.parametrize("row", STEREO, ids=[r["smiles"][:40] for r in STEREO])
def test_wedge_or_hash_is_the_one_rdkit_reads_as_the_smiles(row):
    mol = perceive(read_smiles(row["smiles"]))
    points = [tuple(p) for p in row["points"]]
    for e in row["marks"]:
        hydrogen = tuple(e["hydrogen"]) if e["hydrogen"] else None
        assert e["kind"] in ("wedge", "hash")
        got = mark_kind(mol, points, e["centre"], e["to"], hydrogen)
        if got is None:   # refused: the middle one of three bonds within half a turn
            assert _middle_of_a_half_turn(mol, points, e["centre"], [Mark(e["centre"], "wedge", e["to"])]), e
        else:
            assert got == e["kind"], e


def test_only_4_of_the_places_are_refused_as_the_middle_of_a_half_turn():
    refused = []
    for row in STEREO:
        mol = perceive(read_smiles(row["smiles"]))
        points = [tuple(p) for p in row["points"]]
        for e in row["marks"]:
            hydrogen = tuple(e["hydrogen"]) if e["hydrogen"] else None
            if mark_kind(mol, points, e["centre"], e["to"], hydrogen) is None:
                refused.append((row["smiles"][:12], e["centre"], e["to"]))
    # RDKit reads one of them (strychnine, atom 15 to 14) the opposite way to this package.
    assert len(refused) == 4, refused
    assert ("O=C1C[C@H]2O", 15, 14) in refused


def test_the_middle_of_a_t_is_refused_and_an_end_of_it_is_not():
    mol = perceive(read_smiles("F[C@H](Cl)Br"))
    points = [(-1, 0), (0, 0), (0, 1), (1, 0)]   # F - C - Br in a line, Cl the stem
    assert mark_kind(mol, points, 1, 2) is None
    assert mark_kind(mol, points, 1, 0) is not None
    points = [(-1, 0.2), (0, 0), (0, 1), (1, 0.2)]   # 158 degrees: still within half a turn
    assert mark_kind(mol, points, 1, 2) is None
    points = [(-1, -0.3), (0, 0), (0, 1), (1, -0.3)]   # 213 degrees between F and Br
    assert mark_kind(mol, points, 1, 2) is not None


@pytest.mark.parametrize("row", STEREO, ids=[r["smiles"][:40] for r in STEREO])
def test_the_other_kind_reads_as_the_other_sense(row):
    mol = perceive(read_smiles(row["smiles"]))
    points = [tuple(p) for p in row["points"]]
    for e in row["marks"]:
        hydrogen = tuple(e["hydrogen"]) if e["hydrogen"] else None
        other = "hash" if e["kind"] == "wedge" else "wedge"
        if mark_kind(mol, points, e["centre"], e["to"], hydrogen) is None:
            continue
        got = drawn_chirality(mol, points, e["centre"], [Mark(e["centre"], other, e["to"], hydrogen)])
        assert got is not None and got != mol.atoms[e["centre"]].chirality


def test_the_answer_table_covers_hydrogens_lone_pairs_and_ring_fusions():
    smiles = {r["smiles"] for r in STEREO}
    assert "C[S@](=O)c1ccccc1" in smiles and "[S@@](=O)(C)c1ccccc1" in smiles
    assert "C1CC[C@H]2CCCC[C@@H]2C1" in smiles
    assert sum(1 for r in STEREO for e in r["marks"] if e["to"] is None) >= 40
    assert sum(len(r["marks"]) for r in STEREO) >= 300


STEREO_SMILES = sorted({r["smiles"] for r in STEREO})


@pytest.mark.parametrize("smiles", STEREO_SMILES)
def test_every_stereocentre_gets_one_mark_that_reads_back_right(smiles):
    mol = perceive(read_smiles(smiles))
    points = layout(mol)
    marks = stereo_marks(mol, points)
    centres = [i for i, a in enumerate(mol.atoms) if a.stereocentre]
    assert sorted(m.centre for m in marks) == centres
    bonds = [frozenset((m.centre, m.to)) for m in marks if m.to is not None]
    assert len(bonds) == len(set(bonds))
    for m in marks:
        assert drawn_chirality(mol, points, m.centre, marks) == mol.atoms[m.centre].chirality
        if m.to is not None:
            assert mol.bond_between(m.centre, m.to).order == 1


def test_a_centre_without_a_mark_or_with_bonds_in_a_line_does_not_read():
    mol = perceive(read_smiles("F[C@](Cl)(Br)I"))
    points = [(-1, 0), (0, 0), (1, 0), (0, 1), (0, -1)]
    assert drawn_chirality(mol, points, 1, []) is None
    flat = [(-1, 0), (0, 0), (1, 0), (2, 0), (3, 0)]
    assert drawn_chirality(mol, flat, 1, [Mark(1, "wedge", 0)]) is None
    assert mark_kind(mol, flat, 1, 0) is None


def test_a_drawn_hydrogen_is_a_text_and_a_wedge():
    smiles = "N[C@@H](C)C(=O)O"
    mol = perceive(read_smiles(smiles))
    points = layout(mol)
    row = next(r for r in STEREO if r["smiles"] == smiles)
    # The reference table uses the package's layout, with the drawn hydrogen placed on it.
    assert [tuple(p) for p in row["points"]] == [tuple(p) for p in points]
    e = next(e for e in row["marks"] if e["to"] is None)
    g = draw_molecule(mol, points, STYLE, marks=[Mark(1, e["kind"], None, tuple(e["hydrogen"]))])
    assert [t.name for t in texts_of(g)] == ["Atom N 1", "Atom O 5", "Atom O 6", "Atom H 7"]
    stereo = [s for s in walk([g]) if isinstance(s, (Wedge, Hash))]
    assert [s.name.split()[1] for s in stereo] == ["6"]     # after the 5 bonds
    assert drawing_overlaps([g]) == (0, 0)


# --- mark placement ----------------------------------------------------------------------------

def marked(smiles):
    mol = perceive(read_smiles(smiles))
    points = layout(mol)
    return mol, points, stereo_marks(mol, points)


def test_a_mark_goes_on_a_substituent_and_the_end_of_a_branch_first():
    # menthol: the isopropyl CH, the methyl and the OH; no ring bond and no drawn hydrogen
    _, _, marks = marked(TARGETS["menthol"])
    assert [(m.centre, m.to) for m in marks] == [(3, 1), (6, 7), (9, 10)]
    # leucine: the N (one neighbour) before the CH2 (two) and the COOH carbon (three)
    _, _, marks = marked("CC(C)C[C@H](N)C(=O)O")
    assert [(m.centre, m.to) for m in marks] == [(4, 5)]


def test_the_ring_fusions_of_decalin_get_hydrogens_away_from_the_shared_bond():
    for smiles, same_kind in (("C1CC[C@H]2CCCC[C@@H]2C1", False), ("C1CC[C@H]2CCCC[C@H]2C1", True)):
        mol, points, marks = marked(smiles)
        assert [(m.centre, m.to) for m in marks] == [(3, None), (8, None)]
        # trans: the two hydrogens on opposite faces (a wedge and a hashed wedge); cis: on the same face
        assert (marks[0].kind == marks[1].kind) is same_kind
        for m, other in zip(marks, (8, 3)):
            c, h, o = points[m.centre], m.hydrogen, points[other]
            assert math.dist(c, h) == pytest.approx(OUTSIDE_LENGTHS[0])
            cos = ((h[0] - c[0]) * (o[0] - c[0]) + (h[1] - c[1]) * (o[1] - c[1])) / math.dist(c, h) / math.dist(c, o)
            assert cos == pytest.approx(-1, abs=1e-6)
        _, _, g = drawn(smiles)
        assert sum(t.name.startswith("Atom H") for t in texts_of(g)) == 2


def test_cholesterol_has_hydrogens_at_c8_c9_c14_and_c17_and_no_mark_in_a_ring():
    mol, points, marks = marked(TARGETS["cholesterol"])
    hydrogens = {m.centre: m.hydrogen for m in marks if m.to is None}
    assert sorted(hydrogens) == [8, 11, 15, 16]          # C17, C14, C9, C8
    ring = ring_bonds(mol)
    assert not [m for m in marks if m.to is not None and frozenset((m.centre, m.to)) in ring]
    # C17 and C14 have room outside the rings. The outside position of C9 is taken by the end of
    # the C19 methyl, and that of C8 by the hydrogen of C14: both go inside a ring, and shorter.
    for c, inside, length in ((8, False, OUTSIDE_LENGTHS[0]), (11, False, OUTSIDE_LENGTHS[0]),
                              (15, True, INSIDE_LENGTHS[0]), (16, True, INSIDE_LENGTHS[0])):
        h = hydrogens[c]
        assert any(_inside(h, [points[a] for a in r]) for r in mol.rings) is inside, c
        assert math.dist(points[c], h) == pytest.approx(length), c


def test_a_hydrogen_is_drawn_before_a_wedge_whose_wide_end_is_another_stereocentre():
    mol, _, marks = marked(TARGETS["cholesterol"])
    # C17 (atom 8): its only bond out of the rings goes to C20 (atom 6), a stereocentre
    assert mol.atoms[6].stereocentre and 6 in mol.neighbors(8)
    assert next(m for m in marks if m.centre == 8).to is None
    # quinine: C8 (atom 12) and C9 (atom 10) likewise
    mol, _, marks = marked("COc1ccc2nccc([C@@H](O)[C@@H]3C[C@@H]4CCN3C[C@@H]4C=C)c2c1")
    assert [(m.centre, m.to) for m in marks if m.centre in (10, 12)] == [(10, 11), (12, None)]


def test_a_centre_whose_substituent_is_refused_gets_a_hydrogen_and_not_a_ring_bond():
    smiles = next(s for s in STEREO_SMILES if s.startswith("CC[C@@H]1[C@@]"))   # erythromycin
    mol, points, marks = marked(smiles)
    # atom 14 is in the ring; its methyl (atom 18) is the middle of three bonds within half a turn
    assert mark_kind(mol, points, 14, 18) is None
    assert next(m for m in marks if m.centre == 14).to is None


def test_only_a_centre_with_three_neighbours_and_a_hydrogen_in_brackets_gets_one_drawn():
    for smiles in ("F[C@](Cl)(Br)I", "C[S@](=O)c1ccccc1"):
        mol, points, _ = marked(smiles)
        assert hydrogen_place(mol, points, 1, set()) is None
    mol, points, _ = marked("N[C@@H](C)C(=O)O")
    h = hydrogen_place(mol, points, 1, set())
    assert math.dist(points[1], h) == pytest.approx(OUTSIDE_LENGTHS[0])


def test_a_drawn_hydrogen_keeps_clear_of_those_drawn_before():
    mol, points, marks = marked(TARGETS["cholesterol"])
    c14 = next(m for m in marks if m.centre == 11)
    # On its own, C8 (atom 16) would take its outside position, almost where the hydrogen of C14 is
    # (0.1 bond lengths away; the five-membered ring bends the bonds of C14 slightly).
    assert math.dist(hydrogen_place(mol, points, 16, set()), c14.hydrogen) < 0.2
    before = {m for m in marks if m.to is None and m.centre < 16}
    assert hydrogen_place(mol, points, 16, before) == next(m for m in marks if m.centre == 16).hydrogen


# The quinolizidinium ion and the same rings with other atoms at the fusion: the centre is atom 3.
FUSIONS_WITH_TEXT = [
    ("C1CC[N@@H+]2CCCC[C@@H]2C1", [("N", "normal"), ("+", "sup")]),
    ("C1CC[N@H+]2CCCC[C@@H]2C1", [("N", "normal"), ("+", "sup")]),
    ("C1CC[P@@H+]2CCCC[C@@H]2C1", [("P", "normal"), ("+", "sup")]),
    ("C1CC[B@@H-]2CCCC[C@@H]2C1", [("B", "normal"), ("−", "sup")]),
    ("C1CC[Si@@H]2CCCC[C@@H]2C1", [("Si", "normal")]),
]


@pytest.mark.parametrize("smiles,runs", FUSIONS_WITH_TEXT, ids=[s for s, _ in FUSIONS_WITH_TEXT])
def test_a_hydrogen_drawn_for_a_mark_is_not_said_again_in_the_text_of_its_atom(smiles, runs):
    mol, points, marks = marked(smiles)
    assert [(m.centre, m.to) for m in marks] == [(3, None), (8, None)]
    _, _, group = drawn(smiles)
    # Three labels: the atom without its hydrogen, and the two hydrogens at the wedges.
    said = {t.name: [(run.text, run.kind) for run in t.runs] for t in texts_of(group)}
    element = mol.atoms[3].element
    assert said == {f"Atom {element} 4": runs, "Atom H 11": [("H", "normal")], "Atom H 12": [("H", "normal")]}
    assert drawing_overlaps([group]) == (0, 0)
    # The lines from the centre stop at the bounding box of the label as drawn.
    box = label_box(Label([Run(*r) for r in runs], 0), (0, 0), FONT, STYLE.font_size, STYLE.script_scale)
    text = next(t for t in texts_of(group) if t.name == f"Atom {element} 4")
    assert text.width == pytest.approx(box[2] - box[0])


def test_the_text_keeps_a_hydrogen_that_no_mark_draws():
    # The mark of the N+ goes to its methyl (atom 0): the hydrogen stays in the label.
    smiles = "C[N@@H+]1CCC[C@@H]2CCCC[C@H]12"
    _, _, marks = marked(smiles)
    assert next(m for m in marks if m.centre == 1).to == 0
    _, _, group = drawn(smiles)
    text = next(t for t in texts_of(group) if t.name == "Atom N 2")
    # The ring is on the right of the N, so its H is on the left.
    assert [(run.text, run.kind) for run in text.runs] == [("H", "normal"), ("N", "normal"), ("+", "sup")]
    # With a mark given by hand on a ring bond instead of the drawn hydrogen, the H is back in
    # the label, below the N+ (the outside gap of the ring fusion is below it).
    smiles = "C1CC[N@@H+]2CCCC[C@@H]2C1"
    mol, points, _ = marked(smiles)
    group = draw_molecule(mol, points, STYLE, marks=[Mark(3, "wedge", 2)])
    said = {t.name: "".join(run.text for run in t.runs) for t in texts_of(group)}
    assert said == {"Atom N 4": "N+", "Atom N 4 H": "H"}
    assert said_below(group, "Atom N 4", "Atom N 4 H")


def test_a_line_that_leaves_the_symbol_into_the_h_below_it_stops_at_the_h():
    # The symbol (-1, -1, 1, 1) and its H below (-1, 0.9, 1, 3), y down as on a slide.
    boxes = {0: [(-1, -1, 1, 1), (-1, 0.9, 1, 3)]}
    start, end = clip_ends(((0, 0), (0.5, 5)), (0, 1), boxes)
    assert start == pytest.approx((0.3, 3)) and end == (0.5, 5)
    # Without the H, it stops at the symbol.
    start, _ = clip_ends(((0, 0), (0.5, 5)), (0, 1), {0: boxes[0][:1]})
    assert start == pytest.approx((0.1, 1))
    # A line that passes through the H, with neither end in it, is kept whole (an overlap).
    assert clip_ends(((-3, 2), (3, 2)), (0, 1), boxes) == ((-3, 2), (3, 2))


def test_a_bond_of_the_atom_into_its_h_below_stops_at_the_h_on_the_drawing():
    # P with bonds to the right, left and top, and one 25 degrees off straight down: its H goes
    # below (the position farthest from the bonds), and that bond runs through the H. The layout
    # never produces this (in none of 260 molecules), so the points are given by hand.
    mol = perceive(read_smiles("C[PH](C)(C)C"))
    down = math.radians(-65)
    points = [(-1.0, 0.0), (0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (math.cos(down), math.sin(down))]
    group = draw_molecule(mol, points, STYLE, marks=[])
    said = {t.name: "".join(run.text for run in t.runs) for t in texts_of(group)}
    assert said == {"Atom P 2": "P", "Atom P 2 H": "H"}
    assert said_below(group, "Atom P 2", "Atom P 2 H")
    # The line to the lower carbon is drawn only beyond the H, so it passes through no label.
    assert drawing_overlaps([group]) == (0, 0)
    (line,) = lines_named(group, "Bond 4")
    h = next(t for t in texts_of(group) if t.name == "Atom P 2 H")
    h_middle = h.y + FONT.ascent_height(h.size) - FONT.capital_height(h.size) / 2
    assert min(line.p0[1], line.p1[1]) > h_middle   # it starts below the middle of the H


def test_the_h_above_or_below_has_the_colour_of_its_atom():
    _, _, group = drawn("C1CCNCC1", Style(colour="atoms"))
    colours = {t.name: t.colour for t in texts_of(group)}
    assert colours == {"Atom N 4": Style().element_colours["N"], "Atom N 4 H": Style().element_colours["N"]}


def test_the_lines_and_the_circle_keep_clear_of_the_h_above_or_below():
    for smiles in ("c1cc[nH]c1", "c1ccc2[nH]ccc2c1", "C[NH+](C)C", "C[SiH](C)C", "[cH+]1cccccc1", "CCNCC"):
        for style in (STYLE, Style(aromatic="circle"), Style(colour="atoms_and_bonds")):
            mol, _, group = drawn(smiles, style)
            assert any(t.name.endswith(" H") for t in texts_of(group)), smiles
            assert drawing_overlaps([group]) == (0, 0), (smiles, style)


def said_below(group, name, h_name) -> bool:
    texts = {t.name: t for t in texts_of(group)}
    return texts[h_name].y - texts[name].y == pytest.approx(STYLE.font_size * HYDROGEN_LINE)


def test_a_ring_n_with_its_hydrogen_at_a_mark_is_n_plus_whatever_side_the_ring_is_on():
    # Crinine with a proton on its N (atom 2), which is on the left side of its rings: with
    # its hydrogen the label would be HN+.
    smiles = "C1C[N@H+]2CC3=CC4=C(C=C3[C@]15[C@H]2C[C@H](C=C5)O)OCO4"
    _, _, marks = marked(smiles)
    assert sorted(m.centre for m in marks if m.to is None) == [2, 11]
    _, _, group = drawn(smiles)
    said = {t.name: "".join(run.text for run in t.runs) for t in texts_of(group)}
    assert said == {"Atom N 3": "N+", "Atom O 17": "HO", "Atom O 18": "O", "Atom O 20": "O",
                    "Atom H 21": "H", "Atom H 22": "H"}
    assert drawing_overlaps([group]) == (0, 0)


def test_no_hydrogen_is_drawn_in_the_place_where_the_text_has_its_h():
    # Matrine with a proton on N-1 (atom 7). The gap outside the rings is to the right of the
    # N, where its label has the H: a hydrogen drawn there would sit right after the + of
    # N+, with its wedge hidden between the two. The mark goes on a bond and the label is NH+.
    smiles = "C1C[C@@H]2[C@H]3CCC[N@H+]4[C@H]3[C@@H](CCC4)CN2C(=O)C1"
    mol, points, marks = marked(smiles)
    gap = (points[7][0] + OUTSIDE_LENGTHS[0], points[7][1])
    assert not any(_inside(gap, [points[a] for a in r]) for r in mol.rings)
    assert all(points[j][0] < gap[0] - 0.4 for j in mol.neighbors(7))
    assert next(m for m in marks if m.centre == 7).to is not None
    assert hydrogen_place(mol, points, 7, set()) is None
    _, _, group = drawn(smiles)
    text = next(t for t in texts_of(group) if t.name == "Atom N 8")
    assert [(run.text, run.kind) for run in text.runs] == [("N", "normal"), ("H", "normal"), ("+", "sup")]
    assert sum(t.name.startswith("Atom H") for t in texts_of(group)) == 4
    assert drawing_overlaps([group]) == (0, 0)


def test_what_is_seen_of_the_wedge_to_a_drawn_hydrogen_leaves_out_both_texts():
    # Quinolizidinium: the N+ (atom 3) and the CH (atom 8) of the ring fusion, each with an H
    # directly below it at 0.8 bond lengths. The N+ label takes its half height and margin off the
    # wedge, as the H does: the wedge from the N+ shows too little, the one from the carbon enough.
    mol, points, _ = marked("C1CC[N@@H+]2CCCC[C@@H]2C1")
    scene = _Scene(mol, points, set())
    half = (FONT.capital_height(STYLE.font_size) / 2 + STYLE.label_margin) / STYLE.bond_length
    seen = {}
    for centre in (3, 8):
        c = points[centre]
        reach = scene.reach(centre, (0.0, -1.0))
        assert reach == pytest.approx(half if centre == 3 else 0.0)
        seen[centre] = scene.seen(centre, (c[0], c[1] - 0.8))
        assert seen[centre] == pytest.approx(0.8 - reach - half)
    assert seen[3] < STYLE.shortest_line <= seen[8]


def test_a_bond_with_a_mark_is_not_offered_to_the_centre_at_its_other_end(monkeypatch):
    # In the test molecules no two centres compete for the same bond (a bond between two
    # stereocentres is the last choice for both), so here both centres are offered only that bond.
    from smiles2pptx.depict import stereo

    real = stereo.candidates
    monkeypatch.setattr(stereo, "candidates", lambda mol, points, centre, used: [
        j for j in real(mol, points, centre, used) if isinstance(j, int) and mol.atoms[j].stereocentre])
    _, _, marks = marked("[C@H]([C@H](F)Cl)(F)Cl")
    assert [(m.centre, m.to) for m in marks] == [(0, 1)]     # centre 1 finds the bond taken


def test_over_the_stereo_molecules_no_wide_end_is_a_stereocentre_and_4_marks_are_in_rings():
    hydrogens = ring_marks = 0
    for smiles in STEREO_SMILES:
        mol, points, marks = marked(smiles)
        centres = {i for i, a in enumerate(mol.atoms) if a.stereocentre}
        assert not [m for m in marks if m.to in centres], smiles
        ring = ring_bonds(mol)
        hydrogens += sum(1 for m in marks if m.to is None)
        ring_marks += sum(1 for m in marks if m.to is not None and frozenset((m.centre, m.to)) in ring)
    # The 4 in rings: three centres with four ring bonds (artemisinin, morphine, strychnine),
    # and one in morphine whose hydrogen has no clear position.
    assert (hydrogens, ring_marks) == (30, 4)


@pytest.mark.parametrize("smiles", STEREO_SMILES)
def test_drawn_hydrogens_meet_no_line_or_text(smiles):
    _, _, g = drawn(smiles)
    assert drawing_overlaps([g]) == (0, 0)


# --- the 15 molecules ----------------------------------------------------------------------

@pytest.mark.parametrize("key", sorted(TARGETS))
def test_the_15_molecules_have_no_overlaps_on_the_drawing(key):
    _, _, g = drawn(TARGETS[key])
    assert drawing_overlaps([g]) == (0, 0)
    x0, y0, x1, y1 = g.bbox()
    assert (round(x0, 6), round(y0, 6)) == (0, 0)
    assert outside([g], 960, 540) == 0


def test_the_15_molecules_have_the_shapes_of_their_bonds_and_texts():
    # Lines: one for each single bond, two for a double bond, none for a wedge.
    # Cholesterol: 31 bonds, one double, 8 marks. 4 of the marks go to hydrogens drawn at C8,
    # C9, C14 and C17, so only 4 bonds become wedges (32 - 4 = 28 lines), and there are 4 H
    # labels besides OH.
    expect = {"decanoic": (12, 2, 0), "glycine": (5, 3, 0), "menthol": (8, 1, 3), "nabenzoate": (13, 3, 0),
              "caffeine": (19, 6, 0), "cholesterol": (28, 5, 8), "crown": (18, 6, 0)}
    for key, (lines, texts, stereo) in expect.items():
        _, _, g = drawn(TARGETS[key])
        shapes = list(walk([g]))
        assert sum(isinstance(s, Line) for s in shapes) == lines, key
        assert sum(isinstance(s, Text) for s in shapes) == texts, key
        assert sum(isinstance(s, (Wedge, Hash)) for s in shapes) == stereo, key


# The NH of a ring or a chain, and where its H is written: beside the symbol (HN, NH) or as a
# separate label above or below it.
NH_PLACES = [("O=C1NC(=O)c2ccccc12", "HN"), ("C1=CNC=C1", "below"), ("O=c1cccc[nH]1", "above"),
             ("O=C1C=CC=CN1", "above"), ("c1c[nH]cn1", "below"), ("CCNCC", "above"),
             ("c1ccc2[nH]ccc2c1", "below"), ("C1CCNCC1", "below"), ("c1ccc2[nH]nnc2c1", "below")]


@pytest.mark.parametrize("smiles,place", NH_PLACES, ids=[s for s, _ in NH_PLACES])
def test_the_h_of_an_nh_goes_where_it_is_farthest_from_the_bonds(smiles, place):
    mol, points, group = drawn(smiles)
    n = next(i for i, a in enumerate(mol.atoms) if a.element == "N" and a.h_count == 1)
    texts = {t.name: t for t in texts_of(group)}
    text = texts[f"Atom N {n + 1}"]
    h = texts.get(f"Atom N {n + 1} H")
    if place == "HN":
        assert "".join(run.text for run in text.runs) == "HN" and h is None
    else:
        assert "".join(run.text for run in text.runs) == "N"
        assert [(run.text, run.kind) for run in h.runs] == [("H", "normal")]
        # Directly above or below the N, HYDROGEN_LINE font sizes away, and outside the rings.
        assert h.x + h.width / 2 == pytest.approx(text.x + text.width / 2)
        step = HYDROGEN_LINE * STYLE.font_size * (1 if place == "below" else -1)
        assert h.y - text.y == pytest.approx(step)
        spot = (points[n][0], points[n][1] - step / STYLE.bond_length)
        assert not any(_inside(spot, [points[a] for a in r]) for r in mol.rings)
    assert drawing_overlaps([group]) == (0, 0)
    # Neither bond of the N to a carbon is cut short by the H.
    for j in mol.neighbors(n):
        k = mol.bonds.index(mol.bond_between(n, j))
        if mol.atoms[j].element != "C":
            continue   # a bond between two labels looks short whatever the H does (N-N of a triazole)
        longest = max(math.dist(line.p0, line.p1) for line in lines_named(group, f"Bond {k + 1}"))
        assert longest > STYLE.bond_length / 2


def test_charges_are_superscripts():
    _, _, g = drawn(TARGETS["nabenzoate"])
    runs = {t.name: [(r.text, r.kind) for r in t.runs] for t in texts_of(g)}
    assert ("+", "sup") in runs["Atom Na 10"]
    assert ("−", "sup") in runs["Atom O 3"]
    _, _, g = drawn(TARGETS["mepyridinium"])
    assert [(r.text, r.kind) for r in texts_of(g)[0].runs] == [("N", "normal"), ("+", "sup")]


def test_shapes_are_named_after_their_bond_or_atom():
    mol, _, g = drawn("CC(=O)O")
    names = sorted(s.name for s in walk([g]))
    assert names == ["Atom O 3", "Atom O 4", "Bond 1", "Bond 2", "Bond 2", "Bond 3"]


# --- double bonds ----------------------------------------------------------------------------

def _second_line(g, number):
    lines = lines_named(g, f"Bond {number}")
    assert len(lines) == 2
    return lines


def test_a_ring_double_bond_has_its_second_line_inside_the_ring():
    mol, points, g = drawn("c1ccccc1")
    L = STYLE.bond_length
    ends = [p for s in walk([g]) for p in (s.p0, s.p1)]
    centre = (sum(x for x, _ in ends) / len(ends), sum(y for _, y in ends) / len(ends))
    for k, b in enumerate(mol.bonds):
        if b.order != 2:
            continue
        a, c = _second_line(g, k + 1)
        da = math.dist(centre, ((a.p0[0] + a.p1[0]) / 2, (a.p0[1] + a.p1[1]) / 2))
        dc = math.dist(centre, ((c.p0[0] + c.p1[0]) / 2, (c.p0[1] + c.p1[1]) / 2))
        inner, outer = (a, c) if da < dc else (c, a)
        assert math.dist(inner.p0, inner.p1) == pytest.approx(L * (1 - 2 * STYLE.inner_line_trim))
        assert math.dist(outer.p0, outer.p1) == pytest.approx(L)
        assert abs(da - dc) == pytest.approx(STYLE.double_bond_gap * L)


def test_a_double_bond_to_an_end_atom_is_centred():
    mol, points, g = drawn("CC(=O)C")
    P = [(x, -y) for x, y in points]
    k = next(k for k, b in enumerate(mol.bonds) if b.order == 2)
    assert double_bond_side(mol, k, P) == 0
    a, b = _second_line(g, k + 1)
    L = STYLE.bond_length
    first = lines_named(g, "Bond 1")[0]            # C0-C1: no labels, not clipped
    carbon = (points[1][0] * L + first.p0[0] - points[0][0] * L, -points[1][1] * L + first.p0[1] + points[0][1] * L)
    gap = STYLE.double_bond_gap * L
    assert _line_distance(carbon, a) == pytest.approx(gap / 2)
    assert _line_distance(carbon, b) == pytest.approx(gap / 2)
    assert _line_distance(a.p0, b) == pytest.approx(gap)


def _line_distance(p, line):
    (x0, y0), (x1, y1) = line.p0, line.p1
    return abs((x1 - x0) * (y0 - p[1]) - (x0 - p[0]) * (y1 - y0)) / math.dist(line.p0, line.p1)


def test_a_chain_double_bond_has_its_second_line_on_the_side_of_its_neighbours():
    mol = perceive(read_smiles(r"C/C=C\C"))
    P = [(0.0, 1.0), (0.5, 0.0), (1.5, 0.0), (2.0, 1.0)]   # cis: both neighbours above, y down
    k = next(k for k, b in enumerate(mol.bonds) if b.order == 2)
    assert double_bond_side(mol, k, P) == 1      # the left normal of (0.5,0)->(1.5,0) is +y
    P = [(0.0, -1.0), (0.5, 0.0), (1.5, 0.0), (2.0, -1.0)]
    assert double_bond_side(mol, k, P) == -1


# --- labels and clipping -------------------------------------------------------------------

def test_lines_stop_the_label_margin_short_of_the_text():
    mol, points, g = drawn("OCC(=O)N")
    L = STYLE.bond_length
    for t in texts_of(g):
        baseline = t.y + FONT.ascent_height(t.size)
        box = (t.x, baseline - FONT.capital_height(t.size), t.x + t.width, baseline)
        near = []
        for line in (s for s in walk([g]) if isinstance(s, Line)):
            assert not segment_in_box(line.p0, line.p1, box)
            gap = min(_box_distance(p, box) for p in (line.p0, line.p1))
            if gap < L / 2:
                near.append(gap)
        assert near and all(d == pytest.approx(STYLE.label_margin, abs=1e-6) for d in near)


def _box_distance(p, box):
    # The margin widens the box on every side, so a line ends this far out in x or in y.
    dx = max(box[0] - p[0], 0, p[0] - box[2])
    dy = max(box[1] - p[1], 0, p[1] - box[3])
    return max(dx, dy)


def test_the_capitals_of_a_text_are_centred_on_its_atom():
    mol, points, g = drawn("CO")
    L = STYLE.bond_length
    P = [(x * L, -y * L) for x, y in points]
    line = lines_named(g, "Bond 1")[0]          # from the carbon (no label, not clipped) to the O
    shift = (line.p0[0] - P[0][0], line.p0[1] - P[0][1])
    o = (P[1][0] + shift[0], P[1][1] + shift[1])
    t = texts_of(g)[0]
    baseline = t.y + FONT.ascent_height(t.size)
    assert baseline - FONT.capital_height(t.size) / 2 == pytest.approx(o[1])
    first = [r.text for r in t.runs][0][0]
    assert first == "O" and t.x + FONT.char_width("O") * t.size / 2 == pytest.approx(o[0])


def test_clip_out_keeps_a_segment_that_passes_through_and_drops_one_inside():
    box = (0, 0, 2, 2)
    assert clip_out((-1, 1), (3, 1), box) == ((-1, 1), (3, 1))
    assert clip_out((0.5, 0.5), (1.5, 1.5), box) is None
    assert clip_out((1, 1), (5, 1), box) == ((2, 1), (5, 1))
    assert clip_out((5, 1), (1, 1), box) == ((5, 1), (2, 1))


# --- wedge shapes ------------------------------------------------------------------------------

def test_a_wedge_is_narrow_at_the_centre_and_a_hash_grows_away_from_it():
    # The marks of menthol are on its substituents: a wedge to the isopropyl, hashed wedges to
    # the methyl and to the OH. The one to the OH stops short of the label, and its wide end is
    # as wide as the whole bond would be there (the angle of an unclipped wedge).
    mol, points, g = drawn(TARGETS["menthol"])
    marks = stereo_marks(mol, points)
    L = STYLE.bond_length
    found, clipped = 0, 0
    for s in walk([g]):
        if isinstance(s, Wedge):
            narrow = math.dist(s.points[0], s.points[3])
            wide = math.dist(s.points[1], s.points[2])
            reach = math.dist(_mid(s.points[0], s.points[3]), _mid(s.points[1], s.points[2])) / L
            assert narrow < wide and wide == pytest.approx(STYLE.wedge_width * reach, rel=0.05)
            found += 1
        elif isinstance(s, Hash):
            lengths = [math.dist(a, b) for a, b in s.segments]
            assert all(b >= a - 1e-9 for a, b in zip(lengths, lengths[1:])) and len(lengths) >= 4
            reach = math.dist(_mid(*s.segments[0]), _mid(*s.segments[-1])) / L
            assert lengths[-1] == pytest.approx(STYLE.wedge_width * reach) and lengths[-1] > lengths[0]
            clipped += reach < 0.9
            found += 1
    assert found == len(marks) == 3
    assert clipped == 1


def _mid(a, b):
    return ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)


# --- counts on a drawing ----------------------------------------------------------------------

def test_drawing_overlaps_counts_lines_wedges_and_hashes_through_a_text_and_texts_that_meet():
    size = 14.0
    t = Text(10, 10, 20, FONT.line_height(size), [Run("O")], size)
    baseline = 10 + FONT.ascent_height(size)
    y = baseline - FONT.capital_height(size) / 2          # through the middle of the capitals
    line = Line((0, y), (40, y), 1)
    wedge = Wedge([(0, y - 0.1), (40, y - 2), (40, y + 2), (0, y + 0.1)])
    hashed = Hash([((x, y - 1), (x, y + 1)) for x in (0, 10, 20, 30, 40)], 1)
    assert drawing_overlaps([t, line]) == (1, 0)
    assert drawing_overlaps([t, wedge]) == (1, 0)
    assert drawing_overlaps([t, hashed]) == (1, 0)
    assert drawing_overlaps([t, Line((0, y - 20), (40, y - 20), 1)]) == (0, 0)
    assert drawing_overlaps([t, t.moved(15, 0)]) == (0, 1)
    assert drawing_overlaps([t, t.moved(25, 0)]) == (0, 0)


def test_outside_counts_the_shapes_beyond_each_side_of_the_slide():
    _, _, g = drawn(TARGETS["caffeine"])
    inside = g.centred_at(480, 270)
    n = len(list(walk([g])))
    assert outside([inside], 960, 540) == 0
    x0, y0, x1, y1 = inside.bbox()
    for dx, dy in ((-x1 - 1, 0), (960 - x0 + 1, 0), (0, -y1 - 1), (0, 540 - y0 + 1)):
        assert outside([inside.moved(dx, dy)], 960, 540) == n
    assert 0 < outside([inside.moved(-x0 - 1, 0)], 960, 540) < n
