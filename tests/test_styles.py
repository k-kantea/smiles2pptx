"""Tests for the drawing options of the style: circles in aromatic rings, colours, and
substituents drawn as abbreviations.

- aromatic rings: for 138 SMILES the atoms of Molecule.aromatic_rings are those RDKit finds
  aromatic (tests/data/aromatic_reference.json, made by tools/make_aromatic_reference.py),
  except for three that differ on purpose; a ring written as a Kekule structure counts like
  one written in lower case
- circle: each bond of an aromatic ring is a single line and the ring has a circle; other
  double bonds stay; the circle keeps clear of the labels
- colours: labels in the colours of their elements; with coloured bonds each line, wedge and
  hashed wedge is cut at the middle of its bond
- abbreviations: where a group is found and where not, the molecule that is left, the two
  forms of each label and the side each one extends to
- the 15 molecules in each style: no line or circle through a label, no labels that meet
- the default style uses none of these; the command line and Deck accept the settings
"""
import json
import math
import zipfile
from pathlib import Path

import pytest
from pptx_check import problems

from smiles2pptx import Deck
from smiles2pptx.check import circle_meets_box, drawing_overlaps, outside
from smiles2pptx.chem.abbreviate import abbreviate
from smiles2pptx.chem.perceive import perceive
from smiles2pptx.chem.smiles import read_smiles
from smiles2pptx.cli import main
from smiles2pptx.depict import draw_molecule
from smiles2pptx.depict.bonds import halves, wedge_shapes
from smiles2pptx.depict.labels import atom_label, formula_runs, label_box
from smiles2pptx.depict.stereo import drawn_chirality, stereo_marks
from smiles2pptx.figure import Circle, Hash, Line, Run, Text, Wedge, walk
from smiles2pptx.fonts import get_font
from smiles2pptx.layout import layout
from smiles2pptx.reaction import draw_reaction
from smiles2pptx.style import ABBREVIATIONS, Style
from smiles2pptx.writers.pptx import SLIDE_HEIGHT, SLIDE_WIDTH

DATA = Path(__file__).resolve().parent / "data"
AROMATIC = json.loads((DATA / "aromatic_reference.json").read_text(encoding="utf-8"))["samples"]
TARGETS = {m["id"]: m["smiles"] for m in json.loads((DATA / "test_molecules.json").read_text(encoding="utf-8"))
           if m["target"]}
FONT = get_font("Arial")
BLACK = "000000"
RED, BLUE = Style().element_colours["O"], Style().element_colours["N"]
CIRCLE = Style(aromatic="circle")
ATOMS = Style(colour="atoms")
BONDS = Style(colour="atoms_and_bonds")
ALL = Style(abbreviate="all")
WAYS = {"circle": CIRCLE, "atoms": ATOMS, "atoms_and_bonds": BONDS, "abbreviate": ALL,
        "everything": Style(aromatic="circle", colour="atoms_and_bonds", abbreviate="all")}
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
PPH3 = "c1ccc(cc1)P(c2ccccc2)c3ccccc3"


def read(smiles):
    return perceive(read_smiles(smiles))


def drawn(smiles, style):
    mol = abbreviate(read(smiles), style.groups)
    points = layout(mol, style)
    return mol, points, draw_molecule(mol, points, style)


def shapes(smiles, style, kind=None):
    found = list(walk([drawn(smiles, style)[2]]))
    return found if kind is None else [s for s in found if isinstance(s, kind)]


def aromatic_atoms(mol):
    return sorted({a for r in mol.aromatic_rings for a in mol.rings[r]})


def said(text):
    return "".join(run.text for run in text.runs)


# --- aromatic rings ----------------------------------------------------------------------

@pytest.mark.parametrize("row", AROMATIC, ids=[r["smiles"][:40] for r in AROMATIC])
def test_the_aromatic_atoms_are_those_of_the_table(row):
    assert aromatic_atoms(read(row["smiles"])) == row.get("expect", row["rdkit"])


def test_the_aromatic_table_is_not_all_alike():
    answers = [r.get("expect", r["rdkit"]) for r in AROMATIC]
    assert len(AROMATIC) == 138
    assert sum(1 for a in answers if not a) >= 20 and sum(1 for a in answers if a) >= 100
    assert len({tuple(a) for a in answers}) >= 30
    assert sum(1 for r in AROMATIC if "Kekule" in r["shows"]) >= 25
    assert [r["smiles"] for r in AROMATIC if "expect" in r and r["expect"] != r["rdkit"]] == [
        "C1C=C2C=CC=CC2=C1", "C=C1C=CC=CC1=C", "C1=CC2=CC3=CC=C(N3)C=C4C=CC(=N4)C=C5C=CC(=N5)C=C1N2"]


@pytest.mark.parametrize("lower,kekule,rings", [
    ("c1ccccc1", "C1=CC=CC=C1", 1),
    ("c1ccc2ccccc2c1", "C1=CC=C2C=CC=CC2=C1", 2),
    ("c1ccc2ccccc2c1", "C1=CC2=CC=CC=C2C=C1", 2),
    ("Cn1cnc2c1c(=O)n(C)c(=O)n2C", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", 2),
    ("c1ccc2cccc2cc1", "C1=CC2=CC=CC=CC2=C1", 2),
    ("c1ccn2cccc2c1", "C1=CC2=CC=CN2C=C1", 2),
    ("[cH-]1cccc1", "[CH-]1C=CC=C1", 1),
    ("c1ccc2c(c1)ccc1ccccc12", "C1=CC=C2C(=C1)C=CC1=CC=CC=C21", 3),
])
def test_a_kekule_structure_counts_like_lower_case(lower, kekule, rings):
    assert len(read(lower).aromatic_rings) == len(read(kekule).aromatic_rings) == rings


@pytest.mark.parametrize("smiles,why", [
    ("C1=CC=CC=CC=C1", "8 electrons"),
    ("C1=CC=C1", "4 electrons"),
    ("C1=CCC=C1", "a CH2 in the ring"),
    ("O=C1C=CC(=O)C=C1", "4 electrons: the two C=O give none"),
    ("C=C1C=CC=C1", "a C=C out of the ring"),
    ("C1=CC2=CC=CC2=C1", "8 electrons over the two rings"),
    ("C1CCCCC1", "no double bond"),
])
def test_rings_that_are_not_aromatic(smiles, why):
    assert read(smiles).aromatic_rings == [], why


def test_only_the_bonds_of_aromatic_rings_are_aromatic_bonds():
    assert len(read("Cc1ccccc1").aromatic_bonds()) == 6
    biphenyl = read("c1ccccc1c1ccccc1")
    assert len(biphenyl.bonds) == 13 and len(biphenyl.aromatic_bonds()) == 12
    indane = read("c1ccc2c(c1)CCC2")
    assert len(indane.rings) == 2 and len(indane.aromatic_rings) == 1 and len(indane.aromatic_bonds()) == 6
    assert read("C1=CCCCC1").aromatic_bonds() == set()


# --- circle ------------------------------------------------------------------------------

@pytest.mark.parametrize("smiles", ["c1ccccc1", "C1=CC=CC=C1"])
def test_benzene_with_a_circle_is_six_lines_and_a_circle(smiles):
    found = shapes(smiles, CIRCLE)
    lines = [s for s in found if isinstance(s, Line)]
    circles = [s for s in found if isinstance(s, Circle)]
    assert len(lines) == 6 and len(circles) == 1 and len(found) == 7
    assert sorted(s.name for s in lines) == sorted(f"Bond {k}" for k in range(1, 7))
    ends = {(round(x, 6), round(y, 6)) for s in lines for x, y in (s.p0, s.p1)}
    assert len(ends) == 6   # every line runs from atom to atom
    c = circles[0]
    assert c.name == "Ring 1" and c.width == CIRCLE.line_width and c.colour == BLACK
    assert c.centre == pytest.approx((sum(x for x, _ in ends) / 6, sum(y for _, y in ends) / 6))
    to_bond = CIRCLE.bond_length * math.sqrt(3) / 2
    assert c.radius == pytest.approx(CIRCLE.circle_radius * to_bond)


def test_the_default_way_has_no_circle_and_keeps_the_double_bonds():
    found = shapes("c1ccccc1", Style())
    assert len(found) == 9 and all(isinstance(s, Line) for s in found)


def test_each_aromatic_ring_gets_its_circle_and_other_double_bonds_stay():
    assert len(shapes("c1ccc2ccccc2c1", CIRCLE, Circle)) == 2
    assert len(shapes("c1ccc2ccccc2c1", CIRCLE, Line)) == 11
    assert len(shapes(PPH3, CIRCLE, Circle)) == 3
    assert shapes("C1=CCCCC1", CIRCLE) == shapes("C1=CCCCC1", Style())
    # 2-pyridone: the ring has a circle and the C=O keeps its two lines.
    mol, _, group = drawn("O=c1cccc[nH]1", CIRCLE)
    carbonyl = next(k for k, b in enumerate(mol.bonds) if "O" in (mol.atoms[b.begin].element, mol.atoms[b.end].element))
    lines = [s for s in walk([group]) if isinstance(s, Line)]
    assert len([s for s in lines if s.name == f"Bond {carbonyl + 1}"]) == 2
    assert len(lines) == 6 + 2 and len([s for s in walk([group]) if isinstance(s, Circle)]) == 1
    # indane: only the benzene ring.
    assert len(shapes("c1ccc2c(c1)CCC2", CIRCLE, Circle)) == 1


@pytest.mark.parametrize("smiles", ["C[n+]1ccccc1", "C[n+]1ccn(C)c1", "[O-][n+]1ccccc1", "C[n+]1ccsc1"])
def test_the_circle_is_smaller_where_a_text_reaches_into_the_ring(smiles):
    mol, _, group = drawn(smiles, CIRCLE)
    circle, = [s for s in walk([group]) if isinstance(s, Circle)]
    n = len(mol.rings[0])
    to_bond = CIRCLE.bond_length / (2 * math.tan(math.pi / n))
    assert CIRCLE.circle_min_radius * to_bond - 1e-9 <= circle.radius < CIRCLE.circle_radius * to_bond - 0.5
    assert drawing_overlaps([group]) == (0, 0)


def test_the_circle_meets_no_hydrogen_of_a_ring_atom():
    for smiles in ("O=c1[nH]c(=O)c2[nH]cnc2[nH]1", "c1c[nH]cn1", "C1=CNC=C1", "O=c1cccc[nH]1"):
        assert drawing_overlaps([drawn(smiles, CIRCLE)[2]]) == (0, 0), smiles


def test_a_circle_that_reaches_a_text_is_counted():
    text = Text(10, 0, 10, FONT.line_height(14), [Run("N")], 14)
    box = (10, FONT.ascent_height(14) - FONT.capital_height(14), 20, FONT.ascent_height(14))
    middle = (box[1] + box[3]) / 2
    assert drawing_overlaps([text, Circle((0, middle), 9.9, 1)]) == (0, 0)
    assert drawing_overlaps([text, Circle((0, middle), 10.1, 1)]) == (1, 0)
    assert drawing_overlaps([text, Circle((15, middle), 30, 1)]) == (1, 0)   # the label inside the circle
    assert circle_meets_box((0, 0), 5, (3, 3, 9, 9)) and not circle_meets_box((0, 0), 4, (3, 3, 9, 9))


def test_a_pptx_with_circles_can_be_opened(tmp_path):
    deck = Deck(CIRCLE)
    deck.add_molecule("c1ccc2ccccc2c1")
    deck.save(tmp_path / "a.pptx")
    assert problems(tmp_path / "a.pptx") == []
    with zipfile.ZipFile(tmp_path / "a.pptx") as z:
        slide = z.read("ppt/slides/slide1.xml").decode("utf-8")
    assert slide.count('prst="ellipse"') == 2 and 'name="Ring 2"' in slide


# --- colours -----------------------------------------------------------------------------

def test_the_default_way_is_black():
    for smiles in TARGETS.values():
        assert {s.colour for s in shapes(smiles, Style())} == {BLACK}


def test_atom_texts_have_the_colours_of_their_elements():
    texts = shapes("NCC(=O)O", ATOMS, Text)
    assert {said(t): t.colour for t in texts} == {"H2N": BLUE, "O": RED, "OH": RED}
    assert {s.colour for s in shapes("NCC(=O)O", ATOMS, Line)} == {BLACK}
    assert len(shapes("NCC(=O)O", ATOMS)) == len(shapes("NCC(=O)O", Style()))
    assert {t.colour for t in shapes("[Na+].[Cl-]", ATOMS, Text)} == {BLACK, ATOMS.element_colours["Cl"]}


def test_with_coloured_bonds_a_line_is_cut_at_the_middle_of_its_bond():
    mol, points, group = drawn("CCO", BONDS)
    lines = [s for s in walk([group]) if isinstance(s, Line)]
    assert [(s.name, s.colour) for s in lines] == [("Bond 1", BLACK), ("Bond 2", BLACK), ("Bond 2", RED)]
    black, red = lines[1], lines[2]
    assert black.p1 == pytest.approx(red.p0)
    # The cut is half a bond from the carbon, where the black line starts.
    assert math.dist(black.p0, black.p1) == pytest.approx(BONDS.bond_length / 2)
    text, = [s for s in walk([group]) if isinstance(s, Text)]
    assert text.colour == RED


def test_both_lines_of_a_double_bond_are_cut_and_a_bond_of_one_colour_is_not():
    lines = shapes("CC(C)=O", BONDS, Line)
    assert sorted(s.colour for s in lines if s.name == "Bond 3") == [BLACK, BLACK, RED, RED]
    assert len(lines) == 2 + 4
    peroxide = shapes("COOC", BONDS, Line)
    assert [s.colour for s in peroxide if s.name == "Bond 2"] == [RED]


def test_halves_gives_each_part_the_colour_of_the_nearer_atom():
    p, q = (0.0, 0.0), (10.0, 0.0)
    assert halves((p, q), p, q, ("a", "a")) == [((p, q), "a")]
    assert halves((p, q), p, q, ("a", "b")) == [((p, (5.0, 0.0)), "a"), (((5.0, 0.0), q), "b")]
    assert halves((q, p), p, q, ("a", "b")) == [((q, (5.0, 0.0)), "b"), (((5.0, 0.0), p), "a")]
    # A line beside the bond (of a double bond) is cut level with the middle of the bond.
    assert halves(((1.0, 2.0), (9.0, 2.0)), p, q, ("a", "b"))[0] == (((1.0, 2.0), (5.0, 2.0)), "a")
    # What is left of a line after clipping at its labels, when that lies on one side of the middle.
    assert halves(((6.0, 0.0), (9.0, 0.0)), p, q, ("a", "b")) == [(((6.0, 0.0), (9.0, 0.0)), "b")]
    assert halves(((1.0, 0.0), (4.0, 0.0)), p, q, ("a", "b")) == [(((1.0, 0.0), (4.0, 0.0)), "a")]


def test_a_wedge_and_a_hashed_wedge_are_cut_in_two_colours():
    style = Style()
    centre, far = (0.0, 0.0), (20.0, 0.0)
    seg = (centre, (16.0, 0.0))
    one, = wedge_shapes(3, centre, far, "wedge", seg, style)
    near, rest = wedge_shapes(3, centre, far, "wedge", seg, style, ("a", "b"))
    assert isinstance(one, Wedge) and one.colour == BLACK
    assert (near.colour, rest.colour, near.name, rest.name) == ("a", "b", "Wedge 3", "Wedge 3")
    assert near.bbox()[2] == pytest.approx(10.0) and rest.bbox()[0] == pytest.approx(10.0)
    assert near.points[1] == pytest.approx(rest.points[0]) and near.points[2] == pytest.approx(rest.points[3])
    assert rest.points[1] == pytest.approx(one.points[1]) and rest.points[2] == pytest.approx(one.points[2])
    whole, = wedge_shapes(3, centre, far, "hash", seg, style)
    near, rest = wedge_shapes(3, centre, far, "hash", seg, style, ("a", "b"))
    assert isinstance(near, Hash) and (near.colour, rest.colour) == ("a", "b")
    assert near.segments + rest.segments == whole.segments
    assert max(s[0][0] for s in near.segments) <= 10.0 < min(s[0][0] for s in rest.segments)


def test_the_wedge_to_an_oh_is_black_then_red():
    mol, _, group = drawn("O[C@H](C)CC", BONDS)
    marks = [s for s in walk([group]) if isinstance(s, (Wedge, Hash))]
    assert [s.colour for s in marks] == [BLACK, RED] and len({s.name for s in marks}) == 1
    assert [s.colour for s in shapes("O[C@H](C)CC", ATOMS, (Wedge, Hash))] == [BLACK]


def test_the_wedge_to_a_drawn_hydrogen_follows_the_colour_setting():
    fused = "C1CC[N@@H+]2CCCC[C@@H]2C1"   # the hydrogen of the N is drawn, with a hashed wedge
    marks = shapes(fused, BONDS, (Wedge, Hash))
    assert [s.colour for s in marks if s.name == marks[0].name] == [BLUE, BLACK]
    assert {s.colour for s in shapes(fused, ATOMS, (Wedge, Hash))} == {BLACK}


# --- abbreviations -----------------------------------------------------------------------

def test_the_style_takes_names_all_or_a_list():
    assert Style().abbreviate == () and Style().groups == []
    assert ALL.abbreviate == tuple(ABBREVIATIONS) == ("Ph", "tBu", "CO2H", "OAc")
    assert Style(abbreviate="Ph, CO2H").abbreviate == Style(abbreviate=["Ph", "CO2H", "Ph"]).abbreviate == ("Ph", "CO2H")
    assert Style(abbreviate=("CO2H",)).groups == [("*C(=O)O", "CO2H", "HO2C")]
    assert Style(abbreviate="").abbreviate == ()


@pytest.mark.parametrize("kwargs,says", [
    ({"abbreviate": "Me"}, "'Me'"),
    ({"abbreviate": ("Ph", "Et")}, "'Et'"),
    ({"abbreviate": "X", "abbreviations": {"X": ("C*", "X", "X")}}, "starting with"),
    ({"abbreviate": "X", "abbreviations": {"X": ("X", "X")}}, "starting with"),
])
def test_the_style_refuses_names_and_rows_it_cannot_use(kwargs, says):
    with pytest.raises(ValueError, match=says):
        Style(**kwargs)


@pytest.mark.parametrize("smiles", ["*C*", "*", "*C.C", "*(C)C"])
def test_a_group_must_be_one_group_bonded_to_one_atom(smiles):
    with pytest.raises(ValueError, match="starting with"):
        abbreviate(read("CCO"), [(smiles, "X", "X")])


def test_without_groups_the_molecule_is_given_back():
    mol = read(ASPIRIN)
    assert abbreviate(mol, []) is mol
    plain = read("CC(C)[C@@H]1CC[C@@H](C)C[C@H]1O")
    assert abbreviate(plain, ALL.groups) is plain


def labels(mol):
    return sorted(a.label[0] for a in mol.atoms if a.label)


def test_aspirin_keeps_its_ring_and_two_texts():
    before = read(ASPIRIN)
    mol = abbreviate(before, ALL.groups)
    assert len(before.atoms) == 13 and len(mol.atoms) == 8 and len(mol.bonds) == 8
    assert labels(mol) == ["CO2H", "OAc"]
    assert [a.label for a in mol.atoms if a.label] == [("OAc", "AcO"), ("CO2H", "HO2C")]
    assert mol.rings == [[1, 2, 3, 4, 5, 6]] and mol.aromatic_rings == [0]
    assert all(a.h_count == 0 for a in mol.atoms if a.label)
    assert mol.fragments == [list(range(8))] and mol.smiles == ASPIRIN
    assert len(before.atoms) == 13 and not any(a.label for a in before.atoms)   # the input is unchanged


@pytest.mark.parametrize("smiles,found,atoms", [
    (PPH3, ["Ph", "Ph", "Ph"], 4),
    ("C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=CC=CC=C3", ["Ph", "Ph", "Ph"], 4),   # the rings as Kekule structures
    ("CC(C)(C)c1ccc(O)cc1", ["tBu"], 8),
    ("CCCCCCCCCC(=O)O", ["CO2H"], 10),
    ("OC(=O)CCC(=O)O", ["CO2H", "CO2H"], 4),
    ("O=C([O-])c1ccccc1.[Na+]", ["Ph"], 5),            # the carboxylate has a charge: not CO2H
    ("CC(=O)OC", [], 5),                                  # ester: not CO2H (its OAc would match: see below)
    ("c1ccccc1", [], 6),                                  # no attachment point
    ("Cc1ccccc1C", [], 8),                                # two bonds out of the ring: not Ph
    ("c1ccc(cc1)-c1ccccc1", ["Ph"], 7),                   # biphenyl: one ring stays
    ("c1ccccc1Oc1ccccc1", ["Ph", "Ph"], 3),               # both rings are bonded to the oxygen
    ("CC(C)(C)C(C)(C)C", ["tBu"], 5),                     # one tert-butyl stays
    ("OC(=O)c1ccc(cc1)-c1ccccc1", ["CO2H", "Ph"], 8),
    ("N[C@@H](Cc1ccccc1)C(=O)O", ["CO2H", "Ph"], 5),
])
def test_where_groups_are_found(smiles, found, atoms):
    groups = ALL.groups if smiles != "CC(=O)OC" else Style(abbreviate="CO2H").groups
    mol = abbreviate(read(smiles), groups)
    assert labels(mol) == found and len(mol.atoms) == atoms


def test_the_atoms_of_a_group_have_its_charges_isotopes_and_hydrogens():
    anion = [("*[CH-]C", "X", "X")]
    assert labels(abbreviate(read("CC[CH-]C"), anion)) == ["X"]
    assert labels(abbreviate(read("CC[CH+]C"), anion)) == []
    assert labels(abbreviate(read("CC(C)(C)c1ccc(O)cc1"), ALL.groups)) == ["tBu"]
    assert labels(abbreviate(read("C[13C](C)(C)c1ccc(O)cc1"), ALL.groups)) == []
    assert labels(abbreviate(read("CC(C)(C)c1ccc(O)cc1"), [("*[13C](C)(C)C", "X", "X")])) == []
    assert labels(abbreviate(read("CCC(=O)O"), ALL.groups)) == ["CO2H"]
    assert labels(abbreviate(read("CCC(=O)[O-]"), ALL.groups)) == []
    assert labels(abbreviate(read("CCC(=O)OC"), ALL.groups)) == []
    assert labels(abbreviate(read("CCC(O)O"), ALL.groups)) == []
    assert labels(abbreviate(read("CCC(=O)[O]"), ALL.groups)) == []   # an oxygen without its hydrogen


def test_the_larger_of_two_groups_that_share_atoms_wins():
    acetyl = ("*C(C)=O", "Ac", "Ac")
    for groups in ([acetyl, ABBREVIATIONS["OAc"]], [ABBREVIATIONS["OAc"], acetyl]):
        assert labels(abbreviate(read("CC(=O)Oc1ccccc1C"), groups)) == ["OAc"]
    assert labels(abbreviate(read("CC(=O)Oc1ccccc1C"), [acetyl])) == ["Ac"]
    assert labels(abbreviate(read("CC(=O)c1ccccc1C"), [acetyl, ABBREVIATIONS["OAc"]])) == ["Ac"]
    # Two rows for the same group: the first one in the table is used.
    same = [("*C(=O)O", "CO2H", "HO2C"), ("*C(O)=O", "COOH", "HOOC")]
    assert labels(abbreviate(read("CCC(=O)O"), same)) == ["CO2H"]
    assert labels(abbreviate(read("CCC(=O)O"), same[::-1])) == ["COOH"]


def test_a_group_with_a_stereocentre_or_e_z_is_left_as_it_is():
    group = [("*C(C)(O)F", "X", "X")]
    assert labels(abbreviate(read("CCC(C)(O)F"), group)) == ["X"]
    assert labels(abbreviate(read("CC[C@](C)(O)F"), group)) == []
    propenyl = [("*C=CC", "X", "X")]
    assert labels(abbreviate(read("c1ccccc1C=CC"), propenyl)) == ["X"]
    assert labels(abbreviate(read("c1ccccc1/C=C/C"), propenyl)) == []


def test_what_is_left_is_numbered_again_with_its_e_z_and_stereo():
    mol = abbreviate(read("OC(=O)/C=C/C"), ALL.groups)   # crotonic acid, written from the acid end
    assert [a.element for a in mol.atoms] == ["C", "C", "C", "C"] and labels(mol) == ["CO2H"]
    assert [(b.begin, b.end, b.order) for b in mol.bonds] == [(0, 1, 1), (1, 2, 2), (2, 3, 1)]
    assert mol.bonds[1].cis_trans == (0, 3, "trans")
    p = layout(mol)
    side = lambda a: (p[2][0] - p[1][0]) * (p[a][1] - p[1][1]) - (p[2][1] - p[1][1]) * (p[a][0] - p[1][0])  # noqa: E731
    assert side(0) * side(3) < 0   # E: the two ends on opposite sides of the double bond
    for smiles in ("N[C@@H](Cc1ccccc1)C(=O)O", "N[C@H](Cc1ccccc1)C(=O)O", "O[C@H](C)c1ccccc1",
                   "CC(=O)O[C@@H]1CC[C@@H](C)CC1"):
        mol = abbreviate(read(smiles), ALL.groups)
        points = layout(mol)
        marks = stereo_marks(mol, points)
        centres = [i for i, a in enumerate(mol.atoms) if a.stereocentre]
        assert centres and sorted(m.centre for m in marks) == centres
        for m in marks:
            assert drawn_chirality(mol, points, m.centre, marks) == mol.atoms[m.centre].chirality
    left, right = (abbreviate(read(s), ALL.groups) for s in ("N[C@@H](Cc1ccccc1)C(=O)O", "N[C@H](Cc1ccccc1)C(=O)O"))
    assert [a.chirality for a in left.atoms if a.stereocentre] != [a.chirality for a in right.atoms if a.stereocentre]


def test_the_text_of_a_group_follows_the_side_its_bond_comes_from():
    mol = abbreviate(read("CCC(=O)O"), ALL.groups)   # atoms: C, C, CO2H
    from_left = atom_label(mol, 2, [(0, 0), (1, 0), (2, 0)])
    assert [(r.text, r.kind) for r in from_left.runs] == [("CO", "normal"), ("2", "sub"), ("H", "normal")]
    assert (from_left.symbol_run, from_left.at_end) == (0, False)
    from_right = atom_label(mol, 2, [(4, 0), (3, 0), (2, 0)])
    assert [(r.text, r.kind) for r in from_right.runs] == [("HO", "normal"), ("2", "sub"), ("C", "normal")]
    assert (from_right.symbol_run, from_right.at_end) == (2, True)
    below = atom_label(mol, 2, [(2, -2), (2, -1), (2, 0)])
    assert "".join(r.text for r in below.runs) == "CO2H" and not below.at_end
    # The first letter sits on the atom, or the last one when the bond comes from the right.
    c = FONT.char_width("C") * 14
    box = label_box(from_left, (100, 50), FONT, 14, 0.66)
    assert box[0] == pytest.approx(100 - c / 2)
    assert box[2] - box[0] == pytest.approx(FONT.runs_width(from_left.runs, 14, 0.66))
    box = label_box(from_right, (100, 50), FONT, 14, 0.66)
    assert box[2] == pytest.approx(100 + c / 2)
    assert box[2] - box[0] == pytest.approx(FONT.runs_width(from_right.runs, 14, 0.66))


def test_a_digit_after_a_letter_is_a_subscript():
    assert [(r.text, r.kind) for r in formula_runs("CO2H")] == [("CO", "normal"), ("2", "sub"), ("H", "normal")]
    assert [(r.text, r.kind) for r in formula_runs("tBu")] == [("tBu", "normal")]
    assert [(r.text, r.kind) for r in formula_runs("C12H25")] == [("C", "normal"), ("12", "sub"), ("H", "normal"),
                                                                  ("25", "sub")]


def test_the_texts_of_groups_are_drawn_black_and_named_after_the_group():
    style = Style(colour="atoms_and_bonds", abbreviate="all")
    mol, _, group = drawn(ASPIRIN, style)
    texts = {said(t): t for t in walk([group]) if isinstance(t, Text)}
    assert set(texts) <= {"OAc", "AcO", "CO2H", "HO2C"} and len(texts) == 2
    assert {t.colour for t in texts.values()} == {BLACK}
    assert sorted(t.name for t in texts.values()) == ["Atom CO2H 8", "Atom OAc 1"]
    assert {s.colour for s in walk([group])} == {BLACK}   # the O of OAc does not colour its bond
    assert len(shapes(ASPIRIN, style)) == 9 + 2 + 2      # ring, two bonds, two labels
    assert len(shapes(PPH3, ALL)) == 3 + 4


def test_reactions_and_decks_use_the_settings(tmp_path):
    scheme = draw_reaction("OC(=O)c1ccccc1O.CC(=O)OC(C)=O>>CC(=O)Oc1ccccc1C(=O)O.CC(=O)O", "", "",
                           Style(aromatic="circle", colour="atoms", abbreviate="CO2H,OAc"))
    found = list(walk([scheme]))
    assert sum(isinstance(s, Circle) for s in found) == 2
    assert sum(said(s) in ("CO2H", "HO2C") for s in found if isinstance(s, Text)) == 3
    assert RED in {s.colour for s in found if isinstance(s, Text)}
    deck = Deck(ALL)
    group = deck.add_molecule(PPH3)
    assert sorted(said(s) for s in walk([group]) if isinstance(s, Text)) == ["P", "Ph", "Ph", "Ph"]
    deck.save(tmp_path / "a.pptx")
    assert problems(tmp_path / "a.pptx") == []


# --- the 15 molecules, and the command line ------------------------------------------------------

@pytest.mark.parametrize("way", sorted(WAYS))
def test_the_15_molecules_have_no_overlaps_in_any_way_of_drawing(way):
    style = WAYS[way]
    counts = {"circles": 0, "coloured": 0, "groups": 0}
    for key, smiles in TARGETS.items():
        mol, _, group = drawn(smiles, style)
        placed = group.centred_at(SLIDE_WIDTH / 2, SLIDE_HEIGHT / 2)
        assert drawing_overlaps([placed]) == (0, 0), key
        assert outside([placed], SLIDE_WIDTH, SLIDE_HEIGHT) == 0, key
        found = list(walk([group]))
        counts["circles"] += sum(isinstance(s, Circle) for s in found)
        counts["coloured"] += sum(s.colour != BLACK for s in found)
        counts["groups"] += sum(a.label is not None for a in mol.atoms)
    assert counts == {"circle": {"circles": 9, "coloured": 0, "groups": 0},
                      "atoms": {"circles": 0, "coloured": 38, "groups": 0},
                      "atoms_and_bonds": {"circles": 0, "coloured": 108, "groups": 0},
                      "abbreviate": {"circles": 0, "coloured": 0, "groups": 10},
                      "everything": {"circles": 5, "coloured": 81, "groups": 10}}[way]


def run(args, capsys):
    code = main([str(a) for a in args])
    out = capsys.readouterr()
    return code, out.out, out.err


def slide(path):
    with zipfile.ZipFile(path) as z:
        return z.read("ppt/slides/slide1.xml").decode("utf-8")


def test_the_command_takes_the_three_settings(tmp_path, capsys):
    out = tmp_path / "a.pptx"
    assert run([ASPIRIN, "-o", out], capsys)[0] == 0
    plain = slide(out)
    assert 'prst="ellipse"' not in plain and RED not in plain and "OAc" not in plain
    assert run([ASPIRIN, "--aromatic", "circle", "-o", out], capsys)[0] == 0
    assert slide(out).count('prst="ellipse"') == 1
    for flag in ("--colour", "--color"):
        assert run([ASPIRIN, flag, "atoms_and_bonds", "-o", out], capsys)[0] == 0
        assert RED in slide(out)
    assert run([ASPIRIN, "--abbreviate", "all", "-o", out], capsys)[0] == 0
    assert 'name="Atom OAc 1"' in slide(out) and 'name="Atom CO2H 8"' in slide(out)
    assert run([ASPIRIN, "--abbreviate", "CO2H", "-o", out], capsys)[0] == 0
    assert "Atom OAc" not in slide(out) and "Atom CO2H" in slide(out)
    assert problems(out) == []


@pytest.mark.parametrize("args", [
    ["CCO", "--aromatic", "dots"],
    ["CCO", "--colour", "rainbow"],
    ["CCO", "--abbreviate", "Me"],
])
def test_settings_the_command_does_not_know_stop_with_a_usage_message(args, capsys):
    with pytest.raises(SystemExit) as e:
        main(args)
    assert e.value.code == 2
    assert "usage:" in capsys.readouterr().err
