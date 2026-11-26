"""Tests for reaction schemes (smiles2pptx/reaction.py): reading a reaction SMILES and drawing it.

- reading: three parts split at '>', molecules split at '.', identical SMILES merged into one
  with a count, a ring bond across '.' kept in one molecule, error positions counted in the
  whole reaction
- text by the arrow: a digit right after a Latin letter is a subscript
- drawing: the order of the shapes (molecules, "+", coefficients, the arrow), the centre of
  each molecule on the arrow, the gaps, the arrow length, the agents over the arrow
  and the text above and below it, nested groups, no overlaps and nothing off the slide
"""
import math
import zipfile
import xml.etree.ElementTree as ET

import pytest

from pptx_check import problems
from smiles2pptx.check import drawing_overlaps, outside
from smiles2pptx.chem.smiles import SmilesError
from smiles2pptx.figure import Arrow, Group, Run, Text, walk
from smiles2pptx.fonts import get_font
from smiles2pptx.reaction import condition_runs, draw_reaction, ink_box, read_reaction
from smiles2pptx.style import Style
from smiles2pptx.writers.pptx import SLIDE_HEIGHT, SLIDE_WIDTH, write_pptx

STYLE = Style()
FONT = get_font(STYLE.font)
L = STYLE.bond_length
ASPIRIN = "OC(=O)c1ccccc1O.CC(=O)OC(C)=O>OS(=O)(=O)O>CC(=O)Oc1ccccc1C(=O)O.CC(=O)O"
METHANE = "C.O=O.O=O>>O=C=O.O.O"
ESTER = "CC(=O)O.OCC>>CC(=O)OCC.O"
# Well-known reactions, with the text above and below the arrow.
SCHEMES = [
    (ASPIRIN, "", "85 °C, 15 min"),
    (METHANE, "", ""),
    (ESTER, "H2SO4", "reflux, 2 h"),
    ("CC(=O)O.OCC>OS(=O)(=O)O>CC(=O)OCC.O", "cat.", "reflux, 2 h"),
    ("C=C.[H][H]>[Pd]>CC", "", "rt, 1 atm"),
    ("[H][H].[H][H].O=O>>O.O", "", ""),
    ("N#N.[H][H].[H][H].[H][H]>[Fe]>N.N", "", "400 °C, 200 atm"),
    ("c1ccccc1.BrBr>[Fe](Br)(Br)Br>Brc1ccccc1.Br", "", ""),
    ("OC1OC(CO)C(O)C(O)C1O>>CCO.CCO.O=C=O.O=C=O", "yeast", ""),
    (">>O", "", ""),
    ("O>>", "", ""),
]


def names(group):
    return [c.name for c in group.children]


def row(group):
    """Return the shapes in line with the arrow, in order, and the y of the arrow."""
    arrow = next(c for c in group.children if isinstance(c, Arrow))
    on_line = [c for c in group.children if c.name.split(" ")[0] in ("Reactant", "Product", "Plus", "Coefficient")
               or c is arrow]
    return on_line, arrow.p0[1]


# --- reading ----------------------------------------------------------------------------------

def test_the_same_molecule_is_counted_and_the_parts_are_split():
    r = read_reaction(METHANE)
    assert [(c.smiles, c.count) for c in r.reactants] == [("C", 1), ("O=O", 2)]
    assert r.agents == []
    assert [(c.smiles, c.count) for c in r.products] == [("O=C=O", 1), ("O", 2)]
    for c in r.reactants + r.products:
        assert METHANE[c.position:].startswith(c.smiles)


def test_the_agents_are_the_middle_part():
    r = read_reaction(ASPIRIN)
    assert [c.smiles for c in r.reactants] == ["OC(=O)c1ccccc1O", "CC(=O)OC(C)=O"]
    assert [c.smiles for c in r.agents] == ["OS(=O)(=O)O"]
    assert [c.smiles for c in r.products] == ["CC(=O)Oc1ccccc1C(=O)O", "CC(=O)O"]
    assert all(ASPIRIN[c.position:].startswith(c.smiles) for c in r.reactants + r.agents + r.products)


def test_the_same_text_anywhere_on_a_side_is_counted_in_the_order_first_written():
    r = read_reaction("O.C.O.O>>CO")
    assert [(c.smiles, c.count) for c in r.reactants] == [("O", 3), ("C", 1)]


def test_a_different_text_of_the_same_molecule_is_not_counted_with_it():
    r = read_reaction("CO.OC>>COC.O")
    assert [(c.smiles, c.count) for c in r.reactants] == [("CO", 1), ("OC", 1)]


def test_a_ring_bond_across_a_dot_keeps_the_parts_in_one_molecule():
    r = read_reaction("C1.O.C1>>CC.O")
    assert [(c.smiles, c.count) for c in r.reactants] == [("C1.C1", 1), ("O", 1)]


def test_one_side_may_be_empty():
    assert read_reaction(">>CC").reactants == []
    assert read_reaction("CC>>").products == []


@pytest.mark.parametrize("text, position, reason", [
    ("CCO", 4, "three parts"),          # no '>': reported at the end
    ("C>O", 4, "three parts"),          # one '>'
    ("C>>O>C", 5, "three parts"),       # the third '>'
    ("C>>CQ", 5, "unknown element"),    # not an element; the position counts from the start of the reaction
    ("C(>>C", 2, "unclosed"),           # unclosed '(' in the reactants
    ("C>O(>C", 4, "unclosed"),          # ... and in the agents
    (">O>", 1, "no reactants or products"),  # agents only
    (">>", 1, "no reactants or products"),   # nothing at all
    ("C.>>C", 2, "atom after"),         # '.' with no atom after it
])
def test_a_reaction_that_cannot_be_read_says_where_and_why(text, position, reason):
    with pytest.raises(SmilesError) as e:
        read_reaction(text)
    assert e.value.position == position
    assert e.value.smiles == text
    assert reason in e.value.reason


# --- text by the arrow --------------------------------------------------------------------------

@pytest.mark.parametrize("text, runs", [
    ("H2O", [("H", "normal"), ("2", "sub"), ("O", "normal")]),
    ("H2SO4, Et3N", [("H", "normal"), ("2", "sub"), ("SO", "normal"), ("4", "sub"), (", Et", "normal"),
                     ("3", "sub"), ("N", "normal")]),
    ("C12H22O11", [("C", "normal"), ("12", "sub"), ("H", "normal"), ("22", "sub"), ("O", "normal"),
                   ("11", "sub")]),
    ("85 °C, 15 min", [("85 °C, 15 min", "normal")]),
    ("室温2時間", [("室温2時間", "normal")]),
    ("Pd(PPh3)4", [("Pd(PPh", "normal"), ("3", "sub"), (")4", "normal")]),
    ("", []),
])
def test_a_digit_right_after_a_latin_letter_is_a_subscript(text, runs):
    assert [(r.text, r.kind) for r in condition_runs(text)] == runs


# --- drawing --------------------------------------------------------------------------------------

def test_the_methane_scheme_has_its_shapes_in_order():
    g = draw_reaction(METHANE)
    assert names(g) == ["Reactant 1", "Plus", "Coefficient", "Reactant 2", "Arrow",
                        "Product 1", "Plus", "Coefficient", "Product 2"]
    coefficients = [c for c in g.children if c.name == "Coefficient"]
    assert [c.runs[0].text for c in coefficients] == ["2", "2"]
    assert all(c.size == pytest.approx(STYLE.font_size * STYLE.coefficient_size) for c in coefficients)
    pluses = [c for c in g.children if c.name == "Plus"]
    assert all(c.size == pytest.approx(STYLE.font_size * STYLE.plus_size) for c in pluses)


def test_the_aspirin_scheme_has_the_agent_and_the_text_with_the_arrow():
    g = draw_reaction(ASPIRIN, below="85 °C, 15 min")
    assert names(g) == ["Reactant 1", "Plus", "Reactant 2", "Arrow", "Below arrow", "Reagent 1",
                        "Product 1", "Plus", "Product 2"]


def test_each_molecule_is_a_group_inside_the_scheme_and_the_scheme_starts_at_0_0():
    g = draw_reaction(ASPIRIN, below="85 °C, 15 min", name="Reaction 1")
    assert isinstance(g, Group) and g.name == "Reaction 1"
    groups = [c for c in g.children if isinstance(c, Group)]
    assert [c.name for c in groups] == ["Reactant 1", "Reactant 2", "Reagent 1", "Product 1", "Product 2"]
    assert all(not isinstance(s, Group) for c in groups for s in c.children)
    assert g.bbox()[:2] == pytest.approx((0.0, 0.0))


@pytest.mark.parametrize("smiles, above, below", SCHEMES)
def test_what_is_seen_of_each_molecule_and_sign_is_centred_on_the_arrow(smiles, above, below):
    g = draw_reaction(smiles, above, below)
    shapes, y = row(g)
    for s in shapes:
        if isinstance(s, Arrow):
            continue
        x0, y0, x1, y1 = ink_box([s], FONT)
        assert (y0 + y1) / 2 == pytest.approx(y, abs=1e-6), s.name


@pytest.mark.parametrize("smiles, above, below", SCHEMES)
def test_the_shapes_on_the_line_are_left_to_right_with_the_gap_between(smiles, above, below):
    g = draw_reaction(smiles, above, below)
    shapes, _ = row(g)
    boxes = [ink_box([s], FONT) if not isinstance(s, Arrow) else (s.p0[0], 0, s.p1[0], 0) for s in shapes]
    for a, b in zip(boxes, boxes[1:]):
        assert b[0] - a[2] == pytest.approx(STYLE.reaction_gap * L, abs=1e-6)


def test_the_arrow_is_at_least_its_shortest_length_plus_the_margin():
    arrow = next(c for c in draw_reaction(METHANE).children if isinstance(c, Arrow))
    assert arrow.p1[0] - arrow.p0[0] == pytest.approx((STYLE.arrow_min_length + STYLE.arrow_margin) * L)
    assert arrow.p0[1] == arrow.p1[1] and arrow.p1[0] > arrow.p0[0]
    assert arrow.width == pytest.approx(STYLE.arrow_width * STYLE.line_width)


def test_a_long_text_makes_the_arrow_longer_and_sits_in_its_middle():
    below = "reflux in toluene, 12 h, then rt"
    g = draw_reaction(ESTER, below=below)
    arrow = next(c for c in g.children if isinstance(c, Arrow))
    text = next(c for c in g.children if c.name == "Below arrow")
    width = FONT.runs_width(condition_runs(below), STYLE.font_size, STYLE.script_scale)
    assert text.width == pytest.approx(width)
    assert arrow.p1[0] - arrow.p0[0] == pytest.approx(width + STYLE.arrow_margin * L)
    assert text.x + text.width / 2 == pytest.approx((arrow.p0[0] + arrow.p1[0]) / 2)


def test_the_text_above_sits_just_over_the_arrow_and_the_text_below_just_under_it():
    g = draw_reaction(ESTER, above="H2SO4", below="reflux, 2 h")
    arrow = next(c for c in g.children if isinstance(c, Arrow))
    up = next(c for c in g.children if c.name == "Above arrow")
    down = next(c for c in g.children if c.name == "Below arrow")
    y = arrow.p0[1]
    gap = STYLE.arrow_text_gap * L
    assert up.y + FONT.ascent_height(up.size) == pytest.approx(y - gap)      # baseline
    cap_top = down.y + FONT.ascent_height(down.size) - FONT.capital_height(down.size)
    assert cap_top == pytest.approx(y + gap)
    assert [(r.text, r.kind) for r in up.runs] == [("H", "normal"), ("2", "sub"), ("SO", "normal"), ("4", "sub")]
    assert up.y + up.height < y < down.y


def test_an_agent_sits_over_the_arrow_or_over_the_text_above_it():
    smiles = "CC(=O)O.OCC>OS(=O)(=O)O>CC(=O)OCC.O"
    for above in ("", "cat."):
        g = draw_reaction(smiles, above=above)
        arrow = next(c for c in g.children if isinstance(c, Arrow))
        reagent = next(c for c in g.children if c.name == "Reagent 1")
        x0, _, x1, bottom = ink_box([reagent], FONT)
        if above:
            up = next(c for c in g.children if c.name == "Above arrow")
            under = up.y + FONT.ascent_height(up.size) - FONT.capital_height(up.size)
        else:
            under = arrow.p0[1]
        assert bottom == pytest.approx(under - STYLE.reagent_gap * L)
        assert (x0 + x1) / 2 == pytest.approx((arrow.p0[0] + arrow.p1[0]) / 2)
        assert arrow.p1[0] - arrow.p0[0] >= x1 - x0 + STYLE.arrow_margin * L - 1e-6


def test_two_agents_stand_side_by_side_over_the_arrow():
    g = draw_reaction("CC(=O)O.OCC>OS(=O)(=O)O.c1ccccc1>CC(=O)OCC.O")
    a, b = (ink_box([c], FONT) for c in g.children if c.name.startswith("Reagent"))
    assert b[0] - a[2] == pytest.approx(STYLE.reaction_gap * L)
    assert a[3] == pytest.approx(b[3])


@pytest.mark.parametrize("smiles, above, below", SCHEMES)
def test_the_schemes_have_no_overlaps_and_stay_on_the_slide(smiles, above, below):
    g = draw_reaction(smiles, above, below)
    x0, y0, x1, y1 = g.bbox()
    assert x1 - x0 < SLIDE_WIDTH and y1 - y0 < SLIDE_HEIGHT
    g = g.centred_at(SLIDE_WIDTH / 2, SLIDE_HEIGHT / 2)
    assert drawing_overlaps([g]) == (0, 0)
    assert outside([g], SLIDE_WIDTH, SLIDE_HEIGHT) == 0
    boxes = [c.bbox() for c in g.children]
    for k, a in enumerate(boxes):   # the bounding boxes of the scheme's parts do not meet
        for b in boxes[k + 1:]:
            assert not (a[0] < b[2] and b[0] < a[2] and a[1] < b[3] and b[1] < a[3])


def test_an_arrow_through_a_text_is_counted():
    text = Text(100, 90, 40, 16, [Run("cat.")], 14)
    assert drawing_overlaps([Arrow((80, 100), (160, 100), 1.0)] + [text]) == (1, 0)
    assert drawing_overlaps([Arrow((80, 130), (160, 130), 1.0)] + [text]) == (0, 0)


def test_the_sizes_grow_with_the_scale():
    small, large = (draw_reaction(ESTER, "H2SO4", "reflux, 2 h", style=Style(scale=s)) for s in (1.0, 2.0))
    ws = [g.bbox()[2] - g.bbox()[0] for g in (small, large)]
    hs = [g.bbox()[3] - g.bbox()[1] for g in (small, large)]
    assert ws[1] == pytest.approx(2 * ws[0], rel=1e-6)
    assert hs[1] == pytest.approx(2 * hs[0], rel=1e-6)


def test_a_scheme_is_written_as_nested_groups(tmp_path):
    g = draw_reaction(ASPIRIN, below="85 °C, 15 min", name="Reaction aspirin")
    out = tmp_path / "r.pptx"
    write_pptx(out, [[g.centred_at(SLIDE_WIDTH / 2, SLIDE_HEIGHT / 2)]])
    assert problems(out) == []
    with zipfile.ZipFile(out) as z:
        tree = ET.fromstring(z.read("ppt/slides/slide1.xml"))
    P = "{http://schemas.openxmlformats.org/presentationml/2006/main}"
    tree_root = tree.find(f"{P}cSld/{P}spTree")
    outer = tree_root.findall(f"{P}grpSp")
    assert len(outer) == 1
    assert outer[0].find(f"{P}nvGrpSpPr/{P}cNvPr").get("name") == "Reaction aspirin"
    inner = [x.find(f"{P}nvGrpSpPr/{P}cNvPr").get("name") for x in outer[0].findall(f"{P}grpSp")]
    assert inner == ["Reactant 1", "Reactant 2", "Reagent 1", "Product 1", "Product 2"]
    assert len(outer[0].findall(f"{P}cxnSp")) == 1   # the arrow
