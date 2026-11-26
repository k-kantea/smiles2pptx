"""Tests for bonds that look too short between their atom labels (check.py seen_lengths,
short_lines) and for lengthening such bonds in the layout (layout/lengthen.py).

- the visible length of a bond: the longest of its lines, the two colours of a line added up,
  the centre line of a wedge, the span of the dashes of a hashed wedge, 0 when nothing is drawn
- the 15 molecules in the six styles, and the molecules of the reaction schemes: at least
  Style.shortest_line of every bond is visible
- Br-Br, FeBr3, H2SO4, Cl2, pyridine N-oxide, PPh3 with Ph: short before lengthening
  (shortest_line=0), long enough after it
- only bonds in no ring are lengthened, along their own direction, by moving the side with
  fewer atoms; the other side moves when that side would create an overlap, and neither
  moves when both would
- the 15 molecules get the same layout as without lengthening
- a ring bond between two labels (benzotriazole, tetrazole, the peroxide of artemisinin)
  looks long enough thanks to a narrower margin at its two ends: just the margin it needs,
  and not less than LEAST_MARGIN times the label margin. The other bonds keep the label margin
- the wedge to a hydrogen drawn at a stereocentre looks long enough (N+ at a ring fusion,
  alkaloids, erythromycin), and as long from a labelled centre as from a carbon
"""
import json
import math
from pathlib import Path

import pytest

import smiles2pptx.layout.lengthen as lengthen_module
from smiles2pptx.check import seen_lengths, short_lines
from smiles2pptx.chem.abbreviate import abbreviate
from smiles2pptx.chem.perceive import perceive
from smiles2pptx.chem.smiles import read_smiles
from smiles2pptx.depict import LEAST_MARGIN, draw_molecule
from smiles2pptx.depict.stereo import ring_bonds, stereo_marks
from smiles2pptx.figure import Group, Hash, Line, Text, Wedge, walk
from smiles2pptx.fonts import get_font
from smiles2pptx.layout import layout
from smiles2pptx.layout.lengthen import lengthen
from smiles2pptx.reaction import draw_reaction, read_reaction
from smiles2pptx.style import Style

DATA = Path(__file__).resolve().parent / "data"
TARGETS = {m["id"]: m["smiles"] for m in json.loads((DATA / "test_molecules.json").read_text(encoding="utf-8"))
           if m["target"]}
STYLE = Style()
LEAST = STYLE.shortest_line * STYLE.bond_length
WAYS = {"default": Style(), "circle": Style(aromatic="circle"), "atoms": Style(colour="atoms"),
        "bonds": Style(colour="atoms_and_bonds"), "abbreviate": Style(abbreviate="all"),
        "everything": Style(aromatic="circle", colour="atoms_and_bonds", abbreviate="all")}
REACTIONS = [
    "OC(=O)c1ccccc1O.CC(=O)OC(C)=O>OS(=O)(=O)O>CC(=O)Oc1ccccc1C(=O)O.CC(=O)O",
    "C.O=O.O=O>>O=C=O.O.O",
    "CC(=O)O.OCC>>CC(=O)OCC.O",
    "CC(=O)O.OCC>OS(=O)(=O)O>CC(=O)OCC.O",
    "C=C.[H][H]>[Pd]>CC",
    "N#N.[H][H].[H][H].[H][H]>[Fe]>N.N",
    "c1ccccc1.BrBr>[Fe](Br)(Br)Br>Brc1ccccc1.Br",
    "OC1OC(CO)C(O)C(O)C1O>>CCO.CCO.O=C=O.O=C=O",
    "O=C=O.O=C=O.O>>OC(=O)O.OC(=O)O",
]
# Molecules with bonds between two labelled atoms, outside rings.
SHORT = {"bromine": ("BrBr", ""), "iron bromide": ("[Fe](Br)(Br)Br", ""), "sulfuric acid": ("OS(=O)(=O)O", ""),
         "chlorine": ("ClCl", ""), "pyridine N-oxide": ("[O-][n+]1ccccc1", ""),
         "triphenylphosphine": ("c1ccc(cc1)P(c2ccccc2)c3ccccc3", "all")}

# Molecules with ring bonds between two labelled atoms, and the bonds that look too short with
# the full label margin in the default style (numbered as in the shape names).
RING = {"benzotriazole": ("c1ccc2[nH]nnc2c1", [5]), "tetrazole": ("c1nnn[nH]1", [2, 4]),
        "tetrazolide": ("c1nnn[n-]1", [2, 4]),
        "artemisinin": ("C[C@@H]1CC[C@H]2[C@H](C(=O)O[C@H]3[C@@]24[C@H]1CC[C@](O3)(OO4)C)C", [20])}
# Molecules whose stereocentres have wedges to drawn hydrogens.
DRAWN_H = {
    "strychninium": "O=C1C[C@H]2OCC=C3C[N@@H+]4CC[C@@]56[C@H]4C[C@H]3[C@H]2[C@H]6N1c1ccccc15",
    "crininium": "C1C[N@H+]2CC3=CC4=C(C=C3[C@]15[C@H]2C[C@H](C=C5)O)OCO4",
    "quininium": "COc1ccc2nccc([C@@H](O)[C@@H]3C[C@@H]4CC[N@H+]3C[C@@H]4C=C)c2c1",
    "sparteinium": "C1CC[N@@H+]2C[C@@H]3C[C@H]([C@H]2C1)CN4[C@H]3CCCC4",
    "morphinium": "C[N@H+]1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)O[C@H]3[C@H](C=C4)O",
    "erythromycin": "CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)"
                    "O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)"
                    "(C)O)C)C)O)(C)O",
}


def drawn(smiles, style=STYLE):
    mol = abbreviate(perceive(read_smiles(smiles)), style.groups)
    points = layout(mol, style)
    return mol, points, draw_molecule(mol, points, style)


# --- visible length of a bond --------------------------------------------------------------

def test_what_is_seen_of_a_bond_is_its_longest_line():
    shapes = [
        Line((0, 0), (10, 0), 1, name="Bond 1"),
        # a double bond: the inner line is shorter
        Line((0, 5), (12, 5), 1, name="Bond 2"), Line((2, 7), (8, 7), 1, name="Bond 2"),
        # one line in two colours: 4 + 6
        Line((0, 20), (4, 20), 1, colour="FF0000", name="Bond 3"),
        Line((4, 20), (10, 20), 1, colour="0000FF", name="Bond 3"),
        Text(0, 0, 5, 5, [], 10, name="Atom O 1"),
    ]
    assert seen_lengths(Group(shapes), 4) == pytest.approx([10, 12, 10, 0])
    assert short_lines(Group(shapes), 4, 10.5) == [1, 3, 4]


def test_what_is_seen_of_a_wedge_and_a_hashed_wedge_is_its_centre_line():
    wedge = Wedge([(0, -0.2), (9, -2), (9, 2), (0, 0.2)], name="Wedge 1")
    # a hashed wedge in two colours: two shapes, their dashes 0 to 6 and 8 to 14
    near = Hash([((x, -1), (x, 1)) for x in (0, 2, 4, 6)], 1, name="Hash 2")
    rest = Hash([((x, -1), (x, 1)) for x in (8, 10, 12, 14)], 1, colour="FF0000", name="Hash 2")
    assert seen_lengths(Group([wedge, near, rest]), 2) == pytest.approx([9, 14])


def test_the_wedges_to_drawn_hydrogens_come_after_the_bonds():
    shapes = [Line((0, 0), (10, 0), 1, name="Bond 1"), Wedge([(0, 0), (3, -1), (3, 1), (0, 0)], name="Wedge 3")]
    # one bond; wedge 3 goes to a drawn hydrogen (2 is not drawn)
    assert seen_lengths(Group(shapes), 1) == pytest.approx([10, 3])
    assert short_lines(Group(shapes), 1, 5) == [3]


def test_a_bond_hidden_by_the_texts_is_seen_as_0():
    mol, _, g = drawn("BrBr", Style(shortest_line=0))
    assert len(mol.bonds) == 1
    assert seen_lengths(g, 1)[0] < 0.1 * STYLE.bond_length
    mol, _, g = drawn("[Fe](Br)(Br)Br", Style(shortest_line=0))
    assert sorted(seen_lengths(g, 3))[:2] == [0, 0]


# --- every bond shows at least the minimum -----------------------------------------------------

@pytest.mark.parametrize("way", WAYS)
@pytest.mark.parametrize("key", sorted(TARGETS))
def test_every_bond_of_the_15_molecules_is_seen_long_enough(key, way):
    style = WAYS[way]
    mol, _, g = drawn(TARGETS[key], style)
    assert short_lines(g, len(mol.bonds), style.shortest_line * style.bond_length) == []


@pytest.mark.parametrize("reaction", REACTIONS)
def test_every_bond_of_the_molecules_of_the_reactions_is_seen_long_enough(reaction):
    parts = read_reaction(reaction)
    components = {"Reactant": parts.reactants, "Reagent": parts.agents, "Product": parts.products}
    for style in (WAYS["default"], WAYS["everything"]):
        scheme = draw_reaction(reaction, "", "", style)
        molecules = [g for g in scheme.children if isinstance(g, Group)]
        assert len(molecules) == sum(len(c) for c in components.values())
        for g in molecules:
            kind, number = g.name.split()
            smiles = components[kind][int(number) - 1].smiles
            bonds = len(abbreviate(perceive(read_smiles(smiles)), style.groups).bonds)
            assert short_lines(g, bonds, style.shortest_line * style.bond_length) == [], (reaction, g.name)


@pytest.mark.parametrize("name", SHORT)
def test_bonds_between_two_texts_are_short_before_the_lengthening_and_long_enough_after(name):
    smiles, abbr = SHORT[name]
    mol, _, before = drawn(smiles, Style(abbreviate=abbr, shortest_line=0))
    assert short_lines(before, len(mol.bonds), LEAST) != []
    mol, _, after = drawn(smiles, Style(abbreviate=abbr))
    assert short_lines(after, len(mol.bonds), LEAST) == []


def test_a_bond_leaves_a_text_after_the_h_below_it_too():
    # The symbol (-1, -1, 1, 1) and its H below (-1, 0.9, 1, 3), y down, with the atom at (0, 0).
    parts = [(-1, -1, 1, 1), (-1, 0.9, 1, 3)]
    assert lengthen_module._leave((0, 0), (0, 10), parts) == pytest.approx(0.3)     # down, through the H
    assert lengthen_module._leave((0, 0), (0, -10), parts) == pytest.approx(0.1)    # up, away from it
    assert lengthen_module._leave((0, 0), (10, 0), parts) == pytest.approx(0.1)     # to the side
    # A label below it whose top (2.5) is above the bottom of the H (3): the gap is -0.5.
    assert lengthen_module._gap((0, 0), (0, 3.5), parts, [(-1, 2.5, 1, 4.5)]) == pytest.approx(-0.5)


@pytest.mark.parametrize("scale", [0.7, 1.0, 2.0])
def test_the_least_grows_with_the_scale(scale):
    style = Style(scale=scale)
    mol, _, g = drawn("[Fe](Br)(Br)Br", style)
    seen = seen_lengths(g, len(mol.bonds))
    assert min(seen) >= style.shortest_line * style.bond_length
    assert min(seen) < (style.shortest_line + 0.01) * style.bond_length   # just as long as needed


def test_a_larger_least_lengthens_more():
    for least in (0.3, 0.5):
        style = Style(shortest_line=least)
        mol, points, g = drawn("ClCl", style)
        assert seen_lengths(g, 1)[0] == pytest.approx(least * style.bond_length, abs=1e-4)
    assert math.dist(*drawn("ClCl", Style(shortest_line=0.5))[1]) > math.dist(*drawn("ClCl", Style(shortest_line=0.3))[1])


# --- how the layout lengthens ---------------------------------------------------------------

def test_only_bonds_in_no_ring_are_lengthened_along_their_direction():
    for name in SHORT:
        smiles, abbr = SHORT[name]
        style = Style(abbreviate=abbr)
        mol = abbreviate(perceive(read_smiles(smiles)), style.groups)
        before = layout(mol, Style(abbreviate=abbr, shortest_line=0))
        after = layout(mol, style)
        ring = ring_bonds(mol)
        longer = 0
        for b in mol.bonds:
            d0 = math.dist(before[b.begin], before[b.end])
            d1 = math.dist(after[b.begin], after[b.end])
            if frozenset((b.begin, b.end)) in ring:
                assert d1 == pytest.approx(d0, abs=1e-9), name
            else:
                assert d1 >= d0 - 1e-9, name
                longer += d1 > d0 + 1e-9
            # the direction of every bond is kept
            u0 = ((before[b.end][0] - before[b.begin][0]) / d0, (before[b.end][1] - before[b.begin][1]) / d0)
            u1 = ((after[b.end][0] - after[b.begin][0]) / d1, (after[b.end][1] - after[b.begin][1]) / d1)
            assert u1 == pytest.approx(u0, abs=1e-9), name
        assert longer >= 1, name


def test_the_side_with_fewer_atoms_moves():
    # FeBr3: each Br moves (the Fe with the other two bromines is the larger side), so the Fe stays.
    mol = perceive(read_smiles("[Fe](Br)(Br)Br"))
    style = Style()
    points = {0: (0.0, 0.0), 1: (-0.866, 0.5), 2: (0.866, 0.5), 3: (0.0, -1.0)}
    out = lengthen(mol, points, style)
    assert out[0] == points[0]
    for a in (1, 2, 3):
        assert math.dist(out[0], out[a]) >= 1


class _CountMoved(lengthen_module._Scorer):
    """Scorer that counts each moved atom of `blocked` as an overlap."""
    blocked = set()
    start = {}

    def cost(self, pos, moved, collect=False):
        hard = sum(1 for a in self.blocked if pos[a] != self.start[a])
        return hard, 0.0, []


def test_the_other_side_moves_when_the_smaller_side_would_make_an_overlap(monkeypatch):
    mol = perceive(read_smiles("BrC(C)(C)Br"))   # each C-Br bond splits into one Br and the rest
    points = {0: (-1.0, 0.0), 1: (0.0, 0.0), 2: (0.0, 1.0), 3: (0.0, -1.0), 4: (0.6, 0.0)}
    monkeypatch.setattr(lengthen_module, "_Scorer", _CountMoved)
    _CountMoved.start = dict(points)
    _CountMoved.blocked = {4}
    out = lengthen(mol, points, Style(shortest_line=2.0))
    # Br 0 and C 1 are drawn far enough apart by moving Br 0 (the smaller side); the bond
    # C1-Br4 is too short, and the smaller side (Br 4) may not move: the rest moves instead.
    assert out[4] == points[4]
    assert math.dist(out[1], out[4]) > math.dist(points[1], points[4])


def test_no_side_moves_when_both_would_make_an_overlap(monkeypatch):
    mol = perceive(read_smiles("ClCl"))
    points = {0: (0.0, 0.0), 1: (1.0, 0.0)}
    monkeypatch.setattr(lengthen_module, "_Scorer", _CountMoved)
    _CountMoved.start = dict(points)
    _CountMoved.blocked = {0, 1}
    assert lengthen(mol, points, Style()) == points


def test_no_lengthening_with_a_least_of_0():
    mol = perceive(read_smiles("BrBr"))
    points = {0: (0.0, 0.0), 1: (1.0, 0.0)}
    assert lengthen(mol, points, Style(shortest_line=0)) == points


@pytest.mark.parametrize("key", sorted(TARGETS))
def test_the_15_molecules_are_laid_out_as_without_the_lengthening(key):
    mol = perceive(read_smiles(TARGETS[key]))
    assert layout(mol, Style()) == layout(mol, Style(shortest_line=0))


# --- ring bonds and wedges to drawn hydrogens ------------------------------------------------

def _text_boxes(group, style):
    """Return the bounding boxes of the atom labels (as tall as a capital), without the margin."""
    font = get_font(style.font)
    out = []
    for t in walk([group]):
        if isinstance(t, Text) and t.name.startswith("Atom"):
            baseline = t.y + font.ascent_height(t.size)
            out.append((t.x, baseline - font.capital_height(t.size), t.x + t.width, baseline))
    return out


def _margin_at(p, boxes):
    """Return how far the line end `p` is from the nearest label: the margin widens a box on
    every side, so a line ends this far out in x or in y."""
    return min(max(b[0] - p[0], 0, p[0] - b[2], b[1] - p[1], 0, p[1] - b[3]) for b in boxes)


def _short_with_the_label_margin(smiles, style):
    mol, _, g = drawn(smiles, Style(aromatic=style.aromatic, abbreviate=style.abbreviate, shortest_line=0))
    return short_lines(g, len(mol.bonds), style.shortest_line * style.bond_length)


def test_the_ring_bonds_are_short_with_the_label_margin():
    for name, (smiles, numbers) in RING.items():
        assert _short_with_the_label_margin(smiles, STYLE) == numbers, name
    # With circles, the N=N of benzotriazole is a single line as well, and it is also short.
    assert _short_with_the_label_margin(RING["benzotriazole"][0], Style(aromatic="circle")) == [5, 6]


@pytest.mark.parametrize("way", WAYS)
@pytest.mark.parametrize("name", RING)
def test_a_short_bond_in_a_ring_gets_just_the_margin_it_needs(name, way):
    smiles, _ = RING[name]
    style = WAYS[way]
    short = _short_with_the_label_margin(smiles, style)
    assert short
    least = style.shortest_line * style.bond_length
    mol, _, g = drawn(smiles, style)
    seen = seen_lengths(g, len(mol.bonds))
    for k in short:
        assert least <= seen[k - 1] < least + 0.01, k
    boxes = _text_boxes(g, style)
    m = style.label_margin
    narrowed = 0
    for line in (s for s in walk([g]) if isinstance(s, Line) and s.name.startswith("Bond")):
        k = int(line.name.split()[1])
        for p in (line.p0, line.p1):
            gap = _margin_at(p, boxes)
            if k in short:
                assert gap >= LEAST_MARGIN * m - 1e-6, k
                narrowed += gap < m - 1e-6
            else:
                assert gap >= m - 1e-6, k
    assert narrowed >= 2


def test_the_margin_is_not_narrowed_below_the_least_margin():
    smiles, numbers = RING["benzotriazole"]
    style = Style(shortest_line=0.5)
    mol, _, g = drawn(smiles, style)
    still = short_lines(g, len(mol.bonds), style.shortest_line * style.bond_length)
    assert 5 in still
    boxes = _text_boxes(g, style)
    gaps = [_margin_at(p, boxes) for line in (s for s in walk([g]) if isinstance(s, Line))
            if int(line.name.split()[1]) == 5 for p in (line.p0, line.p1)]
    assert min(gaps) == pytest.approx(LEAST_MARGIN * style.label_margin, abs=1e-6)
    # The line never touches the label: at least a fifth of the margin is left.
    assert min(gaps) >= 0.2 * style.label_margin


@pytest.mark.parametrize("name", DRAWN_H)
def test_the_wedges_to_drawn_hydrogens_are_seen_long_enough(name):
    mol, _, g = drawn(DRAWN_H[name])
    seen = seen_lengths(g, len(mol.bonds))
    assert len(seen) > len(mol.bonds)       # some hydrogens are drawn
    assert min(seen) >= LEAST


@pytest.mark.parametrize("element", ["N@@H+", "P@@H+", "B@@H-", "Si@@H"])
def test_the_wedge_from_a_centre_with_text_is_seen_as_long_as_from_a_carbon(element):
    # A ring fusion with the element at one end and a carbon at the other, both with their H drawn.
    smiles = f"C1CC[{element}]2CCCC[C@@H]2C1"
    mol, points, g = drawn(smiles)
    marks = [m for m in stereo_marks(mol, points) if m.hydrogen is not None]
    assert sorted(mol.atoms[m.centre].element for m in marks) == sorted(["C", element[:element.index("@")]])
    seen = seen_lengths(g, len(mol.bonds))[len(mol.bonds):]
    assert len(seen) == 2 and seen[0] == pytest.approx(seen[1], abs=1e-6)
