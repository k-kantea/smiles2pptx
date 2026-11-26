"""Tests for overlap removal in the layout (layout/untangle.py) and bridged rings (layout/springs.py).

The counts are those of check.py: atoms as points, and the bounding boxes of the atom labels.

- the 15 test molecules and their rewrites by RDKit: both counts are 0
- well-known crowded molecules (all-cis chains, many carboxyl or phosphate groups, a
  macrolide): they have overlaps before untangling and none after, and untangling keeps
  bond lengths, chain angles, regular rings and E/Z
- well-known bridged molecules: no overlaps, bonds within 35 % of the bond length
- the bridged and cage molecules left out of the 15: they are drawn, with no more overlaps
  than measured on 2026-11-20
- every move untangling can make keeps bond lengths, E/Z, straight bonds and ring shapes,
  and a move that bends an angle is taken only when no other move clears as much
"""
import json
import math
from pathlib import Path

import pytest

import smiles2pptx.layout.place as place
from smiles2pptx.check import overlaps
from smiles2pptx.chem.perceive import perceive
from smiles2pptx.chem.smiles import read_smiles
from smiles2pptx.layout import layout
from smiles2pptx.layout.geometry import direction, side
from smiles2pptx.layout.untangle import _Scorer, untangle
from smiles2pptx.style import Style

DATA = Path(__file__).resolve().parent / "data"
MOLECULES = json.loads((DATA / "test_molecules.json").read_text(encoding="utf-8"))
REFERENCE = json.loads((DATA / "perceive_reference.json").read_text(encoding="utf-8"))
TARGET = {m["id"] for m in MOLECULES if m["target"]}
TARGET_SMILES = [m["smiles"] for m in MOLECULES if m["target"]] + [
    r["smiles"] for r in REFERENCE["random"] if r["id"].rsplit("-", 1)[0] in TARGET]

CROWDED = {
    "arachidonic acid": r"CCCCC/C=C\C/C=C\C/C=C\C/C=C\CCCC(=O)O",
    "docosahexaenoic acid": r"CC/C=C\C/C=C\C/C=C\C/C=C\C/C=C\C/C=C\CCC(=O)O",
    "ATP": r"Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O",
    "maleic acid": r"OC(=O)/C=C\C(=O)O",
    "EDTA": r"OC(=O)CN(CC(=O)O)CCN(CC(=O)O)CC(=O)O",
    "citric acid": r"OC(=O)CC(O)(CC(=O)O)C(=O)O",
    "erythromycin": r"CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)"
                    r"O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)"
                    r"(C)O)C)C)O)(C)O",
}
BRIDGED = {
    "norbornane": "C1CC2CCC1C2",
    "bicyclo[2.2.2]octane": "C1CC2CCC1CC2",
    "DABCO": "C1CN2CCN1CC2",
    "tropinone": "CN1C2CCC1CC(=O)C2",
    "atropine": "CN1C2CCC1CC(C2)OC(=O)C(CO)c1ccccc1",
    "quinine": "COc1ccc2nccc([C@@H](O)[C@@H]3C[C@@H]4CCN3C[C@@H]4C=C)c2c1",
    "artemisinin": "C[C@@H]1CC[C@H]2[C@H](C(=O)O[C@H]3[C@@]24[C@H]1CC[C@](O3)(OO4)C)C",
    "strychnine": "O=C1C[C@H]2OCC=C3CN4CC[C@@]56[C@H]4C[C@H]3[C@H]2[C@H]6N1c1ccccc15",
}
# The molecules left out of the 15, and the overlaps (atoms as points) they had on 2026-11-20.
LEFT_OUT = {m["id"]: m["smiles"] for m in MOLECULES if not m["target"]}
LEFT_OUT_OVERLAPS = {"camphor": 1, "adamantane": 2, "morphine": 1, "strychnine": 0, "cubane": 2}


def laid_out(smiles):
    mol = perceive(read_smiles(smiles))
    return mol, layout(mol)


@pytest.mark.parametrize("smiles", TARGET_SMILES)
def test_the_15_molecules_have_no_overlaps(smiles):
    mol, P = laid_out(smiles)
    r = overlaps(mol, P)
    assert (r.points, r.labels) == (0, 0), (r.clash, r.on_bond, r.crossing, r.line_label, r.label_label)


def without_lengthening(monkeypatch):
    # Untangling alone: the bonds that the last step would lengthen (layout/lengthen.py)
    # keep their length, and that step does not clear overlaps on its own.
    monkeypatch.setattr(place, "lengthen", lambda mol, points, style: points)


@pytest.mark.parametrize("name", CROWDED)
def test_crowded_molecules_have_overlaps_before_the_clearing(name, monkeypatch):
    monkeypatch.setattr(place, "untangle", lambda mol, points, style: points)
    without_lengthening(monkeypatch)
    mol, P = laid_out(CROWDED[name])
    r = overlaps(mol, P)
    assert r.points + r.labels > 0


@pytest.mark.parametrize("name", CROWDED)
def test_crowded_molecules_have_no_overlaps_after_the_clearing(name):
    mol, P = laid_out(CROWDED[name])
    r = overlaps(mol, P)
    assert (r.points, r.labels) == (0, 0), (r.clash, r.on_bond, r.crossing, r.line_label, r.label_label)


@pytest.mark.parametrize("name", CROWDED)
def test_the_clearing_keeps_bond_lengths_chain_angles_and_regular_rings(name, monkeypatch):
    without_lengthening(monkeypatch)
    mol, P = laid_out(CROWDED[name])
    for b in mol.bonds:
        assert math.dist(P[b.begin], P[b.end]) == pytest.approx(1, abs=1e-9)
    in_ring = {a for ring in mol.rings for a in ring}
    for a in range(len(mol.atoms)):
        nbrs = mol.neighbors(a)
        if a in in_ring or len(nbrs) < 2:
            continue
        angles = sorted(round(math.degrees(g), 6) for g in _gaps([direction(P[a], P[b]) for b in nbrs]))
        allowed = {2: ([120, 240], [180, 180]), 3: ([120, 120, 120],), 4: ([90, 90, 90, 90],)}[len(nbrs)]
        assert any(angles == pytest.approx(x, abs=1e-6) for x in allowed), (a, angles)
    for ring in mol.rings:
        n = len(ring)
        r = 1 / (2 * math.sin(math.pi / n))
        cx = sum(P[a][0] for a in ring) / n
        cy = sum(P[a][1] for a in ring) / n
        for a in ring:
            assert math.dist(P[a], (cx, cy)) == pytest.approx(r, abs=1e-9)


def test_the_clearing_counts_a_bond_through_the_h_below_a_symbol():
    # Diethylamine with its two bonds from N pointing up: the H goes below the N. The bond C0-C1
    # starts in that H (and misses the N), so untangling has one overlap to clear: the label of
    # the N (atom 2) and that bond. None of the 260 molecules in the checks needs this, so the
    # points are given by hand.
    mol = perceive(read_smiles("CCNCC"))
    pos = {0: (-0.1, -0.7), 1: (-0.87, 0.5), 2: (0.0, 0.0), 3: (0.87, 0.5), 4: (1.74, 0.0)}
    scorer = _Scorer(mol, sorted(pos), Style())
    assert len(scorer.boxes(pos, [2])[2]) == 2
    hard, _, items = scorer.cost(pos, None, collect=True)
    assert (hard, items) == (1, [(2, 0, 1)])


@pytest.mark.parametrize("name", ["arachidonic acid", "docosahexaenoic acid", "maleic acid"])
def test_the_clearing_keeps_the_double_bonds_cis(name):
    mol, P = laid_out(CROWDED[name])
    checked = 0
    for b in mol.bonds:
        if b.cis_trans is None:
            continue
        ref_begin, ref_end, relation = b.cis_trans
        same = side(P[b.begin], P[b.end], P[ref_begin]) == side(P[b.begin], P[b.end], P[ref_end])
        assert same == (relation == "cis")
        checked += 1
    assert checked > 0


@pytest.mark.parametrize("name", BRIDGED)
def test_bridged_molecules_have_no_overlaps(name):
    mol, P = laid_out(BRIDGED[name])
    r = overlaps(mol, P)
    assert (r.points, r.labels) == (0, 0), (r.clash, r.on_bond, r.crossing, r.line_label, r.label_label)


@pytest.mark.parametrize("name", BRIDGED)
def test_bonds_of_bridged_molecules_stay_near_the_bond_length(name):
    mol, P = laid_out(BRIDGED[name])
    worst = max(abs(math.dist(P[b.begin], P[b.end]) - 1) for b in mol.bonds)
    assert worst < 0.35


@pytest.mark.parametrize("name", LEFT_OUT)
def test_left_out_molecules_are_drawn_with_no_more_overlaps_than_before(name):
    mol, P = laid_out(LEFT_OUT[name])
    assert len(P) == len(mol.atoms) and all(math.isfinite(x) and math.isfinite(y) for x, y in P)
    r = overlaps(mol, P)
    assert r.points <= LEFT_OUT_OVERLAPS[name] and r.labels == 0


def test_bridged_layouts_are_the_same_every_time():
    smiles = BRIDGED["quinine"]
    assert laid_out(smiles)[1] == laid_out(smiles)[1]


def _gaps(angles):
    ordered = sorted(a % (2 * math.pi) for a in angles)
    return [(ordered[(i + 1) % len(ordered)] - a) % (2 * math.pi) for i, a in enumerate(ordered)]


# --- every move keeps what it must keep ---------------------------------------------------

MOVE_EXAMPLES = [CROWDED["arachidonic acid"], CROWDED["maleic acid"], CROWDED["citric acid"],
                 r"C/C(C)=C(\C)CC", r"CC/C(C)=C(/CC)C(C)(C)C", r"OC(=O)/C(C)=C/C(=O)O",
                 r"CC(C)(C)/C=C/C(C)(C)C", "CC#CC(C)(C)C", "c1ccccc1C(=O)c1ccccc1"]


@pytest.mark.parametrize("smiles", MOVE_EXAMPLES)
def test_every_move_keeps_bond_lengths_e_z_straight_bonds_and_rings(smiles):
    mol, P = laid_out(smiles)
    atoms = list(range(len(mol.atoms)))
    pos = {a: P[a] for a in atoms}
    scorer = _Scorer(mol, atoms, Style())
    count = 0
    for moved, new, bends in scorer.moves(pos, [tuple(atoms)]):
        trial = dict(pos)
        trial.update(new)
        count += 1
        for b in mol.bonds:
            assert math.dist(trial[b.begin], trial[b.end]) == pytest.approx(1, abs=1e-9)
            if b.cis_trans is not None:
                ref_begin, ref_end, relation = b.cis_trans
                same = side(trial[b.begin], trial[b.end], trial[ref_begin]) == side(trial[b.begin], trial[b.end], trial[ref_end])
                assert same == (relation == "cis"), (moved, bends)
        for a in atoms:
            orders = [b.order for b in mol.bonds if a in (b.begin, b.end)]
            if len(orders) == 2 and (3 in orders or orders.count(2) == 2):
                n1, n2 = mol.neighbors(a)
                angle = abs(math.remainder(direction(trial[a], trial[n1]) - direction(trial[a], trial[n2]), 2 * math.pi))
                assert angle == pytest.approx(math.pi, abs=1e-9)
        for ring in mol.rings:
            for k, a in enumerate(ring):
                b, c = ring[k - 1], ring[(k + 1) % len(ring)]
                before = math.remainder(direction(pos[a], pos[b]) - direction(pos[a], pos[c]), 2 * math.pi)
                after = math.remainder(direction(trial[a], trial[b]) - direction(trial[a], trial[c]), 2 * math.pi)
                assert abs(after) == pytest.approx(abs(before), abs=1e-9)
    assert count > 0


def test_a_flip_is_taken_before_a_bend_that_clears_as_much():
    # Seven carbons on a hexagon: the last one lands on the first. Straightening the chain at
    # the first atom (a bend) and reflecting it in the next bond (a flip) both clear the overlap.
    mol = perceive(read_smiles("CCCCCCC"))
    points = {k: (math.cos(k * math.pi / 3), math.sin(k * math.pi / 3)) for k in range(7)}
    out = untangle(mol, points, Style())
    assert overlaps(mol, [out[a] for a in range(7)]).points == 0
    for a in range(1, 6):
        angle = abs(math.remainder(direction(out[a], out[a - 1]) - direction(out[a], out[a + 1]), 2 * math.pi))
        assert math.degrees(angle) == pytest.approx(120, abs=1e-6)


@pytest.mark.parametrize("name", CROWDED)
def test_crowded_molecules_end_up_as_flat_as_their_turns_allow(name):
    mol, P = laid_out(CROWDED[name])

    def height(turn):
        c, s = math.cos(turn), math.sin(turn)
        ys = [s * x + c * y for x, y in P]
        return max(ys) - min(ys)

    assert height(0) <= min(height(k * math.pi / 6) for k in range(12)) + 1e-6
