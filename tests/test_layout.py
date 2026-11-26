"""Tests for the 2D layout, checked against geometric rules.

The bond length is 1 in the layout. The checks:

- bond lengths: within 5 % of 1 for the 15 test molecules and their rewrites by RDKit
- rings: regular polygons, and the rings of a ring system do not fold onto each other
- angles: 120 degrees in chains (the style's bond angle), 180 at triple bonds and
  allenes, 90 at atoms with four neighbours, and ring substituents halfway round the outside
- chains zigzag (they do not curl), and a double bond with / and \\ has its
  substituents on the sides the SMILES specifies
- the disconnected parts of a SMILES are side by side, and every SMILES in the reference tables
  gets coordinates for every atom
"""
import json
import math
from pathlib import Path

import pytest

from smiles2pptx.chem.perceive import perceive
from smiles2pptx.chem.smiles import read_smiles
from smiles2pptx.layout import layout
from smiles2pptx.style import Style

DATA = Path(__file__).resolve().parent / "data"
MOLECULES = json.loads((DATA / "test_molecules.json").read_text(encoding="utf-8"))
REFERENCE = json.loads((DATA / "perceive_reference.json").read_text(encoding="utf-8"))
TARGET = {m["id"] for m in MOLECULES if m["target"]}

# The 15 test molecules, and the SMILES RDKit wrote for them starting from other atoms.
TARGET_SMILES = [m["smiles"] for m in MOLECULES if m["target"]] + [
    r["smiles"] for r in REFERENCE["random"] if r["id"].rsplit("-", 1)[0] in TARGET]
# Every SMILES in the reference tables, including the 5 bridged and cage molecules.
ALL_SMILES = sorted({r["smiles"] for group in ("molecules", "samples", "random", "perception")
                     for r in REFERENCE[group] if "rdkit_error" not in r})

# Fused and spiro rings beyond the test molecules (well-known ring systems).
RING_EXAMPLES = [
    "C1CC1", "C1CCC1", "C1CCCC1", "C1CCCCCCC1",
    "c1ccc2ccccc2c1",                       # naphthalene
    "c1ccc2cc3ccccc3cc2c1",                 # anthracene
    "c1ccc2c(c1)ccc1ccccc12",               # phenanthrene
    "c1ccc2c(c1)c1ccccc1c1ccccc21",         # triphenylene
    "c1cc2ccc3cccc4ccc(c1)c2c34",           # pyrene
    "c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61",  # coronene
    "c1ccc2[nH]ccc2c1",                     # indole
    "c1ccc2cccc-2cc1",                      # azulene
    "C1CCC2CCCCC2C1",                       # decalin
    "C1CC12CC2",                            # spiropentane
    "C1CCC2(CC1)CCCCC2",                    # spiro[5.5]undecane
    "c1ccc(cc1)-c1ccccc1",                  # biphenyl
    "CC1(C)CCCCC1",                         # 1,1-dimethylcyclohexane
    "c1ccccc1C1(C)CCCCC1",                  # 1-methyl-1-phenylcyclohexane
]
# Chains with triple bonds and cumulated double bonds.
STRAIGHT_EXAMPLES = ["CC#CC", "CC#N", "C=C=C", "O=C=O", "CCC#CC#CCC", "CC=C=CC", "c1ccccc1C#N"]
# Double bonds with / and \.
EZ_EXAMPLES = ["C/C=C/C", "C/C=C\\C", "F/C=C/F", "F/C=C\\F", "C/C=C/C=C/C", "C/C=C\\C=C/C",
               "C(/C)=C/C", "C(\\C)=C/C", "OC(=O)/C=C/C(=O)O", "OC(=O)/C=C\\C(=O)O",
               "CC/C(C)=C(\\C)CC", "c1ccccc1/C=C/c1ccccc1", "c1ccccc1/C=C\\c1ccccc1"]


def placed(smiles, style=None):
    mol = perceive(read_smiles(smiles))
    return mol, layout(mol, style)


def angle_at(P, centre, a, b):
    """Return the angle a-centre-b in degrees."""
    u = (P[a][0] - P[centre][0], P[a][1] - P[centre][1])
    v = (P[b][0] - P[centre][0], P[b][1] - P[centre][1])
    c = (u[0] * v[0] + u[1] * v[1]) / (math.hypot(*u) * math.hypot(*v))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def gaps_at(P, centre, neighbours):
    """Return the angles in degrees between successive bonds around an atom, in angular order."""
    angles = sorted(math.atan2(P[n][1] - P[centre][1], P[n][0] - P[centre][0]) for n in neighbours)
    return [math.degrees((angles[(i + 1) % len(angles)] - a) % (2 * math.pi)) or 360.0
            for i, a in enumerate(angles)]


def in_ring(mol):
    return {a for ring in mol.rings for a in ring}


def straight(mol, i):
    orders = [b.order for b in mol.bonds if i in (b.begin, b.end)]
    return len(orders) == 2 and (3 in orders or orders.count(2) == 2)


def inside(point, polygon):
    x, y = point
    result = False
    for (x0, y0), (x1, y1) in zip(polygon, polygon[1:] + polygon[:1]):
        if (y0 > y) != (y1 > y) and x < x0 + (y - y0) * (x1 - x0) / (y1 - y0):
            result = not result
    return result


# --- bond lengths -------------------------------------------------------------------

def test_target_smiles_cover_the_15_molecules_and_their_rewrites():
    assert len(TARGET) == 15
    assert len(TARGET_SMILES) == 15 + 10 * 15


@pytest.mark.parametrize("smiles", TARGET_SMILES)
def test_every_bond_is_within_5_percent_of_the_bond_length(smiles):
    mol, P = placed(smiles)
    for b in mol.bonds:
        assert abs(math.dist(P[b.begin], P[b.end]) - 1) <= 0.05, (b.begin, b.end)


@pytest.mark.parametrize("smiles", RING_EXAMPLES + STRAIGHT_EXAMPLES + EZ_EXAMPLES)
def test_bonds_of_the_examples_have_the_bond_length(smiles):
    mol, P = placed(smiles)
    for b in mol.bonds:
        assert math.dist(P[b.begin], P[b.end]) == pytest.approx(1, abs=1e-9)


@pytest.mark.parametrize("smiles", ALL_SMILES)
def test_every_atom_of_every_answer_table_smiles_gets_a_point(smiles):
    mol, P = placed(smiles)
    assert len(P) == len(mol.atoms)
    assert all(math.isfinite(x) and math.isfinite(y) for x, y in P)


# --- rings ----------------------------------------------------------------------------

@pytest.mark.parametrize("smiles", TARGET_SMILES + RING_EXAMPLES)
def test_rings_are_regular_polygons(smiles):
    mol, P = placed(smiles)
    for ring in mol.rings:
        n = len(ring)
        cx = sum(P[a][0] for a in ring) / n
        cy = sum(P[a][1] for a in ring) / n
        for a in ring:
            assert math.dist(P[a], (cx, cy)) == pytest.approx(1 / (2 * math.sin(math.pi / n)), abs=1e-9)


@pytest.mark.parametrize("smiles", TARGET_SMILES + RING_EXAMPLES)
def test_atoms_of_a_ring_system_are_at_least_one_bond_apart(smiles):
    mol, P = placed(smiles)
    atoms = sorted(in_ring(mol))
    for i in atoms:
        for j in atoms:
            if i < j and mol.bond_between(i, j) is None:
                assert math.dist(P[i], P[j]) >= 1 - 1e-9, (i, j)


@pytest.mark.parametrize("smiles", TARGET_SMILES + RING_EXAMPLES)
def test_substituents_of_a_ring_share_the_room_outside_evenly(smiles):
    """At a ring atom with two ring neighbours, the bonds out of the ring split the outside
    angle into equal parts (one substituent: halfway round)."""
    mol, P = placed(smiles)
    ring_atoms = in_ring(mol)
    for i in ring_atoms:
        nbrs = mol.neighbors(i)
        ring_nbrs = [n for n in nbrs if any(i in r and n in r for r in mol.rings)]
        outside = [n for n in nbrs if n not in ring_nbrs]
        if len(ring_nbrs) == 2 and outside:
            n = len(next(r for r in mol.rings if i in r))
            inner = 180 * (n - 2) / n
            share = (360 - inner) / (len(outside) + 1)
            assert sorted(gaps_at(P, i, nbrs)) == pytest.approx(sorted([inner] + [share] * (len(outside) + 1)),
                                                                abs=1e-6), i
        for n in outside:
            for ring in mol.rings:
                if i in ring:
                    probe = ((P[i][0] + P[n][0]) / 2, (P[i][1] + P[n][1]) / 2)
                    assert not inside(probe, [P[a] for a in ring]), (i, n)


# --- angles in chains -------------------------------------------------------------------

@pytest.mark.parametrize("smiles", TARGET_SMILES + STRAIGHT_EXAMPLES + EZ_EXAMPLES)
def test_angles_at_chain_atoms(smiles):
    mol, P = placed(smiles)
    ring_atoms = in_ring(mol)
    for i in range(len(mol.atoms)):
        nbrs = mol.neighbors(i)
        if i in ring_atoms or len(nbrs) < 2:
            continue
        if len(nbrs) == 2:
            expected = 180 if straight(mol, i) else 120
            assert angle_at(P, i, *nbrs) == pytest.approx(expected, abs=1e-6), i
        elif len(nbrs) in (3, 4):
            assert gaps_at(P, i, nbrs) == pytest.approx([360 / len(nbrs)] * len(nbrs), abs=1e-6), i


def test_the_style_bond_angle_is_used_in_chains():
    mol, P = placed("CCCCCC", Style(bond_angle=110))
    for i in range(1, 5):
        assert angle_at(P, i, i - 1, i + 1) == pytest.approx(110, abs=1e-6)


ZIGZAGS = [
    ("CCCCCCCCCC(=O)O", list(range(10))),
    ("C" * 16, list(range(16))),
    ("OCCOCCOCCO", list(range(10))),
    ("c1ccccc1OCCCCCC", list(range(5, 13))),
    # The longer branch continues the zigzag: 0-1-2-4-5-... with the methyl 3 to the side.
    ("CCC(C)CCCCCC", [0, 1, 2, 4, 5, 6, 7, 8, 9]),
    ("CC(C)(C)CCCCCC", [0, 1, 4, 5, 6, 7, 8, 9]),
]


@pytest.mark.parametrize("smiles,path", ZIGZAGS)
def test_chains_zigzag_and_do_not_curl(smiles, path):
    """Along the path, atoms three bonds apart are sqrt(7) bond lengths apart, as in a zigzag;
    a chain that turns the same way twice puts them 2 apart. (Through an atom with four
    neighbours the chain goes straight on, which also leaves them sqrt(7) apart.)"""
    mol, P = placed(smiles)
    for a, b in zip(path, path[1:]):
        assert mol.bond_between(a, b) is not None
    for a, d in zip(path, path[3:]):
        assert math.dist(P[a], P[d]) == pytest.approx(math.sqrt(7), abs=1e-9), (a, d)


# --- E/Z ---------------------------------------------------------------------------------

def ez_checks():
    """Return (smiles, begin, end, ref_begin, ref_end, relation) for every E/Z double bond
    outside rings in the examples and the reference-table SMILES."""
    checks = []
    for smiles in ALL_SMILES + EZ_EXAMPLES:
        mol = perceive(read_smiles(smiles))
        for b in mol.bonds:
            if b.cis_trans and not any(b.begin in r and b.end in r for r in mol.rings):
                checks.append((smiles, b.begin, b.end) + tuple(b.cis_trans))
    return checks


EZ_CHECKS = ez_checks()


def test_there_are_double_bonds_with_ez_to_check():
    assert len(EZ_CHECKS) >= 30
    assert {c[-1] for c in EZ_CHECKS} == {"cis", "trans"}


@pytest.mark.parametrize("smiles,begin,end,ref_begin,ref_end,relation", EZ_CHECKS)
def test_double_bonds_have_their_substituents_on_the_sides_the_smiles_says(smiles, begin, end, ref_begin,
                                                                          ref_end, relation):
    mol, P = placed(smiles)

    def side(p):
        (x0, y0), (x1, y1) = P[begin], P[end]
        return (x1 - x0) * (p[1] - y0) - (y1 - y0) * (p[0] - x0)

    s1, s2 = side(P[ref_begin]), side(P[ref_end])
    assert abs(s1) > 0.1 and abs(s2) > 0.1
    assert (s1 * s2 > 0) == (relation == "cis")


def test_crotonic_acid_is_drawn_trans_and_isocrotonic_acid_cis():
    for smiles, same_side in (("C/C=C/C(=O)O", False), ("C/C=C\\C(=O)O", True)):
        mol, P = placed(smiles)

        def side(p):
            (x0, y0), (x1, y1) = P[1], P[2]
            return (x1 - x0) * (p[1] - y0) - (y1 - y0) * (p[0] - x0)

        assert (side(P[0]) * side(P[3]) > 0) == same_side


# --- fragments and the whole drawing ---------------------------------------------------------

def test_the_parts_of_a_smiles_are_side_by_side():
    style = Style()
    mol, P = placed("O=C([O-])c1ccccc1.[Na+]", style)
    first, second = mol.fragments
    right_of_first = max(P[a][0] for a in first)
    assert P[second[0]][0] - right_of_first == pytest.approx(style.fragment_gap, abs=1e-9)
    ys = [P[a][1] for a in first]
    assert P[second[0]][1] == pytest.approx((min(ys) + max(ys)) / 2, abs=1e-9)


def test_three_parts_are_in_the_order_written():
    mol, P = placed("[Na+].CCO.[Cl-]")
    xs = [P[part[0]][0] for part in mol.fragments]
    assert xs == sorted(xs)


def test_the_drawing_is_centred_and_a_long_chain_lies_flat():
    mol, P = placed("CCCCCCCCCC(=O)O")
    xs = [p[0] for p in P]
    ys = [p[1] for p in P]
    assert (min(xs) + max(xs)) / 2 == pytest.approx(0, abs=1e-9)
    assert (min(ys) + max(ys)) / 2 == pytest.approx(0, abs=1e-9)
    assert max(xs) - min(xs) > 4 * (max(ys) - min(ys))


@pytest.mark.parametrize("smiles", ["CC(C)(C)c1ccc(O)cc1", "Oc1ccc(O)cc1", "c1ccc(cc1)-c1ccccc1"])
def test_a_long_molecule_is_turned_to_lie_flat(smiles):
    mol, P = placed(smiles)
    xs = [p[0] for p in P]
    ys = [p[1] for p in P]
    assert max(xs) - min(xs) > max(ys) - min(ys)


def test_the_same_smiles_gives_the_same_layout():
    for smiles in TARGET_SMILES[:15]:
        assert placed(smiles)[1] == placed(smiles)[1]


def test_a_molecule_that_was_not_perceived_is_refused():
    with pytest.raises(ValueError, match="perceive"):
        layout(read_smiles("CCO"))


def test_one_atom():
    assert placed("[Na+]")[1] == [(0.0, 0.0)]
