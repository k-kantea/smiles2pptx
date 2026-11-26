"""Tests for the SMILES reader against the reference table (tests/data/smiles_reference.json).

The table was generated once with RDKit by tools/make_smiles_reference.py; the tests do
not use RDKit. Each SMILES is compared on its atoms (element, aromatic, isotope,
charge, hydrogens in brackets, atom class), its bonds, the bond directions / and \\,
the stereocentres and the number of disconnected parts.
"""
import json
from pathlib import Path

import pytest

from smiles2pptx.chem.elements import ATOMIC_NUMBER, SYMBOLS
from smiles2pptx.chem.molecule import IMPLICIT_H
from smiles2pptx.chem.smiles import SmilesError, read_smiles

DATA = Path(__file__).resolve().parent / "data"
REFERENCE = json.loads((DATA / "smiles_reference.json").read_text(encoding="utf-8"))
FLIP = {"/": "\\", "\\": "/"}


def odd(order):
    """Return True if sorting `order` takes an odd number of swaps."""
    inversions = sum(1 for i in range(len(order)) for j in range(i + 1, len(order)) if order[i] > order[j])
    return inversions % 2 == 1


def by_index(order, symbol):
    """Return the equivalent chirality symbol with the neighbours in ascending order (hydrogen first)."""
    return ("@@" if symbol == "@" else "@") if odd(order) else symbol


def table(smiles):
    """Return the reader's molecule in the format of the reference table."""
    mol = read_smiles(smiles)
    atoms = [[a.atomic_number, a.aromatic, a.isotope, a.charge, a.hydrogens, a.atom_class] for a in mol.atoms]
    bonds = sorted([min(b.begin, b.end), max(b.begin, b.end), b.kind] for b in mol.bonds)
    directional = sorted([b.begin, b.end, b.symbol] if b.begin < b.end else [b.end, b.begin, FLIP[b.symbol]]
                         for b in mol.bonds if b.symbol in FLIP)
    chiral = {str(i): by_index(a.neighbor_order, a.chirality)
              for i, a in enumerate(mol.atoms) if a.chirality in ("@", "@@")}
    return {"smiles": smiles, "atoms": atoms, "bonds": bonds, "directional": directional,
            "chiral": chiral, "fragments": len(mol.fragments)}


def cases(group):
    return [pytest.param(t, id=t.get("id", t["smiles"])) for t in REFERENCE[group]]


@pytest.mark.parametrize("answer", cases("molecules"))
def test_test_molecules_match_the_table(answer):
    assert table(answer["smiles"]) == {k: v for k, v in answer.items() if k != "id"}


@pytest.mark.parametrize("answer", cases("samples"))
def test_notation_samples_match_the_table(answer):
    assert table(answer["smiles"]) == answer


@pytest.mark.parametrize("answer", cases("random"))
def test_reordered_smiles_match_the_table(answer):
    """The same molecules written starting from other atoms and in other orders (by RDKit)."""
    assert table(answer["smiles"]) == {k: v for k, v in answer.items() if k != "id"}


@pytest.mark.parametrize("broken", [pytest.param(b, id=b["problem"]) for b in REFERENCE["broken"]])
def test_broken_smiles_say_where(broken):
    with pytest.raises(SmilesError) as error:
        read_smiles(broken["smiles"])
    assert error.value.position == broken["position"], str(error.value)
    lines = str(error.value).splitlines()
    assert lines[0].startswith(f"cannot parse SMILES at position {broken['position']}: ")
    assert lines[2].index("^") - 2 == broken["position"] - 1  # the caret is under that character


def test_the_table_is_not_all_alike():
    """A table where every entry says the same thing would let a wrong reader pass."""
    entries = REFERENCE["molecules"] + REFERENCE["samples"] + REFERENCE["random"]
    symbols = [s for t in entries for s in t["chiral"].values()]
    directions = [d for t in entries for _, _, d in t["directional"]]
    kinds = {k for t in entries for _, _, k in t["bonds"]}
    assert len({json.dumps(t["atoms"]) for t in entries}) >= 40
    assert symbols.count("@") >= 50 and symbols.count("@@") >= 50
    assert directions.count("/") >= 5 and directions.count("\\") >= 5
    assert kinds == {"single", "double", "triple", "quadruple", "aromatic"}
    assert len({b["position"] for b in REFERENCE["broken"]}) >= 8
    assert sum(not b["rdkit_reads"] for b in REFERENCE["broken"]) >= 10


def test_element_symbols_match_rdkit():
    assert list(SYMBOLS[1:]) == REFERENCE["elements"]
    assert ATOMIC_NUMBER["*"] == 0


@pytest.mark.parametrize("order, symbol, expected", [
    ([0, 1, 2, 3], "@", "@"),
    ([1, 0, 2, 3], "@", "@@"),     # one swap
    ([1, 2, 0, 3], "@@", "@@"),    # three in a cycle: two swaps
    ([3, 0, 1, 2], "@", "@@"),     # four in a cycle: three swaps
    ([0, IMPLICIT_H, 2, 3], "@", "@@"),
])
def test_chirality_by_index(order, symbol, expected):
    assert by_index(order, symbol) == expected


@pytest.mark.parametrize("smiles, atom, order", [
    ("F[C@@H](Cl)Br", 1, [0, IMPLICIT_H, 2, 3]),     # the hydrogen comes right after the preceding atom
    ("[C@@H](F)(Cl)Br", 0, [IMPLICIT_H, 1, 2, 3]),   # no preceding atom: the hydrogen comes first
    ("[C@@H]1(O)CCCCC1", 0, [IMPLICIT_H, 6, 1, 2]),  # the ring closure takes the position of its number
    ("O[C@@H]1CCCCC1", 1, [0, IMPLICIT_H, 6, 2]),
    ("C1CC[C@H](C1)O", 3, [2, IMPLICIT_H, 4, 5]),
])
def test_neighbor_order(smiles, atom, order):
    assert read_smiles(smiles).atoms[atom].neighbor_order == order


def test_rare_stereo_classes_are_read():
    assert read_smiles("F[C@TH1](Cl)(Br)I").atoms[1].chirality == "@"
    assert read_smiles("F[C@TH2](Cl)(Br)I").atoms[1].chirality == "@@"
    assert read_smiles("[Pt@SP1](Cl)(Cl)(N)N").atoms[0].chirality == "@SP1"
    assert read_smiles("[Fe@OH12](F)(F)(F)(F)(F)F").atoms[0].chirality == "@OH12"


def test_ring_bond_numbers_in_all_forms():
    for smiles in ("C1CCCCC1", "C%12CCCCC%12", "C%(123)CCCCC%(123)"):
        mol = read_smiles(smiles)
        assert len(mol.bonds) == 6 and mol.bond_between(0, 5).ring_closure


def test_a_ring_bond_symbol_at_the_closing_end_begins_there():
    bond = read_smiles("C1CCCCCC/C=C/1").bond_between(0, 8)
    assert (bond.begin, bond.end, bond.symbol) == (8, 0, "/")


def test_parts_and_positions():
    mol = read_smiles("CC(=O)[O-].[Na+]")
    assert mol.fragments == [[0, 1, 2, 3], [4]]
    assert [a.position for a in mol.atoms] == [0, 1, 4, 6, 11]
    assert mol.neighbors(1) == [0, 2, 3]
