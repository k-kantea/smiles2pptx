"""Tests for perception against the reference table (tests/data/perceive_reference.json).

The table was generated once with RDKit by tools/make_smiles_reference.py; the tests do
not use RDKit. Each SMILES is compared on its rings (sizes and atoms), the hydrogen count
of every atom, the double bonds of the Kekule structure (how many, and which atoms
they join), its true stereocentres and the E/Z of its double bonds.
"""
import copy
import json
from collections import Counter
from pathlib import Path

import pytest

from smiles2pptx.chem.perceive import PerceiveError, perceive
from smiles2pptx.chem.smiles import read_smiles

DATA = Path(__file__).resolve().parent / "data"
REFERENCE = json.loads((DATA / "perceive_reference.json").read_text(encoding="utf-8"))
GROUPS = ("molecules", "samples", "random", "perception")
ORDER = {"single": 1, "double": 2, "triple": 3, "quadruple": 4}
OTHER = {"cis": "trans", "trans": "cis"}
DIRECTIONAL = ("/", "\\")


def cis_trans(mol, begin, end, ref_begin, ref_end, relation):
    """Return a double bond's E/Z as [lo, hi, ref_lo, ref_hi, relation], where ref_x is the
    lowest-index atom bonded to x other than the partner and relation is "cis" or "trans"
    between ref_lo and ref_hi. `ref_begin` and `ref_end` are the atoms `relation` refers to."""
    refs = {begin: ref_begin, end: ref_end}
    lo, hi = sorted((begin, end))
    out = [lo, hi]
    for atom, partner in ((lo, hi), (hi, lo)):
        lowest = min((k for k in mol.neighbors(atom) if k != partner), default=None)
        if refs[atom] != lowest:
            relation = OTHER[relation]
        out.append(lowest)
    return out + [relation]


def table(smiles):
    """Return the perceived molecule in the format of the reference table."""
    mol = perceive(read_smiles(smiles))
    double = [b for b in mol.bonds if b.order == 2]
    return {"smiles": smiles,
            "rings": sorted(len(r) for r in mol.rings),
            "ring_atoms": sorted({i for r in mol.rings for i in r}),
            "hydrogens": [a.h_count for a in mol.atoms],
            "double_bonds": len(double),
            "double_bond_atoms": sorted({i for b in double for i in (b.begin, b.end)}),
            "stereocentres": [i for i, a in enumerate(mol.atoms) if a.stereocentre],
            "cis_trans": sorted(cis_trans(mol, b.begin, b.end, *b.cis_trans) for b in mol.bonds if b.cis_trans)}


def cases(group):
    return [pytest.param(t, id=t.get("id", t["smiles"])) for t in REFERENCE[group]]


def entries(with_errors=False):
    return [t for g in GROUPS for t in REFERENCE[g] if with_errors or "rdkit_error" not in t]


def check(answer):
    if "rdkit_error" in answer:
        if "kekul" not in answer["rdkit_error"].lower():
            pytest.skip(f"RDKit: {answer['rdkit_error']}")
        with pytest.raises(PerceiveError):
            perceive(read_smiles(answer["smiles"]))
        return
    assert table(answer["smiles"]) == {k: v for k, v in answer.items() if k != "id"}


@pytest.mark.parametrize("answer", cases("molecules"))
def test_test_molecules_match_the_table(answer):
    check(answer)


@pytest.mark.parametrize("answer", cases("samples"))
def test_notation_samples_match_the_table(answer):
    check(answer)


@pytest.mark.parametrize("answer", cases("random"))
def test_reordered_smiles_match_the_table(answer):
    """The same molecules written starting from other atoms and in other orders (by RDKit)."""
    check(answer)


@pytest.mark.parametrize("answer", cases("perception"))
def test_perception_samples_match_the_table(answer):
    """Stereo marks kept and dropped, aromatic rings, valences and rings that cannot be kekulized."""
    check(answer)


def test_the_kekule_structure_is_valid():
    """Aromatic bonds become single or double, other bonds keep their order, and no atom
    takes two double bonds from its aromatic bonds."""
    problems = []
    for t in entries():
        mol = perceive(read_smiles(t["smiles"]))
        from_aromatic = Counter()
        for b in mol.bonds:
            if b.kind == "aromatic":
                if b.order not in (1, 2):
                    problems.append(f"{t['smiles']}: aromatic bond {b.begin}-{b.end} has order {b.order}")
                if b.order == 2:
                    from_aromatic.update((b.begin, b.end))
            elif b.order != ORDER[b.kind]:
                problems.append(f"{t['smiles']}: {b.kind} bond {b.begin}-{b.end} has order {b.order}")
        problems += [f"{t['smiles']}: atom {i} takes {n} double bonds from aromatic bonds"
                     for i, n in sorted(from_aromatic.items()) if n > 1]
    assert not problems, "\n".join(problems[:20])


def snapshot(mol):
    return copy.deepcopy((vars(mol), [vars(a) for a in mol.atoms], [vars(b) for b in mol.bonds]))


def test_perceive_does_not_change_the_molecule_it_is_given():
    changed = []
    for t in entries():
        mol = read_smiles(t["smiles"])
        before = snapshot(mol)
        result = perceive(mol)
        if result is mol or snapshot(mol) != before:
            changed.append(t["smiles"])
    assert not changed, changed[:20]


def marked_on_both_sides(mol, bond):
    """Return True if both atoms of `bond` have another bond written with / or \\."""
    return all(any(b is not bond and b.symbol in DIRECTIONAL for b in mol.bonds if atom in (b.begin, b.end))
               for atom in (bond.begin, bond.end))


def test_the_table_is_not_all_alike():
    """A table where every entry says the same thing would let a wrong perceive pass."""
    found = entries()
    relations = [r[-1] for t in found for r in t["cis_trans"]]
    assert "cis" in relations and "trans" in relations
    assert len({json.dumps(t["stereocentres"]) for t in found}) >= 10
    assert len({json.dumps(t["rings"]) for t in found}) > 3
    assert len({t["double_bonds"] for t in found}) >= 5
    assert {0, 1, 2, 3} <= {h for t in found for h in t["hydrogens"]}
    assert any("kekul" in t["rdkit_error"].lower() for t in entries(with_errors=True) if "rdkit_error" in t)
    dropped_stereo = dropped_cis_trans = 0
    for t in found:
        mol = read_smiles(t["smiles"])
        marked = {i for i, a in enumerate(mol.atoms) if a.chirality in ("@", "@@")}
        dropped_stereo += bool(marked - set(t["stereocentres"]))
        kept = {(r[0], r[1]) for r in t["cis_trans"]}
        dropped_cis_trans += sum(1 for b in mol.bonds if b.kind == "double" and marked_on_both_sides(mol, b)
                                 and (min(b.begin, b.end), max(b.begin, b.end)) not in kept)
    assert dropped_stereo >= 1 and dropped_cis_trans >= 1
