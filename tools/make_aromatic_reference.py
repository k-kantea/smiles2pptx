"""Reference table generator for aromatic rings: tests/data/aromatic_reference.json.

For each SMILES in tests/data/aromatic_samples.json and each test molecule, the table
lists the atoms that RDKit finds aromatic. The tests compare them with the atoms of
Molecule.aromatic_rings. A sample with "expect" in the samples file is one where the
package deliberately differs from RDKit; the table then holds both answers.

Development tool only: it reads the SMILES with RDKit once. Neither the package nor the
tests use RDKit.

Usage:
    python tools/make_aromatic_reference.py
"""
import json
import sys
from pathlib import Path

import rdkit
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "data"


def aromatic_atoms(smiles):
    params = Chem.SmilesParserParams()
    params.removeHs = False   # keep the atom numbering of the SMILES
    mol = Chem.MolFromSmiles(smiles, params)
    if mol is None:
        raise ValueError(f"RDKit cannot read {smiles}")
    return [a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic()]


def main():
    samples = json.loads((DATA / "aromatic_samples.json").read_text(encoding="utf-8"))["samples"]
    molecules = json.loads((DATA / "test_molecules.json").read_text(encoding="utf-8"))
    rows = [{"smiles": m["smiles"], "shows": m["id"]} for m in molecules] + samples
    out = {"made_by": f"tools/make_aromatic_reference.py with RDKit {rdkit.__version__}", "samples": []}
    for row in rows:
        entry = {"smiles": row["smiles"], "shows": row["shows"], "rdkit": aromatic_atoms(row["smiles"])}
        if "expect" in row:
            entry["expect"] = row["expect"]
        out["samples"].append(entry)
    lines = ",\n".join("    " + json.dumps(e, ensure_ascii=False) for e in out["samples"])
    text = '{\n  "made_by": ' + json.dumps(out["made_by"]) + ',\n  "samples": [\n' + lines + "\n  ]\n}\n"
    path = DATA / "aromatic_reference.json"
    path.write_text(text, encoding="utf-8", newline="\n")
    aromatic = sum(1 for e in out["samples"] if e["rdkit"])
    print(f"{path.relative_to(ROOT).as_posix()}: {len(out['samples'])} SMILES, {aromatic} with aromatic atoms, "
          f"{sum('expect' in e for e in out['samples'])} where the package differs from RDKit on purpose")


if __name__ == "__main__":
    sys.exit(main())
