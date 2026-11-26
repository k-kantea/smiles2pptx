"""Reference table generator for wedges: tests/data/stereo_reference.json.

For each SMILES, the table stores the package's layout and, for every stereocentre and every
single bond from it (plus a hydrogen drawn in the widest gap when the centre has one in its
brackets), whether a wedge or a hashed wedge there shows the "@" or "@@" of the SMILES.
RDKit decides: the drawing, with that mark and the marks the package chooses on the other
centres, is written as a molfile, RDKit reads the stereo from the wedges, and its canonical
SMILES must equal that of the input (checked first with the package's own marks). The tests
then check smiles2pptx.depict.stereo against the table, whichever bond the package chooses.

Development tool only: neither the package nor the tests use RDKit.

Usage:
    python tools/make_stereo_reference.py
"""
import json
import math
import sys
from pathlib import Path

from rdkit import Chem, rdBase

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))

from smiles2pptx.chem.perceive import perceive  # noqa: E402
from smiles2pptx.chem.smiles import read_smiles  # noqa: E402
from smiles2pptx.depict.stereo import Mark, stereo_marks  # noqa: E402
from smiles2pptx.layout import layout  # noqa: E402

rdBase.DisableLog("rdApp.*")
DATA = ROOT / "tests" / "data"
EXTRA = [
    "N[C@@H](C)C(=O)O",            # L-alanine
    "C[C@H](N)C(=O)O",             # L-alanine, written from the methyl
    "[C@@H](F)(Cl)Br",             # the centre written first, with a hydrogen
    "F[C@](Cl)(Br)I",              # four atoms, no hydrogen
    "C[S@](=O)c1ccccc1",           # a sulfoxide (a lone pair)
    "[S@@](=O)(C)c1ccccc1",        # a sulfoxide written first
    "C[C@@H]1CCCC[C@H]1C",         # trans-1,2-dimethylcyclohexane
    "C1CC[C@H]2CCCC[C@@H]2C1",     # trans-decalin (ring fusion)
    "C1CC[C@H]2CCCC[C@H]2C1",      # cis-decalin
    "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CC=C4[C@@]3(CC[C@@H](C4)O)C",  # androstenediol (3S,8R,9S,10R,13S,14S,17S)
    "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=CC(=O)CC[C@]34C",         # testosterone (8R,9S,10R,13S,14S,17S)
    "C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc34)[C@@H]1CC[C@@H]2O",            # estradiol
    "CC(C)C[C@H](N)C(=O)O",        # leucine
    "OC[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)C=O",  # glucose, open chain
    "[C@H]([C@H](F)Cl)(F)Cl",      # two centres, each first in the other's list of neighbours
]


def molblock(mol, points, marks):
    hs = [m for m in marks if m.to is None]
    n = len(mol.atoms)
    rows = ["", "  smiles2pptx", "",
            f"{n + len(hs):3d}{len(mol.bonds) + len(hs):3d}  0  0  0  0  0  0  0  0999 V2000"]
    for (x, y), a in zip(points, mol.atoms):
        rows.append(f"{x:10.4f}{y:10.4f}{0:10.4f} {a.element:<3} 0  0  0  0  0  0  0  0  0  0  0  0")
    for m in hs:
        rows.append(f"{m.hydrogen[0]:10.4f}{m.hydrogen[1]:10.4f}{0:10.4f} H   0  0  0  0  0  0  0  0  0  0  0  0")
    on_bond = {frozenset((m.centre, m.to)): m for m in marks if m.to is not None}
    for b in mol.bonds:
        m = on_bond.get(frozenset((b.begin, b.end)))
        if m:
            rows.append(f"{m.centre + 1:3d}{m.to + 1:3d}{b.order:3d}{1 if m.kind == 'wedge' else 6:3d}")
        else:
            rows.append(f"{b.begin + 1:3d}{b.end + 1:3d}{b.order:3d}  0")
    for k, m in enumerate(hs):
        rows.append(f"{m.centre + 1:3d}{n + k + 1:3d}  1{1 if m.kind == 'wedge' else 6:3d}")
    for i, a in enumerate(mol.atoms):
        if a.charge:
            rows.append(f"M  CHG  1 {i + 1:3d} {a.charge:3d}")
    rows.append("M  END")
    return chr(10).join(rows)


def drawing_shows(mol, points, marks, smiles):
    """Return True if RDKit reads the drawing as the same molecule, stereo included."""
    rd = Chem.MolFromMolBlock(molblock(mol, points, marks), removeHs=True)
    return rd is not None and Chem.MolToSmiles(rd) == Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


def widest_gap(points, centre, neighbours):
    c = points[centre]
    angles = sorted(math.atan2(points[j][1] - c[1], points[j][0] - c[0]) for j in neighbours)
    best, at = -1.0, 0.0
    for k, a in enumerate(angles):
        b = angles[(k + 1) % len(angles)] + (2 * math.pi if k == len(angles) - 1 else 0)
        if b - a > best:
            best, at = b - a, (a + b) / 2
    return (c[0] + 0.8 * math.cos(at), c[1] + 0.8 * math.sin(at))


def main():
    from test_overlaps import BRIDGED, CROWDED

    smiles = [m["smiles"] for m in json.loads((DATA / "test_molecules.json").read_text(encoding="utf-8"))]
    smiles += [s for s in list(CROWDED.values()) + list(BRIDGED.values()) if "@" in s] + EXTRA
    rows = []
    for s in dict.fromkeys(smiles):
        mol = perceive(read_smiles(s))
        centres = [i for i, a in enumerate(mol.atoms) if a.stereocentre]
        if not centres:
            continue
        points = layout(mol)
        chosen = stereo_marks(mol, points)
        if not drawing_shows(mol, points, chosen, s):
            raise SystemExit(f"the marks the package chooses do not show {s}")
        entries = []
        for c in centres:
            others = [m for m in chosen if m.centre != c]
            taken = {frozenset((m.centre, m.to)) for m in others if m.to is not None}
            places = [(j, None) for j in mol.neighbors(c)
                      if mol.bond_between(c, j).order == 1 and frozenset((c, j)) not in taken]
            if mol.atoms[c].h_count == 1:
                places.append((None, widest_gap(points, c, mol.neighbors(c))))
            for to, hydrogen in places:
                shown = [k for k in ("wedge", "hash")
                         if drawing_shows(mol, points, others + [Mark(c, k, to, hydrogen)], s)]
                entries.append({"centre": c, "to": to, "hydrogen": hydrogen,
                                "kind": shown[0] if len(shown) == 1 else None})
        rows.append({"smiles": s, "points": [list(p) for p in points], "marks": entries})
    out = DATA / "stereo_reference.json"
    out.write_text(json.dumps(rows, indent=0) + "\n", encoding="utf-8", newline="\n")
    n = sum(len(r["marks"]) for r in rows)
    print(out.name, "molecules", len(rows), "places", n, "unreadable", sum(1 for r in rows for e in r["marks"] if e["kind"] is None))


if __name__ == "__main__":
    main()
