"""Reference table generator for the SMILES reader and perception tests.

  tests/data/smiles_reference.json   what each SMILES states (atoms, bonds, marks as written)
  tests/data/perceive_reference.json what RDKit perceives in it (rings, hydrogens, a Kekule
                                     structure, and which stereo marks are real)

Development tool only: it reads the SMILES with RDKit once. Neither the package nor the
tests use RDKit; the tests compare the reader and perceive() with these tables.

Usage:
    python tools/make_smiles_reference.py
"""
import json
import sys
from pathlib import Path

from rdkit import Chem, rdBase

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "data"
RANDOM_PER_MOLECULE = 10
RANDOM_SEED = 20261120

BOND_KIND = {
    Chem.BondType.SINGLE: "single",
    Chem.BondType.DOUBLE: "double",
    Chem.BondType.TRIPLE: "triple",
    Chem.BondType.QUADRUPLE: "quadruple",
    Chem.BondType.AROMATIC: "aromatic",
}
DIRECTION = {Chem.BondDir.ENDUPRIGHT: "/", Chem.BondDir.ENDDOWNRIGHT: "\\"}
FLIP = {"/": "\\", "\\": "/"}
TETRAHEDRAL = {Chem.ChiralType.CHI_TETRAHEDRAL_CCW: "@", Chem.ChiralType.CHI_TETRAHEDRAL_CW: "@@"}
# The stereo of a double bond, as the relation between its two stereo atoms.
CIS_TRANS = {Chem.BondStereo.STEREOE: "trans", Chem.BondStereo.STEREOTRANS: "trans",
             Chem.BondStereo.STEREOZ: "cis", Chem.BondStereo.STEREOCIS: "cis"}
OTHER = {"cis": "trans", "trans": "cis"}


def as_written(smiles):
    params = Chem.SmilesParserParams()
    params.sanitize = False
    params.removeHs = False
    return Chem.MolFromSmiles(smiles, params)


def odd(order):
    """Return True if sorting `order` takes an odd number of swaps."""
    inversions = sum(1 for i in range(len(order)) for j in range(i + 1, len(order)) if order[i] > order[j])
    return inversions % 2 == 1


def by_index(order, symbol):
    """Return the chirality symbol for the same centre with its neighbours in ascending order
    (a hydrogen, as -1, first)."""
    return ("@@" if symbol == "@" else "@") if odd(order) else symbol


def table(smiles):
    """Return what the SMILES states, in a form that does not depend on how the bonds are stored."""
    mol = as_written(smiles)
    atoms = [[a.GetAtomicNum(), a.GetIsAromatic(), a.GetIsotope() or None, a.GetFormalCharge(),
              a.GetNumExplicitHs() if a.GetNoImplicit() else None, a.GetAtomMapNum() or None]
             for a in mol.GetAtoms()]
    bonds = sorted([min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()), max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
                    BOND_KIND[b.GetBondType()]] for b in mol.GetBonds())
    directional = []
    for b in mol.GetBonds():
        if b.GetBondDir() in DIRECTION:
            i, j, d = b.GetBeginAtomIdx(), b.GetEndAtomIdx(), DIRECTION[b.GetBondDir()]
            directional.append([i, j, d] if i < j else [j, i, FLIP[d]])
    # In RDKit the chiral tag refers to the order of the atom's bonds. Bracket hydrogens
    # become atoms (numbered after the others) so that every neighbour is in that order.
    n = mol.GetNumAtoms()
    mol.UpdatePropertyCache(strict=False)
    with_h = Chem.AddHs(mol, explicitOnly=True)
    chiral = {}
    for a in with_h.GetAtoms():
        if a.GetIdx() < n and a.GetChiralTag() in TETRAHEDRAL:
            order = [b.GetOtherAtomIdx(a.GetIdx()) for b in a.GetBonds()]
            order = [-1 if k >= n else k for k in order]
            chiral[str(a.GetIdx())] = by_index(order, TETRAHEDRAL[a.GetChiralTag()])
    return {"smiles": smiles, "atoms": atoms, "bonds": bonds, "directional": sorted(directional),
            "chiral": chiral, "fragments": len(Chem.GetMolFrags(mol))}


def cis_trans(neighbors, begin, end, ref_begin, ref_end, relation):
    """Return a double bond's E/Z as [lo, hi, ref_lo, ref_hi, relation], where ref_x is the
    lowest-index atom bonded to x other than the partner, and relation ("cis" or "trans") is
    between ref_lo and ref_hi. `ref_begin` and `ref_end` are the atoms the given `relation`
    refers to."""
    refs = {begin: ref_begin, end: ref_end}
    lo, hi = sorted((begin, end))
    out = [lo, hi]
    for atom, partner in ((lo, hi), (hi, lo)):
        lowest = min((k for k in neighbors(atom) if k != partner), default=None)
        if refs[atom] != lowest:
            relation = OTHER[relation]
        out.append(lowest)
    return out + [relation]


def sanitize_error(smiles):
    """Return the first line of RDKit's error message when it cannot sanitize the SMILES."""
    try:
        Chem.SanitizeMol(as_written(smiles))
    except Exception as error:
        return str(error).splitlines()[0]
    raise RuntimeError(f"{smiles}: RDKit did not read it with sanitizing, but sanitizes it alone")


def perceived(smiles):
    """Return what RDKit perceives in the SMILES, with the atoms numbered as written."""
    params = Chem.SmilesParserParams()
    params.removeHs = False  # keep [H] as atoms, so the atom indices follow the SMILES order
    mol = Chem.MolFromSmiles(smiles, params)
    if mol is None:
        return {"smiles": smiles, "rdkit_error": sanitize_error(smiles)}
    assert mol.GetNumAtoms() == as_written(smiles).GetNumAtoms(), smiles

    # The true smallest set of smallest rings (RingInfo holds the symmetrized set, which
    # has more rings than this for cubane).
    rings = [list(r) for r in Chem.GetSSSR(Chem.Mol(mol))]
    ring_atoms = sorted({i for r in rings for i in r})
    assert ring_atoms == [a.GetIdx() for a in mol.GetAtoms() if a.IsInRing()], smiles

    kekule = Chem.Mol(mol)
    Chem.Kekulize(kekule, clearAromaticFlags=True)
    double = [b for b in kekule.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE]

    stereocentres = [a.GetIdx() for a in mol.GetAtoms() if a.GetChiralTag() in TETRAHEDRAL]
    marked = {a.GetIdx() for a in as_written(smiles).GetAtoms() if a.GetChiralTag() in TETRAHEDRAL}
    assert set(stereocentres) <= marked, smiles

    def neighbors(i):
        return [n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors()]

    # Read the same bonds again from the bond directions alone, as cis/trans of the atoms
    # the directions are written on; both readings must agree.
    from_directions = as_written(smiles)
    Chem.SetBondStereoFromDirections(from_directions)
    bonds = []
    for b in mol.GetBonds():
        if b.GetStereo() in CIS_TRANS:
            begin, end = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            entry = cis_trans(neighbors, begin, end, *b.GetStereoAtoms(), CIS_TRANS[b.GetStereo()])
            check = from_directions.GetBondBetweenAtoms(begin, end)
            assert check.GetStereo() in CIS_TRANS, smiles
            a, c = check.GetStereoAtoms()
            if check.GetBeginAtomIdx() != begin:
                a, c = c, a
            assert cis_trans(neighbors, begin, end, a, c, CIS_TRANS[check.GetStereo()]) == entry, smiles
            bonds.append(entry)

    return {"smiles": smiles,
            "rings": sorted(len(r) for r in rings),
            "ring_atoms": ring_atoms,
            "hydrogens": [a.GetTotalNumHs() for a in mol.GetAtoms()],
            "double_bonds": len(double),
            "double_bond_atoms": sorted({i for b in double for i in (b.GetBeginAtomIdx(), b.GetEndAtomIdx())}),
            "stereocentres": stereocentres,
            "cis_trans": sorted(bonds)}


def random_smiles(molecules):
    """Return each test molecule written from other starting atoms and in other atom orders,
    as (id, SMILES) pairs."""
    out = []
    for m in molecules:
        mol = Chem.MolFromSmiles(m["smiles"])
        for k, smiles in enumerate(Chem.MolToRandomSmilesVect(mol, RANDOM_PER_MOLECULE, randomSeed=RANDOM_SEED)):
            out.append((f"{m['id']}-{k}", smiles))
    return out


def write_table(path, out):
    """Write `out` as JSON with one entry per line, so that a change to the table shows up as
    a change to its lines."""
    lines = ["{"]
    keys = list(out)
    for n, key in enumerate(keys):
        value = out[key]
        end = "," if n < len(keys) - 1 else ""
        if isinstance(value, list) and value and isinstance(value[0], dict):
            rows = ",\n".join("    " + json.dumps(v, ensure_ascii=False) for v in value)
            lines.append(f'  "{key}": [\n{rows}\n  ]{end}')
        else:
            lines.append(f'  "{key}": {json.dumps(value, ensure_ascii=False)}{end}')
    lines.append("}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def main():
    molecules = json.loads((DATA / "test_molecules.json").read_text(encoding="utf-8"))
    samples = json.loads((DATA / "smiles_samples.json").read_text(encoding="utf-8"))
    perception = json.loads((DATA / "perceive_samples.json").read_text(encoding="utf-8"))
    random = random_smiles(molecules)
    made_by = f"tools/make_smiles_reference.py with RDKit {rdBase.rdkitVersion}"

    out = {"made_by": made_by,
           "elements": [Chem.GetPeriodicTable().GetElementSymbol(z) for z in range(1, 119)],
           "molecules": [], "samples": [], "random": [], "broken": []}
    for m in molecules:
        out["molecules"].append(dict(table(m["smiles"]), id=m["id"]))
    for s in samples["samples"]:
        out["samples"].append(table(s["smiles"]))
    for name, smiles in random:
        out["random"].append(dict(table(smiles), id=name))
    rdBase.DisableLog("rdApp.error")
    for b in samples["broken"]:
        try:
            rdkit_reads = as_written(b["smiles"]) is not None
        except Exception:
            rdkit_reads = False
        out["broken"].append(dict(b, rdkit_reads=rdkit_reads))
    path = DATA / "smiles_reference.json"
    write_table(path, out)
    chiral = sum(len(t["chiral"]) for g in ("molecules", "samples", "random") for t in out[g])
    directional = sum(len(t["directional"]) for g in ("molecules", "samples", "random") for t in out[g])
    print(f"{path.relative_to(ROOT).as_posix()}: {len(out['molecules'])} molecules, {len(out['samples'])} samples, "
          f"{len(out['random'])} random SMILES, {len(out['broken'])} broken "
          f"({sum(not b['rdkit_reads'] for b in out['broken'])} that RDKit cannot read either); "
          f"{chiral} stereocentres, {directional} bond directions")

    found = {"made_by": made_by,
             "molecules": [dict(perceived(m["smiles"]), id=m["id"]) for m in molecules],
             "samples": [perceived(s["smiles"]) for s in samples["samples"]],
             "random": [dict(perceived(smiles), id=name) for name, smiles in random],
             "perception": [perceived(s["smiles"]) for s in perception["samples"]]}
    path = DATA / "perceive_reference.json"
    write_table(path, found)
    groups = ("molecules", "samples", "random", "perception")
    entries = [t for g in groups for t in found[g]]
    read = [t for t in entries if "rdkit_error" not in t]
    print(f"{path.relative_to(ROOT).as_posix()}: "
          + ", ".join(f"{len(found[g])} {g}" for g in groups) + "; "
          f"{sum(len(t['stereocentres']) for t in read)} stereocentres, "
          f"{sum(len(t['cis_trans']) for t in read)} cis_trans bonds, "
          f"{len(entries) - len(read)} rdkit errors")


if __name__ == "__main__":
    sys.exit(main())
