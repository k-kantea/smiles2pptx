"""Command-line interface: a thin wrapper around Deck.

    smiles2pptx "CC(=O)Oc1ccccc1C(=O)O" -o aspirin.pptx
    smiles2pptx --file molecules.txt -o molecules.pptx
    smiles2pptx --reaction "CC(=O)O.OCC>>CC(=O)OCC.O" --above "H2SO4" --below "reflux, 2 h" -o ester.pptx
    smiles2pptx "CC(=O)Oc1ccccc1C(=O)O" --aromatic circle --colour atoms --abbreviate CO2H,OAc

Each SMILES is drawn on a slide of its own, in this order: the SMILES given as arguments,
then the lines of the files, then --reaction. A SMILES that contains '>' is drawn as a reaction.
"""
from __future__ import annotations

import argparse
import sys
import warnings
from typing import List, Optional, Tuple

from . import __version__
from .api import Deck
from .chem.smiles import SmilesError
from .style import ABBREVIATIONS, AROMATIC_STYLES, COLOUR_STYLES, Style

# (SMILES, group name, source for error messages, text above the arrow, text below the arrow)
Entry = Tuple[str, Optional[str], str, str, str]


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="smiles2pptx",
        description="Draw chemical structures from SMILES as editable PowerPoint shapes, one structure per slide.")
    p.add_argument("smiles", nargs="*", help="SMILES to draw (a SMILES containing '>' is drawn as a reaction)")
    p.add_argument("-f", "--file", action="append", default=[], metavar="FILE",
                   help="read SMILES from FILE, one per line, optionally followed by a name; blank lines and "
                        "lines starting with # are ignored (can be repeated)")
    p.add_argument("-r", "--reaction", metavar="SMILES", help="reaction SMILES to draw (reactants>agents>products)")
    p.add_argument("--above", default="", metavar="TEXT", help="text above the reaction arrow")
    p.add_argument("--below", default="", metavar="TEXT", help="text below the reaction arrow")
    p.add_argument("-o", "--output", default="smiles2pptx.pptx", metavar="FILE",
                   help="output file (default: %(default)s)")
    p.add_argument("--scale", type=float, default=Style().scale,
                   help="drawing scale; 1.0 is the ACS Document 1996 size (default: %(default)s)")
    p.add_argument("--aromatic", choices=AROMATIC_STYLES, default=Style().aromatic,
                   help="draw aromatic rings with alternating double bonds or with a circle (default: %(default)s)")
    p.add_argument("--colour", "--color", choices=COLOUR_STYLES, default=Style().colour,
                   help="black, atom labels in element colours, or labels and bond halves in element colours "
                        "(default: %(default)s)")
    p.add_argument("--abbreviate", default="", metavar="NAMES",
                   help=f"draw these groups as text: names separated by commas, or 'all' "
                        f"(available: {', '.join(ABBREVIATIONS)})")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p


def _read_file(path: str) -> List[Entry]:
    entries = []
    with open(path, encoding="utf-8-sig") as f:
        for number, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(None, 1)
            entries.append((parts[0], parts[1] if len(parts) > 1 else None, f"{path}, line {number}", "", ""))
    return entries


def main(argv: Optional[List[str]] = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    if (args.above or args.below) and not args.reaction:
        parser.error("--above and --below can only be used with --reaction")
    try:
        style = Style(scale=args.scale, aromatic=args.aromatic, colour=args.colour, abbreviate=args.abbreviate)
    except ValueError as e:
        parser.error(str(e))
    entries: List[Entry] = [(s, None, "the command line", "", "") for s in args.smiles]
    for path in args.file:
        try:
            entries += _read_file(path)
        except OSError as e:
            print(f"smiles2pptx: cannot read {path}: {e.strerror}", file=sys.stderr)
            return 1
    if args.reaction:
        entries.append((args.reaction, None, "--reaction", args.above, args.below))
    if not entries:
        parser.error("nothing to draw: give a SMILES, --file or --reaction")
    deck = Deck(style)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        for smiles, name, source, above, below in entries:
            try:
                if ">" in smiles:
                    deck.add_reaction(smiles, above, below, **({"name": name} if name else {}))
                else:
                    deck.add_molecule(smiles, **({"name": name} if name else {}))
            except SmilesError as e:
                print(f"smiles2pptx: {source}: {e}", file=sys.stderr)
                return 1
    for w in caught:
        print(f"smiles2pptx: warning: {w.message}", file=sys.stderr)
    try:
        deck.save(args.output)
    except OSError as e:
        print(f"smiles2pptx: cannot write {args.output}: {e.strerror}", file=sys.stderr)
        return 1
    print(f"{args.output}: {len(deck)} slide{'s' if len(deck) > 1 else ''}")
    return 0
