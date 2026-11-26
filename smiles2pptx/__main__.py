"""Entry point for `python -m smiles2pptx`, which does the same as the smiles2pptx command."""
from .cli import main

if __name__ == "__main__":
    raise SystemExit(main())
