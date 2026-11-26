"""Draw chemical structures from SMILES as editable PowerPoint shapes.

The package uses only the Python standard library.

    from smiles2pptx import Deck
    deck = Deck()
    deck.add_molecule("CC(=O)Oc1ccccc1C(=O)O")
    deck.save("aspirin.pptx")
"""

__version__ = "0.1.0"

from .api import Deck, FitWarning  # noqa: E402
from .chem.smiles import SmilesError  # noqa: E402
from .style import Style  # noqa: E402

__all__ = ["Deck", "FitWarning", "SmilesError", "Style", "__version__"]
