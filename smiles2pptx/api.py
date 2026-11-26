"""Slide decks: Deck puts one molecule or reaction scheme on each slide and writes the
slides to a new pptx file.

    deck = Deck(style=Style(scale=1.4))
    deck.add_molecule("CC(=O)Oc1ccccc1C(=O)O")
    deck.add_reaction("CC(=O)O.OCC>>CC(=O)OCC.O", above="H2SO4", below="reflux, 2 h")
    deck.save("out.pptx")

Each drawing is centred on its slide. A drawing that does not fit inside the margins is
redrawn at a smaller scale so that it fits, and a FitWarning is issued.
"""
from __future__ import annotations

import warnings
from dataclasses import replace
from typing import Callable, List, Optional

from .chem.abbreviate import abbreviate
from .chem.perceive import perceive
from .chem.smiles import read_smiles
from .depict import draw_molecule
from .figure import Group
from .layout import layout
from .reaction import draw_reaction
from .style import Style
from .writers.pptx import SLIDE_HEIGHT, SLIDE_WIDTH, write_pptx

MARGIN = 24.0   # points kept clear at each edge of a slide


class FitWarning(UserWarning):
    """Warning that a drawing was scaled down to fit on its slide."""


class Deck:
    """A presentation with one drawing per slide. `width`, `height` and `margin` are in points."""

    def __init__(self, style: Optional[Style] = None, *, width: float = SLIDE_WIDTH,
                 height: float = SLIDE_HEIGHT, margin: float = MARGIN):
        if margin < 0 or 2 * margin >= min(width, height):
            raise ValueError(f"a margin of {margin} pt leaves no room on a slide of {width} x {height} pt")
        self.style = style or Style()
        self.width, self.height, self.margin = width, height, margin
        self.slides: List[Group] = []

    def __len__(self) -> int:
        return len(self.slides)

    def add_molecule(self, smiles: str, name: str = "Molecule") -> Group:
        """Draw a SMILES (which may have several parts separated by '.') on a new slide.

        Return the drawing as placed on the slide. Raise SmilesError if the SMILES cannot
        be parsed."""
        mol = abbreviate(perceive(read_smiles(smiles)), self.style.groups)
        points = layout(mol, self.style)   # the same at any scale, so a smaller redraw reuses it
        return self._add(lambda style: draw_molecule(mol, points, style, name=name), smiles)

    def add_reaction(self, smiles: str, above: str = "", below: str = "", name: str = "Reaction") -> Group:
        """Draw a reaction SMILES (reactants>agents>products) on a new slide.

        `above` and `below` are the text above and below the arrow. Return the drawing as
        placed on the slide. Raise SmilesError if the SMILES cannot be parsed."""
        return self._add(lambda style: draw_reaction(smiles, above, below, style, name=name), smiles)

    def save(self, path) -> None:
        """Write the slides to a new pptx file (an existing file is replaced)."""
        if not self.slides:
            raise ValueError("deck is empty: add a molecule or reaction first")
        write_pptx(path, [[g] for g in self.slides], font=self.style.font, width=self.width, height=self.height)

    def _add(self, draw: Callable[[Style], Group], smiles: str) -> Group:
        group = draw(self.style)
        room = (self.width - 2 * self.margin, self.height - 2 * self.margin)
        size = _size(group)
        if size[0] > room[0] or size[1] > room[1]:
            # Every length in a drawing is proportional to the scale, so one smaller scale makes
            # it fit; the factor just below 1 keeps rounding from leaving it a hair too large.
            factor = min(room[0] / size[0], room[1] / size[1]) * (1 - 1e-9)
            style = replace(self.style, scale=self.style.scale * factor)
            group = draw(style)
            warnings.warn(f"{smiles}: the drawing ({size[0]:.0f} x {size[1]:.0f} pt) does not fit on the slide "
                          f"({room[0]:.0f} x {room[1]:.0f} pt inside the margins); drawn at scale "
                          f"{style.scale:.2f} instead of {self.style.scale:.2f}", FitWarning, stacklevel=3)
        group = group.centred_at(self.width / 2, self.height / 2)
        self.slides.append(group)
        return group


def _size(group: Group):
    x0, y0, x1, y1 = group.bbox()
    return x1 - x0, y1 - y0
