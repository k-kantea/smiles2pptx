"""Drawing settings, in one place.

The sizes follow the ACS Document 1996 settings. `scale` multiplies all of them
when a drawing is made, because PowerPoint does not change font sizes or line
widths when a group is enlarged. The default scale of 1.4 gives 20.2 pt bonds
and 14 pt labels.

To add an element colour or an abbreviation, add a row to the tables below.

    Style(aromatic="circle")            # a circle in each aromatic ring
    Style(colour="atoms")               # atom labels in element colours
    Style(colour="atoms_and_bonds")     # and the half of each bond next to such an atom
    Style(abbreviate=("Ph", "CO2H"))    # these groups as text; "all" for the whole table
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Tuple

# Sizes at scale 1.0, in points (ACS Document 1996 settings).
BASE_SIZES: Dict[str, float] = {
    "bond_length": 14.4,   # bond length
    "font_size": 10.0,     # atom labels
    "line_width": 0.6,     # bond lines
    "label_margin": 1.6,   # space between a label and the end of a bond
    "wedge_width": 2.8,    # wide end of a wedge bond
    "hash_spacing": 2.5,   # distance between the short lines of a hashed wedge
}

# An H written above or below the element symbol of an atom with two or more bonds (N over H at
# the bottom corner of a ring): distance between the centres of the symbol and the H / label size.
HYDROGEN_LINE = 1.0

# Colours of atom labels when colour is "atoms" or "atoms_and_bonds" (RGB hex).
# Elements that are not listed are black.
ELEMENT_COLOURS: Dict[str, str] = {
    "N": "1F4FC8", "O": "D02626", "S": "A87C00", "P": "D06A00",
    "F": "1B8A3A", "Cl": "1B8A3A", "Br": "9A3A1A", "I": "7A2FA8",
    "B": "B05A5A", "Si": "6E6A3A",
}

# Substituent abbreviations, used only for the names given in Style.abbreviate:
# name -> (the group as a SMILES with "*" for the atom it is bonded to,
#          label when the bond comes from the left, label when the bond comes from the right).
# In the labels a digit right after a letter is a subscript.
ABBREVIATIONS: Dict[str, Tuple[str, str, str]] = {
    "Ph": ("*c1ccccc1", "Ph", "Ph"),
    "tBu": ("*C(C)(C)C", "tBu", "tBu"),
    "CO2H": ("*C(=O)O", "CO2H", "HO2C"),
    "OAc": ("*OC(C)=O", "OAc", "AcO"),
}

AROMATIC_STYLES = ("kekule", "circle")
COLOUR_STYLES = ("black", "atoms", "atoms_and_bonds")


@dataclass(frozen=True)
class Style:
    """Settings for one drawing. Every value has a default; change only what you need."""

    scale: float = 1.4
    font: str = "Arial"
    bond_angle: float = 120.0       # degrees, between two bonds of a chain
    fragment_gap: float = 1.5       # between the parts of a SMILES (split by "."), atom to atom / bond length
    double_bond_gap: float = 0.18   # distance between the two lines / bond length
    inner_line_trim: float = 0.12   # inner line of a ring double bond, shortened at each end / bond length
    shortest_line: float = 0.25     # shortest visible part of a bond between two labels / bond length;
                                    # a bond in no ring is lengthened to reach it (0: never)
    script_scale: float = 0.66      # subscripts and superscripts / label size
    circle_radius: float = 0.62     # aromatic circle radius / distance from the ring centre to a bond
    circle_min_radius: float = 0.4  # the smallest it may become where a label reaches into the ring
    aromatic: str = "kekule"        # "kekule" or "circle"
    colour: str = "black"           # "black", "atoms" or "atoms_and_bonds"
    element_colours: Mapping[str, str] = field(default_factory=lambda: dict(ELEMENT_COLOURS))
    abbreviate: Tuple[str, ...] = ()  # names in `abbreviations` to draw as text; "all" for every one
    abbreviations: Mapping[str, Tuple[str, str, str]] = field(default_factory=lambda: dict(ABBREVIATIONS))
    # Reaction schemes (multiples of the bond length or the label size).
    reaction_gap: float = 0.55        # space between items / bond length
    arrow_min_length: float = 2.5     # shortest arrow / bond length (before the margin)
    arrow_margin: float = 1.2         # added to the arrow length / bond length
    arrow_width: float = 1.3          # arrow line width / bond line width
    plus_size: float = 1.3            # "+" / label size
    coefficient_size: float = 1.2     # "2" in front of a molecule / label size
    arrow_text_gap: float = 0.43      # arrow to the capitals of the text above and below it / bond length
    reagent_gap: float = 0.35         # structures over the arrow to the arrow or the text under them / bond length

    def __post_init__(self):
        if self.scale <= 0:
            raise ValueError(f"scale must be positive, not {self.scale}")
        if self.aromatic not in AROMATIC_STYLES:
            raise ValueError(f"aromatic must be one of {AROMATIC_STYLES}, not {self.aromatic!r}")
        if self.colour not in COLOUR_STYLES:
            raise ValueError(f"colour must be one of {COLOUR_STYLES}, not {self.colour!r}")
        names = self.abbreviate
        if isinstance(names, str):   # "all", or names separated by commas
            names = [n.strip() for n in names.split(",") if n.strip()]
            if names == ["all"]:
                names = list(self.abbreviations)
        names = tuple(dict.fromkeys(names))
        for name in names:
            if name not in self.abbreviations:
                raise ValueError(f"unknown abbreviation {name!r} "
                                 f"(available: {', '.join(self.abbreviations)})")
            row = self.abbreviations[name]
            if len(row) != 3 or not all(isinstance(v, str) and v for v in row) or not row[0].startswith("*"):
                raise ValueError(f"abbreviation {name!r} must be (SMILES starting with '*', "
                                 f"text for a bond from the left, text for a bond from the right), not {row!r}")
        object.__setattr__(self, "abbreviate", names)

    @property
    def bond_length(self) -> float:
        return BASE_SIZES["bond_length"] * self.scale

    @property
    def font_size(self) -> float:
        return BASE_SIZES["font_size"] * self.scale

    @property
    def line_width(self) -> float:
        return BASE_SIZES["line_width"] * self.scale

    @property
    def label_margin(self) -> float:
        return BASE_SIZES["label_margin"] * self.scale

    @property
    def wedge_width(self) -> float:
        return BASE_SIZES["wedge_width"] * self.scale

    @property
    def hash_spacing(self) -> float:
        return BASE_SIZES["hash_spacing"] * self.scale

    @property
    def groups(self) -> List[Tuple[str, str, str]]:
        """The rows of the abbreviations to draw as text, in the order of `abbreviate`."""
        return [tuple(self.abbreviations[name]) for name in self.abbreviate]

    def element_colour(self, symbol: str) -> str:
        if self.colour == "black":
            return "000000"
        return self.element_colours.get(symbol, "000000")
