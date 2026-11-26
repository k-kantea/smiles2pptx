"""Molecule model: atoms and bonds, as the SMILES reader creates them.

The reader records only what the SMILES says. The hydrogens of atoms outside
brackets, the rings, the Kekule structure and which stereo marks are real are
worked out by the next stage (chem/perceive.py), which fills in the fields marked
"perceive" below. The reader leaves them empty.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Set, Tuple

from .elements import ATOMIC_NUMBER

# Placeholder in Atom.neighbor_order for a hydrogen written inside the brackets of an atom ([C@@H]).
IMPLICIT_H = -1


@dataclass
class Atom:
    element: str                      # "C", "Cl", "Na", or "*" for any atom
    aromatic: bool = False            # written in lower case (c, n, [nH])
    bracket: bool = False             # written inside [ ]
    isotope: Optional[int] = None     # [13C] -> 13
    charge: int = 0                   # [O-] -> -1, [Fe+3] -> 3
    hydrogens: Optional[int] = None   # written in brackets ([CH3] -> 3, [C] -> 0); None outside brackets
    chirality: Optional[str] = None   # "@" or "@@"; a rarely used class is kept as written ("@SP1")
    atom_class: Optional[int] = None  # [CH3:1] -> 1
    position: int = 0                 # index of the atom's first character in the SMILES
    # The neighbours in the order the SMILES gives them, which "@" and "@@" refer to:
    # the preceding atom, then the hydrogens in the brackets (IMPLICIT_H), then the ring
    # bonds, the branches and the next atom, in the order they are written.
    neighbor_order: List[int] = field(default_factory=list)
    # perceive: the number of hydrogens that are not separate atoms (as written in brackets;
    # outside brackets, what the normal valence leaves).
    h_count: Optional[int] = None
    # perceive: True if "@" or "@@" is on a real stereocentre (marks on other atoms are ignored).
    stereocentre: bool = False
    # abbreviate: this atom stands for a group and is drawn as text (Ph, CO2H):
    # (label when the bond comes from the left, label when it comes from the right).
    label: Optional[Tuple[str, str]] = None

    @property
    def atomic_number(self) -> int:
        return ATOMIC_NUMBER[self.element]


@dataclass
class Bond:
    begin: int                        # the atom the bond symbol is written after
    end: int
    kind: str = "single"              # single, double, triple, quadruple or aromatic
    symbol: str = ""                  # as written: "", "-", "=", "#", "$", ":", "/" or "\\"
    ring_closure: bool = False        # made by a ring bond number
    position: int = 0                 # index of the character that made the bond
    # perceive: the bond order in the Kekule structure (aromatic bonds become 1 or 2).
    order: Optional[int] = None
    # perceive: on a double bond with E/Z, (a neighbour of begin, a neighbour of end,
    # "cis" or "trans"); the neighbours are the lowest-numbered atoms on each side.
    cis_trans: Optional[Tuple[int, int, str]] = None

    def other(self, atom: int) -> int:
        return self.end if atom == self.begin else self.begin


@dataclass
class Molecule:
    smiles: str
    atoms: List[Atom]
    bonds: List[Bond]
    fragments: List[List[int]]        # atoms of each disconnected part, in the order written
    # perceive: the smallest set of smallest rings, each as its atoms in order around the ring.
    rings: List[List[int]] = field(default_factory=list)
    # perceive: the rings (indices into `rings`) that are aromatic by the 4n+2 rule, however
    # they are written (c1ccccc1 or C1=CC=CC=C1).
    aromatic_rings: List[int] = field(default_factory=list)

    def neighbors(self, atom: int) -> List[int]:
        return [b.other(atom) for b in self.bonds if atom in (b.begin, b.end)]

    def bond_between(self, a: int, b: int) -> Optional[Bond]:
        for bond in self.bonds:
            if {bond.begin, bond.end} == {a, b}:
                return bond
        return None

    def aromatic_bonds(self) -> Set[int]:
        """Return the bonds (indices into `bonds`) that lie in an aromatic ring."""
        pairs = set()
        for r in self.aromatic_rings:
            ring = self.rings[r]
            pairs.update(frozenset((a, b)) for a, b in zip(ring, ring[1:] + ring[:1]))
        return {k for k, b in enumerate(self.bonds) if frozenset((b.begin, b.end)) in pairs}
