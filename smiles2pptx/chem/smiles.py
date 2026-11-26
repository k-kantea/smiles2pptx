"""SMILES reader: turns a SMILES string into a molecule (the first stage).

The reader follows OpenSMILES: atoms in and out of brackets (isotope, element,
chirality, hydrogen count, charge, atom class), the bond symbols - = # $ : / \\,
branches, ring bond numbers (1-9, %nn and %(n...)), lower-case aromatic atoms and
parts separated by '.'. The rarely used chirality classes (@AL, @SP, @TB, @OH)
are read and kept as written; @TH1 and @TH2 are the same as @ and @@.

A SMILES that cannot be parsed raises SmilesError, which gives the position of the
offending character.
"""
from __future__ import annotations

from typing import Dict, List, NoReturn, Optional, Tuple

from .elements import ATOMIC_NUMBER
from .molecule import IMPLICIT_H, Atom, Bond, Molecule

ORGANIC = {"B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"}
AROMATIC_ORGANIC = {"b", "c", "n", "o", "p", "s"}
AROMATIC_IN_BRACKETS = {"b", "c", "n", "o", "p", "s", "se", "as", "te", "si"}
BOND_KIND = {"-": "single", "=": "double", "#": "triple", "$": "quadruple", ":": "aromatic",
             "/": "single", "\\": "single"}
DIRECTIONAL = {"/", "\\"}
FLIP = {"/": "\\", "\\": "/"}
RARE_CHIRALITY = ("TH", "AL", "SP", "TB", "OH")
_OPEN_RING = -2  # placeholder in Atom.neighbor_order for a ring bond that is not closed yet


class SmilesError(ValueError):
    """Error raised for a SMILES that cannot be parsed. `position` counts characters from 1."""

    def __init__(self, smiles: str, index: int, reason: str):
        self.smiles = smiles
        self.index = index
        self.position = index + 1
        self.reason = reason
        super().__init__(f"cannot parse SMILES at position {self.position}: {reason}\n"
                         f"  {smiles}\n  {' ' * index}^")


def read_smiles(smiles: str) -> Molecule:
    """Read one SMILES (which may have several parts separated by '.')."""
    return _Reader(smiles).read()


class _Reader:
    def __init__(self, text: str):
        self.s = text
        self.i = 0
        self.atoms: List[Atom] = []
        self.bonds: List[Bond] = []
        self.prev: Optional[int] = None          # the atom that the next atom or ring bond joins
        self.bond: Optional[Tuple[str, int]] = None  # a pending bond symbol and its index
        self.before_bond: Optional[str] = None   # the token before that bond symbol
        self.branches: List[Tuple[int, int]] = []    # (atom before "(", index of "(")
        # ring bond number -> (atom, bond symbol, index of the symbol or number, slot in neighbor_order)
        self.rings: Dict[int, Tuple[int, str, int, int]] = {}
        self.last: Optional[str] = None          # last token: atom, ring, bond, open, close or dot
        self.last_index = 0
        self.pairs = set()

    def fail(self, index: int, reason: str) -> NoReturn:
        raise SmilesError(self.s, index, reason)

    # --- the whole string -------------------------------------------------------

    def read(self) -> Molecule:
        s = self.s
        while self.i < len(s):
            ch = s[self.i]
            if ch == "[" or ch == "*" or ch.isalpha():
                self.atom()
            elif ch in BOND_KIND:
                self.bond_symbol()
            elif ch == "(":
                self.open_branch()
            elif ch == ")":
                self.close_branch()
            elif ch.isdigit() or ch == "%":
                self.ring_bond()
            elif ch == ".":
                self.dot()
            elif ch.isspace():
                self.fail(self.i, "spaces are not allowed in SMILES")
            else:
                self.fail(self.i, f"invalid character '{ch}'")
        self.finish()
        return Molecule(s, self.atoms, self.bonds, self.fragments())

    def finish(self):
        if not self.atoms:
            self.fail(0, "empty SMILES")
        if self.last == "bond":
            self.fail(self.bond[1], f"missing atom after bond '{self.bond[0]}'")
        if self.last == "dot":
            self.fail(self.last_index, "missing atom after '.'")
        if self.branches:
            self.fail(self.branches[-1][1], "unclosed '('")
        if self.rings:
            number, (_, _, index, _) = min(self.rings.items(), key=lambda item: item[1][2])
            self.fail(index, f"unclosed ring bond {number}")

    def mark(self, kind: str, index: int):
        self.last = kind
        self.last_index = index

    # --- tokens -------------------------------------------------------------------

    def atom(self):
        start = self.i
        atom = self.bracket_atom() if self.s[start] == "[" else self.plain_atom()
        atom.position = start
        index = len(self.atoms)
        self.atoms.append(atom)
        if self.prev is not None:
            symbol = self.bond[0] if self.bond else ""
            self.add_bond(self.prev, index, symbol, start, ring_closure=False)
            self.atoms[self.prev].neighbor_order.append(index)
            atom.neighbor_order.append(self.prev)
        atom.neighbor_order.extend([IMPLICIT_H] * (atom.hydrogens or 0))
        self.prev = index
        self.bond = None
        self.mark("atom", start)

    def plain_atom(self) -> Atom:
        s, i = self.s, self.i
        two = s[i:i + 2]
        if two in ("Cl", "Br"):
            self.i += 2
            return Atom(two)
        ch = s[i]
        if ch == "*":
            self.i += 1
            return Atom("*")
        if ch in ORGANIC:
            self.i += 1
            return Atom(ch)
        if ch in AROMATIC_ORGANIC:
            self.i += 1
            return Atom(ch.upper(), aromatic=True)
        if ch.islower() and i > 0 and s[i - 1:i + 1] in ATOMIC_NUMBER:
            self.fail(i - 1, f"{s[i - 1:i + 1]} must be written in brackets, as [{s[i - 1:i + 1]}]")
        if two in ATOMIC_NUMBER and two[1:].islower():
            self.fail(i, f"{two} must be written in brackets, as [{two}]")
        if ch in ATOMIC_NUMBER:
            self.fail(i, f"{ch} must be written in brackets, as [{ch}]")
        self.fail(i, f"unknown element symbol '{ch}'")

    def bracket_atom(self) -> Atom:
        s, start = self.s, self.i
        j = start + 1

        def at(k):
            if k >= len(s):
                self.fail(start, "unclosed '['")
            return s[k]

        def number(k):
            end = k
            while end < len(s) and s[end].isdigit():
                end += 1
            return (int(s[k:end]) if end > k else None), end

        isotope, j = number(j)
        ch = at(j)
        if ch == "*":
            element, aromatic, j = "*", False, j + 1
        elif ch.isupper():
            two = s[j:j + 2]
            if len(two) == 2 and two[1].islower() and two in ATOMIC_NUMBER:
                element, j = two, j + 2
            elif ch in ATOMIC_NUMBER:
                element, j = ch, j + 1
            else:
                written = two if two[1:].islower() else ch
                self.fail(j, f"unknown element symbol '{written}'")
            aromatic = False
        elif ch.islower():
            two = s[j:j + 2]
            if two in AROMATIC_IN_BRACKETS:
                element, j = two.capitalize(), j + 2
            elif ch in AROMATIC_IN_BRACKETS:
                element, j = ch.upper(), j + 1
            else:
                self.fail(j, f"invalid aromatic symbol '{ch}'")
            aromatic = True
        else:
            self.fail(j, "missing element symbol")

        chirality = None
        if at(j) == "@":
            if at(j + 1) == "@":
                chirality, j = "@@", j + 2
            elif s[j + 1:j + 3] in RARE_CHIRALITY:
                n, end = number(j + 3)
                if n is None:
                    self.fail(j + 3, f"missing number after '@{s[j + 1:j + 3]}'")
                written = s[j:end]
                chirality = {"@TH1": "@", "@TH2": "@@"}.get(written, written)
                j = end
            else:
                chirality, j = "@", j + 1

        hydrogens = 0
        if at(j) == "H":
            n, j = number(j + 1)
            hydrogens = 1 if n is None else n

        charge = 0
        if at(j) in "+-":
            sign = at(j)
            n, end = number(j + 1)
            if n is not None:
                charge, j = n, end
            else:
                k = j
                while k < len(s) and s[k] == sign:
                    k += 1
                charge, j = k - j, k
            if sign == "-":
                charge = -charge

        atom_class = None
        if at(j) == ":":
            atom_class, end = number(j + 1)
            if atom_class is None:
                self.fail(j, "missing atom class number after ':'")
            j = end

        if at(j) != "]":
            self.fail(j, f"unexpected '{s[j]}' in brackets")
        self.i = j + 1
        return Atom(element, aromatic=aromatic, bracket=True, isotope=isotope, charge=charge,
                    hydrogens=hydrogens, chirality=chirality, atom_class=atom_class)

    def bond_symbol(self):
        i = self.i
        if self.last == "bond":
            self.fail(i, "two bond symbols in a row")
        if self.prev is None:
            self.fail(i, f"missing atom before bond '{self.s[i]}'")
        self.bond = (self.s[i], i)
        self.before_bond = self.last
        self.i += 1
        self.mark("bond", i)

    def open_branch(self):
        i = self.i
        if self.prev is None:
            self.fail(i, "missing atom before '('")
        if self.last == "bond":
            self.fail(self.bond[1], "a bond symbol must come after '(', not before it")
        if self.last == "open":
            self.fail(i, "a branch cannot start with '('")
        self.branches.append((self.prev, i))
        self.i += 1
        self.mark("open", i)

    def close_branch(self):
        i = self.i
        if not self.branches:
            self.fail(i, "unmatched ')'")
        if self.last == "open":
            self.fail(i, "empty branch")
        if self.last == "bond":
            self.fail(self.bond[1], f"missing atom after bond '{self.bond[0]}'")
        self.prev = self.branches.pop()[0]
        self.i += 1
        self.mark("close", i)

    def dot(self):
        i = self.i
        if self.last == "bond":
            self.fail(self.bond[1], f"missing atom after bond '{self.bond[0]}'")
        if self.prev is None:
            self.fail(i, "missing atom before '.'")
        if self.branches:
            self.fail(i, "'.' is not allowed inside a branch")
        self.prev = None
        self.i += 1
        self.mark("dot", i)

    def ring_bond(self):
        s, i = self.s, self.i
        after = self.before_bond if self.last == "bond" else self.last
        if self.prev is None or after == "open":
            self.fail(i, "missing atom before ring bond number")
        if s[i] == "%":
            if s[i + 1:i + 2] == "(":
                end = s.find(")", i + 2)
                digits = s[i + 2:end] if end > 0 else ""
                if not digits.isdigit():
                    self.fail(i, "'%(' must be followed by a number and ')'")
                number, self.i = int(digits), end + 1
            else:
                digits = s[i + 1:i + 3]
                if not (len(digits) == 2 and digits.isdigit()):
                    self.fail(i, "'%' must be followed by two digits")
                number, self.i = int(digits), i + 3
        else:
            number, self.i = int(s[i]), i + 1

        symbol, symbol_at = self.bond if self.bond else ("", i)
        atom = self.prev
        if number in self.rings:
            first, first_symbol, first_at, slot = self.rings.pop(number)
            if first == atom:
                self.fail(i, f"ring bond {number} joins an atom to itself")
            begin, end, written = self.ring_pair(first, atom, first_symbol, symbol, symbol_at)
            self.add_bond(begin, end, written, i, ring_closure=True)
            self.atoms[first].neighbor_order[slot] = atom
            self.atoms[atom].neighbor_order.append(first)
        else:
            self.rings[number] = (atom, symbol, symbol_at, len(self.atoms[atom].neighbor_order))
            self.atoms[atom].neighbor_order.append(_OPEN_RING)  # filled in when the ring closes
        self.bond = None
        self.mark("ring", i)

    def ring_pair(self, first: int, second: int, first_symbol: str, second_symbol: str, second_at: int):
        """Return the bond closed by a ring bond number, as (begin, end, symbol).

        A bond symbol can be written at either end; the bond begins at the atom it is written after.
        """
        if first_symbol and second_symbol:
            if first_symbol in DIRECTIONAL and second_symbol in DIRECTIONAL:
                # "/" written at the second atom points from it to the first one.
                if FLIP[second_symbol] != first_symbol:
                    self.fail(second_at, "conflicting directions at the two ends of the ring bond")
            elif first_symbol in DIRECTIONAL or second_symbol in DIRECTIONAL:
                plain = second_symbol if first_symbol in DIRECTIONAL else first_symbol
                if plain != "-":
                    self.fail(second_at, "conflicting bond symbols at the two ends of the ring bond")
                if second_symbol in DIRECTIONAL:
                    return second, first, second_symbol
            elif first_symbol != second_symbol:
                self.fail(second_at, "conflicting bond symbols at the two ends of the ring bond")
            return first, second, first_symbol
        if second_symbol:
            return second, first, second_symbol
        return first, second, first_symbol

    def add_bond(self, begin: int, end: int, symbol: str, position: int, ring_closure: bool):
        pair = frozenset((begin, end))
        if pair in self.pairs:
            self.fail(position, "these two atoms are already bonded")
        self.pairs.add(pair)
        if symbol:
            kind = BOND_KIND[symbol]
        elif self.atoms[begin].aromatic and self.atoms[end].aromatic:
            kind = "aromatic"
        else:
            kind = "single"
        self.bonds.append(Bond(begin, end, kind, symbol, ring_closure, position))

    def fragments(self) -> List[List[int]]:
        parent = list(range(len(self.atoms)))

        def root(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a

        for b in self.bonds:
            parent[root(b.begin)] = root(b.end)
        groups: Dict[int, List[int]] = {}
        for a in range(len(self.atoms)):
            groups.setdefault(root(a), []).append(a)
        return sorted(groups.values(), key=lambda g: g[0])
