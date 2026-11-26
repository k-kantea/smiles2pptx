"""Element symbols, indexed by atomic number (0 is the wildcard atom '*')."""

SYMBOLS = (
    "*",
    "H", "He",
    "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt",
    "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
    "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
)

ATOMIC_NUMBER = {symbol: number for number, symbol in enumerate(SYMBOLS)}

# Normal valences, used to count the hydrogens of atoms written outside brackets and
# to decide which aromatic atoms take a double bond in the Kekule structure. The list
# follows OpenSMILES and RDKit (iodine and phosphorus also have their higher valences).
NORMAL_VALENCE = {
    "B": (3,), "C": (4,), "N": (3, 5), "O": (2,), "F": (1,),
    "Si": (4,), "P": (3, 5, 7), "S": (2, 4, 6), "Cl": (1,),
    "Ge": (4,), "As": (3, 5), "Se": (2, 4, 6), "Br": (1,),
    "Sb": (3, 5), "Te": (2, 4, 6), "I": (1, 3, 5),
}


def normal_valences(element: str, charge: int = 0) -> tuple:
    """Return the normal valences of an atom, or () if it has none.

    A charged atom takes those of the element with the same number of valence electrons
    ([N+] like C, [O+] like N, [C-] like N)."""
    if charge == 0 or element == "*":
        return NORMAL_VALENCE.get(element, ())
    number = ATOMIC_NUMBER[element] - charge
    if not 0 < number < len(SYMBOLS):
        return ()
    return NORMAL_VALENCE.get(SYMBOLS[number], ())
