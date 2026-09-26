"""Element table and Z ranges (shared data of the knowledge layer).

Why it lives here: both the recipe layer (basis-set/ECP matching) and the
analysis layer (structure analysis identifying f-block elements) need it, and
those two layers do not depend on each other -- shared data can only sit in the
knowledge layer.
"""

from __future__ import annotations

_SYMBOLS = (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn "
    "Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd "
    "Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac "
    "Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr"
).split()

ELEMENT_Z: dict[str, int] = {symbol: z for z, symbol in enumerate(_SYMBOLS, start=1)}
LANTHANIDES = range(57, 72)
ACTINIDES = range(89, 104)


class ElementError(ValueError):
    """Unknown element symbol."""


def element_z(symbol: str) -> int:
    """Element symbol -> atomic number; an unknown symbol is an error (no guessing)."""
    key = symbol.strip().capitalize()
    if key not in ELEMENT_Z:
        raise ElementError(
            f"unknown element symbol: {symbol!r}. Next step: check the spelling of the "
            f"element."
        )
    return ELEMENT_Z[key]


def is_f_element(symbol: str) -> bool:
    """Whether this is a lanthanide or actinide element."""
    z = element_z(symbol)
    return z in LANTHANIDES or z in ACTINIDES


def parse_range(text: str) -> range:
    """``"La-Lu"`` / ``"Ce"`` -> a Z range."""
    if "-" not in text:
        z = element_z(text)
        return range(z, z + 1)
    low, _, high = text.partition("-")
    return range(element_z(low), element_z(high) + 1)
