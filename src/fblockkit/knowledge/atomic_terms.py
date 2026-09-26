"""Hund's rules for an l^n shell (shared data of the knowledge layer).

Why it lives here: both the analysis-layer atomic-term check and the
recipe-layer spin-orbit template need the ground term of an l^n configuration,
and the analysis layer is not allowed to be imported by the recipe layer, so the
table and the rule sit with the other shared data.

The values are the textbook ones: S = min(n, 2(2l+1) - n)/2, the tabulated L per
angular momentum (mirrored about the half-filled shell), and J = |L - S| below
half filling, L + S above it.
"""

from __future__ import annotations

from dataclasses import dataclass


class TermError(ValueError):
    """The configuration has no Hund term (with a next step)."""


#: Hund's-rule L for a shell of l = 2 / 3 with n electrons (n = 1..4l+1; mirrored
#: for more than half filling).  Standard textbook values, listed per l because the
#: check is used for d and f shells alike.
HUND_L = {
    2: {1: 2, 2: 3, 3: 3, 4: 2, 5: 0},
    3: {1: 3, 2: 5, 3: 6, 4: 6, 5: 5, 6: 3, 7: 0},
}


@dataclass(frozen=True)
class HundTerm:
    """The Hund's-rule ground term of an l^n configuration."""

    angular: int
    n_electrons: int
    s: float
    l: float
    j: float
    j_maximum: float  # L + S, the value a scalar stretched component carries
    half_filled_or_less: bool


def hund_term(angular: int, n_electrons: int) -> HundTerm:
    """Hund's rules for ``l^n``: S = min(n, 4l+2-n)/2, then the tabulated L and J."""
    capacity = 2 * (2 * angular + 1)
    if not 1 <= n_electrons <= capacity - 1:
        raise TermError(
            f"{n_electrons} electrons do not form a Hund term of an l={angular} shell "
            f"(capacity {capacity}). Next step: check the shell occupation."
        )
    table = HUND_L.get(angular)
    if table is None:
        raise TermError(
            f"no Hund-rule L table for l={angular}. Next step: add it (the values are "
            "textbook; the pattern mirrors about the half-filled shell)."
        )
    less = n_electrons <= 2 * angular + 1
    s = min(n_electrons, capacity - n_electrons) / 2.0
    l_value = table[min(n_electrons, capacity - n_electrons)]
    j = abs(l_value - s) if less else l_value + s
    return HundTerm(
        angular=angular,
        n_electrons=n_electrons,
        s=s,
        l=float(l_value),
        j=j,
        j_maximum=float(l_value) + s,
        half_filled_or_less=less,
    )


__all__ = ["HUND_L", "HundTerm", "TermError", "hund_term"]
