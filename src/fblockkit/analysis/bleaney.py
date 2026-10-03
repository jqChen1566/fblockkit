"""Bleaney's analytic anisotropy comparator (menu 35's optional comparator).

Bleaney's theory (1972) gives the anisotropic part of a lanthanide's molar
magnetic susceptibility in the high-temperature limit, as a product of the
ion's Bleaney constant ``C_J`` and the second-rank axial (rhombic) crystal-
field parameter ``B_0^2`` (``B_2^2``).  This module implements the **modern
SI form** of that relation together with the Bleaney-constant table of the
same source:

    chi_ax = -mu_0 mu_B^2 C_J B_0^2 / (10 (k T)^2)
    chi_rh = -mu_0 mu_B^2 C_J B_2^2 / (30 (k T)^2)

(Chilton, Suturina, Parker, Kuprov et al., Acc. Chem. Res. 2020, eq. (3);
their eq. (2) fixes the companion convention ``chi_ax = 3 chi_z / 2`` and
``chi_rh = (chi_x - chi_y) / 2`` on the traceless tensor, with eigenvalues
ordered ``|chi_x| < |chi_y| < |chi_z|``).  With ``B_0^2`` in joules the
output is the SI molar susceptibility in ``m^3``; the module converts from
wavenumbers (``cm^-1``, the literature's habit) with ``hc x 100``.

The measured anchors (tests pin all of them):

- the **sign structure** against the experimental Delta-chi table of a
  C3-lanthanide tag series (Tb +30.87, Tm -14.35, Yb -6.05, in 1e-32 m^3;
  the 2022 non-canonical-amino-acid tag study): a negative ``C_J`` (Tb) gives
  a positive ``chi_ax`` and vice versa -- all three match;
- the **magnitude**: inverting the experimental Tb tag value through the
  formula implies ``B_0^2 ~ 1.6e2 cm^-1``, inside the typical lanthanide-tag
  range (the tests bound it to 50-400 cm^-1);
- the **T^-2 law** (high-temperature limit) and the direct-value regression
  (Dy, ``B_0^2`` = 1000 cm^-1, 300 K -> 226.512e-32 m^3);
- the **two-table scale bridge**: the classic tabulation (Bleaney's own,
  scaled to ``C_Dy = -100``; as reproduced in the CEST-agent review's
  Table 1) relates to the modern constants by 1.81 +/- 0.03 per ion
  (individual ratios 1.78-1.84 -- inside the classic table's two-
  significant-figure rounding).  Both tables ship here; the comparator uses
  the modern one.  The closed form ``C_J = g_J^2 theta_2 J(J+1)
  (4 J(J+1) - 3)`` (``theta_2`` = ``<J||alpha||J>``, the 2020 review's
  Table 1) reproduces the six modern constants within their printed
  precision, so the 1.81 factor is the two normalisations' scale, not a
  fitted number; the module exposes it as ``c_j_closed_form``.

Boundaries (stated): the formula is the high-temperature limit -- it fails
once the crystal-field splitting is comparable to or larger than ``kT``
(the 2020 review's own warning; the principal axis can even rotate by up to
90 degrees along the series); a same-source published pair (B_0^2, Delta chi)
for one lanthanide system was not located in this round, so the absolute
scale rests on the source's unit-explicit formula plus the anchors above.

The classic framework's relation is now recorded verbatim (Canard et al.,
Inorg. Chem. 2009, eqs. (9)-(16), following Bleaney 1972 and
Golding-Halton 1972): chi_zz - Tr(chi)/3 = 2 N_A C_j B_0^2 and
chi_xx - chi_yy = 2 N_A C_j sqrt(6) B_2^2, with the **absolute**
C_j = -beta^2 (1 + p^j) xi^j / (60 (k T)^2) (beta the Bohr magneton, xi^j a
tabulated coefficient per 4f^n configuration, 1 + p^j the excited-multiplet
term), *tabulated as relative values scaled to C_Dy = -100*.  The numeric
bridge is **closed at the formula level** (2026-10-03): the high-temperature
second-order result for the ground multiplet is the closed form
``g_J^2 theta_2 J(J+1) (4 J(J+1) - 3)``; with the exact-rational ``theta_2``
of the 2020 review's Table 1 it matches the six published modern constants
at their printed resolution (every deviation at or below half of the
value's last printed digit; Tb sits exactly on its rounding boundary), and
it shares the classic table's sign structure for those six ions -- the only
classic-table exceptions, Sm and Eu, are the ions the 1972 source itself
flags for excited-multiplet-dominated behaviour.  The
classic table is therefore the same quantity at two significant figures
under the C_Dy = -100 normalisation, and no separate xi^j tabulation is
needed for the conversion.  One convention item is recorded, not resolved:
the Canard equation's ``2 B_0^2`` entry counts a factor 2 against this
derivation (the two transcriptions' operator normalisation); the reading
notes carry the derivation record.
"""

from __future__ import annotations

from fractions import Fraction

from ..knowledge.models import EVIDENCE_LITERATURE, EVIDENCE_MEASURED, Evidence

__all__ = [
    "BleaneyError",
    "BOHR_MAGNETON",
    "CM1_TO_JOULE",
    "MU_0",
    "K_B",
    "BLEANEY_CJ_MODERN",
    "BLEANEY_CJ_CLASSIC",
    "chi_ax",
    "chi_rh",
    "c_j_closed_form",
    "implied_b02_cm1",
    "comparator_lines",
    "assumptions",
    "evidence",
]

#: CODATA 2018 constants (same values as analysis/pnmr.py).
BOHR_MAGNETON = 9.2740100783e-24  # J/T
MU_0 = 1.25663706212e-6  # N/A^2
K_B = 1.380649e-23  # J/K
#: one wavenumber in joules: h c x 100
CM1_TO_JOULE = 6.62607015e-34 * 2.99792458e8 * 100.0


class BleaneyError(ValueError):
    """The comparator cannot evaluate the request (with a next step)."""


#: The modern Bleaney constants (Chilton et al., Acc. Chem. Res. 2020,
#: eq. (3) and its table): the six ions that source lists.
BLEANEY_CJ_MODERN = {
    "Tb": -158.0,
    "Dy": -181.0,
    "Ho": -71.2,
    "Er": 58.8,
    "Tm": 95.3,
    "Yb": 39.2,
}

#: The classic tabulation (relative values, scaled to C_Dy = -100 by the
#: framework's own convention; as reproduced in the CEST/paracest-agent
#: review's Table 1, 2010) -- the comparator's cross-table check uses it;
#: against the modern absolute constants the per-ion ratio is 1.81 +/- 0.03
#: (documented in the module docstring; a numeric cross-scale regression
#: awaits Bleaney 1972's xi^j tabulation).
BLEANEY_CJ_CLASSIC = {
    "Ce": -6.3,
    "Pr": -11.0,
    "Nd": -4.2,
    "Pm": 2.0,
    "Sm": -0.7,
    "Eu": 4.0,
    "Gd": 0.0,
    "Tb": -86.0,
    "Dy": -100.0,
    "Ho": -39.0,
    "Er": 33.0,
    "Tm": 53.0,
    "Yb": 22.0,
}

#: The closed form's inputs per ion: (J, g_J, theta_2), with theta_2 =
#: <J||alpha||J> as exact rationals from the 2020 review's Table 1
#: ("Equivalence Coefficients for the Low-Energy Terms of Late Ln(III)
#: Ions").  The same coefficients appear as decimals (a_2^J / alpha) in
#: the 1998 bicelle study's Table 1, which traces them to Bleaney 1972 and
#: Abragam & Bleaney 1970, Table 20.
_CJ_CLOSED_FORM = {
    "Tb": (6, Fraction(3, 2), Fraction(-1, 99)),
    "Dy": (Fraction(15, 2), Fraction(4, 3), Fraction(-2, 315)),
    "Ho": (8, Fraction(5, 4), Fraction(-1, 450)),
    "Er": (Fraction(15, 2), Fraction(6, 5), Fraction(4, 1575)),
    "Tm": (6, Fraction(7, 6), Fraction(1, 99)),
    "Yb": (Fraction(7, 2), Fraction(8, 7), Fraction(2, 63)),
}


def _ion_key(ion: str) -> str:
    """Normalise an ion name ("Dy3+", "dy" -> "Dy")."""
    key = ion.strip().rstrip("+-").replace("3", "").strip()
    return key[:1].upper() + key[1:].lower() if key else key


def c_j_closed_form(ion: str) -> float:
    """The Bleaney constant from the closed form ``g_J^2 theta_2 J(J+1)
    (4 J(J+1) - 3)`` -- the high-temperature second-order (van Vleck)
    result for the ground multiplet.

    Exact rational arithmetic throughout, so the well-conditioned ions
    come back bit-exact (Tb -157.5, Ho -71.25, Er 58.752); all six match
    the shipped modern constants within their printed precision (see the
    module docstring).  Ions outside the six raise with a next step.
    """
    key = _ion_key(ion)
    if key not in _CJ_CLOSED_FORM:
        raise BleaneyError(
            f"the ion {ion!r} has no closed-form inputs here (the table "
            f"covers {', '.join(_CJ_CLOSED_FORM)}). Next step: use one of "
            "those ions; for the full series see BLEANEY_CJ_CLASSIC (the "
            "C_Dy = -100 scale)."
        )
    j, g, theta2 = _CJ_CLOSED_FORM[key]
    l = j * (j + 1)
    return float(g * g * theta2 * l * (4 * l - 3))


def _resolve_ion(ion: str) -> float:
    """The modern C_J for an ion name (raises with a next step otherwise)."""
    key = _ion_key(ion)
    if key not in BLEANEY_CJ_MODERN:
        raise BleaneyError(
            f"the ion {ion!r} has no modern Bleaney constant here (the source "
            f"lists {', '.join(BLEANEY_CJ_MODERN)}). Next step: give one of "
            "those ions; the classic table covers the full series but on the "
            "C_Dy = -100 scale (see BLEANEY_CJ_CLASSIC)."
        )
    return BLEANEY_CJ_MODERN[key]


def chi_ax(ion: str, b02_cm1: float, temperature_k: float) -> float:
    """The Bleaney axial anisotropy chi_ax = -mu0 muB^2 C_J B0^2 / (10 (kT)^2).

    ``b02_cm1``: the second-rank axial crystal-field parameter in cm^-1 (the
    literature's convention); the return value is the SI molar susceptibility
    in m^3 (the standard traceless-tensor delta-chi_ax = chi_z - (chi_x +
    chi_y)/2).
    """
    if temperature_k <= 0:
        raise BleaneyError(
            f"the temperature {temperature_k} K is not positive. Next step: "
            "give the measurement temperature (the formula's high-T limit)."
        )
    c_j = _resolve_ion(ion)
    b02_joule = float(b02_cm1) * CM1_TO_JOULE
    return (
        -MU_0
        * BOHR_MAGNETON**2
        * c_j
        * b02_joule
        / (10.0 * (K_B * float(temperature_k)) ** 2)
    )


def chi_rh(ion: str, b22_cm1: float, temperature_k: float) -> float:
    """The Bleaney rhombic anisotropy chi_rh (factor 30; same conventions)."""
    if temperature_k <= 0:
        raise BleaneyError(
            f"the temperature {temperature_k} K is not positive. Next step: "
            "give the measurement temperature."
        )
    c_j = _resolve_ion(ion)
    b22_joule = float(b22_cm1) * CM1_TO_JOULE
    return (
        -MU_0
        * BOHR_MAGNETON**2
        * c_j
        * b22_joule
        / (30.0 * (K_B * float(temperature_k)) ** 2)
    )


def implied_b02_cm1(ion: str, chi_ax_m3: float, temperature_k: float) -> float:
    """Invert the axial relation: the B_0^2 (cm^-1) implied by a chi_ax value."""
    if temperature_k <= 0:
        raise BleaneyError(
            f"the temperature {temperature_k} K is not positive. Next step: "
            "give the measurement temperature (the formula's high-T limit)."
        )
    c_j = _resolve_ion(ion)
    if c_j == 0.0:
        raise BleaneyError(
            f"the ion {ion!r} carries C_J = 0 (the isotropic case); the axial "
            "relation cannot be inverted."
        )
    factor = MU_0 * BOHR_MAGNETON**2 * c_j / (10.0 * (K_B * float(temperature_k)) ** 2)
    return float(chi_ax_m3) / (-factor * CM1_TO_JOULE)


def comparator_lines(ion: str, b02_cm1: float, temperature_k: float) -> tuple[str, ...]:
    """The report block the menu appends when the comparator is requested."""
    value = chi_ax(ion, b02_cm1, temperature_k)
    key = _ion_key(ion)
    lines = [
        "",
        "Bleaney comparator (the high-temperature analytic limit; one-source table):",
        f"  ion {key} (C_J = {BLEANEY_CJ_MODERN[key]:.1f}, modern scale); "
        f"B_0^2 = {float(b02_cm1):g} cm^-1; T = {float(temperature_k):g} K",
        f"  chi_ax = {value / 1e-32:.4f}e-32 m^3 "
        "(SI molar, the standard traceless delta-chi_ax)",
        f"  closed form g_J^2 theta_2 J(J+1)(4J(J+1)-3) = "
        f"{c_j_closed_form(key):.2f} (the shipped constant: "
        f"{BLEANEY_CJ_MODERN[key]:.1f}; both to the source's printed precision)",
    ]
    lines += [f"  - {note}" for note in assumptions()]
    return tuple(lines)


def assumptions() -> tuple[str, ...]:
    """The comparator's assumptions and limits, for the report."""
    return (
        "high-temperature limit: the crystal-field splitting must stay well "
        "below kT (the source's own warning; beyond it the anisotropy and even "
        "the principal axis vary with temperature)",
        "the formula's B_0^2 is the energy-unit axial crystal-field parameter "
        "(entered in cm^-1 here); the modern SI equation and its C_J table come "
        "from one source (Acc. Chem. Res. 2020, eq. (3))",
        "the classic tabulation (scaled to C_Dy = -100) is the same quantity "
        "at two significant figures under its own normalisation; the closed "
        "form g_J^2 theta_2 J(J+1)(4J(J+1)-3) reproduces the modern constants "
        "at their printed precision and is the converter between the scales; "
        "the comparator does not mix the two tables",
        "no same-source published (B_0^2, delta-chi) pair for one lanthanide "
        "system was located in this round; the absolute scale rests on the "
        "unit-explicit source formula plus the sign/magnitude anchors",
    )


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the relation, the tables and the anchors."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "The modern Bleaney relation chi_ax = -mu_0 mu_B^2 C_J B_0^2 / "
                "(10 (kT)^2) and its C_J values (Tb -158, Dy -181, Ho -71.2, "
                "Er +58.8, Tm +95.3, Yb +39.2), with the companion PCS "
                "convention; the classic framework's working equations and the "
                "absolute C_j definition (-beta^2 (1+p) xi^j / (60 (kT)^2), "
                "tabulated relative to C_Dy = -100); the two-table ratio is "
                "1.81 +/- 0.03 per ion; the second-order operator-equivalent "
                "coefficients theta_2 = <J||alpha||J> as exact rationals "
                "(Table 1 of the 2020 review; the same numbers as a_2^J in "
                "the 1998 bicelle study's Table 1, traced there to Bleaney "
                "1972 and Abragam & Bleaney 1970, Table 20)."
            ),
            ref=(
                "Chilton/Suturina/Parker/Kuprov et al., Acc. Chem. Res. 2020, "
                "eqs. (2)-(3) and their table (DOI 10.1021/acs.accounts.0c00275); "
                "Canard/Ryan/Piguet et al., Inorg. Chem. 2009, 48, 2552-2562, "
                "eqs. (9)-(16) (the classic framework and the C_j definition); "
                "the classic table as reproduced in Acta Radiol. 2010 (CEST/"
                "paracest agents review); Prosser/Hwang/Vold, Biophys. J. "
                "1998, 74 (Table 1, p. 2408)"
            ),
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "Anchors: the sign structure matches the experimental C3-tag "
                "Delta-chi table (Tb +30.87, Tm -14.35, Yb -6.05, in 1e-32 m^3); "
                "inverting the Tb value implies B_0^2 ~ 1.6e2 cm^-1 (within the "
                "typical tag range); the T^-2 law is exact in the formula; the "
                "two-table ratio is 1.81 +/- 0.03 per ion; the closed form "
                "g_J^2 theta_2 J(J+1)(4J(J+1)-3) matches the six published "
                "modern constants at their printed resolution (deviations at "
                "or below half of each value's last printed digit, Tb exactly "
                "on its boundary; exact rational arithmetic throughout)."
            ),
            ref="tests/test_bleaney.py; the tag table's 2022 study (see the reading notes)",
        ),
    )
