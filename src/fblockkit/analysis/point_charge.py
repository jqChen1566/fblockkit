"""S2: point-charge crystal-field (CEF) estimate from a structure.

What this module produces
-------------------------
The **exact one-electron matrix elements** of the electrostatic potential of a set
of point charges, over the seven 4f orbitals (``l = 3``, ``m = -3..3``), computed
from a structure alone.  Nothing is fitted and no effective operator is
introduced: the input is charges plus a centre, the output is a lattice sum, a
potential matrix and its eigenvalues.

Model
-----
A point charge ``q_i`` (in units of the elementary charge ``e``) sits at position
``R_i`` relative to the metal centre; one electron (charge ``-e``) at ``r`` feels

    V(r) = -(e^2 / 4 pi eps0) * sum_i q_i / |r - R_i| .                       (1)

For ``r < R_i`` -- the standard assumption of the point-charge model, namely that
the 4f density lies inside the ligand sphere -- every term of (1) has the
multipole expansion

    1 / |r - R_i| = sum_{k>=0} (r^k / R_i^{k+1}) * (4 pi / (2k+1))
                    * sum_{q=-k}^{k} Y_kq*(R_i^) * Y_kq(r^) ,                  (2)

with ``R_i^`` the direction of the ion and ``r^`` the direction of the electron.
Sandwiching (1) between 4f orbitals separates the radial from the angular part.
Only ``k = 0, 2, 4, 6`` survive: the angular integral (a Gaunt coefficient) obeys
the triangle rule ``|3 - 3| <= k <= 3 + 3`` and vanishes for odd ``k``, so for an
f shell the expansion **terminates** and (2) is not a truncation for ``r < R_i``.
The radial part is the expectation value ``<r^k>`` of the 4f radial density.

With the lattice sum defined as (the ``Z_kq`` are the real tesseral harmonics of
the conventions block below; ``q_i`` in units of ``e``)

    L_kq = (4 pi / (2k+1)) * sum_i q_i * Z_kq(R_i^) / R_i^{k+1}
         [units e Angstrom^-(k+1)]                                             (3)

the result is

    V_mm' = -(e^2 / 4 pi eps0) * sum_k <r^k> * sum_q L_kq * G^(kq)_mm' ,
    G^(kq)_mm' = integral dOmega  Y*_3m(Omega) Z_kq(Omega) Y_3m'(Omega) .      (4)

Equation (4) is what :func:`estimate` evaluates exactly (3j symbols, no
quadrature); its correctness for a concrete radial model is checked against a
brute-force 3D numerical integration in ``tests/test_point_charge.py``.

Conventions (must travel with any number quoted from here)
----------------------------------------------------------
* **Orbital basis**: the complex spherical harmonic ``Y_3m`` with the
  Condon-Shortley phase, ordered ``m = -3, ..., +3``; matrix index ``m + 3``.
  All 7x7 matrices here are in that order.
* **Ligand expansion**: the real tesseral harmonics ``Z_kq``,
  ``Z_k0 = Y_k0``, ``Z_kq = sqrt(2)*(-1)^q*Re Y_kq`` for ``q > 0`` (cosine type)
  and ``Z_k,-p = sqrt(2)*(-1)^p*Im Y_kp`` for ``p > 0`` (sine type) -- the
  cosine/sine ``q > 0`` / ``q < 0`` pairing, the same pairing that module A4 uses
  for its operator equivalents.  The two bases are related by a unitary
  transformation, so the complete ``k`` sum in (2) is basis independent: this
  choice makes the lattice sums ``L_kq`` (3) real while the orbital basis stays
  complex, which is exactly the classical crystal-field presentation.
* **Units**: distances in Angstrom, charges in units of ``e``; ``L_kq`` in
  ``e Angstrom^-(k+1)``; ``<r^k>`` in ``Angstrom^k``.  The prefactor is computed
  from the fundamental constants in this module,

      e^2/(4 pi eps0) = 1e8 * e^2 / (4 pi eps0) / (h c)   [cm^-1 Angstrom]
                      = 1.16141e5 cm^-1 Angstrom

  (``e``, ``eps0``: charge and vacuum permittivity; ``h c``: one reciprocal
  centimetre expressed as an energy times a length; ``1e8`` converts ``m^-1`` to
  ``cm^-1 Angstrom^-1``).  The value is about 1.16e5 cm^-1 Angstrom, i.e. two
  unit charges at 1 Angstrom attract with 1.16e5 cm^-1.
* **The radial moments are not shipped**: ``<r^k>`` (k = 2, 4, 6) are specific to
  the ion, the configuration and the relativistic method.  Without them this
  module returns geometry only -- lattice sums in ``e Angstrom^-(k+1)`` and the
  unit-radial-moment potential matrix -- and refuses to print cm^-1 parameters
  (:func:`report` says so explicitly).

Boundaries (what this is *not*)
------------------------------
* The output is a set of **l-shell** parameters ``A_k^q <r^k>`` acting inside the
  4f orbital space.  It is **not** the set of Stevens parameters ``B_k^q`` of the
  ``|J M>`` ground manifold that module A4 fits: turning one into the other needs
  an operator-equivalent projection, which is not performed here.  Numbers from
  the two modules are not comparable parameter by parameter, only through the
  spectra they generate.
* A point-charge estimate is **qualitative**.  Measured deviations between the
  point-charge model and ab initio crystal fields are large enough that the model
  is unusable for quantitative extraction of crystal-field parameters
  (``ungur2017abinitio``); this module is placed at the front of the workflow as
  a sanity check on the geometry and the charge assignment, not as a source of
  parameters (``scheie2021pycrystalfield`` documents the published algorithm
  route that a full treatment follows).

Provenance of the implemented algorithm
---------------------------------------
The published route (build the lattice sums from a structure, evaluate the
one-electron CEF Hamiltonian, diagonalise it) is that of the PyCrystalField
paper; that software is GPL, so only the **published algorithm idea** is reused
here -- no code is linked, copied or transcribed -- and the implementation below
is written from the multipole expansion (2)-(4) with an exact Wigner-3j routine.
"""

from __future__ import annotations

import functools
import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Any, Mapping, Sequence

import numpy as np

from ..knowledge.models import (
    EVIDENCE_LITERATURE,
    EVIDENCE_MEASURED,
    Evidence,
    ReportSection,
)
from .geometry import Atom

# --- fundamental constants (SI 2019 / CODATA 2018 values) --------------------

ELEMENTARY_CHARGE = 1.602176634e-19      # C (exact)
VACUUM_PERMITTIVITY = 8.8541878128e-12   # F/m
PLANCK_CONSTANT = 6.62607015e-34         # J s (exact)
SPEED_OF_LIGHT = 2.99792458e8            # m/s (exact)

#: ``e^2/(4 pi eps0)`` in ``cm^-1 Angstrom``.  Derivation: ``e^2/(4 pi eps0)`` is
#: an energy times a length (J m); dividing by ``h c`` (J m) gives the reciprocal
#: length of the same electrostatic energy, and ``1e8`` converts ``m^-1`` into
#: ``cm^-1 Angstrom^-1`` -- because ``1/r[m] = 1e10 / r[Angstrom]`` and
#: ``1/hc[J cm] = 100/hc[J m]``, hence the factor ``1e10/100 = 1e8``.
#: Value: 1.161410...e5 cm^-1 Angstrom (two unit charges 1 Angstrom apart).
COULOMB_CM1_ANGSTROM = (
    1e8
    * (ELEMENTARY_CHARGE**2 / (4.0 * math.pi * VACUUM_PERMITTIVITY))
    / (PLANCK_CONSTANT * SPEED_OF_LIGHT)
)

#: Ranks that survive for a 4f-4f matrix element: even ``k`` (the 3j ``(3 k 3;
#: 0 0 0)`` vanishes for odd ``k``) with ``k <= 6`` (the triangle rule with
#: ``l = l' = 3``).
F_RANKS = (2, 4, 6)

#: Distance (Angstrom) below which an atom is taken to *be* the centre: it is
#: then excluded, because the point-charge model has no self-interaction term
#: (the metal's own charge does not act on its own electron).
CENTER_TOLERANCE = 1e-8

#: Relative cut-off for :attr:`PCEstimate.nonzero_terms`, measured against the
#: largest ``k >= 2`` lattice sum (the ``k = 0`` monopole term is always non-zero
#: and carries no symmetry information).
NONZERO_RELATIVE_TOLERANCE = 1e-10

_ALLOWED_RADIAL_KEYS = F_RANKS

_EVIDENCE_ALGORITHM = Evidence(
    kind=EVIDENCE_LITERATURE,
    text=(
        "The published route implemented here -- build the crystal-field lattice sums from "
        "a structure, assemble the one-electron CEF Hamiltonian and diagonalise it -- is the "
        "PyCrystalField algorithm route. Only the published algorithm idea is reused: that "
        "software is GPL, so no code is linked, copied or transcribed, and the implementation "
        "in this module is written from the Laplace multipole expansion with an exact "
        "Wigner-3j routine."
    ),
    ref=(
        "Scheie A., J. Appl. Crystallogr., 2021, 54(1), 356-362, "
        "DOI 10.1107/S160057672001554X"
    ),
    bibkey="scheie2021pycrystalfield",
    url="https://doi.org/10.1107/S160057672001554X",
)

_EVIDENCE_CAVEAT = Evidence(
    kind=EVIDENCE_LITERATURE,
    text=(
        "Measured deviations between the point-charge model and ab initio crystal fields are "
        "large enough that the point-charge model is not acceptable for quantitative "
        "crystal-field work: the ab initio route (and its projection onto the |J M> manifold) "
        "is required for parameters, and a point-charge estimate stays a qualitative check of "
        "the geometry and the charge assignment."
    ),
    ref=(
        "Ungur L., Chibotaru L. F., Chem. Eur. J., 2017, 23(15), 3708-3718, "
        "DOI 10.1002/chem.201605102"
    ),
    bibkey="ungur2017abinitio",
    url="https://doi.org/10.1002/chem.201605102",
)

_EVIDENCE_MEASURED = Evidence(
    kind=EVIDENCE_MEASURED,
    text=(
        "Measured while validating this module (see tests/test_point_charge.py): the analytic "
        "multipole result (3)-(4) was checked against a brute-force 3D numerical integration "
        "of equation (1) over a Slater-type 4f radial density (zeta = 6 /Angstrom, cut off at "
        "4 Angstrom, inside the nearest ion at 5.385 Angstrom, so the r < R assumption holds "
        "on the whole grid) and the unit sphere.  The independent quadrature reproduces the "
        "7x7 matrix to 3.6e-14 relative on a 100x48x96 grid (3.6e-14 to 4.9e-14 across grids "
        "from 60x32x64 to 200x96x192), against a matrix whose largest entry is 3.0e4 cm^-1; a "
        "single off-axis charge reproduces it to 3.2e-14.  The Gaunt coefficients were checked "
        "against the same quadrature to 5.8e-16, and the octahedral cage of six equal charges "
        "gives zero k = 2 lattice sums to 4.3e-16 relative to the k = 4 terms while "
        "reproducing the octahedral allowed (k, q) set exactly ((4,0), (4,4), (6,0), (6,4)) "
        "and the 1 + 3 + 3 cubic splitting of the 4f orbitals."
    ),
    ref="reproduced by tests/test_point_charge.py (2026-09-26)",
)


class PointChargeError(ValueError):
    """Invalid point-charge input (structure, charges, centre, radial moments)."""


# --- angular machinery -------------------------------------------------------


def _phase(exponent: int) -> int:
    """``(-1)^exponent`` for integer exponents of either sign."""
    return -1 if exponent % 2 else 1


@functools.lru_cache(maxsize=None)
def wigner_3j(j1: int, j2: int, j3: int, m1: int, m2: int, m3: int) -> float:
    """Wigner 3j symbol ``(j1 j2 j3; m1 m2 m3)`` for integer angular momenta.

    Racah's formula, evaluated exactly: the alternating sum is accumulated as a
    :class:`fractions.Fraction` (exact rational arithmetic on integer
    factorials) and only the final square root is floating point.  Selection
    rules (``m1 + m2 + m3 = 0``, ``|m_i| <= j_i``, triangle inequality) return
    0.0.  Only the ranks used by this module are needed, but the routine is
    general for non-negative integer ``j1, j2, j3``.

    Two of its properties are used as checks in the tests: the orthogonality sum
    ``sum_{m1 m2} (j1 j2 j3; m1 m2 m3)^2 = 1/(2 j3 + 1)`` and the known value
    ``(j j 0; 0 0 0) = (-1)^j / sqrt(2 j + 1)``.
    """
    for limit in (j1, j2, j3):
        if limit < 0:
            raise PointChargeError(
                f"Wigner 3j: negative angular momentum j = {limit}. Next step: pass "
                f"non-negative integer angular momenta (this module uses j <= 6)."
            )
    for value, limit in ((m1, j1), (m2, j2), (m3, j3)):
        if value < -limit or value > limit:
            return 0.0
    if m1 + m2 + m3 != 0:
        return 0.0
    if j3 > j1 + j2 or j3 < abs(j1 - j2):
        return 0.0

    factorial = math.factorial
    delta = Fraction(
        factorial(j1 + j2 - j3) * factorial(j1 - j2 + j3) * factorial(-j1 + j2 + j3),
        factorial(j1 + j2 + j3 + 1),
    )
    weight = (
        factorial(j1 + m1)
        * factorial(j1 - m1)
        * factorial(j2 + m2)
        * factorial(j2 - m2)
        * factorial(j3 + m3)
        * factorial(j3 - m3)
    )
    lower = max(0, j2 - j3 - m1, j1 - j3 + m2)
    upper = min(j1 + j2 - j3, j1 - m1, j2 + m2)
    total = Fraction(0)
    for z in range(lower, upper + 1):
        denominator = (
            factorial(z)
            * factorial(j1 + j2 - j3 - z)
            * factorial(j1 - m1 - z)
            * factorial(j2 + m2 - z)
            * factorial(j3 - j2 + m1 + z)
            * factorial(j3 - j1 - m2 + z)
        )
        total += Fraction(_phase(z), denominator)
    symbol = _phase(j1 - j2 - m3) * total
    return float(symbol * math.sqrt(float(delta * weight)))


def associated_legendre(l: int, m: int, x: Any) -> Any:
    """``P_l^m(x)`` including the Condon-Shortley phase, ``0 <= m <= l``.

    ``x = cos(theta)``; accepts scalars and numpy arrays.  Built from the
    standard recurrences ``P_m^m = (-1)^m (2m-1)!! (1-x^2)^{m/2}`` and
    ``(l-m) P_l^m = (2l-1) x P_{l-1}^m - (l+m-1) P_{l-2}^m`` (both including the
    Condon-Shortley phase), so the harmonic convention is fixed by this routine.
    """
    if m < 0 or m > l:
        raise PointChargeError(
            f"associated Legendre P_{l}^{m} is not defined for these indices "
            f"(need 0 <= m <= l). Next step: use spherical_harmonic(), which handles "
            f"negative m through Y_l,-m = (-1)^m conj(Y_lm)."
        )
    x = np.asarray(x, dtype=float)
    p_mm = np.ones_like(x)
    if m > 0:
        somx2 = np.sqrt(np.maximum(0.0, 1.0 - x * x))
        factor = 1.0
        for _ in range(1, m + 1):
            p_mm = -p_mm * factor * somx2
            factor += 2.0
    if l == m:
        return p_mm
    p_m1 = x * (2 * m + 1) * p_mm
    if l == m + 1:
        return p_m1
    for ll in range(m + 2, l + 1):
        p_ll = ((2 * ll - 1) * x * p_m1 - (ll + m - 1) * p_mm) / (ll - m)
        p_mm, p_m1 = p_m1, p_ll
    return p_m1


def spherical_harmonic(l: int, m: int, theta: Any, phi: Any) -> Any:
    """Complex spherical harmonic ``Y_lm(theta, phi)`` with the Condon-Shortley phase.

    ``Y_lm = sqrt((2l+1)/(4 pi) * (l-|m|)!/(l+|m|)!) P_l^|m|(cos theta) e^{i m phi}``
    for ``m >= 0``, extended to negative ``m`` by ``Y_l,-m = (-1)^m conj(Y_lm)``.
    Accepts scalars and numpy arrays (elementwise).

    This is the convention of the module's orbital basis (``l = 3``), and the
    quadrature validation in the tests uses it as the *definition* of the basis
    while checking everything else (expansion, Gaunt coefficients, prefactor)
    independently.
    """
    if abs(m) > l or l < 0:
        raise PointChargeError(
            f"Y_{l}{m} does not exist (need l >= 0 and |m| <= l). Next step: pass a valid "
            f"(l, m) pair, e.g. (3, -2) for a 4f orbital."
        )
    m_abs = abs(m)
    norm = math.sqrt(
        (2 * l + 1) / (4.0 * math.pi) * math.factorial(l - m_abs) / math.factorial(l + m_abs)
    )
    value = norm * associated_legendre(l, m_abs, np.cos(np.asarray(theta, dtype=float)))
    value = value * np.exp(1j * m * np.asarray(phi, dtype=float))
    if m < 0:
        value = _phase(m_abs) * value
    return value


def real_tesseral(k: int, q: int, theta: Any, phi: Any) -> Any:
    """Real tesseral harmonic ``Z_kq`` (cosine for ``q > 0``, sine for ``q < 0``).

    ``Z_k0 = Y_k0``; ``Z_kq = sqrt(2) (-1)^q Re Y_kq`` (``q > 0``);
    ``Z_k,-p = sqrt(2) (-1)^p Im Y_kp`` (``p > 0``).  For ``l = 1`` this gives
    the familiar ``sqrt(3/4 pi) * (x, y, z)/r``.  The set ``{Z_kq}`` is a
    unitary transformation of ``{Y_kq}``, so
    ``sum_q Z_kq(R^) Z_kq(r^) = sum_q Y_kq*(R^) Y_kq(r^)`` -- the completeness
    identity that lets the ligand expansion be real while the orbital basis
    stays complex.
    """
    if abs(q) > k or k < 0:
        raise PointChargeError(
            f"Z_{k}{q} does not exist (need k >= 0 and |q| <= k). Next step: pass a valid "
            f"(k, q) pair with even k, e.g. (4, 3)."
        )
    if q == 0:
        # Y_k0 is already real (P_k0 is real and e^{i 0 phi} = 1); take the real
        # part so that every Z_kq comes back with a real dtype.
        return np.real(spherical_harmonic(k, 0, theta, phi))
    p = abs(q)
    if q > 0:
        return math.sqrt(2.0) * _phase(q) * np.real(spherical_harmonic(k, q, theta, phi))
    return math.sqrt(2.0) * _phase(p) * np.imag(spherical_harmonic(k, p, theta, phi))


def _tesseral_terms(k: int, q: int) -> tuple[tuple[complex, int], ...]:
    """``Z_kq`` as a linear combination ``sum c_sigma Y_k,sigma``."""
    if q == 0:
        return ((1.0 + 0.0j, 0),)
    p = abs(q)
    root = math.sqrt(2.0) / 2.0
    if q > 0:
        return ((root * _phase(q) + 0.0j, p), (root + 0.0j, -p))
    return ((-1j * root * _phase(p), p), (1j * root, -p))


def gaunt_coefficient(l1: int, m1: int, l2: int, m2: int, l3: int, m3: int) -> float:
    """``integral dOmega Y*_l1m1 Y_l2m2 Y_l3m3`` (complex harmonics, CS phase).

    Gaunt's formula through two Wigner 3j symbols::

        integral = (-1)^m1 * sqrt((2l1+1)(2l2+1)(2l3+1)/(4 pi))
                   * (l1 l2 l3; 0 0 0) * (l1 l2 l3; -m1 m2 m3) .

    The result is real; it is zero unless ``m2 = m1 - m3``, ``k = l2`` is even and
    the three angular momenta form a triangle.
    """
    return (
        _phase(m1)
        * math.sqrt(
            (2 * l1 + 1) * (2 * l2 + 1) * (2 * l3 + 1) / (4.0 * math.pi)
        )
        * wigner_3j(l1, l2, l3, 0, 0, 0)
        * wigner_3j(l1, l2, l3, -m1, m2, m3)
    )


@functools.lru_cache(maxsize=None)
def _angular_matrix(k: int, q: int) -> np.ndarray:
    """``G^(kq)_mm' = integral Y*_3m Z_kq Y_3m' dOmega`` (7x7, index ``m + 3``)."""
    matrix = np.zeros((7, 7), dtype=complex)
    terms = _tesseral_terms(k, q)
    for m_index, m in enumerate(range(-3, 4)):
        for mp_index, mp in enumerate(range(-3, 4)):
            total = 0.0 + 0.0j
            for coefficient, sigma in terms:
                total += coefficient * gaunt_coefficient(3, m, k, sigma, 3, mp)
            matrix[m_index, mp_index] = total
    return matrix


def angular_matrix(k: int, q: int) -> np.ndarray:
    """``G^(kq)_mm' = integral dOmega Y*_3m(Omega) Z_kq(Omega) Y_3m'(Omega)``.

    A fresh 7x7 complex array in the module's basis order (``m = -3..3``,
    index ``m + 3``) for every call; the value is cached internally.
    """
    if k % 2 or k > 6 or k < 0 or abs(q) > k:
        raise PointChargeError(
            f"(k, q) = ({k}, {q}) is not a rank of a 4f-4f matrix element. "
            f"Next step: use even k in {F_RANKS} with |q| <= k -- odd k and k > 6 vanish by "
            f"the 3j triangle rule for l = 3."
        )
    return _angular_matrix(k, q).copy()


# --- inputs ------------------------------------------------------------------


@dataclass(frozen=True)
class _Ion:
    """One point charge of the estimate (element, distance, direction, charge)."""

    element: str
    distance: float
    unit: tuple[float, float, float]
    charge: float


def _element_key(symbol: Any) -> str:
    """Canonical element spelling: stripped, first letter upper case (as parse_xyz)."""
    text = str(symbol).strip()
    if not text:
        raise PointChargeError(
            "an atom carries an empty element symbol. Next step: fix the structure file "
            "(every XYZ coordinate line starts with an element symbol)."
        )
    return text.capitalize()


def _check_center(center: Any) -> tuple[float, float, float]:
    try:
        values = tuple(float(value) for value in center)
    except (TypeError, ValueError):
        raise PointChargeError(
            f"the centre {center!r} is not a triple of numbers. Next step: pass the metal "
            f"site as (x, y, z) in the same units as the atoms (Angstrom)."
        ) from None
    if len(values) != 3:
        raise PointChargeError(
            f"the centre has {len(values)} components, expected 3. Next step: pass the metal "
            f"site as (x, y, z) in the same units as the atoms (Angstrom)."
        )
    if not all(math.isfinite(value) for value in values):
        raise PointChargeError(
            f"the centre {values} contains a non-finite component. Next step: pass a finite "
            f"metal-site position."
        )
    return values


def _check_charges(charges: Any) -> dict[str, float]:
    if charges is None:
        raise PointChargeError(
            "no charges were given. Next step: pass a mapping element -> charge in units of "
            "e, e.g. {'O': -2.0, 'H': 0.4}; every element of the structure must appear in it."
        )
    try:
        items = list(charges.items())
    except AttributeError:
        raise PointChargeError(
            f"the charges {charges!r} are not a mapping. Next step: pass a mapping element -> "
            f"charge in units of e, e.g. {{'O': -2.0, 'H': 0.4}}."
        ) from None
    if not items:
        raise PointChargeError(
            "the charge mapping is empty. Next step: give a charge for every element of the "
            "structure, e.g. {'O': -2.0, 'H': 0.4}."
        )
    resolved: dict[str, float] = {}
    for symbol, value in items:
        key = _element_key(symbol)
        try:
            charge = float(value)
        except (TypeError, ValueError):
            raise PointChargeError(
                f"the charge {value!r} given for {key!r} is not a number. Next step: pass a "
                f"number in units of e, e.g. {{'O': -2.0}}."
            ) from None
        if not math.isfinite(charge):
            raise PointChargeError(
                f"the charge {charge!r} given for {key!r} is not finite. Next step: pass a "
                f"finite charge in units of e."
            )
        resolved[key] = charge
    return resolved


def _collect_ions(
    atoms: Sequence[Atom], center: tuple[float, float, float], charges: Mapping[str, float]
) -> tuple[tuple[_Ion, ...], int]:
    """Ions relative to the centre, plus the number of atoms sitting on the centre."""
    try:
        atom_list = list(atoms)
    except TypeError:
        raise PointChargeError(
            f"the structure {atoms!r} is not a sequence of atoms. Next step: parse the XYZ "
            f"file with fblockkit.analysis.geometry.parse_xyz and pass the atoms here."
        ) from None
    ions: list[_Ion] = []
    used: set[str] = set()
    on_center = 0
    for position, atom in enumerate(atom_list):
        try:
            element = _element_key(atom.element)
            vector = (
                float(atom.x) - center[0],
                float(atom.y) - center[1],
                float(atom.z) - center[2],
            )
        except AttributeError:
            raise PointChargeError(
                f"atom {position} has no element/x/y/z attributes. Next step: pass the atoms "
                f"returned by fblockkit.analysis.geometry.parse_xyz."
            ) from None
        if not all(math.isfinite(component) for component in vector):
            raise PointChargeError(
                f"atom {position} ({element}) has a non-finite coordinate. Next step: fix the "
                f"structure file."
            )
        distance = math.sqrt(sum(component * component for component in vector))
        if distance <= CENTER_TOLERANCE:
            # The atom *is* the centre (the caller normally passes the whole structure,
            # including the metal): the point-charge model has no self-interaction term.
            on_center += 1
            continue
        used.add(element)
        ions.append(
            _Ion(
                element=element,
                distance=distance,
                unit=tuple(component / distance for component in vector),
                charge=0.0,  # filled below, once the missing elements are known
            )
        )
    missing = sorted(element for element in used if element not in charges)
    if missing:
        raise PointChargeError(
            f"no charge given for {', '.join(missing)} (the structure contains "
            f"{', '.join(sorted(used))}). Next step: add every element of the structure to the "
            f"charges mapping, e.g. {{{', '.join(f'{element!r}: -2.0' for element in missing)}}}."
        )
    if not ions:
        raise PointChargeError(
            f"no point charges found: all {len(atom_list)} atom(s) sit on the centre. "
            f"Next step: pass a structure that contains the ligands around the metal site, or "
            f"move the centre onto the metal."
        )
    return (
        tuple(
            _Ion(
                element=ion.element,
                distance=ion.distance,
                unit=ion.unit,
                charge=charges[ion.element],
            )
            for ion in ions
        ),
        on_center,
    )


def _check_radial(radial: Any, ranks: Sequence[int]) -> dict[int, float]:
    """Validate ``radial`` and return it as ``{k: <r^k>}`` including ``<r^0> = 1``."""
    try:
        items = list(radial.items())
    except AttributeError:
        raise PointChargeError(
            f"radial {radial!r} is not a mapping k -> <r^k>. Next step: pass e.g. "
            f"{{2: 0.9, 4: 2.0, 6: 6.0}} in Angstrom^k, or None for a geometry-only estimate."
        ) from None
    required = {k for k in ranks if k != 0}
    resolved: dict[int, float] = {0: 1.0}
    for key, value in items:
        try:
            k = int(key)
        except (TypeError, ValueError):
            raise PointChargeError(
                f"the radial-moment key {key!r} is not an integer rank. Next step: use k in "
                f"{tuple(sorted(required))} (Angstrom^k)."
            ) from None
        if k not in required:
            raise PointChargeError(
                f"radial moment k = {k} is not used by this estimate (the ranks are "
                f"{tuple(sorted(required))}). Next step: pass exactly the moments "
                f"{{2: <r^2>, 4: <r^4>, 6: <r^6>}} in Angstrom^k -- odd k and k > 6 do not "
                f"contribute to a 4f-4f matrix element, and <r^0> = 1 by definition."
            )
        try:
            moment = float(value)
        except (TypeError, ValueError):
            raise PointChargeError(
                f"the radial moment <r^{k}> = {value!r} is not a number. Next step: pass a "
                f"positive number in Angstrom^{k}."
            ) from None
        if not math.isfinite(moment) or moment <= 0.0:
            raise PointChargeError(
                f"the radial moment <r^{k}> = {moment!r} is not positive and finite. "
                f"Next step: pass the expectation value of r^{k} (Angstrom^{k}) of the 4f "
                f"radial density, which is strictly positive."
            )
        resolved[k] = moment
    absent = sorted(required - set(resolved))
    if absent:
        raise PointChargeError(
            f"radial moments {absent} were not given (all of {tuple(sorted(required))} are "
            f"needed for ranks up to k = {ranks[-1]}). Next step: pass the complete mapping "
            f"{{2: <r^2>, 4: <r^4>, 6: <r^6>}} in Angstrom^k, or None to get the geometry-only "
            f"lattice sums."
        )
    return resolved


def _check_kmax(kmax: Any) -> tuple[int, ...]:
    try:
        value = int(kmax)
    except (TypeError, ValueError):
        raise PointChargeError(
            f"kmax = {kmax!r} is not an integer. Next step: pass kmax in (0, 2, 4, 6) -- 6 is "
            f"the largest rank that survives for a 4f-4f matrix element."
        ) from None
    if value % 2 or value < 0 or value > 6:
        raise PointChargeError(
            f"kmax = {kmax!r} is not an allowed rank. Next step: pass kmax in (0, 2, 4, 6); "
            f"odd k and k > 6 give a zero 3j symbol for l = 3, so they cannot contribute."
        )
    return tuple(range(0, value + 1, 2))


# --- result ------------------------------------------------------------------


@dataclass(frozen=True)
class PCEstimate:
    """Point-charge estimate of the 4f crystal-field potential.

    ``lattice_sums`` holds ``L_kq`` (equation (3)) in ``e Angstrom^-(k+1)`` for
    every ``k`` in :attr:`ranks` and every ``q`` in ``-k..k``, keyed by ``(k, q)``
    in ascending order.  ``parameters_cm1`` holds the crystal-field parameters
    ``-(e^2/4 pi eps0) <r^k> L_kq`` in cm^-1, or ``None`` when no radial moments
    were supplied (the refusal is repeated in :attr:`notes`).

    ``potential_matrix`` is the 7x7 matrix of ``V`` over the complex 4f basis
    (``m = -3..3``, index ``m + 3``; Hermitian).  Its scale is set by
    ``potential_units``: cm^-1 when the radial moments were given, otherwise the
    unit-radial-moment matrix in the units of ``lattice_sums``.
    ``potential_eigenvalues`` are its eigenvalues, ascending.

    ``nonzero_terms`` is the ``(k, q)`` pattern with ``k >= 2`` above the relative
    cut-off -- the k = 0 monopole term is excluded from it (it is non-zero for
    any charged environment and carries no symmetry information), while
    ``lattice_sums`` and ``parameters_cm1`` do list it.  ``radial`` echoes the
    moments that were used, including ``<r^0> = 1`` by definition.
    """

    center: tuple[float, float, float]
    lattice_sums: Mapping[tuple[int, int], float]
    parameters_cm1: Mapping[tuple[int, int], float] | None
    potential_matrix: np.ndarray
    potential_eigenvalues: np.ndarray
    potential_units: str
    notes: tuple[str, ...]
    nonzero_terms: tuple[tuple[int, int], ...]
    ranks: tuple[int, ...]
    n_ions: int
    nearest_distance: float
    total_charge: float
    radial: Mapping[int, float] | None

    @property
    def splitting(self) -> float:
        """``max - min`` of the potential eigenvalues (in :attr:`potential_units`)."""
        values = np.asarray(self.potential_eigenvalues, dtype=float)
        return float(values.max() - values.min()) if values.size else 0.0

    @property
    def monopole_shift(self) -> float:
        """The ``k = 0`` (spherical) part of the matrix: ``-L_00 / sqrt(4 pi)`` scaled.

        With radial moments this is the common shift of all seven orbitals in
        cm^-1 (``-(e^2/4 pi eps0) * sum_i q_i / R_i``); without them it is the
        same quantity in the units of ``lattice_sums``.
        """
        scale = COULOMB_CM1_ANGSTROM if self.parameters_cm1 is not None else 1.0
        return float(-scale * self.lattice_sums[(0, 0)] / math.sqrt(4.0 * math.pi))


def estimate(
    atoms: Sequence[Atom],
    charges: Mapping[str, float],
    center: Sequence[float],
    radial: Mapping[int, float] | None = None,
    kmax: int = 6,
) -> PCEstimate:
    """Point-charge crystal-field estimate from a structure (spec S2).

    Parameters
    ----------
    atoms:
        The structure, as returned by
        :func:`fblockkit.analysis.geometry.parse_xyz`.  Atoms sitting on
        ``center`` (the metal itself) are skipped with a note: the model has no
        self-interaction term.
    charges:
        Mapping element symbol -> charge in units of ``e`` (e.g.
        ``{"O": -2.0, "H": 0.4}``).  Every element present in the structure must
        appear; a missing element is an error listing it.
    center:
        ``(x, y, z)`` of the metal site, in the same units as the atoms
        (Angstrom).
    radial:
        Optional mapping ``k -> <r^k>`` in ``Angstrom^k`` for ``k = 2, 4, 6``.
        When given, cm^-1 parameters and a cm^-1 potential matrix are produced.
        When absent, the module refuses to produce cm^-1 numbers (the radial
        expectation values are ion- and method-specific and are not shipped with
        this tool) and returns the geometry-only lattice sums with a note.
    kmax:
        Largest rank kept: 0, 2, 4 or 6 (default 6; the 3j triangle rule makes
        larger and odd ranks identically zero for an f shell).

    Returns
    -------
    PCEstimate
        Lattice sums, parameters (or ``None``), the 7x7 potential matrix, its
        eigenvalues and the notes carrying the assumptions and the conventions.
    """
    ranks = _check_kmax(kmax)
    center_values = _check_center(center)
    charge_map = _check_charges(charges)
    ions, on_center = _collect_ions(atoms, center_values, charge_map)
    moments = _check_radial(radial, ranks) if radial is not None else {k: 1.0 for k in ranks}
    with_radial = radial is not None
    scale = COULOMB_CM1_ANGSTROM if with_radial else 1.0

    # (3): L_kq = (4 pi/(2k+1)) sum_i q_i Z_kq(R_i^)/R_i^(k+1)
    unit = np.array([ion.unit for ion in ions], dtype=float)
    charge_array = np.array([ion.charge for ion in ions], dtype=float)
    distance = np.array([ion.distance for ion in ions], dtype=float)
    theta = np.arccos(np.clip(unit[:, 2], -1.0, 1.0))
    phi = np.arctan2(unit[:, 1], unit[:, 0])
    lattice_sums: dict[tuple[int, int], float] = {}
    for k in ranks:
        weight = charge_array / distance ** (k + 1)
        factor = 4.0 * math.pi / (2 * k + 1)
        for q in range(-k, k + 1):
            values = np.asarray(real_tesseral(k, q, theta, phi), dtype=float)
            lattice_sums[(k, q)] = float(factor * np.sum(weight * values))

    # (4): V_mm' = -(e^2/4 pi eps0) sum_k <r^k> sum_q L_kq G^(kq)_mm'
    matrix = np.zeros((7, 7), dtype=complex)
    for k in ranks:
        for q in range(-k, k + 1):
            value = lattice_sums[(k, q)]
            if value == 0.0:
                continue
            matrix += value * moments[k] * angular_matrix(k, q)
    matrix = -scale * matrix
    hermitian_defect = float(np.abs(matrix - matrix.conj().T).max()) if matrix.size else 0.0
    matrix_scale = float(np.abs(matrix).max()) if matrix.size else 0.0
    matrix = 0.5 * (matrix + matrix.conj().T)
    eigenvalues = np.linalg.eigvalsh(matrix)

    parameters_cm1: dict[tuple[int, int], float] | None = None
    if with_radial:
        parameters_cm1 = {
            (k, q): float(-COULOMB_CM1_ANGSTROM * moments[k] * lattice_sums[(k, q)])
            for k in ranks
            for q in range(-k, k + 1)
        }

    symmetry_scale = max(
        (abs(value) for (k, q), value in lattice_sums.items() if k >= 2),
        default=0.0,
    )
    if symmetry_scale > 0.0:
        nonzero_terms = tuple(
            (k, q)
            for k in ranks
            if k >= 2
            for q in range(-k, k + 1)
            if abs(lattice_sums[(k, q)]) > NONZERO_RELATIVE_TOLERANCE * symmetry_scale
        )
    else:
        nonzero_terms = tuple((k, q) for k in ranks if k >= 2 for q in range(-k, k + 1))

    nearest = min(ion.distance for ion in ions)
    notes: list[str] = []
    notes.append(
        "Assumption: the multipole expansion (Laplace expansion of 1/|r - R_i|) is used for "
        "r < R_i, i.e. the whole 4f density is taken to lie inside the ligand sphere. The "
        "nearest ion of this estimate is at "
        f"{nearest:.3f} Angstrom"
        + (
            f" and sqrt(<r^2>) = {math.sqrt(moments[2]):.3f} Angstrom, a margin of "
            f"{nearest / math.sqrt(moments[2]):.2f} x sqrt(<r^2>)"
            if with_radial and 2 in moments
            else " (no radial moments given, so the margin cannot be checked here)"
        )
        + "; the series terminates exactly at k = 6 for an f shell, so the r < R assumption "
        "is the only approximation in the angular part."
    )
    if with_radial and 2 in moments and nearest < 2.0 * math.sqrt(moments[2]):
        notes.append(
            "Caution: the nearest ion is closer than 2 sqrt(<r^2>), so a non-negligible part "
            "of the 4f density lies outside the ligand sphere and the point-charge expansion "
            "is not trustworthy for this geometry (inference, provisional: the 2 x sqrt(<r^2>) "
            "margin is a working rule, not a published bound)."
        )
    notes.append(
        "Convention: the 4f basis is the complex spherical harmonic Y_3m with the "
        "Condon-Shortley phase, ordered m = -3..3 (matrix index m + 3); the ligand multipole "
        "expansion uses the real tesseral harmonics Z_kq (cosine for q > 0, sine for q < 0). "
        "The lattice sums L_kq are in e Angstrom^-(k+1) (one factor of e is the ligand charge, "
        "the other is the electron charge carried by the prefactor e^2/(4 pi eps0) = "
        f"{COULOMB_CM1_ANGSTROM:.6g} cm^-1 Angstrom)."
    )
    if not with_radial:
        notes.append(
            "No radial expectation values were given, so no cm^-1 crystal-field parameters "
            "are produced: <r^k> (k = 2, 4, 6) is specific to the ion, the configuration and "
            "the relativistic method, and is not shipped with this tool. The lattice sums and "
            "the potential matrix above are geometry only (matrix entries in the units of the "
            "lattice sums, i.e. unit radial moments). Next step: pass radial={2: <r^2>, "
            "4: <r^4>, 6: <r^6>} in Angstrom^k (e.g. from a radial density of your own "
            "calculation) to obtain cm^-1 parameters."
        )
    notes.append(
        "Boundary: these are l-shell parameters A_k^q <r^k> acting inside the 4f orbital "
        "space. They are NOT the Stevens parameters B_k^q of the |J M> ground manifold that "
        "the crystal-field fit (module A4) determines -- the operator-equivalent projection "
        "from the l shell onto the J manifold is not performed here. Compare the two modules "
        "through the spectra they generate, not parameter by parameter."
    )
    if on_center:
        notes.append(
            f"{on_center} atom(s) sit on the centre and were excluded (the point-charge model "
            f"has no self-interaction term: the metal's own charge does not act on its own "
            f"electron)."
        )
    if matrix_scale > 0.0 and hermitian_defect > 1e-10 * matrix_scale:
        notes.append(
            f"The assembled matrix is not Hermitian to rounding (defect "
            f"{hermitian_defect:.3e}, matrix scale {matrix_scale:.3e}); the symmetric part was "
            f"used for the eigenvalues. This should not happen -- check the angular tables."
        )
    notes.append(
        "Qualitative only: measured deviations between point-charge and ab initio crystal "
        "fields make a point-charge estimate unusable for quantitative parameter extraction "
        "(ungur2017abinitio). Use it as a geometry and charge-assignment sanity check."
    )

    return PCEstimate(
        center=center_values,
        lattice_sums=lattice_sums,
        parameters_cm1=parameters_cm1,
        potential_matrix=matrix,
        potential_eigenvalues=eigenvalues,
        potential_units="cm^-1" if with_radial else "e Angstrom^-(k+1)",
        notes=tuple(notes),
        nonzero_terms=nonzero_terms,
        ranks=ranks,
        n_ions=len(ions),
        nearest_distance=float(nearest),
        total_charge=float(sum(ion.charge for ion in ions)),
        radial=dict(moments) if with_radial else None,
    )


# --- report ------------------------------------------------------------------

_SYMMETRY_NOTE = (
    "Symmetry: compare the non-zero (k, q) pattern above with the allowed set of the point "
    "group you confirmed in S1/A4 -- a true three-fold axis keeps q in {0, +-3, +-6} and a "
    "cubic (Oh) site keeps (4, 0), (4, 4), (6, 0), (6, 4) with the whole second rank zero. "
    "More non-zero terms than the point group allows means the geometry or the charge "
    "assignment breaks the symmetry (check the z axis convention first)."
)

_CAVEAT = (
    "Caveat: a point-charge estimate is qualitative. Measured deviations between point-charge "
    "and ab initio crystal fields are large enough that the model is not usable for "
    "quantitative crystal-field parameters (Ungur & Chibotaru, Chem. Eur. J. 2017, 23, "
    "3708-3718, DOI 10.1002/chem.201605102, [ungur2017abinitio]); the quantitative route is "
    "the ab initio one, projected onto the |J M> manifold, whose published algorithm route is "
    "that of Scheie, J. Appl. Crystallogr. 2021, 54, 356-362, DOI 10.1107/S160057672001554X, "
    "[scheie2021pycrystalfield]. Read this section as a check of the geometry and of the "
    "charge assignment, not as a source of parameters."
)


def _format_charges(charges: Mapping[str, float] | None) -> str:
    if not charges:
        return "(none given)"
    try:
        items = sorted(((str(key), float(value)) for key, value in charges.items()))
    except (AttributeError, TypeError, ValueError):
        return str(charges)
    return ", ".join(f"{key} {value:+.3f} e" for key, value in items)


def report(
    estimate: PCEstimate,
    charges: Mapping[str, float],
    radial: Mapping[int, float] | None = None,
) -> ReportSection:
    """Build the S2 report section for one estimate.

    ``charges`` and ``radial`` are the inputs the estimate was built from; they
    are echoed so the section is self-contained.  Every number in the body is
    taken from ``estimate`` and never recomputed here, so a ``radial`` that
    disagrees with the estimate is reported as a note instead of silently
    changing the numbers.
    """
    lines: list[str] = []
    center = estimate.center
    lines.append(
        "V(r) = -(e^2/4 pi eps0) sum_i q_i / |r - R_i|, expanded in spherical harmonics for "
        "r < R_i (the whole 4f density is taken to lie inside the ligand sphere; the expansion "
        "terminates exactly at k = 6 for an f shell)."
    )
    lines.append(
        f"Centre (metal site): ({center[0]:.4f}, {center[1]:.4f}, {center[2]:.4f}) Angstrom"
    )
    lines.append(
        f"Point charges: {estimate.n_ions} ion(s); total charge "
        f"{estimate.total_charge:+.3f} e (ions only, the centre is excluded); nearest ion "
        f"{estimate.nearest_distance:.3f} Angstrom"
    )
    lines.append(f"  Charges used: {_format_charges(charges)}")
    used_moments = estimate.radial
    if used_moments:
        lines.append(
            "  Radial moments used (Angstrom^k): "
            + ", ".join(
                f"<r^{int(k)}> = {float(value):.4f}"
                for k, value in sorted(used_moments.items(), key=lambda item: int(item[0]))
                if int(k) != 0
            )
        )
    else:
        lines.append("  Radial moments used: none (geometry-only estimate)")
    # the numbers below always come from `estimate`; the arguments are only echoed,
    # so a mismatch between them is stated rather than smoothed over
    if radial is not None and not used_moments:
        lines.append(
            "  Note: a radial mapping was passed to report() but the estimate was built "
            "without one -- the geometry-only numbers of the estimate are shown."
        )
    elif radial is None and used_moments:
        lines.append(
            "  Note: report() was called with radial=None; the moments the estimate was built "
            "from are shown."
        )
    lines.append("")

    lines.append(
        "Lattice sums L_kq = (4 pi/(2k+1)) sum_i q_i Z_kq(R_i^)/R_i^(k+1), "
        "units e Angstrom^-(k+1):"
    )
    # Values that are numerically zero (below the relative tolerance of the largest
    # k >= 2 sum) print as 0: their exact digits are trig cancellation noise, and
    # hiding them keeps the report stable across platforms and BLAS builds.
    reference = max(
        (abs(value) for (k, q), value in estimate.lattice_sums.items() if k >= 2),
        default=0.0,
    )
    tolerance = NONZERO_RELATIVE_TOLERANCE * reference

    def _lattice_text(k: int, q: int) -> str:
        value = estimate.lattice_sums[(k, q)]
        if k >= 2 and abs(value) <= tolerance:
            return "+0.000000e+00"
        return f"{value:+.6e}"

    for k in estimate.ranks:
        row = " ".join(f"({k},{q}) {_lattice_text(k, q)}" for q in range(-k, k + 1))
        lines.append(f"  k = {k}: {row}")
    lines.append(
        "  Non-zero symmetry-relevant terms (k >= 2, above "
        f"{NONZERO_RELATIVE_TOLERANCE:.0e} of the largest k >= 2 lattice sum; entries at or "
        "below that line print as 0 in the table above): "
        + (", ".join(f"({k},{q})" for k, q in estimate.nonzero_terms) or "none")
    )
    lines.append("")

    if estimate.parameters_cm1 is not None:
        lines.append(
            "Crystal-field parameters -(e^2/4 pi eps0) <r^k> L_kq, cm^-1 "
            "(k = 0 is the spherical monopole term: a common shift of all seven orbitals, not "
            "a splitting parameter):"
        )
        for k in estimate.ranks:
            row = " ".join(
                f"({k},{q}) {estimate.parameters_cm1[(k, q)]:+.4f}" for q in range(-k, k + 1)
            )
            lines.append(f"  k = {k}: {row}")
    else:
        lines.append(
            "Crystal-field parameters in cm^-1: not produced -- no radial expectation values "
            "<r^k> were given, and they are ion- and method-specific (not shipped with this "
            "tool). The lattice sums above are geometry only."
        )
    lines.append("")

    lines.append(
        "4f orbital energies (eigenvalues of the potential matrix, "
        f"{estimate.potential_units}, ascending):"
    )
    lines.append(
        "  " + " ".join(f"{value:+.4f}" for value in estimate.potential_eigenvalues)
    )
    lines.append(
        f"  Overall splitting (max - min): {estimate.splitting:.4f} {estimate.potential_units}"
    )
    lines.append("")

    lines.append(_SYMMETRY_NOTE)
    lines.append("")
    lines.append(
        "Boundary: these are l-shell parameters acting inside the 4f orbital space; they are "
        "not the Stevens parameters B_k^q of the |J M> manifold fitted by module A4 (no "
        "operator-equivalent projection is performed here)."
    )
    lines.append("")
    lines.append(_CAVEAT)
    lines.append("")
    lines.append("Assumptions and limits carried by this estimate:")
    lines += [f"  - {note}" for note in estimate.notes]
    return ReportSection(title="S2 point-charge crystal-field estimate", body="\n".join(lines))


def evidence() -> tuple[Evidence, ...]:
    """Provenance of this module: two literature anchors and one measurement."""
    return (_EVIDENCE_ALGORITHM, _EVIDENCE_CAVEAT, _EVIDENCE_MEASURED)
