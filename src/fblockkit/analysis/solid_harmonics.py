"""Real solid harmonics and their angular-momentum matrices, derived from polynomials.

Why derive instead of tabulate: the atomic-term check needs the matrices of
L_x, L_y, L_z in the real solid-harmonic basis that ORCA uses for its atomic
orbitals.  Transcribing a real-to-complex transformation table from memory is
exactly the kind of step this project has been burned by before; instead the
matrices are built here from the defining polynomials with exact rational
arithmetic, and every invariant that pins them down is asserted in the tests
(commutation ``[L_i, L_j] = i e_ijk L_k``, ``L^2 = l(l+1) 1``, the ladder
structure of ``L_x +- i L_y``, and the measured ORCA AO order).

Conventions, all pinned to measurements on the fixtures:

- The real solid harmonics are the polynomials

  l=0: 1
  l=1: z; x; y
  l=2: 3z^2-r^2; xz; yz; x^2-y^2; xy
  l=3: 5z^3-3z r^2; x(5z^2-r^2); y(5z^2-r^2); z(x^2-y^2); xyz; x(x^2-3y^2); y(3x^2-y^2)

  normalised to ``<f|f> = 1`` over the unit sphere (the convention that makes
  an AO overlap factorise as radial x angular; measured on the Eu export: the
  f-block overlap is exactly ``R (x) I_7`` with the cross-component blocks at
  2e-16).
- ORCA's order within a shell is the m-indexed one: ``s``; ``pz px py``;
  ``dz2 dxz dyz dx2y2 dxy``; ``f0 f+1 f-1 f+2 f-2 f+3 f-3`` -- i.e. m = 0,
  +1, -1, +2, -2, ... with the standard real combinations (measured on the
  N2 and Eu exports; for d and f alike).  ``component_labels`` spells them the
  way ORCA does.
- ``L = -i (r x grad)``: the matrices returned are Hermitian; in this real
  basis L_z is purely imaginary (antisymmetric) while L_x and L_y are real
  symmetric up to sign conventions that the invariants fix.
"""

from __future__ import annotations

import math
from fractions import Fraction

import numpy as np

__all__ = [
    "component_labels",
    "real_solid_harmonics",
    "angular_momentum_matrices",
]

#: A polynomial keyed by monomial exponents (a, b, c) for x^a y^b z^c.
Polynomial = dict[tuple[int, int, int], Fraction]

# The unnormalised real solid harmonics, in ORCA's m-indexed order per shell.
_RAW: dict[int, tuple[Polynomial, ...]] = {
    0: [{(0, 0, 0): Fraction(1)}],
    1: [
        {(0, 0, 1): Fraction(1)},  # z   (m = 0)
        {(1, 0, 0): Fraction(1)},  # x   (m = +1 position)
        {(0, 1, 0): Fraction(1)},  # y   (m = -1 position)
    ],
    2: [
        # 3z^2-r^2 written out with r^2 = x^2+y^2+z^2: 2z^2-x^2-y^2 (harmonic)
        {(0, 0, 2): Fraction(2), (2, 0, 0): Fraction(-1), (0, 2, 0): Fraction(-1)},
        {(1, 0, 1): Fraction(1)},  # xz
        {(0, 1, 1): Fraction(1)},  # yz
        {(2, 0, 0): Fraction(1), (0, 2, 0): Fraction(-1)},  # x^2-y^2
        {(1, 1, 0): Fraction(1)},  # xy
    ],
    3: [
        # z(5z^2-3r^2) = 2z^3-3x^2z-3y^2z; x(5z^2-r^2) = 4xz^2-x^3-xy^2; cyclic
        {(0, 0, 3): Fraction(2), (2, 0, 1): Fraction(-3), (0, 2, 1): Fraction(-3)},
        {(1, 0, 2): Fraction(4), (3, 0, 0): Fraction(-1), (1, 2, 0): Fraction(-1)},
        {(0, 1, 2): Fraction(4), (2, 1, 0): Fraction(-1), (0, 3, 0): Fraction(-1)},
        {(2, 0, 1): Fraction(1), (0, 2, 1): Fraction(-1)},  # z(x^2-y^2)
        {(1, 1, 1): Fraction(1)},  # xyz
        {(3, 0, 0): Fraction(1), (1, 2, 0): Fraction(-3)},  # x(x^2-3y^2)
        {(2, 1, 0): Fraction(3), (0, 3, 0): Fraction(-1)},  # y(3x^2-y^2)
    ],
}

_COMPONENTS: dict[int, tuple[str, ...]] = {
    0: ("",),
    1: ("z", "x", "y"),
    2: ("z2", "xz", "yz", "x2y2", "xy"),
    3: ("0", "+1", "-1", "+2", "-2", "+3", "-3"),
}


def component_labels(l: int) -> tuple[str, ...]:
    """ORCA's component spelling for a shell, in ORCA's AO order."""
    if l not in _COMPONENTS:
        raise ValueError(
            f"angular momentum l={l} has no recorded ORCA component spelling. "
            "Next step: record it from an export's OrbitalLabels (measured forms are "
            "s, p and f)."
        )
    return _COMPONENTS[l]


# --- exact polynomial arithmetic over the sphere -----------------------------


def _double_factorial(n: int) -> int:
    if n <= 0:
        return 1
    result = 1
    while n > 0:
        result *= n
        n -= 2
    return result


def _angular_integral(a: int, b: int, c: int) -> Fraction:
    """The exact integral of x^a y^b z^c over the unit sphere (0 for odd powers)."""
    if a % 2 or b % 2 or c % 2:
        return Fraction(0)
    numerator = _double_factorial(a - 1) * _double_factorial(b - 1) * _double_factorial(c - 1)
    denominator = _double_factorial(a + b + c + 1)
    return Fraction(4) * Fraction(numerator, denominator)  # 4*pi / (4*pi) cancels


def _inner(left: Polynomial, right: Polynomial) -> Fraction:
    total = Fraction(0)
    for (a1, b1, c1), c1_coeff in left.items():
        for (a2, b2, c2), c2_coeff in right.items():
            total += c1_coeff * c2_coeff * _angular_integral(a1 + a2, b1 + b2, c1 + c2)
    return total


def _derivative(poly: Polynomial, axis: int) -> Polynomial:
    out: Polynomial = {}
    for (a, b, c), coeff in poly.items():
        powers = [a, b, c]
        if powers[axis] == 0:
            continue
        factor = powers[axis]
        powers[axis] -= 1
        key = (powers[0], powers[1], powers[2])
        out[key] = out.get(key, Fraction(0)) + coeff * factor
    return out


def _scale(poly: Polynomial, factor: Fraction) -> Polynomial:
    return {key: coeff * factor for key, coeff in poly.items()}


def _add(left: Polynomial, right: Polynomial) -> Polynomial:
    out = dict(left)
    for key, coeff in right.items():
        out[key] = out.get(key, Fraction(0)) + coeff
    return {key: coeff for key, coeff in out.items() if coeff}


def _multiply_coordinate(poly: Polynomial, axis: int) -> Polynomial:
    out: Polynomial = {}
    for (a, b, c), coeff in poly.items():
        powers = [a, b, c]
        powers[axis] += 1
        key = (powers[0], powers[1], powers[2])
        out[key] = out.get(key, Fraction(0)) + coeff
    return out


def _operator(poly: Polynomial, kind: str) -> Polynomial:
    """Apply the differential part of the angular-momentum operator (no -i factor).

    ``kind`` is ``x``, ``y`` or ``z``; the returned polynomial is the result of
    the real differential operator D_x = y d_z - z d_y (and cyclic), so that
    L_x = -i D_x.
    """
    # explicit and boring beats clever: write the three cases out
    if kind == "x":
        first = _multiply_coordinate(_derivative(poly, 2), 1)  # y * d/dz
        second = _scale(_multiply_coordinate(_derivative(poly, 1), 2), Fraction(-1))  # -z * d/dy
        return _add(first, second)
    if kind == "y":
        first = _multiply_coordinate(_derivative(poly, 0), 2)  # z * d/dx
        second = _scale(_multiply_coordinate(_derivative(poly, 2), 0), Fraction(-1))  # -x * d/dz
        return _add(first, second)
    if kind == "z":
        first = _multiply_coordinate(_derivative(poly, 1), 0)  # x * d/dy
        second = _scale(_multiply_coordinate(_derivative(poly, 0), 1), Fraction(-1))  # -y * d/dx
        return _add(first, second)
    raise ValueError(f"unknown operator axis {kind!r}")  # pragma: no cover


def real_solid_harmonics(l: int) -> tuple[Polynomial, ...]:
    """The 2l+1 *unnormalised* real solid harmonics, in ORCA's order (exact rationals).

    The squared norms are returned by :func:`normalisation_squares`; keeping the
    polynomials unnormalised keeps every step below in exact rational arithmetic
    (the normalisation constants themselves are square roots and therefore
    irrational).
    """
    if l not in _RAW:
        raise ValueError(
            f"l={l} is not recorded. Next step: add its real solid harmonics to _RAW "
            "(they are standard tables; the invariants in the tests will pin them)."
        )
    return tuple(_RAW[l])


def normalisation_squares(l: int) -> tuple[Fraction, ...]:
    """``<p_k|p_k>`` over the unit sphere for each unnormalised harmonic."""
    return tuple(_inner(poly, poly) for poly in real_solid_harmonics(l))


def _gram_is_diagonal(basis: tuple[Polynomial, ...], l: int) -> None:
    """The real solid harmonics of one shell are orthogonal: assert it, do not assume it."""
    for i, left in enumerate(basis):
        for j, right in enumerate(basis):
            if i != j and _inner(left, right) != 0:
                raise ArithmeticError(
                    f"the recorded l={l} harmonics are not orthogonal (<p{i}|p{j}> != 0); "
                    "the table is not a solid-harmonic basis."
                )


def angular_momentum_matrices(l: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The matrices of L_x, L_y, L_z (as -i D) in the *orthonormal* basis of ``l``.

    The projection of ``D_k p_j`` onto the shell is carried out in exact rational
    arithmetic and the residual is verified: a non-zero residual raises instead of
    returning a silently truncated matrix.  The single square root per element
    (the normalisation) is the only step leaving the rationals.
    """
    basis = real_solid_harmonics(l)
    norms = normalisation_squares(l)
    _gram_is_diagonal(basis, l)
    size = len(basis)
    matrices = []
    for kind in ("x", "y", "z"):
        matrix = np.zeros((size, size), dtype=np.complex128)
        for j, poly in enumerate(basis):
            image = _operator(poly, kind)
            residual = image
            for i, other in enumerate(basis):
                # orthonormal-basis matrix element: <p_i|D p_j> / sqrt(N_i N_j)
                raw = _inner(other, image)
                residual = _add(residual, _scale(other, -Fraction(raw, norms[i])))
                matrix[i, j] = -1j * float(raw) / math.sqrt(float(norms[i] * norms[j]))
            if any(residual.values()):
                raise ArithmeticError(
                    f"D_{kind} does not close on the l={l} shell: the expansion residual "
                    f"is non-zero ({residual}). The basis is not a complete shell."
                )
        matrices.append(matrix)
    return tuple(matrices)  # type: ignore[return-value]
