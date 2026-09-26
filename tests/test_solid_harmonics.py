"""Tests for the derived real solid harmonics and their angular-momentum matrices.

These are white-box tests of the derivation: the polynomials are checked to be
harmonic (Laplacian zero) and mutually orthogonal with exact rational
arithmetic, and the matrices are checked against the invariants that pin them
down ([L_i, L_j] = i e_ijk L_k, L^2 = l(l+1), Hermiticity, the m spectrum).
The one thing these tests cannot check is whether ORCA's AO labelling matches
this basis -- that is measured on the fixtures instead (the Eu export's f-block
overlap is exactly R (x) I_7, and the atomic-term check on the f6 fixture is
the end-to-end pin).
"""

from __future__ import annotations

import numpy as np
import pytest

from fblockkit.analysis.solid_harmonics import (
    _add,
    _derivative,
    _inner,
    _operator,
    _scale,
    angular_momentum_matrices,
    component_labels,
    normalisation_squares,
    real_solid_harmonics,
)

SHELLS = (0, 1, 2, 3)


def _laplacian(polynomial):
    second = _add(_derivative(_derivative(polynomial, 0), 0),
                  _add(_derivative(_derivative(polynomial, 1), 1),
                       _derivative(_derivative(polynomial, 2), 2)))
    return second


@pytest.mark.parametrize("l", SHELLS)
def test_the_polynomials_are_harmonic(l):
    # the Laplacian of every basis function must vanish identically; this is the
    # check that caught the r^2 terms written out wrongly during development
    for polynomial in real_solid_harmonics(l):
        assert not _laplacian(polynomial), f"l={l}: {polynomial} is not harmonic"


@pytest.mark.parametrize("l", SHELLS)
def test_the_polynomials_are_mutually_orthogonal(l):
    basis = real_solid_harmonics(l)
    for i, left in enumerate(basis):
        for j, right in enumerate(basis):
            value = _inner(left, right)
            if i == j:
                assert value > 0
            else:
                assert value == 0


@pytest.mark.parametrize("l", SHELLS)
def test_commutation_and_casimir(l):
    lx, ly, lz = angular_momentum_matrices(l)
    size = 2 * l + 1
    assert np.abs(lx @ ly - ly @ lx - 1j * lz).max() < 1e-12
    assert np.abs(ly @ lz - lz @ ly - 1j * lx).max() < 1e-12
    assert np.abs(lz @ lx - lx @ lz - 1j * ly).max() < 1e-12
    assert np.abs(lx @ lx + ly @ ly + lz @ lz - l * (l + 1) * np.eye(size)).max() < 1e-12


@pytest.mark.parametrize("l", SHELLS)
def test_the_matrices_are_hermitian(l):
    for matrix in angular_momentum_matrices(l):
        assert np.abs(matrix - matrix.conj().T).max() < 1e-12


@pytest.mark.parametrize("l", SHELLS)
def test_the_m_spectrum_is_minus_l_to_l(l):
    # the eigenvalues of i L_z in the real basis are the m values; basis-convention free
    _, _, lz = angular_momentum_matrices(l)
    values = np.linalg.eigvalsh(1j * lz)
    assert np.allclose(sorted(values), range(-l, l + 1), atol=1e-12)


def test_the_real_basis_couples_only_plus_minus_m_pairs():
    # measured structure: L_z is purely imaginary, and its only non-zero entries
    # join the (m, -m) partners with magnitude |m|
    _, _, lz = angular_momentum_matrices(3)
    assert np.abs(np.real(lz)).max() < 1e-12
    magnitudes = np.abs(np.imag(lz))
    assert magnitudes[1, 2] == pytest.approx(1.0)  # (+1, -1)
    assert magnitudes[3, 4] == pytest.approx(2.0)  # (+2, -2)
    assert magnitudes[5, 6] == pytest.approx(3.0)  # (+3, -3)
    assert magnitudes[0].max() < 1e-12  # m = 0 does not mix


def test_the_operator_closes_on_every_shell():
    # applying D_x to each basis function must land inside the same shell: the
    # residual of the rational projection is exactly zero
    for l in SHELLS:
        basis = real_solid_harmonics(l)
        norms = normalisation_squares(l)
        for polynomial in basis:
            image = _operator(polynomial, "x")
            residual = image
            for i, other in enumerate(basis):
                residual = _add(residual, _scale(other, -_inner(other, image) / norms[i]))
            assert not residual


def test_component_labels_are_the_measured_orbca_spelling():
    assert component_labels(1) == ("z", "x", "y")
    assert component_labels(2) == ("z2", "xz", "yz", "x2y2", "xy")
    assert component_labels(3) == ("0", "+1", "-1", "+2", "-2", "+3", "-3")


def test_unrecorded_shells_raise_with_a_next_step():
    with pytest.raises(ValueError, match="Next step"):
        real_solid_harmonics(4)
    with pytest.raises(ValueError, match="Next step"):
        component_labels(5)
