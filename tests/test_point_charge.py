"""S2 point-charge tests: the multipole implementation against independent quadrature,
the geometry checks, the refusal paths and the report.

Evidence chain, from the weakest to the strongest:

1. conventions and closed forms pinned by hand: the Wigner 3j values and their
   orthogonality sum, the complex spherical harmonics against their closed forms,
   the real tesseral harmonics against their cartesian forms and against the
   addition theorem;
2. checks that follow from the physics: a single axial charge gives a matrix
   diagonal in m, the trace is the spherical average -7 * sum q_i / R_i, a
   rotation of the whole structure leaves the spectrum invariant, and the
   octahedral cage of six equal charges gives zero k = 2 lattice sums and
   reproduces the octahedral allowed (k, q) set exactly;
3. an **independent numerical quadrature**: the analytic multipole result is
   compared with a brute-force 3D integration of V(r) = -(e^2/4 pi eps0) *
   sum_i q_i / |r - R_i| over a Slater-type 4f radial density and the unit
   sphere.  The radial model is cut off inside the nearest ion, so the r < R
   assumption of the multipole expansion holds on the whole grid and the two
   routes must agree: **measured 3.6e-14 (100x48x96 grid) to 4.9e-14 (fine
   grid) relative**, asserted below at 1e-8 -- six orders below the field's own
   1e-6 target, with margin for platform differences in the quadrature.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from fblockkit.analysis import geometry, point_charge
from fblockkit.analysis.point_charge import (
    COULOMB_CM1_ANGSTROM,
    PointChargeError,
    PCEstimate,
    estimate,
)
from fblockkit.knowledge import sources

# --- the radial model used by the quadrature test ---------------------------
#
# R(r) proportional to r^3 exp(-zeta r) below r_cut and 0 above it.  zeta and
# r_cut are chosen so that (a) the neglected tail of the radial density is
# negligible -- the untruncated weight integral is the closed form 8!/12^9 and
# the truncated one agrees with it to 1.2e-12 relative, all of it the dropped
# tail (asserted below at 1e-9) -- and (b) the support ends well inside the
# nearest ion (4.0 against 5.385 Angstrom), so the multipole expansion is exact,
# not truncated, on the whole quadrature grid.
ZETA = 6.0
R_CUT = 4.0
RADIAL_EXACT_NORM = math.factorial(8) / 12**9  # integral r^8 exp(-12 r) dr on [0, inf)


def _radial(r):
    r = np.asarray(r, dtype=float)
    return np.where(r <= R_CUT, r**3 * np.exp(-ZETA * r), 0.0)


def _radial_grid(nodes: int) -> tuple[np.ndarray, np.ndarray]:
    x, w = np.polynomial.legendre.leggauss(nodes)
    return 0.5 * R_CUT * (x + 1.0), 0.5 * R_CUT * w


def _radial_norm(nodes: int = 400) -> float:
    r, w = _radial_grid(nodes)
    f = _radial(r)
    return float(np.sum(w * r**2 * f**2))


def _radial_moments(nodes: int = 400) -> dict[int, float]:
    r, w = _radial_grid(nodes)
    f = _radial(r)
    norm = float(np.sum(w * r**2 * f**2))
    return {k: float(np.sum(w * r ** (k + 2) * f**2) / norm) for k in (2, 4, 6)}


#: The two ions of the quadrature test: a divalent anion at 6.00 Angstrom and a
#: partly covalent H at 5.385 Angstrom, both inside-outside pair (r < R) for the
#: radial model above.
VALIDATION_IONS = (("O", (0.0, 0.0, 6.0), -2.0), ("H", (3.0, 4.0, 2.0), 0.4))


def _direct_matrix(ions, radial, norm, *, nr=100, nx=48, nphi=96) -> np.ndarray:
    """Equation (1) integrated directly: radial Gauss-Legendre x angular grid.

    Returns ``-e^2/(4 pi eps0) * sum_i q_i integral r^2 dr R(r)^2 dOmega
    Y*_3m (1/|r u - R_i|) Y_3m'`` in cm^-1, i.e. exactly the same object as
    ``PCEstimate.potential_matrix`` but computed without the multipole expansion,
    without Gaunt coefficients and without 3j symbols.  The f-orbital harmonics
    are the module's own ``spherical_harmonic`` (the basis convention is pinned
    separately, by the closed-form tests below); everything else is independent.
    """
    r, wr = _radial_grid(nr)
    xx, wx = np.polynomial.legendre.leggauss(nx)          # xx = cos(theta)
    phi = 2.0 * np.pi * np.arange(nphi) / nphi
    dphi = 2.0 * np.pi / nphi
    theta = np.arccos(xx)[None, :, None]
    weight = np.broadcast_to(wr[:, None, None] * wx[None, :, None] * dphi, (nr, nx, nphi))
    weight = weight * (r[:, None, None] ** 2 * _radial(r)[:, None, None] ** 2 / norm)
    orbitals = np.empty((nr, nx, nphi, 7), dtype=complex)
    for index, m in enumerate(range(-3, 4)):
        orbitals[..., index] = point_charge.spherical_harmonic(3, m, theta, phi[None, None, :])
    flat_orbitals = orbitals.reshape(-1, 7)
    flat_weight = weight.reshape(-1)
    radius = np.broadcast_to(r[:, None, None], (nr, nx, nphi)).reshape(-1)
    unit_x = np.broadcast_to(np.sin(theta) * np.cos(phi)[None, None, :], (nr, nx, nphi)).reshape(-1)
    unit_y = np.broadcast_to(np.sin(theta) * np.sin(phi)[None, None, :], (nr, nx, nphi)).reshape(-1)
    unit_z = np.broadcast_to(xx[None, :, None], (nr, nx, nphi)).reshape(-1)
    total = np.zeros((7, 7), dtype=complex)
    for _element, position, charge in ions:
        distance = np.sqrt(
            (radius * unit_x - position[0]) ** 2
            + (radius * unit_y - position[1]) ** 2
            + (radius * unit_z - position[2]) ** 2
        )
        g = flat_weight * charge / distance
        total += (flat_orbitals.conj() * g[:, None]).T @ flat_orbitals
    return -COULOMB_CM1_ANGSTROM * total


def _cage(distance: float, element: str = "O") -> list[geometry.Atom]:
    return [
        geometry.Atom(element, distance, 0.0, 0.0),
        geometry.Atom(element, -distance, 0.0, 0.0),
        geometry.Atom(element, 0.0, distance, 0.0),
        geometry.Atom(element, 0.0, -distance, 0.0),
        geometry.Atom(element, 0.0, 0.0, distance),
        geometry.Atom(element, 0.0, 0.0, -distance),
    ]


def _rotation(axis, angle: float) -> np.ndarray:
    """Proper rotation matrix (Rodrigues), so no coordinate stays axis aligned."""
    unit = np.asarray(axis, dtype=float)
    unit = unit / np.linalg.norm(unit)
    cross = np.array(
        [
            [0.0, -unit[2], unit[1]],
            [unit[2], 0.0, -unit[0]],
            [-unit[1], unit[0], 0.0],
        ]
    )
    return np.eye(3) + math.sin(angle) * cross + (1.0 - math.cos(angle)) * (cross @ cross)


# --- constants, 3j symbols and harmonics ------------------------------------


def test_coulomb_prefactor_value():
    """e^2/(4 pi eps0) converted from the fundamental constants: ~1.16e5 cm^-1 A."""
    assert COULOMB_CM1_ANGSTROM == pytest.approx(116140.97, abs=0.5)
    # the same number in the more familiar unit: 14.3996 eV Angstrom
    electron_volt = 1.602176634e-19
    assert COULOMB_CM1_ANGSTROM * 1.986445857e-23 / electron_volt == pytest.approx(14.3996, abs=1e-3)


def test_wigner_3j_known_values_and_orthogonality():
    """Racah formula: known closed forms plus the orthogonality sum, exactly."""
    assert point_charge.wigner_3j(1, 1, 0, 0, 0, 0) == pytest.approx(-1.0 / math.sqrt(3.0))
    assert point_charge.wigner_3j(2, 2, 0, 0, 0, 0) == pytest.approx(1.0 / math.sqrt(5.0))
    assert point_charge.wigner_3j(1, 1, 2, 1, -1, 0) == pytest.approx(1.0 / math.sqrt(30.0))
    assert point_charge.wigner_3j(1, 1, 2, 0, 0, 0) == pytest.approx(math.sqrt(2.0 / 15.0))
    # selection rules
    assert point_charge.wigner_3j(1, 1, 2, 1, 1, 0) == 0.0
    assert point_charge.wigner_3j(3, 3, 5, 0, 0, 0) == 0.0  # odd sum -> zero
    assert point_charge.wigner_3j(3, 3, 7, 0, 0, 0) == 0.0  # triangle rule
    # sum over all (m1, m2) triples of (j1 j2 j3; m1 m2 m3)^2 = 1 (each fixed m3
    # contributes 1/(2 j3 + 1)); an exact identity, checked to machine precision
    total = sum(
        point_charge.wigner_3j(3, 3, 6, m1, m2, -m1 - m2) ** 2
        for m1 in range(-3, 4)
        for m2 in range(-3, 4)
    )
    assert total == pytest.approx(1.0, abs=1e-14)


def test_spherical_harmonics_closed_forms():
    """Condon-Shortley complex harmonics against their textbook closed forms."""
    theta, phi = 0.83, -1.27
    assert point_charge.spherical_harmonic(0, 0, theta, phi) == pytest.approx(
        1.0 / math.sqrt(4.0 * math.pi)
    )
    assert point_charge.spherical_harmonic(1, 0, theta, phi) == pytest.approx(
        math.sqrt(3.0 / (4.0 * math.pi)) * math.cos(theta)
    )
    assert point_charge.spherical_harmonic(1, 1, theta, phi) == pytest.approx(
        -math.sqrt(3.0 / (8.0 * math.pi)) * math.sin(theta) * complex(np.cos(phi), np.sin(phi))
    )
    assert point_charge.spherical_harmonic(2, 2, theta, phi) == pytest.approx(
        math.sqrt(15.0 / (32.0 * math.pi)) * math.sin(theta) ** 2 * complex(np.cos(2 * phi), np.sin(2 * phi))
    )
    # Y_l,-m = (-1)^m conj(Y_lm), the phase that fixes the whole basis
    for m in range(1, 4):
        assert point_charge.spherical_harmonic(3, -m, theta, phi) == pytest.approx(
            (-1) ** m * np.conj(point_charge.spherical_harmonic(3, m, theta, phi))
        )


def test_real_tesseral_cartesian_and_addition_theorem():
    """Real harmonics: cartesian forms for l = 1, and the addition theorem."""
    theta, phi = 0.61, 2.4
    scale = math.sqrt(3.0 / (4.0 * math.pi))
    assert point_charge.real_tesseral(1, 1, theta, phi) == pytest.approx(
        scale * math.sin(theta) * math.cos(phi)
    )
    assert point_charge.real_tesseral(1, -1, theta, phi) == pytest.approx(
        scale * math.sin(theta) * math.sin(phi)
    )
    assert point_charge.real_tesseral(1, 0, theta, phi) == pytest.approx(scale * math.cos(theta))

    # sum_q Y*_kq(R) Y_kq(r) = sum_q Z_kq(R) Z_kq(r) = (2k+1)/(4 pi) P_k(cos gamma):
    # the identity that lets the ligand expansion be real while the basis stays complex
    theta2, phi2 = 1.9, -0.4
    cosine = math.cos(theta) * math.cos(theta2) + math.sin(theta) * math.sin(theta2) * math.cos(
        phi - phi2
    )
    legendre = {0: 1.0, 2: (3 * cosine**2 - 1) / 2, 4: (35 * cosine**4 - 30 * cosine**2 + 3) / 8,
                6: (231 * cosine**6 - 315 * cosine**4 + 105 * cosine**2 - 5) / 16}
    for k in (0, 2, 4, 6):
        complex_sum = sum(
            np.conj(point_charge.spherical_harmonic(k, q, theta, phi))
            * point_charge.spherical_harmonic(k, q, theta2, phi2)
            for q in range(-k, k + 1)
        )
        real_sum = sum(
            point_charge.real_tesseral(k, q, theta, phi)
            * point_charge.real_tesseral(k, q, theta2, phi2)
            for q in range(-k, k + 1)
        )
        expected = (2 * k + 1) / (4.0 * math.pi) * legendre[k]
        assert complex_sum.real == pytest.approx(expected, abs=1e-14)
        assert complex_sum.imag == pytest.approx(0.0, abs=1e-14)
        assert real_sum == pytest.approx(expected, abs=1e-14)


def test_gaunt_coefficients_against_numerical_quadrature():
    """Every (k, q, m, m') angular integral against a direct sphere quadrature."""
    xx, wx = np.polynomial.legendre.leggauss(80)
    phi = 2.0 * np.pi * np.arange(160) / 160
    dphi = 2.0 * np.pi / 160
    theta = np.arccos(xx)[:, None]
    weights = wx[:, None] * dphi
    worst = 0.0
    for k in (0, 2, 4, 6):
        for q in range(-k, k + 1):
            for m in range(-3, 4):
                for mp in range(-3, 4):
                    analytic = point_charge.gaunt_coefficient(3, m, k, q, 3, mp)
                    integrand = (
                        np.conj(point_charge.spherical_harmonic(3, m, theta, phi[None, :]))
                        * point_charge.spherical_harmonic(k, q, theta, phi[None, :])
                        * point_charge.spherical_harmonic(3, mp, theta, phi[None, :])
                    )
                    numeric = float(np.sum(weights * integrand.real))
                    imaginary = float(np.sum(weights * integrand.imag))
                    worst = max(worst, abs(analytic - numeric), abs(imaginary))
    # measured 5.8e-16 (float rounding of the quadrature itself)
    assert worst < 1e-12


def test_angular_matrix_against_numerical_quadrature():
    """The 7x7 Z_kq angular matrices, same quadrature, real and imaginary parts."""
    xx, wx = np.polynomial.legendre.leggauss(64)
    phi = 2.0 * np.pi * np.arange(128) / 128
    dphi = 2.0 * np.pi / 128
    theta = np.arccos(xx)[:, None]
    weights = wx[:, None] * dphi
    orbitals = np.stack(
        [
            point_charge.spherical_harmonic(3, m, theta, phi[None, :])
            for m in range(-3, 4)
        ]
    )
    worst = 0.0
    for k in (0, 2, 4, 6):
        for q in range(-k, k + 1):
            matrix = point_charge.angular_matrix(k, q)
            tesseral = point_charge.real_tesseral(k, q, theta, phi[None, :])
            for m_index in range(7):
                for mp_index in range(7):
                    integrand = np.conj(orbitals[m_index]) * tesseral * orbitals[mp_index]
                    numeric = complex(
                        float(np.sum(weights * integrand.real)),
                        float(np.sum(weights * integrand.imag)),
                    )
                    worst = max(worst, abs(matrix[m_index, mp_index] - numeric))
    assert worst < 1e-12


# --- the independent quadrature validation (core of this file) --------------


def test_radial_model_tail_is_negligible():
    """The radial model is cut off inside the nearest ion, and the tail it drops is
    negligible: the truncated weight integral differs from the analytic 8!/12^9 by
    1.2e-12 relative (measured; the closed form of the dropped tail,
    e^-48 * sum_j (8!/(8-j)!) 4^(8-j) / 12^(j+1), is 9.3e-18 = 1.2e-12 of the norm),
    so the r < R condition holds on the whole quadrature grid.  Asserted at 1e-9."""
    truncated = _radial_norm()
    assert abs(truncated - RADIAL_EXACT_NORM) / RADIAL_EXACT_NORM < 1e-9
    for _element, position, _charge in VALIDATION_IONS:
        assert math.dist(position, (0.0, 0.0, 0.0)) > R_CUT


def test_multipole_result_against_direct_quadrature():
    """The analytic matrix versus brute-force 3D integration of equation (1).

    Measured relative deviation: 3.6e-14 on the grid used here (100x48x96) and
    3.6e-14 to 4.9e-14 across grids from 60x32x64 to 200x96x192, against a
    matrix whose largest entry is 3.0e4 cm^-1 -- i.e. agreement at the level of
    the quadrature's own rounding.  Asserted at 1e-8, six orders below the 1e-6
    target, so the check stays green on a different BLAS/quadrature build
    without losing its force.
    """
    moments = _radial_moments()
    atoms = [geometry.Atom(element, *position) for element, position, _q in VALIDATION_IONS]
    charges = {element: charge for element, _p, charge in VALIDATION_IONS}
    result = estimate(atoms, charges, (0.0, 0.0, 0.0), radial=moments)
    direct = _direct_matrix(VALIDATION_IONS, _radial, _radial_norm())
    scale = np.abs(result.potential_matrix).max()
    assert scale > 1e3  # cm^-1, a physically sized 4f potential
    relative = np.abs(direct - result.potential_matrix).max() / scale
    assert relative < 1e-8, f"multipole result deviates by {relative:.3e} relative"


def test_multipole_result_agrees_with_a_lone_charge_too():
    """Same check with a single off-axis charge and an off-origin centre: no symmetry
    of the structure can hide a sign or normalisation error."""
    centre = np.array([0.1, 0.2, -0.1])
    position = np.array([1.1, -4.3, 3.2])
    charge = -0.7
    result = estimate(
        [geometry.Atom("F", *position)], {"F": charge}, tuple(centre), radial=_radial_moments()
    )
    # the quadrature works in the centre's frame, where the ion sits at position - centre
    shifted = (("F", tuple(position - centre), charge),)
    direct = _direct_matrix(shifted, _radial, _radial_norm(), nr=80, nx=40, nphi=80)
    assert np.abs(direct - result.potential_matrix).max() / np.abs(
        result.potential_matrix
    ).max() < 1e-8


# --- physics checks that need no quadrature ---------------------------------


def test_axial_charge_gives_a_matrix_diagonal_in_m():
    """A charge on the z axis leaves only q = 0 lattice sums, so V is diagonal."""
    result = estimate([geometry.Atom("O", 0.0, 0.0, 2.1)], {"O": -2.0}, (0.0, 0.0, 0.0),
                      radial={2: 0.625, 4: 0.5729, 6: 0.7241})
    matrix = result.potential_matrix
    off_diagonal = np.abs(matrix - np.diag(np.diag(matrix))).max()
    assert off_diagonal == pytest.approx(0.0, abs=1e-9)
    assert all(
        abs(result.lattice_sums[(k, q)]) < 1e-16 * abs(result.lattice_sums[(k, 0)])
        for k in (2, 4, 6)
        for q in range(-k, k + 1)
        if q != 0
    )
    # |m| = 3 lies in the xy plane, away from the axial charge: the ordering of the
    # diagonal must follow the orbital shape (m = 0 closest to the charge at +z)
    diagonal = np.diag(matrix).real
    assert diagonal[3] == max(diagonal) and diagonal[0] == min(diagonal)


def test_trace_is_the_spherical_average():
    """sum_m V_mm = -7 * e^2/(4 pi eps0) * sum_i q_i / R_i (Gauss's theorem)."""
    atoms = [geometry.Atom("O", 0.0, 0.0, 2.1), geometry.Atom("H", 1.4, 0.0, 0.0)]
    charges = {"O": -2.0, "H": 0.4}
    result = estimate(atoms, charges, (0.0, 0.0, 0.0), radial={2: 0.625, 4: 0.5729, 6: 0.7241})
    expected = -7.0 * COULOMB_CM1_ANGSTROM * (-2.0 / 2.1 + 0.4 / 1.4)
    assert result.potential_matrix.trace().real == pytest.approx(expected, rel=1e-12)
    assert result.monopole_shift == pytest.approx(-COULOMB_CM1_ANGSTROM * (-2.0 / 2.1 + 0.4 / 1.4), rel=1e-12)


def test_rotation_of_the_whole_structure_leaves_the_spectrum_invariant():
    """Rotating the ligands and the centre leaves the eigenvalues of V invariant."""
    positions = [(1.6, 0.4, -0.9), (0.2, 2.2, 1.1), (-1.9, -0.7, 0.3)]
    centre = (0.3, -0.2, 0.5)
    charges = {"O": -2.0}
    radial = {2: 0.625, 4: 0.5729, 6: 0.7241}
    base = estimate(
        [geometry.Atom("O", *position) for position in positions], charges, centre, radial=radial
    )
    rotation = _rotation([1.0, 2.0, 3.0], 0.7)
    rotated = estimate(
        [geometry.Atom("O", *rotation.dot(position)) for position in positions],
        charges,
        tuple(rotation.dot(centre)),
        radial=radial,
    )
    scale = np.abs(base.potential_eigenvalues).max()
    assert np.abs(
        np.sort(base.potential_eigenvalues) - np.sort(rotated.potential_eigenvalues)
    ).max() < 1e-9 * scale
    # the individual lattice sums do move (they carry the frame), the spectrum does not
    assert max(
        abs(base.lattice_sums[key] - rotated.lattice_sums[key]) for key in base.lattice_sums
    ) > 1e-6


def test_octahedral_cage_has_no_second_rank():
    """Six equal charges in an octahedron: k = 2 vanishes and the allowed set is Oh's."""
    result = estimate(_cage(2.3), {"O": -2.0}, (0.0, 0.0, 0.0),
                      radial={2: 0.625, 4: 0.5729, 6: 0.7241})
    reference = max(
        abs(value) for (k, _q), value in result.lattice_sums.items() if k >= 2
    )
    second_rank = max(abs(result.lattice_sums[(2, q)]) for q in range(-2, 3))
    assert second_rank < 1e-12 * reference
    assert result.nonzero_terms == ((4, 0), (4, 4), (6, 0), (6, 4))
    # 4f in an octahedral field splits as a2u + t1u + t2u (degeneracies 1, 3, 3)
    values = np.sort(result.potential_eigenvalues)
    gaps = np.diff(values)
    assert (gaps > 1e-6).sum() == 2
    assert [int((np.abs(values - level) < 1e-6).sum()) for level in (values[0], values[1], values[4])] == [1, 3, 3]


def test_parameters_and_matrix_are_consistent_with_the_definition():
    """A_k^q = -(prefactor) <r^k> L_kq, and the matrix is Hermitian."""
    moments = {2: 0.625, 4: 0.5729, 6: 0.7241}
    result = estimate(_cage(2.3), {"O": -2.0}, (0.0, 0.0, 0.0), radial=moments)
    assert result.parameters_cm1 is not None
    for (k, q), value in result.parameters_cm1.items():
        assert value == pytest.approx(
            -COULOMB_CM1_ANGSTROM * moments.get(k, 1.0) * result.lattice_sums[(k, q)], rel=1e-12
        )
    assert np.abs(result.potential_matrix - result.potential_matrix.conj().T).max() < 1e-9
    assert result.potential_units == "cm^-1"
    assert np.allclose(
        result.potential_eigenvalues, np.linalg.eigvalsh(result.potential_matrix), atol=1e-9
    )


def test_atoms_on_the_centre_are_excluded_with_a_note():
    """The caller normally passes the whole structure, including the metal."""
    atoms = _cage(2.3) + [geometry.Atom("Gd", 0.0, 0.0, 0.0)]
    result = estimate(atoms, {"O": -2.0, "Gd": 3.0}, (0.0, 0.0, 0.0))
    assert result.n_ions == 6
    assert any("self-interaction" in note for note in result.notes)


# --- the refusal paths ------------------------------------------------------


def test_missing_radial_returns_lattice_sums_only():
    """No radial moments -> geometry only: no cm^-1 parameters, explicit note."""
    result = estimate(_cage(2.3), {"O": -2.0}, (0.0, 0.0, 0.0))
    assert isinstance(result, PCEstimate)
    assert result.parameters_cm1 is None
    assert result.lattice_sums[(0, 0)] != 0.0
    assert result.potential_units == "e Angstrom^-(k+1)"
    assert any("No radial expectation values were given" in note for note in result.notes)
    assert any("not shipped with this tool" in note for note in result.notes)
    # the geometry-only matrix is the unit-radial-moment one: no cm^-1 scale
    assert np.abs(result.potential_matrix).max() < 100.0
    assert np.all(np.isfinite(result.potential_eigenvalues))
    # and the k = 0 monopole is still the spherical average, in the same units
    assert result.monopole_shift == pytest.approx(
        -math.sqrt(4.0 * math.pi) * 6 * (-2.0) / 2.3 / math.sqrt(4.0 * math.pi), rel=1e-12
    )


def test_missing_element_in_charges_is_refused():
    atoms = [geometry.Atom("Gd", 0.0, 0.0, 0.0), geometry.Atom("O", 2.1, 0.0, 0.0)]
    with pytest.raises(PointChargeError) as error:
        estimate(atoms, {"N": -3.0}, (0.0, 0.0, 0.0))
    message = str(error.value)
    assert "no charge given for O" in message
    assert "Next step: " in message
    # the refusal ends with the concrete fix, ready to copy
    assert message.endswith("{'O': -2.0}.")
    # the centre atom (Gd) is not a point charge, so it must not be listed as missing
    assert "Gd" not in message


def test_partial_charge_map_lists_only_the_missing_element():
    atoms = [geometry.Atom("Gd", 0.0, 0.0, 0.0), geometry.Atom("O", 2.1, 0.0, 0.0),
             geometry.Atom("H", 2.9, 0.0, 0.0)]
    with pytest.raises(PointChargeError, match="Next step: ") as error:
        estimate(atoms, {"O": -2.0}, (0.0, 0.0, 0.0))
    message = str(error.value)
    assert "no charge given for H " in message  # H only: O is in the map
    assert message.endswith("{'H': -2.0}.")


def test_empty_charge_map_is_refused():
    with pytest.raises(PointChargeError, match="Next step: "):
        estimate(_cage(2.3), {}, (0.0, 0.0, 0.0))


def test_incomplete_radial_map_is_refused():
    with pytest.raises(PointChargeError, match="Next step: ") as error:
        estimate(_cage(2.3), {"O": -2.0}, (0.0, 0.0, 0.0), radial={2: 0.6, 4: 0.6})
    assert "6" in str(error.value)


def test_odd_or_too_large_radial_rank_is_refused():
    with pytest.raises(PointChargeError, match="Next step: "):
        estimate(_cage(2.3), {"O": -2.0}, (0.0, 0.0, 0.0), radial={2: 0.6, 4: 0.6, 6: 0.6, 8: 0.6})
    with pytest.raises(PointChargeError, match="Next step: "):
        estimate(_cage(2.3), {"O": -2.0}, (0.0, 0.0, 0.0), radial={2: 0.6, 4: 0.6, 6: -0.6})


def test_bad_kmax_and_bad_center_are_refused():
    with pytest.raises(PointChargeError, match="Next step: "):
        estimate(_cage(2.3), {"O": -2.0}, (0.0, 0.0, 0.0), kmax=5)
    with pytest.raises(PointChargeError, match="Next step: "):
        estimate(_cage(2.3), {"O": -2.0}, (0.0, 0.0), kmax=6)
    with pytest.raises(PointChargeError, match="Next step: "):
        estimate([], {"O": -2.0}, (0.0, 0.0, 0.0))


def test_structure_without_ligands_is_refused():
    with pytest.raises(PointChargeError, match="Next step: "):
        estimate([geometry.Atom("Gd", 0.0, 0.0, 0.0)], {"Gd": 3.0}, (0.0, 0.0, 0.0))


def test_kmax_limits_the_ranks_kept():
    result = estimate(_cage(2.3), {"O": -2.0}, (0.0, 0.0, 0.0), kmax=4)
    assert result.ranks == (0, 2, 4)
    assert set(result.lattice_sums) == {
        (k, q) for k in (0, 2, 4) for q in range(-k, k + 1)
    }


# --- evidence, report and the parser-facing path ----------------------------


def test_evidence_bibkeys_exist_in_sources_bib():
    items = point_charge.evidence()
    assert len(items) >= 3
    kinds = {item.kind for item in items}
    assert "literature" in kinds and "measured" in kinds
    for item in items:
        if item.kind != "literature":
            continue
        assert item.bibkey
        entry = sources.get(item.bibkey)  # raises if the key is unknown
        assert entry.field("doi") in item.ref, item.bibkey
    keys = {item.bibkey for item in items if item.kind == "literature"}
    assert keys == {"scheie2021pycrystalfield", "ungur2017abinitio"}


def test_evidence_records_the_licence_boundary():
    """The PyCrystalField route is GPL: only the published algorithm idea is reused."""
    route = next(item for item in point_charge.evidence() if item.bibkey == "scheie2021pycrystalfield")
    assert "GPL" in route.text
    assert "no code is linked" in route.text or "not linked" in route.text


def test_report_end_to_end_on_a_hand_built_structure():
    atoms = _cage(2.3) + [geometry.Atom("Gd", 0.0, 0.0, 0.0)]
    charges = {"O": -2.0, "Gd": 3.0}
    radial = {2: 0.9847, 4: 3.5748, 6: 18.928}
    result = estimate(atoms, charges, (0.0, 0.0, 0.0), radial=radial)
    section = point_charge.report(result, charges, radial)
    assert section.title == "S2 point-charge crystal-field estimate"
    body = section.body
    for fragment in (
        "Lattice sums",
        "Crystal-field parameters",
        "4f orbital energies",
        "O -2.000 e",
        "Symmetry:",
        "ungur2017abinitio",
        "scheie2021pycrystalfield",
        "cm^-1",
        "r < R",
        "not the Stevens parameters",
    ):
        assert fragment in body, fragment
    # every (k, q) of the lattice sums is printed
    for k in (0, 2, 4, 6):
        for q in range(-k, k + 1):
            assert f"({k},{q})" in body


def test_report_flags_a_mismatched_radial_argument():
    """The numbers always come from the estimate; a disagreeing argument is a note."""
    charges = {"O": -2.0}
    radial = {2: 0.625, 4: 0.5729, 6: 0.7241}
    with_moments = estimate(_cage(2.3), charges, (0.0, 0.0, 0.0), radial=radial)
    body = point_charge.report(with_moments, charges, None).body
    assert "report() was called with radial=None" in body
    assert "<r^2> = 0.6250" in body          # the estimate's own moments are shown
    assert "Crystal-field parameters" in body

    without = estimate(_cage(2.3), charges, (0.0, 0.0, 0.0))
    other = point_charge.report(without, charges, radial).body
    assert "the estimate was built without one" in other
    assert "not produced" in other


def test_report_without_radial_states_the_refusal():
    charges = {"O": -2.0}
    result = estimate(_cage(2.3), charges, (0.0, 0.0, 0.0))
    body = point_charge.report(result, charges, None).body
    assert "not produced" in body
    assert "no radial expectation values" in body.lower()
    assert "Boundary:" in body


def test_estimate_accepts_a_parsed_xyz_structure(tmp_path):
    """The intended call path: geometry.parse_xyz -> estimate -> report."""
    xyz = tmp_path / "gd_o6.xyz"
    lines = ["7", "hand-built octahedral cage"]
    lines += [f"O {x} {y} {z}" for x, y, z in
              [(2.3, 0, 0), (-2.3, 0, 0), (0, 2.3, 0), (0, -2.3, 0), (0, 0, 2.3), (0, 0, -2.3)]]
    lines.append("Gd 0 0 0")
    xyz.write_text("\n".join(lines) + "\n", encoding="utf-8")
    atoms = geometry.parse_xyz(xyz)
    centre = atoms[geometry.dominant_center(atoms)]
    result = estimate(atoms, {"O": -2.0, "Gd": 3.0}, (centre.x, centre.y, centre.z),
                      radial={2: 0.9847, 4: 3.5748, 6: 18.928})
    assert result.n_ions == 6
    assert result.splitting > 0.0
    assert point_charge.report(result, {"O": -2.0, "Gd": 3.0}, None).title.startswith("S2")
