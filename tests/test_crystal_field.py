"""Crystal-field (A4) tests: symmetry rules, Stevens operators, linear fit.

Evidence chain exercised here, from the weakest to the strongest:

1. self-checks that follow from the operator definitions (``O_2^0`` closed form,
   Hermiticity, zero trace, ``[J_z, O_k^q] = i q O_k^{-q}``, the textbook D/E
   zero-field-splitting identity);
2. a synthetic round trip on the fixture: diagonalise ``H_CF``, feed the
   eigenpairs back, recover the nine ``B_k^q`` -- plus the measured reason why
   the plain version of that recipe cannot work for a Kramers ion;
3. an **external known answer**: the paper's own two tables must agree with each
   other, i.e. diagonalising the tabulated ``B_k^q`` must reproduce the tabulated
   level energies.  Fixtures: ``fixtures/literature/peng_s1_cf.json`` (C3, Er3+)
   and ``fixtures/literature/peng_s2_s7_oh_cf.json`` (Oh, Dy3+).

Three of those four published columns reproduce to the printing precision of the
table; the fourth is off by a constant factor, which is recorded as a
measurement rather than smoothed over (see the fixture's ``cross_check`` block).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from fblockkit.analysis import crystal_field, geometry
from fblockkit.knowledge import sources
from fblockkit.knowledge.models import EVIDENCE_LITERATURE

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "literature"

#: The two tables are printed to 0.01 cm^-1, so agreement means a deviation
#: below half of the last printed digit.
PRINTING_TOLERANCE = 0.005

#: Half-integer total angular momentum of both fixture ions (Er3+ 4I15/2,
#: Dy3+ 6H15/2): 16 states, eight Kramers doublets, nine C3 parameters.
J_ER = 7.5
J_DY = 7.5


def _fixture(name: str) -> dict:
    return json.loads((FIXTURES / name).read_text(encoding="utf-8"))


def _parameters(fixture: dict, column: str) -> dict[tuple[int, int], float]:
    return {
        (record["k"], record["q"]): record["values"][column]
        for record in fixture["parameters"]["records"]
    }


def _levels(fixture: dict, column: str) -> np.ndarray:
    return np.array([record["values"][column] for record in fixture["levels"]["records"]])


def _doublet_spectrum(parameters, J: float) -> np.ndarray:
    """Ground-manifold levels of ``H_CF``, ground level shifted to zero, Kramers
    doublet members collapsed into one value per doublet."""
    eigenvalues = np.sort(np.linalg.eigvalsh(crystal_field.hamiltonian(parameters, J)))
    return (eigenvalues - eigenvalues[0])[::2]


def _tilted_states(J: float, angles) -> np.ndarray:
    """Deterministic oriented (symmetry-broken) states ``R_z(phi) R_y(theta)|J,J>``.

    The stand-in for the reference method's Haar-random orientation sampling: a
    rotated stretched state is a single-determinant-like state with a
    well-defined orientation, and it is *not* an eigenstate of ``H_CF``.  Built
    here from scratch (ladder operators -> J_y -> matrix exponential) so that it
    does not lean on anything the module under test provides.
    """
    dimension = int(round(2 * J + 1))
    m_values = np.arange(dimension) - J
    j_plus = np.zeros((dimension, dimension), dtype=complex)
    for index in range(dimension - 1):
        j_plus[index + 1, index] = np.sqrt(J * (J + 1) - m_values[index] * (m_values[index] + 1))
    j_y = (j_plus - j_plus.conj().T) / 2j
    values, vectors = np.linalg.eigh(j_y)
    stretched = np.zeros(dimension, dtype=complex)
    stretched[-1] = 1.0
    states = []
    for theta, phi in angles:
        rotate_y = (vectors * np.exp(-1j * theta * values)) @ vectors.conj().T
        rotate_z = np.diag(np.exp(1j * phi * m_values))
        states.append(rotate_z @ rotate_y @ stretched)
    return np.array(states)


# --- symmetry rules ---------------------------------------------------------


def test_allowed_parameters_of_the_three_cases():
    """C3 -> 9, Oh -> 4 (second rank zero), no symmetry -> 27."""
    c3 = crystal_field.allowed_parameters("C3")
    assert len(c3) == 9
    assert c3 == tuple(sorted(c3))  # documented order: k, then q
    assert {k for k, _ in c3} == {2, 4, 6}
    assert {q for _, q in c3} == {0, 3, -3, 6, -6}  # three-fold axis: q multiple of 3
    assert all(abs(q) <= k for k, q in c3)

    oh = crystal_field.allowed_parameters("Oh")
    assert oh == ((4, 0), (4, 4), (6, 0), (6, 4))
    assert all(k != 2 for k, _ in oh)  # no second-rank term survives cubic symmetry

    assert len(crystal_field.allowed_parameters("none")) == 27
    assert crystal_field.allowed_parameters("C1") == crystal_field.allowed_parameters("none")
    assert crystal_field.allowed_parameters("c3") == c3  # spelling is normalised


def test_allowed_parameter_counts_agree_with_geometry_prompt():
    """S1 (geometry) and A4 (crystal field) must not drift apart on this rule."""
    assert geometry.CF_PARAM_COUNTS["C3"] == len(crystal_field.allowed_parameters("C3"))
    assert geometry.CF_PARAM_COUNTS["Oh"] == len(crystal_field.allowed_parameters("Oh"))
    assert geometry.CF_PARAM_COUNTS["none"] == len(crystal_field.allowed_parameters("none"))


def test_allowed_parameters_rejects_unknown_point_group():
    with pytest.raises(crystal_field.CrystalFieldError, match="Next step:"):
        crystal_field.allowed_parameters("D4h")


def test_fixture_nonzero_pattern_is_the_c3_allowed_set():
    """The paper's Table S1 is a worked instance of the symmetry rule."""
    fixture = _fixture("peng_s1_cf.json")
    pattern = tuple((entry[0], entry[1]) for entry in fixture["nonzero_pattern_c3"])
    assert set(pattern) == set(crystal_field.allowed_parameters("C3"))
    assert len(pattern) == len(crystal_field.allowed_parameters("C3"))

    for column in fixture["parameters"]["columns"]:
        non_zero = {key for key, value in _parameters(fixture, column).items() if value != 0.0}
        assert non_zero == set(pattern), column

    oh_fixture = _fixture("peng_s2_s7_oh_cf.json")
    oh_pattern = tuple((entry[0], entry[1]) for entry in oh_fixture["allowed_parameters_oh"])
    assert set(oh_pattern) == set(crystal_field.allowed_parameters("Oh"))
    assert {key for key, value in _parameters(oh_fixture, "HF@HF").items() if value != 0.0} == set(
        oh_pattern
    )


# --- extension Stevens operators --------------------------------------------


@pytest.mark.parametrize("J", [7.5, 6.0, 1.0])
def test_stevens_second_rank_closed_form(J):
    """O_2^0 = 3 J_z^2 - J(J+1) I, exactly (the source table's k=2, q=0 row)."""
    dimension = int(round(2 * J + 1))
    matrix = crystal_field.stevens_matrices(J)[(2, 0)]
    expected = 3.0 * np.diag(np.arange(dimension) - J) ** 2 - J * (J + 1) * np.eye(dimension)
    assert np.array_equal(matrix, expected)


@pytest.mark.parametrize("J", [7.5, 6.0])
def test_stevens_operators_are_hermitian_with_zero_trace(J):
    matrices = crystal_field.stevens_matrices(J)
    assert len(matrices) == 27
    for key, matrix in matrices.items():
        scale = np.abs(matrix).max()
        assert np.abs(matrix - matrix.conj().T).max() <= 1e-12 * scale, key
        assert abs(np.trace(matrix)) <= 1e-12 * scale * matrix.shape[0], key
        # cosine components real symmetric, sine components purely imaginary
        # antisymmetric (the convention this module documents)
        if key[1] > 0:
            assert np.abs(matrix.imag).max() <= 1e-12 * scale, key
        elif key[1] < 0:
            assert np.abs(matrix.real).max() <= 1e-12 * scale, key


def test_stevens_phase_relation_to_the_partner_operator():
    """[J_z, O_k^q] = i q O_k^{-q}: the phase link between the two signs of q.

    This is the fingerprint of the cosine/sine convention.  Note that the
    tesseral operators are *not* eigenoperators of ``ad(J_z)``, which is why the
    relation carries the partner ``O_k^{-q}``; a sign slip in ``c_- = 1/(2i)``
    shows up here immediately.
    """
    J = 7.5
    dimension = int(round(2 * J + 1))
    j_z = np.diag(np.arange(dimension) - J).astype(complex)
    matrices = crystal_field.stevens_matrices(J)
    for (k, q), matrix in matrices.items():
        partner = matrices[(k, -q)] if q != 0 else np.zeros_like(matrix)
        lhs = j_z @ matrix - matrix @ j_z
        assert np.abs(lhs - 1j * q * partner).max() <= 1e-12 * np.abs(matrix).max(), (k, q)


def test_stevens_second_rank_matches_the_textbook_dE_form():
    """Independent algebraic check of the k = 2 operators.

    ``D (J_z^2 - J(J+1)/3) + E (J_x^2 - J_y^2)`` and ``(D/3) O_2^0 + E O_2^2``
    are two spellings of the same operator: ``J_x^2 - J_y^2 = (J_+^2 + J_-^2)/2``
    and ``O_2^2`` is exactly that combination.  Deriving this inside the test
    keeps the check independent of the coefficients tabulated from the source.
    """
    J = 7.5
    dimension = int(round(2 * J + 1))
    m_values = np.arange(dimension) - J
    j_plus = np.zeros((dimension, dimension), dtype=complex)
    for index in range(dimension - 1):
        j_plus[index + 1, index] = np.sqrt(J * (J + 1) - m_values[index] * (m_values[index] + 1))
    j_z = np.diag(m_values).astype(complex)
    j_x = (j_plus + j_plus.conj().T) / 2
    j_y = (j_plus - j_plus.conj().T) / 2j
    matrices = crystal_field.stevens_matrices(J)
    D, E = 1.7, 0.31
    textbook = D * (j_z @ j_z - J * (J + 1) / 3 * np.eye(dimension)) + E * (j_x @ j_x - j_y @ j_y)
    from_module = (D / 3) * matrices[(2, 0)] + E * matrices[(2, 2)]
    assert np.abs(textbook - from_module).max() <= 1e-12 * np.abs(textbook).max()


@pytest.mark.parametrize("J", [7.5, 5.5])
def test_pure_axial_spectrum_is_kramers_paired(J):
    """Half-integer J: a pure B_2^0 spectrum is paired as +-M, every level
    exactly two-fold degenerate (Kramers' theorem).

    The pairing is the M <-> -M one, because ``O_2^0 = 3 J_z^2 - J(J+1)``
    depends on ``M^2`` only.  Note what this does *not* say: the level values
    are not mirror-symmetric about their mean.  They run from ``-J(J+1)`` up to
    ``2J(J+1)``, an asymmetric range, so the set is only centred in the trace
    sense -- which the last assertion checks.
    """
    matrix = crystal_field.stevens_matrices(J)[(2, 0)]
    eigenvalues = np.sort(np.linalg.eigvalsh(matrix))
    assert eigenvalues.size == 2 * J + 1
    for index in range(0, eigenvalues.size, 2):
        assert eigenvalues[index + 1] - eigenvalues[index] < 1e-9
    assert np.unique(np.round(eigenvalues, 9)).size == J + 0.5
    assert abs(eigenvalues.sum()) < 1e-9  # O_2^0 is traceless


@pytest.mark.parametrize("J", [7.5, 6.0])
def test_pure_axial_spectrum_matches_the_closed_form(J):
    """Analytic oracle for a pure B_2^0 term, scale included.

    ``O_2^0 = 3 J_z^2 - J(J+1)`` is diagonal in ``|J M>``, so the whole spectrum
    is ``B_2^0 * (3 M^2 - J(J+1))``.  Note the set is *not* mirror-symmetric
    under ``lambda -> -lambda`` (it runs from ``-B*J(J+1)`` to ``+2 B*J(J+1)``);
    it is only centred, because the operator is traceless.
    """
    B2 = -1.0846603  # the fixture's leading axial parameter
    eigenvalues = np.sort(np.linalg.eigvalsh(crystal_field.hamiltonian({(2, 0): B2}, J)))
    m_values = np.arange(int(round(2 * J + 1))) - J
    assert np.allclose(eigenvalues, np.sort(B2 * (3 * m_values**2 - J * (J + 1))), rtol=1e-12)


def test_pure_axial_spectrum_of_an_integer_J_ion_keeps_one_singlet():
    """Non-Kramers contrast: with integer J the M = 0 state stays unpaired.

    J = 8 is Ho3+ (5I8, the paper's compound 5), whose 17 levels are reported as
    eight quasi-doublets plus one singlet -- the same structure the pure axial
    operator shows here.
    """
    J = 8.0
    eigenvalues = np.sort(np.linalg.eigvalsh(crystal_field.stevens_matrices(J)[(2, 0)]))
    assert eigenvalues.size == 17
    assert eigenvalues[0] == pytest.approx(-J * (J + 1), abs=1e-9)  # the M = 0 state
    assert eigenvalues[1] - eigenvalues[0] > 1e-9  # ... and it is unpaired
    for index in range(1, eigenvalues.size, 2):
        assert eigenvalues[index + 1] - eigenvalues[index] < 1e-9


# --- fitting ----------------------------------------------------------------


def test_round_trip_recovers_the_table_s1_parameters():
    """Round trip on both mean-field columns of the fixture.

    Measured on this fixture: the plain recipe (diagonalise, feed the eigenpairs
    back) is **rank deficient** -- 16 eigenstates of a Kramers manifold give rank
    8 of 10 columns and cond ~ 2e18, because the two members of a doublet have
    equal expectation values for every time-even operator, so each doublet
    contributes one equation instead of two.  The fit still reproduces the
    energies (that is what least squares does with a rank-deficient matrix) while
    the parameters wander.  Adding oriented states -- the reference method's
    sampling device -- restores full rank and the nine B_k^q come back.
    """
    fixture = _fixture("peng_s1_cf.json")
    for column in fixture["parameters"]["columns"][:2]:
        expected = _parameters(fixture, column)
        true = {key: value for key, value in expected.items() if value != 0.0}
        hamiltonian = crystal_field.hamiltonian(true, J_ER)
        levels, vectors = np.linalg.eigh(hamiltonian)
        # rows = coordinate vectors of the eigenstates (columns of `vectors`),
        # i.e. C_{M,i} in Eq. (3).  Note it is NOT conj(eigenvector): the
        # expectation value C^dagger O C is evaluated on the coordinate vector
        # itself, and a conjugated coordinate vector is a different state.
        coefficients = vectors.T

        plain = crystal_field.fit_crystal_field(levels, coefficients, "C3", J_ER)
        assert plain.rank_deficient and plain.ill_conditioned
        assert plain.max_abs_residual < 1e-9  # ... the energies are reproduced
        worst_plain = max(
            abs(plain.parameters[key] - value) / abs(value) for key, value in true.items()
        )
        assert worst_plain > 1e-3  # ... yet the parameters are not determined
        assert any("rank deficient" in note for note in plain.notes)

        extra_coefficients = _tilted_states(
            J_ER, [(0.35 * (index + 1), 2 * np.pi * index / 6) for index in range(6)]
        )
        assert np.allclose(np.linalg.norm(extra_coefficients, axis=1), 1.0)
        extra_levels = np.array([v.conj() @ hamiltonian @ v for v in extra_coefficients]).real

        result = crystal_field.fit_crystal_field(
            np.concatenate([levels, extra_levels]),
            np.vstack([coefficients, extra_coefficients]),
            "C3",
            J_ER,
        )
        assert result.n_parameters == 9 and result.n_states == levels.size + 6
        assert not result.rank_deficient and not result.ill_conditioned
        for key, value in true.items():
            assert result.parameters[key] == pytest.approx(value, rel=1e-8), (column, key)
        assert result.const == pytest.approx(0.0, abs=1e-8)
        assert result.max_abs_residual < 1e-9
        # step 6 of the paper's workflow: spectrum of the fitted Hamiltonian
        assert np.allclose(result.spectrum(), np.sort(levels), atol=1e-8)


def test_kramers_doublet_members_give_identical_design_rows():
    """The mechanism behind the rank deficiency, measured directly.

    Two members of a Kramers doublet have the same expectation value for every
    time-even operator, so their design-matrix rows coincide to machine
    precision -- one equation per doublet, not two.  (Max relative deviation
    1.0e-15 / 2.7e-15 on the two fixture columns.)
    """
    fixture = _fixture("peng_s1_cf.json")
    for column in fixture["parameters"]["columns"][:2]:
        true = {key: value for key, value in _parameters(fixture, column).items() if value != 0.0}
        hamiltonian = crystal_field.hamiltonian(true, J_ER)
        _, vectors = np.linalg.eigh(hamiltonian)
        matrix, _ = crystal_field.design_matrix(vectors.T, "C3", J_ER)
        scale = np.abs(matrix).max()
        for index in range(0, matrix.shape[0], 2):  # the doublets come out adjacent
            assert np.abs(matrix[index] - matrix[index + 1]).max() <= 1e-12 * scale


def test_conjugating_the_sine_parameters_leaves_the_spectrum_unchanged():
    """Why the coordinate-vector convention is easy to get wrong.

    ``C[i, m]`` is the amplitude *itself*, so the row is ``C_i^dagger O C_i``.
    Handing in conjugated coordinates (or reading ``vectors.conj().T`` from an
    eigensolver whose eigenvectors are columns) evaluates a different
    Hamiltonian: conjugation maps ``O_k^q -> O_k^q`` for ``q >= 0`` but
    ``O_k^q -> -O_k^q`` for ``q < 0``.  The two Hamiltonians are different
    operators (here by 1.5e2) yet have *bit-identical* spectra, so a fit on such
    data still reaches tiny residuals -- with the sine-type parameters wrong.
    Measured: max |eig difference| = 0.0.
    """
    fixture = _fixture("peng_s1_cf.json")
    true = {key: value for key, value in _parameters(fixture, "PBE0@HF").items() if value != 0.0}
    conjugated = {key: (value if key[1] >= 0 else -value) for key, value in true.items()}
    direct = crystal_field.hamiltonian(true, J_ER)
    flipped = crystal_field.hamiltonian(conjugated, J_ER)
    assert np.abs(direct - flipped).max() > 1.0  # different operators ...
    assert np.array_equal(np.sort(np.linalg.eigvalsh(direct)), np.sort(np.linalg.eigvalsh(flipped)))

    # ... so the design matrix built from conjugated coordinates is inconsistent
    # with the true Hamiltonian, while the coordinate-vector one is exact.
    levels, vectors = np.linalg.eigh(direct)
    columns = crystal_field.allowed_parameters("C3")
    vector = np.array([true[key] for key in columns] + [0.0])
    exact, _ = crystal_field.design_matrix(vectors.T, "C3", J_ER)
    trapped, _ = crystal_field.design_matrix(vectors.conj().T, "C3", J_ER)
    assert np.abs(exact @ vector - levels).max() < 1e-9
    assert np.abs(trapped @ vector - levels).max() > 1.0


def test_known_answer_diagonalising_table_s1_reproduces_table_s6():
    """External known answer: the paper's two tables must agree with each other.

    Building ``H_CF`` from the nine tabulated PBE0@HF parameters and
    diagonalising gives the eight Kramers doublets of Table S6 to the precision
    to which that table is printed.  This is the strongest check available here:
    it exercises the operator normalisation, the cosine/sine phase convention and
    the sign of every individual B_k^q at once, against an independent published
    calculation rather than against our own reconstruction.
    """
    fixture = _fixture("peng_s1_cf.json")
    tabulated = _levels(fixture, "PBE0@HF")
    computed = _doublet_spectrum(_parameters(fixture, "PBE0@HF"), J_ER)
    assert np.abs(computed - tabulated).max() < PRINTING_TOLERANCE


def test_known_answer_second_system_is_cubic():
    """Same external check for the Oh compound: four parameters, both columns."""
    fixture = _fixture("peng_s2_s7_oh_cf.json")
    for column in fixture["levels"]["columns"]:
        tabulated = _levels(fixture, column)
        computed = _doublet_spectrum(_parameters(fixture, column), J_DY)
        assert np.abs(computed - tabulated).max() < PRINTING_TOLERANCE, column


def test_table_s6_hf_column_is_off_by_a_constant_factor():
    """Recorded inconsistency in the source, locked so it cannot pass unnoticed.

    The HF@HF pair of the paper's two tables disagrees by a constant factor
    (measured 0.76336, uniform over all seven excited doublets to within the
    table's rounding).  Since the PBE0@HF pair of the same tables agrees exactly,
    and both columns of the cubic fixture agree exactly, this is not an
    operator-convention effect -- a wrong convention would not leave a single
    uniform factor behind.  See the fixture's ``cross_check`` block; the
    consequence for users is that only the PBE0@HF pair is usable as-is.
    """
    fixture = _fixture("peng_s1_cf.json")
    computed = _doublet_spectrum(_parameters(fixture, "HF@HF"), J_ER)
    tabulated = _levels(fixture, "HF@HF")
    ratio = computed[1:] / tabulated[1:]
    assert np.allclose(ratio, ratio.mean(), rtol=1e-4)  # uniform, not a fit residual
    assert ratio.mean() == pytest.approx(0.76336, abs=1e-4)


def test_fit_is_exact_when_the_states_are_eigenstates():
    """A non-degenerate sample is an exact known answer for the fitter itself."""
    parameters = {(2, 0): 1.25, (4, 0): -0.04, (4, 3): 0.02, (6, 0): 1e-3}
    all_parameters = {key: 0.0 for key in crystal_field.allowed_parameters("C3")}
    all_parameters.update(parameters)
    hamiltonian = crystal_field.hamiltonian(all_parameters, 6.0, const=-3.5)
    levels, vectors = np.linalg.eigh(hamiltonian)
    result = crystal_field.fit_crystal_field(levels, vectors.T, "C3", 6.0)
    assert result.const == pytest.approx(-3.5, abs=1e-10)
    for key, value in all_parameters.items():
        assert result.parameters[key] == pytest.approx(value, abs=1e-10), key


def test_fit_without_the_constant_column():
    """Cubic case: four columns, so eight eigenstates -- the eight independent
    equations a Kramers manifold provides -- already determine them.  This is the
    same asymmetry the paper reports between its two benchmark compounds: an Oh
    system with four parameters converges with far fewer sampled determinants
    than its C3 system with nine.
    """
    parameters = {(4, 0): -0.5, (4, 4): 0.125, (6, 0): 2e-3, (6, 4): -1e-3}
    hamiltonian = crystal_field.hamiltonian(parameters, J_DY)
    levels, vectors = np.linalg.eigh(hamiltonian)
    result = crystal_field.fit_crystal_field(
        levels, vectors.T, "Oh", J_DY, with_const=False
    )
    assert result.const == 0.0
    assert result.n_parameters == 4
    for key, value in parameters.items():
        assert result.parameters[key] == pytest.approx(value, rel=1e-8), key


def test_design_matrix_shape_and_column_order():
    matrix, columns = crystal_field.design_matrix(
        np.eye(16, dtype=complex), "C3", J_ER, with_const=True
    )
    assert columns == crystal_field.allowed_parameters("C3")
    assert matrix.shape == (16, 10)
    assert np.allclose(matrix[:, -1], 1.0)
    without_const, _ = crystal_field.design_matrix(
        np.eye(16, dtype=complex), "C3", J_ER, with_const=False
    )
    assert without_const.shape == (16, 9)


# --- diagnostics and error paths --------------------------------------------


def test_underdetermined_fit_is_flagged():
    """Fewer sampled states than fitted parameters (C3: 9 B_k^q plus the constant)."""
    coefficients = np.eye(16, dtype=complex)[:4]
    levels = np.array([0.0, 1.0, 2.0, 3.0])
    result = crystal_field.fit_crystal_field(levels, coefficients, "C3", J_ER)
    assert result.underdetermined
    assert any("underdetermined" in note for note in result.notes)


def test_unnormalised_coefficients_are_flagged():
    parameters = {key: 0.0 for key in crystal_field.allowed_parameters("C3")}
    parameters[(2, 0)] = 0.5
    hamiltonian = crystal_field.hamiltonian(parameters, J_ER)
    levels, vectors = np.linalg.eigh(hamiltonian)
    result = crystal_field.fit_crystal_field(levels, 1.1 * vectors.T, "C3", J_ER)
    assert any("norms deviate" in note for note in result.notes)


def test_well_posed_fit_says_so():
    parameters = {key: 0.0 for key in crystal_field.allowed_parameters("C3")}
    parameters[(2, 0)] = 0.5
    hamiltonian = crystal_field.hamiltonian(parameters, J_ER)
    levels, vectors = np.linalg.eigh(hamiltonian)
    extra = _tilted_states(J_ER, [(0.4 * (index + 1), index) for index in range(6)])
    extra_levels = np.array([v.conj() @ hamiltonian @ v for v in extra]).real
    result = crystal_field.fit_crystal_field(
        np.concatenate([levels, extra_levels]),
        np.vstack([vectors.T, extra]),
        "C3",
        J_ER,
    )
    assert result.notes == (
        "well posed: full column rank and cond below the ill-conditioning limit.",
    )


def test_mismatched_coefficient_dimension_reports_next_step():
    with pytest.raises(crystal_field.CrystalFieldError, match="Next step:"):
        crystal_field.fit_crystal_field([0.0] * 4, np.eye(12, dtype=complex)[:4], "C3", J_ER)


def test_mismatched_level_count_reports_next_step():
    with pytest.raises(crystal_field.CrystalFieldError, match="Next step:"):
        crystal_field.fit_crystal_field([0.0] * 5, np.eye(16, dtype=complex), "C3", J_ER)


def test_non_finite_level_reports_next_step():
    levels = np.zeros(16)
    levels[3] = np.nan
    with pytest.raises(crystal_field.CrystalFieldError, match="Next step:"):
        crystal_field.fit_crystal_field(levels, np.eye(16, dtype=complex), "C3", J_ER)


def test_empty_sample_reports_next_step():
    with pytest.raises(crystal_field.CrystalFieldError, match="Next step:"):
        crystal_field.fit_crystal_field([], np.zeros((0, 16), dtype=complex), "C3", J_ER)


def test_non_finite_coefficients_report_next_step():
    coefficients = np.eye(16, dtype=complex)
    coefficients[7, 3] = np.nan
    with pytest.raises(crystal_field.CrystalFieldError, match="Next step:"):
        crystal_field.fit_crystal_field(np.zeros(16), coefficients, "C3", J_ER)


@pytest.mark.parametrize("bad_J", [0.5, 0.0, 7.3, -7.5])
def test_invalid_J_reports_next_step(bad_J):
    with pytest.raises(crystal_field.CrystalFieldError, match="Next step:"):
        crystal_field.stevens_matrices(bad_J)


def test_unknown_operator_index_reports_next_step():
    with pytest.raises(crystal_field.CrystalFieldError, match="Next step:"):
        crystal_field.hamiltonian({(3, 0): 1.0}, J_ER)


# --- provenance -------------------------------------------------------------


def test_evidence_is_wired_to_the_bibliography():
    """Same wiring rule the other layers are held to (see tests/test_sources.py)."""
    items = crystal_field.evidence()
    assert {item.kind for item in items} == {"literature", "manual", "measured"}
    for item in items:
        if item.kind == EVIDENCE_LITERATURE:
            entry = sources.get(item.bibkey)  # raises if the key is unknown
            assert entry.field("doi") in item.ref
    manual = [item for item in items if item.kind == "manual"]
    assert manual and manual[0].url.startswith("https://easyspin.org/")
    assert "Table 1" in manual[0].ref
