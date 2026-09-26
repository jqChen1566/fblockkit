"""Tests for the atomic-term check (analysis/atomic_terms.py).

The anchors are analytic: one electron in any real f component has L = 3 exactly
(<L^2> = 12), a closed-shell pair in one real harmonic has <L^2> = 24 (twice the
one-electron value, the real function carrying no angular momentum of its own),
and every expectation must be invariant under an orbital rotation of the active
space.  Hund's-rule values are pinned against the source's own table entry
(Dy3+, 6H15/2, J = 7.5).

Note that in the *real*-harmonic basis a stretched many-electron state is not a
single determinant -- occupying the components labelled +3 and +2 gives 27, not
the 30 of the M_L = L component -- so no determinant-based stretched anchor is
used here.

The ORCA-convention end of the machinery (that the exported AO block really is
the m-indexed real-harmonic set) is pinned by the fixture test at the bottom,
which needs the f6 export pair.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from fblockkit.analysis.atomic_terms import (
    AtomicTermError,
    analyze_terms,
    evidence,
    hund_term,
    project_shell,
    run,
)
from fblockkit.analysis.entropy_rdm import (
    SpinDensities,
    rotate_densities,
    solve_fci,
    spin_densities,
)
from fblockkit.analysis.solid_harmonics import component_labels
from fblockkit.parsers.fcidump import parse_fcidump
from fblockkit.parsers.orca_json import AoLabel, parse_orca_json

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"


# --- helpers -----------------------------------------------------------------


def _determinant_densities(occupied_up, occupied_down, norb):
    """Spin densities of a determinant with the given occupied orbitals per spin.

    For a determinant ``G[P,Q,R,S] = <a+_P a+_R a_S a_Q> = rho_QP rho_SR -
    rho_SP rho_RQ`` with the spin-orbital 1-RDM (zero between spins).  The
    exchange term therefore survives only inside one spin block: writing it into
    the alpha-beta blocks would be a real error, not a cosmetic one (it made
    <S^2> of a closed-shell pair come out as 0 instead of the correct 6).
    """
    gamma_a = np.zeros((norb, norb))
    gamma_b = np.zeros((norb, norb))
    for orbital in occupied_up:
        gamma_a[orbital, orbital] = 1.0
    for orbital in occupied_down:
        gamma_b[orbital, orbital] = 1.0

    def same_spin(gamma):
        return np.einsum("qp,sr->pqrs", gamma, gamma) - np.einsum("sp,rq->pqrs", gamma, gamma)

    def cross_spin(left, right):
        return np.einsum("qp,sr->pqrs", left, right)

    return SpinDensities(
        gamma_a=gamma_a,
        gamma_b=gamma_b,
        g_aa=same_spin(gamma_a),
        g_ab=cross_spin(gamma_a, gamma_b),
        g_ba=cross_spin(gamma_b, gamma_a),
        g_bb=same_spin(gamma_b),
    )


def _fake_shell_export(components: int, angular: int = 3):
    """One centre, one shell: labels in ORCA's order and an identity overlap."""
    letter = {2: "d", 3: "f"}[angular]
    labels = tuple(
        AoLabel(
            raw=f"0X  1{letter}{name}",
            center=0,
            element="X",
            shell=1,
            angular=letter,
            component=name,
        )
        for name in component_labels(angular)
    )
    return labels, np.eye(components)


# --- Hund's rules -------------------------------------------------------------


def test_hund_terms_match_the_source_table():
    dy = hund_term(3, 9)  # Dy3+, 6H15/2 in the source's Table 1
    assert (dy.s, dy.l, dy.j) == (2.5, 5.0, 7.5)
    assert dy.j_maximum == 7.5
    eu = hund_term(3, 6)  # Eu3+, 7F: less than half filled, J = |L - S| = 0
    assert (eu.s, eu.l, eu.j) == (3.0, 3.0, 0.0)
    assert eu.j_maximum == 6.0
    gd = hund_term(3, 7)  # 8S7/2
    assert (gd.s, gd.l, gd.j) == (3.5, 0.0, 3.5)
    d5 = hund_term(2, 5)  # the half-filled d shell, 6S
    assert (d5.s, d5.l, d5.j) == (2.5, 0.0, 2.5)


def test_impossible_occupations_raise():
    with pytest.raises(AtomicTermError, match="Hund term"):
        hund_term(3, 0)
    with pytest.raises(AtomicTermError, match="Hund term"):
        hund_term(3, 14)


# --- analytic anchors through the whole pipeline ------------------------------


def test_one_electron_has_l_squared_12_whatever_the_component():
    labels, overlap = _fake_shell_export(7)
    for component in range(7):
        densities = _determinant_densities([component], [], 7)
        result = analyze_terms(densities, np.eye(7), overlap, labels, center=0, angular=3)
        assert result.l2 == pytest.approx(12.0, abs=1e-12)
        assert result.l_eff == pytest.approx(3.0, abs=1e-12)
        assert result.s2 == pytest.approx(0.75, abs=1e-12)
        assert result.n_electrons == pytest.approx(1.0, abs=1e-12)


def test_a_closed_shell_pair_has_l_squared_24():
    # two electrons in one real harmonic: <L^2> = 2 <phi|L^2|phi> + 2 sum_alpha
    # <phi|L_alpha|phi>^2 = 2*12 + 0 (a real function has no orbital angular
    # momentum of its own) = 24.  Note that in the *real* basis a stretched
    # many-electron state is NOT a single determinant -- that is why this anchor
    # is the closed-shell pair and not the M_L = L component.
    labels, overlap = _fake_shell_export(7)
    for component in range(7):
        densities = _determinant_densities([component], [component], 7)
        result = analyze_terms(densities, np.eye(7), overlap, labels, center=0, angular=3)
        assert result.l2 == pytest.approx(24.0, abs=1e-12)
        assert result.s2 == pytest.approx(0.0, abs=1e-12)
        assert result.j2 == pytest.approx(24.0, abs=1e-12)


def test_the_quantum_numbers_are_invariant_under_an_orbital_rotation():
    # rotate the active orbitals (and the density objects with them): every
    # expectation is basis independent, which is the invariance criterion the
    # project prefers over external comparisons
    labels, overlap = _fake_shell_export(7)
    densities = _determinant_densities([0, 3], [1], 7)
    rng = np.random.default_rng(11)
    matrix = rng.normal(size=(7, 7))
    rotation, _ = np.linalg.qr(matrix)
    rotated = rotate_densities(densities, rotation)
    plain = analyze_terms(densities, np.eye(7), overlap, labels, center=0, angular=3)
    turned = analyze_terms(rotated, np.eye(7) @ rotation, overlap, labels, center=0, angular=3)
    for name in ("l2", "s2", "ls", "j2", "n_electrons"):
        assert getattr(turned, name) == pytest.approx(getattr(plain, name), abs=1e-10)


def test_a_non_stretched_component_does_not_raise_a_false_alarm():
    # a state with S_eff = L_eff = the Hund values but no stretched M_J must not be
    # reported as a defect: in a scalar calculation J is not yet defined
    labels, overlap = _fake_shell_export(7)
    densities = _determinant_densities([0, 1, 2, 3, 4, 5], [0, 1, 2, 3, 4, 5], 7)
    result = analyze_terms(densities, np.eye(7), overlap, labels, center=0, angular=3)
    # six doubly occupied real f harmonics do not form the Hund term, so this is a
    # NOT-in-the-Hund-term case rather than a J complaint
    assert any("NOT in the Hund term" in verdict for verdict in result.verdicts)
    assert not any("must not be used" in verdict for verdict in result.verdicts)


def test_the_s2_cross_check_refuses_inconsistent_densities():
    labels, overlap = _fake_shell_export(7)
    densities = _determinant_densities([0], [], 7)
    with pytest.raises(AtomicTermError, match="disagree"):
        analyze_terms(
            densities, np.eye(7), overlap, labels, center=0, angular=3, s2_reference=5.0
        )
    # the truthful reference passes and is recorded
    result = analyze_terms(
        densities, np.eye(7), overlap, labels, center=0, angular=3, s2_reference=0.75
    )
    assert any("S^2" in check for check in result.checks)


def test_a_d_shell_projection_is_supported():
    labels, overlap = _fake_shell_export(5, angular=2)
    densities = _determinant_densities([0, 1, 2, 3, 4], [], 5)  # five parallel d electrons
    result = analyze_terms(densities, np.eye(5), overlap, labels, center=0, angular=2)
    assert result.s_eff == pytest.approx(2.5, abs=1e-12)
    assert result.n_electrons == pytest.approx(5.0, abs=1e-12)
    assert result.hund is not None  # d5 is a valid configuration for the table


def test_a_mixed_space_withholds_the_hund_verdict():
    # two of three electrons in the f shell, one in an s-like orbital: the mean
    # character falls below the provisional line
    labels, overlap = _fake_shell_export(7)
    densities = _determinant_densities([0, 1], [0], 7)
    coefficients = np.zeros((7, 7))
    coefficients[:6, :6] = np.eye(6)
    result = analyze_terms(densities, coefficients, overlap, labels, center=0, angular=3)
    body = run(result).body
    assert "withheld" in body
    assert result.hund is None


# --- projection mechanics -----------------------------------------------------


def test_project_shell_reports_the_character_of_each_orbital():
    labels, overlap = _fake_shell_export(7)
    coefficients = np.zeros((7, 2))
    coefficients[:, 0] = np.array([1, 0, 0, 0, 0, 0, 0])
    coefficients[:, 1] = np.array([1, 1, 0, 0, 0, 0, 0]) / math.sqrt(2)
    projection = project_shell(coefficients, overlap, labels, center=0, angular=3)
    assert projection.characters == pytest.approx((1.0, 1.0), abs=1e-12)
    assert projection.n_shells == 1
    assert len(projection.amplitudes[0]) == 7


def test_project_shell_raises_when_the_centre_has_no_such_shell():
    labels, overlap = _fake_shell_export(7)
    with pytest.raises(AtomicTermError, match="no l=3"):
        project_shell(np.eye(7), overlap, labels, center=1, angular=3)


def test_project_shell_raises_when_the_component_spelling_is_broken():
    labels, overlap = _fake_shell_export(7)
    broken = tuple(
        AoLabel(
            raw=label.raw,
            center=label.center,
            element=label.element,
            shell=label.shell,
            angular=label.angular,
            component=("+9" if index == 3 else label.component),
        )
        for index, label in enumerate(labels)
    )
    with pytest.raises(AtomicTermError, match="expected exactly one"):
        project_shell(np.eye(7), overlap, broken, center=0, angular=3)


def test_project_shell_raises_when_the_block_is_not_r_times_identity():
    labels, overlap = _fake_shell_export(7)
    broken = overlap.copy()
    broken[0, 0] = 2.0  # breaks the measured R (x) I structure
    with pytest.raises(AtomicTermError, match="not R"):
        project_shell(np.eye(7), broken, labels, center=0, angular=3)


def test_evidence_carries_the_source_entry():
    text = " ".join(entry.ref + " " + entry.text for entry in evidence())
    assert "10.1021/acs.jpclett.5c02971" in text
    assert "7.5" in text


# --- the real-data pin (needs the f6 fixture pair) ----------------------------

F6_DUMP = FIXTURES / "eu3_f6_casscf.fcidump"
F6_JSON = FIXTURES / "eu3_f6_casscf.canonical.json"


@pytest.mark.skipif(
    not (F6_DUMP.exists() and F6_JSON.exists()),
    reason="the f6 fixture pair is not in the repository yet",
)
def test_the_eu_f6_fixture_is_the_7f_term():
    dump = parse_fcidump(F6_DUMP)
    export = parse_orca_json(F6_JSON)
    state = solve_fci(dump, reference_energy=-10826.6, multiplicity=7)
    densities = spin_densities(state)
    window = range(dump.norb)  # the fixture dump carries only the active block
    coefficients = np.array(export.mo_coefficients).T[:, list(window)]
    result = analyze_terms(
        densities,
        coefficients,
        np.array(export.overlap),
        export.ao_labels,
        center=0,
        angular=3,
        s2_reference=state.s2,
    )
    assert result.l_eff == pytest.approx(3.0, abs=0.02)
    assert result.s_eff == pytest.approx(3.0, abs=0.02)
    assert result.j_eff == pytest.approx(6.0, abs=0.05)  # the stretched 7F component
    assert result.hund is not None and result.hund.l == 3.0
