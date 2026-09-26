"""Tests for the environment spin-polarisation entropy (analysis/environment_spin.py).

The analytic anchors: Delta S_E = 0 for any spin-unpolarised environment (and a
singlet forces that in every basis), and ln 2 for exactly one unpaired electron
confined to the environment block.  The real-data anchor is the N2 CAS(6,6)
fixture, where the singlet spin-flip symmetry is checked to hold in both
available bases.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from fblockkit.analysis.entropy_rdm import (
    SpinDensities,
    rotate_densities,
    rotation_from_coefficients,
    solve_fci,
    spin_densities,
)
from fblockkit.analysis.environment_spin import (
    EnvironmentSpinError,
    analyze,
    environment_spin_entropy,
    evidence,
    loewdin_atom_populations,
    partition_by_centre,
    run,
)
from fblockkit.parsers.fcidump import parse_fcidump
from fblockkit.parsers.orca_json import parse_orca_json

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"
DUMP = FIXTURES / "n2_fcidump.fcidump"
CANONICAL = FIXTURES / "n2_fcidump.canonical.json"
LOCALIZED = FIXTURES / "n2_fcidump.localized.json"
ACTIVE = slice(4, 10)


def _n2_densities() -> SpinDensities:
    dump = parse_fcidump(DUMP)
    state = solve_fci(dump, reference_energy=-108.950671945279, multiplicity=1)
    return spin_densities(state)


def _one_electron(orbital: int, norb: int) -> SpinDensities:
    """A single alpha electron in ``orbital``; every two-body block is zero."""
    zeros = np.zeros((norb, norb, norb, norb))
    gamma_a = np.zeros((norb, norb))
    gamma_a[orbital, orbital] = 1.0
    return SpinDensities(
        gamma_a=gamma_a,
        gamma_b=np.zeros((norb, norb)),
        g_aa=zeros,
        g_ab=zeros.copy(),
        g_ba=zeros.copy(),
        g_bb=zeros.copy(),
    )


# --- the analytic anchors ---------------------------------------------------


def test_a_balanced_environment_gives_zero():
    d = np.diag([0.5, 0.3, 0.2])
    result = environment_spin_entropy(d, d)
    assert result.delta_s == pytest.approx(0.0, abs=1e-14)
    assert result.entropy_total == pytest.approx(result.entropy_alpha + result.entropy_beta)


def test_one_unpaired_electron_in_the_environment_gives_ln2():
    densities = _one_electron(1, norb=2)
    result = environment_spin_entropy(
        densities.gamma_a[1:, 1:], densities.gamma_b[1:, 1:]
    )
    assert result.delta_s == pytest.approx(math.log(2.0), abs=1e-12)
    assert result.electrons_alpha == pytest.approx(1.0)
    assert result.electrons_beta == pytest.approx(0.0)


def test_the_synthetic_two_orbital_value_matches_hand_computation():
    # D_a = diag(0.75, 0.25), D_b = 0 -> D = D_a.  With a fully polarised block the
    # identity Delta S_E = Tr(D) ln 2 holds exactly, so the hand value is ln 2.
    d_a = np.diag([0.75, 0.25])
    d_b = np.zeros((2, 2))
    expected = -2 * (0.375 * math.log(0.375) + 0.125 * math.log(0.125))
    expected -= -(0.75 * math.log(0.75) + 0.25 * math.log(0.25))
    result = environment_spin_entropy(d_a, d_b)
    assert result.delta_s == pytest.approx(expected, abs=1e-12)
    assert result.delta_s == pytest.approx(math.log(2.0), abs=1e-12)


def test_zero_eigenvalues_are_masked_not_logged():
    # an empty environment orbital contributes 0 * ln 0 = 0, never a NaN
    result = environment_spin_entropy(np.diag([1.0, 0.0]), np.diag([0.0, 0.0]))
    assert math.isfinite(result.delta_s)
    assert result.delta_s == pytest.approx(math.log(2.0), abs=1e-12)
    assert result.masked >= 1


# --- the real-data anchor ---------------------------------------------------


def test_a_singlet_state_has_no_polarisation_in_any_basis():
    densities = _n2_densities()
    assert np.abs(densities.gamma_a - densities.gamma_b).max() < 1e-12
    # and after a rotation, which mixes the orbitals
    localized = parse_orca_json(LOCALIZED)
    canonical = parse_orca_json(CANONICAL)
    rotation = rotation_from_coefficients(
        np.array(canonical.mo_coefficients).T,
        np.array(localized.mo_coefficients).T,
        np.array(canonical.overlap),
        active=range(*ACTIVE.indices(28)),
    )
    rotated = rotate_densities(densities, rotation)
    assert np.abs(rotated.gamma_a - rotated.gamma_b).max() < 1e-12
    window = np.arange(2)
    result = environment_spin_entropy(
        rotated.gamma_a[np.ix_(window, window)], rotated.gamma_b[np.ix_(window, window)]
    )
    assert result.delta_s == pytest.approx(0.0, abs=1e-12)


def test_the_atom_assignment_of_the_n2_active_orbitals_is_complete():
    export = parse_orca_json(LOCALIZED)
    coefficients = np.array(export.mo_coefficients)[ACTIVE].T
    populations = loewdin_atom_populations(
        coefficients, np.array(export.overlap), export.ao_labels, len(export.atoms)
    )
    for values in populations:
        assert sum(values) == pytest.approx(1.0, abs=1e-9)
        assert max(values) > 0.5  # every localised orbital belongs to one N
    centres, largest = partition_by_centre(
        coefficients,
        np.array(export.overlap),
        export.ao_labels,
        len(export.atoms),
        cluster_centres=(0,),
    )
    assert set(centres) <= {0, 1}
    assert all(value > 0.5 for value in largest)


def test_the_run_text_reports_the_partition_and_the_reading():
    densities = _one_electron(1, norb=2)
    result = environment_spin_entropy(
        densities.gamma_a[1:, 1:], densities.gamma_b[1:, 1:]
    )
    body = run(
        result,
        environment_orbitals=(1,),
        assignment=((0, 0, 1.0), (1, 1, 1.0)),
        cluster_centres=(0,),
    ).body
    assert "environment" in body and "cluster" in body
    assert "Delta S_E = 0.693147" in body
    assert "carries spin polarisation" in body
    assert "0.007" in body and "2.766" in body


def test_an_empty_environment_is_reported_as_by_construction_not_as_a_pass():
    body = run(
        None,
        environment_orbitals=(),
        assignment=((0, 0, 1.0),),
        cluster_centres=(0,),
    ).body
    assert "empty" in body and "by construction" in body


# --- input validation -------------------------------------------------------


def test_shape_mismatch_is_rejected():
    with pytest.raises(EnvironmentSpinError, match="same orbital set"):
        environment_spin_entropy(np.eye(2), np.eye(3))


def test_a_non_square_block_is_rejected():
    with pytest.raises(EnvironmentSpinError, match="not square"):
        environment_spin_entropy(np.ones((2, 3)), np.ones((2, 3)))


def test_an_empty_block_is_rejected():
    with pytest.raises(EnvironmentSpinError, match="empty"):
        environment_spin_entropy(np.zeros((0, 0)), np.zeros((0, 0)))


def test_a_negative_eigenvalue_is_rejected():
    with pytest.raises(EnvironmentSpinError, match="negative eigenvalue"):
        environment_spin_entropy(np.diag([1.0, -1.0]), np.zeros((2, 2)))


def test_evidence_carries_the_source_anchor():
    entries = evidence()
    text = " ".join(entry.ref + " " + entry.text for entry in entries)
    assert "10.1021/acs.jctc.5c01336" in text
    assert "2.766" in text and "0.007" in text


def test_analyze_on_the_n2_fixture_partitions_and_reads_zero_for_a_singlet():
    export = parse_orca_json(LOCALIZED)
    coefficients = np.array(export.mo_coefficients)[ACTIVE].T
    section = analyze(
        _n2_densities(),
        coefficients,
        np.array(export.overlap),
        export.ao_labels,
        len(export.atoms),
        cluster_centres=(0,),
    )
    assert "Delta S_E" in section.body
    # at least the "two centres" split shows up, and the singlet reading is zero
    assert "environment" in section.body
    if "electrons in the environment block" in section.body:
        assert "Delta S_E = 0.000000" in section.body
    else:
        assert "by construction" in section.body
