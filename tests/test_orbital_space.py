"""Tests for the orbital-space comparison (analysis/orbital_space.py).

The real-data side uses the two N2 exports of the same CAS(6,6) orbital file
(``n2_fcidump.canonical.json`` / ``...localized.json``): the same six active
orbitals before and after an orca_loc IAO-IBO localisation.  The identity limit
-- an orbital rotation inside one space cannot change the space -- is what the
fixture pair measures, and it is the property both diagnostics must have.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from fblockkit.analysis.orbital_space import (
    OrbitalSpaceError,
    compare_spaces,
    evidence,
    run,
)
from fblockkit.parsers.orca_json import parse_orca_json

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"
CANONICAL = FIXTURES / "n2_fcidump.canonical.json"
LOCALIZED = FIXTURES / "n2_fcidump.localized.json"

# orca_loc localized this window only, and it is the CAS(6,6) active space
ACTIVE = slice(4, 10)


def _blocks():
    canonical = parse_orca_json(CANONICAL)
    localized = parse_orca_json(LOCALIZED)
    overlap = np.array(canonical.overlap)
    a = np.array(canonical.mo_coefficients)[ACTIVE].T
    b = np.array(localized.mo_coefficients)[ACTIVE].T
    return a, b, overlap


# --- the identity limit on real data ----------------------------------------


def test_a_space_against_itself_is_perfectly_contained():
    a, _, s = _blocks()
    comparison = compare_spaces(a, a, s)
    assert comparison.n_a == comparison.n_b == 6
    assert len(comparison.singular_values) == 6
    assert max(abs(value - 1.0) for value in comparison.singular_values) < 1e-10
    assert comparison.fraction == pytest.approx(1.0, abs=1e-10)
    assert comparison.deficit == pytest.approx(0.0, abs=1e-10)
    assert comparison.reorthonormalised_a is False
    assert comparison.reorthonormalised_b is False


def test_an_orbital_rotation_inside_the_space_is_invisible():
    # canonical vs localized: same space, different orbitals -- both readings
    # must be blind to the rotation (sigma = 1, "essentially the same space")
    a, b, s = _blocks()
    comparison = compare_spaces(a, b, s)
    assert max(abs(value - 1.0) for value in comparison.singular_values) < 1e-10
    body = run(comparison, label_a="canonical", label_b="localized").body
    assert "essentially the same space" in body
    assert "sigma_F = ||M||_F / sqrt(min(n_A, n_B)) = 1.000000" in body


# --- containment and deficiency ---------------------------------------------


def test_a_subspace_of_a_larger_space_is_contained():
    canonical = parse_orca_json(CANONICAL)
    overlap = np.array(canonical.overlap)
    coefficients = np.array(canonical.mo_coefficients)
    larger = coefficients[4:12].T  # eight orbitals
    smaller = coefficients[4:10].T  # six of them
    comparison = compare_spaces(larger, smaller, overlap)
    assert comparison.n_a == 8 and comparison.n_b == 6
    assert len(comparison.singular_values) == 6
    assert max(abs(value - 1.0) for value in comparison.singular_values) < 1e-10
    assert comparison.fraction == pytest.approx(1.0, abs=1e-10)


def test_a_shifted_window_declares_the_missing_directions():
    canonical = parse_orca_json(CANONICAL)
    overlap = np.array(canonical.overlap)
    coefficients = np.array(canonical.mo_coefficients)
    a = coefficients[4:10].T  # window 4..9
    b = coefficients[5:11].T  # window 5..10: five orbitals in common, one each way
    comparison = compare_spaces(a, b, overlap)
    assert 0.0 <= comparison.fraction < 1.0
    assert 0.0 <= comparison.min_singular < 1.0
    body = run(comparison).body
    assert "misses a substantial part" in body or "contained only approximately" in body
    # the five shared directions stay perfect; only the exchanged pair moves
    assert sum(1 for value in comparison.singular_values if value > 1 - 1e-10) == 5


# --- the orthonormalisation step --------------------------------------------


def test_a_non_orthonormal_block_is_repaired_and_gives_the_same_reading():
    a, b, s = _blocks()
    scaled = a.copy()
    scaled[:, 2] *= 2.5  # still spans the same space, no longer orthonormal
    comparison = compare_spaces(scaled, b, s)
    assert comparison.reorthonormalised_a is True
    assert max(abs(value - 1.0) for value in comparison.singular_values) < 1e-10
    assert comparison.fraction == pytest.approx(1.0, abs=1e-10)
    assert "was not orthonormal in the AO metric" in run(comparison).body


def test_a_linearly_dependent_block_is_rejected():
    a, _, s = _blocks()
    dependent = a.copy()
    dependent[:, 1] = dependent[:, 0]
    with pytest.raises(OrbitalSpaceError, match="linearly dependent"):
        compare_spaces(dependent, a, s)


# --- input validation -------------------------------------------------------


def test_a_dimension_mismatch_is_rejected():
    a, _, s = _blocks()
    with pytest.raises(OrbitalSpaceError, match="must come from the same basis set"):
        compare_spaces(a[:-1], a, s)


def test_a_non_square_overlap_is_rejected():
    a, _, _ = _blocks()
    with pytest.raises(OrbitalSpaceError, match="not square"):
        compare_spaces(a, a, np.eye(3, 4))


def test_an_empty_block_is_rejected():
    a, _, s = _blocks()
    with pytest.raises(OrbitalSpaceError, match="no orbitals"):
        compare_spaces(a[:, :0], a, s)


# --- context ------------------------------------------------------------------


def test_unequal_sizes_skip_the_space_change_reading():
    canonical = parse_orca_json(CANONICAL)
    overlap = np.array(canonical.overlap)
    coefficients = np.array(canonical.mo_coefficients)
    comparison = compare_spaces(coefficients[4:12].T, coefficients[4:10].T, overlap)
    body = run(comparison).body
    assert "space change: not read here" in body
    assert "containment:" in body


def test_evidence_carries_the_literature_anchors():
    entries = evidence()
    assert len(entries) == 3
    text = " ".join(entry.ref + " " + entry.text for entry in entries)
    assert "10.1021/acs.jctc.7b00128" in text  # AVAS paper
    assert "2607.08178" in text  # SA-DMET preprint
    assert any("2.9e-3" in entry.text for entry in entries)
