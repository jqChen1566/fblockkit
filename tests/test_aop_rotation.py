"""Tests for the AOP rotation guess (recipe/aop_rotation.py).

The real-data anchors are the frozen N2 fixtures: the reference is the
converged CASSCF(6,6)/def2-SVP export at 1.600 (n2_cas666_1.600.json), the
target the RHF set of the same geometry (n2_scan_1.600.json + its mkl).  The
engine anchors (the CASSCF validation runs) live in n2_aop_1.600.out /
n2_cas666_1.600.out and are quoted in the fixture README; here the structure
is pinned: the same-geometry rotation reproduces the reference active block
exactly (O_min 1.000000) at containment 0.958 and unitarity at the numerical
floor, the identity case reproduces the reference up to the degenerate
window, and the non-corresponding cross-geometry pair (n2_cas666_1.094.json
at containment 0.022) is refused by the containment gate -- the gate exists
because that same input, before it existed, sent the CASSCF to a wrong
solution (n2_aop_1.600_mismatch.out).
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from fblockkit.parsers.mkl import parse_mkl
from fblockkit.parsers.orca_json import parse_orca_json
from fblockkit.recipe import aop_rotation as aop

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"


def _scan(tag: str):
    return parse_orca_json(FIXTURES / f"n2_scan_{tag}.loc.json")


def _reference():
    return parse_orca_json(FIXTURES / "n2_cas666_1.600.json")


def _target():
    return parse_orca_json(FIXTURES / "n2_scan_1.600.json")


def _rotate_same():
    return aop.rotate_guess(
        _reference(), _target(), [4, 5, 6, 7, 8, 9], n_closed=4
    )


def _achieved(result, reference_export, target_export, reference_rows, target_rows):
    """The achieved overlap block, computed the module's own way: the reference
    coefficients against the guess, in the target's overlap metric."""
    reference = np.array(reference_export.mo_coefficients).T
    coefficients = np.array(result.coefficients)  # n_ao x n_mo
    overlap = np.array(target_export.overlap)
    return reference[:, reference_rows].T @ overlap @ coefficients[:, target_rows]


def test_the_same_geometry_rotation_reproduces_the_reference_exactly():
    """The paper's own setting (a reference active space and the SCF set of the
    same geometry): the built active block reproduces the reference (O_min
    exactly 1), the containment reads the pinned 0.958, and the written set is
    unitary in the target's metric."""
    result = _rotate_same()
    assert result.omin == pytest.approx(1.0, abs=1e-9)
    assert result.containment == pytest.approx(0.958, abs=1e-3)
    assert result.residual < 1e-12
    assert result.n_closed == 4 and result.n_active == 6 and result.n_virtual == 18
    achieved = _achieved(result, _reference(), _target(), [4, 5, 6, 7, 8, 9], [4, 5, 6, 7, 8, 9])
    # the achieved overlap matrix is square-orthogonal (the block IS the reference)
    assert float(np.abs(achieved @ achieved.T - np.eye(6)).max()) < 1e-9


def test_the_identity_rotation_reproduces_the_reference():
    """Reference == target on the localized scan exports: the active block is
    reproduced (up to the arbitrary rotation inside the degenerate window) --
    O_min reads exactly 1 and the written set stays orthonormal."""
    result = aop.rotate_guess(_scan("1.600"), _scan("1.600"), [4, 5, 6], n_closed=4)
    assert result.omin == pytest.approx(1.0, abs=1e-9)
    assert result.residual < 1e-12
    achieved = _achieved(result, _scan("1.600"), _scan("1.600"), [4, 5, 6], [4, 5, 6])
    assert float(np.abs(achieved @ achieved.T - np.eye(3)).max()) < 1e-9


def test_a_non_corresponding_cross_geometry_reference_is_refused():
    """The measured failure behind the containment gate: the CASSCF(6,6)
    reference at 1.094 has active character outside the 1.600 window
    (containment 0.022), which would scramble the closed block -- refused."""
    reference = parse_orca_json(FIXTURES / "n2_cas666_1.094.json")
    with pytest.raises(aop.AopError, match="containment"):
        aop.rotate_guess(reference, _target(), [4, 5, 6, 7, 8, 9], n_closed=4)


def test_the_refusals():
    with pytest.raises(aop.AopError, match="S-Matrix"):
        aop.rotate_guess(
            _reference(), replace(_target(), overlap=None), [4, 5, 6, 7, 8, 9],
            n_closed=4,
        )
    with pytest.raises(aop.AopError, match="basis set and atom order"):
        aop.rotate_guess(
            _reference(), replace(_target(), n_ao=99), [4, 5, 6, 7, 8, 9], n_closed=4
        )
    with pytest.raises(aop.AopError, match="active list is empty"):
        aop.rotate_guess(_reference(), _target(), [], n_closed=4)
    with pytest.raises(aop.AopError, match="outside"):
        aop.rotate_guess(_reference(), _target(), [99], n_closed=4)
    with pytest.raises(aop.AopError, match="closed count"):
        aop.rotate_guess(_reference(), _target(), [4, 5, 6, 7, 8, 9], n_closed=99)
    with pytest.raises(aop.AopError, match="exceeds"):
        aop.rotate_guess(_reference(), _target(), [4, 5, 6, 7, 8, 9], n_closed=25)
    with pytest.raises(aop.AopError, match="carries no MO coefficients"):
        aop.rotate_guess(
            _reference(), replace(_target(), mo_coefficients=()), [4, 5, 6, 7, 8, 9],
            n_closed=4,
        )


def test_the_written_mkl_carries_the_rotation(tmp_path):
    result = _rotate_same()
    template_mkl = parse_mkl(FIXTURES / "n2_scan_1.600.mkl")
    out = tmp_path / "target.aop.fbk.mkl"
    aop.write_guess_mkl(template_mkl, result, out)
    again = parse_mkl(out)
    assert np.abs(np.array(again.coefficients()) - np.array(result.coefficients)).max() < 1e-7
    assert again.atoms == template_mkl.atoms
    assert again.occupations == pytest.approx(result.occupations)


def test_the_engine_validation_anchors_travel_in_the_fixtures():
    """The frozen CASSCF validation outputs carry the measured anchors: the
    aufbau start takes 7 macro-iterations, the AOP start 6, to the same
    solution; the mismatch run wanders to a wrong solution with the
    delocalized-closed warning."""
    aufbau = (FIXTURES / "n2_cas666_1.600.out").read_text(encoding="utf-8")
    aop_run = (FIXTURES / "n2_aop_1.600.out").read_text(encoding="utf-8")
    mismatch = (FIXTURES / "n2_aop_1.600_mismatch.out").read_text(encoding="utf-8")
    assert aufbau.count("MACRO-ITERATION") == 7
    assert "-108.772368734" in aufbau
    assert aop_run.count("MACRO-ITERATION") == 6
    assert "-108.772368704" in aop_run
    assert "delocalized" in mismatch
    assert "-77.092699913" in mismatch


def test_the_report_and_evidence_state_the_protocol_and_the_route():
    body = aop.render(_rotate_same())
    assert "AOP rotation guess" in body
    assert "Eqs. (4)-(11)" in body
    assert "containment" in body and "0.958" in body
    assert "O_min" in body and "1.000000" in body
    assert "orca_2mkl" in body and "moread" in body
    assert "orthonormality residual" in body
    text = " ".join(
        entry.ref + " " + entry.text + " " + entry.url + " " + entry.bibkey
        for entry in aop.evidence()
    )
    assert "5.0058673" in text
    assert "paz2021active" in text
    assert "moread" in text
    assert "0.022" in text
