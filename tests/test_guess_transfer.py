"""Tests for the WASP guess transfer (recipe/guess_transfer.py).

The real-data anchors are the frozen N2 scan: the neighbours at 1.094 and
2.600 Angstrom, the template at 1.600 (its export for the metric and its mkl
for the write).  The distances, the 1/d weights and the Loewdin residual are
pinned exactly; the guess quality is checked only where it is well defined
(the core orbitals, which keep their identity) -- the degenerate bond triad's
per-index overlap is not, because a rotation inside the window is arbitrary.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from fblockkit.parsers.mkl import parse_mkl
from fblockkit.parsers.orca_json import parse_orca_json
from fblockkit.recipe import guess_transfer as gt

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"


def _export(tag: str):
    return parse_orca_json(FIXTURES / f"n2_scan_{tag}.loc.json")


def _result(neighbours=("1.094", "2.600"), template="1.600", **kwargs):
    return gt.interpolate_guess(
        [_export(tag) for tag in neighbours], _export(template), **kwargs
    )


def test_the_weights_follow_the_inverse_distances():
    result = _result()
    by_name = {weight.name: weight for weight in result.weights}
    # RMSD of the two-atom scan: one atom moves by 0.506 / 1.000 Angstrom
    assert by_name["n2_scan_1.094.loc"].distance == pytest.approx(0.357796, abs=1e-5)
    assert by_name["n2_scan_2.600.loc"].distance == pytest.approx(0.707107, abs=1e-5)
    assert by_name["n2_scan_1.094.loc"].weight == pytest.approx(0.664011, abs=1e-5)
    assert by_name["n2_scan_2.600.loc"].weight == pytest.approx(0.335989, abs=1e-5)
    assert result.source == "n2_scan_1.094.loc"


def test_the_mixture_is_orthonormal_in_the_template_metric():
    result = _result()
    template = _export("1.600")
    overlap = np.array(template.overlap)
    coefficients = np.array(result.coefficients)
    identity = coefficients.T @ overlap @ coefficients
    assert np.abs(identity - np.eye(28)).max() < 1e-10
    assert result.residual < 1e-10


def test_the_core_orbitals_transfer_with_high_overlap():
    """The cores keep their identity, so their overlap is a fair quality anchor."""
    result = _result()
    template = _export("1.600")
    overlap = np.array(template.overlap)
    guess = np.array(result.coefficients)
    reference = np.array(template.mo_coefficients).T
    block = guess.T @ overlap @ reference
    assert abs(block[0, 0]) > 0.98
    assert abs(block[2, 2]) > 0.98


def test_delta_restricts_the_neighbourhood():
    result = _result(delta=0.45)
    by_name = {weight.name: weight.weight for weight in result.weights}
    assert by_name["n2_scan_1.094.loc"] == pytest.approx(1.0)
    assert by_name["n2_scan_2.600.loc"] == 0.0
    with pytest.raises(gt.GuessError, match="within delta"):
        _result(delta=0.1)


def test_a_coincident_neighbour_takes_the_whole_weight():
    result = _result(neighbours=("1.600", "1.094"))
    by_name = {weight.name: weight for weight in result.weights}
    assert by_name["n2_scan_1.600.loc"].weight == 1.0
    assert by_name["n2_scan_1.094.loc"].weight == 0.0
    assert any("coincides" in note for note in result.notes)


def test_mismatched_atoms_are_refused():
    neighbour = _export("1.094")
    broken = replace(neighbour, atoms=("N", "O"))
    with pytest.raises(gt.GuessError, match="atom list"):
        gt.interpolate_guess([broken], _export("1.600"))


def test_an_export_without_overlap_is_refused():
    stripped = replace(_export("1.094"), overlap=None)
    with pytest.raises(gt.GuessError, match="S-Matrix"):
        gt.interpolate_guess([stripped], _export("1.600"))


def test_the_written_mkl_carries_the_guess(tmp_path):
    result = _result()
    template_mkl = parse_mkl(FIXTURES / "n2_scan_1.600.mkl")
    out = tmp_path / "target.fbk.mkl"
    gt.write_guess_mkl(template_mkl, result, out)
    again = parse_mkl(out)
    assert np.abs(np.array(again.coefficients()) - np.array(result.coefficients)).max() < 1e-7
    assert again.atoms == template_mkl.atoms
    assert again.occupations == pytest.approx(result.occupations)


def test_the_report_and_evidence_state_the_protocol_and_the_route():
    result = _result()
    body = gt.render(result)
    assert "WASP initial guess" in body
    assert "1/d" in body
    assert "orca_2mkl" in body and "moread" in body
    text = " ".join(
        entry.ref + " " + entry.text + " " + entry.url + " " + entry.bibkey
        for entry in gt.evidence()
    )
    assert "chemrev.5c00866" in text
    assert "wardzala2026multireference" in text
    assert "MOREAD" in text
