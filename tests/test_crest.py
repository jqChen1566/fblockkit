"""Checks of the CREST ensemble reader and the ensemble report (menu 43; plan item 6.1).

Fixtures are a real CREST 3.0.2 quick search (101, 2026-09-30, n-butane on
GFN2; see fixtures/crest/README): crest_conformers.xyz (two conformers,
leading-space atom-count lines, hartree frame comments), crest.energies
(relative kcal/mol) and crest_best.xyz.  Every format convention in the
parser is measured on this run.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from pytest import approx

from fblockkit.parsers.crest import CrestError, read_crest_directory, read_crest_ensemble
from fblockkit.analysis.crest import boltzmann_weights, render

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
CR = FIXTURES / "crest"


def test_the_ensemble_reads():
    ensemble = read_crest_directory(CR)
    assert len(ensemble.frames) == 2
    assert len(ensemble.frames[0].atoms) == 14
    assert ensemble.frames[0].energy_Eh == approx(-13.66512758)
    assert ensemble.frames[1].energy_Eh == approx(-13.66417737)
    assert ensemble.relative_kcal == approx((0.0, 0.596))
    # crest_best.xyz points at the first (lowest) frame
    assert ensemble.best_index == 1
    # the two sources agree within the three-decimal print precision
    assert ensemble.max_rel_deviation_kcal < 0.001


def test_the_directory_can_be_given_as_the_xyz(tmp_path):
    ensemble = read_crest_directory(CR / "crest_conformers.xyz")
    assert len(ensemble.frames) == 2


def test_a_missing_ensemble_reports_next_step(tmp_path):
    with pytest.raises(CrestError, match="crest_conformers.xyz"):
        read_crest_directory(tmp_path)


def test_the_boltzmann_weights():
    ensemble = read_crest_directory(CR)
    weights = boltzmann_weights(ensemble, temperature_K=298.15)
    assert weights[0] + weights[1] == approx(1.0)
    assert weights[0] == approx(0.7322, abs=5e-4)
    # at a very low temperature the ensemble collapses onto the best frame
    cold = boltzmann_weights(ensemble, temperature_K=1.0)
    assert cold[0] == approx(1.0)


def test_the_report_quotes_the_table_and_the_check():
    body = render(read_crest_directory(CR), source="fixtures/crest")
    assert "Conformers: 2" in body
    assert "index 1" in body
    assert "0.596" in body
    assert "Cross-check" in body
    assert "GFN2" in body


def test_a_frame_without_a_numeric_comment_is_refused():
    broken = "2\nnot an energy\nH 0.0 0.0 0.0\nH 0.0 0.0 1.0\n"
    with pytest.raises(CrestError, match="hartree"):
        read_crest_ensemble(broken)
