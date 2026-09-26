"""Tests for the DM-AS selection reader (analysis/dm_selection.py, menu 20).

The real-data anchors are the H2O fixtures: six single-root CASCI candidates
on MP2 orbitals and a PBE0 SCF reference (2.0801 D); the selection picks
(6e, 8o) at 0.0081 D.  The state-averaged CASCI fixture pins the refusal, and
the synthetic readings pin the vector variant and the tie-break.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from fblockkit.analysis import dm_selection as dms
from fblockkit.parsers import parse_auto

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"
CANDIDATES = [("6", 6), ("6", 7), ("6", 8), ("8", 7), ("8", 8), ("10", 8)]


def _reading(nel: str, norb: int):
    path = FIXTURES / f"h2o_dm_casci_e{nel}o{norb}.out"
    result = parse_auto(path)
    magnitude, vector, _ = dms.candidate_dipole(result.sections, str(path))
    return (int(nel), norb, magnitude, vector)


def _reference():
    result = parse_auto(FIXTURES / "h2o_dm_ref_pbe0.out")
    magnitude, vector, _ = dms.reference_dipole(result.sections, "ref")
    return (magnitude, vector, "PBE0/def2-TZVP SCF")


def test_the_selection_picks_the_space_closest_to_the_reference():
    selection = dms.analyze([_reading(*c) for c in CANDIDATES], _reference())
    assert selection.selected == (6, 8)
    assert selection.deviation == pytest.approx(0.0081, abs=5e-4)
    # the ranking is by deviation, best first
    deviations = [row[0] for row in selection.results]
    assert deviations == sorted(deviations)
    by_space = {(row[1], row[2]): row[3] for row in selection.results}
    assert by_space[(6, 8)] == pytest.approx(2.0720, abs=5e-4)
    assert by_space[(6, 6)] == pytest.approx(2.1046, abs=5e-4)


def test_the_reference_reads_the_scf_block_of_an_mp2_output():
    """The prep output carries three blocks (SCF, MP2 x2): the SCF one is the
    KS-DFT-style reference the source used."""
    result = parse_auto(FIXTURES / "h2o_dm_prep_mp2.out")
    magnitude, _, method = dms.reference_dipole(result.sections, "prep")
    assert method == "SCF"
    assert magnitude == pytest.approx(2.16396, abs=5e-5)


def test_a_state_averaged_candidate_is_refused_with_the_rerun_instruction():
    result = parse_auto(FIXTURES / "h2o_dm_casci_e6o6_sa4.out")
    with pytest.raises(dms.DmSelectionError, match="nroots 1"):
        dms.candidate_dipole(result.sections, "sa4")


def test_the_directional_variant_projects_the_vectors():
    readings = [_reading(*c) for c in CANDIDATES]
    reference = _reference()
    selection = dms.analyze(readings, reference, protocol="vgdm")
    # by symmetry of this molecule the vectors are collinear, so the projected
    # deviation equals the magnitude deviation of the collinear component
    assert selection.protocol == "vgdm"
    assert selection.selected == (6, 8)
    with pytest.raises(dms.DmSelectionError, match="vector"):
        dms.analyze(readings, (2.08, None, "value only"), protocol="vgdm")


def test_the_unknown_and_unimplemented_protocols_are_refused_with_reasons():
    readings = [_reading(*c) for c in CANDIDATES]
    with pytest.raises(dms.DmSelectionError, match="per-state dipole"):
        dms.analyze(readings, _reference(), protocol="edm")
    with pytest.raises(dms.DmSelectionError, match="not known"):
        dms.analyze(readings, _reference(), protocol="something")


def test_the_tie_break_prefers_the_smaller_space():
    readings = [
        (6, 6, 2.00, (0.0, 0.0, 2.0)),
        (8, 8, 2.00, (0.0, 0.0, 2.0)),
        (12, 13, 2.00, (0.0, 0.0, 2.0)),
    ]
    selection = dms.analyze(readings, (2.0, (0.0, 0.0, 2.0), "test"))
    assert selection.selected == (6, 6)


def test_no_readings_is_refused():
    with pytest.raises(dms.DmSelectionError, match="no candidate readings"):
        dms.analyze([], _reference())


def test_the_report_states_the_reference_the_ranking_and_the_boundaries():
    selection = dms.analyze([_reading(*c) for c in CANDIDATES], _reference())
    body = dms.render(selection)
    assert "reference: 2.0801 D" in body
    assert "Selected active space: (6e, 8o)" in body
    assert "smallest deviation wins and ties go to the smaller space" in body
    assert "state-averaged" in body  # the refusal boundary is carried
    assert "menu 17" in body  # the scan-continuity pointer


def test_evidence_carries_the_two_papers_and_the_measured_blocks():
    text = " ".join(
        entry.ref + " " + entry.text + " " + entry.url + " " + entry.bibkey
        for entry in dms.evidence()
    )
    assert "2c01128" in text and "6c00473" in text
    assert "kaufold2023dipole" in text and "kaufold2026casci" in text
    assert "State: 0" in text and "2.0801" in text
