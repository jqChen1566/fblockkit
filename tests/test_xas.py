"""Checks of the core-excited-spectra machinery (parsers/rocis_spectra.py +
analysis/xas.py; menu 37; Wave 5.4).

The fixture is a real ORCA 6.1.1 ROCIS run of [FeCl4]2- (x2c-SVPall, ROHF
high-spin d6, NRoots 30, DecomposeFosc, DoSOC): fourteen absorption blocks
(seven plain and seven SOC-corrected; the SOC tables weight fosc by the
initial-state population), 90 excitation-table rows (three spin blocks).

Anchors: the SOC-corrected electric-dipole block holds 2235 state pairs of
which 935 carry fosc > 1e-6; its first non-zero transition is 4-5A -> 5-5A at
715.0726 eV with fosc 0.000892; the largest-gap split lands at 735.602 eV
with clusters [663, 272] and a ratio of 56.355.  The RIXS bookkeeping is
covered by synthetic dictionaries (the measured refusal path of the fixture
run and the intermediate/final-state count of a successful run).
"""

from __future__ import annotations

from pathlib import Path

import pytest
from pytest import approx

from fblockkit.analysis import xas
from fblockkit.diagnosis import references_section
from fblockkit.parsers import parse_auto

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
ROCI = FIXTURES / "rocis"


@pytest.fixture(scope="module")
def parsed():
    return parse_auto(ROCI / "fecl4_xas.out").sections["rocis"]


# --- the parse -----------------------------------------------------------------


def test_the_rocis_blocks_parse_with_both_phases(parsed):
    assert parsed["present"]
    blocks = parsed["spectra"]
    expected = {
        "dipole_length",
        "dipole_velocity",
        "combined_length",
        "combined_velocity",
        "combined_origin_adjusted",
        "combined_origin_independent",
        "combined_origin_independent_velocity",
    }
    assert expected <= set(blocks)
    assert {f"soc_{key}" for key in expected} <= set(blocks)
    assert len(blocks) == 14


def test_the_excitation_table_and_a_first_row(parsed):
    table = parsed["excitation_table"]
    assert len(table) == 90
    assert {row["mult"] for row in table} == {3, 5, 7}
    row = parsed["spectra"]["dipole_length"][0]
    assert (row["i_root"], row["j_root"]) == (0, 1)
    assert row["i_label"] == "5A"
    assert row["ev"] == approx(718.864517)
    assert row["fosc"] == approx(0.00796713, abs=1e-8)
    assert row["cm1"] == approx(5798034.0)


def test_the_soc_block_carries_the_full_state_pairs(parsed):
    rows = parsed["spectra"]["soc_dipole_length"]
    assert len(rows) == 2235
    nonzero = [row for row in rows if row["fosc"] > 1e-6]
    assert len(nonzero) == 935
    first = nonzero[0]
    assert (first["i_root"], first["j_root"]) == (4, 5)
    assert first["ev"] == approx(715.0726)
    assert first["fosc"] == approx(0.000892271, abs=1e-9)


# --- the analysis --------------------------------------------------------------


def test_the_primary_block_prefers_the_soc_corrected_table(parsed):
    key, rows = xas.primary_spectrum(parsed)
    assert key == "soc_dipole_length"
    assert len(rows) == 2235


def test_the_branching_splits_at_the_largest_gap(parsed):
    _, rows = xas.primary_spectrum(parsed)
    nonzero = [row for row in rows if row["fosc"] > 1e-6]
    edge = xas.branching(nonzero)
    assert edge["split_ev"] == approx(735.6025, abs=1e-3)
    low, high = edge["clusters"]
    assert low["n"] == 663 and high["n"] == 272
    assert low["centroid_ev"] == approx(720.93, abs=0.05)
    assert high["centroid_ev"] == approx(742.24, abs=0.05)
    assert edge["ratio_low_high"] == approx(56.355, rel=1e-3)


def test_branching_stays_single_cluster_when_no_gap_reaches_the_threshold():
    rows = [
        {"ev": 700.0 + 0.1 * index, "fosc": 0.01} for index in range(10)
    ]
    edge = xas.branching(rows)
    assert edge["split_ev"] is None
    assert len(edge["clusters"]) == 1 and edge["clusters"][0]["n"] == 10


def test_primary_spectrum_refuses_an_empty_section():
    with pytest.raises(xas.XasError, match="ROCIS"):
        xas.primary_spectrum({"present": False, "spectra": {}})


# --- the render ----------------------------------------------------------------


def test_the_render_reports_the_table_split_and_rixs_status(parsed):
    body = xas.render(parsed, source="fecl4_xas.out", r_stat=2.0)
    assert "primary block: soc_dipole_length" in body
    assert "935 with non-zero fosc" in body
    assert "branching at 735.60 eV" in body
    assert "ratio (low/high) = 56.355" in body
    assert "ratio/stat = 28.177" in body  # 56.355 / 2
    assert "Thole & van der Laan" in body
    # this fixture ran with the RIXS flags off
    assert "RIXS: not requested" in body


def test_the_render_states_the_measured_refusal_and_the_mapspc_recipe():
    data = {
        "present": True,
        "excitation_table": [],
        "spectra": {"dipole_length": [{"i_root": 0, "i_label": "5A", "j_root": 1,
                                       "j_label": "5A", "ev": 718.0, "cm1": 1.0,
                                       "nm": 1.0, "fosc": 0.1, "rest": []}]},
        "riqs_refused": True,
    }
    body = xas.render(data, source="x.out")
    assert "zero intermediate/final states" in body
    assert "6-element window" in body
    data2 = {
        "present": True,
        "excitation_table": [],
        "spectra": data["spectra"],
        "riqs_intermediate": 78,
        "riqs_final": 214,
    }
    body2 = xas.render(data2, source="x.out")
    assert "intermediate states 78" in body2 and "final states 214" in body2
    assert "orca_mapspc" in body2


def test_the_rixs_refusal_fixture_carries_the_engine_warning():
    parsed = parse_auto(ROCI / "fecl4_xas_rixs.out").sections["rocis"]
    assert parsed["present"] and parsed.get("riqs_refused") is True
    body = xas.render(parsed, source="fecl4_xas_rixs.out")
    assert "zero intermediate/final states" in body
    assert "6-element window" in body


def test_the_evidence_is_citable():
    evidence = xas.evidence()
    assert references_section(evidence) is not None
    literature = [item for item in evidence if item.kind == "literature"]
    assert literature and "10.1103/PhysRevB.38.3158" in literature[0].ref
