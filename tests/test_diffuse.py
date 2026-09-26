"""Tests for rule G, the diffuse-orbital (Rydberg) check (analysis/diffuse.py).

Fixture: ``n2_diffuse.out`` -- the N2/def2-SVP CAS(6,6) run repeated with both
prints the check needs (``!PrintBasis`` and the composition table).  The
expected exponents are the ones the run itself prints (1712.84 to 0.1876 for
the s channel, 13.571 to 0.2195 for p), and the expected ranking follows from
them.  ``eu_basis.out`` supplies the f-block face: its four f contractions.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from fblockkit.analysis import diffuse
from fblockkit.analysis.diffuse import DiffuseError, accepts, analyze, evidence, report
from fblockkit.parsers import parse_auto

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"
N2 = FIXTURES / "n2_diffuse.out"
EU = FIXTURES / "eu_basis.out"


def _sections(path):
    return parse_auto(path).sections


# --- the real-data anchor -----------------------------------------------------


def test_the_active_orbitals_are_ranked_by_channel_diffuseness():
    analysis = analyze(_sections(N2))
    assert analysis.active_indices == (4, 5, 6, 7, 8, 9)
    ranked = [(entry["index"], entry["shell"], entry["diffuse"]) for entry in analysis.entries]
    # the s-channel orbitals come first (0.1876 < 0.2195), and the order inside a
    # channel is stable
    assert [item[1] for item in ranked[:2]] == ["s", "s"]
    assert ranked[0][2] == pytest.approx(0.1876459225)
    assert ranked[2][2] == pytest.approx(0.2195434803)
    assert {item[2] for item in ranked} == {0.1876459225, 0.2195434803}


def test_the_channel_weights_and_the_span_are_reported():
    analysis = analyze(_sections(N2))
    by_index = {entry["index"]: entry for entry in analysis.entries}
    # orbital 4 is the sigma bonding orbital: 61% on the N s channel
    assert by_index[4]["weight"] == pytest.approx(61.2, abs=0.1)
    assert by_index[4]["tight"] == pytest.approx(1712.842, abs=0.001)
    # orbitals 5 and 6 are the pi pair: 97.8% on N p
    assert by_index[5]["weight"] == pytest.approx(97.8, abs=0.1)
    assert by_index[5]["tight"] == pytest.approx(13.571, abs=0.001)


def test_the_report_states_the_no_threshold_reading_and_the_remedy():
    body = report(analyze(_sections(N2))).body
    assert "ranking (most diffuse channel first)" in body
    assert "does not transfer" in body
    assert "SA-CASSCF + NEVPT2" in body
    assert "A8" in report(analyze(_sections(N2))).title


def test_the_f_block_fixture_prints_the_f_contractions():
    # the Eu SARC2 basis must show four f contractions, matching the four f shells
    # the export labels carry; the check skips the (absent) ECP blocks
    basis = _sections(EU)["basis"]
    f_shells = [
        min(shell["exponents"])
        for entry in basis["elements"]
        for shell in entry["shells"]
        if shell["angular"] == "f"
    ]
    assert [round(value, 3) for value in sorted(f_shells)] == [0.316, 0.786, 1.956, 4.870]


# --- the module contract -------------------------------------------------------


def test_accepts_only_with_both_prints():
    assert accepts(parse_auto(N2)) is True
    # the plain chain output has the composition table but no basis print
    assert accepts(parse_auto(FIXTURES / "n2_fcidump_step_a.out")) is False


def test_missing_prints_raise_with_a_next_step():
    with pytest.raises(DiffuseError, match="PrintBasis"):
        analyze(_sections(FIXTURES / "n2_fcidump_step_a.out"))
    with pytest.raises(DiffuseError, match="ReducedOrbPopMO_L"):
        analyze(_sections(FIXTURES / "n2_basis.out"))


def test_an_output_without_an_active_space_says_so():
    sections = _sections(N2)
    sections = {
        **sections,
        "orbitals": {**sections["orbitals"], "occupations": [2.0, 2.0] + [0.0] * 26},
    }
    body = report(analyze(sections)).body
    assert "no active orbital" in body


def test_evidence_carries_the_rule_and_the_measured_face():
    text = " ".join(entry.ref + " " + entry.text for entry in evidence())
    assert "2609.13357" in text
    assert "0.1876" in text or "0.1876459225" in text
    assert "Eu" in text
