"""Checks of the Gaussian minimal parser (parsers/gaussian.py; plan item 6.4).

Fixtures are two real G09 Rev D.01 outputs (the gau_orca example set, 2018;
see fixtures/gaussian/README): an H2CO frequency run (one imaginary mode at
-1089.0060 cm-1; the frequency block precedes the closing optimization
activity, so after_geometry is False) and an H2CO TS optimization
(opt(nomicro,calcfc,ts,noeigen); 7 steps; converged).  Both are
external-driver runs -- no SCF Done line, so the SCF facts are absent
(verified absent-field semantics).  The SCF Done regex follows the
published line shape and is exercised here on that shape only; its real
fixture anchor is registered as a follow-up.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from pytest import approx

from fblockkit.parsers import detect_program, facts_from, parse_auto
from fblockkit.parsers.base import ParserError
from fblockkit.parsers.gaussian import GaussianParser, _SCF_DONE_RE
from fblockkit.parsers.orca import OrcaParser

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
GS = FIXTURES / "gaussian"
FREQ = GS / "g09_h2co_freq.out"
TS = GS / "g09_h2co_ts.out"


@pytest.fixture(scope="module")
def freq():
    return parse_auto(FREQ)


@pytest.fixture(scope="module")
def ts():
    return parse_auto(TS)


# --- identification ------------------------------------------------------------


def test_identification_and_negative_control(freq, ts):
    assert freq.program == "gaussian"
    assert ts.program == "gaussian"
    # an ORCA output is not detected as Gaussian, and vice versa
    orca_fixture = FIXTURES / "orca" / "n2_casscf_nevpt2.out"
    assert GaussianParser().detect(orca_fixture) is False
    assert OrcaParser().detect(FREQ) is False
    assert detect_program(FREQ).program == "gaussian"


def test_the_parser_refuses_foreign_files():
    with pytest.raises(ParserError, match="not Gaussian"):
        GaussianParser().parse(FIXTURES / "orca" / "n2_casscf_nevpt2.out")


# --- the facts ------------------------------------------------------------------


def test_the_frequency_probe(freq):
    sections = freq.sections
    assert sections["version"] == "09, Revision D.01"
    assert sections["termination"]["normal"] is True
    assert sections["termination"]["error_terminations"] == 0
    # an external-driver run: no SCF Done line, facts stay absent
    assert sections["gaussian_scf"]["present"] is False
    frequency = sections["gaussian_frequencies"]
    assert frequency["present"] is True
    block = frequency["blocks"][0]
    assert block["n_imaginary"] == 1
    assert block["engine_imaginary_marks"] == [1]
    assert min(block["frequencies"]) == approx(-1089.0060)
    # the frequency block precedes the closing optimization activity
    assert frequency["after_geometry"] is False
    assert sections["gaussian_optimization"]["converged"] is True
    facts = facts_from(freq)
    assert facts["terminated_normally"] is True
    assert facts["frequency_imaginary_count"] == 1
    assert facts["frequency_min_imaginary"] == approx(-1089.0060)
    assert facts["ts_optimization"] is False
    assert "scf_converged" not in facts or facts["scf_converged"] is None


def test_the_ts_probe(ts):
    sections = ts.sections
    optimization = sections["gaussian_optimization"]
    assert optimization["present"] is True
    assert optimization["converged"] is True
    assert optimization["ts"] is True  # the input echo carries opt(...,ts,...)
    assert optimization["steps"] == 7
    assert sections["gaussian_frequencies"]["present"] is False
    assert sections["gaussian_frequencies"]["after_geometry"] is None
    facts = facts_from(ts)
    assert facts["ts_optimization"] is True
    assert facts["optimization_converged"] is True
    assert facts["frequency_present"] is False
    assert facts["terminated_normally"] is True


def test_the_scf_done_line_shape():
    # the regex follows the published shape; a conventional-SCF fixture is
    # registered as a follow-up (the shipped probes are external-driver runs)
    line = " SCF Done:  E(RB3LYP) =  -93.2545622     A.U. after   10 cycles"
    match = _SCF_DONE_RE.search(line)
    assert match is not None
    assert match.group(1) == "RB3LYP"
    assert float(match.group(2)) == approx(-93.2545622)
    assert int(match.group(3)) == 10
