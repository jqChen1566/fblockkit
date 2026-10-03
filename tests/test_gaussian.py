"""Checks of the Gaussian minimal parser (parsers/gaussian.py; plan item 6.4).

Fixtures are two real G09 Rev D.01 outputs (the gau_orca example set, 2018;
see fixtures/gaussian/README): an H2CO frequency run (one imaginary mode at
-1089.0060 cm-1; the frequency block precedes the closing optimization
activity, so after_geometry is False) and an H2CO TS optimization
(opt(nomicro,calcfc,ts,noeigen); 7 steps; converged).  Both are
external-driver runs -- no SCF Done line, so the SCF facts are absent
(verified absent-field semantics).  Two G16 Rev C.01 probes (added
2026-10-02) supply the conventional-SCF anchor (the SCF Done regex is
fixture-verified) and the CASSCF facts (every marker measured; the active
natural occupations come from the final symbolic density matrix's
diagonal).  A third G16 probe (N2 CAS(6,6)/STO-3G, 2026-10-03) covers the
multi-block column pages of that matrix.
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
G16_RHF = GS / "g16_h2o_rhf.log"
G16_CAS = GS / "g16_h2o_cas22.log"
G16_CAS66 = GS / "n2_cas66.log"


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
    # the regex follows the published shape, now fixture-verified by the G16
    # conventional-SCF sample (see the function below)
    line = " SCF Done:  E(RB3LYP) =  -93.2545622     A.U. after   10 cycles"
    match = _SCF_DONE_RE.search(line)
    assert match is not None
    assert match.group(1) == "RB3LYP"
    assert float(match.group(2)) == approx(-93.2545622)
    assert int(match.group(3)) == 10


# --- the G16 conventional-SCF and CASSCF anchors (added 2026-10-02) -----------


def test_the_g16_rhf_anchor():
    """The conventional-SCF sample: the SCF Done regex is fixture-verified."""
    result = parse_auto(G16_RHF)
    assert result.sections["version"] == "16, Revision C.01"
    scf = result.sections["gaussian_scf"]
    assert scf["present"] is True and scf["converged"] is True
    assert scf["cycles"] == 7 and scf["method"] == "RHF"
    assert scf["energy"] == approx(-74.9630631539, abs=1e-10)
    facts = facts_from(result)
    assert facts["scf_converged"] is True
    assert facts["final_energy"] == approx(-74.9630631539, abs=1e-10)
    assert facts["casscf_present"] is False


def test_the_g16_casscf_anchor():
    result = parse_auto(G16_CAS)
    casscf = result.sections["gaussian_casscf"]
    assert casscf["present"] is True
    assert (casscf["active_orbitals"], casscf["active_electrons"]) == (2, 2)
    assert (casscf["core"], casscf["valence"], casscf["virtual"]) == (4, 2, 1)
    assert (casscf["iterations"], casscf["max_iterations"]) == (7, 64)
    assert casscf["energy"] == approx(-74.964316519, abs=1e-9)
    assert casscf["converged"] is True and casscf["converged_via"] is None
    assert casscf["active_occ_min"] == approx(0.00210880, abs=1e-8)
    assert casscf["active_occ_max"] == approx(1.99789, abs=1e-5)
    facts = facts_from(result)
    assert facts["casscf_present"] is True
    assert facts["casscf_converged"] is True
    # G16's marker states no criterion, so the ORCA-only energy-only rule
    # cannot fire on a Gaussian file (nothing guessed)
    assert facts["casscf_converged_via"] is None
    assert facts["final_energy"] == approx(-74.964316519, abs=1e-9)
    # the nearly-empty active orbital sits below the 0.02 documentation line
    assert facts["casscf_active_occ_min"] < 0.02
    # no conventional SCF Done line in this sample; the SCF facts stay absent
    assert result.sections["gaussian_scf"]["present"] is False


def test_the_g16_multiblock_symbolic_density_anchor():
    """The CAS(6,6) probe (2026-10-03) prints the lower triangle in
    five-column pages, so the sixth row's first-page last value is the (6,5)
    element -- a last-value diagonal read fabricates ~1e-18 for the sixth
    orbital.  The pinned occupations discriminate: right read min =
    0.0176416, max = 1.98270."""
    result = parse_auto(G16_CAS66)
    casscf = result.sections["gaussian_casscf"]
    assert casscf["present"] is True
    assert (casscf["active_orbitals"], casscf["active_electrons"]) == (6, 6)
    assert casscf["converged"] is True
    assert casscf["converged_via"] is None
    assert casscf["active_occ_max"] == approx(1.98270, abs=1e-5)
    assert casscf["active_occ_min"] == approx(0.0176416, abs=1e-7)
    facts = facts_from(result)
    assert facts["casscf_converged_via"] is None
    assert facts["casscf_active_occ_min"] < 0.02
