"""Checks of the hyperfine / EFG machinery (parsers/hyperfine.py +
analysis/hyperfine.py; menu 40; Wave 5.7).

Fixtures are real ORCA 6.1.1 EPRNMR runs on CeF3 (2026-09-29): a DFT
(PBE0/x2c-SVPall) probe without and with nuclear parameters, and a
CASSCF+SOC probe.  Anchors: the DFT probe's EFG principal values
V(Tot) = 0.6771983 / 0.7337722 / -1.4109705 a.u. with eta = 0.040096
(derived; the parameter run prints the engine value 0.040223), V(El)+V(Nuc)
reproducing V(Tot) to the printed digits, RHO(0) = 727216.913536801
(all-electron) vs 0.075618828 (small-core ECP, CASSCF probe); the nuclear
V(Nuc) triplet is the same on both probes; the parameter run prints
A(FC) = -11.8696 MHz and A(iso) = -11.8696 MHz; the CASSCF probe switches
the A components off and reports method "CASSCF/ALL STATES AVERAGE".
"""

from __future__ import annotations

from pathlib import Path

import pytest
from pytest import approx

from fblockkit.analysis import hyperfine as hyperfine_analysis
from fblockkit.diagnosis import references_section
from fblockkit.parsers import parse_auto

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
HYP = FIXTURES / "hyperfine"
DFT = HYP / "cef3_epr_dft.out"
NUC = HYP / "cef3_epr_nuc.out"
NUCQ = HYP / "cef3_epr_nuc_q.out"
CASSCF = HYP / "cef3_epr_casscf.out"


@pytest.fixture(scope="module")
def dft():
    return parse_auto(DFT).sections["hyperfine"]


@pytest.fixture(scope="module")
def nuc():
    return parse_auto(NUC).sections["hyperfine"]


@pytest.fixture(scope="module")
def casscf():
    return parse_auto(CASSCF).sections["hyperfine"]


# --- the parse -----------------------------------------------------------------


def test_the_section_and_header_parse(dft):
    assert dft["present"]
    header = dft["header"]
    assert header["method"] == "SCF"
    assert header["multiplicity"] == 2
    assert header["A(iso)"] == 1
    assert header["A(dip)"] == 1
    assert header["A(orb)"] == 0
    assert header["EFG"] == 1
    assert header["Rho(0)"] == 1


def test_the_dft_probe_efg_and_rho(dft):
    (nucleus,) = dft["nuclei"]
    assert nucleus["element"] == "Ce" and nucleus["index"] == 0
    assert nucleus["V_Tot"] == approx([0.6771983, 0.7337722, -1.4109705])
    assert nucleus["V_el"] == approx([0.4892665, 0.5458783, -1.0351447])
    assert nucleus["V_nuc"] == approx([0.1879318, 0.187894, -0.3758258])
    assert nucleus["rho0"] == approx(727216.913536801)
    assert sorted(nucleus["efg_orientation"]) == ["X", "Y", "Z"]
    # the A matrix is printed as zeros without nuclear parameters
    assert nucleus["A_Tot"] == approx([0.0, 0.0, 0.0])
    assert (nucleus["a"] or {}).get("P_MHz_au3") == approx(0.0)


def test_the_efg_metrics_and_sum_check(dft):
    (nucleus,) = dft["nuclei"]
    metrics = hyperfine_analysis.efg_metrics(nucleus)
    assert abs(metrics["Vzz"]) == approx(1.4109705)
    assert metrics["eta"] == approx(0.040096, abs=1e-5)
    summed = [el + nu for el, nu in zip(nucleus["V_el"], nucleus["V_nuc"])]
    assert summed == approx(nucleus["V_Tot"], abs=2e-7)


def test_the_parameter_run_carries_the_a_tensor_and_engine_eta(nuc):
    (nucleus,) = nuc["nuclei"]
    assert nucleus["a"]["I"] == approx(2.5)
    assert nucleus["a"]["P_MHz_au3"] == approx(1.5)
    assert nucleus["A_iso"] == approx(-11.8696)
    assert nucleus["A_Tot"] == approx([-11.5875, -11.9945, -12.0267])
    assert nucleus["A_FC"] == approx([-11.8696, -11.8696, -11.8696])
    assert nucleus["A_SD"] == approx([0.2821, -0.1250, -0.1571])
    assert nucleus["engine_eta"] == approx(0.040223)
    metrics = hyperfine_analysis.efg_metrics(nucleus)
    assert metrics["eta"] == approx(nucleus["engine_eta"], abs=1e-4)


def test_the_casscf_probe_reports_efg_only(casscf):
    header = casscf["header"]
    assert header["method"] == "CASSCF/ALL STATES AVERAGE"
    (nucleus,) = casscf["nuclei"]
    assert nucleus["A_Tot"] is None
    assert nucleus["hfc_flags"] == {"iso": False, "dip": False, "orb": False, "gauge": False}
    assert nucleus["V_Tot"] == approx([-0.0861351, -0.0861479, 0.172283])
    assert nucleus["rho0"] == approx(0.075618828)
    # the nuclear EFG triplet is the same as on the DFT probe (same geometry,
    # different principal-axis orientation)
    dft = parse_auto(DFT).sections["hyperfine"]["nuclei"][0]
    assert sorted(nucleus["V_nuc"]) == approx(sorted(dft["V_nuc"]))


# --- the conversions -----------------------------------------------------------


def test_the_nqcc_conversion_constant():
    # eQVzz/h for Q = 1 barn, Vzz = 1 a.u. -> 234.9648 MHz (CODATA derivation;
    # the literature value is 234.9647)
    assert hyperfine_analysis.nqcc_mhz(1.0, 1.0) == approx(234.9647, abs=2e-3)


def test_the_quadrupole_splitting_chain(dft):
    (nucleus,) = dft["nuclei"]
    metrics = hyperfine_analysis.efg_metrics(nucleus)
    nqcc = hyperfine_analysis.nqcc_mhz(0.5, metrics["Vzz"])
    assert nqcc == approx(-165.7642, rel=1e-4)
    deq = hyperfine_analysis.delta_e_q_mhz(0.5, metrics["Vzz"], metrics["eta"])
    assert deq == approx(82.9043, rel=1e-4)


def test_the_engine_quadrupole_block_matches_the_derived_constant():
    parsed = parse_auto(NUCQ).sections["hyperfine"]
    (nucleus,) = parsed["nuclei"]
    # engine: e**2qQ = -165.764347 MHz for Q = 0.5 barn on the all-electron probe
    assert nucleus["e2qQ_MHz"] == approx(-165.764347, rel=1e-6)
    assert nucleus["quad_tensor"]["Q_barn"] == approx(0.5)
    assert nucleus["quad_tensor"]["I"] == approx(2.5)
    metrics = hyperfine_analysis.efg_metrics(nucleus)
    derived = hyperfine_analysis.nqcc_mhz(0.5, metrics["Vzz"])
    assert nucleus["e2qQ_MHz"] == approx(derived, rel=1e-5)
    body = hyperfine_analysis.render(parsed, source=NUCQ.name)
    assert "engine quadrupole block: e**2qQ = -165.7643 MHz" in body


# --- the render ----------------------------------------------------------------


def test_the_render_reports_both_probe_routes(dft, nuc):
    body_zero = hyperfine_analysis.render(dft, source=DFT.name)
    assert "A tensor (MHz): all zeros" in body_zero
    assert "eta = 0.0401" in body_zero
    assert "Rho(0) = 727217" in body_zero
    assert "basis-domain dependent" in body_zero
    body_nuc = hyperfine_analysis.render(nuc, source=NUC.name, nuclear_Q_barn=0.5)
    assert "A(iso) = -11.8696" in body_nuc
    assert "(engine 0.0402)" in body_nuc
    assert "eQVzz/h = -165.7682 MHz" in body_nuc
    assert "DeltaE_Q = 82.9064 MHz" in body_nuc


def test_the_render_refuses_without_the_section():
    with pytest.raises(hyperfine_analysis.HyperfineError, match="hyperfine"):
        hyperfine_analysis.render({"present": False}, source="x.out")


def test_the_evidence_is_citable():
    evidence = hyperfine_analysis.evidence()
    assert references_section(evidence) is not None
    literature = [item for item in evidence if item.kind == "literature"]
    assert literature and "10.1063/1.5097151" in literature[0].ref
    assert "10.1063/1.1744480" in literature[0].ref
