"""Checks of the AILFT machinery (parsers/ailft.py + analysis/ailft.py; menu 38; Wave 5.5).

The fixture is a real ORCA 6.1.1 run: the Ni(2+) d8 free ion with the AILFT
driver requested through ``ActOrbs dOrbs`` (NEVPT2, three-triplet/fifteen-
singlet roots).  Anchors from the printed output: CASSCF-level Racah
B = 1328.1 / C = 4865.5 / C/B = 3.663, NEVPT2-level B = 1218.0 / C/B = 3.706;
the CASSCF-level RMS is 0.0 (intrinsic) and the NEVPT2-level total 457.4
cm^-1 (blocks 315.9 / 517.9), Pearson 1.000; the SOC constant ZETA_D =
664.14 cm^-1.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from pytest import approx

from fblockkit.analysis import ailft as ailft_analysis
from fblockkit.diagnosis import references_section
from fblockkit.parsers import parse_auto

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
AILFT = FIXTURES / "ailft"
CEF3 = AILFT / "cef3_ailft.out"


@pytest.fixture(scope="module")
def parsed():
    return parse_auto(AILFT / "ni_ailft.out").sections["ailft"]


# --- the parse -----------------------------------------------------------------


def test_the_ailft_blocks_parse_with_both_levels(parsed):
    assert parsed["present"]
    assert parsed["header"]["configuration"] == "d8"
    assert parsed["header"]["ci_blocks"] == 2
    assert parsed["header"]["mo_range"] == (9, 13)
    assert parsed["header"]["center_atom"] == 0
    assert set(parsed["levels"]) == {"casscf", "nevpt2"}


def test_the_casscf_level_parameters(parsed):
    level = parsed["levels"]["casscf"]
    labels = [entry["label"] for entry in level["slater_condon"]]
    assert labels == ["F0dd(from 2el Ints)", "F2dd", "F4dd"]
    assert level["slater_condon"][0]["fixed"] is True
    assert level["slater_condon"][1]["cm1"] == approx(99138.0)
    racah = level["racah"]
    assert racah["B"]["cm1"] == approx(1328.1)
    assert racah["C"]["cm1"] == approx(4865.5)
    assert racah["C_over_B"] == approx(3.663)
    assert level["lft_gbw"] == "ni_ailft.casscf.lft.gbw"


def test_the_nevpt2_level_drops_f0(parsed):
    level = parsed["levels"]["nevpt2"]
    labels = [entry["label"] for entry in level["slater_condon"]]
    assert labels == ["F2dd", "F4dd"]  # F0 / Racah A are not printed at this level
    assert level["racah"]["B"]["cm1"] == approx(1218.0)
    assert level["racah"]["C_over_B"] == approx(3.706)
    assert level["lft_gbw"] == "ni_ailft.nevpt2.lft.gbw"


def test_the_eigenfunctions_and_fit_statistics(parsed):
    split = ailft_analysis.lf_splitting(parsed["levels"]["casscf"])
    assert split is not None and split["spread_cm1"] == approx(0.0)
    assert len(split["energies_cm1"]) == 5
    fit = parsed["fit"]
    assert fit["casscf"]["total_rms_cm1"] == approx(0.0)
    assert fit["nevpt2"]["total_rms_cm1"] == approx(457.4)
    blocks = fit["nevpt2"]["blocks"]
    assert [block["rms_cm1"] for block in blocks] == approx([315.9, 517.9])
    assert len(blocks[0]["roots"]) == 10
    first = blocks[0]["roots"][0]
    assert first["ai_ev"] == approx(-0.002) and first["overlap"] == approx(0.981)
    assert fit["nevpt2"]["pearson"] == approx(1.0)


def test_the_soc_constants(parsed):
    soc = parsed["soc"]
    assert soc["bases"] == ["casscf"]
    assert soc["zeta_cm1"] == approx(664.1)
    assert soc["summary"] == {"ZETA_D": approx(664.14)}
    assert soc["rms_cm1"] == approx(0.0)


# --- the render ----------------------------------------------------------------


def test_the_render_reports_levels_fit_and_soc(parsed):
    body = ailft_analysis.render(parsed, source="ni_ailft.out")
    assert "d8 configuration" in body
    assert "Racah (cm-1): B = 1328.1" in body
    assert "C/B = 3.663" in body
    assert "total RMS = 457.4 cm-1" in body
    assert "ZETA_D = 664.14 cm-1" in body
    assert "beta" not in body  # no free-ion reference passed


def test_the_nephelauxetic_ratios_need_the_caller_reference(parsed):
    body = ailft_analysis.render(
        parsed, source="ni_ailft.out", free_ion_B_cm1=1050.0, free_ion_zeta_cm1=630.0
    )
    assert "beta = B/B0 = 1.265" in body
    assert "relativistic nephelauxetic ratio = 1.054" in body


def test_the_render_refuses_without_an_ailft_section():
    with pytest.raises(ailft_analysis.AilftError, match="AILFT"):
        ailft_analysis.render({"present": False}, source="x.out")


def test_the_f_shell_sample_parses_with_its_own_conventions():
    parsed = parse_auto(CEF3).sections["ailft"]
    assert parsed["header"]["configuration"] == "f1"
    casscf = parsed["levels"]["casscf"]
    # the f-shell VLFT table goes straight to its header (no title line)
    assert casscf["vlft"]["labels"] == ["Orbital", "f0", "f+1", "f-1", "f+2", "f-2", "f+3", "f-3"]
    assert len(casscf["vlft"]["matrix"]) == 7
    # a single f electron has no electron repulsion: only F0 (fixed) and B = 0
    assert [entry["label"] for entry in casscf["slater_condon"]] == ["F0ff"]
    assert casscf["racah"]["B"]["cm1"] == approx(0.0)
    split = ailft_analysis.lf_splitting(casscf)
    assert split["spread_cm1"] == approx(2333.0)
    assert parsed["soc"]["summary"] == {"ZETA_F": approx(638.98)}
    body = ailft_analysis.render(parsed, source="cef3_ailft.out")
    assert "ZETA_F = 638.98 cm-1" in body


def test_the_evidence_is_citable():
    evidence = ailft_analysis.evidence()
    assert references_section(evidence) is not None
    literature = [item for item in evidence if item.kind == "literature"]
    assert literature and "10.1021/acs.jpca.9b11227" in literature[0].ref
    assert "10.1021/acs.inorgchem.7b00642" in literature[0].ref
