"""Checks of the polynuclear-magnetism machinery (parsers/poly_aniso.py +
analysis/poly_aniso.py; menu 39).

The fixture is a real two-center probe of the ORCA 6.1.1 `otool_poly_aniso`
driver (POLY_ANISO v1.0.0): two Co(II) single-ion data files (the menu-36
SINGLE_ANISO probe outputs, copied as aniso_1.input / aniso_2.input), one
Lines-1 pair with J = 0.1 cm^-1 and centered coordinates.  Anchors from the
printed output: g = 2.00000 / 2.00000 / 1.99952 and the spin-orbit spectrum
0 / 0 / 29361.197 / 29361.197 cm^-1 (matching the menu-36 analysis of the
same data); the LINES-1 decomposition weights 32.708 / 93.839 / 11.155 %;
coupled-state relatives 0 / 8.224256e-06 / 0.051286564374 / 0.117103624688
cm^-1; chiT 0.50154932 -> 0.75042950 cm3 K mol-1 over 101 points; the Van
Vleck Z main value 1.500582 at 0.0001 K.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from pytest import approx

from fblockkit.analysis import poly_aniso as poly_analysis
from fblockkit.diagnosis import references_section
from fblockkit.parsers import parse_auto
from fblockkit.parsers.poly_aniso import parse_poly_aniso

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
PROBE = FIXTURES / "poly_aniso" / "two_center_probe.out"


@pytest.fixture(scope="module")
def parsed():
    return parse_auto(PROBE).sections["poly_aniso"]


# --- the parse -----------------------------------------------------------------


def test_the_probe_output_is_recognised_and_parsed(parsed):
    assert parsed["present"]
    assert parsed["counts"] == {"independent": 2, "all": 2}


def test_the_centers_echo_their_single_ion_data(parsed):
    centers = parsed["centers"]
    assert [c["data_file"] for c in centers] == ["aniso_1.input", "aniso_2.input"]
    assert centers[0]["coords"] == approx([0.0, 0.0, 0.0])
    assert centers[1]["coords"] == approx([0.0, 0.0, 3.7])
    for center in centers:
        assert center["n_so_states"] == 4
        assert center["so_spectrum_cm1"] == approx([0.0, 0.0, 29361.197, 29361.197])
        g = center["g_tensor"]
        assert g["X"]["value"] == approx(2.0)
        assert g["Z"]["value"] == approx(1.99952)


def test_the_exchange_block(parsed):
    exchange = parsed["exchange"]
    assert exchange["coupled_states"] == 4
    assert exchange["pairs_total"] == 1
    assert exchange["lines1"] == "INCLUDED"
    assert exchange["lines1_pairs"] == 1
    assert exchange["dipole_dipole"] == "INCLUDED"
    assert exchange["ito"] == "INCLUDED"
    assert exchange["pairs_j"] == [{"pair": 1, "centers": (1, 2), "J_cm1": approx(0.1)}]


def test_the_interaction_decomposition(parsed):
    interaction = parsed["interaction"]
    assert len(interaction) == 8  # 2 pairs-models x 4 terms
    first = interaction[0]
    assert first["model"] == "LINES-1"
    assert first["term"] == "Full Interaction"
    assert first["weight_percent"] == approx(100.0)
    assert first["matrix"][0][0] == approx(0.07447400039159)
    weights = {
        entry["term"]: entry["weight_percent"]
        for entry in interaction
        if entry["model"] == "LINES-1"
    }
    assert weights["Isotropic Term"] == approx(32.708)
    assert weights["Symmetric Term"] == approx(93.839)
    assert weights["Anti-Symmetric Term"] == approx(11.155)


def test_the_coupled_states(parsed):
    states = parsed["coupled_states"]
    assert [s["state"] for s in states] == [1, 2, 3, 4]
    assert [s["relative_cm1"] for s in states] == approx(
        [0.0, 8.224256e-06, 0.051286564374, 0.117103624688]
    )
    assert states[0]["lines_cm1"] == approx(-0.025)


def test_the_population_and_expectation_tables(parsed):
    assert len(parsed["population"]) == 8
    assert parsed["population"][0] == {
        "state": 1,
        "basis_set": 1,
        "weight": approx([0.5, 0.5]),
    }
    assert len(parsed["expectation"]["MS"]) == 8
    assert len(parsed["expectation"]["LJ"]) == 8
    first = parsed["expectation"]["MS"][0]
    assert first["state"] == 1 and first["site"] == 1
    assert first["vec1"] == approx([0.0, 0.0, 1.0])
    assert first["vec2"] == approx([0.0, 0.0, -0.5])


def test_the_chiT_table(parsed):
    chit = parsed["chiT"]
    assert len(chit) == 101
    assert chit[0]["T"] == approx(0.0001)
    assert chit[0]["chiT"] == approx(0.50154932)
    assert chit[-1]["T"] == approx(300.0)
    assert chit[-1]["chiT"] == approx(0.75042950)


def test_the_van_vleck_sequence(parsed):
    vv = parsed["van_vleck"]
    assert len(vv) == 101
    assert vv[0]["T"] == approx(0.0001)
    assert vv[0]["main_values"] == approx([0.001914, 0.002153, 1.500582])
    assert vv[-1]["T"] == approx(300.0)
    assert vv[-1]["main_values"] == approx([0.750350, 0.750442, 0.750496])


def test_the_notes_and_the_absent_section():
    from fblockkit.parsers.base import read_text

    parsed = parse_poly_aniso(read_text(PROBE).splitlines())
    assert any("skipped by the user" in note for note in parsed["notes"])
    assert "finished ok" in parsed["notes"]
    assert parse_poly_aniso(["not a poly aniso output"])["present"] is False


# --- the render ----------------------------------------------------------------


def test_the_render_reports_centers_exchange_and_tables(parsed):
    body = poly_analysis.render(parsed, source=PROBE.name)
    assert "Polynuclear magnetism (POLY_ANISO) report" in body
    assert "centers: 2 independent, 2 in total" in body
    assert "g = 2.0000 / 2.0000 / 1.9995" in body
    assert "J = 0.10000 cm-1" in body
    assert "Isotropic 32.708% / Symmetric 93.839% / Anti-Symmetric 11.155%" in body
    assert "chiT = 0.501549 -> 0.750429 cm3 K mol-1" in body
    assert "1.500582" in body
    assert "never computes or fits them" in body


def test_the_render_refuses_without_a_poly_aniso_section():
    with pytest.raises(poly_analysis.PolyAnisoError, match="POLY_ANISO"):
        poly_analysis.render({"present": False}, source="x.out")


def test_the_evidence_is_citable():
    evidence = poly_analysis.evidence()
    assert references_section(evidence) is not None
    literature = [item for item in evidence if item.kind == "literature"]
    assert literature and "10.1002/chem.201605102" in literature[0].ref
