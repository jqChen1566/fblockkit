"""Checks of the imaginary-mode displacement (recipe/imag_disp.py, menu 30).

The protocol is the NIFREC displacement stage (sum of the imaginary modes or
the most negative one, normalized, coords + vec * disp, rerun with the
frequency check; base_disp 0.1 Angstrom, growing by base_disp).  The ORCA-side
conventions are measured on the F + H2 fixture: the NORMAL MODES block prints
mass-weighted vectors of unit Euclidean norm (the Cartesian pattern is
v_i/sqrt(m_i), masses from the (A.U.) block), and the fixture chain is
engine-validated: the TS's imaginary mode (-90.48 cm**-1) displaced +/-0.1
Angstrom relaxes, both signs, to all-real-frequency minima at
-100.895757 Eh (mirror-image structures -- this mode is the near-linear bend,
so the two signs land on the two mirror images of the same shallow minimum).
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from fblockkit.parsers import parse_auto
from fblockkit.recipe import imag_disp as idisp
from fblockkit.recipe.imag_disp import ImagDispError

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"
FREQ_OUT = FIXTURES / "fhh_optts_freq.out"
BASE_IN = FIXTURES / "inputs" / "fhh_reopt.inp"


@pytest.fixture(scope="module")
def data():
    return idisp.imaginary_modes(parse_auto(FREQ_OUT))


def test_the_parser_reads_the_final_geometry_and_masses():
    result = parse_auto(FREQ_OUT)
    geometry = result.sections["final_geometry"]
    assert geometry["present"]
    assert [atom[0] for atom in geometry["atoms"]] == ["F", "H", "H"]
    assert geometry["masses"] == pytest.approx((18.998, 1.008, 1.008))


def test_the_parser_reads_the_normal_modes():
    modes = parse_auto(FREQ_OUT).sections["normal_modes"]
    assert modes["present"] and modes["n_coord"] == 9
    assert [mode["index"] for mode in modes["modes"]] == list(range(9))
    vector = modes["modes"][6]["vector"]
    # the printed column (the block's first group's second header is 6 7 8)
    assert vector[:4] == pytest.approx((0.005154, 0.005081, -0.000075, -0.549729))
    assert math.sqrt(sum(value * value for value in vector)) == pytest.approx(1.0, abs=1e-6)
    # the H-H stretch (mode 8) is the textbook pattern: the two H z components
    # are opposite, everything else is small
    stretch = modes["modes"][8]["vector"]
    assert stretch[5] == pytest.approx(-stretch[8], abs=1e-3)


def test_the_displacement_vector_de_weights_the_masses(data):
    vector = idisp.displacement_vector(data, selection="lowest")
    assert math.sqrt(sum(value * value for value in vector)) == pytest.approx(1.0, abs=1e-9)
    # the de-weighting convention, pinned numerically: the H1/F amplitude ratio
    # is (v_H/v_F) * sqrt(m_F/m_H) = (0.549729/0.005154) * sqrt(18.998/1.008)
    ratio = vector[3] / vector[0]
    assert ratio == pytest.approx(-463.1, rel=1e-3)
    # one imaginary mode in the fixture, so both selections agree
    assert idisp.displacement_vector(data, selection="sum") == pytest.approx(vector)


def test_the_displacement_and_the_input_edit(data):
    vector = idisp.displacement_vector(data, selection="lowest")
    moved, max_move = idisp.displaced_atoms(data["atoms"], vector, 0.1)
    assert max_move == pytest.approx(0.055, abs=1e-3)
    # the F stays nearly put while the H atoms move (the mode's character)
    assert abs(moved[0][3] - data["atoms"][0][3]) < 1e-4
    assert abs(moved[1][2] - data["atoms"][1][2]) > 0.02
    base = BASE_IN.read_text(encoding="utf-8")
    text = idisp.disp_input(base, moved)
    assert text.count("* xyz 0 2") == 1 and "r2SCAN-3c Opt Freq" in text
    # the block was replaced: the displaced F value is there, the old one is gone
    assert f"{moved[0][1]:.10f}" in text
    assert "-0.000056" not in text
    # a base without the Freq token gains it (the verification is the criterion)
    no_freq = base.replace("! r2SCAN-3c Opt Freq", "! r2SCAN-3c Opt")
    assert "! r2SCAN-3c Opt Freq" in idisp.disp_input(no_freq, moved)
    with pytest.raises(ImagDispError, match="no inline"):
        idisp.disp_input("%maxcore 2000\n! Opt Freq\n", moved)


def test_the_refusals_carry_next_steps(tmp_path):
    with pytest.raises(ImagDispError, match="no imaginary mode"):
        idisp.displacement_vector(
            {"imaginary": (), "vectors": {}, "atoms": (), "masses": ()}, selection="sum"
        )
    with pytest.raises(ImagDispError, match="unknown mode selection"):
        idisp.displacement_vector(
            {"imaginary": (), "vectors": {}, "atoms": (), "masses": ()},
            selection="both",
        )
    # a real ORCA output without a frequency job refuses with the next step
    with pytest.raises(ImagDispError, match="VIBRATIONAL FREQUENCIES"):
        idisp.imaginary_modes(parse_auto(FIXTURES / "n2_hf_clean.out"))


def test_the_fixture_chain_matches_the_generated_inputs(data):
    """The committed restart inputs are what the menu emits, and the engine
    runs carry the cure: both signs reach all-real frequencies."""
    vector = idisp.displacement_vector(data, selection="sum")
    base = BASE_IN.read_text(encoding="utf-8")
    for label, sign in (("disp_p", 1.0), ("disp_m", -1.0)):
        moved, _ = idisp.displaced_atoms(
            data["atoms"], tuple(sign * value for value in vector), 0.1
        )
        produced = idisp.disp_input(base, moved)
        committed = (FIXTURES / "inputs" / f"fhh_reopt.{label}.inp").read_text(
            encoding="utf-8"
        )
        assert produced == committed, label

    ts = parse_auto(FREQ_OUT)
    ts_block = ts.sections["frequencies"]["blocks"][-1]
    assert ts_block["n_imaginary"] == 1
    assert ts_block["min_imaginary"] == pytest.approx(-90.48, abs=0.01)

    results = {}
    for label in ("disp_p", "disp_m"):
        result = parse_auto(FIXTURES / f"fhh_reopt.{label}.out")
        results[label] = result
        # the cure: the rerun's frequency check finds no imaginary mode
        assert result.sections["frequencies"]["blocks"][-1]["n_imaginary"] == 0
    energy_p = results["disp_p"].sections["final_energy"]
    energy_m = results["disp_m"].sections["final_energy"]
    assert energy_p == pytest.approx(-100.895757053, abs=1e-8)
    assert energy_m == pytest.approx(-100.895756991, abs=1e-8)
    # the two signs are mirror images (this mode is the near-linear bend), so
    # they land on the same-energy shallow minimum with flipped x/y
    assert abs(energy_p - energy_m) < 1e-6
    atoms_p = results["disp_p"].sections["final_geometry"]["atoms"]
    atoms_m = results["disp_m"].sections["final_geometry"]["atoms"]
    # mirrored to within the optimisation tolerance (the two signs of the bend)
    assert atoms_p[1][1] == pytest.approx(-atoms_m[1][1], abs=1e-3)
    assert atoms_p[1][2] == pytest.approx(-atoms_m[1][2], abs=1e-3)
