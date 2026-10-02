"""Checks of the magnetocaloric machinery (analysis/magnetocaloric.py; menu 41).

The magnetization fixture is the POLY_ANISO probe with HINT/TMAG
(fixtures/magnetocaloric/poly_mh.out; two Co(II) centers, invented J -- a
format-and-flow probe).  Anchors: the merged table spans 6 temperatures x
71 fields; the Maxwell route gives -DeltaS = 10.8404 J mol-1 K-1 at
T = 1.90 K and 7 T as the probe maximum, decaying with temperature
(9.225 at 2.75 K mid).  The levels route is anchored on the menu-36 Co
spectrum (degenerate ground doublet -> S = R ln 2 = 5.7631) and on a
synthetic non-degenerate ladder exercising the limits.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest
from pytest import approx

from fblockkit.analysis import magnetocaloric as mc
from fblockkit.diagnosis import references_section
from fblockkit.parsers import parse_auto

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
MH = FIXTURES / "magnetocaloric" / "poly_mh.out"
CO = FIXTURES / "single_aniso" / "co_aniso2.out"


@pytest.fixture(scope="module")
def magnet():
    return parse_auto(MH).sections["poly_aniso"]["magnetization"]


# --- the table and the Maxwell route -------------------------------------------


def test_the_magnetization_table_merges(magnet):
    assert magnet["temperatures_K"] == approx([1.8, 2.0, 2.5, 3.0, 4.0, 5.0])
    assert len(magnet["fields_T"]) == 71
    assert magnet["fields_T"][0] == approx(0.0001)
    assert magnet["fields_T"][-1] == approx(7.0)
    assert {len(row) for row in magnet["M_muB"]} == {6}
    assert magnet["M_muB"][0][0] == approx(7.60807e-05)


def test_the_maxwell_entropy_change(magnet):
    rows = mc.maxwell_delta_s(
        magnet["temperatures_K"], magnet["fields_T"], magnet["M_muB"]
    )
    assert [row["T_K"] for row in rows] == approx([1.9, 2.25, 2.75, 3.5, 4.5])
    first, last = rows[0], rows[-1]
    # at 7 T: -DeltaS = 10.8404 (1.9 K) .. 5.821 (4.5 K)
    assert -first["delta_S"][-1] == approx(10.8404, rel=1e-3)
    assert -last["delta_S"][-1] == approx(5.821, rel=1e-3)
    # monotone decay with temperature at fixed field
    at_7T = [-row["delta_S"][-1] for row in rows]
    assert at_7T == sorted(at_7T, reverse=True)
    assert at_7T[0] == approx(max(at_7T), rel=1e-9)


def test_the_maxwell_render(magnet):
    body = mc.render_maxwell(magnet, source=MH.name)
    assert "route: magnetization table (6 temperatures x 71 fields" in body
    assert "maximum: -DeltaS = 10.8404 J mol-1 K-1 at T = 1.90 K, H = 7.000 T" in body
    assert "positive = direct MCE" in body


# --- the levels route ----------------------------------------------------------


def test_the_levels_route_on_the_co_spectrum():
    spectrum = parse_auto(CO).sections["single_aniso"]["segments"][0]["soc_spectrum_cm1"]
    assert spectrum == approx([0.0, 0.0, 29361.197, 29361.197])
    values = mc.entropy_from_levels(spectrum, [1.0, 10.0, 300.0])
    # the degenerate ground doublet pins S = R ln 2 across the range
    for value in values:
        assert value == approx(mc._R_J_MOL_K * math.log(2), rel=1e-6)


def test_the_levels_route_limits_on_a_synthetic_ladder():
    ladder = [0.0, 100.0, 200.0]
    low = mc.entropy_from_levels(ladder, [1.0])[0]
    mid = mc.entropy_from_levels(ladder, [300.0])[0]
    high = mc.entropy_from_levels(ladder, [100000.0])[0]
    assert low == approx(0.0, abs=1e-3)
    assert 0.0 < mid < mc._R_J_MOL_K * math.log(3)
    assert high == approx(mc._R_J_MOL_K * math.log(3), rel=1e-3)
    # monotone in T
    grid = mc.entropy_from_levels(ladder, [10.0, 50.0, 100.0, 300.0])
    assert grid == sorted(grid)


def test_the_levels_render():
    body = mc.render_levels(
        parse_auto(CO).sections["single_aniso"]["segments"][0]["soc_spectrum_cm1"],
        source=CO.name,
    )
    assert "route: spin-orbit levels (4 levels" in body
    assert "high-temperature check: R ln(N) = 11.5263" in body
    assert "upper bound" in body


def test_the_refusals():
    with pytest.raises(mc.MagnetocaloricError, match="level spectrum"):
        mc.render_levels([], source="x.out")
    with pytest.raises(mc.MagnetocaloricError, match="magnetization table"):
        mc.render_maxwell({"temperatures_K": [1.0]}, source="x.out")


def test_the_evidence_is_citable():
    evidence = mc.evidence()
    assert references_section(evidence) is not None
    literature = [item for item in evidence if item.kind == "literature"]
    assert literature and "10.3390/ma13020485" in literature[0].ref
    assert "10.1016/j.jmmm.2019.165933" in literature[0].ref


def test_the_levels_plot_csv_mirrors_the_report_grid():
    """The plot-ready companion of the levels route: the report's default
    temperature grid, the same canonical-ensemble numbers."""
    spectrum = parse_auto(CO).sections["single_aniso"]["segments"][0]["soc_spectrum_cm1"]
    csv_text = mc.levels_plot_csv(spectrum)
    lines = csv_text.splitlines()
    assert lines[0] == "temperature_K,entropy_J_per_mol_per_K"
    assert len(lines) == 1 + 9
    temperatures = [float(line.split(",")[0]) for line in lines[1:]]
    assert temperatures == [1, 2, 5, 10, 20, 50, 100, 200, 300]
    last = lines[-1].split(",")
    assert float(last[1]) == approx(mc.entropy_from_levels(spectrum, [300.0])[0])
    assert csv_text.endswith("\n") and "\r" not in csv_text


def test_the_maxwell_plot_csv_carries_the_full_field_grid(magnet):
    """The plot-ready companion of the Maxwell route: tidy long form over the
    full printed field grid (the report shows five probe fields), the same
    values with the report's sign convention (the column is -DeltaS)."""
    csv_text = mc.maxwell_plot_csv(magnet)
    lines = csv_text.splitlines()
    assert lines[0] == "T_mid_K,field_T,minus_delta_S_J_per_mol_per_K"
    assert len(lines) == 1 + 5 * 71  # 5 temperature pairs x 71 fields
    values = [float(line.split(",")[2]) for line in lines[1:]]
    assert max(values) == approx(10.8404, abs=5e-5)  # the report's probe maximum
    # the first block starts at the first temperature midpoint over the field grid
    first = lines[1].split(",")
    assert float(first[0]) == approx(0.5 * (1.8 + 2.0))
    assert float(first[1]) == approx(0.0001)
    assert csv_text.endswith("\n") and "\r" not in csv_text
