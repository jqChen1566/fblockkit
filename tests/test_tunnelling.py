"""Checks of the tunnelling-relaxation models (analysis/tunnelling.py; menu 46).

The regression spine is the 18-complex literature table
(fixtures/qtm/tau_zeeman_benchmark.yaml): the g values from the JPCL 2018 SI
and log10(tau) from the PCCP 2020 Table 1, two independent data chains of the
same model.  The dipolar side is checked against closed-form analytics (a
single neighbour) plus the published <|E_sf|> calibration (three complexes
from the JPCL SI reproduce the published tau_calc through the 1/k_B density
factor).  The end-to-end chain uses the shipped CO+ SINGLE_ANISO fixture.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest
import yaml
from pytest import approx

from fblockkit.analysis import relaxation
from fblockkit.analysis import tunnelling as tn
from fblockkit.parsers import parse_auto

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
BENCHMARK = FIXTURES / "qtm" / "tau_zeeman_benchmark.yaml"
SINGLE_ANISO = FIXTURES / "single_aniso" / "co_aniso2.out"


# --- the literature regression (Yin & Li equivalent-Zeeman) -------------------


def test_the_zeeman_regression_table():
    """All 18 literature rows reproduce within the print precision."""
    data = yaml.safe_load(BENCHMARK.read_text(encoding="utf-8"))
    assert len(data["rows"]) == 18
    worst = 0.0
    for row in data["rows"]:
        produced = math.log10(tn.zeeman_tau(tuple(row["g"])))
        worst = max(worst, abs(produced - row["log_tau_zeeman"]))
    assert worst < 0.01


def test_the_axial_convention_is_the_largest_g():
    # the g triple is used through gXY^2 (symmetric in the two transverse
    # components), so ordering inside the pair does not matter
    assert tn.axial_parts((19.0, 0.01, 0.02)) == approx((math.hypot(0.01, 0.02), 19.0))
    assert tn.zeeman_tau((19.0, 0.01, 0.02)) == approx(tn.zeeman_tau((0.01, 19.0, 0.02)))


def test_the_rate_and_time_are_reciprocal():
    g = (0.005, 0.008, 19.4)
    rate = tn.zeeman_rate(g)
    assert tn.zeeman_tau(g) == approx(1.0 / (2.0 * rate), rel=1e-12)


def test_a_purely_axial_doublet_has_no_zeeman_channel():
    with pytest.raises(tn.TunnellingError, match="transverse g"):
        tn.zeeman_tau((0.0, 0.0, 20.0))


# --- the published <|E_sf|> calibration (the 1/k_B density factor) ------------


def test_the_e_sf_calibration_reproduces_the_published_times():
    """Three JPCL-SI complexes: <|E_sf|> (K) -> log10(tau) matches the paper."""
    hbar = 1.054571817e-34
    k_b = 1.380649e-23
    cases = (
        ("NAFMIT", 8.671e-9, 3.91),
        ("BAJSIQ", 1.530e-4, -4.59),
        ("OLUJEM", 4.303e-5, -3.48),
    )
    for refcode, e_sf_k, published in cases:
        rate = (2.0 * math.pi * k_b / hbar) * e_sf_k**2
        produced = math.log10(1.0 / (2.0 * rate))
        assert produced == approx(published, abs=0.01), refcode


# --- the dipolar model: closed-form anchor ------------------------------------


def test_the_single_neighbour_analytic_anchor():
    """One 45-degree neighbour: sigma_r equals the hand-derived value exactly,
    and sigma_i is exactly zero (the field has no y component)."""
    d = 3.0
    neighbour = tn.Neighbour(r=(d / math.sqrt(2), 0.0, d / math.sqrt(2)), moment=(0.0, 0.0, 5.0))
    g_center = (0.01, 0.01, 19.0)
    sigma_r2, sigma_i2 = tn.dipolar_sigma((neighbour,), g_center)
    mu_0 = 1.25663706212e-6
    mu_b = 9.2740100783e-24
    r_m = d * 1e-10
    b_x = 1.5 * mu_0 * 5.0 * mu_b / (4.0 * math.pi * r_m**3)
    c = (mu_b / 2.0) * 0.01 * b_x
    assert math.sqrt(sigma_r2) == approx(c, rel=1e-12)
    assert sigma_i2 == 0.0


def test_the_collinear_on_axis_neighbour_cannot_flip():
    """A neighbour on the axis with a collinear moment has a purely axial
    field: no transverse component, no tunnelling channel."""
    neighbour = tn.Neighbour(r=(0.0, 0.0, 6.0), moment=(0.0, 0.0, 9.5))
    sigma_r2, sigma_i2 = tn.dipolar_sigma((neighbour,), (0.01, 0.01, 19.0))
    assert sigma_r2 == 0.0
    assert sigma_i2 == 0.0


def test_dilution_increases_the_time_and_replays():
    """Dilution lengthens tau (fewer active neighbours) and the seeded stream
    makes the table byte-reproducible."""
    rng = __import__("random").Random(3)
    shell = []
    for _ in range(200):
        vector = [rng.gauss(0, 1) for _ in range(3)]
        norm = math.sqrt(sum(v * v for v in vector))
        shell.append(
            tn.Neighbour(r=tuple(10.0 * v / norm for v in vector), moment=(0.0, 0.0, 9.5))
        )
    shell = tuple(shell)
    g_center = (0.01, 0.01, 19.0)
    rows = tn.dilution_medians(shell, g_center, (1.0, 0.3, 0.05), repeats=30)
    assert rows[0]["log_tau_median"] < rows[1]["log_tau_median"] < rows[2]["log_tau_median"]
    again = tn.dilution_medians(shell, g_center, (1.0, 0.3, 0.05), repeats=30)
    assert rows == again


def test_the_neighbour_table_parses_and_refuses():
    table = "# comment\n3  0  0  0 0 9.5\n-3 0 0 0 0 9.5\n\n"
    neighbours = tn.parse_neighbour_table(table)
    assert len(neighbours) == 2
    assert neighbours[0].moment == (0.0, 0.0, 9.5)
    with pytest.raises(tn.TunnellingError, match="six numbers"):
        tn.parse_neighbour_table("1 2 3 4 5\n")
    with pytest.raises(tn.TunnellingError, match="empty"):
        tn.parse_neighbour_table("# only a comment\n")
    with pytest.raises(tn.TunnellingError, match="zero position"):
        tn.parse_neighbour_table("0 0 0 0 0 5\n")


# --- the end-to-end chain on the shipped fixture ------------------------------


def _co_levels():
    result = parse_auto(SINGLE_ANISO)
    segment = result.sections["single_aniso"]["segments"][-1]
    return tn.kd_levels(relaxation.group_metrics(segment))


def test_the_co_fixture_runs_through_the_models():
    levels = _co_levels()
    assert [level.index for level in levels] == [1, 2]
    assert levels[1].energy_cm1 == approx(29361.197)
    tau = tn.zeeman_tau(levels[0].g)
    assert math.log10(tau) == approx(-8.81, abs=0.05)  # the near-isotropic g gives fast QTM
    curve = tn.ueff_curve(levels, temperatures_K=(10.0, 300.0))
    # the excited doublet (29361 cm-1 = 42244 K) is far above the grid: the
    # 10 K weight underflows to exactly zero, the 300 K weight is e^-141
    assert curve[0]["ueff_cm1"] == 0.0
    assert curve[1]["ueff_cm1"] < 1e-50


def test_non_kramers_segments_are_refused():
    rows = [
        {"kramers": False, "energy_cm1": 0.0, "g_principal": (1.0, 2.0, 8.0), "index": 1}
    ]
    with pytest.raises(tn.TunnellingError, match="Kramers"):
        tn.kd_levels(rows)


def test_the_report_carries_the_untouched_window_note():
    levels = _co_levels()
    body = tn.render(
        source="co_aniso2.out",
        levels=levels,
        B_ave_mT=20.0,
        ueff_rows=tn.ueff_curve(levels, temperatures_K=(10.0, 300.0)),
    )
    assert "far above this temperature grid" in body
    assert "42244 K" in body


def test_the_robustness_guards():
    """Robustness round: the physical-range guards refuse with a next step."""
    with pytest.raises(tn.TunnellingError, match="Angstrom"):
        tn.dipolar_tau(
            (tn.Neighbour(r=(0.01, 0.0, 0.0), moment=(0.0, 0.0, 5.0)),),
            (0.01, 0.01, 19.0),
        )
    with pytest.raises(tn.TunnellingError, match="concentration"):
        tn.dilution_medians(
            (tn.Neighbour(r=(5.0, 0.0, 0.0), moment=(0.0, 0.0, 5.0)),),
            (0.01, 0.01, 19.0),
            (1.5,),
        )
    with pytest.raises(tn.TunnellingError, match="repeats"):
        tn.dilution_medians(
            (tn.Neighbour(r=(5.0, 0.0, 0.0), moment=(0.0, 0.0, 5.0)),),
            (0.01, 0.01, 19.0),
            (1.0,),
            repeats=0,
        )
    # a non-finite time renders as an explicit statement, not 'inf s'
    assert "effectively infinite" in tn._format_seconds(float("inf"))
    assert "effectively infinite" in tn._format_seconds(1e310)


def test_the_neighbour_report_block():
    shell = tuple(
        tn.Neighbour(r=(8.0 * math.cos(2 * math.pi * k / 8), 8.0 * math.sin(2 * math.pi * k / 8), 2.0),
                     moment=(0.0, 0.0, 9.5))
        for k in range(8)
    )
    g_center = (0.01, 0.01, 19.0)
    report = {
        "table": "probe.txt",
        "n_neighbours": len(shell),
        "dilution": tn.dilution_medians(shell, g_center, (1.0, 0.25), repeats=12),
        **tn.dipolar_tau(shell, g_center),
    }
    body = tn.render(
        source="co_aniso2.out",
        levels=_co_levels(),
        B_ave_mT=20.0,
        ueff_rows=tn.ueff_curve(_co_levels(), temperatures_K=(300.0,)),
        neighbour_report=report,
    )
    assert "Spin-dipolar model" in body
    assert "Dilution variant" in body
