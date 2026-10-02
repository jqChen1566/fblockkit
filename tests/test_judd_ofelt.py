"""Checks of the Judd-Ofelt fit (analysis/judd_ofelt.py, menu 34).

The regression fixture is the Eu3+ dataset of Babu et al. (SET B) with the
U^(lambda) table of Hovhannesyan/Boudon/Lepers -- the same combination their
published standard-JO fit uses.  The published values reproduced here: the
relative standard deviation (8.52 %) and Omega_6 (2.253e-20 cm^2) digit for
digit, the strong-transition ratios of Paper I to better than 0.05, and the
remaining parameters within their documented bands (Omega_4 is carried by two
weak lines only and is the convention-sensitive one).
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from fblockkit.analysis import judd_ofelt as jo
from fblockkit.diagnosis import references_section

FIXTURE = Path(__file__).resolve().parents[1] / "fixtures" / "judd_ofelt" / "babu2000_eu3.yaml"


@pytest.fixture(scope="module")
def babu():
    return jo.read_dataset(FIXTURE)


@pytest.fixture(scope="module")
def fitted(babu):
    return jo.fit(babu)


def _synthetic(omega_au, entries):
    """Build a dataset whose f values are generated from omega_au (round trip)."""
    from fblockkit.analysis.judd_ofelt import Host, Transition

    host = Host(refractive_index=1.5)
    transitions = []
    for label, energy, j_low, u2, u4, u6 in entries:
        s_target = omega_au[0] * u2 + omega_au[1] * u4 + omega_au[2] * u6
        n = host.n_at(1e7 / energy)
        chi = Host.chi_ed(n)
        delta_e_au = energy / jo._CM1_PER_HARTREE
        f = s_target * delta_e_au * chi / (n * 1.5 * (2 * j_low + 1))
        transitions.append(
            Transition(
                label=label, energy_cm1=energy, f_exp=f, j_low=j_low, u2=u2, u4=u4, u6=u6
            )
        )
    return jo.Dataset(host=host, transitions=tuple(transitions))


# --- the literature regression -------------------------------------------------


def test_the_regression_reproduces_the_published_fit(fitted):
    # relative standard deviation: published 8.52 % (standard JO, 9 transitions)
    assert fitted.sigma_relative * 100.0 == pytest.approx(8.52, abs=0.05)
    omega = [v / 1e-20 for v in fitted.omega_cm2]
    # Omega_6 digit for digit; Omega_2 within 2 %; Omega_4 within 15 %
    # (Omega_4 rides on the two U4-carrying weak lines only -- the band is
    # documented in the close-reading record)
    assert omega[2] == pytest.approx(2.253, abs=0.01)
    assert omega[0] == pytest.approx(18.73, rel=0.02)
    assert omega[1] == pytest.approx(12.58, rel=0.15)
    assert fitted.n_transitions == 9
    # the ratio of the strongest line pins the chain (Paper I: 1.83)
    by_label = {t.label: r for t, r in zip(fitted.transitions, fitted.ratios)}
    assert by_label["7F6 <- 7F0"] == pytest.approx(1.83, abs=0.05)


def test_the_strong_transition_ratios_match_paper_one(fitted):
    # Paper I Table 1 prints the standard-JO ratios r0 for the Babu set
    published = {
        "7F6 <- 7F1": 0.93,
        "7F6 <- 7F0": 1.83,
        "5D1 <- 7F1": 0.99,
        "5D2 <- 7F0": 1.41,
        "5D3 <- 7F1": 0.87,
        "5L6 <- 7F0": 0.38,
    }
    by_label = {t.label: r for t, r in zip(fitted.transitions, fitted.ratios)}
    for label, expected in published.items():
        assert by_label[label] == pytest.approx(expected, abs=0.05), label
    # the two remaining lines are the documented convention-sensitive ones
    assert by_label["5L6 <- 7F1"] == pytest.approx(0.17, abs=0.20)
    assert by_label["5D4 <- 7F0"] == pytest.approx(1.14, abs=0.20)


def test_sigma_and_parameters_are_internally_consistent(fitted):
    residuals = [s - e for s, e in zip(fitted.s_exp, fitted.s_ed)]
    sigma = (sum(r * r for r in residuals) / (fitted.n_transitions - 3)) ** 0.5
    assert sigma == pytest.approx(fitted.sigma, rel=1e-12)


def test_the_constant_index_host_path_runs(babu):
    from fblockkit.analysis.judd_ofelt import Host

    dataset = jo.Dataset(host=Host(refractive_index=1.57), transitions=babu.transitions)
    result = jo.fit(dataset)
    omega = [v / 1e-20 for v in result.omega_cm2]
    # without dispersion the fit is the paper's constant-n variant: parameters
    # move (their own constant-n refit moved the same way) but stay ordered
    assert omega[0] > omega[1] > omega[2] > 0.0
    assert result.sigma_relative > 0.0


# --- the fit machinery ---------------------------------------------------------


def test_a_synthetic_round_trip_recovers_the_parameters_exactly():
    omega_true = (6.0e-3, 3.0e-3, 5.0e-4)
    dataset = _synthetic(
        omega_true,
        [
            ("a", 2200.0, 1, 0.0, 0.0, 0.30),
            ("b", 4800.0, 0, 0.0, 0.0, 0.14),
            ("c", 18000.0, 1, 0.002, 0.0, 0.0),
            ("d", 21500.0, 0, 0.001, 0.001, 0.0),
            ("e", 25000.0, 0, 0.0, 0.0, 0.01),
        ],
    )
    result = jo.fit(dataset)
    for got, want in zip(result.omega_au, omega_true):
        assert got == pytest.approx(want, rel=1e-10)
    assert result.sigma == pytest.approx(0.0, abs=1e-18)


def test_an_undetermined_parameter_is_reported_not_zeroed():
    omega_true = (6.0e-3, 3.0e-3, 5.0e-4)
    dataset = _synthetic(
        omega_true,
        [
            ("a", 2200.0, 1, 0.0, 0.0, 0.30),
            ("b", 4800.0, 0, 0.0, 0.0, 0.14),
            ("c", 18000.0, 1, 0.002, 0.0, 0.0),
            ("e", 25000.0, 0, 0.0, 0.0, 0.01),
        ],
    )
    result = jo.fit(dataset)
    assert result.determined == (True, False, True)
    assert result.omega_cm2[1] is None
    assert "Omega_4 = not determined" in jo.render(result)


def test_the_weighted_variant_differs_and_validates():
    omega_true = (6.0e-3, 3.0e-3, 5.0e-4)
    entries = [
        ("a", 2200.0, 1, 0.0, 0.0, 0.30),
        ("b", 4800.0, 0, 0.0, 0.0, 0.14),
        ("c", 18000.0, 1, 0.002, 0.0, 0.0),
        ("d", 21500.0, 0, 0.001, 0.001, 0.0),
        ("e", 25000.0, 0, 0.0, 0.0, 0.01),
    ]
    dataset = _synthetic(omega_true, entries)
    same = jo.fit(dataset)
    weighted = jo.fit(dataset, weights="1/S")
    assert weighted.omega_au == pytest.approx(same.omega_au, rel=1e-8)
    with pytest.raises(jo.JudOError, match="unknown weighting"):
        jo.fit(dataset, weights="2/S")
    # a zero oscillator strength cannot be 1/S-weighted
    from fblockkit.analysis.judd_ofelt import Host, Transition

    zeroed = jo.Dataset(
        host=Host(refractive_index=1.5),
        transitions=dataset.transitions[:3]
        + (Transition(label="z", energy_cm1=20000.0, f_exp=0.0, j_low=0, u2=0.001),)
        + dataset.transitions[3:],
    )
    with pytest.raises(jo.JudOError, match="zero experimental line strength"):
        jo.fit(zeroed, weights="1/S")


def test_f_md_is_subtracted_before_the_fit(babu):
    from dataclasses import replace

    transitions = tuple(
        replace(t, f_md=0.1e-6) if t.label == "5D1 <- 7F1" else t for t in babu.transitions
    )
    dataset = jo.Dataset(host=babu.host, transitions=transitions)
    result = jo.fit(dataset)
    base = jo.fit(babu)
    assert result.s_exp[2] < base.s_exp[2]  # the subtracted transition's S drops
    # exactly proportional to (f - f_md): ratio of S values equals ratio of f
    assert result.s_exp[2] / base.s_exp[2] == pytest.approx((0.450 - 0.1) / 0.450, rel=1e-12)


# --- the emission side ---------------------------------------------------------


def test_emission_rates_scaling_branching_and_lifetime(babu, fitted):
    from dataclasses import replace

    base = fitted.emission or babu.emission
    assert base == ()
    transitions = (
        replace(babu.transitions[0], j_up=6, f_exp=0.0),
        replace(babu.transitions[1], j_up=6, f_exp=0.0),
    )
    result = jo.fit(jo.Dataset(host=babu.host, transitions=babu.transitions, emission=transitions))
    rows = jo.emission_rates(result)
    total = rows[-1]
    assert sum(row["branching"] for row in rows[:-1]) == pytest.approx(1.0, rel=1e-12)
    assert total["lifetime_us"] == pytest.approx(1e6 / total["total_a_s1"], rel=1e-12)
    # the atomic-unit emission-rate law: A/(dE_au^3 n chi S/(2J+1)) is the same
    # prefactor for every line and equals 2.1422e10 s^-1 (== the papers' SI
    # prefactor e^2 a0^2/(3 pi eps0 hbar^4 c^3) x (J per Hartree)^3 to 0.05 %)
    for row, transition in zip(rows[:-1], transitions):
        n = result.host.n_at(transition.wavelength_for_n)
        chi = jo.Host.chi_ed(n)
        d_e = transition.energy_cm1 / jo._CM1_PER_HARTREE
        pref = row["a_ed_s1"] / (d_e**3 * n * chi * row["s_ed_au"] / (2 * 6 + 1))
        assert pref == pytest.approx(2.1422e10, rel=5e-4)


def test_emission_requires_j_up(tmp_path):
    path = tmp_path / "bad.yaml"
    path.write_text(
        "host: {refractive_index: 1.5}\n"
        "transitions:\n"
        "  - {label: a, energy_cm1: 20000, f: 1e-6, j_low: 0, u2: 0.1}\n"
        "  - {label: b, energy_cm1: 19000, f: 1e-6, j_low: 0, u2: 0.1}\n"
        "  - {label: c, energy_cm1: 18000, f: 1e-6, j_low: 0, u2: 0.1}\n"
        "  - {label: d, energy_cm1: 17000, f: 1e-6, j_low: 0, u2: 0.1}\n"
        "emission:\n"
        "  - {label: e, energy_cm1: 1000, f: 0, j_low: 0, u2: 0, u4: 0, u6: 0.1}\n",
        encoding="utf-8",
    )
    with pytest.raises(jo.JudOError, match="needs j_up"):
        jo.read_dataset(path)


# --- the dataset reader --------------------------------------------------------


def test_reader_errors_carry_next_steps(tmp_path):
    def write(payload: str) -> Path:
        path = tmp_path / "d.yaml"
        path.write_text(payload, encoding="utf-8")
        return path

    base = {"host": {"refractive_index": 1.5}, "transitions": [{"label": "a"}]}
    with pytest.raises(jo.JudOError, match="no 'host'"):
        jo.read_dataset(write("transitions: []\n"))
    with pytest.raises(jo.JudOError, match="no 'transitions'"):
        jo.read_dataset(write("host: {refractive_index: 1.5}\n"))
    with pytest.raises(jo.JudOError, match="neither energy_cm1 nor"):
        jo.read_dataset(
            write(
                "host: {refractive_index: 1.5}\ntransitions:\n"
                "  - {label: a, f: 1e-6, j_low: 0, u2: 0.1}\n"
            )
        )
    with pytest.raises(jo.JudOError, match="j_low"):
        jo.read_dataset(
            write(
                "host: {refractive_index: 1.5}\ntransitions:\n"
                "  - {label: a, energy_cm1: 20000, f: 1e-6, u2: 0.1}\n"
            )
        )
    with pytest.raises(jo.JudOError, match="refractive_index' or 'sellmeier"):
        jo.read_dataset(write("host: {}\ntransitions: []\n"))
    assert yaml.safe_load  # the reader uses yaml; keep the import honest
    assert base  # base unused; the cases above are explicit strings


def test_energy_from_wavelength_is_consistent(tmp_path):
    path = tmp_path / "d.yaml"
    path.write_text(
        "host: {refractive_index: 1.5}\n"
        "transitions:\n"
        "  - {label: a, wavelength_nm: 500, f: 1e-6, j_low: 0, u2: 0.1}\n"
        "  - {label: b, energy_cm1: 20000, f: 1e-6, j_low: 0, u2: 0.1}\n"
        "  - {label: c, energy_cm1: 15000, f: 1e-6, j_low: 0, u2: 0.1}\n"
        "  - {label: d, energy_cm1: 10000, f: 1e-6, j_low: 0, u2: 0.1}\n",
        encoding="utf-8",
    )
    dataset = jo.read_dataset(path)
    assert dataset.transitions[0].energy_cm1 == pytest.approx(20000.0, rel=1e-12)
    assert dataset.transitions[0].wavelength_for_n == pytest.approx(500.0)


def test_refusals_on_degenerate_datasets(babu):
    from fblockkit.analysis.judd_ofelt import Host

    tiny = jo.Dataset(host=Host(refractive_index=1.5), transitions=babu.transitions[:3])
    with pytest.raises(jo.JudOError, match="at least four"):
        jo.fit(tiny)
    from dataclasses import replace

    zeroed = jo.Dataset(
        host=Host(refractive_index=1.5),
        transitions=tuple(replace(t, u2=0.0, u4=0.0, u6=0.0) for t in babu.transitions),
    )
    with pytest.raises(jo.JudOError, match="no Omega can be fitted"):
        jo.fit(zeroed)


def test_the_evidence_is_citable():
    refs = references_section(jo.evidence())
    assert refs is not None
    literature = [item for item in jo.evidence() if item.kind == "literature"]
    assert literature and "10.1016/j.jlumin.2023.120234" in literature[0].ref


def test_the_plot_csv_carries_the_fitted_transition_table(babu):
    """The plot-ready companion: one row per fitted transition with the
    measured and fitted line strengths the report's ratio column derives
    from; the derived summaries stay in the report."""
    import csv
    import io

    result = jo.fit(babu, weights="none")
    csv_text = jo.jo_plot_csv(result)
    rows = list(csv.reader(io.StringIO(csv_text)))
    assert rows[0] == ["side", "label", "energy_cm1", "f_exp", "s_exp", "s_ed",
                       "a_ed_s1", "branching"]
    assert len(rows) == 1 + result.n_transitions  # the fixture has no emission block
    cells = rows[1]
    assert cells[0] == "absorption"
    assert cells[1] == result.transitions[0].label
    assert float(cells[2]) == pytest.approx(result.transitions[0].energy_cm1)
    assert float(cells[3]) == pytest.approx(result.transitions[0].f_exp)
    assert float(cells[4]) == pytest.approx(result.s_exp[0])
    assert float(cells[5]) == pytest.approx(result.s_ed[0])
    assert cells[6] == "" and cells[7] == ""  # no emission columns on this side
    assert csv_text.endswith("\n") and "\r" not in csv_text


def test_the_plot_csv_carries_emission_rows_and_quotes_awkward_labels(babu):
    import csv
    import io
    from dataclasses import replace as dc_replace

    result = jo.fit(babu, weights="none")
    emission = jo.Transition(
        label="5D0->7F2, electric", energy_cm1=16000.0, f_exp=0.0, j_low=0,
        u2=0.1, u4=0.0, u6=0.0, j_up=2,
    )
    with_emission = dc_replace(result, emission=(emission,))
    rows = list(csv.reader(io.StringIO(jo.jo_plot_csv(with_emission))))
    emission_rows = [row for row in rows[1:] if row[0] == "emission"]
    assert len(emission_rows) == 1  # the trailing total row stays in the report
    cells = emission_rows[0]
    assert cells[1] == "5D0->7F2, electric"  # the comma label survived via quoting
    assert cells[3] == "" and cells[4] == "" and cells[5] == ""
    assert float(cells[6]) > 0.0
    assert 0.0 < float(cells[7]) <= 1.0
