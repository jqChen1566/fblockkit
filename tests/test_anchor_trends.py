"""Checks of the anchor trend-consistency module (analysis/anchor_trends.py;
the anchor-QA item of the ML-integration survey).

Every numeric assertion has an independent reference:

- the synthetic series are exact low-order polynomials, for which the local
  polynomial predictor is exact (a degree-two fit over five points of a
  quadratic reproduces it identically), so a clean series must come back at
  floating-point noise, and an injected scale factor has a closed-form
  relative deviation (scaling one point by 1.30 reads 3/13, by 1.08 reads
  2/27 against the covered neighbour values);
- the flagged real cases are the two defects the digitization records and
  settled on the page render -- the pre-correction ``M4`` entries of Pu3+
  and Fm3+ in ANL-89/39 -- injected back into the digitized column, whose
  corrected values are internally smooth at 0.3%;
- the loader is checked against the real fixture files (comment header, row
  count, quoted note cells, skipped fit-metadata groups) and against
  synthetic CSVs for the paths the fixtures do not exercise (an empty value
  cell, a series below the minimum length, an unknown element);
- determinism is a byte comparison of two runs, and the CSV format is
  reconstructed independently in the test (which pins the absence of a
  timestamp).
"""

from __future__ import annotations

import csv
import dataclasses
import math
import shutil
from pathlib import Path

import pytest
from pytest import approx

from fblockkit.analysis.anchor_trends import (
    FILE_PLANS,
    GRADE_OK,
    GRADE_SIGNIFICANT,
    GRADE_SUSPECT,
    FilePlan,
    TrendPoint,
    TrendSeries,
    check_literature,
    check_series,
    ion_atomic_number,
    load_series,
    read_anchor_csv,
    summarize_deviations,
    write_deviation_csv,
)

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "literature"
AN_PARAMETERS = "carnall1989_an3_lacl3_parameters.csv"
LN_PARAMETERS = "carnall1988_laf3_ln_energy_parameters.csv"


def _series(values, *, key="synthetic", file="synthetic.csv", x=None):
    """A TrendSeries over integer coordinates carrying the given values."""
    coordinates = range(len(values)) if x is None else x
    points = tuple(
        TrendPoint(label=f"pt{index}", x=float(coordinate), value=float(value))
        for index, (coordinate, value) in enumerate(zip(coordinates, values))
    )
    return TrendSeries(file=file, key=key, points=points)


def _quadratic(count=11):
    """A smooth series the local degree-two predictor reproduces exactly."""
    return [100.0 + 3.0 * index + 0.5 * index * index for index in range(count)]


def _replace(series, label, value):
    points = tuple(
        dataclasses.replace(point, value=value) if point.label == label else point
        for point in series.points
    )
    return dataclasses.replace(series, points=points)


def _by_label(deviations):
    return {deviation.label: deviation for deviation in deviations}


def _real_series(file_name, key):
    for series in load_series(FIXTURES / file_name):
        if series.key == key:
            return series
    raise AssertionError(f"series {key!r} not found in {file_name}")


# --- the trend check on synthetic data ---------------------------------------------


def test_an_exact_polynomial_series_comes_back_clean():
    deviations = check_series(_series(_quadratic()))
    assert max(deviation.rel_deviation for deviation in deviations) < 1e-9
    assert {deviation.grade for deviation in deviations} == {GRADE_OK}


def test_a_clean_non_monotone_series_comes_back_clean():
    # A downward parabola: smooth, not monotone, still an exact polynomial.
    deviations = check_series(_series([50.0 - (index - 5.0) ** 2 for index in range(11)]))
    assert {deviation.grade for deviation in deviations} == {GRADE_OK}


def test_an_injected_scale_is_graded_with_the_closed_form():
    values = _quadratic()
    series = _replace(_series(values), "pt5", values[5] * 1.30)
    deviations = _by_label(check_series(series))
    assert deviations["pt5"].grade == GRADE_SIGNIFICANT
    # |1.30 v - v| / (1.30 v) = 3/13 for the prediction v of the clean series.
    assert deviations["pt5"].rel_deviation == approx(3.0 / 13.0, rel=1e-6)
    assert all(
        deviation.grade == GRADE_OK
        for label, deviation in deviations.items()
        if label != "pt5"
    )


def test_the_two_review_bands_and_their_parameterisation():
    values = _quadratic()
    series = _replace(_series(values), "pt5", values[5] * 1.08)
    default = _by_label(check_series(series))
    assert default["pt5"].rel_deviation == approx(2.0 / 27.0, rel=1e-6)
    assert default["pt5"].grade == GRADE_SUSPECT
    strict = _by_label(check_series(series, warn=0.01, severe=0.05))
    assert strict["pt5"].grade == GRADE_SIGNIFICANT
    loose = _by_label(check_series(series, warn=0.5, severe=0.9))
    assert loose["pt5"].grade == GRADE_OK


def test_a_gross_outlier_does_not_contaminate_its_neighbours():
    values = _quadratic()
    series = _replace(_series(values), "pt5", values[5] * 3.0)
    deviations = _by_label(check_series(series))
    assert deviations["pt5"].grade == GRADE_SIGNIFICANT
    # The bisquare pass zeroes the displaced point's weight before the
    # predictions are recomputed, so its neighbours keep their clean grades.
    assert all(
        deviation.grade == GRADE_OK
        for label, deviation in deviations.items()
        if label != "pt5"
    )


def test_an_oscillating_series_flags_without_raising():
    values = [10.0 + (3.0 if index % 2 else -3.0) + 0.4 * index for index in range(11)]
    deviations = check_series(_series(values))
    assert len(deviations) == 11
    # The local trend cannot follow an alternating series; the check must
    # report that as deviations rather than fail or invent a trend.
    assert sum(1 for deviation in deviations if deviation.grade != GRADE_OK) >= 3


def test_summarize_deviations_reports_scatter_and_worst_point():
    values = _quadratic()
    series = _replace(_series(values), "pt5", values[5] * 1.30)
    summary = summarize_deviations(check_series(series))
    assert summary.n_points == 11
    assert summary.max_label == "pt5"
    assert summary.max_deviation == approx(3.0 / 13.0, rel=1e-6)
    assert summary.n_significant == 1
    assert summary.n_suspect == 0
    assert summary.median_deviation < 1e-6


def test_summarize_deviations_refuses_empty_and_mixed_collections():
    clean = check_series(_series(_quadratic()))
    with pytest.raises(ValueError, match="Next step:"):
        summarize_deviations(())
    mixed = clean + check_series(_series(_quadratic(), key="other"))
    with pytest.raises(ValueError, match="Next step:"):
        summarize_deviations(mixed)


# --- the two known defects of the ANL-89/39 table ----------------------------------


def test_the_digitized_m4_column_is_internally_smooth():
    deviations = check_series(_real_series(AN_PARAMETERS, "M4"))
    assert max(deviation.rel_deviation for deviation in deviations) < 0.01
    assert {deviation.grade for deviation in deviations} == {GRADE_OK}


def test_the_pre_correction_pu_m4_entry_is_flagged():
    series = _real_series(AN_PARAMETERS, "M4")
    deviations = _by_label(check_series(_replace(series, "Pu3+", 0.388)))
    assert deviations["Pu3+"].grade == GRADE_SUSPECT
    assert deviations["Pu3+"].rel_deviation > 0.10
    assert all(
        deviation.grade == GRADE_OK
        for label, deviation in deviations.items()
        if label != "Pu3+"
    )


def test_the_pre_correction_fm_m4_entry_is_flagged():
    series = _real_series(AN_PARAMETERS, "M4")
    deviations = _by_label(check_series(_replace(series, "Fm3+", 0.312)))
    assert deviations["Fm3+"].grade == GRADE_SIGNIFICANT
    assert deviations["Fm3+"].rel_deviation > 0.40
    assert sum(
        1 for deviation in deviations.values() if deviation.grade == GRADE_SIGNIFICANT
    ) == 1


def test_both_pre_correction_entries_together_are_both_flagged():
    series = _real_series(AN_PARAMETERS, "M4")
    injected = _replace(_replace(series, "Pu3+", 0.388), "Fm3+", 0.312)
    deviations = _by_label(check_series(injected))
    assert deviations["Pu3+"].grade == GRADE_SUSPECT
    assert deviations["Fm3+"].grade == GRADE_SIGNIFICANT
    # The two displaced points sit inside the same column; the check still
    # keeps every other grade out of the significant band.
    assert all(
        deviation.grade != GRADE_SIGNIFICANT
        for label, deviation in deviations.items()
        if label != "Fm3+"
    )


# --- the loader on the real fixture files -------------------------------------------


def test_the_loader_reads_the_comment_header_and_the_parameter_groups():
    table = read_anchor_csv(FIXTURES / LN_PARAMETERS)
    assert table.comments and "Carnall" in table.comments[0]
    assert table.columns[:4] == ("ion", "parameter", "value", "uncertainty")
    assert len(table.rows) == 297
    series = {item.key: item for item in load_series(FIXTURES / LN_PARAMETERS)}
    # Fit metadata is held out of the sweep by the registry's skip list.
    assert "n_levels" not in series
    assert "sigma" not in series
    assert len(series["F2"].points) == 11
    assert [point.label for point in series["F2"].points[:2]] == ["Pr", "Nd"]
    assert series["F2"].points[0].x == 59.0  # Pr
    assert [point.x for point in series["F2"].points] == sorted(
        point.x for point in series["F2"].points
    )
    zeta = series["zeta"]
    assert len(zeta.points) == 13
    assert zeta.points[0].label == "Ce"
    assert zeta.points[-1].label == "Yb"


def test_the_loader_carries_uncertainty_flag_and_note():
    series = {item.key: item for item in load_series(FIXTURES / AN_PARAMETERS)}
    assert "n_levels" not in series and "sigma" not in series
    f2 = {point.label: point for point in series["F2"].points}
    assert f2["Pu3+"].uncertainty == approx(89.0)
    assert f2["Am3+"].uncertainty is None
    assert f2["Am3+"].flag == "fixed"
    assert "extrapolated" in f2["Fm3+"].note
    m4 = {point.label: point for point in series["M4"].points}
    assert m4["Pu3+"].value == approx(0.338)
    assert m4["Fm3+"].value == approx(0.612)


def test_the_wide_table_loads_one_series_per_column_with_source_notes():
    series = {item.key: item for item in load_series(FIXTURES / "carnall1978_laf3_omega_lambda.csv")}
    assert set(series) == {"Omega2", "Omega4", "Omega6"}
    omega2 = {point.label: point for point in series["Omega2"].points}
    assert omega2["Pr3+"].value == approx(0.12)
    assert omega2["Pr3+"].x == 59.0
    assert omega2["Pr3+"].note == "source_ref a"
    assert omega2["Tm3+"].note == "source_ref b,g,h"


def test_empty_value_cells_are_skipped(tmp_path):
    path = tmp_path / "synthetic_parameters.csv"
    path.write_text(
        "# synthetic table\n"
        "ion,parameter,value,uncertainty,flag,note\n"
        "Pr,F2,100,1,,\n"
        "Nd,F2,,2,,\n"
        "Pm,F2,300,,,\n",
        encoding="utf-8",
    )
    series = load_series(path, plan=FilePlan("parameter"))
    f2 = {item.key: item for item in series}["F2"]
    assert [point.label for point in f2.points] == ["Pr", "Pm"]


def test_the_loader_refuses_an_unknown_element(tmp_path):
    path = tmp_path / "unknown_ion.csv"
    path.write_text("ion,parameter,value\nZz3+,F2,1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Next step:"):
        load_series(path, plan=FilePlan("parameter"))


def test_ion_labels_map_to_atomic_numbers():
    assert ion_atomic_number("Pr") == 59
    assert ion_atomic_number("Pr3+") == 59
    assert ion_atomic_number("Eu3+(aq)") == 63
    assert ion_atomic_number("U3+") == 92
    assert ion_atomic_number("No3+") == 102
    with pytest.raises(ValueError, match="Next step:"):
        ion_atomic_number("Zz3+")
    with pytest.raises(ValueError, match="Next step:"):
        ion_atomic_number(59)


def test_load_series_refuses_unknown_names_and_honours_none_plans(tmp_path):
    unknown = tmp_path / "unknown_table.csv"
    unknown.write_text("ion,parameter,value\nPr,X,1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Next step:"):
        load_series(unknown)
    lifetimes = FIXTURES / "carnall1978_laf3_lifetimes.csv"
    assert FILE_PLANS[lifetimes.name].kind == "none"
    assert load_series(lifetimes) == ()
    assert load_series(lifetimes, plan=FilePlan("none", reason="held out")) == ()


def test_the_registry_entries_are_complete():
    assert len(FILE_PLANS) == 13
    for name, plan in FILE_PLANS.items():
        assert name.endswith(".csv")
        if plan.kind == "none":
            assert plan.reason and not plan.columns
        elif plan.kind == "wide":
            assert plan.columns
        else:
            assert plan.kind == "parameter"


# --- reading errors -----------------------------------------------------------------


def test_read_anchor_csv_refuses_a_missing_file(tmp_path):
    with pytest.raises(ValueError, match="Next step:"):
        read_anchor_csv(tmp_path / "absent.csv")


def test_read_anchor_csv_refuses_a_ragged_row(tmp_path):
    path = tmp_path / "ragged.csv"
    path.write_text("# c\nion,parameter,value\nPr,F2,1\nNd,F2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Next step:"):
        read_anchor_csv(path)


def test_a_non_numeric_cell_is_refused(tmp_path):
    path = tmp_path / "bad_value.csv"
    path.write_text("# c\nion,parameter,value\nPr,F2,12x\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Next step:"):
        load_series(path, plan=FilePlan("parameter"))


# --- refusals and validation --------------------------------------------------------


def test_a_series_below_the_minimum_length_is_refused():
    with pytest.raises(ValueError, match="Next step:"):
        check_series(_series([1.0, 2.0, 3.0, 4.0]))


def test_the_controls_are_validated():
    series = _series(_quadratic())
    for controls in (
        {"warn": 0.3, "severe": 0.2},
        {"warn": -0.1},
        {"span": 0.0},
        {"span": 1.5},
        {"degree": 0},
        {"min_neighbors": 2},
        {"min_points": 3},
        {"robust_iterations": 0},
        {"degree": 2.5},
        {"min_points": "5"},
        {"warn": "five"},
    ):
        with pytest.raises(ValueError, match="Next step:"):
            check_series(series, **controls)


def test_non_finite_values_and_a_flat_coordinate_are_refused():
    series = _series(_quadratic())
    broken = dataclasses.replace(
        series,
        points=series.points[:5]
        + (dataclasses.replace(series.points[5], value=float("nan")),)
        + series.points[6:],
    )
    with pytest.raises(ValueError, match="Next step:"):
        check_series(broken)
    flat = TrendSeries(
        file="flat.csv",
        key="K",
        points=tuple(
            TrendPoint(label=f"pt{index}", x=3.0, value=float(index)) for index in range(6)
        ),
    )
    with pytest.raises(ValueError, match="Next step:"):
        check_series(flat)


def test_the_sweep_records_a_short_series_as_a_refusal(tmp_path):
    scratch = tmp_path / "literature"
    shutil.copytree(FIXTURES, scratch)
    # Shorten one registered table (a copy, never the fixture itself).
    (scratch / "carnall1978_laf3_omega_lambda.csv").write_text(
        "# shortened copy for the refusal path\n"
        "ion,Omega2,Omega4,Omega6,source_ref\n"
        "Pr3+,0.12,1.77,4.78,a\n"
        "Nd3+,0.35,2.57,2.50,a\n"
        "Pm3+,0.5,1.9,2.2,b\n",
        encoding="utf-8",
    )
    sweep = check_literature(scratch)
    refusals = {key: reason for _file, key, reason in sweep.refusals()}
    assert set(refusals) == {"Omega2", "Omega4", "Omega6"}
    assert "min_points=5" in refusals["Omega2"]
    assert len(sweep.files) == 13
    assert not any(
        outcome.deviations
        for file_outcome in sweep.files
        if file_outcome.file == "carnall1978_laf3_omega_lambda.csv"
        for outcome in file_outcome.outcomes
    )


def test_the_sweep_refuses_a_missing_registered_file(tmp_path):
    scratch = tmp_path / "literature"
    shutil.copytree(FIXTURES, scratch)
    (scratch / "carnall1978_laf3_lifetimes.csv").unlink()
    with pytest.raises(ValueError, match="Next step:"):
        check_literature(scratch)


# --- the directory sweep and the tidy CSV -------------------------------------------


def test_the_sweep_covers_the_thirteen_registered_tables():
    sweep = check_literature(FIXTURES)
    assert [file_outcome.file for file_outcome in sweep.files] == list(FILE_PLANS)
    for file_outcome in sweep.files:
        if file_outcome.not_analysed:
            assert file_outcome.not_analysed and not file_outcome.outcomes
        else:
            assert file_outcome.outcomes
    for deviation in sweep.deviations():
        assert math.isfinite(deviation.rel_deviation)
        assert 0.0 <= deviation.rel_deviation <= 2.0
        assert deviation.grade in {GRADE_OK, GRADE_SUSPECT, GRADE_SIGNIFICANT}


def test_the_smooth_parameter_families_carry_no_significant_flags():
    # The Racah parameters, the spin-orbit constant and the M family of the
    # actinide table are the smooth columns the check exists for: the
    # digitized values must survive it with a scatter well below the review
    # threshold and no significant grade anywhere.  At most one borderline
    # flag is allowed per family: an endpoint sits on the one-sided
    # prediction of the first (or last) member, which is where a few percent
    # of model bias lives even on a clean column.
    smooth_families = {
        (LN_PARAMETERS, "F2"),
        (LN_PARAMETERS, "F4"),
        (LN_PARAMETERS, "F6"),
        (LN_PARAMETERS, "zeta"),
        (AN_PARAMETERS, "F2"),
        (AN_PARAMETERS, "F4"),
        (AN_PARAMETERS, "F6"),
        (AN_PARAMETERS, "zeta"),
        (AN_PARAMETERS, "M0"),
        (AN_PARAMETERS, "M2"),
        (AN_PARAMETERS, "M4"),
    }
    found = {}
    for file_outcome in check_literature(FIXTURES).files:
        for outcome in file_outcome.outcomes:
            if outcome.deviations:
                found[(file_outcome.file, outcome.series.key)] = outcome.deviations
    for family in smooth_families:
        assert family in found, family
        summary = summarize_deviations(found[family])
        assert summary.n_significant == 0, family
        assert summary.n_suspect <= 1, family
        assert summary.median_deviation < 0.01, family


def test_the_sweep_is_deterministic(tmp_path):
    first = check_literature(FIXTURES)
    second = check_literature(FIXTURES)
    first_csv = tmp_path / "first.csv"
    second_csv = tmp_path / "second.csv"
    write_deviation_csv(first_csv, first.deviations())
    write_deviation_csv(second_csv, second.deviations())
    assert first_csv.read_bytes() == second_csv.read_bytes()


def test_the_deviation_csv_matches_an_independent_reconstruction(tmp_path):
    series = _replace(_series([1.0, 3.0, 5.0, 7.0, 9.0], key="K"), "pt3", 9.9)
    deviations = check_series(series)
    path = write_deviation_csv(tmp_path / "deviations.csv", deviations)
    expected = [
        "file,series,point,x,value,prediction,residual,rel_deviation,grade,"
        "uncertainty,flag,note"
    ]
    for deviation in deviations:
        expected.append(
            f"{deviation.file},{deviation.series},{deviation.label},"
            f"{deviation.x:.10g},{deviation.value:.10g},{deviation.prediction:.10g},"
            f"{deviation.residual:.10g},{deviation.rel_deviation:.10g},"
            f"{deviation.grade},,,"
        )
    text = path.read_text(encoding="utf-8")
    assert text == "\n".join(expected) + "\n"


def test_the_deviation_csv_reads_back_with_the_source_columns(tmp_path):
    deviations = check_series(_real_series(LN_PARAMETERS, "F2"))
    path = write_deviation_csv(tmp_path / "f2.csv", deviations)
    rows = list(csv.DictReader(path.read_text(encoding="utf-8").splitlines()))
    assert len(rows) == len(deviations)
    for row, deviation in zip(rows, deviations):
        assert row["point"] == deviation.label
        assert float(row["value"]) == approx(deviation.value)
        assert float(row["prediction"]) == approx(deviation.prediction)
        assert row["grade"] == deviation.grade
