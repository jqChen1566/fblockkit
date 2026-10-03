"""Anchor trend consistency: a robust leave-one-out sweep of the digitized
literature tables under ``fixtures/literature``.

Several of the anchor CSVs are per-ion parameter tables: one value per
trivalent lanthanide or actinide, in an order fixed by the periodic table
(Carnall 1988 Table 4 for Ln3+:LaF3, Carnall 1989 Table 3 for An3+:LaCl3,
Carnall 1978 Appendix XIV Table 1 for the Judd-Ofelt parameters).  A cell of
such a table that disagrees with its neighbours in the ion sequence is the
kind of fault the digitization's manual spot-checks caught once already: the
pre-correction ``M4`` entries of Pu3+ and Fm3+ in ANL-89/39, settled on the
page render with the invariant ratio ``M4/M0 = 0.385``.  This module turns
that check into a systematic, deterministic sweep:

- :func:`read_anchor_csv` reads a fixture CSV (its ``#`` comment header and
  its data rows);
- :func:`load_series` organises one file into :class:`TrendSeries` following
  :data:`FILE_PLANS`, the per-file registry that also records, for every table
  this check does not apply to, the reason why;
- :func:`check_series` predicts every point from the other points alone
  (leave-one-out) and reports the relative deviation;
- :func:`check_literature` sweeps the whole anchor directory and
  :func:`write_deviation_csv` writes the tidy deviation table.

Method and its limits (the method is a choice, made for n = 8-13 points per
series, and it has known failure modes):

- the predictor is a local polynomial of degree two with tricube weights and a
  nearest-neighbour bandwidth (``neighbors`` nearest other points, the
  farthest at weight zero).  Degree two follows the mild curvature of these
  series without the large endpoint bias a locally linear fit would carry
  (one-sided extrapolation is biased by ``O(h^2 y'')``, which for a convex
  series of length 11-14 is several percent), and the bandwidth keeps the fit
  local enough to follow real structure;
- bisquare reweighting (tuning 4.685, the standard 95%-efficiency value)
  downweights other points whose own leave-one-out residual is large, so a
  single displaced cell cannot drag the prediction of its neighbours and
  produce a second spurious flag.  The pass is iterated
  (``robust_iterations``, default 2): with one displaced point in a short
  series the first residual scale is itself contaminated by that point (its
  neighbours' predictions are pulled with it), and the reweighting has to
  re-run on the improved residuals before the displaced point's weight
  reaches zero;
- leave-one-out is what makes a defect visible at all: a point's own value
  takes no part in its prediction.  A *pair* of consistently displaced cells
  still hides from this check -- it compares each cell with the rest of the
  column, not with the source page;
- the deviation is ``|value - prediction|`` relative to
  ``max(|value|, |prediction|)``, so it is bounded by 2 and never divides by
  zero.  A value above 1 means the prediction lies on the other side of zero
  from the point (a series that crosses zero, or a boundary extrapolation);
  for a point whose value and prediction are of comparable size the metric
  agrees with the ordinary relative error to well below the review
  thresholds;
- the grades are informational review thresholds (defaults: above 5%
  "suspect", above 20% "significant"), parameterised per call.  **A grade is
  not a verdict on the data.**  These series carry genuine non-smooth
  structure -- half-filled-shell effects, the anomalous Sm/Eu Judd-Ofelt
  parameters, parameters fixed or ratio-constrained in the fit, columns
  interpolated in the source, values compiled from different literature
  sources along one row of the table -- and the trend model cannot tell that
  structure apart from a mistyped cell.  A flagged point is a review
  candidate; the dispositions (model limitation, source physics, residual
  digitization doubt) belong in the report beside the numbers, not in the
  grade;
- the first and last points of a series are predicted one-sided, so their
  deviations carry the local model's extrapolation bias on top of any
  mismatch; a borderline flag at an endpoint is a prompt to re-read the page,
  not a finding;
- series of n = 8-13 also mean the check has little power against small
  displacements: a 3% misread cell is indistinguishable from the fit's own
  scatter, and only a gross displacement (the two known cases read 13% and
  48% under this check) stands out reliably.  :func:`summarize_deviations`
  reports a series' own scatter, which is the scale a flagged deviation has
  to be read against: where that scatter is itself several percent, a single
  flag carries no information.

Determinism: the same input gives byte-identical output -- no timestamp, no
randomness, one fixed fitting order and no floating-point accumulation across
calls.  No CSV is ever written back: this module reads the anchors and writes
only its own deviation table.
"""

from __future__ import annotations

import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from .plot_csv import csv_field, plot_number

__all__ = [
    "FILE_PLANS",
    "GRADE_OK",
    "GRADE_SIGNIFICANT",
    "GRADE_SUSPECT",
    "ION_ATOMIC_NUMBERS",
    "AnchorTable",
    "Deviation",
    "DeviationSummary",
    "FileOutcome",
    "FilePlan",
    "SeriesOutcome",
    "Sweep",
    "TrendPoint",
    "TrendSeries",
    "check_literature",
    "check_series",
    "ion_atomic_number",
    "load_series",
    "read_anchor_csv",
    "summarize_deviations",
    "write_deviation_csv",
]

# Review grades (informational; see the module docstring for what they are not).
GRADE_OK = "ok"
GRADE_SUSPECT = "suspect"
GRADE_SIGNIFICANT = "significant"

# Default review thresholds for the bounded relative deviation.
DEFAULT_WARN = 0.05
DEFAULT_SEVERE = 0.20

# Local-fit controls: a degree-two polynomial over the nearest neighbours with
# tricube weights, and iterated bisquare reweighting at the standard tuning.
_DEFAULT_SPAN = 0.6
_DEFAULT_MIN_NEIGHBORS = 6
_DEFAULT_DEGREE = 2
_DEFAULT_MIN_POINTS = 5
_DEFAULT_ROBUST_ITERATIONS = 2
_BISQUARE_TUNING = 4.685
_MAD_TO_SIGMA = 0.6745
# The bisquare scale is floored at this fraction of the series magnitude, so
# that a series whose residuals are pure rounding noise (a synthetic exact
# polynomial, say) keeps unit weights instead of scaling 0/0.
_ROBUST_SCALE_FLOOR = 1e-9

# The ions the sweep places on the atomic-number axis (the f-block sequences
# the anchor tables cover).  Z = 57 + index and Z = 89 + index respectively.
_LANTHANIDE_SYMBOLS = (
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
    "Tm", "Yb", "Lu",
)
_ACTINIDE_SYMBOLS = (
    "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr",
)
ION_ATOMIC_NUMBERS: Mapping[str, int] = {
    **{symbol: 57 + index for index, symbol in enumerate(_LANTHANIDE_SYMBOLS)},
    **{symbol: 89 + index for index, symbol in enumerate(_ACTINIDE_SYMBOLS)},
}

_ION_LABEL = re.compile(r"\s*([A-Z][a-z]?)")


# --- the data model ----------------------------------------------------------------


@dataclass(frozen=True)
class TrendPoint:
    """One cell of a series: a labelled value at an ordering coordinate.

    ``label`` is the point as the source prints it (an ion symbol usually),
    ``x`` its ordering coordinate (the atomic number for the ion series of
    this module, but any monotone coordinate is accepted), ``value`` the
    tabulated number.  ``uncertainty``, ``flag`` and ``note`` carry the
    source's own columns when the table has them (``flag`` is "fixed" or
    "ratio-constrained" in the parameter tables; ``note`` records the
    source's annotations such as interpolation or extrapolation).
    """

    label: str
    x: float
    value: float
    uncertainty: float | None = None
    flag: str = ""
    note: str = ""


@dataclass(frozen=True)
class TrendSeries:
    """One series of the sweep: a named value sequence of one anchor file."""

    file: str
    key: str
    points: tuple[TrendPoint, ...]


@dataclass(frozen=True)
class Deviation:
    """The trend check's result for one point (see the module docstring)."""

    file: str
    series: str
    label: str
    x: float
    value: float
    prediction: float
    residual: float
    rel_deviation: float
    grade: str
    uncertainty: float | None = None
    flag: str = ""
    note: str = ""


@dataclass(frozen=True)
class DeviationSummary:
    """The shape of one series' deviations (see :func:`summarize_deviations`)."""

    file: str
    series: str
    n_points: int
    median_deviation: float
    max_deviation: float
    max_label: str
    n_suspect: int
    n_significant: int


@dataclass(frozen=True)
class AnchorTable:
    """A fixture CSV as read: comment header, column names and string rows."""

    path: Path
    comments: tuple[str, ...]
    columns: tuple[str, ...]
    rows: tuple[Mapping[str, str], ...]
    line_numbers: tuple[int, ...]


@dataclass(frozen=True)
class FilePlan:
    """What the sweep does with one registered anchor file.

    ``kind`` is ``"parameter"`` (a long ``ion,parameter,value`` table grouped
    by parameter), ``"wide"`` (one series per listed value column) or
    ``"none"`` (the file defines no series the trend model applies to, and
    ``reason`` says why).  ``skip`` names parameter groups held out of the
    sweep, and ``columns`` the value columns of a wide table.
    """

    kind: str
    columns: tuple[str, ...] = ()
    skip: tuple[str, ...] = ()
    reason: str = ""


# The registry of the thirteen digitized spectroscopic anchor tables.  The
# order is the order of the sweep and of the report.  A file listed as "none"
# is not an oversight: the trend model's smoothness premise does not hold for
# its rows, and the reason string is what the sweep reports instead.
FILE_PLANS: Mapping[str, FilePlan] = {
    "carnall1968_eu3_aquo_levels_U.csv": FilePlan(
        "none",
        reason=(
            "single ion (Eu3+ aquo): the rows are one J-level ladder with the "
            "U(k) reduced matrix elements, not a scalar series over an ion "
            "coordinate"
        ),
    ),
    "carnall1968_eu3_aquo_parameters.csv": FilePlan(
        "none",
        reason=(
            "single ion in two environments (Eu3+ aquo and Eu3+:LaCl3); the "
            "parameters are not ordered by any coordinate and most appear "
            "once per environment"
        ),
    ),
    "carnall1988_laf3_eu3_cf_levels.csv": FilePlan(
        "none",
        reason=(
            "single ion (Eu3+): the crystal-field components of one J "
            "multiplet form a field-determined ladder, not a smooth sequence"
        ),
    ),
    "carnall1988_laf3_ln_energy_parameters.csv": FilePlan(
        "parameter",
        skip=("n_levels", "sigma"),
        reason="fit metadata (n_levels, sigma) excluded: no trend premise",
    ),
    "gorller1991_eu3_md_strengths.csv": FilePlan(
        "none",
        reason=(
            "single ion; the groups are (symmetry, transition, kind) branches "
            "of one to three component strengths"
        ),
    ),
    "gorller1991_eu3_md_expt_vs_theor.csv": FilePlan(
        "none",
        reason=(
            "experiment-versus-theory comparison: the rows are different "
            "samples of one transition, with no ordering coordinate"
        ),
    ),
    "carnall1978_laf3_omega_lambda.csv": FilePlan(
        "wide", columns=("Omega2", "Omega4", "Omega6")
    ),
    "carnall1978_laf3_eu3_U2.csv": FilePlan(
        "none",
        reason=(
            "single ion; squared unit-tensor matrix elements between level "
            "pairs, sparse by the selection rules (many exact zeros)"
        ),
    ),
    "carnall1989_an3_lacl3_parameters.csv": FilePlan(
        "parameter",
        skip=("n_levels", "sigma"),
        reason="fit metadata (n_levels, sigma) excluded: no trend premise",
    ),
    "carnall1968_ln3_aquo_levels_U.csv": FilePlan(
        "none",
        reason=(
            "level ladders of eight different aquo ions; the J multiplets "
            "differ between ions, so no cross-ion scalar series exists"
        ),
    ),
    "carnall1988_laf3_ho_tb_gd_cf_levels.csv": FilePlan(
        "none",
        reason=(
            "level ladders of three ions (Ho, Tb, Gd): same reason as the "
            "other crystal-field level tables"
        ),
    ),
    "carnall1988_laf3_remaining_cf_levels.csv": FilePlan(
        "none",
        reason=(
            "level ladders of nine ions: same reason as the other "
            "crystal-field level tables"
        ),
    ),
    "carnall1978_laf3_lifetimes.csv": FilePlan(
        "none",
        reason=(
            "the tabulated states differ between ions, so the radiative "
            "lifetimes are state-dominated, not a smooth function of Z; a "
            "trend check would report model inapplicability as flags"
        ),
    ),
}


@dataclass(frozen=True)
class SeriesOutcome:
    """One series of a sweep: checked (deviations) or refused (reason)."""

    series: TrendSeries
    deviations: tuple[Deviation, ...] = ()
    refusal: str = ""


@dataclass(frozen=True)
class FileOutcome:
    """One anchor file of a sweep."""

    file: str
    outcomes: tuple[SeriesOutcome, ...] = ()
    not_analysed: str = ""


@dataclass(frozen=True)
class Sweep:
    """The whole-directory result of :func:`check_literature`."""

    root: Path
    files: tuple[FileOutcome, ...]

    def deviations(self) -> tuple[Deviation, ...]:
        """Every computed deviation, in sweep order."""
        return tuple(
            deviation
            for file_outcome in self.files
            for series_outcome in file_outcome.outcomes
            for deviation in series_outcome.deviations
        )

    def refusals(self) -> tuple[tuple[str, str, str], ...]:
        """``(file, series, reason)`` of every series the check refused."""
        return tuple(
            (file_outcome.file, series_outcome.series.key, series_outcome.refusal)
            for file_outcome in self.files
            for series_outcome in file_outcome.outcomes
            if series_outcome.refusal
        )


# --- reading the fixture CSVs ------------------------------------------------------


def read_anchor_csv(path) -> AnchorTable:
    """Read one anchor CSV: its ``#`` comment header and its data rows.

    Comment lines (the provenance block the digitization writes above every
    table) are returned stripped of the leading ``#``; the first non-comment
    line is the column header, every later non-empty line a row, cells
    stripped of surrounding whitespace (the source prints alignment spaces,
    and the empty cells of an ellipsis or dash stay empty).  Quote characters
    are honoured, so a note cell may contain commas.  Each row keeps the file
    line it came from, for error messages that can be checked against the
    page.
    """
    target = Path(path)
    if not target.is_file():
        raise ValueError(
            f"anchor CSV not found: {target}. Next step: pass a path under "
            f"fixtures/literature (that directory's README lists the tables "
            f"and their sources)."
        )
    comments: list[str] = []
    header: tuple[str, ...] | None = None
    rows: list[Mapping[str, str]] = []
    line_numbers: list[int] = []
    for number, raw in enumerate(target.read_text(encoding="utf-8").splitlines(), start=1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            comments.append(line[1:].strip())
            continue
        cells = [cell.strip() for cell in next(csv.reader([line]))]
        if header is None:
            header = tuple(cells)
            if len(set(header)) != len(header):
                raise ValueError(
                    f"{target.name}:{number}: the header repeats a column "
                    f"name ({header}). Next step: make the column names "
                    f"unique so the rows can be addressed by name."
                )
            continue
        if len(cells) != len(header):
            raise ValueError(
                f"{target.name}:{number}: the row has {len(cells)} cells for "
                f"{len(header)} columns. Next step: fix the row's quoting or "
                f"its cell count so that it matches the header line."
            )
        rows.append(dict(zip(header, cells)))
        line_numbers.append(number)
    if header is None:
        raise ValueError(
            f"{target.name} carries no header line (only comments or "
            f"nothing). Next step: check that the digitized table was copied "
            f"with its column header."
        )
    return AnchorTable(
        path=target,
        comments=tuple(comments),
        columns=header,
        rows=tuple(rows),
        line_numbers=tuple(line_numbers),
    )


def ion_atomic_number(label: str) -> int:
    """The atomic number of an ion label such as ``"Pr"``, ``"U3+"`` or
    ``"Eu3+(aq)"`` (the leading element symbol is read; unknown symbols are
    refused with a Next-step message)."""
    if not isinstance(label, str):
        raise ValueError(
            f"an ion label must be a string, got {label!r}. Next step: pass "
            f"the label as the source prints it (for example 'Pr3+')."
        )
    match = _ION_LABEL.match(label)
    if match is None:
        raise ValueError(
            f"cannot read an element symbol from the ion label {label!r}. "
            f"Next step: pass a label that starts with the element symbol "
            f"(for example 'Pr3+')."
        )
    symbol = match.group(1)
    number = ION_ATOMIC_NUMBERS.get(symbol)
    if number is None:
        raise ValueError(
            f"the ion label {label!r} reads as element {symbol!r}, which is "
            f"not in this module's lanthanide/actinide table. Next step: "
            f"extend ION_ATOMIC_NUMBERS with the element, or pass the "
            f"ordering coordinate explicitly."
        )
    return number


def load_series(path, *, plan: FilePlan | None = None) -> tuple[TrendSeries, ...]:
    """Organise one anchor file into trend series, following :data:`FILE_PLANS`.

    ``plan`` defaults to the registry entry for the file's name; a file with
    no entry is refused rather than guessed at.  Rows whose value cell is
    empty are skipped (that is how the sources' ellipsis and dash cells are
    stored); the points of every series are ordered by ``(x, label)``, so the
    output does not depend on the row order of the file.
    """
    target = Path(path)
    if plan is None:
        plan = FILE_PLANS.get(target.name)
    if plan is None:
        raise ValueError(
            f"no series plan is registered for {target.name!r}. Next step: "
            f"add the file to FILE_PLANS with its table kind, or pass "
            f"plan=FilePlan(...) explicitly."
        )
    if plan.kind == "none":
        return ()
    table = read_anchor_csv(target)
    if plan.kind == "parameter":
        return _parameter_series(table, skip=plan.skip)
    if plan.kind == "wide":
        return _wide_series(table, columns=plan.columns)
    raise ValueError(
        f"unknown plan kind {plan.kind!r} for {target.name!r}. Next step: use "
        f"'parameter', 'wide' or 'none' in the FILE_PLANS entry."
    )


def _parameter_series(table: AnchorTable, *, skip: Sequence[str]) -> tuple[TrendSeries, ...]:
    """Long-format tables: one series per parameter, over the ion sequence."""
    _require_columns(table, ("ion", "parameter", "value"))
    groups: dict[str, list[TrendPoint]] = {}
    for row, line in zip(table.rows, table.line_numbers):
        parameter = row["parameter"]
        if parameter in skip:
            continue
        value_text = row["value"]
        if not value_text:
            continue  # an empty cell is a missing value; the point is dropped
        label = row["ion"]
        if not label:
            raise ValueError(
                f"{table.path.name}:{line}: the ion cell is empty. Next step: "
                f"fill the ion label or drop the row from the digitized table."
            )
        groups.setdefault(parameter, []).append(
            TrendPoint(
                label=label,
                x=ion_atomic_number(label),
                value=_float_cell(value_text, table, "value", line),
                uncertainty=_float_or_none(row.get("uncertainty", ""), table, "uncertainty", line),
                flag=row.get("flag", ""),
                note=row.get("note", ""),
            )
        )
    return tuple(_ordered_series(table, key, points) for key, points in groups.items())


def _wide_series(table: AnchorTable, *, columns: Sequence[str]) -> tuple[TrendSeries, ...]:
    """Wide tables: one series per named value column, over the ion sequence."""
    _require_columns(table, ("ion", *columns))
    groups: dict[str, list[TrendPoint]] = {name: [] for name in columns}
    for row, line in zip(table.rows, table.line_numbers):
        label = row["ion"]
        if not label:
            raise ValueError(
                f"{table.path.name}:{line}: the ion cell is empty. Next step: "
                f"fill the ion label or drop the row from the digitized table."
            )
        source_ref = row.get("source_ref", "")
        note = f"source_ref {source_ref}" if source_ref else ""
        for name in columns:
            value_text = row.get(name, "")
            if not value_text:
                continue  # an empty cell is a missing value; the point is dropped
            groups[name].append(
                TrendPoint(
                    label=label,
                    x=ion_atomic_number(label),
                    value=_float_cell(value_text, table, name, line),
                    note=note,
                )
            )
    return tuple(_ordered_series(table, key, points) for key, points in groups.items())


def _ordered_series(table: AnchorTable, key: str, points: Sequence[TrendPoint]) -> TrendSeries:
    return TrendSeries(
        file=table.path.name,
        key=key,
        points=tuple(sorted(points, key=lambda point: (point.x, point.label))),
    )


def _require_columns(table: AnchorTable, required: Sequence[str]) -> None:
    missing = [name for name in required if name not in table.columns]
    if missing:
        raise ValueError(
            f"{table.path.name} has no column {missing[0]!r}; its columns are "
            f"{list(table.columns)}. Next step: check the digitized header or "
            f"the FILE_PLANS entry for this file."
        )


def _float_cell(text: str, table: AnchorTable, column: str, line: int) -> float:
    try:
        return float(text)
    except ValueError:
        raise ValueError(
            f"{table.path.name}:{line}: the {column} cell reads {text!r}, "
            f"which is not a number. Next step: check that cell against the "
            f"page render and fix the digitized value."
        ) from None


def _float_or_none(text: str, table: AnchorTable, column: str, line: int) -> float | None:
    if not text:
        return None
    return _float_cell(text, table, column, line)


# --- the trend check ---------------------------------------------------------------


def check_series(
    series: TrendSeries,
    *,
    warn: float = DEFAULT_WARN,
    severe: float = DEFAULT_SEVERE,
    span: float = _DEFAULT_SPAN,
    min_neighbors: int = _DEFAULT_MIN_NEIGHBORS,
    degree: int = _DEFAULT_DEGREE,
    min_points: int = _DEFAULT_MIN_POINTS,
    robust: bool = True,
    robust_iterations: int = _DEFAULT_ROBUST_ITERATIONS,
) -> tuple[Deviation, ...]:
    """Predict every point of ``series`` from the other points (leave-one-out)
    and report the relative deviation of every point.

    The predictor is the local polynomial of ``degree`` over the
    ``neighbors`` nearest other points in the scaled ordering coordinate,
    with tricube weights and the farthest neighbour at weight zero
    (``neighbors`` follows from ``span`` as a fraction of the series and the
    ``min_neighbors`` floor).  With ``robust`` a bisquare pass downweights
    points whose own leave-one-out residual is large before the predictions
    are recomputed, and the pass is repeated ``robust_iterations`` times (the
    residual scale of the first pass is still contaminated by a displaced
    point, so one pass alone does not clear it); a residual scale that is
    numerically zero against the series magnitude is floored at that
    magnitude, so a perfectly smooth series keeps unit weights instead of a
    0/0.  ``min_points`` is the shortest series the check accepts: below it
    there are too few distinct neighbours for a leave-one-out prediction to
    mean anything, and the function raises rather than returning noise.

    ``warn`` and ``severe`` grade the bounded relative deviation
    (``|value - prediction| / max(|value|, |prediction|)``); they are review
    thresholds, not error criteria -- see the module docstring.  The
    deviations come back in the series' point order.
    """
    _validate_controls(warn=warn, severe=severe, span=span, min_neighbors=min_neighbors,
                       degree=degree, min_points=min_points,
                       robust_iterations=robust_iterations)
    warn = float(warn)
    severe = float(severe)
    span = float(span)
    if len(series.points) < min_points:
        raise ValueError(
            f"series {series.key!r} of {series.file!r} has "
            f"{len(series.points)} points; the trend check needs at least "
            f"min_points={min_points} (with fewer, leave-one-out prediction "
            f"has no independent neighbours to stand on). Next step: check "
            f"the series with a lower min_points, or leave it out of the "
            f"sweep and say so in the report."
        )
    x = np.array([point.x for point in series.points], dtype=float)
    y = np.array([point.value for point in series.points], dtype=float)
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        bad = series.points[int(np.argmin(np.isfinite(x) & np.isfinite(y)))]
        raise ValueError(
            f"series {series.key!r} of {series.file!r} carries a non-finite "
            f"coordinate or value (first at {bad.label!r}). Next step: drop "
            f"that point before checking, or fix the digitized cell."
        )
    low, high = float(x.min()), float(x.max())
    if high <= low:
        raise ValueError(
            f"every point of series {series.key!r} of {series.file!r} shares "
            f"one coordinate ({low:g}), so no trend can be estimated. Next "
            f"step: pass an ordering coordinate that varies, such as the "
            f"atomic number of each ion."
        )
    scaled = (x - low) / (high - low)
    neighbors = min(x.size - 1, max(min_neighbors, math.ceil(span * (x.size - 1))))
    predictions = _leave_one_out_predictions(
        scaled,
        y,
        degree=degree,
        neighbors=neighbors,
        robust=robust,
        iterations=robust_iterations,
    )
    deviations: list[Deviation] = []
    for point, prediction in zip(series.points, predictions):
        residual = float(point.value - prediction)
        denominator = max(abs(float(point.value)), abs(float(prediction)))
        relative = abs(residual) / denominator if denominator > 0.0 else 0.0
        deviations.append(
            Deviation(
                file=series.file,
                series=series.key,
                label=point.label,
                x=point.x,
                value=point.value,
                prediction=float(prediction),
                residual=residual,
                rel_deviation=relative,
                grade=_grade(relative, warn, severe),
                uncertainty=point.uncertainty,
                flag=point.flag,
                note=point.note,
            )
        )
    return tuple(deviations)


def _validate_controls(
    *, warn: float, severe: float, span: float, min_neighbors: int, degree: int,
    min_points: int, robust_iterations: int,
) -> None:
    for name, value in (("warn", warn), ("severe", severe)):
        try:
            value = float(value)
        except (TypeError, ValueError):
            raise ValueError(
                f"{name} must be a number, got {value!r}. Next step: pass the "
                f"review thresholds as fractions (for example warn=0.05, "
                f"severe=0.20)."
            ) from None
        if not math.isfinite(value) or value < 0.0:
            raise ValueError(
                f"{name} must be a finite non-negative fraction, got "
                f"{value!r}. Next step: pass the review thresholds as "
                f"fractions (for example warn=0.05, severe=0.20)."
            )
    if not float(warn) < float(severe):
        raise ValueError(
            f"the review thresholds must satisfy warn < severe, got "
            f"warn={warn!r} and severe={severe!r}. Next step: pass for "
            f"example warn=0.05, severe=0.20; the two grades are "
            f"informational review bands, not error criteria."
        )
    try:
        span = float(span)
    except (TypeError, ValueError):
        raise ValueError(
            f"span must be a number, got {span!r}. Next step: pass the "
            f"neighbourhood as a fraction of the series (0 < span <= 1, "
            f"default 0.6)."
        ) from None
    if not 0.0 < span <= 1.0:
        raise ValueError(
            f"span must lie in (0, 1], got {span!r}. Next step: pass the "
            f"neighbourhood as a fraction of the series (default 0.6)."
        )
    for name, value, hint in (
        ("degree", degree, "pass the local polynomial degree (default 2)."),
        ("min_neighbors", min_neighbors, "pass the smallest neighbour count a local fit may use (default 6)."),
        ("min_points", min_points, "pass the shortest series the check may run on (default 5)."),
        ("robust_iterations", robust_iterations, "pass the number of bisquare reweighting passes (default 2)."),
    ):
        if isinstance(value, bool) or not isinstance(value, int):
            raise ValueError(
                f"{name} must be an integer, got {value!r}. Next step: {hint}"
            )
    if degree < 1:
        raise ValueError(
            f"degree must be at least 1, got {degree!r}. Next step: pass the "
            f"local polynomial degree (default 2)."
        )
    if min_neighbors < 3:
        raise ValueError(
            f"min_neighbors must be at least 3, got {min_neighbors!r}. Next "
            f"step: pass the smallest neighbour count a local fit may use "
            f"(default 6; three leave a two-point fit after the farthest "
            f"weight vanishes)."
        )
    if min_points < 4:
        raise ValueError(
            f"min_points must be at least 4, got {min_points!r}. Next step: "
            f"pass the shortest series the check may run on (default 5)."
        )
    if robust_iterations < 1:
        raise ValueError(
            f"robust_iterations must be at least 1, got {robust_iterations!r}. "
            f"Next step: pass the number of bisquare reweighting passes "
            f"(default 2), or robust=False to switch the reweighting off "
            f"entirely."
        )


def _leave_one_out_predictions(
    scaled: np.ndarray,
    values: np.ndarray,
    *,
    degree: int,
    neighbors: int,
    robust: bool,
    iterations: int,
) -> np.ndarray:
    """The leave-one-out predictions, with iterated bisquare reweighting."""
    unit = np.ones(scaled.size)
    predictions = _all_local_predictions(scaled, values, unit, degree=degree, neighbors=neighbors)
    if not robust:
        return predictions
    extra = unit
    for _ in range(iterations):
        residual = values - predictions
        scale = float(np.median(np.abs(residual))) / _MAD_TO_SIGMA
        reference = float(np.median(np.abs(values)))
        if reference == 0.0:
            reference = float(np.max(np.abs(values)))
        # The MAD scale collapses when the series is otherwise exact (a
        # synthetic polynomial, say) and only one point is displaced; the
        # floor against the series magnitude keeps such a point downweighted
        # instead of letting a numerically zero scale turn every weight into
        # 0/0.
        floor = _ROBUST_SCALE_FLOOR * reference
        if scale < floor:
            scale = floor
        if scale <= 0.0:
            # Every residual is exactly zero: there is nothing to reweight.
            return predictions
        ratio = residual / (_BISQUARE_TUNING * scale)
        extra = np.where(np.abs(ratio) < 1.0, (1.0 - ratio * ratio) ** 2, 0.0)
        predictions = _all_local_predictions(
            scaled, values, extra, degree=degree, neighbors=neighbors
        )
    return predictions


def _all_local_predictions(
    scaled: np.ndarray,
    values: np.ndarray,
    extra: np.ndarray,
    *,
    degree: int,
    neighbors: int,
) -> np.ndarray:
    return np.array(
        [
            _local_prediction(scaled, values, index, extra, degree=degree, neighbors=neighbors)
            for index in range(scaled.size)
        ]
    )


def _local_prediction(
    scaled: np.ndarray,
    values: np.ndarray,
    index: int,
    extra: np.ndarray,
    *,
    degree: int,
    neighbors: int,
) -> float:
    """The local-polynomial prediction of point ``index`` from the others."""
    distance = np.abs(scaled - scaled[index])
    distance[index] = np.inf
    order = np.argsort(distance, kind="stable")
    bandwidth = distance[order[neighbors - 1]]
    if bandwidth <= 0.0:
        # The nearest neighbours sit on the point's own coordinate; only they
        # carry information, and at full weight.
        weight = np.where(distance == 0.0, extra, 0.0)
    else:
        ratio = distance / bandwidth
        weight = np.where(ratio < 1.0, (1.0 - ratio**3) ** 3, 0.0) * extra
    weight[index] = 0.0
    active = weight > 0.0
    count = int(active.sum())
    if count == 0:  # unreachable with neighbors >= 1; kept as a guard
        return float(values[index])
    offsets = scaled[active] - scaled[index]
    active_values = values[active]
    active_weight = weight[active]
    distinct = np.unique(offsets).size
    if count < 2 or distinct < 2:
        return float(np.sum(active_weight * active_values) / np.sum(active_weight))
    local_degree = min(degree, distinct - 1)
    design = np.vstack([offsets**power for power in range(local_degree + 1)]).T
    root = np.sqrt(active_weight)
    coefficients, *_ = np.linalg.lstsq(design * root[:, None], active_values * root, rcond=None)
    return float(coefficients[0])


def summarize_deviations(deviations) -> DeviationSummary:
    """Summarise one series' deviations: its own scatter and its worst point.

    ``median_deviation`` is the robust scatter of the series' leave-one-out
    deviations -- the scale a single flagged point has to be read against.
    Where that scatter is itself several percent, the series varies
    point-to-point as much as any displacement a review threshold would flag,
    and no single flag separates a mistyped cell from genuine roughness (the
    module docstring's model-limitation case).  The deviations must all
    belong to one ``(file, series)`` pair; an empty or mixed collection is
    refused, not averaged over.
    """
    rows = tuple(deviations)
    if not rows:
        raise ValueError(
            "there are no deviations to summarise. Next step: run "
            "check_series on the series first, and pass its result here."
        )
    keys = {(row.file, row.series) for row in rows}
    if len(keys) != 1:
        raise ValueError(
            f"the deviations cover {len(keys)} (file, series) pairs, and a "
            f"summary is per series. Next step: group the deviations by "
            f"their (file, series) fields and summarise one group at a time."
        )
    worst = max(rows, key=lambda row: row.rel_deviation)
    file, series = next(iter(keys))
    return DeviationSummary(
        file=file,
        series=series,
        n_points=len(rows),
        median_deviation=float(np.median([row.rel_deviation for row in rows])),
        max_deviation=float(worst.rel_deviation),
        max_label=worst.label,
        n_suspect=sum(row.grade == GRADE_SUSPECT for row in rows),
        n_significant=sum(row.grade == GRADE_SIGNIFICANT for row in rows),
    )


def _grade(relative: float, warn: float, severe: float) -> str:
    if relative > severe:
        return GRADE_SIGNIFICANT
    if relative > warn:
        return GRADE_SUSPECT
    return GRADE_OK


# --- the directory sweep -----------------------------------------------------------


def check_literature(
    root,
    *,
    warn: float = DEFAULT_WARN,
    severe: float = DEFAULT_SEVERE,
    span: float = _DEFAULT_SPAN,
    min_neighbors: int = _DEFAULT_MIN_NEIGHBORS,
    degree: int = _DEFAULT_DEGREE,
    min_points: int = _DEFAULT_MIN_POINTS,
    robust: bool = True,
    robust_iterations: int = _DEFAULT_ROBUST_ITERATIONS,
) -> Sweep:
    """Run the trend check over every registered anchor file under ``root``.

    Every file of :data:`FILE_PLANS` runs in registry order; a registered file
    that is missing raises (the sweep is a check of a directory that is
    supposed to be complete).  Series shorter than ``min_points`` are recorded
    as refusals, not checked; files whose plan is ``"none"`` are recorded with
    the plan's reason.  The keyword controls are those of :func:`check_series`.
    """
    directory = Path(root)
    files: list[FileOutcome] = []
    for name, plan in FILE_PLANS.items():
        target = directory / name
        if not target.is_file():
            raise ValueError(
                f"registered anchor file not found: {target}. Next step: "
                f"check the fixtures/literature directory (the README lists "
                f"the thirteen tables and their sources)."
            )
        if plan.kind == "none":
            files.append(FileOutcome(file=name, not_analysed=plan.reason))
            continue
        outcomes: list[SeriesOutcome] = []
        for series in load_series(target, plan=plan):
            if len(series.points) < min_points:
                outcomes.append(
                    SeriesOutcome(
                        series=series,
                        refusal=(
                            f"{len(series.points)} points, below "
                            f"min_points={min_points}: too short for a "
                            f"leave-one-out trend"
                        ),
                    )
                )
                continue
            outcomes.append(
                SeriesOutcome(
                    series=series,
                    deviations=check_series(
                        series,
                        warn=warn,
                        severe=severe,
                        span=span,
                        min_neighbors=min_neighbors,
                        degree=degree,
                        min_points=min_points,
                        robust=robust,
                        robust_iterations=robust_iterations,
                    ),
                )
            )
        files.append(FileOutcome(file=name, outcomes=tuple(outcomes)))
    return Sweep(root=directory, files=tuple(files))


# --- the tidy deviation table ------------------------------------------------------

_DEVIATION_COLUMNS = (
    "file",
    "series",
    "point",
    "x",
    "value",
    "prediction",
    "residual",
    "rel_deviation",
    "grade",
    "uncertainty",
    "flag",
    "note",
)


def write_deviation_csv(path, deviations) -> Path:
    """Write the deviation table as a tidy CSV file (the plot-ready style).

    One header line with underscored column names, one row per deviation at
    ten significant digits (:func:`plot_csv.plot_number`), text fields quoted
    only where quoting is required (:func:`plot_csv.csv_field`).  The rows are
    written in the order given; ``check_literature`` produces them in
    sweep order, and the loader orders each series by coordinate, so the same
    input gives a byte-identical file -- there is no timestamp.  An absent
    uncertainty is an empty cell.  The parent directory must already exist.
    """
    target = Path(path)
    lines = [",".join(_DEVIATION_COLUMNS)]
    for deviation in deviations:
        cells = (
            deviation.file,
            deviation.series,
            deviation.label,
            plot_number(deviation.x),
            plot_number(deviation.value),
            plot_number(deviation.prediction),
            plot_number(deviation.residual),
            plot_number(deviation.rel_deviation),
            deviation.grade,
            "" if deviation.uncertainty is None else plot_number(deviation.uncertainty),
            deviation.flag,
            deviation.note,
        )
        lines.append(",".join(csv_field(cell) for cell in cells))
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target
