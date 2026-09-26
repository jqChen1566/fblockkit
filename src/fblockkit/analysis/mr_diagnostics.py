"""A3: combined multi-reference (MR) character diagnosis -- one labelled panel.

Purpose: one ORCA output carries at most a few MR-character indicators. This module
collects the ones the output has, prints them in one table together with the published
line each is judged against, and gives a combined verdict -- so the user can judge
whether a single-reference treatment is adequate. Analysis only: numbers, thresholds and
verdicts; the module never writes to the user's files.

Interface (the analysis-layer convention, see ``analysis/__init__.py``):
``accepts(result) -> bool``, ``run(result) -> ReportSection``, ``evidence()`` -> the
provenance of every threshold and statement used here.

Indicators (parser section keys; see ``parsers/orca.py``)

1. T1 diagnostic -- ``sections["cc"]["t1"]``, printed by a coupled-cluster (CCSD) job in
   the "COUPLED CLUSTER ENERGY" block, next to the singles norm it is derived from
   (``sections["cc"]["singles_norm"]``). Line: 0.02 for closed-shell references, the
   **commonly used screening convention** (Lee & Taylor 1989) -- a screening convention,
   not a hard law. T1 and the singles norm are one indicator, not two: measured on the
   fixtures, singles_norm / T1 is sqrt(10) = 3.162278 for both N2 outputs and
   sqrt(14) = 3.741657 for F2, i.e. the square root of the correlated-electron count.
   ORCA prints T1 only, not D1: the string "D1" occurs nowhere in the three CCSD
   fixtures (measured 2026-09-26; the jobs and their T1 values are recorded in the
   fixture table, ``fixtures/orca/README.md``), so this panel offers no D1-based line.

2. Natural-orbital occupations -- ``sections["casscf"]["active_occupations"]``, the
   N(occ)= line of a CASSCF output (a tuple of floats, possibly empty). Two indicators
   are read from it, both free of invented thresholds:

2a. fractional active orbitals: the count of active orbitals with a *fractional*
    occupation, defined as strictly inside the 0.02-1.98 window **and** not integer-
    valued within ``INTEGER_TOLERANCE``. The window is ORCA's auto-ICE active-space
    convention (manual Sec. 3.14), borrowed here as the definition of a fractional
    orbital; no published threshold exists for the count itself (the review that would
    be the place for one gives none). Rationale for the integer test (this module's
    convention, provisional): a single configuration has integer natural occupations
    (0, 1 or 2), so an occupation inside the window that is away from an integer means
    the active space is not one configuration -- while an occupation of exactly 1 is a
    single unpaired electron and is not fractional.

2b. ``max(s_bound)``, the largest single-orbital entropy **bound**, computed with the A2
    machinery (``analysis.entropy`` -- one convention for the whole toolkit). The bound
    is the maximum-entropy completion of the four occupation weights for a spin-summed
    occupation; it never under-estimates the true four-state entropy, so the 0.14 line
    (Stein & Reiher 2016) transfers only in the **exclusion direction**: a bound at or
    below the line excludes multireference character on this metric, while a bound above
    it decides nothing (the true entropy needs the 2-RDM route; see the A2 module and
    the A6 section). Wardzala et al. 2026 quotes the same 0.14 line as the
    **M-diagnostic** line: a second, different use of the number, stated as such in the
    body. The Zs(1) > 0.2 line means something else again and is not used here.

Combined verdict (three values): "MR character indicated" when at least one available
indicator points at MR character (T1 above its screening line, or a non-zero fractional
count); "MR character not indicated" when no available indicator points at it -- the
entropy bound then acts as an exclusion where it can (a bound at or below the line
excludes multireference character on that metric; a bound above the line is indecisive
and is reported as such, never as a positive signal); "cannot be judged (missing data)"
when no indicator is available at all. The body always names which indicators fired,
which were indecisive and which were unavailable: one available indicator never stands
in silently for the panel. An output with no indicator at all (``accepts`` False) still
returns that third verdict instead of raising -- the report pipeline gates on
``accepts``, so this path serves direct callers.

Calibration range (caveat carried in the body): both lines were calibrated on main-group
and 3d systems; the review that quotes the 0.14 M-diagnostic line contains **no
lanthanide or actinide example** (its systems are main-group molecules and 3d/4d/5d
transition-metal complexes). For an f-block system the lines are therefore a default
starting value that must be verified, not a verdict.

Measured on the fixtures (ORCA 6.1.1; the values are pinned in
``tests/test_mr_diagnostics.py``):

- ``n2_ccsd.out`` T1 = 0.013031252 and ``f2_ccsd.out`` T1 = 0.011538126 (both below the
  0.02 line), ``n2_stretch_ccsd.out`` T1 = 0.045766132 (above) -- the two sides of the
  screening;
- ``n2_casscf_nevpt2.out`` N(occ) = (1.99236, 1.70922, 1.70922, 0.29360, 0.29360,
  0.00200): 4 fractional orbitals, max(s_bound) = 0.834236 (above 0.14 -- not decidable,
  a loose bound);
- ``generated_ce3_sarc2.out`` (f1 Ce3+) N(occ) = (1, 0, 0, 0, 0, 0, 0): **0** fractional
  orbitals, and max(s_bound) = 1.386294 = ln 4 -- the open-shell single-occupation limit
  of the bound, not multireference character. The panel therefore does not read this as
  a positive signal: the entropy bound is an exclusion test and the verdict rests on the
  occupation count (and, if present, T1);
- ``fblock_dft_la_complex.out`` (a plain DFT single point) carries neither indicator, so
  ``accepts`` is False.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from ..knowledge.models import (
    EVIDENCE_LITERATURE,
    EVIDENCE_MANUAL,
    EVIDENCE_MEASURED,
    Evidence,
    ParseResult,
    ReportSection,
)
from .entropy import MULTIREFERENCE_THRESHOLD, single_orbital_entropy_bound

TITLE = "A3 multi-reference character"

T1_SCREENING_THRESHOLD = 0.02  # closed-shell screening convention (Lee & Taylor 1989)
FRACTIONAL_WINDOW_LOW = 0.02   # active-orbital window of ORCA's auto-ICE convention
FRACTIONAL_WINDOW_HIGH = 1.98  # (borrowed as the definition of a fractional orbital)
INTEGER_TOLERANCE = 1e-3       # an occupation within this of 0/1/2 counts as integer-valued

VERDICT_INDICATED = "MR character indicated"
VERDICT_NOT_INDICATED = "MR character not indicated"
VERDICT_UNKNOWN = "cannot be judged (missing data)"

_LABEL_T1 = "T1 diagnostic"
_LABEL_FRACTIONAL = "fractional active orbitals"
_LABEL_ENTROPY = "max single-orbital entropy bound"

_STATUS_ABOVE = "above the line"
_STATUS_AT_OR_BELOW = "at or below the line"
_STATUS_PRESENT = "present"
_STATUS_NONE = "none"
_STATUS_UNAVAILABLE = "unavailable"
_STATUS_EXCLUDES = "at or below: excludes"
_STATUS_INDECISIVE = "above: indecisive"

_NEEDS_CCSD = (
    "T1 needs a coupled-cluster (CCSD) output: ORCA prints it inside the "
    "\"COUPLED CLUSTER ENERGY\" block."
)
_NEEDS_CASSCF = (
    "the occupation-based indicators need a CASSCF output with an N(occ)= line."
)

_EVIDENCE_T1 = Evidence(
    kind=EVIDENCE_LITERATURE,
    text=(
        "T1 diagnostic: the norm of the coupled-cluster single-excitation amplitudes, "
        "normalised by the number of correlated electrons. Values above about 0.02 for a "
        "closed-shell reference are the commonly used screening convention for whether a "
        "single-reference treatment is adequate; it is a screening convention, not a hard "
        "law, and a value below the line does not prove single-reference adequacy."
    ),
    ref="Lee T. J., Taylor P. R., Int. J. Quantum Chem., 1989, 36(S23), 199-207, DOI 10.1002/qua.560360824",
    url="https://doi.org/10.1002/qua.560360824",
    bibkey="lee1989diagnostic",
)

_EVIDENCE_ENTROPY = Evidence(
    kind=EVIDENCE_LITERATURE,
    text=(
        "Multireference criterion: max s1 > 0.14, with s1 the four-state single-orbital "
        "entropy (Stein & Reiher 2016). This panel applies the line to the A2 bound one-"
        "sidedly: a bound at or below 0.14 excludes multireference character on this "
        "metric, a bound above it decides nothing because the bound is loose."
    ),
    ref="Stein C. J., Reiher M., J. Chem. Theory Comput., 2016, 12(4), 1760-1771, DOI 10.1021/acs.jctc.6b00156",
    url="https://doi.org/10.1021/acs.jctc.6b00156",
    bibkey="stein2016automated",
)

_EVIDENCE_M_DIAGNOSTIC = Evidence(
    kind=EVIDENCE_LITERATURE,
    text=(
        "The same 0.14 line is quoted by this review as the M-diagnostic line: M > 0.14 "
        "marks significant multiconfigurational character (a CCSD(T)/CASSCF pre-screening "
        "for whether a multireference treatment is needed). The number is shared between "
        "two different diagnostics; the Zs(1) > 0.2 line means something else again and is "
        "not used by this panel."
    ),
    ref=(
        "Wardzala J. J. et al., Chem. Rev., 2026, 126(8), 4592-4618, "
        "DOI 10.1021/acs.chemrev.5c00866 (Sec. 2.3.2)"
    ),
    url="https://doi.org/10.1021/acs.chemrev.5c00866",
    bibkey="wardzala2026multireference",
)

_EVIDENCE_FBLOCK_SCOPE = Evidence(
    kind=EVIDENCE_LITERATURE,
    text=(
        "Calibration range of these thresholds: the review's systems are main-group "
        "molecules and 3d/4d/5d transition-metal complexes -- it contains no lanthanide or "
        "actinide example at all. For an f-block system the lines quoted by this panel are "
        "a default starting value that must be verified, not a verdict."
    ),
    ref=(
        "Wardzala J. J. et al., Chem. Rev., 2026, 126(8), 4592-4618, "
        "DOI 10.1021/acs.chemrev.5c00866; scope noted in this group's reading record "
        "文献细读/细读_自动活性空间综述_ChemRev2026.md (2026-09-25, fact 14)"
    ),
    url="https://doi.org/10.1021/acs.chemrev.5c00866",
    bibkey="wardzala2026multireference",
)

_EVIDENCE_WINDOW = Evidence(
    kind=EVIDENCE_MANUAL,
    text=(
        "\"All orbitals between occupation number say 1.98 down to 0.02 will be included "
        "in the active space.\" (the auto-ICE selection convention; this panel borrows the "
        "0.02/1.98 interval as its definition of a fractionally occupied active orbital)"
    ),
    ref="ORCA 6.1 manual Sec. 3.14 (ICE-CI and auto-ICE)",
    url="https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/iceci.html",
)

_EVIDENCE_MEASURED = Evidence(
    kind=EVIDENCE_MEASURED,
    text=(
        "Measured on the ORCA 6.1.1 fixtures (2026-09-26; the jobs are listed in the "
        "fixture table, fixtures/orca/README.md): the CCSD outputs print a T1 diagnostic "
        "and its singles norm, and no D1 diagnostic (the string \"D1\" occurs nowhere in "
        "n2_ccsd.out, f2_ccsd.out or n2_stretch_ccsd.out); T1 = 0.013031252 (N2), "
        "0.011538126 (F2), 0.045766132 (stretched N2); singles_norm / T1 = sqrt(10) = "
        "3.162278 on both N2 outputs and sqrt(14) = 3.741657 on F2 (one indicator, not "
        "two). CASSCF: N2 CAS(6,6) N(occ) = (1.99236, 1.70922, 1.70922, 0.29360, 0.29360, "
        "0.00200) gives 4 fractional orbitals and max(s_bound) = 0.834236 (indecisive); "
        "Ce3+ 4f1 N(occ) = (1, 0, 0, 0, 0, 0, 0) gives 0 fractional orbitals and "
        "max(s_bound) = 1.386294 = ln 4, the open-shell single-occupation limit of the "
        "bound (indecisive, not a positive signal)."
    ),
    ref=(
        "fixtures/orca/ (ORCA 6.1.1 outputs; the D1 measurement is also recorded in "
        "src/fblockkit/parsers/orca.py)"
    ),
)


@dataclass(frozen=True)
class Indicator:
    """One MR-character indicator of the panel, as found in one output."""

    key: str          # "t1" / "fractional_occupations" / "max_entropy"
    label: str        # the name printed in the table
    available: bool   # False = this output does not carry the indicator
    fired: bool       # True = the indicator points at MR character (meaningful when available)
    value_text: str   # the measured value, formatted for the table
    line_text: str    # the published line, formatted for the table
    status: str       # the table's verdict cell
    summary: str      # compact form for the combined verdict's fired / not-fired lists
    reason: str       # the per-indicator verdict sentence
    requirement: str = ""  # what the indicator needs from the output (when unavailable)


def _num(value: float) -> str:
    """Format a measurement for the report: the digits as printed, no rounding gloss."""
    return f"{value:.9g}"


def _is_fractional(occupation: float) -> bool:
    """Strictly inside the 0.02-1.98 window and not integer-valued (``INTEGER_TOLERANCE``)."""
    if not FRACTIONAL_WINDOW_LOW < occupation < FRACTIONAL_WINDOW_HIGH:
        return False
    return abs(occupation - round(occupation)) > INTEGER_TOLERANCE


def fractional_orbitals(occupations: Sequence[float]) -> tuple[float, ...]:
    """The active orbitals with a fractional occupation, in input order.

    Fractional = strictly inside the 0.02-1.98 window (ORCA's auto-ICE active-space
    convention, manual Sec. 3.14) **and** not within ``INTEGER_TOLERANCE`` of an integer.
    The integer test is this module's convention (provisional): a single configuration
    has integer natural occupations (0, 1 or 2), so an occupation inside the window away
    from an integer means the active space is not one configuration; an occupation of
    exactly 1 is a single unpaired electron and is not fractional. The tolerance is a
    numerical guard, larger than the printed precision of the N(occ)= line (5 decimals),
    so residual noise on a near-integer occupation is not reported as a fractional orbital.
    """
    return tuple(n for n in occupations if _is_fractional(n))


def max_entropy_bound(occupations: Sequence[float]) -> float:
    """max(s_bound) over the active orbitals (0.0 when there is no occupation).

    The bound itself comes from the A2 module -- one convention for the whole toolkit;
    the formula and its one-sided use are defined there.
    """
    return max((single_orbital_entropy_bound(n) for n in occupations), default=0.0)


def _t1_indicator(t1: float | None, singles_norm: float | None) -> Indicator:
    if t1 is None:
        return Indicator(
            key="t1",
            label=_LABEL_T1,
            available=False,
            fired=False,
            value_text="n/a",
            line_text=_num(T1_SCREENING_THRESHOLD),
            status=_STATUS_UNAVAILABLE,
            summary="n/a",
            reason=(
                f"not available: {_NEEDS_CCSD} Next step: run a CCSD single point at the "
                "same geometry and rerun this panel on that output."
            ),
            requirement="needs a CCSD output",
        )
    line = _num(T1_SCREENING_THRESHOLD)
    norm = (
        f" (singles norm {_num(singles_norm)}; T1 is derived from it, so the two are one "
        "indicator, not two)"
        if singles_norm is not None
        else ""
    )
    fired = t1 > T1_SCREENING_THRESHOLD
    convention = (
        f"The {line} line is the commonly used screening convention for closed-shell "
        "references (Lee & Taylor 1989), not a hard law. ORCA prints T1 only, not D1 "
        "(measured on the fixtures)."
    )
    if fired:
        reason = (
            f"value {_num(t1)}{norm}: above the {line} screening line for closed-shell "
            "references -> this indicator flags the single-reference treatment. "
            + convention
        )
        summary = f"{_num(t1)} > {line}"
    else:
        reason = (
            f"value {_num(t1)}{norm}: at or below the {line} screening line for "
            "closed-shell references -> this indicator gives no evidence of multi-"
            "reference character; a value below the line is a screening result, not proof "
            "of single-reference adequacy. " + convention
        )
        summary = f"{_num(t1)} <= {line}"
    return Indicator(
        key="t1",
        label=_LABEL_T1,
        available=True,
        fired=fired,
        value_text=_num(t1),
        line_text=line,
        status=_STATUS_ABOVE if fired else _STATUS_AT_OR_BELOW,
        summary=summary,
        reason=reason,
    )


def _fractional_indicator(occupations: tuple[float, ...]) -> Indicator:
    if not occupations:
        return Indicator(
            key="fractional_occupations",
            label=_LABEL_FRACTIONAL,
            available=False,
            fired=False,
            value_text="n/a",
            line_text=f"window {_num(FRACTIONAL_WINDOW_LOW)}-{_num(FRACTIONAL_WINDOW_HIGH)}",
            status=_STATUS_UNAVAILABLE,
            summary="n/a",
            reason=(
                f"not available: {_NEEDS_CASSCF} Next step: run a CASSCF job with the "
                "target active space and rerun this panel on that output."
            ),
            requirement="needs a CASSCF output with an N(occ)= line",
        )
    window = f"{_num(FRACTIONAL_WINDOW_LOW)}-{_num(FRACTIONAL_WINDOW_HIGH)}"
    fractional = fractional_orbitals(occupations)
    count = len(fractional)
    value_text = f"{count} of {len(occupations)} in {window}"
    fired = count > 0
    if fired:
        reason = (
            f"{count} of {len(occupations)} active orbitals have a fractional occupation "
            f"(inside {window} and not integer-valued) -> the active-space wavefunction is "
            "not one configuration. The window is ORCA's active-space convention "
            "(auto-ICE) and is used here as the definition of a fractional orbital; the "
            "count itself has no published threshold (the review gives none), so the "
            "panel reports the count, not a line."
        )
    else:
        reason = (
            f"no active orbital has a fractional occupation (all {len(occupations)} "
            f"occupations are integer-valued: 0, 1 or 2) -> this indicator gives no "
            "evidence of multi-reference character. The window is ORCA's active-space "
            "convention (auto-ICE); the count itself has no published threshold. An "
            "occupation of exactly 1 (a single unpaired electron) is not fractional."
        )
    return Indicator(
        key="fractional_occupations",
        label=_LABEL_FRACTIONAL,
        available=True,
        fired=fired,
        value_text=value_text,
        line_text=f"window {window}",
        status=_STATUS_PRESENT if fired else _STATUS_NONE,
        summary=f"{count} of {len(occupations)} fractional",
        reason=reason,
    )


def _entropy_indicator(occupations: tuple[float, ...]) -> Indicator:
    if not occupations:
        return Indicator(
            key="max_entropy",
            label=_LABEL_ENTROPY,
            available=False,
            fired=False,
            value_text="n/a",
            line_text=_num(MULTIREFERENCE_THRESHOLD),
            status=_STATUS_UNAVAILABLE,
            summary="n/a",
            reason=(
                f"not available: {_NEEDS_CASSCF} Next step: run a CASSCF job with the "
                "target active space and rerun this panel on that output."
            ),
            requirement="needs a CASSCF output with an N(occ)= line",
        )
    line = _num(MULTIREFERENCE_THRESHOLD)
    max_s = max_entropy_bound(occupations)
    excludes = max_s <= MULTIREFERENCE_THRESHOLD
    convention = (
        "Convention: s_bound is the maximum-entropy completion of the four occupation "
        "weights for a spin-summed occupation (A2 module), so it never under-estimates "
        "the true four-state entropy; it is tight only for uncorrelated spin channels. An "
        "orbital with occupation 1 gives s_bound = ln 4 = 1.386294 by construction. The "
        "0.14 line (Stein & Reiher 2016) applies to the true four-state entropy; Wardzala "
        "et al. 2026 quotes the same number as the M-diagnostic line. The Zs(1) > 0.2 "
        "line means something else again and is not used here."
    )
    if excludes:
        reason = (
            f"max(s_bound) = {max_s:.6f} <= {line}: the true four-state entropy is also "
            "below the Stein & Reiher line, so this metric excludes a multireference "
            "reading of the active space. " + convention
        )
    else:
        reason = (
            f"max(s_bound) = {max_s:.6f} > {line}: indecisive -- the bound is above the "
            "line, which is not evidence of multireference character. The true four-state "
            "entropy needs the 2-RDM route (orca_2json RDM2_aa/ab/bb plus orca_loc "
            "orbitals; see the A6 section). " + convention
        )
    return Indicator(
        key="max_entropy",
        label=_LABEL_ENTROPY,
        available=True,
        fired=False,  # the bound can exclude; it can never confirm MR character
        value_text=f"{max_s:.6f}",
        line_text=f"{line} (exclusion)",
        status=_STATUS_EXCLUDES if excludes else _STATUS_INDECISIVE,
        summary=f"{max_s:.6f} {'<=' if excludes else '>'} {line}"
        + (" (excludes)" if excludes else " (indecisive)"),
        reason=reason,
    )


def collect(result: ParseResult) -> tuple[Indicator, ...]:
    """The panel's indicators for one parse result, in table order (always three items)."""
    cc = result.sections.get("cc") or {}
    casscf = result.sections.get("casscf") or {}
    occupations = tuple(casscf.get("active_occupations") or ())
    return (
        _t1_indicator(cc.get("t1"), cc.get("singles_norm")),
        _fractional_indicator(occupations),
        _entropy_indicator(occupations),
    )


def combined_verdict(indicators: Sequence[Indicator]) -> str:
    """The panel verdict from the indicators that are available.

    One of "MR character indicated" / "MR character not indicated" /
    "cannot be judged (missing data)"; an empty panel (no indicator available) gives the
    third value. A verdict is drawn from the available indicators only -- a missing one is
    named by the caller, never silently replaced.
    """
    available = tuple(item for item in indicators if item.available)
    if not available:
        return VERDICT_UNKNOWN
    if any(item.fired for item in available):
        return VERDICT_INDICATED
    return VERDICT_NOT_INDICATED


def accepts(result: ParseResult) -> bool:
    """True when the output carries at least one indicator this panel uses."""
    cc = result.sections.get("cc") or {}
    casscf = result.sections.get("casscf") or {}
    return bool(cc.get("t1") is not None or casscf.get("active_occupations"))


def _row(label: str, value: str, line: str, status: str) -> str:
    return f"  {label:<34}{value:<24}{line:<20}{status}"


def _caveats(indicators: tuple[Indicator, ...], available: tuple[Indicator, ...]) -> list[str]:
    lines = [
        "Caveats (carried with the verdict):",
        "- Threshold calibration range: the T1 0.02 line and the max(s1) 0.14 line were "
        "calibrated on main-group and 3d systems; the review that quotes the 0.14 "
        "M-diagnostic line contains no lanthanide or actinide example at all (Wardzala et "
        "al. 2026). For an f-block system these are a default starting value that must be "
        "verified, not a verdict.",
        "- Two diagnostics share the number 0.14: max(s1) > 0.14 (Stein & Reiher 2016) and "
        "the M diagnostic > 0.14 (Wardzala et al. 2026). The Zs(1) > 0.2 line means "
        "something else again and is not used by this panel.",
    ]
    if any(item.key == "t1" and item.available for item in indicators):
        lines.append(
            "- A T1 at or below the line is a screening result, not proof of "
            "single-reference adequacy; ORCA prints T1 only, not D1 (measured on the "
            "fixtures), so no D1-based line is offered."
        )
    if any(
        item.key in ("fractional_occupations", "max_entropy") and item.available
        for item in indicators
    ):
        lines.append(
            "- Entropy bound convention: s_bound is the maximum-entropy completion of the "
            "four occupation weights for a spin-summed occupation (A2 module); it never "
            "under-estimates the true four-state entropy and is used here as an exclusion "
            "test only. An orbital with occupation 1 gives s_bound = ln 4 = 1.386294 by "
            "construction (measured on the Ce3+ fixture: N(occ) = (1, 0, 0, 0, 0, 0, 0), "
            "0 fractional orbitals, max(s_bound) = 1.386294 -> indecisive). The true "
            "four-state entropy needs the 2-RDM route (orca_2json RDM2 + orca_loc)."
        )
    if len(available) < len(indicators):
        lines.append(
            "- The verdict rests on the available indicators only; the ones missing from "
            "this output are named above and are not replaced by another indicator."
        )
    return lines


def run(result: ParseResult) -> ReportSection:
    """Build the A3 report section: indicator table, per-indicator verdicts, combined
    verdict and caveats."""
    indicators = collect(result)
    available = tuple(item for item in indicators if item.available)
    fired = tuple(item for item in available if item.fired)
    not_fired = tuple(item for item in available if not item.fired)
    missing = tuple(item for item in indicators if not item.available)

    if not available:
        return ReportSection(
            title=TITLE,
            body="\n".join(
                [
                    "Refused: this output carries no multi-reference indicator this panel "
                    "can use.",
                    f"- {_LABEL_T1}: not available (it is printed only by a coupled-cluster job).",
                    f"- {_LABEL_FRACTIONAL} and {_LABEL_ENTROPY}: not available (they need "
                    "a CASSCF output with an N(occ)= line).",
                    "",
                    f"Combined verdict: {VERDICT_UNKNOWN}",
                    "  unavailable: " + "; ".join(item.label for item in indicators),
                    "",
                    "Next step: run a CCSD single point for the T1 indicator, or a CASSCF "
                    "job with the target active space for the occupation indicators, then "
                    "rerun this panel on that output.",
                ]
            ),
        )

    lines = [
        "One panel over the multi-reference (MR) character indicators this output "
        "carries; every line is labelled with its source. The verdict below is drawn from "
        "the available indicators only.",
        "",
        f"Indicator table ({len(indicators)} indicators; line = the published criterion, "
        "n/a = not available in this output):",
        _row("indicator", "value", "line", "verdict"),
        *(
            _row(item.label, item.value_text, item.line_text, item.status)
            for item in indicators
        ),
        "",
        "Per-indicator verdicts:",
        *(f"- {item.label}: {item.reason}" for item in indicators),
        "",
        f"Combined verdict: {combined_verdict(indicators)}",
        "  fired: "
        + ("; ".join(f"{item.label} ({item.summary})" for item in fired) or "(none)"),
        "  not fired: "
        + ("; ".join(f"{item.label} ({item.summary})" for item in not_fired) or "(none)"),
        "  unavailable: "
        + ("; ".join(f"{item.label} ({item.requirement})" for item in missing) or "(none)"),
    ]
    if missing:
        lines.append(
            f"Panel completeness: {len(available)} of {len(indicators)} indicators "
            "available; the verdict rests on the available indicator(s) only and the "
            "missing one(s) are named here, never replaced silently."
        )
    lines += ["", *_caveats(indicators, available)]
    return ReportSection(title=TITLE, body="\n".join(lines))


def evidence() -> tuple[Evidence, ...]:
    """This panel's provenance (thresholds, calibration range, fixture measurements)."""
    return (
        _EVIDENCE_T1,
        _EVIDENCE_ENTROPY,
        _EVIDENCE_M_DIAGNOSTIC,
        _EVIDENCE_FBLOCK_SCOPE,
        _EVIDENCE_WINDOW,
        _EVIDENCE_MEASURED,
    )
