"""A10: the magnetic-doublet criterion ``g_T * theta_3 = 20`` (an empirical metric).

For a mononuclear Dy(III) single-molecule magnet the relaxation barrier is read
off the Kramers-doublet ladder: a doublet that still has an excited state above
it supports excitation (a step of the barrier), while a doublet whose transverse
g-value is large enough opens ground-state quantum tunnelling (QTM) and caps the
barrier.  The source's empirical metric for that reading is

    g_T = (g1 + g2 + g3 sin(theta3)) / 3,        line:  g_T * theta3 = 20

with ``g1 <= g2 <= g3`` the three principal g values of the doublet and
``theta3`` the angle (in **degrees**) between this doublet's ``g3`` axis and the
ground doublet's ``g3`` axis.  Below the line the doublet can still act as a
step of the barrier; above it, QTM is opened.

This module is pure post-processing: the criterion needs nothing but the
per-doublet g tensors, so it runs from a small JSON table (menu 16).  The table
gives, per doublet, the three principal values and either ``theta3`` directly
(as published tables give it) or the ``g3`` axis itself, in which case the angle
to the reference doublet's axis is computed here.  The angle is folded to
[0, 90] degrees because a principal axis's sign is arbitrary.

What the criterion is and is not
--------------------------------
The source is explicit that this is an empirical metric, not a derived formula:
in principle the ``g1``/``g2`` weights should be reduced as ``g3`` rotates into
the plane, but that information is not available from the historic calculations
it was calibrated on.  The line's physical background is the material's internal
dipolar field, which varies with packing and dilution -- so the line is not a
universal constant.  The calibration set is 19 mononuclear Dy(III) systems
(2016-2025, one consistent methodology: 9-in-7 active space, ANO-RCC three-tier
basis, SA-CASSCF-SO in (Open)Molcas), and it contains exactly one known outlier
([Dy(Cp^ttt)2]+, whose QTM doublet falls deep in the safe zone) -- a failure mode
this module cannot detect from the g tensors alone, and says so in its output.

The calibration numbers are regression-checked (all 38 tabulated g_T values are
recomputed from the raw columns of the source's Table S2) in
``tests/test_magnetic_doublets.py`` against
``fixtures/literature/chilton_s2_dy19.json``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ..knowledge.models import EVIDENCE_LITERATURE, EVIDENCE_MEASURED, Evidence, ReportSection

__all__ = [
    "MagneticError",
    "Doublet",
    "DoubletTable",
    "DoubletVerdict",
    "CriterionReport",
    "g_transverse",
    "theta3_degrees",
    "parse_doublets",
    "analyze",
    "render",
    "run",
    "evidence",
]

#: The criterion line: below it the doublet supports excitation above it, at or
#: above it the doublet facilitates quantum tunnelling (the source's own line).
CRITERION_LINE = 20.0

#: The calibration envelope, recomputed from the source's Table S2 columns by the
#: regression test: the largest "safe" product and the smallest "definite QTM"
#: product.  The line sits in this gap.
SAFE_SIDE_MAX = 15.31
QTM_SIDE_MIN = 30.25

#: The source's recorded outlier (the only doublet whose side of the line
#: disagrees with its assigned role).
KNOWN_OUTLIER = "[Dy(Cp^ttt)2]+"

VERDICT_SUPPORTS = "supports excitation"
VERDICT_QTM = "facilitates QTM"


class MagneticError(ValueError):
    """The doublet table cannot be read or the criterion cannot be applied."""


@dataclass(frozen=True)
class Doublet:
    """One Kramers doublet: principal g values and its angle to the reference axis.

    ``theta3_deg`` and ``axis3`` are the two ways to say where the doublet's
    ``g3`` axis points: an angle already measured (published tables) or the axis
    itself (from a SINGLE_ANISO-style export), the angle then being computed
    against the reference doublet's axis.
    """

    label: str
    g1: float
    g2: float
    g3: float
    theta3_deg: float | None = None
    axis3: tuple[float, float, float] | None = None
    energy: float | None = None


@dataclass(frozen=True)
class DoubletTable:
    """The parsed input table."""

    system: str
    doublets: tuple[Doublet, ...]
    reference: int  # index of the reference doublet (its g3 axis is the quantisation axis)


@dataclass(frozen=True)
class DoubletVerdict:
    """The criterion's reading of one doublet."""

    index: int
    label: str
    energy: float | None
    g1: float
    g2: float
    g3: float
    theta3_deg: float
    g_transverse: float
    product: float
    verdict: str


@dataclass(frozen=True)
class CriterionReport:
    """The verdicts plus the criteria and boundaries the output carries."""

    verdicts: tuple[DoubletVerdict, ...]
    reference: int
    computed_from_axes: bool  # at least one angle was measured against the reference axis
    checks: tuple[str, ...]
    notes: tuple[str, ...]


# --- the criterion arithmetic -------------------------------------------------


def g_transverse(g1: float, g2: float, g3: float, theta3_deg: float) -> float:
    """``g_T = (g1 + g2 + g3 sin(theta3)) / 3`` -- theta3 in DEGREES.

    Degrees is a measured fact of the calibration table: recomputing its g_T
    values with radians disagrees in the second digit (the regression test pins
    this down on all 38 tabulated doublets).
    """
    return (g1 + g2 + g3 * math.sin(math.radians(theta3_deg))) / 3.0


def theta3_degrees(
    axis_a: tuple[float, float, float], axis_b: tuple[float, float, float]
) -> float:
    """The angle (degrees, folded to [0, 90]) between two principal axes.

    The fold is required because an eigenvector's sign is arbitrary: the source's
    own theta3 values run to 89.9 degrees, i.e. it reports the unsigned angle.
    """
    norm_a = math.sqrt(sum(value * value for value in axis_a))
    norm_b = math.sqrt(sum(value * value for value in axis_b))
    if norm_a <= 0.0 or norm_b <= 0.0:
        raise MagneticError(
            "a g3 axis has zero length, so no angle can be measured. Next step: give "
            "the unit axis of the largest principal g value."
        )
    cosine = sum(a * b for a, b in zip(axis_a, axis_b)) / (norm_a * norm_b)
    cosine = min(1.0, max(-1.0, cosine))
    return math.degrees(math.acos(abs(cosine)))


# --- the input table ----------------------------------------------------------


def _number(value, what: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise MagneticError(
            f"{what} is {value!r}, which is not a number. Next step: give a number."
        ) from exc
    if not math.isfinite(result):
        raise MagneticError(
            f"{what} is not finite ({value!r}). Next step: give a finite number."
        )
    return result


def _vector(value, what: str) -> tuple[float, float, float]:
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        raise MagneticError(
            f"{what} is {value!r}; three numbers are needed. Next step: give the axis "
            "as [x, y, z]."
        )
    return tuple(_number(entry, what) for entry in value)  # type: ignore[return-value]


def _coerce_g(row: dict, where: str) -> tuple[float, float, float]:
    if "g" not in row:
        raise MagneticError(
            f"doublet {where} has no 'g' entry. Next step: give the three principal "
            "g values as \"g\": [g1, g2, g3]."
        )
    values = row["g"]
    if not isinstance(values, (list, tuple)) or len(values) != 3:
        raise MagneticError(
            f"doublet {where} has 'g' = {values!r}; three values are needed. Next step: "
            "give the three principal g values as a list of three numbers."
        )
    numbers = tuple(_number(entry, f"doublet {where} g value") for entry in values)
    # the source's formula puts the LARGEST value in the sin term; sorting here
    # makes the input order irrelevant (the g3 axis, when given, is defined as the
    # axis of the largest value regardless)
    return tuple(sorted(numbers))  # type: ignore[return-value]


def parse_doublets(payload) -> DoubletTable:
    """Read the doublet table from the JSON payload (menu 16's input schema).

    Schema::

        {
          "system": "optional label",
          "reference": 0,            # optional index or label; default 0
          "doublets": [
            {"label": "ground", "g": [0.41, 0.44, 9.04], "energy": 0.0,
             "axis3": [0.0, 0.0, 1.0]},
            {"label": "1st excited", "g": [0.41, 0.44, 9.04], "theta3": 9.52}
          ]
        }

    Per doublet, exactly one of ``theta3`` (degrees, as published tables give it)
    and ``axis3`` (the unit axis of the largest principal value) says where the
    doublet's g3 axis points; the reference doublet may omit both (its theta3 is
    0 by definition when it provides the axis others are measured against).
    """
    if not isinstance(payload, dict):
        raise MagneticError(
            "the input is not a JSON object. Next step: give an object with a "
            "'doublets' list (see the user guide for the schema)."
        )
    rows = payload.get("doublets")
    if not isinstance(rows, list) or not rows:
        raise MagneticError(
            "the input has no non-empty 'doublets' list. Next step: add the "
            "per-doublet entries, each with 'g' and 'theta3' or 'axis3'."
        )
    doublets: list[Doublet] = []
    for position, row in enumerate(rows):
        if not isinstance(row, dict):
            raise MagneticError(
                f"doublet {position} is {row!r}, not an object. Next step: give each "
                "doublet as an object."
            )
        label = row.get("label")
        where = f"{position} ({label!r})" if label else str(position)
        g1, g2, g3 = _coerce_g(row, where)
        theta3 = row.get("theta3")
        if theta3 is not None:
            theta3 = _number(theta3, f"doublet {where} theta3")
            if not 0.0 <= theta3 <= 180.0:
                raise MagneticError(
                    f"doublet {where} has theta3 = {theta3}, outside [0, 180] degrees. "
                    "Next step: check the units -- the criterion's angle is in degrees."
                )
        axis3 = row.get("axis3")
        if axis3 is not None:
            axis3 = _vector(axis3, f"doublet {where} axis3")
        energy = row.get("energy")
        if energy is not None:
            energy = _number(energy, f"doublet {where} energy")
        doublets.append(
            Doublet(
                label=str(row.get("label") or f"doublet {position}"),
                g1=g1,
                g2=g2,
                g3=g3,
                theta3_deg=theta3,
                axis3=axis3,
                energy=energy,
            )
        )
    reference = payload.get("reference", 0)
    if isinstance(reference, str):
        match = [index for index, item in enumerate(doublets) if item.label == reference]
        if len(match) != 1:
            raise MagneticError(
                f"the reference label {reference!r} matches {len(match)} doublet(s). "
                "Next step: give the reference as an index, or make the label unique."
            )
        reference_index = match[0]
    else:
        reference_index = int(_number(reference, "reference index"))
        if not 0 <= reference_index < len(doublets):
            raise MagneticError(
                f"the reference index {reference_index} is outside the table's "
                f"{len(doublets)} doublet(s). Next step: check the reference."
            )
    return DoubletTable(
        system=str(payload.get("system") or "the given table"),
        doublets=tuple(doublets),
        reference=reference_index,
    )


# --- the analysis -------------------------------------------------------------


def analyze(table: DoubletTable) -> CriterionReport:
    """Apply the criterion to every doublet (the reference's theta3 is 0)."""
    needs_axis = [
        index
        for index, doublet in enumerate(table.doublets)
        if index != table.reference and doublet.theta3_deg is None
    ]
    reference = table.doublets[table.reference]
    if needs_axis:
        if reference.axis3 is None:
            raise MagneticError(
                f"doublet(s) {needs_axis} give no theta3, so their angle must be "
                "measured against the reference doublet's g3 axis -- but the reference "
                "doublet has no 'axis3'. Next step: give the reference's axis3, or give "
                "theta3 directly on every doublet."
            )
        if reference.theta3_deg not in (None, 0.0):
            raise MagneticError(
                f"the reference doublet carries theta3 = {reference.theta3_deg}, but "
                "its g3 axis is the quantisation axis the other doublets are measured "
                "against, so its theta3 is 0 by definition. Next step: drop the "
                "reference's theta3 (or point 'reference' at another doublet)."
            )
    verdicts: list[DoubletVerdict] = []
    computed_from_axes = False
    for index, doublet in enumerate(table.doublets):
        if doublet.theta3_deg is not None:
            theta3 = doublet.theta3_deg
        elif index == table.reference:
            theta3 = 0.0
        elif doublet.axis3 is not None:
            theta3 = theta3_degrees(doublet.axis3, reference.axis3)  # type: ignore[arg-type]
            computed_from_axes = True
        else:
            raise MagneticError(
                f"doublet {index} has neither 'theta3' nor 'axis3'. Next step: give one "
                "of the two (theta3 in degrees, or the unit axis of its largest g value)."
            )
        value = g_transverse(doublet.g1, doublet.g2, doublet.g3, theta3)
        product = value * theta3
        verdicts.append(
            DoubletVerdict(
                index=index,
                label=doublet.label,
                energy=doublet.energy,
                g1=doublet.g1,
                g2=doublet.g2,
                g3=doublet.g3,
                theta3_deg=theta3,
                g_transverse=value,
                product=product,
                verdict=VERDICT_SUPPORTS if product < CRITERION_LINE else VERDICT_QTM,
            )
        )
    reference_line = (
        f"the reference doublet is the table's entry {table.reference} "
        f"('{reference.label}'): its g3 axis is the quantisation axis the computed "
        "angles are measured against, and its own theta3 is 0 by definition"
        if computed_from_axes
        else f"every theta3 is used as given (as a published table measures it, "
        f"against the ground doublet it refers to); entry {table.reference} "
        f"('{reference.label}') is marked as the table's reference row -- set "
        "'reference' to mark another entry"
    )
    checks = (
        f"the criterion: g_T = (g1 + g2 + g3 sin(theta3)) / 3 with theta3 in DEGREES, "
        f"compared with g_T * theta3 = {CRITERION_LINE:g}: below the line the doublet "
        "can still act as a step of the barrier (supports excitation above it), at or "
        "above it the doublet facilitates quantum tunnelling",
        "the three principal values are sorted internally, the largest taking the sin "
        "term (the source's convention); a 'axis3' entry is the axis of that largest "
        "value",
        reference_line,
        "an angle computed from axes is folded to [0, 90] degrees (a principal axis's "
        "sign is arbitrary); theta3 given in the table is used as given",
        f"calibration context: on the source's 19-system table the line sits in the "
        f"gap between the largest 'safe' product ({SAFE_SIDE_MAX:.2f}) and the smallest "
        f"'definite QTM' product ({QTM_SIDE_MIN:.2f})",
    )
    notes = (
        "The metric is empirical, not derived: the source states that the g1/g2 weights "
        "should in principle be reduced as g3 rotates into the plane, but those historic "
        "calculations do not carry that information. The line's physical background is "
        "the material's internal dipolar field, which varies with packing and dilution, "
        "so the line is not a universal constant.",
        "Applicability domain: the line was calibrated on 19 mononuclear Dy(III) "
        "complexes (2016-2025, one consistent methodology -- 9-in-7 active space, "
        "ANO-RCC three-tier basis, SA-CASSCF-SO in (Open)Molcas). Other ions or "
        "nuclearities would need their own calibration.",
        f"Known failure mode the criterion cannot see: in the source's set exactly one "
        f"doublet is misread -- {KNOWN_OUTLIER}'s QTM doublet falls deep in the safe "
        "zone (product 2.52), and the source offers no explanation. A 'safe' verdict "
        "here is therefore the criterion's reading, not a proof of barrier behaviour.",
    )
    return CriterionReport(
        verdicts=tuple(verdicts),
        reference=table.reference,
        computed_from_axes=computed_from_axes,
        checks=checks,
        notes=notes,
    )


# --- output -------------------------------------------------------------------


def render(report: CriterionReport) -> str:
    """The doublet table with each doublet's verdict."""
    lines = [
        "Magnetic-doublet criterion (g_T = (g1 + g2 + g3 sin(theta3)) / 3, degrees; "
        "line g_T * theta3 = 20):",
        f"  {'#':>1}{'':>2} {'doublet':<34} {'E(cm-1)':>8}  {'g1':>6}  {'g2':>6}  "
        f"{'g3':>6}  {'theta3(deg)':>11}  {'g_T':>7}  {'g_T*theta3':>10}  reading",
    ]
    for verdict in report.verdicts:
        energy = f"{verdict.energy:>8.1f}" if verdict.energy is not None else "       -"
        marker = " *" if verdict.index == report.reference else "  "
        lines.append(
            f"  {verdict.index:>1}{marker} {verdict.label:<34} {energy}  "
            f"{verdict.g1:>6.2f}  {verdict.g2:>6.2f}  {verdict.g3:>6.2f}  "
            f"{verdict.theta3_deg:>11.2f}  {verdict.g_transverse:>7.3f}  "
            f"{verdict.product:>10.2f}  {verdict.verdict}"
        )
    lines.append(
        "  (* = the reference doublet, quantisation axis for the computed angles)"
        if report.computed_from_axes
        else "  (* = the reference row; every theta3 was given in the table, no axis was used)"
    )
    lines += ["", "Criteria and boundaries:"]
    for item in report.checks + report.notes:
        lines.append(f"  - {item}")
    return "\n".join(lines)


def run(table: DoubletTable) -> ReportSection:
    """The analyser entry point for menu 16."""
    return ReportSection(
        title="A10 magnetic-doublet criterion (g_T * theta_3)",
        body=render(analyze(table)),
    )


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the criterion and of its regression check."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "The empirical metric g_T = (g1 + g2 + g3 sin(theta3)) / 3 with the "
                "line g_T * theta_3 = 20 (theta3 in degrees), calibrated on 19 "
                "high-performance monometallic Dy(III) single-molecule magnets "
                "(2016-2025, one methodology: 9-in-7 active space, ANO-RCC "
                "three-tier basis, SA-CASSCF-SO in (Open)Molcas; the source's seven "
                "items). The source states the metric is empirical (the g1/g2 "
                "plane-out correction is not available) and records exactly one "
                "outlier, [Dy(Cp^ttt)2]+, whose QTM doublet falls in the safe zone."
            ),
            ref=(
                "Supporting Information, section 4.2-4.3 (the 19-system compilation "
                "and Table S2; the criterion text above it)"
            ),
            url="https://doi.org/10.1039/d5cs00493d",
            bibkey="chilton2025abinitio",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "Regression check of the criterion on the source's own table: all 38 "
                "tabulated g_T values (19 systems x a safe and a definite-QTM "
                "doublet) are recomputed from the raw g/theta3 columns and agree with "
                "the table within its three-decimal rounding (the largest "
                "disagreement is 8e-4); the separation the source claims is "
                "reproduced -- the largest safe product is 15.31, the smallest "
                "definite-QTM product is 30.25, and the single doublet on the wrong "
                "side of the line is the source's recorded outlier."
            ),
            ref="tests/test_magnetic_doublets.py; fixtures/literature/chilton_s2_dy19.json",
        ),
    )
