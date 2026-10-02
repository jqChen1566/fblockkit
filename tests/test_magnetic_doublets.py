"""Regression checks of the g_T criterion (analysis/magnetic_doublets.py).

The fixture is the source's own calibration table (19 mononuclear Dy(III) SMMs,
a "safe" and a "definite QTM" doublet each, with the raw g values and theta3
alongside the tabulated g_T).  The tests recompute every tabulated value, and
pin the two facts the criterion stands on: the line sits in the measured gap,
and exactly one doublet disagrees with it -- the source's recorded outlier.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from fblockkit.analysis import magnetic_doublets
from fblockkit.analysis.magnetic_doublets import (
    CRITERION_LINE,
    KNOWN_OUTLIER,
    QTM_SIDE_MIN,
    SAFE_SIDE_MAX,
    MagneticError,
    g_transverse,
    parse_doublets,
    theta3_degrees,
)

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "literature"
TABLE = FIXTURES / "chilton_s2_dy19.json"


def _fixture() -> dict:
    return json.loads(TABLE.read_text(encoding="utf-8"))


def _doublets():
    """Yield (index, role, row-dict) for every tabulated doublet."""
    for row in _fixture()["rows"]:
        for role in ("safe", "qtm"):
            yield row["index"], role, row[role]


# --- the calibration table ----------------------------------------------------


def test_the_tabulated_g_T_values_are_reproduced():
    """Every one of the 38 tabulated g_T values is recomputed from its raw columns.

    Tolerance 0.001: the table rounds g_T to three decimals and the raw columns
    carry two decimals, so the largest disagreement measures 8e-4 (recorded in
    the fixture's cross_check block).
    """
    worst = 0.0
    for index, role, entry in _doublets():
        recomputed = g_transverse(entry["g1"], entry["g2"], entry["g3"], entry["theta3_deg"])
        worst = max(worst, abs(recomputed - entry["g_T"]))
    assert worst < 0.001, f"worst disagreement {worst}"


def test_the_line_sits_in_the_measured_gap():
    """The safe-side maximum and the QTM-side minimum straddle the line 20."""
    products = {
        (index, role): g_transverse(
            entry["g1"], entry["g2"], entry["g3"], entry["theta3_deg"]
        )
        * entry["theta3_deg"]
        for index, role, entry in _doublets()
    }
    outlier = _fixture()["criterion"]["known_outlier_index"]
    safe_max = max(
        value for (index, role), value in products.items() if role == "safe"
    )
    qtm_min = min(
        value
        for (index, role), value in products.items()
        if role == "qtm" and index != outlier
    )
    assert safe_max == pytest.approx(15.306, abs=1e-2)
    assert qtm_min == pytest.approx(30.245, abs=1e-2)
    # the two constants the module prints are these measured numbers
    assert safe_max == pytest.approx(SAFE_SIDE_MAX, abs=1e-2)
    assert qtm_min == pytest.approx(QTM_SIDE_MIN, abs=1e-2)
    assert SAFE_SIDE_MAX < CRITERION_LINE < QTM_SIDE_MIN


def test_exactly_one_doublet_disagrees_with_the_line():
    """The only doublet on the wrong side of the line is the recorded outlier."""
    outlier = _fixture()["criterion"]["known_outlier_index"]
    disagreements = []
    for index, role, entry in _doublets():
        product = (
            g_transverse(entry["g1"], entry["g2"], entry["g3"], entry["theta3_deg"])
            * entry["theta3_deg"]
        )
        expected_side = "safe" if product < CRITERION_LINE else "qtm"
        if expected_side != role:
            disagreements.append((index, role))
    assert disagreements == [(outlier, "qtm")]
    assert KNOWN_OUTLIER.startswith("[Dy(Cp^ttt)2]+")


def test_the_two_worked_examples_of_the_close_reading():
    """The close reading recomputed two rows by hand; both are reproduced here."""
    rows = {row["index"]: row for row in _fixture()["rows"]}
    first = rows[1]["safe"]
    assert g_transverse(first["g1"], first["g2"], first["g3"], first["theta3_deg"]) == pytest.approx(0.206, abs=1e-3)
    second = rows[2]["qtm"]
    assert g_transverse(second["g1"], second["g2"], second["g3"], second["theta3_deg"]) == pytest.approx(0.741, abs=1e-3)


def test_degrees_not_radians():
    """The formula's angle is in degrees -- a measured fact of the table.

    With radians the reproduction collapses on most rows: that is the negative
    control for the unit (the close reading confirmed the units the same way).
    """
    good = 0
    bad_radians = 0
    for _, _, entry in _doublets():
        degrees = g_transverse(entry["g1"], entry["g2"], entry["g3"], entry["theta3_deg"])
        radians = (
            entry["g1"]
            + entry["g2"]
            + entry["g3"] * math.sin(entry["theta3_deg"])
        ) / 3.0
        if abs(degrees - entry["g_T"]) < 1e-3:
            good += 1
        if abs(radians - entry["g_T"]) > 1e-3:
            bad_radians += 1
    assert good == 38
    assert bad_radians >= 30


# --- the axes path ------------------------------------------------------------


def test_theta3_from_axes_and_the_sign_fold():
    z = (0.0, 0.0, 1.0)
    tilted = (math.sin(math.radians(30.0)), 0.0, math.cos(math.radians(30.0)))
    assert theta3_degrees(tilted, z) == pytest.approx(30.0, abs=1e-9)
    # a principal axis's sign is arbitrary: the fold reads the unsigned angle
    assert theta3_degrees(z, (0.0, 0.0, -1.0)) == pytest.approx(0.0, abs=1e-9)
    mirrored = (-math.sin(math.radians(30.0)), 0.0, math.cos(math.radians(30.0)))
    assert theta3_degrees(mirrored, z) == pytest.approx(30.0, abs=1e-9)
    # the fold also caps at 90 degrees (the source's values run to 89.9)
    assert theta3_degrees((1.0, 0.0, 0.0), (0.0, 1.0, 0.0)) == pytest.approx(90.0, abs=1e-9)


def test_a_zero_axis_is_refused():
    with pytest.raises(MagneticError, match="zero length"):
        theta3_degrees((0.0, 0.0, 0.0), (0.0, 0.0, 1.0))


def test_analyze_computes_the_angle_against_the_reference_axis():
    payload = {
        "system": "synthetic ladder",
        "doublets": [
            {"label": "ground", "g": [0.0, 0.0, 20.0], "axis3": [0.0, 0.0, 1.0]},
            {
                "label": "excited",
                "g": [1.0, 2.0, 10.0],
                "axis3": [math.sin(math.radians(40.0)), 0.0, math.cos(math.radians(40.0))],
            },
            {"label": "table row", "g": [0.41, 0.44, 9.04], "theta3": 9.52},
        ],
    }
    report = magnetic_doublets.analyze(parse_doublets(payload))
    verdicts = {item.label: item for item in report.verdicts}
    assert verdicts["ground"].theta3_deg == pytest.approx(0.0)
    assert verdicts["excited"].theta3_deg == pytest.approx(40.0, abs=1e-9)
    assert verdicts["table row"].theta3_deg == pytest.approx(9.52)
    assert verdicts["table row"].g_transverse == pytest.approx(0.782, abs=1e-3)
    assert verdicts["table row"].verdict == magnetic_doublets.VERDICT_SUPPORTS


def test_analyze_reads_the_source_row_eleven():
    """The manual example: the two tabulated doublets of [Dy(BC4Ph5)2]-."""
    payload = {
        "system": "[Dy(BC4Ph5)2]- (the source's Table S2 row 7)",
        "doublets": [
            {"label": "safe doublet", "g": [0.41, 0.44, 9.04], "theta3": 9.52},
            {"label": "definite QTM doublet", "g": [3.45, 3.73, 5.96], "theta3": 17.58},
        ],
    }
    report = magnetic_doublets.analyze(parse_doublets(payload))
    safe, qtm = report.verdicts
    assert safe.verdict == magnetic_doublets.VERDICT_SUPPORTS
    assert qtm.verdict == magnetic_doublets.VERDICT_QTM
    assert qtm.product == pytest.approx(52.62, abs=5e-2)
    # the reference defaults to the first doublet
    assert report.reference == 0


def test_the_reference_axis_is_required_when_a_theta3_is_missing():
    payload = {
        "doublets": [
            {"label": "ground", "g": [0.0, 0.0, 20.0]},
            {"label": "excited", "g": [1.0, 2.0, 10.0], "axis3": [0.0, 1.0, 0.0]},
        ]
    }
    with pytest.raises(MagneticError, match="no 'axis3'"):
        magnetic_doublets.analyze(parse_doublets(payload))


def test_a_reference_with_a_nonzero_angle_is_refused_when_axes_are_used():
    payload = {
        "doublets": [
            {"label": "ground", "g": [0.0, 0.0, 20.0], "theta3": 5.0, "axis3": [0.0, 0.0, 1.0]},
            {"label": "excited", "g": [1.0, 2.0, 10.0], "axis3": [0.0, 1.0, 0.0]},
        ]
    }
    with pytest.raises(MagneticError, match="0 by definition"):
        magnetic_doublets.analyze(parse_doublets(payload))


def test_a_doublet_with_neither_angle_nor_axis_is_refused():
    payload = {
        "doublets": [
            {"label": "ground", "g": [0.0, 0.0, 20.0], "axis3": [0.0, 0.0, 1.0]},
            {"label": "excited", "g": [1.0, 2.0, 10.0]},
        ]
    }
    with pytest.raises(MagneticError, match="neither 'theta3' nor 'axis3'"):
        magnetic_doublets.analyze(parse_doublets(payload))


# --- input refusals -----------------------------------------------------------


def test_the_three_g_values_are_required():
    with pytest.raises(MagneticError, match="three values"):
        parse_doublets({"doublets": [{"g": [1.0, 2.0], "theta3": 5.0}]})
    with pytest.raises(MagneticError, match="no 'g' entry"):
        parse_doublets({"doublets": [{"theta3": 5.0}]})


def test_an_empty_table_is_refused():
    with pytest.raises(MagneticError, match="'doublets' list"):
        parse_doublets({"doublets": []})
    with pytest.raises(MagneticError, match="'doublets' list"):
        parse_doublets({})


def test_an_out_of_range_angle_is_refused():
    with pytest.raises(MagneticError, match="degrees"):
        parse_doublets({"doublets": [{"g": [1.0, 2.0, 3.0], "theta3": 190.0}]})


def test_the_reference_may_be_given_by_label():
    payload = {
        "reference": "ground state",
        "doublets": [
            {"label": "excited", "g": [1.0, 2.0, 3.0], "theta3": 5.0},
            {"label": "ground state", "g": [0.0, 0.0, 20.0], "axis3": [0.0, 0.0, 1.0]},
        ],
    }
    table = parse_doublets(payload)
    assert table.reference == 1
    report = magnetic_doublets.analyze(table)
    assert report.verdicts[1].theta3_deg == 0.0


# --- the report ---------------------------------------------------------------


def test_the_report_states_the_criteria_and_the_boundaries():
    payload = {
        "doublets": [
            {"label": "safe doublet", "g": [0.41, 0.44, 9.04], "theta3": 9.52},
            {"label": "definite QTM doublet", "g": [3.45, 3.73, 5.96], "theta3": 17.58},
        ]
    }
    body = magnetic_doublets.render(magnetic_doublets.analyze(parse_doublets(payload)))
    assert "DEGREES" in body
    assert "empirical, not derived" in body
    assert "19 mononuclear Dy(III)" in body
    assert KNOWN_OUTLIER in body
    assert "15.31" in body and "30.25" in body
    assert magnetic_doublets.VERDICT_QTM in body and magnetic_doublets.VERDICT_SUPPORTS in body


def test_the_report_says_whether_an_axis_was_used():
    """With every theta3 given, the report must not claim the reference's is 0."""
    table_only = magnetic_doublets.analyze(
        parse_doublets(
            {
                "doublets": [
                    {"label": "a", "g": [0.41, 0.44, 9.04], "theta3": 9.52},
                    {"label": "b", "g": [3.45, 3.73, 5.96], "theta3": 17.58},
                ]
            }
        )
    )
    body = magnetic_doublets.render(table_only)
    assert "no angle needed" in body
    assert "0 by definition" not in body
    with_axes = magnetic_doublets.analyze(
        parse_doublets(
            {
                "doublets": [
                    {"label": "ground", "g": [0.0, 0.0, 20.0], "axis3": [0.0, 0.0, 1.0]},
                    {"label": "excited", "g": [1.0, 2.0, 10.0], "axis3": [0.0, 1.0, 0.0]},
                ]
            }
        )
    )
    body = magnetic_doublets.render(with_axes)
    assert "quantisation axis for the computed angles" in body
    assert "0 by definition" in body


def test_evidence_carries_the_source_and_the_reproduction():
    text = " ".join(
        entry.ref + " " + entry.text + " " + entry.url + " " + entry.bibkey
        for entry in magnetic_doublets.evidence()
    )
    assert "10.1039/d5cs00493d" in text
    assert "chilton2025abinitio" in text
    assert "38" in text and "15.31" in text and "30.25" in text


def test_the_dy_acac_fixture_feeds_menu_16_end_to_end():
    """The ligand-field-bearing multi-doublet case: Dy(acac)3(H2O)2 -- the
    classic mononuclear Dy(III) SMM (fixtures/single_aniso/dy_acac.*, the
    molecular unit of its crystal structure).  The engine ladder is
    cross-checked here against the all-electron CASSCF-SO column of the
    JCTC 2025 benchmark (the same complex as its "1Dy"): the ordering and
    splittings reproduce with a uniform ~7% underestimation from the
    def2-ECP basis against the benchmark's ANO-RCC-DKH one (measured
    mean |delta| = 19.7 cm-1, max 33.5 cm-1; pinned with margin)."""
    from fblockkit.parsers import parse_auto
    from fblockkit.parsers import single_aniso as orca_single_aniso

    fixture = Path(__file__).resolve().parents[1] / "fixtures" / "single_aniso"
    parsed = parse_auto(fixture / "dy_acac.out")
    segments = parsed.sections["single_aniso"]["segments"]
    payload = orca_single_aniso.doublets_payload(segments, system="dy_acac.out")
    rows = payload["doublets"]
    assert len(rows) == 8
    ground = rows[0]
    assert ground["g"][2] == pytest.approx(19.5276, abs=5e-4)  # measured pin
    assert ground["g"][0] < 0.02 and ground["g"][1] < 0.02  # strong Ising
    assert rows[2]["g"][2] == pytest.approx(11.0709, abs=5e-4)
    computed = [row["energy"] for row in rows]
    assert computed == sorted(computed)
    assert computed[1] == pytest.approx(139.335, abs=1e-3)
    assert computed[7] == pytest.approx(496.655, abs=1e-3)
    jctc = (0.0, 153.2, 228.3, 281.2, 313.3, 406.1, 465.0, 530.2)
    deltas = [abs(value - reference) for value, reference in zip(computed, jctc)]
    assert sum(deltas) / len(deltas) < 25.0  # measured 19.7
    assert max(deltas) < 40.0  # measured 33.5
    table = parse_doublets(payload)
    report = magnetic_doublets.analyze(table)
    readings = [verdict.verdict for verdict in report.verdicts]
    assert readings[0] == magnetic_doublets.VERDICT_SUPPORTS
    assert readings[1] == magnetic_doublets.VERDICT_SUPPORTS
    assert all(
        reading == magnetic_doublets.VERDICT_QTM for reading in readings[2:]
    )
