"""Regression checks of the Bleaney comparator (analysis/bleaney.py).

The anchors: the sign structure against the experimental C3-tag Delta-chi
table (Tb +30.87, Tm -14.35, Yb -6.05, in 1e-32 m^3; the 2022 tag study),
the Tb value's implied B_0^2 (~156 cm^-1), the T^-2 law, the direct Dy
value, the two-table scale bridge (1.81 +/- 0.03 per ion), and the closed
form g_J^2 theta_2 J(J+1)(4J(J+1)-3) against the published modern
constants at their printed precision.
"""

from __future__ import annotations

import math

import pytest

from fblockkit.analysis import bleaney
from fblockkit.analysis.bleaney import BleaneyError
from fblockkit.knowledge.models import EVIDENCE_LITERATURE, EVIDENCE_MEASURED


def test_the_direct_value_and_the_t_squared_law():
    value = bleaney.chi_ax("Dy", 1000.0, 300.0)
    assert value / 1e-32 == pytest.approx(226.5120047, abs=1e-6)
    half = bleaney.chi_ax("Dy", 1000.0, 150.0)
    assert half / value == pytest.approx(4.0, abs=1e-12)


def test_the_sign_structure_matches_the_published_tag_table():
    """The formula's signs (set by C_J's signs) reproduce every published tag
    value's sign when inverted through a positive B_0^2."""
    for ion, published in (("Tb", 30.87), ("Tm", -14.35), ("Yb", -6.05)):
        implied = bleaney.implied_b02_cm1(ion, published * 1e-32, 300.0)
        assert implied > 0.0  # a physical axial field; sign comes from C_J
        rebuilt = bleaney.chi_ax(ion, implied, 300.0)
        assert math.copysign(1.0, rebuilt) == math.copysign(1.0, published)


def test_the_tb_magnitude_anchor():
    implied = bleaney.implied_b02_cm1("Tb", 30.87e-32, 300.0)
    assert implied == pytest.approx(156.123, abs=0.01)
    assert 50.0 < implied < 400.0  # the typical lanthanide-tag range


def test_the_two_table_scale_bridge():
    # the classic table is the framework's own relative normalisation:
    # C_Dy = -100 by definition; the closed form supplies the principled
    # converter between the two scales (see the module docstring)
    assert bleaney.BLEANEY_CJ_CLASSIC["Dy"] == -100.0
    ratios = [
        bleaney.BLEANEY_CJ_MODERN[ion] / bleaney.BLEANEY_CJ_CLASSIC[ion]
        for ion in bleaney.BLEANEY_CJ_MODERN
    ]
    assert all(1.75 < ratio < 1.85 for ratio in ratios)
    assert sum(ratios) / len(ratios) == pytest.approx(1.806, abs=0.005)


def test_the_closed_form_reproduces_the_published_constants():
    """C_J = g_J^2 theta_2 J(J+1)(4J(J+1)-3) with exact rationals: the
    well-conditioned ions come back bit-exact, and all six sit inside the
    published values' 3-s.f. printing band (max deviation 0.5, Tb)."""
    closed = {ion: bleaney.c_j_closed_form(ion) for ion in bleaney.BLEANEY_CJ_MODERN}
    assert closed["Tb"] == pytest.approx(-157.5, abs=1e-12)
    assert closed["Dy"] == pytest.approx(-3808 / 21, abs=1e-12)
    assert closed["Ho"] == pytest.approx(-71.25, abs=1e-12)
    assert closed["Er"] == pytest.approx(58.752, abs=1e-12)
    assert closed["Tm"] == pytest.approx(18865 / 198, abs=1e-12)
    assert closed["Yb"] == pytest.approx(1920 / 49, abs=1e-12)
    for ion, published in bleaney.BLEANEY_CJ_MODERN.items():
        assert abs(closed[ion] - published) <= 0.51  # the printed 3 s.f. band
    with pytest.raises(BleaneyError, match="no closed-form inputs"):
        bleaney.c_j_closed_form("Gd")


def test_the_comparator_block_carries_the_closed_form_line():
    body = "\n".join(bleaney.comparator_lines("Dy", 1000.0, 300.0))
    assert "closed form" in body
    assert "-181.33" in body  # the 2-dp closed-form value
    assert "C_J = -181.0" in body  # the shipped constant


def test_the_rhombic_relation_and_the_refusals():
    assert bleaney.chi_rh("Dy", 200.0, 300.0) / 1e-32 == pytest.approx(
        15.1008003, abs=1e-6
    )
    # both relations are linear in B and share the prefactor structure
    assert bleaney.chi_rh("Dy", 3.0, 300.0) == pytest.approx(
        bleaney.chi_ax("Dy", 1.0, 300.0) / 10.0, rel=1e-12
    )
    with pytest.raises(BleaneyError, match="no modern Bleaney constant"):
        bleaney.chi_ax("Gd", 100.0, 300.0)
    with pytest.raises(BleaneyError, match="not positive"):
        bleaney.chi_ax("Dy", 100.0, 0.0)
    assert "kT" in " ".join(bleaney.assumptions())
    kinds = {item.kind for item in bleaney.evidence()}
    assert EVIDENCE_LITERATURE in kinds and EVIDENCE_MEASURED in kinds
