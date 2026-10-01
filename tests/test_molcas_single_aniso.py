"""Tests for the OpenMolcas SINGLE_ANISO reader (parsers/molcas_single_aniso.py).

The real-data anchor is the Dy(III) 4f9 single-ion chain
(fixtures/openmolcas/dy_smoke.*): one pseudospin multiplet of effective
S~ = 1/2 built from the two degenerate lowest spin-orbit states, with the
measured g values and main axes, and the sign check of the g-value product.

The last test closes the loop into menu 16: the payload built from the parsed
multiplets is fed through the criterion machinery unchanged, so the engine
output and the hand-written JSON table enter by one door.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from fblockkit.analysis import magnetic_doublets
from fblockkit.parsers.molcas_single_aniso import (
    MolcasAnisoError,
    doublets_payload,
    parse_single_aniso,
    parse_susceptibility,
)

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "openmolcas"
OUT = FIXTURES / "dy_smoke.out"


@pytest.fixture(scope="module")
def text() -> str:
    return OUT.read_text(encoding="utf-8")


def test_the_dy_multiplet_parses_with_all_measured_fields(text):
    multiplets = parse_single_aniso(text)
    assert len(multiplets) == 1
    m = multiplets[0]
    assert m.index == 1
    assert m.effective_spin == "1/2"
    assert m.so_states == ((1, 0.0), (2, 0.0))
    assert m.tunnelling_splitting_cm1 == pytest.approx(0.0, abs=1e-12)
    assert m.g_values[0] == pytest.approx(0.91315132997786)
    assert m.g_values[1] == pytest.approx(1.99057019632358)
    assert m.g_values[2] == pytest.approx(13.14431673558191)
    # the main magnetic axes, in the initial Cartesian frame
    assert m.axes[2][0] == pytest.approx(0.28993600214377)
    assert m.axes[2][1] == pytest.approx(-0.09211736592909)
    assert m.axes[2][2] == pytest.approx(0.95260249084031)
    assert m.check_sign == pytest.approx(-19.908976)
    assert m.sign_product == -1


def test_the_payload_orders_g_ascending_and_carries_the_axis(text):
    payload = doublets_payload(parse_single_aniso(text), system="Dy 4f9 smoke")
    assert payload["reference"] == 0
    assert len(payload["doublets"]) == 1
    row = payload["doublets"][0]
    # g ascending (menu 16's g1 <= g2 <= g3 convention); already ascending here
    assert row["g"][0] == pytest.approx(0.91315132997786)
    assert row["g"][2] == pytest.approx(13.14431673558191)
    # the axis travels with the largest g value (gZ's Zm axis)
    assert row["axis3"][2] == pytest.approx(0.95260249084031)
    assert row["energy"] == pytest.approx(0.0)


def test_the_payload_feeds_menu_16_end_to_end(text):
    payload = doublets_payload(parse_single_aniso(text))
    table = magnetic_doublets.parse_doublets(payload)
    report = magnetic_doublets.analyze(table)
    assert report.reference == 0
    assert not report.computed_from_axes  # a single doublet: theta3 is 0 by definition
    verdict = report.verdicts[0]
    assert verdict.theta3_deg == 0.0
    expected = (0.91315132997786 + 1.99057019632358) / 3.0
    assert verdict.g_transverse == pytest.approx(expected)


def test_a_text_without_the_block_is_refused_with_a_next_step():
    with pytest.raises(MolcasAnisoError, match="Next step"):
        parse_single_aniso("hello world, no pseudospin here")


def test_the_susceptibility_table_parses(text):
    table = parse_susceptibility(text)
    assert table.n_points == 301
    assert table.temperature_range_k == (0.0, 300.0)
    assert len(table.points) == 301
    first, last = table.points[0], table.points[-1]
    assert first.temperature_k == pytest.approx(0.0001)
    assert first.chi_t_cm3k_mol == pytest.approx(14.18870104)
    assert last.temperature_k == pytest.approx(300.0)
    assert last.chi_t_cm3k_mol == pytest.approx(14.21767223)
    assert last.chi_cm3_mol == pytest.approx(0.0473922)


def test_a_text_without_the_susceptibility_section_is_refused():
    with pytest.raises(MolcasAnisoError, match="Next step"):
        parse_susceptibility("no susceptibility section in here")


def test_an_incomplete_g_table_is_refused(text):
    # drop the gZ line; the block header survives, so the partial table must be
    # refused rather than read with a missing component
    gz_line = next(
        line for line in text.splitlines() if line.strip().startswith("gZ =")
    )
    mangled = text.replace(gz_line, "")
    with pytest.raises(MolcasAnisoError, match="incomplete g-tensor"):
        parse_single_aniso(mangled)
