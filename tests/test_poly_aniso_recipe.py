"""Tests for the POLY_ANISO input writer (recipe/poly_aniso.py).

The byte anchor is the fixture's two-centre probe plan: the generated input
that reproduces ``fixtures/poly_aniso/two_center_probe.out`` byte for byte
when run with the fixture's ``aniso_1/aniso_2.input`` files (measured on
6.1.1, recorded in the module evidence and the fixture README).
"""

from __future__ import annotations

import pytest

from fblockkit.recipe import poly_aniso

PROBE_TEXT = """&POLY_ANISO

NNEQ
  2  T
  1  1
  2  2

PAIR
  1
  1 2 0.1

COOR
  0 0 0
  0 0 3.7

TINT
  0 300 101

End of Input
"""


def _probe_plan(**overrides):
    arguments = dict(
        equivalent_centres=[1, 1],
        spin_orbit_states=[2, 2],
        pairs=[(1, 2, 0.1)],
        coordinates=[(0.0, 0.0, 0.0), (0.0, 0.0, 3.7)],
        temperature_grid=(0.0, 300.0, 101),
    )
    arguments.update(overrides)
    return poly_aniso.build_input(**arguments)


def test_the_probe_plan_writes_the_frozen_bytes():
    """The plan behind the acceptance run: the exact input text."""
    plan = _probe_plan()
    assert plan.text == PROBE_TEXT
    assert plan.total_centres == 2
    assert plan.exchange_basis == 4


def test_optional_blocks_can_be_omitted():
    plan = poly_aniso.build_input(
        equivalent_centres=[2], spin_orbit_states=[2], pairs=[(1, 2, -0.5)]
    )
    assert "COOR" not in plan.text and "TINT" not in plan.text
    assert "PAIR\n  1\n  1 2 -0.5\n" in plan.text
    assert plan.text.startswith("&POLY_ANISO\n") and plan.text.endswith("End of Input\n")


def test_the_refusals():
    with pytest.raises(poly_aniso.PolyAnisoPlanError, match="one number per centre type"):
        poly_aniso.build_input(equivalent_centres=[1, 1], spin_orbit_states=[2], pairs=[(1, 2, 0.1)])
    with pytest.raises(poly_aniso.PolyAnisoPlanError, match="up to 6"):
        _probe_plan(equivalent_centres=[1] * 7, spin_orbit_states=[2] * 7, pairs=[(1, 2, 0.1)])
    with pytest.raises(poly_aniso.PolyAnisoPlanError, match="at least one centre"):
        _probe_plan(equivalent_centres=[0, 1])
    with pytest.raises(poly_aniso.PolyAnisoPlanError, match="at least two"):
        _probe_plan(spin_orbit_states=[1, 2])
    with pytest.raises(poly_aniso.PolyAnisoPlanError, match="no coupled pairs"):
        _probe_plan(pairs=[])
    with pytest.raises(poly_aniso.PolyAnisoPlanError, match="outside the cluster"):
        _probe_plan(pairs=[(1, 3, 0.1)])
    with pytest.raises(poly_aniso.PolyAnisoPlanError, match="with itself"):
        _probe_plan(pairs=[(1, 1, 0.1)])
    with pytest.raises(poly_aniso.PolyAnisoPlanError, match="appears twice"):
        _probe_plan(pairs=[(1, 2, 0.1), (2, 1, 0.2)])
    with pytest.raises(poly_aniso.PolyAnisoPlanError, match="one 'x y z' line"):
        _probe_plan(coordinates=[(0.0, 0.0, 0.0)])
    with pytest.raises(poly_aniso.PolyAnisoPlanError, match="not increasing"):
        _probe_plan(temperature_grid=(300.0, 0.0, 101))
    with pytest.raises(poly_aniso.PolyAnisoPlanError, match="not increasing"):
        _probe_plan(temperature_grid=(0.0, 300.0, 1))


def test_the_plan_render_states_the_checklist_and_boundaries():
    body = poly_aniso.render_plan(_probe_plan(), path="poly_aniso.input")
    assert "centre types: 2" in body
    assert "exchange basis size 4" in body
    assert "otool_poly_aniso" in body
    assert "aniso_1.input" in body
    assert "J values and the coordinates are the caller's" in body
    assert "LIN3/LIN9" in body


def test_evidence_states_the_format_source_and_the_acceptance():
    text = " ".join(
        entry.ref + " " + entry.text + " " + entry.url + " " + entry.bibkey
        for entry in poly_aniso.evidence()
    )
    assert "7.18" in text
    assert "aniso_1.input" in text
    assert "byte for byte" in text
