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

SYMM_PROBE_TEXT = """&POLY_ANISO

NNEQ
  1  T
  2
  2

SYMM
  2
  1  0  0
  0  1  0
  0  0  1
  -1  0  0
  0  -1  0
  0  0  -1

PAIR
  1
  1 2 0.1

End of Input
"""

LIN3_PROBE_TEXT = """&POLY_ANISO

NNEQ
  2  T
  1  1
  2  2

LIN3
  1
  1 2 0.1 0.1 0.1

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
        equivalent_centres=[1, 1], spin_orbit_states=[2, 2], pairs=[(1, 2, -0.5)]
    )
    assert "COOR" not in plan.text and "TINT" not in plan.text
    assert "SYMM" not in plan.text
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


def test_a_wrong_pair_arity_is_refused_per_model():
    with pytest.raises(poly_aniso.PolyAnisoPlanError, match=r"reads 'i j J' \(three"):
        _probe_plan(pairs=[(1, 2, 0.1, 0.2)])
    with pytest.raises(poly_aniso.PolyAnisoPlanError, match=r"reads 'i j Jx Jy Jz' \(five"):
        _probe_plan(pairs=[(1, 2, 0.1)], pair_model="lin3")


def test_the_lin9_form_is_refused_with_the_measured_engine_defect():
    with pytest.raises(poly_aniso.PolyAnisoPlanError, match="input_process.f90"):
        _probe_plan(pair_model="lin9")
    with pytest.raises(poly_aniso.PolyAnisoPlanError, match="unknown pair model"):
        _probe_plan(pair_model="lin1")


def test_the_missing_symm_is_refused_with_the_measured_reason():
    """A type with several equivalent centres needs its rotation matrices:
    the driver's own check exists (measured) but it still exits 0, so the
    refusal is client-side."""
    with pytest.raises(poly_aniso.PolyAnisoPlanError, match="SYMM is mandatory"):
        poly_aniso.build_input(
            equivalent_centres=[2], spin_orbit_states=[2], pairs=[(1, 2, 0.1)]
        )


def test_the_symm_shape_is_validated():
    identity = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
    inversion = tuple(-value for value in identity)
    with pytest.raises(poly_aniso.PolyAnisoPlanError, match="matrix list"):
        _probe_plan(symmetry=[(identity, identity)])  # 2 types, only 1 list
    with pytest.raises(poly_aniso.PolyAnisoPlanError, match=r"rotation matrix\(es\)"):
        _probe_plan(
            equivalent_centres=[2], spin_orbit_states=[2], pairs=[(1, 2, 0.1)],
            symmetry=[(identity,)],  # 2 centres, 1 matrix
        )
    with pytest.raises(poly_aniso.PolyAnisoPlanError, match="nine"):
        _probe_plan(
            equivalent_centres=[1, 1], pairs=[(1, 2, 0.1)],
            symmetry=[((1.0, 0.0, 0.0),), (identity,)],
        )
    # the accepted shape: identity + inversion for the two-centre type
    plan = poly_aniso.build_input(
        equivalent_centres=[2], spin_orbit_states=[2], pairs=[(1, 2, 0.1)],
        symmetry=[(identity, inversion)],
    )
    assert "SYMM\n  2\n  1  0  0\n  0  1  0\n  0  0  1\n  -1  0  0\n  0  -1  0\n  0  0  -1\n" in plan.text


def test_the_symm_variant_writes_the_probe_bytes():
    """The input frozen as fixtures/poly_aniso/symm_probe.polyinp: the
    one-type two-centre cluster with identity+inversion matrices, accepted
    by the driver end to end (rc = 0)."""
    identity = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
    inversion = tuple(-value for value in identity)
    plan = poly_aniso.build_input(
        equivalent_centres=[2], spin_orbit_states=[2], pairs=[(1, 2, 0.1)],
        symmetry=[(identity, inversion)],
    )
    assert plan.text == SYMM_PROBE_TEXT
    assert plan.symmetry == ((identity, inversion),)
    assert plan.pair_model == "lines"


def test_the_lin3_variant_writes_the_probe_bytes():
    """The input frozen as fixtures/poly_aniso/lin3_probe.polyinp: the same
    two-type cluster as the byte anchor, with the axis-diagonal pair form."""
    plan = poly_aniso.build_input(
        equivalent_centres=[1, 1], spin_orbit_states=[2, 2],
        pairs=[(1, 2, 0.1, 0.1, 0.1)], pair_model="lin3",
    )
    assert plan.text == LIN3_PROBE_TEXT
    assert plan.pair_model == "lin3"


def test_the_plan_render_states_the_checklist_and_boundaries():
    body = poly_aniso.render_plan(_probe_plan(), path="poly_aniso.input")
    assert "centre types: 2" in body
    assert "exchange basis size 4" in body
    assert "otool_poly_aniso" in body
    assert "aniso_1.input" in body
    assert "J values, the coordinates and the SYMM rotation matrices are the caller's" in body
    assert "SYMM is mandatory" in body
    assert "LIN9 form is a documented termination" in body
    lin3_body = poly_aniso.render_plan(
        poly_aniso.build_input(
            equivalent_centres=[1, 1], spin_orbit_states=[2, 2],
            pairs=[(1, 2, 0.1, 0.2, 0.3)], pair_model="lin3",
        ),
        path="poly_aniso.input",
    )
    assert "LIN3 axis-diagonal" in lin3_body and "[0.1 0.2 0.3]" in lin3_body


def test_evidence_states_the_format_source_and_the_acceptance():
    text = " ".join(
        entry.ref + " " + entry.text + " " + entry.url + " " + entry.bibkey
        for entry in poly_aniso.evidence()
    )
    assert "7.18" in text
    assert "aniso_1.input" in text
    assert "byte for byte" in text
    assert "SYMM is mandatory" in text
    assert "LIN9" in text
