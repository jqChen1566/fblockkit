"""Tests for the relativistic tier tree and the two-component SO-CASSCF template.

The anchors: the tier tree's rows carry the source's own numbers (the Nd3+
increment of 1.6-4.9 cm^-1 for the f block against ~25% of the splitting for
light molecules, the QED scale, the amfX2C warning), and the template's state
count is *computed* -- the ground-term manifold (2S+1)(2L+1) from the Hund term
-- which for Nd3+ must give the source's 52 states.
"""

from __future__ import annotations

import io

import pytest

from fblockkit.recipe import plan_relativistic, plan_so_casscf
from fblockkit.recipe.relativistic import RelativisticError
from fblockkit.ui import Session
from fblockkit.ui.handlers import HANDLERS


# --- the tier tree ------------------------------------------------------------


def test_the_tree_picks_the_tier_by_need():
    scalar = plan_relativistic(needs_soc=False, f_block=True)
    assert scalar.tier.startswith("SFX2C-1e")
    light = plan_relativistic(needs_soc=True, f_block=False)
    assert light.tier == "X2Ccorr"
    heavy = plan_relativistic(needs_soc=True, f_block=True)
    assert heavy.tier.startswith("X2CAMF")
    core = plan_relativistic(needs_soc=True, f_block=False, core_spectroscopy=True)
    assert "QED" in core.tier


def test_the_f_block_case_is_the_avoid_over_configuration_one():
    plan = plan_relativistic(needs_soc=True, f_block=True)
    text = " ".join(plan.notes)
    assert "1.6-4.9 cm^-1" in text
    assert "avoid-over-configuration" in text


def test_amfx2c_is_ruled_out_in_every_row():
    for needs_soc in (True, False):
        plan = plan_relativistic(needs_soc=needs_soc, f_block=True)
        assert any("amfX2C" in note for note in plan.notes)


def test_non_boolean_flags_are_refused():
    with pytest.raises(RelativisticError, match="booleans"):
        plan_relativistic(needs_soc="yes", f_block=True)


# --- the SO-CASSCF template ---------------------------------------------------


def test_the_state_window_is_the_whole_ground_term_manifold():
    nd = plan_so_casscf("Nd")
    assert (nd.n_f_electrons, nd.s, nd.l_value, nd.n_states) == (3, 1.5, 6.0, 52)
    # the same rule sizes every Ln3+ centre
    dy = plan_so_casscf("Dy")
    assert (dy.n_f_electrons, dy.s, dy.l_value, dy.n_states) == (9, 2.5, 5.0, 66)
    er = plan_so_casscf("Er")
    assert (er.n_f_electrons, er.n_states) == (11, 52)


def test_the_template_states_the_basis_tiers_and_the_sensitivity_check():
    nd = plan_so_casscf("Nd")
    text = " ".join(nd.basis_notes + nd.numerical_notes + nd.outputs + nd.checks)
    assert "Dyall VTZ-SO" in text
    assert "Cholesky" in text
    assert "Kramers doublet" in text
    assert "centre of gravity" in text
    assert "4349.4" in text and "1850.3" in text  # the make-or-break example
    assert "414I" in text or "4I9/2" in text or "4I" in text  # the J-multiplet grouping


def test_a_non_f_block_element_is_refused_with_a_next_step():
    with pytest.raises(RelativisticError, match="not an f-block element"):
        plan_so_casscf("Fe")
    with pytest.raises(RelativisticError, match="not an element symbol"):
        plan_so_casscf("Zz")
    with pytest.raises(RelativisticError, match="empty f shell"):
        plan_so_casscf("La")


# --- the menu-3 wiring --------------------------------------------------------


def _run_menu_three(tmp_path, elements, targets, spin="2"):
    xyz = tmp_path / "atom.xyz"
    xyz.write_text("1\nstructure\nCe 0 0 0\n", encoding="utf-8")
    out = io.StringIO()
    Session(
        lines=[
            "3", elements, "0", spin, "3", targets,
            str(xyz), "1", "1,7,2,1", "default", "", "0",
        ],
        out=out,
    ).run(HANDLERS)
    return out.getvalue()


def test_menu_three_adds_the_tier_when_the_profile_needs_it(tmp_path):
    text = _run_menu_three(tmp_path, "Ce,O,H", "magnetic")
    assert "Relativistic tier:" in text
    assert "X2CAMF" in text or "X2Ccorr" in text


def test_menu_three_stays_quiet_for_a_light_energy_profile(tmp_path):
    text = _run_menu_three(tmp_path, "C,H", "energy", spin="1")
    assert "Relativistic tier:" not in text
