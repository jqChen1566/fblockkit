"""Tests for the DMET embedding recipe (recipe/dmet.py) and its menu-6 hook.

The recipe is guidance: the anchors are the source's own settings (the loose
double criterion, the Loewdin split with the f-centre impurity, R = Delta S_E,
SA-CASSCF in the cluster, SOMF spin-orbit, the ANO-RCC basis tiers) and the
project's settled convention for the core-contribution accounting.  The f-count
derivation is pinned against the lanthanide series.
"""

from __future__ import annotations

import io

import pytest

from fblockkit.recipe import plan_dmet, render_dmet
from fblockkit.recipe.dmet import DmetError
from fblockkit.ui import Session
from fblockkit.ui.handlers import HANDLERS


# --- the f-count derivation ---------------------------------------------------


def test_the_f_count_follows_the_tripositive_ion():
    assert plan_dmet("Dy").n_f_electrons == 9
    assert plan_dmet("Er").n_f_electrons == 11
    assert plan_dmet("Ce").n_f_electrons == 1
    assert plan_dmet("Yb").n_f_electrons == 13
    assert plan_dmet("U").n_f_electrons == 3  # actinide: Z - 89
    # the symbol is normalised the way the rest of the toolkit normalises it
    assert plan_dmet("dy").element == "Dy"


def test_a_non_f_block_element_is_refused():
    with pytest.raises(DmetError, match="not an f-block element"):
        plan_dmet("Fe")
    with pytest.raises(DmetError, match="not an element symbol"):
        plan_dmet("Xx")


def test_an_empty_f_shell_is_refused():
    with pytest.raises(DmetError, match="empty f shell"):
        plan_dmet("La")


# --- the content anchors ------------------------------------------------------


def test_the_steps_carry_the_source_settings():
    plan = plan_dmet("Dy")
    text = render_dmet(plan)
    assert "delta E < 1e-4 AND |g| < 0.01" in text  # the loose pair, both criteria
    assert "R = Delta S_E" in text
    assert "CAS(9e,7o)" in text  # the f-count-derived cluster CAS
    assert "SOMF" in text
    assert "Loewdin-orthogon" in text or "Loewdin-orthogonality" in text
    assert "1e-13" in text  # the bath threshold
    for tier in ("ANO-RCC-VTZP", "ANO-RCC-VDZP", "Cholesky"):
        assert tier in text, tier


def test_the_gates_quantify_the_failure_mode():
    text = render_dmet(plan_dmet("Dy"))
    assert "0.007" in text and "2.766" in text  # the source's Delta S_E anchors
    assert "42x" in text and "357.7" in text and "8.6" in text  # the MAE cost
    assert "orbital gradient" in text  # the double criterion


def test_the_boundaries_state_the_accounting_decision_and_the_licence_gate():
    plan = plan_dmet("Dy")
    text = render_dmet(plan)
    assert "Core-contribution accounting" in text
    assert "spin-orbit" in text  # the two accountings are named
    assert "convention-free" in text  # the explicit Fock form
    assert "licence" in text  # liblan is not bundled
    assert "One-shot DMET" in text


def test_evidence_carries_the_source():
    plan = plan_dmet("Dy")
    text = " ".join(entry.ref + " " + entry.text for entry in plan.evidence)
    assert "10.1021/acs.jctc.5c01336" in text
    assert "0.6-7.8" in text


# --- the menu-6 hook -----------------------------------------------------------


def test_menu_six_prints_the_recipe_for_liblan():
    out = io.StringIO()
    Session(lines=["6", "liblan", "0"], out=out).run(HANDLERS)
    text = out.getvalue()
    assert "liblan" in text
    assert "DMET embedding recipe for Dy" in text
    assert "(f-count)e,7o" in text


def test_menu_six_leaves_other_tools_without_a_recipe():
    out = io.StringIO()
    Session(lines=["6", "openmolcas", "0"], out=out).run(HANDLERS)
    text = out.getvalue()
    assert "OpenMolcas" in text
    assert "DMET embedding recipe" not in text
