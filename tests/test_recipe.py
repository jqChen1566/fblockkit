"""Recipe-layer tests: G1 basis-set/ECP recommendation, G4 convergence template,
ORCA input rendering.

Discipline: the basis table is as strict as the rule table -- an entry without provenance
must be rejected; a generated input must be ASCII and directly runnable (a missing
auxiliary basis for TRAH, say, must raise rather than produce a broken input).
"""

from __future__ import annotations

import pytest

from fblockkit.knowledge.loader import RuleError, load_rules
from fblockkit.knowledge.models import SystemProfile
from fblockkit.recipe import (
    BasisDataError,
    facts_from_profile,
    load_basis_entries,
    plan_convergence,
    recommend,
    recommend_basis_ecp,
    recommend_with_advice,
    render_orca_input,
    run_guidance,
)
from fblockkit.recipe.render import RenderError


def _profile(**kwargs) -> SystemProfile:
    base = dict(elements=("Ce", "O", "H"), charge=-2, f_count=1, targets=("energy",))
    base.update(kwargs)
    return SystemProfile(**base)


# --- G1 data table -----------------------------------------------------------


def test_basis_table_loads_with_evidence():
    entries = load_basis_entries()
    assert len(entries) >= 7
    ids = [entry.id for entry in entries]
    assert len(ids) == len(set(ids))
    for entry in entries:
        assert entry.evidence, entry.id
        assert entry.note.strip(), entry.id


def test_basis_entry_without_evidence_rejected(tmp_path):
    bad = tmp_path / "bad.yaml"
    bad.write_text(
        "- id: X1\n  kind: recommend\n  elements: La-Lu\n"
        "  basis: SARC2-DKH-QZVP\n  note: n\n  evidence: []\n",
        encoding="utf-8",
    )
    with pytest.raises(BasisDataError, match="provenance"):
        load_basis_entries(bad)


def test_basis_entry_unknown_element_rejected(tmp_path):
    bad = tmp_path / "bad.yaml"
    bad.write_text(
        "- id: X1\n  kind: recommend\n  elements: Xx-Yy\n  basis: b\n  note: n\n"
        "  evidence:\n    - {kind: manual, text: t, ref: r}\n",
        encoding="utf-8",
    )
    with pytest.raises(BasisDataError, match="unknown element"):
        load_basis_entries(bad)


def test_recommend_basis_for_ce3():
    advice = recommend_basis_ecp(_profile())
    ids = {entry.id for entry in advice.recommendations}
    assert advice.applicable
    assert "BE-SARC2-LN" in ids  # first choice for lanthanide wavefunction methods
    assert "BE-DEF2-LN" in ids
    assert "BE-5F-IN-CORE" not in ids  # Pa-Lr only


def test_recommend_basis_for_pu3_includes_5f_in_core_with_refusals():
    advice = recommend_basis_ecp(
        _profile(elements=("Pu", "Cl"), metal_valence=3, f_count=5)
    )
    rec_ids = {entry.id for entry in advice.recommendations}
    assert {"BE-SARC-DKH-AN", "BE-ANO-RCC", "BE-STUTTGART-ECP"} <= rec_ids
    assert {entry.id for entry in advice.refusals} == {"BE-5F-IN-CORE"}
    refusals = advice.refusal_texts()
    assert any("f–f" in text for text in refusals)
    assert any("SOC" in text for text in refusals)
    assert {entry.id for entry in advice.cautions} == {"BE-DEF2-AN-CAUTION"}


def test_recommend_basis_pu4_excludes_5f_in_core():
    advice = recommend_basis_ecp(
        _profile(elements=("Pu", "Cl"), metal_valence=4, f_count=4)
    )
    assert "BE-5F-IN-CORE" not in {entry.id for entry in advice.refusals}


def test_recommend_basis_without_f_element():
    advice = recommend_basis_ecp(SystemProfile(elements=("C", "H"), targets=("energy",)))
    assert advice.applicable is False
    assert advice.recommendations == ()
    assert any("contains no lanthanide/actinide" in note for note in advice.notes)


# --- recipe recommendation ---------------------------------------------------


def test_facts_from_profile():
    facts = facts_from_profile(_profile(targets=("energy", "geometry")))
    assert facts["f_block"] is True
    assert facts["geometry_task"] is True
    assert facts["soc_task"] is False
    facts = facts_from_profile(_profile(elements=("Gd",), targets=("magnetic",), f_count=7))
    assert facts["soc_task"] is True


def test_recommend_energy_profile():
    rec = recommend(_profile())
    assert "CASSCF + SC-NEVPT2" in rec.method_chain
    assert rec.key_settings  # the standing guidance from A3 / A4 and so on
    assert rec.evidence
    assert rec.basis_ecp


def test_recommend_geometry_profile_refuses_pt2_in_loop():
    rec = recommend(_profile(targets=("energy", "geometry")))
    assert "CASSCF (geometry/frequency reference)" in rec.method_chain
    assert any("optimisation loop" in text for text in rec.refusals)


def test_recommend_rejects_diagnosis_rules():
    diagnosis = next(rule for rule in load_rules() if rule.kind == "diagnosis")
    with pytest.raises(ValueError, match="recipe rules"):
        recommend(_profile(), rules=[diagnosis])


def test_diagnosis_rule_with_method_rejected(tmp_path):
    (tmp_path / "bad.yaml").write_text(
        "- id: X1\n  kind: diagnosis\n  severity: warn\n  method: NEVPT2\n  title: t\n"
        "  condition: {all: []}\n  action: a\n"
        "  evidence:\n    - {kind: measured, text: t, ref: r}\n",
        encoding="utf-8",
    )
    with pytest.raises(RuleError, match="method"):
        load_rules(tmp_path)


# --- G4 convergence template -------------------------------------------------


def test_plan_convergence_default_first():
    plan = plan_convergence(difficulty="default")
    assert plan.simple_keywords == ()
    assert plan.needs_auxiliary is False
    assert any("default" in note for note in plan.notes)


def test_plan_convergence_difficult_uses_trah():
    plan = plan_convergence(difficulty="difficult", pt2=True)
    assert "TRAH" in plan.simple_keywords
    assert plan.needs_auxiliary is True
    assert any("NEVPT2" in line for line in plan.scf_block)
    assert any("orbital convergence" in note for note in plan.notes)


def test_plan_convergence_rejects_unknown_difficulty():
    with pytest.raises(ValueError, match="difficulty"):
        plan_convergence(difficulty="whatever")


# --- input rendering ---------------------------------------------------------


def _render(**kwargs) -> str:
    rec, advice = recommend_with_advice(_profile(elements=("Ce", "O", "H")))
    entry = next(e for e in advice.recommendations if e.id == "BE-SARC2-LN")
    params = dict(
        geometry="Ce 0.0 0.0 0.0\nO 0.0 0.0 2.0",
        charge=1,
        mult=2,
        basis_entry=entry,
        method_keywords=["TightSCF", "NEVPT2"],
        casscf={"nel": 1, "norb": 7, "mult": 2, "nroots": 1},
        convergence=plan_convergence(pt2=True),
    )
    params.update(kwargs)
    return render_orca_input(rec, **params)


def test_render_orca_input_basic():
    text = _render()
    text.encode("ascii")  # the input file must be ASCII
    assert "! DKH2 SARC2-DKH-QZVP TightSCF NEVPT2" in text
    assert "%casscf" in text and "nel 1" in text and "norb 7" in text
    assert "* xyz 1 2" in text
    assert "Ce 0.0 0.0 0.0" in text
    assert text.rstrip().endswith("*")


def test_render_requires_auxiliary_for_trah():
    with pytest.raises(RenderError, match="auxiliary"):
        _render(convergence=plan_convergence(difficulty="difficult"))


def test_render_with_auxiliary_emits_trah():
    text = _render(
        convergence=plan_convergence(difficulty="difficult"),
        auxiliary="SARC2-DKH-QZVP/JK",
    )
    assert "TRAH" in text
    assert "SARC2-DKH-QZVP/JK" in text


def test_render_rejects_bad_geometry():
    with pytest.raises(RenderError, match="element symbol"):
        _render(geometry="0.0 0.0 0.0")


def test_render_rejects_incomplete_casscf():
    with pytest.raises(RenderError, match="nroots"):
        _render(casscf={"nel": 1, "norb": 7, "mult": 2})


def test_run_guidance_contains_diagnostics_and_refusals():
    rec = recommend(_profile(elements=("Pu", "Cl"), metal_valence=3, f_count=5))
    text = run_guidance(rec, plan_convergence())
    assert "diagnosis layer" in text
    assert "5f" in text or "SOC" in text  # the refusals travel with it
