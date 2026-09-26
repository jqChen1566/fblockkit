"""Diagnosis-layer tests: rule evaluation, cross-level consistency, report rendering.

Fixtures: fixtures/orca/ (real ORCA 6.1.1 output) and fixtures/literature/pucl3_s18.json
(literature Table S18; provenance in that directory's README).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from fblockkit.diagnosis import (
    CrossLevelError,
    DiagnosisError,
    build_report,
    cross_level_check,
    diagnose,
    diagnose_facts,
    load_records,
    render_markdown,
    to_markdown,
)
from fblockkit.knowledge.loader import RuleError, load_rules
from fblockkit.knowledge.models import SEVERITIES
from fblockkit.parsers import parse_auto

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"
PUCL3 = Path(__file__).resolve().parents[1] / "fixtures" / "literature" / "pucl3_s18.json"


def _findings(name: str) -> dict[str, str]:
    """Fixture -> {rule_id: severity}."""
    return {f.rule_id: f.severity for f in diagnose(parse_auto(FIXTURES / name))}


# --- shape of the rule table -------------------------------------------------


def test_all_diagnosis_rules_have_severity():
    rules = load_rules()
    diagnosis = [r for r in rules if r.kind == "diagnosis"]
    assert diagnosis
    for rule in diagnosis:
        assert rule.severity in SEVERITIES, rule.id
    for rule in rules:
        if rule.kind == "recipe":
            assert not rule.severity, rule.id


def test_diagnosis_rule_without_severity_rejected(tmp_path):
    (tmp_path / "bad.yaml").write_text(
        "- id: X1\n  kind: diagnosis\n  title: t\n"
        "  condition: {all: []}\n  action: a\n"
        "  evidence:\n    - {kind: measured, text: t, ref: r}\n",
        encoding="utf-8",
    )
    with pytest.raises(RuleError, match="severity"):
        load_rules(tmp_path)


def test_recipe_rule_with_severity_rejected(tmp_path):
    (tmp_path / "bad.yaml").write_text(
        "- id: X1\n  kind: recipe\n  title: t\n  severity: warn\n"
        "  condition: {all: []}\n  action: a\n"
        "  evidence:\n    - {kind: measured, text: t, ref: r}\n",
        encoding="utf-8",
    )
    with pytest.raises(RuleError, match="severity"):
        load_rules(tmp_path)


def test_engine_rejects_recipe_rules():
    recipe = next(rule for rule in load_rules() if rule.kind == "recipe")
    with pytest.raises(DiagnosisError, match="diagnosis rules"):
        diagnose_facts({}, rules=[recipe])


# --- criterion behaviour on the fixtures -------------------------------------


def test_caspt2_run_findings():
    findings = _findings("n2_caspt2.out")
    # the D3 mandatory-check note (info) and the CASSCF energy-criterion note
    assert findings["DG-CASPT2-WEIGHTS"] == "info"
    assert findings["DG-CASSCF-ENERGY-ONLY-CONVERGENCE"] == "warn"
    # this job has weight 0.9347 >= 0.9 and smallest denominator 0.2025 >= 0.01,
    # so neither threshold rule may fire
    assert "DG-CASPT2-LOW-WEIGHT" not in findings
    assert "DG-CASPT2-SMALL-DENOMINATOR" not in findings
    # the smallest active occupation 0.00208 < 0.02 -> a nearly empty orbital
    assert findings["DG-ACTIVE-SPACE-NEAR-EMPTY"] == "warn"


def test_caspt2_weight_thresholds_fire():
    """The threshold rules are checked with synthetic facts (the fixture's own weight is healthy)."""
    low_weight = diagnose_facts(
        {"caspt2_present": True, "caspt2_min_reference_weight": 0.85}
    )
    assert {f.rule_id for f in low_weight} >= {"DG-CASPT2-LOW-WEIGHT"}
    small_denom = diagnose_facts(
        {"caspt2_present": True, "caspt2_min_denominator": 0.005}
    )
    assert {f.rule_id for f in small_denom} >= {"DG-CASPT2-SMALL-DENOMINATOR"}


def test_nevpt2_run_findings():
    findings = _findings("n2_casscf_nevpt2.out")
    # all class contributions are negative (V1_i/Vm1_a) -> no sign of an intruder state
    assert "DG-NEVPT2-FALSE-INTRUDER" not in findings
    assert findings["DG-CASSCF-ENERGY-ONLY-CONVERGENCE"] == "warn"


def test_nevpt2_false_intruder_fires_on_positive_class():
    findings = diagnose_facts({"nevpt2_max_hole_particle": 0.012})
    assert {f.rule_id for f in findings} >= {"DG-NEVPT2-FALSE-INTRUDER"}


def test_casscf_not_converged_has_refusal():
    findings = diagnose_facts({"casscf_present": True, "casscf_converged": False})
    by_id = {f.rule_id: f for f in findings}
    assert by_id["DG-CASSCF-NOT-CONVERGED"].severity == "error"
    assert by_id["DG-CASSCF-NOT-CONVERGED"].refusals


def test_soc_run_findings():
    findings = _findings("co_plus_soc.out")
    assert findings["DG-SOC-STATE-IDENTITY"] == "info"
    # this job converged by the gradient criterion -> the energy-criterion rule must not fire
    assert "DG-CASSCF-ENERGY-ONLY-CONVERGENCE" not in findings


def test_clean_dft_run_has_no_findings():
    assert diagnose(parse_auto(FIXTURES / "fblock_dft_la_complex.out")) == ()


def test_aborted_runs_reported():
    assert _findings("scf_noconv.out")["DG-RUN-ABORTED"] == "error"
    assert _findings("fblock_dft_gd_crash.out")["DG-RUN-ABORTED"] == "error"
    # an unconverged (aborted) run must not fire the pseudo-convergence rule
    # (its precondition is "convergence was reported")
    assert "DG-SCF-PSEUDO-CONVERGENCE" not in _findings("scf_noconv.out")


def test_pseudo_convergence_fires_only_when_converged():
    findings = diagnose_facts({"scf_converged": True, "scf_cycles": 2})
    assert {f.rule_id for f in findings} >= {"DG-SCF-PSEUDO-CONVERGENCE"}


def test_generated_input_roundtrip_diagnosis():
    """Closed-loop fixture (generated by the recipe layer -> a real ORCA run): only the
    "nearly empty orbital" note should fire.

    For an f1 system, putting the whole f shell into the window is a routine choice in
    magnetic/spectroscopic studies, so this conclusion is a warn-level "please confirm"
    (convention in the "known-behaviour notes" of fixtures/orca/README.md).
    """
    findings = _findings("generated_ce3_sarc2.out")
    assert set(findings) == {"DG-ACTIVE-SPACE-NEAR-EMPTY"}


# --- cross-level solution consistency (literature fixture PuCl3 Table S18) ---


def test_pucl3_fixture_matches_literature_spans():
    records = load_records(PUCL3)
    assert len(records) == 21
    hf = [r["energies"]["HF"] for r in records]
    ccsd = [r["energies"]["CCSD(T)"] for r in records]
    # in step with the full-table span quoted in the literature review §4.6
    assert min(hf) == pytest.approx(-1930.87769)
    assert max(hf) == pytest.approx(-1930.81091)
    assert min(ccsd) == pytest.approx(-1932.97951)
    assert max(ccsd) == pytest.approx(-1932.93063)


def test_cross_level_check_fires_on_pucl3():
    findings = cross_level_check(load_records(PUCL3))
    by_id = {f.rule_id: f for f in findings}
    assert set(by_id) == {"XL-CROSS-LEVEL-BEST", "XL-CROSS-LEVEL-TOP-N"}
    # the two levels disagree on the best solution: HF favours No.3, CCSD(T) favours
    # No.1, a penalty of about 0.06 kcal/mol
    assert "No.3" in by_id["XL-CROSS-LEVEL-BEST"].message
    assert "No.1" in by_id["XL-CROSS-LEVEL-BEST"].message
    assert "0.06 kcal/mol" in by_id["XL-CROSS-LEVEL-BEST"].message
    # the second-lowest HF solution No.10 (-1930.87697) drops out of the top three at
    # the CCSD(T) level
    assert "No.10" in by_id["XL-CROSS-LEVEL-TOP-N"].message
    # the provenance is literature (with its DOI)
    assert by_id["XL-CROSS-LEVEL-BEST"].evidence[0].kind == "literature"
    assert "10.1021/acs.jctc.4c01189" in by_id["XL-CROSS-LEVEL-BEST"].evidence[0].ref


def test_cross_level_check_silent_when_order_preserved():
    records = [
        {"label": "A", "energies": {"HF": -10.0, "CCSD(T)": -20.0}},
        {"label": "B", "energies": {"HF": -9.0, "CCSD(T)": -19.0}},
        {"label": "C", "energies": {"HF": -8.0, "CCSD(T)": -18.0}},
    ]
    assert cross_level_check(records) == ()


def test_cross_level_check_rejects_missing_level():
    with pytest.raises(CrossLevelError, match="has no"):
        cross_level_check([{"label": "A", "energies": {"HF": -1.0}}])


# --- report rendering --------------------------------------------------------


def test_report_rendering_from_fixture():
    findings = diagnose(parse_auto(FIXTURES / "n2_caspt2.out"))
    report = build_report(findings, subject="n2_caspt2.out")
    markdown = to_markdown(report)
    assert "# fBlockKit report" in markdown
    assert "Subject: n2_caspt2.out" in markdown
    assert "[Info]" in markdown and "[Warning]" in markdown
    assert "DG-CASPT2-WEIGHTS" in markdown
    assert "## Provenance" in markdown
    assert "ORCA 6.1 manual §3.17" in markdown
    # severity ordering: warning before info
    assert markdown.index("[Warning]") < markdown.index("[Info]")


def test_render_markdown_empty_findings():
    text = render_markdown((), subject="clean.out")
    assert "No diagnostic rule matched" in text


def test_report_combines_analysis_sections():
    """Analysis and diagnosis combine into one report: the analysis sections sit after
    the summary and before the findings."""
    from fblockkit.analysis import run_all

    result = parse_auto(FIXTURES / "n2_casscf_orbcomp.out")
    sections = run_all(result)
    # A1 (the composition table is in this fixture) + A2 (entropy) + A3 (the MR panel
    # accepts any CASSCF output with active occupations)
    assert len(sections) == 3
    markdown = to_markdown(build_report(diagnose(result), sections=sections, subject="n2"))
    assert markdown.index("## Summary") < markdown.index("## A2") < markdown.index("## Findings")
    assert "## A1" in markdown and "## A3" in markdown
