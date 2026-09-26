"""Knowledge-layer tests: condition evaluation, rule loading and validation
(including the rejection paths).

Discipline: a positive case must be paired with a negative one -- a validator that
rejects a rule without provenance matters as much as one that accepts a valid rule.
"""

from __future__ import annotations

import pytest

from fblockkit.knowledge.conditions import (
    AllOf,
    AnyOf,
    ConditionError,
    Predicate,
    evaluate,
    from_mapping,
)
from fblockkit.knowledge.loader import RuleError, evaluate_rules, load_rules
from fblockkit.knowledge.models import EVIDENCE_KINDS

# --- condition evaluation ---------------------------------------------------


def test_predicate_operators():
    facts = {"n": 5, "kind": "f_block", "targets": ["energy", "magnetic"]}
    assert evaluate(Predicate("n", "eq", 5), facts)
    assert evaluate(Predicate("n", "ne", 4), facts)
    assert evaluate(Predicate("n", "gt", 4), facts)
    assert evaluate(Predicate("n", "ge", 5), facts)
    assert evaluate(Predicate("n", "lt", 6), facts)
    assert evaluate(Predicate("n", "le", 5), facts)
    assert evaluate(Predicate("kind", "in", ["f_block", "d_block"]), facts)
    assert evaluate(Predicate("kind", "not_in", ["main_group"]), facts)
    assert evaluate(Predicate("targets", "contains", "energy"), facts)
    assert not evaluate(Predicate("targets", "contains", "spectra"), facts)


def test_missing_field_is_none_not_error():
    assert not evaluate(Predicate("absent", "eq", 1), {})
    assert evaluate(Predicate("absent", "eq", None), {})


def test_all_any_semantics():
    t = Predicate("x", "eq", 1)
    f = Predicate("x", "eq", 2)
    assert evaluate(AllOf((t, t)), {"x": 1})
    assert not evaluate(AllOf((t, f)), {"x": 1})
    assert evaluate(AnyOf((f, t)), {"x": 1})
    assert not evaluate(AnyOf((f, f)), {"x": 1})
    # an empty conjunction is true and an empty disjunction false (used by unconditional rules)
    assert evaluate(AllOf(()), {})
    assert not evaluate(AnyOf(()), {})


def test_from_mapping_forms_and_rejections():
    node = from_mapping({"field": "a", "op": "eq", "value": 1})
    assert isinstance(node, Predicate)
    node = from_mapping({"all": [{"field": "a", "op": "eq", "value": 1}]})
    assert isinstance(node, AllOf)
    with pytest.raises(ConditionError, match="unknown operator"):
        from_mapping({"field": "a", "op": "zzz", "value": 1})
    with pytest.raises(ConditionError, match="condition must be"):
        from_mapping({"foo": 1})


# --- rule loading and validation --------------------------------------------


def test_load_shipped_rules_valid():
    rules = load_rules()
    assert len(rules) >= 12
    ids = [rule.id for rule in rules]
    assert len(ids) == len(set(ids))
    for rule in rules:
        assert rule.evidence, f"{rule.id} has no provenance"
        for item in rule.evidence:
            assert item.kind in EVIDENCE_KINDS
            assert item.text and item.ref


def test_rule_without_evidence_rejected(tmp_path):
    bad = tmp_path / "bad.yaml"
    bad.write_text(
        "- id: X1\n  kind: recipe\n  title: t\n"
        "  condition: {all: []}\n  action: a\n  evidence: []\n",
        encoding="utf-8",
    )
    with pytest.raises(RuleError, match="provenance"):
        load_rules(tmp_path)


def test_duplicate_id_rejected(tmp_path):
    entry = (
        "- id: X1\n  kind: recipe\n  title: t\n"
        "  condition: {all: []}\n  action: a\n"
        "  evidence:\n    - {kind: measured, text: t, ref: r}\n"
    )
    (tmp_path / "a.yaml").write_text(entry, encoding="utf-8")
    (tmp_path / "b.yaml").write_text(entry, encoding="utf-8")
    with pytest.raises(RuleError, match="duplicate"):
        load_rules(tmp_path)


def test_bad_condition_rejected(tmp_path):
    (tmp_path / "bad.yaml").write_text(
        "- id: X1\n  kind: recipe\n  title: t\n"
        "  condition: {field: a, op: zzz, value: 1}\n  action: a\n"
        "  evidence:\n    - {kind: measured, text: t, ref: r}\n",
        encoding="utf-8",
    )
    with pytest.raises(RuleError, match="invalid condition"):
        load_rules(tmp_path)


# --- behaviour of the shipped rules -----------------------------------------


def test_a2_fires_for_f_block_energy():
    rules = load_rules()
    hits = {r.id for r in evaluate_rules(rules, {"f_block": True, "targets": ["energy"]})}
    assert "A2-f-block-energy-nevpt2" in hits
    hits = {r.id for r in evaluate_rules(rules, {"f_block": False, "targets": ["energy"]})}
    assert "A2-f-block-energy-nevpt2" not in hits


def test_b3_fires_at_mid_size():
    rules = load_rules()
    hits = {r.id for r in evaluate_rules(rules, {"active_space_orbitals": 20})}
    assert "B3-large-active-space" in hits
    hits = {r.id for r in evaluate_rules(rules, {"active_space_orbitals": 8})}
    assert "B3-large-active-space" not in hits


def test_diagnosis_rule_fires():
    rules = load_rules()
    hits = {r.id for r in evaluate_rules(rules, {"caspt2_present": True})}
    assert "DG-CASPT2-WEIGHTS" in hits


# --- evaluation of the newer fact fields (vocabulary in rules/README.md) -----


def test_diagnosis_fact_fields_evaluate():
    """Evaluation cases for the fields the parser layer added (scf_converged /
    casscf_* / caspt2_min_*)."""
    facts = {
        "scf_converged": False,
        "casscf_present": True,
        "casscf_converged": True,
        "caspt2_min_reference_weight": 0.86,
        "caspt2_min_denominator": 0.20,
    }
    assert evaluate(Predicate("scf_converged", "eq", False), facts)
    assert evaluate(Predicate("casscf_converged", "eq", True), facts)
    assert evaluate(Predicate("caspt2_min_reference_weight", "lt", 0.9), facts)
    assert evaluate(Predicate("caspt2_min_denominator", "le", 0.25), facts)


def test_frequency_fact_fields_evaluate():
    """Evaluation cases for the frequency / optimisation fields (vocabulary in
    rules/README.md), including the negative-direction contrasts."""
    facts = {
        "ts_optimization": True,
        "optimization_converged": False,
        "frequency_present": True,
        "frequency_after_geometry": True,
        "frequency_imaginary_count": 1,
        "frequency_min_imaginary": -90.48,
    }
    assert evaluate(Predicate("ts_optimization", "eq", True), facts)
    assert evaluate(Predicate("optimization_converged", "eq", False), facts)
    assert evaluate(Predicate("frequency_imaginary_count", "ge", 1), facts)
    assert evaluate(Predicate("frequency_min_imaginary", "gt", -50.0), facts) is False
    # a missing field stays False for ordering and equality comparisons alike
    assert evaluate(Predicate("frequency_min_imaginary", "gt", -50.0), {}) is False
    assert evaluate(Predicate("frequency_imaginary_count", "ge", 1), {}) is False


def test_solution_branch_fact_fields_evaluate():
    """Evaluation cases for the composition-derived fields (vocabulary in
    rules/README.md): the f-block element presence and the largest f weight over the
    active orbitals."""
    facts = {"f_block_element_present": True, "active_f_weight_max": 0.0}
    assert evaluate(Predicate("f_block_element_present", "eq", True), facts)
    assert evaluate(Predicate("active_f_weight_max", "lt", 10.0), facts)
    # without the composition table both fields are absent: neither comparison fires
    assert evaluate(Predicate("f_block_element_present", "eq", True), {}) is False
    assert evaluate(Predicate("active_f_weight_max", "lt", 10.0), {}) is False


# --- regression: two defects that were fixed --------------------------------


def test_missing_field_ordering_ops_return_false():
    """Regression: an ordering comparison on a missing field must return False rather
    than raising TypeError."""
    for op in ("gt", "ge", "lt", "le"):
        assert evaluate(Predicate("absent", op, 1), {}) is False


def test_load_rules_reports_all_problems(tmp_path):
    """Regression: with more than one bad file, all problems must be reported at once
    rather than raising on the first one."""
    bad_rule = (
        "- id: X1\n  kind: recipe\n  title: t\n"
        "  condition: {all: []}\n  action: a\n  evidence: []\n"
    )
    (tmp_path / "a.yaml").write_text(bad_rule, encoding="utf-8")
    (tmp_path / "b.yaml").write_text(bad_rule, encoding="utf-8")
    with pytest.raises(RuleError) as excinfo:
        load_rules(tmp_path)
    message = str(excinfo.value)
    assert "a.yaml" in message and "b.yaml" in message
    assert "2 problem(s)" in message
