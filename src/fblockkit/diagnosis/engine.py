"""Diagnosis engine: rule evaluation and conclusion generation
(architecture design v0.1 §3 L2'').

The input is a parser-layer result (or an equivalent fact table) and the output
is a sequence of Findings graded by severity. This layer only evaluates
deterministically and never modifies any of the user's files; every conclusion
carries its provenance (from the rule's evidence).

Discipline:

- only diagnosis rules (kind=diagnosis) are accepted; a recipe rule mixed in is
  an error rather than being skipped silently;
- the rule order is the conclusion order (auditable: conclusions correspond to
  the rule files line by line).
"""

from __future__ import annotations

from typing import Iterable, Mapping

from ..knowledge.loader import evaluate_rules, load_rules
from ..knowledge.models import RULE_DIAGNOSIS, Finding, ParseResult, Rule
from ..parsers import facts_from


class DiagnosisError(ValueError):
    """Invalid input for the diagnosis layer."""


def finding_from_rule(rule: Rule) -> Finding:
    """Translate a matched diagnosis rule into one Finding (conclusion = title,
    action = action)."""
    return Finding(
        severity=rule.severity,
        message=rule.title,
        evidence=rule.evidence,
        suggested_fix=rule.action,
        refusals=rule.refusal,
        rule_id=rule.id,
    )


def diagnose_facts(
    facts: Mapping[str, object], rules: Iterable[Rule] | None = None
) -> tuple[Finding, ...]:
    """Evaluate against one fact table (handy for checking hypotheses that were never
    written to disk; also the test entry point)."""
    if rules is None:
        rules = (rule for rule in load_rules() if rule.kind == RULE_DIAGNOSIS)
    materialized = tuple(rules)
    for rule in materialized:
        if rule.kind != RULE_DIAGNOSIS:
            raise DiagnosisError(
                f"the diagnosis engine accepts diagnosis rules only, got {rule.id!r} "
                f"(kind={rule.kind!r}). Next step: filter by kind before passing them in."
            )
    return tuple(finding_from_rule(rule) for rule in evaluate_rules(materialized, facts))


def diagnose(result: ParseResult, rules: Iterable[Rule] | None = None) -> tuple[Finding, ...]:
    """Check up on the parse result of one program output."""
    return diagnose_facts(facts_from(result), rules)
