"""Recipe-layer main entry: system profile -> rule matching + G1 basis advice ->
Recommendation.

Flow (architecture design v0.1 §3 L2):

1. ``facts_from_profile``: derive the recipe rules' fact fields from a SystemProfile;
2. rule matching (rules of kind=recipe; a diagnosis rule mixed in is an error);
3. G1 basis-set / ECP advice (``basis_ecp.recommend_basis_ecp``);
4. assemble a Recommendation: method chain / basis sets / key settings / evidence /
   refusals / warnings.

Refusal semantics: matched refusal-type rules and non-applicable tiers are written
into ``Recommendation.refusals`` rather than being silently downgraded.
"""

from __future__ import annotations

from typing import Iterable

from ..knowledge.loader import evaluate_rules, load_rules
from ..knowledge.models import RULE_RECIPE, Recommendation, Rule, SystemProfile
from .basis_ecp import (
    BasisAdvice,
    BasisEntry,
    is_f_element,
    recommend_basis_ecp,
)

__all__ = [
    "facts_from_profile",
    "recommend",
    "recommend_with_advice",
]


def facts_from_profile(profile: SystemProfile) -> dict[str, object]:
    """Derive the recipe rules' fact fields (vocabulary in knowledge/rules/README.md).

    Conventions:
    - ``f_block``: the element list holds an Ln/An, or the user reported an f electron count;
    - ``soc_task``: the targets include magnetic or spectra (in the f block such targets
      go through the SOC chain);
    - ``active_space_orbitals``: supplied by the active-space template (G2, a later step) --
      it is not given here, so the B-group size rules do not fire at the profile stage
      (a missing field evaluates to None).
    """
    targets = tuple(str(t) for t in profile.targets)
    return {
        "f_block": profile.f_count > 0 or any(is_f_element(s) for s in profile.elements),
        "targets": list(targets),
        "geometry_task": "geometry" in targets,
        "soc_task": ("magnetic" in targets) or ("spectra" in targets),
        "open_shells": profile.open_shells,
    }


def recommend_with_advice(
    profile: SystemProfile,
    rules: Iterable[Rule] | None = None,
    entries: Iterable[BasisEntry] | None = None,
) -> tuple[Recommendation, BasisAdvice]:
    """As ``recommend``, but also returns the structured G1 advice (the renderer needs
    the selected entry)."""
    if rules is None:
        rules = tuple(rule for rule in load_rules() if rule.kind == RULE_RECIPE)
    materialized = tuple(rules)
    for rule in materialized:
        if rule.kind != RULE_RECIPE:
            raise ValueError(
                f"the recipe layer accepts recipe rules only, got {rule.id!r} "
                f"(kind={rule.kind!r}). Next step: filter by kind before passing them in."
            )
    facts = facts_from_profile(profile)
    matched = evaluate_rules(materialized, facts)
    advice = recommend_basis_ecp(profile, tuple(entries) if entries is not None else None)

    method_chain = tuple(dict.fromkeys(rule.method for rule in matched if rule.method))
    key_settings = tuple(rule.action for rule in matched)
    refusals = tuple(text for rule in matched for text in rule.refusal)
    refusals += advice.refusal_texts()
    warnings = tuple(f"Note: {entry.note.strip()}" for entry in advice.cautions)
    warnings += advice.notes
    evidence = tuple(item for rule in matched for item in rule.evidence) + advice.evidence
    recommendation = Recommendation(
        method_chain=method_chain,
        basis_ecp=advice.lines(),
        active_space="",
        key_settings=key_settings,
        evidence=evidence,
        refusals=refusals,
        warnings=warnings,
    )
    return recommendation, advice


def recommend(
    profile: SystemProfile,
    rules: Iterable[Rule] | None = None,
    entries: Iterable[BasisEntry] | None = None,
) -> Recommendation:
    """System profile -> recommendation (with basis advice, key settings and provenance)."""
    return recommend_with_advice(profile, rules, entries)[0]
