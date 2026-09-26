"""Loading and validating knowledge data (rules for now).

Validation is strict on purpose: a rule without provenance is rejected at
load time (architecture design v0.1, section 2). Error messages say what is
wrong, what the value was, and what to do -- the same discipline the CLI
follows for user-facing errors.

Loading collects *all* problems across files and entries before raising, so
rule curators can fix a batch in one pass instead of one file per run.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import yaml

from .conditions import ConditionError, evaluate, from_mapping
from .models import (
    CONFIDENCE_CONFIRMED,
    CONFIDENCE_PROVISIONAL,
    EVIDENCE_KINDS,
    EVIDENCE_LITERATURE,
    RULE_DIAGNOSIS,
    RULE_KINDS,
    SEVERITIES,
    Evidence,
    Rule,
)
from .sources import BibDataError, index_by_key

DEFAULT_RULES_DIR = Path(__file__).resolve().parent / "rules"

_BIB_CACHE: dict[str, object] | None = None


def bib_keys() -> dict[str, object]:
    """Citation keys of sources.bib (cached; literature evidence must point here)."""
    global _BIB_CACHE
    if _BIB_CACHE is None:
        try:
            _BIB_CACHE = index_by_key()
        except BibDataError as exc:
            raise RuleError(f"sources.bib cannot be loaded: {exc}") from exc
    return _BIB_CACHE


class RuleError(ValueError):
    """A rule file failed validation."""


def _load_yaml(path: Path) -> Any:
    try:
        with path.open("r", encoding="utf-8") as handle:
            return yaml.safe_load(handle)
    except yaml.YAMLError as exc:
        raise RuleError(f"{path.name}: YAML parse failed: {exc}") from exc


def _build_evidence(entries: Any, where: str) -> tuple[Evidence, ...]:
    if not isinstance(entries, Sequence) or isinstance(entries, (str, bytes)) or not entries:
        raise RuleError(
            f"{where}: evidence is missing or empty -- every rule must carry its "
            f"provenance (manual/literature/measured)"
        )
    built = []
    for index, entry in enumerate(entries):
        if not isinstance(entry, Mapping):
            raise RuleError(f"{where}: evidence[{index}] must be a mapping")
        kind = entry.get("kind")
        if kind not in EVIDENCE_KINDS:
            raise RuleError(
                f"{where}: evidence[{index}].kind={kind!r} is not one of {EVIDENCE_KINDS}"
            )
        text = entry.get("text")
        ref = entry.get("ref")
        if not isinstance(text, str) or not text.strip():
            raise RuleError(f"{where}: evidence[{index}].text is missing")
        if not isinstance(ref, str) or not ref.strip():
            raise RuleError(
                f"{where}: evidence[{index}].ref is missing (section number / DOI / "
                f"measurement id)"
            )
        bibkey = str(entry.get("bibkey", "")).strip()
        if kind == EVIDENCE_LITERATURE:
            if not bibkey:
                raise RuleError(
                    f"{where}: evidence[{index}] is literature but has no bibkey "
                    "(a literature reference must point into sources.bib so the output "
                    "can print the full citation)"
                )
            if bibkey not in bib_keys():
                raise RuleError(
                    f"{where}: evidence[{index}].bibkey={bibkey!r} is not in sources.bib"
                )
        built.append(
            Evidence(
                kind=kind,
                text=text,
                ref=ref,
                url=str(entry.get("url", "")),
                bibkey=bibkey,
            )
        )
    return tuple(built)


def _build_rule(entry: Any, where: str) -> Rule:
    if not isinstance(entry, Mapping):
        raise RuleError(f"{where}: a rule must be a mapping")
    rule_id = entry.get("id")
    if not isinstance(rule_id, str) or not rule_id.strip():
        raise RuleError(f"{where}: id is missing")
    where = f"{where}[{rule_id}]"
    kind = entry.get("kind")
    if kind not in RULE_KINDS:
        raise RuleError(f"{where}: kind={kind!r} is not one of {RULE_KINDS}")
    title = entry.get("title")
    if not isinstance(title, str) or not title.strip():
        raise RuleError(f"{where}: title is missing")
    action = entry.get("action")
    if not isinstance(action, str) or not action.strip():
        raise RuleError(f"{where}: action is missing")
    try:
        condition = from_mapping(entry.get("condition"))
    except ConditionError as exc:
        raise RuleError(f"{where}: invalid condition: {exc}") from exc
    confidence = entry.get("confidence", CONFIDENCE_CONFIRMED)
    if confidence not in (CONFIDENCE_CONFIRMED, CONFIDENCE_PROVISIONAL):
        raise RuleError(
            f"{where}: confidence={confidence!r} must be {CONFIDENCE_CONFIRMED} or "
            f"{CONFIDENCE_PROVISIONAL}"
        )
    refusal = entry.get("refusal", [])
    if not isinstance(refusal, Sequence) or isinstance(refusal, (str, bytes)):
        raise RuleError(f"{where}: refusal must be a list")
    severity = entry.get("severity", "")
    if kind == RULE_DIAGNOSIS:
        if severity not in SEVERITIES:
            raise RuleError(
                f"{where}: a diagnosis rule must give severity in {SEVERITIES}, "
                f"currently {severity!r}"
            )
    elif severity:
        raise RuleError(f"{where}: a recipe rule must not carry severity ({severity!r})")
    method = str(entry.get("method", ""))
    if kind == RULE_DIAGNOSIS and method:
        raise RuleError(f"{where}: a diagnosis rule must not carry method ({method!r})")
    return Rule(
        id=rule_id,
        kind=kind,
        title=title,
        condition=condition,
        action=action,
        evidence=_build_evidence(entry.get("evidence"), where),
        refusal=tuple(str(item) for item in refusal),
        confidence=confidence,
        severity=severity,
        method=method,
    )


def load_rules(directory: Path | str = DEFAULT_RULES_DIR) -> tuple[Rule, ...]:
    """Load and validate every ``*.yaml`` rule file in ``directory``.

    All problems are collected and reported together (see module docstring).
    """
    path = Path(directory)
    if not path.is_dir():
        raise RuleError(f"rule directory does not exist: {path}")
    rules: list[Rule] = []
    seen: dict[str, str] = {}
    errors: list[str] = []
    for yaml_file in sorted(path.glob("*.yaml")):
        try:
            data = _load_yaml(yaml_file)
        except RuleError as exc:
            errors.append(str(exc))
            continue
        if data is None:
            continue
        if not isinstance(data, Sequence) or isinstance(data, (str, bytes)):
            errors.append(f"{yaml_file.name}: the top level must be a list of rules")
            continue
        for entry in data:
            try:
                rule = _build_rule(entry, yaml_file.name)
            except RuleError as exc:
                errors.append(str(exc))
                continue
            if rule.id in seen:
                errors.append(
                    f"duplicate rule id: {rule.id!r} (already seen in {seen[rule.id]}, "
                    f"now again in {yaml_file.name})"
                )
                continue
            seen[rule.id] = yaml_file.name
            rules.append(rule)
    if errors:
        joined = "\n".join(f"- {message}" for message in errors)
        raise RuleError(f"rule loading failed ({len(errors)} problem(s)):\n{joined}")
    return tuple(rules)


def evaluate_rules(
    rules: Iterable[Rule], facts: Mapping[str, object]
) -> tuple[Rule, ...]:
    """Return the rules whose condition matches ``facts`` (order preserved)."""
    return tuple(rule for rule in rules if evaluate(rule.condition, facts))
