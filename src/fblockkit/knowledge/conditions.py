"""Structured conditions for rule matching.

A condition is a small declarative tree evaluated against a flat mapping of
facts (field -> value). Deliberately simple and auditable: fixed operator
vocabulary, AND/OR combination, no eval, no code strings.

Semantics:
- ``AllOf([])`` is True (empty conjunction).
- ``AnyOf([])`` is False (empty disjunction).
- A missing field evaluates to None: equality tests compare against None,
  while ordering comparisons (gt/ge/lt/le) return False instead of raising.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence, Union


@dataclass(frozen=True)
class Predicate:
    """A single field/op/value test."""

    field: str
    op: str
    value: Any


@dataclass(frozen=True)
class AllOf:
    children: tuple["Condition", ...]


@dataclass(frozen=True)
class AnyOf:
    children: tuple["Condition", ...]


Condition = Union[Predicate, AllOf, AnyOf]

OPERATORS = ("eq", "ne", "in", "not_in", "contains", "gt", "ge", "lt", "le")
_ORDER_OPS = ("gt", "ge", "lt", "le")


class ConditionError(ValueError):
    """Raised for malformed condition trees."""


def _compare(op: str, fact: Any, value: Any) -> bool:
    if op in _ORDER_OPS and fact is None:
        # Documented semantics: a missing fact never raises -- ordering
        # comparisons on it are simply False (equality ops still
        # compare against None).
        return False
    if op == "eq":
        return fact == value
    if op == "ne":
        return fact != value
    if op == "in":
        return fact in value
    if op == "not_in":
        return fact not in value
    if op == "contains":
        return value in fact
    if op == "gt":
        return fact > value
    if op == "ge":
        return fact >= value
    if op == "lt":
        return fact < value
    if op == "le":
        return fact <= value
    raise ConditionError(f"unknown operator: {op!r} (allowed: {', '.join(OPERATORS)})")


def evaluate(condition: Condition, facts: Mapping[str, Any]) -> bool:
    """Evaluate a condition tree against a facts mapping.

    A missing field evaluates to None rather than raising -- rule authors get
    a clean ``False`` instead of a crash when a fact is unavailable (ordering
    comparisons return False; ``eq None`` still matches, by design).
    """
    if isinstance(condition, AllOf):
        return all(evaluate(child, facts) for child in condition.children)
    if isinstance(condition, AnyOf):
        return any(evaluate(child, facts) for child in condition.children)
    if isinstance(condition, Predicate):
        fact = facts.get(condition.field)
        return _compare(condition.op, fact, condition.value)
    raise ConditionError(f"not a condition node: {condition!r}")


def from_mapping(entry: Any) -> Condition:
    """Build a Condition from its YAML representation.

    Accepted forms::

        {field: f_block, op: eq, value: true}
        {all: [<condition>, ...]}
        {any: [<condition>, ...]}
    """
    if not isinstance(entry, Mapping):
        raise ConditionError(f"condition must be a mapping, got {type(entry).__name__}")
    keys = set(entry)
    if keys == {"all"}:
        children = entry["all"]
        if not isinstance(children, Sequence):
            raise ConditionError("'all' must be a sequence of conditions")
        return AllOf(tuple(from_mapping(child) for child in children))
    if keys == {"any"}:
        children = entry["any"]
        if not isinstance(children, Sequence):
            raise ConditionError("'any' must be a sequence of conditions")
        return AnyOf(tuple(from_mapping(child) for child in children))
    if keys == {"field", "op", "value"}:
        op = entry["op"]
        if op not in OPERATORS:
            raise ConditionError(
                f"unknown operator: {op!r} (allowed: {', '.join(OPERATORS)})"
            )
        if not isinstance(entry["field"], str) or not entry["field"]:
            raise ConditionError("'field' must be a non-empty string")
        return Predicate(field=entry["field"], op=op, value=entry["value"])
    raise ConditionError(
        "condition must be {field, op, value} or {all: [...]} or {any: [...]}; "
        f"got keys: {sorted(keys)}"
    )
