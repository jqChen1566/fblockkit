"""Tool-index layer (T1/T2): registration, search and onboarding guidance for the
A/B/C tiers of external tools.

Data: ``knowledge/toolindex/tools.yaml`` (every entry carries relation and
provenance; loading validates strictly). The meaning of relation (in step with the
architecture design): absorb = already built in / directly absorbable; interface =
we interface with it (generate its input, read its output); index = registered only
(a compute-kernel type, possibly with an algorithm borrowed from it).

Discipline: the license field follows the official source; where this group's
documents do not state it, it is recorded as "to be checked" rather than guessed.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import yaml

from ..knowledge.loader import bib_keys
from ..knowledge.models import (
    EVIDENCE_KINDS,
    EVIDENCE_LITERATURE,
    RELATION_ABSORB,
    RELATION_INDEX,
    RELATION_INTERFACE,
    Evidence,
)

DEFAULT_TOOLS_FILE = Path(__file__).resolve().parent.parent / "knowledge" / "toolindex" / "tools.yaml"

RELATIONS = (RELATION_ABSORB, RELATION_INTERFACE, RELATION_INDEX)
RELATION_LABELS = {
    RELATION_ABSORB: "A absorb directly",
    RELATION_INTERFACE: "B interface",
    RELATION_INDEX: "C registered (index)",
}

# relation is the "strategy" (how we mean to use it) and status the "present state"
# (whether it is used yet) -- the two must be recorded separately, otherwise
# "registered" gets read as "used".
STATUS_ACTIVE = "active"
STATUS_PLANNED = "planned"
STATUSES = (STATUS_ACTIVE, STATUS_PLANNED)
STATUS_LABELS = {STATUS_ACTIVE: "in use", STATUS_PLANNED: "registered (not in use)"}


class ToolIndexError(ValueError):
    """The tool-index data is invalid."""


@dataclass(frozen=True)
class ToolRecord:
    id: str
    name: str
    purpose: str
    relation: str
    status: str
    license: str
    source: str
    note: str
    evidence: tuple[Evidence, ...] = ()

    def line(self) -> str:
        label = RELATION_LABELS.get(self.relation, self.relation)
        state = STATUS_LABELS.get(self.status, self.status)
        return f"{self.name} ({label}; {state}; license: {self.license})"


def _build(entry: object, where: str) -> ToolRecord:
    if not isinstance(entry, dict):
        raise ToolIndexError(f"{where}: an entry must be a mapping")
    tool_id = entry.get("id")
    if not isinstance(tool_id, str) or not tool_id.strip():
        raise ToolIndexError(f"{where}: id is missing")
    where = f"{where}[{tool_id}]"
    for field in ("name", "purpose", "license", "source"):
        value = entry.get(field)
        if not isinstance(value, str) or not value.strip():
            raise ToolIndexError(f"{where}: {field} is missing")
    relation = entry.get("relation")
    if relation not in RELATIONS:
        raise ToolIndexError(f"{where}: relation={relation!r} is not one of {RELATIONS}")
    status = entry.get("status")
    if status not in STATUSES:
        raise ToolIndexError(
            f"{where}: status={status!r} is not one of {STATUSES} "
            f"(relation records the strategy and status the present state; both are required)"
        )
    raw_evidence = entry.get("evidence")
    if not isinstance(raw_evidence, list) or not raw_evidence:
        raise ToolIndexError(
            f"{where}: evidence is missing or empty -- every entry must carry its provenance"
        )
    evidence = []
    for index, item in enumerate(raw_evidence):
        if not isinstance(item, dict) or item.get("kind") not in EVIDENCE_KINDS:
            raise ToolIndexError(f"{where}: evidence[{index}] is invalid")
        if not str(item.get("text", "")).strip() or not str(item.get("ref", "")).strip():
            raise ToolIndexError(f"{where}: evidence[{index}] is missing text or ref")
        bibkey = str(item.get("bibkey", "")).strip()
        if item["kind"] == EVIDENCE_LITERATURE:
            if not bibkey:
                raise ToolIndexError(
                    f"{where}: evidence[{index}] is literature but has no bibkey (it must "
                    "point into sources.bib)"
                )
            if bibkey not in bib_keys():
                raise ToolIndexError(
                    f"{where}: evidence[{index}].bibkey={bibkey!r} is not in sources.bib"
                )
        evidence.append(
            Evidence(
                kind=item["kind"],
                text=str(item["text"]).strip(),
                ref=str(item["ref"]).strip(),
                url=str(item.get("url", "")),
                bibkey=bibkey,
            )
        )
    return ToolRecord(
        id=tool_id,
        name=entry["name"],
        purpose=entry["purpose"],
        relation=relation,
        status=status,
        license=entry["license"],
        source=entry["source"],
        note=str(entry.get("note", "")).strip(),
        evidence=tuple(evidence),
    )


def load_tools(path: Path | str = DEFAULT_TOOLS_FILE) -> tuple[ToolRecord, ...]:
    """Load and validate the tool index; all problems are reported at once."""
    p = Path(path)
    if not p.is_file():
        raise ToolIndexError(f"tool index file does not exist: {p}")
    try:
        data = yaml.safe_load(p.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise ToolIndexError(f"{p.name}: YAML parse failed: {exc}") from exc
    if not isinstance(data, list):
        raise ToolIndexError(f"{p.name}: the top level must be a list of entries")
    records, seen, errors = [], set(), []
    for entry in data:
        try:
            record = _build(entry, p.name)
        except ToolIndexError as exc:
            errors.append(str(exc))
            continue
        if record.id in seen:
            errors.append(f"duplicate entry id: {record.id!r} ({p.name})")
            continue
        seen.add(record.id)
        records.append(record)
    if errors:
        raise ToolIndexError(
            "tool index loading failed (%d problem(s)):\n%s"
            % (len(errors), "\n".join(f"- {e}" for e in errors))
        )
    return tuple(records)


def search(query: str, records: Sequence[ToolRecord] | None = None) -> tuple[ToolRecord, ...]:
    """Keyword search: every space-separated term must hit (name/purpose/note/source),
    case-insensitive."""
    if records is None:
        records = load_tools()
    terms = [term.lower() for term in re.split(r"\s+", query.strip()) if term]
    if not terms:
        return ()
    hits = []
    for record in records:
        haystack = " ".join(
            (record.name, record.purpose, record.note, record.source, record.relation)
        ).lower()
        if all(term in haystack for term in terms):
            hits.append(record)
    hits.sort(key=lambda record: (terms[0] not in record.name.lower(), record.name.lower()))
    return tuple(hits)


def guide(tool_id: str, records: Sequence[ToolRecord] | None = None) -> str:
    """T2: onboarding notes for one tool (install/license/how to obtain + provenance)."""
    if records is None:
        records = load_tools()
    for record in records:
        if record.id == tool_id:
            lines = [
                f"{record.name} -- {RELATION_LABELS.get(record.relation, record.relation)}",
                f"Purpose: {record.purpose}",
                f"Status: {STATUS_LABELS.get(record.status, record.status)}",
                f"License: {record.license}",
                f"Source: {record.source}",
            ]
            if record.note:
                lines.append(f"Note: {record.note}")
            lines.append("Provenance:")
            for item in record.evidence:
                lines.append(f"  - [{item.kind}] {item.text} ({item.ref})")
            return "\n".join(lines)
    raise ToolIndexError(
        f"the index has no tool with id={tool_id!r}. Next step: look keywords up with "
        f"search, or check the spelling of the id."
    )


__all__ = [
    "RELATIONS",
    "RELATION_LABELS",
    "ToolIndexError",
    "ToolRecord",
    "guide",
    "load_tools",
    "search",
]
