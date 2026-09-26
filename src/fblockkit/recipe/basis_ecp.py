"""G1: Ln/An basis-set / ECP recommendation (the data side).

Data file ``knowledge/basis_ecp.yaml``: every entry carries the element range
covered, the applicable valences and targets, the ORCA keywords and its
provenance; loading validates strictly (no provenance means rejection, the same
discipline as the rule table).

Matching dimensions: element (Z range) x valence (all accepted when unspecified)
x target (all accepted when unspecified). The output falls into three tiers:
recommend (usable) / caution (boundary reminder) / refuse (unusable or not
automatically generatable).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import yaml

from ..knowledge.elements import (
    ELEMENT_Z,
    ElementError,
    element_z,
    is_f_element,
    parse_range as _parse_range,
)
from ..knowledge.loader import bib_keys
from ..knowledge.models import EVIDENCE_KINDS, EVIDENCE_LITERATURE, Evidence, SystemProfile

DEFAULT_BASIS_FILE = Path(__file__).resolve().parent.parent / "knowledge" / "basis_ecp.yaml"

KINDS = ("recommend", "caution", "refuse")


class BasisDataError(ValueError):
    """The basis-set data file is invalid."""


@dataclass(frozen=True)
class BasisEntry:
    """One basis-set / ECP suggestion (range covered + keywords + provenance)."""

    id: str
    kind: str
    elements: range
    valence: tuple[int, ...]
    targets: tuple[str, ...]
    note: str
    basis: str = ""
    ecp: str = ""
    hamiltonian: str = ""
    auxiliary: str = ""
    refusals: tuple[str, ...] = ()
    evidence: tuple[Evidence, ...] = ()

    def applies(self, profile: SystemProfile, z_values: Sequence[int]) -> bool:
        if not any(z in self.elements for z in z_values):
            return False
        if self.valence and profile.metal_valence and profile.metal_valence not in self.valence:
            return False
        if self.targets and profile.targets and not set(self.targets) & set(profile.targets):
            return False
        return True

    def line(self) -> str:
        """Single-line description (ASCII: it may end up in an ORCA input comment, so
        no full-width punctuation)."""
        head = self.basis or self.ecp
        details = []
        if self.hamiltonian:
            details.append(self.hamiltonian)
        if self.auxiliary:
            details.append(f"aux {self.auxiliary}")
        if self.ecp and self.basis:
            details.append(self.ecp)
        return f"{head} ({'; '.join(details)})" if details else head


@dataclass(frozen=True)
class BasisAdvice:
    """The summary output of one match."""

    applicable: bool  # does the system hold f-block elements (otherwise this table does not apply)
    recommendations: tuple[BasisEntry, ...] = ()
    cautions: tuple[BasisEntry, ...] = ()
    refusals: tuple[BasisEntry, ...] = ()
    evidence: tuple[Evidence, ...] = ()
    notes: tuple[str, ...] = ()

    def lines(self) -> tuple[str, ...]:
        """Text lines for Recommendation.basis_ecp (recommendations first, boundary
        reminders after)."""
        out = [entry.line() for entry in self.recommendations]
        out += [f"Note: {entry.note.strip()}" for entry in self.cautions]
        out += [f"Not applicable: {entry.note.strip()}" for entry in self.refusals]
        return tuple(out)

    def refusal_texts(self) -> tuple[str, ...]:
        return tuple(text for entry in self.refusals for text in entry.refusals)


def _build_entry(entry: Any, where: str) -> BasisEntry:
    if not isinstance(entry, Mapping):
        raise BasisDataError(f"{where}: an entry must be a mapping")
    entry_id = entry.get("id")
    if not isinstance(entry_id, str) or not entry_id.strip():
        raise BasisDataError(f"{where}: id is missing")
    where = f"{where}[{entry_id}]"
    kind = entry.get("kind")
    if kind not in KINDS:
        raise BasisDataError(f"{where}: kind={kind!r} is not one of {KINDS}")
    elements_text = entry.get("elements", "")
    try:
        elements = _parse_range(str(elements_text)) if elements_text else range(0)
    except (ElementError, BasisDataError) as exc:
        raise BasisDataError(f"{where}: invalid elements: {exc}") from exc
    if not elements and kind != "caution":
        raise BasisDataError(
            f"{where}: elements is missing (the element range covered, e.g. La-Lu)"
        )
    note = entry.get("note")
    if not isinstance(note, str) or not note.strip():
        raise BasisDataError(f"{where}: note is missing")
    basis = str(entry.get("basis", ""))
    ecp = str(entry.get("ecp", ""))
    if kind == "recommend" and not (basis or ecp):
        raise BasisDataError(f"{where}: a recommend entry must give basis or ecp")
    raw_evidence = entry.get("evidence")
    if not isinstance(raw_evidence, Sequence) or isinstance(raw_evidence, (str, bytes)) or not raw_evidence:
        raise BasisDataError(
            f"{where}: evidence is missing or empty -- every entry must carry its provenance"
        )
    evidence = []
    for index, item in enumerate(raw_evidence):
        if not isinstance(item, Mapping) or item.get("kind") not in EVIDENCE_KINDS:
            raise BasisDataError(
                f"{where}: evidence[{index}] is invalid (kind must be one of {EVIDENCE_KINDS})"
            )
        text, ref = item.get("text"), item.get("ref")
        if not (isinstance(text, str) and text.strip() and isinstance(ref, str) and ref.strip()):
            raise BasisDataError(f"{where}: evidence[{index}] is missing text or ref")
        bibkey = str(item.get("bibkey", "")).strip()
        if item["kind"] == EVIDENCE_LITERATURE:
            if not bibkey:
                raise BasisDataError(
                    f"{where}: evidence[{index}] is literature but has no bibkey (it must "
                    "point into sources.bib)"
                )
            if bibkey not in bib_keys():
                raise BasisDataError(
                    f"{where}: evidence[{index}].bibkey={bibkey!r} is not in sources.bib"
                )
        evidence.append(
            Evidence(
                kind=item["kind"],
                text=text,
                ref=ref,
                url=str(item.get("url", "")),
                bibkey=bibkey,
            )
        )
    valence = tuple(int(v) for v in entry.get("valence", ()) or ())
    targets = tuple(str(t) for t in entry.get("targets", ()) or ())
    refusals = tuple(str(r) for r in entry.get("refusals", ()) or ())
    return BasisEntry(
        id=entry_id,
        kind=kind,
        elements=elements,
        valence=valence,
        targets=targets,
        note=note,
        basis=basis,
        ecp=ecp,
        hamiltonian=str(entry.get("hamiltonian", "")),
        auxiliary=str(entry.get("auxiliary", "")),
        refusals=refusals,
        evidence=tuple(evidence),
    )


def load_basis_entries(path: Path | str = DEFAULT_BASIS_FILE) -> tuple[BasisEntry, ...]:
    """Load and validate the basis table; all problems are reported at once."""
    p = Path(path)
    if not p.is_file():
        raise BasisDataError(f"basis data file does not exist: {p}")
    try:
        data = yaml.safe_load(p.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise BasisDataError(f"{p.name}: YAML parse failed: {exc}") from exc
    if not isinstance(data, Sequence) or isinstance(data, (str, bytes)):
        raise BasisDataError(f"{p.name}: the top level must be a list of entries")
    entries, seen, errors = [], {}, []
    for item in data:
        try:
            entry = _build_entry(item, p.name)
        except BasisDataError as exc:
            errors.append(str(exc))
            continue
        if entry.id in seen:
            errors.append(f"duplicate entry id: {entry.id!r} ({p.name})")
            continue
        seen[entry.id] = True
        entries.append(entry)
    if errors:
        raise BasisDataError(
            "basis table loading failed (%d problem(s)):\n%s"
            % (len(errors), "\n".join(f"- {e}" for e in errors))
        )
    return tuple(entries)


def recommend_basis_ecp(
    profile: SystemProfile, entries: Sequence[BasisEntry] | None = None
) -> BasisAdvice:
    """Match by element range x valence x target; with no f-block element in the system
    it returns applicable=False explicitly."""
    if entries is None:
        entries = load_basis_entries()
    f_symbols = [symbol for symbol in profile.elements if is_f_element(symbol)]
    if not f_symbols:
        return BasisAdvice(
            applicable=False,
            notes=(
                "The system (by the given elements) contains no lanthanide/actinide -- "
                "this table is aimed at Ln/An calculations and does not apply.",
            ),
        )
    z_values = [element_z(symbol) for symbol in f_symbols]
    matched = [entry for entry in entries if entry.applies(profile, z_values)]
    notes: list[str] = []
    if profile.metal_valence:
        notes.append(
            f"Matched at metal valence {profile.metal_valence} (entries that do not list "
            f"a valence are treated as the general tier)."
        )
    else:
        notes.append(
            "Metal valence not specified: valence-specific tiers (5f-in-core being "
            "trivalent only, for instance) were left out; give the valence and query again."
        )
    return BasisAdvice(
        applicable=True,
        recommendations=tuple(e for e in matched if e.kind == "recommend"),
        cautions=tuple(e for e in matched if e.kind == "caution"),
        refusals=tuple(e for e in matched if e.kind == "refuse"),
        evidence=tuple(ev for entry in matched for ev in entry.evidence),
        notes=tuple(notes),
    )
