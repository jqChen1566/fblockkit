"""Core data models (architecture design v0.1, section 2).

All models are immutable dataclasses. Provenance is mandatory by design:
every user-facing statement must carry where it comes from -- a manual quote,
a literature reference, or a measurement made and recorded by this group.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping

from .conditions import Condition

# --- closed vocabularies -------------------------------------------------

EVIDENCE_MANUAL = "manual"
EVIDENCE_LITERATURE = "literature"
EVIDENCE_MEASURED = "measured"
EVIDENCE_KINDS = (EVIDENCE_MANUAL, EVIDENCE_LITERATURE, EVIDENCE_MEASURED)

RULE_RECIPE = "recipe"
RULE_DIAGNOSIS = "diagnosis"
RULE_KINDS = (RULE_RECIPE, RULE_DIAGNOSIS)

CONFIDENCE_CONFIRMED = "confirmed"      # backed directly by manual/literature
CONFIDENCE_PROVISIONAL = "provisional"  # our inference, pending verification

SEVERITIES = ("info", "warn", "error", "refuse")

RELATION_ABSORB = "absorb"      # tier A: implemented inside fBlockKit
RELATION_INTERFACE = "interface"  # tier B: generate its input / read its output
RELATION_INDEX = "index"        # tier C: registered for orientation only


# --- models --------------------------------------------------------------

@dataclass(frozen=True)
class Evidence:
    """One provenance record attached to a rule, finding or recommendation."""

    kind: str  # manual / literature / measured
    text: str
    ref: str
    url: str = ""
    bibkey: str = ""  # required for kind == "literature"; must exist in sources.bib


@dataclass(frozen=True)
class SystemProfile:
    """What the user has told us about the system under study."""

    elements: tuple[str, ...] = ()
    charge: int = 0
    spin: int = 0  # 2S
    open_shells: int = 0
    f_count: int = 0
    metal_valence: int = 0  # f-block metal valence (III=3 / IV=4; 0 = not specified)
    targets: tuple[str, ...] = ()  # energy / geometry / excited / magnetic / spectra
    budget: str = "standard"       # screening / standard / high


@dataclass(frozen=True)
class Rule:
    """An auditable recipe or diagnosis rule (architecture design v0.1, 6.1)."""

    id: str
    kind: str  # recipe / diagnosis
    title: str
    condition: Condition
    action: str
    evidence: tuple[Evidence, ...]
    refusal: tuple[str, ...] = ()
    confidence: str = CONFIDENCE_CONFIRMED
    severity: str = ""  # required for diagnosis rules (info/warn/error/refuse); forbidden for recipe rules
    method: str = ""  # recipe rules only: the method name this rule contributes (collected into Recommendation.method_chain)


@dataclass(frozen=True)
class Recommendation:
    """Output of the recipe layer."""

    method_chain: tuple[str, ...] = ()
    basis_ecp: tuple[str, ...] = ()
    active_space: str = ""
    key_settings: tuple[str, ...] = ()
    evidence: tuple[Evidence, ...] = ()
    refusals: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()


@dataclass(frozen=True)
class ParseResult:
    """Structured view of an output file produced by the parsers layer."""

    program: str
    path: str
    sections: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class Finding:
    """One diagnostic conclusion."""

    severity: str  # info / warn / error / refuse
    message: str
    evidence: tuple[Evidence, ...] = ()
    suggested_fix: str = ""
    fix_template: str = ""
    refusals: tuple[str, ...] = ()  # the rule's "must not" statements (refusal semantics)
    rule_id: str = ""  # id of the rule that produced this finding (auditable)


@dataclass(frozen=True)
class ToolEntry:
    """One entry of the external-tool index (A/B/C tiers)."""

    name: str
    purpose: str
    license: str
    source: str
    relation: str  # absorb / interface / index
    status: str = "active"
    notes: str = ""


@dataclass(frozen=True)
class ReportSection:
    title: str
    body: str


@dataclass(frozen=True)
class Report:
    title: str
    sections: tuple[ReportSection, ...] = ()
    findings: tuple[Finding, ...] = ()
    provenance: tuple[Evidence, ...] = ()
