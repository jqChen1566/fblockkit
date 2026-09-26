"""Diagnosis layer: output check-up, cross-level consistency check and report
rendering (architecture design v0.1 §3 L2'').

Public interface:

- ``diagnose(result)`` / ``diagnose_facts(facts)``: evaluate against a parse result / a
  fact table and produce Findings;
- ``cross_level_check(records)``: cross-level ordering consistency of multiple solutions
  (literature fixture PuCl3);
- ``triage(result)`` / ``propose_fixes(input_text, findings)``: SCF rescue triage and
  corrected-input generation;
- ``build_report`` / ``to_markdown`` / ``render_markdown``: report assembly and rendering.

This layer never modifies the user's files; any "corrected input" is a separately
generated new file (``propose_fixes`` returns the text of new files).
"""

from .cross_level import CrossLevelError, cross_level_check, load_records
from .engine import DiagnosisError, diagnose, diagnose_facts, finding_from_rule
from .report import (
    EVIDENCE_LABELS,
    REFERENCES_TITLE,
    SEVERITY_LABELS,
    build_report,
    references_section,
    render_markdown,
    to_markdown,
)
from .scf_rescue import FixProposal, ScfRescueError, propose_fixes, triage

__all__ = [
    "CrossLevelError",
    "DiagnosisError",
    "EVIDENCE_LABELS",
    "FixProposal",
    "REFERENCES_TITLE",
    "SEVERITY_LABELS",
    "ScfRescueError",
    "build_report",
    "cross_level_check",
    "diagnose",
    "diagnose_facts",
    "finding_from_rule",
    "load_records",
    "propose_fixes",
    "references_section",
    "render_markdown",
    "to_markdown",
    "triage",
]
