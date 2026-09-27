"""The common handler group."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from ...diagnosis import SEVERITY_LABELS
from ...parsers import parse_auto
from ..session import Session
def _plain(value: Any) -> Any:
    if isinstance(value, tuple):
        return [_plain(item) for item in value]
    if isinstance(value, dict):
        return {key: _plain(item) for key, item in value.items()}
    return value


def _write_report_files(path: Path, markdown: str, payload: dict[str, Any]) -> tuple[Path, Path]:
    md_path = path.with_name(path.name + ".fbk.md")
    json_path = path.with_name(path.name + ".fbk.json")
    md_path.write_text(markdown, encoding="utf-8")
    json_path.write_text(
        json.dumps(_plain(payload), ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return md_path, json_path


def _report_payload(report, extra: dict[str, Any]) -> dict[str, Any]:
    return {
        "title": report.title,
        "sections": [{"title": section.title, "body": section.body} for section in report.sections],
        "findings": [
            {
                "severity": finding.severity,
                "message": finding.message,
                "rule_id": finding.rule_id,
                "suggested_fix": finding.suggested_fix,
                "refusals": list(finding.refusals),
            }
            for finding in report.findings
        ],
        "provenance": [
            {"kind": item.kind, "text": item.text, "ref": item.ref, "url": item.url}
            for item in report.provenance
        ],
        **extra,
    }


def _summarize_findings(session: Session, findings) -> None:
    if not findings:
        session.say("Diagnosis: no diagnostic rule matched.")
        return
    counts: dict[str, int] = {}
    for finding in findings:
        counts[finding.severity] = counts.get(finding.severity, 0) + 1
    text = "; ".join(
        f"{SEVERITY_LABELS.get(severity, severity)} {count}" for severity, count in counts.items()
    )
    session.say(f"Diagnosis: {text} (details in the report file).")


def casscf_reference(session: Session, output_text: str) -> tuple[bool, float | None]:
    """The printed CASSCF energy of an optional output file (the CI-root anchor).

    Returns ``(ok, energy)``; with no output requested it returns ``(True, None)``
    and with a missing or unconverged CASSCF section it says the shared guidance
    and returns ``(False, None)`` so the caller stops -- the flow the
    FCIDUMP-driven menus (24/25/26) share.  A malformed file raises the parser's
    ``ParserError`` into the caller's own error path.
    """
    if not output_text:
        return True, None
    result = parse_auto(Path(output_text))
    casscf = result.sections.get("casscf", {})
    if not casscf.get("present") or not casscf.get("converged"):
        session.say(
            "The output carries no converged CASSCF section to match the CI "
            "root against. Next step: give the converged run's output, or leave "
            "the output prompt empty to use the lowest root of the FCIDUMP's "
            "Ms sector."
        )
        return False, None
    return True, casscf.get("energy")
