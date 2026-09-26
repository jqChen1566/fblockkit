"""Diagnosis report: Report assembly and Markdown rendering.

Report discipline: every conclusion carries its provenance; provenance is
summarised at the end of the document under three headings (manual /
literature / measured by this group). The output is Markdown text (the user
interface also saves the same content as JSON, which is a later step).
"""

from __future__ import annotations

from typing import Iterable, Sequence

from ..knowledge import sources
from ..knowledge.models import Evidence, Finding, Report, ReportSection

# The references block always renders last (after the provenance summary).
REFERENCES_TITLE = "References"

SEVERITY_LABELS = {
    "refuse": "Refuse",
    "error": "Error",
    "warn": "Warning",
    "info": "Info",
}
# Display order: most severe first, info last
_SEVERITY_ORDER = ("refuse", "error", "warn", "info")

EVIDENCE_LABELS = {
    "manual": "Manual",
    "literature": "Literature",
    "measured": "Measured",
}


def _severity_rank(finding: Finding) -> int:
    return _SEVERITY_ORDER.index(finding.severity) if finding.severity in _SEVERITY_ORDER else len(_SEVERITY_ORDER)


def references_section(provenance: Sequence[Evidence]) -> ReportSection | None:
    """Complete citations + paste-ready BibTeX for every literature item cited.

    Returns None when nothing cited a bibliography entry.
    """
    keys: list[str] = []
    for item in provenance:
        if item.bibkey and item.bibkey not in keys:
            keys.append(item.bibkey)
    if not keys:
        return None
    lines = ["Complete citations:"]
    for key in keys:
        lines.append(f"- [{key}] {sources.format_citation(sources.get(key))}")
    lines += ["", "BibTeX (paste-ready):"]
    for key in keys:
        lines += ["", "```bibtex", sources.to_bibtex(sources.get(key)), "```"]
    return ReportSection(title=REFERENCES_TITLE, body="\n".join(lines))


def build_report(
    findings: Iterable[Finding],
    title: str = "fBlockKit report",
    subject: str = "",
    sections: Iterable[ReportSection] = (),
    extra_evidence: Iterable[Evidence] = (),
) -> Report:
    """Assemble a Report from a sequence of Findings (sorted by severity, with
    deduplicated provenance).

    ``sections`` merges in the analysis layer's report sections (A1/A2/S1 and so
    on), giving a complete "analysis + diagnosis" report; ``extra_evidence``
    collects the analysis modules' provenance (A2's threshold reference, say) so
    that the References block covers every citation.
    """
    ordered = sorted(findings, key=_severity_rank)
    provenance: list[Evidence] = []
    seen: set[tuple[str, str, str]] = set()
    for finding in ordered:
        for item in finding.evidence:
            key = (item.kind, item.text, item.ref)
            if key not in seen:
                seen.add(key)
                provenance.append(item)
    for item in extra_evidence:
        key = (item.kind, item.text, item.ref)
        if key not in seen:
            seen.add(key)
            provenance.append(item)
    counts = {severity: 0 for severity in _SEVERITY_ORDER}
    for finding in ordered:
        counts[finding.severity] = counts.get(finding.severity, 0) + 1
    statistics = "; ".join(
        f"{SEVERITY_LABELS[severity]} {counts[severity]}"
        for severity in _SEVERITY_ORDER
        if counts.get(severity)
    )
    body_lines = []
    if subject:
        body_lines.append(f"Subject: {subject}")
    body_lines.append(f"Counts: {statistics if statistics else 'no findings'}")
    all_sections = [ReportSection(title="Summary", body="\n".join(body_lines)), *sections]
    refs = references_section(provenance)
    if refs is not None:
        all_sections.append(refs)
    return Report(
        title=title,
        sections=tuple(all_sections),
        findings=tuple(ordered),
        provenance=tuple(provenance),
    )


def to_markdown(report: Report) -> str:
    """Render as Markdown."""
    lines: list[str] = [f"# {report.title}", ""]
    late = [section for section in report.sections if section.title == REFERENCES_TITLE]
    for section in report.sections:
        if section.title == REFERENCES_TITLE:
            continue
        lines += [f"## {section.title}", "", section.body, ""]
    lines += ["## Findings", ""]
    if not report.findings:
        lines += ["No diagnostic rule matched.", ""]
    for index, finding in enumerate(report.findings, start=1):
        label = SEVERITY_LABELS.get(finding.severity, finding.severity)
        lines.append(f"### {index}. [{label}] {finding.message}")
        lines.append("")
        if finding.suggested_fix:
            lines.append(f"- Suggested action: {finding.suggested_fix}")
        for refusal in finding.refusals:
            lines.append(f"- Must not: {refusal}")
        if finding.rule_id:
            lines.append(f"- Rule: {finding.rule_id}")
        if finding.evidence:
            lines.append("- Evidence:")
            for item in finding.evidence:
                kind = EVIDENCE_LABELS.get(item.kind, item.kind)
                lines.append(f"  - [{kind}] {item.text}")
                ref = f"    - Source: {item.ref}"
                if item.url:
                    ref += f" ({item.url})"
                lines.append(ref)
        lines.append("")
    lines += ["## Provenance", ""]
    if not report.provenance:
        lines += ["(none)", ""]
    for item in report.provenance:
        kind = EVIDENCE_LABELS.get(item.kind, item.kind)
        lines.append(f"- [{kind}] {item.ref}")
    lines.append("")
    for section in late:  # references render last, in full
        lines += [f"## {section.title}", "", section.body, ""]
    return "\n".join(lines)


def render_markdown(
    findings: Sequence[Finding], title: str = "fBlockKit report", subject: str = ""
) -> str:
    """One-step form: a sequence of Findings -> Markdown text."""
    return to_markdown(build_report(findings, title=title, subject=subject))
