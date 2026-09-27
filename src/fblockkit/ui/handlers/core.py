"""The core handler group."""

from __future__ import annotations

from pathlib import Path

from .common import _write_report_files, _report_payload, _summarize_findings
from ...analysis import evidence_for, geometry as geometry_analysis, run_all
from ...diagnosis import build_report, diagnose, references_section, to_markdown
from ...parsers import ParserError, parse_auto
from ..session import Session
# --- 1 check-up and characterisation ----------------------------------------


def report_output(session: Session) -> None:
    path_text = session.ask("ORCA output file path")
    if not path_text:
        session.say("Cancelled (no path given).")
        return
    path = Path(path_text)
    try:
        result = parse_auto(path)
    except ParserError as exc:
        session.say(f"Parse failed: {exc}")
        return
    sections = run_all(result)
    findings = diagnose(result)
    # analysis modules carry their own evidence (literature thresholds etc.) -
    # merge it so the References block covers analysis citations as well
    analysis_evidence = evidence_for(result)
    report = build_report(
        findings,
        subject=str(path),
        sections=sections,
        extra_evidence=analysis_evidence,
    )
    payload = _report_payload(
        report,
        {
            "program": result.program,
            "input_path": str(path),
            "version": result.sections.get("version", ""),
            "terminated_normally": result.sections.get("terminated_normally"),
            "final_energy": result.sections.get("final_energy"),
        },
    )
    md_path, json_path = _write_report_files(path, to_markdown(report), payload)
    session.say(f"Report written: {md_path}")
    session.say(f"Data written: {json_path}")
    session.say(f"Analysis sections: {len(sections)}.")
    _summarize_findings(session, findings)


# --- 2 coordination geometry ------------------------------------------------


def geometry_report(session: Session) -> None:
    path_text = session.ask("XYZ structure file path")
    if not path_text:
        session.say("Cancelled (no path given).")
        return
    point_group = session.ask("Point group (optional)")
    path = Path(path_text)
    try:
        section = geometry_analysis.run(path, cf_point_group=point_group or None)
    except geometry_analysis.StructureError as exc:
        session.say(f"Structure read failed: {exc}")
        return
    session.say(section.body)
    body = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(geometry_analysis.evidence())
    if refs is not None:
        body += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = path.with_name(path.name + ".fbk.md")
    md_path.write_text(body, encoding="utf-8")
    session.say(f"Report written: {md_path}")


# --- 8 / 0 ------------------------------------------------------------------


def save_script(session: Session) -> None:
    # take the snapshot first: the saving interaction itself (the path line) does not
    # enter the script, otherwise replaying it would carry one line too many
    snapshot = session.snapshot()
    path_text = session.ask("Script save path", default="fbk_session.txt")
    if path_text is None:
        session.say("Cancelled (input finished).")
        return
    target = Path(path_text)
    target.write_text("\n".join(snapshot) + "\n", encoding="utf-8")
    session.say(f"Saved {len(snapshot)} input line(s) to: {target}")




def quit_session(session: Session) -> None:
    session.say("Goodbye.")
    session.stop = True


