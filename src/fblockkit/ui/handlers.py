"""Menu handler functions: every handler takes the session and produces text and files.

Discipline: the user's input files are only ever read; every product is a new file
(``*.fbk.md`` / ``*.fbk.json`` / ``*.fbk.inp``) and never overwrites the original; the
output contains no unstable content such as timestamps (so a script replays
byte-identically).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from ..analysis import atomic_terms, avas, cf_declaration, crystal_field, point_charge
from ..analysis import entropy_rdm, environment_spin
from ..analysis import apc, ass1st, dm_selection, magnetic_doublets, orbital_mapping
from ..analysis import orbital_portrait, qicas
from ..recipe import ass1st as ass1st_recipe
from ..recipe import dm_batch, guess_transfer
from ..analysis import orbital_space as orbital_space_analysis
from ..analysis import evidence_for, run_all
from ..analysis import geometry as geometry_analysis
from ..diagnosis import SEVERITY_LABELS, build_report, diagnose, references_section, to_markdown
from ..diagnosis import CrossLevelError, cross_level_check, load_records
from ..diagnosis import ScfRescueError, propose_fixes
from ..diagnosis import triage as scf_triage
from ..knowledge.models import ReportSection, SystemProfile
from ..parsers import ParserError, parse_auto
from ..parsers.fcidump import parse_fcidump
from ..parsers.mkl import parse_mkl
from ..parsers.orca_json import parse_orca_json
from ..recipe import (
    BasisAdvice,
    BasisDataError,
    RenderError,
    plan_convergence,
    recommend_basis_ecp,
    recommend_with_advice,
    render_orca_input,
    run_guidance,
)
from ..toolindex import ToolIndexError
from ..toolindex import guide as tool_guide_text
from ..toolindex import search as tool_search_fn
from .session import Session


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


# --- 3 generate an input ----------------------------------------------------


def _parse_profile(session: Session) -> SystemProfile | None:
    elements_text = session.ask("Element list (comma separated, e.g. La,O,H)")
    if not elements_text:
        session.say("Cancelled (no element given).")
        return None
    charge = session.ask("System charge", default="0")
    multiplicity = session.ask("Spin multiplicity (2S+1)", default="1")
    valence = session.ask("f-block metal valence")
    targets_text = session.ask(
        "Targets (comma separated: energy/geometry/excited/magnetic/spectra)", default="energy"
    )
    try:
        return SystemProfile(
            elements=tuple(item.strip() for item in elements_text.split(",") if item.strip()),
            charge=int(charge or 0),
            spin=int(multiplicity or 1) - 1,
            metal_valence=int(valence) if valence else 0,
            targets=tuple(item.strip() for item in (targets_text or "").split(",") if item.strip()),
        )
    except ValueError:
        session.say("Numeric parsing failed (charge/multiplicity/valence must be integers). Cancelled.")
        return None


def _print_advice(session: Session, advice: BasisAdvice) -> None:
    if not advice.applicable:
        for note in advice.notes:
            session.say(f"- {note}")
        return
    session.say(f"Basis-set / ECP recommendation ({len(advice.recommendations)} tier(s)):")
    for index, entry in enumerate(advice.recommendations, start=1):
        session.say(f"  {index}. {entry.line()}")
        session.say(f"     {entry.note.strip()}")
    for entry in advice.cautions:
        session.say(f"Note: {entry.note.strip()}")
    for entry in advice.refusals:
        session.say(f"Not applicable: {entry.note.strip()}")
        for refusal in entry.refusals:
            session.say(f"  - {refusal}")
    for note in advice.notes:
        session.say(f"- {note}")


def generate_input(session: Session) -> None:
    profile = _parse_profile(session)
    if profile is None:
        return
    geometry_path = session.ask("Structure file (XYZ) path")
    if not geometry_path:
        session.say("Cancelled (no structure file given -- generating an input needs a geometry).")
        return
    try:
        recommendation, advice = recommend_with_advice(profile)
    except (BasisDataError, ValueError) as exc:
        session.say(f"Recipe layer failed: {exc}")
        return
    session.say(f"Method chain: {' / '.join(recommendation.method_chain) or '(undecided)'}")
    _print_advice(session, advice)

    basis_entry = None
    if advice.recommendations:
        choice = session.ask(f"Choose a basis tier (1-{len(advice.recommendations)})", default="1")
        try:
            index = int(choice or "1") - 1
            basis_entry = advice.recommendations[index]
        except (ValueError, IndexError):
            session.say(f"Invalid choice ({choice!r}) -- falling back to tier 1.")
            basis_entry = advice.recommendations[0]

    # G2: offer a starting active space before asking (suggestion only)
    from ..knowledge.elements import is_f_element
    from ..recipe.active_space import ActiveSpaceError, suggest_active_space

    for element in profile.elements:
        if not is_f_element(element):
            continue
        try:
            suggestions = suggest_active_space(element, profile.metal_valence or 3)
        except ActiveSpaceError:
            break
        session.say(f"Starting active-space suggestions (G2, suggestions only; {element}):")
        for suggestion in suggestions:
            mark = "confirmed" if suggestion.confidence == "confirmed" else "provisional"
            session.say(f"  {suggestion.label}: {suggestion.casscf_line()} [{mark}]")
            session.say(f"    {suggestion.rationale}")
        break

    casscf = None
    casscf_text = session.ask("Active space nel,norb,mult,nroots")
    if casscf_text:
        fields = [item.strip() for item in casscf_text.split(",")]
        if len(fields) != 4 or not all(item.isdigit() for item in fields):
            session.say(f"Active-space format is wrong ({casscf_text!r}) -- no %casscf block will be generated this time.")
        else:
            nel, norb, mult, nroots = (int(item) for item in fields)
            casscf = {"nel": nel, "norb": norb, "mult": mult, "nroots": nroots}

    difficulty = session.ask("Convergence tier (default = defaults first; difficult = use TRAH)", default="default")
    if difficulty not in ("default", "difficult"):
        session.say(f"Invalid convergence tier ({difficulty!r}) -- using default.")
        difficulty = "default"
    auxiliary = ""
    if difficulty == "difficult":
        auxiliary = session.ask("The /C-type auxiliary basis TRAH needs") or ""
        if not auxiliary:
            session.say("No /C auxiliary basis given -- the manual does not allow TRAH without one, switching to the default convergence tier.")
            difficulty = "default"
    maxcore_text = session.ask(
        "MaxCore per process in MB (Enter = 2000; a measured f-block TRAH-CASSCF "
        "needed 9345 MB)",
        default="2000",
    )
    try:
        maxcore = int(maxcore_text or "2000")
        if maxcore <= 0:
            raise ValueError
    except ValueError:
        session.say(f"Invalid MaxCore ({maxcore_text!r}) -- using 2000 MB.")
        maxcore = 2000
    plan = plan_convergence(
        difficulty=difficulty,
        pt2=any("NEVPT2" in m for m in recommendation.method_chain),
        maxcore=maxcore,
    )
    keywords = ["TightSCF"]
    if any("NEVPT2" in m for m in recommendation.method_chain):
        keywords.append("NEVPT2")

    geometry_text = Path(geometry_path).read_text(encoding="utf-8", errors="replace")
    # the first two XYZ lines (atom count and comment) are not part of the coordinate block
    geometry_body = "\n".join(
        line for line in geometry_text.splitlines()[2:] if line.strip()
    )
    try:
        text = render_orca_input(
            recommendation,
            geometry_body,
            profile.charge,
            profile.spin + 1,
            basis_entry=basis_entry,
            method_keywords=keywords,
            casscf=casscf,
            convergence=plan,
            auxiliary=auxiliary,
            maxcore=maxcore,
        )
    except RenderError as exc:
        session.say(f"Rendering failed: {exc}")
        return
    target = Path(geometry_path).with_name(Path(geometry_path).stem + ".fbk.inp")
    target.write_text(text, encoding="ascii")
    session.say(f"Input file written: {target}")
    session.say(
        run_guidance(
            recommendation, plan, basis_entry, _relativistic_for(profile)
        )
    )


# --- 4 basis query ----------------------------------------------------------


def _relativistic_for(profile):
    """The relativistic tier the profile's needs imply (None when nothing is asked)."""
    from ..knowledge.elements import ElementError, is_f_element
    from ..recipe import plan_relativistic

    f_block = False
    for symbol in profile.elements:
        try:
            f_block = f_block or is_f_element(symbol)
        except ElementError:
            continue
    targets = {target.strip().lower() for target in profile.targets}
    needs_soc = bool(targets & {"magnetic", "spectra"})
    if not (f_block or needs_soc):
        return None
    return plan_relativistic(needs_soc=needs_soc, f_block=f_block)


def basis_query(session: Session) -> None:
    profile = _parse_profile(session)
    if profile is None:
        return
    _print_advice(session, recommend_basis_ecp(profile))


# --- 5/6 tool index ---------------------------------------------------------


def tool_search(session: Session) -> None:
    query = session.ask("Search keywords (several words, space separated)")
    if not query:
        session.say("Cancelled (no keyword given).")
        return
    try:
        hits = tool_search_fn(query)
    except ToolIndexError as exc:
        session.say(f"Index loading failed: {exc}")
        return
    if not hits:
        session.say(f"No hits (keywords: {query}).")
        return
    for record in hits:
        session.say(f"- {record.line()}")
        session.say(f"  {record.purpose}")
        if record.note:
            session.say(f"  Note: {record.note.strip()}")


def tool_guide(session: Session) -> None:
    tool_id = session.ask("Tool index id (e.g. openmolcas; search with menu 5)")
    if not tool_id:
        session.say("Cancelled (no id given).")
        return
    try:
        session.say(tool_guide_text(tool_id))
    except ToolIndexError as exc:
        session.say(f"{exc}")
        return
    recipe = _tool_recipe(tool_id)
    if recipe is not None:
        session.say(recipe)


def _tool_recipe(tool_id: str) -> str | None:
    """The workflow recipe for the tools that have one (menu 6 prints it with the guide)."""
    if tool_id.strip().lower() == "liblan":
        from ..recipe import plan_dmet, render_dmet

        plan = plan_dmet("Dy")
        return (
            render_dmet(plan)
            + "\n\n(printed for Dy as the example: the cluster CAS size follows the f "
            "count of the centre, (f-count)e,7o for any Ln3+/An3+ ion)"
        )
    return None


# --- 7 cross-level consistency ----------------------------------------------


def cross_level(session: Session) -> None:
    """Compare solutions across theory levels (A7) from a records JSON file."""
    path_text = session.ask("Records JSON path (each item carries label and energies; see section 7 of the guide)")
    if not path_text:
        session.say("Cancelled (no path given).")
        return
    try:
        records = load_records(path_text)
    except CrossLevelError as exc:
        session.say(f"Read failed: {exc}")
        return
    findings = cross_level_check(records)
    if not findings:
        session.say("The two levels agree: no sign of a cross-level inconsistency.")
        return
    session.say(f"{len(findings)} cross-level risk(s) found:")
    for finding in findings:
        session.say(f"[{SEVERITY_LABELS.get(finding.severity, finding.severity)}] {finding.message}")
        session.say(f"  Action: {finding.suggested_fix}")
        for item in finding.evidence:
            session.say(f"  Evidence: [{item.kind}] {item.ref}")


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


# --- 11 point-charge crystal-field estimate ---------------------------------


def point_charge_estimate(session: Session) -> None:
    """S2: point-charge crystal-field estimate from a structure (XYZ + charges)."""
    path_text = session.ask("XYZ structure file path")
    if not path_text:
        session.say("Cancelled (no path given).")
        return
    path = Path(path_text)
    try:
        atoms = geometry_analysis.parse_xyz(path)
    except geometry_analysis.StructureError as exc:
        session.say(f"Structure read failed: {exc}")
        return
    centre = atoms[geometry_analysis.dominant_center(atoms)]
    session.say(
        f"Centre: {centre.element} at ({centre.x}, {centre.y}, {centre.z}) "
        "(the dominant centre of the structure)"
    )
    charges_text = session.ask("Point charges per element (e.g. O=-2,H=0.4)")
    charges: dict[str, float] = {}
    for chunk in (charges_text or "").split(","):
        if not chunk.strip():
            continue
        symbol, separator, value = chunk.partition("=")
        if not separator:
            session.say(f"Charge entry {chunk.strip()!r} is not element=charge; skipped.")
            continue
        try:
            charges[symbol.strip()] = float(value)
        except ValueError:
            session.say(f"Charge value {value.strip()!r} is not a number; skipped.")
    radial_text = session.ask(
        "Radial moments r2,r4,r6 in Angstrom^k (Enter = geometry-only output)", default=""
    )
    radial: dict[int, float] | None = None
    if radial_text:
        values = [item.strip() for item in radial_text.split(",")]
        if len(values) != 3:
            session.say("Expected three values (r2,r4,r6); giving geometry-only output.")
        else:
            try:
                radial = {k: float(v) for k, v in zip((2, 4, 6), values)}
            except ValueError:
                session.say("Radial moments must be numbers; giving geometry-only output.")
                radial = None
    try:
        estimate = point_charge.estimate(
            atoms, charges, (centre.x, centre.y, centre.z), radial
        )
        section = point_charge.report(estimate, charges, radial)
    except point_charge.PointChargeError as exc:
        session.say(f"Point-charge estimate refused: {exc}")
        return
    session.say(section.body)
    body = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(point_charge.evidence())
    if refs is not None:
        body += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = path.with_name(path.name + ".fbk.md")
    md_path.write_text(body, encoding="utf-8")
    session.say(f"Report written: {md_path}")


def quit_session(session: Session) -> None:
    session.say("Goodbye.")
    session.stop = True


# --- 9 SCF rescue -----------------------------------------------------------


def scf_rescue(session: Session) -> None:
    """Triage the SCF of one output file and, optionally, write corrected inputs."""
    output_text = session.ask("ORCA output file path")
    if not output_text:
        session.say("Cancelled (no path given).")
        return
    try:
        result = parse_auto(output_text)
    except ParserError as exc:
        session.say(f"Parse failed: {exc}")
        return
    block = result.sections.get("scf", {}).get("convergence_block") or ()
    if block:
        # the verbatim view: the check mode decides which rows are enforced, so the
        # informational rows are shown exactly as printed and never interpreted here
        session.say(
            "SCF CONVERGENCE block, verbatim from the output (printed for reference; "
            "the check mode decides which rows are actually enforced):"
        )
        for line in block:
            session.say(line)
        session.say("")
    try:
        findings = scf_triage(result)
    except ScfRescueError as exc:
        session.say(f"SCF triage failed: {exc}")
        return
    if not findings:
        session.say("SCF looks healthy: no triage finding.")
        return
    for finding in findings:
        label = SEVERITY_LABELS.get(finding.severity, finding.severity)
        session.say(f"[{label}] {finding.message}")
        session.say(f"  action: {finding.suggested_fix}")
        for refusal in finding.refusals:
            session.say(f"  do not: {refusal}")
    input_text = session.ask("Matching input file path (Enter = skip corrected inputs)")
    if not input_text:
        return
    source = Path(input_text)
    try:
        original = source.read_text(encoding="utf-8", errors="replace")
        proposals = propose_fixes(original, findings)
    except (OSError, ScfRescueError) as exc:
        session.say(f"Corrected-input generation failed: {exc}")
        return
    if not proposals:
        session.say("No corrected input proposed (the input already contains the recommended settings).")
        return
    for proposal in proposals:
        target = source.with_name(f"{source.stem}.fix_{proposal.name}.inp")
        target.write_text(proposal.content, encoding="utf-8")
        session.say(f"Corrected input written: {target}")
        session.say(f"  rationale: {proposal.rationale}")
        for change in proposal.changes:
            session.say(f"  change: {change}")


# --- 10 crystal-field fit ---------------------------------------------------


def _coefficient_value(entry) -> complex:
    """A coefficient may be a real number or a [re, im] pair."""
    if isinstance(entry, (int, float)):
        return complex(entry)
    if isinstance(entry, (list, tuple)) and len(entry) == 2:
        return complex(entry[0], entry[1])
    raise ValueError(f"coefficient {entry!r} must be a number or a [re, im] pair")


def crystal_field_fit(session: Session) -> None:
    """Fit B_k^q from a levels + coefficients JSON file and write a report."""
    path_text = session.ask("CF input JSON path (levels + coefficients; see the user guide)")
    if not path_text:
        session.say("Cancelled (no path given).")
        return
    path = Path(path_text)
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        point_group = str(data["point_group"])
        j_value = float(data["J"])
        levels = [float(value) for value in data["levels"]]
        coefficients = [
            [_coefficient_value(entry) for entry in row] for row in data["coefficients"]
        ]
        result = crystal_field.fit_crystal_field(levels, coefficients, point_group, j_value)
    except (OSError, KeyError, TypeError, ValueError) as exc:
        session.say(
            f"CF fit failed: {exc}. Next step: provide a JSON with keys "
            "point_group, J, levels (list) and coefficients (list of rows of "
            "numbers or [re, im] pairs)."
        )
        return

    lines = [
        f"Point group {result.point_group}, J = {result.j}",
        f"Sampled states: {result.n_states}; fitted parameters: {result.n_parameters}",
        f"Condition number: {result.condition_number:.3g}; rank {result.rank}",
        f"max |residual| = {result.max_abs_residual:.4g} (caller's energy unit)",
        "",
        "B_k^q:",
    ]
    for (k, q), value in sorted(result.parameters.items()):
        lines.append(f"  ({k},{q:>2})  {value: .10g}")
    lines.append(f"  const    {result.const: .10g}")
    for note in result.notes:
        lines.append(f"note: {note}")
    body = "\n".join(lines)
    session.say(body)

    # A5: projection-basis declaration check on the same JSON (the five items any
    # B_k^q set must record before it may be compared with another scheme)
    a5_payload: dict[str, Any] = {
        "parameters": {f"{k},{q}": value for (k, q), value in result.parameters.items()},
        "declaration": data.get("declaration") or {},
    }
    if isinstance(data.get("compare"), dict):
        a5_payload["compare"] = data["compare"]
    a5_section = None
    try:
        a5_section = cf_declaration.run(a5_payload)
    except cf_declaration.DeclarationError as exc:
        session.say(f"Projection-basis declaration check refused: {exc}")

    report_lines = f"## A4 crystal-field fit\n\n{body}\n"
    if a5_section is not None:
        report_lines += f"\n## {a5_section.title}\n\n{a5_section.body}\n"
    refs = references_section(crystal_field.evidence() + cf_declaration.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = path.with_name(path.name + ".fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")
    if a5_section is not None:
        declared = sum(1 for item in cf_declaration.check_declaration(a5_payload["declaration"]) if item.present)
        session.say(
            f"Projection-basis declaration: {declared} of "
            f"{len(cf_declaration.DECLARATION_FIELDS)} item(s) recorded (section A5 in the report)."
        )
        if "compare" in a5_payload:
            refusals = [
                item
                for item in cf_declaration.compare_parameters(
                    {
                        "parameters": a5_payload["parameters"],
                        "declaration": a5_payload["declaration"],
                    },
                    a5_payload["compare"],
                )
                if item.verdict == cf_declaration.VERDICT_REFUSED
            ]
            session.say(
                f"Parameter comparison against the bundled set: "
                f"{len(refusals)} (k, q) refused without a matching declaration."
            )


def exact_entropy(session: Session) -> None:
    """Menu 12: the exact four-state entropy from a CASSCF output + its FCIDUMP.

    The route and its cross-checks live in ``analysis.entropy_rdm``; this handler
    only collects the four paths, infers (or accepts) the active-orbital window,
    and writes the report next to the FCIDUMP.
    """
    output_text = session.ask(
        "Converged CASSCF output path (the run whose orbitals were dumped)"
    )
    if not output_text:
        session.say("Cancelled (no output path given).")
        return
    fcidump_text = session.ask("FCIDUMP path (written by the !FCIDUMP dump run)")
    if not fcidump_text:
        session.say("Cancelled (no FCIDUMP path given).")
        return
    canonical_text = session.ask(
        "orca_2json export of the canonical gbw (Enter = skip the localized-basis step)"
    )
    localized_text = ""
    if canonical_text:
        localized_text = session.ask(
            "orca_2json export after orca_loc (the localized step needs both exports)"
        )
    window_text = session.ask(
        "Active window 'first last' in ORCA's 0-based orbital numbering (Enter = infer)"
    )
    cluster_text = (
        session.ask(
            "Cluster centres for the environment-spin partition (comma-separated atom "
            "indices; Enter = the f-block centre, if any)"
        )
        if localized_text
        else ""
    )
    try:
        result = parse_auto(Path(output_text))
        casscf = result.sections.get("casscf", {})
        if not casscf.get("present"):
            session.say(
                "This output has no CASSCF section. Next step: use the converged "
                "CASSCF run whose orbitals the FCIDUMP was dumped from."
            )
            return
        if not casscf.get("converged"):
            session.say(
                "The CASSCF did not report convergence, so the dumped Hamiltonian "
                "describes intermediate orbitals and no engine number can validate "
                "the reconstruction. Next step: converge the CASSCF (see the "
                "convergence guidance) and dump again."
            )
            return
        energy = casscf.get("energy")
        if energy is None:
            session.say(
                "The CASSCF section carries no final energy to check against. "
                "Next step: use the converged run's output file."
            )
            return
        states = casscf.get("states") or ()
        multiplicity = states[0].get("mult") if states else None
        printed_occ = tuple(casscf.get("active_occupations") or ())
        table_occ = tuple(result.sections.get("orbitals", {}).get("occupations") or ())
        dump = parse_fcidump(Path(fcidump_text))
        if window_text:
            try:
                first, last = (
                    int(token) for token in window_text.replace(",", " ").split()
                )
            except ValueError:
                session.say(
                    "The window must be two integers, e.g. '4 9'. Cancelled."
                )
                return
            window: tuple[int, ...] | None = tuple(range(first, last + 1))
        else:
            window, note = entropy_rdm.infer_active_window(
                table_occ, dump.norb, reference_occupations=printed_occ or None
            )
            session.say(f"Active window inferred: {list(window)} ({note}).")
        canonical = (
            parse_orca_json(Path(canonical_text)) if canonical_text else None
        )
        localized = (
            parse_orca_json(Path(localized_text)) if localized_text else None
        )
        analysis = entropy_rdm.analyze(
            dump,
            reference_energy=float(energy),
            multiplicity=multiplicity,
            reference_occupations=printed_occ or None,
            canonical=canonical,
            localized=localized,
            active_window=window,
        )
        section = entropy_rdm.run(analysis)
        extras: list[tuple[Any, tuple[Any, ...]]] = []
        if canonical is not None and analysis.active_window is not None:
            term_sections = atomic_terms.analyze_export(
                analysis.densities,
                analysis.state.s2,
                canonical,
                analysis.active_window,
            )
            if term_sections:
                extras.append((term_sections, atomic_terms.evidence()))
        if localized is not None and analysis.active_window is not None:
            cluster_centres = _cluster_centres(cluster_text, canonical)
            if cluster_centres is None:
                session.say(
                    "Environment-spin partition skipped: no cluster centre was given and "
                    "the molecule has no f-block element. Next step: give the cluster "
                    "centre index when the environment-spin entropy is wanted."
                )
            else:
                rotation = entropy_rdm.rotation_from_coefficients(
                    np.array(canonical.mo_coefficients).T,
                    np.array(localized.mo_coefficients).T,
                    np.array(canonical.overlap),
                    active=analysis.active_window,
                )
                rotated = entropy_rdm.rotate_densities(analysis.densities, rotation)
                local_coefficients = np.array(localized.mo_coefficients).T[
                    :, list(analysis.active_window)
                ]
                extras.append(
                    (
                        (
                            environment_spin.analyze(
                                rotated,
                                local_coefficients,
                                np.array(localized.overlap),
                                localized.ao_labels,
                                len(localized.atoms),
                                cluster_centres=cluster_centres,
                            ),
                        ),
                        environment_spin.evidence(),
                    )
                )
    except (ParserError, entropy_rdm.EntropyRdmError, OSError, ValueError) as exc:
        session.say(f"Exact entropy analysis failed: {exc}")
        return

    session.say(section.body)
    report_lines = f"## {section.title}\n\n{section.body}\n"
    evidence_entries = list(entropy_rdm.evidence())
    for sections, entries in extras:
        for extra in sections:
            session.say(extra.body)
            report_lines += f"\n## {extra.title}\n\n{extra.body}\n"
        evidence_entries.extend(entries)
    refs = references_section(tuple(evidence_entries))
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = Path(fcidump_text).with_name(Path(fcidump_text).name + ".fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")


def _cluster_centres(text: str, canonical) -> tuple[int, ...] | None:
    """Parse the cluster-centre prompt; None when there is nothing to partition on."""
    from ..knowledge.elements import is_f_element

    if text:
        try:
            return tuple(int(token) for token in text.replace(",", " ").split())
        except ValueError:
            return None
    if canonical is not None:
        detected = tuple(
            index for index, symbol in enumerate(canonical.atoms) if is_f_element(symbol)
        )
        if detected:
            return detected
    return None


# --- 13 orbital-space comparison ---------------------------------------------


def orbital_space(session: Session) -> None:
    """Menu 13: sigma_F and the space-change SVD from two orca_2json exports.

    Both readings come from the singular values of C_A^T S C_B; the handler only
    collects the two exports and their orbital windows (ORCA's 0-based numbering,
    the same convention as menu 12) and writes the report next to the first
    export.
    """
    first_text = session.ask("First orca_2json export (space A) path")
    if not first_text:
        session.say("Cancelled (no export path given).")
        return
    first_window_text = session.ask(
        "Window 'first last' of space A in ORCA's 0-based orbital numbering "
        "(Enter = all orbitals)"
    )
    second_text = session.ask("Second orca_2json export (space B) path")
    if not second_text:
        session.say("Cancelled (no second export path given).")
        return
    second_window_text = session.ask(
        "Window 'first last' of space B (Enter = all orbitals)"
    )
    try:
        first = parse_orca_json(Path(first_text))
        second = parse_orca_json(Path(second_text))
        overlap = np.array(first.overlap)
        if first.overlap is None or second.overlap is None:
            session.say(
                "One of the exports carries no S-Matrix, so the overlap of the two "
                "spaces cannot be formed. Next step: re-export both gbw files with "
                "orca_2json (the S-Matrix is part of the default export)."
            )
            return
        coefficients_a = np.array(first.mo_coefficients).T[
            :, _window(first_window_text, first.n_mo)
        ]
        coefficients_b = np.array(second.mo_coefficients).T[
            :, _window(second_window_text, second.n_mo)
        ]
        if first.n_ao != second.n_ao or not np.allclose(
            overlap, np.array(second.overlap), atol=1e-8
        ):
            session.say(
                "The two exports do not share one basis: their AO dimensions or overlap "
                "matrices differ. Next step: compare spaces of the same system in the "
                "same basis set (e.g. before and after orca_loc)."
            )
            return
        comparison = orbital_space_analysis.compare_spaces(coefficients_a, coefficients_b, overlap)
        window_a = _window(first_window_text, first.n_mo)
        window_b = _window(second_window_text, second.n_mo)
        section = orbital_space_analysis.run(
            comparison,
            label_a=_space_label(first, first_window_text),
            label_b=_space_label(second, second_window_text),
            jaccard_index=orbital_space_analysis.jaccard(window_a, window_b),
        )
    except (ParserError, orbital_space_analysis.OrbitalSpaceError, OSError, ValueError) as exc:
        session.say(f"Orbital-space comparison failed: {exc}")
        return

    session.say(section.body)
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(orbital_space_analysis.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = Path(first_text).with_name(Path(first_text).name + ".fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")


def _window(text: str, n_mo: int) -> list[int]:
    """A 'first last' window in ORCA's 0-based numbering; empty text means everything."""
    if not text:
        return list(range(n_mo))
    first, last = (int(token) for token in text.replace(",", " ").split())
    return list(range(first, last + 1))


def _space_label(export, window_text: str) -> str:
    span = window_text.strip() or f"0..{export.n_mo - 1}"
    return f"{Path(export.base_name).name} [{span}]"



# --- 14 AVAS target projection -----------------------------------------------


def avas_target(session: Session) -> None:
    """Menu 14: the AVAS projection of a target AO shell onto an orbital export."""
    export_text = session.ask("orca_2json export path (any gbw export)")
    if not export_text:
        session.say("Cancelled (no export path given).")
        return
    centre_text = session.ask(
        "Target centre (atom index, 0-based; Enter = the f-block element, when there is one)"
    )
    angular_text = session.ask(
        "Target angular momentum as a letter (s/p/d/f/g; Enter = f)", default="f"
    )
    shells_text = session.ask(
        "Target shell number(s), comma separated as ORCA labels them (Enter = every shell "
        "of that angular momentum)"
    )
    threshold_text = session.ask(
        f"Truncation threshold (Enter = {avas.DEFAULT_THRESHOLD:g}; the source's range is 0.05-0.1)"
    )
    option_text = session.ask(
        "Open-shell option (2 = alpha orbitals only; Enter = 3, keep every singly "
        "occupied orbital)"
    )
    letters = {letter: value for value, letter in _ANGULAR_LETTERS().items()}
    try:
        export = parse_orca_json(Path(export_text))
        centre = int(centre_text) if centre_text else None
        key = (angular_text or "f").strip().lower()
        if key.isdigit():
            angular = int(key)
        elif key in letters:
            angular = letters[key]
        else:
            session.say(f"{angular_text!r} is not an angular-momentum letter (s/p/d/f/g).")
            return
        shells = (
            {int(token) for token in shells_text.replace(",", " ").split()}
            if shells_text
            else None
        )
        threshold = float(threshold_text) if threshold_text else avas.DEFAULT_THRESHOLD
        option = int(option_text) if option_text else avas.DEFAULT_OPTION
        analysis = avas.analyze(
            export,
            centre=centre,
            angular=angular,
            shells=shells,
            threshold=threshold,
            option=option,
        )
        section = avas.run(analysis)
    except (ParserError, avas.AvasError, OSError, ValueError) as exc:
        session.say(f"AVAS projection failed: {exc}")
        return

    session.say(section.body)
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(avas.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = Path(export_text).with_name(Path(export_text).name + ".avas.fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")


def _ANGULAR_LETTERS() -> dict:
    from ..parsers.orca_json import ANGULAR_LETTERS

    return dict(ANGULAR_LETTERS)



# --- 15 orbital portrait -----------------------------------------------------


def orbital_portrait_report(session: Session) -> None:
    """Menu 15: the deterministic descriptor panel of an orbital set (an export in)."""
    export_text = session.ask("orca_2json export path")
    if not export_text:
        session.say("Cancelled (no export path given).")
        return
    window_text = session.ask(
        "Window 'first last' in ORCA's 0-based numbering (Enter = the orbitals with "
        "fractional occupations)"
    )
    try:
        export = parse_orca_json(Path(export_text))
        window = None
        if window_text:
            first, last = (int(token) for token in window_text.replace(",", " ").split())
            window = list(range(first, last + 1))
        section = orbital_portrait.run(export, window=window)
    except (ParserError, orbital_portrait.PortraitError, OSError, ValueError) as exc:
        session.say(f"Orbital portrait failed: {exc}")
        return
    session.say(section.body)
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(orbital_portrait.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = Path(export_text).with_name(Path(export_text).name + ".portrait.fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")


# --- 16 magnetic doublets ----------------------------------------------------


def magnetic_doublets_report(session: Session) -> None:
    """Menu 16: the g_T * theta_3 criterion over a Kramers-doublet table (JSON in)."""
    path_text = session.ask(
        "Kramers-doublet table JSON path (per doublet: the three g values and theta3 "
        "or the g3 axis; see the user guide)"
    )
    if not path_text:
        session.say("Cancelled (no path given).")
        return
    path = Path(path_text)
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        table = magnetic_doublets.parse_doublets(payload)
        section = magnetic_doublets.run(table)
    except (OSError, ValueError) as exc:
        session.say(
            f"Magnetic-doublet reading failed: {exc}"
        )
        return
    session.say(section.body)
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(magnetic_doublets.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = path.with_name(path.name + ".magnetic.fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")


# --- 17 cross-structure orbital mapping --------------------------------------


def orbital_mapping_report(session: Session) -> None:
    """Menu 17: map the orbitals of a structure series and make a selection consistent."""
    path_text = session.ask(
        "Structures manifest JSON path (a list of structures with their exports; see the "
        "user guide)"
    )
    if not path_text:
        session.say("Cancelled (no path given).")
        return
    path = Path(path_text)
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        entries = payload["structures"]
        exports = []
        names = []
        structures = []
        for entry in entries:
            export_path = Path(str(entry["export"]))
            if not export_path.is_absolute():
                export_path = path.parent / export_path
            export = parse_orca_json(export_path)
            name = str(entry.get("name") or export_path.name)
            exports.append(export)
            names.append(name)
            structures.append(orbital_mapping.descriptors_from_export(export, name))
        tau = payload.get("tau")
        selections = payload.get("selections")
        # the optional active-space overlap check runs over adjacent pairs
        active_overlap = None
        actives = payload.get("active")
        if actives:
            active_overlap = orbital_mapping.active_overlap_series(exports, actives, names)
        section = orbital_mapping.run(
            structures,
            tau=float(tau) if tau is not None else None,
            selections=selections,
            active_overlap=active_overlap,
        )
    except (OSError, KeyError, TypeError, ValueError, ParserError) as exc:
        session.say(f"Orbital mapping failed: {exc}")
        return
    session.say(section.body)
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(orbital_mapping.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = path.with_name(path.name + ".mapping.fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")


# --- 18 WASP guess transfer --------------------------------------------------


def wasp_guess(session: Session) -> None:
    """Menu 18: interpolate neighbour orbitals into a gbw-ready mkl (WASP)."""
    path_text = session.ask(
        "Series manifest JSON path (neighbours and the template mkl; see the user guide)"
    )
    if not path_text:
        session.say("Cancelled (no path given).")
        return
    path = Path(path_text)

    def _resolve(entry) -> Path:
        candidate = Path(str(entry))
        return candidate if candidate.is_absolute() else path.parent / candidate

    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        neighbours = [
            parse_orca_json(_resolve(entry["export"])) for entry in payload["structures"]
        ]
        template_entry = payload["template"]
        template_export = parse_orca_json(_resolve(template_entry["export"]))
        template_mkl_path = _resolve(template_entry["mkl"])
        template_mkl = parse_mkl(template_mkl_path)
        delta = payload.get("delta")
        result = guess_transfer.interpolate_guess(
            neighbours,
            template_export,
            delta=float(delta) if delta is not None else None,
        )
    except (OSError, KeyError, TypeError, ValueError, ParserError) as exc:
        session.say(f"Guess interpolation failed: {exc}")
        return
    body = guess_transfer.render(result)
    session.say(body)
    guess_path = template_mkl_path.with_name(template_mkl_path.stem + ".fbk.mkl")
    try:
        guess_transfer.write_guess_mkl(template_mkl, result, guess_path)
    except (OSError, guess_transfer.GuessError) as exc:
        session.say(f"Writing the guess mkl failed: {exc}")
        return
    session.say(f"Guess written: {guess_path}")
    session.say(
        f"Next step: run ``orca_2mkl {guess_path.stem} -gbw`` next to that file, then "
        f"start the target calculation with ``!moread`` and %moinp \"{guess_path.stem}.gbw\"."
    )
    section = ReportSection(
        title="G5 WASP initial guess (interpolated orbitals for a geometry series)",
        body=body,
    )
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(guess_transfer.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = path.with_name(path.name + ".guess.fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")


# --- 19 dipole-moment candidate batch ----------------------------------------


def dm_batch_generate(session: Session) -> None:
    """Menu 19: write the DM-AS candidate batch (prep + CASCI inputs + manifest)."""
    xyz_text = session.ask("Structure file (XYZ) path")
    if not xyz_text:
        session.say("Cancelled (no structure given).")
        return
    path = Path(xyz_text)
    try:
        structure = path.read_text(encoding="utf-8")
    except OSError as exc:
        session.say(f"Reading the structure failed: {exc}")
        return
    charge_text = session.ask("Charge (Enter = 0)")
    mult_text = session.ask("Multiplicity (Enter = 1; the protocol is built for singlets)")
    basis = session.ask("Basis keyword (Enter = def2-TZVP)") or "def2-TZVP"
    prep = session.ask(
        "Orbital preparation: Enter = MP2 natural orbitals (the source's choice); 'hf' = "
        "plain HF orbitals"
    ) or "mp2"
    max_norb_text = session.ask("Largest active space in orbitals (Enter = 14)")
    extras_text = session.ask(
        "Include the PASS+ extras (the 4-electron and two-virtual rows; Enter = no)"
    )
    try:
        batch = dm_batch.plan_batch(
            xyz_text=structure,
            charge=int(charge_text) if charge_text else 0,
            multiplicity=int(mult_text) if mult_text else 1,
            basis=basis.strip(),
            prep=prep.strip().lower(),
            include_pass_plus=bool(extras_text.strip()),
            max_norb=int(max_norb_text) if max_norb_text else 14,
        )
    except (ValueError, dm_batch.DmBatchError) as exc:
        session.say(f"Batch generation failed: {exc}")
        return
    directory = path.with_name(path.stem + ".fbk.dm")
    try:
        directory.mkdir(exist_ok=True)
        (directory / batch.prep_name).write_text(batch.prep_input(), encoding="utf-8")
        for candidate in batch.candidates:
            (directory / batch.candidate_name(candidate)).write_text(
                batch.candidate_input(candidate), encoding="utf-8"
            )
        (directory / batch.script_name).write_text(batch.render_script(), encoding="utf-8")
        (directory / batch.manifest_name).write_text(
            dm_batch.manifest_text(batch), encoding="utf-8"
        )
    except OSError as exc:
        session.say(f"Writing the batch failed: {exc}")
        return
    session.say(f"DM-AS candidate batch written: {directory}")
    session.say(
        f"  {len(batch.candidates)} candidate input(s) + {batch.prep_name} "
        f"({batch.prep} orbitals) + {batch.script_name} + {batch.manifest_name}"
    )
    session.say("  each candidate: !NoIter moread + %moinp the prep gbw + "
                "%casscf nel/norb/nroots 1 (the form whose S0 dipole ORCA prints)")
    session.say(f"Next step: run the script on the cluster ({batch.script_name}), fill the "
                "reference output into the manifest, then use menu 20.")


# --- 20 dipole-moment selection ----------------------------------------------


def dm_select(session: Session) -> None:
    """Menu 20: rank the candidate spaces by dipole-moment deviation (DM-AS)."""
    path_text = session.ask("Batch manifest JSON path (from menu 19, reference filled in)")
    if not path_text:
        session.say("Cancelled (no manifest given).")
        return
    path = Path(path_text)

    def _resolve(entry) -> Path:
        candidate = Path(str(entry))
        return candidate if candidate.is_absolute() else path.parent / candidate

    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        readings = []
        for entry in payload["candidates"]:
            result = parse_auto(_resolve(entry["output"]))
            magnitude, vector, _ = dm_selection.candidate_dipole(
                result.sections, str(entry["output"])
            )
            readings.append((int(entry["nel"]), int(entry["norb"]), magnitude, vector))
        reference_entry = payload.get("reference")
        if not isinstance(reference_entry, dict):
            raise dm_selection.DmSelectionError(
                "the manifest carries no reference block. Next step: give "
                '{"output": "dft.out"} or {"magnitude_debye": ..., "source": ...}.'
            )
        if "output" in reference_entry:
            result = parse_auto(_resolve(reference_entry["output"]))
            magnitude, vector, method = dm_selection.reference_dipole(
                result.sections, str(reference_entry["output"])
            )
            source = f"{reference_entry['output']} ({method} block)"
        elif "magnitude_debye" in reference_entry:
            magnitude = float(reference_entry["magnitude_debye"])
            vector = reference_entry.get("vector")
            source = str(reference_entry.get("source") or "supplied value")
        else:
            raise dm_selection.DmSelectionError(
                "the reference block carries neither 'output' nor 'magnitude_debye'. "
                'Next step: give one of the two.'
            )
        selection = dm_selection.analyze(
            readings,
            (magnitude, tuple(vector) if vector is not None else None, source),
            protocol=str(payload.get("protocol") or "gdm"),
        )
    except (OSError, KeyError, TypeError, ValueError, ParserError) as exc:
        session.say(f"Dipole-moment selection failed: {exc}")
        return
    body = dm_selection.render(selection)
    session.say(body)
    section = ReportSection(
        title="A12 dipole-moment active-space selection",
        body=body,
    )
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(dm_selection.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = path.with_name(path.name + ".dm_select.fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")


def apc_ranking(session: Session) -> None:
    """Menu 21: rank orbitals by the approximate pair coefficient (APC) and select the
    active space at a CSF cap (an orca_2json export in)."""
    path_text = session.ask(
        "orca_2json export path (needs the FockMatrix K-Matrix block; the user guide "
        "has the exact request)"
    )
    if not path_text:
        session.say("Cancelled (no export path given).")
        return
    variant_text = session.ask("Ranking variant: apc or apcx (Enter = apc)").strip().lower()
    if variant_text not in ("", "apc", "apcx"):
        session.say(f"Cancelled (unknown variant {variant_text!r}; use apc or apcx).")
        return
    window_text = session.ask(
        "Candidate window: the N lowest virtuals in energy (Enter = 23, the source's "
        "general scheme)"
    ).strip()
    cap_text = session.ask(
        "CSF cap: max(7,6), max(8,8), max(10,10), max(12,12) or an integer "
        "(Enter = max(8,8))"
    ).strip()
    delta_text = session.ask(
        "Model gap: energies or fock (Enter = energies; choose fock for localized orbitals)"
    ).strip().lower()
    if delta_text not in ("", "energies", "fock"):
        session.say(f"Cancelled (unknown model-gap source {delta_text!r}).")
        return
    try:
        window_size = int(window_text) if window_text else apc.DEFAULT_WINDOW
    except ValueError:
        session.say(f"Cancelled (the window must be an integer, got {window_text!r}).")
        return
    if not cap_text:
        cap, cap_label = apc.CSF_CAPS["max(8,8)"], "max(8,8)"
    elif cap_text in apc.CSF_CAPS:
        cap, cap_label = apc.CSF_CAPS[cap_text], cap_text
    else:
        try:
            cap, cap_label = int(cap_text), cap_text
        except ValueError:
            session.say(
                f"Cancelled (unknown cap {cap_text!r}; use a preset or an integer)."
            )
            return
    path = Path(path_text)
    try:
        export = parse_orca_json(path)
        section = apc.run(
            export,
            variant="APCX" if variant_text == "apcx" else "APC",
            window_size=window_size,
            delta=delta_text or "energies",
            cap=cap,
            cap_label=cap_label,
        )
    except (ParserError, apc.ApcError, OSError, ValueError) as exc:
        session.say(f"APC ranking failed: {exc}")
        return
    session.say(section.body)
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(apc.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = path.with_name(path.name + ".apc.fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")


def _ass1st_coordinates(atoms) -> tuple[tuple[str, float, float, float], ...]:
    return tuple((atom.element, atom.x, atom.y, atom.z) for atom in atoms)


def ass1st_start(session: Session) -> None:
    """Menu 22: write the ASS1ST round-1 input (structure + initial space in)."""
    xyz_text = session.ask("Structure file (XYZ) path")
    if not xyz_text:
        session.say("Cancelled (no structure given).")
        return
    charge_text = session.ask("Charge (Enter = 0)").strip()
    mult_text = session.ask("Multiplicity (Enter = 1)").strip()
    space_text = session.ask(
        "Initial active space nel,norb (small but chemically reasonable; the user "
        "guide has the source's advice)"
    ).strip()
    states_text = session.ask("States to average, nroots (Enter = 1)").strip()
    keywords_text = session.ask("Method/basis keywords (Enter = RHF def2-SVP TightSCF)").strip()
    maxcore_text = session.ask("MaxCore in MB (Enter = 2000)").strip()
    path = Path(xyz_text)
    try:
        charge = int(charge_text) if charge_text else 0
        multiplicity = int(mult_text) if mult_text else 1
        tokens = [token for token in space_text.replace(",", " ").split() if token]
        if len(tokens) != 2:
            raise ass1st_recipe.Ass1stInputError(
                f"the initial active space {space_text!r} is not 'nel,norb'."
            )
        n_electrons, n_orbitals = (int(token) for token in tokens)
        n_states = int(states_text) if states_text else 1
        maxcore = int(maxcore_text) if maxcore_text else 2000
        atoms = geometry_analysis.parse_xyz(path)
        text = ass1st_recipe.round_one_input(
            _ass1st_coordinates(atoms),
            charge=charge,
            multiplicity=multiplicity,
            n_electrons=n_electrons,
            n_orbitals=n_orbitals,
            n_states=n_states,
            keywords=keywords_text or "RHF def2-SVP TightSCF",
            maxcore=maxcore,
        )
    except (
        ass1st_recipe.Ass1stInputError,
        geometry_analysis.StructureError,
        OSError,
        ValueError,
    ) as exc:
        session.say(f"ASS1ST round-1 input failed: {exc}")
        return
    stem = ass1st_recipe.round_one_stem(path.stem)
    inp_path = path.with_name(f"{stem}.inp")
    conf_path = path.with_name(f"{stem}.json.conf")
    inp_path.write_text(text, encoding="utf-8")
    conf_path.write_text(ass1st_recipe.export_conf(), encoding="utf-8")
    session.say(text)
    session.say(
        f"Written: {inp_path} (the round-1 input) and {conf_path} (the export "
        "request). Next steps:\n"
        f"  1. run it:        orca {inp_path}\n"
        f"  2. export it:     orca_2json {stem}.gbw   (the .json.conf is already "
        "beside the .gbw)\n"
        f"  3. select:        menu 23 with {stem}.json -> the quasi-NOON tables and "
        "the next round's input"
    )


def ass1st_round(session: Session) -> None:
    """Menu 23: analyse one ASS1ST round and write the next round's input."""
    path_text = session.ask(
        "orca_2json export path (a CASSCF + FIC-NEVPT2 round; the user guide has the "
        "exact request)"
    )
    if not path_text:
        session.say("Cancelled (no export path given).")
        return
    band_text = session.ask(
        "NOON band: T (band [T, 2-T]) or T_ext,T_int (Enter = 0.05)"
    ).strip()
    weights_text = session.ask(
        "State weights, comma-separated (Enter = equal; used only for state-averaged "
        "rounds)"
    ).strip()
    prev_text = session.ask(
        "Previously visited spaces 'ne,no' space-separated (Enter = none; enables the "
        "cycle warning)"
    ).strip()
    keywords_text = session.ask(
        "Method/basis keywords for the next round (Enter = RHF def2-SVP TightSCF)"
    ).strip()
    path = Path(path_text)
    try:
        export = parse_orca_json(path)
        band = ass1st.parse_band(band_text) if band_text else (
            ass1st.ASS1ST_RECOMMENDED_THRESHOLD,
            2.0 - ass1st.ASS1ST_RECOMMENDED_THRESHOLD,
        )
        weights = None
        if weights_text:
            tokens = [token for token in weights_text.replace(",", " ").split() if token]
            weights = tuple(float(token) for token in tokens)
        previous = []
        for token in prev_text.split():
            pair = token.replace(",", " ").split()
            if len(pair) != 2:
                raise ass1st.Ass1stError(
                    f"the visited space {token!r} is not 'ne,no'."
                )
            previous.append((int(pair[0]), int(pair[1])))
        round_ = ass1st.analyze_round(
            export, weights=weights, band=band, previous_spaces=tuple(previous)
        )
    except (ParserError, ass1st.Ass1stError, OSError, ValueError) as exc:
        session.say(f"ASS1ST round failed: {exc}")
        return
    section = ass1st.run(export, weights=weights, band=band, previous_spaces=tuple(previous))
    session.say(section.body)
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(ass1st.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = path.with_name(path.name + ".ass1st.fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")
    if round_.suggestion.self_consistent:
        session.say(
            f"Self-consistent at ({round_.partition.n_electrons}e, "
            f"{round_.partition.n_orbitals}o): no next round is needed. Use this "
            "active space for the production calculation."
        )
        return
    suggestion = round_.suggestion
    stem = ass1st_recipe.next_round_stem(export.base_name)
    coordinates = tuple(
        (element, *position) for element, position in zip(export.atoms, export.coordinates)
    )
    text = ass1st_recipe.next_round_input(
        coordinates,
        charge=export.charge,
        multiplicity=export.multiplicity,
        n_electrons=suggestion.n_electrons,
        n_orbitals=suggestion.n_orbitals,
        n_states=round_.n_states,
        keywords=keywords_text or "RHF def2-SVP TightSCF",
    )
    inp_path = path.with_name(f"{stem}.inp")
    conf_path = path.with_name(f"{stem}.json.conf")
    inp_path.write_text(text, encoding="utf-8")
    conf_path.write_text(ass1st_recipe.export_conf(), encoding="utf-8")
    session.say(
        f"Next round written: {inp_path} at ({suggestion.n_electrons}e, "
        f"{suggestion.n_orbitals}o), with {conf_path}. Run it, export "
        f"{stem}.json, and feed that back to this menu; give the visited spaces "
        "(including this round's) to get the cycle warning right."
    )


def qicas_optimize(session: Session) -> None:
    """Menu 24: optimize an active space by the QICAS F_QI minimization (an FCIDUMP in)."""
    fcidump_text = session.ask(
        "FCIDUMP path (a converged CASSCF dump; the window QICAS optimizes within)"
    )
    if not fcidump_text:
        session.say("Cancelled (no FCIDUMP path given).")
        return
    output_text = session.ask(
        "CASSCF output path (matches the CI root to the printed energy; Enter = skip)"
    ).strip()
    space_text = session.ask(
        "Target active space nel,norb within the window (e.g. 4,4; the source's "
        "(N_CAS, D_CAS))"
    ).strip()
    pairs_text = session.ask(
        "Rotation set: touch (every pair touching a non-active orbital, the source's "
        "chemical-accuracy choice) or exclusive (active/non-active only, its economical "
        "variant) (Enter = touch)"
    ).strip().lower()
    if pairs_text not in ("", "touch", "exclusive"):
        session.say(f"Cancelled (unknown rotation set {pairs_text!r}).")
        return
    try:
        reference = None
        if output_text:
            result = parse_auto(Path(output_text))
            casscf = result.sections.get("casscf", {})
            if not casscf.get("present") or not casscf.get("converged"):
                session.say(
                    "The output carries no converged CASSCF section to match the CI "
                    "root against. Next step: give the converged run's output, or leave "
                    "the output prompt empty to use the lowest root of the FCIDUMP's "
                    "Ms sector."
                )
                return
            reference = casscf.get("energy")
        tokens = [token for token in space_text.replace(",", " ").split() if token]
        if len(tokens) != 2:
            raise qicas.QicasError(
                f"the target space {space_text!r} is not 'nel,norb'."
            )
        n_cas, n_act = (int(token) for token in tokens)
        dump = parse_fcidump(Path(fcidump_text))
        section = qicas.run(
            dump,
            n_cas=n_cas,
            n_active_orbitals=n_act,
            pairs_mode=pairs_text or "touch",
            reference_energy=reference,
        )
    except (ParserError, qicas.QicasError, OSError, ValueError) as exc:
        session.say(f"QICAS optimization failed: {exc}")
        return
    session.say(section.body)
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(qicas.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = Path(fcidump_text).with_name(Path(fcidump_text).name + ".qicas.fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")


HANDLERS = {
    "report_output": report_output,
    "geometry_report": geometry_report,
    "generate_input": generate_input,
    "basis_query": basis_query,
    "tool_search": tool_search,
    "tool_guide": tool_guide,
    "cross_level": cross_level,
    "save_script": save_script,
    "scf_rescue": scf_rescue,
    "crystal_field_fit": crystal_field_fit,
    "point_charge_estimate": point_charge_estimate,
    "exact_entropy": exact_entropy,
    "orbital_space": orbital_space,
    "avas_target": avas_target,
    "orbital_portrait": orbital_portrait_report,
    "magnetic_doublets": magnetic_doublets_report,
    "orbital_mapping": orbital_mapping_report,
    "wasp_guess": wasp_guess,
    "dm_batch": dm_batch_generate,
    "dm_select": dm_select,
    "apc_ranking": apc_ranking,
    "ass1st_start": ass1st_start,
    "ass1st_round": ass1st_round,
    "qicas_optimize": qicas_optimize,
    "quit": quit_session,
}
