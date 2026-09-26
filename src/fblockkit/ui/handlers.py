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

from ..analysis import cf_declaration, crystal_field, point_charge
from ..analysis import entropy_rdm
from ..analysis import evidence_for, run_all
from ..analysis import geometry as geometry_analysis
from ..diagnosis import SEVERITY_LABELS, build_report, diagnose, references_section, to_markdown
from ..diagnosis import CrossLevelError, cross_level_check, load_records
from ..diagnosis import ScfRescueError, propose_fixes
from ..diagnosis import triage as scf_triage
from ..knowledge.models import SystemProfile
from ..parsers import ParserError, parse_auto
from ..parsers.fcidump import parse_fcidump
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
    session.say(run_guidance(recommendation, plan, basis_entry))


# --- 4 basis query ----------------------------------------------------------


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
    except (ParserError, entropy_rdm.EntropyRdmError, OSError, ValueError) as exc:
        session.say(f"Exact entropy analysis failed: {exc}")
        return

    session.say(section.body)
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(entropy_rdm.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = Path(fcidump_text).with_name(Path(fcidump_text).name + ".fbk.md")
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
    "quit": quit_session,
}
