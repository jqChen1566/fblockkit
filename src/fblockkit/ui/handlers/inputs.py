"""The inputs handler group."""

from __future__ import annotations

from pathlib import Path

from ...diagnosis import CrossLevelError, SEVERITY_LABELS, ScfRescueError, cross_level_check, load_records, propose_fixes, triage as scf_triage
from ...knowledge.models import SystemProfile
from ...parsers import ParserError, parse_auto
from ...analysis import geometry as geometry_analysis
from ...recipe import deltascf as deltascf_recipe
from ...recipe import (
    BasisAdvice,
    BasisDataError,
    RenderError,
    plan_convergence,
    recommend_basis_ecp,
    recommend_with_advice,
    render_orca_input,
    run_guidance,
)
from ...toolindex import ToolIndexError, guide as tool_guide_text, search as tool_search_fn
from ..session import Session
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
    from ...knowledge.elements import is_f_element
    from ...recipe.active_space import ActiveSpaceError, suggest_active_space

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
    from ...knowledge.elements import ElementError, is_f_element
    from ...recipe import plan_relativistic

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
        from ...recipe import plan_dmet, render_dmet

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


def deltascf_generate(session: Session) -> None:
    """Menu 27: write a DeltaSCF / MOM excited-state SCF input."""
    xyz_text = session.ask("Structure file (XYZ) path")
    if not xyz_text:
        session.say("Cancelled (no structure given).")
        return
    charge_text = session.ask("Charge of the reference system (Enter = 0)").strip()
    mult_text = session.ask("Multiplicity (Enter = 1)").strip()
    spec_text = session.ask(
        "Occupation spec: an ALPHACONF list (e.g. '0,1' = HOMO->LUMO, '0,0,1' = "
        "HOMO->LUMO+1), or 'ionize N' for IONIZEALPHA"
    ).strip()
    beta_text = session.ask(
        "BETACONF list (Enter = none; give one only with an ALPHACONF spec)"
    ).strip()
    keywords_text = session.ask(
        "Method/basis keywords (Enter = PBE0 def2-TZVP UHF; every source example "
        "uses UHF)"
    ).strip()
    mom_text = session.ask(
        "MOM metric: mom (regular), pmom (the projection-operator variant), or imom "
        "(KeepInitialRef) (Enter = mom)"
    ).strip().lower()
    tactics_text = session.ask(
        "Hard-case tactics: none, freeze (FreezeAndRelease) or gmf (Enter = none)"
    ).strip().lower()
    gbw_text = session.ask(
        "Ground-state gbw path for the MORead start (Enter = skip)"
    ).strip()
    path = Path(xyz_text)
    try:
        charge = int(charge_text) if charge_text else 0
        multiplicity = int(mult_text) if mult_text else 1
        alpha_conf = beta_conf = ionize = None
        tokens = spec_text.lower().split()
        if len(tokens) == 2 and tokens[0] == "ionize":
            ionize = int(tokens[1])
        elif tokens:
            alpha_conf = spec_text
        else:
            raise deltascf_recipe.DeltaScfError("no occupation spec was given.")
        if beta_text:
            beta_conf = beta_text
        atoms = geometry_analysis.parse_xyz(path)
        coordinates = tuple((atom.element, atom.x, atom.y, atom.z) for atom in atoms)
        text = deltascf_recipe.deltascf_input(
            coordinates,
            charge=charge,
            multiplicity=multiplicity,
            alpha_conf=alpha_conf,
            beta_conf=beta_conf,
            ionize_alpha=ionize,
            keywords=keywords_text or "PBE0 def2-TZVP UHF",
            mom=mom_text or "mom",
            tactics=tactics_text or "none",
            gs_gbw=gbw_text or None,
        )
    except (
        deltascf_recipe.DeltaScfError,
        geometry_analysis.StructureError,
        OSError,
        ValueError,
    ) as exc:
        session.say(f"DeltaSCF input failed: {exc}")
        return
    inp_path = path.with_name(f"{path.stem}.dscf.inp")
    inp_path.write_text(text, encoding="utf-8")
    session.say(text)
    session.say(f"Written: {inp_path}. Run it with ORCA, then check the result:")
    for line in deltascf_recipe.run_guidance_lines():
        session.say(f"  - {line}")


def ras_ormas_generate(session: Session) -> None:
    """Menu 28: write a RAS / ORMAS (generalized active space) input."""
    from ...recipe import ras_ormas as ras_ormas_recipe

    xyz_text = session.ask("Structure file (XYZ) path")
    if not xyz_text:
        session.say("Cancelled (no structure given).")
        return
    charge_text = session.ask("Charge of the system (Enter = 0)").strip()
    mult_text_input = session.ask("Multiplicity (Enter = 1)").strip()
    route_text = session.ask(
        "Route: casscf (orbital-optimized RASSCF/ORMAS-SCF) or rasci (CI-only) "
        "(Enter = casscf)"
    ).strip().lower()
    route = route_text or "casscf"
    space_text = session.ask(
        "Partition mask type: ras (three RAS blocks) or ormas (any number of "
        "sub-spaces) (Enter = ormas)"
    ).strip().lower()
    space = space_text or "ormas"
    if space == "ras":
        mask_text = session.ask(
            "RAS mask: '<nel>:<n1> <h1>/<n2>/<n3> <p3>' "
            "(e.g. '6:2 2/2/2 2', the manual's own form)"
        ).strip()
    else:
        mask_text = session.ask(
            "ORMAS mask: '<nel>: <m1> <min1> <max1>, <m2> <min2> <max2>, ...' "
            "(e.g. '6: 2 0 4, 2 0 4, 2 0 4')"
        ).strip()
    state_mult_text = session.ask(
        "Multiplicities (comma list; Enter = the structure's multiplicity)"
    ).strip()
    nroots_text = session.ask("Roots per multiplicity (Enter = 1)").strip() or "1"
    cistep_text = exc_text = ""
    if route == "rasci":
        cistep_text = session.ask(
            "CIStep: accci, csfci, detci, treecsf (Enter = the module's default)"
        ).strip().lower()
        exc_text = session.ask("ExcLevel on top of the references (Enter = the module's default)").strip()
    keywords_text = session.ask(
        "Method/basis keywords (Enter = RHF def2-SVP; the manual's examples use a "
        "plain SCF line)"
    ).strip()
    path = Path(xyz_text)
    try:
        charge = int(charge_text or "0")
        multiplicity = int(mult_text_input or "1")
        atoms = geometry_analysis.parse_xyz(path)
        coordinates = tuple((atom.element, atom.x, atom.y, atom.z) for atom in atoms)
        text = ras_ormas_recipe.ras_ormas_input(
            coordinates,
            charge=charge,
            multiplicity=multiplicity,
            space=space,
            mask=mask_text,
            route=route,
            mult=state_mult_text or None,
            nroots=nroots_text,
            cistep=cistep_text or None,
            exc_level=int(exc_text) if exc_text else None,
            keywords=keywords_text or "RHF def2-SVP",
        )
    except (
        ras_ormas_recipe.RasOrmasError,
        geometry_analysis.StructureError,
        OSError,
        ValueError,
    ) as exc:
        session.say(f"RAS/ORMAS input failed: {exc}")
        return
    suffix = ".rasormas.inp" if route == "casscf" else ".rasci.inp"
    inp_path = path.with_name(f"{path.stem}{suffix}")
    inp_path.write_text(text, encoding="utf-8")
    session.say(text)
    session.say(f"Written: {inp_path}. Run it with ORCA, then check the result:")
    for line in ras_ormas_recipe.run_guidance_lines(route):
        session.say(f"  - {line}")
