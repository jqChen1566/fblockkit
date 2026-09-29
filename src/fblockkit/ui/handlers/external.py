"""The external-program handler group (the B-layer menus).

What belongs here: the menus whose counterpart is another program rather than
the engine -- writing the input that program consumes and reading the output
it leaves behind.  Four program families so far: pysisyphus (menus 32/33,
the PES-exploration pair), xTB (menu 42, the pre-screening capture), CREST
(menu 43, the conformer ensemble) and MOKIT/automr (menus 44/45, the
generated input and the run report).
"""

from __future__ import annotations

from pathlib import Path

from ...analysis import crest as crest_analysis
from ...analysis import geometry as geometry_analysis
from ...analysis import mokit as mokit_analysis
from ...analysis import pysisyphus_run as pysisyphus_analysis
from ...analysis import xtb as xtb_analysis
from ...diagnosis import references_section
from ...knowledge.models import ReportSection
from ...parsers import ParserError, parse_auto
from ...parsers.crest import CrestError, read_crest_directory
from ...parsers.mokit import MokitError, read_mokit_run_file
from ...parsers.pysisyphus import PysisyphusError
from ...parsers.xtb import XtbError, read_xtb_run_file
from ...recipe import mokit as mokit_recipe
from ...recipe import pysisyphus as pysisyphus_recipe
from ...recipe.mokit import MokitInputError
from ..session import Session


def _atoms_from_structure_source(path: Path) -> tuple[tuple[str, float, float, float], ...]:
    """Atoms from an XYZ file, an ORCA input (inline block) or an ORCA output."""
    suffix = path.suffix.lower()
    if suffix == ".xyz":
        return tuple(
            (atom.element, atom.x, atom.y, atom.z)
            for atom in geometry_analysis.parse_xyz(path)
        )
    if suffix == ".out":
        result = parse_auto(path)
        section = result.sections.get("final_geometry", {})
        if not section.get("present"):
            raise PysisyphusError(
                "the ORCA output carries no final CARTESIAN COORDINATES block. Next "
                "step: give the output of a job that reached a geometry (or an XYZ "
                "file)."
            )
        return tuple(section["atoms"])
    if suffix == ".inp":
        text = path.read_text(encoding="utf-8", errors="replace")
        return pysisyphus_recipe.atoms_from_orca_input(text)
    raise PysisyphusError(
        f"unsupported structure source {path.name!r}: give an XYZ file (.xyz), an "
        "ORCA input with inline coordinates (.inp) or an ORCA output (.out)."
    )


def pysisyphus_generate(session: Session) -> None:
    """Menu 32: write a pysisyphus run input (YAML + structure) for a chosen method."""
    source_text = session.ask(
        "Structure source path (XYZ / ORCA input / ORCA output)"
    )
    if not source_text:
        session.say("Cancelled (no structure given).")
        return
    source = Path(source_text)
    keywords = session.ask(
        "ORCA keywords for every call (any simple line; Enter = HF def2-SVP)"
    ) or "HF def2-SVP"
    charge_text = session.ask("Charge", default="0")
    mult_text = session.ask("Multiplicity (2S+1)", default="1")
    job = (
        session.ask(
            "Job: min (minimum optimisation) or ts (transition-state search) "
            "(Enter = min)",
            default="min",
        )
        .strip()
        .lower()
    )
    default_thresh = "baker" if job == "ts" else "gau"
    thresh = session.ask(
        "Convergence threshold (gau_loose/gau/gau_tight/gau_vtight/baker)",
        default=default_thresh,
    ).strip()
    pal_text = session.ask("Cores per ORCA call (pal)", default="1")
    mem_text = session.ask("Memory per core in MB (mem)", default="1500")
    blocks = session.ask(
        "ORCA block string, optional (e.g. '%scf maxiter 300 end'; Enter = none)",
        default="",
    ).strip()
    hessian_init = ""
    rx_modes = ""
    if job == "ts":
        hessian_init = session.ask(
            "Model Hessian (fischer/lindh/simple/swart/unit; Enter = fischer; 'calc' "
            "is refused: the ORCA-6 Hessian parse stop makes it unusable)",
            default="fischer",
        ).strip()
        rx_modes = session.ask(
            "Optional rx_modes (YAML flow, e.g. [[[[DIHEDRAL, 2, 0, 1, 3], 1]]]; "
            "Enter = none)",
            default="",
        ).strip()
    max_cycles_text = session.ask(
        "Cycle limit (Enter = the program's own default 150)", default=""
    ).strip()
    try:
        atoms = _atoms_from_structure_source(source)
        yaml_text = pysisyphus_recipe.input_yaml(
            xyz_fn=f"{source.stem}.pysisyphus.xyz",
            keywords=keywords,
            charge=int(charge_text),
            mult=int(mult_text),
            job=job,
            pal=int(pal_text),
            mem=int(mem_text),
            thresh=thresh,
            blocks=blocks or None,
            hessian_init=hessian_init or None,
            rx_modes=rx_modes or None,
            max_cycles=int(max_cycles_text) if max_cycles_text else None,
        )
    except (PysisyphusError, ParserError, OSError, ValueError) as exc:
        session.say(f"pysisyphus input generation failed: {exc}")
        return
    xyz_path = source.with_name(f"{source.stem}.pysisyphus.xyz")
    yaml_path = source.with_name(f"{source.stem}.pysisyphus.yaml")
    xyz_path.write_text(
        pysisyphus_recipe.structure_xyz(
            atoms, comment=f"generated by fBlockKit (menu 32) from {source.name}"
        ),
        encoding="utf-8",
    )
    yaml_path.write_text(yaml_text, encoding="utf-8")
    body = pysisyphus_recipe.render(
        xyz_name=xyz_path.name,
        yaml_name=yaml_path.name,
        n_atoms=len(atoms),
        job=job,
        keywords=keywords,
        charge=int(charge_text),
        mult=int(mult_text),
        pal=int(pal_text),
        mem=int(mem_text),
        thresh=thresh,
        hessian_init=hessian_init or None,
    )
    session.say(body)
    session.say(f"Written: {xyz_path}")
    session.say(f"Written: {yaml_path}")
    session.say("Next steps:")
    for line in pysisyphus_recipe.run_guidance_lines():
        session.say(f"  - {line}")
    section = ReportSection(
        title="pysisyphus input (B-layer generation)",
        body=body + "\n\n" + "\n".join(f"- {line}" for line in pysisyphus_recipe.run_guidance_lines()),
    )
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(pysisyphus_recipe.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = source.with_name(f"{source.stem}.pysisyphus.fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")


def pysisyphus_report(session: Session) -> None:
    """Menu 33: read a pysisyphus run back (cross-checks + boundary survey)."""
    path_text = session.ask(
        "pysisyphus run directory, or its console capture file (e.g. run.out)"
    )
    if not path_text:
        session.say("Cancelled (no path given).")
        return
    try:
        data = pysisyphus_analysis.read_run(path_text)
    except (PysisyphusError, pysisyphus_analysis.PysisyphusRunError, OSError) as exc:
        session.say(f"pysisyphus run read failed: {exc}")
        return
    body = pysisyphus_analysis.render(data)
    session.say(body)
    section = ReportSection(title="pysisyphus run report (B-layer reading)", body=body)
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(pysisyphus_analysis.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = data.directory / f"{data.log_path.stem}.pysisyphus.fbk.md"
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")


def _write_report(md_path: Path, title: str, body: str, evidence) -> None:
    """Write one ``*.fbk.md`` companion (title + body + references), the house shape."""
    section = ReportSection(title=title, body=body)
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(evidence)
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path.write_text(report_lines, encoding="utf-8")


def xtb_report(session: Session) -> None:
    """Menu 42: read an xTB capture back (run facts + pre-screening notes)."""
    path_text = session.ask(
        "xTB capture file path (the redirected xtb stdout/stderr, e.g. xtb.out)"
    )
    if not path_text:
        session.say("Cancelled (no path given).")
        return
    path = Path(path_text)
    try:
        run = read_xtb_run_file(path)
    except (XtbError, OSError) as exc:
        session.say(f"xTB run read failed: {exc}")
        return
    body = xtb_analysis.render(run, source=path.name)
    session.say(body)
    md_path = path.with_name(f"{path.stem}.xtb.fbk.md")
    _write_report(md_path, "xTB run report (B-layer reading)", body, xtb_analysis.evidence())
    session.say(f"Report written: {md_path}")


def crest_report(session: Session) -> None:
    """Menu 43: read a CREST conformer ensemble (sorted table + weights)."""
    path_text = session.ask(
        "CREST run directory, or its crest_conformers.xyz file"
    )
    if not path_text:
        session.say("Cancelled (no path given).")
        return
    path = Path(path_text)
    try:
        ensemble = read_crest_directory(path)
    except (CrestError, OSError) as exc:
        session.say(f"CREST ensemble read failed: {exc}")
        return
    directory = path if path.is_dir() else path.parent
    body = crest_analysis.render(ensemble, source=directory.name)
    session.say(body)
    md_path = directory / "crest_ensemble.fbk.md"
    _write_report(md_path, "CREST conformer ensemble (B-layer reading)", body, crest_analysis.evidence())
    session.say(f"Report written: {md_path}")


def mokit_generate(session: Session) -> None:
    """Menu 44: write an automr input (Gaussian-style .gjf + mokit{} block)."""
    source_text = session.ask("Structure source path (XYZ / ORCA input / ORCA output)")
    if not source_text:
        session.say("Cancelled (no structure given).")
        return
    source = Path(source_text)
    try:
        atoms = _atoms_from_structure_source(source)
    except (PysisyphusError, ParserError, OSError) as exc:
        session.say(f"Structure read failed: {exc}")
        return
    method = session.ask(
        "Method (Enter = CASSCF; add the size to pin it, e.g. CASSCF(6,6); NEVPT2 ok)",
        default="CASSCF",
    )
    basis = session.ask("Basis set", default="cc-pVDZ")
    options = session.ask(
        "mokit{} options, comma-separated (Enter = none; e.g. GVB_prog=Gaussian)",
        default="",
    )
    charge_text = session.ask("Charge", default="0")
    mult_text = session.ask("Multiplicity (2S+1)", default="1")
    mem_text = session.ask("Memory in GB (%mem)", default="4")
    nproc_text = session.ask("Cores (%nprocshared)", default="2")
    try:
        gjf_text = mokit_recipe.input_gjf(
            atoms=atoms,
            method=method,
            basis=basis,
            mokit_options=options,
            charge=int(charge_text),
            mult=int(mult_text),
            memory_gb=int(mem_text),
            nproc=int(nproc_text),
        )
    except (MokitInputError, ValueError) as exc:
        session.say(f"Input generation failed: {exc}")
        return
    gjf_path = source.with_name(f"{source.stem}_automr.gjf")
    gjf_path.write_text(gjf_text, encoding="utf-8")
    body = mokit_recipe.render(
        gjf_text=gjf_text,
        gjf_name=gjf_path.name,
        n_atoms=len(atoms),
        method=method.strip() or "CASSCF",
        basis=basis.strip() or "cc-pVDZ",
    )
    session.say(body)
    session.say(f"Written: {gjf_path}")
    md_path = source.with_name(f"{source.stem}_automr.fbk.md")
    _write_report(md_path, "MOKIT automr input (B-layer generation)", body, mokit_recipe.evidence())
    session.say(f"Report written: {md_path}")


def mokit_report(session: Session) -> None:
    """Menu 45: read a MOKIT/automr run back (stages, energy chain, termination)."""
    path_text = session.ask("automr output file path (the captured 'automr x.gjf' stdout)")
    if not path_text:
        session.say("Cancelled (no path given).")
        return
    path = Path(path_text)
    try:
        run = read_mokit_run_file(path)
    except (MokitError, OSError) as exc:
        session.say(f"automr run read failed: {exc}")
        return
    body = mokit_analysis.render(run, source=path.name)
    session.say(body)
    md_path = path.with_name(f"{path.stem}.mokit.fbk.md")
    _write_report(md_path, "MOKIT automr run report (B-layer reading)", body, mokit_analysis.evidence())
    session.say(f"Report written: {md_path}")
