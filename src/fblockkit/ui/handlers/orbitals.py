"""The orbitals handler group."""

from __future__ import annotations

import json
import numpy as np
from pathlib import Path
from typing import Any

from ...analysis import atomic_terms, avas, entropy_rdm, environment_spin, magnetic_doublets, orbital_mapping, orbital_portrait, orbital_space as orbital_space_analysis
from ...diagnosis import references_section
from ...knowledge.models import ReportSection
from ...parsers import ParserError, parse_auto
from ...parsers.fcidump import parse_fcidump
from ...parsers.mkl import parse_mkl
from ...parsers.molcas_single_aniso import (
    MolcasAnisoError,
    doublets_payload,
    parse_single_aniso,
    parse_susceptibility,
)
from ...parsers.orca_json import parse_orca_json
from ...recipe import aop_rotation, guess_transfer
from ..session import Session
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
    from ...knowledge.elements import is_f_element

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
    from ...parsers.orca_json import ANGULAR_LETTERS

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
    """Menu 16: the g_T * theta_3 criterion over a Kramers-doublet table.

    Three input routes.  The hand-written JSON table (per doublet: the three
    g values and theta3 or the g3 axis; see the user guide), the raw text
    output of an OpenMolcas SINGLE_ANISO run, or an ORCA output whose
    %casscf ANISO block ran -- the g tensors of every parsed pseudospin
    multiplet from either engine are converted into the same table, so the
    engine printouts and the manual table enter the criterion by one door.
    """
    path_text = session.ask(
        "Kramers-doublet table JSON path, an OpenMolcas SINGLE_ANISO output path, "
        "or an ORCA output with the ANISO block (per doublet: the three g values "
        "and theta3 or the g3 axis; see the user guide)"
    )
    if not path_text:
        session.say("Cancelled (no path given).")
        return
    path = Path(path_text)
    chi_line = ""
    try:
        text = path.read_text(encoding="utf-8")
        if text.lstrip().startswith("{"):
            payload = json.loads(text)
        elif "CALCULATION OF PSEUDOSPIN HAMILTONIAN TENSORS" in text:
            payload = doublets_payload(parse_single_aniso(text))
            try:
                susceptibility = parse_susceptibility(text)
                last = susceptibility.points[-1]
                chi_line = (
                    f"\nThe output's own susceptibility table: "
                    f"{susceptibility.n_points} temperature points; "
                    f"chi*T(300 K) = {last.chi_t_cm3k_mol:.4f} cm3*K/mol "
                    "(the free-ion 6H15/2 Curie value is about 14.2).\n"
                )
            except MolcasAnisoError:
                chi_line = ""
        else:
            from ...parsers import single_aniso as orca_single_aniso

            parsed = None
            try:
                parsed = parse_auto(path)
            except (ParserError, OSError):
                parsed = None
            segments = (
                ((parsed.sections.get("single_aniso") or {}).get("segments") or [])
                if parsed is not None
                else []
            )
            if not segments:
                session.say(
                    "The file is neither a JSON doublet table, an OpenMolcas "
                    "SINGLE_ANISO output, nor an ORCA output with a SINGLE_ANISO "
                    "section. Next step: give a JSON object with a 'doublets' "
                    "list, an OpenMolcas SINGLE_ANISO output, or an ORCA output "
                    "whose %casscf ANISO block ran (MLTP must be given)."
                )
                return
            payload = orca_single_aniso.doublets_payload(segments, system=path.name)
        table = magnetic_doublets.parse_doublets(payload)
        section = magnetic_doublets.run(table)
    except (OSError, ValueError) as exc:
        session.say(
            f"Magnetic-doublet reading failed: {exc}"
        )
        return
    session.say(section.body + chi_line)
    report_lines = f"## {section.title}\n\n{section.body}\n{chi_line}"
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


# --- 47 AOP rotation guess ----------------------------------------------------


def aop_rotation_guess(session: Session) -> None:
    """Menu 47: rotate a target's orbitals onto a reference active space (two-step SVD)."""
    path_text = session.ask(
        "Rotation manifest JSON path (a reference active space and the target's "
        "export + mkl; see the user guide)"
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
        reference_entry = payload["reference"]
        reference = parse_orca_json(_resolve(reference_entry["export"]))
        active = list(reference_entry["active"])
        target_entry = payload["target"]
        target = parse_orca_json(_resolve(target_entry["export"]))
        template_mkl_path = _resolve(target_entry["mkl"])
        template_mkl = parse_mkl(template_mkl_path)
        result = aop_rotation.rotate_guess(
            reference, target, active, n_closed=int(payload["closed"])
        )
    except (OSError, KeyError, TypeError, ValueError, ParserError) as exc:
        session.say(f"AOP rotation failed: {exc}")
        return
    body = aop_rotation.render(result)
    session.say(body)
    guess_path = template_mkl_path.with_name(template_mkl_path.stem + ".aop.fbk.mkl")
    try:
        aop_rotation.write_guess_mkl(template_mkl, result, guess_path)
    except (OSError, aop_rotation.AopError) as exc:
        session.say(f"Writing the guess mkl failed: {exc}")
        return
    session.say(f"Guess written: {guess_path}")
    session.say(
        f"Next step: run ``orca_2mkl {guess_path.stem} -gbw`` next to that file, then "
        f"start the target calculation with ``!moread`` and %moinp \"{guess_path.stem}.gbw\"."
    )
    section = ReportSection(
        title="G7 AOP rotation guess (a target orbital set aligned onto a reference active space)",
        body=body,
    )
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(aop_rotation.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = path.with_name(path.name + ".aop.fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")


def perturb_batch(session: Session) -> None:
    """Menu 29: perturb a converged reference's orbitals into K restart inputs."""
    from ...recipe import perturb_guess

    mkl_text = session.ask(
        "Reference mkl path (``orca_2mkl <base> -mkl`` of the converged run)"
    )
    if not mkl_text:
        session.say("Cancelled (no reference mkl given).")
        return
    input_text = session.ask("Base ORCA input path (the input that produced the reference run)")
    if not input_text:
        session.say("Cancelled (no base input given).")
        return
    starts_text = session.ask("Number of perturbed starts (Enter = 3)").strip()
    seed_text = session.ask(
        "Random seed; it fixes every pair and angle, so the batch replays "
        "byte-identically (Enter = 20260927)"
    ).strip()
    pairs_text = session.ask("Pairs per start (Enter = 10, the source's value)").strip()
    window_text = session.ask("Orbital window per side (Enter = 15, the source's value)").strip()
    mkl_path = Path(mkl_text)
    base_path = Path(input_text)
    try:
        starts = int(starts_text or "3")
        seed = int(seed_text or "20260927")
        n_pairs = int(pairs_text or "10")
        window = int(window_text or "15")
        if starts < 1:
            raise perturb_guess.PerturbError("the number of starts must be at least one.")
        reference = parse_mkl(mkl_path)
        base_text = base_path.read_text(encoding="utf-8", errors="replace")
        geometry_note = perturb_guess.check_geometry_match(base_text, reference)
        written: list[tuple[Path, Path, str]] = []
        records: list[str] = []
        for index in range(1, starts + 1):
            start = perturb_guess.perturb_mkl(
                reference, seed=seed + index - 1, n_pairs=n_pairs, window=window
            )
            stem = f"{mkl_path.stem}.p{index}.fbk"
            mkl_out = mkl_path.with_name(stem + ".mkl")
            mkl_out.write_text(start.mkl.render(), encoding="utf-8")
            variant = perturb_guess.mo_read_variant(base_text, stem + ".gbw")
            inp_out = base_path.with_name(f"{base_path.stem}.p{index}.inp")
            inp_out.write_text(variant, encoding="utf-8")
            records.append(perturb_guess.render(start, seed=seed + index - 1, index=index))
            written.append((mkl_out, inp_out, stem))
    except (perturb_guess.PerturbError, ParserError, OSError, ValueError) as exc:
        session.say(f"Perturbation batch failed: {exc}")
        return
    session.say(
        f"Reference: {mkl_path.name} ({reference.n_mo} orbitals, "
        f"{'unrestricted' if reference.unrestricted else 'restricted'}); "
        f"{starts} perturbed start(s), seed {seed}."
    )
    if geometry_note:
        session.say(f"Note: {geometry_note}")
    for record in records:
        session.say(record)
    for mkl_out, inp_out, stem in written:
        session.say(f"Written: {mkl_out}")
        session.say(f"Written: {inp_out}")
    session.say("Next steps:")
    session.say(
        "  - convert each orbital file next to itself: "
        f"``orca_2mkl {written[0][2]} -gbw`` (the inputs already point at "
        f"{written[0][2]}.gbw);"
    )
    session.say(
        "  - run the inputs and compare the converged energies with the reference's: "
        "a start finding a LOWER energy heals a wrongly converged reference; the "
        "same (or a higher) energy is no information -- the source is explicit that "
        "the test cannot guarantee detection"
    )
    session.say(
        "  - the perturbed columns are mixtures, not eigenfunctions: the occupations "
        "and orbital energies in the file are the reference's, and ORCA "
        "re-determines everything after reading the guess"
    )
    session.say(
        "  - complementary check: a stability analysis detects unstable solutions "
        "(saddle points) but cannot distinguish local from global minima (the "
        "source's own boundary)"
    )
    body = "\n".join(records) + "\n\nNext steps:\n" + "\n".join(
        "  - " + line
        for line in (
            f"convert with ``orca_2mkl <name> -gbw`` and run: {', '.join(w[2] for w in written)}",
            "a lower converged energy than the reference's identifies wrong convergence",
            "the perturbation cannot guarantee detection (the source's boundary)",
        )
    )
    section = ReportSection(
        title="perturbed multistart batch (randomized occupied-virtual mixing)",
        body=body,
    )
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(perturb_guess.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = base_path.with_name(base_path.name + ".perturb.fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")




def state_data_report(session: Session) -> None:
    """Menu 31: per-state CASSCF data from the property file."""
    from ...analysis import state_data as state_data_analysis
    from ...parsers.orca_property import parse_property

    path_text = session.ask("ORCA property file path (``<base>.property.txt``)")
    if not path_text:
        session.say("Cancelled (no path given).")
        return
    path = Path(path_text)
    try:
        sections = parse_property(path)
        table = state_data_analysis.state_table(sections)
        lines = state_data_analysis.transitions(sections)
    except (ParserError, state_data_analysis.StateDataError, OSError) as exc:
        session.say(f"State data failed: {exc}")
        return
    body = state_data_analysis.render(table, lines, source=path.name)
    session.say(body)
    section = ReportSection(title="CASSCF state data (property file)", body=body)
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(state_data_analysis.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = path.with_name(path.name + ".states.fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")


def judd_ofelt_fit(session: Session) -> None:
    """Menu 34: the standard Judd-Ofelt fit from a transition dataset."""
    from ...analysis import judd_ofelt as judd_ofelt_analysis

    path_text = session.ask("Judd-Ofelt dataset path (YAML; see the formats chapter)")
    if not path_text:
        session.say("Cancelled (no path given).")
        return
    weighting = session.ask(
        "Least squares: unweighted (Enter) or the normalized 1/S variant (type 1/S)"
    ).strip()
    weights = "1/S" if weighting == "1/S" else "none"
    path = Path(path_text)
    try:
        dataset = judd_ofelt_analysis.read_dataset(path)
        result = judd_ofelt_analysis.fit(dataset, weights=weights)
        emission = judd_ofelt_analysis.emission_rates(result)
    except (judd_ofelt_analysis.JudOError, OSError) as exc:
        session.say(f"Judd-Ofelt fit failed: {exc}")
        return
    body = judd_ofelt_analysis.render(result, source=path.name)
    if emission:
        body += "\n" + judd_ofelt_analysis.render_emission(emission)
    session.say(body)
    section = ReportSection(title="Judd-Ofelt intensity parameters", body=body)
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(judd_ofelt_analysis.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = path.with_name(path.name + ".jo.fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")
    csv_path = path.with_name(path.name + ".jo.fbk.csv")
    csv_path.write_text(judd_ofelt_analysis.jo_plot_csv(result), encoding="utf-8")
    session.say(
        f"Plot data written: {csv_path} (transition table; columns "
        "side,label,energy_cm1,f_exp,s_exp,s_ed,a_ed_s1,branching)"
    )


def ailft_report(session: Session) -> None:
    """Menu 38: ab initio ligand-field analysis from an ORCA AILFT output."""
    from ...analysis import ailft as ailft_analysis

    path_text = session.ask("ORCA output path (a CASSCF run with the AILFT driver)")
    if not path_text:
        session.say("Cancelled (no path given).")
        return
    path = Path(path_text)
    try:
        parsed = parse_auto(path)
    except (ParserError, OSError) as exc:
        session.say(f"Ligand-field analysis failed: {exc}")
        return
    data = parsed.sections.get("ailft") or {}
    ref_text = session.ask(
        "Free-ion reference from the built-in table: element+charge "
        "(e.g. Nd3+; Enter = skip/manual)"
    ).strip()
    reference: tuple[str, dict] | None = None
    if ref_text:
        from ...knowledge import ailft_references as free_ion_table

        entry = free_ion_table.lookup(ref_text)
        if entry is None:
            session.say(
                f"{ref_text!r} is not in the built-in table. It covers the "
                "trivalent lanthanide and actinide series plus the measured "
                "Ni(2+) probe (the manual lists the ions); enter the "
                "references manually below, or skip."
            )
        else:
            reference = (ref_text, entry)
    b0_text = session.ask("Free-ion Racah B (cm-1) for the nephelauxetic ratio (Enter = skip)").strip()
    zeta_text = session.ask("Free-ion SOC constant zeta0 (cm-1) (Enter = skip)").strip()
    def _opt(text: str) -> float | None:
        if not text:
            return None
        try:
            return float(text)
        except ValueError:
            session.say(f"Not a number: {text!r}; skipping this reference.")
            return None
    try:
        casscf_status = parsed.sections.get("casscf") or {}
        convergence: tuple[bool | None, str | None] = (
            (casscf_status.get("converged"), casscf_status.get("converged_via"))
            if casscf_status.get("present")
            else (None, None)
        )
        body = ailft_analysis.render(
            data,
            source=path.name,
            free_ion_B_cm1=_opt(b0_text),
            free_ion_zeta_cm1=_opt(zeta_text),
            convergence=convergence,
            free_ion_reference=reference,
        )
    except ailft_analysis.AilftError as exc:
        session.say(f"Ligand-field analysis failed: {exc}")
        return
    session.say(body)
    section = ReportSection(title="Ab initio ligand-field analysis", body=body)
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(ailft_analysis.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = path.with_name(path.name + ".ailft.fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")


def xas_report(session: Session) -> None:
    """Menu 37: core-excited spectra -- read a ROCIS output or write inputs.

    Mode 1 reads the ROCIS absorption blocks, reports the transition table
    and the edge branching ratio, and states the RIXS bookkeeping of the
    run; mode 2 writes a ROCIS input that carries the off-resonance XES
    request (the plain RIXS channel is its carrier); mode 3 writes the
    two-step CAS-CI core-excited XAS inputs; mode 4 writes their RAS-CI
    XES counterpart (the emission side of the same two-step shape).
    """
    from ...analysis import xas as xas_analysis

    mode = session.ask(
        "What do you need? (1) read a ROCIS output (Enter), (2) write a ROCIS "
        "XES input, (3) write a CAS-CI core-excited XAS input, (4) write a "
        "RAS-CI core-excited XES input"
    ).strip()
    if mode == "2":
        _rocis_xes_write(session)
        return
    if mode == "3":
        _casci_xas_write(session)
        return
    if mode == "4":
        _casci_xes_write(session)
        return
    path_text = session.ask("ORCA output path (a ROCIS core-excited-spectra run)")
    if not path_text:
        session.say("Cancelled (no path given).")
        return
    path = Path(path_text)
    try:
        parsed = parse_auto(path)
    except (ParserError, OSError) as exc:
        session.say(f"Core-excited spectra report failed: {exc}")
        return
    data = parsed.sections.get("rocis") or {}
    stat_text = session.ask(
        "Statistical branching ratio for the comparison (Enter = report the ratio only)"
    ).strip()
    r_stat = None
    if stat_text:
        try:
            r_stat = float(stat_text)
        except ValueError:
            session.say(f"Not a number: {stat_text!r}; reporting the ratio only.")
    try:
        body = xas_analysis.render(data, source=path.name, r_stat=r_stat)
    except xas_analysis.XasError as exc:
        session.say(f"Core-excited spectra report failed: {exc}")
        return
    session.say(body)
    section = ReportSection(title="Core-excited spectra (XAS/RIXS)", body=body)
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(xas_analysis.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = path.with_name(path.name + ".xas.fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")
    csv_path = path.with_name(path.name + ".xas.fbk.csv")
    csv_path.write_text(xas_analysis.xas_plot_csv(data), encoding="utf-8")
    session.say(
        f"Plot data written: {csv_path} (primary-block transitions; columns "
        "i_root,i_label,j_root,j_label,energy_eV,fosc)"
    )


def _rocis_xes_write(session: Session) -> None:
    """Menu 37, mode 2: the ROCIS XES input writer (off-resonance emission)."""
    from ...analysis import geometry as geometry_analysis
    from ...recipe import rocis_xes as xes_recipe

    xyz_text = session.ask("Structure file (XYZ) path")
    if not xyz_text:
        session.say("Cancelled (no structure given).")
        return
    charge_text = session.ask("Charge of the system (Enter = 0)").strip()
    multiplicity_text = session.ask(
        "Multiplicity (Enter = 1; it also becomes the ROCIS ReferenceMult)"
    ).strip()
    rohf_text = session.ask(
        "Unpaired electrons for the ROHF high-spin preparation (Enter = no "
        "%scf block)"
    ).strip()
    element_text = session.ask(
        "XASelems: the 0-based position of the core element among the atoms "
        "(Enter = 0)"
    ).strip()
    nroots_text = session.ask(
        "NRoots (Enter = 30; enough roots are needed for the plain RIXS "
        "channel, the carrier of the XES table -- measured: 10 skipped it, "
        "30 covered it)"
    ).strip()
    window_text = session.ask(
        "OrbWin: six integers 'd1s,d1e,d2s,d2e,accs,acce' -- the two "
        "spin-orbit-split core ranges, then a wide acceptor (the probe's own "
        "numbering: '6,6,7,8,0,2000')"
    )
    rixssoc_text = session.ask(
        "Include the SOC-corrected RIXS channel? (y/N; it stores large "
        "transition-density files -- measured past 36 GB at NRoots 30 on the "
        "probe)"
    ).strip().lower()
    elastic_text = session.ask("Include the elastic line? (Enter = yes)").strip().lower()
    keywords_text = session.ask(
        "Method/basis keywords (Enter = x2c x2c-SVPall AutoAux TightSCF)"
    ).strip()
    nprocs_text = session.ask("Parallel processes (Enter = 8)").strip()
    maxcore_text = session.ask("MaxCore in MB (Enter = 4000)").strip()
    path = Path(xyz_text)
    try:
        atoms = geometry_analysis.parse_xyz(path)
        coordinates = tuple((atom.element, atom.x, atom.y, atom.z) for atom in atoms)
        window = tuple(int(token) for token in window_text.replace(",", " ").split())
        plan = xes_recipe.build_input(
            coordinates,
            charge=int(charge_text or "0"),
            multiplicity=int(multiplicity_text or "1"),
            xas_element=int(element_text or "0"),
            nroots=int(nroots_text or "30"),
            window=window,
            do_rixssoc=rixssoc_text in ("y", "yes"),
            do_elastic=elastic_text not in ("n", "no"),
            rohf_electrons=int(rohf_text) if rohf_text else None,
            keywords=keywords_text or xes_recipe.DEFAULT_KEYWORDS,
            nprocs=int(nprocs_text or "8"),
            maxcore=int(maxcore_text or "4000"),
        )
    except (
        xes_recipe.RocisXesError,
        geometry_analysis.StructureError,
        OSError,
        ValueError,
    ) as exc:
        session.say(f"ROCIS XES input failed: {exc}")
        return
    inp_path = path.with_name(f"{path.stem}.xes.inp")
    inp_path.write_text(plan.text, encoding="utf-8")
    body = xes_recipe.render_plan(plan, path=inp_path)
    session.say(body)
    section = ReportSection(
        title="ROCIS XES input (off-resonance X-ray emission)", body=body
    )
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(xes_recipe.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = inp_path.with_name(inp_path.name + ".fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")


def _casci_xas_write(session: Session) -> None:
    """Menu 37, mode 3: the two-step CAS-CI core-excited XAS input writer."""
    from ...analysis import geometry as geometry_analysis
    from ...recipe import casci_xas as casci_recipe

    xyz_text = session.ask("Structure file (XYZ) path")
    if not xyz_text:
        session.say("Cancelled (no structure given).")
        return
    charge_text = session.ask("Charge of the system (Enter = 0)").strip()
    multiplicity_text = session.ask("Multiplicity (Enter = 1)").strip()
    space_text = session.ask(
        "Valence active space of step 1 as 'nel norb' (e.g. '6 5' for a d6 shell)"
    )
    nroots_text = session.ask("Step-1 roots (Enter = 5)").strip()
    cores_text = session.ask(
        "Core orbital indices to rotate in (0-based, from the step-1 output's "
        "orbital table; e.g. '6 7 8' for the Fe 2p near -700 eV)"
    )
    step2_mult_text = session.ask(
        "Step-2 multiplicities (Enter = the ground multiplicity and ground-2)"
    ).strip()
    step2_nroots_text = session.ask("Step-2 roots per multiplicity (Enter = 20,20)").strip()
    keywords_text = session.ask(
        "Method/basis keywords (Enter = def2-SVP def2-SVP/C TightSCF)"
    ).strip()
    nprocs_text = session.ask("Parallel processes (Enter = 8)").strip()
    maxcore_text = session.ask("MaxCore in MB (Enter = 4000)").strip()
    path = Path(xyz_text)
    try:
        atoms = geometry_analysis.parse_xyz(path)
        coordinates = tuple((atom.element, atom.x, atom.y, atom.z) for atom in atoms)
        space_tokens = space_text.replace(",", " ").split()
        plan = casci_recipe.build_inputs(
            coordinates,
            charge=int(charge_text or "0"),
            multiplicity=int(multiplicity_text or "1"),
            valence_nel=int(space_tokens[0]),
            valence_norb=int(space_tokens[1]),
            step1_nroots=int(nroots_text or "5"),
            core_orbitals=[int(token) for token in cores_text.replace(",", " ").split()],
            step2_mult=step2_mult_text or None,
            step2_nroots=step2_nroots_text or None,
            keywords=keywords_text or casci_recipe.DEFAULT_KEYWORDS,
            maxcore=int(maxcore_text or "4000"),
            nprocs=int(nprocs_text or "8"),
            gbw_name=f"{path.stem}.casci_xas.step1.gbw",
        )
    except (
        casci_recipe.CasciXasError,
        geometry_analysis.StructureError,
        OSError,
        ValueError,
        IndexError,
    ) as exc:
        session.say(f"CAS-CI XAS input failed: {exc}")
        return
    step1_path = path.with_name(f"{path.stem}.casci_xas.step1.inp")
    step2_path = path.with_name(f"{path.stem}.casci_xas.step2.inp")
    step1_path.write_text(plan.step1_text, encoding="utf-8")
    step2_path.write_text(plan.step2_text, encoding="utf-8")
    body = casci_recipe.render_plan(plan, step1_path=step1_path, step2_path=step2_path)
    session.say(body)
    section = ReportSection(
        title="CAS-CI core-excited XAS input (two-step protocol)", body=body
    )
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(casci_recipe.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = step1_path.with_name(f"{path.stem}.casci_xas.fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")


def _casci_xes_write(session: Session) -> None:
    """Menu 37, mode 4: the two-step RAS-CI core-excited XES input writer.

    The emission side of the CAS-CI/RAS-CI two-step protocol (manual
    section 3.13.19): the same valence SA-CASSCF, then the core orbitals
    (e.g. a metal 1s and 3p for K-beta emission) rotated into the window
    and one CAS-CI iteration over the singly-core-excited space (the
    ``refs ras(Nel: NRAS1 1 / NRAS2 / 0 0)`` restriction), with the
    ``XESSOC``/``XASMOs`` rel request.
    """
    from ...analysis import geometry as geometry_analysis
    from ...recipe import casci_xas as casci_recipe

    xyz_text = session.ask("Structure file (XYZ) path")
    if not xyz_text:
        session.say("Cancelled (no structure given).")
        return
    charge_text = session.ask("Charge of the system (Enter = 0)").strip()
    multiplicity_text = session.ask("Multiplicity (Enter = 1)").strip()
    space_text = session.ask(
        "Valence active space of step 1 as 'nel norb' (e.g. '6 5' for a d6 shell)"
    )
    nroots_text = session.ask("Step-1 roots (Enter = 5)").strip()
    cores_text = session.ask(
        "Core orbital indices to rotate in (0-based, from the step-1 output's "
        "orbital table; for a K-beta emission the metal 1s and 3p, e.g. "
        "'0 26 27 28' for Fe)"
    )
    xasmo_text = session.ask(
        "XASMOs -- the global index of the rotated 1s MO (Enter = the window "
        "head: the lowest-index core moves to the leading slot)"
    ).strip()
    step2_mult_text = session.ask(
        "Step-2 multiplicities (Enter = the ground multiplicity and ground-2)"
    ).strip()
    step2_nroots_text = session.ask(
        "Step-2 roots per multiplicity (Enter = 1000,1000: the manual's large "
        "number; the engine adjusts it to the CSF count of the restricted space)"
    ).strip()
    keywords_text = session.ask(
        "Method/basis keywords (Enter = def2-SVP def2-SVP/C TightSCF)"
    ).strip()
    nprocs_text = session.ask("Parallel processes (Enter = 8)").strip()
    maxcore_text = session.ask("MaxCore in MB (Enter = 4000)").strip()
    path = Path(xyz_text)
    try:
        atoms = geometry_analysis.parse_xyz(path)
        coordinates = tuple((atom.element, atom.x, atom.y, atom.z) for atom in atoms)
        space_tokens = space_text.replace(",", " ").split()
        plan = casci_recipe.build_inputs(
            coordinates,
            charge=int(charge_text or "0"),
            multiplicity=int(multiplicity_text or "1"),
            valence_nel=int(space_tokens[0]),
            valence_norb=int(space_tokens[1]),
            step1_nroots=int(nroots_text or "5"),
            core_orbitals=[int(token) for token in cores_text.replace(",", " ").split()],
            step2_mult=step2_mult_text or None,
            step2_nroots=step2_nroots_text or None,
            keywords=keywords_text or casci_recipe.DEFAULT_KEYWORDS,
            maxcore=int(maxcore_text or "4000"),
            mode="xes",
            xas_mo=int(xasmo_text) if xasmo_text else None,
            nprocs=int(nprocs_text or "8"),
            gbw_name=f"{path.stem}.casci_xes.step1.gbw",
        )
    except (
        casci_recipe.CasciXasError,
        geometry_analysis.StructureError,
        OSError,
        ValueError,
        IndexError,
    ) as exc:
        session.say(f"RAS-CI XES input failed: {exc}")
        return
    step1_path = path.with_name(f"{path.stem}.casci_xes.step1.inp")
    step2_path = path.with_name(f"{path.stem}.casci_xes.step2.inp")
    step1_path.write_text(plan.step1_text, encoding="utf-8")
    step2_path.write_text(plan.step2_text, encoding="utf-8")
    body = casci_recipe.render_plan(plan, step1_path=step1_path, step2_path=step2_path)
    session.say(body)
    section = ReportSection(
        title="RAS-CI core-excited XES input (two-step protocol)", body=body
    )
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(casci_recipe.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = step1_path.with_name(f"{path.stem}.casci_xes.fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")


def relaxation_report(session: Session) -> None:
    """Menu 36: magnetic relaxation / QTM metrics from an ORCA output.

    Reads either the SINGLE_ANISO embedded section (the ANISO sub-block of
    %casscf) or the Orca_Magrelax section; a file with both gets both parts.
    """
    from ...analysis import relaxation as relaxation_analysis

    path_text = session.ask(
        "ORCA output path (with a SINGLE_ANISO (ANISO sub-block) or MAGRELAX section)"
    )
    if not path_text:
        session.say("Cancelled (no path given).")
        return
    path = Path(path_text)
    try:
        parsed = parse_auto(path)
    except (ParserError, OSError) as exc:
        session.say(f"Magnetic relaxation report failed: {exc}")
        return
    blocks: list[str] = []
    single_aniso = parsed.sections.get("single_aniso") or {}
    segments = single_aniso.get("segments") or []
    for index, segment in enumerate(segments):
        blocks.append(
            relaxation_analysis.render_single_aniso(
                segment, segment_index=index, n_segments=len(segments), source=path.name
            )
        )
    magrelax_missing = ""
    try:
        magrelax_data = relaxation_analysis.read_magrelax(path)
    except relaxation_analysis.RelaxationError as exc:
        magrelax_data = None
        magrelax_missing = str(exc)
    if magrelax_data is not None:
        blocks.append(relaxation_analysis.render_magrelax(magrelax_data, path.name))
    if not blocks:
        session.say(f"No magnetic relaxation data found. {magrelax_missing}")
        return
    body = "\n\n".join(blocks)
    session.say(body)
    section = ReportSection(title="Magnetic relaxation and QTM", body=body)
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(relaxation_analysis.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = path.with_name(path.name + ".relax.fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")
    if magrelax_data is not None and magrelax_data.tau_s:
        csv_path = path.with_name(path.name + ".relax.fbk.csv")
        csv_path.write_text(
            relaxation_analysis.magrelax_plot_csv(magrelax_data), encoding="utf-8"
        )
        session.say(
            f"Plot data written: {csv_path} (tau(T) table; "
            "columns temperature_K,rate_per_s,tau_s)"
        )


def tunnelling_report(session: Session) -> None:
    """Menu 46: tunnelling-relaxation prediction from a SINGLE_ANISO output.

    The equivalent-Zeeman model (Yin & Li 2020) runs on the per-doublet g
    values and energies alone; the spin-dipolar model (Aravena 2018/2026,
    dilution variant Llanos & Aravena 2019) takes an optional neighbour
    table (dx dy dz mx my mz per line, in the central g-frame).
    """
    from ...analysis import relaxation as relaxation_analysis
    from ...analysis import tunnelling as tunnelling_analysis
    from ...analysis.tunnelling import TunnellingError

    path_text = session.ask("ORCA output path (with a SINGLE_ANISO section)")
    if not path_text:
        session.say("Cancelled (no path given).")
        return
    path = Path(path_text)
    try:
        parsed = parse_auto(path)
    except (ParserError, OSError) as exc:
        session.say(f"Tunnelling prediction failed: {exc}")
        return
    segments = (parsed.sections.get("single_aniso") or {}).get("segments") or []
    if not segments:
        session.say(
            "No SINGLE_ANISO section found in this output. Next step: give an "
            "output of a %casscf run with the ANISO sub-block (menu 36's input)."
        )
        return
    b_text = session.ask(
        "B_ave in mT (the empirical internal-field scale; Enter = 20.0)",
        default="20.0",
    )
    try:
        b_ave = float(b_text)
    except ValueError:
        session.say("B_ave must be a number (mT).")
        return
    table_text = session.ask(
        "Neighbour table path for the spin-dipolar model (Enter = skip)"
    )
    neighbour_report = None
    try:
        segment = segments[-1]
        rows = relaxation_analysis.group_metrics(segment)
        levels = tunnelling_analysis.kd_levels(rows)
        ueff_rows = tunnelling_analysis.ueff_curve(levels, B_ave_mT=b_ave)
        if table_text:
            table_path = Path(table_text)
            neighbours = tunnelling_analysis.parse_neighbour_table(
                table_path.read_text(encoding="utf-8", errors="replace")
            )
            neighbour_report = {
                "table": table_path.name,
                "n_neighbours": len(neighbours),
                **tunnelling_analysis.dipolar_tau(neighbours, levels[0].g),
            }
            dilution_text = session.ask(
                "Dilution concentrations, comma-separated (Enter = skip)"
            )
            if dilution_text.strip():
                try:
                    concentrations = tuple(
                        float(part) for part in dilution_text.split(",") if part.strip()
                    )
                except ValueError:
                    session.say("The dilution concentrations must be numbers; skipping them.")
                    concentrations = ()
                if concentrations:
                    neighbour_report["dilution"] = tunnelling_analysis.dilution_medians(
                        neighbours, levels[0].g, concentrations
                    )
    except (TunnellingError, OSError) as exc:
        session.say(f"Tunnelling prediction failed: {exc}")
        return
    body = tunnelling_analysis.render(
        source=f"{path.name} (segment {len(segments)} of {len(segments)})",
        levels=levels,
        B_ave_mT=b_ave,
        ueff_rows=ueff_rows,
        neighbour_report=neighbour_report,
    )
    session.say(body)
    section = ReportSection(title="Quantum-tunnelling relaxation prediction", body=body)
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(tunnelling_analysis.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = path.with_name(path.name + ".qtm.fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")


def pnmr_report(session: Session) -> None:
    """Menu 35: pseudocontact shifts from a susceptibility tensor + a structure."""
    from ...analysis import pnmr as pnmr_analysis

    path_text = session.ask("pNMR run file path (YAML; see the formats chapter)")
    if not path_text:
        session.say("Cancelled (no path given).")
        return
    path = Path(path_text)
    try:
        data = pnmr_analysis.read_run(path)
    except (pnmr_analysis.PnmrError, OSError) as exc:
        session.say(f"pNMR run failed: {exc}")
        return
    body = pnmr_analysis.render(data)
    ion_text = session.ask(
        "Bleaney comparator: the Ln(III) ion (Tb, Dy, Ho, Er, Tm, Yb; Enter = skip)"
    ).strip()
    if ion_text:
        b02_text = session.ask(
            "B_0^2 in cm^-1 (the axial crystal-field parameter; Enter = skip)"
        ).strip()
        if b02_text:
            temperature_text = session.ask(
                "Temperature in K for the comparator (Enter = 300)"
            ).strip()
            from ...analysis import bleaney as bleaney_analysis

            try:
                b02_value = float(b02_text)
                temperature_value = float(temperature_text or "300")
            except ValueError:
                session.say(
                    "Bleaney comparator skipped: the B_0^2 and temperature "
                    "entries must be numbers (B_0^2 in cm^-1). Next step: give "
                    "numeric values, or leave the ion question empty to skip "
                    "the comparator."
                )
            else:
                try:
                    block = bleaney_analysis.comparator_lines(
                        ion_text, b02_value, temperature_value
                    )
                except bleaney_analysis.BleaneyError as exc:
                    session.say(f"Bleaney comparator skipped: {exc}")
                else:
                    body = body + "\n" + "\n".join(block)
    session.say(body)
    section = ReportSection(title="pNMR pseudocontact shifts", body=body)
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(pnmr_analysis.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = path.with_name(path.name + ".pnmr.fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")


def poly_aniso_report(session: Session) -> None:
    """Menu 39: polynuclear magnetism -- read a POLY_ANISO output or write the input.

    Mode 1 reads the cluster report of the ORCA POLY_ANISO driver (the
    per-center single-ion data, the exchange decomposition, the coupled
    states, the chiT(T) table and the Van Vleck susceptibility tensors);
    mode 2 writes the driver's input plus the cluster checklist.
    """
    from ...analysis import poly_aniso as poly_analysis

    mode = session.ask(
        "What do you need? (1) read a poly_aniso.output report, (2) write a "
        "POLY_ANISO input + checklist (Enter = 1)"
    ).strip()
    if mode == "2":
        _poly_aniso_write(session)
        return
    path_text = session.ask("poly_aniso.output path (from otool_poly_aniso)")
    if not path_text:
        session.say("Cancelled (no path given).")
        return
    path = Path(path_text)
    try:
        parsed = parse_auto(path)
    except (ParserError, OSError) as exc:
        session.say(f"Polynuclear magnetism report failed: {exc}")
        return
    data = parsed.sections.get("poly_aniso") or {}
    try:
        body = poly_analysis.render(data, source=path.name)
    except poly_analysis.PolyAnisoError as exc:
        session.say(f"Polynuclear magnetism report failed: {exc}")
        return
    session.say(body)
    section = ReportSection(title="Polynuclear magnetism (POLY_ANISO)", body=body)
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(poly_analysis.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = path.with_name(path.name + ".polyaniso.fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")


def _poly_aniso_write(session: Session) -> None:
    """Menu 39, mode 2: the POLY_ANISO input writer (cluster description -> input)."""
    from ...recipe import poly_aniso as poly_recipe

    path_text = session.ask(
        "Output path for the generated input (Enter = ./poly_aniso.input)"
    )
    path = Path(path_text) if path_text else Path("poly_aniso.input")

    def _integers(text: str, count: int | None, what: str):
        try:
            values = [int(token) for token in text.replace(",", " ").split()]
        except ValueError:
            session.say(f"{what} must be integers, space-separated. Cancelled.")
            return None
        if not values or (count is not None and len(values) != count):
            session.say(
                f"{what} needs {count if count else 'at least one'} integer(s). "
                "Cancelled."
            )
            return None
        return values

    types = _integers(
        session.ask("Number of non-equivalent centre types (1-6)"), 1, "The type count"
    )
    if types is None:
        return
    n_types = types[0]
    centres = _integers(
        session.ask(f"Equivalent centres per type ({n_types} integers)"),
        n_types,
        "The centre-count line",
    )
    if centres is None:
        return
    states = _integers(
        session.ask(
            f"Low-lying spin-orbit functions per centre type ({n_types} integers)"
        ),
        n_types,
        "The spin-orbit-count line",
    )
    if states is None:
        return

    def _vector(text: str):
        try:
            values = [float(token) for token in text.replace(",", " ").split()]
        except ValueError:
            return None
        return values if len(values) == 3 else None

    # the SYMM block is mandatory whenever a type carries more than one
    # equivalent centre (the driver's own measured check; it aborts with a
    # serious-error banner yet still exits 0, so the writer refuses instead)
    _identity = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
    symmetry = None
    if any(count > 1 for count in centres):
        symmetry = []
        for type_index, count in enumerate(centres, start=1):
            if count == 1:
                symmetry.append((_identity,))
                continue
            matrices = []
            for matrix_index in range(1, count + 1):
                first_row = session.ask(
                    f"Type {type_index}, rotation matrix {matrix_index} of {count}: "
                    "row 1 'r11 r12 r13' (Enter = the identity matrix)"
                )
                if not first_row:
                    matrices.append(_identity)
                    continue
                rows = [_vector(first_row)]
                for row_index in (2, 3):
                    rows.append(
                        _vector(
                            session.ask(
                                f"Type {type_index}, matrix {matrix_index}, row "
                                f"{row_index} 'r{row_index}1 r{row_index}2 r{row_index}3'"
                            )
                        )
                    )
                if any(row is None for row in rows):
                    session.say("A matrix row must be three numbers. Cancelled.")
                    return
                matrices.append(tuple(value for row in rows for value in row))
            symmetry.append(tuple(matrices))

    coordinates = None
    first = session.ask(
        "Coordinates of type 1, 'x y z' in Angstrom (Enter = skip the COOR block)"
    )
    if first:
        row = _vector(first)
        if row is None:
            session.say("The coordinates must be three numbers. Cancelled.")
            return
        coordinates = [row]
        for index in range(2, n_types + 1):
            row = _vector(session.ask(f"Coordinates of type {index}, 'x y z'"))
            if row is None:
                session.say("The coordinates must be three numbers. Cancelled.")
                return
            coordinates.append(row)

    pair_model = session.ask(
        "Pair model: Enter = Lines isotropic 'i j J' / type lin3 = axis-diagonal "
        "'i j Jx Jy Jz'"
    ).strip().lower() or "lines"
    if pair_model not in ("lines", "lin3", "lin9"):
        session.say(
            f"Unknown pair model {pair_model!r}; the writer renders 'lines' or "
            "'lin3'. Cancelled."
        )
        return
    pairs = []
    while True:
        label = "'i j Jx Jy Jz'" if pair_model == "lin3" else "'i j J'"
        j_text = "J values in cm-1" if pair_model == "lin3" else "J in cm-1"
        prompt = (
            f"Coupled pair {label} ({j_text}; Enter = done)"
            if pairs
            else f"Coupled pair {label} ({j_text}; at least one pair)"
        )
        line = session.ask(prompt)
        if not line:
            break
        try:
            tokens = line.replace(",", " ").split()
            if pair_model == "lin3":
                if len(tokens) != 5:
                    raise ValueError
                pairs.append(
                    (int(tokens[0]), int(tokens[1]), float(tokens[2]),
                     float(tokens[3]), float(tokens[4]))
                )
            else:
                if len(tokens) != 3:
                    raise ValueError
                first_site, second_site, coupling = (
                    int(tokens[0]), int(tokens[1]), float(tokens[2])
                )
                pairs.append((first_site, second_site, coupling))
        except (IndexError, ValueError):
            expected = "'i j Jx Jy Jz'" if pair_model == "lin3" else "'i j J'"
            session.say(
                f"A pair line must read {expected} (the two site indices and the "
                "J value(s)). Cancelled."
            )
            return

    grid = None
    grid_text = session.ask(
        "Susceptibility grid 't_min t_max n_points' (Enter = skip the TINT block)"
    )
    if grid_text:
        try:
            tokens = grid_text.replace(",", " ").split()
            grid = (float(tokens[0]), float(tokens[1]), int(tokens[2]))
        except (IndexError, ValueError):
            session.say("The grid must read 't_min t_max n_points'. Cancelled.")
            return
    try:
        plan = poly_recipe.build_input(
            equivalent_centres=centres,
            spin_orbit_states=states,
            pairs=pairs,
            coordinates=coordinates,
            temperature_grid=grid,
            symmetry=symmetry,
            pair_model=pair_model,
        )
    except poly_recipe.PolyAnisoPlanError as exc:
        session.say(f"Writing the POLY_ANISO input failed: {exc}")
        return
    path.write_text(plan.text, encoding="utf-8")
    body = poly_recipe.render_plan(plan, path=path)
    session.say(body)
    section = ReportSection(
        title="Polynuclear magnetism input (POLY_ANISO exchange cluster)", body=body
    )
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(poly_recipe.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = path.with_name(path.name + ".fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")


def hyperfine_report(session: Session) -> None:
    """Menu 40: hyperfine / EFG report from an ORCA EPRNMR run.

    Reads the electric and magnetic hyperfine structure section: the A
    tensor per nucleus, the EFG principal values with the electron/nuclear
    decomposition, and Rho(0); optionally converts the EFG to the nuclear
    quadrupole coupling constant and the first-order Mossbauer splitting
    with a caller-supplied nuclear quadrupole moment.
    """
    from ...analysis import hyperfine as hyperfine_analysis

    path_text = session.ask("ORCA output path (an EPRNMR run with a Nuclei list)")
    if not path_text:
        session.say("Cancelled (no path given).")
        return
    path = Path(path_text)
    q_text = session.ask(
        "Nuclear quadrupole moment Q (barn) for the quadrupole-splitting conversion (Enter = skip)"
    ).strip()
    q_value: float | None = None
    if q_text:
        try:
            q_value = float(q_text)
        except ValueError:
            session.say(f"Not a number: {q_text!r}; skipping the conversion.")
    try:
        parsed = parse_auto(path)
    except (ParserError, OSError) as exc:
        session.say(f"Hyperfine report failed: {exc}")
        return
    data = parsed.sections.get("hyperfine") or {}
    try:
        body = hyperfine_analysis.render(data, source=path.name, nuclear_Q_barn=q_value)
    except hyperfine_analysis.HyperfineError as exc:
        session.say(f"Hyperfine report failed: {exc}")
        return
    session.say(body)
    section = ReportSection(title="Hyperfine / EFG report", body=body)
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(hyperfine_analysis.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = path.with_name(path.name + ".hyperfine.fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")


def magnetocaloric_report(session: Session) -> None:
    """Menu 41: magnetic entropy and the magnetocaloric effect.

    Routes by content: a POLY_ANISO output with the HINT/TMAG magnetization
    table goes through the Maxwell relation; a SINGLE_ANISO output (or the
    per-center spectra of a POLY_ANISO run) goes through the partition
    function of the spin-orbit levels.
    """
    from ...analysis import magnetocaloric as mc_analysis

    path_text = session.ask("ORCA output path (SINGLE_ANISO levels or a POLY_ANISO M(H) run)")
    if not path_text:
        session.say("Cancelled (no path given).")
        return
    path = Path(path_text)
    try:
        parsed = parse_auto(path)
    except (ParserError, OSError) as exc:
        session.say(f"Magnetocaloric report failed: {exc}")
        return
    poly = parsed.sections.get("poly_aniso") or {}
    single = parsed.sections.get("single_aniso") or {}
    try:
        if poly.get("magnetization"):
            body = mc_analysis.render_maxwell(poly["magnetization"], source=path.name)
            plot_csv = mc_analysis.maxwell_plot_csv(poly["magnetization"])
            plot_kind = (
                "Maxwell -DeltaS(T, H) table; columns "
                "T_mid_K,field_T,minus_delta_S_J_per_mol_per_K"
            )
        else:
            spectrum: list[float] = []
            for segment in single.get("segments") or []:
                if segment.get("soc_spectrum_cm1"):
                    spectrum = segment["soc_spectrum_cm1"]
                    break
            if not spectrum:
                for center in poly.get("centers") or []:
                    if center.get("so_spectrum_cm1"):
                        spectrum = center["so_spectrum_cm1"]
                        break
            body = mc_analysis.render_levels(spectrum, source=path.name)
            plot_csv = mc_analysis.levels_plot_csv(spectrum)
            plot_kind = "S(T) table; columns temperature_K,entropy_J_per_mol_per_K"
    except mc_analysis.MagnetocaloricError as exc:
        session.say(f"Magnetocaloric report failed: {exc}")
        return
    session.say(body)
    section = ReportSection(title="Magnetic entropy / magnetocaloric report", body=body)
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(mc_analysis.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = path.with_name(path.name + ".mce.fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")
    csv_path = path.with_name(path.name + ".mce.fbk.csv")
    csv_path.write_text(plot_csv, encoding="utf-8")
    session.say(f"Plot data written: {csv_path} ({plot_kind})")


def state_tracking_report(session: Session) -> None:
    """Menu 49: track one state across a sequence of runs (density-matrix walk).

    The measure and its declared substitutions live in
    ``analysis.state_tracking``; this handler collects the ordered exports,
    reads each run's per-root energies from its sibling ``.out`` when present,
    and writes the lineage report next to the first export.
    """
    from ...analysis import state_tracking as tracking

    runs_text = session.ask(
        "Run exports in tracking order (orca_2json paths, comma-separated; "
        "first the run that holds the target state)"
    )
    if not runs_text:
        session.say("Cancelled (no paths given).")
        return
    paths = [Path(item.strip()) for item in runs_text.split(",") if item.strip()]
    if len(paths) < 2:
        session.say(
            "State tracking failed: at least two runs are needed. Next step: "
            "give the exports in tracking order, first the target state's run."
        )
        return
    root_text = session.ask(
        "Tracking target: 0-based root index in the first export (Enter = 0)"
    ).strip() or "0"
    try:
        target_root = int(root_text)
    except ValueError:
        session.say(f"Not an integer root index: {root_text!r}; cancelled.")
        return

    def _energies(path: Path) -> dict[tuple[int, int], float]:
        sibling = path.with_suffix(".out")
        if not sibling.exists():
            return {}
        parsed = parse_auto(sibling)
        table: dict[tuple[int, int], float] = {}
        casscf = parsed.sections.get("casscf") or {}
        for block in casscf.get("states", ()):
            for root in block["roots"]:
                table[(block["mult"], root["root"])] = root["energy"]
        return table

    runs = []
    for path in paths:
        try:
            export = parse_orca_json(path)
            runs.append(
                tracking.build_run(export, energies=_energies(path), name=path.name)
            )
        except (ParserError, tracking.StateTrackingError, OSError) as exc:
            session.say(f"State tracking failed on {path.name}: {exc}")
            return
    try:
        result = tracking.track(runs, target_root)
    except tracking.StateTrackingError as exc:
        session.say(f"State tracking failed: {exc}")
        return
    body = tracking.render(result, source=paths[0].name)
    session.say(body)
    section = ReportSection(title="Cross-run state tracking", body=body)
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(tracking.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = paths[0].with_name(paths[0].name + ".track.fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")
