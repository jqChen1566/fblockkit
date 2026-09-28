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
from ...parsers.orca_json import parse_orca_json
from ...recipe import guess_transfer
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
        title="4.1 perturbed multistart batch (randomized occupied-virtual mixing)",
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
    session.say(body)
    section = ReportSection(title="pNMR pseudocontact shifts", body=body)
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(pnmr_analysis.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = path.with_name(path.name + ".pnmr.fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")
