"""The selection handler group."""

from __future__ import annotations

import json
from pathlib import Path

from ...analysis import aegiss, apc, ass1st, dm_selection, geometry as geometry_analysis, pios, qicas, tnass
from ...diagnosis import references_section
from ...knowledge.models import ReportSection
from ...parsers import ParserError, parse_auto
from ...parsers.fcidump import parse_fcidump
from ...parsers.mkl import parse_mkl
from ...parsers.orca_json import parse_orca_json
from ...recipe import ass1st as ass1st_recipe, dm_batch
from .common import casscf_reference
from ..session import Session
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
    apc_n_text = session.ask(
        "Balanced variant: set aside the N highest-entropy virtuals per entropy "
        "evaluation (Enter = 0, the classic scheme; 2 is the source's recommendation "
        "for large molecules)"
    ).strip()
    try:
        apc_n = int(apc_n_text) if apc_n_text else 0
    except ValueError:
        session.say(f"Cancelled (APC-N must be an integer, got {apc_n_text!r}).")
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
            apc_n=apc_n,
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
        ok, reference = casscf_reference(session, output_text)
        if not ok:
            return
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


def aegiss_select(session: Session) -> None:
    """Menu 25: the AEGISS selection -- entropy screening + an AO projection."""
    export_text = session.ask(
        "orca_2json export path (the gbw whose orbitals were dumped; needs S and the "
        "orbital labels)"
    )
    if not export_text:
        session.say("Cancelled (no export path given).")
        return
    fcidump_text = session.ask("FCIDUMP path (the same run's !FCIDUMP dump)")
    if not fcidump_text:
        session.say("Cancelled (no FCIDUMP path given).")
        return
    output_text = session.ask(
        "CASSCF output path (matches the CI root to the printed energy; Enter = skip)"
    ).strip()
    label_text = session.ask(
        "AO label: element + angular momentum (+ optional shell/component), e.g. "
        "'C pz' or 'Fe d'"
    ).strip()
    if not label_text:
        session.say("Cancelled (no AO label given).")
        return
    tau_text = session.ask(
        "Entropy fraction tau (line = tau * S_max; Enter = 0.1, the source's default)"
    ).strip()
    epsilon_text = session.ask(
        "Projection threshold on the weight (Enter = 0.5, the source's value)"
    ).strip()
    try:
        ok, reference = casscf_reference(session, output_text)
        if not ok:
            return
        export = parse_orca_json(Path(export_text))
        dump = parse_fcidump(Path(fcidump_text))
        section = aegiss.run(
            export,
            dump,
            label=label_text,
            tau=float(tau_text) if tau_text else aegiss.AEGISS_DEFAULT_TAU,
            epsilon=float(epsilon_text) if epsilon_text else aegiss.AEGISS_DEFAULT_EPSILON,
            reference_energy=reference,
        )
    except (ParserError, aegiss.AegissError, OSError, ValueError) as exc:
        session.say(f"AEGISS selection failed: {exc}")
        return
    session.say(section.body)
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(aegiss.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = Path(export_text).with_name(Path(export_text).name + ".aegiss.fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")


def tnass_select(session: Session) -> None:
    """Menu 26: the TNASS subset selection (Renyi-2 bipartition entropy)."""
    fcidump_text = session.ask(
        "FCIDUMP path (a converged CASSCF dump; the subset is selected within it)"
    )
    if not fcidump_text:
        session.say("Cancelled (no FCIDUMP path given).")
        return
    output_text = session.ask(
        "CASSCF output path (matches the CI root to the printed energy; Enter = skip)"
    ).strip()
    size_text = session.ask("Target size: the number of active spatial orbitals").strip()
    method_text = session.ask(
        "Method: greedy, 'block K' (e.g. 'block 2'), or brute (Enter = greedy)"
    ).strip().lower()
    try:
        ok, reference = casscf_reference(session, output_text)
        if not ok:
            return
        n_target = int(size_text)
        tokens = method_text.split()
        method, block_size = "greedy", None
        if tokens == ["brute"]:
            method = "brute"
        elif len(tokens) == 2 and tokens[0] == "block":
            method, block_size = "block", int(tokens[1])
        elif tokens not in ([], ["greedy"]):
            raise tnass.TnassError(
                f"unknown method {method_text!r}; use 'greedy', 'block K', or 'brute'."
            )
        dump = parse_fcidump(Path(fcidump_text))
        section = tnass.run(
            dump,
            n_target=n_target,
            method=method,
            block_size=block_size,
            reference_energy=reference,
        )
    except (ParserError, tnass.TnassError, OSError, ValueError) as exc:
        session.say(f"TNASS selection failed: {exc}")
        return
    session.say(section.body)
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(tnass.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = Path(fcidump_text).with_name(Path(fcidump_text).name + ".tnass.fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")




# --- 48 PiOS pi-orbital active space -----------------------------------------


def pios_select(session: Session) -> None:
    """Menu 48: build the pi-orbital active space of a conjugated system (PiOS)."""
    path_text = session.ask(
        "PiOS manifest JSON path (an export, the pi-system atoms and the mkl; see "
        "the user guide)"
    )
    if not path_text:
        session.say("Cancelled (no manifest given).")
        return
    path = Path(path_text)

    def _resolve(entry) -> Path:
        candidate = Path(str(entry))
        return candidate if candidate.is_absolute() else path.parent / candidate

    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        export = parse_orca_json(_resolve(payload["export"]))
        atoms = [int(index) for index in payload["atoms"]]
        template_mkl_path = _resolve(payload["mkl"])
        template_mkl = parse_mkl(template_mkl_path)
        result = pios.select_pi_space(
            export,
            atoms,
            charge=int(payload["charge"]) if "charge" in payload else None,
            contributions=payload.get("contributions"),
            pi_electrons=(
                int(payload["pi_electrons"]) if "pi_electrons" in payload else None
            ),
        )
    except (OSError, KeyError, TypeError, ValueError, ParserError) as exc:
        session.say(f"PiOS selection failed: {exc}")
        return
    body = pios.render(result)
    session.say(body)
    mkl_path = template_mkl_path.with_name(template_mkl_path.stem + ".pios.fbk.mkl")
    try:
        pios.write_mkl(template_mkl, result, mkl_path)
    except (OSError, pios.PiosError) as exc:
        session.say(f"Writing the pi-space mkl failed: {exc}")
        return
    session.say(f"Pi space written: {mkl_path}")
    n_inactive = sum(1 for value in export.mo_occupations if value > 1.99) - result.n_occupied_pi
    session.say(
        f"Next step: run ``orca_2mkl {mkl_path.stem} -gbw`` next to that file, then "
        f"start the MCSCF with ``!moread`` and %moinp \"{mkl_path.stem}.gbw\" — the "
        f"written partition is {n_inactive} inactive, then the "
        f"{len(result.indices)} pi orbitals as the active window."
    )
    section = ReportSection(
        title="A13 PiOS pi-orbital active space (oriented-p projection)", body=body
    )
    report_lines = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(pios.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = path.with_name(path.name + ".pios.fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")
