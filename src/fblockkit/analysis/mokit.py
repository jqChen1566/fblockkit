"""the MOKIT/automr run report (menu 45).

Renders one automr capture: the program-path survey, the run settings and
merged mokit{} options, the final strategy table, the stage sequence, the
energy chain, the automatically determined active space, the radical-index
tables and the termination state -- plus the reading notes that state what
was and was not read (the .fch side products are listed and parsed; see
``parsers/mokit_fch.py`` for their measured layout).
Every anchor is measured (``parsers/mokit.py``; MOKIT 1.2.8, 2026-09-30).
"""

from __future__ import annotations

from ..knowledge.models import EVIDENCE_MEASURED, Evidence
from ..parsers.mokit import MokitRun
from ..parsers.mokit_fch import stage_label

__all__ = ["render", "evidence"]

_RADICAL_LABELS = {
    "biradical": "biradical character y0",
    "tetraradical": "tetraradical character y1",
    "yamaguchi_unpaired": "Yamaguchi unpaired electrons",
    "head_gordon_min": "Head-Gordon unpaired (min)",
    "head_gordon_squared": "Head-Gordon unpaired (squared)",
}

_STAGE_TITLES = {
    "do_hf": "HF",
    "get_paired_LMO": "UNO rotation",
    "do_gvb": "GVB",
    "do_cas": "CASSCF/CASCI",
}


def render(run: MokitRun, *, source: str, fch_files: tuple = ()) -> str:
    """The report body (plain text, English, no timestamps of the day).

    ``fch_files``: the .fch side products read next to the output (an iterable
    of :class:`fblockkit.parsers.mokit_fch.FchFile`, or of ready-made problem
    strings for files that could not be read).
    """
    lines: list[str] = [
        "MOKIT automr run report",
        f"Source: {source}",
        "",
        f"Version: {run.version} ({run.built})",
    ]
    if run.program_paths:
        found = ", ".join(
            f"{name.split('_path')[0]}={value}"
            for name, value in sorted(run.program_paths.items())
        )
        lines.append(f"Program paths: {found}")
    if run.method_basis:
        lines.append(
            f"Run settings: memory {run.memory}, nproc {run.nproc}, method/basis "
            f"{run.method_basis}"
        )
    if run.keywords:
        lines.append(f"MOKIT options (merged): {run.keywords}")
    else:
        lines.append("MOKIT options (merged): (none)")
    if run.strategy_number is not None:
        active = [name for name, value in run.strategy_flags.items() if value]
        lines.append(
            f"Strategy (No. {run.strategy_number}) active flags: "
            + (", ".join(active) if active else "(none)")
        )
    if run.stages:
        titled = [_STAGE_TITLES.get(stage, stage) for stage in run.stages]
        lines.append("Stages: " + " -> ".join(titled))

    if run.energies:
        lines.append("")
        lines.append("Energy chain:")
        for label, value in run.energies:
            extra = ""
            for s2_label, s2_value in run.s2_values:
                if s2_label == label:
                    extra = f"  <S**2> = {s2_value:.3f}"
                    break
            lines.append(f"  {label:<7} {value:.8f} Eh{extra}")

    if run.gvb_order is not None:
        lines.append(
            f"GVB: {run.gvb_order} pair(s), program {run.gvb_program}"
        )
    if run.active_space is not None:
        electrons, orbitals = run.active_space
        lines.append(
            f"Active space (automatically determined): CAS({electrons}e,{orbitals}o), "
            f"program {run.casscf_program}"
        )

    for position, group in enumerate(run.radical_groups, start=1):
        lines.append("")
        lines.append(f"Radical index (table {position} of {len(run.radical_groups)}):")
        for name, value in group:
            lines.append(f"  {_RADICAL_LABELS.get(name, name)}: {value:.3f}")

    if run.side_products:
        lines.append("")
        lines.append("Side products echoed by the run (audit trail, not parsed here):")
        for name in run.side_products:
            lines.append(f"  {name}")

    if fch_files:
        lines.append("")
        lines.append(
            f".fch side products read next to the output ({len(fch_files)} file(s)):"
        )
        for item in fch_files:
            if isinstance(item, str):
                lines.append(f"  - {item}")
                continue
            spin = (
                f"{item.n_alpha}alpha/{item.n_beta}beta"
                if item.n_alpha is not None and item.n_beta is not None
                else "spin counts not carried"
            )
            mo_note = (
                f"{len(item.alpha_energies)} alpha MOs, no beta block"
                if item.beta_energies is None
                else f"{len(item.alpha_energies)} alpha + "
                f"{len(item.beta_energies)} beta MOs"
            )
            lines.append(
                f"  - {item.name} -- {stage_label(item.name)} (name-derived): "
                f"nbf {item.nbf}; {len(item.atomic_numbers)} atoms; charge "
                f"{item.charge}, mult {item.multiplicity}, {item.n_electrons} e "
                f"({spin}); {mo_note}"
            )
            facts: list[str] = []
            if item.scf_energy is not None:
                facts.append(f"SCF {item.scf_energy:.8f} Eh")
            if (
                item.total_energy is not None
                and item.total_energy != item.scf_energy
            ):
                facts.append(f"total {item.total_energy:.8f} Eh")
            if item.dipole_au is not None:
                facts.append(
                    "dipole (as stored) "
                    f"({item.dipole_au[0]:.3f}, {item.dipole_au[1]:.3f}, "
                    f"{item.dipole_au[2]:.3f})"
                )
            if facts:
                lines.append("      " + "; ".join(facts))
        lines.append(
            "    Notes: the stage notes are derived from the file names (MOKIT's "
            "naming); coordinates and the dipole vector are reported as stored "
            "(the file's own units; the coordinates are in Bohr); the format "
            "carries no occupation numbers."
        )

    lines.append("")
    if run.terminated:
        lines.append(
            "Termination: Normal termination of AutoMR -- the workflow ran to its end."
        )
    else:
        lines.append(
            "Termination: no 'Normal termination of AutoMR' line found -- the run "
            "stopped early; treat every number above as partial."
        )

    lines.append("")
    lines.append(
        "Reading notes: the active space on the CASSCF line is the automatically "
        "determined selection (GVB natural-orbital occupations above the 0.02 "
        "threshold); writing CASSCF(n,m) in the route line pins the size instead. "
        "The energy chain shows the flow RHF/UHF (the lower one is kept) -> GVB -> "
        "CASCI/CASSCF, so consecutive entries are different wave-function levels, "
        "not an error. GVB runs need a backend program (GAMESS by default; Gaussian "
        "and QChem are the alternates -- PySCF is not a GVB backend); the CASSCF "
        "stage defaults to PySCF. The natural-orbital .fch side products are read "
        "back when they sit next to the output (the section above; the reader is "
        "parsers/mokit_fch.py, whose anchor policy is in its docstring); "
        "dynamic-correlation stages (CASPT2/NEVPT2/DMRG) appear in the strategy "
        "table only when requested."
    )
    return "\n".join(lines)


def evidence() -> tuple[Evidence, ...]:
    return (
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "Fixture anchors measured on MOKIT 1.2.8 (101, 2026-09-30; H2O "
                "CASSCF/cc-pVDZ with mokit{GVB_prog=Gaussian} on Gaussian 16 + "
                "PySCF): the banner and version line, the program-path block "
                "(gms_path NOT FOUND), the merged-options line, both strategy-table "
                "printings (No. 0 / No. 1), the stage sequence do_hf -> "
                "get_paired_LMO -> do_gvb -> do_cas, the energy chain "
                "RHF/UHF/UHF2/GVB/CASCI/CASSCF, the CAS(4e,4o) selection, the three "
                "Radical-index tables and the closing 'Normal termination of AutoMR' "
                "line."
            ),
            ref="fixtures/mokit/ (README records the run command and the environment)",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "The citation form is the project's own registration (MOKIT has no "
                "program paper; the official citation is 'Jingxiang Zou, Molecular "
                "Orbital Kit (MOKIT)', printed by the README): the tool index entry "
                "records the same."
            ),
            ref="调研/合成_cjq12_全景与方案_20260925.md section 1 (MOKIT entry)",
            url="https://gitlab.com/jxzou/mokit",
        ),
    )
