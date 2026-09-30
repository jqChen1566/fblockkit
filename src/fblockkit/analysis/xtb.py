"""the xTB pre-screening report (menu 42).

Renders one xTB run capture into the run facts (task kinds, the energy
chain markers, the frequency set with its imaginary count, the
thermochemistry summary) and the reading notes that fix the position of
this engine in the protocol chain.  The full anchor list and every format
convention live in ``parsers/xtb.py`` (measured, xTB 6.7.1).
"""

from __future__ import annotations

from ..knowledge.models import EVIDENCE_LITERATURE, EVIDENCE_MEASURED, Evidence
from ..parsers.xtb import XtbRun

__all__ = ["render", "evidence"]

_MAX_FREQUENCY_VALUES = 60


def render(run: XtbRun, *, source: str) -> str:
    """The report body (plain text, English, no timestamps of the day)."""
    lines: list[str] = [
        "xTB run report",
        f"Source: {source}",
        "",
        f"Program: xTB {run.version} (build {run.build}, compiled on {run.compiled_on})",
        f"Task: {_task_text(run.tasks)}",
        "",
    ]
    if run.energy_Eh is not None:
        lines.append(f"Total energy: {run.energy_Eh:.12f} Eh")
    if run.gradient_norm is not None:
        lines.append(f"Gradient norm: {run.gradient_norm:.12f} Eh/alpha")
    if run.homo_lumo_gap_eV is not None:
        lines.append(f"HOMO-LUMO gap: {run.homo_lumo_gap_eV:.6f} eV")
    if run.scf_cycles:
        rendered = ", ".join(str(value) for value in run.scf_cycles)
        lines.append(
            f"SCF: converged in {run.scf_cycles[-1]} iteration(s) "
            f"({len(run.scf_cycles)} convergence marker(s): {rendered})"
        )
    if run.opt_converged is True:
        lines.append(f"Geometry optimisation: converged after {run.opt_cycles} iteration(s)")
    elif run.opt_converged is False:
        lines.append(
            f"Geometry optimisation: FAILED to converge in {run.opt_cycles} iteration(s) "
            "-- the run stopped at its --cycles limit; the last geometry is a partial "
            "step, not a minimum"
        )
    lines.extend(_frequency_lines(run))
    lines.extend(_thermo_lines(run))
    if run.finished is not None:
        lines.append(f"Run finished: {run.finished}")
    lines.append("")
    lines.append(
        "Reading notes: xTB (GFN2-xTB here) is the Tier-1 pre-screening level of this "
        "group's protocol chain -- its structures and energies rank and seed, they are "
        "not carried to a higher level as results (a measured example: a GFN2 "
        "transition state re-optimised at r2SCAN-3c showed seven imaginary modes, so a "
        "pre-screening structure is a starting point only). GFN2 covers all fifteen "
        "lanthanides; the actinides have no semi-empirical coverage in this stack "
        "(GFN-FF fails silently for them), so pre-screening actinide systems is out of "
        "scope for this engine. The 'normal termination of xtb' line goes to stderr -- "
        "keep it by redirecting 2>&1 when capturing the run."
    )
    return "\n".join(lines)


def _task_text(tasks: tuple[str, ...]) -> str:
    names = {
        "single_point": "single point",
        "optimization": "geometry optimisation",
        "frequencies": "frequencies",
    }
    return " + ".join(names.get(task, task) for task in tasks)


def _frequency_lines(run: XtbRun) -> list[str]:
    if not run.frequency_groups:
        return []
    values = run.frequencies
    imaginary = [value for value in values if value < 0]
    lines = [
        "",
        f"Frequencies: {len(values)} mode(s) in the harmonic set; "
        f"imaginary: {len(imaginary)}",
    ]
    if run.engine_n_imaginary is not None:
        lines.append(
            f"  Engine counters: {run.engine_n_frequencies} vibration(s), "
            f"{run.engine_n_imaginary} imaginary"
        )
    shown = values[:_MAX_FREQUENCY_VALUES]
    rendered = "  ".join(f"{value:.2f}" for value in shown)
    if len(values) > _MAX_FREQUENCY_VALUES:
        rendered += f"  ... ({len(values) - _MAX_FREQUENCY_VALUES} more)"
    lines.append(f"  {rendered}")
    if run.frequency_groups and len(run.frequency_groups) > 1:
        lines.append(
            f"  (the block is printed {len(run.frequency_groups)} times; the identical "
            "copies are kept as measured)"
        )
    if imaginary:
        lines.append(
            f"  Most negative mode: {min(imaginary):.2f} cm-1"
            + (
                f"; engine: found {run.imag_found} significant imaginary frequency"
                if run.imag_found is not None
                else ""
            )
        )
        lines.append(
            "  An imaginary mode marks a saddle or a non-stationary geometry at this "
            "level; a transition-state candidate must be re-verified at the target level."
        )
    return lines


def _thermo_lines(run: XtbRun) -> list[str]:
    if not run.thermo:
        return []
    names = (
        ("free_energy_Eh", "total free energy"),
        ("total_energy_Eh", "total energy"),
        ("zpe_Eh", "zero point energy"),
    )
    lines = ["", "Thermochemistry (from the engine's printed table):"]
    for key, label in names:
        if key in run.thermo:
            lines.append(f"  {label}: {run.thermo[key]:.12f} Eh")
    return lines


def evidence() -> tuple[Evidence, ...]:
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "The xTB program citation printed by the engine itself (banner): program "
                "paper of the engine whose capture this menu reads."
            ),
            ref="Bannwarth C., Caldeweyher E., Ehlert S., Hansen A., Pracht P., Seibert J., Spicher S., Grimme S., WIREs Comput. Mol. Sci., 2021, 11, e01493, DOI 10.1002/wcms.1493",
            bibkey="bannwarth2021xtb",
        ),
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "The GFN2-xTB parametrisation paper (printed in the engine banner for "
                "the GFN2 level used by the fixtures)."
            ),
            ref="Bannwarth C., Ehlert S., Grimme S., J. Chem. Theory Comput., 2019, 15, 1652-1671, DOI 10.1021/acs.jctc.8b01176",
            bibkey="bannwarth2019gfn2",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "Fixture anchors measured on xTB 6.7.1 (101, 2026-09-30): the version "
                "banner, the result box, both optimisation markers (converged after 5; "
                "FAILED TO CONVERGE with --cycles 2), the doubly printed frequency "
                "block (water: nine values, imaginary 0; planar ammonia: imaginary 1, "
                "most negative -1337.66 cm-1) and the closing marker."
            ),
            ref="fixtures/xtb/ (README records the run commands)",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "The cross-level boundary is this group's own record: a GFN2 transition "
                "state re-optimised at r2SCAN-3c carried seven imaginary modes (SCINE "
                "stack assessment, cjq11); lanthanide coverage of GFN2 and the silent "
                "failure of GFN-FF on actinides are recorded in the same assessment."
            ),
            ref="SCINE capability assessment (this group, 2026-09-20/21)",
        ),
    )
