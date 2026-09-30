"""the B-layer evaluation of pysisyphus -- read one of its runs (menu 33).

The B-layer relation is "generate its input / read its output" (no code is
linked in; pysisyphus is GPL-3.0 and stays an external program).  This module
is the reading half: it takes a run directory -- or the console capture alone
-- and reports the run against itself, with the cross-checks the measured
format layer allows:

- the trajectory frame count against the console's cycle count;
- the last frame's energy against the calculator call that produced it (the
  run's final call when it converged there, the second-to-last one when it
  stopped -- a stopped run evaluates one extra geometry, the closing state,
  after its last cycle; measured);
- the closing state (``Final summary``) against the **last** calculator call;
- the closing geometry file against the closing state;
- for a converged run, the closing geometry's coordinates against the closing
  frame (a stopped run's closing geometry is the extra evaluation and has no
  frame; measured).

The calculator calls are the per-cycle ORCA outputs
(``qm_calcs/calculator_*.orca.out``), parsed by this toolkit's ORCA parser --
so the run's own energies are re-derived from ORCA's files rather than trusted
to the console text.

The evaluation's boundary survey (all measured on pysisyphus 1.0.0 driving
ORCA 6.1.1, 2026-09-28; the fixtures under ``fixtures/pysisyphus/``):

- **Gradient-driven work is exact.**  Minimum optimisation and -- with a model
  Hessian -- TS optimisation run end to end; the TS fixture reproduces the
  energy asserted in pysisyphus's own example set (Baker test-set case 11,
  ``-154.050455732882``) to 7e-8.
- **A quantum Hessian crashes the pairing.**  ORCA 6 writes an extra
  ``$multiplicity`` block in its ``.hess`` file, and pysisyphus's pyparsing
  grammar for that file (measured identical in release 1.0.0 and on master)
  stops at it: ``ParseException: Expected W:(0-9), found '$'`` -- so
  ``hessian_init: calc`` (the ``tsopt`` default), ``do_hess`` and every
  frequency-dependent step of the run end in
  ``RunAfterCalculationFailedException``.  The affected jobs need a model
  Hessian (``hessian_init: fischer``) and an outside frequency check (a plain
  ORCA ``Freq`` run; the menu-30 chain is the toolkit's own re-optimisation
  route).  This is a defect of the interface, not of the science: ORCA's own
  ``.engrad`` sidecar parses fine, which is why everything gradient-driven is
  reliable.
- ``optimization.h5`` carries the full per-cycle history but is HDF5; h5py is
  not a dependency of this toolkit, so the reader takes the same numbers from
  the text artifacts instead.

The crash classification is part of the report: a run that ended in the
post-processing is named as such, and the ``$multiplicity`` signature is
recognised (with the backup directory searched for the offending ``.hess``
when it is still present).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from ..knowledge.models import EVIDENCE_LITERATURE, EVIDENCE_MEASURED, Evidence
from ..parsers import ParserError, parse_auto
from ..parsers.pysisyphus import (
    Frame,
    PysisyphusError,
    find_run_log,
    read_log,
    read_run_record,
    read_trj,
    read_xyz,
)

__all__ = [
    "PysisyphusRunError",
    "Check",
    "RunData",
    "read_run",
    "render",
    "evidence",
]

_ENERGY_TOL = 5e-9  # the artifacts print energies at 8 decimals
_COORD_TOL = 1e-8
_ORCA_TOL = 1e-6  # ORCA prints more decimals than the console summary


class PysisyphusRunError(ValueError):
    """The run directory/capture cannot be read as a pysisyphus run (with a next step)."""


@dataclass(frozen=True)
class Check:
    """One cross-check between the run's own artifacts."""

    name: str
    ok: bool | None  # None = not applicable (an artifact was not given)
    detail: str


@dataclass(frozen=True)
class OrcaCall:
    """One calculator call of the run (the fields the checks need)."""

    path: Path
    energy: float | None
    terminated: bool | None


@dataclass(frozen=True)
class RunData:
    """One pysisyphus run, read and cross-checked."""

    log_path: Path
    directory: Path
    facts: dict
    record: dict | None = None
    frames: tuple[Frame, ...] | None = None
    final_frame: Frame | None = None
    orca_paths: tuple[Path, ...] = ()
    orca_calls: tuple[OrcaCall, ...] = ()  # the last two calls, parsed
    crash_kind: str | None = None
    crash_hess: Path | None = None
    checks: tuple[Check, ...] = field(default_factory=tuple)


def _trajectory_path(directory: Path, kind: str | None) -> Path | None:
    """The run's trajectory: the kind-matching name, else the only non-COS .trj."""
    preferred = "ts_optimization.trj" if kind and "TS" in kind.upper() else "optimization.trj"
    candidate = directory / preferred
    if candidate.is_file():
        return candidate
    others = [
        path
        for path in sorted(directory.glob("*.trj"))
        if not any(tag in path.name for tag in ("cos_hei", "current", "cycle_", "image_"))
    ]
    if len(others) == 1:
        return others[0]
    return None


def _final_geometry_path(directory: Path, facts: dict) -> Path | None:
    """The closing geometry: the console's own file name, else the unique candidate."""
    if (name := facts.get("final_geometry_fn")) is not None:
        candidate = directory / name
        if candidate.is_file():
            return candidate
    others = sorted(directory.glob("*final_geometry.xyz"))
    if len(others) == 1:
        return others[0]
    return None


def _orca_outputs(directory: Path) -> tuple[Path, ...]:
    """The ``qm_calcs/calculator_*.orca.out`` files, ordered by call index."""
    qm_calcs = directory / "qm_calcs"
    if not qm_calcs.is_dir():
        return ()
    keyed: list[tuple[tuple[int, int], Path]] = []
    for path in qm_calcs.glob("calculator_*.orca.out"):
        parts = path.stem.split(".")
        try:
            key = (int(parts[0].rsplit("_", 1)[1]), int(parts[1]))
        except (IndexError, ValueError):
            continue
        keyed.append((key, path))
    return tuple(path for _, path in sorted(keyed))


def _find_crashed_backup(directory: Path) -> Path | None:
    """The ``crashed_*`` backup directory of this run, when it is still around."""
    candidates = sorted(path for path in directory.glob("crashed_*") if path.is_dir())
    return candidates[0] if len(candidates) == 1 else None


def _classify_crash(text: str, exceptions: tuple[str, ...]) -> str:
    """Name the crash: the measured ORCA-6 Hessian parse stop, or the generic one."""
    if any("ParseException" in line for line in exceptions) or "parse_hess_file" in text:
        return "hessian-parse"
    return "postprocessing"


def _frame_geometry_matches(one: Frame, other: Frame, tol: float = _COORD_TOL) -> bool:
    if len(one.atoms) != len(other.atoms):
        return False
    for (symbol_a, *xyz_a), (symbol_b, *xyz_b) in zip(one.atoms, other.atoms):
        if symbol_a.lower() != symbol_b.lower():
            return False
        if any(abs(a - b) > tol for a, b in zip(xyz_a, xyz_b)):
            return False
    return True


def read_run(path: str | Path) -> RunData:
    """Read a pysisyphus run: the given capture file, or the run directory.

    The console capture is required (it carries the outcome marker, the cycle
    table and the closing summary); the trajectory, ``RUN.yaml``, the closing
    geometry and the calculator outputs are picked up from the directory when
    present and cross-checked.
    """
    source = Path(path)
    log_path = find_run_log(source)
    directory = log_path.parent
    try:
        text = log_path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise PysisyphusRunError(f"the capture cannot be read: {exc}") from exc
    facts = read_log(text)

    record = None
    record_path = directory / "RUN.yaml"
    if record_path.is_file():
        try:
            record = read_run_record(record_path.read_text(encoding="utf-8"))
        except PysisyphusError:
            record = None  # a foreign RUN.yaml is not part of the run being read

    frames = None
    if (trj_path := _trajectory_path(directory, facts["kind"])) is not None:
        try:
            frames = read_trj(trj_path.read_text(encoding="utf-8", errors="replace"))
        except PysisyphusError as exc:
            raise PysisyphusRunError(
                f"the trajectory {trj_path.name} does not parse: {exc}"
            ) from exc

    final_frame = None
    final_path = _final_geometry_path(directory, facts)
    if final_path is not None:
        try:
            final_frame = read_xyz(final_path.read_text(encoding="utf-8", errors="replace"))
        except PysisyphusError as exc:
            raise PysisyphusRunError(
                f"the closing geometry {final_path.name} does not parse: {exc}"
            ) from exc

    orca_paths = _orca_outputs(directory)
    orca_calls: list[OrcaCall] = []
    for candidate in orca_paths[-2:]:
        try:
            result = parse_auto(candidate)
        except ParserError:
            continue
        orca_calls.append(
            OrcaCall(
                path=candidate,
                energy=result.sections.get("final_energy"),
                terminated=result.sections.get("terminated_normally"),
            )
        )
    last_call = orca_calls[-1] if orca_calls else None
    converged_before_crash = bool(
        facts["crash"] is not None and facts["crash"]["converged_before_crash"]
    )
    if facts["outcome"] == "converged" or converged_before_crash:
        cycle_call = last_call  # the run converged at its last cycle's geometry
    elif len(orca_calls) > 1:
        cycle_call = orca_calls[-2]  # the extra call is the closing evaluation
    else:
        cycle_call = None

    crash_kind = None
    crash_hess = None
    if facts["crash"] is not None:
        crash_kind = _classify_crash(text, facts["crash"]["exceptions"])
        backup = _find_crashed_backup(directory)
        if backup is not None:
            hess_files = sorted(backup.glob("*.hess"))
            if hess_files:
                crash_hess = hess_files[0]

    checks: list[Check] = []
    cycles = facts["cycles"]
    if frames is not None:
        checks.append(
            Check(
                name="trajectory length vs cycle count",
                ok=len(frames) == len(cycles),
                detail=f"{len(frames)} frame(s) vs {len(cycles)} cycle row(s)",
            )
        )
    # The last cycle's geometry is the geometry of the last cycle-producing
    # calculator call: the last call itself when the run converged there, and
    # the second-to-last call when the run stopped (the stopped run evaluates
    # one extra geometry -- the closing state -- after the last cycle;
    # measured on the fixture: 3 cycle rows against 4 calls, the extra call's
    # energy matching the closing summary).
    if frames and cycle_call is not None and frames[-1].energy is not None:
        energy = cycle_call.energy
        ok = energy is not None and abs(frames[-1].energy - energy) <= _ORCA_TOL
        checks.append(
            Check(
                name="closing frame energy vs its calculator call",
                ok=ok,
                detail=f"frame {frames[-1].energy:.8f} Eh vs {cycle_call.path.name}"
                + (f" {energy:.9f} Eh" if energy is not None else " (no energy parsed)"),
            )
        )
    else:
        checks.append(
            Check(
                name="closing frame energy vs its calculator call",
                ok=None,
                detail=(
                    "the run directory holds no trajectory file"
                    if frames is None
                    else "no per-cycle calculator call could be identified next to the "
                    "capture"
                ),
            )
        )
    if last_call is not None and "energy" in facts["final"]:
        summary = facts["final"]["energy"]
        energy = last_call.energy
        ok = energy is not None and abs(energy - summary) <= _ORCA_TOL
        checks.append(
            Check(
                name="closing state (final summary) vs the last calculator call",
                ok=ok,
                detail=f"summary {summary:.8f} Eh vs {last_call.path.name}"
                + (f" {energy:.9f} Eh" if energy is not None else " (no energy parsed)"),
            )
        )
    else:
        checks.append(
            Check(
                name="closing state (final summary) vs the last calculator call",
                ok=None,
                detail="no qm_calcs/calculator_*.orca.out found next to the capture",
            )
        )
    if final_frame is not None and "energy" in facts["final"]:
        summary = facts["final"]["energy"]
        energy = final_frame.energy
        ok = energy is not None and abs(energy - summary) <= _ENERGY_TOL
        checks.append(
            Check(
                name="closing geometry file vs the closing state",
                ok=ok,
                detail=f"geometry {energy!r} vs summary {summary:.8f} Eh",
            )
        )
    if frames and final_frame is not None:
        if facts["outcome"] == "converged" or converged_before_crash:
            ok = _frame_geometry_matches(frames[-1], final_frame)
            detail = "atoms and coordinates within 1e-8 Angstrom"
        else:
            ok = None
            detail = (
                "the stopped run's closing geometry is the extra evaluation after "
                "the last cycle, so it has no trajectory frame (measured)"
            )
        checks.append(
            Check(
                name="closing geometry coordinates vs the closing frame",
                ok=ok,
                detail=detail,
            )
        )

    return RunData(
        log_path=log_path,
        directory=directory,
        facts=facts,
        record=record,
        frames=frames,
        final_frame=final_frame,
        orca_paths=orca_paths,
        orca_calls=tuple(orca_calls),
        crash_kind=crash_kind,
        crash_hess=crash_hess,
        checks=tuple(checks),
    )


_CHECK_MARK = {True: "ok", False: "DIFFERS", None: "n/a"}


def _cycle_table(facts: dict) -> list[str]:
    """The cycle table as a compact ASCII block (stars preserved)."""
    if not facts["cycles"]:
        return ["  (the capture carries no cycle table)"]
    header = (
        "  cycle   d(energy)  max|force|  rms|force|   max|step|   rms|step|"
        "  s/cycle"
    )
    lines = [header]

    def cell(value: float | None, starred: bool) -> str:
        text = "nan" if value is None else f"{value:.6f}"
        return f"{text}{'*' if starred else ' '}"

    for row in facts["cycles"]:
        seconds = "-" if row.seconds is None else f"{row.seconds:.3f}"
        lines.append(
            f"  {row.cycle:5d}   "
            f"{cell(row.d_energy, 'd_energy' in row.starred):>12} "
            f"{cell(row.max_force, 'max_force' in row.starred):>11} "
            f"{cell(row.rms_force, 'rms_force' in row.starred):>11} "
            f"{cell(row.max_step, 'max_step' in row.starred):>11} "
            f"{cell(row.rms_step, 'rms_step' in row.starred):>11} "
            f"{seconds:>8}"
        )
    return lines


def render(data: RunData) -> str:
    """The menu's report: the run's own facts, the cross-checks, the boundaries."""
    facts = data.facts
    out: list[str] = [f"pysisyphus run: {data.directory}"]
    if facts["version"]:
        program = f"pysisyphus {facts['version']}"
        if facts["python"]:
            program += f" (Python {facts['python']})"
        if facts["executed_at"]:
            program += f"; executed at {facts['executed_at']} on '{facts['host']}'"
        out.append(f"  program: {program}")
    job = f"  job: {facts['kind']}"
    if facts["system"]:
        job += f"; system: {facts['system']} ({facts['n_atoms']} atoms)"
    if facts["coord_type"]:
        job += f"; coordinate system: {facts['coord_type']} ({facts['n_coord']})"
    out.append(job)
    calculator = "  calculator:"
    if facts["calculator"]:
        calculator += f" {facts['calculator']} ({facts['calculator_id']})"
    if facts["charge"] is not None:
        calculator += f"; charge {facts['charge']}, multiplicity {facts['multiplicity']}"
    if facts["optimizer"]:
        calculator += f"; optimizer: {facts['optimizer']}"
    out.append(calculator)
    if facts["thresholds"]:
        labels = {
            "max_force": "max(|force|)",
            "rms_force": "rms(force)",
            "max_step": "max(|step|)",
            "rms_step": "rms(step)",
        }
        parts = []
        for name, label in labels.items():
            if name not in facts["thresholds"]:
                continue
            value, overachieved = facts["thresholds"][name]
            text = f"{label} <= {value:g}"
            if overachieved is not None:
                text += f" (overachieved {overachieved:g})"
            parts.append(text)
        out.append("  convergence thresholds: " + ", ".join(parts))
    if data.record is not None:
        calc = data.record.get("calc") or {}
        opt = data.record.get("opt") or data.record.get("tsopt") or {}
        record_line = (
            f"  run record (RUN.yaml): {calc.get('type', '?')} "
            f"'{calc.get('keywords', '?')}', pal {calc.get('pal', '?')}, "
            f"mem {calc.get('mem', '?')} MB"
        )
        if "thresh" in opt:
            record_line += f", thresh {opt['thresh']}"
        if "max_cycles" in opt:
            record_line += f", max_cycles {opt['max_cycles']}"
        if data.record.get("version"):
            record_line += f" (record version {data.record['version']})"
        out.append(record_line)
    out.append("")
    out += _cycle_table(facts)
    out.append("")
    if facts["crash"] is not None:
        crash = facts["crash"]
        out.append("outcome: the run crashed in its post-processing.")
        if crash["converged_before_crash"]:
            out.append(
                "  the optimisation part itself had already printed its converged marker;"
                " the failure came after it."
            )
        if data.crash_kind == "hessian-parse":
            out.append(
                "  class: the ORCA-6 Hessian parse stop -- ORCA 6.1.1 writes a "
                "'$multiplicity' block into its .hess file and this pysisyphus "
                "interface's grammar stops at it (measured identical in 1.0.0 and on "
                "master): ParseException at that block, then "
                "RunAfterCalculationFailedException."
            )
            if data.crash_hess is not None:
                hess_text = data.crash_hess.read_text(encoding="utf-8", errors="replace")
                carries = "$multiplicity" in hess_text
                out.append(
                    f"  the backup {data.crash_hess} "
                    + ("carries the '$multiplicity' block (the signature is confirmed)."
                       if carries else
                       "does not carry '$multiplicity' (the signature is not confirmed).")
                )
            out.append(
                "  next steps: use a model Hessian for the TS search "
                "(hessian_init: fischer -- the tsopt default is 'calc') and verify the "
                "frequencies outside pysisyphus (a plain ORCA Freq run on the closing "
                "geometry; menu 30 is the toolkit's own re-optimisation route)."
            )
        else:
            out.append(
                "  class: a post-processing exception; next step: inspect the "
                "crashed_* backup directory copied next to the capture."
            )
        if crash["exceptions"]:
            out.append(f"  exception: {crash['exceptions'][-1]}")
    else:
        outcome = {
            "converged": "converged (the program's own marker)",
            "cycles exceeded": (
                "the number of cycles was exceeded -- not converged, even where the "
                "closing force values already look small (the marker is the ground "
                "truth, measured on the stopping fixture)"
            ),
            "unknown": (
                "no closing marker recognised (the capture may be truncated); "
                "next step: give the full console capture of the run."
            ),
        }.get(facts["outcome"], facts["outcome"])
        out.append(f"outcome: {outcome}")
    final = facts["final"]
    if "energy" in final:
        line = f"  final: energy {final['energy']:.8f} Eh"
        if "max_forces_internal" in final:
            line += f"; max|force| {final['max_forces_internal']:.6f} (internal)"
        if "max_forces_cartesian" in final:
            line += f" / {final['max_forces_cartesian']:.6f} (cartesian)"
        out.append(line)
    artifacts = []
    if data.frames is not None:
        artifacts.append(f"trajectory ({len(data.frames)} frames)")
    if data.final_frame is not None:
        artifacts.append("closing geometry")
    if data.record is not None:
        artifacts.append("run record RUN.yaml")
    if data.orca_paths:
        artifacts.append(
            f"{len(data.orca_paths)} calculator output(s), last {data.orca_paths[-1].name}"
        )
    if (data.directory / "optimization.h5").is_file():
        artifacts.append("structured history optimization.h5 (not read: HDF5)")
    if artifacts:
        out.append("  artifacts: " + "; ".join(artifacts))
    if data.checks:
        out.append("")
        out.append("cross-checks:")
        for check in data.checks:
            out.append(f"  - {check.name}: {_CHECK_MARK[check.ok]} ({check.detail})")
    if data.orca_calls and data.orca_calls[-1].terminated is False:
        out.append(
            "  note: the last calculator ORCA output does not carry ORCA's normal "
            "termination banner; next step: inspect that output directly (menu 1)."
        )
    out += [
        "",
        "boundaries (the survey of the ORCA 6.1.1 pairing):",
        "  - gradient-driven work is exact through ORCA's .engrad sidecar (the checks "
        "above); a TS search is reliable with a model Hessian;",
        "  - frequency-dependent steps inside pysisyphus (hessian_init: calc, do_hess) "
        "fail on the '$multiplicity' block -- use a model Hessian and check "
        "frequencies outside pysisyphus;",
        "  - the structured history (optimization.h5) is HDF5 and is not read here "
        "(h5py is not a dependency); the text artifacts carry the same numbers;",
        "  - qm_calcs/cur_out is a symlink into the scratch directory, which "
        "pysisyphus cleans after every call, so it dangles after the run; the "
        "per-call copies (calculator_*.orca.*) are the durable ones.",
    ]
    return "\n".join(out)


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the interface survey and of the measured artifact formats."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "pysisyphus is a GPL-3.0 Python suite for potential-energy-surface "
                "exploration (optimisations, TS searches, IRC, chain-of-states) that "
                "drives external quantum-chemistry programs; the citation is the "
                "program's own (printed in every run's banner)."
            ),
            ref=(
                "Steinmetzer J., Kupfer S., Graefe S., Int. J. Quantum Chem. 2021, "
                "121(3), e26390, DOI 10.1002/qua.26390 (software: GPL-3.0, "
                "github.com/eljost/pysisyphus, release 1.0.0 of 2025-08-31)"
            ),
            url="https://doi.org/10.1002/qua.26390",
            bibkey="steinmetzer2021pysisyphus",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "Measured on pysisyphus 1.0.0 + ORCA 6.1.1 (fixtures/pysisyphus/, "
                "2026-09-28): the console/trj/geometry/RUN.yaml formats above; the "
                "gradient route is exact (H2O/def2-SVP optimisation end to end, and "
                "the Baker case-11 TS at HF/3-21G reaches -154.05045566 Eh, "
                "reproducing pysisyphus's own asserted -154.050455732882 to 7e-8); "
                "and the quantum-Hessian route crashes with ParseException at ORCA "
                "6's '$multiplicity' Hessian block (grammar identical in release "
                "1.0.0 and on master), so hessian_init: calc and do_hess are "
                "unusable in this pairing."
            ),
            ref="tests/test_pysisyphus.py; fixtures/pysisyphus/ (h2o_opt, butadiene_ts, h2o_stop3, hess_crash)",
        ),
    )
