"""the run artifacts of pysisyphus (formats; the reader behind menu 33).

pysisyphus (Steinmetzer, Kupfer, Graefe, Int. J. Quantum Chem. 2021, 121,
e26390; GPL-3.0) explores potential-energy surfaces by orchestrating external
calculators -- here ORCA -- and writes, next to the input YAML, a small set of
run artifacts.  Every convention below is **measured** on runs of pysisyphus
1.0.0 driving ORCA 6.1.1 (fixtures ``fixtures/pysisyphus/``, 2026-09-28):

- the **console capture** (the reader's primary source): the ``# RUNNING ... #``
  banner, the system/calculator/charge block, the convergence thresholds (each
  with its overachieved value in parentheses), the per-cycle table
  ``cycle  Δ(energy)  max(|force|)  rms(force)  max(|step|)  rms(step)  s/cycle``
  whose converged entries carry a trailing ``*``, the closing marker
  (``Converged!`` or ``Number of cycles exceeded!``), the ``Final summary``
  block (internal/cartesian force norms and the energy, both in hartree), the
  ``Wrote final, hopefully optimized, geometry to '<file>'`` line, and -- on a
  crash -- the traceback with the ``crashed_<name>`` backup path;
- ``optimization.trj`` (``ts_optimization.trj`` for TS runs): concatenated XYZ
  frames, one per cycle; the comment line carries the energy in the
  ``f"{energy: >20.8f} , "`` form (``Geometry.comment``; an unevaluated geometry
  carries an empty comment instead);
- ``final_geometry.xyz`` (``ts_final_geometry.xyz``; TS runs also copy it to
  ``ts_opt.xyz``): the closing geometry, same comment convention, **without** a
  trailing newline;
- ``RUN.yaml``: the normalised input record (all keys expanded, unused sections
  ``null``) plus a ``version`` key;
- ``qm_calcs/calculator_<i>.<cycle>.orca.{inp,out,engrad,...}``: every ORCA
  call of the run (the console points into this directory), readable by the
  ORCA parser; ``qm_calcs/cur_out`` is a symlink into the (already cleaned)
  scratch directory, so it does not survive the run.

This module is the format layer only: it turns those files into plain data.
The run-level semantics -- cross-checks between the artifacts, the crash
classification, the report -- live in ``analysis/pysisyphus_run.py``.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path

__all__ = [
    "PysisyphusError",
    "Frame",
    "CycleRow",
    "read_xyz",
    "read_trj",
    "read_log",
    "read_run_record",
    "find_run_log",
]


class PysisyphusError(ValueError):
    """A pysisyphus artifact is missing or of an unrecognised layout (with a next step)."""


@dataclass(frozen=True)
class Frame:
    """One geometry of a trajectory (or a standalone XYZ): energy and atoms.

    ``comment`` is the comment line with the energy prefix removed -- the
    measured layout is ``"<energy> , <rest>"`` and pysisyphus itself strips the
    energy prefix when it re-reads its own files, so the field is carried here
    for symmetry with that behaviour.
    """

    energy: float | None
    comment: str
    atoms: tuple[tuple[str, float, float, float], ...]


@dataclass(frozen=True)
class CycleRow:
    """One row of the console cycle table (``None`` where the entry was ``nan``)."""

    cycle: int
    d_energy: float | None
    max_force: float | None
    rms_force: float | None
    max_step: float | None
    rms_step: float | None
    seconds: float | None
    starred: frozenset[str]


#: the energy prefix of a frame comment: f"{energy: >20.8f} , " (Geometry.comment)
_COMMENT_ENERGY_RE = re.compile(r"^\s*([-+]?\d+\.\d+)\s*,(\s*)(.*)$")

_BANNER_RE = re.compile(r"#\s*RUNNING\s+(.+?)\s*#")
_VERSION_RE = re.compile(r"Version\s+(\S+)\s+\(Python\s+([^,]+),")
_EXECUTED_RE = re.compile(r"Executed at\s+(.+?)\s+on\s+'([^']+)'")
_SYSTEM_RE = re.compile(r"Input geometry:\s*Geometry\((.+?),\s*(\d+)\s+atoms?\)")
_COORD_RE = re.compile(r"Coordinate system:\s*(\S+)\s*\n\s*Coordinate number:\s*(\d+)")
_CALC_RE = re.compile(r"Calculator:\s*(\w+)\(([\w.-]+)\)")
_CHARGE_RE = re.compile(r"Charge:\s*(-?\d+)")
_MULT_RE = re.compile(r"Multiplicity:\s*(\d+)")
_OPTIMIZER_RE = re.compile(r"Optimizer:\s*(\S+)")
_THRESHOLD_RE = re.compile(
    r"(max|rms)\(\|?(force|step)\|?\)\s*<=\s*([0-9.]+)(?:,\s*\(([0-9.]+)\))?"
)
_FINAL_GEOM_RE = re.compile(
    r"Wrote final, hopefully optimized, geometry to '([^']+)'"
)
_DURATION_RE = re.compile(r"pysisyphus run took\s+([\d:]+)\s+h")
_CITATION_RE = re.compile(r"https://doi\.org/\S+")
_FINAL_SUMMARY_RE = re.compile(
    r"(max\(forces, internal\)|rms\(forces, internal\)|max\(forces,cartesian\)|"
    r"rms\(forces,cartesian\)|energy):\s*([-+\d.]+)"
)
_CRASHED_RE = re.compile(r"RunAfterCalculationFailedException")
_COPIED_RE = re.compile(r"Copied contents of\s*'([^']+)'\s*to\s*'([^']+)'")
_EXCEPTION_LINE_RE = re.compile(r"^\w[\w.]*(?:Error|Exception)\b.*$", re.MULTILINE)

_TABLE_HEADER_RE = re.compile(r"^\s*cycle\s+.*Δ\(energy\).*$")
_TABLE_ROW_FIELDS = 7  # cycle + six numeric columns


def _xyz_frames(text: str) -> tuple[Frame, ...]:
    """The (one or more) XYZ frames in a text, comment lines included.

    The measured files have no trailing newline on the last atom line, so the
    walk is line-based and does not rely on blank-line separators.
    """
    lines = text.splitlines()
    frames: list[Frame] = []
    index = 0
    while index < len(lines):
        line = lines[index].strip()
        if not line:
            index += 1
            continue
        try:
            natoms = int(line.split()[0])
        except ValueError:
            raise PysisyphusError(
                f"the line {lines[index]!r} is not an atom count; the file is not the "
                "XYZ/trajectory layout pysisyphus writes. Next step: give a "
                "pysisyphus artifact (.trj / *_final_geometry.xyz)."
            ) from None
        if index + 1 >= len(lines):
            raise PysisyphusError(
                "the file ends after an atom count (no comment/atoms). Next step: "
                "check whether the file was truncated."
            )
        comment_line = lines[index + 1]
        energy = None
        comment = comment_line.strip()
        if (match := _COMMENT_ENERGY_RE.match(comment_line)) is not None:
            energy = float(match.group(1))
            comment = match.group(3).strip()
        atoms: list[tuple[str, float, float, float]] = []
        for atom_line in lines[index + 2 : index + 2 + natoms]:
            fields = atom_line.split()
            if len(fields) < 4:
                raise PysisyphusError(
                    f"the coordinate line {atom_line!r} carries fewer than four fields."
                )
            atoms.append((fields[0], float(fields[1]), float(fields[2]), float(fields[3])))
        if len(atoms) != natoms:
            raise PysisyphusError(
                f"the frame declares {natoms} atoms but holds {len(atoms)} coordinate "
                "lines (the file looks truncated)."
            )
        frames.append(Frame(energy=energy, comment=comment, atoms=tuple(atoms)))
        index += 2 + natoms
    if not frames:
        raise PysisyphusError(
            "the file holds no XYZ frame. Next step: give a pysisyphus trajectory "
            "(.trj) or a final-geometry file."
        )
    return tuple(frames)


def read_xyz(text: str) -> Frame:
    """A single XYZ structure in the pysisyphus convention (with the energy comment)."""
    frames = _xyz_frames(text)
    if len(frames) != 1:
        raise PysisyphusError(
            f"the file holds {len(frames)} frames; a single structure was expected. "
            "Next step: use the trajectory reader for multi-frame files."
        )
    return frames[0]


def read_trj(text: str) -> tuple[Frame, ...]:
    """All frames of a pysisyphus trajectory, in file order (one per cycle)."""
    return _xyz_frames(text)


def _parse_cycles(lines: list[str]) -> tuple[CycleRow, ...]:
    """The cycle table: the rows between the header and the first non-row line.

    Long tables restart after every ten rows with a dashes separator and no
    repeated header (measured on the TS fixture: rows 0-9, dashes, rows 10-19),
    so a dashes line inside the table is skipped rather than treated as the end.
    """
    start = next(
        (index for index, line in enumerate(lines) if _TABLE_HEADER_RE.match(line)), None
    )
    if start is None:
        return ()
    rows: list[CycleRow] = []
    for line in lines[start + 1 :]:
        fields = line.split()
        if len(fields) == 1 and set(fields[0]) == {"-"}:
            continue  # the dashes lines (under the header and between blocks)
        if len(fields) != _TABLE_ROW_FIELDS or not fields[0].isdigit():
            if rows:
                break
            continue
        starred = set()
        values: list[float | None] = []
        names = ("d_energy", "max_force", "rms_force", "max_step", "rms_step", "seconds")
        for name, token in zip(names, fields[1:]):
            raw = token[:-1] if token.endswith("*") else token
            if token.endswith("*"):
                starred.add(name)
            try:
                value = float(raw)
            except ValueError:
                raise PysisyphusError(
                    f"the cycle-table entry {token!r} is not a number; the table layout "
                    "differs from the one this reader was measured on."
                ) from None
            values.append(None if math.isnan(value) else value)
        rows.append(CycleRow(cycle=int(fields[0]), starred=frozenset(starred), **dict(zip(names, values))))
    return tuple(rows)


def read_log(text: str) -> dict:
    """The facts of a pysisyphus console capture (one run per capture).

    The ``RUNNING`` banner is required -- it is what tells a pysisyphus capture
    apart from any other log.  The outcome is classified from the program's own
    closing markers (``Converged!`` / ``Number of cycles exceeded!``), never
    from the force values: the measured non-converged fixture stops with
    forces already below the thresholds, and the marker is the ground truth.
    """
    if (banner := _BANNER_RE.search(text)) is None:
        raise PysisyphusError(
            "the file carries no '# RUNNING ... #' banner, so it is not a pysisyphus "
            "console capture. Next step: give the capture of the run, e.g. the file "
            "you get from ``pysis main.yaml > pysis.out``."
        )
    lines = text.splitlines()
    facts: dict = {
        "version": None,
        "python": None,
        "executed_at": None,
        "host": None,
        "kind": banner.group(1).strip(),
        "system": None,
        "n_atoms": None,
        "coord_type": None,
        "n_coord": None,
        "calculator": None,
        "calculator_id": None,
        "charge": None,
        "multiplicity": None,
        "optimizer": None,
        "thresholds": {},
        "cycles": _parse_cycles(lines),
        "final": {},
        "final_geometry_fn": None,
        "duration": None,
        "citation": None,
        "crash": None,
    }
    if (match := _VERSION_RE.search(text)) is not None:
        facts["version"] = match.group(1)
        facts["python"] = match.group(2).strip()
    if (match := _EXECUTED_RE.search(text)) is not None:
        facts["executed_at"] = match.group(1)
        facts["host"] = match.group(2)
    if (match := _SYSTEM_RE.search(text)) is not None:
        facts["system"] = match.group(1)
        facts["n_atoms"] = int(match.group(2))
    if (match := _COORD_RE.search(text)) is not None:
        facts["coord_type"] = match.group(1)
        facts["n_coord"] = int(match.group(2))
    if (match := _CALC_RE.search(text)) is not None:
        facts["calculator"] = match.group(1)
        facts["calculator_id"] = match.group(2)
    if (match := _CHARGE_RE.search(text)) is not None:
        facts["charge"] = int(match.group(1))
    if (match := _MULT_RE.search(text)) is not None:
        facts["multiplicity"] = int(match.group(1))
    if (match := _OPTIMIZER_RE.search(text)) is not None:
        facts["optimizer"] = match.group(1)
    for match in _THRESHOLD_RE.finditer(text):
        kind = f"{match.group(1)}_{match.group(2)}"
        overachieved = float(match.group(4)) if match.group(4) else None
        facts["thresholds"][kind] = (float(match.group(3)), overachieved)
    summary_tail = text.split("Final summary:", 1)[-1]
    for match in _FINAL_SUMMARY_RE.finditer(summary_tail):
        label = match.group(1)
        key = {
            "max(forces, internal)": "max_forces_internal",
            "rms(forces, internal)": "rms_forces_internal",
            "max(forces,cartesian)": "max_forces_cartesian",
            "rms(forces,cartesian)": "rms_forces_cartesian",
            "energy": "energy",
        }[label]
        facts["final"][key] = float(match.group(2))
    if (match := _FINAL_GEOM_RE.search(text)) is not None:
        facts["final_geometry_fn"] = match.group(1)
    if (match := _DURATION_RE.search(text)) is not None:
        facts["duration"] = match.group(1)
    if (match := _CITATION_RE.search(text)) is not None:
        facts["citation"] = match.group(0)
    if _CRASHED_RE.search(text):
        crash: dict = {"backup": None, "source": None, "exceptions": ()}
        if (match := _COPIED_RE.search(text)) is not None:
            crash["source"] = match.group(1)
            crash["backup"] = match.group(2)
        crash["exceptions"] = tuple(
            line.strip() for line in _EXCEPTION_LINE_RE.findall(text)
        )
        crash["converged_before_crash"] = "Converged!" in text
        facts["crash"] = crash
        facts["outcome"] = "crashed"
    elif "Converged!" in text:
        facts["outcome"] = "converged"
    elif "Number of cycles exceeded!" in text:
        facts["outcome"] = "cycles exceeded"
    else:
        facts["outcome"] = "unknown"
    return facts


def read_run_record(text: str) -> dict:
    """The normalised run record (``RUN.yaml``) as a plain dict.

    The file is written by pysisyphus itself (the input with every key
    expanded, unused sections set to ``null``, plus ``version``), so the YAML
    is loaded as-is; the reader only requires the top level to be a mapping.
    """
    import yaml

    try:
        record = yaml.safe_load(text)
    except yaml.YAMLError as exc:
        raise PysisyphusError(
            f"RUN.yaml does not parse as YAML ({exc}). Next step: give the run "
            "directory of a pysisyphus run (the file is written by the program)."
        ) from exc
    if not isinstance(record, dict):
        raise PysisyphusError(
            "RUN.yaml's top level is not a mapping; the file is not the run record "
            "pysisyphus writes."
        )
    return record


def find_run_log(path: Path) -> Path:
    """The console capture: the given file, or the unique candidate in a directory.

    In a directory the candidates are the top-level ``*.log``/``*.out``/``*.txt``
    files whose content carries the ``RUNNING`` banner; the run's own logging
    files (``pysisyphus.log`` etc.) do not carry it, so the scan is unambiguous
    for a standard run directory.
    """
    if path.is_file():
        return path
    if not path.is_dir():
        raise PysisyphusError(
            f"{path} is neither a file nor a directory. Next step: check the path."
        )
    candidates = []
    for pattern in ("*.log", "*.out", "*.txt"):
        for candidate in sorted(path.glob(pattern)):
            try:
                head = candidate.read_text(encoding="utf-8", errors="replace")
            except OSError:
                continue
            if _BANNER_RE.search(head):
                candidates.append(candidate)
    if not candidates:
        raise PysisyphusError(
            f"no file under {path} carries the pysisyphus console banner. Next step: "
            "give the console capture itself (e.g. ``pysis main.yaml > pysis.out``) "
            "or the run directory that holds it."
        )
    if len(candidates) > 1:
        names = ", ".join(candidate.name for candidate in candidates)
        raise PysisyphusError(
            f"several files carry the console banner ({names}); the reader expects one "
            "run per capture. Next step: give the capture file directly."
        )
    return candidates[0]
