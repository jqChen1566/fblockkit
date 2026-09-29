"""Gaussian (G09/G16) output: the minimal fact set (plan item 6.4).

Termination, SCF, frequency and optimization facts mapped onto the
diagnosis engine's shared fact vocabulary (``knowledge/rules/README.md``);
CASSCF facts are deliberately not emitted yet -- no Gaussian CASSCF sample
is in hand and nothing is guessed (registered as a follow-up).

Fixture provenance (measured 2026-09-29): the two shipped probes are real
G09 Rev D.01 outputs (the gau_orca example set, 2018).  Both are
external-driver runs -- the energy came from an attached program, so they
carry **no ``SCF Done`` line**; they anchor the termination / frequency /
optimization facts and the absent-field semantics.  The ``SCF Done`` parse
follows the published line shape ("SCF Done:  E(RB3LYP) =  -... A.U. after
N cycles", stable across G09/G16); its fixture anchor is pending a
conventional-SCF sample (registered -- do not treat the regex as
fixture-verified).  The optimization-cycle marker is "Step number  N out of
a maximum of  M" (measured); frequency blocks are collected from the
"Frequencies --" lines (three values per line, continuation lines carry
their own prefix).

Menu 1 (the check-up report) consumes this through ``parse_auto``: the
identification banner is "Entering Gaussian System"; the diagnosis rules
keyed on shared facts (termination, SCF, frequencies, optimization) run
unchanged, while ORCA-specific fields stay absent and their rules do not
fire.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from ..knowledge.models import ParseResult
from .base import ParserError, read_text, register

PROGRAM = "gaussian"

_BANNER = "Entering Gaussian System"
_VERSION_RE = re.compile(r"Gaussian (\d+), Revision ([A-Z0-9.]+)")
_TERMINATION_RE = re.compile(r"Normal termination of Gaussian (\d+) at (.+)\.\s*$")
_ERROR_TERMINATION = "Error termination"
# published line shape; see the module docstring for the fixture caveat
_SCF_DONE_RE = re.compile(
    r"SCF Done:\s+E\((\w+)\)\s+=\s+(-?\d+\.\d+)\s+A\.U\. after\s+(\d+) cycles"
)
_SCF_FAILURE = "Convergence failure"
_FREQ_LINE_RE = re.compile(r"^\s*Frequencies --\s+(.*\S)\s*$")
_IMAG_MARK_RE = re.compile(r"(\d+) imaginary frequencies? \(negative Signs\)")
_BERNY_RE = re.compile(r"Berny optimization")
_STEP_RE = re.compile(r"Step number\s+(\d+)\s+out of a maximum of\s+(\d+)")
_STATIONARY_RE = re.compile(r"Stationary point found")
_OPT_COMPLETED_RE = re.compile(r"Optimization completed")
_INPUT_HEADER_RE = re.compile(r"^\s*#\s*(.+)$")


class GaussianParser:
    """The minimal Gaussian facts; heavier menus stay ORCA-only."""

    program = PROGRAM

    def detect(self, path: str | Path) -> bool:
        """Read the head of the file only (no full parse)."""
        text = read_text(path)
        head = "\n".join(text.splitlines()[:120])
        return _BANNER in head

    def parse(self, path: str | Path) -> ParseResult:
        text = read_text(path)
        lines = text.splitlines()
        if _BANNER not in "\n".join(lines[:120]):
            raise ParserError(
                f"{path} is not Gaussian output (no 'Entering Gaussian System' "
                "banner found)."
            )
        version = ""
        if (match := _VERSION_RE.search(text)) is not None:
            version = f"{match.group(1)}, Revision {match.group(2)}"
        sections: dict[str, Any] = {
            "version": version,
            "termination": _parse_termination(text),
            "gaussian_scf": _parse_scf(lines),
            "gaussian_frequencies": _parse_frequencies(lines),
            "gaussian_optimization": _parse_optimization(lines),
        }
        # the after-geometry question, same semantics as the ORCA parser:
        # with an optimisation, a frequency block printed before the last
        # cycle is not a verification of the optimised structure; without
        # one the question does not arise
        frequencies = sections["gaussian_frequencies"]
        optimization = sections["gaussian_optimization"]
        if frequencies["present"]:
            last_cycle = optimization["last_cycle_line"]
            frequencies["after_geometry"] = (
                last_cycle is None or frequencies["last_line"] > last_cycle
            )
        else:
            frequencies["after_geometry"] = None
        return ParseResult(program=PROGRAM, path=str(path), sections=sections)


def _parse_termination(text: str) -> dict[str, Any]:
    found = _TERMINATION_RE.findall(text)
    errors = text.count(_ERROR_TERMINATION)
    return {
        "normal": bool(found) and errors == 0,
        "line": found[-1][1].strip() if found else None,
        "error_terminations": errors,
    }


def _parse_scf(lines: list[str]) -> dict[str, Any]:
    """The last SCF Done line (published shape; fixture anchor registered)."""
    converged: bool | None = None
    energy: float | None = None
    cycles: int | None = None
    method: str | None = None
    for line in lines:
        if _SCF_FAILURE in line:
            converged = False
            continue
        if (match := _SCF_DONE_RE.search(line)) is not None:
            method = match.group(1)
            energy = float(match.group(2))
            cycles = int(match.group(3))
            converged = True
    return {
        "present": energy is not None,
        "converged": converged,
        "energy": energy,
        "cycles": cycles,
        "method": method,
    }


def _parse_frequencies(lines: list[str]) -> dict[str, Any]:
    frequencies: list[float] = []
    marks: list[int] = []
    last_line = 0
    for index, line in enumerate(lines, start=1):
        if (match := _FREQ_LINE_RE.match(line)) is not None:
            frequencies.extend(float(value) for value in match.group(1).split())
            last_line = index
            continue
        if (match := _IMAG_MARK_RE.search(line)) is not None:
            marks.append(int(match.group(1)))
    if not frequencies:
        return {"present": False, "blocks": [], "after_geometry": None}
    block = {
        "frequencies": frequencies,
        "n_imaginary": sum(1 for value in frequencies if value < 0),
        "engine_imaginary_marks": marks,
        "last_block_line": last_line,
    }
    return {
        "present": True,
        "blocks": [block],
        "last_line": last_line,
        "after_geometry": None,  # filled by the caller
    }


def _parse_optimization(lines: list[str]) -> dict[str, Any]:
    present = False
    converged = False
    steps = 0
    last_cycle_line: int | None = None
    ts = False
    header_seen = False
    for index, line in enumerate(lines, start=1):
        if not header_seen and (match := _INPUT_HEADER_RE.match(line)) is not None:
            header_seen = True
            if re.search(r"\bts\b", match.group(1), re.IGNORECASE):
                ts = True
            continue
        if _BERNY_RE.search(line):
            present = True
            last_cycle_line = index
            continue
        if (match := _STEP_RE.search(line)) is not None:
            present = True
            steps = int(match.group(1))
            last_cycle_line = index
            continue
        if _OPT_COMPLETED_RE.search(line) or _STATIONARY_RE.search(line):
            converged = True
            last_cycle_line = index
    return {
        "present": present,
        "converged": converged if present else None,
        "steps": steps,
        "ts": ts,
        "last_cycle_line": last_cycle_line,
    }


def facts(result: ParseResult) -> dict[str, Any]:
    """Map a Gaussian parse result onto the shared fact fields."""
    sections = result.sections
    termination = sections.get("termination") or {}
    scf = sections.get("gaussian_scf") or {}
    frequencies = sections.get("gaussian_frequencies") or {}
    optimization = sections.get("gaussian_optimization") or {}
    blocks = frequencies.get("blocks") or []
    last_block = blocks[-1] if blocks else None
    minimum = None
    if last_block is not None:
        imaginary = [value for value in last_block["frequencies"] if value < 0]
        minimum = min(imaginary) if imaginary else None
    return {
        "terminated_normally": termination.get("normal"),
        "scf_converged": scf.get("converged"),
        "scf_cycles": scf.get("cycles"),
        "ts_optimization": optimization.get("ts", False),
        "optimization_converged": optimization.get("converged"),
        "frequency_present": frequencies.get("present"),
        "frequency_after_geometry": frequencies.get("after_geometry"),
        "frequency_imaginary_count": (
            last_block["n_imaginary"] if last_block is not None else None
        ),
        "frequency_min_imaginary": minimum,
    }


register(GaussianParser())
