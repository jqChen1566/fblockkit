"""ORCA output parser (architecture design v0.1 §3 L1).

Coverage (v0.1, deep): program version, normal-termination marker, SCF
convergence block (cycle count, energy history, and -- since v0.2 -- the
header-driven DIIS/SOSCF iteration tables and the SCF convergence summary with
the ConvCheckMode), orbital energy/occupation table, CASSCF convergence and
states, NEVPT2/CASPT2 results (reference weight and smallest energy denominator
-- the fact source for the D3 criteria), SOC markers, the coupled-cluster T1
diagnostic, vibrational frequencies (imaginary-mode count), the local spin
analysis block (fragment expectation values), geometry optimisation and
transition-state context, errors and warnings.

Every regular expression was checked against real ORCA 6.1.1 output (fixtures in
``fixtures/orca/``; their provenance and generation are recorded in
``fixtures/orca/README.md``). The measured fixtures, not manual examples, define
the format -- manual examples can come from older versions (the CASPT2 reference
weight line style has changed, for instance).

Conventions:

- ``scf.energy`` is the ``Total Energy`` of the ``TOTAL SCF ENERGY`` block;
  ``final_energy`` is the closing ``FINAL SINGLE POINT ENERGY`` (the final value
  after post-DFT processing; the two agree for a plain HF/DFT single point and
  the latter governs when post-processing is present).
- ``orbitals`` takes the **last** ``ORBITAL ENERGIES`` table: in a multi-method
  chain the last table reflects the final orbitals (e.g. CASSCF natural-orbital
  occupations) rather than the initial SCF.
- The CASSCF ``energy`` is the ``Final CASSCF energy`` (state-averaged
  convention); individual state energies are listed separately in ``states``
  (state-specific convention). The two conventions must not be mixed (see the
  diagnosis-layer criteria).
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any

from ..knowledge.elements import ElementError, is_f_element
from ..knowledge.models import ParseResult
from .base import ParserError, float_or_none, read_text, register

PROGRAM = "orca"

# --- identification and global markers --------------------------------------

_BANNER_RE = re.compile(r"\*\s+O\s+R\s+C\s+A\s+\*")
_VERSION_RE = re.compile(r"Program Version\s+(\d+(?:\.\d+)*)")
TERMINATED_MARK = "****ORCA TERMINATED NORMALLY****"

# --- SCF --------------------------------------------------------------------

_SCF_CONVERGED_RE = re.compile(r"SCF CONVERGED AFTER\s+(\d+)\s+CYCLES")
_SCF_NOT_CONVERGED_RE = re.compile(r"SCF NOT CONVERGED AFTER\s+(\d+)\s+CYCLES")
# The DIIS and SOSCF iteration tables share the first two columns (index, energy),
# while the third (Delta-E/RMSDP) is in scientific notation -- that is the
# fingerprint, so density-matrix lines or basis-set contraction tables are not
# matched by mistake.
_SCF_ITER_RE = re.compile(
    r"^\s{0,10}(\d{1,5})\s+(-?\d+\.\d{6,})\s+[-+]?\d+\.\d+e[+-]\d{2}\s"
)
_TOTAL_SCF_RE = re.compile(r"^Total Energy\s+:\s+(-?\d+\.\d+)\s+Eh", re.MULTILINE)
_FINAL_ENERGY_RE = re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)")

# --- ORBITAL ENERGIES table -------------------------------------------------

_ORB_HEADER = "NO   OCC          E(Eh)"
_ORB_ROW_RE = re.compile(
    r"^\s*(\d+)\s+(\d+\.\d{4})\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$"
)

# --- CASSCF -----------------------------------------------------------------

_CASSCF_BANNER = "ORCA-CASSCF"
_MACRO_ITER_RE = re.compile(r"^MACRO-ITERATION\s+(\d+):")
_CASSCF_ITER_RE = re.compile(r"E\(CAS\)=\s+(-?\d+\.\d+)\s+Eh")
# CASSCF convergence has two marker kinds (both measured in the fixtures):
#   ---- THE CAS-SCF ENERGY   HAS CONVERGED ----   (energy criterion)
#   ---- THE CAS-SCF GRADIENT HAS CONVERGED ----   (gradient criterion)
# They do not mean the same thing (our group's lesson: the convergence criterion
# is the gradient, not the energy), hence converged_via is recorded separately.
_CASSCF_CONVERGED_RE = re.compile(r"----\s*THE CAS-SCF (ENERGY|GRADIENT)\s+HAS CONVERGED")
# Active-orbital natural occupations (printed once per macro iteration; take the last)
_N_OCC_RE = re.compile(r"N\(occ\)=\s+(.*\S)\s*$")
_FINAL_CASSCF_RE = re.compile(r"^Final CASSCF energy\s+:\s+(-?\d+\.\d+)\s+Eh", re.MULTILINE)
_BLOCK_RE = re.compile(r"CAS-SCF STATES FOR BLOCK\s+(\d+)\s+MULT=\s+(\d+)\s+NROOTS=\s+(\d+)")
_ROOT_RE = re.compile(
    r"^ROOT\s+(\d+):\s+E=\s+(-?\d+\.\d+)\s+Eh"
    r"(?:\s+(-?\d+\.\d+)\s+eV\s+(-?\d+\.\d+)\s+cm\*\*-1)?"
)

# --- perturbation layer (NEVPT2 / CASPT2 share the skeleton) ----------------

_PT2_SECTION_RE = {
    "nevpt2": re.compile(r"<\s*NEVPT2\s*>"),
    "caspt2": re.compile(r"<\s*CASPT2\s*>"),
}
_PT2_END_RE = {
    "nevpt2": re.compile(r"TIMINGS NEVPT2"),
    "caspt2": re.compile(r"TIMINGS CASPT2"),
}
_PT2_STATE_RE = re.compile(r"^\s+MULT\s+(\d+),\s+ROOT\s+(\d+)\s*$")
_PT2_DE_RE = re.compile(r"Total Energy Correction\s+:\s+dE\s+=\s+(-?\d+\.\d+)")
_PT2_E0_RE = re.compile(r"Reference\s+Energy\s+:\s+E0\s+=\s+(-?\d+\.\d+)")
_PT2_W0_RE = re.compile(r"Reference\s+Weight\s+:\s+W0\s+=\s+(-?\d+\.\d+)")
# Reference weight printed on its own during the iteration phase (no W0 marker;
# appears after each state converges in CASPT2)
_PT2_WEIGHT_ITER_RE = re.compile(r"^\s*Reference Weight\s+=\s+(-?\d+\.\d+)\s*$")
_PT2_ETOT_RE = re.compile(r"Total Energy \(E0\+dE\)\s+:\s+E\s+=\s+(-?\d+\.\d+)")
_SMALLEST_DENOM_RE = re.compile(r"smallest energy denominator\s+(\w+)\s+=\s+(-?\d+\.\d+)")
# Per-class correlation contributions (V0_ijab/IJAB, V1_i/ITUV, Vm1_a/TUVA, ...;
# two layouts exist, in the computing section and in the results section)
_PT2_CLASS_RE = re.compile(r"Class\s+(\w+)\s*:\s*dE\s*=\s*(-?\d+\.\d+)")

# --- per-MO composition table (Loewdin ORBITAL-COMPOSITIONS) ----------------
# Requires %output Print[P_ReducedOrbPopMO_L] 1 (manual §2.9, Print keyword);
# the table body is a repeating block: MO index row / energy row / occupation row
# / per-AO percentage rows.
_ORBCOMP_MARKER = "LOEWDIN ORBITAL-COMPOSITIONS"
_ORBCOMP_SECOND = "LOEWDIN REDUCED ACTIVE MOs"  # second table duplicates the first; stop there
_ORBCOMP_HEADER_RE = re.compile(r"^\s*\d+(?:\s+\d+)+\s*$")
_ORBCOMP_AO_ROW_RE = re.compile(r"^\s*(\d+)\s+([A-Z][a-z]?)\s+(\S+)\s+(.*?)\s*$")
_NUMBER_RE = re.compile(r"-?\d+\.\d+")


def _floats(text: str) -> list[float]:
    return [float(token) for token in _NUMBER_RE.findall(text)]


def _parse_orbital_composition(lines: list[str]) -> dict[str, Any]:
    """Parse the per-MO composition table; present=False when it was not printed
    (not an error -- it is an optional print)."""
    start = None
    end = len(lines)
    for i, line in enumerate(lines):
        if start is None and _ORBCOMP_MARKER in line:
            start = i + 1
        elif start is not None and _ORBCOMP_SECOND in line:
            end = i
            break
    if start is None:
        return {"present": False, "orbitals": ()}

    orbitals: list[dict[str, Any]] = []
    state = "seek"
    indices: list[int] = []
    entries: list[dict[str, Any]] = []
    weights: list[dict[tuple[int, str, str], float]] = []

    def _flush() -> None:
        for entry, shell_map in zip(entries, weights):
            entry["shells"] = tuple(
                {"atom": atom, "element": element, "shell": shell, "weight": weight}
                for (atom, element, shell), weight in sorted(
                    shell_map.items(), key=lambda item: item[1], reverse=True
                )
            )
            orbitals.append(entry)

    for line in lines[start:end]:
        stripped = line.strip()
        if stripped and set(stripped) <= {"-", " "}:
            continue  # divider line (dividers inside a block contain inner spaces:
                      # '--------  --------')
        if state == "seek":
            if _ORBCOMP_HEADER_RE.match(line):
                indices = [int(token) for token in stripped.split()]
                entries = [
                    {"index": index, "energy": None, "occupation": None, "shells": ()}
                    for index in indices
                ]
                weights = [dict() for _ in indices]
                state = "energy"
            continue
        if state == "energy":
            values = _floats(stripped)
            if len(values) == len(indices):
                for entry, value in zip(entries, values):
                    entry["energy"] = value
                state = "occupation"
            continue
        if state == "occupation":
            values = _floats(stripped)
            if len(values) == len(indices):
                for entry, value in zip(entries, values):
                    entry["occupation"] = value
                state = "rows"
            continue
        match = _ORBCOMP_AO_ROW_RE.match(line) if stripped else None
        if match is not None:
            atom = int(match.group(1))
            element = match.group(2)
            shell = match.group(3)[0]  # s / p / d / f / g (first char of the AO label)
            values = _floats(match.group(4))
            for shell_map, value in zip(weights, values):
                key = (atom, element, shell)
                shell_map[key] = shell_map.get(key, 0.0) + value
            continue
        _flush()
        entries, weights, indices = [], [], []
        state = "seek"
    if state == "rows":
        _flush()
    return {"present": True, "orbitals": tuple(orbitals)}


# --- SOC markers ------------------------------------------------------------
# How the marker set was fixed (measured on the fixtures, see
# fixtures/orca/README.md):
# - a bare "SOC" is unusable: the author list at the top of an ORCA output
#   contains "Local ZFS, SOC";
# - "QDPT", "SOC integrals" and "SOC RMEs" are unusable: they also appear in
#   NEVPT2/CASPT2 output (the 'MRCI SOC BLOCK INPUT' comment and the timings
#   lines of the perturbation module), which would give false positives;
# - the markers below appear only when SOC is really being computed (checked one
#   by one against all five fixtures).
_SOC_MARKERS_RE = re.compile(
    r"Calculating SOCInts|Doing QDPT with ONLY SOC|NONZERO SOC MATRIX ELEMENTS|"
    r"SOC CORRECTED|SOC MATRIX|Lowest eigenvalue of the SOC matrix"
)

# --- coupled cluster diagnostics (T1) ---------------------------------------
# ORCA prints the T1 diagnostic (Lee & Taylor) inside the "COUPLED CLUSTER
# ENERGY" block, together with the singles norm it is derived from; no D1
# diagnostic is printed (measured on the CCSD fixtures, 2026-09-26).
_CC_SECTION_RE = re.compile(r"COUPLED CLUSTER ENERGY")
_T1_RE = re.compile(r"^T1 diagnostic\s+\.\.\.\s+(-?\d+\.\d+)")
_SINGLES_NORM_RE = re.compile(r"^Singles Norm <S\|S>\*\*1/2\s+\.\.\.\s+(-?\d+\.\d+)")


def _parse_cc(lines: list[str]) -> dict[str, Any]:
    """The coupled-cluster energy block: the T1 diagnostic and the singles norm.

    Last occurrence wins (a job may chain several CC treatments).
    """
    present = False
    t1: float | None = None
    singles_norm: float | None = None
    for line in lines:
        if _CC_SECTION_RE.search(line):
            present = True
        elif (m := _T1_RE.match(line)) is not None:
            t1 = float(m.group(1))
        elif (m := _SINGLES_NORM_RE.match(line)) is not None:
            singles_norm = float(m.group(1))
    return {"present": present or t1 is not None, "t1": t1, "singles_norm": singles_norm}


# --- vibrational frequencies -------------------------------------------------
# Measured format (fixtures: h2o_freq_min, nh3_planar_freq, h2o_linear_freq):
# a "VIBRATIONAL FREQUENCIES" header, the scaling-factor line, then one row per
# mode "     6:    -724.56 cm**-1  ***imaginary mode***" (imaginary modes are
# negative and carry the marker), the block ends at "NORMAL MODES".  A job may
# print more than one block (e.g. an initial Hessian before a TS optimisation
# and the final verification afterwards) -- the last block describes the final
# structure, which is the convention used here.
_FREQ_HEADER = "VIBRATIONAL FREQUENCIES"
_FREQ_END = "NORMAL MODES"
_FREQ_SCALING_RE = re.compile(r"Scaling factor for frequencies\s*=\s*(-?\d+\.\d+)")
_FREQ_ROW_RE = re.compile(r"^\s*(\d+):\s+(-?\d+\.\d+)\s+cm\*\*-1(\s+\*\*\*imaginary mode\*\*\*)?\s*$")


def _parse_frequency_block(lines: list[str], start: int) -> dict[str, Any]:
    scaling: float | None = None
    modes: list[dict[str, Any]] = []
    for line in lines[start + 1:]:
        if line.strip() == _FREQ_END:
            break
        if (m := _FREQ_SCALING_RE.search(line)) is not None:
            scaling = float(m.group(1))
            continue
        if (m := _FREQ_ROW_RE.match(line)) is not None:
            value = float(m.group(2))
            modes.append(
                {
                    "index": int(m.group(1)),
                    "wavenumber": value,
                    "imaginary": bool(m.group(3)) or value < 0.0,
                }
            )
        elif modes and line.strip():
            break  # left the block body
    imaginary = tuple(m for m in modes if m["imaginary"])
    return {
        "scaling": scaling,
        "modes": tuple(modes),
        "n_modes": len(modes),
        "imaginary": imaginary,
        "n_imaginary": len(imaginary),
        "min_wavenumber": min((m["wavenumber"] for m in modes), default=None),
        "min_imaginary": min((m["wavenumber"] for m in imaginary), default=None),
        "start_line": start,
    }


def _parse_frequencies(lines: list[str]) -> dict[str, Any]:
    starts = [i for i, line in enumerate(lines) if line.strip() == _FREQ_HEADER]
    if not starts:
        return {"present": False, "blocks": (), "last_block_line": None}
    blocks = tuple(_parse_frequency_block(lines, start) for start in starts)
    return {"present": True, "blocks": blocks, "last_block_line": starts[-1]}


# --- local spin analysis -----------------------------------------------------
# Measured format (fixtures: n2_stretch_local_spin, n2_stretch_casscf_local_spin;
# ORCA 6.1 manual §5.1.10): the block appears when the input divides the
# molecule into fragments (element symbol with a parenthesised fragment index,
# e.g. "N(1)"); it prints "<SA*SB>" (fragment-pair expectation values, lower
# triangle) and a "<SzA>  Seff(A)" table ("n.a." for <SzA> when the total state
# is a singlet).  A CASSCF job prints one block per root plus a state-average
# block; the UHF job prints several blocks (initial guess and final density).
_LOCAL_SPIN_MARKER = "LOCAL SPIN ANALYSIS"
_LS_FRAGS_RE = re.compile(r"^Number of fragments\s+=\s+(\d+)")
_LS_ATOMS_RE = re.compile(r"^Number of atoms\s+=\s+(\d+)")
_LS_BASIS_RE = re.compile(r"^Number of basis functions\s+=\s+(\d+)")
_LS_STATE_RE = re.compile(r"^State to be analyzed\s+=\s+(.+?)\s*$")
_LS_BLOCK_RE = re.compile(r"^State belongs to block\s+=\s+(.+?)\s*$")
_LS_MULT_RE = re.compile(r"^Multiplicity\s+=\s+(\d+)")
_LS_SAB_HEADER = "<SA*SB>"
_LS_SZ_HEADER = "<SzA>"
_LS_ROW_RE = re.compile(r"^\s*(\d+)\s*:\s+(.*?)\s*$")
_LS_NUMBER_RE = re.compile(r"-?\d+\.\d+")
# atom lines in the echoed input, carrying the fragment index: "|  7> N(1) 0 0 0"
_ECHO_ATOM_RE = re.compile(r"^\|\s*\d+>\s*([A-Z][a-z]?)\((\d+)\)")


def _ls_rows(chunk: list[str], header: str) -> dict[int, str]:
    """The ``index: values`` rows that follow a table header inside one block."""
    index = next((i for i, line in enumerate(chunk) if header in line), None)
    if index is None:
        return {}
    rows: dict[int, str] = {}
    for line in chunk[index + 1:]:
        stripped = line.strip()
        if not stripped or set(stripped) <= {"-"}:
            if rows:
                break  # a blank line after the rows ends the table
            continue
        if (m := _LS_ROW_RE.match(line)) is None:
            if rows:
                break
            continue
        rows[int(m.group(1))] = m.group(2)
    return rows


def _parse_local_spin_block(chunk: list[str]) -> dict[str, Any] | None:
    fields: dict[str, Any] = {
        "n_fragments": None,
        "n_atoms": None,
        "n_basis_functions": None,
        "state": None,
        "state_label": "",
        "block_label": "",
        "multiplicity": None,
        "sab": (),
        "sz": (),
        "sz_na": False,
        "seff": (),
    }
    for line in chunk:
        stripped = line.strip()
        if (m := _LS_FRAGS_RE.match(stripped)) is not None:
            fields["n_fragments"] = int(m.group(1))
        elif (m := _LS_ATOMS_RE.match(stripped)) is not None:
            fields["n_atoms"] = int(m.group(1))
        elif (m := _LS_BASIS_RE.match(stripped)) is not None:
            fields["n_basis_functions"] = int(m.group(1))
        elif (m := _LS_STATE_RE.match(stripped)) is not None:
            label = m.group(1)
            fields["state_label"] = label
            if label.isdigit():
                fields["state"] = int(label)
        elif (m := _LS_BLOCK_RE.match(stripped)) is not None:
            fields["block_label"] = m.group(1)
        elif (m := _LS_MULT_RE.match(stripped)) is not None:
            fields["multiplicity"] = int(m.group(1))

    sab_rows = _ls_rows(chunk, _LS_SAB_HEADER)
    if not sab_rows:
        return None
    size = max(sab_rows)
    matrix = [[0.0] * size for _ in range(size)]
    for row_index, text in sab_rows.items():
        values = [float(token) for token in _LS_NUMBER_RE.findall(text)]
        for column, value in enumerate(values, start=1):
            matrix[row_index - 1][column - 1] = value
            matrix[column - 1][row_index - 1] = value
    fields["sab"] = tuple(tuple(row) for row in matrix)

    sz: list[float | None] = []
    seff: list[float] = []
    for _, text in sorted(_ls_rows(chunk, _LS_SZ_HEADER).items()):
        if "n.a." in text:
            sz.append(None)
        else:
            values = _LS_NUMBER_RE.findall(text)
            sz.append(float(values[0]) if values else None)
        numbers = _LS_NUMBER_RE.findall(text)
        if numbers:
            seff.append(float(numbers[-1]))
    fields["sz"] = tuple(sz)
    fields["seff"] = tuple(seff)
    fields["sz_na"] = bool(sz) and all(value is None for value in sz)
    if fields["n_fragments"] is None:
        fields["n_fragments"] = size
    return fields


def _parse_local_spin(lines: list[str]) -> dict[str, Any]:
    starts = [i for i, line in enumerate(lines) if _LOCAL_SPIN_MARKER in line]
    if not starts:
        return {"present": False, "blocks": (), "fragment_elements": ()}
    blocks = []
    for position, start in enumerate(starts):
        end = starts[position + 1] if position + 1 < len(starts) else len(lines)
        block = _parse_local_spin_block(lines[start:end])
        if block is not None:
            blocks.append(block)
    fragments: dict[int, list[str]] = {}
    for line in lines:
        if (m := _ECHO_ATOM_RE.match(line)) is not None:
            fragments.setdefault(int(m.group(2)), []).append(m.group(1))
    size = max(fragments, default=0)
    fragment_elements = tuple(tuple(fragments.get(index, ())) for index in range(1, size + 1))
    return {"present": True, "blocks": tuple(blocks), "fragment_elements": fragment_elements}


# --- geometry optimisation / transition-state context ------------------------
# "Following TS mode number" appears only in transition-state optimisations
# (measured on fhh_optts_nofreq); the echoed input keyword lines are the
# fallback (comments in the echo start with "#", so only lines beginning with
# "!" count).  The optimisation status text is "THE OPTIMIZATION HAS CONVERGED"
# (h2o_freq_min) or "The optimization has not yet converged" per cycle; the last
# occurrence wins.
_TS_MODE_RE = re.compile(r"Following TS mode number")
_ECHO_KEYWORD_RE = re.compile(r"^\|\s*\d+>\s*!(.*)$")
_GEOM_CYCLE_RE = re.compile(r"GEOMETRY OPTIMIZATION CYCLE\s+\d+")
_OPT_CONVERGED_RE = re.compile(r"THE OPTIMIZATION HAS CONVERGED", re.IGNORECASE)
_OPT_NOT_CONVERGED_RE = re.compile(r"THE OPTIMIZATION HAS NOT YET CONVERGED", re.IGNORECASE)


def _parse_optimization(lines: list[str]) -> dict[str, Any]:
    cycle_lines = [i for i, line in enumerate(lines) if _GEOM_CYCLE_RE.search(line)]
    ts = any(_TS_MODE_RE.search(line) for line in lines)
    if not ts:
        for line in lines:
            if (m := _ECHO_KEYWORD_RE.match(line)) is not None and re.search(
                r"\bOptTS\b", m.group(1)
            ):
                ts = True
                break
    converged: bool | None = None
    for line in lines:
        if _OPT_CONVERGED_RE.search(line) is not None:
            converged = True
        elif _OPT_NOT_CONVERGED_RE.search(line) is not None:
            converged = False
    return {
        "ts": ts,
        "has_cycles": bool(cycle_lines),
        "converged": converged if cycle_lines else None,
        "last_cycle_line": cycle_lines[-1] if cycle_lines else None,
    }


# --- errors and warnings ----------------------------------------------------

_ERROR_MARKERS = (
    # case-insensitive: measured on the Yb fixtures, ORCA writes "Aborting the run"
    # where the LEANSCF abort writes "aborting the run"
    re.compile(r"ORCA finished by error termination", re.IGNORECASE),
    re.compile(r"mpirun noticed that process rank", re.IGNORECASE),
    re.compile(r"Segmentation fault", re.IGNORECASE),
    re.compile(r"aborting the run", re.IGNORECASE),
    re.compile(r"Error \(ORCA", re.IGNORECASE),
    # measured on generated_yb3_sarc2*.out: the unconverged-wavefunction abort and the
    # out-of-memory abort (both are explicit ORCA abort reasons worth reporting)
    re.compile(r"IS NOT FULLY CONVERGED", re.IGNORECASE),
    re.compile(r"OUT OF MEMORY", re.IGNORECASE),
)
_WARNING_LINE_RE = re.compile(r"^\s*(?:WARNING|Warning)\s*:")


def _slice(text: str, start: re.Pattern[str], end: re.Pattern[str]) -> str:
    """Text between the first hit of start and the first hit of end after it; with
    no end hit, up to the end of the text."""
    m = start.search(text)
    if not m:
        return ""
    tail = text[m.end():]
    m2 = end.search(tail)
    return tail[: m2.start()] if m2 else tail


def _parse_scf(lines: list[str]) -> dict[str, Any]:
    converged: bool | None = None
    cycles: int | None = None
    energies: list[float] = []
    for line in lines:
        if (m := _SCF_CONVERGED_RE.search(line)) is not None:
            converged, cycles = True, int(m.group(1))
        elif (m := _SCF_NOT_CONVERGED_RE.search(line)) is not None:
            converged, cycles = False, int(m.group(1))
        elif (m := _SCF_ITER_RE.match(line)) is not None:
            value = float_or_none(m.group(2))
            if value is not None:
                energies.append(value)
    text = "\n".join(lines)
    energy = None
    if (m := _TOTAL_SCF_RE.search(text)) is not None:
        energy = float_or_none(m.group(1))
    return {
        "converged": converged,
        "cycles": cycles,
        "energy": energy,
        "energies": tuple(energies),
        # the header-driven tables and the convergence summary (see above); read here
        # so that the diagnosis layer works on the parse result alone
        **_scan_scf_tables(lines),
    }


# --- SCF convergence tables (the reader lives here since v0.2) ---------------
# The DIIS/SOSCF iteration table and the "SCF CONVERGENCE" summary are parsed here so
# that the diagnosis layer can triage an SCF from the parse result alone, without
# re-reading the file.  Both readers are header-driven and were developed against
# fixtures/orca/ (see the fixture table): a row is attributed to the table whose header
# most recently preceded it, and it is accepted only when its first field continues the
# cycle count -- ORCA prints many other numeric tables, and the TRAH table ("Iter.") is
# not an "Iteration" table at all.
_ITER_TABLE_RE = re.compile(r"^\s*Iteration\s+Energy \(Eh\)")
_TABLE_DIVIDER_CHARS = frozenset("- ")
_TURN_ON_DIIS = "***Turning on AO-DIIS***"
_DIIS_RESET = "****Resetting DIIS****"
_CONVERGER_SWITCH_RE = re.compile(r"\*\*\*\s*Initializing\s+(\w+)\s*\*\*\*")
#: AutoTRAH announces itself in prose rather than with an "Initializing" banner.
_AUTO_TRAH_MARKER = "Leaving SCF to start the TRAH-SCF procedure"
_SCF_SOLVER_BANNER = "ORCA LEAN-SCF"
_CRITERIA_TITLE = "SCF CONVERGENCE"
_CRITERION_RE = re.compile(
    r"^\s*Last\s+(?P<name>.+?)\s*\.\.\.\s*"
    r"(?P<value>[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\s+"
    r"Tolerance\s*:\s*(?P<tolerance>[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\s*$"
)
#: "Convergence Check Mode ConvCheckMode   .... Total+1el-Energy"
_MODE_LABEL_RE = re.compile(r"Convergence Check Mode\s+ConvCheckMode\s*\.*\s*(?P<label>\S.*?)\s*$")
#: ORCA echoes the input file at the top of the output; explicit settings and the
#: !ExtremeSCF keyword can be seen there.
_ECHO_LINE_RE = re.compile(r"^\s*\|\s*\d+>")
_MODE_SETTING_RE = re.compile(r"ConvCheckMode\s+(\d)", re.IGNORECASE)
_EXTREME_KEYWORD_RE = re.compile(r"\bExtremeSCF\b", re.IGNORECASE)

#: ConvCheckMode values (ORCA 6.1 manual §2.6.1): 0 = every criterion has to be
#: satisfied, 1 = one criterion is enough, 2 = the total-energy and one-electron-energy
#: changes.  Mode 2 is ORCA's default and is what every standard preset sets, so it is
#: the assumption when an output says nothing about the mode.
CONVCHECK_ALL = 0
CONVCHECK_ONE_IS_ENOUGH = 1
CONVCHECK_ENERGY = 2
DEFAULT_CONVCHECK_MODE = CONVCHECK_ENERGY

#: Printed mode labels.  ORCA does not document these labels in the manual; the mapping
#: is measured on the fixtures (every output there that prints the SCF settings shows it)
#: and it agrees with the manual's wording for mode 2.  An unrecognized label falls back
#: to the default mode, and the consumer is told so through ``check_mode_source``.
_CONVCHECK_LABELS = {"total+1el-energy": CONVCHECK_ENERGY}


def _table_columns(header: str) -> tuple[str, ...]:
    """Column names of an SCF iteration table.

    A unit token such as ``(Eh)`` continues the previous name, so the number of names
    equals the number of data fields in a row.
    """
    names: list[str] = []
    for token in header.split():
        if token.startswith("(") and names:
            names[-1] = f"{names[-1]} {token}"
        else:
            names.append(token)
    return tuple(names)


def _resolve_check_mode(echo_text: str, mode_label: str) -> tuple[int, str]:
    """Which ConvCheckMode the run used, and where that was read from.

    An explicit setting in the echoed input wins (it is what ORCA was told to do); the
    printed label comes next; otherwise the documented default is assumed and the source
    string says so, so a consumer never claims more than it knows.
    """
    setting = _MODE_SETTING_RE.search(echo_text)
    if setting is not None and setting.group(1) in "012":
        return int(setting.group(1)), f"the echoed input setting 'ConvCheckMode {setting.group(1)}'"
    if _EXTREME_KEYWORD_RE.search(echo_text) is not None:
        return CONVCHECK_ALL, "the echoed input keyword '!ExtremeSCF' (mode 0)"
    if mode_label:
        known = _CONVCHECK_LABELS.get(mode_label.lower())
        if known is not None:
            return known, f"the printed mode label {mode_label!r}"
        return DEFAULT_CONVCHECK_MODE, (
            f"assumed default: the printed mode label {mode_label!r} is not one we know"
        )
    return DEFAULT_CONVCHECK_MODE, "assumed default (the output does not print the mode)"


def _scan_scf_tables(lines: list[str]) -> dict[str, Any]:
    """Read the SCF iteration tables and the convergence summary.

    Returns the measured signals (everything optional): whether the LEAN-SCF solver
    banner was seen, the DIIS-error history, the error at the AO-DIIS switch, DIIS
    resets, converger switches (SOSCF/TRAH), the largest absolute energy step, the
    printed convergence criteria with their tolerances, and the ConvCheckMode with the
    source it was read from.
    """
    solver_seen = False
    columns: tuple[str, ...] = ()
    diis_rows: list[tuple[int, float]] = []
    energy_steps: list[float] = []
    diis_at_switch: float | None = None
    diis_switch_cycle: int | None = None
    resets = 0
    switches: list[str] = []
    criteria: list[tuple[str, float, float]] = []
    expected_index = 0
    in_criteria = False
    mode_label = ""
    echo: list[str] = []

    for line in lines:
        stripped = line.strip()
        if _ECHO_LINE_RE.match(line) is not None:
            echo.append(line.split(">", 1)[1])
        if (label := _MODE_LABEL_RE.search(line)) is not None:
            mode_label = label.group("label")
        if in_criteria:
            match = _CRITERION_RE.match(line)
            if match is not None:
                value = float_or_none(match.group("value"))
                tolerance = float_or_none(match.group("tolerance"))
                if value is not None and tolerance is not None:
                    criteria.append((match.group("name").strip(), value, tolerance))
                continue
            if not criteria and (not stripped or set(stripped) <= _TABLE_DIVIDER_CHARS):
                continue  # divider and blank line between the title and the table
            in_criteria = False  # block ended; fall through to the other readers

        if _SCF_SOLVER_BANNER in line:
            solver_seen = True
        if _ITER_TABLE_RE.match(line) is not None:
            columns = _table_columns(line)
            expected_index = 0
            continue
        if columns and stripped:
            fields = line.split()
            try:
                index = int(fields[0])
            except ValueError:
                index = None
            values = [float_or_none(token) for token in fields[1:]]
            if (
                index is not None
                and len(fields) == len(columns)
                and index == expected_index + 1
                and all(value is not None and math.isfinite(value) for value in values)
            ):
                expected_index = index
                if "DIISErr" in columns:
                    diis_rows.append((index, values[columns.index("DIISErr") - 1]))
                    if "Delta-E" in columns:
                        energy_steps.append(abs(values[columns.index("Delta-E") - 1]))
                continue

        if _TURN_ON_DIIS in line and diis_at_switch is None and diis_rows:
            diis_at_switch = diis_rows[-1][1]
            diis_switch_cycle = diis_rows[-1][0]
        if _DIIS_RESET in line:
            resets += 1
        if (match := _CONVERGER_SWITCH_RE.search(line)) is not None:
            switches.append(match.group(1).upper())
        if _AUTO_TRAH_MARKER in line:
            switches.append("TRAH")
        if stripped == _CRITERIA_TITLE:
            in_criteria = True

    check_mode, mode_source = _resolve_check_mode(" ".join(echo), mode_label)
    return {
        "solver_seen": solver_seen,
        "diis_rows": tuple(diis_rows),
        "diis_error_at_switch": diis_at_switch,
        "diis_switch_cycle": diis_switch_cycle,
        "diis_resets": resets,
        "converger_switches": tuple(switches),
        "max_abs_energy_step": max(energy_steps) if energy_steps else None,
        "criteria": tuple(criteria),
        "check_mode": check_mode,
        "check_mode_source": mode_source,
    }


def _parse_orbitals(lines: list[str]) -> dict[str, Any]:
    """Take the last ORBITAL ENERGIES table (convention in the module docstring)."""
    last_at = None
    for i, line in enumerate(lines):
        if _ORB_HEADER in line:
            last_at = i
    if last_at is None:
        return {"energies": (), "occupations": ()}
    energies: list[float] = []
    occupations: list[float] = []
    for line in lines[last_at + 1:]:
        if (m := _ORB_ROW_RE.match(line)) is not None:
            occupations.append(float(m.group(2)))
            energies.append(float(m.group(3)))
        elif energies and line.strip():
            # the table has ended (there is no blank line inside it, so any other
            # content means we have left the table body)
            break
    return {"energies": tuple(energies), "occupations": tuple(occupations)}


def _parse_casscf(lines: list[str]) -> dict[str, Any]:
    present = any(_CASSCF_BANNER in line for line in lines)
    if not present:
        return {"present": False}

    converged: bool | None = None
    converged_via: str | None = None
    macro_iterations = 0
    energy: float | None = None
    states: list[dict[str, Any]] = []
    active_occupations: tuple[float, ...] = ()
    current: dict[str, Any] | None = None
    for line in lines:
        if _MACRO_ITER_RE.match(line) is not None:
            macro_iterations += 1
        if (m := _CASSCF_CONVERGED_RE.search(line)) is not None:
            converged = True
            converged_via = m.group(1).lower()
        if (m := _N_OCC_RE.search(line)) is not None:
            values = tuple(
                v for v in (float_or_none(t) for t in m.group(1).split()) if v is not None
            )
            if values:
                active_occupations = values
        if (m := _FINAL_CASSCF_RE.match(line)) is not None:
            energy = float_or_none(m.group(1))
        if (m := _BLOCK_RE.search(line)) is not None:
            current = {"block": int(m.group(1)), "mult": int(m.group(2)), "nroots": int(m.group(3)), "roots": []}
            states.append(current)
        elif (m := _ROOT_RE.match(line)) is not None and current is not None:
            current["roots"].append(
                {
                    "root": int(m.group(1)),
                    "energy": float_or_none(m.group(2)),
                    "de_ev": float_or_none(m.group(3)) if m.group(3) else None,
                    "de_cm1": float_or_none(m.group(4)) if m.group(4) else None,
                }
            )
    return {
        "present": True,
        "converged": converged,
        "converged_via": converged_via,
        "macro_iterations": macro_iterations or None,
        "energy": energy,
        "states": tuple(states),
        "active_occupations": active_occupations,
    }


def _parse_pt2(lines: list[str], kind: str) -> dict[str, Any]:
    """Parse a NEVPT2 / CASPT2 section; both share the block structure and CASPT2
    additionally carries the weight and the denominators."""
    text = "\n".join(lines)
    section = _slice(text, _PT2_SECTION_RE[kind], _PT2_END_RE[kind])
    if not section:
        return {"present": False}

    states: dict[tuple[int, int], dict[str, Any]] = {}
    current: dict[str, Any] | None = None
    denominators: dict[str, float] = {}
    for line in section.splitlines():
        if (m := _PT2_STATE_RE.match(line)) is not None:
            key = (int(m.group(1)), int(m.group(2)))
            current = states.setdefault(
                key,
                {"mult": key[0], "root": key[1], "denominators": {}, "classes": {}},
            )
            continue
        if current is None:
            continue
        if (m := _PT2_DE_RE.search(line)) is not None:
            current["de"] = float_or_none(m.group(1))
        elif (m := _PT2_E0_RE.search(line)) is not None:
            current["e0"] = float_or_none(m.group(1))
        elif (m := _PT2_W0_RE.search(line)) is not None:
            current["weight"] = float_or_none(m.group(1))
        elif (m := _PT2_ETOT_RE.search(line)) is not None:
            current["energy"] = float_or_none(m.group(1))
        elif (m := _PT2_WEIGHT_ITER_RE.match(line)) is not None:
            current.setdefault("weight", float_or_none(m.group(1)))
        elif (m := _PT2_CLASS_RE.search(line)) is not None:
            value = float_or_none(m.group(2))
            if value is not None:
                current["classes"][m.group(1)] = value
        elif kind == "caspt2" and (m := _SMALLEST_DENOM_RE.search(line)) is not None:
            denominators[m.group(1)] = float(m.group(2))
            current["denominators"][m.group(1)] = float(m.group(2))

    out: dict[str, Any] = {"present": True, "states": tuple(states.values())}
    if kind == "caspt2":
        weights = [s["weight"] for s in states.values() if s.get("weight") is not None]
        out["min_reference_weight"] = min(weights) if weights else None
        all_denoms = [v for s in states.values() for v in s["denominators"].values()]
        out["min_denominator"] = min(all_denoms) if all_denoms else None
    return out


def _parse_errors_warnings(lines: list[str]) -> tuple[tuple[str, ...], tuple[str, ...]]:
    errors: list[str] = []
    warnings: list[str] = []
    for line in lines:
        for pattern in _ERROR_MARKERS:
            if pattern.search(line):
                errors.append(line.strip())
                break
        else:
            if _WARNING_LINE_RE.match(line) is not None:
                warnings.append(line.strip())
    return tuple(errors), tuple(warnings)


class OrcaParser:
    """ORCA output parser (program name ``orca``)."""

    program = PROGRAM

    def detect(self, path: str | Path) -> bool:
        """Read the first lines of the file to decide whether it is ORCA output
        (no full parse)."""
        text = read_text(path)
        head = "\n".join(text.splitlines()[:120])
        return bool(_BANNER_RE.search(head)) and "Program Version" in head

    def parse(self, path: str | Path) -> ParseResult:
        text = read_text(path)
        lines = text.splitlines()
        if not _BANNER_RE.search("\n".join(lines[:120])):
            raise ParserError(
                f"{path} is not ORCA output (no ORCA banner found). "
                "Next step: confirm the file comes from ORCA; for another program's "
                "output use the matching parser."
            )
        version = ""
        if (m := _VERSION_RE.search(text)) is not None:
            version = m.group(1)
        errors, warnings = _parse_errors_warnings(lines)
        final_energy = None
        if (m := _FINAL_ENERGY_RE.search(text)) is not None:
            final_energy = float_or_none(m.group(1))
        optimization = _parse_optimization(lines)
        frequencies = _parse_frequencies(lines)
        # does the last frequency block describe the final structure?  With a
        # geometry optimisation, a block printed before the last cycle is not a
        # verification of the optimised geometry (e.g. the initial Hessian of a
        # TS optimisation); without one the question does not arise.
        if frequencies["present"]:
            last_cycle = optimization["last_cycle_line"]
            frequencies["after_geometry"] = (
                last_cycle is None or frequencies["last_block_line"] > last_cycle
            )
        else:
            frequencies["after_geometry"] = None
        sections: dict[str, Any] = {
            "version": version,
            "terminated_normally": TERMINATED_MARK in text,
            "final_energy": final_energy,
            "scf": _parse_scf(lines),
            "orbitals": _parse_orbitals(lines),
            "orbital_composition": _parse_orbital_composition(lines),
            "casscf": _parse_casscf(lines),
            "nevpt2": _parse_pt2(lines, "nevpt2"),
            "caspt2": _parse_pt2(lines, "caspt2"),
            "cc": _parse_cc(lines),
            "frequencies": frequencies,
            "local_spin": _parse_local_spin(lines),
            "optimization": optimization,
            "soc_present": bool(_SOC_MARKERS_RE.search(text)),
            "errors": errors,
            "warnings": warnings,
        }
        return ParseResult(program=PROGRAM, path=str(path), sections=sections)


def facts(result: ParseResult) -> dict[str, Any]:
    """Map an ORCA parse result onto the diagnosis rules' fact fields.

    The field vocabulary is in ``knowledge/rules/README.md`` ("fact field
    vocabulary"); a new field must be registered there first. Missing items are
    left out of the fact table (a condition's missing field means None, as the
    README states).
    """
    s = result.sections
    scf = s.get("scf", {})
    casscf = s.get("casscf", {})
    nevpt2 = s.get("nevpt2", {})
    caspt2 = s.get("caspt2", {})
    frequencies = s.get("frequencies", {})
    optimization = s.get("optimization", {})
    last_freq_block = frequencies["blocks"][-1] if frequencies.get("blocks") else None
    active_occ = casscf.get("active_occupations") or ()
    # NEVPT2 1-hole (V1_i/ITUV) and 1-particle (Vm1_a/TUVA) class contributions:
    # manual §3.16.6 states that their being positive or large relative to IJAB is
    # a sign of (false) intruder states.
    hole_particle = [
        value
        for state in nevpt2.get("states", ())
        for name, value in state.get("classes", {}).items()
        if name in ("V1_i", "Vm1_a")
    ]
    # Composition-derived signals for the "f-block system, but no f character in the
    # active orbitals" check (a d-type solution branch; see the rule file): the element
    # presence comes from the composition table's own atom labels, and the f weight of
    # an orbital is the sum over its f-shell rows. Active orbitals are those with a
    # fractional occupation (the same 0.02-1.98 window the other checks use).
    composition = s.get("orbital_composition", {})
    orbitals = composition.get("orbitals") or ()
    f_block_element_present: bool | None = None
    active_f_weight_max: float | None = None
    if orbitals:
        elements = {shell["element"] for orbital in orbitals for shell in orbital["shells"]}
        f_block_element_present = False
        for element in sorted(elements):
            try:
                f_block_element_present = f_block_element_present or is_f_element(element)
            except ElementError:
                continue  # an unexpected label never turns into a verdict
        f_weights = [
            sum(shell["weight"] for shell in orbital["shells"] if shell["shell"] == "f")
            for orbital in orbitals
            if orbital["occupation"] is not None
            and 0.02 < orbital["occupation"] < 1.98
        ]
        active_f_weight_max = max(f_weights) if f_weights else None

    candidates: dict[str, Any] = {
        "terminated_normally": s.get("terminated_normally"),
        "scf_cycles": scf.get("cycles"),
        "scf_converged": scf.get("converged"),
        "casscf_present": casscf.get("present"),
        "casscf_converged": casscf.get("converged"),
        "casscf_converged_via": casscf.get("converged_via"),
        "casscf_active_occ_min": min(active_occ) if active_occ else None,
        "casscf_active_occ_max": max(active_occ) if active_occ else None,
        "nevpt2_max_hole_particle": max(hole_particle) if hole_particle else None,
        "caspt2_present": caspt2.get("present"),
        "caspt2_min_reference_weight": caspt2.get("min_reference_weight"),
        "caspt2_min_denominator": caspt2.get("min_denominator"),
        "soc_present": s.get("soc_present"),
        "f_block_element_present": f_block_element_present,
        "active_f_weight_max": active_f_weight_max,
        "ts_optimization": optimization.get("ts", False),
        "optimization_converged": optimization.get("converged"),
        "frequency_present": frequencies.get("present"),
        "frequency_after_geometry": frequencies.get("after_geometry"),
        "frequency_imaginary_count": (
            last_freq_block["n_imaginary"] if last_freq_block is not None else None
        ),
        "frequency_min_imaginary": (
            last_freq_block["min_imaginary"] if last_freq_block is not None else None
        ),
    }
    return {key: value for key, value in candidates.items() if value is not None}


register(OrcaParser())
