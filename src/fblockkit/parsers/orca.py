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
from . import ailft as _ailft
from . import rocis_spectra as _rocis_spectra
from . import single_aniso as _single_aniso
from .base import ParserError, float_or_none, read_text, register

PROGRAM = "orca"

# --- identification and global markers --------------------------------------

_BANNER_RE = re.compile(r"\*\s+O\s+R\s+C\s+A\s+\*")
_VERSION_RE = re.compile(r"Program Version\s+(\d+(?:\.\d+)*)")
TERMINATED_MARK = "****ORCA TERMINATED NORMALLY****"

# --- SCF --------------------------------------------------------------------

# The verdict banners: measured forms "*           SCF CONVERGED AFTER  21 CYCLES *"
# and "*       CASSCF NOT CONVERGED AFTER 783 CYCLES *" (the SCF inside CASSCF is
# preceded by an uppercase letter, which is what the lookbehind blocks -- measured
# on generated_yb3_sarc2_trah.out, where the CASSCF macro-iteration verdict must not
# be read as an SCF verdict; it is captured in the casscf section instead).
_SCF_CONVERGED_RE = re.compile(r"(?<![A-Z])SCF CONVERGED AFTER\s+(\d+)\s+CYCLES")
_SCF_NOT_CONVERGED_RE = re.compile(r"(?<![A-Z])SCF NOT CONVERGED AFTER\s+(\d+)\s+CYCLES")
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
# The occupation column may carry a sign: ORCA prints tiny negative natural
# occupations as "-0.0000" (measured, Eu3+ fixture 2026-09-26) -- without the
# optional sign the row fails to match and the table parse stops there.
_ORB_ROW_RE = re.compile(
    r"^\s*(\d+)\s+(-?\d+\.\d{4})\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$"
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
# The module-level failure banner (measured on generated_yb3_sarc2_trah.out:
# "*       CASSCF NOT CONVERGED AFTER 783 CYCLES       *") -- the flip side of
# the SCF lookbehind above: this verdict belongs to the casscf section.
_CASSCF_NOT_CONVERGED_RE = re.compile(r"CASSCF NOT CONVERGED AFTER\s+(\d+)\s+CYCLES")
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
    weights: list[dict[tuple[int, str, str, int], float]] = []

    def _flush() -> None:
        for entry, shell_map in zip(entries, weights):
            entry["shells"] = tuple(
                {
                    "atom": atom,
                    "element": element,
                    "shell": shell,
                    "shell_index": shell_index,
                    "weight": weight,
                }
                for (atom, element, shell, shell_index), weight in sorted(
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
            label = match.group(3)
            shell = label[0]  # s / p / d / f / g (first char of the AO label)
            # the AO label also carries the shell counter ("1s", "2p", "1dz2"): it
            # indexes the shells of that (element, angular momentum) channel in the
            # order the basis prints them, which is what the basis reader below and
            # the diffuse-orbital check need
            # this table labels AOs by letter only ("s", "pz", "dz2"), so the shell
            # counter is unknown here (None); the orca_2json exports carry the
            # numbered labels ("1s", "2p") and their readers give real indices
            shell_index = None
            values = _floats(match.group(4))
            for shell_map, value in zip(weights, values):
                key = (atom, element, shell, shell_index)
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

# --- printed basis set (BASIS SET IN INPUT FORMAT) ---------------------------
# Written only when the run asks for it (!PrintBasis): the block lists, per element,
# the shells in ORCA's own input format -- one header "<L> <nprim>" per contracted
# shell and one "index exponent coefficient(s)" row per primitive.  The order inside
# one angular momentum runs tight to diffuse (measured on the N2 def2-SVP fixture),
# but the diffuse-orbital check takes the minimum exponent rather than relying on it.

_BASIS_MARKER = "BASIS SET IN INPUT FORMAT"
_BASIS_SHELL_RE = re.compile(r"^\s*([SPDFGHI])\s+(\d+)\s*$")
_BASIS_PRIMITIVE_RE = re.compile(r"^\s*\d+\s+([-+\d.eEdD]+(?:\s+[-+\d.eEdD]+)+)\s*$")


def _parse_basis(lines: list[str]) -> dict[str, Any]:
    """Per-element shells with their primitive exponents (requires !PrintBasis)."""
    start = None
    for index, line in enumerate(lines):
        if _BASIS_MARKER in line:
            start = index + 1
            break
    if start is None:
        return {"present": False, "elements": (), "ecp_blocks": 0}
    elements: list[dict[str, Any]] = []
    element: str | None = None
    shells: list[dict[str, Any]] = []
    shell: dict[str, Any] | None = None
    remaining = 0
    ecp_blocks = 0
    in_ecp = False

    def _close_shell() -> None:
        nonlocal shell
        if shell is not None:
            shells.append(shell)
            shell = None

    def _close_element() -> None:
        nonlocal element
        _close_shell()
        if element is not None:
            elements.append({"element": element, "shells": tuple(shells)})
            element = None
        shells.clear()

    for line in lines[start:]:
        stripped = line.strip()
        if in_ecp:
            if "end;" in stripped or stripped == "end":
                in_ecp = False
            continue
        if stripped.startswith("NewGTO"):
            _close_element()
            tokens = stripped.split()
            if len(tokens) < 2:
                raise ParserError(
                    "a NewGTO line in the printed basis carries no element symbol. "
                    "Next step: report the output; the basis reader cannot name the "
                    "element."
                )
            element = tokens[1].capitalize()
            shell = None
            remaining = 0
            continue
        if stripped.startswith("NewECP"):
            # an ECP block sits next to the basis for the elements that have one; it
            # carries no Gaussian exponents, so it is counted and skipped
            _close_element()
            ecp_blocks += 1
            in_ecp = True
            continue
        if stripped.startswith("#"):
            continue
        if stripped and set(stripped) <= {"-", " "}:
            continue  # the divider under the section marker
        if stripped in ("end;", "end"):
            _close_element()
            continue
        if not stripped:
            continue
        match = _BASIS_SHELL_RE.match(line)
        if match is not None:
            _close_shell()
            shell = {"angular": match.group(1).lower(), "exponents": [], "coefficients": []}
            remaining = int(match.group(2))
            continue
        match = _BASIS_PRIMITIVE_RE.match(line)
        if match is not None:
            if shell is None or remaining <= 0:
                raise ParserError(
                    f"a primitive row appears outside any shell in the printed basis "
                    f"({stripped!r}). Next step: report the output; the basis block does "
                    "not follow the expected shape."
                )
            numbers = [
                float(token.replace("D", "E").replace("d", "e"))
                for token in match.group(1).split()
            ]
            shell["exponents"].append(numbers[0])
            shell["coefficients"].append(tuple(numbers[1:]))
            remaining -= 1
            continue
        # anything else ends the block (the startup banner follows it)
        break
    _close_element()
    if not elements:
        return {"present": False, "elements": (), "ecp_blocks": ecp_blocks}
    for entry in elements:
        for item in entry["shells"]:
            item["exponents"] = tuple(item["exponents"])
            item["coefficients"] = tuple(item["coefficients"])
    return {"present": True, "elements": tuple(elements), "ecp_blocks": ecp_blocks}


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


# --- final geometry and normal modes -----------------------------------------
# Measured format (fixture: fhh_optts_freq and the whole set): the last
# "CARTESIAN COORDINATES (ANGSTROEM)" block carries the final structure
# ("<symbol> <x> <y> <z>"); the "(A.U.)" block repeats it in bohr and prints,
# per row, "<no> <symbol> <ZA> <frag> <mass> <x> <y> <z>" -- the atomic masses
# are read from there so any mass de-weighting uses the run's own values.
# The "NORMAL MODES" block prints the mass-weighted Hessian's eigenvectors in
# groups of up to six: a header line of mode indices, then one row per
# Cartesian coordinate ("<i> <v0> <v1> ...").  The block's own header says the
# vectors are "Cartesian displacements weighted by ... 1/sqrt(m[i])"; their
# Euclidean norm is 1, so the Cartesian displacement pattern is v_i /
# sqrt(m_i) (the reading the recipe layer uses).  The last block describes the
# final structure (the frequencies' convention).
_FINAL_GEO_HEADER = "CARTESIAN COORDINATES (ANGSTROEM)"
_FINAL_MASS_HEADER = "CARTESIAN COORDINATES (A.U.)"
_NM_HEADER = "NORMAL MODES"
_NM_END = "IR SPECTRUM"
_GEO_ROW_RE = re.compile(r"^\s*([A-Za-z]{1,2})\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$")
_MASS_ROW_RE = re.compile(
    r"^\s*\d+\s+([A-Za-z]{1,2})\s+-?\d+\.\d+\s+\d+\s+(-?\d+\.\d+)\s+"
)
_NM_LABEL_ROW_RE = re.compile(r"^\s*\d+(\s+-?\d+\.\d+|\s+\d+\.\d+e[+-]\d+)*\s*$")


def _parse_final_geometry(lines: list[str]) -> dict[str, Any]:
    starts = [i for i, line in enumerate(lines) if line.strip() == _FINAL_GEO_HEADER]
    if not starts:
        return {"present": False, "atoms": (), "masses": ()}
    start = starts[-1]
    atoms: list[tuple[str, float, float, float]] = []
    for line in lines[start + 2:]:
        if (m := _GEO_ROW_RE.match(line)) is not None:
            atoms.append((m.group(1), float(m.group(2)), float(m.group(3)), float(m.group(4))))
            continue
        if atoms:
            break  # the block body has ended
    masses: list[float] = []
    mass_starts = [i for i, line in enumerate(lines) if line.strip() == _FINAL_MASS_HEADER]
    if mass_starts:
        for line in lines[mass_starts[-1] + 2:]:
            if (m := _MASS_ROW_RE.match(line)) is not None:
                masses.append(float(m.group(2)))
                continue
            if masses:
                break
    if len(masses) != len(atoms):
        masses = []  # the (A.U.) block does not mirror the (ANGSTROEM) one
    return {"present": True, "atoms": tuple(atoms), "masses": tuple(masses)}


def _parse_normal_modes(lines: list[str]) -> dict[str, Any]:
    starts = [i for i, line in enumerate(lines) if line.strip() == _NM_HEADER]
    if not starts:
        return {"present": False, "n_coord": 0, "modes": ()}
    start = starts[-1]
    columns: dict[int, dict[int, float]] = {}
    pending: list[int] = []
    for line in lines[start + 1:]:
        stripped = line.strip()
        if stripped == _NM_END:
            break
        tokens = stripped.split()
        if not tokens:
            continue
        if stripped.startswith("-----") and any(columns.values()):
            break  # the block's trailing rule (the title rule comes before data)
        if all(token.isdigit() for token in tokens):
            # a header line: the mode indices of the next column group
            pending = [int(token) for token in tokens]
            for index in pending:
                columns.setdefault(index, {})
            continue
        if not tokens[0].isdigit() or not pending:
            continue
        coordinate = int(tokens[0])
        values = [float(token) for token in tokens[1:]]
        if len(values) != len(pending):
            continue
        for index, value in zip(pending, values):
            columns[index][coordinate] = value
    modes = []
    for index in sorted(columns):
        column = columns[index]
        if not column:
            continue
        vector = tuple(column.get(i, 0.0) for i in range(max(column) + 1))
        modes.append({"index": index, "vector": vector})
    n_coord = max((len(mode["vector"]) for mode in modes), default=0)
    return {"present": True, "n_coord": n_coord, "modes": tuple(modes)}


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

# --- dipole-moment blocks (SCF / CASSCF / CASCI, one format) ------------------
# Measured on ORCA 6.1.1 (2026-09-27): the SCF block has no State line, the
# CASSCF/CASCI block prints "State: <n>" (single root: 0, relaxed density; a
# state-averaged run prints one block with "Method: CASSCF/ALL STATES AVERAGE",
# "State: -1" and the unrelaxed density -- per-state dipoles are not printed).
_DIPOLE_TITLE = "DIPOLE MOMENT"
_DIPOLE_INT_RE = {
    "state": re.compile(r"^State\s*:\s*(-?\d+)\s*$"),
    "multiplicity": re.compile(r"^Multiplicity\s*:\s*(-?\d+)\s*$"),
    "irrep": re.compile(r"^Irrep\s*:\s*(-?\d+)\s*$"),
}
_DIPOLE_TOTAL_RE = re.compile(
    r"^Total Dipole Moment\s*:\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*$"
)
_DIPOLE_AU_RE = re.compile(r"^Magnitude \(a\.u\.\)\s*:\s*([-\d.eE+]+)\s*$")
_DIPOLE_DEBYE_RE = re.compile(r"^Magnitude \(Debye\)\s*:\s*([-\d.eE+]+)\s*$")


def _parse_dipole(lines: list[str]) -> dict[str, Any]:
    """Every DIPOLE MOMENT block of the output, in file order.

    One format serves SCF, CASSCF and CASCI runs (measured); a block is kept
    only when it carries the total vector and both magnitudes, so the property
    summary lines (``Dipole moment ... YES``) do not produce phantom blocks.
    """
    blocks = []
    for index, line in enumerate(lines):
        if line.strip() != _DIPOLE_TITLE:
            continue
        block = _parse_dipole_block(lines[index : index + 30])
        if block is not None:
            blocks.append(block)
    return {"present": bool(blocks), "blocks": tuple(blocks)}


def _parse_dipole_block(chunk: list[str]) -> dict[str, Any] | None:
    block: dict[str, Any] = {}
    for line in chunk:
        stripped = line.strip()
        if stripped.startswith("Method") and ":" in stripped:
            block["method"] = stripped.split(":", 1)[1].strip()
        elif stripped.startswith("Level") and ":" in stripped:
            block["level"] = stripped.split(":", 1)[1].strip()
        elif stripped.startswith("Energy") and ":" in stripped:
            block["energy"] = float(stripped.split(":", 1)[1].strip().split()[0])
        else:
            for key, pattern in _DIPOLE_INT_RE.items():
                if (match := pattern.match(stripped)) is not None:
                    block[key] = int(match.group(1))
                    break
            else:
                if (match := _DIPOLE_TOTAL_RE.match(stripped)) is not None:
                    block["total"] = tuple(float(match.group(i)) for i in (1, 2, 3))
                elif (match := _DIPOLE_AU_RE.match(stripped)) is not None:
                    block["magnitude_au"] = float(match.group(1))
                elif (match := _DIPOLE_DEBYE_RE.match(stripped)) is not None:
                    block["magnitude_debye"] = float(match.group(1))
    if "total" not in block or "magnitude_debye" not in block:
        return None
    return block
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
    diis_blocks: list[tuple[tuple[int, float], ...]] = []
    diis_block_switches: list[tuple[int | None, float | None]] = []
    block: list[tuple[int, float]] = []
    block_switch: tuple[int, float] | None = None
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
            # a new table ends the previous block; measured: a geometry
            # optimisation prints one table per SCF cycle and the cycle count
            # restarts at 1, so the blocks are per-geometry-step
            if block:
                diis_blocks.append(tuple(block))
                diis_block_switches.append(block_switch or (None, None))
                block = []
                block_switch = None
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
                    row = (index, values[columns.index("DIISErr") - 1])
                    diis_rows.append(row)
                    block.append(row)
                    if "Delta-E" in columns:
                        energy_steps.append(abs(values[columns.index("Delta-E") - 1]))
                continue

        if _TURN_ON_DIIS in line and diis_at_switch is None and diis_rows:
            diis_at_switch = diis_rows[-1][1]
            diis_switch_cycle = diis_rows[-1][0]
        if _TURN_ON_DIIS in line and block_switch is None and block:
            block_switch = block[-1]
        if _DIIS_RESET in line:
            resets += 1
        if (match := _CONVERGER_SWITCH_RE.search(line)) is not None:
            switches.append(match.group(1).upper())
        if _AUTO_TRAH_MARKER in line:
            switches.append("TRAH")
        if stripped == _CRITERIA_TITLE:
            in_criteria = True

    if block:
        diis_blocks.append(tuple(block))
        diis_block_switches.append(block_switch or (None, None))
    check_mode, mode_source = _resolve_check_mode(" ".join(echo), mode_label)
    return {
        "solver_seen": solver_seen,
        "diis_rows": tuple(diis_rows),
        "diis_blocks": tuple(diis_blocks),
        "diis_block_switches": tuple(diis_block_switches),
        "diis_error_at_switch": diis_at_switch,
        "diis_switch_cycle": diis_switch_cycle,
        "diis_resets": resets,
        "converger_switches": tuple(switches),
        "max_abs_energy_step": max(energy_steps) if energy_steps else None,
        "criteria": tuple(criteria),
        "check_mode": check_mode,
        "check_mode_source": mode_source,
        "convergence_block": _raw_convergence_block(lines),
    }


def _raw_convergence_block(lines: list[str]) -> tuple[str, ...]:
    """The last SCF CONVERGENCE block, verbatim (title plus the criterion rows).

    The check mode decides which of these rows are actually enforced; the verbatim
    view exists so the informational rows (density, DIIS error, ...) can be read
    as printed, without any interpretation imposed on them.
    """
    start = None
    for i, line in enumerate(lines):
        if line.strip() == _CRITERIA_TITLE:
            start = i
    if start is None:
        return ()
    block: list[str] = []
    row_seen = False
    for line in lines[start:]:
        if line.strip() == _CRITERIA_TITLE and not block:
            block.append(line)
            continue
        stripped = line.strip()
        if _CRITERION_RE.match(line) is not None:
            block.append(line)
            row_seen = True
            continue
        if not row_seen and (not stripped or set(stripped) <= _TABLE_DIVIDER_CHARS):
            block.append(line)  # divider / blank line between the title and the table
            continue
        break
    return tuple(block)


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


# --- 5.2: the QDPT magnetic-property blocks (EPR g/D tensors and chi) ---------
# Measured on ORCA 6.1.1 (fixtures co_plus_qdpt.*, o2_qdpt2.*; see the fixtures
# README for the run conditions):
#
# - the effective-Hamiltonian g-matrix prints as
#   "ELECTRONIC G-MATRIX FROM EFFECTIVE HAMILTONIAN" -> "Spin multiplicity = n"
#   -> "total g-matrix:" + a 3x3 table, then "g-factors: g1 g2 g3 iso = ...";
#   other g-like tables exist in the same output (the Kramers-pair Zeeman
#   matrices, "ELECTRONIC G-MATRIX: S contribution") and are NOT this tensor;
# - the zero-field splitting prints as several "ZERO-FIELD SPLITTING" blocks,
#   each with a variant line ("2ND ORDER SOC CONTRIBUTION", "EFFECTIVE
#   HAMILTONIAN SOC CONTRIBUTION", "... SOC and SSC CONTRIBUTION") and a set
#   of "Raw matrix (cm-1)" + "Eigenvalues (traceless)" + "Eigenvectors" +
#   Euler angles + "D = ..." + "E/D = ..." tables;
# - with DoSusceptibility, "SOC CORRECTED MAGNETIZATION AND/OR SUSCEPTIBILITY"
#   carries, per temperature, "TEMPERATURE/K: x" + "Tensor in molecular frame
#   (cm3*K/mol)" + a 3x3 table.
_FLOAT3_RE = re.compile(r"^\s*([-+0-9.eEdD]+)\s+([-+0-9.eEdD]+)\s+([-+0-9.eEdD]+)\s*$")
_G_ROW_RE = re.compile(
    r"^\s*([0-2])\s+([-+0-9.eEdD]+)\s+([-+0-9.eEdD]+)\s+([-+0-9.eEdD]+)\s*$"
)
_G_FACTORS_RE = re.compile(
    r"g-factors:\s*\n\s*([-+0-9.eEdD]+)\s+([-+0-9.eEdD]+)\s+([-+0-9.eEdD]+)"
    r"\s+iso\s*=\s*([-+0-9.eEdD]+)"
)
_SPIN_MULT_RE = re.compile(r"Spin multiplicity\s*=\s*(\d+)")
_D_LINE_RE = re.compile(r"^D\s*=\s*([-+0-9.eEdD]+)\s+cm-1")
_ED_LINE_RE = re.compile(r"^E/D\s*=\s*([-+0-9.eEdD]+)")
_TEMPERATURE_RE = re.compile(r"TEMPERATURE/K:\s*([-+0-9.eEdD]+)")
_COLUMN_HEADER_RE = re.compile(r"^\s*[0-2]\s+[0-2]\s+[0-2]\s*$")


def _read_matrix(lines: list[str], start: int, *, indexed: bool = False):
    """Three consecutive 3-number lines (optionally with a leading row index)."""
    rows: list[tuple[float, ...]] = []
    index = start
    while index < len(lines) and len(rows) < 3:
        line = lines[index]
        if not line.strip():
            index += 1
            continue
        if not rows and _COLUMN_HEADER_RE.match(line) is not None:
            index += 1  # the "0 1 2" column-id line above an indexed table
            continue
        match = _G_ROW_RE.match(line) if indexed else _FLOAT3_RE.match(line)
        if match is None:
            return None
        values = tuple(float_or_none(value) for value in match.groups()[-3:])
        if any(value is None for value in values):
            return None
        rows.append(values)  # type: ignore[arg-type]
        index += 1
    return tuple(rows) if len(rows) == 3 else None


def _parse_epr(lines: list[str]) -> dict[str, Any]:
    """The effective-Hamiltonian g-matrix (the QDPT EPR block)."""
    result: dict[str, Any] = {
        "present": False,
        "g_matrix": None,
        "g_factors": None,
        "iso": None,
        "multiplicity": None,
        "n_blocks": 0,
    }
    text = "\n".join(lines)
    if "ELECTRONIC G-MATRIX FROM EFFECTIVE HAMILTONIAN" not in text:
        return result
    result["present"] = True
    for index, line in enumerate(lines):
        if "ELECTRONIC G-MATRIX FROM EFFECTIVE HAMILTONIAN" in line:
            result["n_blocks"] += 1
            if result["g_matrix"] is not None:
                continue  # keep the first block (the ground multiplet), count the rest
            for probe in range(index + 1, min(index + 20, len(lines))):
                if (match := _SPIN_MULT_RE.search(lines[probe])) is not None:
                    result["multiplicity"] = int(match.group(1))
                if "total g-matrix:" in lines[probe]:
                    matrix = _read_matrix(lines, probe + 1, indexed=True)
                    if matrix is not None:
                        result["g_matrix"] = matrix
                    break
    if (match := _G_FACTORS_RE.search(text)) is not None:
        result["g_factors"] = tuple(float(match.group(i)) for i in (1, 2, 3))
        result["iso"] = float(match.group(4))
    return result


#: preference order when several ZFS variants are present (the measured variants)
ZFS_PREFERENCE = (
    "effective hamiltonian soc and ssc",
    "2nd order soc and ssc",
    "effective hamiltonian soc",
    "2nd order soc",
)


def _parse_zfs(lines: list[str]) -> dict[str, Any]:
    """All zero-field-splitting blocks (variant, matrices, D and E/D)."""
    blocks: list[dict[str, Any]] = []
    header_indices = [
        index for index, line in enumerate(lines) if line.strip() == "ZERO-FIELD SPLITTING"
    ]
    for header in header_indices:
        variant = None
        for probe in range(header + 1, min(header + 4, len(lines))):
            if lines[probe].strip():
                variant = lines[probe].strip()
                break
        d_line = d_value = ed_value = None
        matrix = eigenvalues = eigenvectors = euler = None
        for probe in range(header + 1, min(header + 120, len(lines))):
            line = lines[probe]
            if line.strip() == "ZERO-FIELD SPLITTING" and probe != header:
                break
            if (match := _D_LINE_RE.match(line.strip())) is not None:
                d_line = probe
                d_value = float_or_none(match.group(1))
                break
        if d_line is None:
            continue  # a header without its D line (an aborted block): skip
        for probe in range(header + 1, min(d_line + 3, len(lines))):
            stripped = lines[probe].strip()
            if stripped.startswith("Raw matrix (cm-1)"):
                matrix = _read_matrix(lines, probe + 1)
            elif stripped.startswith("Eigenvalues (traceless)"):
                if (match := _FLOAT3_RE.match(lines[probe + 1])) is not None:
                    eigenvalues = tuple(float_or_none(v) for v in match.groups())
            elif stripped.startswith("Eigenvectors:"):
                eigenvectors = _read_matrix(lines, probe + 1)
            elif stripped.startswith("Euler angles"):
                if (match := _FLOAT3_RE.match(lines[probe + 1])) is not None:
                    euler = tuple(float_or_none(v) for v in match.groups())
            elif (match := _ED_LINE_RE.match(stripped)) is not None:
                ed_value = float_or_none(match.group(1))
        blocks.append(
            {
                "variant": variant,
                "raw_matrix": matrix,
                "eigenvalues_traceless": eigenvalues,
                "eigenvectors": eigenvectors,
                "euler_deg": euler,
                "D_cm1": d_value,
                "E_over_D": ed_value,
            }
        )
    preferred = None
    for wanted in ZFS_PREFERENCE:
        for block in blocks:
            if wanted in (block["variant"] or "").strip().lower():
                preferred = block
                break
        if preferred is not None:
            break
    return {"present": bool(blocks), "blocks": tuple(blocks), "preferred": preferred}


def _parse_susceptibility(lines: list[str]) -> dict[str, Any]:
    """The temperature-dependent molar susceptibility tensors (cm3*K/mol)."""
    temperatures: list[float] = []
    tensors: list[tuple] = []
    for index, line in enumerate(lines):
        if (match := _TEMPERATURE_RE.search(line)) is None:
            continue
        for probe in range(index + 1, min(index + 6, len(lines))):
            if "Tensor in molecular frame" in lines[probe]:
                matrix = _read_matrix(lines, probe + 1)
                if matrix is not None:
                    temperatures.append(float(match.group(1)))
                    tensors.append(matrix)
                break
    return {
        "present": bool(tensors),
        "temperatures_k": tuple(temperatures),
        "tensors_cm3k_mol": tuple(tensors),
    }


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
        if (m := _CASSCF_NOT_CONVERGED_RE.search(line)) is not None:
            converged = False
            converged_via = f"not converged after {int(m.group(1))} cycles"

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
            "basis": _parse_basis(lines),
            "casscf": _parse_casscf(lines),
            "nevpt2": _parse_pt2(lines, "nevpt2"),
            "caspt2": _parse_pt2(lines, "caspt2"),
            "cc": _parse_cc(lines),
            "frequencies": frequencies,
            "final_geometry": _parse_final_geometry(lines),
            "normal_modes": _parse_normal_modes(lines),
            "local_spin": _parse_local_spin(lines),
            "epr": _parse_epr(lines),
            "zfs": _parse_zfs(lines),
            "susceptibility": _parse_susceptibility(lines),
            "single_aniso": _single_aniso.parse_segments(lines),
            "rocis": _rocis_spectra.parse_rocis(lines),
            "ailft": _ailft.parse_ailft(lines),
            "optimization": optimization,
            "dipole": _parse_dipole(lines),
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
