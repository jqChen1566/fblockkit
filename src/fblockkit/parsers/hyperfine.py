"""Hyperfine / EFG section of an ORCA EPRNMR run (the Mossbauer-facing block).

``%eprnmr`` with a ``Nuclei`` list prints per nucleus the electric and
magnetic hyperfine structure (manual section 7.51.3; the relevant flags are
``aiso``/``adip``/``aorb`` for the A tensor and ``fgrad``/``rho`` for the
electric field gradient and the density at the nucleus).  Two measured
interface facts:

- when ``Nuclei`` is given, **the coordinates block must precede %eprnmr**
  (otherwise the engine stops with "[EPRNMR] block: nuclear properties are
  requested but no coordinates have been read!");
- the A tensor in MHz needs the per-nucleus nuclear parameters (I, P: the
  ``PPP``/``III`` keywords); without them the engine prints
  ``Isotope=0 I=0.0 P=0.0000`` and an all-zero A matrix (measured).

Measured section shape (ORCA 6.1.1, probes on CeF3, 2026-09-29): a header
with the per-property nucleus counts, then per nucleus a block with

- ``Nucleus <i><El>: A : Isotope=.. I=.. P=.. MHz/au**3`` and
  ``Q : Isotope=.. I=.. Q=.. barn``;
- ``HFC: iso=.. dip=.. orb=.. gauge=..`` and ``EFG: fgrad=.. rho=..``;
- ``Total HFC matrix (all values in MHz)``: a 3x3 matrix printed as three
  bare rows, then the A(FC)/A(SD) rows, a dashed rule, the A(Tot) row with
  a trailing ``A(iso)=`` element, and an ``Orientation:`` matrix;
- ``Raw EFG matrix (all values in a.u.**-3)``: the 3x3 matrix, the
  ``V(Tot)`` principal-value row, the ``Orientation:`` matrix, and the
  ``Contributions:`` rows ``V(El )`` / ``V(Nuc)`` followed by
  ``RHO(0)= .. a.u.**-3``.

The value of RHO(0) depends on the basis-set domain: an all-electron basis
carries the core density (order 10^5 a.u.**-3 for Ce), a small-core ECP
basis reports the valence density only (order 10^-1) -- same label, do not
compare across domains (measured: 727216.9 vs 0.075618828 on the same
geometry).
"""

from __future__ import annotations

import re
from typing import Any

__all__ = ["parse_hyperfine"]

_SECTION = "ELECTRIC AND MAGNETIC HYPERFINE STRUCTURE"
_FLOAT = r"[-+]?\d+\.\d+(?:[eE][-+]?\d+)?"
_COUNT_RE = re.compile(r"Number of nuclei to compute (\w+\([^)]*\)|\w+)\s*:\s*(\d+)")
_METHOD_RE = re.compile(r"^\s*Method\s*:\s*(.+?)\s*$")
_MULT_RE = re.compile(r"^\s*Multiplicity\s*:\s*(\d+)\s*$")
_NUCLEUS_RE = re.compile(r"^\s*Nucleus\s+(\d+)([A-Za-z]{1,2}):\s*A\s*:\s*(.+)$")
_Q_RE = re.compile(r"^\s*Q\s*:\s*(.+)$")
_FLAGS_RE = re.compile(r"^\s*(HFC|EFG):\s*(.+)$")
_ISOTOPE_RE = re.compile(
    r"Isotope=\s*(-?\d+)\s+I=\s*([\d.]+)\s+P=\s*(-?[\d.]+)\s*MHz/au\*\*3"
)
_QUAD_RE = re.compile(r"Isotope=\s*(-?\d+)\s+I=\s*([\d.]+)\s+Q=\s*(-?[\d.]+)\s*barn")
_ROW3_RE = re.compile(rf"^\s*({_FLOAT})\s+({_FLOAT})\s+({_FLOAT})\s*$")
_ABC_RE = re.compile(rf"^\s*A\((FC|SD)\)\s+({_FLOAT})\s+({_FLOAT})\s+({_FLOAT})\s*$")
_ATOT_RE = re.compile(
    rf"^\s*A\(Tot\)\s+({_FLOAT})\s+({_FLOAT})\s+({_FLOAT})\s+A\(iso\)=\s*({_FLOAT})\s*$"
)
_VTOT_RE = re.compile(
    rf"^\s*V\(Tot\)\s+({_FLOAT})\s+({_FLOAT})\s+({_FLOAT})\s*$"
)
_VSRC_RE = re.compile(
    rf"^\s*V\((El\s*|Nuc)\)\s+({_FLOAT})\s+({_FLOAT})\s+({_FLOAT})\s*$"
)
_RHO_RE = re.compile(rf"RHO\(0\)=\s*({_FLOAT})\s*a\.u\.\*\*-3")
# the engine's own quadrupole-tensor block (appears when Q/I parameters are
# given; measured 2026-09-29): "e**2qQ = .. MHz", "e**2qQ/(4I*(2I-1))= .. MHz"
# and an "eta = .." line, with a NOTE giving the diagonal representation.
_QTENSOR_RE = re.compile(
    rf"Quadrupole tensor eigenvalues \(in MHz;Q=\s*({_FLOAT})\s*I=\s*({_FLOAT})\)"
)
_EQ2Q_DIV_RE = re.compile(rf"e\*\*2qQ/\(4I\*\(2I-1\)\)=\s*({_FLOAT})\s*MHz")
_EQ2Q_RE = re.compile(rf"e\*\*2qQ\s*=\s*({_FLOAT})\s*MHz")
_ETA_RE = re.compile(rf"^\s*eta\s*=\s*({_FLOAT})\s*$")
_ORIENT_RE = re.compile(rf"^\s*([XYZ])\s+({_FLOAT})\s+({_FLOAT})\s+({_FLOAT})\s*$")
_HFC_TITLE = "Total HFC matrix"
_EFG_TITLE = "Raw EFG matrix"


def parse_hyperfine(lines: list[str]) -> dict[str, Any]:
    """The per-nucleus hyperfine and EFG blocks; values verbatim."""
    text = "\n".join(lines)
    if _SECTION not in text:
        return {"present": False, "header": {}, "nuclei": []}
    out: dict[str, Any] = {"present": True, "header": {}, "nuclei": []}
    nucleus: dict[str, Any] | None = None
    table: str | None = None  # "hfc" or "efg"
    row_buffer: list[list[float]] = []
    orient_buffer: dict[str, list[float]] = {}

    def _close_nucleus() -> None:
        nonlocal row_buffer, orient_buffer
        row_buffer = []
        orient_buffer = {}

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if (match := _COUNT_RE.search(line)) is not None:
            out["header"][match.group(1)] = int(match.group(2))
            continue
        if (match := _METHOD_RE.match(line)) is not None and "method" not in out["header"]:
            out["header"]["method"] = match.group(1)
            continue
        if (match := _MULT_RE.match(line)) is not None and "multiplicity" not in out["header"]:
            out["header"]["multiplicity"] = int(match.group(1))
            continue
        if (match := _NUCLEUS_RE.match(line)) is not None:
            _close_nucleus()
            nucleus = {
                "index": int(match.group(1)),
                "element": match.group(2),
                "a": {},
                "q": {},
                "hfc_flags": {},
                "efg_flags": {},
                "hfc_matrix": [],
                "A_FC": None,
                "A_SD": None,
                "A_Tot": None,
                "A_iso": None,
                "hfc_orientation": {},
                "efg_raw": [],
                "V_Tot": None,
                "efg_orientation": {},
                "V_el": None,
                "V_nuc": None,
                "rho0": None,
            }
            if (m2 := _ISOTOPE_RE.search(match.group(3))) is not None:
                nucleus["a"] = {
                    "isotope": int(m2.group(1)),
                    "I": float(m2.group(2)),
                    "P_MHz_au3": float(m2.group(3)),
                }
            out["nuclei"].append(nucleus)
            table = None
            row_buffer = []
            orient_buffer = {}
            continue
        if nucleus is None:
            continue
        if (match := _Q_RE.match(line)) is not None:
            if (m2 := _QUAD_RE.search(match.group(1))) is not None:
                nucleus["q"] = {
                    "isotope": int(m2.group(1)),
                    "I": float(m2.group(2)),
                    "Q_barn": float(m2.group(3)),
                }
            continue
        if (match := _FLAGS_RE.match(line)) is not None:
            flags = {}
            for key, value in re.findall(r"(\w+)\s*=\s*(YES|NO)", match.group(2)):
                flags[key.lower()] = value == "YES"
            if match.group(1) == "HFC":
                nucleus["hfc_flags"] = flags
            else:
                nucleus["efg_flags"] = flags
            continue
        if _HFC_TITLE in line:
            table = "hfc"
            row_buffer = []
            orient_buffer = {}
            continue
        if _EFG_TITLE in line:
            table = "efg"
            row_buffer = []
            orient_buffer = {}
            continue
        if table == "hfc":
            if (match := _ABC_RE.match(line)) is not None:
                key = "A_" + match.group(1)
                nucleus[key] = [float(match.group(i)) for i in (2, 3, 4)]
                continue
            if (match := _ATOT_RE.match(line)) is not None:
                nucleus["A_Tot"] = [float(match.group(i)) for i in (1, 2, 3)]
                nucleus["A_iso"] = float(match.group(4))
                continue
            if (match := _ORIENT_RE.match(line)) is not None:
                orient_buffer[match.group(1)] = [float(match.group(i)) for i in (2, 3, 4)]
                if len(orient_buffer) == 3:
                    nucleus["hfc_orientation"] = orient_buffer
                continue
            if (match := _ROW3_RE.match(line)) is not None:
                row_buffer.append([float(match.group(i)) for i in (1, 2, 3)])
                if len(row_buffer) == 3:
                    nucleus["hfc_matrix"] = row_buffer
                continue
            continue
        if table == "efg":
            if (match := _VTOT_RE.match(line)) is not None:
                nucleus["V_Tot"] = [float(match.group(i)) for i in (1, 2, 3)]
                continue
            if (match := _VSRC_RE.match(line)) is not None:
                key = "V_el" if match.group(1).strip() == "El" else "V_nuc"
                nucleus[key] = [float(match.group(i)) for i in (2, 3, 4)]
                continue
            if (match := _RHO_RE.search(line)) is not None:
                nucleus["rho0"] = float(match.group(1))
                continue
            if (match := _QTENSOR_RE.search(line)) is not None:
                nucleus["quad_tensor"] = {
                    "Q_barn": float(match.group(1)),
                    "I": float(match.group(2)),
                }
                continue
            if (match := _EQ2Q_DIV_RE.search(line)) is not None:
                nucleus["e2qQ_div4I"] = float(match.group(1))
                continue
            if (match := _EQ2Q_RE.search(line)) is not None:
                nucleus["e2qQ_MHz"] = float(match.group(1))
                continue
            if (match := _ETA_RE.match(line)) is not None:
                nucleus["engine_eta"] = float(match.group(1))
                continue
            if (match := _ORIENT_RE.match(line)) is not None:
                orient_buffer[match.group(1)] = [float(match.group(i)) for i in (2, 3, 4)]
                if len(orient_buffer) == 3:
                    nucleus["efg_orientation"] = orient_buffer
                continue
            if (match := _ROW3_RE.match(line)) is not None:
                row_buffer.append([float(match.group(i)) for i in (1, 2, 3)])
                if len(row_buffer) == 3:
                    nucleus["efg_raw"] = row_buffer
                continue
            continue
    _close_nucleus()
    for item in out["nuclei"]:
        item.pop("_pending", None)
    return out
