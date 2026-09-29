"""ROCIS embedded spectra (the ROCIS module's absorption-spectrum blocks).

The ORCA ROCIS module (manual section 5.7; ``%rocis`` input block) prints, on
top of its excitation-energy table, one absorption-spectrum block per
requested oscillator-strength model.  Measured on ORCA 6.1.1 (fixtures
``rocis/``; [FeCl4]2- L-edge probe, 2026-09-29):

- ``ROCIS-EXCITATION SPECTRA`` section header, then ``ABSORPTION SPECTRUM VIA
  TRANSITION ELECTRIC DIPOLE MOMENTS`` (fosc(D2), D2, DX/DY/DZ) and ``... VIA
  TRANSITION VELOCITY DIPOLE MOMENTS`` (fosc(P2), P2, PX/PY/PZ);
- with ``DecomposeFosc`` up to five ``COMBINED ELECTRIC DIPOLE + MAGNETIC
  DIPOLE + ELECTRIC QUADRUPOLE`` variants (length / velocity / origin
  adjusted / origin independent) whose columns are fosc(D2) fosc(M2)
  fosc(Q2), a total and a component/total ratio;
- with ``DoSOC`` the SOC-corrected blocks repeat the same tables after the
  spin-orbit step; their fosc column is **weighted by the initial-state
  population** (column header ``fosc(D2) (*population)``) -- the reader keeps
  the raw value and records that semantic.

Transition labels are ``<root>-<mult><letter>`` (e.g. ``0-5A -> 1-5A``, the
SOC-corrected tables print ``5.0A``); the reader keeps the label verbatim
and parses the two root indices.  All energies are kept as printed (eV,
cm^-1, nm); oscillator strengths in atomic units as printed.
"""

from __future__ import annotations

import re
from typing import Any

__all__ = ["parse_rocis"]

_EXCITATION_ROW_RE = re.compile(
    r"^\s*(\d+)\s+(\d+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*$"
)
_SPEC_ROW_RE = re.compile(
    r"^\s*(\d+)-([\w.]+)\s*->\s*(\d+)-([\w.]+)\s+"
    r"([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.eE+]+)\s+(.*)$"
)
_BLOCK_TITLES = (
    ("dipole_length", "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"),
    ("dipole_velocity", "ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS"),
    ("combined_length", "ABSORPTION SPECTRUM COMBINED ELECTRIC DIPOLE + MAGNETIC DIPOLE + ELECTRIC QUADRUPOLE SPECTRUM"),
    ("combined_velocity", "ABSORPTION SPECTRUM COMBINED ELECTRIC DIPOLE + MAGNETIC DIPOLE + ELECTRIC QUADRUPOLE SPECTRUM (Velocity)"),
    ("combined_origin_adjusted", "ABSORPTION SPECTRUM COMBINED ELECTRIC DIPOLE + MAGNETIC DIPOLE + ELECTRIC QUADRUPOLE SPECTRUM (origin adjusted)"),
    ("combined_origin_independent", "ABSORPTION SPECTRUM COMBINED ELECTRIC DIPOLE + MAGNETIC DIPOLE + ELECTRIC QUADRUPOLE SPECTRUM (Origin Independent, Length)"),
    ("combined_origin_independent_velocity", "ABSORPTION SPECTRUM COMBINED ELECTRIC DIPOLE + MAGNETIC DIPOLE + ELECTRIC QUADRUPOLE SPECTRUM (Origin Independent, Velocity)"),
)


def _classify(title: str) -> str | None:
    """Block key for a spectrum-block title (exact-prefix match, longest first)."""
    for key, marker in sorted(_BLOCK_TITLES, key=lambda item: -len(item[1])):
        if title.strip().startswith(marker):
            return key
    return None


def _power_two(rest: str) -> bool:
    """D2/P2 blocks have three trailing vector components after the strength."""
    return len(rest.split()) >= 3


def parse_rocis(lines: list[str]) -> dict[str, Any]:
    """Excitation-energy table plus every absorption-spectrum block."""
    result: dict[str, Any] = {
        "present": False,
        "excitation_table": [],
        "spectra": {},
    }
    text = "\n".join(lines)
    if "ROCIS-EXCITATION SPECTRA" not in text:
        return result  # not a ROCIS run (the marker is printed by every one)
    # the excitation-energy table sits between its header and a blank/dashed line
    in_table = False
    for line in lines:
        if line.strip() == "Excitation energies":
            in_table = True
            continue
        if in_table:
            if (match := _EXCITATION_ROW_RE.match(line)) is not None:
                result["excitation_table"].append(
                    {
                        "root": int(match.group(1)),
                        "mult": int(match.group(2)),
                        "eh": float(match.group(3)),
                        "cm1": float(match.group(4)),
                        "ev": float(match.group(5)),
                    }
                )
            elif result["excitation_table"] and line.strip() and not line.strip()[0].isdigit():
                in_table = False
    # spectrum blocks: title line, then rows until the block's dashed footer
    index = 0
    soc_phase = False
    while index < len(lines):
        line = lines[index]
        title = line.strip()
        if title.startswith("SOC CORRECTED"):
            soc_phase = True
            title = title[len("SOC CORRECTED") :].strip()
        key = _classify(title) if title.startswith("ABSORPTION SPECTRUM") else None
        if key is not None:
            block_key = ("soc_" if soc_phase else "") + key
            rows: list[dict[str, Any]] = []
            probe = index + 1
            saw_data = False
            while probe < len(lines):
                candidate = lines[probe]
                if candidate.strip().startswith("-----") and saw_data:
                    break
                if (match := _SPEC_ROW_RE.match(candidate)) is not None:
                    saw_data = True
                    tail = match.group(9)
                    tokens = tail.split()
                    row: dict[str, Any] = {
                        "i_root": int(match.group(1)),
                        "i_label": match.group(2),
                        "j_root": int(match.group(3)),
                        "j_label": match.group(4),
                        "ev": float(match.group(5)),
                        "cm1": float(match.group(6)),
                        "nm": float(match.group(7)),
                        "fosc": float(match.group(8)),
                        "rest": [float(token) for token in tokens] if tokens else [],
                    }
                    rows.append(row)
                probe += 1
                if len(rows) > 4000:  # a runaway guard, not a format fact
                    break
            if rows:
                result["spectra"][block_key] = rows
                result["present"] = True
            index = probe
        else:
            index += 1
    # RIXS bookkeeping (the WARNING path is the measured refusal mode)
    for line in lines:
        if "zero number of intermediate and/or final states" in line:
            result["riqs_refused"] = True
        elif "Intermediate States:" in line:
            result["riqs_intermediate"] = int(line.split(":")[1])
        elif "Final States:" in line and "intermediate" not in line:
            result.setdefault("riqs_final", int(line.split(":")[1]))
    return result
