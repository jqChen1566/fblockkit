"""SINGLE_ANISO embedded output (the ORCA interface to otool_single_aniso).

The ORCA ``ANISO`` sub-block inside ``%casscf`` (manual §5.32) calls the
SINGLE_ANISO utility on the SOC data file and prints its output in full.  The
format below was measured on ORCA 6.1.1 (fixtures ``single_aniso/``):

- **One segment per spin-orbit pass.**  With ``DoSSC true`` the output
  contains *two* complete segments (the first from the SOC-only QDPT data,
  the second after the spin-spin coupling pass; each segment ends with the
  ``SINGLE_ANISO finished sucessfully`` banner).  Parse segments in order and
  compare their spectra before choosing one.
- **``MLTP`` must be given explicitly**: without it the header reports
  ``Computation of the EPR g-tensor ... NOT INCLUDED`` and no pseudospin
  group block is printed at all (measured; the manual's "default = ground
  term multiplicity" wording notwithstanding).
- Each pseudospin group block carries the g-tensor table (principal values +
  main magnetic axes), the magnetic-moment ITO decomposition (not parsed
  here), the ZFS matrix and D-tensor tables when the pseudospin dimension
  is > 2 (a Kramers doublet prints no ZFS/D block) and the angular-moment
  expectation values along the main magnetic axes.
- The ``AB INITIO BLOCKING BARRIER`` block (``UBAR true``) prints the Zeeman
  eigenstates (Kramers: ``1+ 1- E`` columns; integer pseudospin: ``1+ 0 1- E``)
  and the magnetic-moment matrix elements between them, intra-group
  (opposite magnetisation) and inter-group ("neighbouring multiplets"),
  averaged per the printed formula ``(|mu_X| + |mu_Y| + |mu_Z|)/3``.  The
  ``in cases with even number of electrons ... check the tunnelling
  splitting instead`` sentence is a fixed template warning: it is printed
  for odd-electron systems as well and is recorded, not interpreted.

Numbers are kept in the file's own units: cm^-1 for energies, Bohr
magnetons for the matrix elements.
"""

from __future__ import annotations

import re
from fractions import Fraction
from typing import Any

__all__ = ["parse_segments", "doublets_payload"]

_FLOAT = r"-?\d+\.\d+(?:E[+-]?\d+)?"
_NUMBER = rf"(-?[\d.]+(?:E[+-]?\d+)?)"

_START_MARK = "SINGLE_ANISO Program"
_END_MARK = "SINGLE_ANISO finished sucessfully"

_MULT_RE = re.compile(
    r"CALCULATION OF PSEUDOSPIN HAMILTONIAN TENSORS FOR THE MULTIPLET\s+(\d+)\s*"
    rf"\(\s*effective S\s*=\s*([\d.]+?)(?:/(\d+))?\s*\)"
)
_SPIN_ORBIT_STATE_RE = re.compile(rf"spin-orbit state \d+\. energy\(\d+\) =\s*{_NUMBER} cm-1")
_G_ROW_RE = re.compile(
    rf"g([XYZ])\s*=\s*{_NUMBER}\s*\|\s*([XYZ])m\s*\|\s*{_NUMBER}\s+{_NUMBER}\s+{_NUMBER}\s*\|"
)
_SIGN_RE = re.compile(r"The sign of the product gX \* gY \* gZ for multiplet (\d+):\s*([<>])\s*0")
_D_ROW_RE = re.compile(
    rf"D([xyz])\s*=\s*{_NUMBER}\s*\|\s*([XYZ])a\s*\|\s*{_NUMBER}\s+{_NUMBER}\s+{_NUMBER}\s*\|\s*"
    rf"{_NUMBER}\s+{_NUMBER}\s+{_NUMBER}\s*\|"
)
_AM_RE = re.compile(
    rf"<1\|\s*mu_([XYZ])\s*\|1>\s*\|\s*{_NUMBER}\s*\|\s*{_NUMBER}\s*\|\s*{_NUMBER}\s*\|"
)
_MATRIX_ELEMENT_RE = re.compile(
    rf"^\s*(\d+)\.\s*\|\s*<\s*(\d+)\.([0-9][+-]?)\s*\|\s*mu_([XYZ])\s*\|\s*"
    rf"(\d+)\.([0-9][+-]?)\s*>\s*\|\s*"
    rf"(-?[\d.]+(?:E[+-]?\d+)?)\s+(-?[\d.]+(?:E[+-]?\d+)?)\s*\|\s*"
    rf"(-?[\d.]+(?:E[+-]?\d+)?)?\s*\|"
)
_ZEEMAN_ROW_RE = re.compile(r"^\s*\d+\.\s*\|.*\|\s*$")
_EVEN_NOTE = "in cases with even number of electrons"
_SPECTRUM_COLON_RE = re.compile(rf"^\s*:\s*((?:{_FLOAT}\s*)+)$")


def _parse_spin(match: re.Match[str]) -> float:
    if match.group(3) is not None:
        return float(match.group(2)) / int(match.group(3))
    return float(match.group(2))


def _collect_spectrum(lines: list[str], marker: str) -> list[float]:
    """Floats on the marker line after its colon plus the following colon lines."""
    for index, line in enumerate(lines):
        if marker not in line:
            continue
        values: list[float] = []
        if ":" in line:
            values.extend(float(token) for token in line.split(":", 1)[1].split() if token)
        for probe in range(index + 1, len(lines)):
            if (match := _SPECTRUM_COLON_RE.match(lines[probe])) is None:
                break
            values.extend(float(token) for token in match.group(1).split())
        return values
    return []


def _parse_header(lines: list[str]) -> dict[str, Any]:
    header: dict[str, Any] = {
        "datafile": None,
        "n_soc_states": None,
        "n_spin_free_states": None,
        "g_tensor_included": None,
        "ubar_included": None,
    }
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("Data file for this center"):
            header["datafile"] = stripped.rsplit(":", 1)[1].strip()
        elif stripped.startswith("Number of spin-orbit states"):
            header["n_soc_states"] = int(stripped.rsplit(":", 1)[1])
        elif stripped.startswith("Number of spin-free states"):
            header["n_spin_free_states"] = int(stripped.rsplit(":", 1)[1])
        elif stripped.startswith("Computation of the EPR g-tensor"):
            header["g_tensor_included"] = "NOT INCLUDED" not in stripped
        elif stripped.startswith("Estimation of blocking barrier"):
            header["ubar_included"] = "INCLUDED" in stripped and "NOT" not in stripped
    return header


def _parse_group(lines: list[str], start: int, stop: int, index: int, spin: float) -> dict[str, Any]:
    group: dict[str, Any] = {
        "index": index,
        "spin": spin,
        "soc_state_energies_cm1": [],
        "g_values": None,
        "g_axes": None,
        "g_sign_product": None,
        "d_values": None,
        "d_axes_main": None,
        "d_axes_cartesian": None,
        "d_cm1": None,
        "e_cm1": None,
        "angular_moments_z": None,
    }
    g_values: dict[str, float] = {}
    g_axes: dict[str, tuple[float, float, float]] = {}
    d_values: dict[str, float] = {}
    d_axes_main: dict[str, tuple[float, float, float]] = {}
    d_axes_cartesian: dict[str, tuple[float, float, float]] = {}
    for position in range(start, stop):
        line = lines[position]
        if (match := _SPIN_ORBIT_STATE_RE.search(line)) is not None:
            group["soc_state_energies_cm1"].append(float(match.group(1)))
        elif (match := _G_ROW_RE.search(line)) is not None:
            g_values[match.group(1)] = float(match.group(2))
            g_axes[match.group(3)] = (float(match.group(4)), float(match.group(5)), float(match.group(6)))
        elif (match := _SIGN_RE.search(line)) is not None:
            group["g_sign_product"] = -1 if match.group(2) == "<" else 1
        elif (match := _D_ROW_RE.search(line)) is not None:
            d_values[match.group(1)] = float(match.group(2))
            d_axes_main[match.group(3)] = (float(match.group(4)), float(match.group(5)), float(match.group(6)))
            d_axes_cartesian[match.group(3)] = (float(match.group(7)), float(match.group(8)), float(match.group(9)))
        elif (
            line.strip().startswith("D =")
            and position >= 1
            and "Anisotropy parameters" in lines[position - 1]
        ):
            group["d_cm1"] = float(line.split("=", 1)[1])
        elif (
            line.strip().startswith("E =")
            and group["d_cm1"] is not None
            and any("Anisotropy parameters" in probe for probe in lines[max(0, position - 2) : position])
        ):
            group["e_cm1"] = float(line.split("=", 1)[1])
        elif (match := _AM_RE.search(line)) is not None:
            if match.group(1) == "Z":
                group["angular_moments_z"] = tuple(float(match.group(i)) for i in (2, 3, 4))
    if g_values:
        group["g_values"] = (g_values.get("X"), g_values.get("Y"), g_values.get("Z"))
        group["g_axes"] = tuple(g_axes.get(axis) for axis in ("X", "Y", "Z"))
    if d_values:
        group["d_values"] = (d_values.get("x"), d_values.get("y"), d_values.get("z"))
        group["d_axes_main"] = tuple(d_axes_main.get(axis) for axis in ("X", "Y", "Z"))
        group["d_axes_cartesian"] = tuple(d_axes_cartesian.get(axis) for axis in ("X", "Y", "Z"))
    return group


def _parse_ubar(lines: list[str], start: int, stop: int) -> dict[str, Any]:
    ubar: dict[str, Any] = {
        "present": False,
        "zeeman_states": [],
        "matrix_elements": [],
        "even_electron_note": False,
    }
    for position in range(start, stop):
        line = lines[position]
        if _EVEN_NOTE in line:
            ubar["even_electron_note"] = True
        if (match := _MATRIX_ELEMENT_RE.match(line)) is not None:
            average = match.group(9)
            ubar["matrix_elements"].append(
                {
                    "source_mult": int(match.group(1)),
                    "source_state": match.group(3),
                    "component": match.group(4),
                    "target_mult": int(match.group(5)),
                    "target_state": match.group(6),
                    "delta_mult": int(match.group(5)) - int(match.group(1)),
                    "real": float(match.group(7)),
                    "imag": float(match.group(8)),
                    "average": float(average) if average else None,
                }
            )
            ubar["present"] = True
        elif _ZEEMAN_ROW_RE.match(line) is not None and "<" not in line:
            cells = [cell.strip() for cell in line.split("|") if cell.strip()]
            if len(cells) >= 3 and all(re.fullmatch(_FLOAT, cell) for cell in cells[1:]):
                numeric = [float(cell) for cell in cells[1:]]
                entry = {"mult": int(cells[0].rstrip("."))}
                if len(numeric) == 3:  # Kramers: 1+ 1- E
                    entry.update(m_plus=numeric[0], m_zero=None, m_minus=numeric[1], energy_cm1=numeric[2])
                    entry["columns"] = 3
                elif len(numeric) == 4:  # integer pseudospin: 1+ 0 1- E
                    entry.update(
                        m_plus=numeric[0], m_zero=numeric[1], m_minus=numeric[2], energy_cm1=numeric[3]
                    )
                    entry["columns"] = 4
                else:
                    continue  # an unexpected column count: leave the row out
                ubar["zeeman_states"].append(entry)
                ubar["present"] = True
    return ubar


def _segment_bounds(lines: list[str]) -> list[tuple[int, int]]:
    bounds: list[tuple[int, int]] = []
    start = None
    for index, line in enumerate(lines):
        if _START_MARK in line and start is None:
            start = index
        elif _END_MARK in line and start is not None:
            bounds.append((start, index + 1))
            start = None
    return bounds


def parse_segments(lines: list[str]) -> dict[str, Any]:
    """All SINGLE_ANISO segments of an ORCA output, in print order."""
    segments: list[dict[str, Any]] = []
    for start, stop in _segment_bounds(lines):
        chunk = lines[start:stop]
        segment: dict[str, Any] = _parse_header(chunk)
        segment["soc_spectrum_cm1"] = _collect_spectrum(chunk, "Spin-orbit energy spectra")
        segment["spin_free_spectrum_cm1"] = _collect_spectrum(chunk, "Spin-free energy spectra")
        # group blocks: from each MULTIPLET header to the next block boundary
        group_heads: list[tuple[int, int, float]] = []
        for offset, line in enumerate(chunk):
            if (match := _MULT_RE.search(line)) is not None:
                group_heads.append((offset, int(match.group(1)), _parse_spin(match)))
        boundaries = [head[0] for head in group_heads]
        for position, (offset, index, spin) in enumerate(group_heads):
            if position + 1 < len(group_heads):
                stop_at = group_heads[position + 1][0]
            else:
                stop_at = len(chunk)
                for probe in range(offset, len(chunk)):
                    if "AB INITIO BLOCKING BARRIER" in chunk[probe]:
                        stop_at = probe
                        break
            segment.setdefault("groups", []).append(
                _parse_group(chunk, offset, stop_at, index, spin)
            )
        ubar_start = None
        for offset, line in enumerate(chunk):
            if "AB INITIO BLOCKING BARRIER" in line:
                ubar_start = offset
                break
        segment["ubar"] = (
            _parse_ubar(chunk, ubar_start, len(chunk)) if ubar_start is not None else None
        )
        text = "\n".join(chunk)
        segment["chi_present"] = "CALCULATION OF THE MAGNETIC SUSCEPTIBILITY" in text
        segment["magnetization_present"] = "CALCULATION OF THE MOLAR MAGNETIZATION" in text
        segment.setdefault("groups", [])
        segments.append(segment)
    return {"present": bool(segments), "segments": segments}


def _spin_text(spin: float) -> str:
    fraction = Fraction(spin).limit_denominator(2)
    if fraction.denominator == 1:
        return str(fraction.numerator)
    return f"{fraction.numerator}/{fraction.denominator}"


def doublets_payload(
    segments: list[dict[str, Any]],
    *,
    segment_index: int = 0,
    system: str | None = None,
) -> dict:
    """Build the menu-16 JSON table from one segment's pseudospin groups.

    Rows mirror the OpenMolcas payload (parsers/molcas_single_aniso.py): the
    g values ascending, the axis of the largest g value, and the group's
    lowest spin-orbit energy (cm^-1).  Groups without a g table (a run whose
    MLTP did not cover them) are skipped; a segment with no g-carrying
    group is refused.
    """
    if not segments or not 0 <= segment_index < len(segments):
        raise ValueError(
            f"no SINGLE_ANISO segment at index {segment_index}. Next step: check "
            "that the ORCA %casscf ANISO block ran (MLTP must be given "
            "explicitly)."
        )
    rows = []
    for group in segments[segment_index].get("groups") or []:
        g_values = group.get("g_values")
        g_axes = group.get("g_axes")
        if not g_values or not g_axes or not all(
            value is not None for value in (*g_values, *g_axes)
        ):
            continue
        triples = sorted(zip(g_values, g_axes), key=lambda item: item[0])
        row = {
            "label": (
                f"multiplet {group['index']} (effective S = {_spin_text(group['spin'])})"
            ),
            "g": [value for value, _ in triples],
            "axis3": list(triples[-1][1]),
        }
        energies = group.get("soc_state_energies_cm1") or []
        if energies:
            row["energy"] = min(energies)
        rows.append(row)
    if not rows:
        raise ValueError(
            "the SINGLE_ANISO segment carries no pseudospin group with a g "
            "tensor; a run without MLTP prints no g analysis. Next step: give "
            "MLTP explicitly in the %casscf ANISO block and rerun."
        )
    return {
        "system": system or "ORCA SINGLE_ANISO output",
        "reference": 0,
        "doublets": rows,
    }
