"""POLY_ANISO output (the ORCA polynuclear-magnetism driver).

``otool_poly_aniso`` (POLY_ANISO v1.0.0, Ungur & Chibotaru) is shipped with
ORCA and is called independently of it::

    otool_poly_aniso < poly_aniso.input > poly_aniso.output

The input starts with ``&POLY_ANISO``; the per-center ab initio data files
are the ``<job>.CASSCF.anisofile`` files produced by the SINGLE_ANISO runs,
copied under the names ``aniso_1.input``, ``aniso_2.input``, ... (the names
are mandatory -- a renamed file aborts the run in ``inquire_key_presence``).
Exchange constants are user input (PAIR/LIN1, LIN3, LIN9), not computed by
the program.

Measured output structure (ORCA 6.1.1, ``otool_poly_aniso`` v1.0.0 compiled
2025-12-04; two-center probe, 2026-09-29, fixtures ``poly_aniso/``):

- per-center blocks (``Center type N:`` -- note: the sample prints
  ``Center type 1`` for every data-file block, so the blocks are keyed by
  their ``Data file for this center`` field and by file order);
- the exchange block (coupled-state count, pair count, the wrapped
  ``Magnetic dipole-dipole coupling for all / interacting metal pairs`` and
  ``Decomposition of magnetic interaction into / irreducible tensor
  operators (ITO)`` labels, the per-pair J list);
- the run configuration echo (temperature interval, susceptibility tensors,
  plots);
- FIRST-ORDER ANISOTROPIC MAGNETIC COUPLING: per pair and per model
  (LINES-1, DIPOLE-DIPOLE), the full 3x3 interaction matrix plus its
  Isotropic / Symmetric / Anti-Symmetric decomposition with weight % (the
  matrix rows are 12-decimal E-notation; each matrix is split across the
  neighbouring table rows);
- the coupled-states table (Lines / dipole / total, absolute and relative);
- POPULATION ANALYSIS and the EXPECTATION VALUES tables;
- the chiT table (T, statistical sum Z, chiT for zJ=0, chiT, chi, 1/chi);
- the VAN VLECK SUSCEPTIBILITY TENSOR sequence (one 3x3 + main values +
  main axes per printed temperature);
- closing notes (skipped torque/magnetization, HAPPY LANDING,
  ``finished sucessfully``).
"""

from __future__ import annotations

import re
from typing import Any

__all__ = ["parse_poly_aniso"]

_CENTER_RE = re.compile(r"^\s*Center type (\d+):\s*$")
_FIELD_RE = re.compile(r"^\s*(.+?)\s+-{3,}\s*:\s*(.*?)\s*$")
_G_ROW_RE = re.compile(r"g([XYZ])\s*=\s*(-?[\d.]+)\s+(\w+)\s*:\s*([-0-9. ]+)$")
# interface-wrap labels: the value sits on the line after the label
_WRAPPED_FIELDS = (
    "Magnetic dipole-dipole coupling for all",
    "Decomposition of magnetic interaction into",
)
_PAIR_J_RE = re.compile(r"pair\s+(\d+)\s*:\s*centers\s+(\d+)\s*--\s*(\d+)\s*:\s*(-?[\d.]+)")
# a 12-decimal E-notation float, as in the interaction matrices
_E12 = r"-?\d\.\d{12}E[+-]\d{2}"
_MATRIX_ROW_RE = re.compile(rf"({_E12})\s+({_E12})\s+({_E12})")
_LABEL_ROW_RE = re.compile(
    r"^\s*(\d+)\s+(\d+)\s*--\s*(\d+)\s*\|\s*([\w-]+)\s*\|\s*(.+?)\s*\|\s*([\d.]+)\s*\|"
)
_CHIT_ROW_RE = re.compile(
    r"^\s*\|\s*([\d.]+)\s*\|\s*([\d.eE+-]+)\s*\|\s*([\d.eE+-]+)\s*\|\s*"
    r"([\d.eE+-]+)\s*\|\s*([\d.eE+-]+)\s*\|\s*([\d.eE+-]+)\s*\|\s*$"
)
_COUPLED_ROW_RE = re.compile(
    r"^\s*(\d+)\s*\|\s*(-?[\d.]+)\s*\|\s*(-?[\d.]+)\s*\|\s*(-?[\d.]+)\s*\|\s*(-?[\d.]+)\s*\|\s*$"
)
_VV_ROW_RE = re.compile(r"^\s*([\d.]*)\s*\|\s*([xyz])\s*\|(.*)$")


def _bare_floats(tokens: list[str]) -> bool:
    for token in tokens:
        try:
            float(token)
        except ValueError:
            return False
    return bool(tokens)


def parse_poly_aniso(lines: list[str]) -> dict[str, Any]:
    """Section-wise reader; every field keeps the printed value verbatim."""
    text = "\n".join(lines)
    if "POLY_ANISO Program" not in text:
        return {
            "present": False,
            "counts": {},
            "centers": [],
            "exchange": {},
            "config": {},
            "interaction": [],
            "coupled_states": [],
            "population": [],
            "expectation": {"MS": [], "LJ": []},
            "chiT": [],
            "van_vleck": [],
            "notes": [],
        }
    out: dict[str, Any] = {
        "present": True,
        "counts": {},
        "centers": [],
        "exchange": {"pairs_j": []},
        "config": {},
        "interaction": [],
        "coupled_states": [],
        "population": [],
        "expectation": {"MS": [], "LJ": []},
        "chiT": [],
        "van_vleck": [],
        "notes": [],
    }
    section: str | None = None
    center: dict[str, Any] | None = None
    pending_spectrum: list[float] = []
    pending_wrap: str | None = None
    expect_table: str | None = None
    last_expect_state: int | None = None
    vv_block: list[dict[str, Any]] = []
    matrix_rows: list[list[float]] = []
    labels: list[dict[str, Any]] = []

    def _flush_spectrum() -> None:
        nonlocal pending_spectrum
        if center is not None and pending_spectrum:
            center.setdefault("so_spectrum_cm1", []).extend(pending_spectrum)
        pending_spectrum = []

    def _flush_matrix() -> None:
        while len(matrix_rows) >= 3:
            triple = matrix_rows[:3]
            del matrix_rows[:3]
            entry = {"matrix": triple}
            if len(labels) > len(out["interaction"]):
                entry.update(labels[len(out["interaction"])])
            out["interaction"].append(entry)

    def _flush_vv() -> None:
        while len(vv_block) >= 3:
            rows = vv_block[:3]
            del vv_block[:3]
            entry = {
                "T": None,
                "tensor": [row["row"] for row in rows],
                "main_values": [row["main_value"] for row in rows],
                "main_axes": [row["main_axis"] for row in rows],
            }
            for row in rows:
                if row["T"] is not None:
                    entry["T"] = row["T"]
            out["van_vleck"].append(entry)

    for line in lines:
        stripped = line.strip()
        # ---- section switches ------------------------------------------------
        if "FIRST-ORDER ANISOTROPIC MAGNETIC COUPLING" in line:
            section = "interaction"
            continue
        if line.lstrip().startswith("Coupled|"):
            section = "coupled_states"
            continue
        if "POPULATION ANALYSIS" in line:
            section = "population"
            continue
        if stripped == "EXPECTATION VALUES":
            section = "expectation"
            expect_table = None
            continue
        if "CALCULATION OF THE MAGNETIC SUSCEPTIBILITY" in line:
            section = "chiT"
            continue
        if "VAN VLECK SUSCEPTIBILITY TENSOR" in line:
            section = "van_vleck"
            continue
        if "HAPPY" in line and "LANDING" in line:
            out["notes"].append("finished ok")
            continue
        if stripped.endswith("skipped by the user"):
            out["notes"].append(stripped)
            continue
        # ---- wrapped labels (value on the following line) --------------------
        if pending_wrap is not None:
            line = pending_wrap + " " + stripped
            pending_wrap = None
        elif stripped.startswith(_WRAPPED_FIELDS):
            pending_wrap = stripped
            continue
        # ---- centers ---------------------------------------------------------
        if (match := _CENTER_RE.match(line)) is not None:
            _flush_spectrum()
            center = {"type": int(match.group(1))}
            out["centers"].append(center)
            continue
        if (match := _G_ROW_RE.search(line)) is not None and center is not None:
            center.setdefault("g_tensor", {})[match.group(1)] = {
                "value": float(match.group(2)),
                "pseudospin": match.group(3),
                "main_axis": [float(value) for value in match.group(4).split()],
            }
            continue
        # ---- key : value fields ----------------------------------------------
        # two measured variants are not padded with dashes: the coordinates
        # line ("... in Angstrom  :  x y z") and the LINES-1 pair count
        if center is not None and section is None and "Cartesian coordinates" in line:
            tokens = line.split(":")[-1].split()
            if len(tokens) == 3 and _bare_floats(tokens):
                center["coords"] = [float(item) for item in tokens]
                continue
        if "Number of pairs for LINES-1 model" in line and ":" in line:
            out["exchange"]["lines1_pairs"] = int(line.split(":")[-1])
            continue
        if (match := _FIELD_RE.match(line)) is not None:
            key, value = match.group(1), match.group(2)
            if key.startswith("Number of independent centers"):
                out["counts"]["independent"] = int(value)
            elif key.startswith("Number of all centers"):
                out["counts"]["all"] = int(value)
            elif key.startswith("Temperature interval"):
                out["config"]["temperature_interval"] = value
            elif key.startswith("Magnetic susceptibility tensors"):
                out["config"]["susceptibility"] = value
            elif key.startswith("Automatic generation of PLOTS"):
                out["config"]["plots"] = value
            elif key.startswith("Number of exchange coupled states"):
                out["exchange"]["coupled_states"] = int(value)
            elif key.startswith("Number of interacting pairs"):
                out["exchange"]["pairs_total"] = int(value)
            elif key.startswith("Magnetic dipole-dipole coupling"):
                out["exchange"]["dipole_dipole"] = value
            elif key.startswith("Lines-1 model"):
                out["exchange"]["lines1"] = value
            elif key.startswith("Number of pairs for LINES-1"):
                out["exchange"]["lines1_pairs"] = int(value)
            elif key.startswith("Decomposition of magnetic interaction"):
                out["exchange"]["ito"] = value
            elif center is not None and section is None:
                if key.startswith("Data file for this center"):
                    center["data_file"] = value
                elif key.startswith("Cartesian coordinates"):
                    center["coords"] = [float(item) for item in value.split()]
                elif key.startswith("Number of centers of this kind"):
                    center["n_equivalent"] = int(value)
                elif key.startswith("Local basis for exchange interaction"):
                    center["exchange_basis"] = int(value)
                elif key.startswith("Number of spin-orbit states"):
                    center["n_so_states"] = int(value)
                elif key.startswith("Number of spin-free states"):
                    center["n_spin_free"] = int(value)
                elif key.startswith("Spin-orbit energy spectra"):
                    pending_spectrum = [float(item) for item in value.split()]
            continue
        if (match := _PAIR_J_RE.search(line)) is not None:
            out["exchange"]["pairs_j"].append(
                {
                    "pair": int(match.group(1)),
                    "centers": (int(match.group(2)), int(match.group(3))),
                    "J_cm1": float(match.group(4)),
                }
            )
            continue
        # ---- spectrum continuation -------------------------------------------
        if pending_spectrum and stripped and _bare_floats(stripped.split()):
            pending_spectrum.extend(float(item) for item in stripped.split())
            continue
        _flush_spectrum()
        # ---- interaction matrices --------------------------------------------
        if section == "interaction":
            if (match := _LABEL_ROW_RE.match(line)) is not None:
                labels.append(
                    {
                        "pair": (int(match.group(2)), int(match.group(3))),
                        "model": match.group(4),
                        "term": match.group(5),
                        "weight_percent": float(match.group(6)),
                    }
                )
                if (row := _MATRIX_ROW_RE.search(line)) is not None:
                    matrix_rows.append([float(row.group(i)) for i in (1, 2, 3)])
                    _flush_matrix()
                continue
            if (row := _MATRIX_ROW_RE.search(line)) is not None:
                matrix_rows.append([float(row.group(i)) for i in (1, 2, 3)])
                _flush_matrix()
                continue
            continue
        # ---- coupled states ---------------------------------------------------
        if section == "coupled_states":
            if (match := _COUPLED_ROW_RE.match(line)) is not None:
                out["coupled_states"].append(
                    {
                        "state": int(match.group(1)),
                        "lines_cm1": float(match.group(2)),
                        "dipole_cm1": float(match.group(3)),
                        "total_cm1": float(match.group(4)),
                        "relative_cm1": float(match.group(5)),
                    }
                )
            continue
        # ---- population -------------------------------------------------------
        if section == "population":
            parts = [item.strip() for item in line.split("|")]
            if len(parts) >= 4 and parts[0].isdigit() and parts[1].isdigit():
                out["population"].append(
                    {
                        "state": int(parts[0]),
                        "basis_set": int(parts[1]),
                        "weight": [_float_or_none(item) for item in parts[2:-1]],
                    }
                )
            continue
        # ---- expectation values -----------------------------------------------
        if section == "expectation":
            if "MAGNETIC MOMENT (M=-L-2S)" in line:
                expect_table = "MS"
                continue
            if "ORBITAL MOMENT (L)" in line:
                expect_table = "LJ"
                continue
            parts = [item.strip() for item in line.split("|")]
            if len(parts) >= 6 and parts[0].isdigit():
                last_expect_state = int(parts[0])
            if len(parts) >= 6 and _float_or_none(parts[1]) is not None and expect_table is not None:
                state = last_expect_state
                record = {
                    "state": state,
                    "site": int(parts[1]),
                    "vec1": [float(item) for item in parts[2].split()],
                    "norm1": _float_or_none(parts[3]),
                    "vec2": [float(item) for item in parts[4].split()],
                    "norm2": _float_or_none(parts[5]),
                }
                out["expectation"][expect_table].append(record)
            continue
        # ---- chiT -------------------------------------------------------------
        if section == "chiT" and (match := _CHIT_ROW_RE.match(line)) is not None:
            out["chiT"].append(
                {
                    "T": float(match.group(1)),
                    "Z": float(match.group(2)),
                    "chiT_zJ0": float(match.group(3)),
                    "chiT": float(match.group(4)),
                    "chi": float(match.group(5)),
                    "inv_chi": float(match.group(6)),
                }
            )
            continue
        # ---- van Vleck tensors --------------------------------------------------
        if section == "van_vleck" and (match := _VV_ROW_RE.match(line)) is not None:
            tensor_tokens = match.group(3).split("|")
            if len(tensor_tokens) <= 2:
                continue
            row_tokens = tensor_tokens[0].split()
            if len(row_tokens) != 3 or not _bare_floats(row_tokens):
                continue
            main_tokens = tensor_tokens[1].strip()
            axis_tokens = tensor_tokens[2].split() if len(tensor_tokens) > 2 else []
            vv_block.append(
                {
                    "T": _float_or_none(match.group(1)),
                    "row": [float(item) for item in row_tokens],
                    "main_value": _float_or_none(main_tokens.split(":")[-1]),
                    "main_axis": [float(item) for item in axis_tokens] if _bare_floats(axis_tokens) else [],
                }
            )
            _flush_vv()
            continue
    _flush_spectrum()
    _flush_matrix()
    return out


def _float_or_none(text: str) -> float | None:
    try:
        return float(text)
    except (TypeError, ValueError):
        return None
