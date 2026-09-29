"""AILFT embedded output (the ORCA ab initio ligand-field theory module).

The CASSCF module's AILFT driver (manual section 3.13.16; requested through
``ActOrbs``/``LFTCase``) prints, after the fit, one parameter block per
theoretical level (CASSCF, and NEVPT2 when requested):

- ``AILFT MATRIX ELEMENTS (<level>)`` -- the ligand-field one-electron
  matrix VLFT (a.u.), the Slater-Condon parameters (F0/F2/F4 for d shells,
  F0/F2/F4/F6 for f shells; three units each; "(fixed)" marks F0 from the
  raw two-electron integrals), the Racah parameters (A, B, C and C/B) and
  the ligand-field one-electron eigenfunctions, plus the ``*.lft.gbw``
  file name;
- ``COMPARISON OF AB INITIO AND LIGAND FIELD RESULTS`` -- per-block
  AI-root/LF-root tables (energies in eV, overlap S, Delta), the per-block
  and total RMS errors, the 95% confidence intervals and Pearson's
  correlation coefficient;
- ``SPIN ORBIT COUPLING (based on CASSCF orbitals)`` -- the AI-SOC and
  LF-SOC matrix elements, the SOC fit (a, b) and the fitted SOC constant
  zeta, summarised in the ``-----SOC-CONSTANTS-----`` block (ZETA_D for d
  shells, ZETA_F for f shells).

Measured on ORCA 6.1.1 (fixtures ``ailft/``; Ni(2+) d8 free ion and the
CeF3 f1 probe, 2026-09-29).  Only the printed values are kept; the reader
never refits.
"""

from __future__ import annotations

import re
from typing import Any

__all__ = ["parse_ailft"]

_HEADER_CFG_RE = re.compile(r"^\s*(\w+) configuration\s*$")
_HEADER_BLOCKS_RE = re.compile(r"^\s*(\d+) CI blocks\s*$")
_HEADER_MOS_RE = re.compile(r"^\s*MOs (\d+) to (\d+)\s*$")
_CENTER_RE = re.compile(r"Metal/Atom center is atom (\d+)")
_LEVEL_RE = re.compile(r"AILFT MATRIX ELEMENTS \((\w+)\)")
_VLFT_ROW_RE = re.compile(r"^\s*(\S+)\s+((?:-?\d+\.\d+\s+){1,8}-?\d+\.\d+)\s*$")
_SCP_RE = re.compile(
    r"^\s*(.+?)\s+=\s+(-?[\d.]+) a\.u\.\s+=\s+(-?[\d.]+) eV\s+=\s+(-?[\d.]+) cm\*\*-1( \(fixed\))?\s*$"
)
_CB_RE = re.compile(r"^\s*C/B\s+=\s+(-?[\d.]+)\s*$")
_EIGEN_ROW_RE = re.compile(r"^\s*(\d+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+((?:-?\d+\.\d+\s+){1,8}-?\d+\.\d+)\s*$")
_LFT_GBW_RE = re.compile(r"Ligand field orbitals were stored in (\S+)")
_BLOCK_RE = re.compile(r"^\s*Block\s+(\d+)\s*$")
_ROOT_RE = re.compile(
    r"^\s*AI-Root\s+(\d+):\s+E\(AI\)=\s+(-?[\d.]+) eV\s+->\s+LF-Root\s+(\d+):\s+(-?[\d.]+) eV\s+"
    r"S=\s+([\d.]+)\s+Delta=\s+(-?[\d.]+) eV\s*$"
)
_RMS_RE = re.compile(r"RMS error for this block = \s+(-?[\d.]+) eV = \s+(-?[\d.]+) cm\*\*-1")
_TOTAL_RMS_RE = re.compile(r"Total RMS error g= \s+(-?[\d.]+) eV = \s+(-?[\d.]+) cm\*\*-1")
_PEARSON_RE = re.compile(r"Pearson's correlation coefficient = \s+(-?[\d.]+)")
_REF_AI_RE = re.compile(r"Reference energy AI\s+= \s+(-?[\d.]+) au")
_SOC_ZETA_RE = re.compile(r"SOC constant zeta = \s+(-?[\d.]+) eV = \s+(-?[\d.]+) cm\*\*-1")
_ZETA_SUMMARY_RE = re.compile(r"^(ZETA_[A-Z]) = (-?[\d.]+)\s*$")
_SOC_HEADER_RE = re.compile(r"SPIN ORBIT COUPLING \(based on (\w+) orbitals\)")
_SOC_RMS_RE = re.compile(r"RMS error of nonzero matrix elements = \s+(-?[\d.]+) cm\*\*-1")


def _numbers(text: str) -> list[float]:
    return [float(token) for token in text.split()]


def parse_ailft(lines: list[str]) -> dict[str, Any]:
    """The AILFT parameter blocks, fit statistics and SOC constants."""
    text = "\n".join(lines)
    if "AB INITIO LIGAND FIELD THEORY" not in text:
        return {"present": False, "header": {}, "levels": {}, "fit": {}, "soc": None}
    header: dict[str, Any] = {}
    levels: dict[str, dict[str, Any]] = {}
    fit: dict[str, dict[str, Any]] = {}
    soc: dict[str, Any] | None = None
    current: str | None = None
    table: str | None = None  # "vlft" or "eigenfunctions"
    racah_mode = False
    in_fit = False
    in_soc = False
    block: dict[str, Any] | None = None
    for line in lines:
        stripped = line.strip()
        if (match := _HEADER_CFG_RE.match(line)) is not None and "configuration" in line:
            header.setdefault("configuration", match.group(1))
            continue
        if (match := _HEADER_BLOCKS_RE.match(line)) is not None:
            header.setdefault("ci_blocks", int(match.group(1)))
            continue
        if (match := _HEADER_MOS_RE.match(line)) is not None and "MOs" in line:
            header.setdefault("mo_range", (int(match.group(1)), int(match.group(2))))
            continue
        if (match := _CENTER_RE.search(line)) is not None:
            header.setdefault("center_atom", int(match.group(1)))
            continue
        if (match := _LEVEL_RE.search(line)) is not None:
            current = match.group(1).lower()
            levels.setdefault(
                current,
                {
                    "vlft": None,
                    "slater_condon": [],
                    "racah": {},
                    "eigenfunctions": None,
                    "lft_gbw": None,
                },
            )
            in_fit = False
            table = None
            continue
        if current is None:
            continue
        if (match := _LFT_GBW_RE.search(line)) is not None:
            levels[current]["lft_gbw"] = match.group(1)
            continue
        if "Ligand field one-electron matrix VLFT" in line:
            table = "vlft"
            continue
        if "ligand field one electron eigenfunctions" in line:
            table = "eigenfunctions"
            continue
        if line.lstrip().startswith("Orbital"):
            # the d-shell VLFT table announces itself above; the f-shell table
            # goes straight to its header; unrelated log lines also start with
            # "Orbital" (e.g. "Orbital improvement steps"), so the header is
            # recognised by its second token: "Energy" for the eigenfunction
            # table, an orbital label for the VLFT tables
            tokens = stripped.split()
            second = tokens[1] if len(tokens) > 1 else ""
            if second == "Energy":
                table = "eigenfunctions"
            elif second and (
                (second[0] in "df" and second[1:2].isdigit())
                or second in ("dz2", "dxz", "dyz", "dxy", "dx2-y2", "dx2y2")
            ):
                table = "vlft"
            else:
                continue  # an unrelated line that merely starts with "Orbital"
            # fall through to the table handlers below
        if "Racah Parameters" in line:
            racah_mode = True
            continue
        if racah_mode and stripped.startswith("----") and levels[current]["racah"]:
            racah_mode = False
            continue
        if not in_fit and (match := _SCP_RE.match(line)) is not None:
            entry = {
                "label": match.group(1).strip(),
                "au": float(match.group(2)),
                "ev": float(match.group(3)),
                "cm1": float(match.group(4)),
                "fixed": match.group(5) is not None,
            }
            if racah_mode:
                levels[current]["racah"][entry["label"]] = entry
            else:
                levels[current]["slater_condon"].append(entry)
            continue
        if (match := _CB_RE.match(line)) is not None:
            levels[current]["racah"]["C_over_B"] = float(match.group(1))
            continue
        if table == "vlft":
            if line.lstrip().startswith("Orbital"):
                if levels[current].get("vlft") is None:
                    levels[current]["vlft"] = {"labels": [], "matrix": []}
                if not levels[current]["vlft"]["matrix"]:
                    levels[current]["vlft"]["labels"] = stripped.split()
                continue
            if (match := _VLFT_ROW_RE.match(line)) is not None:
                vlft = levels[current]["vlft"] or {"labels": None, "matrix": []}
                vlft["matrix"].append(
                    [match.group(1)] + _numbers(match.group(2))
                )
                levels[current]["vlft"] = vlft
                continue
            if stripped and levels[current].get("vlft") and levels[current]["vlft"]["matrix"]:
                table = None
        if table == "eigenfunctions":
            if line.lstrip().startswith("Orbital"):
                levels[current]["eigenfunctions"] = []
                continue
            if (match := _EIGEN_ROW_RE.match(line)) is not None and levels[current]["eigenfunctions"] is not None:
                levels[current]["eigenfunctions"].append(
                    {
                        "index": int(match.group(1)),
                        "ev": float(match.group(2)),
                        "cm1": float(match.group(3)),
                        "coeffs": _numbers(match.group(4)),
                    }
                )
                continue
            if stripped and levels[current]["eigenfunctions"]:
                table = None
        if (match := _REF_AI_RE.search(line)) is not None:
            entry = fit.setdefault(current, {"blocks": [], "total_rms_ev": None, "total_rms_cm1": None, "pearson": None, "reference_ai_au": None})
            entry["reference_ai_au"] = float(match.group(1))
            continue
        if "COMPARISON OF AB INITIO AND LIGAND FIELD RESULTS" in line:
            in_fit = True
            fit.setdefault(current, {"blocks": [], "total_rms_ev": None, "total_rms_cm1": None, "pearson": None, "reference_ai_au": None})
            block = None
            continue
        if in_fit and current in fit:
            entry = fit[current]
            if _SCP_RE.match(line) is not None:
                continue  # the confidence-interval rows reuse the SCP line shape
            if (match := _BLOCK_RE.match(line)) is not None:
                block = {"index": int(match.group(1)), "roots": [], "rms_ev": None, "rms_cm1": None}
                entry["blocks"].append(block)
                continue
            if block is not None and (match := _ROOT_RE.match(line)) is not None:
                block["roots"].append(
                    {
                        "ai_root": int(match.group(1)),
                        "ai_ev": float(match.group(2)),
                        "lf_root": int(match.group(3)),
                        "lf_ev": float(match.group(4)),
                        "overlap": float(match.group(5)),
                        "delta_ev": float(match.group(6)),
                    }
                )
                continue
            if block is not None and (match := _RMS_RE.search(line)) is not None:
                block["rms_ev"] = float(match.group(1))
                block["rms_cm1"] = float(match.group(2))
                continue
            if (match := _TOTAL_RMS_RE.search(line)) is not None:
                entry["total_rms_ev"] = float(match.group(1))
                entry["total_rms_cm1"] = float(match.group(2))
                continue
            if (match := _PEARSON_RE.search(line)) is not None:
                entry["pearson"] = float(match.group(1))
                in_fit = False
                continue
        if (match := _SOC_HEADER_RE.search(line)) is not None:
            in_soc = True
            soc = soc or {"bases": [], "zeta_ev": None, "zeta_cm1": None, "summary": {}, "rms_cm1": None}
            soc["bases"].append(match.group(1).lower())
            continue
        if in_soc and soc is not None:
            if (match := _SOC_ZETA_RE.search(line)) is not None:
                soc["zeta_ev"] = float(match.group(1))
                soc["zeta_cm1"] = float(match.group(2))
                continue
            if (match := _SOC_RMS_RE.search(line)) is not None:
                soc["rms_cm1"] = float(match.group(1))
                continue
            if (match := _ZETA_SUMMARY_RE.match(stripped)) is not None:
                soc["summary"][match.group(1)] = float(match.group(2))
                continue
    return {"present": True, "header": header, "levels": levels, "fit": fit, "soc": soc}
