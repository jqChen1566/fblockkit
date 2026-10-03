"""The MOKIT/automr run output (plan item 6.1; the reader behind the automr menu).

MOKIT (Jingxiang Zou, Molecular Orbital Kit; Apache-2.0) drives a black-box
multireference workflow -- HF, UNO rotation, GVB pairs, CASCI/CASSCF, and
the dynamic-correlation family -- through one ``automr x.gjf`` call.  This
module reads the run output back; every convention below is **measured** on
MOKIT 1.2.8 (fixtures ``fixtures/mokit/``, 2026-09-30; H2O CASSCF probe with
``mokit{GVB_prog=Gaussian}`` on Gaussian 16 + PySCF):

- the identification banner ``Output of AutoMR of MOKIT`` and the version
  line ``Version: 1.2.8 (2026-Aug-18)``;
- the program-path block (``MOKIT_ROOT``, ``gau_path``, ``gms_path``,
  ``orca_path`` ...; unset entries read ``NOT FOUND``) and the run settings
  line ``memory = 4GB, nproc = 2, method/basis = casscf/cc-pvdz``;
- the merged ``mokit{}`` keywords (``gvb_prog=gaussian``);
- the strategy table(s): ``No. Strategy = N`` plus flag rows
  (``UNO = T   GVB = T   CASSCF = T``); the table is reprinted whenever the
  strategy updates, so the **last** values are the effective ones;
- the stage markers ``Enter subroutine <name>...`` in execution order;
- the GVB/CASSCF stage lines ``GVB(4) using program gaussian`` and
  ``CASSCF(4e,4o) using program pyscf`` (the automatically determined active
  space is the one on the CASSCF line);
- the energy chain rows ``E(UHF) = -75.82292560 a.u., <S**2>= 1.066`` /
  ``E(GVB) = ...`` / ``E(CASCI) = ...`` / ``E(CASSCF) = ...`` in order;
- the ``Radical index`` tables (printed once per stage; the two heading
  forms ``biradical character   (1-2t/(1+t^2)) y0=`` and ``biradical
  character y0 = n_LUNO =`` both occur);
- the closing line ``Normal termination of AutoMR at <timestamp>``.

Side products (``*_rhf.fch``, ``*_uno_asrot.fch``, ``*2gvb4_s.fch``,
``*_CASSCF_NO.fch``, stage ``.dat`` files) are echoed by this reader from the
run text; the ``.fch`` files themselves are parsed by
``parsers/mokit_fch.py`` (menu 45's side-product section).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

__all__ = [
    "MokitError",
    "MokitRun",
    "read_mokit_run",
    "read_mokit_run_file",
]


class MokitError(ValueError):
    """Not a MOKIT automr output (or the banner is missing)."""


@dataclass(frozen=True)
class MokitRun:
    """The measured content of one automr stdout capture."""

    version: str
    built: str
    program_paths: dict[str, str]
    memory: str
    nproc: str
    method_basis: str
    keywords: str
    strategy_number: int | None
    strategy_flags: dict[str, bool]
    stages: tuple[str, ...]
    energies: tuple[tuple[str, float], ...]
    s2_values: tuple[tuple[str, float], ...]
    gvb_order: int | None
    gvb_program: str | None
    active_space: tuple[int, int] | None   # (electrons, orbitals) from the CASSCF line
    casscf_program: str | None
    radical_groups: tuple[tuple[tuple[str, float], ...], ...]
    terminated: bool
    terminated_at: str | None
    side_products: tuple[str, ...] = field(default=())


_BANNER = "Output of AutoMR of MOKIT"
_VERSION_RE = re.compile(r"Version:\s*(\S+)\s+\((\S+)\)")
_PATH_RE = re.compile(r"^(\w+_path)\s*=\s*(\S.*?)\s*$")
_SETTINGS_RE = re.compile(r"memory =\s*(\S+?), nproc =\s*(\d+), method/basis = (\S+)")
_KEYWORDS_HEADER = "Keywords in MOKIT{} are merged and shown as follows:"
_STRATEGY_RE = re.compile(r"No\. Strategy = (\d+)")
_FLAG_RE = re.compile(r"([A-Za-z_][A-Za-z0-9_]*)\s*=\s*([TF])\b")
_STAGE_RE = re.compile(r"Enter subroutine (\w+)")
_GVB_RE = re.compile(r"GVB\((\d+)\) using program (\w+)")
_CASSCF_RE = re.compile(r"CASSCF\((\d+)e,(\d+)o\) using program (\w+)")
_ENERGY_RE = re.compile(r"E\((\w+)\)\s*=\s*(-?\d+\.\d+) a\.u\.")
_S2_RE = re.compile(r"E\((\w+)\)\s*=\s*-?\d+\.\d+ a\.u\., <S\*\*2>=\s*([\d.]+)")
_RADICAL_TITLE_RE = re.compile(r"^-+\s*Radical index\s*-+$")
_RADICAL_ROWS = (
    ("biradical", re.compile(r"biradical character.*?y0\s*=.*?([\d.]+)\s*$")),
    ("tetraradical", re.compile(r"tetraradical character.*?y1\s*=.*?([\d.]+)\s*$")),
    ("yamaguchi_unpaired", re.compile(r"Yamaguchi's unpaired electrons.*?:\s*([\d.]+)\s*$")),
    ("head_gordon_min", re.compile(r"Head-Gordon's unpaired electrons\(sum_n min.*?:\s*([\d.]+)\s*$")),
    ("head_gordon_squared", re.compile(r"Head-Gordon's unpaired electrons\(sum_n \(n\(2-n\)\)\^2.*?:\s*([\d.]+)\s*$")),
)
_TERMINATION_RE = re.compile(r"Normal termination of AutoMR at (.+?)\s*$")
_SIDE_SUFFIXES = (".fch", ".fchk", ".dat", ".gjf")


def read_mokit_run(text: str) -> MokitRun:
    """Parse one automr output; raise MokitError without the identification banner."""
    if _BANNER not in text:
        raise MokitError(
            "no 'Output of AutoMR of MOKIT' banner found: this does not look like a "
            "MOKIT automr output. Next step: give the .out file captured from "
            "'automr x.gjf > x.out'."
        )
    lines = text.splitlines()
    version, built = "", ""
    if (found := _VERSION_RE.search(text)) is not None:
        version, built = found.group(1), found.group(2)

    program_paths: dict[str, str] = {}
    for line in lines:
        if (found := _PATH_RE.match(line)) is not None:
            program_paths[found.group(1)] = found.group(2)
    memory = nproc = method_basis = ""
    if (found := _SETTINGS_RE.search(text)) is not None:
        memory, nproc, method_basis = found.group(1), found.group(2), found.group(3)

    keywords = ""
    for position, line in enumerate(lines):
        if line.strip() == _KEYWORDS_HEADER:
            for following in lines[position + 1 :]:
                if following.strip():
                    keywords = following.strip()
                    break
            break

    strategy_number: int | None = None
    strategy_flags: dict[str, bool] = {}
    for line in lines:
        if (found := _STRATEGY_RE.search(line)) is not None:
            strategy_number = int(found.group(1))
        for name, value in _FLAG_RE.findall(line):
            strategy_flags[name] = value == "T"

    stages = tuple(_STAGE_RE.findall(text))

    energies: list[tuple[str, float]] = []
    for label, value in _ENERGY_RE.findall(text):
        energies.append((label, float(value)))
    s2_values = tuple((label, float(value)) for label, value in _S2_RE.findall(text))

    gvb_order: int | None = None
    gvb_program: str | None = None
    if (found := _GVB_RE.search(text)) is not None:
        gvb_order = int(found.group(1))
        gvb_program = found.group(2)
    active_space: tuple[int, int] | None = None
    casscf_program: str | None = None
    if (found := _CASSCF_RE.search(text)) is not None:
        active_space = (int(found.group(1)), int(found.group(2)))
        casscf_program = found.group(3)

    radical_groups = _radical_groups(lines)

    terminated = False
    terminated_at: str | None = None
    if (found := _TERMINATION_RE.search(text)) is not None:
        terminated = True
        terminated_at = found.group(1)

    return MokitRun(
        version=version,
        built=built,
        program_paths=program_paths,
        memory=memory,
        nproc=nproc,
        method_basis=method_basis,
        keywords=keywords,
        strategy_number=strategy_number,
        strategy_flags=strategy_flags,
        stages=stages,
        energies=tuple(energies),
        s2_values=s2_values,
        gvb_order=gvb_order,
        gvb_program=gvb_program,
        active_space=active_space,
        casscf_program=casscf_program,
        radical_groups=radical_groups,
        terminated=terminated,
        terminated_at=terminated_at,
        side_products=_side_products(text),
    )


def _radical_groups(lines: list[str]) -> tuple[tuple[tuple[str, float], ...], ...]:
    """The Radical index tables, one group per heading block."""
    groups: list[tuple[tuple[str, float], ...]] = []
    position = 0
    while position < len(lines):
        if _RADICAL_TITLE_RE.match(lines[position].strip()):
            entries: list[tuple[str, float]] = []
            for following in lines[position + 1 :]:
                if following.strip() == "":
                    break
                for name, pattern in _RADICAL_ROWS:
                    if (found := pattern.search(following)) is not None:
                        entries.append((name, float(found.group(1))))
            if entries:
                groups.append(tuple(entries))
        position += 1
    return tuple(groups)


def _side_products(text: str) -> tuple[str, ...]:
    """The side-product file names echoed by the run (audit trail), in first-seen order.

    The echo lines are the ``$<utility> <args...>`` audit trail; file names may sit
    anywhere in the argument list (``$dat2fch x.dat x.fch -gvb 4``), so every token
    of every echoed command is inspected.
    """
    seen: list[str] = []
    for line in text.splitlines():
        strip = line.strip()
        if not strip.startswith("$"):
            continue
        for token in strip.split()[1:]:
            if token.endswith(_SIDE_SUFFIXES) and token not in seen:
                seen.append(token)
    return tuple(seen)


def read_mokit_run_file(path: str | Path) -> MokitRun:
    """Read a capture file; raise MokitError if it is missing."""
    source = Path(path)
    if not source.is_file():
        raise MokitError(f"{source} does not exist. Next step: give the automr output path.")
    return read_mokit_run(source.read_text(encoding="utf-8", errors="replace"))
