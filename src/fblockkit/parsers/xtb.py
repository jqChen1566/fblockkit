"""The xTB run capture (plan item 6.1; the reader behind the pre-screening menu).

xTB (Bannwarth, Caldeweyher, Ehlert, Hansen, Pracht, Seibert, Spicher,
Grimme, WIREs Comput. Mol. Sci. 2021, 11, e01493; GFN2-xTB: Bannwarth,
Ehlert, Grimme, J. Chem. Theory Comput. 2019, 15, 1652) is the Tier-1
pre-screening engine of this group's protocol chain.  This module reads the
captured stdout of a run back into plain data; every convention below is
**measured** on xTB 6.7.1 (fixtures ``fixtures/xtb/``, 2026-09-30):

- the banner line ``* xtb version 6.7.1 (edcfbbe) compiled by '...' on
  2024-07-22`` (the identification anchor);
- the result box printed by every run kind: ``TOTAL ENERGY ... Eh`` /
  ``GRADIENT NORM ... Eh/alpha`` / ``HOMO-LUMO GAP ... eV`` (the last value
  is the one of the final evaluation);
- the SCF marker ``*** convergence criteria satisfied after N iterations
  ***`` (a run may satisfy it several times; the entries are collected);
- the geometry-optimisation markers: ``*** GEOMETRY OPTIMIZATION CONVERGED
  AFTER N ITERATIONS ***`` or -- with ``--cycles N`` exhausted --
  ``*** FAILED TO CONVERGE GEOMETRY OPTIMIZATION IN N ITERATIONS ***``;
- the frequency block(s): the title line ``vibrational frequencies (cm-1)``
  (``projected vibrational frequencies`` after an optimisation) followed by
  ``eigval :`` rows, six values each, running over the line group; the same
  values are printed **twice** (once from the Hessian step, once inside the
  thermochemistry section) -- they are kept as groups, and the report uses
  the first group;
- the thermochemistry summary: the ``# frequencies`` / ``# imaginary freq.``
  setup counters, the ``found N significant imaginary frequency`` line
  (with its ``imag cut-off``), and the ``:: total free energy ... Eh ::``
  family of boxed lines;
- the closing marker ``* finished run on YYYY/MM/DD at HH:MM:SS.sss``
  (stdout); note that ``normal termination of xtb`` goes to **stderr** and
  is therefore not part of a stdout-only capture.

Mode note: with a linear starting geometry the two rotational zero modes
are reported among the low eigenvalues and no imaginary bend appears in the
harmonic picture (measured on linear H2O: five zero modes and a positive
degenerate bend); this module only reports what the file carries.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

__all__ = [
    "XtbError",
    "XtbRun",
    "read_xtb_run",
    "read_xtb_run_file",
]


class XtbError(ValueError):
    """Not an xTB run capture (or a capture without the version banner)."""


@dataclass(frozen=True)
class XtbRun:
    """The measured content of one xtb stdout capture."""

    version: str
    build: str
    compiled_by: str
    compiled_on: str
    tasks: tuple[str, ...]           # any of: single_point / optimization / frequencies
    energy_Eh: float | None
    gradient_norm: float | None
    homo_lumo_gap_eV: float | None
    scf_cycles: tuple[int, ...]      # one entry per SCF convergence marker
    opt_converged: bool | None       # None = no optimisation marker in the capture
    opt_cycles: int | None
    frequency_groups: tuple[tuple[float, ...], ...]
    engine_n_frequencies: int | None  # the thermochemistry setup counter
    engine_n_imaginary: int | None
    imag_found: int | None            # the "found N significant imaginary" line
    thermo: dict[str, float] = field(default_factory=dict)
    finished: str | None = None       # "YYYY/MM/DD HH:MM:SS.sss"

    @property
    def frequencies(self) -> tuple[float, ...]:
        """The first frequency group (the Hessian step; the later copy repeats it)."""
        return self.frequency_groups[0] if self.frequency_groups else ()

    @property
    def n_imaginary(self) -> int:
        """Imaginary modes counted from the first frequency group."""
        return sum(1 for value in self.frequencies if value < 0)


_BANNER_RE = re.compile(
    r"\* xtb version (\S+) \(([0-9a-f]+)\) compiled by '([^']*)' on ([\d-]+)"
)
_TOTAL_ENERGY_RE = re.compile(r"TOTAL ENERGY\s+(-?\d+\.\d+) Eh")
_GRADIENT_RE = re.compile(r"GRADIENT NORM\s+(-?\d+\.\d+) Eh/")
_GAP_RE = re.compile(r"HOMO-LUMO GAP\s+(-?\d+\.\d+) eV")
_SCF_RE = re.compile(r"\*\*\* convergence criteria satisfied after (\d+) iterations \*\*\*")
_OPT_CONVERGED_RE = re.compile(
    r"\*\*\* GEOMETRY OPTIMIZATION CONVERGED AFTER (\d+) ITERATIONS \*\*\*"
)
_OPT_FAILED_RE = re.compile(
    r"\*\*\* FAILED TO CONVERGE GEOMETRY OPTIMIZATION IN (\d+) ITERATIONS \*\*\*"
)
_FREQ_TITLE_RE = re.compile(r"(?:projected )?vibrational frequencies \(cm")
_EIGVAL_RE = re.compile(r"^eigval :\s*(.+?)\s*$")
_N_FREQ_RE = re.compile(r":\s*# frequencies\s+(\d+)")
_N_IMAG_RE = re.compile(r":\s*# imaginary freq\.\s+(\d+)")
_IMAG_FOUND_RE = re.compile(r"found\s+(\d+)\s+significant imaginary frequency")
_THERMO_LINE_RE = re.compile(r"::\s+(.+?)\s+(-?\d+\.\d+) Eh\s+::")
_FINISHED_RE = re.compile(r"\* finished run on ([\d/]+) at ([\d:.]+)")

_THERMO_KEYS = {
    "total free energy": "free_energy_Eh",
    "total energy": "total_energy_Eh",
    "zero point energy": "zpe_Eh",
}


def read_xtb_run(text: str) -> XtbRun:
    """Parse one xtb stdout capture; raise XtbError without the version banner."""
    match = _BANNER_RE.search(text)
    if match is None:
        raise XtbError(
            "no '* xtb version ...' banner found: this does not look like an xTB "
            "stdout capture. Next step: redirect the xtb stdout (and stderr) into a "
            "file and give that file."
        )
    lines = text.splitlines()

    # the result box repeats per geometry step; the last entry is the final one
    energies = [float(value) for value in _TOTAL_ENERGY_RE.findall(text)]
    gradients = [float(value) for value in _GRADIENT_RE.findall(text)]
    gaps = [float(value) for value in _GAP_RE.findall(text)]
    scf_cycles = tuple(int(value) for value in _SCF_RE.findall(text))

    opt_converged: bool | None = None
    opt_cycles: int | None = None
    for line in lines:
        if (found := _OPT_CONVERGED_RE.search(line)) is not None:
            opt_converged = True
            opt_cycles = int(found.group(1))
        elif (found := _OPT_FAILED_RE.search(line)) is not None:
            opt_converged = False
            opt_cycles = int(found.group(1))

    frequency_groups = _frequency_groups(lines)

    tasks: list[str] = []
    if opt_converged is not None:
        tasks.append("optimization")
    if frequency_groups:
        tasks.append("frequencies")
    if not tasks:
        tasks.append("single_point")

    engine_n_freq: int | None = None
    engine_n_imag: int | None = None
    imag_found: int | None = None
    thermo: dict[str, float] = {}
    for line in lines:
        if (found := _N_FREQ_RE.search(line)) is not None and engine_n_freq is None:
            engine_n_freq = int(found.group(1))
        if (found := _N_IMAG_RE.search(line)) is not None and engine_n_imag is None:
            engine_n_imag = int(found.group(1))
        if (found := _IMAG_FOUND_RE.search(line)) is not None:
            imag_found = int(found.group(1))
        if (found := _THERMO_LINE_RE.search(line)) is not None:
            key = _THERMO_KEYS.get(found.group(1).strip())
            if key is not None:
                thermo[key] = float(found.group(2))

    finished: str | None = None
    if (found := _FINISHED_RE.search(text)) is not None:
        finished = f"{found.group(1)} {found.group(2)}"

    return XtbRun(
        version=match.group(1),
        build=match.group(2),
        compiled_by=match.group(3),
        compiled_on=match.group(4),
        tasks=tuple(tasks),
        energy_Eh=energies[-1] if energies else None,
        gradient_norm=gradients[-1] if gradients else None,
        homo_lumo_gap_eV=gaps[-1] if gaps else None,
        scf_cycles=scf_cycles,
        opt_converged=opt_converged,
        opt_cycles=opt_cycles,
        frequency_groups=frequency_groups,
        engine_n_frequencies=engine_n_freq,
        engine_n_imaginary=engine_n_imag,
        imag_found=imag_found,
        thermo=thermo,
        finished=finished,
    )


def _frequency_groups(lines: list[str]) -> tuple[tuple[float, ...], ...]:
    """Consecutive ``eigval :`` rows form one group; a non-row line closes it.

    A frequency run prints the same values twice (the Hessian step and the
    thermochemistry section), so the capture usually carries two identical
    groups; both are kept as measured.
    """
    groups: list[tuple[float, ...]] = []
    current: list[float] = []
    for line in lines:
        if (match := _EIGVAL_RE.match(line)) is not None:
            current.extend(float(value) for value in match.group(1).split())
            continue
        if current:
            groups.append(tuple(current))
            current = []
    if current:
        groups.append(tuple(current))
    return tuple(groups)


def read_xtb_run_file(path: str | Path) -> XtbRun:
    """Read a capture file; raise XtbError if it is missing."""
    source = Path(path)
    if not source.is_file():
        raise XtbError(f"{source} does not exist. Next step: give the capture file path.")
    return read_xtb_run(source.read_text(encoding="utf-8", errors="replace"))
