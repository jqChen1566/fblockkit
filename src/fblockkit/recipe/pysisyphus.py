"""generating a pysisyphus input (menu 32).

The B-layer relation to pysisyphus is "generate its input / read its output":
the program (GPL-3.0) stays external and orchestrates the engine runs itself,
while this toolkit writes the YAML it consumes and reads the artifacts it
leaves behind (menu 33).  The YAML schema below is the one pysisyphus's own
example set uses and the one the measured runs of 2026-09-28 accept
(pysisyphus 1.0.0, ORCA 6.1.1 as the calculator):

- ``geom``: ``type`` (redund / cart / tric) and ``fn`` (an XYZ file);
- ``calc``: ``type`` (``orca5``; ``orca`` is the same interface under its
  legacy name), ``keywords`` (the ``!`` line of every ORCA call; the engine
  adds ``engrad`` itself for gradients), ``charge``, ``mult``, ``pal``,
  ``mem`` (MB per core) and optionally ``blocks`` (a quoted block string);
- ``opt``: ``type`` (``rfo``), ``thresh`` and optionally ``max_cycles``;
- ``tsopt``: ``type`` (``rsprfo``), ``thresh``, ``hessian_init`` and
  optionally ``rx_modes`` (a YAML flow sequence selecting the mode to follow,
  e.g. ``[[[[DIHEDRAL, 2, 0, 1, 3], 1]]]``).

Two measured boundaries shape the generator:

- **The quantum-Hessian route is refused.**  With ``hessian_init: calc`` (the
  ``tsopt`` default) or ``do_hess: true`` the run crashes on ORCA 6's
  ``.hess`` file: its ``$multiplicity`` block is not in the grammar
  pysisyphus ships (measured identical in release 1.0.0 and on master), so
  the post-processing raises.  The generator therefore writes a **model
  Hessian** (``hessian_init: fischer`` by default; ``lindh``/``simple``/
  ``swart``/``unit`` are the other model choices) and refuses ``calc`` with
  the reason; the frequency verification happens outside pysisyphus (a plain
  ORCA ``Freq`` run -- menu 30 is the toolkit's own route back).
- **The threshold names are pysisyphus's own** (``gau_loose``, ``gau``,
  ``gau_tight``, ``gau_vtight``, ``baker``; the dummy ``never`` is refused
  because it sets a billion cycles and disables dumping).

The job kinds offered are the two this toolkit's rescue chains use: a minimum
optimisation (``opt``) and a transition-state search (``tsopt``).  The file
written here is plain YAML; nothing in it is engine-specific beyond the
``keywords`` string, which is passed through verbatim.
"""

from __future__ import annotations

import re
from typing import Sequence

from ..knowledge.models import EVIDENCE_LITERATURE, EVIDENCE_MEASURED, Evidence

__all__ = [
    "PysisyphusError",
    "PYSISYPHUS_THRESHOLDS",
    "MODEL_HESSIANS",
    "atoms_from_orca_input",
    "structure_xyz",
    "input_yaml",
    "render",
    "run_guidance_lines",
    "evidence",
]

#: the threshold set of pysisyphus's own optimizers (its Thresh literal, minus
#: the dummy "never", which sets max_cycles to a billion and disables dumping)
PYSISYPHUS_THRESHOLDS = ("gau_loose", "gau", "gau_tight", "gau_vtight", "baker")

#: the model (non-quantum) Hessians of pysisyphus's guess_hessians module
MODEL_HESSIANS = ("fischer", "lindh", "simple", "swart", "unit")

_COORD_TYPES = ("redund", "cart", "tric")
_JOB_KEYS = {"min": "opt", "ts": "tsopt"}
_DEFAULT_OPTIMIZER = {"min": "rfo", "ts": "rsprfo"}
_DEFAULT_THRESHOLD = {"min": "gau", "ts": "baker"}


class PysisyphusError(ValueError):
    """The pysisyphus input cannot be built from the given values (with a next step)."""


def _check_ascii(text: str, label: str) -> str:
    if not text:
        raise PysisyphusError(f"the {label} is empty. Next step: give a value.")
    try:
        text.encode("ascii")
    except UnicodeEncodeError as exc:
        raise PysisyphusError(
            f"the {label} is not pure ASCII ({exc}); pysisyphus inputs are."
        ) from exc
    if "\n" in text or "\r" in text:
        raise PysisyphusError(f"the {label} carries a line break; it must be one line.")
    return text


def _yaml_scalar(text: str) -> str:
    """A plain scalar when that is safe, a double-quoted one otherwise."""
    if re.fullmatch(r"[A-Za-z0-9_.\-+/*() ]+", text) and text.strip() == text:
        return text
    return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'


_XYZ_HEADER_RE = re.compile(r"^\s*\*\s*xyz\s+[-\d]", re.IGNORECASE)


def atoms_from_orca_input(text: str) -> tuple[tuple[str, float, float, float], ...]:
    """The atoms of an ORCA input's inline ``* xyz`` block.

    The same measured convention as the menu-30 reader: the header line is
    ``* xyz <charge> <mult>`` and the block ends at the next line starting
    with ``*``.
    """
    lines = text.splitlines()
    start = next((i for i, line in enumerate(lines) if _XYZ_HEADER_RE.match(line)), None)
    if start is None:
        raise PysisyphusError(
            "the file carries no inline '* xyz' block, so no structure can be taken "
            "from it. Next step: give an XYZ file, an ORCA input with inline "
            "coordinates, or an ORCA output."
        )
    end = next(
        (i for i in range(start + 1, len(lines)) if lines[i].strip().startswith("*")),
        None,
    )
    if end is None:
        raise PysisyphusError("the '* xyz' block has no closing '*'.")
    atoms: list[tuple[str, float, float, float]] = []
    for line in lines[start + 1 : end]:
        fields = line.split()
        if not fields:
            continue
        if len(fields) < 4:
            raise PysisyphusError(f"invalid coordinate line in the '* xyz' block: {line!r}")
        atoms.append((fields[0], float(fields[1]), float(fields[2]), float(fields[3])))
    if not atoms:
        raise PysisyphusError("the '* xyz' block holds no atoms.")
    return tuple(atoms)


def structure_xyz(atoms: Sequence[tuple[str, float, float, float]], *, comment: str = "") -> str:
    """A standard XYZ block (the format pysisyphus reads for ``geom.fn``)."""
    if not atoms:
        raise PysisyphusError(
            "the structure carries no atoms. Next step: give a structure file with "
            "coordinates."
        )
    lines = [str(len(atoms)), comment]
    for symbol, x, y, z in atoms:
        _check_ascii(symbol, "element symbol")
        lines.append(f"{symbol:<3s} {x:>15.10f} {y:>15.10f} {z:>15.10f}")
    return "\n".join(lines) + "\n"


def input_yaml(
    *,
    xyz_fn: str,
    keywords: str,
    charge: int = 0,
    mult: int = 1,
    job: str = "min",
    calc_type: str = "orca5",
    pal: int = 1,
    mem: int = 1500,
    thresh: str | None = None,
    coord_type: str = "redund",
    max_cycles: int | None = None,
    hessian_init: str | None = None,
    rx_modes: str | None = None,
    blocks: str | None = None,
) -> str:
    """The YAML input of one pysisyphus job (minimum optimisation or TS search).

    ``job='min'`` writes an ``opt`` block (rfo); ``job='ts'`` writes a
    ``tsopt`` block (rsprfo) with a model Hessian --- ``hessian_init: calc``
    is refused because the ORCA-6 Hessian parse stop makes it unusable in
    this pairing (measured; see the module docstring), and ``do_hess`` is not
    offered for the same reason.
    """
    if job not in _JOB_KEYS:
        raise PysisyphusError(
            f"unknown job {job!r}; this generator writes 'min' (minimum optimisation) "
            "or 'ts' (transition-state search)."
        )
    keywords = _check_ascii(keywords, "keywords")
    if not isinstance(charge, int):
        raise PysisyphusError("the charge must be an integer.")
    if not isinstance(mult, int) or mult < 1:
        raise PysisyphusError("the multiplicity must be a positive integer.")
    if not isinstance(pal, int) or pal < 1:
        raise PysisyphusError("pal (cores per calculation) must be a positive integer.")
    if not isinstance(mem, int) or mem < 100:
        raise PysisyphusError(
            "mem (MB per core) must be an integer of at least 100; next step: give "
            "the per-core memory ORCA should use."
        )
    if coord_type not in _COORD_TYPES:
        raise PysisyphusError(
            f"unknown coordinate system {coord_type!r}; use one of "
            f"{', '.join(_COORD_TYPES)}."
        )
    threshold = thresh or _DEFAULT_THRESHOLD[job]
    if threshold not in PYSISYPHUS_THRESHOLDS:
        raise PysisyphusError(
            f"unknown threshold {threshold!r}; pysisyphus's optimizer thresholds are "
            f"{', '.join(PYSISYPHUS_THRESHOLDS)} (the dummy 'never' is not offered: it "
            "sets a billion cycles and disables the dump files)."
        )
    hessian = hessian_init or ("fischer" if job == "ts" else None)
    if job == "ts":
        if hessian == "calc":
            raise PysisyphusError(
                "hessian_init 'calc' is refused: the quantum Hessian route crashes "
                "with ORCA 6 (its .hess '$multiplicity' block is not in this "
                "interface's grammar -- measured on release 1.0.0 and master), and "
                "'calc' is the tsopt default, so it must be named explicitly. Next "
                "step: keep a model Hessian (fischer/lindh/simple/swart/unit) and "
                "verify the frequencies with a plain ORCA Freq run."
            )
        if hessian not in MODEL_HESSIANS:
            raise PysisyphusError(
                f"unknown hessian_init {hessian!r}; the model Hessians are "
                f"{', '.join(MODEL_HESSIANS)} (a file path is possible in pysisyphus "
                "but is not written by this menu)."
            )
        if rx_modes is not None:
            _check_ascii(rx_modes, "rx_modes")
            if not (rx_modes.startswith("[") and rx_modes.endswith("]")):
                raise PysisyphusError(
                    "rx_modes must be a YAML flow sequence written literally, e.g. "
                    "[[[[DIHEDRAL, 2, 0, 1, 3], 1]]] (it selects the mode to follow)."
                )
    elif rx_modes is not None:
        raise PysisyphusError("rx_modes belongs to a TS search (job='ts').")
    if blocks is not None:
        _check_ascii(blocks, "blocks")
    if max_cycles is not None and (not isinstance(max_cycles, int) or max_cycles < 1):
        raise PysisyphusError("max_cycles must be a positive integer.")

    lines = [
        "# pysisyphus input, generated by fBlockKit (menu 32).",
        "# Run it with:  pysis <this file>   (keep the console:  pysis <this file> > run.out)",
        "geom:",
        f" type: {coord_type}",
        f" fn: {xyz_fn}",
        "calc:",
        f" type: {calc_type}",
        f" keywords: {_yaml_scalar(keywords)}",
        f" charge: {charge}",
        f" mult: {mult}",
        f" pal: {pal}",
        f" mem: {mem}",
    ]
    if blocks is not None:
        lines.append(f' blocks: "{blocks}"')
    lines += [
        f"{_JOB_KEYS[job]}:",
        f" type: {_DEFAULT_OPTIMIZER[job]}",
        f" thresh: {threshold}",
    ]
    if max_cycles is not None:
        lines.append(f" max_cycles: {max_cycles}")
    if job == "ts":
        lines.append(
            f" # model Hessian (the quantum-Hessian route crashes on ORCA 6's .hess; "
            f"measured)"
        )
        lines.append(f" hessian_init: {hessian}")
        if rx_modes is not None:
            lines.append(f" rx_modes: {rx_modes}")
        lines.append(" do_hess: false  # frequency verification stays outside pysisyphus")
    text = "\n".join(lines) + "\n"
    return text


def render(*, xyz_name: str, yaml_name: str, n_atoms: int, job: str, keywords: str,
           charge: int, mult: int, pal: int, mem: int, thresh: str | None,
           hessian_init: str | None) -> str:
    """The menu's record: what was written and what the two files carry."""
    kind = (
        "minimum optimisation (opt/rfo)"
        if job == "min"
        else "transition-state search (tsopt/rsprfo, model Hessian)"
    )
    lines = [
        f"pysisyphus input for {n_atoms} atom(s)",
        f"  job: {kind}",
        f"  calculator: orca5, keywords '{keywords}', charge {charge}, multiplicity {mult}, "
        f"pal {pal}, mem {mem} MB per core",
        f"  threshold: {thresh or _DEFAULT_THRESHOLD[job]}",
    ]
    if job == "ts":
        lines.append(f"  hessian: model ({hessian_init or 'fischer'})")
    lines.append(f"  written: {xyz_name} (structure, XYZ) and {yaml_name} (the input)")
    return "\n".join(lines)


def run_guidance_lines() -> tuple[str, ...]:
    """The run-side guidance (environment, scratch, the reading-back chain)."""
    return (
        "run it on a machine with pysisyphus installed and the engine configured: "
        "the program reads the engine command from ~/.pysisyphusrc (an INI file: "
        "[orca5] cmd=/path/to/orca), or from $PATH when the section is absent",
        "give pysisyphus a real scratch disk before starting (it stages every "
        "calculator call in a temporary directory under $TMPDIR and copies the "
        "results back into qm_calcs/ at the end of each call)",
        "keep the console:  pysis <file> > run.out  -- that capture is the run's "
        "record (the cycle table, the outcome marker, the closing summary), and it "
        "is what menu 33 reads back",
        "frequency verification stays outside pysisyphus in this pairing: the "
        "quantum-Hessian route (hessian_init: calc / do_hess) crashes on ORCA 6's "
        ".hess format (measured), so run a plain ORCA Freq job on the closing "
        "geometry -- menu 30 is the toolkit's own route back from an imaginary mode",
        "bring the whole run directory back to menu 33: it cross-checks the console, "
        "the trajectory, the closing geometry and the per-cycle ORCA outputs against "
        "each other",
    )


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the YAML schema and of the two measured boundaries."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "pysisyphus reads its jobs from YAML files (geom/calc/opt/tsopt "
                "sections) and drives external calculators; the citation is the "
                "program's own, printed in every run's banner."
            ),
            ref=(
                "Steinmetzer J., Kupfer S., Graefe S., Int. J. Quantum Chem. 2021, "
                "121(3), e26390, DOI 10.1002/qua.26390 (software: GPL-3.0, "
                "github.com/eljost/pysisyphus, release 1.0.0 of 2025-08-31)"
            ),
            url="https://doi.org/10.1002/qua.26390",
            bibkey="steinmetzer2021pysisyphus",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "Measured on pysisyphus 1.0.0 + ORCA 6.1.1 (fixtures/pysisyphus/, "
                "2026-09-28): the generated schema runs end to end for a minimum "
                "optimisation (H2O/def2-SVP) and for a TS search with a model "
                "Hessian (Baker case-11 butadiene at HF/3-21G, converging to "
                "pysisyphus's own asserted energy to 7e-8), while hessian_init: calc "
                "and do_hess crash in the .hess parse (the '$multiplicity' block "
                "ORCA 6 writes is not in the interface's grammar; identical in "
                "release 1.0.0 and on master)."
            ),
            ref="tests/test_pysisyphus.py; fixtures/pysisyphus/",
        ),
    )
