"""the ASS1ST round inputs (menus 22 and 23).

One ASS1ST round is a CASSCF calculation whose perturbation-theory density is
kept for the next selection step.  Both menus write the same input shape,
measured to run on ORCA 6.1.1:

    ! <keywords> FIC-NEVPT2 KeepDens
    %maxcore <MB>
    %casscf
      nel <n>
      norb <m>
      mult <M>
      nroots <R>
      PTSettings
        Density Unrelaxed
        NatOrbs true
      end
    end
    * xyz <charge> <M>
    ...

The pieces are each load-bearing (measured): ``FIC-NEVPT2`` because ORCA
answers the unrelaxed density only for the FIC ansatz (SC is refused); the
``PTSettings`` block because ``Density Unrelaxed`` is what makes the density
the perturbation-theory one and ``NatOrbs true`` keeps the natural-orbital
basis; ``KeepDens`` because the density is read from the run's
``<base>.densities`` sidecar.  Menu 22 names the first round ``<stem>.r1``;
menu 23 strips a trailing ``.rN`` and writes ``.r{N+1}`` so a chain keeps the
structure's name.  Every text is ASCII (the renderer asserts it, as ORCA
inputs must be).

The companion export request (written next to each input as
``<base>.json.conf``) is the one the selection reads:
``{"MOCoefficients": true, "1elIntegrals": ["S"], "Densities": ["all"]}``.
"""

from __future__ import annotations

import re
from pathlib import Path

from ..parsers.mkl import MklError, MklFile

__all__ = [
    "Ass1stInputError",
    "EXPORT_REQUEST",
    "export_conf",
    "round_one_input",
    "next_round_input",
    "next_round_stem",
    "round_one_stem",
    "write_qno_mkl",
]

#: The ``<base>.json.conf`` body the selection round needs (measured).
EXPORT_REQUEST = '{ "MOCoefficients": true, "1elIntegrals": ["S"], "Densities": ["all"] }\n'

_ROUND_SUFFIX = re.compile(r"\.r(\d+)$")


class Ass1stInputError(ValueError):
    """Invalid round parameters; the message carries the next step."""


def _validate(n_electrons: int, n_orbitals: int, multiplicity: int, n_states: int) -> None:
    if n_orbitals < 1:
        raise Ass1stInputError(f"an active space of {n_orbitals} orbitals does not exist.")
    if not 0 <= n_electrons <= 2 * n_orbitals:
        raise Ass1stInputError(
            f"{n_electrons} electrons do not fit in {n_orbitals} active orbitals."
        )
    if multiplicity < 1:
        raise Ass1stInputError(f"multiplicity {multiplicity} is not positive.")
    if n_states < 1:
        raise Ass1stInputError(f"{n_states} states is not a positive count.")


def _render(
    *,
    keywords: str,
    maxcore: int,
    charge: int,
    multiplicity: int,
    n_electrons: int,
    n_orbitals: int,
    n_states: int,
    coordinates: tuple[tuple[str, float, float, float], ...],
) -> str:
    lines = [
        f"! {keywords} FIC-NEVPT2 KeepDens",
        f"%maxcore {maxcore}",
        "%casscf",
        f"  nel {n_electrons}",
        f"  norb {n_orbitals}",
        f"  mult {multiplicity}",
        f"  nroots {n_states}",
        "  PTSettings",
        "    Density Unrelaxed",
        "    NatOrbs true",
        "  end",
        "end",
        f"* xyz {charge} {multiplicity}",
    ]
    for element, x, y, z in coordinates:
        lines.append(f"{element:<2} {x:>16.10f} {y:>16.10f} {z:>16.10f}")
    lines.append("*")
    text = "\n".join(lines) + "\n"
    try:
        text.encode("ascii")
    except UnicodeEncodeError as exc:  # pragma: no cover - guarded at the callers
        raise Ass1stInputError(
            f"the generated input is not ASCII ({exc}); ORCA inputs must be pure ASCII."
        ) from exc
    return text


def round_one_input(
    coordinates: tuple[tuple[str, float, float, float], ...],
    *,
    charge: int,
    multiplicity: int,
    n_electrons: int,
    n_orbitals: int,
    n_states: int = 1,
    keywords: str = "RHF def2-SVP TightSCF",
    maxcore: int = 2000,
) -> str:
    """The round-1 input: a small chemically reasonable starting space."""
    _validate(n_electrons, n_orbitals, multiplicity, n_states)
    return _render(
        keywords=keywords,
        maxcore=maxcore,
        charge=charge,
        multiplicity=multiplicity,
        n_electrons=n_electrons,
        n_orbitals=n_orbitals,
        n_states=n_states,
        coordinates=coordinates,
    )


def next_round_input(
    coordinates: tuple[tuple[str, float, float, float], ...],
    *,
    charge: int,
    multiplicity: int,
    n_electrons: int,
    n_orbitals: int,
    n_states: int = 1,
    keywords: str = "RHF def2-SVP TightSCF",
    maxcore: int = 2000,
) -> str:
    """The next round's input: the same shape at the suggested space."""
    return round_one_input(
        coordinates,
        charge=charge,
        multiplicity=multiplicity,
        n_electrons=n_electrons,
        n_orbitals=n_orbitals,
        n_states=n_states,
        keywords=keywords,
        maxcore=maxcore,
    )


def round_one_stem(structure_stem: str) -> str:
    """The round-1 base name for a structure file stem (``n2.xyz`` -> ``n2.r1``)."""
    return f"{_ROUND_SUFFIX.sub('', structure_stem)}.r1"


def next_round_stem(export_base_name: str) -> str:
    """The next round's base name: strip a trailing ``.rN`` and increment.

    ``n2.r1`` -> ``n2.r2``; a base without a round suffix starts at ``.r2``
    (its own round is taken as round 1).
    """
    match = _ROUND_SUFFIX.search(export_base_name)
    if match is None:
        return f"{export_base_name}.r2"
    return f"{export_base_name[:match.start()]}.r{int(match.group(1)) + 1}"


def export_conf() -> str:
    """The ``<base>.json.conf`` body for the round's export request."""
    return EXPORT_REQUEST


def write_qno_mkl(
    template_mkl: MklFile,
    coefficients,
    occupations,
    path,
) -> None:
    """Write a quasi-natural orbital set into a copy of the template mkl.

    ``coefficients`` (AO x MO) and ``occupations`` come from the analysis layer's
    quasi-natural export (``analysis.ass1st.quasi_natural_export``), already
    ordered so the engine's by-orbital-order window takes the inactive prefix
    and then the active window.  The template supplies the geometry and basis
    metadata; energies are left at the template's values because the
    quasi-natural orbitals carry no energies of their own (the engine
    recomputes them on read).
    """
    try:
        updated = template_mkl.with_orbitals(
            coefficients,
            occupations=[float(value) for value in occupations],
        )
    except MklError as exc:
        raise Ass1stInputError(
            f"the template mkl does not accept the quasi-natural coefficients: {exc}"
        ) from exc
    Path(path).write_text(updated.render(), encoding="utf-8")
