"""4.2: CASSCF per-state data from the property file (menu 31).

What a single output can and cannot say about state tracking (the availability
survey of Wave 4.2, measured on ORCA 6.1.1):

- **CI vectors are not persistable.**  The CASSCF module keeps them in run-time
  temporaries (ORCA's TCIVectorStorage restart files), and the ``.cis`` file the
  JSON exporter looks for belongs to the CIS/STEOM modules -- ``orca_2json``
  prints ``Warning (ORCA_JSON): STEOM CIS file could not be found: <base>.cis``
  on a CASSCF job.  The CI vectors cannot be asked for with a print flag.
- **State characters are printed twice**: the initial ``INITIAL CI STATE
  CHECK`` and the final ``CAS-SCF STATES FOR BLOCK`` blocks carry the dominant
  CSF occupations per root (the ``.out`` side; the two snapshots are taken
  before and after the orbital optimisation).
- **Per-state energies are structured in the property file**: ``$CAS_SCF_Energies``
  carries ``totalEnergy``/``Mult``/``Irrep`` per state (and the block/root
  indices), byte-matching the final ``ROOT n:`` lines of the output (this
  module's measured cross-check).
- **Per-state observables**: the state-averaged dipole in the property file is
  one vector (x, y, z; ``State -1``), not a per-root table -- the earlier
  menu-19/20 conclusion stands, with the addition that per-root dipoles ARE
  available from single-root runs (the menu-19 batch route), and per-root
  *densities* from the FIC-NEVPT2 sidecar (the menu-23 chain).
- Consequently, identity tracking along a series needs per-state observables
  from such runs; it is registered as a candidate increment, not claimed here.

What this module does: read ``$CAS_SCF_Energies`` (the per-state table) and the
``CASSCF_Absorption_Spectrum`` section (per-transition state pairs,
multiplicities and excitation energies) out of a parsed property file, with the
measured cross-checks (the energies against the output, the eV/cm-1 conversion
against itself) and honest boundaries for the columns the manual leaves
unnamed.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..knowledge.models import EVIDENCE_MANUAL, EVIDENCE_MEASURED, Evidence
from ..parsers.orca_property import PropertySection

__all__ = [
    "StateDataError",
    "StateEntry",
    "Transition",
    "state_table",
    "transitions",
    "render",
    "evidence",
]

#: eV -> cm-1 (the 2018 CODATA value; used only as an internal consistency gate).
EV_TO_CM1 = 8065.544

_PROPERTY_FILE_URL = (
    "https://www.faccts.de/docs/orca/6.1/manual/contents/utilitiesvisualization/"
    "property_file_list.html"
)


class StateDataError(ValueError):
    """The property file lacks the sections this analysis needs (with a next step)."""


@dataclass(frozen=True)
class StateEntry:
    """One CASSCF state: its indices, multiplicity, irrep and energy."""

    index: int  # row in the property array (0-based)
    block: int | None
    root: int | None
    mult: int | None
    irrep: int | None
    energy_eh: float


@dataclass(frozen=True)
class Transition:
    """One absorption line: the state pair and its excitation energy."""

    initial: int
    initial_irrep: int
    final: int
    final_irrep: int
    mult_initial: int | None
    mult_final: int | None
    de_ev: float
    de_cm1: float
    extra: tuple[float, ...]  # the columns the manual leaves unnamed


def _section(sections: tuple[PropertySection, ...], name: str) -> PropertySection:
    for section in sections:
        if section.name == name:
            return section
    raise StateDataError(
        f"the property file carries no ${name} section. Next step: give the property "
        "file of a CASSCF run (``<base>.property.txt``, written automatically)."
    )


def _column(section: PropertySection, field: str, *, required: bool = True):
    item = section.field(field)
    if item is None:
        if required:
            raise StateDataError(
                f"the ${section.name} section carries no '{field}' field; the schema in "
                "the manual's property-file appendix lists it, so this file is of a "
                "different layout."
            )
        return None
    return item.value


def state_table(sections: tuple[PropertySection, ...]) -> dict:
    """The per-state table of ``$CAS_SCF_Energies`` (with the run's own summary)."""
    section = _section(sections, "CAS_SCF_Energies")
    energies = tuple(row[0] for row in _column(section, "totalEnergy"))
    mults = _column(section, "Mult", required=False)
    irreps = _column(section, "Irrep", required=False)
    blocks = _column(section, "Block", required=False)
    roots = _column(section, "Root", required=False)
    n_states = len(energies)
    for name, column in (
        ("Mult", mults),
        ("Irrep", irreps),
        ("Block", blocks),
        ("Root", roots),
    ):
        if column is not None and len(column) != n_states:
            raise StateDataError(
                f"the ${section.name} '{name}' array holds {len(column)} entries while "
                f"totalEnergy holds {n_states}."
            )
    states = tuple(
        StateEntry(
            index=index,
            block=blocks[index][0] if blocks is not None else None,
            root=roots[index][0] if roots is not None else None,
            mult=mults[index][0] if mults is not None else None,
            irrep=irreps[index][0] if irreps is not None else None,
            energy_eh=energy,
        )
        for index, energy in enumerate(energies)
    )
    return {
        "method": _column(section, "Method", required=False),
        "final_energy": _column(section, "finalEnergy", required=False),
        "n_active_electrons": _column(section, "numOfActiveEl", required=False),
        "n_active_orbitals": _column(section, "numOfActiveOrbs", required=False),
        "states": states,
    }


def transitions(sections: tuple[PropertySection, ...]) -> tuple[Transition, ...]:
    """The ``CASSCF_Absorption_Spectrum`` transitions (state pairs + energies).

    The schema (manual's property-file appendix) names ``States`` (initial
    state, initial irrep, final state, final irrep) and ``Multiplicities``, and
    lists ``ExcitationEnergies`` without column names; measured, the first two
    columns are the excitation energy in eV and in cm-1 (they reproduce the
    ``ROOT n:`` lines' ``eV``/``cm**-1`` values, this module's cross-check).
    The remaining columns are carried through as ``extra`` without
    interpretation.
    """
    section = _section(sections, "CASSCF_Absorption_Spectrum")
    pairs = _column(section, "States")
    multiplicities = _column(section, "Multiplicities", required=False)
    energies = _column(section, "ExcitationEnergies")
    if len(pairs) != len(energies) or (
        multiplicities is not None and len(multiplicities) != len(energies)
    ):
        raise StateDataError(
            "the CASSCF_Absorption_Spectrum arrays disagree about the transition count "
            f"({len(pairs)} state rows, {len(energies)} energy rows)."
        )
    collected = []
    for index, (pair, row) in enumerate(zip(pairs, energies)):
        if len(row) < 2:
            raise StateDataError(
                f"transition {index} carries {len(row)} value(s); the measured layout "
                "starts with the excitation energy in eV and in cm-1."
            )
        de_ev, de_cm1 = row[0], row[1]
        if de_ev > 1.0 and abs(de_ev * EV_TO_CM1 - de_cm1) > 0.5:
            raise StateDataError(
                f"transition {index} is inconsistent with itself: {de_ev:g} eV x "
                f"{EV_TO_CM1:g} = {de_ev * EV_TO_CM1:.1f} != {de_cm1:.1f} cm-1, so the "
                "first two columns are not the eV/cm-1 pair the measured layout has."
            )
        collected.append(
            Transition(
                initial=pair[0],
                initial_irrep=pair[1],
                final=pair[2],
                final_irrep=pair[3],
                mult_initial=multiplicities[index][0] if multiplicities is not None else None,
                mult_final=multiplicities[index][1] if multiplicities is not None else None,
                de_ev=de_ev,
                de_cm1=de_cm1,
                extra=tuple(row[2:]),
            )
        )
    return tuple(collected)


def render(table: dict, lines: tuple[Transition, ...], *, source: str) -> str:
    """The menu's report: the per-state table, the transitions, the boundaries."""
    out = [
        f"CASSCF state data from {source}",
        f"  method: {table['method']}; final (SA) energy: "
        f"{table['final_energy']:.9f} Eh; active space: "
        f"({table['n_active_electrons']}e, {table['n_active_orbitals']}o)",
        "",
        "  state  block  root  mult  irrep        energy (Eh)      (relative, eV)",
    ]
    ground = min(state.energy_eh for state in table["states"]) if table["states"] else 0.0
    for state in table["states"]:
        block = "-" if state.block is None else str(state.block)
        root = "-" if state.root is None else str(state.root)
        mult = "-" if state.mult is None else str(state.mult)
        irrep = "-" if state.irrep is None else str(state.irrep)
        relative = (state.energy_eh - ground) * 27.211386
        out.append(
            f"  {state.index:5d}  {block:>5}  {root:>4}  {mult:>4}  {irrep:>5}  "
            f"{state.energy_eh:>16.9f}  {relative:>14.4f}"
        )
    if lines:
        out += [
            "",
            "  absorption transitions (electric-dipole route, density: see the section):",
            "  initial -> final (irrep)   mult (i->f)      dE (eV)      dE (cm**-1)",
        ]
        for line in lines:
            out.append(
                f"  {line.initial} -> {line.final} ({line.initial_irrep}->{line.final_irrep})   "
                f"{line.mult_initial}->{line.mult_final}        "
                f"{line.de_ev:>10.4f}   {line.de_cm1:>12.1f}"
            )
        if any(value != 0.0 for value in lines[0].extra):
            unnamed = ", ".join(f"{value:.6g}" for value in lines[0].extra)
            out.append(
                "  note: this transition carries further columns the manual's schema "
                f"leaves unnamed ({unnamed}); they are printed in the report file, not "
                "interpreted."
            )
    out += [
        "",
        "Boundaries (the Wave-4.2 availability survey):",
        "  - the CI vectors are not persistable (run-time temporaries only; the .cis "
        "file belongs to the CIS/STEOM modules), so state identity along a series "
        "cannot be read from one output;",
        "  - the .out side prints the per-root dominant CSF occupations twice (initial "
        "state check and final states block), and this property file adds the "
        "structured per-state energies and transitions;",
        "  - per-state observables for identity tracking: single-root runs carry the "
        "per-state dipole (the menu-19 route), the FIC-NEVPT2 sidecar carries per-root "
        "densities (the menu-23 chain) -- a tracker across a series is registered as a "
        "candidate increment.",
    ]
    return "\n".join(out)


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the property-file schema and of the measured cross-checks."""
    return (
        Evidence(
            kind=EVIDENCE_MANUAL,
            text=(
                "The property file's section schema (CAS_SCF_Energies with "
                "totalEnergy/Mult/Irrep per state; CASSCF_Absorption_Spectrum with "
                "States/Multiplicities/ExcitationEnergies; Absorption_Spectrum field "
                "names) is the manual's own appendix list."
            ),
            ref="ORCA 6.1 manual, property-file list (appendix to the utilities chapter)",
            url=_PROPERTY_FILE_URL,
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "Measured on ORCA 6.1.1 (fixtures n2_sa.*, h2o_sa.*): the property "
                "file's per-state totalEnergy values reproduce the output's final "
                "'ROOT n: E=' lines to their printed precision (N2: -108.9756667043 / "
                "-108.5872679717 / -108.5582885454 Eh), the excitation energies' first "
                "two columns are the eV/cm-1 pair (they reproduce the output's "
                "eV/cm**-1 values and satisfy the 8065.544 conversion), and the "
                "state-averaged dipole is one x/y/z vector (State -1) rather than a "
                "per-root table."
            ),
            ref="tests/test_orca_property.py, tests/test_state_data.py; fixtures/orca/n2_sa.*",
        ),
    )
