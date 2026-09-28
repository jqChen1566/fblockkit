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
multiplicities, excitation energies and the electric-dipole data) out of a
parsed property file, with the measured cross-checks (the energies against the
output, the eV/cm-1 conversion against itself, the oscillator strength against
(2/3) dE_au D2, D2 against the dipole components' squares) and honest
boundaries.

The absorption column layout was pinned in Wave 5.1 on purpose-built ORCA
probes with **nonzero** values (``fixtures/orca/h2o_absp.*``; the earlier N2
fixture has symmetry-forbidden zeros only): 11 columns =
[eV, cm-1, nm, fosc, D2, then the complex dipole components as (re, im)
pairs DX, DY, DZ]; a SOC run prints the two copies of the section
(``&RelCorrection`` 1 and 2) and the higher one is read.  The associated
magnetic-dipole data live in ``$CASSCF_ECD_Spectrum``.  This is the data the
Judd-Ofelt menu (34) points at for per-transition intensities.
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
#: eV per Hartree (CODATA 2018; the f = (2/3) dE_au D2 gate).
EV_PER_HARTREE = 27.211386245988

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
    """One absorption line: the state pair, its energy and its dipole data.

    The measured layout of ``ExcitationEnergies`` (Wave 5.1, the H2O probe
    fixtures ``h2o_absp.*`` with nonzero values): column 0/1 = energy in eV
    and cm-1, 2 = wavelength (nm), 3 = oscillator strength fosc (D2 gauge),
    4 = D2 (au^2), then **six numbers = the complex dipole components as
    (real, imaginary) pairs**: DX(5,6), DY(7,8), DZ(9,10) -- measured by
    aligning the columns with the output's ``ABSORPTION SPECTRUM`` block
    (DX of the strong line lands on column 5, DY on 7, DZ on 9; the
    imaginary parts vanish for real wave functions).  The associated
    magnetic-dipole data live in the separate ``$CASSCF_ECD_Spectrum``
    section.  The two internal gates below (f = (2/3) dE_au D2 and
    D2 = |DX|^2+|DY|^2+|DZ|^2) reproduce the printed numbers to their
    precision.
    """

    initial: int
    initial_irrep: int
    final: int
    final_irrep: int
    mult_initial: int | None
    mult_final: int | None
    de_ev: float
    de_cm1: float
    wavelength_nm: float | None = None
    fosc: float | None = None
    d2_au: float | None = None
    dipoles: tuple[complex, complex, complex] | None = None
    extra: tuple[float, ...] = ()  # any columns beyond the measured eleven


def _rel_correction(section: PropertySection) -> int:
    item = section.field("RelCorrection")
    if item is None or not isinstance(item.value, (int, float)):
        return -1
    return int(item.value)


def _section(sections: tuple[PropertySection, ...], name: str) -> PropertySection:
    matches = [section for section in sections if section.name == name]
    if not matches:
        raise StateDataError(
            f"the property file carries no ${name} section. Next step: give the "
            "property file of a CASSCF run (``<base>.property.txt``, written "
            "automatically)."
        )
    # Measured (h2o_absp_soc fixture): a SOC run prints the absorption section
    # twice -- the pre-SOC copy (&RelCorrection 1) and the SOC-corrected one
    # (&RelCorrection 2).  The highest relativistic treatment is the one the
    # user asked for, so that copy is read (deterministic: the first of the
    # equally-ranked ones, in file order).
    return max(matches, key=_rel_correction)


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
    """The ``CASSCF_Absorption_Spectrum`` transitions (state pairs, energies, dipoles).

    The measured column layout is documented on :class:`Transition`: the first
    two columns are the excitation energy in eV and cm-1 (they reproduce the
    ``ROOT n:`` lines' ``eV``/``cm**-1`` values), then wavelength, oscillator
    strength, dipole strength, and the six numbers of the complex dipole
    components.  Two measured internal gates run per row: the oscillator
    strength against (2/3) dE_au D2, and D2 against the dipole components
    (both reproduce the printed values to their precision).  A SOC run prints
    the section twice; the copy with the highest ``&RelCorrection`` is read
    (:func:`_section`).
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
        wavelength = row[2] if len(row) > 2 else None
        fosc = row[3] if len(row) > 3 else None
        d2 = row[4] if len(row) > 4 else None
        dipoles = None
        if len(row) >= 11:
            dxr, dxi, dyr, dyi, dzr, dzi = row[5:11]
            dipoles = (complex(dxr, dxi), complex(dyr, dyi), complex(dzr, dzi))
        if wavelength is not None and de_cm1 > 0.0:
            expected = 1e7 / de_cm1
            if abs(expected - wavelength) > 1e-4 * wavelength:
                raise StateDataError(
                    f"transition {index}: wavelength {wavelength:g} nm does not match "
                    f"the energy ({expected:.4g} nm); the column is not the measured "
                    "wavelength."
                )
        if fosc is not None and d2 is not None and d2 > 0.0:
            expected = (2.0 / 3.0) * (de_ev / EV_PER_HARTREE) * d2
            if abs(expected - fosc) > 1e-3 * fosc:
                raise StateDataError(
                    f"transition {index}: fosc {fosc:.6g} does not satisfy "
                    f"f = (2/3) dE_au D2 ({expected:.6g}); the column order is not "
                    "the measured one."
                )
        if d2 is not None and dipoles is not None and d2 > 0.0:
            strength = sum(abs(component) ** 2 for component in dipoles)
            if abs(strength - d2) > 1e-3 * d2:
                raise StateDataError(
                    f"transition {index}: D2 {d2:.6g} does not equal the dipole "
                    f"components' sum |DX|^2+|DY|^2+|DZ|^2 ({strength:.6g})."
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
                wavelength_nm=wavelength,
                fosc=fosc,
                d2_au=d2,
                dipoles=dipoles,
                extra=tuple(row[11:]),
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
            "  initial -> final (irrep)   mult (i->f)      dE (eV)      dE (cm**-1)   fosc",
        ]
        for line in lines:
            fosc = "-" if line.fosc is None else f"{line.fosc:>10.4e}"
            out.append(
                f"  {line.initial} -> {line.final} ({line.initial_irrep}->{line.final_irrep})   "
                f"{line.mult_initial}->{line.mult_final}        "
                f"{line.de_ev:>10.4f}   {line.de_cm1:>12.1f}   {fosc}"
            )
        if any(value != 0.0 for value in lines[0].extra):
            unnamed = ", ".join(f"{value:.6g}" for value in lines[0].extra)
            out.append(
                "  note: this transition carries further columns beyond the measured "
                f"layout ({unnamed}); they are printed in the report file, not "
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
        "candidate increment;",
        "  - the absorption columns are the measured layout [eV, cm-1, nm, fosc, D2, "
        "DX(re,im), DY(re,im), DZ(re,im)] (nonzero probes, Wave 5.1); a SOC run's "
        "duplicate section is read at its highest &RelCorrection, and the "
        "magnetic-dipole data live in the separate ECD section (not shown here).",
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
                "Measured on ORCA 6.1.1 (fixtures n2_sa.*, h2o_sa.*, h2o_absp.*): the "
                "property file's per-state totalEnergy values reproduce the output's "
                "final 'ROOT n: E=' lines to their printed precision (N2: "
                "-108.9756667043 / -108.5872679717 / -108.5582885454 Eh), the "
                "absorption section's 11 columns are [eV, cm-1, nm, fosc, D2, DX(re,im), "
                "DY(re,im), DZ(re,im)] (pinned on the nonzero H2O probe: the "
                "oscillator strengths satisfy f = (2/3) dE_au D2 and D2 = the squared "
                "components), a SOC run's duplicate section is resolved by "
                "&RelCorrection (highest wins), and the state-averaged dipole is one "
                "x/y/z vector (State -1) rather than a per-root table."
            ),
            ref="tests/test_state_data.py, tests/test_orca_property.py; fixtures/orca/n2_sa.*, h2o_absp.*",
        ),
    )
