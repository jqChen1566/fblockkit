"""ASS1ST active-space construction from NEVPT2 quasi-natural occupations.

The source scheme (Khedkar & Roemelt, JCTC 2019 and 2020) builds an active
space bottom-up: from a small but chemically reasonable initial CASSCF, a
perturbation-theory density is made and its internal/internal and
external/external blocks are diagonalized *separately*, giving two sets of
quasi-natural orbitals; every orbital whose occupation falls inside a threshold
band is taken as strongly correlated, and the active space is grown (an active
orbital whose occupation has drifted to ~2 or ~0 is reassigned out again --
shrinking is allowed).  The construction is applied round after round until
nothing changes.

This module is the *selection round*: it reads one export of a
CASSCF + perturbation-theory run (menus 22 generates the round-1 input; menu 23
runs this analysis and writes the next round's input), prints the block
quasi-occupation tables with the band marked, and proposes the next active
space with the source's electron/orbital bookkeeping:

    add an internal quasi-orbital (occupation ~2):  +2 electrons, +1 orbital
    add an external quasi-orbital (occupation ~0):  +0 electrons, +1 orbital
    reassign an active orbital to internal:         -2 electrons, -1 orbital
    reassign an active orbital to external:         -0 electrons, -1 orbital

(checked against the source's worked example, iron porphyrin (16,13) ->
(10,10) by removing three ~doubly occupied active orbitals).

Measured interface (ORCA 6.1.1, fixtures ``n2_ass1st.json`` /
``n2_ass1st_sa.json``)
--------------------------------------------------------------------------

- ORCA carries the *unrelaxed* perturbation-theory density only for its
  FIC-NEVPT2 ansatz (measured: SC-NEVPT2 plus ``Density Unrelaxed`` is refused
  with "Unrelaxed densities are only available for FIC ansatz").  The source
  works with the SC-NEVPT2 first-order density; the substitution FIC for SC is
  carried as an honest deviation (same object family -- the density of the
  zeroth-plus-first-order wavefunction -- but a different contraction, so
  counts near the thresholds may shift; no literature comparison exists for the
  exchange).
- The density lands in the run's ``<base>.densities`` sidecar under
  ``Tdens-CASNEV.mult.<M>.root.<R>.p`` (one entry per state; the companion
  ``Tdens-CAS`` entry is its CASSCF reference density), and ``orca_2json``
  with ``"Densities": ["all"]`` exports them; the measured convention is
  ``D_AO = S D_json S`` (validated: the CASSCF reference density is recovered
  exactly diagonal in its natural-orbital basis with the printed occupations,
  couplings 1e-15; the NEVPT2 density's whole-space eigenvalues reproduce
  ORCA's printed "Natural Orbital Occupation Numbers" block to 5e-9).
- ORCA's printed NOON block is the *whole-space* naturalization, not the
  per-orbital diagonal and not the block spectra; this module therefore does
  the block diagonalization itself, which is the source's construction.  On
  the N2 fixture the internal block is already nearly diagonal, while the
  external block is not (off-diagonals as large as the occupations), so the
  block step is what makes the external side meaningful at all.
- State averaging: one density per state; this module averages the MO-basis
  blocks with the given weights (the source's eq. 23) -- equal weights by
  default, which is ORCA's own default (its output prints ``ROOT=n WEIGHT=``).
- The band is ``[T_ext, T_int]`` with ``T_int = 2 - T_ext`` for a single
  threshold ``T`` (the source's convention; T = 0.05 conservative, 0.03
  relaxed), or two independent lines (the 2020 paper uses pairs like
  0.03/1.96).  The reassignment test uses the same band against the CASSCF
  natural occupations of the current round.

The measured regression anchors are the N2/def2-SVP CAS(6,6) rounds: the
internal block spectra 2.00000, 2.00000, 1.98196, 1.97680 and the suggestion
(6,6) -> (4,4) (the pi/pi* quartet -- drop sigma-2p to internal, sigma*-2p to
external), reproduced by both the single-root and the 3-singlet state-averaged
fixture; tests/test_ass1st.py pins them.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np

from ..knowledge.models import EVIDENCE_LITERATURE, EVIDENCE_MEASURED, Evidence, ReportSection
from ..parsers.orca_json import OrcaJson

__all__ = [
    "ASS1ST_RECOMMENDED_THRESHOLD",
    "Ass1stError",
    "QuasiOrbital",
    "Partition",
    "Suggestion",
    "Ass1stRound",
    "parse_band",
    "analyze_round",
    "render",
    "run",
    "evidence",
]

#: The source's conservative single threshold (its relaxed value is 0.03).
ASS1ST_RECOMMENDED_THRESHOLD = 0.05

#: Occupation tolerance for the internal (2) / external (0) partition: the
#: export prints integer occupations exactly, fractional ones to full precision.
_OCCUPATION_TOLERANCE = 1e-6

# Measured density entry naming: Tdens-CASNEV.mult.<M>.root.<R>.p
_DENSITY_RE = re.compile(r"^Tdens-CASNEV\.mult\.(\d+)\.root\.(\d+)\.p$")


class Ass1stError(ValueError):
    """A defect that stops the round; the message carries the next step."""


# --- models -------------------------------------------------------------------


@dataclass(frozen=True)
class QuasiOrbital:
    """One quasi-natural orbital of a block: its occupation and band membership."""

    occupation: float
    in_band: bool


@dataclass(frozen=True)
class Partition:
    """The current round's orbital partition and active-space electron count."""

    internal: tuple[int, ...]
    active: tuple[int, ...]
    external: tuple[int, ...]
    n_electrons: int  # sum of the active natural occupations, rounded

    @property
    def n_orbitals(self) -> int:
        return len(self.active)


@dataclass(frozen=True)
class Suggestion:
    """The source's next-space proposal (empty additions = self-consistent)."""

    add_internal: tuple[QuasiOrbital, ...]
    add_external: tuple[QuasiOrbital, ...]
    reassign_internal: tuple[int, ...]
    reassign_external: tuple[int, ...]
    n_electrons: int
    n_orbitals: int

    @property
    def self_consistent(self) -> bool:
        return not (
            self.add_internal or self.add_external
            or self.reassign_internal or self.reassign_external
        )


@dataclass(frozen=True)
class Ass1stRound:
    """One analysed round, ready to render (and to build the next input from)."""

    base_name: str
    n_states: int
    weights: tuple[float, ...]
    ext_line: float
    int_line: float
    state_averaged: bool
    partition: Partition
    active_occupations: tuple[tuple[int, float], ...]
    internal_block: tuple[QuasiOrbital, ...]
    external_block: tuple[QuasiOrbital, ...]
    suggestion: Suggestion
    multiplicity: int
    cycle_note: str | None


# --- parsing ------------------------------------------------------------------


def parse_band(spec: str) -> tuple[float, float]:
    """Parse a threshold band: 'T' -> (T, 2-T); 'T_ext,T_int' -> two lines."""
    tokens = [token for token in spec.replace(",", " ").split() if token]
    if not tokens:
        raise Ass1stError(
            "the threshold band is empty. Next step: give one value (band [T, 2-T], "
            "the source's convention) or two independent lines 'T_ext,T_int'."
        )
    try:
        values = [float(token) for token in tokens]
    except ValueError as exc:
        raise Ass1stError(
            f"the threshold band {spec!r} is not numeric. Next step: give one value or "
            "'T_ext,T_int'."
        ) from exc
    if len(values) == 1:
        ext_line, int_line = values[0], 2.0 - values[0]
    elif len(values) == 2:
        ext_line, int_line = values
    else:
        raise Ass1stError(
            f"the threshold band {spec!r} has {len(values)} values; give one or two."
        )
    if not 0.0 < ext_line < int_line < 2.0:
        raise Ass1stError(
            f"the band ({ext_line}, {int_line}) is not a sensible occupation window: "
            "the external line must sit between 0 and the internal line, which must sit "
            "below 2."
        )
    return ext_line, int_line


# --- analysis -----------------------------------------------------------------


def _nevt2_densities(export: OrcaJson) -> tuple[tuple[int, str], ...]:
    """The NEVPT2 unrelaxed density entries, ordered by (multiplicity, root)."""
    if export.densities is None:
        raise Ass1stError(
            "the export carries no Densities block. Next step: re-export with orca_2json, "
            'using \'{"MOCoefficients": true, "1elIntegrals": ["S"], "Densities": ["all"]}\' '
            "as <base>.json.conf, from a run that wrote its density sidecar "
            "(the generated round-1 input already carries KeepDens)."
        )
    entries: list[tuple[int, int, str]] = []
    for title, _ in export.densities:
        match = _DENSITY_RE.match(title)
        if match is not None:
            entries.append((int(match.group(1)), int(match.group(2)), title))
    if not entries:
        raise Ass1stError(
            "the export carries no Tdens-CASNEV.* density (the NEVPT2 unrelaxed density). "
            "Next step: run the round with '!FIC-NEVPT2 KeepDens' and the %casscf "
            "PTSettings block 'Density Unrelaxed' + 'NatOrbs true' (menu 22 writes it), "
            "then re-export."
        )
    multiplicities = {entry[0] for entry in entries}
    if len(multiplicities) != 1:
        raise Ass1stError(
            f"the export mixes densities of multiplicities {sorted(multiplicities)}; "
            "one round analyses one spin multiplicity."
        )
    entries.sort()
    return tuple((root, title) for _, root, title in entries)


def _validate_export(export: OrcaJson) -> None:
    if export.overlap is None:
        raise Ass1stError(
            "the export carries no S-Matrix, which the density convention needs "
            "(measured: D_AO = S D_json S). Next step: re-export with "
            '\'{"MOCoefficients": true, "1elIntegrals": ["S"], "Densities": ["all"]}\'.'
        )


def analyze_round(
    export: OrcaJson,
    *,
    weights: tuple[float, ...] | None = None,
    band: tuple[float, float] = (ASS1ST_RECOMMENDED_THRESHOLD, 2.0 - ASS1ST_RECOMMENDED_THRESHOLD),
    previous_spaces: tuple[tuple[int, int], ...] = (),
) -> Ass1stRound:
    """Analyse one round: block quasi-occupations, the band, and the next space."""
    _validate_export(export)
    entries = _nevt2_densities(export)
    n_states = len(entries)
    if weights is None:
        weights = tuple(1.0 / n_states for _ in range(n_states))
    if len(weights) != n_states:
        raise Ass1stError(
            f"{len(weights)} state weights for {n_states} density entries; give one per "
            "state (or none for equal weights)."
        )
    if any(weight <= 0.0 for weight in weights):
        raise Ass1stError("state weights must be positive.")
    total = sum(weights)
    weights = tuple(weight / total for weight in weights)

    coefficients = np.asarray(export.mo_coefficients)
    overlap = np.asarray(export.overlap)
    occupancy = np.asarray(export.mo_occupations)
    record = dict(export.densities or ())
    density_mo = np.zeros((export.n_mo, export.n_mo))
    for weight, (_, title) in zip(weights, entries):
        density_ao = np.asarray(record[title])
        density_mo += weight * (coefficients @ (overlap @ density_ao @ overlap) @ coefficients.T)

    internal = tuple(i for i, occ in enumerate(occupancy) if abs(occ - 2.0) <= _OCCUPATION_TOLERANCE)
    external = tuple(i for i, occ in enumerate(occupancy) if abs(occ) <= _OCCUPATION_TOLERANCE)
    active = tuple(
        i for i, occ in enumerate(occupancy)
        if _OCCUPATION_TOLERANCE < occ < 2.0 - _OCCUPATION_TOLERANCE
    )
    if not internal:
        raise Ass1stError(
            "the export has no internal (doubly occupied) orbital; a CASSCF export "
            "always has one. Next step: check that this is the round's own export."
        )
    if not active:
        raise Ass1stError(
            "the export has no fractionally occupied orbital, so it is not a CASSCF "
            "round with an active space. Next step: feed the export of the round's "
            "CASSCF calculation."
        )
    ext_line, int_line = band
    if not 0.0 < ext_line < int_line < 2.0:
        raise Ass1stError(
            f"the band ({ext_line}, {int_line}) is not a sensible occupation window."
        )

    def block(values: np.ndarray, line: float, keep_above: bool) -> tuple[QuasiOrbital, ...]:
        spectra = sorted(np.linalg.eigvalsh(values), reverse=True)
        return tuple(
            QuasiOrbital(occupation=float(v), in_band=(v >= line) if keep_above else (v <= line))
            for v in spectra
        )

    internal_block = block(density_mo[np.ix_(internal, internal)], int_line, keep_above=False)
    external_block = block(density_mo[np.ix_(external, external)], ext_line, keep_above=True)

    reassign_internal = tuple(i for i in active if occupancy[i] >= int_line)
    reassign_external = tuple(i for i in active if occupancy[i] <= ext_line)
    add_internal = tuple(orb for orb in internal_block if orb.in_band)
    add_external = tuple(orb for orb in external_block if orb.in_band)

    n_electrons = int(round(float(occupancy[list(active)].sum())))
    n_orbitals = len(active)
    next_electrons = (
        n_electrons + 2 * len(add_internal) - 2 * len(reassign_internal)
    )
    next_orbitals = (
        n_orbitals + len(add_internal) + len(add_external)
        - len(reassign_internal) - len(reassign_external)
    )
    suggestion = Suggestion(
        add_internal=add_internal,
        add_external=add_external,
        reassign_internal=reassign_internal,
        reassign_external=reassign_external,
        n_electrons=next_electrons,
        n_orbitals=next_orbitals,
    )
    if not suggestion.self_consistent:
        if next_orbitals < 2 or next_electrons < 0 or next_electrons > 2 * next_orbitals:
            raise Ass1stError(
                f"the bookkeeping suggests the invalid space ({next_electrons}e, "
                f"{next_orbitals}o); relax the band or choose the space manually "
                "(the source leaves tie-breaking to chemical judgment)."
            )
        cycle_note = None
        if (next_electrons, next_orbitals) in previous_spaces:
            cycle_note = (
                f"the next space ({next_electrons}e, {next_orbitals}o) was already "
                "visited on this chain: the source warns that growing and shrinking can "
                "cycle between spaces -- pick one of the visited spaces by chemical "
                "judgment (its own recommendation) or adjust the band."
            )
    else:
        cycle_note = None

    return Ass1stRound(
        base_name=export.base_name,
        n_states=n_states,
        weights=weights,
        ext_line=ext_line,
        int_line=int_line,
        state_averaged=n_states > 1,
        partition=Partition(
            internal=internal, active=active, external=external, n_electrons=n_electrons
        ),
        active_occupations=tuple((i, float(occupancy[i])) for i in active),
        internal_block=internal_block,
        external_block=external_block,
        suggestion=suggestion,
        multiplicity=export.multiplicity,
        cycle_note=cycle_note,
    )


# --- output -------------------------------------------------------------------


def render(round_: Ass1stRound) -> str:
    """The two block spectra, the band, and the next-space suggestion."""
    partition, suggestion = round_.partition, round_.suggestion
    lines = [
        "ASS1ST selection round (NEVPT2 quasi-natural occupation numbers):",
        f"  system: {round_.base_name}; {partition.n_orbitals} active orbitals "
        f"({partition.n_electrons}e) against {len(partition.internal)} internal and "
        f"{len(partition.external)} external; "
        + (
            f"state-averaged over {round_.n_states} states, weights "
            + ", ".join(f"{w:.3f}" for w in round_.weights)
            if round_.state_averaged
            else "single state"
        ),
        f"  band: external line {round_.ext_line:g}, internal line {round_.int_line:g} "
        "(the source's quasi-NOON window; 2 - T for a single threshold T)",
        "",
        f"  {'internal block':<16}{'occupation':>12}  in band",
    ]
    for quasi in round_.internal_block:
        lines.append(
            f"  {'':<16}{quasi.occupation:>12.5f}  {'yes' if quasi.in_band else ''}"
        )
    lines.append(f"  {'external block':<16}{'occupation':>12}  in band")
    for quasi in round_.external_block:
        lines.append(
            f"  {'':<16}{quasi.occupation:>12.5f}  {'yes' if quasi.in_band else ''}"
        )
    lines += ["", "Active orbitals (CASSCF natural occupations of this round):"]
    for index, occupation in round_.active_occupations:
        lines.append(f"  MO {index:>3}  occupation {occupation:.5f}")
    lines += ["", "Suggestion:"]
    if suggestion.self_consistent:
        lines.append(
            f"  self-consistent: no orbital crosses the band in either block and no "
            f"active orbital drifted out; keep ({partition.n_electrons}e, "
            f"{partition.n_orbitals}o). If a previous round suggested the same space, "
            "the chain is finished."
        )
    else:
        if suggestion.add_internal:
            values = ", ".join(f"{orb.occupation:.5f}" for orb in suggestion.add_internal)
            lines.append(f"  add to active, from internal: {len(suggestion.add_internal)} "
                         f"orbital(s) ({values})")
        if suggestion.add_external:
            values = ", ".join(f"{orb.occupation:.5f}" for orb in suggestion.add_external)
            lines.append(f"  add to active, from external: {len(suggestion.add_external)} "
                         f"orbital(s) ({values})")
        if suggestion.reassign_internal:
            lines.append("  reassign active to internal (occupation at the top of the "
                         "band): " + ", ".join(str(i) for i in suggestion.reassign_internal))
        if suggestion.reassign_external:
            lines.append("  reassign active to external (occupation at the bottom of the "
                         "band): " + ", ".join(str(i) for i in suggestion.reassign_external))
        lines.append(
            f"  next space: ({suggestion.n_electrons}e, {suggestion.n_orbitals}o) "
            f"-- from ({partition.n_electrons}e, {partition.n_orbitals}o)"
        )
    if round_.cycle_note:
        lines.append(f"  cycle warning: {round_.cycle_note}")
    lines += ["", "Boundaries and checks:"]
    for item in (
        "the density is ORCA's FIC-NEVPT2 unrelaxed density; the source works with the "
        "SC-NEVPT2 first-order density -- same object family (the 0+1-order density), "
        "different contraction, so counts near the thresholds may shift (carried as an "
        "honest substitution; no literature comparison exists)",
        "the source's own caveats: the outcome depends on the initial space (different "
        "sensible initials converge together only with conservative thresholds), and "
        "growing plus shrinking can cycle between two spaces -- break ties by chemical "
        "judgment",
        "every block table is this tool's own block diagonalization of the exported "
        "density (the source's construction); ORCA's printed NOON list is the "
        "whole-space naturalization and is not used",
        f"multiplicity {round_.multiplicity}; the next round's input is generated "
        "next to this report by the menu",
    ):
        lines.append(f"  - {item}")
    return "\n".join(lines)


def run(export: OrcaJson, **kwargs) -> ReportSection:
    """The analyser entry point for menu 23."""
    return ReportSection(
        title="3.3 ASS1ST selection round",
        body=render(analyze_round(export, **kwargs)),
    )


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the scheme, its state-averaged extension and the interface."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "ASS1ST: bottom-up active-space construction from the quasi-natural "
                "occupation numbers of an SC-NEVPT2 first-order density -- separate "
                "diagonalization of the internal/internal and external/external density "
                "blocks, a threshold band [T, 2-T] (0.05 conservative, 0.03 relaxed, "
                "separate lines allowed), growth plus reassignment (shrinking), iteration "
                "to self-consistency, and the recorded warnings (dependence on the "
                "initial space; possible cycling between two spaces)."
            ),
            ref=(
                "J. Chem. Theory Comput. 2019, 15, 3522-3536, sections 2.1-2.5 (the "
                "procedure and Figure 3); the electron/orbital bookkeeping matches the "
                "iron-porphyrin worked example"
            ),
            url="https://doi.org/10.1021/acs.jctc.8b01293",
            bibkey="khedkar2019ass1st",
        ),
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "SA-ASS1ST: per-state density blocks averaged with the state weights "
                "(eq. 23) for excited-state targets; the 2020 extension also removes the "
                "4-RDM terms by replacing the active-space ionisation energies and "
                "electron affinities with an averaged active-orbital energy and a "
                "0.25 Eh level shift."
            ),
            ref=(
                "J. Chem. Theory Comput. 2020, 16, 4993-5005, sections 2.4-2.5 "
                "(state averaging and the IPEA-type approximation)"
            ),
            url="https://doi.org/10.1021/acs.jctc.0c00332",
            bibkey="khedkar2020sa",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "ORCA interface, measured on 6.1.1: the unrelaxed density exists only for "
                "FIC-NEVPT2 (SC is refused: 'Unrelaxed densities are only available for "
                "FIC ansatz'); the density appears in the sidecar as "
                "Tdens-CASNEV.mult.<M>.root.<R>.p per state and exports with "
                "Densities:['all'] under the convention D_AO = S D_json S (validated "
                "against the CASSCF reference density: exactly diagonal in its "
                "natural-orbital basis, couplings 1e-15; the NEVPT2 density's whole-space "
                "eigenvalues reproduce the printed NOON block to 5e-9). Regression "
                "anchors: N2/def2-SVP CAS(6,6) internal block 2.00000, 2.00000, 1.98196, "
                "1.97680 and the (6,6) -> (4,4) suggestion, from both the single-root and "
                "the three-singlet state-averaged fixture."
            ),
            ref="tests/test_ass1st.py; fixtures/orca/n2_ass1st.json, n2_ass1st_sa.json",
        ),
    )
