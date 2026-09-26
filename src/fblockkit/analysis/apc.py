"""3.1: ranked-orbital active-space selection with approximate pair coefficients.

The ranked-orbital framework (King & Gagliardi, JCTC 2021) separates *how the
orbitals are made* from *how they are ranked by importance*: every candidate
orbital (all doubly occupied ones plus a window of virtuals) gets an importance
number, and the active space is the top of that ranking, trimmed by a size cap.
The cap is a maximum number of configuration state functions,

    N_CSF = C(L, alpha) C(L, beta) - C(L, alpha+1) C(L, beta-1)      (eq. 2)

with L orbitals and (alpha, beta) the spin counts at S = S_z = 0 (even electron
count) or 1/2 (odd).  The source lists the familiar reference spaces
max(7,6) -> 490, max(8,8) -> 1764, max(10,10) -> 19404, max(12,12) -> 226512
CSFs (its sentence writes max(N_elec, N_orbs); (7e, 6o) evaluates to 210 under
eq. 2 while 490 is a 7-orbital, 6- or 8-electron space -- the notation is
internally ambiguous, the operative constraint is N_CSF <= cap, and the preset
names are carried verbatim).  The selection drops the lowest-ranked orbital
repeatedly until N_CSF <= cap, skipping a drop that would leave fewer than one
occupied or fewer than two unoccupied orbitals in the space (the source's
solver-stability floor).  Entropies that agree to 1e-12 count as equal and rank
by ORCA's orbital index: physically degenerate partners differ only in the last
floating-point digits, which can depend on the BLAS the K-matrix transform ran
under, and without the rule a report could differ between machines.

The ranking implemented here is the source's *approximate pair coefficient*
(APC).  Each doubly occupied orbital i is paired with each window virtual a
through the two-configuration model |20> + c |02>, whose exact pair coefficient
is c = -(12|12)/(Delta + sqrt((12|12)^2 + Delta^2)); APC replaces the exchange
integral by the exchange-matrix diagonal and the model gap by the orbital-energy
difference (or the Fock diagonal, for non-canonical orbitals),

    c_ia = -x_a / (D_ia + sqrt(x_a^2 + D_ia^2)),   x_a = 0.5 K_aa,  D_ia = eps_a - eps_i
                                                                     (eqs. 18, 19)

and the one-orbital entropies follow from the binary entropies of the gathered
coefficients,

    S^i = H( sum_a c_ia^2 / (1 + sum_a c_ia^2) ),   S^a = H( sum_i c_ia^2 / (1 + sum_i c_ia^2) )   (eqs. 12, 13)

with H the binary entropy.  Variant APCX keeps the exact exchange integrals
(ia|ia) instead of the 0.5 K_aa approximation -- cheaper for us than for the
source, because the same ``orca_2json`` export that carries the window can carry
both -- and is offered for comparison; the source finds APC the better ranking
for canonical HF orbitals (a fortunate cancellation of error: the diagonal-sum
approximation overestimates, the pairwise-only model underestimates).

Measured interface (ORCA 6.1.1, fixtures ``h2_apc.json`` / ``n2_apc.json``)
---------------------------------------------------------------------------

- The source's ``0.5 K_aa`` equals ``-diag(C K C^T)[a]`` for the exported
  ``K-Matrix`` block: the two agree as ``sum_i (a i | a i)`` over the doubly
  occupied orbitals, verified against the ``MO_IAJB`` entries to 4e-12 on N2
  and to 7e-16 on the H2 two-configuration model (where eq. 18 is exact).
- ``D_ia`` defaults to the exported orbital energies; ``delta="fock"`` uses the
  diagonal of the full Fock in the export's own orbital basis,
  ``F = H + J + K`` (measured: the exported ``F-Matrix`` block equals J + K and
  excludes the core Hamiltonian; ``diag(C (H + J + K) C^T)`` reproduces the
  orbital energies of the canonical N2 export to 6e-10, so the two sources are
  interchangeable there and only the Fock diagonal is meaningful for localized
  sets).
- The candidate window is the ``window_size`` virtuals lowest in energy
  (default 23, the source's general-scheme choice).  For APCX the export's
  ``MO_IAJB`` window must cover occupied x candidate virtuals; the refusal names
  the exact eight-integer ``OrbWin`` line to use (measured: the four-integer
  form is rejected by ``orca_2json``).

Scope and boundaries
--------------------

- Closed-shell RHF exports only (occupations 0/2, single spin block).  The
  source's rule for singly occupied orbitals (assign the maximum approximated
  entropy) and the UNO(HS) variant need multi-spin exports and are not covered.
- The source's own deficiencies travel with the output: APC systematically
  overestimates the doubly-occupied orbital entropies (R^2 = 0.64, MAE 0.0240
  against DMRG; virtuals 0.83/0.0064; ranking precision about 88%), and it is
  expected to perform worse in much larger systems and where the HF determinant
  is a poor approximation.  Both are printed with every report.

The regression anchors are the H2 two-configuration model (the exact pair
coefficient from the FCIDUMP 2x2 CI equals the eq. 14 closed form to 2e-17 and
the engine energy to 12 digits; eq. 18 holds to 7e-16) and the N2/def2-SVP
export, whose max(10,10) selection lands exactly on the source's (10,10) =
19404 entry in tests/test_apc.py.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from ..knowledge.models import EVIDENCE_LITERATURE, EVIDENCE_MEASURED, Evidence, ReportSection
from ..parsers.orca_json import OrcaJson

__all__ = [
    "ApcError",
    "CSF_CAPS",
    "DEFAULT_WINDOW",
    "Candidate",
    "Ranking",
    "Selection",
    "ApcReport",
    "n_csf",
    "pair_coefficients",
    "rank_orbitals",
    "select",
    "analyze",
    "render",
    "run",
    "evidence",
]

#: The source's general-scheme candidate window: all doubly occupied orbitals
#: plus the lowest 23 virtuals in energy.
DEFAULT_WINDOW = 23

#: The source's named reference caps, carried verbatim (the label sentence in
#: the source is ambiguous between (N_elec, N_orbs) and (N_orbs, N_elec); the
#: numbers are what its own figures use, and the operative constraint is the
#: CSF count).
CSF_CAPS: dict[str, int] = {
    "max(7,6)": 490,
    "max(8,8)": 1764,
    "max(10,10)": 19404,
    "max(12,12)": 226512,
}

#: Occupations are read from the export (which prints them to four decimals);
#: anything farther than this from 0 or 2 is not a closed-shell RHF orbital.
OCCUPATION_TOLERANCE = 1e-3

#: Ranking resolution: entropies agreeing to this many decimals count as equal and
#: order by index, so the rank order and the drop order do not depend on the BLAS.
_RANK_RESOLUTION = 12

_ROLE_OCCUPIED = "doubly occupied"
_ROLE_VIRTUAL = "virtual"

_VARIANT_APC = "APC"
_VARIANT_APCX = "APCX"

_DELTA_ENERGIES = "orbital energies"
_DELTA_FOCK = "Fock diagonal"

# The export requests that make each variant runnable; these exact lines are
# what the refusals hand back as the next step.
_APC_CONF = '{"MOCoefficients": true, "FockMatrix": ["K"]}'
_FOCK_CONF = '{"MOCoefficients": true, "1elIntegrals": ["H"], "FockMatrix": ["J", "K"]}'


class ApcError(ValueError):
    """A defect that stops the ranking; the message carries the next step."""


# --- models ------------------------------------------------------------------


@dataclass(frozen=True)
class Candidate:
    """One ranked orbital: its ORCA index, its role, its approximated entropy."""

    index: int
    role: str
    entropy: float


@dataclass(frozen=True)
class Ranking:
    """The ranked candidate orbitals of one export."""

    base_name: str
    variant: str
    delta_source: str
    occupied: tuple[int, ...]
    virtuals: tuple[int, ...]
    candidates: tuple[Candidate, ...]
    window: tuple[int, int]
    n_virtual_available: int


@dataclass(frozen=True)
class Selection:
    """The ranked-orbital selection at one CSF cap."""

    cap: int
    cap_label: str
    n_electrons: int
    n_orbitals: int
    n_csfs: int
    cap_reached: bool
    members: tuple[int, ...]
    dropped: tuple[int, ...]
    floor_skips: tuple[int, ...]


@dataclass(frozen=True)
class ApcReport:
    """One ranking plus its selection, ready to render."""

    ranking: Ranking
    selection: Selection
    occupied_energy: tuple[tuple[int, float], ...]
    window_energies: tuple[tuple[int, float], ...]


# --- equations ---------------------------------------------------------------


def n_csf(n_electrons: int, n_orbitals: int) -> int:
    """Equation (2): the CSF count at S = S_z = 0 (even) or 1/2 (odd)."""
    if n_orbitals < 0 or n_electrons < 0 or n_electrons > 2 * n_orbitals:
        raise ApcError(
            f"cannot count the CSFs of {n_electrons} electrons in {n_orbitals} orbitals: "
            "the space does not exist."
        )
    alpha = (n_electrons + 1) // 2
    beta = n_electrons // 2
    total = math.comb(n_orbitals, alpha) * math.comb(n_orbitals, beta)
    if alpha + 1 <= n_orbitals and beta >= 1:
        total -= math.comb(n_orbitals, alpha + 1) * math.comb(n_orbitals, beta - 1)
    return total


def _binary_entropy(p: float) -> float:
    """-p ln p - (1-p) ln(1-p), with the exact limits at p = 0 and 1."""
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -(p * math.log(p) + (1.0 - p) * math.log(1.0 - p))


# --- extraction --------------------------------------------------------------


def _require_closed_shell_rff(export: OrcaJson) -> None:
    """The v1 scope: closed-shell RHF exports with a single spin block."""
    if export.hftyp != "RHF":
        raise ApcError(
            f"the export reports HFTyp {export.hftyp!r}; this ranking covers closed-shell "
            "RHF exports only (the source's singly-occupied-orbital rule and its UNO "
            "variant need multi-spin exports). Next step: run the SCF as RHF and re-export "
            "with orca_2json."
        )
    for index, occupation in enumerate(export.mo_occupations):
        if min(abs(occupation), abs(occupation - 2.0)) > OCCUPATION_TOLERANCE:
            raise ApcError(
                f"orbital {index} has occupation {occupation:g}, which is neither 0 nor 2: "
                "the export is not a closed-shell RHF orbital set. Next step: run the SCF "
                "as RHF and re-export with orca_2json."
            )


def _exchange_diagonal(export: OrcaJson) -> list[float]:
    """The source's 0.5 K_aa, in ORCA's export sign: -diag(C K C^T) per orbital."""
    if export.exchange is None:
        raise ApcError(
            "the export carries no K-Matrix block (FockMatrix request). Next step: re-export "
            f"with {{{_APC_CONF}}} written as <base>.json.conf next to the .gbw."
        )
    if len(export.exchange) != 1:
        raise ApcError(
            f"the K-Matrix block carries {len(export.exchange)} spin matrices; this ranking "
            "covers single-spin (RHF) exports. Next step: run the SCF as RHF and re-export."
        )
    coefficients = np.asarray(export.mo_coefficients)
    exchange = np.asarray(export.exchange[0])
    in_mo = coefficients @ exchange @ coefficients.T
    return [-float(value) for value in np.diag(in_mo)]


def _fock_diagonal(export: OrcaJson) -> list[float]:
    """The diagonal of the full Fock in the export's own orbital basis.

    Measured on the canonical exports: the F-Matrix block equals J + K without
    the core Hamiltonian, so the full Fock is H + J + K; its MO diagonal
    reproduces the canonical orbital energies to 6e-10 (asserted in the tests).
    """
    for name, block in (
        ("H-Matrix", export.hamiltonian),
        ("J-Matrix", export.coulomb),
        ("K-Matrix", export.exchange),
    ):
        if block is None:
            raise ApcError(
                f"the export carries no {name} block, so the Fock diagonal cannot be built. "
                f"Next step: re-export with {{{_FOCK_CONF}}} written as <base>.json.conf "
                "next to the .gbw."
            )
    for name, block in (("J-Matrix", export.coulomb), ("K-Matrix", export.exchange)):
        assert block is not None  # checked above
        if len(block) != 1:
            raise ApcError(
                f"the {name} block carries {len(block)} spin matrices; this ranking covers "
                "single-spin (RHF) exports. Next step: run the SCF as RHF and re-export."
            )
    coefficients = np.asarray(export.mo_coefficients)
    hamiltonian = np.asarray(export.hamiltonian)
    assert export.coulomb is not None and export.exchange is not None
    fock = np.asarray(export.coulomb[0]) + np.asarray(export.exchange[0])
    in_mo = coefficients @ (hamiltonian + fock) @ coefficients.T
    return [float(value) for value in np.diag(in_mo)]


def _exact_pair_integrals(
    export: OrcaJson, occupied: list[int], virtuals: list[int]
) -> dict[tuple[int, int], float]:
    """The exact (ia|ia) entries per pair, from the export's MO_IAJB block."""
    if export.mo_iajb is None:
        first_virtual = virtuals[0] if virtuals else 0
        last_virtual = virtuals[-1] if virtuals else 0
        window = [occupied[0], occupied[-1], first_virtual, last_virtual, 0, 0, 0, 0]
        raise ApcError(
            "the APCX variant needs the MO_IAJB two-electron block and the export carries "
            "none. Next step: re-export with a configuration such as "
            '{ "MOCoefficients": true, "2elIntegrals": ["MO_IAJB"], "OrbWin": '
            f"{window} }} (the input window is written as eight integers, the second all "
            "zeros; the four-integer form is rejected by orca_2json)."
        )
    diagonal = {
        (i, a): value for i, j, a, b, value in export.mo_iajb if i == j and a == b
    }
    missing = [
        (i, a) for i in occupied for a in virtuals if (i, a) not in diagonal
    ]
    if missing:
        first_virtual = virtuals[0] if virtuals else 0
        last_virtual = virtuals[-1] if virtuals else 0
        window = [occupied[0], occupied[-1], first_virtual, last_virtual, 0, 0, 0, 0]
        preview = ", ".join(f"({i},{a})" for i, a in missing[:3])
        raise ApcError(
            f"the export's MO_IAJB window does not cover {len(missing)} occupied-virtual "
            f"pairs this ranking needs (for example {preview}). Next step: re-export with "
            f'\'{{ "2elIntegrals": ["MO_IAJB"], "OrbWin": {window} }}\' (eight integers, the '
            "second window all zeros); or choose a smaller virtual window."
        )
    return diagonal


# --- ranking and selection ----------------------------------------------------


def _prepare(
    export: OrcaJson,
    *,
    variant: str,
    window_size: int,
    delta: str,
) -> tuple[list[int], list[int], list[float], dict[tuple[int, int], float], str]:
    """Validate the export and gather the working quantities for the ranking."""
    _require_closed_shell_rff(export)
    if variant not in (_VARIANT_APC, _VARIANT_APCX):
        raise ApcError(f"unknown ranking variant {variant!r}; use APC or APCX.")
    if window_size < 1:
        raise ApcError(f"window_size {window_size} leaves no virtual orbital to pair with.")
    occupied = [index for index, occ in enumerate(export.mo_occupations) if occ > 1.0]
    virtual_all = [index for index, occ in enumerate(export.mo_occupations) if occ <= 1.0]
    if not occupied:
        raise ApcError("the export has no doubly occupied orbital; nothing to rank.")
    if not virtual_all:
        raise ApcError(
            "the export has no virtual orbital; the pair coefficients need at least one "
            "(the source's stability floor of two virtuals applies while dropping)."
        )
    virtual_sorted = sorted(virtual_all, key=lambda index: export.mo_energies[index])
    virtuals = sorted(virtual_sorted[:window_size])

    if delta == "energies":
        energies = list(export.mo_energies)
        delta_source = _DELTA_ENERGIES
    elif delta == "fock":
        energies = _fock_diagonal(export)
        delta_source = _DELTA_FOCK
    else:
        raise ApcError(f"unknown delta source {delta!r}; use 'energies' or 'fock'.")

    coefficients: dict[tuple[int, int], float] = {}
    if variant == _VARIANT_APC:
        # 0.5 K_aa is a property of the virtual orbital alone.
        diagonal = _exchange_diagonal(export)
        for i in occupied:
            for a in virtuals:
                numerator = diagonal[a]
                gap = energies[a] - energies[i]
                coefficients[(i, a)] = -numerator / (gap + math.hypot(gap, numerator))
    else:
        # The APCX numerator is the pair's own (ia|ia).
        exact = _exact_pair_integrals(export, occupied, virtuals)
        for (i, a), numerator in exact.items():
            gap = energies[a] - energies[i]
            coefficients[(i, a)] = -numerator / (gap + math.hypot(gap, numerator))
    return occupied, virtuals, energies, coefficients, delta_source


def pair_coefficients(
    export: OrcaJson,
    *,
    variant: str = _VARIANT_APC,
    window_size: int = DEFAULT_WINDOW,
    delta: str = "energies",
) -> tuple[tuple[int, int, float], ...]:
    """The pair coefficients c_ia (eq. 19) as ``(i, a, c)`` rows, rank order aside.

    Exposed for the regression anchors (the H2 two-configuration model compares
    them against the exact CI coefficient from the FCIDUMP).
    """
    occupied, virtuals, _, coefficients, _ = _prepare(
        export, variant=variant, window_size=window_size, delta=delta
    )
    return tuple(
        (i, a, coefficients[(i, a)]) for i in occupied for a in virtuals
    )


def rank_orbitals(
    export: OrcaJson,
    *,
    variant: str = _VARIANT_APC,
    window_size: int = DEFAULT_WINDOW,
    delta: str = "energies",
) -> Ranking:
    """Rank the candidate orbitals of one export by APC (or APCX) entropy.

    Candidates are all doubly occupied orbitals plus the ``window_size``
    virtuals lowest in energy.  ``delta`` = ``"energies"`` uses the exported
    orbital energies (canonical orbitals); ``"fock"`` uses the Fock diagonal in
    the export's own basis (non-canonical or localized sets).
    """
    occupied, virtuals, _, coefficients, delta_source = _prepare(
        export, variant=variant, window_size=window_size, delta=delta
    )
    n_virtual_available = sum(1 for occ in export.mo_occupations if occ <= 1.0)

    candidates: list[Candidate] = []
    for i in occupied:
        total = sum(coefficients[(i, a)] ** 2 for a in virtuals)
        p = total / (1.0 + total)
        candidates.append(Candidate(index=i, role=_ROLE_OCCUPIED, entropy=_binary_entropy(p)))
    for a in virtuals:
        total = sum(coefficients[(i, a)] ** 2 for i in occupied)
        p = total / (1.0 + total)
        candidates.append(Candidate(index=a, role=_ROLE_VIRTUAL, entropy=_binary_entropy(p)))
    # Entropies that agree to 1e-12 rank by index.  Physically degenerate pairs
    # (the pi_u pair of a linear molecule, for instance) differ only in the last
    # floating-point digits, which may depend on the BLAS the K-matrix transform
    # ran under; without this rule the drop order -- and, at a cap boundary, the
    # selected set -- could differ between machines.
    candidates.sort(key=lambda item: (-round(item.entropy, _RANK_RESOLUTION), item.index))

    return Ranking(
        base_name=export.base_name,
        variant=variant,
        delta_source=delta_source,
        occupied=tuple(occupied),
        virtuals=tuple(virtuals),
        candidates=tuple(candidates),
        window=(virtuals[0], virtuals[-1]),
        n_virtual_available=n_virtual_available,
    )


def select(ranking: Ranking, *, cap: int, cap_label: str | None = None) -> Selection:
    """The source's drop procedure: shed the lowest-ranked orbital until the
    CSF count fits the cap, skipping drops that would break the space floor."""
    if cap < 1:
        raise ApcError(f"a CSF cap of {cap} is not a positive count.")
    occupied = set(ranking.occupied)
    virtuals = set(ranking.virtuals)
    # Ascending importance: the drop order (ties resolved by index, the same
    # resolution rule the ranking used).
    order = [
        candidate.index
        for candidate in sorted(
            ranking.candidates, key=lambda item: (round(item.entropy, _RANK_RESOLUTION), item.index)
        )
    ]
    dropped: list[int] = []
    dropped_set: set[int] = set()
    floor_skips: list[int] = []

    def current() -> tuple[int, int, int]:
        n_occ = sum(1 for index in occupied if index not in dropped_set)
        n_virt = sum(1 for index in virtuals if index not in dropped_set)
        return n_occ, n_virt, n_csf(2 * n_occ, n_occ + n_virt)

    n_occ, n_virt, csfs = current()
    while csfs > cap:
        for index in order:
            if index in dropped_set:
                continue
            if index in occupied and n_occ - 1 < 1:
                if index not in floor_skips:
                    floor_skips.append(index)
                continue
            if index in virtuals and n_virt - 1 < 2:
                if index not in floor_skips:
                    floor_skips.append(index)
                continue
            dropped.append(index)
            dropped_set.add(index)
            n_occ, n_virt, csfs = current()
            break
        else:
            break  # nothing droppable: the cap is unreachable under the floor
    members = tuple(
        sorted(index for index in (occupied | virtuals) if index not in dropped_set)
    )
    return Selection(
        cap=cap,
        cap_label=cap_label if cap_label is not None else str(cap),
        n_electrons=2 * n_occ,
        n_orbitals=n_occ + n_virt,
        n_csfs=csfs,
        cap_reached=csfs <= cap,
        members=members,
        dropped=tuple(dropped),
        floor_skips=tuple(floor_skips),
    )


def analyze(
    export: OrcaJson,
    *,
    variant: str = _VARIANT_APC,
    window_size: int = DEFAULT_WINDOW,
    delta: str = "energies",
    cap: int = CSF_CAPS["max(8,8)"],
    cap_label: str | None = None,
) -> ApcReport:
    """Rank and select in one pass (the menu's entry point)."""
    ranking = rank_orbitals(export, variant=variant, window_size=window_size, delta=delta)
    selection = select(ranking, cap=cap, cap_label=cap_label)
    occupied_energy = tuple((i, export.mo_energies[i]) for i in ranking.occupied)
    window_energies = tuple((a, export.mo_energies[a]) for a in ranking.virtuals)
    return ApcReport(
        ranking=ranking,
        selection=selection,
        occupied_energy=occupied_energy,
        window_energies=window_energies,
    )


# --- output -------------------------------------------------------------------


def render(report: ApcReport) -> str:
    """The ranked table, the selection and the boundaries."""
    ranking, selection = report.ranking, report.selection
    lines = [
        f"Ranked-orbital active-space selection ({ranking.variant}), "
        f"{ranking.delta_source} as the model gap:",
        f"  system: {ranking.base_name}: doubly occupied {len(ranking.occupied)}, virtual "
        f"window {ranking.window[0]}-{ranking.window[1]} ({len(ranking.virtuals)} of "
        f"{ranking.n_virtual_available} virtuals, lowest in energy).",
        "",
        f"  {'rank':>4}  {'MO':>3}  {'role':<15}  {'S':>8}",
    ]
    for rank, candidate in enumerate(ranking.candidates, start=1):
        lines.append(
            f"  {rank:>4}  {candidate.index:>3}  {candidate.role:<15}  "
            f"{candidate.entropy:>8.4f}"
        )
    lines += ["", "Selection:"]
    reached = (
        f"N_CSF = {selection.n_csfs} <= {selection.cap}"
        if selection.cap_reached
        else (
            f"N_CSF = {selection.n_csfs} > {selection.cap}: the cap cannot be reached "
            "without breaking the source's two-virtual floor"
        )
    )
    lines.append(
        f"  cap {selection.cap_label} = {selection.cap} CSFs; selected "
        f"({selection.n_electrons}e, {selection.n_orbitals}o), {reached}"
    )
    lines.append(
        "  active orbitals: "
        + ", ".join(str(index) for index in selection.members)
        if selection.members
        else "  active orbitals: none"
    )
    lines.append(
        "  dropped (in drop order): "
        + (", ".join(str(index) for index in selection.dropped) if selection.dropped else "none")
    )
    if selection.floor_skips:
        lines.append(
            "  kept by the stability floor (a drop that would leave fewer than one "
            "occupied or fewer than two unoccupied orbitals): "
            + ", ".join(str(index) for index in selection.floor_skips)
        )
    lines += ["", "Boundaries and checks:"]
    for item in (
        "closed-shell RHF exports only; the source's singly-occupied-orbital rule "
        "(assign the maximum approximated entropy) and its UNO variant are not covered",
        "APC systematically overestimates the doubly-occupied orbital entropies "
        "(R^2 = 0.64, MAE 0.0240 against DMRG; virtuals 0.83/0.0064; ranking precision "
        "about 88% on the source's set)",
        "the source expects the scheme to perform worse in much larger systems and where "
        "the HF determinant is a poor approximation; the ranking is a screening device, "
        "not a converged answer",
        "the entropies and the ranking depend on the candidate window; the window is "
        "printed above, and N_CSF is checked with eq. (2) at every drop",
    ):
        lines.append(f"  - {item}")
    return "\n".join(lines)


def run(
    export: OrcaJson,
    *,
    variant: str = _VARIANT_APC,
    window_size: int = DEFAULT_WINDOW,
    delta: str = "energies",
    cap: int = CSF_CAPS["max(8,8)"],
    cap_label: str | None = None,
) -> ReportSection:
    """The analyser entry point for menu 21."""
    return ReportSection(
        title=f"3.1 ranked-orbital active-space selection ({variant})",
        body=render(
            analyze(
                export,
                variant=variant,
                window_size=window_size,
                delta=delta,
                cap=cap,
                cap_label=cap_label,
            )
        ),
    )


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the framework, the approximation and the measured interface."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "The ranked-orbital scheme (rank by importance, trim to a CSF cap with "
                "eq. 2, keep at least one occupied and two unoccupied orbitals) and the "
                "approximate pair coefficient: the two-configuration pair coefficient "
                "c = -(12|12)/(Delta + sqrt((12|12)^2 + Delta^2)) with (ia|ia) ~ 0.5 K_aa "
                "and Delta ~ eps_a - eps_i (eqs. 9-19), the one-orbital entropies of "
                "eqs. 12-13, and the APCX exact-integral variant. The source's own error "
                "record: APC overestimates doubly occupied orbitals (R^2 = 0.64, "
                "MAE 0.0240 vs DMRG; virtuals 0.83/0.0064; precision about 88%), performs "
                "worse than APC for HF orbitals as APCX, and is expected to degrade in "
                "much larger systems or with a poor HF reference."
            ),
            ref=(
                "J. Chem. Theory Comput. 2021, 17, 2817-2831, sections 3.2 and 3.12 "
                "(eqs. 2, 9-19 and the performance discussion)"
            ),
            url="https://doi.org/10.1021/acs.jctc.1c00037",
            bibkey="king2021ranked",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "Interface conventions measured on ORCA 6.1.1 exports: the exported "
                "K-Matrix equals the exchange contribution as it enters the Fock, so the "
                "source's 0.5 K_aa is -diag(C K C^T)[a] = sum_i (a i | a i) (matched to "
                "4e-12 against the MO_IAJB entries on the N2/def2-SVP export, and to "
                "7e-16 on the H2/STO-3G two-configuration model, where eq. 18 is exact); "
                "diag(C (H + J + K) C^T) reproduces the orbital energies to 6e-10 (the "
                "exported F-Matrix block excludes the core Hamiltonian); the MO_IAJB "
                "entries are (i, j, a, b, (ia|jb)) with the internal indices first, and "
                "the input OrbWin line must be written as eight integers. Regression "
                "anchors: the H2 exact pair coefficient from the FCIDUMP 2x2 CI matches "
                "the engine's CASSCF energy to 12 digits and the eq. 14 closed form to "
                "2e-17; the N2 max(10,10) selection lands exactly on the source's "
                "(10,10) = 19404 table entry."
            ),
            ref="tests/test_apc.py; fixtures/orca/h2_apc.json; fixtures/orca/n2_apc.json",
        ),
    )
