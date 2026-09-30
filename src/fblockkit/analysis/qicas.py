"""QICAS -- quantum-information-assisted active-space optimization.

The out-of-CAS correlation of a CAS(N_CAS, D_CAS) scheme in an ordered orbital
basis B is the sum of the single-orbital (four-state) entropies over the
*non-active* orbitals,

    F_QI(B) = sum_{i not in A} S(rho_i)                                  (eq. 2)

and the source's central result (their Theorem 1) bounds the CASCI error from
above by it: ``E_CASCI(B) - E_FCI <= Delta E_max / ln(4) * F_QI(B)``.  F_QI
depends only on the 1- and 2-RDM -- not on the Hamiltonian -- so the active
space can be *optimized* by rotating the orbitals to minimize F_QI, which is
what QICAS (Ding, Knecht & Schilling, JPCL 2023) does.  The minimization is
the source's gradient-free Jacobi sweep (their Appendix B): for every orbital
pair touching a non-active orbital, scan the rotation angle in [0, pi) with a
1e-2 grid, refine the best angle at 1e-4, accept the rotation when the
improvement exceeds eps1, rotate the density objects, and repeat the sweep for
up to N_cycle cycles (or until a cycle improves by less than 10*eps1).  Every
non-active orbital whose occupancy settles above (below) 1 is read as closed
(virtual).  The source's optional size-selection variant (Appendix C, minimize
the *total* orbital entropy and read the plateau of the resulting threshold
diagram) is not implemented here.

This tool's variant, stated honestly
------------------------------------

The source drives QICAS with a low-bond-dimension DMRG ground state over the
full orbital space (or a large subset of it).  This toolbox has no DMRG: its
(1-, 2-RDM) come from the *exact* determinant CI of the FCIDUMP route (the
same machinery as the four-state entropy, menu 12), which covers the FCIDUMP's
orbital window.  QICAS is therefore applied *on the window* -- the source's
own subset application, with exact instead of approximate RDMs inside it, and
the partition chosen as follows: the dumped orbital order is re-sorted by
natural occupation, the k = (N - N_CAS)/2 highest occupations become the
closed block, the lowest become the virtual block, and the rest is active
(the source fixes the blocks as index ranges of an ordered basis; the
occupation re-sorting is this tool's documented choice of that order).  The
window restriction means F_QI sums over the window's non-active orbitals only
-- it is not comparable with the source's full-space numbers.

What the report verifies on real data
-------------------------------------

- the exact window FCI energy of the CI solution (cross-checkable against the
  run's printed CASSCF energy, the menu-12 anchor);
- F_QI and the per-orbital entropy profile before and after the optimization;
- the CASCI(N_CAS, D_CAS) energy in the initial and the optimized basis (the
  optimized CASCI is the source's headline deliverable: it approaches the
  window FCI);
- the Theorem-1 bound ``Delta E <= Delta E_max / ln(4) * F_QI`` evaluated with
  the numbers (Delta E_max from the CI spectrum's span; the Ms sector's span
  is a lower bound of the full Delta E_max, so a passing check is stricter
  than the theorem guarantees, never weaker).

The regression anchors are the N2/def2-SVP CAS(6,6) FCIDUMP fixture at the
(4,4)-in-(6,6) and (2,4)-in-(6,6) targets; tests/test_qicas.py pins the
measured numbers.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from ..knowledge.models import EVIDENCE_LITERATURE, EVIDENCE_MEASURED, Evidence, ReportSection
from ..parsers.fcidump import Fcidump
from . import entropy_rdm
from .entropy_rdm import EntropyRdmError, SpinDensities, orbital_entropy, spin_densities

__all__ = [
    "QicasError",
    "CSpace",
    "RotationReport",
    "QicasResult",
    "partition_by_occupation",
    "out_of_cas_correlation",
    "minimize_rotations",
    "rotate_integrals",
    "casci_energy",
    "analyze",
    "render",
    "run",
    "evidence",
]

#: The source's acceptance threshold for one Jacobi rotation (its C2 value;
#: the Cr2 run used 1e-7).
QICAS_EPSILON = 1e-8

#: The source's cap on full sweeps through the pair list.
QICAS_MAX_CYCLES = 200


class QicasError(ValueError):
    """The QICAS run cannot proceed on the given data (with a next step)."""


# --- models -------------------------------------------------------------------


@dataclass(frozen=True)
class CSpace:
    """The closed / active / virtual partition of the window (index tuples)."""

    closed: tuple[int, ...]
    active: tuple[int, ...]
    virtual: tuple[int, ...]

    @property
    def n_orbitals(self) -> int:
        return len(self.active)


@dataclass(frozen=True)
class RotationReport:
    """One optimization run's bookkeeping."""

    pairs_mode: str
    accepted: int
    cycles: int
    fqi_initial: float
    fqi_final: float
    epsilon1: float
    matrix: tuple[tuple[float, ...], ...]


@dataclass(frozen=True)
class QicasResult:
    """Everything the report renders (and the tests pin)."""

    n_window: int
    n_window_electrons: int
    space: CSpace
    space_final: CSpace
    occupations: tuple[float, ...]
    entropy_initial: tuple[float, ...]
    entropy_final: tuple[float, ...]
    rotation: RotationReport
    energy_fci: float
    energy_casci_initial: float
    energy_casci_final: float
    gap_to_fci: float
    theorem_bound: float
    final_occupations: tuple[float, ...]
    source_energy_checked: bool


# --- the partition ------------------------------------------------------------


def partition_by_occupation(
    occupations: tuple[float, ...], n_cas: int, n_active_orbitals: int
) -> CSpace:
    """Order the window by occupation and cut closed / active / virtual blocks."""
    n_window = len(occupations)
    if (n_window - n_cas) % 2 != 0:
        raise QicasError(
            f"the window holds {n_window} electrons and the target active space takes "
            f"{n_cas}; the difference must be even (closed orbitals hold pairs)."
        )
    k_closed = (n_window - n_cas) // 2
    n_virtual = n_window - k_closed - n_active_orbitals
    if k_closed < 0 or n_virtual < 0 or n_active_orbitals < 1:
        raise QicasError(
            f"the target ({n_cas}e, {n_active_orbitals}o) does not fit the window: it "
            f"needs {k_closed} closed and {n_virtual} virtual orbitals beside the "
            f"{n_active_orbitals} active ones. Next step: choose a target with "
            f"0 <= closed, 0 <= virtual (window {n_window} orbitals, "
            f"{n_window} electrons)."
        )
    order = sorted(range(n_window), key=lambda i: (-occupations[i], i))
    closed = tuple(sorted(order[:k_closed]))
    active = tuple(sorted(order[k_closed:k_closed + n_active_orbitals]))
    virtual = tuple(sorted(order[k_closed + n_active_orbitals:]))
    return CSpace(closed=closed, active=active, virtual=virtual)


def out_of_cas_correlation(dens: SpinDensities, nonactive: tuple[int, ...]) -> float:
    """F_QI: the sum of four-state entropies over the non-active orbitals."""
    return float(sum(orbital_entropy(dens, i) for i in nonactive))


# --- the Jacobi optimization --------------------------------------------------


def _weights_entropy(n_a: float, n_b: float, p2: float) -> float:
    """The four-state entropy from (n_up, n_down, double occupancy).

    Mirrors :func:`entropy_rdm.orbital_weights` (round-off snapping and the
    inconsistency check); a unit test asserts the equivalence on slices.
    """
    weights = (1.0 - n_a - n_b + p2, n_a - p2, n_b - p2, p2)
    total = 0.0
    for value in weights:
        if value < 0.0:
            if value < -1e-8:
                raise QicasError(
                    f"a negative four-state weight ({value:.3e}) appeared during the "
                    "rotation scan: the density objects are inconsistent."
                )
            value = 0.0
        if value > 1e-14:
            total -= value * math.log(value)
    return float(total)


def _pair_entropies(dens: SpinDensities, i: int, j: int, theta: float) -> tuple[float, float]:
    """The entropies of the two orbitals after rotating the (i, j) pair by theta.

    Exact and O(1): a rotation mixing only i and j changes the single-orbital
    densities of i and j alone, and the new values follow from the 2x2 slices
    of the density objects.
    """
    c, s = math.cos(theta), math.sin(theta)
    t = np.array([[c, -s], [s, c]])
    idx = np.array([i, j])
    g_a = dens.gamma_a[np.ix_(idx, idx)]
    g_b = dens.gamma_b[np.ix_(idx, idx)]
    g_2 = dens.g_ab[np.ix_(idx, idx, idx, idx)]
    n_a = np.einsum("pn,qn,pq->n", t, t, g_a)
    n_b = np.einsum("pn,qn,pq->n", t, t, g_b)
    p2 = np.einsum("pn,qn,rn,sn,pqrs->n", t, t, t, t, g_2)
    return (
        _weights_entropy(float(n_a[0]), float(n_b[0]), float(p2[0])),
        _weights_entropy(float(n_a[1]), float(n_b[1]), float(p2[1])),
    )


def _rotation_two(n: int, i: int, j: int, theta: float) -> np.ndarray:
    c, s = math.cos(theta), math.sin(theta)
    t = np.eye(n)
    t[i, i], t[i, j] = c, -s
    t[j, i], t[j, j] = s, c
    return t


def _best_angle(dens: SpinDensities, i: int, j: int, weights: tuple[float, float]) -> tuple[float, float]:
    """Scan theta for one pair; return (best theta, improvement over theta = 0).

    The scan mirrors the source's Appendix B: a 1e-2 grid over [0, pi), then a
    1e-4 refinement in the winner's neighbourhood.
    """
    base = _pair_entropies(dens, i, j, 0.0)
    base_value = weights[0] * base[0] + weights[1] * base[1]

    def value(theta: float) -> float:
        pair = _pair_entropies(dens, i, j, theta)
        return weights[0] * pair[0] + weights[1] * pair[1]

    grid = np.arange(0.0, math.pi, 1e-2)
    values = np.array([value(theta) for theta in grid])
    best = float(grid[int(np.argmin(values))])
    fine = np.arange(best - 1e-2, best + 1e-2 + 1e-6, 1e-4)
    fine = fine[(fine >= 0.0) & (fine < math.pi)]
    if fine.size:
        fine_values = np.array([value(theta) for theta in fine])
        candidate = float(fine[int(np.argmin(fine_values))])
        if value(candidate) < value(best):
            best = candidate
    return best, base_value - value(best)


def minimize_rotations(
    dens: SpinDensities,
    nonactive: tuple[int, ...],
    *,
    pairs_mode: str = "touch",
    epsilon1: float = QICAS_EPSILON,
    max_cycles: int = QICAS_MAX_CYCLES,
    seed: int = 0,
) -> tuple[SpinDensities, RotationReport]:
    """Minimize F_QI by the source's Jacobi sweep over orbital rotations.

    ``pairs_mode``: ``"touch"`` rotates every pair touching a non-active
    orbital (the source's C2 choice); ``"exclusive"`` only active/non-active
    pairs (their Cr2 economy variant; F_QI is *not* invariant under rotations
    within the non-active space, which is why the default includes them).
    The pair order is shuffled per cycle with a fixed seed, so a report is
    reproducible (the source's shuffle is stochastic; ours is not).
    """
    if pairs_mode not in ("touch", "exclusive"):
        raise QicasError(f"unknown rotation set {pairs_mode!r}; use touch or exclusive.")
    n = dens.norb
    nonactive_set = set(nonactive)
    if pairs_mode == "touch":
        allowed = lambda i, j: (i in nonactive_set) or (j in nonactive_set)  # noqa: E731
    else:
        allowed = lambda i, j: (i in nonactive_set) != (j in nonactive_set)  # noqa: E731
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n) if allowed(i, j)]

    fqi_initial = out_of_cas_correlation(dens, nonactive)
    u = np.eye(n)
    accepted = 0
    cycles = 0
    if not pairs or not nonactive:
        report = RotationReport(
            pairs_mode=pairs_mode,
            accepted=0,
            cycles=0,
            fqi_initial=fqi_initial,
            fqi_final=fqi_initial,
            epsilon1=epsilon1,
            matrix=tuple(tuple(float(v) for v in row) for row in u),
        )
        return dens, report

    rng = np.random.default_rng(seed)
    current = dens
    fqi = fqi_initial
    improvement_stop = 10.0 * epsilon1
    for cycle in range(max_cycles):
        cycles = cycle + 1
        fqi_at_start = fqi
        for k in rng.permutation(len(pairs)):
            i, j = pairs[int(k)]
            weights = (1.0 if i in nonactive_set else 0.0, 1.0 if j in nonactive_set else 0.0)
            theta, improvement = _best_angle(current, i, j, weights)
            if improvement > epsilon1:
                t = _rotation_two(n, i, j, theta)
                current = entropy_rdm.rotate_densities(current, t)
                # In the ``u[old, new]`` convention sequential rotations compose
                # as ``u_total = t_1 @ t_2 @ ...`` (the new transform multiplies
                # from the right); a unit test pins this against a replay of
                # rotate_densities with the accumulated matrix.
                u = u @ t
                fqi -= improvement
                accepted += 1
        if fqi_at_start - fqi < improvement_stop:
            break
    report = RotationReport(
        pairs_mode=pairs_mode,
        accepted=accepted,
        cycles=cycles,
        fqi_initial=fqi_initial,
        fqi_final=out_of_cas_correlation(current, nonactive),
        epsilon1=epsilon1,
        matrix=tuple(tuple(float(v) for v in row) for row in u),
    )
    return current, report


# --- integrals and the CASCI energy -------------------------------------------


def _dense_eri(dump: Fcidump) -> np.ndarray:
    n = dump.norb
    g = np.zeros((n, n, n, n))
    for (p, q, r, s), value in dump.g.items():
        g[p, q, r, s] = value
    return g


def rotate_integrals(dump: Fcidump, u: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rotate the FCIDUMP integrals to the new basis (``u[old, new]``)."""
    h = np.asarray(dump.h, dtype=np.float64)
    g = _dense_eri(dump)
    h_new = u.T @ h @ u
    g_new = np.einsum("pi,qj,rk,sl,pqrs->ijkl", u, u, u, u, g)
    return h_new, g_new


def casci_energy(
    dump: Fcidump,
    u: np.ndarray,
    space: CSpace,
    *,
    reference_energy: float | None = None,
) -> float:
    """The CASCI(N_CAS, D_CAS) energy in the basis rotated by ``u``.

    Closed orbitals are folded into an effective one-electron operator (and a
    core energy); the active block is solved exactly with the same determinant
    solver the four-state entropy uses.
    """
    h, g = rotate_integrals(dump, u)
    closed = list(space.closed)
    active = list(space.active)
    n_act = len(active)
    k_closed = len(closed)
    n_cas = dump.nelec - 2 * k_closed

    # Fock folding: h_eff = h + sum_i [2 (pq|ii) - (pi|iq)] over the closed set
    # (chemist notation; the same convention the CASSCF Hamiltonian uses).
    h_eff = np.array(h)
    for i in closed:
        h_eff += 2.0 * g[:, :, i, i]
        h_eff -= g[:, i, i, :]
    ecore = float(dump.ecore + 2.0 * sum(h[i, i] for i in closed))
    for i in closed:
        for j in closed:
            ecore += 2.0 * float(g[i, i, j, j]) - float(g[i, j, j, i])

    h_block = tuple(
        tuple(float(h_eff[p, q]) for q in active) for p in active
    )
    g_block: dict[tuple[int, int, int, int], float] = {}
    for a_i, p in enumerate(active):
        for a_j, q in enumerate(active):
            for a_k, r in enumerate(active):
                for a_l, s in enumerate(active):
                    g_block[(a_i, a_j, a_k, a_l)] = float(g[p, q, r, s])
    effective = Fcidump(
        norb=n_act,
        nelec=n_cas,
        ms2=dump.ms2,
        ecore=ecore,
        h=h_block,
        g=g_block,
        orbsym=(1,) * n_act,
        isym=dump.isym,
    )
    state = entropy_rdm.solve_fci(effective, reference_energy=reference_energy)
    return float(state.energy_total)


# --- the analysis -------------------------------------------------------------


def analyze(
    dump: Fcidump,
    *,
    n_cas: int,
    n_active_orbitals: int,
    pairs_mode: str = "touch",
    epsilon1: float = QICAS_EPSILON,
    max_cycles: int = QICAS_MAX_CYCLES,
    seed: int = 0,
    reference_energy: float | None = None,
) -> QicasResult:
    """One QICAS run: ground-state RDMs, the partition, and the optimization."""
    if n_cas < 0 or n_cas > dump.nelec:
        raise QicasError(
            f"the target active space takes {n_cas} electrons from the window's "
            f"{dump.nelec}. Next step: choose 0 <= N_CAS <= N."
        )
    try:
        state = entropy_rdm.solve_fci(dump, reference_energy=reference_energy)
        dens = spin_densities(state)
    except EntropyRdmError as exc:
        raise QicasError(str(exc)) from exc
    occupations = tuple(
        float(dens.gamma_a[i, i] + dens.gamma_b[i, i]) for i in range(dump.norb)
    )
    space = partition_by_occupation(occupations, n_cas, n_active_orbitals)
    nonactive = tuple(sorted(space.closed + space.virtual))

    entropy_before = entropy_rdm.entropy_spectrum(dens)
    fqi_before = out_of_cas_correlation(dens, nonactive)
    energy_casci_initial = casci_energy(dump, np.eye(dump.norb), space)
    energy_fci = float(state.energy_total)

    optimized, rotation = minimize_rotations(
        dens,
        nonactive,
        pairs_mode=pairs_mode,
        epsilon1=epsilon1,
        max_cycles=max_cycles,
        seed=seed,
    )
    entropy_after = entropy_rdm.entropy_spectrum(optimized)
    final_occupations = tuple(
        float(optimized.gamma_a[i, i] + optimized.gamma_b[i, i]) for i in range(dump.norb)
    )
    # The source's reading of the optimized basis: a non-active orbital counts as
    # closed (virtual) if its occupancy is larger (smaller) than 1.  The
    # optimizer may swap the closed and virtual slot contents (a rotation that
    # lowers F_QI), so the final partition is re-derived here and the CASCI
    # check uses it; its shape may differ from the requested target.
    closed_final = tuple(
        sorted(i for i in nonactive if final_occupations[i] > 1.0)
    )
    virtual_final = tuple(
        sorted(i for i in nonactive if final_occupations[i] <= 1.0)
    )
    space_final = CSpace(closed=closed_final, active=space.active, virtual=virtual_final)
    energy_casci_final = casci_energy(dump, np.asarray(rotation.matrix), space_final)

    delta_after = energy_casci_final - energy_fci
    delta_max = float(state.spectrum[-1] - state.spectrum[0]) if state.spectrum else 0.0
    bound = delta_max / math.log(4.0) * rotation.fqi_final

    return QicasResult(
        n_window=dump.norb,
        n_window_electrons=dump.nelec,
        space=space,
        space_final=space_final,
        occupations=occupations,
        entropy_initial=entropy_before,
        entropy_final=entropy_after,
        rotation=rotation,
        energy_fci=energy_fci,
        energy_casci_initial=energy_casci_initial,
        energy_casci_final=energy_casci_final,
        gap_to_fci=delta_after,
        theorem_bound=bound,
        final_occupations=final_occupations,
        source_energy_checked=reference_energy is not None,
    )


# --- output -------------------------------------------------------------------


def render(result: QicasResult) -> str:
    """The entropy profiles, the rotation bookkeeping and the energy checks."""
    space = result.space
    lines = [
        "QICAS orbital optimization (out-of-CAS correlation over orbital rotations):",
        f"  window: {result.n_window} orbitals, {result.n_window_electrons} electrons; "
        f"target ({result.n_window_electrons - 2 * len(space.closed)}e, "
        f"{space.n_orbitals}o) with "
        f"{len(space.closed)} closed, {len(space.active)} active, "
        f"{len(space.virtual)} virtual",
        f"  requested partition (by natural occupation): closed {list(space.closed)}, "
        f"active {list(space.active)}, virtual {list(space.virtual)}",
        "",
        f"  {'orbital':>7}  {'occupation':>10}  {'S initial':>9}  {'S optimized':>11}",
    ]
    for i in range(result.n_window):
        lines.append(
            f"  {i:>7}  {result.occupations[i]:>10.5f}  "
            f"{result.entropy_initial[i]:>9.4f}  {result.entropy_final[i]:>11.4f}"
        )
    rotation = result.rotation
    lines += [
        "",
        f"  F_QI (non-active entropies): {rotation.fqi_initial:.4f} -> "
        f"{rotation.fqi_final:.4f}  ({rotation.accepted} rotations accepted in "
        f"{rotation.cycles} cycle(s), set '{rotation.pairs_mode}', eps1 = "
        f"{rotation.epsilon1:g})",
        "",
        "Energy checks (window: exact CI; CASCI in the partition):",
        f"  E(window FCI)          = {result.energy_fci:.9f} Eh",
        f"  E(CASCI, initial)      = {result.energy_casci_initial:.9f} Eh  "
        f"(+{result.energy_casci_initial - result.energy_fci:.6e} vs FCI)",
        f"  E(CASCI, optimized)    = {result.energy_casci_final:.9f} Eh  "
        f"(+{result.gap_to_fci:.6e} vs FCI)",
        f"  Theorem-1 check: gap {result.gap_to_fci:.3e} <= dE_max/ln4 * F_QI = "
        f"{result.theorem_bound:.3e} (holds)"
        if result.gap_to_fci <= result.theorem_bound
        else (
            f"  Theorem-1 check: the gap {result.gap_to_fci:.3e} exceeds the bound "
            f"{result.theorem_bound:.3e} -- the source's proof assumes the reference "
            "state overlap deficit is below 1/2, which a target this far from the "
            "window ground state can leave; the inequality is then not informative "
            "for this partition (it is reported, not hidden)"
        ),
    ]
    final = result.space_final
    lines.append(
        "  final partition (occupancy reading of the optimized basis): closed "
        f"{list(final.closed)}, active {list(final.active)}, virtual "
        f"{list(final.virtual)}"
        + (
            "  [same shape as requested]"
            if len(final.closed) == len(space.closed)
            and len(final.virtual) == len(space.virtual)
            else (
                f"  [the optimizer re-shuffled the non-active slots: the final space "
                f"reads as ({result.n_window_electrons - 2 * len(final.closed)}e, "
                f"{len(final.active)}o)]"
            )
        )
    )
    lines += ["", "Boundaries and checks:"]
    for item in (
        "the RDMs are the *exact* ground state of the FCIDUMP window (the four-state "
        "entropy route); the source drives QICAS with a low-bond-dimension DMRG ground "
        "state over the full space -- this is the source's own subset application, so "
        "F_QI here sums over the window's non-active orbitals only and is not "
        "comparable with the source's full-space values",
        "the partition is this tool's documented ordering choice (the dumped orbitals "
        "re-sorted by occupation); the source fixes closed/active/virtual as index "
        "ranges of an ordered basis",
        "the optimized rotation is reported as a matrix and used for the CASCI check; "
        "writing it back into a .gbw (the mkl route of menu 18) is not implemented",
        "the source's size-selection variant (minimize the total orbital entropy and "
        "read the plateau of the threshold diagram; their Appendix C) is not "
        "implemented",
        f"the run's printed CASSCF energy cross-check: "
        + ("applied (the CI root was matched to it)" if result.source_energy_checked
           else "not applied (no output file was given; the lowest root of the "
                "FCIDUMP's Ms sector was used)"),
    ):
        lines.append(f"  - {item}")
    return "\n".join(lines)


def run(dump: Fcidump, **kwargs) -> ReportSection:
    """The analyser entry point for menu 24."""
    return ReportSection(
        title="3.4 QICAS active-space optimization",
        body=render(analyze(dump, **kwargs)),
    )


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the scheme, its theorem and the measured interface."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "QICAS: the out-of-CAS correlation F_QI(B) = sum of non-active orbital "
                "entropies; Theorem 1 (Delta E <= Delta E_max / ln4 * F_QI) makes it a "
                "valid cost function; the gradient-free Jacobi scan over orbital "
                "rotations (Appendix B: 1e-2 grid, 1e-4 refinement, eps1 acceptance, "
                "N_cycle sweeps; 'touch' pairs for chemical accuracy, active/non-active "
                "pairs for the economical variant) and the occupancy classification of "
                "the non-active orbitals. Source headline: C2/cc-pVDZ CAS(8,8) reaches "
                "the CASSCF energy within chemical accuracy (F_QI 0.64 -> 0.38 over the "
                "full space; 39/41 dissociation points within 1.6 mHa); Cr2 CASSCF "
                "converges in 1-2 macro iterations from the QICAS orbitals."
            ),
            ref=(
                "J. Phys. Chem. Lett. 2023, 14, 11022-11029, eqs. (2)-(8), Theorem 1, "
                "and Appendices A-C"
            ),
            url="https://doi.org/10.1021/acs.jpclett.3c02536",
            bibkey="ding2023qicas",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "Implementation variant and regression anchors on the N2/def2-SVP "
                "CAS(6,6) FCIDUMP: the RDMs are exact (four-state entropy route), the "
                "partition is occupation-ordered, the F_QI candidate evaluation during "
                "the Jacobi scan is an exact 2-slice contraction (checked against a "
                "full density rotation in the tests), and the Theorem-1 inequality is "
                "verified with the run's own numbers on the (4,4)-in-(6,6) and "
                "(2,4)-in-(6,6) targets."
            ),
            ref="tests/test_qicas.py; fixtures/orca/n2_fcidump.fcidump",
        ),
    )
