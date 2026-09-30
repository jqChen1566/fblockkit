"""TNASS -- subset selection by the Renyi-2 entropy of the bipartition.

TNASS (Mingare, Heuzé & Coveney, arXiv 2026) selects the active space as the
subset ``A`` of spatial orbitals that **maximizes the Renyi-2 entropy**
``S2(A) = -log Tr(rho_A^2)`` of the bipartition between ``A`` and its
complement -- a genuinely multi-orbital entanglement measure that goes beyond
the single-orbital entropy of the AutoCAS family.  Their entanglement feature
(Eq. 1: ``|EF> = sum_b e^{-S2(b)} |b>``, the amplitudes being subsystem
purities) carries ``S2`` for every partition of the full system and is built
as a tensor network from a low-bond-dimension DMRG state; the selection
algorithms on top of that oracle are pure combinatorics: exact brute force
over ``C(N, n)`` subsets, the greedy method of their Eq. (11) (start empty,
repeatedly add the orbital that maximizes ``S2(A u {i})``), and the block
greedy variant that mixes single-orbital-entropy initialization with greedy
blocks of size ``k`` (``k = 1`` reducing to the single-orbital ranking,
``k = n`` to plain greedy).

This tool's variant, stated honestly
------------------------------------

There is no tensor network here: the ``S2`` oracle is computed *exactly* over
the FCIDUMP window from the same determinant CI the four-state-entropy route
uses.  For a spatial subset ``A``, the reduced state is assembled by grouping
the CI vector by the occupation pattern of ``A`` (the environment index being
the complement's pattern), and

    ``Tr(rho_A^2) = Tr((V V^T)^2)``,  ``V[sigma, env] = c(sigma, env)``

so the purity -- and ``S2`` -- follows from an exact contraction of the wave
function; the fermionic signs are carried by the determinant list itself
(nothing is reordered).  The source works over the full orbital space with an
approximate MPS whose useful bond dimensions are 4-6 (their own scaling
table); the exact route is the counterpart restricted to the window the
FCIDUMP covers, and it removes the bond-dimension truncation *inside* that
window.

The regression anchors: the one-orbital ``S2`` reproduces
``-log(sum_k w_k^2)`` from the four-state weights (the menu-12 densities) to
1e-12; the two-orbital ``S2`` reproduces the purity assembled independently
from the spin-resolved 2-RDM; the full-window subset has ``S2 = 0`` (a pure
state); and on the benzene pi platform the greedy four-orbital selection lands
on the strongly entangled quartet of the manifold (the exact values are pinned
in the tests).
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass

import numpy as np

from ..knowledge.models import EVIDENCE_LITERATURE, EVIDENCE_MEASURED, Evidence, ReportSection
from ..parsers.fcidump import Fcidump
from . import entropy_rdm
from .entropy_rdm import EntropyRdmError, FciState

__all__ = [
    "TnassError",
    "BRUTE_FORCE_CAP",
    "GreedyStep",
    "TnassResult",
    "renyi2_entropy",
    "greedy_selection",
    "block_greedy_selection",
    "brute_selection",
    "analyze",
    "render",
    "run",
    "evidence",
]

#: The exact combinatorial search is capped at this many subsets (the source's
#: brute force scales as C(N, n); beyond the cap the greedy family is the route).
BRUTE_FORCE_CAP = 20000


class TnassError(ValueError):
    """The selection cannot proceed on the given data (with a next step)."""


@dataclass(frozen=True)
class GreedyStep:
    """One selection step: which orbital joined, and the subset's S2."""

    added: int
    subset: tuple[int, ...]
    s2: float


@dataclass(frozen=True)
class TnassResult:
    """The selection trace, ready to render."""

    method: str
    n_target: int
    block_size: int | None
    steps: tuple[GreedyStep, ...]
    single_s2: tuple[float, ...]  # S2 of every one-orbital subset (the seed ranking)
    selected: tuple[int, ...]
    n_electrons: int
    energy_fci: float


# --- the exact S2 oracle ------------------------------------------------------


def renyi2_entropy(state: FciState, subset: tuple[int, ...]) -> float:
    """S2(A) = -log Tr(rho_A^2) for a spatial subset, exactly from the CI vector."""
    norb = state.norb
    if not subset:
        raise TnassError("the empty subset has no bipartition entropy.")
    for index in subset:
        if not 0 <= index < norb:
            raise TnassError(
                f"orbital index {index} is outside the window's {norb} orbitals."
            )
    if len(subset) == norb:
        return 0.0  # the whole system is pure against its empty complement
    mask_a = 0
    for index in subset:
        mask_a |= 1 << index
    mask_env = ((1 << norb) - 1) & ~mask_a

    # Group the CI vector by (A-pattern, environment-pattern).  The determinant
    # list already carries the fermionic signs, so grouping is exact.
    groups: dict[tuple[int, int], dict[tuple[int, int], float]] = {}
    for coefficient, (alpha, beta) in zip(state.coefficients, state.determinants):
        if coefficient == 0.0:
            continue
        sigma = (alpha & mask_a, beta & mask_a)
        env = (alpha & mask_env, beta & mask_env)
        groups.setdefault(sigma, {})[env] = float(coefficient)

    sigmas = sorted(groups)
    envs = sorted({env for group in groups.values() for env in group})
    env_index = {env: i for i, env in enumerate(envs)}
    matrix = np.zeros((len(sigmas), len(envs)))
    for i, sigma in enumerate(sigmas):
        for env, value in groups[sigma].items():
            matrix[i, env_index[env]] = value
    if matrix.shape[0] <= matrix.shape[1]:
        product = matrix @ matrix.T
    else:
        product = matrix.T @ matrix
    purity = float((product**2).sum())
    if purity <= 0.0:
        raise TnassError("a zero purity appeared; the CI vector is not normalised.")
    return -math.log(purity)


def _solve(dump: Fcidump, reference_energy: float | None) -> FciState:
    try:
        return entropy_rdm.solve_fci(dump, reference_energy=reference_energy)
    except EntropyRdmError as exc:
        raise TnassError(str(exc)) from exc


# --- the selection algorithms -------------------------------------------------


def greedy_selection(
    state: FciState, n_target: int, *, seed: tuple[int, ...] = ()
) -> list[GreedyStep]:
    """The source's Eq. (11): grow the subset by the maximal S2 increase."""
    norb = state.norb
    if not 1 <= n_target <= norb:
        raise TnassError(
            f"the target size {n_target} is outside 1..{norb} (the window's orbitals)."
        )
    subset = tuple(seed)
    if any(not 0 <= index < norb for index in subset):
        raise TnassError("the seed contains an orbital outside the window.")
    steps: list[GreedyStep] = []
    while len(subset) < n_target:
        best_index, best_value = None, -math.inf
        for candidate in range(norb):
            if candidate in subset:
                continue
            value = renyi2_entropy(state, tuple(sorted(subset + (candidate,))))
            if value > best_value + 1e-12:
                best_index, best_value = candidate, value
        subset = tuple(sorted(subset + (best_index,)))  # type: ignore[arg-type]
        steps.append(GreedyStep(added=int(best_index), subset=subset, s2=best_value))
    return steps


def block_greedy_selection(
    state: FciState, n_target: int, block_size: int
) -> list[GreedyStep]:
    """The source's block greedy: entropy-initialised blocks of greedy growth.

    ``block_size = 1`` reduces to the single-orbital-S2 ranking; a block of
    size ``k`` starts from the highest-S2 unchosen orbital and grows greedily
    for ``k`` orbitals in total, then restarts.
    """
    norb = state.norb
    if not 1 <= n_target <= norb:
        raise TnassError(
            f"the target size {n_target} is outside 1..{norb} (the window's orbitals)."
        )
    if not 1 <= block_size <= max(1, n_target):
        raise TnassError(
            f"the block size {block_size} must lie in 1..{n_target}."
        )
    singles = {
        index: renyi2_entropy(state, (index,)) for index in range(norb)
    }
    chosen: tuple[int, ...] = ()
    steps: list[GreedyStep] = []
    while len(chosen) < n_target:
        remaining = [i for i in range(norb) if i not in chosen]
        block_start = max(remaining, key=lambda i: (singles[i], -i))
        block = (block_start,)
        chosen = tuple(sorted(chosen + block))
        steps.append(GreedyStep(added=block_start, subset=chosen, s2=renyi2_entropy(state, chosen)))
        for _ in range(block_size - 1):
            if len(chosen) >= n_target:
                break
            best_index, best_value = None, -math.inf
            for candidate in range(norb):
                if candidate in chosen:
                    continue
                value = renyi2_entropy(state, tuple(sorted(chosen + (candidate,))))
                if value > best_value + 1e-12:
                    best_index, best_value = candidate, value
            chosen = tuple(sorted(chosen + (best_index,)))  # type: ignore[arg-type]
            steps.append(GreedyStep(added=int(best_index), subset=chosen, s2=best_value))
    return steps


def brute_selection(state: FciState, n_target: int) -> list[GreedyStep]:
    """The source's exact combinatorial search: the global S2 maximum.

    Enumerates all ``C(N, n)`` subsets (capped); the trace has a single step
    carrying the winning subset and its S2.
    """
    norb = state.norb
    if not 1 <= n_target <= norb:
        raise TnassError(
            f"the target size {n_target} is outside 1..{norb} (the window's orbitals)."
        )
    count = math.comb(norb, n_target)
    if count > BRUTE_FORCE_CAP:
        raise TnassError(
            f"the exact search would evaluate C({norb}, {n_target}) = {count} subsets, "
            f"above the cap of {BRUTE_FORCE_CAP}. Next step: use the greedy method (or "
            "a smaller window)."
        )
    best: tuple[int, ...] | None = None
    best_value = -math.inf
    for combo in itertools.combinations(range(norb), n_target):
        value = renyi2_entropy(state, combo)
        if value > best_value + 1e-12:
            best, best_value = combo, value
    assert best is not None
    return [GreedyStep(added=-1, subset=best, s2=best_value)]


# --- the workflow -------------------------------------------------------------


def analyze(
    dump: Fcidump,
    *,
    n_target: int,
    method: str = "greedy",
    block_size: int | None = None,
    reference_energy: float | None = None,
) -> TnassResult:
    """One TNASS selection: exact S2 subsets over the window.

    ``method``: ``"greedy"`` (Eq. 11), ``"block"`` (with ``block_size``), or
    ``"brute"`` (the exact combinatorial search, capped).
    """
    state = _solve(dump, reference_energy)
    singles = tuple(renyi2_entropy(state, (index,)) for index in range(dump.norb))
    if method == "greedy" and block_size is None:
        label = "greedy"
        steps = greedy_selection(state, n_target)
    elif method == "block":
        if block_size is None:
            raise TnassError("the block method needs a block size.")
        label = f"block greedy (k={block_size})"
        steps = block_greedy_selection(state, n_target, block_size)
    elif method == "brute" and block_size is None:
        label = "brute force (the global maximum)"
        steps = brute_selection(state, n_target)
    else:
        raise TnassError(
            f"unknown method {method!r}; use 'greedy', 'block' (with a block size) or "
            "'brute'."
        )
    selected = steps[-1].subset if steps else ()
    densities = entropy_rdm.spin_densities(state)
    occupations = [
        float(densities.gamma_a[i, i] + densities.gamma_b[i, i]) for i in range(dump.norb)
    ]
    n_electrons = int(round(sum(occupations[i] for i in selected)))
    return TnassResult(
        method=label,
        n_target=n_target,
        block_size=block_size,
        steps=tuple(steps),
        single_s2=singles,
        selected=selected,
        n_electrons=n_electrons,
        energy_fci=float(state.energy_total),
    )


# --- output -------------------------------------------------------------------


def render(result: TnassResult) -> str:
    """The selection trace and the resulting space."""
    lines = [
        "TNASS active-space selection (Renyi-2 entropy of the bipartition):",
        f"  method: {result.method}; target size {result.n_target} of the window's "
        "spatial orbitals",
        "",
        f"  {'orbital':>7}  {'S2(single)':>10}  (the seed ranking)",
    ]
    for index, value in enumerate(result.single_s2):
        lines.append(f"  {index:>7}  {value:>10.5f}")
    lines += ["", f"  {'step':>4}  {'added':>5}  {'subset':<24}  {'S2(A)':>8}"]
    for step_number, step in enumerate(result.steps, start=1):
        subset = ", ".join(str(i) for i in step.subset)
        added = f"{step.added:>5}" if step.added >= 0 else f"{'global':>5}"
        lines.append(f"  {step_number:>4}  {added}  {subset:<24}  {step.s2:>8.5f}")
    lines += [
        "",
        f"  selected: orbitals {list(result.selected)} "
        f"({result.n_electrons}e, {len(result.selected)}o)",
        f"  exact window FCI energy: {result.energy_fci:.9f} Eh",
        "",
        "Boundaries and checks:",
    ]
    for item in (
        "the S2 oracle is exact over the FCIDUMP window (the determinant CI of the "
        "four-state-entropy route); the source builds it as a tensor-network "
        "entanglement feature over the full orbital space with an approximate MPS "
        "(its useful bond dimensions are 4-6) -- there is no bond-dimension "
        "truncation here, but the window boundary is the restriction instead",
        "the selection maximizes the bipartition entanglement with the window "
        "complement; the source's best-k variant (choose k by the CASCI energy) is "
        "not implemented -- run the menu at several target sizes and compare",
        "the subset is a spatial-orbital set (both spins travel together), matching "
        "the source's selection domain",
        "delivering the space to a CASSCF needs the orbital-order machinery "
        "(menu 22's boundary note)",
    ):
        lines.append(f"  - {item}")
    return "\n".join(lines)


def run(dump: Fcidump, **kwargs) -> ReportSection:
    """The analyser entry point for menu 26."""
    return ReportSection(
        title="TNASS selection (Renyi-2 bipartition)",
        body=render(analyze(dump, **kwargs)),
    )


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the method and the exact-oracle variant."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "TNASS: select the spatial-orbital subset A maximizing the Renyi-2 "
                "entropy S2(A) = -log Tr(rho_A^2) of the A-versus-complement "
                "bipartition; the entanglement feature |EF> = sum_b e^{-S2(b)}|b> "
                "carries every partition's purity (built from a low-bond-dimension DMRG "
                "state as an MPS of bond dimension chi^4); the selection algorithms are "
                "the exact brute force over C(N, n), the greedy method of Eq. (11), and "
                "the block greedy family (k = 1: the single-orbital ranking; k = n: "
                "greedy).  The source reports lower ground-state energies and more "
                "accurate dipoles than single-orbital-entropy or HOMO/LUMO-window "
                "selection at equal active-space sizes."
            ),
            ref=(
                "arXiv:2608.03645v1, sections 2.3 (Eqs. 1-2), 3.1-3.3 "
                "(the EF construction and the selection algorithms)"
            ),
            url="https://doi.org/10.48550/arXiv.2608.03645",
            bibkey="mingare2026tnass",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "The exact S2 oracle and its anchors: the one-orbital S2 reproduces "
                "-log(sum w^2) from the four-state weights to 1e-12; the two-orbital S2 "
                "reproduces the purity assembled independently from the spin-resolved "
                "2-RDM; the full-window subset has S2 = 0; and the greedy four-orbital "
                "selection on the benzene pi platform lands on the strongly entangled "
                "quartet (values pinned in the tests)."
            ),
            ref="tests/test_tnass.py; fixtures/orca/benzene.fcidump",
        ),
    )
