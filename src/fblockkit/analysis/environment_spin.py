"""Environment spin-polarisation entropy: the Delta S_E criterion of Ai et al. (JCTC 2025).

The quantity (their Eq. 9)

    Delta S_E = -2 Tr[ (D_E/2) ln(D_E/2) ] + Tr[ D_E^a ln D_E^a ] + Tr[ D_E^b ln D_E^b ]

with ``D_E = D_E^alpha + D_E^beta`` the spin-resolved *truncated* (sub-block)
density matrices of the environment vanishing exactly when ``D_E^alpha =
D_E^beta``: it measures spin polarisation that sits in the environment orbitals
rather than on the metal centre.  In the source's workflow it is the R-DIIS
residual, i.e. the criterion for "the spin polarisation is all on the metal";
their 1Dy case separates a wrong solution (2.766) from the correct one (0.007),
and using the wrong one costs a factor 42 in the fitted crystal-field MAE
(357.7 cm^-1 vs 8.6 cm^-1).

Where the numbers come from here: the exact-route density objects
(:mod:`fblockkit.analysis.entropy_rdm`) supply the spin-resolved 1-RDMs in the
localised active basis, and the environment is the set of active orbitals whose
centre is not the cluster centre -- the centre assignment is the largest Löwdin
population over an atom's AO block, the same "largest weight" rule the toolkit
uses elsewhere for composition decisions.

Boundary of the criterion in this route (stated because it decides which
questions the number can answer):

- the inactive orbitals of a CASSCF wave function are doubly occupied by
  construction, so they contribute zero spin polarisation and the environment
  here is the *active* environment; a mean-field solution's ligand-core spin
  polarisation (the source works with ROHF densities) is not visible to this
  route.  For that use the local-spin analysis of an SCF output (menu 1,
  ``analysis.local_spin``), which reports per-fragment <S_z> and <S_A S_B> for
  exactly that data;
- with a pure metal-centred active space the environment is empty and
  Delta S_E = 0 by construction -- reported as such rather than as a pass;
- for a closed-shell (S = 0) state the singlet symmetry forces
  ``gamma^alpha = gamma^beta`` in *every* orbital basis, so the criterion is
  blind by construction (a clean negative control on the N2 fixture).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from ..knowledge.models import (
    EVIDENCE_LITERATURE,
    EVIDENCE_MEASURED,
    Evidence,
    ReportSection,
)

__all__ = [
    "EnvironmentSpinError",
    "analyze",
    "EnvironmentSpinEntropy",
    "environment_spin_entropy",
    "loewdin_atom_populations",
    "partition_by_centre",
    "run",
    "evidence",
]


class EnvironmentSpinError(ValueError):
    """The environment-spin analysis cannot run on the given data (with a next step)."""

#: Default masking threshold for eigenvalues in the logarithms.
DEFAULT_TOLERANCE = 1e-12

#: Provisional reading line: the source's correct solution measured 0.007 and its
#: wrong solution 2.766; the line below sits an order of magnitude above the
#: correct-solution scale and eight orders above the exact route's noise.
PROVISIONAL_LINE = 0.05


@dataclass(frozen=True)
class EnvironmentSpinEntropy:
    """The three entropy terms of Eq. 9 and the resulting Delta S_E."""

    eigenvalues_total: tuple[float, ...]
    eigenvalues_alpha: tuple[float, ...]
    eigenvalues_beta: tuple[float, ...]
    entropy_total: float
    entropy_alpha: float
    entropy_beta: float
    delta_s: float
    electrons_alpha: float
    electrons_beta: float
    masked: int


def _entropy(values: np.ndarray, tolerance: float) -> tuple[float, int]:
    """-sum v ln v over the eigenvalues above the mask, and how many were masked."""
    alive = values > tolerance
    masked = int(values.size - alive.sum())
    positive = values[alive]
    return float(-np.sum(positive * np.log(positive))), masked


def _eigenvalues(matrix: np.ndarray, name: str, tolerance: float) -> np.ndarray:
    values = np.linalg.eigvalsh(np.asarray(matrix, dtype=np.float64))
    if values.min() < -1e-8:
        raise EnvironmentSpinError(
            f"the {name} block has a negative eigenvalue ({values.min():.3e}): a "
            "sub-block of a density matrix cannot do that. Next step: check that the "
            "block and the density objects belong to the same state and basis."
        )
    return np.clip(values, 0.0, None)


def environment_spin_entropy(
    d_alpha: np.ndarray,
    d_beta: np.ndarray,
    *,
    tolerance: float = DEFAULT_TOLERANCE,
) -> EnvironmentSpinEntropy:
    """Eq. 9 on the environment sub-blocks of the two spin densities.

    Eigenvalues at or below ``tolerance`` are masked out of the logarithms (the
    x ln x -> 0 limit; the masking is explicit here for the same reason the
    project masks it in every other log-of-density step).
    """
    a = np.asarray(d_alpha, dtype=np.float64)
    b = np.asarray(d_beta, dtype=np.float64)
    if a.shape != b.shape:
        raise EnvironmentSpinError(
            f"the alpha and beta environment blocks have shapes {a.shape} and {b.shape}. "
            "Next step: both must come from the same orbital set."
        )
    if a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise EnvironmentSpinError(
            f"the environment blocks are {a.shape}, not square. Next step: pass the "
            "sub-block of the density matrices over the environment orbitals."
        )
    if a.shape[0] == 0:
        raise EnvironmentSpinError(
            "the environment block is empty. Next step: this state has no environment "
            "orbitals -- report zero polarisation by construction instead of calling "
            "this function."
        )
    values_a = _eigenvalues(a, "alpha", tolerance)
    values_b = _eigenvalues(b, "beta", tolerance)
    values_total = _eigenvalues((a + b) / 2.0, "spin-summed (halved)", tolerance)
    # -2 Tr[(D/2) ln(D/2)]: with the eigenvalues of D/2 (the halved block)
    entropy_total, _ = _entropy(values_total, tolerance)
    entropy_total *= 2.0
    entropy_alpha, masked_a = _entropy(values_a, tolerance)
    entropy_beta, masked_b = _entropy(values_b, tolerance)
    delta_s = entropy_total - entropy_alpha - entropy_beta
    # the criterion's zero is exact; the exact route leaves round-off of the order of
    # 1e-16 on it, which is reported as a clean zero (measured: -5.6e-17 on the N2
    # singlet, whose value must be 0 by symmetry)
    if abs(delta_s) < 1e-12:
        delta_s = 0.0
    return EnvironmentSpinEntropy(
        eigenvalues_total=tuple(float(v) * 2.0 for v in values_total),
        eigenvalues_alpha=tuple(float(v) for v in values_a),
        eigenvalues_beta=tuple(float(v) for v in values_b),
        entropy_total=entropy_total,
        entropy_alpha=entropy_alpha,
        entropy_beta=entropy_beta,
        delta_s=delta_s,
        electrons_alpha=float(np.trace(a)),
        electrons_beta=float(np.trace(b)),
        masked=masked_a + masked_b,
    )


def loewdin_atom_populations(
    coefficients: np.ndarray, overlap: np.ndarray, labels, n_atoms: int
) -> tuple[tuple[float, ...], ...]:
    """Per-atom Löwdin population of every orbital given (orbitals in the columns).

    ``p_A(i) = sum_{mu in A} |(S^{1/2} C)_mu,i|^2`` -- the Löwdin (symmetrically
    orthogonalised) atomic population, which is non-negative and sums to exactly
    one over the atoms for a normalised orbital.  The naive same-atom block sum
    ``sum_{mu,nu in A} C S C`` must NOT be used: it drops the cross-atom entries
    and does not sum to 1 (measured: 1.056 on the N2 fixture).
    """
    coef = np.asarray(coefficients, dtype=np.float64)
    s = np.asarray(overlap, dtype=np.float64)
    if coef.ndim != 2 or coef.shape[0] != s.shape[0]:
        raise EnvironmentSpinError(
            f"the coefficient block is {coef.shape} against an overlap matrix of "
            f"{s.shape}. Next step: pass the AO coefficient matrix of the same export."
        )
    blocks = [[] for _ in range(n_atoms)]
    for index, label in enumerate(labels):
        if not 0 <= label.center < n_atoms:
            raise EnvironmentSpinError(
                f"an AO label names centre {label.center} but the export lists only "
                f"{n_atoms} atom(s). Next step: check the export."
            )
        blocks[label.center].append(index)
    values, vectors = np.linalg.eigh(s)
    if values.min() <= 1e-10:
        raise EnvironmentSpinError(
            f"the AO overlap matrix is singular (smallest eigenvalue {values.min():.3e}). "
            "Next step: check the basis set of the export."
        )
    root = vectors @ np.diag(values**0.5) @ vectors.T
    symmetrised = root @ coef
    populations = []
    for orbital in range(coef.shape[1]):
        column = symmetrised[:, orbital]
        populations.append(
            tuple(float((column[indices] ** 2).sum()) if indices else 0.0 for indices in blocks)
        )
    return tuple(populations)


def partition_by_centre(
    coefficients: np.ndarray,
    overlap: np.ndarray,
    labels,
    n_atoms: int,
    *,
    cluster_centres: tuple[int, ...],
) -> tuple[tuple[int, ...], tuple[float, ...]]:
    """Assign every orbital to its largest-population centre; return (centres, populations).

    The largest-weight rule (never a ratio or an absolute threshold) is the one
    the toolkit's composition checks converged on; ``cluster_centres`` names the
    centres treated as the cluster, everything else is environment.
    """
    populations = loewdin_atom_populations(coefficients, overlap, labels, n_atoms)
    centres = tuple(
        int(np.argmax(values)) if max(values) > 0 else -1 for values in populations
    )
    return centres, tuple(max(values) for values in populations)


def analyze(
    densities,
    coefficients: np.ndarray,
    overlap: np.ndarray,
    labels,
    n_atoms: int,
    *,
    cluster_centres: tuple[int, ...],
) -> ReportSection:
    """Partition the active orbitals, build the environment blocks and render the section.

    ``densities`` are the spin-resolved 1-RDMs in the basis the ``coefficients``
    belong to (for the localised-environment reading: both rotated to the
    orca_loc basis first).
    """
    centres, populations = partition_by_centre(
        coefficients, overlap, labels, n_atoms, cluster_centres=cluster_centres
    )
    assignment = tuple(
        (index, centre, population)
        for index, (centre, population) in enumerate(zip(centres, populations))
    )
    environment = tuple(
        index for index, centre in enumerate(centres) if centre not in cluster_centres
    )
    if not environment:
        # an empty environment is a statement, not a computation: report the
        # partition and the by-construction zero
        result = None
    else:
        block = np.ix_(environment, environment)
        result = environment_spin_entropy(
            densities.gamma_a[block], densities.gamma_b[block]
        )
    return run(
        result,
        environment_orbitals=environment,
        assignment=assignment,
        cluster_centres=cluster_centres,
    )


def run(
    result: EnvironmentSpinEntropy | None,
    *,
    environment_orbitals: tuple[int, ...],
    assignment: tuple[tuple[int, int, float], ...],
    cluster_centres: tuple[int, ...],
) -> ReportSection:
    """Render the Delta S_E section of the exact-state report.

    ``assignment`` lists (orbital index, assigned centre, population) per active
    orbital; ``cluster_centres`` are the centres treated as the cluster.
    ``result`` may be ``None``: that is the empty-environment case, which is a
    statement ("zero by construction"), not a computation.
    """
    lines = [
        "Partition of the active orbitals (largest Löwdin population per centre; "
        "the centres treated as the cluster: "
        + ", ".join(str(c) for c in cluster_centres)
        + "):",
    ]
    for orbital, centre, population in assignment:
        role = "cluster" if centre in cluster_centres else "environment"
        lines.append(
            f"  orbital {orbital}: centre {centre} (population {population:.4f}) -> {role}"
        )
    lines.append("")
    if not environment_orbitals:
        lines += [
            "Environment: none -- every active orbital is assigned to a cluster centre.",
            "",
            "Reading:",
            "  - the environment is empty in this active space: Delta S_E = 0 by "
            "construction, which says nothing about the metal centre's own polarisation",
        ]
    else:
        assert result is not None
        lines += [
            f"Environment: {len(environment_orbitals)} active orbital(s) "
            f"{list(environment_orbitals)}; electrons in the environment block: "
            f"{result.electrons_alpha:.6f} alpha + {result.electrons_beta:.6f} beta",
            f"  eigenvalues of D_E (alpha+beta): "
            + " ".join(f"{value:.6f}" for value in result.eigenvalues_total),
            f"  Delta S_E = -2 Tr[(D/2) ln(D/2)] + Tr[D_a ln D_a] + Tr[D_b ln D_b]",
            f"            = {result.entropy_total:.6f} - {result.entropy_alpha:.6f} "
            f"- {result.entropy_beta:.6f} = {result.delta_s:.6f}"
            + (
                f"  ({result.masked} eigenvalue(s) masked at the tolerance)"
                if result.masked
                else ""
            ),
            "",
            "Reading (provisional line; the source's anchors are 0.007 for its correct "
            "solution and 2.766 for its wrong one, and the criterion's zero is exact):",
        ]
    if environment_orbitals and result is not None:
        if result.delta_s <= PROVISIONAL_LINE:
            lines.append(
                f"  - the environment carries no appreciable spin polarisation "
                f"(Delta S_E = {result.delta_s:.6f})"
            )
        else:
            lines.append(
                f"  - the environment carries spin polarisation (Delta S_E = "
                f"{result.delta_s:.6f}, above the provisional line): for a single-ion "
                "target this is the source's failure mode -- check whether those "
                "electrons belong there (a ligand-centred open shell) or the solution "
                "is wrong"
            )
    lines += [
        "",
        "Boundary of this route: the inactive orbitals of the CASSCF wave function are "
        "doubly occupied and carry no spin by construction, so this is the *active* "
        "environment. A mean-field solution's ligand spin polarisation is diagnosed by "
        "the local-spin analysis of the SCF output (menu 1) instead.",
    ]
    return ReportSection(title="Environment spin-polarisation entropy (Delta S_E)", body="\n".join(lines))


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the Delta S_E criterion."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "Delta S_E (Eq. 9) is the environment spin-polarisation entropy "
                "-2 Tr[(D_E/2) ln(D_E/2)] + Tr[D_E^a ln D_E^a] + Tr[D_E^b ln D_E^b]; it "
                "vanishes exactly when the environment's alpha and beta densities are "
                "equal, i.e. when all spin polarisation sits on the metal. The source "
                "uses it as the R-DIIS residual and measures 2.766 (wrong solution) vs "
                "0.007 (correct solution) for a 1Dy complex, where the wrong solution "
                "costs a factor 42 in the crystal-field MAE (357.7 vs 8.6 cm^-1)."
            ),
            ref=(
                "Ai et al., J. Chem. Theory Comput., 2025, 21, 9631-9640, "
                "DOI 10.1021/acs.jctc.5c01336 (Eq. 9), code github.com/IrisA144/"
                "liblan_preview"
            ),
            bibkey="ai2025density",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "Exact-route limits (tests/test_environment_spin.py): a closed-shell "
                "state gives Delta S_E = 0 in every orbital basis (the S = 0 spin-flip "
                "symmetry forces gamma_a = gamma_b -- measured on the N2 CAS(6,6) "
                "fixture), one unpaired electron wholly inside the environment block "
                "gives ln 2 to machine precision, and a synthetic two-orbital "
                "environment reproduces the hand-computed entropy of its eigenvalues."
            ),
            ref="tests/test_environment_spin.py; fixtures/orca/README.md",
        ),
    )
