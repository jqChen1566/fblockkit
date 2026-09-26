"""Orbital-space comparison: the subspace fraction ``sigma_F`` and the space-change SVD.

Two published, zero-external-reference diagnostics share one computation -- the
singular values of the overlap matrix between two orbital sets, both
orthonormalised internally first (the DMET source states the result is
independent of the orthonormalisation scheme):

- ``sigma_F`` (Guan & Jiang, state-averaged DMET, Eq. 10): the RMS singular
  value of ``M = C_A^T S C_B``, normalised by ``sqrt(min(|A|, |B|))``.  It
  answers "how much of space B is contained in space A"; the source reports
  that the ordering of the deficit ``1 - sigma_F`` tracks the ordering of the
  energy error, so its intended use is *ranking* candidate active spaces
  against a fuller reference space.  A large ``sigma_F`` cannot exclude that A
  is too large -- the check is one-sided by construction;
- ``S_change`` (Sayfutyarova, Sun, Chan & Knizia, Eq. 16): for two sets of the
  same size (the initial and final active space of an orbital optimisation),
  singular values near 1 mean the space did not change, and a singular value
  near 0 means one initial active orbital was replaced by an unrelated one.
  The source's reading: the initial space lacked an element needed to represent
  the correlated wave function, so it must be changed or enlarged.

Both diagnosed quantities need only two coefficient blocks and the AO overlap
matrix of one basis -- exactly what two ``orca_2json`` exports of the same gbw
(or of two gbw files of the same system and basis) provide.  The absolute
reading bands used below are provisional; they are anchored on the sources'
own measured numbers (recorded in the report and in ``evidence``), and the
ranking use -- comparing several candidate spaces for one system -- does not
depend on them.
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

__all__ = ["SpaceComparison", "compare_spaces", "jaccard", "run", "evidence"]


class OrbitalSpaceError(ValueError):
    """The comparison cannot run on the given coefficients (with a next step)."""


#: |C^T S C - 1| tolerance for treating a coefficient block as orthonormal.
ORTHONORMAL_TOLERANCE = 1e-6

#: Provisional reading bands (see the module docstring for their anchors).
DEFICIT_SMALL = 1e-4
DEFICIT_LARGE = 1e-2
CHANGE_ESSENTIALLY_UNCHANGED = 0.9
CHANGE_REPLACED = 0.5


def jaccard(window_a, window_b) -> float:
    """The Jaccard index of two orbital windows (|A and B| / |A or B|).

    The evaluation protocol of the RLEASE source reports it next to the energy
    error to locate where a disagreement comes from: a low Jaccard index means
    the two calculations picked different orbitals, a high one with different
    energies means the same space behaved differently.
    """
    set_a = {int(index) for index in window_a}
    set_b = {int(index) for index in window_b}
    union = set_a | set_b
    if not union:
        raise OrbitalSpaceError(
            "both windows are empty, so the Jaccard index is undefined. Next step: give "
            "at least one orbital in a window."
        )
    return len(set_a & set_b) / len(union)


@dataclass(frozen=True)
class SpaceComparison:
    """Singular values of ``C_A^T S C_B`` and the two readings built from them.

    ``singular_values`` is descending, of length ``min(n_a, n_b)``.  Both
    coefficient blocks were orthonormalised in the AO metric before the overlap
    was formed; ``reorthonormalised_a`` / ``_b`` record whether that was
    necessary (``False`` = the block already satisfied ``C^T S C = 1``).
    """

    n_a: int
    n_b: int
    singular_values: tuple[float, ...]
    fraction: float
    deficit: float
    min_singular: float
    reorthonormalised_a: bool
    reorthonormalised_b: bool
    note: str


def _as_block(matrix: np.ndarray, overlap: np.ndarray, name: str) -> np.ndarray:
    block = np.asarray(matrix, dtype=np.float64)
    if block.ndim != 2:
        raise OrbitalSpaceError(
            f"{name} has {block.ndim} dimensions, not 2. Next step: pass an AO "
            "coefficient block with the orbitals in the columns."
        )
    n_ao = overlap.shape[0]
    if block.shape[0] != n_ao:
        raise OrbitalSpaceError(
            f"{name} has {block.shape[0]} AO rows but the overlap matrix is {n_ao} x "
            f"{n_ao}. Next step: the two exports and the overlap matrix must come from "
            "the same basis set."
        )
    if block.shape[1] == 0:
        raise OrbitalSpaceError(
            f"{name} has no orbitals. Next step: give a window with at least one orbital."
        )
    return block


def _orthonormalise(block: np.ndarray, overlap: np.ndarray) -> tuple[np.ndarray, bool]:
    """Return the block made orthonormal in the AO metric, and whether that was needed."""
    gram = block.T @ overlap @ block
    deviation = float(np.abs(gram - np.eye(gram.shape[0])).max())
    if deviation <= ORTHONORMAL_TOLERANCE:
        return block, False
    values, vectors = np.linalg.eigh(gram)
    if values.min() <= 1e-12:
        raise OrbitalSpaceError(
            f"the coefficient block is linearly dependent (smallest metric eigenvalue "
            f"{values.min():.3e}). Next step: remove the redundant orbital from the "
            "window -- the subspace comparison is defined on linearly independent sets."
        )
    inverse_sqrt = vectors @ np.diag(values**-0.5) @ vectors.T
    return block @ inverse_sqrt, True


def compare_spaces(
    block_a: np.ndarray,
    block_b: np.ndarray,
    overlap: np.ndarray,
    *,
    note: str = "",
) -> SpaceComparison:
    """Compute the singular values of ``C_A^T S C_B`` and both readings.

    ``block_a`` / ``block_b`` are AO coefficient matrices with the orbitals in
    the columns; ``overlap`` is the AO overlap matrix of the same basis.
    """
    s = np.asarray(overlap, dtype=np.float64)
    if s.ndim != 2 or s.shape[0] != s.shape[1]:
        raise OrbitalSpaceError(
            f"the overlap matrix is {s.shape}, not square. Next step: export the "
            "S-Matrix together with the MO coefficients (same orca_2json run)."
        )
    a = _as_block(block_a, s, "space A")
    b = _as_block(block_b, s, "space B")
    a, reorth_a = _orthonormalise(a, s)
    b, reorth_b = _orthonormalise(b, s)
    overlap_matrix = a.T @ s @ b
    singular = np.linalg.svd(overlap_matrix, compute_uv=False)
    # Numerical hygiene: singular values of a product of orthonormal blocks are
    # bounded by 1; clip the round-off above it so sqrt() below stays real.
    singular = np.clip(singular, 0.0, 1.0)
    n_min = min(a.shape[1], b.shape[1])
    fraction = float(math.sqrt(float((singular**2).sum()) / n_min))
    return SpaceComparison(
        n_a=int(a.shape[1]),
        n_b=int(b.shape[1]),
        singular_values=tuple(float(value) for value in singular),
        fraction=fraction,
        deficit=1.0 - fraction,
        min_singular=float(singular.min()),
        reorthonormalised_a=reorth_a,
        reorthonormalised_b=reorth_b,
        note=note,
    )


def _definite_note(comparison: SpaceComparison) -> str:
    if comparison.deficit <= DEFICIT_SMALL:
        return "the smaller space is contained in the larger to numerical precision"
    if comparison.deficit <= DEFICIT_LARGE:
        return (
            "the smaller space is contained only approximately: the deficit is small "
            "but well above the numerical level"
        )
    return (
        "the smaller space misses a substantial part of the larger one; the listed "
        "smallest singular values name the directions that are missing"
    )


def run(
    comparison: SpaceComparison,
    *,
    label_a: str = "A",
    label_b: str = "B",
    jaccard_index: float | None = None,
) -> ReportSection:
    """Render the comparison report (both readings, with their provisional bands)."""
    singulars = " ".join(f"{value:.4f}" for value in comparison.singular_values)
    lines = [
        f"Space {label_a}: {comparison.n_a} orbitals; space {label_b}: "
        f"{comparison.n_b} orbitals.",
    ]
    if comparison.reorthonormalised_a or comparison.reorthonormalised_b:
        which = ", ".join(
            name
            for name, flag in (
                (label_a, comparison.reorthonormalised_a),
                (label_b, comparison.reorthonormalised_b),
            )
            if flag
        )
        lines.append(
            f"  ({which} was not orthonormal in the AO metric and was orthonormalised "
            "with its own overlap first; the source states the result is independent "
            "of the orthonormalisation scheme)"
        )
    lines += [
        "",
        f"Singular values of C_{label_a}^T S C_{label_b} (descending, "
        f"{len(comparison.singular_values)} values):",
        f"  {singulars}",
        f"  sigma_F = ||M||_F / sqrt(min(n_A, n_B)) = {comparison.fraction:.6f}",
        f"  deficit 1 - sigma_F = {comparison.deficit:.3e}; smallest singular value = "
        f"{comparison.min_singular:.6f}",
        "",
        "Reading (provisional bands; the sigma_F source -- Guan & Jiang, Eq. (10), "
        "arXiv:2607.08178 as recorded in the project notes, not re-verified -- and the S_change "
        "source report 1 - sigma_F = "
        "2.9e-3/6.8e-5 and a 0.65-0.99 SVD range respectively, and the stated use is "
        "ranking candidate spaces, not an absolute pass mark):",
        f"  - containment: {_definite_note(comparison)} "
        f"(deficit vs the 1e-4 / 1e-2 provisional lines)",
    ]
    if comparison.n_a == comparison.n_b:
        if comparison.min_singular >= CHANGE_ESSENTIALLY_UNCHANGED:
            verdict = "the two spaces are essentially the same space"
        elif comparison.min_singular >= CHANGE_REPLACED:
            verdict = (
                "the spaces are related by a moderate rotation; no initial active "
                "orbital was replaced outright"
            )
        else:
            verdict = (
                "an initial active orbital was largely replaced during the "
                "optimisation -- the initial space lacked an element (AVAS paper, "
                "Eq. 16 and its reading)"
            )
        lines.append(
            f"  - space change (equal sizes; the source's measured range for its "
            f"own runs is 0.65-0.99): {verdict} (smallest singular value vs the "
            f"0.9 / 0.5 provisional lines)"
        )
    else:
        lines.append(
            "  - space change: not read here -- the source's criterion compares two "
            "sets of the same size (initial vs final active space)"
        )
    if jaccard_index is not None:
        lines.append(
            f"  - Jaccard index of the two windows = {jaccard_index:.4f}: reported next "
            "to the overlap so that a disagreement can be split into 'different "
            "orbitals chosen' (low index) against 'same space, different behaviour' "
            "(high index) -- the evaluation protocol of the descriptor-panel source"
        )
    if comparison.note:
        lines += ["", comparison.note]
    return ReportSection(
        title="Orbital-space comparison (sigma_F / space-change SVD)",
        body="\n".join(lines),
    )


def evidence() -> tuple[Evidence, ...]:
    """Provenance of both diagnostics (rule/evidence discipline of the project)."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "The subspace fraction sigma_F = ||M||_F / sqrt(min(|A|,|B|)) with "
                "M_ij = <psi_i^A|psi_j^B> measures how much of one orbital space is "
                "contained in another; each set is orthonormalised first and the "
                "result is independent of the orthonormalisation scheme. The source "
                "reports that the ordering of 1 - sigma_F tracks the ordering of the "
                "energy error (its own values: 2.9e-3 for a loose starting space, "
                "6.8e-5 after a localisation)."
            ),
            ref=(
                "Guan Z.-B., Jiang H., 'State-Averaged Density Matrix Embedding Theory "
                "for Local Excitations', Eq. (10); arXiv:2607.08178 (identifier as "
                "recorded in the project's reading notes; not re-verified against arXiv)"
            ),
        ),
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "The space-change criterion S_change = (C_act^final)^T S C_act^initial: "
                "singular values near 1 mean the active space hardly changed, and a "
                "singular value near 0 means one initial active orbital was replaced by "
                "an unrelated orbital -- the initial space lacked an element needed for "
                "the strong correlation, so it must be modified or enlarged."
            ),
            ref=(
                "Sayfutyarova E. R., Sun Q., Chan G. K.-L., Knizia G., J. Chem. Theory "
                "Comput., 2017, 13, 4063-4078, DOI 10.1021/acs.jctc.7b00128 (Eq. 16; "
                "measured range 0.65-0.99)"
            ),
            bibkey="sayfutyarova2017automated",
            url="https://doi.org/10.1021/acs.jctc.7b00128",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "Identity limit on real data (fixtures n2_fcidump.canonical.json and "
                "n2_fcidump.localized.json, the same six active orbitals of the same "
                "CAS(6,6) before and after an orca_loc IAO-IBO localisation): the six "
                "singular values are 1.0000 to better than 1e-10, so sigma_F = 1 and "
                "the space-change reading is 'unchanged' -- the check is blind to an "
                "orbital rotation inside the same space, as it must be."
            ),
            ref="tests/test_orbital_space.py; fixtures/orca/README.md",
        ),
    )
