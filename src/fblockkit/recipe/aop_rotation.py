"""G7: the AOP rotation guess -- a two-step SVD rotation onto a reference active space.

A MCSCF calculation starts wherever its initial guess sits, and the wrong
basin means the wrong active space.  The source (Paz, Baleeva & Glover 2021)
rotates the current orbital set onto a *reference* active space by the SVD of
their mutual overlap and uses the rotated orbitals as the initial guess.  The
rotation runs in two steps so that the closed and the virtual subspaces do
not mix into the active one (its Eqs. (4)-(11)):

1. the lowest ``n_closed + n_active`` orbitals of the target, ``C_1``
   (Eq. (6)), are rotated by the right singular vectors of
   ``S_MO1 = C_ref,act S C_1^T`` (Eq. (4)); the first ``n_active`` rows of
   the rotated block are the closed-active combinations that align best with
   the reference active orbitals;
2. ``C_2`` (Eq. (7)) concatenates those first ``n_active`` rows with the
   untouched virtual block, and a second SVD of
   ``S_MO2 = C_ref,act S C_2^T`` (Eq. (5)) rotates ``C_2`` so that its first
   ``n_active`` rows align best;
3. the final set is ``C' = C'_clsd + C'_act + C'_virt`` (Eqs. (8)-(11)): the
   closed block from step 1, the active and virtual blocks from step 2.  The
   result is a unitary rotation of the target's own orbital set -- the module
   verifies that in the target's overlap metric and reports the residual.

The alignment is read through ``O_min``, the smallest singular value of the
overlap SVD -- the source's diagnostic (its Eq. (3)): near 1 every aligned
direction is good, a small value flags at least one poor one.  The source's
working criterion ``O_min >= 0.85`` (calibrated on the bimodal distribution
of its 3000-configuration dataset, its Fig. S6) is quoted for its setting:
there it is evaluated on *converged* orbitals, here on the guess.

What the module measures and refuses
------------------------------------

- The reference and the target must share the basis set and the atom order
  (their coefficient matrices must live in the same AO frame); the overlap
  block uses the *target's* S matrix -- exact when the shared atoms sit at
  the same geometry (the source's own application), the source's small-step
  approximation otherwise (the same boundary menu 17's overlap block
  carries); distant geometries are a demonstration, not a calibrated reading.
- **Containment gate**: the reference active space must be represented in
  the target's lowest ``n_closed + n_active`` orbitals at least to the
  source's own alignment line (0.85).  Below it the two-step construction
  cannot keep a clean closed block -- the leftover target directions it
  would assign to the closed manifold are then not the closed orbitals, and
  the resulting guess can send the CASSCF to a wrong solution (measured;
  see the module's evidence).  The gate is this project's; the source's
  settings have containment by construction.
- The active list comes from the reference (0-based indices of its active
  orbitals); the target's partition (``n_closed`` closed orbitals, the same
  ``n_active``, the rest virtual) is declared by the caller and checked
  against the target's orbital count.
- The written file comes from the target's own orbital set: the mkl template
  is the target's, and the occupation/energy tags written with the guess are
  the target's (the rotation preserves the set, not the individual
  occupancies).

The division of labour matches the WASP route's: this module builds the
guess, ORCA re-optimises the orbitals.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..knowledge.models import EVIDENCE_LITERATURE, EVIDENCE_MEASURED, Evidence
from ..parsers.mkl import MklError, MklFile
from .guess_transfer import write_guess_mkl as _write_guess_mkl

__all__ = [
    "AopError",
    "AopResult",
    "rotate_guess",
    "write_guess_mkl",
    "render",
    "evidence",
]


class AopError(ValueError):
    """The rotation cannot be built from the given inputs (with a next step)."""


#: Gate on the reference's representation in the target's closed+active window
#: (the smallest singular value of the first SVD).  Below it the two-step
#: construction cannot keep a clean closed block: the leftover target
#: directions it would assign to the closed manifold are then not the closed
#: orbitals.  Set at the source's own alignment line (0.85), which is there
#: calibrated on *achieved* alignments; here it gates the *containment* step.
#: Measured: a same-geometry reference reads 0.958 and proceeds; a
#: non-corresponding cross-geometry reference reads 0.022 and sent the CASSCF
#: to a wrong solution before this gate existed.
CONTAINMENT_FLOOR = 0.85


#: Singular values within this window count as degenerate: the right singular
#: vectors of that cluster are re-based on the target's own orbitals, so the
#: written file does not inherit the BLAS build's arbitrary choice inside a
#: degenerate singular subspace (measured: without this, the N2 pi pair swapped
#: between machines).
SINGULAR_WINDOW = 1e-9


def _stable_right_singular(matrix):
    """The right singular vectors, degenerate clusters re-based deterministically.

    The SVD's right singular vectors are unique only up to a rotation inside a
    cluster of degenerate singular values, and that choice differs between
    LAPACK builds.  The cluster's *span* is stable, so each degenerate cluster
    is re-expressed on the target's own orbitals: the cluster's projector is
    applied to the unit vectors of the coefficient space in index order and
    Gram-Schmidt-orthonormalised (twice; Euclidean metric -- the right vectors
    live in the Euclidean coefficient space); the signs follow the 2^-j
    weighted sum.  A cluster whose projected basis falls short refuses, because
    the rotation would then depend on the diagonaliser's own choice.
    """
    _, singular, vh = np.linalg.svd(matrix)
    # the trailing rows beyond the singular-value count form the (degenerate)
    # null-space block of a non-square matrix -- always re-based
    boundaries = [(len(singular), vh.shape[0])] if vh.shape[0] > len(singular) else []
    start = 0
    for position in range(1, len(singular) + 1):
        if position < len(singular) and singular[position - 1] - singular[position] <= SINGULAR_WINDOW:
            continue
        if position - start > 1:
            boundaries.append((start, position))
        start = position
    for start, position in boundaries:
        if position - start > 1:
            cluster = vh[start:position]  # (m, n), orthonormal rows
            projector = cluster.T @ cluster
            basis = []
            for unit in np.eye(projector.shape[0]):
                vector = projector @ unit
                for _ in range(2):
                    for previous in basis:
                        vector = vector - float(previous @ vector) * previous
                norm = float(np.sqrt(max(float(vector @ vector), 0.0)))
                if norm > 1e-6:
                    basis.append(vector / norm)
                    if len(basis) == cluster.shape[0]:
                        break
            if len(basis) != cluster.shape[0]:
                raise AopError(
                    "a degenerate singular cluster could not be re-based on the "
                    "target's orbitals (the projected basis fell short); the "
                    "rotation would depend on the linear-algebra build's own "
                    "choice there. Next step: report this input with its exports."
                )
            fixed = np.array(basis)
            weights = 2.0 ** -np.arange(fixed.shape[1])
            for row in range(fixed.shape[0]):
                if float(fixed[row] @ weights) < 0.0:
                    fixed[row] = -fixed[row]
            vh[start:position] = fixed
    # every row's sign: the 2^-j weighted sum decides (a singular vector is
    # unique up to its sign, and the sign of a noise-crossing component can
    # flip between builds)
    weights = 2.0 ** -np.arange(vh.shape[1])
    for row in range(vh.shape[0]):
        if float(vh[row] @ weights) < 0.0:
            vh[row] = -vh[row]
    return singular, vh


@dataclass(frozen=True)
class AopResult:
    """The rotated guess and everything the report needs to say about it."""

    coefficients: tuple[tuple[float, ...], ...]  # n_ao x n_mo, export row order
    reference: str
    target: str
    n_closed: int
    n_active: int
    n_virtual: int
    containment: float  # the reference's representation in the target window (>= floor)
    omin: float  # the built active block's smallest singular value vs the reference
    residual: float  # max |C' S C'^T - I| of the written set
    occupations: tuple[float, ...]
    energies: tuple[float, ...]
    checks: tuple[str, ...]
    notes: tuple[str, ...]


def rotate_guess(reference, target, active, *, n_closed: int) -> AopResult:
    """Rotate the target's orbitals onto the reference active space (two-step SVD).

    ``reference``/``target``: ``orca_2json`` exports of the same system and
    basis (same atoms, same atom order); ``active``: the 0-based indices of
    the reference's active orbitals; ``n_closed``: the target's number of
    closed-shell orbitals (the active block follows, then the virtuals).
    """
    if target.overlap is None:
        raise AopError(
            f"the target export '{target.base_name}' carries no S-Matrix, so the "
            "rotation has no metric. Next step: re-export with the S-Matrix (the "
            "orca_2json default plus the 1elIntegrals configuration)."
        )
    if not reference.mo_coefficients or not target.mo_coefficients:
        raise AopError(
            "one of the exports carries no MO coefficients, so there is no orbital "
            "set to rotate. Next step: re-export both files with "
            '"MOCoefficients": true.'
        )
    if tuple(reference.atoms) != tuple(target.atoms) or reference.n_ao != target.n_ao:
        raise AopError(
            f"the reference '{reference.base_name}' and the target '{target.base_name}' "
            "do not share the basis set and atom order (their coefficient matrices "
            "must live in the same AO frame). Next step: use a reference of the same "
            "system in the same basis and atom order."
        )
    indices = [int(index) for index in active]
    if not indices:
        raise AopError(
            "the reference active list is empty. Next step: give the 0-based indices "
            "of the reference's active orbitals."
        )
    c_ref = np.array(reference.mo_coefficients, dtype=float)
    c_target = np.array(target.mo_coefficients, dtype=float)
    overlap = np.array(target.overlap, dtype=float)
    for index in indices:
        if not 0 <= index < c_ref.shape[0]:
            raise AopError(
                f"the reference active index {index} is outside the export "
                f"'{reference.base_name}'s {c_ref.shape[0]} orbitals. Next step: check "
                "the active list."
            )
    n_active = len(indices)
    n_mo = c_target.shape[0]
    if not 0 <= n_closed < n_mo:
        raise AopError(
            f"the closed count {n_closed} is outside the target's {n_mo} orbitals. "
            "Next step: give the target's number of closed-shell orbitals."
        )
    if n_closed + n_active > n_mo:
        raise AopError(
            f"closed + active = {n_closed} + {n_active} = {n_closed + n_active} exceeds "
            f"the target's {n_mo} orbitals. Next step: check the closed count and the "
            "reference active list (the partition must fit the target's orbital count)."
        )
    c_ref_active = c_ref[indices, :]
    # step 1: rotate the lowest n_closed + n_active target orbitals (Eqs. (4), (6))
    c1 = c_target[: n_closed + n_active, :]
    singular1, vh1 = _stable_right_singular(c_ref_active @ overlap @ c1.T)
    containment = float(singular1.min())
    if containment < CONTAINMENT_FLOOR:
        raise AopError(
            "the reference active space is not sufficiently represented in the "
            f"target's lowest {n_closed + n_active} orbitals (containment "
            f"{containment:.3f} < {CONTAINMENT_FLOOR:g}), so the two-step rotation "
            "cannot keep a clean closed block: the leftover target directions it "
            "would assign to the closed manifold are then not the closed orbitals. "
            "Next step: use a reference whose active space corresponds to the "
            "target's closed+active window (e.g. CASSCF orbitals pushed to a "
            "close-by geometry), or extend the window if the target partition "
            "allows it."
        )
    c1_rotated = vh1 @ c1
    # step 2: the rotated active block + the untouched virtuals (Eqs. (5), (7))
    c2 = np.vstack((c1_rotated[:n_active, :], c_target[n_closed + n_active :, :]))
    _, vh2 = _stable_right_singular(c_ref_active @ overlap @ c2.T)
    c2_rotated = vh2 @ c2
    # assembly (Eqs. (8)-(11)): closed from step 1, active and virtual from step 2
    c_prime = np.vstack(
        (c1_rotated[n_active:, :], c2_rotated[:n_active, :], c2_rotated[n_active:, :])
    )
    achieved = c_ref_active @ overlap @ c_prime[n_closed : n_closed + n_active, :].T
    omin = float(np.linalg.svd(achieved, compute_uv=False).min())
    residual = float(np.abs(c_prime @ overlap @ c_prime.T - np.eye(n_mo)).max())
    n_virtual = n_mo - n_closed - n_active
    checks = (
        "the two-step rotation follows the source's Eqs. (4)-(11): the first SVD "
        "rotates the lowest closed+active target orbitals (Eq. (6)), the second "
        "concatenates the rotated active block with the untouched virtuals (Eq. (7)) "
        "and rotates again; the final set takes the closed block from the first step "
        "and the active and virtual blocks from the second (Eqs. (8)-(11))",
        "the overlap block uses the *target's* S matrix -- the source's own setting "
        "(chromophore orbitals inside a same-geometry condensed-phase basis); between "
        "distant geometries it is the small-step approximation the source uses along "
        "its interpolations, so the alignment there is a demonstration, not a "
        "calibrated reading",
        f"the reference's representation in the target's closed+active window "
        f"(containment) is {containment:.6f}: the smallest singular value of the "
        "reference-to-window overlap; the gate requires the source's own alignment "
        f"line ({CONTAINMENT_FLOOR:g}), below which the closed block the construction "
        "would build is not the closed manifold (measured: a non-corresponding "
        "cross-geometry reference reads 0.022 and sent the CASSCF to a wrong "
        "solution)",
        f"O_min = {omin:.6f} is the smallest singular value of the built active "
        "block's overlap with the reference -- the source's alignment diagnostic (its "
        "Eq. (3)); its working criterion O_min >= 0.85 is calibrated on its own "
        "condensed-phase dataset and read there off *converged* orbitals, while this "
        "value belongs to the guess (measured: with a same-geometry reference the "
        "built block reproduces the reference exactly, O_min = 1.000000)",
        f"the guess is a unitary rotation of the target's own orbital set; the "
        f"orthonormality residual max |C' S C'^T - I| in the target's metric is "
        f"{residual:.2e}, and the occupation/energy tags written with it are the "
        "target's",
    )
    notes = (
        "Within a degenerate reference window the individual orbital assignment is "
        "arbitrary (any rotation among degenerate reference orbitals aligns equally "
        "well); the SVD fixes the subspace, which is what the O_min diagnostic reads",
        "The written file is a Molekel mkl: run ``orca_2mkl <name>.fbk -gbw`` to turn "
        "it into a ``.gbw`` and read it with ``!moread`` + ``%moinp`` (the route "
        "measured for menu 18).",
        "This module builds the guess; ORCA re-optimises the orbitals. The rotated "
        "set is not a converged orbital set and is not claimed to be one.",
    )
    return AopResult(
        coefficients=tuple(tuple(float(value) for value in row) for row in c_prime.T),
        reference=reference.base_name,
        target=target.base_name,
        n_closed=int(n_closed),
        n_active=int(n_active),
        n_virtual=int(n_virtual),
        containment=containment,
        omin=omin,
        residual=residual,
        occupations=tuple(float(value) for value in target.mo_occupations),
        energies=tuple(float(value) for value in target.mo_energies),
        checks=checks,
        notes=notes,
    )


def write_guess_mkl(template_mkl: MklFile, result: AopResult, path) -> None:
    """Write the rotated guess into a copy of the target's mkl (the shared route)."""
    try:
        _write_guess_mkl(template_mkl, result, path)
    except MklError as exc:  # the template and the coefficients disagree about the basis
        raise AopError(
            f"the target mkl does not accept the rotated coefficients: {exc}"
        ) from exc


def render(result: AopResult) -> str:
    """The partition, the diagnostics and the boundaries, as the menu prints them."""
    lines = [
        "AOP rotation guess (two-step SVD rotation of the target's orbitals onto a "
        "reference active space):",
        f"  reference: {result.reference}",
        f"  target:    {result.target}",
        f"  partition at the target: closed {result.n_closed}, active "
        f"{result.n_active}, virtual {result.n_virtual}",
        "",
        f"  reference containment in the window: {result.containment:.6f}",
        f"  O_min of the built active block:    {result.omin:.6f}",
        f"  orthonormality residual:            {result.residual:.2e}",
        "",
        "Criteria and boundaries:",
    ]
    for item in result.checks + result.notes:
        lines.append(f"  - {item}")
    return "\n".join(lines)


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the rotation protocol and of the write-back route."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "AOP-MCSCF: rotate the current orbital set onto a reference active "
                "space by the SVD of the MO overlap and use the rotated orbitals as "
                "the MCSCF initial guess.  The two-step procedure (first the lowest "
                "closed+active orbitals, then the rotated active block with the "
                "untouched virtuals) keeps the closed and virtual subspaces from "
                "mixing into the active one (Eqs. (4)-(11)); O_min, the smallest "
                "singular value of the overlap SVD (Eq. (3)), is the alignment "
                "diagnostic, with O_min >= 0.85 as the working criterion calibrated "
                "on the bimodal distribution of the source's 3000-configuration "
                "dataset.  The source reports the method converging over 90% of "
                "3000 thermally sampled condensed-phase configurations where an "
                "aufbau initial guess failed for every one of them, and compares "
                "favourably with the swaps-based alternative (its Table S3)."
            ),
            ref=(
                "Paz A. S. P., Baleeva N. S., Glover W. J., J. Chem. Phys., 2021, 155, "
                "071103, DOI 10.1063/5.0058673 (Eqs. (1)-(11); the O_min criterion in "
                "section III)"
            ),
            url="https://doi.org/10.1063/5.0058673",
            bibkey="paz2021active",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "The ORCA route and the guess quality, measured on 6.1.1 "
                "(N2/def2-SVP, fixtures/orca/n2_cas666_1.600.*): with the reference "
                "the converged CASSCF(6,6) active space at the target's own geometry "
                "and the target the RHF orbital set, the rotation reproduces the "
                "reference active block exactly (O_min 1.000000; containment 0.958; "
                "orthonormality residual 9.8e-15 in the target's metric), and a "
                "CASSCF(6,6) started from the written guess (``!NoIter moread`` + "
                "``%casscf MaxIter``) converges to the reference solution in 6 "
                "macro-iterations where the aufbau start takes 7 (energies "
                "-108.772368705 vs -108.772368734 Eh). The containment gate came out "
                "of a measured failure: a non-corresponding cross-geometry reference "
                "(CASSCF(6,6) at 1.094 against the RHF set at 1.600; containment "
                "0.022) produced a guess whose closed block ORCA flagged as "
                "delocalized, and the CASSCF wandered to a wrong solution (-77.09 Eh "
                "after 30 macro-iterations) -- the gate now refuses that input."
            ),
            ref="tests/test_aop_rotation.py; fixtures/orca/n2_cas666_1.600.*",
        ),
    )
