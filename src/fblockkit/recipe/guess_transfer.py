"""G5: the WASP initial guess for a geometry series -- write-back through the mkl.

Keeping one active space consistent along a path (menu 17's mapping) is one
half of the cross-structure problem; the other half is *starting* every
structure's calculation from an orbital set that already lies in the right
basin.  The source protocol (the automatic-active-space review's WASP, its
Eqs. (11)-(13)) builds the guess for a new geometry from the already-computed
neighbours:

    C = sum_beta  w_beta C_beta / (sum_beta w_beta),
    w_beta = 1 / d(G, G_beta),

over the neighbourhood ``d(G, G_beta) <= delta`` (``d`` the RMSD between the
structures), followed by an orthonormalisation of the mixture -- which is what
this module does, and then writes the result into a **gbw-ready Molekel mkl**
(``parsers/mkl.py``): ``orca_2mkl <name>.fbk -gbw`` turns it into a ``.gbw``
that ORCA accepts through ``!moread`` + ``%moinp`` (measured on 6.1.1: the
round trip preserves coefficients to the print precision and a cross-geometry
file is accepted as ``INITIAL GUESS: MOREAD``).

What the module measures and refuses
------------------------------------

- The structures must share the basis set and the atom order (the same check
  the mapping makes); the mixture is built from the exports' AO-major
  coefficient matrices, which is meaningful exactly under that condition.
  Rotated structures are out of scope: the coefficients live in each
  structure's own AO frame, so the source's RMSD metric is used without
  alignment (documented, not silently aligned).
- The orthonormalisation is Loewdin's, in the **target geometry's** overlap
  metric (from the template export); the orthogonalisation residual is
  reported.  A near-singular mixture (a neighbour that duplicates another's
  orbitals) is refused with the residual named.
- When one neighbour sits exactly on the target geometry its weight dominates
  everything (1/d diverges): the mixture degenerates to that neighbour, and
  the report says so.
- The occupations and orbital energies written into the mkl come from the
  heaviest neighbour (the orbitals are a mixture; the tags have to come from
  somewhere and the nearest converged structure is the defensible choice).

The division of labour follows the source: this module builds the guess, ORCA
re-optimises the orbitals; nothing here claims the mixture is a converged
orbital set.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..knowledge.models import EVIDENCE_LITERATURE, EVIDENCE_MEASURED, Evidence
from ..parsers.mkl import MklError, MklFile

__all__ = [
    "GuessError",
    "GuessResult",
    "interpolate_guess",
    "write_guess_mkl",
    "render",
    "evidence",
]

#: Below this distance a neighbour counts as "at the target geometry" and the
#: mixture degenerates to it (the 1/d weight diverges otherwise).
COINCIDENT_TOLERANCE = 1e-6

#: Eigenvalues of the overlap metric below this floor make the orthonormalisation
#: meaningless (the mixture has linearly dependent orbitals).
SINGULAR_FLOOR = 1e-8


class GuessError(ValueError):
    """The guess cannot be built from the given structures (with a next step)."""


@dataclass(frozen=True)
class GuessWeight:
    """One neighbour's distance and normalised weight."""

    name: str
    distance: float  # RMSD in Angstrom
    weight: float  # normalised; 0 when outside the neighbourhood


@dataclass(frozen=True)
class GuessResult:
    """The interpolated guess and everything the report needs to say about it."""

    coefficients: tuple[tuple[float, ...], ...]  # n_ao x n_mo, export row order
    weights: tuple[GuessWeight, ...]
    delta: float | None
    residual: float  # max |c^T S c - I| after the Loewdin orthonormalisation
    source: str  # the neighbour the occupations/energies are taken from
    occupations: tuple[float, ...]
    energies: tuple[float, ...]
    checks: tuple[str, ...]
    notes: tuple[str, ...]


def _matrix(export) -> np.ndarray:
    """The AO-major coefficient matrix of an export (n_ao x n_mo)."""
    return np.array(export.mo_coefficients, dtype=float).T


def _coordinates(export) -> np.ndarray:
    return np.array(export.coordinates, dtype=float)


def rmsd(a, b) -> float:
    """The root-mean-square distance between two structures of the same atom order."""
    difference = _coordinates(a) - _coordinates(b)
    return float(np.sqrt((difference**2).sum(axis=1).mean()))


def interpolate_guess(neighbours, template, *, delta: float | None = None) -> GuessResult:
    """Build the WASP guess for the template geometry from the neighbour exports.

    ``neighbours``: exports of already-computed structures; ``template``: the
    export of the target geometry (its overlap matrix orthonormalises the
    mixture).  ``delta``: the neighbourhood radius in Angstrom (``None``: every
    given neighbour counts).
    """
    if not neighbours:
        raise GuessError(
            "no neighbour structures were given. Next step: give the exports of the "
            "already-computed structures of the series."
        )
    for export in (template, *neighbours):
        if export.overlap is None:
            raise GuessError(
                f"the export '{export.base_name}' carries no S-Matrix, so the mixture "
                "cannot be orthonormalised in the target geometry's metric. Next step: "
                "re-export with the S-Matrix (the orca_2json default plus the "
                "1elIntegrals configuration)."
            )
        if not export.coordinates:
            raise GuessError(
                f"the export '{export.base_name}' carries no coordinates, so the "
                "structure distances are undefined. Next step: re-export with the "
                "coordinates (the orca_2json default)."
            )
    template_atoms = tuple(template.atoms)
    template_nao = template.n_ao
    template_nmo = template.n_mo
    for export in neighbours:
        if tuple(export.atoms) != template_atoms:
            raise GuessError(
                f"the neighbour '{export.base_name}' has a different atom list than the "
                "template, so the coefficient matrices are not comparable. Next step: "
                "interpolate between structures of the same system (same atoms, same "
                "order, same basis)."
            )
        if export.n_ao != template_nao or export.n_mo != template_nmo:
            raise GuessError(
                f"the neighbour '{export.base_name}' has {export.n_ao} AOs x "
                f"{export.n_mo} MOs while the template has {template_nao} x "
                f"{template_nmo}. Next step: interpolate within one and the same basis "
                "set and electron count."
            )
    distances = [(export, rmsd(export, template)) for export in neighbours]
    if delta is not None:
        kept = [(export, distance) for export, distance in distances if distance <= delta]
        if not kept:
            raise GuessError(
                f"no neighbour lies within delta = {delta:g} Angstrom of the template "
                f"(closest: {min(distance for _, distance in distances):g}). Next step: "
                "increase delta or give neighbours of the target geometry."
            )
    else:
        kept = distances
    nearest = min(kept, key=lambda item: item[1])
    if nearest[1] < COINCIDENT_TOLERANCE:
        weights = [
            GuessWeight(
                name=export.base_name,
                distance=distance,
                weight=1.0 if export is nearest[0] else 0.0,
            )
            for export, distance in distances
        ]
        mixture = _matrix(nearest[0])
        notes = (
            "A neighbour coincides with the template geometry (distance below "
            f"{COINCIDENT_TOLERANCE:g} Angstrom), so its weight dominates the 1/d "
            "mixture and the guess is that structure's orbital set as it stands.",
        )
    else:
        raw = np.array([1.0 / distance for _, distance in kept])
        raw /= raw.sum()
        weight_by_id = {id(export): float(value) for (export, _), value in zip(kept, raw)}
        weights = [
            GuessWeight(
                name=export.base_name,
                distance=distance,
                weight=weight_by_id.get(id(export), 0.0),
            )
            for export, distance in distances
        ]
        mixture = sum(
            weight * _matrix(export)
            for (export, _), weight in zip(kept, raw)
        )
        notes = ()
    overlap = np.array(template.overlap, dtype=float)
    metric = mixture.T @ overlap @ mixture
    eigenvalues, eigenvectors = np.linalg.eigh(metric)
    if eigenvalues.min() <= SINGULAR_FLOOR:
        raise GuessError(
            "the weighted mixture is (near-)singular: its overlap metric has a "
            f"eigenvalue {eigenvalues.min():.3g} <= {SINGULAR_FLOOR:g}, so the "
            "orthonormalisation is meaningless. Next step: check the neighbours (a "
            "duplicated structure or a broken export) or raise delta."
        )
    inverse_root = eigenvectors @ np.diag(eigenvalues ** -0.5) @ eigenvectors.T
    orthonormal = mixture @ inverse_root
    residual = float(np.abs(orthonormal.T @ overlap @ orthonormal - np.eye(orthonormal.shape[1])).max())
    source_export = min(kept, key=lambda item: item[1])[0]
    checks = (
        f"the mixture C = sum_beta w_beta C_beta / sum w_beta with w_beta = 1/d, over the "
        f"neighbourhood delta = {f'{delta:g} Angstrom' if delta is not None else 'all given neighbours'} "
        "(the source's Eqs. (11)-(13)); d is the RMSD without alignment (the "
        "coefficients live in each structure's own AO frame)",
        f"the mixture is orthonormalised in the target geometry's overlap metric "
        f"(Loewdin); the residual max |C^T S C - I| is {residual:.2e}",
        f"the occupations and orbital energies written with the guess come from the "
        f"nearest neighbour ('{source_export.base_name}'); the orbitals themselves are "
        "the mixture",
    )
    all_notes = notes + (
        "The written file is a Molekel mkl: run ``orca_2mkl <name>.fbk -gbw`` to turn "
        "it into a ``.gbw`` and read it with ``!moread`` + ``%moinp`` (measured: ORCA "
        "accepts the converted file as INITIAL GUESS: MOREAD).",
        "This module builds the guess; ORCA re-optimises the orbitals. The mixture is "
        "not a converged orbital set and is not claimed to be one.",
    )
    return GuessResult(
        coefficients=tuple(tuple(float(value) for value in row) for row in orthonormal),
        weights=tuple(weights),
        delta=delta,
        residual=residual,
        source=source_export.base_name,
        occupations=tuple(float(value) for value in source_export.mo_occupations),
        energies=tuple(float(value) for value in source_export.mo_energies),
        checks=checks,
        notes=all_notes,
    )


def write_guess_mkl(template_mkl: MklFile, result: GuessResult, path) -> None:
    """Write the guess into a copy of the template mkl (the gbw-ready artifact)."""
    try:
        updated = template_mkl.with_orbitals(
            result.coefficients,
            occupations=result.occupations,
            energies=result.energies,
        )
    except MklError as exc:  # the template and the exports disagree about the basis
        raise GuessError(
            f"the template mkl does not accept the interpolated coefficients: {exc}"
        ) from exc
    Path(path).write_text(updated.render(), encoding="utf-8")


def render(result: GuessResult) -> str:
    """The weights, the residual and the boundaries, as the menu prints them."""
    lines = [
        "WASP initial guess for a geometry series (weighted orbital interpolation, "
        "1/d over the neighbourhood):",
        "  neighbour            d (Angstrom)      weight",
    ]
    for weight in result.weights:
        lines.append(f"  {weight.name:<20} {weight.distance:>12.6f} {weight.weight:>12.6f}")
    lines += [
        "",
        f"  orthonormalisation residual: {result.residual:.2e}",
        f"  occupations/energies from:   {result.source}",
        "",
        "Criteria and boundaries:",
    ]
    for item in result.checks + result.notes:
        lines.append(f"  - {item}")
    return "\n".join(lines)



def evidence() -> tuple[Evidence, ...]:
    """Provenance of the WASP protocol and of the write-back route."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "WASP: the initial guess for a new geometry of a series as the "
                "1/d-weighted interpolation of the neighbouring structures' orbital "
                "coefficients (Eqs. (11)-(13) of the source), orthonormalised before "
                "use; the source presents it as the route that keeps the active space "
                "consistent across a geometry series, transferable to CASSCF / CASPT2 / "
                "NEVPT2 because all of them hinge on a reliable initial guess."
            ),
            ref=(
                "the automatic-active-space review (WASP section, Eqs. (11)-(13); "
                "the project's close reading 细读_自动活性空间综述_ChemRev2026.md section 3.7)"
            ),
            url="https://doi.org/10.1021/acs.chemrev.5c00866",
            bibkey="wardzala2026multireference",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "The ORCA write-back route this module needs, measured on 6.1.1 "
                "(2026-09-27): ``orca_2mkl <base> -mkl`` produces a plain-text Molekel "
                "file (geometry in Angstrom, verbatim basis, grouped coefficients, "
                "occupations); after editing, ``orca_2mkl <new> -gbw`` converts back and "
                "ORCA accepts the result as INITIAL GUESS: MOREAD. Round-trip fidelity "
                "5e-8 (the file's print precision); the mkl's p-shell component order "
                "was measured to be (y, z, x) against the export's (z, x, y) and is "
                "translated by the reader/writer."
            ),
            ref="tests/test_mkl.py, tests/test_guess_transfer.py; fixtures/orca/n2_scan_*",
        ),
    )
