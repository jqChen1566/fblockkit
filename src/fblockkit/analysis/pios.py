"""PiOS: the pi-orbital active space for conjugated systems (Sayfutyarova & Hammes-Schiffer 2019).

The protocol builds the active space of a conjugated pi-system from a
single-reference wavefunction, in the paper's four steps:

1. **the plane** (its Eqs. (1)-(3)): from the positions of the pi-system
   atoms, the centroid, the (unweighted) inertia tensor, and its smallest
   eigenvector ``n`` -- the direction of the local ``p_z'`` orbitals,
   ``p'_z,A = n_x p_x,A + n_y p_y,A + n_z p_z,A``;
2. **the electron count** (its Section IIC): bonds from the covalent-radius
   criterion ``r_ij < 1.3 (R_i + R_j)``, the per-atom pi contributions
   (an sp2 carbon counts one; N/P count one on two sigma bonds and two on
   three), adjustable per atom and for the total charge;
3. **the projection and selection** (its Eqs. (11)-(15)): with the oriented
   p'_z columns collected in ``O``,
   ``S_AO,pi = O^T S O``, ``X = O^T S C_occ``,
   ``S_occ,pi = X^T S_AO,pi^-1 X``; its eigenvectors rotated into the
   occupied block, and the ``N_pi,occ = N_pi,e / 2`` orbitals with the
   largest eigenvalues selected (the virtual block the same way);
4. **the semicanonicalisation** (its Eq. (16)): the selected blocks are
   rotated to diagonalise the Fock matrix within themselves, giving
   energy-ordered pi orbitals.

The orbitals the module writes are a full orthonormal set in ORCA's
partition order -- inactive occupied | pi occupied | pi virtual | inactive
virtual -- so a ``%casscf`` run with ``nel = N_pi,e``, ``norb = |M|`` finds
the pi space as its active window (the same write-back route as menu 47).

What this ORCA route changes, and what it measures
--------------------------------------------------

- The source builds the oriented orbitals in an auxiliary minimal basis
  (MINAO) and projects through IAOs (its Eqs. (6)-(8)); that needs the
  cross-basis overlap ``S_12``, which ORCA does not export.  This module
  uses the calculation's own basis instead -- each pi atom's *valence p
  shell* (the outermost p shell of that atom, measured to be the last one
  in the export's label order) -- the same substitution menu 14's AVAS
  route carries.  The projection eigenvalues are the paper's own validity
  measure for the choice of atoms and electron count, and the report prints
  the selection spectrum (selected vs excluded) for that reading.
- The Fock matrix for Eq. (16) is assembled as ``H + J + K`` from the
  export's blocks (``1elIntegrals: ["H", "S"]`` and ``FockMatrix: ["J",
  "K"]``), cross-checked in the tests against the canonical orbital
  energies.
- The written basis of every block is a Gram-Schmidt of the *parents*
  projected onto it (the pi block from the parents with the largest pi
  character, the inactive block from the rest, in index order), so it
  never inherits the diagonaliser's BLAS-dependent choice inside a
  degenerate cluster: under a 1e-14 input perturbation the written
  coefficients move at the 1e-12 level (measured; the examples' captures
  are compared across machines).

Scope: one pi-system per call (the paper's multiple-fragment extension is
registered); closed-shell SCF exports (the source's closed-shell Slater
determinant); approximately planar atom sets (the report prints the maximum
out-of-plane deviation and refuses beyond :data:`PLANARITY_CEILING`).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..knowledge.models import (
    EVIDENCE_LITERATURE,
    EVIDENCE_MEASURED,
    Evidence,
    ReportSection,
)
from ..parsers.mkl import MklError, MklFile
from ..parsers.orca_json import OrcaJson

__all__ = [
    "PiosError",
    "PiosResult",
    "select_pi_space",
    "write_mkl",
    "render",
    "run",
    "evidence",
]

#: Single-bond covalent radii in Angstrom (Cordero et al., Dalton Trans. 2008,
#: 2832-2838 -- the source's Ref. 58).  Elements outside this table are
#: refused: the sigma-bond count of *every* atom of the molecule enters the
#: connectivity, so no radius means no degree.
COVALENT_RADII = {
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39,
}

#: The source's bonding criterion: bonded when r_ij < 1.3 (R_i + R_j).
BOND_FACTOR = 1.3

#: Beyond this maximum out-of-plane deviation (Angstrom) the atom set is not
#: the source's "approximately planar" case and the pi direction is refused.
PLANARITY_CEILING = 0.5

class PiosError(ValueError):
    """The pi space cannot be built from the given inputs (with a next step)."""


@dataclass(frozen=True)
class PiosResult:
    """The pi space and everything the report needs to say about it."""

    base_name: str
    elements: tuple[str, ...]
    indices: tuple[int, ...]  # 0-based atom indices as given
    normal: tuple[float, float, float]
    centroid: tuple[float, float, float]
    max_out_of_plane: float
    degrees: tuple[int, ...]
    contributions: tuple[int, ...]
    charge: int
    pi_electrons: int
    n_occupied_pi: int
    n_virtual_pi: int
    occ_spectrum: tuple[float, ...]  # all occupied projection eigenvalues, desc
    vir_spectrum: tuple[float, ...]
    occ_energies: tuple[float, ...]  # semicanonical pi energies, energy order
    vir_energies: tuple[float, ...]
    occ_parents: tuple[tuple[tuple[int, float], ...], ...]  # per pi orbital
    vir_parents: tuple[tuple[tuple[int, float], ...], ...]
    coefficients: tuple[tuple[float, ...], ...]  # n_ao x n_mo, export row order
    occupations: tuple[float, ...]
    energies: tuple[float, ...]
    residual: float  # max |C'^T S C' - I|
    checks: tuple[str, ...]
    notes: tuple[str, ...]


def _element_contribution(element: str, degree: int) -> int:
    """The source's per-atom pi contribution (Section IIC)."""
    if element == "C":
        return 1
    if element in ("N", "P"):
        if degree == 2:
            return 1
        if degree == 3:
            return 2
        raise PiosError(
            f"the {element} atom is connected by {degree} sigma bonds; the source's "
            "rule covers two (one pi electron) and three (two pi electrons). Next "
            "step: give this atom's contribution explicitly in the manifest."
        )
    raise PiosError(
        f"the source's automatic rule covers C, N and P; '{element}' is outside it. "
        "Next step: give the atom's contribution explicitly in the manifest."
    )


def _degrees(export: OrcaJson, indices: tuple[int, ...]) -> tuple[int, ...]:
    """The sigma-bond count of every pi atom, from the covalent-radius criterion."""
    coordinates = np.asarray(export.coordinates, dtype=float)
    elements = tuple(export.atoms)
    for element in elements:
        if element not in COVALENT_RADII:
            raise PiosError(
                f"the covalent radius of '{element}' is not in the table, so the "
                "connectivity cannot be built. Next step: use a system of the "
                "tabulated main-group elements."
            )
    positions = coordinates[list(indices)]
    radii = np.array([COVALENT_RADII[elements[index]] for index in indices])
    degrees = []
    for position, (atom, radius) in enumerate(zip(indices, radii)):
        limit = BOND_FACTOR * (radius + np.array(
            [COVALENT_RADII[elements[other]] for other in range(len(elements))]
        ))
        delta = coordinates - positions[position]
        distances = np.linalg.norm(delta, axis=1)
        bonded = (distances < limit) & (np.arange(len(elements)) != atom)
        degrees.append(int(bonded.sum()))
    return tuple(degrees)


def _normal(export: OrcaJson, indices: tuple[int, ...]):
    """The centroid, the pi direction and the planarity measure (Eqs. (1)-(3))."""
    coordinates = np.asarray(export.coordinates, dtype=float)[list(indices)]
    centroid = coordinates.mean(axis=0)
    delta = coordinates - centroid
    tensor = delta.T @ delta  # the unweighted inertia tensor of Eq. (2)
    eigenvalues, eigenvectors = np.linalg.eigh(tensor)
    if eigenvalues[1] - eigenvalues[0] < 1e-6 * (eigenvalues[2] + 1e-12):
        raise PiosError(
            "the atom set has (near-)degenerate smallest moments of inertia "
            "(collinear or rotationally symmetric), so the pi direction is "
            "undefined. Next step: check the atom list (the pi system must span "
            "a plane)."
        )
    normal = eigenvectors[:, 0]
    # canonical orientation: first non-negligible component positive
    for value in normal:
        if abs(value) > 1e-9:
            if value < 0.0:
                normal = -normal
            break
    deviation = float(np.abs(delta @ normal).max())
    if deviation > PLANARITY_CEILING:
        raise PiosError(
            f"the atom set is not approximately planar (maximum out-of-plane "
            f"deviation {deviation:.3f} Angstrom > {PLANARITY_CEILING:g}). Next "
            "step: check the atom list."
        )
    return centroid, normal, deviation


def _oriented_matrix(export: OrcaJson, indices: tuple[int, ...], normal: np.ndarray):
    """The oriented p'_z columns ``O`` (n_ao x |M|) on each atom's valence p shell."""
    labels = export.ao_labels
    if labels is None:
        raise PiosError(
            "the export carries no AO labels, so the p shells cannot be located. "
            "Next step: re-export with the orca_2json default (the labels come with "
            "the MOCoefficients request)."
        )
    columns = np.zeros((export.n_ao, len(indices)))
    shells_used = []
    for column, atom in enumerate(indices):
        shell_numbers = sorted(
            {label.shell for label in labels if label.center == atom and label.angular == "p"}
        )
        if not shell_numbers:
            raise PiosError(
                f"the atom {atom} ({export.atoms[atom]}) carries no p shell in the "
                "export, so it cannot contribute a pi orbital. Next step: check the "
                "atom list."
            )
        shell = shell_numbers[-1]  # the outermost (valence) p shell
        shells_used.append(shell)
        for label_index, label in enumerate(labels):
            if label.center != atom or label.angular != "p" or label.shell != shell:
                continue
            component = {"x": normal[0], "y": normal[1], "z": normal[2]}.get(label.component)
            if component is None:
                raise PiosError(
                    f"the p component '{label.component}' of atom {atom} is not a "
                    "Cartesian direction. Next step: check the export's labels."
                )
            columns[label_index, column] = component
    return columns, tuple(shells_used)


def _fock_matrix(export: OrcaJson) -> np.ndarray:
    """The full Fock matrix ``H + J + K`` in the AO basis (measured composition)."""
    for name, block in (
        ("H-Matrix", export.hamiltonian),
        ("J-Matrix", export.coulomb),
        ("K-Matrix", export.exchange),
    ):
        if block is None:
            raise PiosError(
                f"the export carries no {name} block, so the Fock matrix for the "
                "semicanonicalisation cannot be built. Next step: re-export with "
                '{"MOCoefficients": true, "1elIntegrals": ["H", "S"], "FockMatrix": '
                '["J", "K"]} written as <base>.json.conf next to the .gbw.'
            )
        if name != "H-Matrix" and len(block) != 1:
            raise PiosError(
                f"the {name} block carries {len(block)} spin matrices; this menu "
                "covers closed-shell (RHF) exports. Next step: run the SCF as RHF "
                "and re-export."
            )
    hamiltonian = np.asarray(export.hamiltonian, dtype=float)
    coulomb = np.asarray(export.coulomb[0], dtype=float)
    exchange = np.asarray(export.exchange[0], dtype=float)
    return hamiltonian + coulomb + exchange


#: The pi-block parents must stand out from the next candidate by this margin:
#: below it the partition of parents between the pi and the inactive block is
#: numerically ambiguous, and the menu refuses instead of guessing.
PARENT_MARGIN = 1e-6


def _parent_norms(projector, parents, overlap):
    """The S-norm of every parent SCF orbital projected onto a subspace."""
    return np.array(
        [
            float(
                np.sqrt(
                    max(float((projector @ parent) @ overlap @ (projector @ parent)), 0.0)
                )
            )
            for parent in parents.T
        ]
    )


def _gram_schmidt(indices, projections, overlap, count, *, block):
    """S-metric Gram-Schmidt of the given projected parents (twice-orthogonalised).

    The accepted count must reach ``count`` exactly; a shortfall means the
    block is not spanned by the parent set (or a direction is numerically
    consumed), and the menu refuses instead of writing an unstable basis.
    """
    basis = []
    for index in indices:
        if len(basis) == count:
            break
        vector = projections[int(index)]
        for _ in range(2):  # twice is enough: the second pass recovers orthogonality
            for previous in basis:
                vector = vector - float(previous @ overlap @ vector) * previous
        norm = float(np.sqrt(max(float(vector @ overlap @ vector), 0.0)))
        if norm > 1e-6:
            basis.append(vector / norm)
    if len(basis) != count:
        raise PiosError(
            f"the deterministic basis of the {block} block could not be built "
            "(the parent Gram-Schmidt fell short). Next step: report this input "
            "with its export; the block spans directions outside the parent set."
        )
    result = np.array(basis).T
    weights_ao = 2.0 ** -np.arange(result.shape[0])
    for column in range(result.shape[1]):
        if float(weights_ao @ result[:, column]) < 0.0:
            result[:, column] *= -1
    return result


def _block_bases(projector, parents, overlap, count, *, block):
    """The pi block and the inactive block of one manifold from the parents.

    Both bases are S-metric Gram-Schmidt constructions of the parents projected
    onto the block (twice-orthogonalised, near-unity residuals because the
    parents of a conjugated system's pi space are its canonical pi orbitals,
    which are S-orthogonal).  The ancestors are partitioned first, and the
    partition is the only order-sensitive step: the ``count`` parents with the
    largest projection norm form the pi block (ordered by descending norm,
    index for ties), the rest form the inactive block (processed in index
    order).  The partition refuses when the last selected and the first
    rejected parent lie within :data:`PARENT_MARGIN` -- there the two blocks
    are numerically ambiguous.  Under a 1e-13 input perturbation the written
    bases move at the 1e-11 level (measured); the diagonaliser's own
    within-degeneracy choice is not inherited because the subspaces enter only
    through their projectors, which are eigenvector sums and stable.
    """
    norms = _parent_norms(projector, parents, overlap)
    # descending projection norm; within 1e-9 the parent index decides (the
    # norms of symmetry-equivalent parents are equal up to noise, and an order
    # swap between such near-ties moves the written basis by the tie size only)
    by_norm = np.argsort(-norms, kind="stable")
    clusters: list[list[int]] = []
    for index in by_norm:
        if clusters and abs(norms[clusters[-1][-1]] - norms[index]) <= 1e-9:
            clusters[-1].append(int(index))
        else:
            clusters.append([int(index)])
    order = [index for cluster in clusters for index in sorted(cluster)]
    selected = order[:count]
    chosen = set(selected)
    if norms[order[count - 1]] - norms[order[count]] < PARENT_MARGIN:
        raise PiosError(
            f"the {block} block is not cleanly separated in the parent set: the "
            f"last selected parent reads {norms[order[count - 1]]:.6f} and the "
            f"first rejected {norms[order[count]]:.6f} (margin "
            f"{PARENT_MARGIN:g}). Next step: check the atom set and the electron "
            "count; the partition between the pi and inactive parents is ambiguous "
            "here."
        )
    pi_projections = [projector @ parent for parent in parents.T]
    rest_projections = [
        parent - projected for parent, projected in zip(parents.T, pi_projections)
    ]
    pi_basis = _gram_schmidt(
        selected, pi_projections, overlap, count, block=f"pi ({block})"
    )
    rejected = [int(index) for index in range(len(norms)) if int(index) not in chosen]
    rest_basis = _gram_schmidt(
        rejected, rest_projections, overlap, len(rejected), block=f"inactive ({block})"
    )
    return pi_basis, rest_basis


def _parent_weights(orbital, parents, overlap, count=2):
    """The ``count`` parent SCF orbitals carrying the largest weight in an orbital."""
    weights = np.abs(parents.T @ overlap @ orbital) ** 2
    order = np.argsort(-weights)[:count]
    return tuple((int(index), float(weights[index])) for index in order)


def select_pi_space(
    export: OrcaJson,
    indices,
    *,
    charge: int | None = None,
    contributions=None,
    pi_electrons: int | None = None,
) -> PiosResult:
    """Build the pi-orbital active space of the atom set ``indices`` (0-based)."""
    if export.overlap is None:
        raise PiosError(
            f"the export '{export.base_name}' carries no S-Matrix. Next step: "
            "re-export with the S-Matrix (the orca_2json default)."
        )
    if not export.mo_coefficients:
        raise PiosError(
            "the export carries no MO coefficients. Next step: re-export with "
            '{"MOCoefficients": true}.'
        )
    atoms = tuple(int(index) for index in indices)
    if len(atoms) < 3:
        raise PiosError(
            f"the pi system has {len(atoms)} atom(s); a pi plane needs at least "
            "three. Next step: check the atom list."
        )
    if len(set(atoms)) != len(atoms):
        raise PiosError("the atom list repeats an index. Next step: give each atom once.")
    for atom in atoms:
        if not 0 <= atom < len(export.atoms):
            raise PiosError(
                f"the atom index {atom} is outside the export's {len(export.atoms)} "
                "atoms. Next step: check the atom list (0-based, the export's order)."
            )
    overlap = np.asarray(export.overlap, dtype=float)
    coefficients = np.asarray(export.mo_coefficients, dtype=float).T  # AO-major
    occupations = np.asarray(export.mo_occupations, dtype=float)
    occupied = [index for index, value in enumerate(occupations) if value > 1.99]
    virtuals = [index for index, value in enumerate(occupations) if value < 0.01]
    if len(occupied) + len(virtuals) != len(occupations):
        raise PiosError(
            f"the export '{export.base_name}' has fractional occupations (an active "
            "space), but the source's protocol starts from a closed-shell Slater "
            "determinant. Next step: export the SCF run (RHF) the pi space is built "
            "from."
        )
    if len(occupied) == 0:
        raise PiosError(
            "the export carries no doubly occupied orbital. Next step: check the "
            "occupations (a closed-shell SCF export is expected)."
        )
    fock = _fock_matrix(export)
    centroid, normal, deviation = _normal(export, atoms)
    oriented, shells = _oriented_matrix(export, atoms, normal)
    pi_overlap = oriented.T @ overlap @ oriented
    eigenvalues = np.linalg.eigvalsh(pi_overlap)
    if eigenvalues.min() < 1e-8:
        raise PiosError(
            "the oriented p orbitals are linearly dependent (the smallest overlap "
            f"eigenvalue is {eigenvalues.min():.2e}). Next step: check the atom "
            "list (this arises for duplicated or coincident atoms)."
        )
    inverse_overlap = np.linalg.inv(pi_overlap)
    degrees = _degrees(export, atoms)
    if contributions is None:
        counted = tuple(_element_contribution(export.atoms[atom], degree)
                        for atom, degree in zip(atoms, degrees))
    else:
        counted = tuple(int(value) for value in contributions)
        if len(counted) != len(atoms):
            raise PiosError(
                f"the manifest gives {len(counted)} contributions for {len(atoms)} "
                "atoms. Next step: give one contribution per atom, in the same order."
            )
    effective_charge = int(export.charge if charge is None else charge)
    if pi_electrons is None:
        total = sum(counted) - effective_charge
    else:
        total = int(pi_electrons)
    if total <= 0 or total > 2 * len(atoms) or total % 2 != 0:
        raise PiosError(
            f"the pi electron count {total} is not an even number in "
            f"1..{2 * len(atoms)}. Next step: check the contributions and the charge "
            "(or give it explicitly in the manifest)."
        )
    n_occupied = total // 2
    n_virtual = len(atoms) - n_occupied
    if n_occupied > len(occupied) or n_virtual > len(virtuals):
        raise PiosError(
            f"the pi space needs {n_occupied} occupied and {n_virtual} virtual "
            f"orbitals while the export has {len(occupied)} and {len(virtuals)}. "
            "Next step: check the electron count and the atom list."
        )

    def block_spectrum(block_indices):
        block = coefficients[:, block_indices]
        x = oriented.T @ overlap @ block
        spectrum = x.T @ inverse_overlap @ x
        values, vectors = np.linalg.eigh(spectrum)
        values = np.clip(values, 0.0, None)  # a projection spectrum is non-negative
        order = np.argsort(-values)
        return block @ vectors[:, order], values[order], x  # x: pi AO rows

    occ_rotated, occ_values, _ = block_spectrum(occupied)
    vir_rotated, vir_values, _ = block_spectrum(virtuals)
    occ_selected = occ_rotated[:, :n_occupied]
    vir_selected = vir_rotated[:, :n_virtual]
    occ_inactive = occ_rotated[:, n_occupied:]
    vir_inactive = vir_rotated[:, n_virtual:]
    occ_parents = coefficients[:, occupied]
    vir_parents = coefficients[:, virtuals]

    def projector(columns):
        """The S-orthogonal projector onto a subspace given S-orthonormal columns."""
        return (columns @ columns.T) @ overlap

    occ_pi, occ_rest = _block_bases(
        projector(occ_selected), occ_parents, overlap, n_occupied, block="occupied"
    )
    vir_pi, vir_rest = _block_bases(
        projector(vir_selected), vir_parents, overlap, n_virtual, block="virtual"
    )
    # the paper's semicanonical energies: the Fock eigenvalues of each selected
    # block (basis-independent values); the per-orbital tags written into the mkl
    # are the Fock Rayleigh quotients of the written orbitals themselves
    occ_energies = tuple(
        float(value) for value in np.linalg.eigvalsh(occ_pi.T @ fock @ occ_pi)
    )
    vir_energies = tuple(
        float(value) for value in np.linalg.eigvalsh(vir_pi.T @ fock @ vir_pi)
    )

    full = np.hstack([occ_rest, occ_pi, vir_pi, vir_rest])
    residual = float(np.abs(full.T @ overlap @ full - np.eye(full.shape[1])).max())
    occupations = tuple([2.0] * len(occupied) + [0.0] * len(virtuals))
    energies = tuple(float(column @ fock @ column) for column in full.T)
    occ_blocks = tuple(
        _parent_weights(occ_pi[:, column], occ_parents, overlap)
        for column in range(occ_pi.shape[1])
    )
    vir_blocks = tuple(
        _parent_weights(vir_pi[:, column], vir_parents, overlap)
        for column in range(vir_pi.shape[1])
    )
    def rayleigh(orbital):
        return float(np.linalg.norm(oriented.T @ overlap @ orbital))
    checks = (
        "the construction follows the source's Eqs. (1)-(16): the plane from the "
        "inertia tensor, the oriented valence p orbitals, the occupied/virtual "
        "projection spectra, the selection of N_pi,occ = N_pi,e / 2 and its virtual "
        "partner, and the Fock semicanonicalisation of the selected blocks",
        f"the oriented orbitals live on each atom's valence p shell (the outermost "
        f"p shell of the export; used shells: {', '.join(str(shell) for shell in shells)}); "
        "the source builds them in an auxiliary MINAO basis through IAOs, which "
        "needs the cross-basis overlap ORCA does not export -- this is the "
        "calculation's own basis, the same substitution menu 14 carries",
        f"the Fock matrix is assembled as H + J + K from the export's blocks; the "
        f"written set is orthonormal to {residual:.1e} in the export's overlap",
        "the written basis of every block is the Gram-Schmidt of the parent SCF "
        "orbitals projected onto it (pi character decides the partition; index "
        "order within each pool), so the files do not inherit the "
        "diagonaliser's BLAS-dependent choice of a degenerate cluster's basis "
        "(a determinism convention, not part of the source)",
    )
    notes = (
        "The written file is a Molekel mkl in ORCA's partition order -- inactive "
        "occupied | pi occupied | pi virtual | inactive virtual -- so a %casscf "
        f"run with nel {total}, norb {len(atoms)} takes the pi space as its active "
        "window: run ``orca_2mkl <name>.fbk -gbw`` and start with ``!moread`` + "
        "``%moinp`` (the menu 47 route)",
        "Scope: one pi system per call and closed-shell exports; the projection "
        "eigenvalues in the report are the source's own validity measure for the "
        "chosen atoms and electron count (a selected eigenvalue below an excluded "
        "one means the choices deserve a second look)",
    )
    return PiosResult(
        base_name=export.base_name,
        elements=tuple(export.atoms[atom] for atom in atoms),
        indices=atoms,
        normal=tuple(float(value) for value in normal),
        centroid=tuple(float(value) for value in centroid),
        max_out_of_plane=deviation,
        degrees=degrees,
        contributions=counted,
        charge=effective_charge,
        pi_electrons=total,
        n_occupied_pi=n_occupied,
        n_virtual_pi=n_virtual,
        occ_spectrum=tuple(float(value) for value in occ_values),
        vir_spectrum=tuple(float(value) for value in vir_values),
        occ_energies=occ_energies,
        vir_energies=vir_energies,
        occ_parents=occ_blocks,
        vir_parents=vir_blocks,
        coefficients=tuple(tuple(float(value) for value in row) for row in full),
        occupations=occupations,
        energies=tuple(float(value) for value in energies),
        residual=residual,
        checks=checks,
        notes=notes,
    )


def write_mkl(template_mkl: MklFile, result: PiosResult, path) -> None:
    """Write the pi-space orbitals into a copy of the system's mkl.

    The same write-back route menu 18 uses (the recipe layer's writer): the
    template supplies the basis and the labels, the coefficients carry the
    orca_2json row order, and ``orca_2mkl <name> -gbw`` turns the result into
    a gbw ORCA reads through ``!moread``.
    """
    try:
        updated = template_mkl.with_orbitals(
            result.coefficients,
            occupations=result.occupations,
            energies=result.energies,
        )
    except MklError as exc:  # the template and the coefficients disagree about the basis
        raise PiosError(
            f"the mkl does not accept the pi-space coefficients: {exc}"
        ) from exc
    Path(path).write_text(updated.render(), encoding="utf-8")


def render(result: PiosResult) -> str:
    """The plane, the atoms, the spectra and the boundaries, as the menu prints them."""
    lines = [
        "PiOS pi-orbital active space (Huckel-style construction on the oriented "
        "valence p orbitals):",
        f"  system: {result.base_name}",
        "  atom    element   sigma-bonds   pi electrons",
    ]
    for position, atom in enumerate(result.indices):
        lines.append(
            f"  {atom:>4}    {result.elements[position]:<7}   "
            f"{result.degrees[position]:>6}       {result.contributions[position]:>6}"
        )
    lines += [
        "",
        f"  pi plane normal:       ({result.normal[0]:+.6f}, {result.normal[1]:+.6f}, "
        f"{result.normal[2]:+.6f})",
        f"  centroid:              ({result.centroid[0]:+.6f}, {result.centroid[1]:+.6f}, "
        f"{result.centroid[2]:+.6f})",
        f"  max out-of-plane:      {result.max_out_of_plane:.6f} Angstrom",
        f"  pi electrons:          {result.pi_electrons} "
        f"(charge {result.charge:+d}; {sum(result.contributions)} from the atoms)",
        f"  active space:          CAS({result.pi_electrons}e, {len(result.indices)}o) "
        f"= {result.n_occupied_pi} occupied + {result.n_virtual_pi} virtual pi orbitals",
        "",
        "Occupied projection spectrum (selected first):",
        "   orbital      sigma",
    ]
    for position, value in enumerate(result.occ_spectrum):
        mark = "  <- selected" if position < result.n_occupied_pi else ""
        lines.append(f"   {position:>7}  {value:>9.6f}{mark}")
    lines += ["", "Virtual projection spectrum (selected first):", "   orbital      sigma"]
    for position, value in enumerate(result.vir_spectrum):
        mark = "  <- selected" if position < result.n_virtual_pi else ""
        lines.append(f"   {position:>7}  {value:>9.6f}{mark}")
    lines += ["", "Selected pi orbitals (semicanonical energies, parent SCF orbitals):"]
    for label, energies, parents in (
        ("occ", result.occ_energies, result.occ_parents),
        ("vir", result.vir_energies, result.vir_parents),
    ):
        for position, (energy, weight) in enumerate(zip(energies, parents)):
            roots = ", ".join(f"{index} ({value:.3f})" for index, value in weight)
            lines.append(f"   {label} {position}:  E = {energy:+.6f} Eh   parents: {roots}")
    lines += ["", "Criteria and boundaries:"]
    for item in result.checks + result.notes:
        lines.append(f"  - {item}")
    return "\n".join(lines)


def run(export: OrcaJson, indices, **kwargs) -> ReportSection:
    """The analyser entry point for menu 48."""
    result = select_pi_space(export, indices, **kwargs)
    return ReportSection(
        title="A13 PiOS pi-orbital active space (oriented-p projection)",
        body=render(result),
    )


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the PiOS protocol and of its ORCA route."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "PiOS: build the active space of a conjugated pi-system from a "
                "single-reference wavefunction -- the pi plane from the positions "
                "of the pi atoms (Eqs. (1)-(3)), the pi electron count from atomic "
                "connectivity with covalent radii (Section IIC), the occupied and "
                "virtual pi orbitals from the eigenvectors of the projected overlap "
                "S_occ,pi = X^T S_AO,pi^-1 X (Eqs. (11)-(15)), and the Fock "
                "semicanonicalisation of the selected blocks (Eq. (16)).  The "
                "source's own validity measure is the projection spectrum, and its "
                "benzene benchmark reports the SVD eigenvalues between the initial "
                "pi space and the CASSCF-optimised active space as 0.9708, 0.9709, "
                "0.9875, 0.9998, 0.9998, 1.0 -- the pi orbitals are near-optimal "
                "from the outset."
            ),
            ref=(
                "Sayfutyarova E. R., Hammes-Schiffer S., J. Chem. Theory Comput., "
                "2019, 15, 1679-1689, DOI 10.1021/acs.jctc.8b01196 (Eqs. (1)-(16); "
                "the benzene benchmark in Section III A)"
            ),
            url="https://doi.org/10.1021/acs.jctc.8b01196",
            bibkey="sayfutyarova2019constructing",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "The ORCA route, measured on 6.1.1 (fixtures/orca/benzene_rhf.*): "
                "the RHF/cc-pVDZ export (S, H, J, K blocks) gives the six pi "
                "orbitals of benzene -- the three occupied selected at projection "
                "eigenvalues 0.7789 / 0.7649 / 0.7649 with the next-largest "
                "excluded at 0.000, the three virtual at 1.000 / 1.000 / 1.000 "
                "(next excluded 0.235); the occupied-side deficit is the d "
                "polarisation a pure-p projection cannot carry and the source's "
                "IAO route does -- the documented substitution of this menu.  The "
                "projector trace is exact on both sides (2.309 + 3.691 = 6, the "
                "dimension of the p_z space), and the Fock diagonal of the written "
                "orbitals reproduces the canonical orbital energies (the same "
                "cross-check as menu 21's Fock assembly).  End to end: a "
                "CASSCF(6,6) started from the written set (``!NoIter moread`` + "
                "``%casscf MaxIter``) converges to the fixture's solution "
                "(-230.793818903 vs -230.793818898 Eh) in 13 macro-iterations "
                "where the aufbau start takes 7, and the SVD eigenvalues between "
                "the guess's active space and the converged one read 0.9999 / "
                "0.9999 / 0.9999 / 0.7573 / 0.7573 / 0.6558 -- the virtual gap "
                "equals the converged pi* orbitals' own out-of-plane character "
                "(0.7573 / 0.6558, measured directly), i.e. exactly the "
                "polarisation the pure-p route cannot carry (the source reports "
                "0.9708-1.0 with its IAO/aug-cc-pVTZ setup)."
            ),
            ref="tests/test_pios.py; fixtures/orca/benzene_rhf.*, benzene_pios.*",
        ),
    )
