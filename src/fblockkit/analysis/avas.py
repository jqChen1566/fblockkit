"""AVAS: an active space built from target atomic orbitals (Sayfutyarova et al. 2017).

The construction (the paper's Eqs. 2-9, 13):

1. choose a set A of target AOs (the paper uses a minimal auxiliary basis; here
   the target is a subset of the calculation's own AOs -- the "non-minimal ANO"
   route the paper describes in its outlook, which is also the f-block route of
   record because the explicit list does not go through the minimal basis at
   all: the ecosystem's minimal-basis coverage has no entry for Eu, measured
   through the engine's own localization tooling);
2. the projector onto ``span(A)`` in the AO metric, applied to a coefficient
   vector ``c``, is ``c -> i sigma^{-1} i^T S c`` with ``i`` the embedding of the
   A rows and ``sigma = S[AA]`` the target overlap;
3. the overlap spectra of that projector over the occupied and the virtual
   orbital blocks are ``S^A_occ = (S C_occ)_A^T sigma^{-1} (S C_occ)_A`` and the
   same for the virtual block; at most ``|A|`` eigenvalues are non-zero because
   ``span(A)`` has that dimension;
4. rotating the blocks by the eigenvectors splits them into the part that
   overlaps the target (into the active space) and the part that does not (core
   or truly virtual), and the truncation threshold (0.05-0.1 in the source, the
   only user-chosen number) drops the negligible overlaps.

What this module reports is therefore the *spectra* and the resulting
``(n_electrons, n_orbitals)`` together with the paper's readings -- the virtual
side decides whether AVAS suits the system at all (all-small virtual overlaps
mean the antibonding target character is ligand-centred, the [CuCl4]2- case,
and the remedy is to add ligand AOs to the target), and the occupied side's low
eigenvalues name the orbitals that are strongly mixed with the environment.
Two boundaries are stated in the report: the source's own "both criteria can
falsify but not confirm", and the ORCA-side one -- ORCA's CASSCF takes its
active space by *orbital order*, not by an index list, so the deliverable here
is a verdict on a window plus the (n_el, n_orb) for the input; the ready block
for ORCA's built-in AVAS is emitted with the f-block caveat.

Open shells (the paper's Section 2.5): option 2 (use the alpha orbitals only)
can move an unoccupied beta orbital into the core and make CASCI lie above the
variational HF energy; option 3 keeps every singly occupied orbital in the
active space and is the default recommended here for f-block targets.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..parsers.orca_json import ANGULAR_LETTERS
from .solid_harmonics import component_labels
from ..knowledge.models import (
    EVIDENCE_LITERATURE,
    EVIDENCE_MEASURED,
    Evidence,
    ReportSection,
)

__all__ = [
    "AvasError",
    "AvasAnalysis",
    "target_aos",
    "projection_spectrum",
    "analyze",
    "run",
    "evidence",
]


class AvasError(ValueError):
    """The AVAS analysis cannot run on the given data (with a next step)."""

#: Occupation splitting for the occupied / virtual blocks of the export.
OCCUPIED_ABOVE = 0.5

#: Default truncation threshold (the source's 0.05-0.1 range).
DEFAULT_THRESHOLD = 0.1

#: Default for open-shell targets: option 3 keeps every singly occupied orbital.
DEFAULT_OPTION = 3


@dataclass(frozen=True)
class AvasAnalysis:
    """The projection spectra and the resulting active-space recommendation."""

    centre: int
    element: str
    angular: int
    shells: tuple[int, ...]
    target_aos: tuple[int, ...]
    threshold: float
    option: int
    occupied_indices: tuple[int, ...]
    virtual_indices: tuple[int, ...]
    partial_indices: tuple[int, ...]
    singly_occupied: tuple[int, ...]
    occupied_overlaps: tuple[float, ...]
    virtual_overlaps: tuple[float, ...]
    active_occupied: tuple[int, ...]
    active_virtual: tuple[int, ...]
    n_electrons: float
    n_orbitals: int
    checks: tuple[str, ...]
    verdicts: tuple[str, ...]


def target_aos(labels, centre: int, angular: int, shells=None) -> tuple[int, ...]:
    """Indices of the target AOs: one centre, one angular momentum, the given shells.

    ``shells=None`` takes every shell of that angular momentum on the centre.
    """
    letter = ANGULAR_LETTERS.get(angular)
    if letter is None:
        raise AvasError(
            f"angular momentum l={angular} has no recorded label letter. Next step: "
            "extend ANGULAR_LETTERS in the export reader."
        )
    if labels is None:
        raise AvasError(
            "the export carries no AO labels, so the target AOs cannot be named. Next "
            "step: re-export with the orbital labels (the orca_2json default)."
        )
    selected = [
        index
        for index, label in enumerate(labels)
        if label.center == centre
        and label.angular == letter
        and (shells is None or label.shell in shells)
    ]
    if not selected:
        wanted = "any shell" if shells is None else f"shell(s) {sorted(shells)}"
        raise AvasError(
            f"centre {centre} has no l={angular} ({letter}) AO in {wanted}. Next step: "
            "check the centre index, the angular momentum and the shell numbers against "
            "the export's AO labels."
        )
    return tuple(selected)


def _sigma_inverse(overlap: np.ndarray, target: tuple[int, ...]) -> np.ndarray:
    block = overlap[np.ix_(target, target)]
    values, vectors = np.linalg.eigh(block)
    if values.min() <= 1e-10:
        raise AvasError(
            f"the target overlap matrix is singular (smallest eigenvalue "
            f"{values.min():.3e}): the selected target AOs are linearly dependent. "
            "Next step: remove the redundant shell from the target set."
        )
    return vectors @ np.diag(values**-1) @ vectors.T


def projection_spectrum(
    coefficients: np.ndarray, overlap: np.ndarray, target: tuple[int, ...]
) -> tuple[float, ...]:
    """Eigenvalues of the target projector over the block of orbitals given.

    ``coefficients`` has the orbitals in the columns (the block is either the
    occupied or the virtual set).  The eigenvalues lie in [0, 1]; at most
    ``len(target)`` of them are non-zero.
    """
    coef = np.asarray(coefficients, dtype=np.float64)
    s = np.asarray(overlap, dtype=np.float64)
    if coef.ndim != 2 or coef.shape[0] != s.shape[0]:
        raise AvasError(
            f"the coefficient block is {coef.shape} against an overlap matrix of "
            f"{s.shape}. Next step: pass AO coefficients (orbitals in the columns) of "
            "the same export."
        )
    if coef.shape[1] == 0:
        return ()
    sigma_inverse = _sigma_inverse(s, target)
    projected = (s @ coef)[list(target), :]  # (S C)_A
    matrix = projected.T @ sigma_inverse @ projected
    matrix = 0.5 * (matrix + matrix.T)
    values = np.linalg.eigvalsh(matrix)
    return tuple(float(value) for value in np.clip(values, 0.0, 1.0)[::-1])


def _occupation_blocks(occupations) -> tuple[list[int], list[int], list[int]]:
    occupied, virtual, partial = [], [], []
    for index, value in enumerate(occupations):
        if value > OCCUPIED_ABOVE:
            occupied.append(index)
            if value < 1.5:
                partial.append(index)
        elif value > 1e-6:
            partial.append(index)
            virtual.append(index)
        else:
            virtual.append(index)
    return occupied, virtual, partial


def analyze(
    export,
    *,
    centre: int | None = None,
    angular: int = 3,
    shells=None,
    threshold: float = DEFAULT_THRESHOLD,
    option: int = DEFAULT_OPTION,
) -> AvasAnalysis:
    """Run the AVAS projection on an orbital export and recommend an active space.

    ``centre=None`` picks the first f-block element of the export (the same
    element test the other analyses use); pass an index to override.
    """
    from ..knowledge.elements import is_f_element

    if export.overlap is None:
        raise AvasError(
            "the export carries no overlap matrix, so no projection can be formed. "
            "Next step: re-export with the S-Matrix (the orca_2json default)."
        )
    if centre is None:
        centres = [
            index for index, symbol in enumerate(export.atoms) if is_f_element(symbol)
        ]
        if not centres:
            raise AvasError(
                "the export has no f-block element to target. Next step: give the centre "
                "index explicitly."
            )
        centre = centres[0]
    if not 0 <= centre < len(export.atoms):
        raise AvasError(
            f"centre {centre} is outside the export's {len(export.atoms)} atom(s). Next "
            "step: check the index."
        )
    if option not in (2, 3):
        raise AvasError(
            f"option {option} is not 2 or 3 (the paper's open-shell choices). Next step: "
            "use 3 (keep every singly occupied orbital; the default for f-block "
            "targets) or 2 (alpha orbitals only; only when a high-spin state is "
            "excluded)."
        )
    target = target_aos(export.ao_labels, centre, angular, shells)
    coefficients = np.array(export.mo_coefficients, dtype=np.float64).T
    overlap = np.array(export.overlap, dtype=np.float64)
    occupied, virtual, partial = _occupation_blocks(export.mo_occupations)
    occupied_overlaps = projection_spectrum(
        coefficients[:, occupied], overlap, target
    )
    virtual_overlaps = projection_spectrum(coefficients[:, virtual], overlap, target)

    active_occupied = tuple(
        index
        for index, value in zip(occupied, occupied_overlaps)
        if value >= threshold
    )
    active_virtual = tuple(
        index for index, value in zip(virtual, virtual_overlaps) if value >= threshold
    )
    singly_occupied = tuple(index for index in occupied if index in partial)
    if option == 3:
        # every singly occupied orbital stays in the active space, whether or not it
        # overlaps the target
        active_occupied = tuple(sorted(set(active_occupied) | set(singly_occupied)))

    n_electrons = 0.0
    for index in active_occupied:
        n_electrons += (
            1.0 if index in singly_occupied else float(export.mo_occupations[index])
        )
    for index in active_virtual:
        n_electrons += float(export.mo_occupations[index])
    n_orbitals = len(active_occupied) + len(active_virtual)

    checks = [
        f"target: {len(target)} l={angular} AO(s) of {export.atoms[centre]} at centre "
        f"{centre}, shell(s) {sorted({export.ao_labels[i].shell for i in target})}",
        f"occupied block: {len(occupied)} orbital(s) (occupation > {OCCUPIED_ABOVE}); "
        f"virtual block: {len(virtual)}",
    ]
    if partial:
        checks.append(
            f"the export carries {len(partial)} fractional occupation(s) "
            f"({', '.join(f'{export.mo_occupations[i]:.4f}' for i in partial[:6])}"
            + (", ..." if len(partial) > 6 else "")
            + "): this is a correlated export (CASSCF), so the occupied/virtual split "
            "is a convention, not the SCF aufbau. For a pure AVAS construction use an "
            "SCF export"
        )
    if singly_occupied:
        checks.append(
            f"singly occupied orbital(s) (occupation in 0.5-1.5): {list(singly_occupied)}"
        )
    verdicts = [
        f"occupied side: {len(active_occupied)} orbital(s) above the {threshold:g} "
        f"threshold; the largest overlap(s) "
        + " ".join(f"{value:.4f}" for value in occupied_overlaps[:4])
        + (" ..." if len(occupied_overlaps) > 4 else "")
        + " -- values near 1 are essentially the target AO itself, and low values among "
        "the kept ones name orbitals strongly mixed with their environment"
    ]
    above = [value for value in virtual_overlaps if value >= threshold]
    if above:
        verdicts.append(
            f"virtual side: {len(active_virtual)} orbital(s) above the threshold "
            f"(largest {max(virtual_overlaps):.4f}) -- these are the target-centred "
            "antibonding combinations that belong in the active space"
        )
    else:
        verdicts.append(
            f"virtual side: no orbital above the threshold (largest "
            f"{max(virtual_overlaps, default=0.0):.4f}) -- the antibonding target "
            "character is ligand-centred, so this target alone gives an active space "
            "without virtual orbitals (the source's [CuCl4]2- case). Next step: add the "
            "ligand AOs to the target set, or read the deficiency as 'AVAS with this "
            "target is not suited to the system'"
        )
    if option == 3 and singly_occupied:
        verdicts.append(
            "option 3 is applied: every singly occupied orbital stays in the active "
            "space regardless of its target overlap (the source's safeguard against the "
            "option-2 failure mode, in which an unoccupied beta orbital is forced into "
            "the core and CASCI can end up above the variational HF energy)"
        )
    verdicts.append(
        f"recommended active space: ({n_electrons:g} electrons, {n_orbitals} orbitals) = "
        f"{len(active_occupied)} occupied + {len(active_virtual)} virtual above the "
        "threshold (the electron count uses each active orbital's own occupation from "
        "the export, so a correlated export contributes its fractional occupations)"
    )
    verdicts.append(
        "both AVAS criteria are of the falsify-only kind (the source states it): a "
        "reading here that shows nothing wrong does not establish that the space is the "
        "best one -- the space-change SVD and the 'CASCI below the variational HF "
        "energy' check can still falsify it"
    )
    return AvasAnalysis(
        centre=centre,
        element=export.atoms[centre],
        angular=angular,
        shells=tuple(sorted({export.ao_labels[i].shell for i in target})),
        target_aos=target,
        threshold=threshold,
        option=option,
        occupied_indices=tuple(occupied),
        virtual_indices=tuple(virtual),
        partial_indices=tuple(partial),
        singly_occupied=singly_occupied,
        occupied_overlaps=occupied_overlaps,
        virtual_overlaps=virtual_overlaps,
        active_occupied=active_occupied,
        active_virtual=active_virtual,
        n_electrons=n_electrons,
        n_orbitals=n_orbitals,
        checks=tuple(checks),
        verdicts=tuple(verdicts),
    )


def _orca_block(analysis: AvasAnalysis) -> list[str]:
    """The measured ORCA ``%scf avas`` block for this target (one entry per component).

    The block grammar was measured on ORCA 6.1.1 (see fixtures/orca/README.md): each
    target function is one (shell, l, m_l, centre) entry, and the four lists run in
    parallel.  An empty list means the component spelling for this angular momentum
    is not recorded and no block is emitted.
    """
    components = component_labels(analysis.angular)
    if not components or not analysis.shells:
        return []
    shells = [shell for shell in analysis.shells for _ in components]
    m_values = [component for _ in analysis.shells for component in components]
    return [
        "%scf",
        "  avas",
        "    system",
        "      shell  " + ", ".join(str(shell) for shell in shells),
        "      l      " + ", ".join(str(analysis.angular) for _ in shells),
        "      m_l    " + ", ".join(m_values),
        "      center " + ", ".join(str(analysis.centre) for _ in shells),
        "    end",
        "  end",
        "end",
    ]


def run(analysis: AvasAnalysis) -> ReportSection:
    """Render the AVAS section."""
    lines = ["Target and blocks:"]
    for check in analysis.checks:
        lines.append(f"  - {check}")
    lines += [
        "",
        "Occupied-side target overlaps (descending):",
        "  " + " ".join(f"{value:.4f}" for value in analysis.occupied_overlaps),
        "Virtual-side target overlaps (descending):",
        "  " + " ".join(f"{value:.4f}" for value in analysis.virtual_overlaps),
        "",
        f"Recommendation ({analysis.threshold:g} truncation, option {analysis.option}): "
        f"({analysis.n_electrons:g} electrons, {analysis.n_orbitals} orbitals) -- "
        f"occupied {list(analysis.active_occupied)}, virtual "
        f"{list(analysis.active_virtual)} (0-based orbital indices of the export)",
        "",
        "Reading:",
    ]
    for verdict in analysis.verdicts:
        lines.append(f"  - {verdict}")
    block = _orca_block(analysis)
    if block:
        lines += [
            "",
            "For the same target, ORCA's built-in AVAS block (when its own minimal basis "
            "covers the element):",
            *[f"  {line}" for line in block],
        ]
    lines += [
        "",
        "Boundary of this route: the target set here is a subset of the calculation's own "
        "AOs (the source's non-minimal-ANO variant), which is what makes the projection "
        "computable from an export alone; the explicit list does not go through the "
        "minimal basis at all, and the ecosystem's minimal-basis coverage has no entry "
        "for Eu (measured through the engine's own localization tooling). And "
        "ORCA's CASSCF takes its active space by orbital order, not by an index list: the "
        "(n_el, n_orb) above sizes the input and the window, it does not name the "
        "orbitals to a restart.",
    ]
    return ReportSection(title="AVAS target projection (active-space construction)", body="\n".join(lines))


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the AVAS construction."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "AVAS projects a target AO set onto the SCF occupied and virtual blocks "
                "and rotates each block so that its overlap with the target is diagonal; "
                "the orbitals with non-zero overlap enter the active space, at most twice "
                "the target size. The truncation threshold (0.05-0.1) is the only "
                "user-chosen number besides the target and the SCF state. On open shells "
                "the source's option 2 (alpha orbitals only) can force an unoccupied "
                "beta orbital into the core -- CASCI can then lie above the variational "
                "HF energy; option 3 keeps every singly occupied orbital and is the safe "
                "choice for f-block targets. Both quality criteria of the source are "
                "falsify-only: the SVD overlap (0.65-0.99 measured) and the stability of "
                "properties against the space size."
            ),
            ref=(
                "Sayfutyarova E. R., Sun Q., Chan G. K.-L., Knizia G., J. Chem. Theory "
                "Comput., 2017, 13, 4063-4078, DOI 10.1021/acs.jctc.7b00128 (Eqs. 2-10; "
                "Sections 2.3-2.5, 4.1-4.3, 5.1)"
            ),
            bibkey="sayfutyarova2017automated",
            url="https://doi.org/10.1021/acs.jctc.7b00128",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "The projector identity is checked on the fixtures: for a target that is "
                "the whole AO set the occupied spectrum is exactly 1 in every entry, the "
                "eigenvalues stay in [0, 1], and the sum of the occupied-side "
                "eigenvalues equals the trace of the projected overlap (the target "
                "population of the occupied block). The engine's minimal-basis "
                "localization tooling has no entry for Eu (\"The minimal basis set is "
                "not defined for element Eu\", fixtures/orca/README.md), and the "
                "target-set route here avoids any minimal basis by construction, which "
                "is why it is the f-block route."
            ),
            ref="tests/test_avas.py; fixtures/orca/README.md",
        ),
    )
