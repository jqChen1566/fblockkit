"""Atomic-term check: the effective L, S and J of a computed state, against Hund's rules.

The check comes from the crystal-field workflow of Peng et al. (JPCL 2025): before
any mean-field density is used to fit a crystal field, its effective quantum
numbers must belong to the Hund ground term of the ion -- they report J = 7.5
(the Hund value) only for the Hartree-Fock density, while M06/PBE0/B3LYP/PBE/
r2SCAN give 7.3 down to 5.9, i.e. a density that must not be used.  The rule is
stated there as "J = L + S (the maximum J)": in a scalar (spin-orbit-free)
calculation the computed component of a Hund term is the stretched one, whose
<J^2> works out to J = L + S exactly; the physical Hund third-rule J (|L - S|
for a less-than-half-filled shell) is only reached once spin-orbit coupling is
included, and both numbers are reported here so the distinction is explicit.

How the numbers are obtained: the active-space CI state is exact in this toolkit
(``analysis.entropy_rdm``), so the expectations of L^2, S^2 and L.S are computed
from its 1- and 2-RDMs together with the angular-momentum operators expressed in
the active-orbital basis.  Those operators are built by projecting the active
orbitals onto the f AOs of the chosen centre (the measured structure of an ORCA
export: the f block is ``(shell, component)`` ordered and its overlap is exactly
``R (x) I``, so the orthonormal f basis and the pure angular matrices separate --
see ``analysis.solid_harmonics``).  The f character of each active orbital is
reported as the quality measure of that projection, and the Hund verdict is only
issued when the projection captures the active space essentially completely.

Scope notes (each one measured or from the source):

- the check describes the *computed* state: a state outside the Hund term is
  reported as such, and a state whose active space is not one clean l shell at all
  (the f-block "wrong-solution" failure mode in which the electrons sit in d/s/p
  orbitals and the f shell is empty) is reported with its measured l characters
  instead of quantum numbers;
- <S^2> computed here from the densities is cross-checked against the value the
  CI solver already reports from its own S^2 matrix (a second, independent code
  path), and the check is part of the report.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..parsers.orca_json import ANGULAR_LETTERS
from ..knowledge.models import (
    EVIDENCE_LITERATURE,
    EVIDENCE_MEASURED,
    Evidence,
    ReportSection,
)
from .entropy_rdm import (
    EntropyRdmError,
    SpinDensities,
    expectation_product,
    spin_orbital_matrix,
)
from .solid_harmonics import angular_momentum_matrices, component_labels

__all__ = [
    "AtomicTermError",
    "ShellProjection",
    "AtomicTermAnalysis",
    "HundTerm",
    "project_shell",
    "analyze_terms",
    "hund_term",
    "run",
    "evidence",
]


class AtomicTermError(ValueError):
    """The atomic-term check cannot run on the given data (with a next step)."""


#: |<f|f> - 1| tolerance when verifying the measured (shell, component) block structure.
_BLOCK_TOLERANCE = 1e-8

#: Provisional quality line: below this mean f character the Hund verdict is withheld.
CHARACTER_LINE = 0.9

#: Hund's-rule L for a shell of l = 2 / 3 with n electrons (n = 1..4l+1; mirrored
#: for more than half filling).  Standard textbook values, listed per l because the
#: check is used for d and f shells alike.
_HUND_L = {
    2: {1: 2, 2: 3, 3: 3, 4: 2, 5: 0},
    3: {1: 3, 2: 5, 3: 6, 4: 6, 5: 5, 6: 3, 7: 0},
}


@dataclass(frozen=True)
class HundTerm:
    """The Hund's-rule ground term of an l^n configuration."""

    angular: int
    n_electrons: int
    s: float
    l: float
    j: float
    j_maximum: float  # L + S, the value a scalar stretched component carries
    half_filled_or_less: bool


def hund_term(angular: int, n_electrons: int) -> HundTerm:
    """Hund's rules for ``l^n``: S = min(n, 4l+2-n)/2, then the tabulated L and J."""
    capacity = 2 * (2 * angular + 1)
    if not 1 <= n_electrons <= capacity - 1:
        raise AtomicTermError(
            f"{n_electrons} electrons do not form a Hund term of an l={angular} shell "
            f"(capacity {capacity}). Next step: check the shell occupation."
        )
    table = _HUND_L.get(angular)
    if table is None:
        raise AtomicTermError(
            f"no Hund-rule L table for l={angular}. Next step: add it (the values are "
            "textbook; the pattern mirrors about the half-filled shell)."
        )
    less = n_electrons <= 2 * angular + 1
    s = min(n_electrons, capacity - n_electrons) / 2.0
    l_value = table[min(n_electrons, capacity - n_electrons)]
    j = abs(l_value - s) if less else l_value + s
    return HundTerm(
        angular=angular,
        n_electrons=n_electrons,
        s=s,
        l=float(l_value),
        j=j,
        j_maximum=float(l_value) + s,
        half_filled_or_less=less,
    )


@dataclass(frozen=True)
class ShellProjection:
    """The active orbitals' coefficients in the orthonormal l-shell basis of one centre.

    ``amplitudes[i]`` is active orbital ``i``'s Löwdin amplitude vector over the
    orthonormal shell basis, flat with the (shell, component) index (shell
    major); ``characters[i] = sum |amplitude|^2`` is that orbital's l character
    (a projection norm in [0, 1]); ``radial_overlap`` is the measured shell x
    shell overlap R.  The amplitudes are kept *per shell* rather than summed
    over shells: both the character (a sum of squares) and the operator matrices
    (``L = sum_s d_s^T A d_s``, because the operator is diagonal in the shell
    index) need the shell structure.
    """

    center: int
    element: str
    angular: int
    n_shells: int
    component_labels: tuple[str, ...]
    amplitudes: tuple[tuple[float, ...], ...]
    characters: tuple[float, ...]
    radial_overlap: tuple[tuple[float, ...], ...]

    @property
    def n_components(self) -> int:
        return 2 * self.angular + 1


def _component_order(labels, angular: int, center: int) -> tuple[list[int], list[int], str]:
    """Index and shell number of the AOs of ``(center, angular)`` in component order."""
    wanted = ANGULAR_LETTERS.get(angular)
    if wanted is None:
        raise AtomicTermError(
            f"angular momentum l={angular} has no ORCA label letter recorded. "
            "Next step: extend ANGULAR_LETTERS in the export reader."
        )
    components = component_labels(angular)
    selected = [
        (index, label)
        for index, label in enumerate(labels)
        if label.center == center and label.angular == wanted
    ]
    if not selected:
        raise AtomicTermError(
            f"the export has no l={angular} ({wanted}) AO on centre {center}. Next step: "
            "check the centre index and the shell you are analysing."
        )
    if len(selected) % len(components):
        raise AtomicTermError(
            f"centre {center} has {len(selected)} l={angular} AOs, not a multiple of "
            f"{len(components)}. Next step: check the AO labelling."
        )
    order: list[int] = []
    shells: list[int] = []
    element = ""
    for shell_number in sorted({label.shell for _, label in selected}):
        for component in components:
            matches = [
                index
                for index, label in selected
                if label.shell == shell_number and label.component == component
            ]
            if len(matches) != 1:
                raise AtomicTermError(
                    f"the l={angular} shell {shell_number} of centre {center} has "
                    f"{len(matches)} AOs labelled {component!r}, expected exactly one. "
                    "Next step: the AO labelling does not follow the measured mean "
                    "component spelling -- report the labels so the reader can be extended."
                )
            order.append(matches[0])
            shells.append(shell_number)
            element = selected[0][1].element
    return order, shells, element


def project_shell(
    coefficients: np.ndarray,
    overlap: np.ndarray,
    labels,
    *,
    center: int,
    angular: int = 3,
) -> ShellProjection:
    """Project active orbitals onto the orthonormal l-shell basis of one centre.

    ``coefficients`` is an AO coefficient matrix with the active orbitals in the
    columns; ``overlap`` the AO overlap matrix; ``labels`` the parsed AO labels
    of the same export.  The measured structure (the l block is ordered by
    (shell, component) and its overlap is ``R (x) I``) is verified, not assumed.
    """
    coef = np.asarray(coefficients, dtype=np.float64)
    s = np.asarray(overlap, dtype=np.float64)
    if coef.ndim != 2 or coef.shape[0] != s.shape[0]:
        raise AtomicTermError(
            f"the coefficient block is {coef.shape} against an overlap matrix of "
            f"{s.shape}. Next step: pass the AO coefficient matrix (orbitals in the "
            "columns) of the same export as the overlap."
        )
    order, shells, element = _component_order(labels, angular, center)
    n_components = 2 * angular + 1
    n_shells = len(order) // n_components
    block = s[np.ix_(order, order)]
    # measured structure: block[(s1, c1), (s2, c2)] = R[s1, s2] * delta(c1, c2), i.e.
    # every shell pair carries a multiple of the component identity (the radial
    # overlap between different shells is NOT zero -- only the cross-component
    # entries are)
    radial = np.zeros((n_shells, n_shells))
    for first in range(n_shells):
        rows = slice(first * n_components, (first + 1) * n_components)
        for second in range(n_shells):
            cols = slice(second * n_components, (second + 1) * n_components)
            piece = block[rows, cols]
            target = piece[0, 0] * np.eye(n_components)
            deviation = float(np.abs(piece - target).max())
            if deviation > _BLOCK_TOLERANCE:
                raise AtomicTermError(
                    f"the l={angular} block of centre {center} is not R (x) I "
                    f"(deviation {deviation:.2e} in the shell pair ({first}, {second})). "
                    "Next step: the AO order or normalisation of this export differs from "
                    "the measured convention -- report it."
                )
            radial[first, second] = piece[0, 0]
    if np.abs(radial - radial.T).max() > _BLOCK_TOLERANCE:
        raise AtomicTermError(
            f"the radial overlap matrix of the l={angular} shells of centre {center} is "
            "not symmetric. Next step: report the export -- this is not the measured "
            "convention."
        )
    values, vectors = np.linalg.eigh(radial)
    if values.min() <= 1e-10:
        raise AtomicTermError(
            f"the l={angular} shells of centre {center} are linearly dependent "
            f"(smallest radial eigenvalue {values.min():.3e}). Next step: check the "
            "basis set of the export."
        )
    root = vectors @ np.diag(values**0.5) @ vectors.T
    # rows of `chosen` are (shell, component) ordered; the Loewdin amplitudes are
    # d[(s, c), i] = sum_s' sqrt(R)[s, s'] C[(s', c), i] -- exactly S_f^{1/2} C_f
    chosen = coef[order, :].reshape(n_shells, n_components, -1)
    amplitudes = np.einsum("st,tci->sci", root, chosen, optimize=True)
    flat = amplitudes.reshape(n_shells * n_components, -1).T  # (n_active, shell*component)
    characters = tuple(float(value) for value in np.einsum("ik,ik->i", flat, flat))
    return ShellProjection(
        center=center,
        element=element,
        angular=angular,
        n_shells=n_shells,
        component_labels=component_labels(angular),
        amplitudes=tuple(tuple(float(v) for v in row) for row in flat),
        characters=characters,
        radial_overlap=tuple(tuple(float(v) for v in row) for row in radial),
    )


@dataclass(frozen=True)
class AtomicTermAnalysis:
    """Effective quantum numbers of one computed state in one projected shell."""

    projection: ShellProjection
    l2: float
    s2: float
    ls: float
    j2: float
    l_eff: float
    s_eff: float
    j_eff: float
    n_electrons: float
    hund: HundTerm | None
    checks: tuple[str, ...]
    verdicts: tuple[str, ...]


def _effective(j2: float) -> float:
    """The effective quantum number from ``<J^2> = J (J+1)``."""
    return (-1.0 + (1.0 + 4.0 * max(j2, 0.0)) ** 0.5) / 2.0


def _clean(value: float) -> float:
    """Normalise a signed zero for clean reporting (-0.0 + 0.0 is +0.0)."""
    return float(value) + 0.0


def _spin_matrices(norb: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The spin-orbital matrices of S_x, S_y, S_z in the active basis (2n x 2n)."""
    spatial = np.eye(norb)
    half = 0.5
    sx = np.zeros((2 * norb, 2 * norb), dtype=np.complex128)
    sy = np.zeros_like(sx)
    sz = np.zeros_like(sx)
    sx[:norb, norb:] = half * spatial
    sx[norb:, :norb] = half * spatial
    sy[:norb, norb:] = -0.5j * spatial
    sy[norb:, :norb] = 0.5j * spatial
    sz[:norb, :norb] = half * spatial
    sz[norb:, norb:] = -half * spatial
    return sx, sy, sz


def analyze_terms(
    densities: SpinDensities,
    coefficients: np.ndarray,
    overlap: np.ndarray,
    labels,
    *,
    center: int,
    angular: int = 3,
    s2_reference: float | None = None,
) -> AtomicTermAnalysis:
    """Compute <L^2>, <S^2>, <L.S>, <J^2> and the effective quantum numbers.

    ``s2_reference`` is the value the CI solver reports from its own S^2 matrix;
    the cross-check against the density-based value is recorded in ``checks``.
    """
    projection = project_shell(
        coefficients, overlap, labels, center=center, angular=angular
    )
    norb = densities.norb
    if len(projection.amplitudes) != norb:
        raise AtomicTermError(
            f"the coefficient block has {len(projection.amplitudes)} active orbitals but the "
            f"state has {norb}. Next step: pass the same active window that produced the "
            "state."
        )
    # (n_shells, n_components, n_active): L = sum_s d_s^T A d_s because the operator
    # is diagonal in the shell index
    n_shells = projection.n_shells
    n_components = projection.n_components
    shells = (
        np.array(projection.amplitudes)
        .reshape(norb, n_shells, n_components)
        .transpose(1, 2, 0)
    )
    l_matrices = angular_momentum_matrices(angular)
    l_effective_series = [
        sum(block.conj().T @ matrix @ block for block in shells) for matrix in l_matrices
    ]
    spin_matrices = _spin_matrices(norb)

    def spin_free(matrix):
        return spin_orbital_matrix(matrix, norb)

    l_squared = sum(
        expectation_product(spin_free(matrix), spin_free(matrix), densities).real
        for matrix in l_effective_series
    )
    s_squared = sum(
        expectation_product(matrix, matrix, densities).real for matrix in spin_matrices
    )
    l_dot_s = sum(
        expectation_product(spin_free(l_mat), s_mat, densities).real
        for l_mat, s_mat in zip(l_effective_series, spin_matrices)
    )
    j_squared = l_squared + s_squared + 2.0 * l_dot_s

    occupation = np.diag(densities.gamma_a + densities.gamma_b)
    n_electrons = float(sum(c * o for c, o in zip(projection.characters, occupation)))
    checks = []
    if s2_reference is not None:
        deviation = abs(s_squared - s2_reference)
        checks.append(
            f"<S^2> from the density objects = {s_squared:.10f} vs the CI solver's own "
            f"S^2 matrix value {s2_reference:.10f} (deviation {deviation:.2e})"
        )
        if deviation > 1e-6:
            raise AtomicTermError(
                f"the two independent <S^2> values disagree by {deviation:.2e}: the "
                "density objects and the CI state are inconsistent. Next step: rebuild "
                "the densities from this state."
            )
    mean_character = float(np.mean(projection.characters))
    checks.append(
        f"active-orbital l={angular} character: mean {mean_character:.4f}, minimum "
        f"{min(projection.characters):.4f} over {norb} orbital(s)"
    )
    verdicts = []
    hund: HundTerm | None = None
    if mean_character < CHARACTER_LINE:
        verdicts.append(
            f"the active space is not a clean l={angular} shell (mean character "
            f"{mean_character:.3f} is below the {CHARACTER_LINE} line): the quantum "
            "numbers are not formed and the Hund verdict is withheld. For an f-block "
            "target this is the signature of a wrong solution in which the f shell is "
            "not the valence shell of the calculation"
        )
    else:
        rounded = int(round(n_electrons))
        if abs(n_electrons - rounded) > 0.15:
            verdicts.append(
                f"the projected shell holds {n_electrons:.3f} electrons, which does not "
                "name a configuration: no Hund verdict"
            )
        else:
            hund = hund_term(angular, rounded)
            verdicts.append(
                f"Hund's rules for l={angular}, n={rounded}: S = {hund.s:g}, "
                f"L = {hund.l:g}, J = {hund.j:g} "
                f"({'|L-S|' if hund.half_filled_or_less else 'L+S'})"
            )
            for name, value, target in (
                ("S", _effective(s_squared), hund.s),
                ("L", _effective(l_squared), hund.l),
            ):
                if abs(value - target) > 0.1:
                    verdicts.append(
                        f"NOT in the Hund term: {name}_eff = {value:.3f} against the "
                        f"Hund value {target:g}"
                    )
            if all(
                abs(value - target) <= 0.1
                for value, target in (
                    (_effective(s_squared), hund.s),
                    (_effective(l_squared), hund.l),
                )
            ):
                j_eff = _effective(j_squared)
                # J is only meaningful for a J eigenstate.  A scalar (spin-orbit-free)
                # state is one (M_L, M_S) component: the stretched one gives J = L+S,
                # the fully antiparallel one gives |L-S|, and any other component has
                # no definite J at all -- none of those three is a defect, which is
                # why the source's discriminating check is the L and S comparison
                # above (their r2SCAN entry moves L to 4.6, not J alone).
                if abs(j_eff - hund.j_maximum) <= 0.1:
                    verdicts.append(
                        f"the computed component is the maximal-J one: J_eff = "
                        f"{j_eff:.3f} = L + S = {hund.j_maximum:g} (the value a scalar "
                        f"stretched density carries; the Hund third-rule J = {hund.j:g} "
                        "is reached after spin-orbit coupling)"
                    )
                elif abs(j_eff - abs(hund.l - hund.s)) <= 0.1:
                    verdicts.append(
                        f"the computed component is the minimal-J one: J_eff = "
                        f"{j_eff:.3f} = |L - S| "
                        f"({'the Hund third-rule value' if hund.j == abs(hund.l - hund.s) else 'a higher multiplet'})"
                    )
                else:
                    verdicts.append(
                        f"J_eff = {j_eff:.3f} is a non-stretched component of the Hund "
                        "term (M_J is neither maximal nor minimal): in a scalar "
                        "calculation J is not yet defined -- spin-orbit coupling decides "
                        "the multiplet, and this is not a defect"
                    )
    return AtomicTermAnalysis(
        projection=projection,
        l2=l_squared,
        s2=s_squared,
        ls=l_dot_s,
        j2=j_squared,
        l_eff=_effective(l_squared),
        s_eff=_effective(s_squared),
        j_eff=_effective(j_squared),
        n_electrons=n_electrons,
        hund=hund,
        checks=tuple(checks),
        verdicts=tuple(verdicts),
    )


def analyze_export(
    densities: SpinDensities,
    s2_reference: float,
    export,
    window,
    *,
    centres: tuple[int, ...] | None = None,
) -> tuple[ReportSection, ...]:
    """One atomic-term section per f-block centre of an orbital export.

    ``centres=None`` auto-detects the lanthanide/actinide atoms of the export
    (the same element test the d-solution diagnosis uses).  A centre whose shell
    cannot be projected (no such AOs in the basis) yields an abstention section
    with the reason, never a silently missing section.
    """
    from ..knowledge.elements import is_f_element

    if export.ao_labels is None or export.overlap is None:
        raise AtomicTermError(
            "the export carries no AO labels or no overlap matrix, so no centre can be "
            "projected. Next step: re-export with the S-Matrix and the orbital labels "
            "(the defaults of orca_2json)."
        )
    if centres is None:
        centres = tuple(
            index for index, symbol in enumerate(export.atoms) if is_f_element(symbol)
        )
    if not centres:
        return ()
    coefficients = np.array(export.mo_coefficients).T[:, list(window)]
    overlap = np.array(export.overlap)
    sections = []
    for centre in centres:
        try:
            result = analyze_terms(
                densities,
                coefficients,
                overlap,
                export.ao_labels,
                center=centre,
                s2_reference=s2_reference,
            )
        except AtomicTermError as exc:
            sections.append(
                ReportSection(
                    title="Atomic-term check (L, S, J against Hund's rules)",
                    body=(
                        f"Centre {centre} ({export.atoms[centre]}): not checked. "
                        f"{exc}"
                    ),
                )
            )
            continue
        sections.append(run(result))
    return tuple(sections)


def run(analysis: AtomicTermAnalysis, *, method_note: str = "") -> ReportSection:
    """Render the atomic-term section of the exact-state report."""
    projection = analysis.projection
    lines = [
        f"Target: the l={projection.angular} shell "
        f"({projection.component_labels[0]!r}..{projection.component_labels[-1]!r}) of "
        f"{projection.element} at centre {projection.center}; {projection.n_shells} "
        "contained shell(s) in the basis",
        "",
        "Quality of the projection (the check only describes the projected part):",
    ]
    for check in analysis.checks:
        lines.append(f"  - {check}")
    lines += [
        "",
        f"<L^2> = {_clean(analysis.l2):.6f}  -> L_eff = {_clean(analysis.l_eff):.4f}",
        f"<S^2> = {_clean(analysis.s2):.6f}  -> S_eff = {_clean(analysis.s_eff):.4f}",
        f"<L.S> = {_clean(analysis.ls):.6f}",
        f"<J^2> = {_clean(analysis.j2):.6f}  -> J_eff = {_clean(analysis.j_eff):.4f}",
        f"projected shell occupation = {_clean(analysis.n_electrons):.4f} electron(s)",
        "",
        "Verdict:",
    ]
    for verdict in analysis.verdicts:
        lines.append(f"  - {verdict}")
    if method_note:
        lines += ["", method_note]
    return ReportSection(title="Atomic-term check (L, S, J against Hund's rules)", body="\n".join(lines))


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the atomic-term check."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "Before a mean-field density is used for crystal-field fitting, its "
                "effective quantum numbers must belong to the Hund ground term: the "
                "source reports J = 7.5 (the Hund value for Dy3+, 6H15/2) for the HF "
                "density against 7.3/7.1/7.1/7.1/6.7/5.9 for M06/PBE0/B97-D/B3LYP/PBE/"
                "r2SCAN, with the rule stated as J = L + S (the maximum J)."
            ),
            ref=(
                "Peng et al., J. Phys. Chem. Lett., 2025, 16, 12312-12320, "
                "DOI 10.1021/acs.jpclett.5c02971 (arXiv:2505.16905v2), Table 1"
            ),
            bibkey="peng2025accurate",
            url="https://doi.org/10.1021/acs.jpclett.5c02971",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "The <S^2> used here is computed from the spin-resolved density objects "
                "with the two-operator expectation formula, an independent code path "
                "from the CI solver's own S^2 matrix; on the N2 (singlet) and Eu3+ "
                "(septet) fixtures the two values agree to 1e-15 and 0.0 respectively. "
                "The angular-momentum matrices themselves are derived from the solid "
                "harmonic polynomials and pinned by [L_i, L_j] = i e_ijk L_k, "
                "L^2 = l(l+1) and the m = -l..l spectrum (tests/test_solid_harmonics.py)."
            ),
            ref="tests/test_solid_harmonics.py, tests/test_atomic_terms.py",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "The projection convention was measured on the Eu export: the f block of "
                "an ORCA basis is ordered by (shell, component) with the component "
                "spelling 0, +1, -1, +2, -2, +3, -3, and its overlap matrix is exactly "
                "R (x) I with the cross-component blocks at 2e-16 -- so the orthonormal "
                "shell basis separates into radial and angular parts."
            ),
            ref="measured 2026-09-26 on the Eu (SARC2-DKH-QZVP) export; fixtures/orca/README.md",
        ),
    )
