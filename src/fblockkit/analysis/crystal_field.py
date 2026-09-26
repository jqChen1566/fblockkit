"""A4: crystal-field (CF) parameters ``B_k^q`` from sampled states -- a linear fit.

Model (reference paper, Eq. (3)): an effective operator inside the lowest
``|J M>`` manifold, written with extended Stevens operator equivalents,

    H_CF(J) = sum_{k=2,4,6} sum_{q=-k..k} B_k^q O_k^q(J) ,                 (2)
    E_i     = const + sum_{J M M'} C*_{J M,i} C_{J M',i} B_k^q [O_k^q]_{M M'}

so the sampled energies ``E_i`` are **linear** in the ``B_k^q``.  Any state
inside the manifold may be sampled -- it does not have to be an eigenstate;
that is the mechanism by which mean-field-cost states determine a
multireference-quality effective Hamiltonian.

This module is pure numpy: sample energies plus projection coefficients plus
the point group in, ``B_k^q`` plus a reconstructed spectrum out.  No engine
dependency, no file parsing.

Operator-equivalent source (do not re-derive from memory)
---------------------------------------------------------
Read on 2026-09-25 from EasySpin documentation, page "Stevens operators",
Table 1 "Extended Stevens operators O_k^q",
https://easyspin.org/documentation/stevensoperators.html
(equation image ``eqn/stevensoperators14.png``).  Transcribed rows (``+-`` is
the two signs of ``q``; the bare "1" in rows is the unit operator)::

    k=2  q=0     3 S_z^2 - s
    k=2  q=+-1   c_+- [S_z, S_+ +- S_-]_+
    k=2  q=+-2   c_+- (S_+^2 +- S_-^2)
    k=4  q=0     35 S_z^4 - (30 s - 25) S_z^2 + (3 s^2 - 6 s)
    k=4  q=+-1   c_+- [7 S_z^3 - (3 s + 1) S_z, S_+ +- S_-]_+
    k=4  q=+-2   c_+- [7 S_z^2 - (s + 5), S_+^2 +- S_-^2]_+
    k=4  q=+-3   c_+- [S_z, S_+^3 +- S_-^3]_+
    k=4  q=+-4   c_+- (S_+^4 +- S_-^4)
    k=6  q=0     231 S_z^6 - (315 s - 735) S_z^4 + (105 s^2 - 525 s + 294) S_z^2
                 - (5 s^3 - 40 s^2 + 60 s)
    k=6  q=+-1   c_+- [33 S_z^5 - (30 s - 15) S_z^3 + (5 s^2 - 10 s + 12) S_z,
                       S_+ +- S_-]_+
    k=6  q=+-2   c_+- [33 S_z^4 - (18 s + 123) S_z^2 + (s^2 + 10 s + 102),
                       S_+^2 +- S_-^2]_+
    k=6  q=+-3   c_+- [11 S_z^3 - (3 s + 59) S_z, S_+^3 +- S_-^3]_+
    k=6  q=+-4   c_+- [11 S_z^2 - (s + 38), S_+^4 +- S_-^4]_+
    k=6  q=+-5   c_+- [S_z, S_+^5 +- S_-^5]_+
    k=6  q=+-6   c_+- (S_+^6 +- S_-^6)

with, quoted from the same table's footnote, ``[A,B]_+ = (AB + BA)/2``,
``s = S(S+1)``, ``c_+ = 1/2``, ``c_- = 1/(2i)``.  The page is written for spin
operators; the operator equivalents used here are the same polynomials acting
in the ``|J M>`` manifold (that is, ``S -> J``, ``s = J(J+1)``).  EasySpin's
page states "each of the operators ... is hermitian, and the associated
coefficients ... are always real-valued" and warns that "conventions for the
normalization and the phase of these operators vary widely".

Convention stated explicitly
----------------------------
* Basis: ``|J M>`` with ``M = -J, -J+1, ..., +J`` (ascending; index ``m`` maps
  to ``M = m - J``).  ``J_+ |J M> = sqrt(J(J+1) - M(M+1)) |J M+1>``.
* Extended Stevens operators of the Rudowicz/Ryabov lineage -- the same
  lineage cited by the reference paper -- in the cosine/sine (tesseral) form:
  ``q > 0`` gives the cosine component, ``q < 0`` the sine component, and
  ``q = 0`` the axial one.
* Each ``O_k^q`` is **individually Hermitian**; ``O_k^q`` and ``O_k^{-q}`` are
  *independent* operators (they are not related by ``O_k^q^dagger =
  (-1)^q O_k^{-q}``).  The structural relation measured for this convention is

      [J_z, O_k^q] = i q O_k^{-q},

  which fixes the relative phase of the two signs (verified to 1e-15 relative
  in the tests).  Consequences: the fitted ``B_k^q`` are real, and a negative
  ``q`` is a separate parameter -- which is why the reference paper's Table S1
  lists e.g. both ``(4,3)`` and ``(4,-3)``.
* Units: the model is unit-agnostic; the reference paper's ``B_k^q`` and level
  energies are in cm^-1.
* Comparability: numbers computed under a different operator convention,
  manifold, unit or z-axis definition are **not** comparable.  Any output of
  this module should be reported together with the four items above (operator
  convention, projection manifold, units, z-axis).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from ..knowledge.models import (
    EVIDENCE_LITERATURE,
    EVIDENCE_MANUAL,
    EVIDENCE_MEASURED,
    Evidence,
)

#: Ranks kept: even ``k`` only (time-reversal symmetry) and ``k <= 6`` (4f
#: shell, ``l = 3``; the reference paper states higher ranks are conventionally
#: ignored).
K_RANKS = (2, 4, 6)

#: Condition-number above which the least-squares problem is reported as
#: ill-conditioned.  Measured discrimination on the reference fixture: a
#: well-posed sampling of the C3 problem gives cond ~ 7e5, while sampling only
#: Kramers-degenerate eigenstates gives cond ~ 2e18 (see the tests).
ILL_CONDITION_LIMIT = 1e8

#: Deviation of a sampled state's coefficient norm from 1 above which the fit is
#: flagged: Eq. (3) presumes normalised projection coefficients.
NORM_TOLERANCE = 1e-6

_ALLOWED_PARAMETERS: Mapping[str, tuple[tuple[int, int], ...]] = {
    # A true three-fold axis quantises |q| to multiples of 3 (q = 0, +-3, +-6).
    "C3": ((2, 0), (4, -3), (4, 0), (4, 3), (6, -6), (6, -3), (6, 0), (6, 3), (6, 6)),
    # Cubic with a four-fold axis as z; no second-rank term survives.
    "OH": ((4, 0), (4, 4), (6, 0), (6, 4)),
    # No symmetry at all: all 27 parameters are independent.
    "C1": tuple(
        (k, q) for k in K_RANKS for q in range(-k, k + 1)
    ),
}

#: Accepted spellings of the point groups (normalised: lower case, no spaces,
#: no underscores).
_POINT_GROUP_ALIASES: Mapping[str, str] = {
    "c3": "C3",
    "s6": "C3",
    "d3": "C3",
    "oh": "OH",
    "c1": "C1",
    "none": "C1",
    "nosymmetry": "C1",
}


class CrystalFieldError(ValueError):
    """Invalid crystal-field fitting input."""


_EVIDENCE_SYMMETRY = Evidence(
    kind=EVIDENCE_LITERATURE,
    text=(
        "Non-zero CF parameters are fixed by the point group: a true C3 axis leaves 9 "
        "parameters (q in {0, +-3, +-6}, even k), cubic symmetry with a four-fold axis "
        "leaves 4 -- (4,0), (4,4), (6,0), (6,4), with the whole second rank zero -- and a "
        "fit that needs all 27 means the frame carries no true symmetry (then the z-axis "
        "definition must be recorded). Table S1 of this paper lists the nine non-zero "
        "B_k^q of its C3 compound, Table S2 block A the four of its Oh compound."
    ),
    ref=(
        "Peng L., Liu S., Zhang X., Chen X., Li C., Ung S. F., Cheng H.-P., Chan G. K.-L., "
        "J. Phys. Chem. Lett., 2025, 16(47), 12312-12320, DOI 10.1021/acs.jpclett.5c02971"
    ),
    bibkey="peng2025accurate",
    url="https://doi.org/10.1021/acs.jpclett.5c02971",
)

_EVIDENCE_OPERATORS = Evidence(
    kind=EVIDENCE_MANUAL,
    text=(
        "Extended Stevens operator equivalents O_k^q (k = 2, 4, 6) taken from Table 1 "
        "'Extended Stevens operators O_k^q' of the EasySpin documentation page 'Stevens "
        "operators'; the page defines [A,B]_+ = (AB + BA)/2, s = S(S+1), c_+ = 1/2, "
        "c_- = 1/(2i), states that every operator of the table is Hermitian with real "
        "coefficients, and warns that normalisation and phase conventions vary widely. The "
        "same table lists O_k^{-q} separately from O_k^q (cosine/sine tesseral pair), which "
        "is the convention of the reference paper's Rudowicz/Ryabov lineage."
    ),
    ref="EasySpin documentation, 'Stevens operators', Table 1 (read 2026-09-25)",
    url="https://easyspin.org/documentation/stevensoperators.html",
)

_EVIDENCE_KRAMERS = Evidence(
    kind=EVIDENCE_MEASURED,
    text=(
        "Measured while validating this module (see tests/test_crystal_field.py): sampling "
        "only the eigenstates of a half-integer-J (Kramers) CF Hamiltonian gives a "
        "rank-deficient design matrix -- for the J = 15/2 C3 fixture the 16 eigenstates "
        "yield rank 8 of 10 columns and cond 2.4e18 (HF@HF column) / 5.5e17 (PBE0@HF), "
        "because the two members of each Kramers doublet have equal expectation values for "
        "every time-even operator, so each doublet contributes one equation instead of two; "
        "the two rows of a doublet agree to 1e-15 relative. Adding oriented "
        "(non-symmetry-adapted) sampled states restores rank 10, drops cond to 7e5, and "
        "recovers the input B_k^q to 3.6e-13 (HF@HF) / 6.6e-13 (PBE0@HF) relative. This is "
        "why the reference method samples randomly oriented determinants rather than "
        "eigenstates."
    ),
    ref="reproduced by tests/test_crystal_field.py (2026-09-25 fixture peng_s1_cf.json)",
)


def allowed_parameters(point_group: str) -> tuple[tuple[int, int], ...]:
    """``(k, q)`` pairs allowed by the point group, ordered by ``k`` then ``q``.

    ``"C3"`` (a true three-fold axis) -> 9, ``"Oh"`` (cubic) -> 4, ``"C1"`` /
    ``"none"`` (no symmetry) -> 27.  The tuple order defines the column order of
    the design matrix, so it is part of the fit output.
    """
    key = str(point_group).strip().lower().replace("_", "").replace(" ", "")
    canonical = _POINT_GROUP_ALIASES.get(key)
    if canonical is None:
        raise CrystalFieldError(
            f"unknown point group {point_group!r}: cannot tell which B_k^q are allowed. "
            f"Next step: pass one of 'C3' (three-fold axis), 'Oh' (cubic) or 'C1'/'none' "
            f"(no symmetry); for any other point group work out the allowed (k, q) from "
            f"the site symmetry and extend the table in this module."
        )
    return _ALLOWED_PARAMETERS[canonical]


def _check_J(J: float) -> int:
    """Validate ``J`` and return the manifold dimension ``2J + 1``."""
    try:
        value = float(J)
    except (TypeError, ValueError):
        raise CrystalFieldError(
            f"J must be a number, got {J!r}. Next step: pass the total angular momentum "
            f"of the manifold, e.g. 7.5 for Er3+ (4I15/2)."
        ) from None
    dimension = 2.0 * value + 1.0
    if abs(dimension - round(dimension)) > 1e-9 or value < 1.0:
        raise CrystalFieldError(
            f"J = {J!r} does not define a |J M> manifold of at least 3 states (2J+1 must be "
            f"a positive integer >= 3; the k = 2 operator equivalents vanish for J = 1/2). "
            f"Next step: pass a physical total angular momentum such as 6 (Tb3+, 7F6) or "
            f"7.5 (Er3+, 4I15/2)."
        )
    return int(round(dimension))


def _ladder_operators(J: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """``(M, J_z, J_+, J_-)`` in the ascending ``|J M>`` basis."""
    dimension = _check_J(J)
    m_values = np.arange(dimension, dtype=float) - float(J)
    j_z = np.diag(m_values).astype(complex)
    j_plus = np.zeros((dimension, dimension), dtype=complex)
    for index in range(dimension - 1):
        # J_+ |M> = sqrt(J(J+1) - M(M+1)) |M+1>; clip rounding noise at the top
        # of the ladder, where the radicand is exactly zero.
        radicand = J * (J + 1.0) - m_values[index] * (m_values[index] + 1.0)
        j_plus[index + 1, index] = math.sqrt(max(0.0, radicand))
    return m_values, j_z, j_plus, j_plus.conj().T


def _symmetrized(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """``[A, B]_+ = (AB + BA)/2`` as defined in the source table."""
    return (a @ b + b @ a) / 2.0


def _z_part(k: int, q_abs: int, s: float, j_z: np.ndarray, eye: np.ndarray) -> np.ndarray:
    """The ``J_z``-only factor of row ``(k, |q|)`` of the source table.

    Row by row from the transcription in the module docstring; ``q_abs == 0``
    rows are the whole operator and ``q_abs == k`` rows have factor ``1``.
    """
    j_z2 = j_z @ j_z
    if q_abs == 0:
        if k == 2:
            return 3.0 * j_z2 - s * eye
        if k == 4:
            return 35.0 * np.linalg.matrix_power(j_z, 4) - (30.0 * s - 25.0) * j_z2 + (
                3.0 * s * s - 6.0 * s
            ) * eye
        if k == 6:
            return (
                231.0 * np.linalg.matrix_power(j_z, 6)
                - (315.0 * s - 735.0) * np.linalg.matrix_power(j_z, 4)
                + (105.0 * s * s - 525.0 * s + 294.0) * j_z2
                - (5.0 * s**3 - 40.0 * s * s + 60.0 * s) * eye
            )
    if (k, q_abs) == (2, 1):
        return j_z
    if (k, q_abs) == (2, 2):
        return eye
    if (k, q_abs) == (4, 1):
        return 7.0 * np.linalg.matrix_power(j_z, 3) - (3.0 * s + 1.0) * j_z
    if (k, q_abs) == (4, 2):
        return 7.0 * j_z2 - (s + 5.0) * eye
    if (k, q_abs) == (4, 3):
        return j_z
    if (k, q_abs) == (4, 4):
        return eye
    if (k, q_abs) == (6, 1):
        return (
            33.0 * np.linalg.matrix_power(j_z, 5)
            - (30.0 * s - 15.0) * np.linalg.matrix_power(j_z, 3)
            + (5.0 * s * s - 10.0 * s + 12.0) * j_z
        )
    if (k, q_abs) == (6, 2):
        return (
            33.0 * np.linalg.matrix_power(j_z, 4)
            - (18.0 * s + 123.0) * j_z2
            + (s * s + 10.0 * s + 102.0) * eye
        )
    if (k, q_abs) == (6, 3):
        return 11.0 * np.linalg.matrix_power(j_z, 3) - (3.0 * s + 59.0) * j_z
    if (k, q_abs) == (6, 4):
        return 11.0 * j_z2 - (s + 38.0) * eye
    if (k, q_abs) == (6, 5):
        return j_z
    if (k, q_abs) == (6, 6):
        return eye
    raise CrystalFieldError(
        f"no operator equivalent implemented for (k, q) = ({k}, {q_abs}). "
        f"Next step: use k in {K_RANKS} with |q| <= k."
    )


def stevens_matrices(J: float) -> dict[tuple[int, int], np.ndarray]:
    """``[O_k^q(J)]_{M M'}`` for k in {2, 4, 6} and all q in ``-k..k``.

    Returns a dict keyed by ``(k, q)``; each value is a complex ``(2J+1,
    (2J+1))`` matrix in the ascending ``|J M>`` basis.  Cosine components
    (``q > 0``) are real symmetric, sine components (``q < 0``) purely
    imaginary antisymmetric, and all of them are Hermitian -- see the module
    docstring for the source and the convention.
    """
    _, j_z, j_plus, j_minus = _ladder_operators(J)
    eye = np.eye(j_z.shape[0], dtype=complex)
    s = float(J) * (float(J) + 1.0)

    matrices: dict[tuple[int, int], np.ndarray] = {}
    for k in K_RANKS:
        matrices[(k, 0)] = _z_part(k, 0, s, j_z, eye)
        for q_abs in range(1, k + 1):
            z_part = _z_part(k, q_abs, s, j_z, eye)
            up = np.linalg.matrix_power(j_plus, q_abs)
            down = np.linalg.matrix_power(j_minus, q_abs)
            matrices[(k, q_abs)] = 0.5 * _symmetrized(z_part, up + down)
            matrices[(k, -q_abs)] = (1.0 / 2.0j) * _symmetrized(z_part, up - down)
    return matrices


def hamiltonian(
    parameters: Mapping[tuple[int, int], float], J: float, *, const: float = 0.0
) -> np.ndarray:
    """``const * 1 + sum B_k^q O_k^q`` as a ``(2J+1, 2J+1)`` Hermitian matrix."""
    dimension = _check_J(J)
    matrices = stevens_matrices(J)
    total = const * np.eye(dimension, dtype=complex)
    for key, value in parameters.items():
        matrix = matrices.get((int(key[0]), int(key[1])))
        if matrix is None:
            raise CrystalFieldError(
                f"(k, q) = {tuple(key)} is outside this module's operator set. "
                f"Next step: use k in {K_RANKS} with |q| <= k."
            )
        total = total + float(value) * matrix
    return total


def _design(
    coefficients: Sequence[Sequence[complex]] | np.ndarray,
    columns: tuple[tuple[int, int], ...],
    J: float,
    *,
    with_const: bool,
) -> tuple[np.ndarray, float]:
    """Least-squares matrix of Eq. (3): one row per sampled state.

    Column ``(k, q)`` holds ``Re C_i^dagger O_k^q C_i``.  Also returns the
    largest imaginary part dropped in the process: for a Hermitian operator and
    a state in the declared basis it is zero up to rounding, so a large value
    means the coefficients are not in the documented ``|J M>`` basis.
    """
    dimension = _check_J(J)
    state = np.asarray(coefficients, dtype=complex)
    if state.ndim != 2 or state.shape[1] != dimension:
        raise CrystalFieldError(
            f"coefficient array has shape {state.shape}, expected (n_states, {dimension}) "
            f"for J = {J} (columns ordered as M = -J ... +J). "
            f"Next step: pass the projection coefficients C[i, m] of each sampled state onto "
            f"|J M>, with m = M + J; check the M ordering of the producing code."
        )
    if state.shape[0] == 0:
        raise CrystalFieldError(
            "no sampled states given (0 rows). Next step: pass at least as many sampled "
            "states as there are allowed B_k^q plus one for the constant."
        )
    if not np.all(np.isfinite(state)):
        raise CrystalFieldError(
            "coefficients contain non-finite values. Next step: drop or fix the failed "
            "projected states before fitting (a state whose projection did not converge is "
            "not a sample)."
        )
    matrices = stevens_matrices(J)
    block = np.empty((state.shape[0], len(columns)), dtype=float)
    imaginary = 0.0
    for position, key in enumerate(columns):
        value = np.einsum("im,mn,in->i", state.conj(), matrices[key], state)
        block[:, position] = value.real
        imaginary = max(imaginary, float(np.abs(value.imag).max()))
    if with_const:
        block = np.column_stack([block, np.ones(state.shape[0])])
    return block, imaginary


def design_matrix(
    coefficients: Sequence[Sequence[complex]] | np.ndarray,
    point_group: str,
    J: float,
    *,
    with_const: bool = True,
) -> tuple[np.ndarray, tuple[tuple[int, int], ...]]:
    """Eq. (3) matrix, one row per sampled state, and its column keys.

    Columns are the allowed ``(k, q)`` pairs in :func:`allowed_parameters`
    order, plus a trailing all-ones column when ``with_const`` is set.  The
    least-squares matrix itself is what ``numpy.linalg.cond`` and
    ``numpy.linalg.matrix_rank`` should be applied to.
    """
    columns = allowed_parameters(point_group)
    matrix, _ = _design(coefficients, columns, J, with_const=with_const)
    return matrix, columns


@dataclass(frozen=True)
class CFResult:
    """Outcome of a CF fit (all energies and parameters in the caller's units).

    ``reconstructed_levels[i]`` is the model energy of sampled state ``i``
    (``const + sum B_k^q C_i^dagger O_k^q C_i``) and ``residuals[i]`` is
    ``levels[i] - reconstructed_levels[i]``, both per sampled state -- *not* the
    level structure of the fitted Hamiltonian, which is :meth:`spectrum`.
    Sequence fields are tuples so the record stays immutable.
    """

    parameters: Mapping[tuple[int, int], float]
    const: float
    reconstructed_levels: tuple[float, ...]
    residuals: tuple[float, ...]
    condition_number: float
    rank: int
    n_states: int
    n_parameters: int
    underdetermined: bool
    rank_deficient: bool
    ill_conditioned: bool
    notes: tuple[str, ...]
    j: float
    point_group: str

    @property
    def max_abs_residual(self) -> float:
        return max((abs(value) for value in self.residuals), default=0.0)

    def spectrum(self, *, ascending: bool = True) -> np.ndarray:
        """Eigenvalues of the fitted ``H_CF`` (paper step 6).

        For a fit on sampled eigenstates these reproduce
        :attr:`reconstructed_levels`; for oriented (non-eigenstate) sampling
        they need not -- the reconstruction is per state, the spectrum is the
        level structure a measurement or a multireference calculation would be
        compared against.
        """
        values = np.linalg.eigvalsh(hamiltonian(self.parameters, self.j, const=self.const))
        return np.sort(values) if ascending else values


def fit_crystal_field(
    levels: Sequence[float] | np.ndarray,
    coefficients: Sequence[Sequence[complex]] | np.ndarray,
    point_group: str,
    J: float,
    *,
    with_const: bool = True,
) -> CFResult:
    """Fit ``B_k^q`` and the constant from sampled energies and coefficients.

    ``levels[i]`` is the energy of sampled state ``i`` and ``coefficients[i, m]``
    its amplitude on ``|J M>`` with ``M = m - J``, i.e. the state is
    ``sum_m coefficients[i, m] |J M>`` and the row is ``C_i^dagger O C_i``.
    Pass the amplitudes themselves, not their complex conjugates: conjugating
    the coordinates conjugates the operator expectations, which flips the sign
    of every ``q < 0`` parameter while leaving the spectrum unchanged, so such a
    fit reaches a small residual with wrong sine-type parameters.
    Coefficients are assumed normalised (Eq. (3) presumes it, and a deviation is
    reported in the notes).  The solution is ``numpy.linalg.lstsq`` over the
    symmetry-allowed columns; diagnostics flag a rank-deficient, underdetermined
    or ill-conditioned system rather than returning silently meaningless
    parameters.
    """
    energy = np.asarray(levels, dtype=float).ravel()
    columns = allowed_parameters(point_group)
    matrix, imaginary = _design(coefficients, columns, J, with_const=with_const)
    if matrix.shape[0] != energy.size:
        raise CrystalFieldError(
            f"got {energy.size} level energies but {matrix.shape[0]} sampled states. "
            f"Next step: pass one energy per state, in the same order as the coefficient "
            f"rows."
        )
    if not np.all(np.isfinite(energy)):
        raise CrystalFieldError(
            "level energies contain non-finite values. Next step: drop or fix the failed "
            "sampled states before fitting (a non-converged mean-field state is not a "
            "sample)."
        )

    solution, _, rank, _ = np.linalg.lstsq(matrix, energy, rcond=None)
    parameters = {key: float(solution[index]) for index, key in enumerate(columns)}
    const = float(solution[-1]) if with_const else 0.0

    predicted = matrix @ solution
    with np.errstate(divide="ignore", invalid="ignore"):
        condition = float(np.linalg.cond(matrix))
    n_columns = matrix.shape[1]
    underdetermined = matrix.shape[0] < n_columns
    rank_deficient = int(rank) < n_columns
    ill_conditioned = not np.isfinite(condition) or condition > ILL_CONDITION_LIMIT

    notes: list[str] = []
    if underdetermined:
        notes.append(
            f"underdetermined: {matrix.shape[0]} sampled states < {n_columns} fitted "
            f"parameters -- the least-squares solution is one of many; only the "
            f"reconstructed energies carry information."
        )
    if rank_deficient:
        notes.append(
            f"rank deficient: only {int(rank)} of {n_columns} columns are independent. "
            f"Degenerate partners - the two members of a Kramers doublet for half-integer J, "
            f"or a time-reversal-related +-M pair for integer J - have identical expectation "
            f"values for every time-even operator, so sampling eigenstates gives one "
            f"independent equation per degenerate pair instead of two. Supply sampled "
            f"(non-eigen-) states to restore full rank."
        )
    if ill_conditioned:
        notes.append(
            f"ill-conditioned design matrix (cond = {condition:.3e} > "
            f"{ILL_CONDITION_LIMIT:.0e}): small energy errors can move B_k^q a long way. "
            f"Sample more (and differently oriented) states before trusting the parameters."
        )
    norms = np.linalg.norm(np.asarray(coefficients, dtype=complex), axis=1)
    worst_norm = float(np.abs(norms - 1.0).max()) if norms.size else 0.0
    if worst_norm > NORM_TOLERANCE:
        notes.append(
            f"coefficient norms deviate from 1 by up to {worst_norm:.3e}: Eq. (3) presumes "
            f"normalised projection coefficients, so the fitted B_k^q carry this bias."
        )
    if imaginary > 1e-8 * max(1.0, float(np.abs(matrix).max())):
        notes.append(
            f"discarded imaginary parts up to {imaginary:.3e}: the expectation values should "
            f"be real for states in the documented |J M> basis. Check the M ordering and the "
            f"normalisation of the coefficients before reading the parameters."
        )
    if not notes:
        notes.append(
            "well posed: full column rank and cond below the ill-conditioning limit."
        )

    return CFResult(
        parameters=parameters,
        const=const,
        reconstructed_levels=tuple(float(value) for value in predicted),
        residuals=tuple(float(value) for value in (energy - predicted)),
        condition_number=condition,
        rank=int(rank),
        n_states=int(matrix.shape[0]),
        n_parameters=len(columns),
        underdetermined=underdetermined,
        rank_deficient=rank_deficient,
        ill_conditioned=ill_conditioned,
        notes=tuple(notes),
        j=float(J),
        point_group=str(point_group),
    )


def evidence() -> tuple[Evidence, ...]:
    """Provenance of this module's rules (see architecture design v0.1 section 2)."""
    return (_EVIDENCE_SYMMETRY, _EVIDENCE_OPERATORS, _EVIDENCE_KRAMERS)
