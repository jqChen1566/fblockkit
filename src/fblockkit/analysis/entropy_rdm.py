"""Exact four-state single-orbital entropy, rebuilt from an FCIDUMP.

Why this module exists
----------------------
The single-orbital entropy of the autoCAS protocol (Stein & Reiher; see the A2
module :mod:`fblockkit.analysis.entropy` for the definition and the references)
is fixed by three numbers per orbital -- the alpha occupation ``n_a``, the beta
occupation ``n_b`` and the double-occupancy probability ``p2``::

    w = [1 - n_a - n_b + p2,  n_a - p2,  n_b - p2,  p2]
    s = -sum w_j ln w_j

An ORCA output prints only the *spin-summed* natural occupations, so an output
file alone can never produce this quantity (that is why A2 reports a rigorous
*upper bound* instead).  ORCA's own 2-RDM export (``orca_2json`` RDM2_aa /
RDM2_ab / RDM2_bb from ``.RDM2`` files) is the documented route, but on ORCA
6.1.1 it is unreachable for CAS-type methods: requesting the 2-body density
(``%autoci Density2 ...``) either does nothing (unrelaxed), aborts with
"sigma equations for the chosen method are not implemented" (FIC-NEVPT2), or
crashes in the response solver; this was measured over six FIC methods and three
density kinds (2026-09-26, records in the project session log).

The route used here: rebuild the wavefunction from data ORCA *does* write.
ORCA's ``!FCIDUMP`` keyword dumps the converged active-space Hamiltonian
(effective one-electron integrals, chemist-notation two-electron integrals and
the core energy) into a standard FCIDUMP.  This module solves that small CAS-CI
with its own determinant solver (exact, deterministic), rebuilds the
spin-resolved density objects, optionally rotates them into the localised
orbital basis, and evaluates the four-state entropy.  Nothing is taken from
undocumented files.

Why the three numbers are enough (fixed-particle-number argument)
-----------------------------------------------------------------
For a state of fixed total particle number N, the one-orbital reduced density
matrix in the four Fock states of one spatial orbital is exactly *diagonal* in
any orbital basis:

- the |0> <-> |up-down> coherence needs a rest configuration shared by a
  zero-occupancy and a double-occupancy bra/ket, which would have to carry both
  N_rest = N and N_rest = N - 2 -- impossible;
- the |up> <-> |down> coherence is the spin-flip expectation <a+_up a_down>,
  which vanishes for any S_z eigenstate (the CI is solved in a fixed S_z
  sector, so the eigenstates are S_z eigenstates).

Its diagonal is exactly the four weights above, so the eigenvalue entropy is
the population entropy and the three numbers (n_a, n_b, p2) are complete.

What is validated, and how
--------------------------
The reconstruction is only as good as its cross-checks, and the module computes
them instead of trusting itself:

- the CI energy eigenvalue must match the CASSCF energy printed in the run
  output (selected by ``reference_energy`` when given); measured agreement on
  the N2 fixture is at the 1e-13 level;
- the natural occupations of the reconstructed 1-RDM must match the ``N(occ)=``
  line of the output (the caller compares them; :func:`natural_occupations`
  supplies the numbers);
- :func:`energy_from_densities` recomputes the active-space energy from the
  density objects alone -- a convention test for the tensor blocks that also
  makes the rotation test possible (the energy is invariant under an orbital
  rotation of the integrals and the densities together);
- the four-state weights of a single determinant collapse to one unit weight
  (zero entropy), and the weights always sum to one.

The rotation from canonical to localised active orbitals is built by
:func:`rotation_from_coefficients` from the two coefficient matrices and the AO
overlap matrix (the ``orca_2json`` export of the canonical and the ``orca_loc``
gbw files).  Entropies/plateau/threshold language is only applied to the
localised-basis spectrum, in step with the A2 preconditions.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np

from ..knowledge.models import (
    EVIDENCE_LITERATURE,
    EVIDENCE_MANUAL,
    EVIDENCE_MEASURED,
    Evidence,
    ReportSection,
)
from ..parsers.fcidump import Fcidump

#: Deterministic dense CI: the solver stores a dense Hamiltonian, so the
#: determinant count is capped (the cap is about memory/time, not about CI
#: theory -- the f-block use cases sit far below it).
DEFAULT_MAX_DETERMINANTS = 3000

#: |<S^2> - S(S+1)| tolerance when selecting a state by multiplicity.
S2_TOLERANCE = 0.05


class EntropyRdmError(ValueError):
    """The exact-entropy route cannot run on the given data (with a next step)."""


@dataclass(frozen=True)
class FciState:
    """One CI eigenstate of the FCIDUMP Hamiltonian (active space only)."""

    norb: int
    nelec: int
    ms2: int
    na: int
    nb: int
    energy_active: float
    ecore: float
    s2: float
    root_index: int
    determinants: tuple[tuple[int, int], ...]
    coefficients: tuple[float, ...]

    @property
    def energy_total(self) -> float:
        return self.energy_active + self.ecore


@dataclass(frozen=True)
class SpinDensities:
    """Spin-resolved density objects, in the E-operator convention.

    All arrays are index-major ``[p, q]`` / ``[p, q, r, s]`` with the operator
    definitions::

        gamma_a[p, q] = <a+_p,alpha a_q,alpha>
        gamma_b[p, q] = <a+_p,beta  a_q,beta>
        g_aa[p, q, r, s] = <a+_p,alpha a+_r,alpha a_s,alpha a_q,alpha>
        g_ab[p, q, r, s] = <a+_p,alpha a+_r,beta  a_s,beta  a_q,alpha>
        g_ba[p, q, r, s] = <a+_p,beta  a+_r,alpha a_s,alpha a_q,beta>
        g_bb[p, q, r, s] = <a+_p,beta  a+_r,beta  a_s,beta  a_q,beta>

    The diagonal ``g_ab[i, i, i, i]`` is the double-occupancy probability
    ``<n_i,alpha n_i,beta>`` of orbital ``i``.
    """

    gamma_a: np.ndarray
    gamma_b: np.ndarray
    g_aa: np.ndarray
    g_ab: np.ndarray
    g_ba: np.ndarray
    g_bb: np.ndarray

    @property
    def norb(self) -> int:
        return int(self.gamma_a.shape[0])


# --------------------------------------------------------------------------
# determinant machinery
# --------------------------------------------------------------------------

_ANN = "a"  # apply an annihilation operator
_CRE = "c"  # apply a creation operator


def _strings(norb: int, nocc: int) -> tuple[int, ...]:
    return tuple(sum(1 << i for i in combo) for combo in itertools.combinations(range(norb), nocc))


def _annihilate(string: int, k: int) -> tuple[int, float] | None:
    if not (string >> k) & 1:
        return None
    sign = -1.0 if bin(string & ((1 << k) - 1)).count("1") % 2 else 1.0
    return string ^ (1 << k), sign


def _create(string: int, k: int) -> tuple[int, float] | None:
    if (string >> k) & 1:
        return None
    sign = -1.0 if bin(string & ((1 << k) - 1)).count("1") % 2 else 1.0
    return string | (1 << k), sign


class _Tables:
    """Vectorised operator application over the determinant list.

    ``apply(chain)`` evaluates an operator chain on *every* determinant at
    once.  Single operators act on the alpha/beta occupation strings (an
    intermediate of a chain is generally *not* a determinant, so the string
    level is the right one); only after the last operator are the two strings
    mapped back to a determinant index.  Both the Hamiltonian build and every
    expectation value go through this.  The per-operator string tables are
    dense in the string value, which bounds the orbital count this route
    handles (the determinant cap binds first for every realistic window).
    """

    _MAX_STRING_ORBITALS = 16

    def __init__(self, dets: Sequence[tuple[int, int]], norb: int) -> None:
        if norb > self._MAX_STRING_ORBITALS:
            raise EntropyRdmError(
                f"the deterministic solver handles at most "
                f"{self._MAX_STRING_ORBITALS} active orbitals (got {norb}). "
                "Next step: shrink the active space, or use the A2 bound analysis "
                "for this run."
            )
        self.norb = norb
        self.ndet = len(dets)
        a_strings = sorted({d[0] for d in dets})
        b_strings = sorted({d[1] for d in dets})
        self._nb_strings = len(b_strings)
        size = 1 << norb
        a_index = np.full(size, -1, dtype=np.int64)
        b_index = np.full(size, -1, dtype=np.int64)
        for i, s in enumerate(a_strings):
            a_index[s] = i
        for i, s in enumerate(b_strings):
            b_index[s] = i
        self._a_index, self._b_index = a_index, b_index
        self._a0 = np.array([d[0] for d in dets], dtype=np.int64)
        self._b0 = np.array([d[1] for d in dets], dtype=np.int64)
        self._string_tables: dict[tuple[str, int, int], tuple[np.ndarray, np.ndarray]] = {}
        for spin in (0, 1):
            for kind in (_ANN, _CRE):
                op = _annihilate if kind == _ANN else _create
                for k in range(norb):
                    new = np.zeros(size, dtype=np.int64)
                    sign = np.zeros(size, dtype=np.int8)
                    for s in range(size):
                        result = op(s, k)
                        if result is not None:
                            new[s], sign[s] = result[0], int(result[1])
                    self._string_tables[(kind, spin, k)] = (new, sign.astype(np.float64))

    def apply(self, chain: Iterable[tuple[str, int, int]]) -> tuple[np.ndarray, np.ndarray]:
        """Apply an operator chain (in application order: first element acts first).

        Returns ``(target_index, sign)`` arrays over the determinant list;
        dead branches carry index -1 and sign 0.
        """
        a = self._a0
        b = self._b0
        sign = np.ones(self.ndet, dtype=np.float64)
        for kind, spin, k in chain:
            new, s = self._string_tables[(kind, spin, k)]
            current = a if spin == 0 else b
            sign = sign * s[current]
            if spin == 0:
                a = new[current]
            else:
                b = new[current]
        ai = self._a_index[a]
        bi = self._b_index[b]
        invalid = (ai < 0) | (bi < 0)
        target = np.where(invalid, -1, ai * self._nb_strings + bi)
        sign = np.where(invalid, 0.0, sign)
        return target, sign


def _chain_2e(spin_p: int, p: int, spin_r: int, r: int, spin_s: int, s: int, spin_q: int, q: int):
    """Application order for <a+_p a+_r a_s a_q> (rightmost operator first)."""
    return [
        (_ANN, spin_q, q),
        (_ANN, spin_s, s),
        (_CRE, spin_r, r),
        (_CRE, spin_p, p),
    ]


# --------------------------------------------------------------------------
# CI solver
# --------------------------------------------------------------------------


def dispatch_determinants(norb: int, na: int, nb: int) -> tuple[tuple[int, int], ...]:
    """All (alpha-string, beta-string) determinant pairs of an (na, nb) sector.

    Strings are ordered by their bit value so that the determinant index is
    ``a_position * n_beta_strings + b_position`` with both positions in the
    same order the lookup tables use (itertools order is *not* bit order).
    """
    a_strings = sorted(_strings(norb, na))
    b_strings = sorted(_strings(norb, nb))
    return tuple((a, b) for a in a_strings for b in b_strings)


def solve_fci(
    dump: Fcidump,
    *,
    reference_energy: float | None = None,
    multiplicity: int | None = None,
    max_determinants: int = DEFAULT_MAX_DETERMINANTS,
    energy_tolerance: float = 1e-6,
) -> FciState:
    """Solve the CAS-CI of an FCIDUMP and return the selected eigenstate.

    Selection: with ``reference_energy`` (the CASSCF energy printed by the run
    that produced the orbitals) pick the eigenvalue closest to it -- at
    converged orbitals that eigenvalue *is* the CASSCF energy, which makes the
    match stronger than any heuristic.  Otherwise pick the lowest root whose
    <S^2> matches ``multiplicity`` (default: the FCIDUMP's own MS2 sector).
    """
    norb, nelec, ms2 = dump.norb, dump.nelec, dump.ms2
    if (nelec - ms2) % 2 != 0:
        raise EntropyRdmError(
            f"NELEC={nelec} and MS2={ms2} have different parity. Next step: the "
            "FCIDUMP header is inconsistent; re-dump it."
        )
    na, nb = (nelec + ms2) // 2, (nelec - ms2) // 2
    if not (0 <= nb <= na <= norb):
        raise EntropyRdmError(
            f"the FCIDUMP electron count does not fit the orbital count "
            f"(NELEC={nelec}, MS2={ms2}, NORB={norb}). Next step: check the dump."
        )
    dets = dispatch_determinants(norb, na, nb)
    if len(dets) > max_determinants:
        raise EntropyRdmError(
            f"the CAS-CI has {len(dets)} determinants in the Ms sector, above the "
            f"deterministic solver's cap of {max_determinants}. Next step: use the "
            "exact route on a smaller active space (the f-block windows it targets "
            "sit far below the cap), or keep the A2 bound analysis for this run."
        )
    tables = _Tables(dets, norb)
    ndet = len(dets)

    h = np.asarray(dump.h, dtype=np.float64)
    hamiltonian = np.zeros((ndet, ndet), dtype=np.float64)
    rows = np.arange(ndet, dtype=np.int64)
    for p in range(norb):
        for q in range(norb):
            v = h[p, q]
            if v == 0.0:
                continue
            for spin in (0, 1):
                idx, sign = tables.apply([(_ANN, spin, q), (_CRE, spin, p)])
                alive = sign != 0.0
                np.add.at(hamiltonian, (idx[alive], rows[alive]), v * sign[alive])
    for (p, q, r, s), v in dump.g.items():
        if v == 0.0:
            continue
        for spin_p, spin_r in ((0, 0), (0, 1), (1, 0), (1, 1)):
            chain = _chain_2e(spin_p, p, spin_r, r, spin_r, s, spin_p, q)
            idx, sign = tables.apply(chain)
            alive = sign != 0.0
            np.add.at(hamiltonian, (idx[alive], rows[alive]), 0.5 * v * sign[alive])
    if not np.allclose(hamiltonian, hamiltonian.T, atol=1e-10):
        raise EntropyRdmError(
            "internal error: the CI Hamiltonian came out non-symmetric. Next step: "
            "this is a defect of the solver; keep the fixture and report it."
        )

    s2_matrix = _spin_squared_matrix(tables, norb, na, nb, ndet, rows)
    eigenvalues, eigenvectors = np.linalg.eigh(hamiltonian)

    root = _select_root(
        eigenvalues,
        eigenvectors,
        s2_matrix,
        reference_active=(
            None if reference_energy is None else reference_energy - dump.ecore
        ),
        multiplicity=multiplicity,
        energy_tolerance=energy_tolerance,
    )
    vector = eigenvectors[:, root]
    s2 = float(vector @ s2_matrix @ vector)
    if reference_energy is not None and multiplicity is not None:
        # the energy-selected state must also carry the requested multiplicity;
        # a mismatch means the two inputs disagree and must not be papered over
        target = (multiplicity**2 - 1) / 4.0
        if abs(s2 - target) > S2_TOLERANCE:
            raise EntropyRdmError(
                f"the CI root matching the reference energy ({eigenvalues[root]:.9f} Eh) "
                f"has <S^2> = {s2:.3f}, not the {target:.3f} of multiplicity "
                f"{multiplicity}. Next step: check that the run output and the FCIDUMP "
                "belong to the same job (same active space and multiplicity)."
            )
    return FciState(
        norb=norb,
        nelec=nelec,
        ms2=ms2,
        na=na,
        nb=nb,
        energy_active=float(eigenvalues[root]),
        ecore=dump.ecore,
        s2=s2,
        root_index=int(root),
        determinants=dets,
        coefficients=tuple(float(c) for c in vector),
    )


def _spin_squared_matrix(tables: _Tables, norb: int, na: int, nb: int, ndet: int, rows) -> np.ndarray:
    """<S^2> = S_z^2 + S_z + sum_pq a+_p,beta a_p,alpha a+_q,alpha a_q,beta (matrix form)."""
    ms = 0.5 * (na - nb)
    matrix = np.eye(ndet, dtype=np.float64) * (ms * ms + ms)
    for p in range(norb):
        for q in range(norb):
            chain = [(_ANN, 1, q), (_CRE, 0, q), (_ANN, 0, p), (_CRE, 1, p)]
            idx, sign = tables.apply(chain)
            alive = sign != 0.0
            np.add.at(matrix, (idx[alive], rows[alive]), sign[alive])
    return matrix


def _select_root(
    eigenvalues,
    eigenvectors,
    s2_matrix,
    *,
    reference_active,
    multiplicity,
    energy_tolerance,
):
    if reference_active is not None:
        best = int(np.argmin(np.abs(eigenvalues - reference_active)))
        if abs(eigenvalues[best] - reference_active) > energy_tolerance:
            raise EntropyRdmError(
                f"no CI eigenvalue matches the reference CASSCF energy within "
                f"{energy_tolerance} Eh (closest active-space value: "
                f"{eigenvalues[best]:.9f} at root {best}; reference minus core "
                f"energy: {reference_active:.9f}). Next step: check that the "
                "FCIDUMP was dumped for *these* converged orbitals (same job, same "
                "active space); a mismatch means the two files are unrelated."
            )
        return best
    if multiplicity is None:
        return 0
    target = (multiplicity**2 - 1) / 4.0
    for root in range(len(eigenvalues)):
        vector = eigenvectors[:, root]
        s2 = float(vector @ s2_matrix @ vector)
        if abs(s2 - target) <= S2_TOLERANCE:
            return root
    raise EntropyRdmError(
        f"no root in the FCIDUMP's M_s sector has <S^2> matching multiplicity "
        f"{multiplicity} (target {target:.2f}). Next step: check the requested "
        "multiplicity against the active-space electron count."
    )


# --------------------------------------------------------------------------
# density objects
# --------------------------------------------------------------------------


def _expect(tables: _Tables, coefficients: np.ndarray, chain) -> float:
    idx, sign = tables.apply(chain)
    idx = np.where(idx < 0, 0, idx)
    return float(np.dot(coefficients * sign, coefficients[idx]))


def spin_densities(state: FciState) -> SpinDensities:
    """Rebuild the spin-resolved one- and two-particle density objects.

    Definition labels follow :class:`SpinDensities`; the blocks are computed
    directly over the CI vector with the same vectorised tables that built the
    Hamiltonian.
    """
    norb = state.norb
    tables = _Tables(state.determinants, norb)
    c = np.asarray(state.coefficients, dtype=np.float64)
    gamma_a = np.zeros((norb, norb))
    gamma_b = np.zeros((norb, norb))
    g_aa = np.zeros((norb, norb, norb, norb))
    g_ab = np.zeros((norb, norb, norb, norb))
    g_ba = np.zeros((norb, norb, norb, norb))
    g_bb = np.zeros((norb, norb, norb, norb))
    for p in range(norb):
        for q in range(norb):
            gamma_a[p, q] = _expect(tables, c, [(_ANN, 0, q), (_CRE, 0, p)])
            gamma_b[p, q] = _expect(tables, c, [(_ANN, 1, q), (_CRE, 1, p)])
            for r in range(norb):
                for s in range(norb):
                    g_aa[p, q, r, s] = _expect(
                        tables, c, _chain_2e(0, p, 0, r, 0, s, 0, q)
                    )
                    g_ab[p, q, r, s] = _expect(
                        tables, c, _chain_2e(0, p, 1, r, 1, s, 0, q)
                    )
                    g_ba[p, q, r, s] = _expect(
                        tables, c, _chain_2e(1, p, 0, r, 0, s, 1, q)
                    )
                    g_bb[p, q, r, s] = _expect(
                        tables, c, _chain_2e(1, p, 1, r, 1, s, 1, q)
                    )
    return SpinDensities(
        gamma_a=gamma_a, gamma_b=gamma_b, g_aa=g_aa, g_ab=g_ab, g_ba=g_ba, g_bb=g_bb
    )


def natural_occupations(dens: SpinDensities) -> tuple[float, ...]:
    """Eigenvalues of the spin-summed 1-RDM (descending) -- the N(occ) cross-check."""
    total = dens.gamma_a + dens.gamma_b
    values = np.linalg.eigvalsh((total + total.T) / 2.0)
    return tuple(float(v) for v in values[::-1])


def spin_occupations(dens: SpinDensities) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Diagonal (in the current basis) alpha and beta occupations, in orbital order."""
    return (
        tuple(float(v) for v in np.diag(dens.gamma_a)),
        tuple(float(v) for v in np.diag(dens.gamma_b)),
    )


def energy_from_densities(dump: Fcidump, dens: SpinDensities) -> float:
    """Rebuild the active-space energy from the density objects alone.

    Convention test for the tensor blocks (it must reproduce the CI energy) and
    the vehicle for the rotation-invariance test: rotating the integrals and the
    density objects together must leave this number unchanged.
    """
    gamma = dens.gamma_a + dens.gamma_b
    norb = dens.norb
    energy = 0.0
    h = np.asarray(dump.h, dtype=np.float64)
    for p in range(norb):
        for q in range(norb):
            energy += h[p, q] * gamma[q, p]
    # E_pq E_rs = delta_qr E_ps + sum_sigma,tau a+_p,sigma a+_r,tau a_s,tau a_q,sigma
    # (the anticommutator step), so in H2 = 1/2 sum (pq|rs)(E_pq E_rs - delta_qr E_ps)
    # the delta terms cancel exactly and only the four spin blocks remain
    two_body = dens.g_aa + dens.g_ab + dens.g_ba + dens.g_bb
    for (p, q, r, s), v in dump.g.items():
        energy += 0.5 * v * two_body[p, q, r, s]
    return energy


# --------------------------------------------------------------------------
# rotation to the localised basis
# --------------------------------------------------------------------------


def rotation_from_coefficients(
    canonical: np.ndarray,
    localized: np.ndarray,
    overlap: np.ndarray,
    active: Sequence[int],
) -> np.ndarray:
    """The active-block orbital rotation ``u[canonical, localized]``.

    ``canonical`` / ``localized`` are AO coefficient matrices with the active
    orbitals in the columns given by ``active``; ``overlap`` is the AO overlap
    matrix S.  For orthonormal orbitals ``u = C_can^T S C_loc`` is orthogonal,
    which callers can (and tests do) assert.
    """
    can = np.asarray(canonical, dtype=np.float64)[:, list(active)]
    loc = np.asarray(localized, dtype=np.float64)[:, list(active)]
    s = np.asarray(overlap, dtype=np.float64)
    if can.shape != loc.shape:
        raise EntropyRdmError(
            f"coefficient blocks differ in shape: {can.shape} vs {loc.shape}. "
            "Next step: the two orca_2json exports must come from the same gbw file "
            "(same basis, same number of orbitals)."
        )
    if s.shape != (can.shape[0], can.shape[0]):
        raise EntropyRdmError(
            f"the overlap matrix is {s.shape} but the coefficient block has "
            f"{can.shape[0]} AO rows. Next step: export the S-Matrix together with "
            "the MO coefficients (same orca_2json run)."
        )
    return can.T @ s @ loc


def rotate_densities(dens: SpinDensities, u: np.ndarray) -> SpinDensities:
    """Rotate every density object to the new orbital basis.

    ``u[old, new]`` with the operator convention ``a+_new,i = sum_p u[p, i] a+_p``;
    expectation values of the new operators in the same state are then exact
    contractions, e.g. ``gamma' = u^T gamma u`` and a four-index transform for
    each two-particle block.
    """
    u = np.asarray(u, dtype=np.float64)
    n = dens.norb
    if u.shape != (n, n):
        raise EntropyRdmError(
            f"rotation matrix is {u.shape} but the active space has {n} orbitals. "
            "Next step: build the rotation from the active columns of both "
            "coefficient matrices."
        )
    return SpinDensities(
        gamma_a=u.T @ dens.gamma_a @ u,
        gamma_b=u.T @ dens.gamma_b @ u,
        g_aa=np.einsum("pi,qj,rk,sl,pqrs->ijkl", u, u, u, u, dens.g_aa),
        g_ab=np.einsum("pi,qj,rk,sl,pqrs->ijkl", u, u, u, u, dens.g_ab),
        g_ba=np.einsum("pi,qj,rk,sl,pqrs->ijkl", u, u, u, u, dens.g_ba),
        g_bb=np.einsum("pi,qj,rk,sl,pqrs->ijkl", u, u, u, u, dens.g_bb),
    )


# --------------------------------------------------------------------------
# the four-state entropy
# --------------------------------------------------------------------------


def orbital_weights(dens: SpinDensities, i: int) -> tuple[float, float, float, float]:
    """The four occupation weights of orbital ``i``: (empty, up, down, double).

    Negative round-off (|w| < 1e-12) is snapped to zero; a genuinely negative
    weight signals an inconsistent density set and raises.
    """
    n_a = float(dens.gamma_a[i, i])
    n_b = float(dens.gamma_b[i, i])
    p2 = float(dens.g_ab[i, i, i, i])
    w = (1.0 - n_a - n_b + p2, n_a - p2, n_b - p2, p2)
    snapped = []
    for value in w:
        if value < 0.0:
            if value < -1e-8:
                raise EntropyRdmError(
                    f"negative four-state weight {value:.3e} for orbital {i}: the "
                    "density objects are inconsistent. Next step: verify that the "
                    "CI solution really is the intended state (check the energy "
                    "against the run output)."
                )
            value = 0.0
        snapped.append(float(value))
    return tuple(snapped)  # type: ignore[return-value]


def orbital_entropy(dens: SpinDensities, i: int) -> float:
    """The four-state single-orbital entropy of orbital ``i`` (exact, any basis)."""
    entropy = 0.0
    for w in orbital_weights(dens, i):
        if w > 1e-14:
            entropy -= w * math.log(w)
    return entropy + 0.0  # normalise a signed zero to +0.0 for clean reporting


def entropy_spectrum(dens: SpinDensities) -> tuple[float, ...]:
    """The four-state entropy of every active orbital, in orbital order."""
    return tuple(orbital_entropy(dens, i) for i in range(dens.norb))


# --------------------------------------------------------------------------
# orchestration: the analysis record, the active window, the report section
# --------------------------------------------------------------------------


def infer_active_window(
    occupations: Sequence[float],
    norb_expected: int,
    *,
    reference_occupations: Sequence[float] | None = None,
) -> tuple[tuple[int, ...], str]:
    """Locate the active-orbital window inside the ORCA occupation table.

    ``occupations`` is the OCC column of the last ORBITAL ENERGIES table (0-based
    orbital order).  The window is the contiguous span that carries every
    fractionally-occupied orbital; when the span is shorter than the FCIDUMP's
    orbital count (an active orbital can print occupation 0.0000), it is grown
    towards the neighbour farther from a fully occupied core orbital.  With
    ``reference_occupations`` (the ``N(occ)=`` line) the chosen window is
    validated against those occupations -- a shifted window of the right size
    would otherwise pass unnoticed.  Returns ``(window, note)``; the caller can
    always override the window by hand.
    """
    fractional = [
        i for i, occ in enumerate(occupations) if 1e-4 < occ < 1.9999
    ]
    if not fractional:
        raise EntropyRdmError(
            "no fractionally-occupied orbital in the occupation table. Next step: "
            "check that the output is the converged CASSCF run of this active space."
        )
    low, high = fractional[0], fractional[-1]
    note = "span of the fractionally-occupied orbitals"
    while high - low + 1 < norb_expected:
        # an active orbital can print 0.0000 (or 2.0000); grow towards the
        # neighbour whose occupation is farther from 2.0 -- a core orbital, the
        # wrong side, carries 2.0, while an empty active orbital carries 0.0
        below = occupations[low - 1] if low > 0 else None
        above = occupations[high + 1] if high < len(occupations) - 1 else None
        take_below = None
        if below is not None and above is None:
            take_below = True
        elif above is not None and below is None:
            take_below = False
        elif below is not None and above is not None:
            take_below = below < above
        if take_below is None:
            break
        if take_below:
            low -= 1
            note = (
                "span extended downwards (the active space contains an orbital "
                "printed as 0.0000 or 2.0000)"
            )
        else:
            high += 1
            note = (
                "span extended upwards (the active space contains an orbital "
                "printed as 0.0000 or 2.0000)"
            )
    window = tuple(range(low, high + 1))
    if len(window) != norb_expected:
        raise EntropyRdmError(
            f"the inferred active window has {len(window)} orbital(s) but the FCIDUMP "
            f"has {norb_expected}. Next step: give the window explicitly (the "
            "first and last active orbital, as numbered in ORCA's orbital table)."
        )
    if reference_occupations is not None:
        chosen = sorted((occupations[i] for i in window), reverse=True)
        printed = sorted(reference_occupations, reverse=True)
        if len(chosen) == len(printed):
            deviation = max(abs(a - b) for a, b in zip(chosen, printed))
            if deviation > 5e-4:
                raise EntropyRdmError(
                    f"the inferred active window's occupations do not match the printed "
                    f"N(occ)= line (max |deviation| = {deviation:.2e}). Next step: give "
                    "the window explicitly (first and last active orbital, as numbered "
                    "in ORCA's orbital table)."
                )
    return window, note


@dataclass(frozen=True)
class ExactEntropyAnalysis:
    """The complete record of one exact-entropy analysis (cross-checks included)."""

    state: FciState
    densities: SpinDensities
    canonical_spectrum: tuple[float, ...]
    occupation_deviation: float | None
    checks: tuple[str, ...]
    localized_spectrum: tuple[float, ...] | None = None
    localized_weights: tuple[tuple[float, float, float, float], ...] | None = None
    active_window: tuple[int, ...] | None = None


def analyze(
    dump: Fcidump,
    *,
    reference_energy: float,
    multiplicity: int | None = None,
    reference_occupations: Sequence[float] | None = None,
    canonical=None,
    localized=None,
    active_window: Sequence[int] | None = None,
) -> ExactEntropyAnalysis:
    """Run the whole exact route and collect the engine cross-checks.

    ``canonical`` / ``localized`` are :class:`~fblockkit.parsers.orca_json.OrcaJson`
    exports of the gbw before and after ``orca_loc``; both must be given (with
    ``active_window``) for the localized-basis spectrum.
    """
    state = solve_fci(
        dump, reference_energy=reference_energy, multiplicity=multiplicity
    )
    densities = spin_densities(state)
    canonical_spectrum = entropy_spectrum(densities)
    checks = [
        f"CI energy + core energy = {state.energy_total:.12f} Eh vs the run's "
        f"{reference_energy:.12f} Eh (deviation {abs(state.energy_total - reference_energy):.2e}; "
        "the solver refuses a match beyond 1e-6)"
    ]
    occupation_deviation: float | None = None
    if reference_occupations is not None:
        computed = sorted(natural_occupations(densities), reverse=True)
        printed = sorted(reference_occupations, reverse=True)
        if len(computed) != len(printed):
            raise EntropyRdmError(
                f"the output printed {len(printed)} active occupation(s) but the FCIDUMP "
                f"has {len(computed)} active orbital(s). Next step: check that the "
                "output and the FCIDUMP belong to the same active space."
            )
        occupation_deviation = max(
            abs(a - b) for a, b in zip(computed, printed)
        )
        checks.append(
            f"natural occupations vs the printed N(occ): max |deviation| = "
            f"{occupation_deviation:.2e} (the print precision is 1e-5)"
        )
    localized_spectrum = None
    localized_weights = None
    window: tuple[int, ...] | None = None
    if canonical is not None or localized is not None:
        if canonical is None or localized is None:
            raise EntropyRdmError(
                "the localized spectrum needs BOTH orca_2json exports (the canonical "
                "and the orca_loc one). Next step: export the two gbw files with "
                "orca_2json and give both paths."
            )
        if active_window is None:
            raise EntropyRdmError(
                "the localized-basis rotation needs the active-orbital window. "
                "Next step: give the window explicitly (first and last active "
                "orbital, as numbered in ORCA's orbital table)."
            )
        window = tuple(int(i) for i in active_window)
        if len(window) != dump.norb:
            raise EntropyRdmError(
                f"the active window has {len(window)} orbital(s) but the FCIDUMP has "
                f"{dump.norb}. Next step: check the window against the orbital table."
            )
        rotation = rotation_from_coefficients(
            np.array(canonical.mo_coefficients).T,
            np.array(localized.mo_coefficients).T,
            np.array(canonical.overlap),
            active=window,
        )
        orthogonality = float(np.abs(rotation.T @ rotation - np.eye(dump.norb)).max())
        checks.append(
            f"localization rotation orthogonality: max |u^T u - I| = {orthogonality:.2e} "
            "(two exports of the same gbw, one localized window)"
        )
        localized_densities = rotate_densities(densities, rotation)
        localized_spectrum = entropy_spectrum(localized_densities)
        localized_weights = tuple(
            orbital_weights(localized_densities, i) for i in range(dump.norb)
        )
    return ExactEntropyAnalysis(
        state=state,
        densities=densities,
        canonical_spectrum=canonical_spectrum,
        occupation_deviation=occupation_deviation,
        checks=tuple(checks),
        localized_spectrum=localized_spectrum,
        localized_weights=localized_weights,
        active_window=window,
    )


def _largest_gap(values: Sequence[float]) -> tuple[int, float]:
    if len(values) < 2:
        return len(values), 0.0
    gaps = [
        (values[i] - values[i + 1]) / values[i] if values[i] > 0 else 0.0
        for i in range(len(values) - 1)
    ]
    index = max(range(len(gaps)), key=lambda i: gaps[i])
    return index + 1, gaps[index]


def run(analysis: ExactEntropyAnalysis, *, threshold: float = 0.14) -> ReportSection:
    """Render the A2-complement report section (exact route, cross-checks included)."""
    state = analysis.state
    lines = [
        f"Active space: {state.norb} orbitals / {state.nelec} electrons "
        f"(M_s sector with n_alpha={state.na}, n_beta={state.nb}); "
        f"determinants solved: {len(state.determinants)}; selected root "
        f"{state.root_index} with <S^2> = {state.s2:.4f}",
        "",
        "Cross-checks against the engine (the reconstruction is trusted only with them):",
    ]
    for check in analysis.checks:
        lines.append(f"  - {check}")
    lines += [
        "",
        "Four-state single-orbital entropy (canonical orbital basis of the dump; "
        "informational -- the 0.14 line and the plateau reading are defined in a "
        "localized basis):",
        "  " + " ".join(f"{value:.4f}" for value in analysis.canonical_spectrum),
    ]
    if analysis.localized_spectrum is not None:
        window = analysis.active_window
        gap_split, gap = _largest_gap(
            sorted(analysis.localized_spectrum, reverse=True)
        )
        above = [
            window[i]
            for i, value in enumerate(analysis.localized_spectrum)
            if value > threshold
        ]
        lines += [
            "",
            f"Four-state single-orbital entropy (localized basis: IAO-IBO rotation of "
            f"the active window {list(window)} from the user's orca_loc):",
            "  " + " ".join(f"{value:.4f}" for value in analysis.localized_spectrum),
            f"  max s1 = {max(analysis.localized_spectrum):.4f} vs the Stein & Reiher "
            f"line {threshold}: " + (
                f"candidate orbitals (original numbering): {above}"
                if above
                else "no orbital above the line"
            ),
        ]
        if len(analysis.localized_spectrum) >= 2:
            if gap >= 0.25:
                lines.append(
                    f"  Plateau/cliff: the largest relative gap ({gap:.1%}) sits after "
                    f"item {gap_split} of the descending order."
                )
            else:
                lines.append(
                    f"  No significant plateau in this window (largest relative gap "
                    f"{gap:.1%}); this is a property of the window, not of a solver."
                )
        lines += [
            "",
            "Four weights per localized orbital (empty, up, down, double):",
        ]
        for i, weights in enumerate(analysis.localized_weights or ()):
            orbital = window[i] if window is not None else i
            lines.append(
                f"  orbital {orbital}: " + " ".join(f"{w:.4f}" for w in weights)
            )
    lines += [
        "",
        "Preconditions carried with these numbers:",
        "- The four-state density is diagonal for a fixed-particle-number state (an "
        "exact statement, not an approximation), so the three numbers n_up, n_down and "
        "p2 are complete -- and they exist only in a reconstruction like this one; a "
        "printed N(occ) line cannot supply them.",
        "- Read the threshold against the localized-basis spectrum only; the canonical "
        "spectrum is reported for orientation.",
        "- The solver is exact within the dumped active space; the usual active-space "
        "caveats (size, composition, solution branch) transfer unchanged.",
    ]
    return ReportSection(
        title="A2x exact four-state single-orbital entropy (FCIDUMP route)",
        body="\n".join(lines),
    )


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the exact route (rule/evidence discipline of the project)."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "The single-orbital entropy is defined from the four occupation states "
                "of an orbital (empty / up / down / up-down), s = -sum w ln w; the "
                "0.14 multireference line and the localized-basis precondition are the "
                "autoCAS protocol's."
            ),
            ref=(
                "Stein C. J., Reiher M., J. Chem. Theory Comput., 2016, 12(4), "
                "1760-1771, DOI 10.1021/acs.jctc.6b00156"
            ),
            bibkey="stein2016automated",
            url="https://doi.org/10.1021/acs.jctc.6b00156",
        ),
        Evidence(
            kind=EVIDENCE_MANUAL,
            text=(
                "ORCA writes the active-space Hamiltonian to a standard FCIDUMP with "
                "the !FCIDUMP keyword; the dump run reports 'IS NOT FULLY CONVERGED' "
                "by design (it saves the integrals after its first macro-iteration). "
                "The manual documents the .RDM2 export route for 2-RDMs via "
                "orca_2json (section 9.3); on 6.1.1 that route was measured to be "
                "unreachable for CAS-type methods (six FIC methods and three density "
                "requests; failure modes recorded in this module's docstring)."
            ),
            ref="ORCA 6.1 manual sections 9.3 (orca_2json) and the CASSCF module; measured 2026-09-26",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "Route validation (2026-09-26, fixture n2_fcidump.*): the CI energy "
                "reproduces the printed CASSCF energy at 4e-13 Eh, the natural "
                "occupations reproduce the N(occ)= line, and the energy rebuilt from "
                "the density objects is invariant under an orbital rotation applied "
                "to the integrals and the densities together (tests/test_entropy_rdm.py)."
            ),
            ref="tests/test_entropy_rdm.py; fixtures/orca/README.md (FCIDUMP-route section)",
        ),
    )


__all__ = [
    "DEFAULT_MAX_DETERMINANTS",
    "EntropyRdmError",
    "ExactEntropyAnalysis",
    "FciState",
    "SpinDensities",
    "analyze",
    "dispatch_determinants",
    "energy_from_densities",
    "entropy_spectrum",
    "evidence",
    "infer_active_window",
    "natural_occupations",
    "orbital_entropy",
    "orbital_weights",
    "rotate_densities",
    "rotation_from_coefficients",
    "run",
    "solve_fci",
    "spin_densities",
    "spin_occupations",
]
