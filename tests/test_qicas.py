"""Regression checks of the QICAS optimization (analysis/qicas.py).

The fixture is the N2/def2-SVP CAS(6,6) FCIDUMP (the four-state-entropy route's
own file, already validated against the run's CASSCF energy).  Two targets are
pinned:

- **(4,4)-in-(6,6)** (one closed, four active, one virtual): the dumped
  natural orbitals are already F_QI-minimal for this split -- F_QI = 0.030170
  with zero accepted rotations -- and the CASCI in that partition sits
  +5.331e-3 Eh above the window FCI, inside the Theorem-1 bound (1.034e-1).
- **(2,4)-in-(6,6)** (a deliberately awkward request that splits the pi pair):
  F_QI starts at 0.248094 (CASCI +6.697e-2 off), the Jacobi sweep accepts two
  rotations taking F_QI to the same 0.030170 fixed point, and the occupancy
  reading of the optimized basis reclassifies the space to the (4e, 4o)
  arrangement -- the CASCI gap improves to +5.331e-3, equal to the (4,4)
  target's.  That is the source's story on real data, with the run's own
  numbers.

Infrastructure checks pin two implementation laws: the accumulated rotation
matrix must replay the density rotations exactly (the ``u[old, new]``
composition law -- the module's first bug), and the O(1) pair-slice entropy
evaluation used inside the Jacobi scan must equal a full density rotation.
"""

from __future__ import annotations

import numpy as np
import pytest

from fblockkit.analysis import entropy_rdm, qicas
from fblockkit.analysis.qicas import QicasError, partition_by_occupation
from fblockkit.parsers.fcidump import parse_fcidump
from fblockkit.parsers.orca_json import parse_orca_json

from pathlib import Path

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"
FCIDUMP = FIXTURES / "n2_fcidump.fcidump"

REFERENCE = -108.950671945  # the run's printed CASSCF energy (the menu-12 anchor)


def _dump():
    return parse_fcidump(FCIDUMP)


def _dens():
    return entropy_rdm.spin_densities(entropy_rdm.solve_fci(_dump(), reference_energy=REFERENCE))


# --- infrastructure laws -------------------------------------------------------


def test_partition_orders_by_occupation_and_validates():
    occupations = (0.02, 1.9, 0.05, 1.95, 0.03, 1.93)
    # descending: 3 (1.95), 5 (1.93), 1 (1.9), then the near-empty three
    space = partition_by_occupation(occupations, n_cas=2, n_active_orbitals=4)
    # k = (6 - 2)/2 = 2 closed = the two largest occupations; 0 virtual left
    assert space.closed == (3, 5)
    assert space.active == (0, 1, 2, 4)
    assert space.virtual == ()
    # a virtual-containing request: k = 1 closed, four active, one virtual
    mixed = partition_by_occupation(occupations, n_cas=4, n_active_orbitals=4)
    assert mixed.closed == (3,)
    assert mixed.active == (1, 2, 4, 5)
    assert mixed.virtual == (0,)
    with pytest.raises(QicasError, match="must be even"):
        partition_by_occupation(occupations, n_cas=3, n_active_orbitals=2)
    with pytest.raises(QicasError, match="does not fit"):
        partition_by_occupation(occupations, n_cas=2, n_active_orbitals=6)


def test_the_accumulated_rotation_replays_the_density_rotations():
    """The composition law of the u[old, new] convention: the matrix returned by
    the optimizer, applied once to the starting densities, must reproduce the
    optimized ones exactly (pins the u = u @ t accumulation)."""
    dens = _dens()
    optimized, report = qicas.minimize_rotations(dens, (0, 2))
    u = np.asarray(report.matrix)
    assert report.accepted >= 1
    replayed = entropy_rdm.rotate_densities(dens, u)
    assert np.abs(replayed.gamma_a - optimized.gamma_a).max() < 1e-12
    assert np.abs(replayed.g_ab - optimized.g_ab).max() < 1e-10
    assert np.abs(u @ u.T - np.eye(dens.norb)).max() < 1e-10


def test_the_pair_slice_evaluation_equals_a_full_rotation():
    """The exact O(1) candidate evaluation inside the Jacobi scan, against a full
    density rotation for the same angle."""
    dens = _dens()
    for i, j, theta in ((0, 2, 0.31), (1, 5, 1.2), (3, 0, 2.7)):
        t = qicas._rotation_two(dens.norb, i, j, theta)
        rotated = entropy_rdm.rotate_densities(dens, t)
        fast = qicas._pair_entropies(dens, i, j, theta)
        slow = (
            entropy_rdm.orbital_entropy(rotated, i),
            entropy_rdm.orbital_entropy(rotated, j),
        )
        assert fast == pytest.approx(slow, abs=1e-12)


def test_the_folding_reproduces_the_window_fci_when_nothing_is_closed():
    """CASCI with an empty closed set is the window FCI itself; this pins the
    effective-Hamiltonian folding (core energy + Fock contraction signs)."""
    dump = _dump()
    space = qicas.CSpace(closed=(), active=(0, 1, 2, 3, 4, 5), virtual=())
    energy = qicas.casci_energy(dump, np.eye(6), space)
    assert energy == pytest.approx(REFERENCE, abs=1e-9)


# --- the (4,4) target ----------------------------------------------------------


def test_the_4x4_target_is_already_at_the_fixed_point():
    result = qicas.analyze(
        _dump(), n_cas=4, n_active_orbitals=4, reference_energy=REFERENCE
    )
    assert result.space.closed == (0,)
    assert result.space.active == (1, 2, 3, 4)
    assert result.space.virtual == (5,)
    assert result.rotation.accepted == 0
    assert result.rotation.fqi_initial == pytest.approx(0.030170, abs=1e-4)
    assert result.gap_to_fci == pytest.approx(5.331e-3, abs=1e-5)
    assert result.theorem_bound == pytest.approx(1.034e-1, abs=1e-4)
    assert result.gap_to_fci <= result.theorem_bound
    assert result.space_final == result.space


# --- the (2,4) target: the source's story on real data -------------------------


def test_the_2x4_target_is_repaired_to_the_4x4_space():
    result = qicas.analyze(
        _dump(), n_cas=2, n_active_orbitals=4, reference_energy=REFERENCE
    )
    # the request splits the pi pair (closed {sigma-2p, one pi}); the sweep
    # finds the F_QI fixed point and the occupancy reading reclassifies
    assert result.rotation.fqi_initial == pytest.approx(0.248094, abs=1e-4)
    assert result.rotation.fqi_final == pytest.approx(0.030170, abs=1e-4)
    assert result.rotation.accepted >= 2
    assert result.energy_casci_initial - result.energy_fci == pytest.approx(6.697e-2, abs=1e-4)
    assert result.gap_to_fci == pytest.approx(5.331e-3, abs=1e-5)
    assert result.gap_to_fci <= result.theorem_bound
    assert result.space_final.closed == (2,)
    assert result.space_final.active == (1, 3, 4, 5)
    assert result.space_final.virtual == (0,)


def test_the_exclusive_rotation_set_runs_and_stays_above_the_touch_set():
    """The economical variant (active/non-active pairs only) is a documented
    restriction of the default; on this window it must not beat the touch set's
    final F_QI."""
    touch = qicas.analyze(_dump(), n_cas=2, n_active_orbitals=4, reference_energy=REFERENCE)
    exclusive = qicas.analyze(
        _dump(),
        n_cas=2,
        n_active_orbitals=4,
        reference_energy=REFERENCE,
        pairs_mode="exclusive",
    )
    assert exclusive.rotation.fqi_final >= touch.rotation.fqi_final - 1e-9


# --- output --------------------------------------------------------------------


def test_the_report_prints_the_partitions_and_the_theorem_check():
    result = qicas.analyze(
        _dump(), n_cas=4, n_active_orbitals=4, reference_energy=REFERENCE
    )
    body = qicas.render(result)
    assert "QICAS orbital optimization" in body
    assert "requested partition" in body and "final partition" in body
    assert "Theorem-1 check" in body and "(holds)" in body
    assert "subset application" in body  # the honest window restriction


def test_evidence_carries_the_source():
    bibkeys = {item.bibkey for item in qicas.evidence() if item.bibkey}
    assert bibkeys == {"ding2023qicas"}


# --- the optimized-basis export -------------------------------------------------


def test_the_export_rotates_the_window_in_place():
    result = qicas.analyze(
        _dump(), n_cas=2, n_active_orbitals=4, reference_energy=REFERENCE
    )
    export = parse_orca_json(FIXTURES / "n2_fcidump.canonical.json")
    qe = qicas.quasi_export(result, export)
    assert qe.n_closed == 1 and qe.n_active == 4
    overlap = np.asarray(export.overlap)
    gram = qe.coefficients.T @ overlap @ qe.coefficients
    assert np.abs(gram - np.eye(export.n_mo)).max() < 1e-10
    # the column order is preserved: outside the window the occupations are the
    # export's own, inside it they are exactly the optimized ones, and the
    # rotation did move the window contents
    occ0 = np.asarray(export.mo_occupations, dtype=float)
    assert np.allclose(qe.occupations[:4], occ0[:4])
    assert np.allclose(
        qe.occupations[4:10], np.asarray(result.final_occupations), atol=1e-9
    )
    assert not np.allclose(qe.occupations[4:10], occ0[4:10])


def test_the_export_is_deterministic():
    result = qicas.analyze(
        _dump(), n_cas=2, n_active_orbitals=4, reference_energy=REFERENCE
    )
    export = parse_orca_json(FIXTURES / "n2_fcidump.canonical.json")
    first = qicas.quasi_export(result, export)
    second = qicas.quasi_export(result, export)
    assert np.array_equal(first.coefficients, second.coefficients)
    assert np.array_equal(first.occupations, second.occupations)


def test_the_write_back_round_trips_through_the_mkl(tmp_path):
    from fblockkit.parsers.mkl import parse_mkl
    from fblockkit.recipe import ass1st as ass1st_recipe

    result = qicas.analyze(
        _dump(), n_cas=2, n_active_orbitals=4, reference_energy=REFERENCE
    )
    export = parse_orca_json(FIXTURES / "n2_fcidump.canonical.json")
    qe = qicas.quasi_export(result, export)
    template = parse_mkl(FIXTURES / "n2_scan_1.600.mkl")
    out = tmp_path / "qicas.mkl"
    ass1st_recipe.write_qno_mkl(template, qe.coefficients, qe.occupations, out)
    back = parse_mkl(out)
    coefficients = np.asarray(back.coefficients())
    assert coefficients.shape == qe.coefficients.shape
    assert np.abs(coefficients - qe.coefficients).max() < 1e-6
