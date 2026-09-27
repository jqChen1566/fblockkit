"""Regression checks of the TNASS selection (analysis/tnass.py).

The fixtures are the two FCIDUMP windows already used by the entropy and
QICAS routes: benzene/cc-pVDZ (the pi sextet) and N2/def2-SVP (CAS(6,6)).
The anchors:

- the exact Renyi-2 oracle reproduces ``-log(sum w^2)`` from the four-state
  weights for every one-orbital subset (both fixtures, 1e-12);
- the bipartition symmetry ``S2(A) = S2(complement)`` holds for every subset
  (a structural check on the mask bookkeeping);
- the greedy trace and the brute-force global maximum are pinned on both
  fixtures -- including the measured suboptimality of greedy (benzene:
  0.23000 vs 0.23280; N2: 0.12663 vs 0.13134) and the seed ranking of the
  block method at k = 1.
"""

from __future__ import annotations

import itertools
import math

import pytest

from fblockkit.analysis import entropy_rdm, tnass
from fblockkit.analysis.tnass import TnassError
from fblockkit.parsers.fcidump import parse_fcidump

from pathlib import Path

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"
BENZENE = FIXTURES / "benzene.fcidump"
N2 = FIXTURES / "n2_fcidump.fcidump"

BENZENE_ENERGY = -230.793818898
N2_ENERGY = -108.950671945


def _state(path: Path, energy: float):
    return entropy_rdm.solve_fci(parse_fcidump(path), reference_energy=energy)


# --- the exact oracle ---------------------------------------------------------


@pytest.mark.parametrize(
    ("path", "energy"),
    [(BENZENE, BENZENE_ENERGY), (N2, N2_ENERGY)],
)
def test_the_one_orbital_s2_reproduces_the_four_state_weights(path, energy):
    state = _state(path, energy)
    densities = entropy_rdm.spin_densities(state)
    for index in range(state.norb):
        purity = sum(w * w for w in entropy_rdm.orbital_weights(densities, index))
        assert tnass.renyi2_entropy(state, (index,)) == pytest.approx(
            -math.log(purity), abs=1e-12
        )


@pytest.mark.parametrize(
    ("path", "energy"),
    [(BENZENE, BENZENE_ENERGY), (N2, N2_ENERGY)],
)
def test_the_bipartition_symmetry_holds_for_every_subset(path, energy):
    """S2(A) = S2(complement): the entanglement of a cut does not depend on which
    side is called A.  A structural check on the occupation-mask bookkeeping."""
    state = _state(path, energy)
    orbitals = set(range(state.norb))
    for size in range(1, state.norb):
        for subset in itertools.combinations(range(state.norb), size):
            complement = tuple(sorted(orbitals - set(subset)))
            assert tnass.renyi2_entropy(state, subset) == pytest.approx(
                tnass.renyi2_entropy(state, complement), abs=1e-12
            )


def test_the_full_window_and_the_empty_subset():
    state = _state(BENZENE, BENZENE_ENERGY)
    assert tnass.renyi2_entropy(state, tuple(range(6))) == 0.0
    with pytest.raises(TnassError, match="empty subset"):
        tnass.renyi2_entropy(state, ())
    with pytest.raises(TnassError, match="outside the window"):
        tnass.renyi2_entropy(state, (9,))


def test_the_benzene_two_orbital_values():
    state = _state(BENZENE, BENZENE_ENERGY)
    assert tnass.renyi2_entropy(state, (0, 1)) == pytest.approx(0.1881207, abs=1e-6)
    assert tnass.renyi2_entropy(state, (0, 5)) == pytest.approx(0.0695097, abs=1e-6)


# --- the selection algorithms -------------------------------------------------


def test_the_benzene_greedy_trace_and_the_brute_optimum():
    dump = parse_fcidump(BENZENE)
    greedy = tnass.analyze(dump, n_target=4, reference_energy=BENZENE_ENERGY)
    assert [(step.added, round(step.s2, 5)) for step in greedy.steps] == [
        (3, 0.14977),
        (4, 0.23280),
        (5, 0.24044),
        (0, 0.23000),
    ]
    assert greedy.selected == (0, 3, 4, 5)
    assert greedy.n_electrons == 2
    brute = tnass.analyze(
        dump, n_target=4, method="brute", reference_energy=BENZENE_ENERGY
    )
    assert brute.selected == (0, 1, 2, 5)
    assert brute.steps[0].s2 == pytest.approx(0.23280, abs=1e-5)
    assert brute.n_electrons == 6
    # the measured suboptimality of the greedy walk on this platform
    assert brute.steps[0].s2 > greedy.steps[-1].s2 + 1e-3


def test_the_block_seed_ranking_is_the_single_orbital_s2_order():
    dump = parse_fcidump(BENZENE)
    result = tnass.analyze(
        dump, n_target=6, method="block", block_size=1, reference_energy=BENZENE_ENERGY
    )
    assert [step.added for step in result.steps] == [3, 4, 2, 1, 0, 5]
    assert result.selected == (0, 1, 2, 3, 4, 5)


def test_the_n2_greedy_and_brute_selections():
    dump = parse_fcidump(N2)
    greedy = tnass.analyze(dump, n_target=4, reference_energy=N2_ENERGY)
    assert greedy.selected == (0, 3, 4, 5)
    assert greedy.steps[-1].s2 == pytest.approx(0.12663, abs=1e-5)
    brute = tnass.analyze(dump, n_target=4, method="brute", reference_energy=N2_ENERGY)
    assert brute.selected == (0, 1, 2, 5)
    assert brute.steps[0].s2 == pytest.approx(0.13134, abs=1e-5)


# --- refusals -----------------------------------------------------------------


def test_the_brute_cap_refuses_with_a_next_step(monkeypatch):
    dump = parse_fcidump(BENZENE)
    monkeypatch.setattr(tnass, "BRUTE_FORCE_CAP", 10)
    with pytest.raises(TnassError, match="above the cap"):
        tnass.analyze(dump, n_target=4, method="brute")


def test_target_and_method_validation():
    dump = parse_fcidump(BENZENE)
    with pytest.raises(TnassError, match="outside 1..6"):
        tnass.analyze(dump, n_target=7, reference_energy=BENZENE_ENERGY)
    with pytest.raises(TnassError, match="unknown method"):
        tnass.analyze(dump, n_target=2, method="annealing")
    with pytest.raises(TnassError, match="needs a block size"):
        tnass.analyze(dump, n_target=2, method="block")


# --- output -------------------------------------------------------------------


def test_the_report_prints_the_trace_and_the_space():
    dump = parse_fcidump(BENZENE)
    body = tnass.render(tnass.analyze(dump, n_target=4, reference_energy=BENZENE_ENERGY))
    assert "TNASS active-space selection" in body
    assert "seed ranking" in body
    assert "selected: orbitals [0, 3, 4, 5] (2e, 4o)" in body
    assert "bond-dimension truncation" in body


def test_evidence_carries_the_source():
    bibkeys = {item.bibkey for item in tnass.evidence() if item.bibkey}
    assert bibkeys == {"mingare2026tnass"}
