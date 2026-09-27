"""Checks of the perturbed multistart batch (recipe/perturb_guess.py, menu 29).

The algorithm is the source's section IV.C (Vaucher & Reiher, JCTC 2017):
random occupied-virtual pair mixings, Eqs. (2)-(3), of a converged reference's
orbitals, written back through the measured mkl route.  The fixture chain is
engine-validated on the CH4 dissociation trap (UKS PBE0/def2-SVP, one C-H at
2.6 Angstrom): the guess propagated from the equilibrium orbitals keeps the
restricted solution (-40.111230225721 Eh), ORCA's default guess lands on
another restricted one (-40.183584827774 Eh), and the three menu-generated
perturbed starts all reach the broken-symmetry solution (-40.2416 Eh,
<S**2> 0.971943) -- a 0.130 Eh healing and a further 0.058 Eh below the
default guess.
"""

from __future__ import annotations

import random
from pathlib import Path

import numpy as np
import pytest

from fblockkit.parsers import parse_auto
from fblockkit.parsers.mkl import parse_mkl
from fblockkit.recipe import perturb_guess as pg
from fblockkit.recipe.perturb_guess import PerturbError, PerturbationPlan, PairMix

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"
REF_MKL = FIXTURES / "ch4_diss_prop.mkl"
REF_INPUT = FIXTURES / "inputs" / "ch4_diss_prop.inp"


def test_plan_selects_the_source_windows_and_draws_pairs():
    occupations = [1.0] * 20 + [0.0] * 30
    energies = [float(i) for i in range(50)]
    plan = pg.plan_perturbation(occupations, energies, n_pairs=10, window=15, rng=random.Random(1))
    assert plan.window_occupied == tuple(range(5, 20))
    assert plan.window_virtual == tuple(range(20, 35))
    assert len(plan.pairs) == 10
    for pair in plan.pairs:
        assert pair.occupied in plan.window_occupied
        assert pair.virtual in plan.window_virtual
        assert 0.0 <= pair.angle_deg < 90.0  # uniform in [0, 90)


def test_the_seed_reproduces_every_draw():
    occupations = [1.0, 1.0, 0.0, 0.0, 0.0]
    energies = [0.0, 1.0, 2.0, 3.0, 4.0]
    first = pg.plan_perturbation(occupations, energies, rng=random.Random(7))
    second = pg.plan_perturbation(occupations, energies, rng=random.Random(7))
    third = pg.plan_perturbation(occupations, energies, rng=random.Random(8))
    assert first == second
    assert first != third


def test_plan_refusals_carry_next_steps():
    with pytest.raises(PerturbError, match="at least one"):
        pg.plan_perturbation([1.0, 0.0], [0.0, 1.0], n_pairs=0, rng=random.Random(1))
    with pytest.raises(PerturbError, match="occupied-virtual separation"):
        pg.plan_perturbation([1.0, 1.0], [0.0, 1.0], rng=random.Random(1))
    with pytest.raises(PerturbError, match="disagree"):
        pg.plan_perturbation([1.0, 0.0], [0.0], rng=random.Random(1))


def test_mixing_is_a_column_rotation():
    # an orthonormal 4x4 matrix; a 90-degree pair swaps the columns (with the
    # source's sign convention), a 0-degree pair changes nothing
    identity = tuple(tuple(1.0 if i == j else 0.0 for j in range(4)) for i in range(4))
    swap = PerturbationPlan(
        pairs=(PairMix(occupied=0, virtual=2, angle_deg=90.0),),
        window_occupied=(0,),
        window_virtual=(2,),
    )
    mixed = np.array(pg.mix_orbitals(identity, swap))
    assert np.allclose(mixed[:, 0], [0.0, 0.0, 1.0, 0.0])
    assert np.allclose(mixed[:, 2], [-1.0, 0.0, 0.0, 0.0])
    noop = PerturbationPlan(
        pairs=(PairMix(occupied=0, virtual=2, angle_deg=0.0),),
        window_occupied=(0,),
        window_virtual=(2,),
    )
    assert np.allclose(np.array(pg.mix_orbitals(identity, noop)), identity)
    # every mixing is orthogonal, so the set stays orthonormal
    assert np.allclose(mixed.T @ mixed, np.eye(4))
    with pytest.raises(PerturbError, match="outside"):
        pg.mix_orbitals(identity, PerturbationPlan((PairMix(0, 9, 30.0),), (0,), (9,)))


def test_perturb_mkl_keeps_the_reference_tags_and_reproduces():
    reference = parse_mkl(REF_MKL)
    first = pg.perturb_mkl(reference, seed=20260927, n_pairs=10, window=15)
    second = pg.perturb_mkl(reference, seed=20260927, n_pairs=10, window=15)
    assert first.mkl.render() == second.mkl.render()  # byte-stable
    assert first.mkl.occupations == reference.occupations
    assert first.mkl.beta_occupations == reference.beta_occupations
    assert [g.energies for g in first.mkl.groups] == [g.energies for g in reference.groups]
    assert first.beta is not None  # an unrestricted reference perturbs both spins
    # the coefficients actually changed on both spins
    assert np.abs(
        np.array(first.mkl.coefficients()) - np.array(reference.coefficients())
    ).max() > 1e-3
    assert np.abs(
        np.array(first.mkl.beta_coefficients()) - np.array(reference.beta_coefficients())
    ).max() > 1e-3
    alpha_only = pg.perturb_mkl(reference, seed=20260927, perturb_beta=False)
    assert alpha_only.beta is None
    assert alpha_only.mkl.beta_coefficients() == reference.beta_coefficients()


def test_mo_read_variant_edits_the_input():
    base = (
        "! UKS PBE0 def2-SVP TightSCF\n"
        "%maxcore 2000\n"
        "* xyz 0 1\n"
        "H 0 0 0\n"
        "H 0 0 1.0\n"
        "*\n"
    )
    variant = pg.mo_read_variant(base, "ref.p1.fbk.gbw")
    assert "! UKS PBE0 def2-SVP TightSCF MORead" in variant
    assert '%moinp "ref.p1.fbk.gbw"' in variant
    # an existing %moinp is replaced, not duplicated
    again = pg.mo_read_variant(variant, "ref.p2.fbk.gbw")
    assert again.count("%moinp") == 1
    assert '%moinp "ref.p2.fbk.gbw"' in again
    assert again.count("MORead") == 1  # not doubled
    with pytest.raises(PerturbError, match="no simple"):
        pg.mo_read_variant("%maxcore 2000\n* xyz 0 1\nH 0 0 0\n*\n", "x.gbw")


def test_the_geometry_cross_check():
    reference = parse_mkl(REF_MKL)
    base = REF_INPUT.read_text(encoding="utf-8")
    assert pg.check_geometry_match(base, reference) is None  # agrees in every digit
    shifted = base.replace("2.600000", "2.700000")
    with pytest.raises(PerturbError, match="must run at the reference geometry"):
        pg.check_geometry_match(shifted, reference)
    external = base.split("* xyz")[0] + '* xyzfile 0 1 geom.xyz\n'
    note = pg.check_geometry_match(external, reference)
    assert note is not None and "not cross-checked" in note


def test_the_fixture_chain_matches_the_generated_inputs():
    """The committed starts are exactly what the menu emits (same seed), and
    the run outputs carry the healing: trap < fresh < perturbed, spin flipped."""
    reference = parse_mkl(REF_MKL)
    base = REF_INPUT.read_text(encoding="utf-8")
    for index in (1, 2, 3):
        start = pg.perturb_mkl(reference, seed=20260927 + index - 1, n_pairs=10, window=15)
        stem = f"ch4_diss_prop.p{index}.fbk"
        produced = pg.mo_read_variant(base, stem + ".gbw")
        committed = (FIXTURES / "inputs" / f"ch4_diss_prop.p{index}.inp").read_text(
            encoding="utf-8"
        )
        assert produced == committed, index

    def energy(name: str) -> float:
        return parse_auto(FIXTURES / f"{name}.out").sections["final_energy"]

    trap = energy("ch4_diss_prop")
    fresh = energy("ch4_diss_fresh")
    assert trap == pytest.approx(-40.111230225721, abs=1e-9)
    assert fresh == pytest.approx(-40.183584827774, abs=1e-9)
    assert fresh < trap  # the propagated guess sits above even the default one
    for index in (1, 2, 3):
        healed = energy(f"ch4_diss_prop.p{index}")
        assert healed == pytest.approx(-40.2416020, abs=1e-5)
        assert healed < fresh  # the perturbed starts reach the broken-symmetry solution
