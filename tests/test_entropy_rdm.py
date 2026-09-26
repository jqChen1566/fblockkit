"""Exact four-state entropy from an FCIDUMP: solver, densities, rotation.

The reference values are the ones ORCA itself printed for the N2 fixture
(fixtures/orca/n2_fcidump_step_a.out):

- CASSCF energy  -108.950671945279 Eh
- N(occ) = 1.99345 1.93634 1.93634 0.06599 0.06599 0.00190

The module docstring of fblockkit.analysis.entropy_rdm states the contract: a
reconstruction is only trusted when these two engine numbers come back, plus
the internal conventions (energy from the density objects, rotation invariance).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from fblockkit.analysis.entropy_rdm import (
    EntropyRdmError,
    energy_from_densities,
    entropy_spectrum,
    infer_active_window,
    natural_occupations,
    orbital_weights,
    rotate_densities,
    rotation_from_coefficients,
    solve_fci,
    spin_densities,
)
from fblockkit.parsers.fcidump import Fcidump, parse_fcidump

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"
DUMP = FIXTURES / "n2_fcidump.fcidump"

REFERENCE_ENERGY = -108.950671945279
REFERENCE_OCCUPATIONS = (1.99345, 1.93634, 1.93634, 0.06599, 0.06599, 0.00190)
# canonical-basis four-state entropies from the route-validation run (2026-09-26);
# orbitals 1/2 and 3/4 are degenerate pairs, so the test only fixes the multiset
REFERENCE_ENTROPY_MULTISET = (0.0078, 0.0223, 0.2258, 0.2258, 0.2302, 0.2302)


@pytest.fixture(scope="module")
def dump() -> Fcidump:
    return parse_fcidump(DUMP)


@pytest.fixture(scope="module")
def state(dump):
    return solve_fci(dump, reference_energy=REFERENCE_ENERGY, multiplicity=1)


@pytest.fixture(scope="module")
def densities(state):
    return spin_densities(state)


# --- the two engine cross-checks ---------------------------------------------


def test_ci_energy_reproduces_the_casscf_energy(state):
    assert state.energy_active + state.ecore == pytest.approx(REFERENCE_ENERGY, abs=1e-9)


def test_natural_occupations_match_the_printed_n_occ(densities):
    computed = natural_occupations(densities)
    assert computed == pytest.approx(REFERENCE_OCCUPATIONS, abs=1e-4)


def test_spin_state_is_a_singlet(state):
    assert state.s2 == pytest.approx(0.0, abs=1e-8)


def test_selection_paths_agree(dump, state):
    by_multiplicity = solve_fci(dump, multiplicity=1)
    assert by_multiplicity.energy_active == pytest.approx(state.energy_active, abs=1e-10)


# --- internal conventions -----------------------------------------------------


def test_energy_from_densities_reproduces_the_ci_energy(dump, state, densities):
    rebuilt = energy_from_densities(dump, densities)
    assert rebuilt == pytest.approx(state.energy_active, abs=1e-8)


def test_four_state_weights_sum_to_one(densities):
    for i in range(densities.norb):
        assert sum(orbital_weights(densities, i)) == pytest.approx(1.0, abs=1e-8)


def test_entropy_spectrum_matches_the_recorded_values(densities):
    computed = sorted(entropy_spectrum(densities))
    assert computed == pytest.approx(sorted(REFERENCE_ENTROPY_MULTISET), abs=1e-4)


def test_rotation_preserves_the_energy(dump, densities):
    """The strongest convention test: rotate the integrals and the densities
    together -- the energy expression must be invariant (it exercises gamma,
    the four spin blocks and every index convention at once)."""
    rng = np.random.default_rng(7)
    q, _ = np.linalg.qr(rng.standard_normal((dump.norb, dump.norb)))
    rotated_dens = rotate_densities(densities, q)
    n = dump.norb
    h = np.asarray(dump.h)
    h_rot = q.T @ h @ q
    g_dense = np.zeros((n, n, n, n))
    for key, value in dump.g.items():
        g_dense[key] = value
    g_rot = np.einsum("pi,qj,rk,sl,pqrs->ijkl", q, q, q, q, g_dense)
    # the rotated tensor is dense: every index quadruple must be carried over
    # (zero-in-the-old-basis entries are not zero after the rotation)
    rotated_dump = Fcidump(
        norb=n,
        nelec=dump.nelec,
        ms2=dump.ms2,
        ecore=dump.ecore,
        h=tuple(tuple(row) for row in h_rot),
        g={
            (i, j, k, l): float(g_rot[i, j, k, l])
            for i in range(n)
            for j in range(n)
            for k in range(n)
            for l in range(n)
        },
        orbsym=dump.orbsym,
        isym=dump.isym,
    )
    assert energy_from_densities(rotated_dump, rotated_dens) == pytest.approx(
        energy_from_densities(dump, densities), abs=1e-8
    )


def test_identity_rotation_is_a_no_op(densities):
    same = rotate_densities(densities, np.eye(densities.norb))
    assert entropy_spectrum(same) == pytest.approx(entropy_spectrum(densities), abs=1e-12)


def test_rotation_from_coefficients_is_orthogonal(dump):
    """The rotation built from the two real exports (canonical and IAO-IBO
    localised) is an orthogonal matrix -- validated against the coefficients and
    the AO overlap in the orca_2json fixtures."""
    canonical = json.loads((FIXTURES / "n2_fcidump.canonical.json").read_text())
    localized = json.loads((FIXTURES / "n2_fcidump.localized.json").read_text())
    c_can = np.array([mo["MOCoefficients"] for mo in canonical["Molecule"]["MolecularOrbitals"]["MOs"]]).T
    c_loc = np.array([mo["MOCoefficients"] for mo in localized["Molecule"]["MolecularOrbitals"]["MOs"]]).T
    s = np.array(canonical["Molecule"]["S-Matrix"])
    u = rotation_from_coefficients(c_can, c_loc, s, active=[4, 5, 6, 7, 8, 9])
    assert u.T @ u == pytest.approx(np.eye(6), abs=1e-10)


def test_localized_spectrum_uses_the_same_weights(dump, densities):
    """Rotating with the real localisation gives a spectrum of the same
    multiset-only character: weights stay in [0, 1] and each orbital's entropy
    stays below ln 4."""
    canonical = json.loads((FIXTURES / "n2_fcidump.canonical.json").read_text())
    localized = json.loads((FIXTURES / "n2_fcidump.localized.json").read_text())
    c_can = np.array([mo["MOCoefficients"] for mo in canonical["Molecule"]["MolecularOrbitals"]["MOs"]]).T
    c_loc = np.array([mo["MOCoefficients"] for mo in localized["Molecule"]["MolecularOrbitals"]["MOs"]]).T
    s = np.array(canonical["Molecule"]["S-Matrix"])
    u = rotation_from_coefficients(c_can, c_loc, s, active=[4, 5, 6, 7, 8, 9])
    loc_dens = rotate_densities(densities, u)
    spectrum = entropy_spectrum(loc_dens)
    assert len(spectrum) == 6
    for i, value in enumerate(spectrum):
        assert 0.0 <= value <= np.log(4.0) + 1e-12, i
        assert sum(orbital_weights(loc_dens, i)) == pytest.approx(1.0, abs=1e-8)
    # the total correlation content is basis-invariant: sum of occupations is N
    assert float(np.trace(loc_dens.gamma_a + loc_dens.gamma_b)) == pytest.approx(6.0, abs=1e-8)
    # and the localised basis must not leave the entropy spectrum identical to
    # the canonical one (IAO-IBO mixes the bonding/antibonding pairs)
    assert spectrum != pytest.approx(entropy_spectrum(densities), abs=1e-6)


# --- degenerate limits and refusals ------------------------------------------


def _tiny_dump(tmp_path, text: str) -> Fcidump:
    path = tmp_path / "FCIDUMP"
    path.write_text(text, encoding="utf-8")
    return parse_fcidump(path)


def test_single_determinant_has_zero_entropy(tmp_path):
    """A closed-shell determinant: filling the lower orbital twice gives zero
    entropy on every orbital (the classic sanity limit)."""
    dump = _tiny_dump(
        tmp_path,
        "&FCI NORB=2,NELEC=2,MS2=0\n/\n"
        "  -1.0 1 1 0 0\n   0.5 2 2 0 0\n   0.0 0 0 0 0\n",
    )
    state = solve_fci(dump)
    densities = spin_densities(state)
    assert entropy_spectrum(densities) == pytest.approx((0.0, 0.0), abs=1e-10)
    assert natural_occupations(densities) == pytest.approx((2.0, 0.0), abs=1e-10)


def test_reference_energy_mismatch_refuses(dump):
    with pytest.raises(EntropyRdmError, match="Next step"):
        solve_fci(dump, reference_energy=-100.0)


def test_impossible_multiplicity_refuses(dump):
    with pytest.raises(EntropyRdmError, match="Next step"):
        solve_fci(dump, multiplicity=9)


def test_determinant_cap_refuses(dump):
    with pytest.raises(EntropyRdmError, match="cap"):
        solve_fci(dump, reference_energy=REFERENCE_ENERGY, max_determinants=100)


# --- orchestration: analysis record, window inference, report section ---------


def test_infer_active_window_plain_span():
    occupations = [2.0] * 4 + list(REFERENCE_OCCUPATIONS) + [0.0] * 18
    window, _ = infer_active_window(occupations, 6)
    assert window == (4, 5, 6, 7, 8, 9)


def test_infer_active_window_with_empty_active_orbital():
    """The f-block shape: one active orbital prints occupation 0.0000; the
    extension must grow away from the 2.0 core neighbour."""
    occupations = [2.0] * 37 + [1.0] + [0.0] * 6 + [0.0] * 20
    window, note = infer_active_window(occupations, 7)
    assert window == tuple(range(37, 44))
    assert "extended upwards" in note


def test_infer_active_window_with_interior_zero():
    occupations = [2.0] * 4 + [1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0] + [0.0] * 20
    window, _ = infer_active_window(occupations, 7)
    assert window == tuple(range(4, 11))


def test_infer_active_window_refuses_on_count_mismatch():
    occupations = [2.0] * 4 + list(REFERENCE_OCCUPATIONS) + [0.0] * 18
    with pytest.raises(EntropyRdmError, match="Next step"):
        # 28 orbitals in the table can never form a 30-orbital window
        infer_active_window(occupations, 30)


def test_infer_active_window_validates_against_the_printed_n_occ():
    """A shifted window of the right size must be caught: the occupations of the
    chosen orbitals are compared with the N(occ)= line."""
    correct = [2.0] * 4 + list(REFERENCE_OCCUPATIONS) + [0.0] * 18
    window, _ = infer_active_window(
        correct, 6, reference_occupations=REFERENCE_OCCUPATIONS
    )
    assert window == (4, 5, 6, 7, 8, 9)
    shifted = (
        [2.0] * 5 + [1.93634, 1.93634, 0.06599, 0.06599, 0.00190, 0.0] + [0.0] * 17
    )
    with pytest.raises(EntropyRdmError, match="N\\(occ\\)"):
        infer_active_window(shifted, 6, reference_occupations=REFERENCE_OCCUPATIONS)


def _exports():
    from fblockkit.parsers.orca_json import parse_orca_json

    return (
        parse_orca_json(FIXTURES / "n2_fcidump.canonical.json"),
        parse_orca_json(FIXTURES / "n2_fcidump.localized.json"),
    )


def test_analyze_collects_the_cross_checks_and_both_spectra(dump):
    from fblockkit.analysis.entropy_rdm import analyze

    canonical, localized = _exports()
    analysis = analyze(
        dump,
        reference_energy=REFERENCE_ENERGY,
        multiplicity=1,
        reference_occupations=REFERENCE_OCCUPATIONS,
        canonical=canonical,
        localized=localized,
        active_window=(4, 5, 6, 7, 8, 9),
    )
    assert analysis.occupation_deviation is not None
    assert analysis.occupation_deviation < 1e-4
    assert len(analysis.checks) >= 3
    assert sorted(analysis.canonical_spectrum) == pytest.approx(
        sorted(REFERENCE_ENTROPY_MULTISET), abs=1e-4
    )
    assert analysis.localized_spectrum is not None
    assert len(analysis.localized_spectrum) == 6
    assert analysis.active_window == (4, 5, 6, 7, 8, 9)
    for weights in analysis.localized_weights or ():
        assert sum(weights) == pytest.approx(1.0, abs=1e-8)


def test_analyze_refuses_a_lone_export(dump):
    from fblockkit.analysis.entropy_rdm import analyze

    canonical, _ = _exports()
    with pytest.raises(EntropyRdmError, match="BOTH"):
        analyze(
            dump,
            reference_energy=REFERENCE_ENERGY,
            canonical=canonical,
            active_window=(4, 5, 6, 7, 8, 9),
        )


def test_run_section_states_the_route_and_the_checks(dump):
    from fblockkit.analysis.entropy_rdm import analyze, run as entropy_run

    canonical, localized = _exports()
    analysis = analyze(
        dump,
        reference_energy=REFERENCE_ENERGY,
        multiplicity=1,
        reference_occupations=REFERENCE_OCCUPATIONS,
        canonical=canonical,
        localized=localized,
        active_window=(4, 5, 6, 7, 8, 9),
    )
    section = entropy_run(analysis)
    assert "FCIDUMP route" in section.title
    assert "Cross-checks against the engine" in section.body
    assert "-108.950671945" in section.body
    assert "IAO-IBO" in section.body
    assert "0.14" in section.body
    assert "Plateau" in section.body or "plateau" in section.body


def test_evidence_carries_the_three_provenance_kinds():
    from fblockkit.analysis.entropy_rdm import evidence

    kinds = {item.kind for item in evidence()}
    assert kinds == {"literature", "manual", "measured"}
    for item in evidence():
        assert item.text and item.ref


# --- the f-block chain (Eu3+) -------------------------------------------------


EU_DUMP = FIXTURES / "eu3_fcidump.fcidump"
EU_REFERENCE_ENERGY = -10826.615512593075


def test_eu3_high_spin_chain():
    """The f-block fixture: Eu3+ 4f6 in the MS2=6 sector (7 determinants).
    The maximal-weight 7F state is a single determinant, so the reconstruction
    must reproduce the printed energy, give <S^2> = 12 exactly and a zero
    entropy spectrum -- a clean negative control next to the N2 positive."""
    dump = parse_fcidump(EU_DUMP)
    assert (dump.norb, dump.nelec, dump.ms2) == (7, 6, 6)
    state = solve_fci(dump, reference_energy=EU_REFERENCE_ENERGY, multiplicity=7)
    assert state.energy_active + state.ecore == pytest.approx(
        EU_REFERENCE_ENERGY, abs=1e-9
    )
    assert len(state.determinants) == 7
    assert state.s2 == pytest.approx(12.0, abs=1e-8)
    densities = spin_densities(state)
    assert natural_occupations(densities) == pytest.approx(
        (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0), abs=1e-8
    )
    assert entropy_spectrum(densities) == pytest.approx((0.0,) * 7, abs=1e-10)
