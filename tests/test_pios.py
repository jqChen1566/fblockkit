"""Tests for the PiOS pi-orbital active space (analysis/pios.py).

The real-data anchors are the RHF/cc-pVDZ benzene export
(fixtures/orca/benzene_rhf.json, the S/H/J/K chain measured on 6.1.1): the pi
plane is the molecular plane, the six carbons contribute one electron each
(CAS(6e, 6o)), the occupied projection spectrum selects 0.7789 / 0.7649 /
0.7649 with the next eigenvalue excluded at 0.000 and the virtual side
1.000 / 1.000 / 1.000 (next 0.235); the projector trace closes exactly at six
(2.309 + 3.691); and the written orbitals' Fock energies reproduce the
canonical occupied/virtual energies as sets (the menu-21 cross-check).  The
CASSCF(6,6) engine anchors live in fixtures/orca/benzene_pios*.out.
"""

from __future__ import annotations

import numpy as np
import pytest

from dataclasses import replace
from pathlib import Path

from fblockkit.analysis import pios
from fblockkit.parsers.mkl import parse_mkl
from fblockkit.parsers.orca_json import parse_orca_json

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"
ATOMS = [0, 1, 2, 3, 4, 5]


def _export():
    return parse_orca_json(FIXTURES / "benzene_rhf.json")


def _select():
    return pios.select_pi_space(_export(), ATOMS)


def test_the_plane_and_the_atom_rules():
    result = _select()
    assert np.allclose(result.normal, (0.0, 0.0, 1.0), atol=1e-8)
    assert result.max_out_of_plane < 1e-8
    assert result.elements == ("C",) * 6
    assert result.degrees == (3,) * 6
    assert result.contributions == (1,) * 6
    assert result.pi_electrons == 6
    assert result.n_occupied_pi == 3 and result.n_virtual_pi == 3
    assert result.charge == 0


def test_the_projection_spectra_and_the_trace_invariant():
    """The selection spectra, and the projector trace closing at the dimension
    of the p_z space (the internal consistency of the whole construction)."""
    result = _select()
    assert result.occ_spectrum[:4] == pytest.approx(
        (0.778934, 0.764871, 0.764871, 0.0), abs=1e-5
    )
    assert result.vir_spectrum[:5] == pytest.approx(
        (1.0, 1.0, 1.0, 0.235129, 0.235129), abs=1e-5
    )
    total = sum(result.occ_spectrum) + sum(result.vir_spectrum)
    assert total == pytest.approx(6.0, abs=1e-6)
    # the selection gaps are clean on both sides (the source's validity reading)
    assert result.occ_spectrum[2] > 0.75 > 0.01 > result.occ_spectrum[3]
    assert result.vir_spectrum[2] > 0.99 > 0.5 > result.vir_spectrum[3]


def test_the_selected_energies_and_the_fock_cross_check():
    """The semicanonical energies of the selected blocks (anchors), the Fock
    assembly checked against the canonical orbital energies, and the trace
    invariant: the written set is complete and S-orthonormal, so its block
    eigenvalues sum to the canonical energy sum."""
    export = _export()
    result = _select()
    assert result.occ_energies == pytest.approx((-0.498453, -0.333508, -0.333508), abs=1e-5)
    assert result.vir_energies == pytest.approx((0.281786, 0.281786, 0.470139), abs=1e-5)
    # the Fock assembly: f = H + J + K reproduces the canonical energies on the
    # diagonal of the MO-basis representation (measured 1.1e-7)
    f = pios._fock_matrix(export)
    coefficients = np.asarray(export.mo_coefficients).T
    overlap = np.asarray(export.overlap)
    f_mo = coefficients.T @ f @ coefficients
    assert np.abs(np.diag(f_mo) - np.asarray(export.mo_energies)).max() < 1e-6
    # the trace invariant: the written block spans the occupied (virtual) space,
    # so its Fock Rayleigh sum equals the canonical occupied (virtual) energy sum
    n_occ = sum(1 for value in export.mo_occupations if value > 1.99)
    assert sum(result.energies[:n_occ]) == pytest.approx(
        sum(export.mo_energies[:n_occ]), abs=1e-5
    )
    assert sum(result.energies[n_occ:]) == pytest.approx(
        sum(export.mo_energies[n_occ:]), abs=1e-5
    )
    assert result.residual < 1e-10


def test_the_parent_weights_follow_the_canonical_orbitals():
    result = _select()
    # the a2u member is uniquely pinned; the e1g pair is degenerate and its two
    # members may swap (either is equally valid -- the pair spans the same space)
    assert [index for index, _ in result.occ_parents[0]] == [16, 17]
    assert result.occ_parents[0][0][1] == pytest.approx(1.0, abs=1e-6)
    assert {result.occ_parents[1][0][0], result.occ_parents[2][0][0]} == {19, 20}


def test_the_written_mkl_carries_the_pi_space(tmp_path):
    result = _select()
    template_mkl = parse_mkl(FIXTURES / "benzene_rhf.mkl")
    out = tmp_path / "benzene_pios.fbk.mkl"
    pios.write_mkl(template_mkl, result, out)
    again = parse_mkl(out)
    assert np.abs(np.array(again.coefficients()) - np.array(result.coefficients)).max() < 1e-7
    assert again.occupations == pytest.approx(result.occupations)


def test_the_element_rules():
    assert pios._element_contribution("C", 2) == 1
    assert pios._element_contribution("C", 4) == 1
    assert pios._element_contribution("N", 2) == 1
    assert pios._element_contribution("N", 3) == 2
    assert pios._element_contribution("P", 3) == 2
    with pytest.raises(pios.PiosError, match="two .* and three"):
        pios._element_contribution("N", 1)
    with pytest.raises(pios.PiosError, match="covers C, N and P"):
        pios._element_contribution("S", 2)


def test_the_refusals():
    export = _export()
    with pytest.raises(pios.PiosError, match="S-Matrix"):
        pios.select_pi_space(replace(export, overlap=None), ATOMS)
    with pytest.raises(pios.PiosError, match="at least three"):
        pios.select_pi_space(export, [0, 1])
    with pytest.raises(pios.PiosError, match="outside"):
        pios.select_pi_space(export, [0, 1, 2, 3, 4, 99])
    with pytest.raises(pios.PiosError, match="repeats"):
        pios.select_pi_space(export, [0, 0, 1, 2, 3, 4])
    with pytest.raises(pios.PiosError, match="H-Matrix"):
        pios.select_pi_space(replace(export, hamiltonian=None), ATOMS)
    with pytest.raises(pios.PiosError, match="one contribution per atom"):
        pios.select_pi_space(export, ATOMS, contributions=[1, 1, 1])
    with pytest.raises(pios.PiosError, match="even number"):
        pios.select_pi_space(export, ATOMS, pi_electrons=7)
    with pytest.raises(pios.PiosError, match="even number"):
        pios.select_pi_space(export, ATOMS, pi_electrons=14)
    # a fractional-occupation (CASSCF) export is refused with the next step
    casscf = parse_orca_json(FIXTURES / "benzene.json")
    with pytest.raises(pios.PiosError, match="fractional occupations"):
        pios.select_pi_space(casscf, ATOMS)
    # a lifted atom leaves the planar regime (the fitted-plane deviation of a
    # single lifted atom is ~0.47x the lift: 1.2 Angstrom reads ~0.57)
    coordinates = [list(row) for row in export.coordinates]
    coordinates[2][2] = 1.2
    lifted = replace(export, coordinates=tuple(tuple(row) for row in coordinates))
    with pytest.raises(pios.PiosError, match="not approximately planar"):
        pios.select_pi_space(lifted, ATOMS)


def test_the_engine_validation_anchors():
    """The frozen CASSCF(6,6) run started from the written pi space: same
    solution as the aufbau fixture, and the SVD between the two active spaces
    -- occupied side essentially exact, virtual side carrying exactly the
    converged pi*'s own out-of-plane character."""
    out = (FIXTURES / "benzene_pios.out").read_text(encoding="utf-8")
    assert out.count("MACRO-ITERATION") == 13
    assert "-230.793818903" in out
    assert "delocalized" not in out
    converged = parse_orca_json(FIXTURES / "benzene_pios.json")
    active = [
        index
        for index, value in enumerate(converged.mo_occupations)
        if 0.001 < value < 1.999
    ]
    assert active == [18, 19, 20, 21, 22, 23]
    guess = np.asarray(_select().coefficients)[:, 18:24]
    coefficients = np.asarray(converged.mo_coefficients).T
    overlap = np.asarray(_export().overlap)
    singular = np.linalg.svd(guess.T @ overlap @ coefficients[:, active], compute_uv=False)
    assert singular == pytest.approx(
        (0.9999, 0.9999, 0.9999, 0.7573, 0.7573, 0.6558), abs=1e-4
    )


def test_the_report_and_evidence_state_the_protocol_and_the_route():
    body = pios.render(_select())
    assert "PiOS pi-orbital active space" in body
    assert "CAS(6e, 6o)" in body
    assert "0.778934" in body and "1.000000" in body
    assert "valence p shell" in body
    assert "moread" in body
    text = " ".join(
        entry.ref + " " + entry.text + " " + entry.url + " " + entry.bibkey
        for entry in pios.evidence()
    )
    assert "8b01196" in text
    assert "sayfutyarova2019constructing" in text
    assert "0.7789" in text
