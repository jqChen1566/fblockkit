"""Regression checks of the APC ranked-orbital selection (analysis/apc.py).

Three anchors carry this implementation:

1. the source's own CSF table (max(7,6) = 490, max(8,8) = 1764,
   max(10,10) = 19404, max(12,12) = 226512) recomputed from its eq. 2, together
   with the recorded notation ambiguity ((7e, 6o) evaluates to 210; 490 is a
   7-orbital space with 6 or 8 electrons);
2. the H2/STO-3G two-configuration model (fixtures ``h2_apc.json`` and
   ``h2_apc.fcidump``, where the model is the *whole* space): the exported
   exchange diagonal equals the exact (12|12) to 7e-16 (eq. 18 is exact here),
   the FCIDUMP 2x2 CI reproduces the engine's CASSCF energy to 12 digits, the
   eq. 14 closed form equals the eigenvector coefficient to 2e-17, and the
   eq. 19 pair coefficient measures the Koopmans approximation honestly
   (c_APC = -0.072125 against the exact -0.113263, entropy 0.0324 against
   0.0679 -- this STO-3G case is the model's worst corner, not its typical
   accuracy);
3. the N2/def2-SVP export (``n2_apc.json``): the exchange identity
   -diag(C K C^T)[a] = sum_i (a i | a i) holds to 1e-9 against the MO_IAJB
   entries, the Fock diagonal reproduces the orbital energies to 1e-8, the
   pi/pi* quartet ranks first, and the max(10,10) selection lands exactly on
   the source's (10,10) = 19404 table entry.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pytest

from fblockkit.analysis import apc
from fblockkit.analysis.apc import ApcError, n_csf
from fblockkit.parsers import parse_auto
from fblockkit.parsers.fcidump import parse_fcidump
from fblockkit.parsers.orca_json import parse_orca_json

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"
N2 = FIXTURES / "n2_apc.json"
H2 = FIXTURES / "h2_apc.json"
H2_DUMP = FIXTURES / "h2_apc.fcidump"
H2_CASSCF_OUT = FIXTURES / "h2_apc_casscf_a.out"


def _variant(tmp_path: Path, mutate: Callable[[dict[str, Any]], None], source: Path = N2) -> Path:
    document = json.loads(source.read_text(encoding="utf-8"))
    mutate(document["Molecule"])
    target = tmp_path / "variant.json"
    target.write_text(json.dumps(document), encoding="utf-8")
    return target


def _exchange_diagonal(export) -> np.ndarray:
    """-diag(C K C^T): the source's 0.5 K_aa in ORCA's export sign."""
    coefficients = np.array(export.mo_coefficients)
    exchange = np.array(export.exchange[0])
    return -np.diag(coefficients @ exchange @ coefficients.T)


def _binary_entropy(p: float) -> float:
    return -(p * math.log(p) + (1.0 - p) * math.log(1.0 - p))


# --- eq. 2 and the presets ----------------------------------------------------


def test_the_sources_csf_table_is_reproduced():
    """The four catalogue values recomputed from eq. 2, plus the notation record:
    the source's sentence writes max(N_elec, N_orbs) while its own 490 is a
    7-orbital, 6- or 8-electron space ((7e, 6o) evaluates to 210)."""
    assert n_csf(6, 7) == 490
    assert n_csf(8, 7) == 490
    assert n_csf(7, 6) == 210  # the literal (N_elec, N_orbs) reading, for the record
    assert n_csf(8, 8) == 1764
    assert n_csf(10, 10) == 19404
    assert n_csf(12, 12) == 226512
    assert n_csf(10, 7) == 196  # the shape N2 selects at the max(7,6) level
    assert apc.CSF_CAPS == {
        "max(7,6)": 490,
        "max(8,8)": 1764,
        "max(10,10)": 19404,
        "max(12,12)": 226512,
    }


# --- the H2 two-configuration anchor -----------------------------------------


def test_eq_18_is_exact_on_the_two_configuration_model():
    """0.5 K_aa equals the exact (12|12) here: the exported exchange diagonal and
    the MO_IAJB entry (and the FCIDUMP integral) agree to 1e-12."""
    export = parse_orca_json(H2)
    x = _exchange_diagonal(export)
    assert x[1] == pytest.approx(0.18121046221962572, abs=1e-12)
    assert export.mo_iajb is not None
    assert x[1] == pytest.approx(export.mo_iajb[0][4], abs=1e-12)
    dump = parse_fcidump(H2_DUMP)
    assert x[1] == pytest.approx(dump.g[(0, 1, 0, 1)], abs=1e-12)


def test_the_exact_model_chain_reproduces_the_engine_energy():
    """The FCIDUMP 2x2 CI (with the core energy) reproduces the CASSCF energy the
    engine printed, the eq. 14 closed form reproduces the CI eigenvector
    coefficient, and the eq. 19 pair coefficient from the module matches the
    measured value."""
    dump = parse_fcidump(H2_DUMP)
    e_20 = 2.0 * dump.h[0][0] + dump.g[(0, 0, 0, 0)]
    e_02 = 2.0 * dump.h[1][1] + dump.g[(1, 1, 1, 1)]
    coupling = dump.g[(0, 1, 0, 1)]
    eigenvalues, eigenvectors = np.linalg.eigh(
        np.array([[e_20, coupling], [coupling, e_02]])
    )
    energy = eigenvalues[0] + dump.ecore
    printed = parse_auto(H2_CASSCF_OUT).sections["final_energy"]
    assert printed is not None
    assert energy == pytest.approx(printed, abs=1e-10)

    c_exact = eigenvectors[1, 0] / eigenvectors[0, 0]
    delta = (e_02 - e_20) / 2.0
    c_closed = -coupling / (delta + math.hypot(delta, coupling))
    assert c_exact == pytest.approx(c_closed, abs=1e-12)

    ((i, a, c_apc),) = apc.pair_coefficients(parse_orca_json(H2), window_size=1)
    assert (i, a) == (0, 1)
    assert c_apc == pytest.approx(-0.0721245883402186, abs=1e-9)
    # the honest record of the approximation on this extreme model: the exact
    # entropy from c_exact against the APC entropy from c_apc
    s_exact = _binary_entropy(c_exact**2 / (1.0 + c_exact**2))
    report = apc.analyze(parse_orca_json(H2), cap=apc.CSF_CAPS["max(7,6)"], cap_label="max(7,6)")
    assert s_exact == pytest.approx(0.0679216490, abs=1e-9)
    assert report.ranking.candidates[0].entropy == pytest.approx(0.0324025405, abs=1e-9)
    assert report.selection.members == (0, 1)
    assert report.selection.n_csfs == 3


# --- the N2 export ------------------------------------------------------------


def test_the_exchange_identity_holds_against_the_window_entries():
    """-diag(C K C^T)[a] = sum_i (a i | a i) over the doubly occupied i, for every
    virtual of the exported window (measured worst case 4e-12)."""
    export = parse_orca_json(N2)
    x = _exchange_diagonal(export)
    ranking = apc.rank_orbitals(export, window_size=10)
    assert ranking.window == (7, 16)
    assert ranking.occupied == (0, 1, 2, 3, 4, 5, 6)
    assert export.mo_iajb is not None
    diagonal = {(i, a): value for i, j, a, b, value in export.mo_iajb if i == j and a == b}
    for a in ranking.virtuals:
        total = sum(diagonal[(i, a)] for i in ranking.occupied)
        assert x[a] == pytest.approx(total, abs=1e-9)


def test_the_fock_diagonal_reproduces_the_orbital_energies():
    """The full Fock in the export's own basis is H + J + K (the exported F block
    equals J + K and excludes H); its MO diagonal reproduces the canonical
    orbital energies, so the two delta sources are interchangeable on a
    canonical export."""
    export = parse_orca_json(N2)
    coefficients = np.array(export.mo_coefficients)
    fock = (
        np.array(export.hamiltonian)
        + np.array(export.coulomb[0])
        + np.array(export.exchange[0])
    )
    in_mo = coefficients @ fock @ coefficients.T
    assert np.abs(np.diag(in_mo) - np.array(export.mo_energies)).max() < 1e-8
    by_energies = apc.rank_orbitals(export)
    by_fock = apc.rank_orbitals(export, delta="fock")
    entropies_energies = {c.index: c.entropy for c in by_energies.candidates}
    entropies_fock = {c.index: c.entropy for c in by_fock.candidates}
    assert entropies_energies.keys() == entropies_fock.keys()
    worst = max(
        abs(entropies_energies[index] - entropies_fock[index])
        for index in entropies_energies
    )
    # Measured 1.5e-11: the two sources differ by the H + J + K reconstruction
    # precision only.  (Near-degenerate virtual pairs may swap in the rank
    # *order* at that scale; the entropies per orbital are what is asserted.)
    assert worst < 1e-8


def test_the_ranking_puts_the_pi_system_first():
    """N2's strongest correlation is the pi/pi* quartet; both variants rank it in the
    top six, and the exact-integral variant puts it first outright."""
    export = parse_orca_json(N2)
    apc_ranking = apc.rank_orbitals(export, window_size=10)
    assert {candidate.index for candidate in apc_ranking.candidates[:4]} == {5, 6, 7, 8}
    # the two degenerate pairs order by index (the ranking-resolution rule)
    assert [candidate.index for candidate in apc_ranking.candidates[:2]] == [7, 8]
    assert apc_ranking.candidates[0].entropy == pytest.approx(0.3816, abs=5e-4)
    apcx_ranking = apc.rank_orbitals(export, variant="APCX", window_size=10)
    assert [candidate.index for candidate in apcx_ranking.candidates[:2]] == [7, 8]
    assert {candidate.index for candidate in apcx_ranking.candidates[:4]} == {5, 6, 7, 8}


def test_the_two_variants_agree_on_the_top_set_and_the_magnitude_direction():
    """The source's documented behaviour: APC (diagonal-sum approximation)
    overestimates relative to the exact integrals of APCX, while the top-six
    *set* is the same -- the ranking is what the selection consumes."""
    export = parse_orca_json(N2)
    apc_ranking = apc.rank_orbitals(export, window_size=10)
    apcx_ranking = apc.rank_orbitals(export, variant="APCX", window_size=10)
    assert {c.index for c in apc_ranking.candidates[:6]} == {
        c.index for c in apcx_ranking.candidates[:6]
    } == {3, 4, 5, 6, 7, 8}
    by_index = {candidate.index: candidate.entropy for candidate in apcx_ranking.candidates}
    for candidate in apc_ranking.candidates[:6]:
        assert candidate.entropy > by_index[candidate.index] > 0.0


def test_the_selection_reproduces_the_source_levels():
    """The four preset caps on N2: (10e, 7o) at max(7,6) -- one of the shapes the
    source reports for that level -- and the max(10,10) space with the source's
    own N_CSF value, 19404."""
    export = parse_orca_json(N2)
    expected = {
        "max(7,6)": ((10, 7), 196, (2, 3, 4, 5, 6, 7, 8)),
        "max(8,8)": ((10, 8), 1176, (2, 3, 4, 5, 6, 7, 8, 17)),
        "max(10,10)": ((10, 10), 19404, (2, 3, 4, 5, 6, 7, 8, 9, 10, 17)),
        "max(12,12)": ((10, 12), 169884, (2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 17)),
    }
    for label, (shape, csfs, members) in expected.items():
        selection = apc.analyze(
            export, cap=apc.CSF_CAPS[label], cap_label=label
        ).selection
        assert (selection.n_electrons, selection.n_orbitals) == shape
        assert selection.n_csfs == csfs
        assert selection.members == members
        assert selection.cap_reached is True


def test_degenerate_entropies_rank_and_drop_by_index():
    """Entropies agreeing to 1e-12 count as equal and order by index: the rule that
    keeps the rank order, the drop order and the selection independent of the BLAS
    the K-matrix transform ran under.  The synthetic pair differs by 5e-16, with the
    *higher* index holding the smaller value -- so a plain ascending sort would drop
    index 2, and only the tie rule drops index 1."""
    ranking = apc.Ranking(
        base_name="synthetic",
        variant="APC",
        delta_source="orbital energies",
        occupied=(0,),
        virtuals=(1, 2, 3),
        candidates=(
            apc.Candidate(index=2, role="virtual", entropy=0.25),
            apc.Candidate(index=1, role="virtual", entropy=0.25 + 5e-16),
            apc.Candidate(index=3, role="virtual", entropy=0.30),
            apc.Candidate(index=0, role="doubly occupied", entropy=0.5),
        ),
        window=(1, 3),
        n_virtual_available=3,
    )
    selection = apc.select(ranking, cap=6)
    # (2e, 4o) = 10 CSFs > 6: the first drop is the tie pair's lower index
    assert selection.dropped == (1,)
    assert selection.members == (0, 2, 3)
    assert selection.n_csfs == 6
    assert selection.cap_reached is True


def test_the_space_floor_blocks_an_impossible_drop():
    """H2 has one occupied and one virtual orbital; a cap below its 3 CSFs cannot be
    reached without breaking the floor, and the selection says so instead of
    returning an unreasonable space."""
    report = apc.analyze(parse_orca_json(H2), cap=2)
    assert report.selection.cap_reached is False
    assert report.selection.n_csfs == 3
    assert set(report.selection.floor_skips) == {0, 1}
    assert report.selection.members == (0, 1)


# --- refusals -----------------------------------------------------------------


def test_a_missing_exchange_block_names_the_exact_request(tmp_path):
    variant = _variant(tmp_path, lambda molecule: molecule.pop("K-Matrix"))
    with pytest.raises(ApcError, match=r'FockMatrix.*\["K"\]') as excinfo:
        apc.rank_orbitals(parse_orca_json(variant))
    assert "K-Matrix" in str(excinfo.value)


def test_apcx_without_the_window_block_names_the_eight_integer_line(tmp_path):
    def mutate(molecule: dict[str, Any]) -> None:
        molecule.pop("2elIntegrals")

    variant = _variant(tmp_path, mutate)
    with pytest.raises(ApcError, match="eight integers") as excinfo:
        apc.rank_orbitals(parse_orca_json(variant), variant="APCX")
    assert "OrbWin" in str(excinfo.value)


def test_apcx_beyond_the_exported_window_is_refused_with_the_corrected_line():
    """The export covers 10 virtuals; asking APCX for the default 23 must refuse and
    name the eight-integer OrbWin that covers the request."""
    export = parse_orca_json(N2)
    with pytest.raises(ApcError, match=r"OrbWin.*\[0, 6, 7, 27, 0, 0, 0, 0\]"):
        apc.rank_orbitals(export, variant="APCX")


def test_an_open_shell_export_is_refused(tmp_path):
    variant = _variant(tmp_path, lambda molecule: molecule.__setitem__("HFTyp", "ROHF"))
    with pytest.raises(ApcError, match="closed-shell"):
        apc.rank_orbitals(parse_orca_json(variant))


def test_a_non_integer_occupation_is_refused(tmp_path):
    def mutate(molecule: dict[str, Any]) -> None:
        molecule["MolecularOrbitals"]["MOs"][5]["Occupancy"] = 1.0

    variant = _variant(tmp_path, mutate)
    with pytest.raises(ApcError, match="neither 0 nor 2"):
        apc.rank_orbitals(parse_orca_json(variant))


def test_a_multi_spin_k_block_is_refused(tmp_path):
    def mutate(molecule: dict[str, Any]) -> None:
        molecule["K-Matrix"].append(molecule["K-Matrix"][0])

    variant = _variant(tmp_path, mutate)
    with pytest.raises(ApcError, match="spin matrices"):
        apc.rank_orbitals(parse_orca_json(variant))


# --- output -------------------------------------------------------------------


def test_the_report_prints_the_selection_and_the_boundaries():
    export = parse_orca_json(N2)
    body = apc.render(
        apc.analyze(export, cap=apc.CSF_CAPS["max(10,10)"], cap_label="max(10,10)")
    )
    assert "Ranked-orbital active-space selection (APC)" in body
    assert "cap max(10,10) = 19404 CSFs; selected (10e, 10o)" in body
    assert "active orbitals: 2, 3, 4, 5, 6, 7, 8, 9, 10, 17" in body
    assert "screening device" in body
    assert "R^2 = 0.64" in body


def test_evidence_carries_the_source_and_the_measured_record():
    kinds = {item.kind for item in apc.evidence()}
    assert kinds == {"literature", "measured"}
    bibkeys = {item.bibkey for item in apc.evidence() if item.bibkey}
    assert bibkeys == {"king2021ranked"}
