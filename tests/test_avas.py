"""Tests for the AVAS target projection (analysis/avas.py).

The exact anchors: with the target equal to the whole AO set the projector is the
identity, so every occupied overlap is exactly 1; the eigenvalues stay in
[0, 1]; at most |A| of them are non-zero; and summed over a complete orbital set
they equal |A| (the trace identity, derived from C C^T = S^-1).  The real-data
anchor is the N2 CAS(6,6) export: targeting the two nitrogen 2p shells
reproduces the textbook (6 electrons, 6 orbitals) active space.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from fblockkit.analysis.avas import (
    AvasError,
    analyze,
    evidence,
    projection_spectrum,
    run,
    target_aos,
)
from fblockkit.parsers.orca_json import AoLabel, OrcaJson, parse_orca_json

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"
CANONICAL = FIXTURES / "n2_fcidump.canonical.json"


def _export():
    export = parse_orca_json(CANONICAL)
    coefficients = np.array(export.mo_coefficients).T
    return export, coefficients, np.array(export.overlap)


# --- the exact identities -----------------------------------------------------


def test_the_whole_ao_set_is_the_identity_projector():
    export, coefficients, overlap = _export()
    target = tuple(range(export.n_ao))
    occupied = [i for i, occ in enumerate(export.mo_occupations) if occ > 0.5]
    spectrum = projection_spectrum(coefficients[:, occupied], overlap, target)
    assert len(spectrum) == len(occupied)
    assert max(abs(value - 1.0) for value in spectrum) < 1e-10


def test_the_eigenvalues_are_bounded_and_few():
    export, coefficients, overlap = _export()
    target = target_aos(export.ao_labels, 0, 1, {2})
    spectrum = projection_spectrum(coefficients, overlap, target)
    assert all(0.0 <= value <= 1.0 for value in spectrum)
    nonzero = [value for value in spectrum if value > 1e-10]
    assert len(nonzero) <= len(target)


def test_the_sum_over_a_complete_orbital_set_is_the_target_size():
    # C C^T = S^-1 for a complete orthonormal MO set, so the projected overlaps
    # sum to |A| exactly -- an independent check of the sigma^-1 construction
    export, coefficients, overlap = _export()
    target = target_aos(export.ao_labels, 0, 1, {2}) + target_aos(export.ao_labels, 1, 1, {2})
    spectrum = projection_spectrum(coefficients, overlap, target)
    assert sum(spectrum) == pytest.approx(len(target), abs=1e-10)


# --- the real-data anchor -----------------------------------------------------


def test_the_n2_2p_target_reproduces_the_textbook_cas66():
    export, _, _ = _export()
    result = analyze(export, centre=0, angular=1, shells={2}, threshold=0.1)
    assert result.n_orbitals == 6
    assert result.n_electrons == pytest.approx(6.0, abs=0.2)
    assert result.active_virtual  # the pi* and sigma* combinations are found
    assert len(result.active_occupied) == 3
    # the three occupied p-derived orbitals are the sigma and the pi pair
    assert [round(value, 4) for value in result.occupied_overlaps[:3]] == [
        0.6963,
        0.5595,
        0.5595,
    ]


def test_the_report_names_the_fractional_occupation_convention():
    export, _, _ = _export()
    body = run(analyze(export, centre=0, angular=1, shells={2})).body
    assert "fractional occupation" in body
    assert "recommended active space" in body
    assert "(6" in body  # the (n_el, n_orb) line
    assert "falsify" in body  # the source's own limit
    assert "orbital order" in body  # the ORCA-side boundary


def test_the_target_selection_is_by_centre_shell_and_angular_momentum():
    export = parse_orca_json(CANONICAL)
    assert target_aos(export.ao_labels, 0, 1, {2}) == (6, 7, 8)
    assert target_aos(export.ao_labels, 0, 1, {1}) == (3, 4, 5)
    assert len(target_aos(export.ao_labels, 0, 2)) == 5


# --- errors and the f-block default ------------------------------------------


def test_a_missing_shell_raises_with_a_next_step():
    export = parse_orca_json(CANONICAL)
    with pytest.raises(AvasError, match="Next step"):
        target_aos(export.ao_labels, 0, 3)  # no f shell in def2-SVP


def test_no_f_block_element_raises_and_an_explicit_centre_works():
    export = parse_orca_json(CANONICAL)
    with pytest.raises(AvasError, match="no f-block element"):
        analyze(export)  # N2 has no f-block centre to auto-detect
    result = analyze(export, centre=1, angular=1)  # every p shell of the second N
    assert result.centre == 1 and result.element == "N"
    assert result.shells == (1, 2)


def test_an_unknown_option_is_refused():
    export = parse_orca_json(CANONICAL)
    with pytest.raises(AvasError, match="option 2 or 3|not 2 or 3"):
        analyze(export, centre=0, angular=1, option=1)


# --- open shells ---------------------------------------------------------------


def _open_shell_export() -> OrcaJson:
    """A minimal synthetic export: one doubly, one singly occupied, one empty MO."""
    labels = tuple(
        AoLabel(raw=f"0X  1{name}", center=0, element="X", shell=1, angular="p", component=name)
        for name in ("z", "x", "y")
    )
    overlap = np.eye(3)
    # orbital 0 = pz (doubly), orbital 1 = px (singly), orbital 2 = py (empty)
    coefficients = (
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
    )
    return OrcaJson(
        base_name="synthetic",
        charge=1,
        multiplicity=2,
        hftyp="ROHF",
        point_group="C1",
        atoms=("X",),
        coordinates=((0.0, 0.0, 0.0),),
        n_mo=3,
        n_ao=3,
        mo_coefficients=coefficients,
        mo_occupations=(2.0, 1.0, 0.0),
        mo_energies=(-1.0, -0.5, -0.1),
        overlap=tuple(tuple(float(v) for v in row) for row in overlap),
        ao_labels=labels,
    )


def test_option_three_keeps_every_singly_occupied_orbital():
    export = _open_shell_export()
    # target only the empty py: the singly occupied px is kept anyway by option 3
    result = analyze(
        export, centre=0, angular=1, shells={1}, threshold=0.1, option=3
    )
    assert 1 in result.active_occupied
    assert result.singly_occupied == (1,)
    assert any("option 3" in verdict for verdict in result.verdicts)
    # with option 2 the same orbital is not forced in -- but here it is a target
    # orbital anyway, so the space is identical; the note is what differs
    option_two = analyze(export, centre=0, angular=1, shells={1}, threshold=0.1, option=2)
    assert option_two.active_occupied == result.active_occupied
    assert not any("option 3" in verdict for verdict in option_two.verdicts)


def test_an_orca_block_is_emitted_for_f_targets():
    export = _open_shell_export()
    body = run(analyze(export, centre=0, angular=1, shells={1})).body
    assert "%scf" in body and "m_l" in body and "center 0" in body


def test_evidence_carries_the_source():
    text = " ".join(entry.ref + " " + entry.text for entry in evidence())
    assert "10.1021/acs.jctc.7b00128" in text
    assert "0.05-0.1" in text or "0.05--0.1" in text
