"""Tests for the Molekel MKL reader/writer (parsers/mkl.py).

The fixture is the ORCA-written mkl of the N2/def2-SVP RHF at 1.600 Angstrom
(with the same run's canonical ``orca_2json`` export next to it): the parse is
validated against ORCA's own export -- the reader must reproduce the gbw
convention coefficients to the mkl's print precision (5e-8), which also pins
the measured p-shell row order -- and the writer is validated by a full parse
-> edit -> render -> parse cycle.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from fblockkit.parsers.mkl import MklError, parse_mkl
from fblockkit.parsers.orca_json import parse_orca_json

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"
MKL = FIXTURES / "n2_scan_1.600.mkl"
EXPORT = FIXTURES / "n2_scan_1.600.json"  # the same run's canonical orbitals


def test_the_geometry_and_blocks_are_read():
    mkl = parse_mkl(MKL)
    assert mkl.atoms == ((7, 0.0, 0.0, 0.0), (7, 0.0, 0.0, 1.6))
    assert (mkl.n_ao, mkl.n_mo) == (28, 28)
    assert [len(group.labels) for group in mkl.groups] == [5, 5, 5, 5, 5, 3]
    assert sum(mkl.occupations) == pytest.approx(14.0, abs=1e-9)
    assert mkl.charge_mult == (0, 1)


def test_the_coefficients_reproduce_orcas_own_export():
    """The reader returns gbw-convention rows: identical to ORCA's export.

    This pins two measured facts at once: the group layout (label line,
    energy line, one row per AO) and the p-shell component reordering
    (the mkl stores (py, pz, px) where the export lists (pz, px, py)).
    """
    mkl = parse_mkl(MKL)
    export = parse_orca_json(EXPORT)
    coefficients = np.array(mkl.coefficients())
    reference = np.array(export.mo_coefficients).T
    assert np.abs(coefficients - reference).max() < 1e-7
    overlap = np.array(export.overlap)
    identity = coefficients.T @ overlap @ coefficients
    assert np.abs(identity - np.eye(28)).max() < 1e-6


def test_the_writer_round_trips_an_edit(tmp_path):
    mkl = parse_mkl(MKL)
    export = parse_orca_json(EXPORT)
    replacement = np.array(export.mo_coefficients).T  # same orbitals back in
    updated = mkl.with_orbitals(replacement)
    path = tmp_path / "edited.mkl"
    path.write_text(updated.render(), encoding="utf-8")
    again = parse_mkl(path)
    assert np.abs(np.array(again.coefficients()) - replacement).max() < 1e-7
    assert again.atoms == mkl.atoms
    assert again.occupations == pytest.approx(mkl.occupations)
    # the rendered file keeps the block terminators the converter expects
    text = path.read_text(encoding="utf-8")
    assert text.count("$END") == 6  # one per section, as the measured file
    assert text.count(" $END") == 2  # the coefficient and occupation terminators
    assert text.count("$COEFF_ALPHA") == 1 and text.count("$OCC_ALPHA") == 1


def test_a_wrong_shaped_replacement_is_refused():
    mkl = parse_mkl(MKL)
    with pytest.raises(MklError, match="template holds"):
        mkl.with_orbitals([[0.0] * 5])


def test_a_file_without_the_sections_is_refused(tmp_path):
    path = tmp_path / "not.mkl"
    path.write_text("hello world\n", encoding="utf-8")
    with pytest.raises(MklError, match="measured MKL sections"):
        parse_mkl(path)


# --- the unrestricted (beta) extension: the UKS CH4 reference of menu 29 ---------

UKS_MKL = FIXTURES / "ch4_diss_prop.mkl"  # the menu-29 reference (UKS PBE0/def2-SVP)


def test_the_beta_block_of_an_unrestricted_file_is_read():
    mkl = parse_mkl(UKS_MKL)
    assert mkl.unrestricted
    assert (mkl.n_ao, mkl.n_mo) == (34, 34)
    assert len(mkl.beta_coefficients()) == 34
    assert len(mkl.beta_occupations) == mkl.n_mo
    assert sum(mkl.occupations) == pytest.approx(5.0, abs=1e-9)
    assert sum(mkl.beta_occupations) == pytest.approx(5.0, abs=1e-9)


def test_the_beta_blocks_survive_a_render_parse_cycle(tmp_path):
    mkl = parse_mkl(UKS_MKL)
    path = tmp_path / "uks.mkl"
    path.write_text(mkl.render(), encoding="utf-8")
    text = path.read_text(encoding="utf-8")
    assert text.count("$COEFF_ALPHA") == 1 and text.count("$COEFF_BETA") == 1
    assert text.count("$OCC_ALPHA") == 1 and text.count("$OCC_BETA") == 1
    again = parse_mkl(path)
    assert again.render() == text  # stable
    assert np.abs(
        np.array(again.beta_coefficients()) - np.array(mkl.beta_coefficients())
    ).max() < 1e-7
    assert again.beta_occupations == pytest.approx(mkl.beta_occupations)


def test_a_restricted_file_refuses_a_beta_replacement():
    mkl = parse_mkl(MKL)
    assert not mkl.unrestricted
    assert mkl.beta_coefficients() == ()
    with pytest.raises(MklError, match="restricted file"):
        mkl.with_orbitals(mkl.coefficients(), beta_coefficients=mkl.coefficients())


def test_a_beta_replacement_reaches_the_written_file(tmp_path):
    mkl = parse_mkl(UKS_MKL)
    alpha = [list(row) for row in mkl.coefficients()]
    beta = [list(row) for row in mkl.beta_coefficients()]
    # a transposition of two alpha columns and a sign flip on one beta column
    for row in alpha:
        row[1], row[2] = row[2], row[1]
    for row in beta:
        row[0] = -row[0]
    updated = mkl.with_orbitals(alpha, beta_coefficients=beta)
    path = tmp_path / "edited_uks.mkl"
    path.write_text(updated.render(), encoding="utf-8")
    again = parse_mkl(path)
    assert np.abs(np.array(again.coefficients()) - np.array(alpha)).max() < 1e-7
    assert np.abs(np.array(again.beta_coefficients()) - np.array(beta)).max() < 1e-7
