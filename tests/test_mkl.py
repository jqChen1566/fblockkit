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
