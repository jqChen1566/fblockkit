"""FCIDUMP reader tests.

Fixture: fixtures/orca/n2_fcidump.fcidump -- the active-space Hamiltonian of the
N2 CASSCF(6,6) reference (ORCA 6.1.1, the two-step !FCIDUMP recipe); see
fixtures/orca/README.md for provenance.  The numeric values double as the
reference for the CI tests in tests/test_entropy_rdm.py.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from fblockkit.parsers.base import ParserError
from fblockkit.parsers.fcidump import parse_fcidump

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"
DUMP = FIXTURES / "n2_fcidump.fcidump"


@pytest.fixture(scope="module")
def dump():
    return parse_fcidump(DUMP)


def test_header(dump):
    assert dump.norb == 6
    assert dump.nelec == 6
    assert dump.ms2 == 0
    assert dump.ecore == pytest.approx(-97.5043826450, abs=1e-9)
    assert dump.orbsym == (1,) * 6
    assert dump.isym == 1


def test_one_electron_block_is_symmetric(dump):
    for p in range(dump.norb):
        for q in range(dump.norb):
            assert dump.h[p][q] == dump.h[q][p]
    # the Fermi-contact-like diagonal of orbital 6 (the weakly occupied sigma*)
    assert dump.h[5][5] == pytest.approx(-1.52850478652576526883e00, abs=1e-12)


def test_two_electron_values_and_symmetry(dump):
    # three body lines of the fixture, with their 8-fold permutation families
    assert dump.g[(0, 0, 0, 0)] == pytest.approx(5.26622545380221795952e-01, abs=1e-15)
    value_1010 = 3.96641108080904650213e-02  # line "2 1 2 1"
    for perm in [(1, 0, 1, 0), (0, 1, 1, 0), (1, 0, 0, 1), (0, 1, 0, 1)]:
        assert dump.g[perm] == pytest.approx(value_1010, abs=1e-15)
    value_1100 = 4.80049958301277002182e-01  # line "2 2 1 1"
    assert dump.g[(1, 1, 0, 0)] == pytest.approx(value_1100, abs=1e-15)
    assert dump.g[(0, 0, 1, 1)] == pytest.approx(value_1100, abs=1e-15)
    # every stored key is a valid 4-tuple of 0-based indices
    assert all(0 <= i < dump.norb for key in dump.g for i in key)


def _write(tmp_path, text: str) -> Path:
    path = tmp_path / "FCIDUMP"
    path.write_text(text, encoding="utf-8")
    return path


def test_missing_file(tmp_path):
    with pytest.raises(ParserError, match="Next step"):
        parse_fcidump(tmp_path / "nope.fcidump")


def test_not_a_fcidump(tmp_path):
    with pytest.raises(ParserError, match="Next step"):
        parse_fcidump(FIXTURES / "n2_fcidump_step_a.out")


def test_truncated_dump_missing_core_line(tmp_path):
    text = (DUMP.read_text(encoding="utf-8").splitlines())[:-1]
    path = _write(tmp_path, "\n".join(text))
    with pytest.raises(ParserError, match="core-energy"):
        parse_fcidump(path)


def test_malformed_body_line(tmp_path):
    path = _write(
        tmp_path,
        "&FCI NORB=2,NELEC=2,MS2=0\n/\n  1.0 1 1 1\n  0.0 0 0 0 0\n",
    )
    with pytest.raises(ParserError, match="5 fields"):
        parse_fcidump(path)


def test_index_out_of_range(tmp_path):
    path = _write(
        tmp_path,
        "&FCI NORB=2,NELEC=2,MS2=0\n/\n  1.0 1 3 1 1\n  0.0 0 0 0 0\n",
    )
    with pytest.raises(ParserError, match="out of range"):
        parse_fcidump(path)


def test_header_without_header_line(tmp_path):
    path = _write(tmp_path, "not an fcidump at all\n/\n0.0 0 0 0 0\n")
    with pytest.raises(ParserError, match="&FCI"):
        parse_fcidump(path)
