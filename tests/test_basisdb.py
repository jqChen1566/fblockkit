"""Checks of the basis-library query layer (knowledge/basisdb.py; plan item 6.2).

The fixture is a four-row miniature baked from the real deployed library
(101, ~/projects/orca_basis_sets; BSE v0.12 snapshot): def2-SVP and
def2-TZVP (both cover Ce), 3-21G (does not) -- see fixtures/basisdb/README.
Anchors from the real rows: def2-SVP ahlrichs / rel=none / built-in / ECP /
1058 contracted functions / 30 refs.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from pytest import approx

from fblockkit.knowledge import basisdb

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
MINI = FIXTURES / "basisdb" / "mini_basis.db"


def test_the_element_query_finds_the_covering_sets():
    data = basisdb.query_element(MINI, "Ce")
    names = [entry["name"] for entry in data["orbital"]]
    assert names == ["def2-svp", "def2-tzvp"]
    first = data["orbital"][0]
    assert first["family"] == "ahlrichs"
    assert first["relativistic"] == "none"
    assert first["orca_builtin"] is True
    assert first["has_ecp"] is True
    assert first["size"] == 1058
    assert first["ref_count"] == 30
    assert "Ce" in first["elements"]
    assert data["auxiliary_counts"] == {}


def test_the_render_reports_the_library_leg(tmp_path):
    body = basisdb.render_query(MINI, "Ce")
    assert "Basis-library query: Ce at" in body
    assert "2 entries: 2 orbital + none auxiliary" in body
    assert "[built-in] def2-svp" in body
    assert "1058 fn  refs 30  ECP" in body
    assert "GTOName" in body


def test_the_render_truncates_long_lists():
    body = basisdb.render_query(MINI, "Ce", max_rows=1)
    assert "... (1 more; the full listing is INDEX.md in the library)" in body


def test_the_refusals():
    with pytest.raises(basisdb.BasisdbError, match="no entry covering U"):
        basisdb.render_query(MINI, "U")
    with pytest.raises(basisdb.BasisdbError, match="not an element symbol"):
        basisdb.query_element(MINI, "ce")
    with pytest.raises(basisdb.BasisdbError, match="not an element symbol"):
        basisdb.query_element(MINI, "Carbon")


def test_the_library_lookup_order(monkeypatch, tmp_path):
    import shutil

    # isolate from any real deployment on this host (101 carries one)
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    monkeypatch.delenv("FBK_BASISDB", raising=False)
    # found via an explicit file path (the file name needs not be basis.db)
    assert basisdb.find_library(MINI) == MINI
    assert basisdb.find_library(tmp_path / "missing") is None
    # $FBK_BASISDB may be the index file itself ...
    monkeypatch.setenv("FBK_BASISDB", str(MINI))
    assert basisdb.find_library() == MINI
    # ... or a directory holding a basis.db
    library_dir = tmp_path / "lib"
    library_dir.mkdir()
    shutil.copy(MINI, library_dir / "basis.db")
    monkeypatch.setenv("FBK_BASISDB", str(library_dir))
    assert basisdb.find_library() == library_dir / "basis.db"
    monkeypatch.delenv("FBK_BASISDB")
    assert basisdb.find_library() is None
