"""Analysis-layer tests: A1 orbital composition, A2 entropy spectrum, S1 coordination
geometry.

Fixtures: fixtures/orca/*.out (including two samples generated specifically for A1, with
per-MO composition tables).
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from fblockkit.analysis import composition, entropy, geometry
from fblockkit.parsers import parse_auto

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"


def _result(name: str):
    return parse_auto(FIXTURES / name)


# --- A1 orbital composition --------------------------------------------------


def test_composition_accepts_and_classifies():
    result = _result("n2_casscf_orbcomp.out")
    assert composition.accepts(result)
    rows = composition.orbital_rows(result)
    assert len(rows) == 28
    # MO 0 is the sigma bonding orbital of N=N: about half s on each N
    first = rows[0]
    assert first.occupation == pytest.approx(2.0)
    assert first.dominant["shell"] == "s"
    assert first.dominant["weight"] == pytest.approx(49.9, abs=0.2)
    # the shell weights of one orbital add up to about 100% (Loewdin percentages)
    assert sum(item["weight"] for item in first.shells) == pytest.approx(100.0, abs=1.0)


def test_composition_ranking_partitions_by_occupation():
    """Regression (this group's EuF lesson): a ranking must be partitioned by occupation
    first; virtual orbitals must not leak into the occupied ranking."""
    rows = composition.orbital_rows(_result("ce3_orbcomp.out"))
    occupied = composition.shell_ranking(rows, "f", partition="occupied", top_n=5)
    assert occupied
    assert all(row.occupation > 0.02 for row in occupied)
    # the top three f weights in the whole space are all virtual (occ = 0) -- they must be
    # excluded once partitioned
    top_overall = composition.shell_ranking(rows, "f", top_n=3)
    assert all(row.occupation <= 0.02 for row in top_overall)
    assert all(row.shell_weight("f") > 99.0 for row in top_overall)
    occupied_indices = {row.index for row in occupied}
    assert not occupied_indices & {row.index for row in top_overall}


def test_composition_section_reports_active_orbital_shell():
    """A measured fact of the Ce3+ fixture: the singly occupied active orbital is 100% Ce-d
    (every f orbital is virtual).

    That is, this calculation landed on the d1 solution (4f not occupied) -- A1 correctly
    reports its dominant shell as d; for an f1 target the initial guess must be changed or
    multiple starting points used (another instance of the multiple-solution problem; see
    the fixture README notes).
    """
    section = composition.run(_result("ce3_orbcomp.out"), shell="f", top_n=3)
    assert section.title.startswith("A1")
    assert "Dominant shells of the active orbitals" in section.body
    assert "MO   27  occ 1.0000  Ce1 d 100.0%" in section.body
    assert composition.evidence()[0].kind == "measured"


def test_composition_rejects_output_without_table():
    with pytest.raises(ValueError, match="Next step: "):
        composition.run(_result("fblock_dft_la_complex.out"))


# --- A2 entropy spectrum -----------------------------------------------------


def test_entropy_bound_formula_and_limits():
    assert entropy.single_orbital_entropy_bound(0.0) == 0.0
    assert entropy.single_orbital_entropy_bound(2.0) == 0.0
    # the open-shell single-occupation limit of the bound is ln 4
    assert entropy.single_orbital_entropy_bound(1.0) == pytest.approx(math.log(4.0))
    with pytest.raises(ValueError, match="out of range"):
        entropy.single_orbital_entropy_bound(2.5)


def test_entropy_bound_is_an_upper_bound_over_all_completions():
    """The defining property (module docstring): for every admissible completion of the
    four occupation weights with spin-summed occupation n -- alpha/beta occupations
    n_a + n_b = n, double occupancy p2 in [max(0, n-1), min(n_a, n_b)] -- the true
    four-state entropy -sum w ln w is at most s_bound(n); the spin-independent completion
    n_a = n_b = n/2, p2 = n_a * n_b attains the bound exactly."""
    import numpy as np

    for n in (0.05, 0.3, 0.7, 1.0, 1.4, 1.9):
        bound = entropy.single_orbital_entropy_bound(n)
        # the uncorrelated completion attains the bound exactly
        n_a = n_b = n / 2.0
        p2 = n_a * n_b
        w = np.array([1 - 2 * n_a + p2, n_a - p2, n_b - p2, p2])
        s_uncorrelated = float(-(w[w > 0] * np.log(w[w > 0])).sum())
        assert s_uncorrelated == pytest.approx(bound, rel=1e-12)
        # every admissible completion stays at or below it
        for n_alpha in np.linspace(0.0, min(1.0, n), 41):
            n_beta = n - n_alpha
            if n_beta < 0.0 or n_beta > 1.0:
                continue
            for p2 in np.linspace(max(0.0, n - 1.0), min(n_alpha, n_beta), 41):
                weights = np.array(
                    [1 - n_alpha - n_beta + p2, n_alpha - p2, n_beta - p2, p2]
                )
                positive = weights[weights > 0]
                s_true = float(-(positive * np.log(positive)).sum())
                assert s_true <= bound + 1e-12


def test_entropy_bound_excludes_for_a_closed_shell_single_determinant():
    """The one-sided reading: an all-integer (0/2) occupation spectrum gives bound 0, so
    the exclusion direction fires on a synthetic closed-shell case."""
    from fblockkit.knowledge.models import ParseResult

    result = ParseResult(
        program="orca",
        path="synthetic",
        sections={"casscf": {"active_occupations": (2.0, 2.0, 2.0, 0.0, 0.0)}},
    )
    section = entropy.run(result)
    assert "max(s_bound) = 0.0000" in section.body
    assert "Exclusion test" in section.body and "no orbital shows significant" in section.body


def test_entropy_section_on_f1_atom():
    section = entropy.run(_result("generated_ce3_sarc2.out"))
    assert "max(s_bound) = 1.3863" in section.body  # f1 (occupation 1) -> ln 4
    assert "NOT decidable" in section.body
    assert "localized-orbital basis" in section.body


def test_entropy_section_on_n2_casscf():
    section = entropy.run(_result("n2_casscf_nevpt2.out"))
    assert section.title.startswith("A2")
    assert "Active space: 6 orbitals" in section.body
    # the orbital with occupation 0.00208 has a small but non-zero entropy
    assert "0.0" in section.body


def test_entropy_rejects_non_casscf_output():
    with pytest.raises(ValueError, match="Next step: "):
        entropy.run(_result("fblock_dft_la_complex.out"))


# --- S1 coordination geometry ------------------------------------------------


_OCTAHEDRON = """7
comment
Ce  0.0 0.0 0.0
O   2.4 0.0 0.0
O  -2.4 0.0 0.0
O   0.0 2.4 0.0
O   0.0 -2.4 0.0
O   0.0 0.0 2.4
O   0.0 0.0 -2.4
"""

_LINEAR = """2
comment
Ce 0.0 0.0 0.0
O  0.0 0.0 2.3
"""


def test_geometry_parses_and_finds_center():
    atoms = geometry.parse_xyz(_OCTAHEDRON)
    assert len(atoms) == 7
    assert atoms[0].element == "Ce"
    assert geometry.dominant_center(atoms) == 0
    shell, _ = geometry.coordination_shell(atoms, 0)
    assert len(shell) == 6


def test_geometry_octahedron():
    section = geometry.run(_OCTAHEDRON)
    assert "Coordination shell: 6 ligands" in section.body
    assert "nearly isotropic" in section.body
    # three inversion pairs of ligands
    assert "Inversion pairs: 3/6" in section.body
    assert "C3 → 9" in section.body and "Oh → 4" in section.body


def test_geometry_square_plane_is_planar():
    xyz = """5
comment
Ce 0.0 0.0 0.0
O  2.0 0.0 0.0
O -2.0 0.0 0.0
O  0.0 2.0 0.0
O  0.0 -2.0 0.0
"""
    section = geometry.run(xyz)
    assert "nearly planar distribution" in section.body


def test_geometry_linear_two_coordinate():
    section = geometry.run(_LINEAR)
    assert "nearly linear distribution" in section.body


def test_geometry_rejects_bad_input(tmp_path):
    bad = tmp_path / "bad.xyz"
    bad.write_text("3\ncomment\nCe 0 0 0\n", encoding="utf-8")
    with pytest.raises(geometry.StructureError, match="Next step: "):
        geometry.parse_xyz(bad)


def test_geometry_cf_count_for_known_point_group():
    section = geometry.run(_OCTAHEDRON, cf_point_group="Oh")
    assert "Point group confirmed by the user, Oh: 4 B_k^q allowed" in section.body
