"""Tests for the cross-structure orbital mapping (analysis/orbital_mapping.py).

The real-data anchors are the frozen N2 scan (1.094 / 1.600 / 2.600 Angstrom,
RHF/def2-SVP, occupied block localized by IAO-BOYS, virtual block by
Foster-Boys, exports carrying the T-Matrix): the two 1s cores map uniquely
across the three structures, the three equivalent bond orbitals form one
degenerate class, the two one-centre hybrids are the non-matchable set at the
source's tau = 0.5, and the tau curve shows the source's plateau behaviour.
The map machinery itself (the set-valued m_LK condition, the union rule) is
pinned with small synthetic descriptor tables.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from fblockkit.analysis import orbital_mapping as om
from fblockkit.parsers.orca_json import parse_orca_json

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"
SCAN = ("1.094", "1.600", "2.600")


def _structures():
    return [
        om.descriptors_from_export(parse_orca_json(FIXTURES / f"n2_scan_{r}.loc.json"), f"r={r}")
        for r in SCAN
    ]


def _orbitals(*rows):
    """Small hand-built descriptors: (index, occupation, kinetic, {shell: population})."""
    descriptors = []
    for index, occupation, kinetic, shells in rows:
        keys = sorted(shells)
        descriptors.append(
            om.OrbitalDescriptor(
                index=index,
                occupation=occupation,
                kinetic=kinetic,
                shells=tuple((om.ShellRef(*key), shells[key]) for key in keys),
            )
        )
    return descriptors


def _structure(name, *rows):
    return om.StructureDescription(name=name, orbitals=tuple(_orbitals(*rows)))


# --- the real scan ------------------------------------------------------------


def test_the_descriptors_read_the_t_matrix():
    structures = _structures()
    reference = structures[0]
    assert len(reference.orbitals) == 28
    cores = [orbital for orbital in reference.orbitals if orbital.occupation > 1.9 and orbital.kinetic > 20.0]
    assert len(cores) == 2
    assert cores[0].kinetic == pytest.approx(22.6155, abs=1e-3)
    # the core orbitals are 1s on one centre each
    for orbital in cores:
        dominant = orbital.dominant(1)[0]
        assert dominant[0].shell == 1 and dominant[0].angular == "s"
        assert dominant[1] == pytest.approx(1.0, abs=0.01)


def test_the_cores_map_uniquely_and_the_hybrids_do_not_match():
    result = om.analyze(_structures())
    occupied = result.classes["occupied"]
    # two single-member classes = the two cores; one 3x3 class = the bond triad
    singles = [klass for klass in occupied if len(klass) == 3]
    assert len(singles) == 2
    for klass in singles:
        assert sorted(orbital for _, orbital in klass) in ([0, 0, 0], [2, 2, 2])
    triads = [klass for klass in occupied if len(klass) == 9]
    assert len(triads) == 1
    for structure, orbital in triads[0]:
        assert orbital in (4, 5, 6)
    # the non-matchable set: the two one-centre hybrids in every structure
    assert result.no_match["occupied"] == ((1, 3), (1, 3), (1, 3))


def test_the_tau_curve_shows_the_source_plateaus():
    result = om.analyze(_structures())
    curve = {tau: (occ, virt) for tau, occ, virt in result.tau_curve}
    assert curve[0.05][0] == 21  # nothing but nothing: all occupied unmatched
    assert curve[0.2][0] == 15  # the whole valence unmatched
    assert curve[0.5][0] == 6  # the hybrids only
    assert result.tau_occupied == 0.5 and result.tau_virtual == 0.5
    assert not result.tau_given


def test_the_classes_form_a_partition():
    result = om.analyze(_structures())
    for name, expected in (("occupied", 21), ("virtual", 63)):
        members = [pair for klass in result.classes[name] for pair in klass]
        assert len(members) == expected
        assert len(set(members)) == len(members)


def test_a_selection_becomes_a_consistent_active_space():
    structures = _structures()
    # selecting one bond orbital of one structure pulls the whole degenerate
    # triad into every structure's consistent space -- and nothing else
    result = om.analyze(structures, selections={"r=1.094": [4]})
    assert result.consistent[0] == (4, 5, 6)
    assert result.consistent[1] == (4, 5, 6)
    assert result.consistent[2] == (4, 5, 6)
    # selecting a non-matchable orbital pulls the non-matchable block
    result = om.analyze(structures, selections={"r=2.600": [3]})
    assert result.consistent[0] == (1, 3)
    assert result.consistent[2] == (1, 3)


def test_an_explicit_tau_is_used_as_given():
    result = om.analyze(_structures(), tau=0.3)
    assert result.tau_given and result.tau_occupied == 0.3 and result.tau_virtual == 0.3
    assert len(result.tau_curve) == 1
    assert result.no_match["occupied"] == ((1, 3, 4, 5, 6),) * 3


# --- the map machinery on synthetic tables ------------------------------------


SHELL_A = (0, "N", 1, "s")
SHELL_B = (0, "N", 2, "s")


def test_the_m_lk_condition_requires_the_whole_set_in_both_directions():
    # structure L: orbitals {0,1} match each other (a degenerate pair) and match
    # structure K's {0,1}; orbital 2 matches nothing
    a = _structure(
        "L",
        (0, 2.0, 1.0, {SHELL_A: 1.0}),
        (1, 2.0, 1.01, {SHELL_A: 1.0}),
        (2, 2.0, 9.0, {SHELL_B: 1.0}),
    )
    b = _structure(
        "K",
        (0, 2.0, 1.02, {SHELL_A: 1.0}),
        (1, 2.0, 1.03, {SHELL_A: 1.0}),
        (2, 2.0, 9.5, {SHELL_B: 1.0}),
    )
    result = om.analyze([a, b], tau=0.1)
    assert result.classes["occupied"] == (((0, 0), (0, 1), (1, 0), (1, 1)), ((0, 2), (1, 2)))
    # structure 2's orbital 2 is non-matchable (9.0 vs 9.5: 0.5 > tau)
    assert result.no_match["occupied"] == ((2,), (2,))


def test_a_one_sided_match_falls_into_the_non_matchable_block():
    # L's {0,1} form a self-set, but its members map onto different K sets
    # (L0 -> {K0}, L1 -> {K0,K1}): the set-to-set condition fails for some
    # structure pair, so the whole set is non-matchable -- the source's "for
    # some structure combination" clause
    a = _structure(
        "L",
        (0, 2.0, 1.00, {SHELL_A: 1.0}),
        (1, 2.0, 1.05, {SHELL_A: 1.0}),
    )
    b = _structure(
        "K",
        (0, 2.0, 1.06, {SHELL_A: 1.0}),
        (1, 2.0, 1.11, {SHELL_A: 1.0}),
    )
    result = om.analyze([a, b], tau=0.1)
    assert result.no_match["occupied"] == ((0, 1), (0, 1))


# --- refusals -----------------------------------------------------------------


def test_a_structure_without_the_kinetic_matrix_is_refused():
    export = parse_orca_json(FIXTURES / "n2_fcidump.localized.json")
    assert export.kinetic is None
    with pytest.raises(om.MappingError, match="T-Matrix"):
        om.descriptors_from_export(export, "no-T")


def test_fractional_occupations_are_refused():
    export = parse_orca_json(FIXTURES / "n2_scan_1.094.loc.json")
    from dataclasses import replace

    fractional = replace(export, mo_occupations=(1.0,) + export.mo_occupations[1:])
    with pytest.raises(om.MappingError, match="fractional"):
        om.descriptors_from_export(fractional, "cas")


def test_different_shell_sets_are_refused():
    a = _structure("L", (0, 2.0, 1.0, {SHELL_A: 1.0}))
    b = _structure("K", (0, 2.0, 1.0, {(1, "O", 1, "s"): 1.0}))
    with pytest.raises(om.MappingError, match="shell set"):
        om.analyze([a, b], tau=0.1)


def test_a_single_structure_is_refused():
    a = _structure("L", (0, 2.0, 1.0, {SHELL_A: 1.0}))
    with pytest.raises(om.MappingError, match="at least two"):
        om.analyze([a], tau=0.1)


def test_different_orbital_counts_are_refused():
    a = _structure("L", (0, 2.0, 1.0, {SHELL_A: 1.0}), (1, 2.0, 1.0, {SHELL_A: 1.0}))
    b = _structure("K", (0, 2.0, 1.0, {SHELL_A: 1.0}))
    with pytest.raises(om.MappingError, match="orbital"):
        om.analyze([a, b], tau=0.1)


# --- the report and evidence --------------------------------------------------


def test_the_report_states_the_maps_the_curve_and_the_boundaries():
    structures = _structures()
    body = om.render(om.analyze(structures, selections={"r=1.094": [4]}), structures)
    assert "Cross-structure orbital mapping" in body
    assert "Tau curve" in body
    assert "non-matchable" in body
    assert "Consistent active space" in body
    assert "Loewdin populations" in body  # the documented population substitution
    assert "IAO-BOYS" in body  # the measured localization boundary


def test_evidence_carries_the_source_and_the_route():
    text = " ".join(
        entry.ref + " " + entry.text + " " + entry.url + " " + entry.bibkey
        for entry in om.evidence()
    )
    assert "10.1021/acs.jpclett.2c03905" in text
    assert "bensberg2023corresponding" in text
    assert "5e-7" in text or "5e-07" in text


# --- the active-space overlap check (AOP) -------------------------------------


def _scan(tag: str):
    from fblockkit.parsers.orca_json import parse_orca_json

    return parse_orca_json(FIXTURES / f"n2_scan_{tag}.loc.json")


def test_the_active_overlap_reads_small_and_large_steps():
    """The AOP scalar on the real scan: cores survive everything, the bond
    triad survives the 0.01-Angstrom step and degrades across 0.5 Angstrom."""
    triad = [4, 5, 6]
    assert om.active_overlap_determinant(_scan("1.600"), _scan("1.610"), triad) == pytest.approx(
        0.9955, abs=1e-3
    )
    assert om.active_overlap_determinant(_scan("1.094"), _scan("1.600"), triad) == pytest.approx(
        0.7172, abs=1e-3
    )
    assert om.active_overlap_determinant(_scan("1.094"), _scan("2.600"), [0, 2]) == pytest.approx(
        1.0007, abs=1e-3
    )


def test_an_active_inactive_exchange_drives_the_determinant_to_zero():
    """Swapping one active orbital with an inactive one is the source's failure
    mode -- while a rotation inside the degenerate window leaves it unchanged."""
    from dataclasses import replace

    a = _scan("1.600")
    b = _scan("1.610")
    rows = [list(row) for row in b.mo_coefficients]
    rows[4], rows[20] = rows[20], rows[4]
    exchanged = replace(b, mo_coefficients=tuple(tuple(row) for row in rows))
    assert om.active_overlap_determinant(a, exchanged, [4, 5, 6]) < 1e-3
    rows = [list(row) for row in b.mo_coefficients]
    rows[4], rows[5] = rows[5], rows[4]
    rotated = replace(b, mo_coefficients=tuple(tuple(row) for row in rows))
    assert om.active_overlap_determinant(a, rotated, [4, 5, 6]) == pytest.approx(0.9955, abs=1e-3)


def test_the_series_runs_over_adjacent_pairs():
    exports = [_scan(tag) for tag in ("1.600", "1.610", "2.600")]
    names = ["r=1.600", "r=1.610", "r=2.600"]
    actives = {name: [4, 5, 6] for name in names}
    rows = om.active_overlap_series(exports, actives, names)
    assert [row[:2] for row in rows] == [("r=1.600", "r=1.610"), ("r=1.610", "r=2.600")]
    assert rows[0][2] > 0.99 and rows[1][2] < 0.9


def test_the_aop_refusals():
    exports = [_scan("1.600"), _scan("1.610")]
    with pytest.raises(om.MappingError, match="no active-space list"):
        om.active_overlap_series(exports, {}, ["a", "b"])
    wrong = {"a": [4, 5, 6], "b": [4, 5]}
    with pytest.raises(om.MappingError, match="equally sized"):
        om.active_overlap_series(exports, wrong, ["a", "b"])
    with pytest.raises(om.MappingError, match="outside"):
        om.active_overlap_determinant(exports[0], exports[1], [99])


def test_the_report_carries_the_aop_block_when_given():
    structures = _structures()
    rows = om.active_overlap_series(
        [_scan(tag) for tag in SCAN],
        {f"r={tag}": [4, 5, 6] for tag in SCAN},
        [f"r={tag}" for tag in SCAN],
    )
    section = om.run(structures, active_overlap=rows)
    assert "Active-space overlap" in section.body
    assert "preserved" in section.body
