"""Regression checks of the AEGISS selection (analysis/aegiss.py).

The fixtures are the benzene/cc-pVDZ pi platform of the source's own example:
``benzene.json`` (the export of the CASSCF gbw) and ``benzene.fcidump`` (the
same run's dump; the run was prepared with an AVAS-prepped pi window, and its
printed CASSCF energy -230.793818898 agrees with the source's -230.793770 to
5e-5).  The anchors:

- the exact window FCI energy reproduces the engine's print to eleven digits;
- the exact four-state entropies [0.174, 0.342] all clear both the source's
  10 % and its benzene 20 % screen lines;
- the projection on ``C pz`` (the whole pz family, 12 target functions) gives
  weights [3.087, 1.867, 1.867, 0.992, 0.992, 0.810] -- all above the 0.5
  threshold -- recovering the textbook (6e, 6o) pi space;
- the measured weight-definition deviation: the source's signed row sum
  cancels for nodal pi orbitals (pinned here via the a2u-only survival on the
  shell-2 label), which is why the projection norm is used;
- a sigma orbital measures exactly 0.0000 on the same target (negative
  control).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable

import pytest

from fblockkit.analysis import aegiss
from fblockkit.analysis.aegiss import AegissError, parse_ao_group
from fblockkit.parsers.fcidump import parse_fcidump
from fblockkit.parsers.orca_json import parse_orca_json

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"
EXPORT = FIXTURES / "benzene.json"
DUMP = FIXTURES / "benzene.fcidump"

ENGINE_ENERGY = -230.793818898


def _export():
    return parse_orca_json(EXPORT)


def _dump():
    return parse_fcidump(DUMP)


def _variant(tmp_path: Path, mutate: Callable[[dict[str, Any]], None]) -> Path:
    document = json.loads(EXPORT.read_text(encoding="utf-8"))
    mutate(document["Molecule"])
    target = tmp_path / "variant.json"
    target.write_text(json.dumps(document), encoding="utf-8")
    return target


# --- the AO label -------------------------------------------------------------


def test_parse_ao_group_reads_the_label_grammar():
    group = parse_ao_group("C 2pz")
    assert (group.element, group.shell, group.angular, group.component) == ("C", 2, "p", "z")
    group = parse_ao_group("Fe d")
    assert (group.element, group.shell, group.angular, group.component) == ("Fe", None, "d", None)
    group = parse_ao_group("N p")
    assert group.component is None and group.shell is None
    with pytest.raises(AegissError, match="does not split"):
        parse_ao_group("C")
    with pytest.raises(AegissError, match="does not parse"):
        parse_ao_group("C 2q3")


def test_matching_aos_resolves_shell_and_component():
    labels = _export().ao_labels
    assert len(aegiss.matching_aos(labels, parse_ao_group("C pz"))) == 12
    assert len(aegiss.matching_aos(labels, parse_ao_group("C 2pz"))) == 6
    assert len(aegiss.matching_aos(labels, parse_ao_group("C 1pz"))) == 6
    assert len(aegiss.matching_aos(labels, parse_ao_group("C p"))) == 36  # 6 C x 2 shells x 3
    assert len(aegiss.matching_aos(labels, parse_ao_group("H s"))) == 12  # 6 H x 2 s shells
    assert aegiss.matching_aos(labels, parse_ao_group("Xx d")) == ()


# --- the two screens ----------------------------------------------------------


def test_entropy_screen_lines_and_refusals():
    entropies = (0.18129, 0.33934, 0.33934, 0.34171, 0.34171, 0.17426)
    kept, line = aegiss.entropy_screen(entropies, 0.1)
    assert kept == (0, 1, 2, 3, 4, 5)  # every pi orbital clears the 10 % line
    assert line == pytest.approx(0.1 * 0.34171, abs=1e-9)
    kept, line = aegiss.entropy_screen(entropies, 0.6)
    assert kept == (1, 2, 3, 4)  # a2u and b2g fall below 0.6 * S_max
    with pytest.raises(AegissError, match="not in"):
        aegiss.entropy_screen(entropies, 0.0)
    with pytest.raises(AegissError, match="zero"):
        aegiss.entropy_screen((0.0, 0.0), 0.1)


def test_projection_weights_measure_pi_and_reject_sigma():
    export = _export()
    window = (18, 19, 20, 21, 22, 23)
    aos = aegiss.matching_aos(export.ao_labels, parse_ao_group("C pz"))
    weights = aegiss.projection_weights(export, window, aos)
    assert weights == pytest.approx(
        (3.0874, 1.8672, 1.8672, 0.9920, 0.9920, 0.8096), abs=1e-3
    )
    # negative control: a sigma orbital carries exactly zero pz content
    sigma = aegiss.projection_weights(export, (17,), aos)
    assert abs(sigma[0]) < 1e-9
    # the shell-resolved label resolves the pi manifold into its subsets
    aos2 = aegiss.matching_aos(export.ao_labels, parse_ao_group("C 2pz"))
    weights2 = aegiss.projection_weights(export, window, aos2)
    assert weights2 == pytest.approx(
        (2.0872, 0.9491, 0.9491, 0.1873, 0.1873, 0.0473), abs=1e-3
    )


# --- the workflow -------------------------------------------------------------


def test_the_benzene_selection_recovers_the_pi_sextet():
    result = aegiss.analyze(_export(), _dump(), label="C pz", reference_energy=ENGINE_ENERGY)
    assert result.window == (18, 19, 20, 21, 22, 23)
    assert result.kept == (0, 1, 2, 3, 4, 5)
    assert result.selected == (0, 1, 2, 3, 4, 5)
    assert (result.n_electrons, result.n_orbitals) == (6, 6)
    assert result.n_group_aos == 12
    assert result.energy_fci == pytest.approx(ENGINE_ENERGY, abs=1e-9)
    assert result.engine_energy == ENGINE_ENERGY
    assert result.entropies == pytest.approx(
        (0.18129, 0.33934, 0.33934, 0.34171, 0.34171, 0.17426), abs=1e-4
    )


def test_the_shell_resolved_label_narrows_the_space():
    """Measured resolution effect: with the shell index given, the projection
    norm keeps only the strongly overlapping subset at epsilon = 0.5."""
    result = aegiss.analyze(_export(), _dump(), label="C 2pz", reference_energy=ENGINE_ENERGY)
    assert result.n_group_aos == 6
    assert result.selected == (0, 1, 2)
    assert (result.n_electrons, result.n_orbitals) == (6, 3)


def test_the_twenty_percent_line_keeps_the_full_set():
    result = aegiss.analyze(
        _export(), _dump(), label="C pz", tau=0.2, reference_energy=ENGINE_ENERGY
    )
    assert result.kept == (0, 1, 2, 3, 4, 5)
    assert (result.n_electrons, result.n_orbitals) == (6, 6)


# --- refusals -----------------------------------------------------------------


def test_a_missing_labels_block_is_refused(tmp_path):
    variant = _variant(tmp_path, lambda m: m["MolecularOrbitals"].pop("OrbitalLabels"))
    with pytest.raises(AegissError, match="OrbitalLabels"):
        aegiss.analyze(parse_orca_json(variant), _dump(), label="C pz")


def test_a_missing_overlap_is_refused(tmp_path):
    variant = _variant(tmp_path, lambda m: m.pop("S-Matrix"))
    with pytest.raises(AegissError, match="S-Matrix"):
        aegiss.analyze(parse_orca_json(variant), _dump(), label="C pz")


def test_an_unmatched_label_is_refused():
    with pytest.raises(AegissError, match="matches no orbital label"):
        aegiss.analyze(_export(), _dump(), label="Xx d")


def test_a_window_mismatch_is_refused():
    """An export without fractional occupations (an RHF export) cannot pair with a
    CASSCF dump: the window sizes disagree and the analysis refuses."""
    with pytest.raises(AegissError, match="do not belong to the same run"):
        aegiss.analyze(
            parse_orca_json(FIXTURES / "n2_apc.json"), _dump(), label="C pz"
        )


# --- output -------------------------------------------------------------------


def test_the_report_prints_the_screens_and_the_final_space():
    body = aegiss.render(
        aegiss.analyze(_export(), _dump(), label="C pz", reference_energy=ENGINE_ENERGY)
    )
    assert "AEGISS active-space selection" in body
    assert "entropy line" in body and "projection threshold" in body
    assert "final active space: (6e, 6o)" in body
    assert "signed row sum" in body  # the deviation is stated


def test_evidence_carries_the_source():
    bibkeys = {item.bibkey for item in aegiss.evidence() if item.bibkey}
    assert bibkeys == {"tarocco2026aegiss"}
