"""Tests for the cross-run state tracker (analysis/state_tracking.py, menu 49).

The fixture sequence is the N2/def2-SVP pair already in the tree: the
state-averaged CAS(6,6) export (three roots) and the single-root CASSCF
export (one root) of the same geometry.  Pinned anchors:

- the 1-RDM convention ``Gamma = C (S D_json S) C^T``: traces equal the
  electron count (14) and the same-run self-difference vanishes (the
  rotation of a run onto itself is the identity);
- the within-run Frobenius differences: D(0,1) = 0.90035, D(0,2) = 0.88714 and
  D(1,2) = 0.01325 (the near-degenerate pi pair -- the metric's hard case);
- the cross-run difference of the same physical state measures 0.05212 after
  the rotation ``M = C_cand S C_target^T`` (orthogonal to 2e-14), while the
  wrong starting states measure 0.907 / 0.894 -- the tracking decision, with
  the margin reported beside it;
- the energies come from the sibling ``.out`` CAS-SCF STATES block (the
  parser's measured section) and ``W0 = (E_t - E_k)^2`` enters Q when both
  energies are known; an incomplete energy set falls back to the density
  ordering;
- refusals: a second run at another geometry (perturbed coordinates), with
  another element (the atom gate), with another basis, a missing per-root
  density sidecar, a missing S matrix, and a target root that does not exist.
"""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pytest

from fblockkit.analysis import state_tracking
from fblockkit.analysis.state_tracking import StateTrackingError, build_run, track
from fblockkit.parsers.base import parse_auto
from fblockkit.parsers.orca_json import parse_orca_json

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"
SA = FIXTURES / "n2_ass1st_sa.json"
SA_OUT = FIXTURES / "n2_ass1st_sa.out"
SINGLE = FIXTURES / "n2_ass1st.json"
SINGLE_OUT = FIXTURES / "n2_ass1st.out"

# the block's own ROOT lines (not the initial-state snapshot):
SA_ENERGIES = {
    (1, 0): -108.9409927681,
    (1, 1): -108.5567547967,
    (1, 2): -108.5267602967,
}
SINGLE_ENERGY = -108.9506719453


def _energies(path: Path) -> dict[tuple[int, int], float]:
    """Per-root energies from the sibling .out, as the handler assembles them."""
    parsed = parse_auto(path)
    table: dict[tuple[int, int], float] = {}
    for block in parsed.sections["casscf"]["states"]:
        for root in block["roots"]:
            table[(block["mult"], root["root"])] = root["energy"]
    return table


def _variant(tmp_path: Path, mutate: Callable[[dict[str, Any]], None]) -> Path:
    document = json.loads(SA.read_text(encoding="utf-8"))
    mutate(document["Molecule"])
    target = tmp_path / "variant.json"
    target.write_text(json.dumps(document), encoding="utf-8")
    return target


# --- the run builder ----------------------------------------------------------


def test_the_builder_reads_the_three_roots_with_the_electron_count():
    run = build_run(parse_orca_json(SA), name="sa")
    assert run.root_labels() == ((1, 0), (1, 1), (1, 2))
    assert (run.n_mo, run.n_ao) == (28, 28)
    for root in run.roots:
        gamma = np.asarray(root.gamma)
        assert np.trace(gamma) == pytest.approx(14.0, abs=1e-9)  # N2: 14 electrons
        assert np.abs(gamma - gamma.T).max() < 1e-10
        assert root.energy is None  # no energies were handed in


def test_the_builder_refuses_without_densities(tmp_path):
    def mutate(molecule: dict[str, Any]) -> None:
        for key in [k for k in molecule["Densities"] if k.startswith("Tdens-CAS.")]:
            molecule["Densities"].pop(key)

    variant = _variant(tmp_path, mutate)
    with pytest.raises(StateTrackingError, match="Tdens-CAS"):
        build_run(parse_orca_json(variant), name="variant")


def test_the_builder_refuses_without_the_overlap(tmp_path):
    variant = _variant(tmp_path, lambda molecule: molecule.pop("S-Matrix"))
    with pytest.raises(StateTrackingError, match="S matrix"):
        build_run(parse_orca_json(variant), name="variant")


# --- the walk on one run (self-rotation) --------------------------------------


def test_the_same_run_scores_its_own_roots():
    """Target root 0 against its own run: the self-difference vanishes and the
    other roots separate by the measured Frobenius values."""
    run = build_run(parse_orca_json(SA), name="sa")
    result = track([run, run], target_root=0)
    step = result.steps[0]
    scores = {score.root: score for score in step.candidates}
    assert scores[0].d == pytest.approx(0.0, abs=1e-10)
    assert scores[1].d == pytest.approx(0.90035, abs=1e-4)
    assert scores[2].d == pytest.approx(0.88714, abs=1e-4)
    assert scores[0].chosen and not scores[1].chosen and not scores[2].chosen
    # the runner-up is root 2 (0.88714); the margin is its difference
    assert step.margin == pytest.approx(0.88714, abs=1e-4)
    assert not step.energy_used


def test_the_pi_pair_separates_but_with_a_small_margin():
    """Target root 1: the near-degenerate pi partner (D = 0.01325) is told
    apart from root 0 (0.90035), with the small margin reported."""
    run = build_run(parse_orca_json(SA), name="sa")
    result = track([run, run], target_root=1)
    step = result.steps[0]
    scores = {score.root: score for score in step.candidates}
    assert scores[1].chosen
    assert scores[2].d == pytest.approx(0.01325, abs=1e-4)
    assert step.margin == pytest.approx(0.01325, abs=1e-4)


# --- the cross-run walk -------------------------------------------------------


def test_the_rotation_between_runs_is_orthogonal():
    sa = build_run(parse_orca_json(SA), name="sa")
    sg = build_run(parse_orca_json(SINGLE), name="single")
    c_target = np.asarray(sa.coefficients)
    c_cand = np.asarray(sg.coefficients)
    overlap = np.asarray(sg.overlap)
    rotation = c_cand @ overlap @ c_target.T
    assert np.abs(rotation @ rotation.T - np.eye(rotation.shape[0])).max() < 1e-10
    # one geometry, one S matrix (the fixtures' own gate)
    assert np.abs(np.asarray(sa.overlap) - overlap).max() < 1e-12


def test_the_same_state_tracks_across_runs_and_the_wrong_ones_separate():
    sa = build_run(parse_orca_json(SA), name="sa")
    sg = build_run(parse_orca_json(SINGLE), name="single")
    result = track([sa, sg], target_root=0)
    step = result.steps[0]
    assert len(step.candidates) == 1
    score = step.candidates[0]
    assert score.root == 0 and score.chosen
    assert score.d == pytest.approx(0.05212, abs=5e-5)
    assert step.margin is None  # a single candidate
    wrong = track([sa, sg], target_root=2)
    assert wrong.steps[0].candidates[0].d == pytest.approx(0.89401, abs=1e-3)


def test_the_energies_come_from_the_sibling_out_and_make_w0():
    assert _energies(SA_OUT) == pytest.approx(SA_ENERGIES)
    assert _energies(SINGLE_OUT)[(1, 0)] == pytest.approx(SINGLE_ENERGY)
    sa = build_run(parse_orca_json(SA), energies=_energies(SA_OUT), name="sa")
    sg = build_run(parse_orca_json(SINGLE), energies=_energies(SINGLE_OUT), name="single")
    result = track([sa, sg], target_root=0)
    step = result.steps[0]
    assert step.energy_used
    score = step.candidates[0]
    expected_w0 = (SA_ENERGIES[(1, 0)] - SINGLE_ENERGY) ** 2
    assert score.w0 == pytest.approx(expected_w0, rel=1e-6)
    assert score.q == pytest.approx(score.w0 + score.d)


def test_an_incomplete_energy_set_falls_back_to_the_density_ordering():
    sa = build_run(parse_orca_json(SA), name="sa")  # target energy unknown
    sg = build_run(parse_orca_json(SINGLE), energies=_energies(SINGLE_OUT))
    result = track([sa, sg], target_root=0)
    step = result.steps[0]
    assert not step.energy_used
    assert step.candidates[0].w0 is None
    assert step.candidates[0].q == step.candidates[0].d


# --- the average-size walk (SA 3 roots -> SA 4 roots) -------------------------


def test_the_sequence_walk_across_the_average_size():
    """SA(3 roots) -> SA(4 roots): the ground and first two excited states
    match with small differences; the printed-degenerate pair (roots 2/3 in
    the SA4 run) collapses the margin, the metric's hard case on record."""
    sa = build_run(parse_orca_json(SA), energies=_energies(SA_OUT), name="sa")
    sa4 = build_run(
        parse_orca_json(FIXTURES / "n2_ass1st_sa4.json"),
        energies=_energies(FIXTURES / "n2_ass1st_sa4.out"),
        name="sa4",
    )
    assert sa4.root_labels() == ((1, 0), (1, 1), (1, 2), (1, 3))
    zero = track([sa, sa4], target_root=0)
    scores = {score.root: score for score in zero.steps[0].candidates}
    assert scores[0].chosen
    assert scores[0].d == pytest.approx(7.2159e-03, abs=1e-5)
    one = track([sa, sa4], target_root=1)
    scores = {score.root: score for score in one.steps[0].candidates}
    assert scores[1].chosen
    assert scores[1].d == pytest.approx(1.2558e-02, abs=1e-5)
    assert one.steps[0].margin == pytest.approx(6.68e-03, abs=1e-4)
    two = track([sa, sa4], target_root=2)
    scores = {score.root: score for score in two.steps[0].candidates}
    assert scores[2].d == pytest.approx(scores[3].d, abs=1e-12)
    assert two.steps[0].margin == pytest.approx(0.0, abs=1e-12)
    body = state_tracking.render(two, source="n2_ass1st_sa.json")
    assert "nearly degenerate in this metric" in body


# --- refusals -----------------------------------------------------------------


def test_a_perturbed_geometry_is_refused(tmp_path):
    def mutate(molecule: dict[str, Any]) -> None:
        molecule["Atoms"][1]["Coords"] = [0.0, 0.0, 1.095]

    variant = _variant(tmp_path, mutate)
    sa = build_run(parse_orca_json(SA), name="sa")
    other = build_run(parse_orca_json(variant), name="shifted")
    with pytest.raises(StateTrackingError, match="Angstrom"):
        track([sa, other], target_root=0)


def test_another_element_is_refused(tmp_path):
    def mutate(molecule: dict[str, Any]) -> None:
        molecule["Atoms"][0]["ElementLabel"] = "C"

    variant = _variant(tmp_path, mutate)
    sa = build_run(parse_orca_json(SA), name="sa")
    other = build_run(parse_orca_json(variant), name="carbon")
    with pytest.raises(StateTrackingError, match="atoms"):
        track([sa, other], target_root=0)


def test_another_ao_dimension_is_refused():
    sa = build_run(parse_orca_json(SA), name="sa")
    other = replace(sa, name="other", n_ao=40, n_mo=40)
    with pytest.raises(StateTrackingError, match="AO space"):
        track([sa, other], target_root=0)


def test_a_missing_target_root_names_the_available_ones():
    sa = build_run(parse_orca_json(SA), name="sa")
    sg = build_run(parse_orca_json(SINGLE), name="single")
    with pytest.raises(StateTrackingError, match="\\(mult 1, root 2\\)"):
        track([sa, sg], target_root=7)


def test_a_single_run_is_refused():
    sa = build_run(parse_orca_json(SA), name="sa")
    with pytest.raises(StateTrackingError, match="at least two runs"):
        track([sa], target_root=0)


# --- output -------------------------------------------------------------------


def test_the_report_prints_the_lineage_and_the_declarations():
    sa = build_run(parse_orca_json(SA), energies=_energies(SA_OUT), name="n2_ass1st_sa")
    sg = build_run(parse_orca_json(SINGLE), energies=_energies(SINGLE_OUT), name="n2_ass1st")
    result = track([sa, sg], target_root=0)
    body = state_tracking.render(result, source="n2_ass1st_sa.json")
    assert "Cross-run state tracking" in body
    assert "Step 1 -> 2 (n2_ass1st)" in body
    assert "[chosen]" in body
    assert "a single candidate" in body
    assert "Boundaries" in body and "W1" in body
    assert "outside" in body  # the cross-geometry scope note
    # a multi-candidate step prints the margin instead
    multi = track([sa, sa], target_root=0)
    body_multi = state_tracking.render(multi, source="n2_ass1st_sa.json")
    assert "margin over the runner-up" in body_multi


def test_evidence_carries_the_source_paper():
    bibkeys = {item.bibkey for item in state_tracking.evidence() if item.bibkey}
    assert bibkeys == {"tran2019tracking"}
