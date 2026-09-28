"""Checks of the pysisyphus B-layer pair (menus 32/33; Wave 4.6).

The fixtures are real runs of pysisyphus 1.0.0 driving ORCA 6.1.1 on server
101 (2026-09-28; see fixtures/pysisyphus/README.md for the provenance and the
measured format conventions).  The four runs cover the whole outcome space the
reader must classify: a converged minimum optimisation, a converged TS search
whose cycle table restarts after ten rows, a run stopped by its cycle limit
(with the extra closing evaluation), and a run that crashed in its final-
Hessian post-processing -- the measured ORCA-6 ``$multiplicity`` parse stop.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from fblockkit.analysis import pysisyphus_run as prun
from fblockkit.diagnosis import references_section
from fblockkit.parsers import pysisyphus as pyparse
from fblockkit.parsers.pysisyphus import PysisyphusError as ArtifactError
from fblockkit.recipe.pysisyphus import PysisyphusError as InputError
from fblockkit.recipe import pysisyphus as recipe

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "pysisyphus"
H2O_OPT = FIXTURES / "h2o_opt"
TS = FIXTURES / "butadiene_ts"
STOP3 = FIXTURES / "h2o_stop3"
CRASH = FIXTURES / "hess_crash"


# --- the format layer --------------------------------------------------------


def test_the_console_facts_of_the_converged_run():
    facts = pyparse.read_log((H2O_OPT / "run_stdout.log").read_text(encoding="utf-8"))
    assert facts["kind"] == "OPTIMIZATION"
    assert facts["version"] == "1.0.0"
    assert facts["system"] == "H₂O" and facts["n_atoms"] == 3
    assert facts["coord_type"] == "redund" and facts["n_coord"] == 3
    assert (facts["calculator"], facts["calculator_id"]) == ("ORCA", "calculator_000")
    assert (facts["charge"], facts["multiplicity"]) == (0, 1)
    assert facts["optimizer"] == "rfo"
    assert facts["outcome"] == "converged"
    assert facts["final"]["energy"] == pytest.approx(-75.96133849, abs=1e-12)
    assert facts["final_geometry_fn"] == "final_geometry.xyz"
    assert facts["citation"] == "https://doi.org/10.1002/qua.26390"
    # the four thresholds, with the overachieved values on the force lines only
    assert facts["thresholds"]["max_force"] == (0.00045, 9e-05)
    assert facts["thresholds"]["rms_force"] == (0.0003, 6e-05)
    assert facts["thresholds"]["max_step"] == (0.0018, None)
    assert facts["thresholds"]["rms_step"] == (0.0012, None)
    assert len(facts["cycles"]) == 4
    # the first row is the nan* energy change of the starting cycle
    assert facts["cycles"][0].d_energy is None
    assert "d_energy" in facts["cycles"][0].starred
    # the closing row is starred on every column (the run converged there)
    assert facts["cycles"][-1].starred == {
        "d_energy",
        "max_force",
        "rms_force",
        "max_step",
        "rms_step",
    }


def test_the_cycle_table_restarts_after_ten_rows():
    # measured on the TS run: rows 0-9, a dashes separator, rows 10-19
    facts = pyparse.read_log((TS / "run_stdout.log").read_text(encoding="utf-8"))
    assert facts["kind"] == "TS-OPTIMIZATION"
    assert facts["optimizer"] == "rsprfo"
    assert [row.cycle for row in facts["cycles"]] == list(range(20))
    assert facts["cycles"][-1].d_energy == pytest.approx(-1e-06, abs=1e-12)


def test_the_stopped_run_reports_its_own_marker():
    facts = pyparse.read_log((STOP3 / "run_stdout.log").read_text(encoding="utf-8"))
    assert facts["outcome"] == "cycles exceeded"
    assert len(facts["cycles"]) == 3
    # the closing forces already sit below the thresholds; the marker is the
    # outcome, so the reader must not "converge" on the force values
    assert facts["final"]["max_forces_internal"] <= facts["thresholds"]["max_force"][0]


def test_the_crash_capture_carries_the_backup_and_the_exceptions():
    facts = pyparse.read_log((CRASH / "crash_stdout.log").read_text(encoding="utf-8"))
    assert facts["outcome"] == "crashed"
    assert facts["crash"]["converged_before_crash"] is True
    assert facts["crash"]["backup"].endswith("crashed_calculator_000")
    assert any("ParseException" in line for line in facts["crash"]["exceptions"])
    assert any(
        "RunAfterCalculationFailedException" in line
        for line in facts["crash"]["exceptions"]
    )


def test_the_trajectory_frames_and_the_comment_energy():
    frames = pyparse.read_trj((H2O_OPT / "optimization.trj").read_text(encoding="utf-8"))
    assert len(frames) == 4
    assert frames[0].energy == pytest.approx(-75.96073698, abs=1e-12)
    assert frames[-1].energy == pytest.approx(-75.96133849, abs=1e-12)
    assert frames[-1].atoms[0] == ("O", -0.00316071, -0.01122882, 0.0)
    assert frames[-1].comment == ""
    closing = pyparse.read_xyz(
        (H2O_OPT / "final_geometry.xyz").read_text(encoding="utf-8")
    )
    assert closing.energy == frames[-1].energy
    assert closing.atoms == frames[-1].atoms


def test_the_run_record_is_the_program_s_own_yaml():
    record = pyparse.read_run_record((TS / "RUN.yaml").read_text(encoding="utf-8"))
    assert record["version"] == "1.0.0"
    assert record["calc"]["type"] == "orca5"
    assert record["calc"]["keywords"] == "hf 3-21g"
    assert record["tsopt"]["type"] == "rsprfo"
    assert record["tsopt"]["hessian_init"] == "fischer"
    assert record["tsopt"]["thresh"] == "baker"


def test_refusals_carry_next_steps(tmp_path):
    with pytest.raises(ArtifactError, match="RUNNING"):
        pyparse.read_log("this is not a pysisyphus console capture\n")
    with pytest.raises(ArtifactError, match="atom count"):
        pyparse.read_trj("not a trajectory at all\n")
    # a directory without a console capture is refused by name
    with pytest.raises(ArtifactError, match="console banner"):
        pyparse.find_run_log(tmp_path)
    # two captures in one directory: the reader wants one run per capture
    (tmp_path / "a.out").write_text("x\n# RUNNING OPTIMIZATION #\n", encoding="utf-8")
    (tmp_path / "b.out").write_text("y\n# RUNNING OPTIMIZATION #\n", encoding="utf-8")
    with pytest.raises(ArtifactError, match="one run per capture"):
        pyparse.find_run_log(tmp_path)


# --- the analysis layer ------------------------------------------------------


def test_the_converged_run_cross_checks_all_clean():
    data = prun.read_run(H2O_OPT)
    assert [check.ok for check in data.checks] == [True] * 5
    assert data.orca_calls[-1].energy == pytest.approx(-75.961338488, abs=1e-9)
    assert data.orca_calls[-1].terminated is True


def test_the_ts_run_reproduces_the_upstream_assertion():
    data = prun.read_run(TS)
    assert [check.ok for check in data.checks] == [True] * 5
    # pysisyphus's own example asserts -154.050455732882 for this case; this run
    # (a model Hessian instead of the quantum one) lands within 1e-7
    assert data.facts["final"]["energy"] == pytest.approx(-154.050455732882, abs=1e-7)


def test_the_stopped_run_binds_its_closing_state_to_the_extra_call():
    data = prun.read_run(STOP3)
    outcomes = {check.name: check for check in data.checks}
    # the closing state is the extra evaluation after the last cycle ...
    assert outcomes["closing state (final summary) vs the last calculator call"].ok
    assert data.orca_calls[-1].path.name == "calculator_000.003.orca.out"
    # ... and its geometry has no trajectory frame (measured, said so honestly)
    geometry_check = outcomes["closing geometry coordinates vs the closing frame"]
    assert geometry_check.ok is None
    assert "extra evaluation" in geometry_check.detail


def test_the_crash_is_classified_with_the_multiplicity_signature():
    data = prun.read_run(CRASH)
    assert data.crash_kind == "hessian-parse"
    assert data.crash_hess is not None
    assert "$multiplicity" in data.crash_hess.read_text(encoding="utf-8")
    body = prun.render(data)
    assert "$multiplicity" in body
    assert "signature is confirmed" in body
    assert "hessian_init: fischer" in body
    # the optimisation part itself was converged before the crash
    assert "already printed its converged marker" in body


def test_the_report_is_deterministic_and_names_the_artifacts():
    data = prun.read_run(H2O_OPT)
    assert prun.render(data) == prun.render(data)
    body = prun.render(data)
    assert "trajectory (4 frames)" in body
    assert "optimization.h5" in body  # named, and said not to be read
    assert "4 frame(s) vs 4 cycle row(s)" in body and ": ok" in body


def test_the_evidence_is_complete_and_citable():
    for evidence in (prun.evidence(), recipe.evidence()):
        refs = references_section(evidence)
        assert refs is not None
    literature = [item for item in prun.evidence() if item.kind == "literature"]
    assert literature and "10.1002/qua.26390" in literature[0].ref


# --- the generator -----------------------------------------------------------


def _load(yaml_text: str) -> dict:
    parsed = yaml.safe_load(yaml_text)
    assert isinstance(parsed, dict)
    return parsed


def test_the_minimum_input_matches_the_measured_schema():
    text = recipe.input_yaml(xyz_fn="h2o.pysisyphus.xyz", keywords="HF def2-SVP", pal=2)
    record = _load(text)
    assert record["geom"] == {"type": "redund", "fn": "h2o.pysisyphus.xyz"}
    assert record["calc"]["type"] == "orca5"
    assert record["calc"]["keywords"] == "HF def2-SVP"
    assert record["calc"]["pal"] == 2 and record["calc"]["mem"] == 1500
    assert record["opt"]["type"] == "rfo" and record["opt"]["thresh"] == "gau"
    assert "tsopt" not in record


def test_the_ts_input_carries_a_model_hessian_and_the_rx_modes():
    text = recipe.input_yaml(
        xyz_fn="x.xyz",
        keywords="HF 3-21G",
        job="ts",
        rx_modes="[[[[DIHEDRAL, 2, 0, 1, 3], 1]]]",
    )
    record = _load(text)
    assert record["tsopt"]["type"] == "rsprfo"
    assert record["tsopt"]["thresh"] == "baker"
    assert record["tsopt"]["hessian_init"] == "fischer"
    assert record["tsopt"]["rx_modes"] == [[[["DIHEDRAL", 2, 0, 1, 3], 1]]]
    assert record["tsopt"]["do_hess"] is False
    assert "# model Hessian" in text


def test_the_quantum_hessian_is_refused_with_the_measured_reason():
    with pytest.raises(InputError, match="multiplicity"):
        recipe.input_yaml(xyz_fn="x.xyz", keywords="HF def2-SVP", job="ts",
                          hessian_init="calc")


def test_the_generator_refusals_carry_next_steps():
    cases = (
        (dict(xyz_fn="x.xyz", keywords=""), "empty"),
        (dict(xyz_fn="x.xyz", keywords="PBE0\nD4"), "line break"),
        (dict(xyz_fn="x.xyz", keywords="PBE0", job="irc"), "unknown job"),
        (dict(xyz_fn="x.xyz", keywords="PBE0", thresh="never"), "unknown threshold"),
        (dict(xyz_fn="x.xyz", keywords="PBE0", coord_type="zmat"), "coordinate system"),
        (dict(xyz_fn="x.xyz", keywords="PBE0", pal=0), "positive integer"),
        (dict(xyz_fn="x.xyz", keywords="PBE0", mem=10), "at least 100"),
        (dict(xyz_fn="x.xyz", keywords="PBE0", charge=0.5), "integer"),
        (dict(xyz_fn="x.xyz", keywords="PBE0", mult=0), "positive integer"),
        (dict(xyz_fn="x.xyz", keywords="PBE0", max_cycles=0), "positive integer"),
        (dict(xyz_fn="x.xyz", keywords="PBE0", job="ts", hessian_init="xtb"), "unknown hessian_init"),
        (dict(xyz_fn="x.xyz", keywords="PBE0", job="ts", rx_modes="DIHEDRAL 2 0 1 3"), "flow sequence"),
        (dict(xyz_fn="x.xyz", keywords="PBE0", rx_modes="[[BOND, 0, 1]]"), "belongs to a TS"),
        (dict(xyz_fn="x.xyz", keywords="Ångström"), "ASCII"),
    )
    for kwargs, fragment in cases:
        with pytest.raises(InputError, match=fragment):
            recipe.input_yaml(**kwargs)


def test_the_structure_writer_and_the_orca_input_reader():
    atoms = (("O", 0.0, 0.0, 0.0), ("H", 0.96, 0.0, 0.0))
    text = recipe.structure_xyz(atoms, comment="test")
    lines = text.splitlines()
    assert lines[0] == "2" and lines[1] == "test"
    assert lines[2].split()[0] == "O"
    with pytest.raises(InputError, match="no atoms"):
        recipe.structure_xyz(())
    block = "*xyz -1 2\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\n*\n"
    assert recipe.atoms_from_orca_input(block)[1][0] == "H"
    with pytest.raises(InputError, match="no inline"):
        recipe.atoms_from_orca_input("!HF def2-SVP\n")


def test_the_guidance_names_the_measured_boundary():
    text = "\n".join(recipe.run_guidance_lines())
    assert "pysisyphusrc" in text
    assert "TMPDIR" in text
    assert "menu 33" in text
    assert "outside pysisyphus" in text
