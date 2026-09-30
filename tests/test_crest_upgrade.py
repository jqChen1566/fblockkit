"""Tests for the CREST-ensemble upgrade writer (recipe/crest_upgrade.py).

The real-data anchors are the n-butane CREST fixture ensemble (two
conformers, anti/gauche): the generated first input is pinned byte for byte
and is the one whose engine run (``fixtures/crest/upgrade_acceptance.out``)
converged to -158.392455872472 Eh in five cycles (recorded in the module
evidence and the fixture README).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from fblockkit.parsers.crest import read_crest_directory
from fblockkit.recipe import crest_upgrade

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "crest"

EXPECTED_HEAD = """! Opt r2SCAN-3c
%maxcore 2000
* xyz 0 1
"""


def _ensemble():
    return read_crest_directory(FIXTURES)


def test_the_upgrade_batch_writes_every_conformer(tmp_path):
    plan = crest_upgrade.write_upgrade_inputs(_ensemble(), tmp_path / "upgrade")
    assert [item.name for item in plan.inputs] == ["conf_01.opt.inp", "conf_02.opt.inp"]
    assert plan.method == "r2SCAN-3c" and plan.charge == 0 and plan.multiplicity == 1
    assert plan.inputs[0].relative_kcal == pytest.approx(0.000, abs=1e-3)
    assert plan.inputs[1].relative_kcal == pytest.approx(0.596, abs=1e-3)
    assert plan.inputs[0].energy_Eh == pytest.approx(-13.66512758, abs=1e-8)
    for item in plan.inputs:
        written = (tmp_path / "upgrade" / item.name).read_text(encoding="utf-8")
        assert written == item.text
        assert written.startswith(EXPECTED_HEAD)
        assert written.rstrip().endswith("*")
        body = written.splitlines()
        assert any(line.startswith("C ") for line in body)
        assert any(line.startswith("H ") for line in body)


def test_the_count_takes_the_lowest_conformers(tmp_path):
    plan = crest_upgrade.write_upgrade_inputs(_ensemble(), tmp_path / "upgrade", count=1)
    assert [item.name for item in plan.inputs] == ["conf_01.opt.inp"]
    assert not (tmp_path / "upgrade" / "conf_02.opt.inp").exists()


def test_the_method_and_charge_travel_into_the_input(tmp_path):
    plan = crest_upgrade.write_upgrade_inputs(
        _ensemble(), tmp_path / "upgrade", method="PBE0 D4 def2-TZVP", charge=-1,
        multiplicity=2,
    )
    assert plan.inputs[0].text.startswith("! Opt PBE0 D4 def2-TZVP\n")
    assert "* xyz -1 2" in plan.inputs[0].text


def test_the_refusals(tmp_path):
    ensemble = _ensemble()
    with pytest.raises(crest_upgrade.UpgradeError, match="method line is empty"):
        crest_upgrade.write_upgrade_inputs(ensemble, tmp_path, method="  ")
    with pytest.raises(crest_upgrade.UpgradeError, match="outside the ensemble"):
        crest_upgrade.write_upgrade_inputs(ensemble, tmp_path, count=3)
    with pytest.raises(crest_upgrade.UpgradeError, match="outside the ensemble"):
        crest_upgrade.write_upgrade_inputs(ensemble, tmp_path, count=0)


def test_the_plan_render_states_the_batch_and_boundaries(tmp_path):
    plan = crest_upgrade.write_upgrade_inputs(_ensemble(), tmp_path / "upgrade")
    body = crest_upgrade.render_plan(plan)
    assert "conf_01.opt.inp" in body and "conf_02.opt.inp" in body
    assert "orca conf_01.opt.inp" in body
    assert "pre-screen ensemble's, used as-is" in body
    assert "registered as the next increment" in body


def test_evidence_states_the_acceptance_run():
    text = " ".join(entry.text + " " + entry.ref for entry in crest_upgrade.evidence())
    assert "HURRAY" in text
    assert "-158.392455872472" in text
    assert "upgrade_acceptance.out" in text
