"""Tests for the DM-AS candidate batch generator (recipe/dm_batch.py, menu 19).

The batch is text: the tests pin the PASS/PASS+ candidate families against the
source's constraints, the generated inputs against the measured ORCA route
(MP2 NatOrbs preparation; !NoIter + moread CASCI with nroots 1), and the
manifest/script against what menu 20 expects.
"""

from __future__ import annotations

import json

import pytest

from fblockkit.recipe import dm_batch
from fblockkit.recipe.dm_batch import DmBatchError, pass_candidates, plan_batch

XYZ = "3\nwater\nO 0.0 0.0 0.11779\nH 0.0 0.757 -0.47116\nH 0.0 -0.757 -0.47116\n"


def _batch(**overrides):
    settings = dict(
        xyz_text=XYZ,
        charge=0,
        multiplicity=1,
        basis="def2-TZVP",
    )
    settings.update(overrides)
    return plan_batch(**settings)


def _pairs(candidates):
    return {(candidate.nel, candidate.norb) for candidate in candidates}


def test_the_pass_family_follows_the_source_constraints():
    candidates = pass_candidates()
    assert len(candidates) == 35  # n_o runs over 9 + 8 + 7 + 6 + 5 values
    pairs = _pairs(candidates)
    assert (6, 6) in pairs  # n_o = n_e/2 + 3 is the first allowed
    assert (6, 9) in pairs and (14, 14) in pairs
    assert not any(nel == 4 for nel, _ in pairs)
    assert (6, 5) not in pairs  # n_o = n_e/2 + 2 is excluded from PASS


def test_the_pass_plus_extras_add_the_flagged_rows():
    added = _pairs(pass_candidates(include_pass_plus=True)) - _pairs(pass_candidates())
    assert (4, 4) in added and (4, 14) in added  # the 4-electron rows
    assert (6, 5) in added and (14, 9) in added  # the two-virtual rows
    assert len(added) == 11 + 5


def test_the_preparation_input_carries_the_mp2_natural_orbitals():
    batch = _batch()
    text = batch.prep_input()
    assert "! RHF def2-TZVP" in text
    assert "NatOrbs true" in text and "Density relaxed" in text
    assert "* xyz 0 1" in text
    hf = _batch(prep="hf").prep_input()
    assert "NatOrbs" not in hf


def test_the_candidate_inputs_are_the_measured_casci_route():
    batch = _batch()
    candidate = batch.candidates[0]
    text = batch.candidate_input(candidate)
    assert "! RHF def2-TZVP NoIter moread" in text
    assert '%moinp "prep.gbw"' in text
    assert f" nel {candidate.nel}" in text and f" norb {candidate.norb}" in text
    assert " nroots 1" in text and " mult 1" in text


def test_the_script_runs_the_preparation_first_and_every_candidate():
    batch = _batch()
    script = batch.render_script()
    assert script.index("prep.inp") < script.index("cand_e6o6.inp")
    for candidate in batch.candidates:
        assert batch.candidate_name(candidate) in script
    assert "menu 20" in script


def test_the_manifest_has_one_entry_per_candidate_and_a_placeholder_reference():
    batch = _batch()
    payload = json.loads(dm_batch.manifest_text(batch))
    assert len(payload["candidates"]) == len(batch.candidates)
    assert payload["candidates"][0]["output"] == "cand_e6o6.out"
    assert "reference" in payload and payload["protocol"] == "gdm"


def test_the_refusals_name_the_next_step():
    with pytest.raises(DmBatchError, match="singlet"):
        _batch(multiplicity=2)
    with pytest.raises(DmBatchError, match="charged"):
        _batch(charge=-1)
    with pytest.raises(DmBatchError, match="preparation level"):
        _batch(prep="ccsd")
    with pytest.raises(DmBatchError, match="no PASS candidate"):
        _batch(max_norb=5)


def test_evidence_carries_the_two_papers_and_the_measured_route():
    text = " ".join(
        entry.ref + " " + entry.text + " " + entry.url + " " + entry.bibkey
        for entry in dm_batch.evidence()
    )
    assert "2c01128" in text and "6c00473" in text
    assert "kaufold2026casci" in text
    assert "NatOrbs" in text and "MaxMacroIter" in text
