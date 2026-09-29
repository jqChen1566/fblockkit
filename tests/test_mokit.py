"""Checks of the MOKIT/automr reader and input generator (menus 44/45; plan item 6.1).

Fixtures are real MOKIT 1.2.8 artifacts (101, 2026-09-30; see
fixtures/mokit/README): an H2O CASSCF probe driven with
``mokit{GVB_prog=Gaussian}`` (Gaussian 16 for GVB, PySCF for CASSCF) and an
N2 probe whose input was produced by this module's generator (the
generator's own acceptance run).  The reader's anchors and the generator's
template conventions are measured on these files.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from pytest import approx

from fblockkit.parsers.mokit import MokitError, read_mokit_run_file
from fblockkit.recipe import mokit as mokit_recipe
from fblockkit.recipe.mokit import MokitInputError

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
MO = FIXTURES / "mokit"


def test_the_h2o_probe():
    run = read_mokit_run_file(MO / "h2o_gvb_automr.out")
    assert run.version == "1.2.8"
    assert run.built == "2026-Aug-18"
    assert run.program_paths["gms_path"] == "NOT FOUND"
    assert run.program_paths["gau_path"].endswith("/g16")
    assert (run.memory, run.nproc) == ("4GB", "2")
    assert run.method_basis == "casscf/cc-pvdz"
    assert run.keywords == "gvb_prog=gaussian"
    # the strategy table is reprinted on updates; the last number is effective
    assert run.strategy_number == 1
    assert run.strategy_flags["UNO"] is True
    assert run.strategy_flags["GVB"] is True
    assert run.strategy_flags["CASSCF"] is True
    assert run.strategy_flags["CASPT2"] is False
    assert run.stages == ("do_hf", "get_paired_LMO", "do_gvb", "do_cas")
    chain = dict(run.energies)
    assert chain["UHF"] == approx(-75.82292560)
    assert chain["GVB"] == approx(-75.91806514)
    assert chain["CASSCF"] == approx(-75.90823047)
    assert dict(run.s2_values)["UHF"] == approx(1.066)
    assert (run.gvb_order, run.gvb_program) == (4, "gaussian")
    assert run.active_space == (4, 4)
    assert run.casscf_program == "pyscf"
    # three radical-index tables (after UNO, GVB and CASSCF)
    assert len(run.radical_groups) == 3
    assert run.radical_groups[0][0] == ("biradical", approx(0.109))
    assert run.radical_groups[2][0] == ("biradical", approx(0.165))
    assert run.terminated is True
    assert run.terminated_at == "Wed Sep 30 00:09:49 2026"
    assert "h2o_gvb_uhf_uno_asrot2gvb4_s.fch" in run.side_products


def test_the_generated_h2o_probe():
    """The generator's own acceptance run: its input was produced by menu 44."""
    run = read_mokit_run_file(MO / "h2o_generated_automr.out")
    assert run.terminated is True
    assert run.active_space == (4, 4)  # automatically determined
    assert run.casscf_program == "pyscf"
    assert run.gvb_order == 4


def test_foreign_files_are_refused():
    with pytest.raises(MokitError, match="AutoMR"):
        read_mokit_run_file(FIXTURES / "orca" / "n2_casscf_nevpt2.out")


def test_missing_file_reports_next_step(tmp_path):
    with pytest.raises(MokitError, match="does not exist"):
        read_mokit_run_file(tmp_path / "nope.out")


# --- the generator -------------------------------------------------------------


def test_the_generator_reproduces_its_fixture():
    """The shipped generated input is byte-identical to a fresh generation."""
    atoms = (
        ("O", -0.23497692, 0.90193619, -0.068688),
        ("H", 1.26502308, 0.90193619, -0.068688),
        ("H", -0.73568721, 2.31589843, -0.068688),
    )
    produced = mokit_recipe.input_gjf(
        atoms=atoms,
        method="CASSCF",
        basis="cc-pVDZ",
        mokit_options="GVB_prog=Gaussian",
        charge=0,
        mult=1,
    )
    assert produced == (MO / "h2o_generated.gjf").read_text(encoding="utf-8")


def test_the_generator_normalises_gvb_programs():
    text = mokit_recipe.input_gjf(
        atoms=(("N", 0.0, 0.0, 0.0), ("N", 0.0, 0.0, 1.1)),
        method="CASSCF",
        basis="cc-pVDZ",
        mokit_options="gvb_prog=gamess,ist=1",
        charge=0,
        mult=1,
    )
    assert "mokit{GVB_prog=GAMESS,ist=1}" in text


def test_the_generator_refuses_foreign_gvb_backends():
    with pytest.raises(MokitInputError, match="GAMESS"):
        mokit_recipe.input_gjf(
            atoms=(("N", 0.0, 0.0, 0.0),),
            method="CASSCF",
            basis="cc-pVDZ",
            mokit_options="GVB_prog=PySCF",
            charge=0,
            mult=1,
        )


def test_the_generator_refuses_malformed_options():
    with pytest.raises(MokitInputError, match="key=value"):
        mokit_recipe.input_gjf(
            atoms=(("N", 0.0, 0.0, 0.0),),
            method="CASSCF",
            basis="cc-pVDZ",
            mokit_options="this is not an option",
            charge=0,
            mult=1,
        )


def test_the_generator_refuses_method_with_basis():
    with pytest.raises(MokitInputError, match="basis"):
        mokit_recipe.input_gjf(
            atoms=(("N", 0.0, 0.0, 0.0),),
            method="CASSCF/cc-pVDZ",
            basis="cc-pVDZ",
            mokit_options="",
            charge=0,
            mult=1,
        )
