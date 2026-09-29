"""Checks of the xTB capture reader (parsers/xtb.py; plan item 6.1, menu 42).

Fixtures are real xTB 6.7.1 captures (101, 2026-09-30; see
fixtures/xtb/README): water single point, water opt+hess (nine modes, zero
imaginary), planar ammonia --hess (one imaginary mode, -1337.66 cm-1) and
water with --cycles 2 (optimisation failed to converge).  Every anchor in
the parser is measured on these files.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from pytest import approx

from fblockkit.parsers.xtb import XtbError, read_xtb_run_file

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
XT = FIXTURES / "xtb"


def test_the_opt_hess_capture():
    run = read_xtb_run_file(XT / "h2o_ohess.out")
    assert run.version == "6.7.1"
    assert run.build == "edcfbbe"
    assert run.compiled_on == "2024-07-22"
    assert run.tasks == ("optimization", "frequencies")
    assert run.energy_Eh == approx(-5.070544013094)
    assert run.gradient_norm == approx(0.000932444846)
    assert run.homo_lumo_gap_eV == approx(14.397675893256)
    assert run.opt_converged is True
    assert run.opt_cycles == 5
    # the SCF marker occurs once per SCF run; the capture carries two
    assert run.scf_cycles == (3, 3)
    # the frequency block is printed twice; both groups are kept as measured
    assert [len(group) for group in run.frequency_groups] == [9, 9]
    assert run.frequencies == run.frequency_groups[1]
    assert run.n_imaginary == 0
    assert run.engine_n_frequencies == 3
    assert run.engine_n_imaginary == 0
    # the three real modes of water
    real = [value for value in run.frequencies if abs(value) > 1.0]
    assert real == approx([1538.60, 3642.93, 3657.89])
    assert run.thermo["free_energy_Eh"] == approx(-5.068023573707)
    assert run.thermo["zpe_Eh"] == approx(0.020137684918)
    assert run.finished == "2026/09/29 23:33:55.899"


def test_the_imaginary_mode_capture():
    run = read_xtb_run_file(XT / "nh3_planar_hess.out")
    assert run.tasks == ("frequencies",)
    assert run.n_imaginary == 1
    assert run.engine_n_imaginary == 1
    assert run.imag_found == 1
    imaginary = [value for value in run.frequencies if value < 0]
    assert min(imaginary) == approx(-1337.66)


def test_the_failed_optimisation_capture():
    run = read_xtb_run_file(XT / "h2o_noconv.out")
    assert run.tasks == ("optimization",)
    assert run.opt_converged is False
    assert run.opt_cycles == 2
    assert run.energy_Eh == approx(-5.070336253370)


def test_the_single_point_capture():
    run = read_xtb_run_file(XT / "h2o_sp.out")
    assert run.tasks == ("single_point",)
    assert run.frequency_groups == ()
    assert run.n_imaginary == 0
    assert run.energy_Eh == approx(-5.070425872940)
    assert run.gradient_norm == approx(0.024265751947)


def test_foreign_files_are_refused():
    with pytest.raises(XtbError, match="xtb version"):
        read_xtb_run_file(FIXTURES / "orca" / "n2_casscf_nevpt2.out")


def test_missing_file_reports_next_step(tmp_path):
    with pytest.raises(XtbError, match="does not exist"):
        read_xtb_run_file(tmp_path / "nope.out")
