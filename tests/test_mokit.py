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

from fblockkit.analysis import mokit as mokit_analysis
from fblockkit.parsers import ParserError
from fblockkit.parsers.mokit import MokitError, read_mokit_run_file
from fblockkit.parsers.mokit_fch import read_fch, stage_label
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


# --- the .fch side products ---------------------------------------------------


def test_the_fch_reader_pins_the_rhf_stage():
    """The layout anchor: the parsed coefficients (column-major) reconstruct the
    file's own stored density triangle (row-major lower) to 1e-8."""
    rhf = read_fch(MO / "h2o_gvb_rhf.fch")
    assert rhf.nbf == 24 and rhf.charge == 0 and rhf.multiplicity == 1
    assert (rhf.n_alpha, rhf.n_beta, rhf.n_electrons) == (5, 5, 10)
    assert rhf.alpha_energies[0] == approx(-20.6296358, abs=1e-6)
    assert rhf.coordinates_bohr[0][0] == approx(-0.444042024, abs=1e-9)
    assert rhf.beta_energies is None and rhf.beta_coefficients is None
    assert rhf.scf_energy == approx(-75.78429272840782, abs=1e-12)
    nocc = rhf.n_electrons // 2
    pairs = [(p, q) for p in range(rhf.nbf) for q in range(p + 1)]
    assert rhf.total_scf_density is not None
    worst = max(
        abs(
            2.0
            * sum(
                rhf.alpha_coefficients[p][m] * rhf.alpha_coefficients[q][m]
                for m in range(nocc)
            )
            - rhf.total_scf_density[k]
        )
        for k, (p, q) in enumerate(pairs)
    )
    assert worst < 1e-8


def test_the_fch_reader_pins_the_transformed_stage():
    no = read_fch(MO / "h2o_gvb_uhf_gvb4_CASSCF_NO.fch")
    assert no.total_energy == approx(-75.91806513607675, abs=1e-12)
    assert no.beta_energies is None
    assert no.dipole_au is not None
    assert no.dipole_au[0] == approx(0.307792233, abs=1e-9)
    assert len(no.alpha_energies) == no.nbf
    assert stage_label(no.name) == "the CASSCF natural orbitals"
    assert stage_label("h2o_gvb_uhf_uno_asrot2gvb4.fch").startswith(
        "the UNO active-space rotation"
    )
    assert stage_label("mystery.fch") == "mystery"


def test_the_fch_reader_refuses(tmp_path):
    with pytest.raises(ParserError, match="does not exist"):
        read_fch(tmp_path / "missing.fch")
    junk = tmp_path / "junk.fch"
    junk.write_text("just a line\n", encoding="utf-8")
    with pytest.raises(ParserError, match="title/level"):
        read_fch(junk)
    incomplete = tmp_path / "incomplete.fch"
    incomplete.write_text("title\nlevel\nSomething  I  1\n", encoding="utf-8")
    with pytest.raises(ParserError, match="Number of basis functions"):
        read_fch(incomplete)


def test_the_report_renders_the_fch_section():
    run = read_mokit_run_file(MO / "h2o_gvb_automr.out")
    files = (
        read_fch(MO / "h2o_gvb_rhf.fch"),
        read_fch(MO / "h2o_gvb_uhf_uno_asrot2gvb4.fch"),
    )
    body = mokit_analysis.render(run, source="h2o_gvb_automr.out", fch_files=files)
    assert ".fch side products read next to the output (2 file(s)):" in body
    assert "name-derived" in body and "no beta block" in body
    assert "-75.78429273 Eh" in body
    body2 = mokit_analysis.render(
        run, source="x.out", fch_files=("broken.fch: not read (short file)",)
    )
    assert "broken.fch: not read" in body2


def test_a_truncated_array_is_refused(tmp_path):
    """Dropping one value from an array must raise -- the reader stops at the
    next section instead of swallowing its numbers as data (the overread
    that silently fabricated a coefficient with the next header's count,
    found by the 2026-10-03 review)."""
    source = MO / "h2o_gvb_rhf.fch"
    lines = source.read_text(encoding="utf-8").splitlines(keepends=True)
    start = next(
        i for i, ln in enumerate(lines) if ln.startswith("Alpha Orbital Energies")
    )
    end = next(
        i
        for i in range(start + 1, len(lines))
        if lines[i].strip() and lines[i].lstrip()[0] not in "0123456789+-."
    )
    tokens = lines[end - 1].split()
    lines[end - 1] = " ".join(tokens[:-1]) + "\n"
    broken = tmp_path / "truncated.fch"
    broken.write_text("".join(lines), encoding="utf-8")
    with pytest.raises(
        ParserError, match="Alpha Orbital Energies.*were read"
    ):
        read_fch(broken)
