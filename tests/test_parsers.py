"""Parser-layer tests: field extraction from real ORCA 6.1.1 output fixtures
(fixtures/orca/).

All asserted values are measured from the fixtures (provenance and generation in
fixtures/orca/README.md). Two regression cases cover real defects the fixtures
caught during development: the SCF iteration regex matching numeric tables by
mistake, and the CASSCF convergence marker accepting only the ENERGY criterion.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from fblockkit.parsers import ParserError, available, detect_program, facts_from, parse_auto
from fblockkit.parsers.orca import OrcaParser

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"

ALL_FIXTURES = (
    "n2_casscf_nevpt2.out",
    "n2_caspt2.out",
    "co_plus_soc.out",
    "co_plus_soc_min.out",
    "scf_noconv.out",
    "fblock_dft_la_complex.out",
    "fblock_dft_gd_crash.out",
    "generated_ce3_sarc2.out",
    "n2_hf_clean.out",
    "n2_casscf_orbcomp.out",
    "ce3_orbcomp.out",
    "n2_ccsd.out",
    "f2_ccsd.out",
    "n2_stretch_ccsd.out",
    "h2o_freq_min.out",
    "nh3_planar_optts_freq.out",
    "nh3_planar_freq.out",
    "h2o_linear_freq.out",
    "fhh_optts_nofreq.out",
    "fhh_optts_freq.out",
    "n2_stretch_local_spin.out",
    "n2_stretch_casscf_local_spin.out",
    "generated_yb3_sarc2.out",
    "generated_yb3_sarc2_trah.out",
)


def _sections(name: str) -> dict:
    return parse_auto(FIXTURES / name).sections


# --- identification and generic behaviour -----------------------------------


def test_registry_has_orca():
    assert [p.program for p in available()] == ["orca"]


def test_all_fixtures_parse_as_orca():
    for name in ALL_FIXTURES:
        result = parse_auto(FIXTURES / name)
        assert result.program == "orca"
        assert result.sections["version"] == "6.1.1"


def test_detect_rejects_other_files():
    assert detect_program(FIXTURES / "inputs" / "n2_caspt2.inp") is None
    with pytest.raises(ParserError, match="Next step: "):
        parse_auto(FIXTURES / "inputs" / "n2_caspt2.inp")


def test_missing_file_error_has_next_action(tmp_path):
    with pytest.raises(ParserError, match="Next step: "):
        parse_auto(tmp_path / "does_not_exist.out")


def test_garbage_content_rejected(tmp_path):
    bogus = tmp_path / "bogus.out"
    bogus.write_text("this is not the output of any program\n", encoding="utf-8")
    with pytest.raises(ParserError, match="Next step: "):
        OrcaParser().parse(bogus)


# --- SCF --------------------------------------------------------------------


def test_scf_block_of_completed_dft():
    scf = _sections("fblock_dft_la_complex.out")["scf"]
    assert scf["converged"] is True
    assert scf["cycles"] == 64
    assert len(scf["energies"]) == 64
    assert scf["energy"] == pytest.approx(-1241.3910657529527)


def test_scf_not_converged_path():
    sections = _sections("scf_noconv.out")
    assert sections["scf"]["converged"] is False
    assert sections["scf"]["cycles"] == 2
    assert sections["terminated_normally"] is False
    assert sections["final_energy"] is None
    assert sections["errors"]


def test_casscf_only_run_has_no_scf_section():
    """An N2 job with HFTyp=CASSCF has no separate SCF section in the output --
    None is the correct value, not a parse gap."""
    scf = _sections("n2_casscf_nevpt2.out")["scf"]
    assert scf["converged"] is None
    assert scf["cycles"] is None
    assert scf["energies"] == ()


def test_scf_history_ignores_density_matrix_rows():
    """Regression: the SCF iteration-row fingerprint (third column in scientific
    notation) must not match numeric tables such as density-matrix rows."""
    scf = _sections("n2_casscf_nevpt2.out")["scf"]
    assert scf["energies"] == ()


def test_scf_tables_are_parsed():
    """The DIIS iteration table and the SCF convergence summary are parser output since
    v0.2 (the diagnosis layer's SCF triage works on the parse result alone)."""
    scf = _sections("n2_hf_clean.out")["scf"]
    assert scf["solver_seen"] is True
    assert scf["diis_rows"] == ((1, 0.0941), (2, 0.0682), (3, 0.0471), (4, 0.0336))
    assert scf["diis_switch_cycle"] == 2
    assert scf["diis_error_at_switch"] == pytest.approx(0.0682)
    assert scf["converger_switches"] == ("SOSCF",)
    assert len(scf["criteria"]) == 6
    names = [name for name, _value, _tol in scf["criteria"]]
    assert "Energy change" in names and "Orbital Gradient" in names
    assert scf["check_mode"] == 2
    assert "Total+1el-Energy" in scf["check_mode_source"]


def test_scf_tables_absent_without_an_scf_section():
    """A CASSCF-type job prints no SCF tables: the fields stay empty and the mode is the
    documented default, with the source saying that it was assumed."""
    scf = _sections("co_plus_soc.out")["scf"]
    assert scf["solver_seen"] is False
    assert scf["diis_rows"] == ()
    assert scf["criteria"] == ()
    assert scf["check_mode"] == 2
    assert "assumed default" in scf["check_mode_source"]


def test_scf_diis_resets_are_counted():
    scf = _sections("fblock_dft_la_complex.out")["scf"]
    assert scf["diis_resets"] == 1
    assert len(scf["diis_rows"]) == 29


# --- orbital table ----------------------------------------------------------


def test_orbital_table_taken_from_last_block():
    orbitals = _sections("n2_casscf_nevpt2.out")["orbitals"]
    assert len(orbitals["energies"]) == 28
    # Natural-orbital occupations after CASSCF are non-integer (measured:
    # 1.9924 / 1.7092 / 0.2936 / 0.0020)
    assert orbitals["occupations"][4] == pytest.approx(1.9924)
    assert orbitals["occupations"][7] == pytest.approx(0.2936)
    assert orbitals["occupations"][9] == pytest.approx(0.0020)
    assert orbitals["energies"][0] == pytest.approx(-15.657211)


def test_orbital_table_absent_when_run_aborted():
    assert _sections("scf_noconv.out")["orbitals"]["energies"] == ()


# --- CASSCF -----------------------------------------------------------------


def test_casscf_energy_convergence_marker():
    casscf = _sections("n2_casscf_nevpt2.out")["casscf"]
    assert casscf["present"] is True
    assert casscf["converged"] is True
    assert casscf["converged_via"] == "energy"
    assert casscf["macro_iterations"] == 9
    assert casscf["energy"] == pytest.approx(-108.749354438)
    block = casscf["states"][0]
    assert (block["mult"], block["nroots"]) == (1, 2)
    roots = block["roots"]
    assert roots[0]["energy"] == pytest.approx(-108.9436489636)
    assert roots[1]["energy"] == pytest.approx(-108.5550599116)
    assert roots[1]["de_ev"] == pytest.approx(10.574)


def test_casscf_gradient_convergence_marker():
    """Regression: the convergence marker can be the GRADIENT criterion (measured
    in the CO+ SOC job)."""
    casscf = _sections("co_plus_soc.out")["casscf"]
    assert casscf["converged"] is True
    assert casscf["converged_via"] == "gradient"
    assert [(b["mult"], b["nroots"]) for b in casscf["states"]] == [(4, 3), (2, 3)]


def test_casscf_absent_in_pure_dft():
    assert _sections("fblock_dft_la_complex.out")["casscf"]["present"] is False


def test_generated_input_output_roundtrip():
    """A real run of a recipe-layer generated input (closed-loop fixture):
    termination, convergence and active occupations all parse."""
    sections = _sections("generated_ce3_sarc2.out")
    assert sections["terminated_normally"] is True
    assert sections["final_energy"] == pytest.approx(-8851.275311962716)
    casscf = sections["casscf"]
    assert casscf["converged"] is True
    assert casscf["converged_via"] == "gradient"
    assert casscf["active_occupations"] == (1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)


# --- NEVPT2 / CASPT2 --------------------------------------------------------


def test_nevpt2_states_and_energies():
    nevpt2 = _sections("n2_casscf_nevpt2.out")["nevpt2"]
    assert nevpt2["present"] is True
    assert len(nevpt2["states"]) == 2
    root0 = nevpt2["states"][0]
    assert (root0["mult"], root0["root"]) == (1, 0)
    assert root0["de"] == pytest.approx(-0.20632901479979)
    assert root0["e0"] == pytest.approx(-108.94364896355802)
    assert root0["energy"] == pytest.approx(-109.14997797835781)


def test_caspt2_weights_and_denominators():
    caspt2 = _sections("n2_caspt2.out")["caspt2"]
    assert caspt2["present"] is True
    assert len(caspt2["states"]) == 2
    root0, root1 = caspt2["states"]
    assert root0["weight"] == pytest.approx(0.94568724878961)
    assert root1["weight"] == pytest.approx(0.93467343170849)
    # 8 smallest-energy-denominator classes per state
    # (IJAB/ITAB/IJTA/TUAB/IJTU/ITAU/TUVA/ITUV)
    assert len(root0["denominators"]) == 8
    assert root1["denominators"]["ITUV"] == pytest.approx(0.202530004)
    assert caspt2["min_reference_weight"] == pytest.approx(0.93467343170849)
    assert caspt2["min_denominator"] == pytest.approx(0.202530004)


def test_nevpt2_run_has_no_caspt2_section():
    assert _sections("n2_casscf_nevpt2.out")["caspt2"]["present"] is False


# --- SOC markers ------------------------------------------------------------


def test_soc_present_only_in_soc_run():
    expected = {
        "n2_casscf_nevpt2.out": False,
        "n2_caspt2.out": False,
        "co_plus_soc.out": True,
        "co_plus_soc_min.out": True,  # PrintLevel 1: the lower print bound at which SOC is detected
        "scf_noconv.out": False,
        "fblock_dft_la_complex.out": False,
        "fblock_dft_gd_crash.out": False,
    }
    for name, flag in expected.items():
        assert _sections(name)["soc_present"] is flag, name


# --- abnormal termination ---------------------------------------------------


def test_crashed_run_reports_errors():
    sections = _sections("fblock_dft_gd_crash.out")
    assert sections["terminated_normally"] is False
    assert sections["final_energy"] is None
    assert any("Segmentation fault" in e for e in sections["errors"])


# --- coupled-cluster T1 diagnostic (A3 input) -------------------------------


def test_t1_diagnostic_parsed():
    cc = _sections("n2_ccsd.out")["cc"]
    assert cc["present"] is True
    assert cc["t1"] == pytest.approx(0.013031252)
    assert cc["singles_norm"] == pytest.approx(0.041208436)
    assert _sections("f2_ccsd.out")["cc"]["t1"] == pytest.approx(0.011538126)


def test_t1_above_threshold_fixture():
    """Stretched N2 (2.0 A): T1 = 0.0458, above the 0.02 screening line (the
    flagged side of the A3 multi-reference screening)."""
    cc = _sections("n2_stretch_ccsd.out")["cc"]
    assert cc["t1"] == pytest.approx(0.045766132)


def test_cc_absent_in_plain_jobs():
    cc = _sections("n2_hf_clean.out")["cc"]
    assert cc["present"] is False
    assert cc["t1"] is None


# --- vibrational frequencies (D4 input) -------------------------------------


def test_frequency_minimum_has_no_imaginary():
    frequencies = _sections("h2o_freq_min.out")["frequencies"]
    assert frequencies["present"] is True
    block = frequencies["blocks"][-1]
    assert block["n_modes"] == 9  # 6 translations/rotations + 3 real modes
    assert block["n_imaginary"] == 0
    assert block["min_imaginary"] is None
    assert block["modes"][6]["wavenumber"] == pytest.approx(1653.25)
    assert block["modes"][-1]["wavenumber"] == pytest.approx(3932.73)
    assert frequencies["after_geometry"] is True


def test_frequency_single_imaginary_mode():
    block = _sections("nh3_planar_freq.out")["frequencies"]["blocks"][-1]
    assert block["n_modes"] == 12
    assert block["n_imaginary"] == 1
    assert block["min_imaginary"] == pytest.approx(-724.56)
    assert block["imaginary"][0]["index"] == 6


def test_frequency_two_imaginary_modes():
    block = _sections("h2o_linear_freq.out")["frequencies"]["blocks"][-1]
    assert block["n_imaginary"] == 2
    assert block["min_imaginary"] == pytest.approx(-1493.31)


def test_frequency_absent_in_non_frequency_job():
    frequencies = _sections("n2_hf_clean.out")["frequencies"]
    assert frequencies["present"] is False
    assert frequencies["after_geometry"] is None


# --- local spin analysis (A6 input) -----------------------------------------


def test_local_spin_uhf_blocks():
    local_spin = _sections("n2_stretch_local_spin.out")["local_spin"]
    assert local_spin["present"] is True
    assert local_spin["fragment_elements"] == (("N",), ("N",))
    assert len(local_spin["blocks"]) == 3  # initial guess / after SCF / final duplicate
    last = local_spin["blocks"][-1]
    assert last["n_fragments"] == 2
    assert last["sab"][0][0] == pytest.approx(3.6304)
    assert last["sab"][0][1] == pytest.approx(-2.1835)
    assert last["sz"] == (pytest.approx(1.4371), pytest.approx(-1.4371))
    assert last["seff"] == (pytest.approx(1.4699), pytest.approx(1.4699))
    assert last["sz_na"] is False


def test_local_spin_casscf_per_root_blocks():
    local_spin = _sections("n2_stretch_casscf_local_spin.out")["local_spin"]
    assert len(local_spin["blocks"]) == 3  # state average + two roots
    roots = [b for b in local_spin["blocks"] if b["state"] is not None]
    assert [b["state"] for b in roots] == [0, 1]
    assert roots[0]["seff"][0] == pytest.approx(1.4110)
    assert roots[1]["seff"][0] == pytest.approx(0.5003)
    assert roots[0]["sz_na"] is True  # a singlet: <SzA> is n.a. by definition
    average = [b for b in local_spin["blocks"] if b["state"] is None][0]
    assert average["state_label"].startswith("State average")
    assert average["n_fragments"] == 2


def test_local_spin_absent_without_fragments():
    assert _sections("n2_hf_clean.out")["local_spin"]["present"] is False


# --- geometry optimisation / transition state (D4 input) --------------------


def test_optimization_converged_minimum():
    optimization = _sections("h2o_freq_min.out")["optimization"]
    assert optimization["ts"] is False
    assert optimization["has_cycles"] is True
    assert optimization["converged"] is True


def test_optimization_absent_without_cycles():
    optimization = _sections("nh3_planar_freq.out")["optimization"]
    assert optimization["has_cycles"] is False
    assert optimization["converged"] is None


def test_optts_marker_in_output_text():
    optimization = _sections("fhh_optts_nofreq.out")["optimization"]
    assert optimization["ts"] is True
    assert optimization["converged"] is True  # HURRAY ... THE OPTIMIZATION HAS CONVERGED


def test_optts_not_converged_hits_iteration_limit():
    """Planar NH3 OptimTS(+Freq) that reached the iteration limit without
    converging; ORCA also skipped the requested frequency step (measured)."""
    sections = _sections("nh3_planar_optts_freq.out")
    optimization = sections["optimization"]
    assert optimization["ts"] is True
    assert optimization["converged"] is False
    assert sections["frequencies"]["present"] is False


def test_frequency_facts_for_verified_ts():
    facts = facts_from(parse_auto(FIXTURES / "fhh_optts_freq.out"))
    assert facts["ts_optimization"] is True
    assert facts["optimization_converged"] is True
    assert facts["frequency_present"] is True
    assert facts["frequency_after_geometry"] is True
    assert facts["frequency_imaginary_count"] == 1
    assert facts["frequency_min_imaginary"] == pytest.approx(-90.48)


def test_frequency_facts_for_ts_without_verification():
    facts = facts_from(parse_auto(FIXTURES / "fhh_optts_nofreq.out"))
    assert facts["ts_optimization"] is True
    assert facts["frequency_present"] is False
    assert "frequency_imaginary_count" not in facts
    assert "frequency_min_imaginary" not in facts


def test_frequency_facts_for_clean_minimum():
    facts = facts_from(parse_auto(FIXTURES / "h2o_freq_min.out"))
    assert facts["frequency_imaginary_count"] == 0  # 0 is a valid fact, kept
    assert "frequency_min_imaginary" not in facts


# --- solution-branch facts (composition table) -------------------------------


def test_solution_branch_facts_from_the_composition_table():
    """Both Ce fixtures landed on the d1 branch (the singly occupied active orbital is
    100% Ce-d, so the largest f weight over the active orbitals is 0.0%); the N2 fixture
    carries no f-block element, and a run without the composition table carries neither
    fact."""
    ce = facts_from(parse_auto(FIXTURES / "ce3_orbcomp.out"))
    assert ce["f_block_element_present"] is True
    assert ce["active_f_weight_max"] == pytest.approx(0.0)
    generated = facts_from(parse_auto(FIXTURES / "generated_ce3_sarc2.out"))
    assert generated["active_f_weight_max"] == pytest.approx(0.0)
    n2 = facts_from(parse_auto(FIXTURES / "n2_casscf_orbcomp.out"))
    assert n2["f_block_element_present"] is False
    plain = facts_from(parse_auto(FIXTURES / "n2_hf_clean.out"))
    assert "f_block_element_present" not in plain
    assert "active_f_weight_max" not in plain


# --- aborted-CASSCF fixtures (the Yb end-to-end arc) ------------------------


def test_aborted_casscf_fixtures_capture_their_abort_reasons():
    """Measured: an unconverged-wavefunction abort and an out-of-memory abort are both
    reported in the parse result's errors (the error markers are case-insensitive --
    ORCA writes "Aborting the run" here and "aborting the run" in the LEANSCF abort)."""
    default = _sections("generated_yb3_sarc2.out")
    assert default["terminated_normally"] is False
    assert any("IS NOT FULLY CONVERGED" in e for e in default["errors"])
    assert any("Aborting the run" in e for e in default["errors"])
    assert default["casscf"]["present"] is True
    assert default["casscf"]["converged"] is None
    trah = _sections("generated_yb3_sarc2_trah.out")
    assert trah["terminated_normally"] is False
    assert any("OUT OF MEMORY" in e for e in trah["errors"])


# --- fact-field mapping -----------------------------------------------------


def test_facts_for_caspt2_run():
    facts = facts_from(parse_auto(FIXTURES / "n2_caspt2.out"))
    assert facts["caspt2_present"] is True
    assert facts["casscf_converged"] is True
    assert facts["caspt2_min_reference_weight"] == pytest.approx(0.93467343170849)
    # False is a valid fact (kept); None (unknown) does not enter the fact table
    assert facts["soc_present"] is False
    assert "scf_cycles" not in facts  # a CASSCF job has no separate SCF section


def test_facts_drop_unknown_fields():
    facts = facts_from(parse_auto(FIXTURES / "scf_noconv.out"))
    assert facts["scf_cycles"] == 2
    assert facts["scf_converged"] is False
    assert "casscf_converged" not in facts


# --- the Eu3+ fixture's "-0.0000" occupation (regression) --------------------


def test_orbital_table_survives_negative_zero_occupation(tmp_path):
    """ORCA prints tiny negative natural occupations as '-0.0000' (measured on
    the Eu3+ fixture, 2026-09-26).  The row must match, or the table parse stops
    there and every later orbital is lost."""
    text = (
        "                     *   *   *   *\n"
        "                     * O   R   C   A *\n"
        "                     *   *   *   *\n"
        "ORBITAL ENERGIES\n"
        "----------------\n"
        "\n"
        "  NO   OCC          E(Eh)            E(eV) \n"
        "  32   1.0000      -1.329467       -36.1766 \n"
        "  33   -0.0000      -0.495719       -13.4892 \n"
        "  34   0.0000      -2.508227       -68.2523 \n"
        "\n"
        "trailing text\n"
    )
    path = tmp_path / "eu_fragment.out"
    path.write_text(text, encoding="utf-8")
    result = OrcaParser().parse(path)
    occupations = result.sections["orbitals"]["occupations"]
    assert tuple(occupations) == pytest.approx((1.0, -0.0, 0.0))
    energies = result.sections["orbitals"]["energies"]
    assert tuple(energies) == pytest.approx((-1.329467, -0.495719, -2.508227))
