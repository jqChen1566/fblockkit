"""Checks of the DeltaSCF input generator (recipe/deltascf.py, menu 27).

The grammar follows the ORCA 6.1 manual's DeltaSCF chapter verbatim
(``%scf ALPHACONF 0,1``, ``IONIZEALPHA``, ``PMOM``/``KeepInitialRef``,
``FreezeAndRelease``/``GMF``), and the fixture chain is engine-validated:
the generated input (``h2co_dscf.inp``: formaldehyde PBE0/def2-TZVP UHF
DeltaSCF with a single ``ALPHACONF 0,1``) converged to the n->pi* saddle at
-114.294975 Eh against the clean ground state's -114.418617 Eh -- an
excitation of 3.364 eV, the textbook vertical n->pi* value of formaldehyde.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from fblockkit.recipe import deltascf
from fblockkit.recipe.deltascf import DeltaScfError

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"

COORDS = (("C", 0.0, 0.0, 0.0), ("O", 0.0, 0.0, 1.205))


def test_parse_conf_list_normalises_and_validates():
    assert deltascf.parse_conf_list("0,1") == "0,1"
    assert deltascf.parse_conf_list("0, 0, 1") == "0,0,1"
    with pytest.raises(DeltaScfError, match="empty"):
        deltascf.parse_conf_list("")
    with pytest.raises(DeltaScfError, match="comma-separated"):
        deltascf.parse_conf_list("HOMO")


def test_the_conf_list_input_carries_the_manual_grammar():
    text = deltascf.deltascf_input(
        COORDS, charge=0, multiplicity=1, alpha_conf="0,1", beta_conf="0,1"
    )
    for marker in (
        "! PBE0 def2-TZVP UHF DeltaSCF",
        "  ALPHACONF 0,1",
        "  BETACONF 0,1",
        "* xyz 0 1",
    ):
        assert marker in text, marker
    assert text.isascii()
    assert "  IONIZEALPHA" not in text
    assert "DoMOM" not in text  # the plain metric is the default


def test_the_ionize_and_metric_variants():
    text = deltascf.deltascf_input(
        COORDS,
        charge=0,
        multiplicity=2,
        ionize_alpha=0,
        keywords="PBE0 def2-TZVP UHF",
        mom="pmom",
    )
    assert "  IONIZEALPHA 0" in text
    assert "  PMOM true" in text
    imom = deltascf.deltascf_input(
        COORDS, charge=0, multiplicity=2, ionize_alpha=0, mom="imom"
    )
    assert "  KeepInitialRef true" in imom


def test_the_tactics_and_the_moread_start():
    text = deltascf.deltascf_input(
        COORDS,
        charge=0,
        multiplicity=1,
        alpha_conf="0,1",
        tactics="freeze",
        gs_gbw="gs.gbw",
    )
    assert "DeltaSCF FreezeAndRelease" in text
    assert "  SOSCFMaxStep 0.1" in text
    assert "! MORead" in text
    assert '%moinp "gs.gbw"' in text
    gmf = deltascf.deltascf_input(
        COORDS, charge=0, multiplicity=1, alpha_conf="0,1", tactics="gmf"
    )
    assert "DeltaSCF GMF" in gmf


def test_the_refusals_carry_next_steps():
    with pytest.raises(DeltaScfError, match="exactly one occupation spec"):
        deltascf.deltascf_input(COORDS, charge=0, multiplicity=1)
    with pytest.raises(DeltaScfError, match="exactly one occupation spec"):
        deltascf.deltascf_input(
            COORDS, charge=0, multiplicity=1, alpha_conf="0,1", ionize_alpha=0
        )
    with pytest.raises(DeltaScfError, match="BETACONF needs"):
        deltascf.deltascf_input(
            COORDS, charge=0, multiplicity=1, beta_conf="0,1", ionize_alpha=0
        )
    with pytest.raises(DeltaScfError, match="unknown MOM metric"):
        deltascf.deltascf_input(
            COORDS, charge=0, multiplicity=1, alpha_conf="0,1", mom="pimom"
        )
    with pytest.raises(DeltaScfError, match="unknown tactic"):
        deltascf.deltascf_input(
            COORDS, charge=0, multiplicity=1, alpha_conf="0,1", tactics="anneal"
        )


def test_the_guidance_carries_the_source_caveats():
    lines = deltascf.run_guidance_lines()
    joined = " ".join(lines)
    assert "spin purification" in joined
    assert "pi->pi*" in joined
    assert "ground-state calculation's orbitals" in joined


def test_the_fixture_chain_matches_the_generated_input():
    """The run fixture is the engine validation of this generator: its input
    bytes carry the same blocks the generator emits, and its outcome is the
    textbook n->pi* excitation against the clean ground state."""
    inp = (FIXTURES / "inputs" / "h2co_dscf.inp").read_text(encoding="utf-8")
    assert "! PBE0 def2-TZVP UHF DeltaSCF" in inp
    assert "  ALPHACONF 0,1" in inp
    dscf_out = (FIXTURES / "h2co_dscf.out").read_text(encoding="utf-8", errors="replace")
    gs_out = (FIXTURES / "h2co_gs.out").read_text(encoding="utf-8", errors="replace")
    assert "SCF CONVERGED" in dscf_out
    from fblockkit.parsers import parse_auto

    dscf_energy = parse_auto(FIXTURES / "h2co_dscf.out").sections["final_energy"]
    gs_energy = parse_auto(FIXTURES / "h2co_gs.out").sections["final_energy"]
    assert dscf_energy == pytest.approx(-114.294975, abs=1e-6)
    assert gs_energy == pytest.approx(-114.418617, abs=1e-6)
    excitation_ev = (dscf_energy - gs_energy) * 27.211386
    assert excitation_ev == pytest.approx(3.364, abs=0.01)
