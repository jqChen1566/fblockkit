"""Checks of the RAS/ORMAS input generator (recipe/ras_ormas.py, menu 28).

The masks follow the ORCA 6.1 manual's model-space sections verbatim
(``RAS(Nel: NRAS1 MaxHoles / NRAS2 / NRAS3 MaxParticles)`` and
``ORMAS(nel: m1 min1 max1, m2 min2 max2, ...)``), and the fixture chain is
engine-validated: on N2/def2-SVP at 1.10 Angstrom the unconstrained ORMAS
mask reproduces the plain CASSCF(6,6) energy to all 12 printed digits
(-108.989034756374), the restricted RAS raises it (-108.985623236851), and
the CI-only %rasci route sits above both (-108.921051085808 with RAS, the
superset ORMAS space at -108.921178474942).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from fblockkit.parsers import parse_auto
from fblockkit.recipe import ras_ormas
from fblockkit.recipe.ras_ormas import RasOrmasError

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"

COORDS = (("N", 0.0, 0.0, 0.0), ("N", 0.0, 0.0, 1.10))


def test_parse_ras_mask_normalises_the_manual_forms():
    assert ras_ormas.parse_ras_mask("6:2 2/2/2 2") == (6, 2, 2, 2, 2, 2)
    assert ras_ormas.parse_ras_mask("11:3 1/5/0 0") == (11, 3, 1, 5, 0, 0)
    # the manual's own RIXS example spells a space before the RAS3 pair
    assert ras_ormas.parse_ras_mask("12:4 1/5/ 0 0") == (12, 4, 1, 5, 0, 0)
    assert ras_ormas.format_ras_mask((6, 2, 2, 2, 2, 2)) == "RAS(6:2 2/2/2 2)"


def test_parse_ras_mask_refusals():
    with pytest.raises(RasOrmasError, match="no ':'"):
        ras_ormas.parse_ras_mask("2 2/2/2 2")
    with pytest.raises(RasOrmasError, match="slash"):
        ras_ormas.parse_ras_mask("6:2 2/2 2")
    with pytest.raises(RasOrmasError, match="field counts"):
        ras_ormas.parse_ras_mask("6:2 2 2/2/2 2")
    with pytest.raises(RasOrmasError, match="non-integer"):
        ras_ormas.parse_ras_mask("6:2 2/x/2 2")
    with pytest.raises(RasOrmasError, match="negative"):
        ras_ormas.parse_ras_mask("-2:2 2/2/2 2")


def test_parse_ormas_mask_takes_commas_and_slashes():
    nel, subspaces = ras_ormas.parse_ormas_mask("6: 2 0 4, 2 0 4, 2 0 4")
    assert nel == 6
    assert subspaces == ((2, 0, 4), (2, 0, 4), (2, 0, 4))
    assert ras_ormas.parse_ormas_mask("6: 2 0 4/2 0 4/2 0 4") == (nel, subspaces)
    assert (
        ras_ormas.format_ormas_mask(nel, subspaces)
        == "ORMAS(6: 2 0 4, 2 0 4, 2 0 4)"
    )


def test_parse_ormas_mask_refusals():
    with pytest.raises(RasOrmasError, match="no ':'"):
        ras_ormas.parse_ormas_mask("2 0 4")
    with pytest.raises(RasOrmasError, match="three integers"):
        ras_ormas.parse_ormas_mask("6: 2 0")
    with pytest.raises(RasOrmasError, match="orbital"):
        ras_ormas.parse_ormas_mask("6: 0 0 0")
    with pytest.raises(RasOrmasError, match="2\\*morb"):
        ras_ormas.parse_ormas_mask("6: 2 0 5")
    with pytest.raises(RasOrmasError, match="minimum above"):
        ras_ormas.parse_ormas_mask("6: 2 3 2")
    # the engine's own wording for the over-full mask
    with pytest.raises(RasOrmasError, match="smaller than the number of active electrons"):
        ras_ormas.parse_ormas_mask("6: 2 0 1, 2 0 1")
    with pytest.raises(RasOrmasError, match="cannot carry"):
        ras_ormas.parse_ormas_mask("1: 2 3 4")
    with pytest.raises(RasOrmasError, match="25"):
        ras_ormas.parse_ormas_mask("6: " + ", ".join(["1 0 2"] * 26))


def test_the_casscf_input_carries_the_manual_grammar():
    text = ras_ormas.ras_ormas_input(
        COORDS, charge=0, multiplicity=1, space="ras", mask="6:2 2/2/2 2"
    )
    for marker in (
        "! RHF def2-SVP",
        "%casscf",
        "  nel 6",
        "  norb 6",
        "  refs",
        "    RAS(6:2 2/2/2 2)",
        "  end",  # the refs sub-block's own end (measured: required)
        "end",
        "* xyz 0 1",
    ):
        assert marker in text, marker
    assert text.isascii()
    # the mask's numbers appear exactly once each in the nel/norb lines
    assert "  nel 8" not in text


def test_the_ormas_input_derives_the_numbers_from_the_mask():
    text = ras_ormas.ras_ormas_input(
        COORDS, charge=0, multiplicity=1, space="ormas",
        mask="6: 2 0 4, 2 0 4, 2 0 4",
    )
    assert "  nel 6" in text
    assert "  norb 6" in text
    assert "    ORMAS(6: 2 0 4, 2 0 4, 2 0 4)" in text


def test_the_rasci_variant_and_its_options():
    text = ras_ormas.ras_ormas_input(
        COORDS, charge=0, multiplicity=1, space="ras", mask="6:2 2/2/2 2",
        route="rasci", cistep="accci", exc_level=0,
        mult="1,3", nroots="5,5",
    )
    for marker in ("%rasci", "  refs", "    RAS(6:2 2/2/2 2)", "  mult 1,3",
                   "  nroots 5,5", "  cistep accci", "  ExcLevel 0"):
        assert marker in text, marker
    assert "%casscf" not in text
    plain = ras_ormas.ras_ormas_input(
        COORDS, charge=0, multiplicity=1, space="ras", mask="6:2 2/2/2 2",
        route="rasci",
    )
    assert "cistep" not in plain
    assert "ExcLevel" not in plain


def test_the_refusals_carry_next_steps():
    with pytest.raises(RasOrmasError, match="unknown partition type"):
        ras_ormas.ras_ormas_input(
            COORDS, charge=0, multiplicity=1, space="gas", mask="6:2 2/2/2 2"
        )
    with pytest.raises(RasOrmasError, match="unknown route"):
        ras_ormas.ras_ormas_input(
            COORDS, charge=0, multiplicity=1, space="ras", mask="6:2 2/2/2 2",
            route="mcscf",
        )
    # the engine's own RAS orbital-sum refusal, pre-empted with its wording
    with pytest.raises(RasOrmasError, match="must sum up to NORB"):
        ras_ormas.ras_ormas_input(
            COORDS, charge=0, multiplicity=1, space="ras", mask="6:2 2/2/2 2",
            norb=8,
        )
    with pytest.raises(RasOrmasError, match="active-electron count"):
        ras_ormas.ras_ormas_input(
            COORDS, charge=0, multiplicity=1, space="ras", mask="6:2 2/2/2 2",
            nel=8,
        )
    # the ORMAS override, turned into a refusal
    with pytest.raises(RasOrmasError, match="overwrites the nel/norb lines"):
        ras_ormas.ras_ormas_input(
            COORDS, charge=0, multiplicity=1, space="ormas",
            mask="6: 2 0 4, 2 0 4, 2 0 4", norb=8,
        )
    with pytest.raises(RasOrmasError, match="overwrites the nel line"):
        ras_ormas.ras_ormas_input(
            COORDS, charge=0, multiplicity=1, space="ormas",
            mask="6: 2 0 4, 2 0 4, 2 0 4", nel=4,
        )
    with pytest.raises(RasOrmasError, match="unknown CIStep"):
        ras_ormas.ras_ormas_input(
            COORDS, charge=0, multiplicity=1, space="ras", mask="6:2 2/2/2 2",
            route="rasci", cistep="davidson",
        )
    with pytest.raises(RasOrmasError, match="ExcLevel"):
        ras_ormas.ras_ormas_input(
            COORDS, charge=0, multiplicity=1, space="ras", mask="6:2 2/2/2 2",
            route="rasci", exc_level=-1,
        )
    with pytest.raises(RasOrmasError, match="nroots"):
        ras_ormas.ras_ormas_input(
            COORDS, charge=0, multiplicity=1, space="ras", mask="6:2 2/2/2 2",
            mult="1,3", nroots="5,5,5",
        )


def test_the_guidance_carries_the_source_facts():
    casscf_lines = " ".join(ras_ormas.run_guidance_lines("casscf"))
    assert "active-active rotation" in casscf_lines
    assert "reproduces the full CAS" in casscf_lines
    assert "overwrites nel/norb" in casscf_lines
    rasci_lines = " ".join(ras_ormas.run_guidance_lines("rasci"))
    assert "standalone CI" in rasci_lines
    assert "frozen core" in rasci_lines


def test_the_fixture_chain_matches_the_generated_inputs():
    """The run fixtures are the engine validation of this generator: their
    committed input bytes are exactly what the generator emits for the same
    arguments, and their energies carry the three anchors."""
    cases = {
        "n2_ras": dict(space="ras", mask="6:2 2/2/2 2"),
        "n2_ormas": dict(space="ormas", mask="6: 2 0 4, 2 0 4, 2 0 4"),
        "n2_rasci": dict(
            space="ras", mask="6:2 2/2/2 2", route="rasci",
            cistep="accci", exc_level=0,
        ),
        "n2_rasci_ormas": dict(
            space="ormas", mask="6: 2 0 4, 2 0 4, 2 0 4", route="rasci"
        ),
    }
    for name, kwargs in cases.items():
        produced = ras_ormas.ras_ormas_input(
            COORDS, charge=0, multiplicity=1, **kwargs
        )
        committed = (FIXTURES / "inputs" / f"{name}.inp").read_text(encoding="utf-8")
        assert produced == committed, name

    def energy(name: str) -> float:
        return parse_auto(FIXTURES / f"{name}.out").sections["final_energy"]

    cas = energy("n2_cas666_ref")
    assert cas == pytest.approx(-108.989034756374, abs=1e-9)
    # the identity anchor: the unconstrained ORMAS mask IS the full CAS
    assert energy("n2_ormas") == cas
    # the restricted RAS raises the MCSCF energy
    assert energy("n2_ras") > cas
    # the CI-only route sits above the orbital-optimized one, and its
    # unrestricted (ORMAS) space reaches below its restricted (RAS) space
    assert energy("n2_rasci") > energy("n2_ras")
    assert energy("n2_rasci_ormas") < energy("n2_rasci")
