"""A6 local spin analysis tests: block selection, per-fragment tables, the <SA*SB> share
formula, the singlet refusal, the spin-purity check and the Delta S_E boundary.

Fixtures: fixtures/orca/n2_stretch_local_spin.out (UHF, three blocks, fragments N(1)/N(2))
and fixtures/orca/n2_stretch_casscf_local_spin.out (CASSCF, a state-average block plus two
root blocks); both are read through the parsers layer, so the format is not re-implemented
here.  The metal/ligand path uses a synthetic ParseResult, because no f-block local spin
output exists in the fixtures yet -- the pinned section format it uses is covered by
tests/test_parsers.py.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from fblockkit.analysis import local_spin
from fblockkit.knowledge import sources
from fblockkit.knowledge.models import (
    EVIDENCE_LITERATURE,
    EVIDENCE_MANUAL,
    EVIDENCE_MEASURED,
    ParseResult,
)
from fblockkit.parsers import parse_auto

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"


def _result(name: str) -> ParseResult:
    return parse_auto(FIXTURES / name)


def _synthetic(
    *,
    sab=((8.75, -0.20, -0.20), (-0.20, 0.75, 0.05), (-0.20, 0.05, 0.75)),
    sz=(2.65, -0.075, -0.075),
    seff=(2.5, 0.5, 0.5),
    multiplicity=None,
    fragment_elements=(("Dy",), ("Cl", "Cl"), ("Cl", "Cl")),
) -> ParseResult:
    """A hand-built local spin section: a Dy(III) centre with two chloride fragments.

    The pinned format is the one the parser produces (see tests/test_parsers.py); the
    numbers are chosen so that the share formula can be checked by hand:
    <S^2> = 8.75 + 0.75 + 0.75 - 4 * 0.20 + 2 * 0.05 = 9.55.
    """
    block = {
        "state": None,
        "state_label": "",
        "block_label": "",
        "multiplicity": multiplicity,
        "n_fragments": len(sab),
        "n_atoms": 7,
        "n_basis_functions": 120,
        "sab": sab,
        "sz": sz,
        "sz_na": all(value is None for value in sz),
        "seff": seff,
    }
    return ParseResult(
        program="orca",
        path="synthetic",
        sections={
            "local_spin": {
                "present": True,
                "blocks": (block,),
                "fragment_elements": fragment_elements,
            }
        },
    )


# --- (a) the UHF fixture: last block only ------------------------------------


def test_accepts_and_selects_the_last_scf_block():
    """The SCF job prints three blocks (initial guess / after SCF / final duplicate); only
    the last one belongs to the converged wavefunction."""
    result = _result("n2_stretch_local_spin.out")
    assert local_spin.accepts(result) is True
    selected = local_spin.relevant_blocks(result)
    assert len(selected) == 1
    label, block = selected[0]
    assert label == "SCF final wavefunction (block 3 of 3)"
    assert local_spin.total_s2(block["sab"]) == pytest.approx(2.8938)


def test_uhf_tables_and_pinned_values():
    result = _result("n2_stretch_local_spin.out")
    section = local_spin.run(result)
    assert section.title == "A6 local spin analysis"

    # the final block's own numbers (fixtures/orca/README.md)
    assert "1.4371" in section.body and "-1.4371" in section.body
    assert "+1.4371" in section.body
    assert "1.4699" in section.body
    assert "3.6304" in section.body and "-2.1835" in section.body
    # both fragments are declared, with their element symbols
    assert "  1: N" in section.body and "  2: N" in section.body
    # the earlier blocks are listed so that the choice of block is visible
    assert "1.5231" in section.body  # initial-guess block, unused
    assert "block 3 of 3" in section.body
    # no f-block element among the declared fragments: tables still printed, with the note
    assert "f-block fragment(s): none" in section.body


def test_uhf_sz_sum_is_zero_for_the_ms_zero_state():
    """Sum of <SzA> = +1.4371 - 1.4371 = 0: the local projections add up to the total M_S."""
    result = _result("n2_stretch_local_spin.out")
    body = local_spin.run(result).body
    assert "Sum of <SzA> = +0.0000" in body
    assert "+1.4371 / -1.4371" in body


def test_uhf_seff_and_share_are_pinned():
    result = _result("n2_stretch_local_spin.out")
    _, block = local_spin.relevant_blocks(result)[0]
    rows = local_spin.fragment_rows(block, (("N",), ("N",)))
    assert [row.number for row in rows] == [1, 2]
    assert rows[0].seff == pytest.approx(1.4699)
    assert rows[0].sz == pytest.approx(1.4371)
    assert rows[1].sz == pytest.approx(-1.4371)
    assert rows[0].elements == ("N",)
    assert rows[0].f_block is False
    # <S^2> = 3.6304 + 3.6304 - 2 * 2.1835 = 2.8938; each fragment carries half of it
    assert local_spin.total_s2(block["sab"]) == pytest.approx(2.8938)
    assert local_spin.spin_share(block["sab"], 0) == pytest.approx(0.5)
    assert local_spin.spin_share(block["sab"], 1) == pytest.approx(0.5)
    body = local_spin.run(result).body
    assert "share(M) = (sum over all fragments B of <S_M*S_B>) / <S^2>_total" in body
    assert "0.5000" in body


def test_uhf_final_block_is_spin_contaminated():
    """The final UHF block: <S^2> = 2.8938 against S(S+1) = 0 of the singlet its printed
    M_S = 0 admits -- spin contamination of the broken-symmetry solution."""
    body = local_spin.run(_result("n2_stretch_local_spin.out")).body
    assert "Total <S^2> = sum over all fragment pairs (A, B) of <SA*SB> = 2.8938" in body
    assert "Spin purity" in body
    assert "a singlet, S = 0" in body
    assert "the deviation is +2.8938" in body
    assert "spin-contaminated" in body
    assert "2.8938 is not S(S+1) of any spin eigenstate" in body
    # the output's own broken-symmetry warning is used, so the wording stays factual
    assert "broken-symmetry" in body


# --- (b) the CASSCF fixture: state average + roots ---------------------------


def test_casscf_blocks_state_average_first_then_roots():
    result = _result("n2_stretch_casscf_local_spin.out")
    selected = local_spin.relevant_blocks(result)
    assert len(selected) == 3
    labels = [label for label, _ in selected]
    assert labels[0].startswith("state average")
    assert labels[1].startswith("root 0")
    assert labels[2].startswith("root 1")
    assert "State average densities" in labels[0]
    assert "multiplicity 1" in labels[0]
    body = local_spin.run(result).body
    assert "state-specific analysis of all roots" in body
    assert "State to be analyzed = 0" in body and "State to be analyzed = 1" in body
    for value in ("1.0252", "1.4110", "0.5003", "3.4021", "0.7506"):
        assert value in body, value


def test_casscf_singlet_prints_na_and_refuses_the_share():
    """Every root is a singlet: <SzA> is "n.a." and <S^2> = 0, so the share is 0/0."""
    result = _result("n2_stretch_casscf_local_spin.out")
    for _, block in local_spin.relevant_blocks(result):
        assert block["sz_na"] is True
        assert local_spin.total_s2(block["sab"]) == pytest.approx(0.0, abs=1e-9)
        assert local_spin.spin_share(block["sab"], 0) is None
    body = local_spin.run(result).body
    assert "n.a. (singlet)" in body
    assert "all <SzA> values are zero by definition" in body
    assert "undefined for a total singlet" in body
    assert "Sum of <SzA>: not available" in body
    # the local-spin magnitudes are reported instead of the refused share
    assert "Local-spin magnitude: Seff = 1.4110" in body
    assert "Local-spin magnitude: Seff = 0.5003" in body


def test_casscf_blocks_are_spin_pure():
    """The CASSCF blocks match S(S+1) = 0 of their declared multiplicity, i.e. they are
    spin-pure -- no contamination line."""
    body = local_spin.run(_result("n2_stretch_casscf_local_spin.out")).body
    assert "against S(S+1) = 0.0000 for the declared multiplicity 1 (S = 0)" in body
    assert "this block is spin-pure" in body
    assert "spin-contaminated" not in body
    assert "spin contamination" not in body


# --- (c) the metal path (synthetic f-block case) -----------------------------


def test_metal_fragment_path_and_share_formula():
    result = _synthetic()
    assert local_spin.accepts(result) is True
    block = result.sections["local_spin"]["blocks"][0]
    sab = block["sab"]
    # <S^2> = 8.75 + 0.75 + 0.75 - 4 * 0.20 + 2 * 0.05 = 9.55
    assert local_spin.total_s2(sab) == pytest.approx(9.55)
    # share(Dy) = (8.75 - 0.20 - 0.20) / 9.55, share(Cl) = (0.75 - 0.20 + 0.05) / 9.55
    assert local_spin.spin_share(sab, 0) == pytest.approx(8.35 / 9.55)
    assert local_spin.spin_share(sab, 1) == pytest.approx(0.60 / 9.55)
    assert local_spin.spin_share(sab, 2) == pytest.approx(0.60 / 9.55)
    # the shares partition the total <S^2> exactly
    assert sum(local_spin.spin_share(sab, index) for index in range(3)) == pytest.approx(1.0)
    # fragment identification through the knowledge layer's element table
    assert local_spin.f_block_fragments((("Dy",), ("Cl", "Cl"), ("Cl", "Cl"))) == (1,)

    body = local_spin.run(result).body
    assert "f-block fragment(s): 1 (Dy)" in body
    assert "metal fragment(s) Seff = 2.5000 (fragment 1, Dy)" in body
    assert "The metal carries the larger local spin (mean 2.5000 vs 0.5000" in body
    assert "share(M) = (sum over all fragments B of <S_M*S_B>) / <S^2>_total" in body
    assert "share(fragment 1, Dy) = 8.3500 / 9.5500 = 0.8743" in body
    assert "(f-block)" in body
    assert "The shares sum to 1.0000" in body
    assert "Sum of <SzA> = +2.5000" in body
    # no multiplicity is declared, but the printed M_S = 5/2 bites
    assert "M_S = +2.5000" in body
    assert "S = 5/2" in body


def test_metal_share_is_undefined_for_a_synthetic_singlet():
    sab = ((0.75, -0.75), (-0.75, 0.75))  # <S^2> = 0.75 - 0.75 - 0.75 + 0.75 = 0
    assert local_spin.total_s2(sab) == pytest.approx(0.0)
    assert local_spin.spin_share(sab, 0) is None
    result = _synthetic(
        sab=sab,
        sz=(0.0, 0.0),
        seff=(0.5, 0.5),
        fragment_elements=(("Dy",), ("Cl", "Cl")),
    )
    body = local_spin.run(result).body
    assert "undefined for a total singlet" in body
    assert "0/0" in body
    assert "Local-spin magnitude: metal fragment(s) Seff = 0.5000" in body


def test_synthetic_missing_fragment_elements_skips_the_metal_reading():
    """Without the echoed atom lines the fragments cannot be named: the tables stay, the
    metal/ligand reading is refused and says why."""
    result = _synthetic(fragment_elements=())
    body = local_spin.run(result).body
    assert "f-block fragment(s): not determined" in body
    assert "Next step: " in body
    assert "share(fragment 1" in body  # the share is still reported


def test_partial_na_is_not_labelled_a_singlet():
    """The "(singlet)" tag describes the block, not one row: a fragment that prints n.a.
    inside a spin-polarised block is not a singlet value."""
    body = local_spin.run(_synthetic(sz=(None, 1.5, -1.5))).body
    assert "n.a. (singlet)" not in body
    assert "n.a." in body
    assert "over 2 of 3 fragments" in body


def test_fragment_list_size_mismatch_is_refused():
    result = _synthetic(fragment_elements=(("Dy",), ("Cl", "Cl")))
    with pytest.raises(local_spin.LocalSpinError, match="Next step: "):
        local_spin.run(result)
    with pytest.raises(local_spin.LocalSpinError, match="Next step: "):
        local_spin.spin_share(((1.0,), (1.0,)), 2)


def test_unknown_element_symbols_are_noted_not_guessed():
    result = _synthetic(fragment_elements=(("Dy",), ("Xx",), ("Cl", "Cl")))
    body = local_spin.run(result).body
    assert "Not classified" in body
    assert "Xx" in body
    assert local_spin.unknown_symbols((("Dy",), ("Xx",))) == ("Xx",)


# --- (d) an output without the block ----------------------------------------


def test_accepts_false_without_the_block():
    result = _result("n2_hf_clean.out")
    assert result.sections["local_spin"]["present"] is False
    assert local_spin.accepts(result) is False
    assert local_spin.relevant_blocks(result) == ()
    with pytest.raises(local_spin.LocalSpinError, match="Next step: "):
        local_spin.run(result)


# --- evidence ----------------------------------------------------------------


def test_evidence_sections_and_bibkeys():
    items = local_spin.evidence()
    assert {item.kind for item in items} == {
        EVIDENCE_MANUAL,
        EVIDENCE_LITERATURE,
        EVIDENCE_MEASURED,
    }
    literature = [item for item in items if item.kind == EVIDENCE_LITERATURE]
    assert {item.bibkey for item in literature} == {
        "ai2025density",
        "herrmann2005comparative",
    }
    for item in literature:
        assert item.bibkey, "literature evidence without bibkey"
        entry = sources.get(item.bibkey)  # raises if the key is not in sources.bib
        assert entry.field("doi") in item.ref, f"{item.bibkey}: ref must carry the DOI"


def test_manual_evidence_quotes_the_section_and_its_sentences():
    manual = [item for item in local_spin.evidence() if item.kind == EVIDENCE_MANUAL]
    assert len(manual) == 1
    item = manual[0]
    assert item.url == (
        "https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/"
        "population.html"
    )
    assert "5.1.10" in item.text
    assert "implemented in the SCF and CASSCF modules" in item.text
    assert "All that is required is to divide the molecule into fragments" in item.text
    assert "if nroots > 1, the printing will contain the state-specific analysis of all roots" in item.text
    assert "for a singlet state all <SzA> values are zero by definition" in item.text
    assert "Seff(A)" in item.text and "<SA*SB>" in item.text


def test_measured_evidence_pins_the_fixture_values():
    measured = [item for item in local_spin.evidence() if item.kind == EVIDENCE_MEASURED]
    assert len(measured) == 1
    text = measured[0].text
    for value in ("-2.1835", "+2.1641", "3.6304", "1.4699", "2.8938", "1.4110", "0.5003"):
        assert value in text, value


# --- the boundary paragraph --------------------------------------------------


def test_body_carries_the_delta_s_e_boundary():
    """The surrogate statement travels with the report: A6 is not Delta S_E, the real
    quantity needs orca_2json + orca_loc, and the scale it stands in for is quoted."""
    for result in (
        _result("n2_stretch_local_spin.out"),
        _result("n2_stretch_casscf_local_spin.out"),
        _synthetic(),
    ):
        body = local_spin.run(result).body
        assert "Boundary of this analysis" in body
        assert "NOT the environment spin polarisation entropy Delta S_E" in body
        assert "orca_2json" in body and "RDM2_aa" in body and "orca_loc" in body
        assert "IAO-IBO" in body
        assert "2.766" in body and "0.007" in body
        assert "8.6" in body and "357.7" in body
        assert "discriminating power is limited" in body
        # Eq. (9) is spelled out for the reader
        assert "Tr[D_E^alpha ln D_E^alpha]" in body
