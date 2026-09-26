"""Regression checks of the ASS1ST selection round (analysis/ass1st.py, recipe/ass1st.py).

The fixtures are two N2/def2-SVP CAS(6,6) rounds of ORCA 6.1.1
(``n2_ass1st`` single-root, ``n2_ass1st_sa`` averaged over three singlet
states), a round-2 rerun of the suggested (4e, 4o) space (``n2_ass1st.r2``),
and the printed NOON block inside each ``.out``.  The anchors:

1. the density convention ``D_AO = S D_json S``: the CASSCF reference density
   (``Tdens-CAS...``) is recovered *exactly diagonal* in the exported
   natural-orbital basis with the printed occupations (couplings at 1e-15);
2. ORCA's printed "Natural Orbital Occupation Numbers" block equals the
   whole-space eigenvalues of the recovered NEVPT2 density (print precision);
3. the single-root round: internal block 2.00000, 2.00000, 1.98196, 1.97680,
   external top 0.01566, and the suggestion (6,6) -> (4,4) by reassigning the
   sigma-2p (1.99345) to internal and the sigma*-2p (0.00190) to external;
4. the state-averaged round reproduces the same suggestion, and unequal
   weights move the block spectra (weight normalization);
5. the round-2 rerun is the recorded chain caveat: ORCA's default window put
   the (4,4) space on a *different-shaped* solution (occupations
   1.99226/1.93756/0.06480/0.00537), and the next suggestion would shrink to
   (2,2) -- the measured evidence for the "check the round's active orbitals"
   boundary.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pytest

from fblockkit.analysis import ass1st
from fblockkit.analysis.ass1st import Ass1stError, parse_band
from fblockkit.parsers.orca_json import parse_orca_json
from fblockkit.recipe import ass1st as recipe

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"
ST = FIXTURES / "n2_ass1st.json"
ST_OUT = FIXTURES / "n2_ass1st.out"
SA = FIXTURES / "n2_ass1st_sa.json"
R2 = FIXTURES / "n2_ass1st.r2.json"


def _variant(tmp_path: Path, mutate: Callable[[dict[str, Any]], None]) -> Path:
    document = json.loads(ST.read_text(encoding="utf-8"))
    mutate(document["Molecule"])
    target = tmp_path / "variant.json"
    target.write_text(json.dumps(document), encoding="utf-8")
    return target


def _printed_noons(path: Path) -> list[float]:
    text = path.read_text(encoding="utf-8", errors="replace")
    block = re.search(r"Natural Orbital Occupation Numbers:\n((?:N\[\s*\d+\] =.*\n)+)", text)
    assert block is not None, "the fixture has no printed NOON block"
    return [float(m.group(1)) for m in re.finditer(r"=\s*([\d.]+)", block.group(1))]


# --- the density convention ---------------------------------------------------


def test_the_casscf_reference_density_is_recovered_exactly_diagonal():
    """The measured convention D_AO = S D_json S, pinned by the one density whose
    answer is known: the CASSCF reference density is diagonal in the exported
    natural-orbital basis, with the printed occupations."""
    export = parse_orca_json(ST)
    record = dict(export.densities)
    coefficients = np.array(export.mo_coefficients)
    overlap = np.array(export.overlap)
    density = coefficients @ (overlap @ np.array(record["Tdens-CAS.mult.1.root.0.p"]) @ overlap) @ coefficients.T
    diagonal = np.diag(density)
    assert np.abs(diagonal - np.array(export.mo_occupations)).max() < 1e-9
    off = np.abs(density - np.diag(diagonal)).max()
    assert off < 1e-10


def test_the_nevpt2_density_reproduces_orcas_printed_noons():
    """The whole-space eigenvalues of the recovered Tdens-CASNEV density equal the
    "Natural Orbital Occupation Numbers" ORCA printed (to its print precision) --
    which also pins that the printed block is the whole-space naturalization,
    while this module does the source's block diagonalization itself."""
    export = parse_orca_json(ST)
    record = dict(export.densities)
    coefficients = np.array(export.mo_coefficients)
    overlap = np.array(export.overlap)
    density = coefficients @ (overlap @ np.array(record["Tdens-CASNEV.mult.1.root.0.p"]) @ overlap) @ coefficients.T
    eigenvalues = np.sort(np.linalg.eigvalsh(density))[::-1]
    printed = np.array(_printed_noons(ST_OUT))
    assert eigenvalues.shape == printed.shape
    assert np.abs(eigenvalues - printed).max() < 1e-7


# --- the single-root round ----------------------------------------------------


def test_parse_band_accepts_one_or_two_lines():
    assert parse_band("0.05") == (0.05, 1.95)
    assert parse_band("0.03,1.96") == (0.03, 1.96)
    assert parse_band("0.03 1.96") == (0.03, 1.96)
    with pytest.raises(Ass1stError, match="not a sensible occupation window"):
        parse_band("1.5")
    with pytest.raises(Ass1stError, match="one or two"):
        parse_band("0.1 0.2 0.3")


def test_the_single_root_round_reproduces_the_anchors():
    round_ = ass1st.analyze_round(parse_orca_json(ST))
    assert round_.partition.internal == (0, 1, 2, 3)
    assert round_.partition.active == (4, 5, 6, 7, 8, 9)
    assert (round_.partition.n_electrons, round_.partition.n_orbitals) == (6, 6)
    internal = [orb.occupation for orb in round_.internal_block]
    assert internal == pytest.approx([2.0, 2.0, 1.98196, 1.97680], abs=5e-5)
    assert round_.external_block[0].occupation == pytest.approx(0.01566, abs=5e-5)
    suggestion = round_.suggestion
    assert suggestion.reassign_internal == (4,)  # sigma-2p at 1.99345
    assert suggestion.reassign_external == (9,)  # sigma*-2p at 0.00190
    assert not suggestion.add_internal and not suggestion.add_external
    assert (suggestion.n_electrons, suggestion.n_orbitals) == (4, 4)
    assert not suggestion.self_consistent


def test_the_band_marks_the_drifting_orbitals():
    """The active orbital occupations are printed with the round; the band's lines
    decide the reassignments (here the same suggestion at the 0.03/1.96 pair)."""
    round_ = ass1st.analyze_round(parse_orca_json(ST), band=(0.03, 1.96))
    assert round_.suggestion.reassign_internal == (4,)
    assert round_.suggestion.reassign_external == (9,)
    assert (round_.suggestion.n_electrons, round_.suggestion.n_orbitals) == (4, 4)


def test_the_cycle_warning_fires_on_a_visited_space():
    round_ = ass1st.analyze_round(parse_orca_json(ST), previous_spaces=((4, 4),))
    assert round_.cycle_note is not None
    assert "(4e, 4o)" in round_.cycle_note
    assert ass1st.analyze_round(parse_orca_json(ST)).cycle_note is None


# --- the state-averaged round -------------------------------------------------


def test_the_state_averaged_round_reproduces_the_same_suggestion():
    round_ = ass1st.analyze_round(parse_orca_json(SA))
    assert round_.state_averaged and round_.n_states == 3
    assert round_.weights == pytest.approx((1 / 3, 1 / 3, 1 / 3))
    internal = [orb.occupation for orb in round_.internal_block]
    assert internal == pytest.approx([2.0, 2.0, 1.98195, 1.97374], abs=5e-5)
    assert (round_.suggestion.n_electrons, round_.suggestion.n_orbitals) == (4, 4)


def test_unequal_weights_are_normalized_and_move_the_spectra():
    equal = ass1st.analyze_round(parse_orca_json(SA))
    tilted = ass1st.analyze_round(parse_orca_json(SA), weights=(2.0, 1.0, 1.0))
    assert sum(tilted.weights) == pytest.approx(1.0)
    assert tilted.weights[0] == pytest.approx(0.5)
    assert tilted.internal_block[2].occupation != pytest.approx(
        equal.internal_block[2].occupation, abs=1e-7
    )
    with pytest.raises(Ass1stError, match="one per state"):
        ass1st.analyze_round(parse_orca_json(SA), weights=(0.5, 0.5))


# --- the round-2 rerun (the recorded chain caveat) ----------------------------


def test_the_round_two_rerun_shows_the_chain_caveat():
    """Documented measured behaviour: the suggested (4,4) space, rerun from the
    default guess, converged to a different-shaped solution -- the next
    suggestion would shrink to (2,2).  The analysis reports whatever the round
    actually is; the check that the round's active orbitals are the intended
    set is the caller's (the report prints their CASSCF occupations)."""
    round_ = ass1st.analyze_round(parse_orca_json(R2))
    assert (round_.partition.n_electrons, round_.partition.n_orbitals) == (4, 4)
    active = dict(round_.active_occupations)
    assert active[5] == pytest.approx(1.99226, abs=5e-5)
    assert active[8] == pytest.approx(0.00537, abs=5e-5)
    assert round_.suggestion.reassign_internal == (5,)
    assert round_.suggestion.reassign_external == (8,)
    assert (round_.suggestion.n_electrons, round_.suggestion.n_orbitals) == (2, 2)


# --- refusals -----------------------------------------------------------------


def test_a_missing_densities_block_names_the_request(tmp_path):
    variant = _variant(tmp_path, lambda molecule: molecule.pop("Densities"))
    with pytest.raises(Ass1stError, match="Densities") as excinfo:
        ass1st.analyze_round(parse_orca_json(variant))
    assert "orca_2json" in str(excinfo.value)


def test_a_missing_nevpt2_density_names_the_run_settings(tmp_path):
    def mutate(molecule: dict[str, Any]) -> None:
        for name in [k for k in molecule["Densities"] if "CASNEV" in k]:
            molecule["Densities"][name.replace("CASNEV", "CASRENAMED")] = molecule["Densities"].pop(name)

    variant = _variant(tmp_path, mutate)
    with pytest.raises(Ass1stError, match="FIC-NEVPT2 KeepDens"):
        ass1st.analyze_round(parse_orca_json(variant))


def test_a_missing_overlap_names_the_request(tmp_path):
    variant = _variant(tmp_path, lambda molecule: molecule.pop("S-Matrix"))
    with pytest.raises(Ass1stError, match="S-Matrix"):
        ass1st.analyze_round(parse_orca_json(variant))


# --- the round inputs (recipe) ------------------------------------------------


def test_round_one_input_carries_the_measured_blocks():
    text = recipe.round_one_input(
        (("N", 0.0, 0.0, 0.0), ("N", 0.0, 0.0, 1.094)),
        charge=0,
        multiplicity=1,
        n_electrons=6,
        n_orbitals=6,
        n_states=3,
    )
    for marker in (
        "! RHF def2-SVP TightSCF FIC-NEVPT2 KeepDens",
        "%maxcore 2000",
        "  nel 6",
        "  norb 6",
        "  nroots 3",
        "  PTSettings",
        "    Density Unrelaxed",
        "    NatOrbs true",
        "* xyz 0 1",
        "N      0.0000000000",
    ):
        assert marker in text, marker
    assert text.isascii()
    assert recipe.export_conf().startswith('{ "MOCoefficients": true')
    with pytest.raises(recipe.Ass1stInputError, match="do not fit"):
        recipe.round_one_input((("N", 0, 0, 0),), charge=0, multiplicity=1,
                              n_electrons=5, n_orbitals=2)


def test_the_round_stems_chain():
    assert recipe.round_one_stem("n2") == "n2.r1"
    assert recipe.round_one_stem("n2.r1") == "n2.r1"
    assert recipe.next_round_stem("n2.r1") == "n2.r2"
    assert recipe.next_round_stem("n2_ass1st.r3") == "n2_ass1st.r4"
    assert recipe.next_round_stem("n2") == "n2.r2"


# --- output -------------------------------------------------------------------


def test_the_report_prints_the_tables_and_the_suggestion():
    body = ass1st.render(ass1st.analyze_round(parse_orca_json(ST)))
    assert "ASS1ST selection round" in body
    assert "internal block" in body and "external block" in body
    assert "next space: (4e, 4o) -- from (6e, 6o)" in body
    assert "SC-NEVPT2 first-order density" in body  # the substitution is stated
    assert "whole-space naturalization" in body


def test_evidence_carries_the_two_papers():
    bibkeys = {item.bibkey for item in ass1st.evidence() if item.bibkey}
    assert bibkeys == {"khedkar2019ass1st", "khedkar2020sa"}
