"""A3 combined multi-reference character panel: T1 + active-space occupation indicators.

Fixtures: fixtures/orca/*.out (real ORCA 6.1.1 outputs, read-only). The T1 values, the
occupation counts and the entropies pinned here were computed from those files once and
are the panel's regression values; the entropy convention is the A2 module's
(spin-summed occupation n, x = n/2).
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from fblockkit.analysis import entropy, mr_diagnostics
from fblockkit.knowledge import sources
from fblockkit.knowledge.models import EVIDENCE_KINDS, EVIDENCE_LITERATURE, ParseResult
from fblockkit.parsers import parse_auto

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"


def _result(name: str):
    return parse_auto(FIXTURES / name)


def _indicator(result, key: str):
    return next(item for item in mr_diagnostics.collect(result) if item.key == key)


# --- (a) T1 diagnostic: above and below the 0.02 screening line ----------------


@pytest.mark.parametrize(
    ("name", "value", "fired"),
    [
        ("n2_ccsd.out", "0.013031252", False),
        ("f2_ccsd.out", "0.011538126", False),
        ("n2_stretch_ccsd.out", "0.045766132", True),
    ],
)
def test_t1_verdict_on_ccsd_fixtures(name, value, fired):
    result = _result(name)
    assert mr_diagnostics.accepts(result)
    t1 = _indicator(result, "t1")
    assert t1.available and t1.fired is fired
    assert t1.value_text == value
    section = mr_diagnostics.run(result)
    assert section.title == "A3 multi-reference character"
    assert value in section.body
    expected = (
        mr_diagnostics.VERDICT_INDICATED if fired else mr_diagnostics.VERDICT_NOT_INDICATED
    )
    assert f"Combined verdict: {expected}" in section.body
    # the reference line is stated as a screening convention, and the singles norm is
    # not silently promoted to a second indicator
    assert "screening convention" in section.body
    assert "not a hard law" in section.body
    assert "ORCA prints T1 only, not D1" in section.body


def test_t1_and_singles_norm_are_one_indicator():
    """Measured on the fixtures: singles_norm / T1 is the square root of the correlated
    electron count (10 for both N2 jobs, 14 for F2), so the two printed numbers are not
    independent indicators."""
    for name, n_correlated in (
        ("n2_ccsd.out", 10),
        ("f2_ccsd.out", 14),
        ("n2_stretch_ccsd.out", 10),
    ):
        cc = _result(name).sections["cc"]
        assert cc["singles_norm"] / cc["t1"] == pytest.approx(
            math.sqrt(n_correlated), abs=1e-6
        )


# --- (b) occupation-based indicators on the CASSCF fixture ---------------------


def test_occupation_indicators_on_n2_casscf():
    result = _result("n2_casscf_nevpt2.out")
    assert mr_diagnostics.accepts(result)
    occupations = result.sections["casscf"]["active_occupations"]
    assert occupations == pytest.approx(
        (1.99236, 1.70922, 1.70922, 0.29360, 0.29360, 0.00200)
    )
    # pinned measured values (recomputed from the fixture once)
    spectrum = entropy.entropy_bound_spectrum(occupations)
    assert spectrum == pytest.approx(
        (
            0.050161126086058244,
            0.8292570416590108,
            0.8292570416590108,
            0.83423593952425,
            0.83423593952425,
            0.015814510224464173,
        )
    )
    assert mr_diagnostics.max_entropy_bound(occupations) == pytest.approx(
        0.83423593952425
    )
    assert mr_diagnostics.fractional_orbitals(occupations) == pytest.approx(
        (1.70922, 1.70922, 0.29360, 0.29360)
    )
    assert len(mr_diagnostics.fractional_orbitals(occupations)) == 4

    section = mr_diagnostics.run(result)
    assert "4 of 6 in 0.02-1.98" in section.body
    assert "0.834236" in section.body
    assert f"Combined verdict: {mr_diagnostics.VERDICT_INDICATED}" in section.body
    assert "fired: fractional active orbitals (4 of 6 fractional)" in section.body
    # the entropy bound is an exclusion test and is indecisive here: it appears under
    # "not fired", never as a positive signal
    assert "max single-orbital entropy bound (0.834236 > 0.14 (indecisive))" in section.body
    # T1 is missing on a plain CASSCF output and the panel says what it would need
    assert "T1 needs a coupled-cluster (CCSD) output" in section.body
    assert "Panel completeness: 2 of 3 indicators" in section.body


def test_fractional_window_excludes_integer_occupations():
    """The definition of a fractional occupation: inside 0.02-1.98 and not integer-valued
    (0, 1 or 2 -- a single configuration's natural occupations)."""
    assert mr_diagnostics.fractional_orbitals(()) == ()
    assert mr_diagnostics.fractional_orbitals((2.0, 1.0, 0.0)) == ()
    assert mr_diagnostics.fractional_orbitals((1.9999, 0.0000)) == ()  # outside the window
    assert mr_diagnostics.fractional_orbitals((1.5000, 0.5000)) == pytest.approx(
        (1.5000, 0.5000)
    )
    # the tolerance is a numerical guard: a near-integer occupation inside the window
    # counts as integer-valued, an occupation further away does not
    assert mr_diagnostics.fractional_orbitals((1.0005,)) == ()
    assert mr_diagnostics.fractional_orbitals((1.0020,)) == pytest.approx((1.0020,))


# --- (c) the f1 edge case: no fractional orbital, half-filled entropy ----------


def test_f1_edge_case_ce3():
    """Ce3+ 4f1, whole f shell in the window: N(occ) = (1, 0, 0, 0, 0, 0, 0).

    The fractional count is 0 (the single occupation is exactly 1), and max(s_bound)
    is ln 4 = 1.386294 -- the open-shell single-occupation limit of the bound, which is
    indecisive (not a positive signal). The verdict rests on the count: no indicator
    points at multireference character.
    """
    result = _result("generated_ce3_sarc2.out")
    assert mr_diagnostics.accepts(result)
    occupations = result.sections["casscf"]["active_occupations"]
    assert occupations == pytest.approx((1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
    assert mr_diagnostics.fractional_orbitals(occupations) == ()
    spectrum = entropy.entropy_bound_spectrum(occupations)
    assert sum(1 for value in spectrum if value == 0.0) == 6
    assert spectrum[0] == pytest.approx(math.log(4.0))
    assert mr_diagnostics.max_entropy_bound(occupations) == pytest.approx(
        math.log(4.0)
    )

    section = mr_diagnostics.run(result)
    assert "0 of 7 in 0.02-1.98" in section.body
    assert "1.386294" in section.body
    assert "not fired: fractional active orbitals (0 of 7 fractional)" in section.body
    assert "(indecisive)" in section.body
    assert "ln 4" in section.body
    assert f"Combined verdict: {mr_diagnostics.VERDICT_NOT_INDICATED}" in section.body


# --- (d) a plain DFT output carries no indicator ------------------------------


def test_accepts_is_false_for_a_plain_dft_output():
    result = _result("fblock_dft_la_complex.out")
    assert mr_diagnostics.accepts(result) is False
    assert all(not item.available for item in mr_diagnostics.collect(result))


# --- (e) provenance: evidence and bibkeys -------------------------------------


def test_evidence_is_non_empty_and_every_bibkey_resolves():
    items = mr_diagnostics.evidence()
    assert items
    table = sources.index_by_key()
    literature = [item for item in items if item.kind == EVIDENCE_LITERATURE]
    assert literature
    for item in items:
        assert item.kind in EVIDENCE_KINDS
        assert item.text and item.ref
        if item.kind != EVIDENCE_LITERATURE:
            continue
        assert item.bibkey in table, item.bibkey
        # the reference text carries the complete citation (the DOI of the bib entry)
        assert sources.get(item.bibkey).field("doi") in item.ref
    # the thresholds and the statement set of this panel are all covered
    keys = {item.bibkey for item in literature}
    assert {"lee1989diagnostic", "stein2016automated", "wardzala2026multireference"} <= keys
    assert any(item.kind == "manual" for item in items)  # the 0.02-1.98 window
    assert any(item.kind == "measured" for item in items)  # the fixture measurements


def test_literature_evidence_of_the_panel_is_reported_with_its_bibkey():
    """The 0.14 line is used twice; both uses are carried, and the second one cites the
    review that calibrates it on main-group/3d systems only."""
    texts = " ".join(item.text for item in mr_diagnostics.evidence())
    assert "M-diagnostic" in texts
    assert "no lanthanide or actinide example" in texts


# --- (f) cannot be judged / missing indicators --------------------------------


def test_cannot_be_judged_when_no_indicator_is_available():
    result = _result("fblock_dft_la_complex.out")
    section = mr_diagnostics.run(result)
    assert f"Combined verdict: {mr_diagnostics.VERDICT_UNKNOWN}" in section.body
    assert "Refused" in section.body
    assert "Next step: " in section.body
    assert mr_diagnostics.combined_verdict(()) == mr_diagnostics.VERDICT_UNKNOWN

    empty = ParseResult(
        program="orca",
        path="(none)",
        sections={
            "cc": {"present": False, "t1": None, "singles_norm": None},
            "casscf": {"present": False, "active_occupations": ()},
        },
    )
    assert mr_diagnostics.accepts(empty) is False
    assert mr_diagnostics.combined_verdict(mr_diagnostics.collect(empty)) == (
        mr_diagnostics.VERDICT_UNKNOWN
    )


def test_a_partial_panel_names_the_missing_indicator():
    """One available indicator must not stand in silently for the whole panel: the body
    states which indicator is missing and what it needs."""
    ccsd_only = mr_diagnostics.run(_result("n2_ccsd.out"))
    assert "unavailable: fractional active orbitals (needs a CASSCF output with an N(occ)= line)" in ccsd_only.body
    assert "Panel completeness: 1 of 3 indicators" in ccsd_only.body
    assert "the occupation-based indicators need a CASSCF output with an N(occ)= line" in ccsd_only.body

    casscf_only = mr_diagnostics.run(_result("n2_casscf_nevpt2.out"))
    assert "unavailable: T1 diagnostic (needs a CCSD output)" in casscf_only.body
    assert "T1 needs a coupled-cluster (CCSD) output" in casscf_only.body
