"""G2 tests: f-electron counting, starting-space suggestions, verification protocol."""

from __future__ import annotations

import pytest

from fblockkit.recipe.active_space import (
    ActiveSpaceError,
    f_electron_count,
    suggest_active_space,
    verification_evidence,
    verification_protocol,
)


def test_f_electron_count_known_ions():
    # standard f^n configurations for well-known ions
    assert f_electron_count("Ce", 3) == 1
    assert f_electron_count("Ce", 4) == 0
    assert f_electron_count("Gd", 3) == 7
    assert f_electron_count("Eu", 2) == 7
    assert f_electron_count("Dy", 3) == 9
    assert f_electron_count("U", 3) == 3
    assert f_electron_count("U", 4) == 2
    assert f_electron_count("Pu", 3) == 5
    assert f_electron_count("Th", 4) == 0


def test_f_electron_count_rejects_non_f_and_bad_valence():
    with pytest.raises(ActiveSpaceError, match="Next step"):
        f_electron_count("Fe", 3)
    with pytest.raises(ActiveSpaceError, match="Next step"):
        f_electron_count("Dy", 5)


def test_minimal_suggestion_matches_documented_precedent():
    minimal, double = suggest_active_space("Dy", 3)
    assert (minimal.n_electrons, minimal.n_orbitals) == (9, 7)
    assert minimal.casscf_line() == "CAS(9e,7o)"
    assert minimal.confidence == "confirmed"
    assert double.confidence == "provisional"
    assert double.n_orbitals == 12


def test_actinide_iv_double_shell_carries_valence_requirement():
    _, double = suggest_active_space("U", 4)
    text = " ".join(double.caveats)
    assert "5f and 6d must both be in the valence space" in text
    assert any(item.bibkey == "lu2025normconserving" for item in double.evidence)


def test_regular_configuration_caveat_present():
    minimal, _ = suggest_active_space("Pu", 3)
    assert any("6d" in caveat for caveat in minimal.caveats)


def test_verification_protocol_is_honest_about_s_change():
    lines = "\n".join(verification_protocol())
    assert "0.02-1.98" in lines
    assert "0.14" in lines
    assert "NOT implemented in v0.1" in lines  # S_change: stated, not faked
    assert any(item.kind == "measured" for item in verification_evidence())
    assert any(item.bibkey for item in verification_evidence())
