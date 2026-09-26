"""Documentation/menu isomorphism check (architecture design §3 L3: section number =
menu path number).

This is one of the acceptance items of step 8: changing the menu requires changing the
documentation and vice versa -- checked mechanically by a test rather than promised by
hand.
"""

from __future__ import annotations

import re
from pathlib import Path

from fblockkit.ui import load_menu

DOCS = Path(__file__).resolve().parents[1] / "docs" / "USER_GUIDE.md"


def test_doc_headings_match_menu_numbers():
    text = DOCS.read_text(encoding="utf-8")
    documented = {match.group(1) for match in re.finditer(r"^## (\d+)\s", text, re.M)}
    menu_numbers = {item.number for item in load_menu()}
    assert documented == menu_numbers, (
        f"documentation and menu are not isomorphic: the documentation has "
        f"{documented - menu_numbers} extra, the menu has {menu_numbers - documented} "
        f"extra (changing the menu requires changing the documentation)"
    )


def test_doc_mentions_each_menu_title_topic():
    """Every menu item's key topic words must appear in the documentation (so a title
    cannot be emptied out)."""
    text = DOCS.read_text(encoding="utf-8")
    for keyword in (
        "Check-up and characterisation",
        "Coordination geometry",
        "Generate an ORCA input",
        "Basis-set / ECP recommendation",
        "Search the tool index",
        "Tool onboarding notes",
        "Cross-level solution consistency",
        "Save this session as a script",
        "SCF rescue",
        "Crystal-field fit",
        "Exact four-state entropy",
        "Orbital-space comparison",
        "AVAS target projection",
    ):
        assert keyword in text, keyword


def test_doc_states_the_script_replay_rule():
    # normalize whitespace so rewrapping the paragraph does not break the check
    text = " ".join(DOCS.read_text(encoding="utf-8").split())
    assert "an empty line means" in text
    assert "fblockkit run" in text
