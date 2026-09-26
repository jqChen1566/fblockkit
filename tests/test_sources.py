"""Reference database tests: BibTeX completeness and citation wiring.

Discipline under test: every `literature` evidence item (rules, basis tables,
tool index, analysis/diagnosis constants) must carry a bibkey that exists in
sources.bib, and report output must print complete citations plus paste-ready
BibTeX - so no user ever has to reconstruct a reference by hand.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from fblockkit.analysis import composition, entropy, geometry
from fblockkit.diagnosis import build_report, diagnose, references_section, to_markdown
from fblockkit.diagnosis.cross_level import _PUCL3_EVIDENCE
from fblockkit.knowledge import sources
from fblockkit.knowledge.loader import RuleError, load_rules
from fblockkit.knowledge.models import EVIDENCE_LITERATURE, Evidence
from fblockkit.parsers import parse_auto
from fblockkit.recipe import load_basis_entries
from fblockkit.recipe.basis_ecp import BasisDataError
from fblockkit.toolindex import ToolIndexError, load_tools

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"


def test_sources_bib_is_complete():
    entries = sources.load_bib()
    assert len(entries) >= 8
    for entry in entries:
        for field in sources.REQUIRED_FIELDS:
            assert entry.field(field), f"{entry.key} missing {field}"
        assert entry.field("doi").startswith("10.")


def test_sources_bib_rejects_incomplete_entry(tmp_path):
    bad = tmp_path / "sources.bib"
    bad.write_text(
        "@article{x1,\n  author = {A, B},\n  title = {T},\n  journal = {J},\n}\n",
        encoding="utf-8",
    )
    with pytest.raises(sources.BibDataError, match="missing required field"):
        sources.load_bib(bad)


def test_unknown_bibkey_reports_next_step():
    with pytest.raises(sources.BibDataError, match="Next step"):
        sources.get("no-such-key")


def _literature_evidence_everywhere():
    from fblockkit.recipe import active_space

    for rule in load_rules():
        yield from ((f"rule {rule.id}", item) for item in rule.evidence)
    for entry in load_basis_entries():
        yield from ((f"basis {entry.id}", item) for item in entry.evidence)
    for record in load_tools():
        yield from ((f"tool {record.id}", item) for item in record.evidence)
    for module in (composition, entropy, geometry):
        yield from ((f"module {module.__name__}", item) for item in module.evidence())
    yield ("cross_level", _PUCL3_EVIDENCE)
    for element, valence in (("Dy", 3), ("U", 4)):
        for suggestion in active_space.suggest_active_space(element, valence):
            yield from (
                (f"active_space {element}{valence}", item) for item in suggestion.evidence
            )
    yield from (
        ("active_space verification", item) for item in active_space.verification_evidence()
    )


def test_every_literature_evidence_resolves_to_bib():
    total = 0
    for where, item in _literature_evidence_everywhere():
        if item.kind != EVIDENCE_LITERATURE:
            continue
        total += 1
        assert item.bibkey, f"{where}: literature evidence without bibkey"
        entry = sources.get(item.bibkey)  # raises if unknown
        assert entry.field("doi") in item.ref, (
            f"{where}: ref text must carry the DOI {entry.field('doi')} (complete citation)"
        )
    assert total >= 8, "expected literature evidence in rules/basis/tools/modules"


def test_rule_with_literature_but_no_bibkey_rejected(tmp_path):
    (tmp_path / "bad.yaml").write_text(
        "- id: X1\n  kind: recipe\n  title: t\n  condition: {all: []}\n  action: a\n"
        "  evidence:\n    - {kind: literature, text: t, ref: r}\n",
        encoding="utf-8",
    )
    with pytest.raises(RuleError, match="bibkey"):
        load_rules(tmp_path)


def test_rule_with_unknown_bibkey_rejected(tmp_path):
    (tmp_path / "bad.yaml").write_text(
        "- id: X1\n  kind: recipe\n  title: t\n  condition: {all: []}\n  action: a\n"
        "  evidence:\n    - {kind: literature, text: t, ref: r, bibkey: nope}\n",
        encoding="utf-8",
    )
    with pytest.raises(RuleError, match="sources.bib"):
        load_rules(tmp_path)


def test_basis_entry_with_literature_needs_bibkey(tmp_path):
    bad = tmp_path / "bad.yaml"
    bad.write_text(
        "- id: X1\n  kind: recommend\n  elements: La-Lu\n  basis: b\n  note: n\n"
        "  evidence:\n    - {kind: literature, text: t, ref: r}\n",
        encoding="utf-8",
    )
    with pytest.raises(BasisDataError, match="bibkey"):
        load_basis_entries(bad)


def test_tool_entry_with_literature_needs_bibkey(tmp_path):
    bad = tmp_path / "tools.yaml"
    bad.write_text(
        "- id: x\n  name: X\n  purpose: p\n  relation: index\n  status: planned\n"
        "  license: MIT\n  source: s\n  evidence:\n    - {kind: literature, text: t, ref: r}\n",
        encoding="utf-8",
    )
    with pytest.raises(ToolIndexError, match="bibkey"):
        load_tools(bad)


def test_report_references_are_complete_and_paste_ready():
    result = parse_auto(FIXTURES / "generated_ce3_sarc2.out")
    evidence = tuple(
        item
        for module in (composition, entropy)
        if module.accepts(result)
        for item in module.evidence()
    )
    markdown = to_markdown(
        build_report(diagnose(result), subject="generated_ce3_sarc2.out", extra_evidence=evidence)
    )
    assert "## References" in markdown
    index = markdown.index("## References")
    assert markdown.index("## Findings") < index  # references come last
    for key in ("stein2016automated", "stein2019autocas"):
        assert f"[{key}]" in markdown
    assert "10.1021/acs.jctc.6b00156" in markdown
    assert "```bibtex" in markdown and "@article{stein2016automated," in markdown


def test_references_section_absent_without_literature():
    assert references_section([Evidence(kind="measured", text="t", ref="r")]) is None
