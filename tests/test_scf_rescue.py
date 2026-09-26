"""SCF rescue triage and corrected-input proposals (feature D2).

Real fixtures: ``fixtures/orca/`` (ORCA 6.1.1 outputs, sources in that directory's
README). Three paths have no fixture -- a clean short run, a pseudo-convergence and a
run whose criteria are above tolerance -- so they are built synthetically here; every
synthetic output keeps the layout of the real ones (banner, DIIS table, verdict,
``SCF CONVERGENCE`` block).
"""

from __future__ import annotations

import difflib
import shutil
from pathlib import Path

import pytest

from fblockkit.diagnosis import (
    SEVERITY_LABELS,
    ScfRescueError,
    propose_fixes,
    render_markdown,
    triage,
)
from fblockkit.diagnosis.scf_rescue import (
    FIX_PRESCF,
    FIX_SLOWCONV,
    RULE_ABORTED_NO_VERDICT,
    RULE_CRITERIA_UNMET,
    RULE_DIIS_REBOUND,
    RULE_LONG_CONVERGENCE,
    RULE_NOT_CONVERGED,
    RULE_PSEUDO_CONVERGENCE,
)
from fblockkit.knowledge.models import ParseResult
from fblockkit.parsers import parse_auto

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"
SCF_NOCONV_INP = FIXTURES / "inputs" / "scf_noconv.inp"

TRIAGE_IDS = (
    RULE_NOT_CONVERGED,
    RULE_ABORTED_NO_VERDICT,
    RULE_PSEUDO_CONVERGENCE,
    RULE_LONG_CONVERGENCE,
    RULE_DIIS_REBOUND,
    RULE_CRITERIA_UNMET,
)


def _findings(name: str) -> dict[str, object]:
    return {finding.rule_id: finding for finding in triage(parse_auto(FIXTURES / name))}


def _proposals_for(name: str):
    """Triage a fixture output and propose fixes for its original input."""
    findings = triage(parse_auto(FIXTURES / name))
    return propose_fixes(SCF_NOCONV_INP.read_text(encoding="utf-8"), findings)


# --- synthetic outputs ------------------------------------------------------

_BANNER = (
    "                                 *****************\n"
    "                                 * O   R   C   A *\n"
    "                                 *****************\n"
    "Program Version 6.1.1 - RELEASE\n"
)
_SCF_TABLE_HEADER = (
    "                                      ORCA LEAN-SCF\n"
    "                              memory conserving SCF solver\n"
    "----------------------------------------D-I-I-S-----------------------------------------\n"
    "Iteration    Energy (Eh)           Delta-E    RMSDP     MaxDP     DIISErr   Damp  Time(sec)\n"
    "-----------------------------------------------------------------------------------------\n"
)


def _row(cycle: int, energy: float, delta_e: float, rmsdp: float, maxdp: float, diis: float) -> str:
    return (
        f"{cycle:5d}    {energy:.16f}    {delta_e:+.2e}  {rmsdp:.2e}  {maxdp:.2e}  "
        f"{diis:.2e}  0.700   0.1"
    )


def _criteria_block(rows: list[tuple[str, float, float]]) -> str:
    lines = ["---------------", "SCF CONVERGENCE", "---------------", ""]
    for name, value, tolerance in rows:
        lines.append(
            f"  Last {name:<24} ...    {value:.4e}  Tolerance :   {tolerance:.4e}"
        )
    return "\n".join(lines)


def _echo_block(input_lines: list[str]) -> str:
    """The input echo ORCA writes at the top of an output ('|  7>   MaxIter 3')."""
    return "\n".join(
        f"|{index:3d}> {line}" for index, line in enumerate(input_lines, start=1)
    )


def _output(
    *,
    cycles: int,
    diis_errors: list[float],
    criteria: list[tuple[str, float, float]],
    mode_label: str | None = "Total+1el-Energy",
    echo: list[str] | None = None,
) -> str:
    """A minimal ORCA-like output: banner, input echo, one DIIS table, verdict, criteria."""
    rows = []
    energy = -108.90
    for index, error in enumerate(diis_errors, start=1):
        energy -= 0.01 / index
        rows.append(
            _row(index, energy, -0.01 / index, 1e-3 / index, 1e-2 / index, error)
        )
        if index == 2:
            rows.append("                               ***Turning on AO-DIIS***")
    settings = (
        f" Convergence Check Mode ConvCheckMode   .... {mode_label}\n"
        if mode_label is not None
        else ""
    )
    echo_block = _echo_block(echo) + "\n" if echo else ""
    verdict = f"*           SCF CONVERGED AFTER  {cycles} CYCLES          *"
    return (
        _BANNER
        + echo_block
        + settings
        + _SCF_TABLE_HEADER
        + "\n".join(rows)
        + "\n\n"
        + verdict
        + "\n\n"
        + _criteria_block(criteria)
        + "\n\n                             ****ORCA TERMINATED NORMALLY****\n"
    )


_MET = [
    ("Energy change", 5.0e-10, 1.0e-08),
    ("MAX-Density change", 8.0e-09, 1.0e-07),
    ("RMS-Density change", 3.0e-10, 5.0e-09),
    ("DIIS Error", 5.0e-09, 5.0e-07),
    ("Orbital Gradient", 2.0e-07, 1.0e-05),
    ("Orbital Rotation", 1.0e-06, 1.0e-05),
]

#: A calm DIIS error series: rising start-up cycles, monotone fall after AO-DIIS is
#: switched on (cycle 2), so no rebound can be reported from it.
_CALM = [9.0e-02, 5.0e-02, 1.0e-04, 1.0e-06, 1.0e-08]


def _write(tmp_path: Path, name: str, text: str) -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


# --- triage on the real fixtures --------------------------------------------


def test_not_converged_run_is_triaged_with_its_numbers():
    by_id = _findings("scf_noconv.out")
    finding = by_id[RULE_NOT_CONVERGED]
    assert finding.severity == "error"
    assert "SCF NOT CONVERGED AFTER 2 CYCLES" in finding.message
    assert "4.710e-02" in finding.message  # last DIIS error, read from the table
    assert finding.evidence and all(item.ref for item in finding.evidence)
    assert any("MaxIter" in refusal for refusal in finding.refusals)


def test_not_converged_does_not_fire_pseudo_or_rebound_rules():
    """The DIIS error fell monotonically in scf_noconv.out: nothing else may fire."""
    assert set(_findings("scf_noconv.out")) == {RULE_NOT_CONVERGED}


def test_long_scf_run_is_triaged():
    by_id = _findings("fblock_dft_la_complex.out")
    finding = by_id[RULE_LONG_CONVERGENCE]
    assert finding.severity == "warn"
    assert "64 cycles" in finding.message
    assert "SOSCF" in finding.message  # measured converger switch
    # A rising *start-up* DIIS error (cycle 1: 3.75e-01) is not a rebound: the peak is
    # only taken after AO-DIIS was switched on, where this run falls monotonically.
    assert RULE_DIIS_REBOUND not in by_id


def test_clean_hf_run_stays_silent():
    """Real negative control: an 8-cycle HF run that ORCA accepted.

    Three printed rows sit above their tolerances (MAX-Density 5.8573e-05 vs 1.0000e-05,
    RMS-Density 8.9188e-06 vs 1.0000e-06, DIIS Error 1.8058e-03 vs 1.0000e-06), but the
    default ConvCheckMode=2 enforces the energy change only, and that one is met
    (2.0488e-07 vs 1.0000e-06). Nothing may fire here.
    """
    assert triage(parse_auto(FIXTURES / "n2_hf_clean.out")) == ()


def test_informational_rows_above_tolerance_do_not_fire():
    """Same lesson on the 64-cycle La fixture: its energy change is met (0.75 x tolerance)."""
    by_id = _findings("fblock_dft_la_complex.out")
    assert set(by_id) == {RULE_LONG_CONVERGENCE}


def test_crashed_scf_is_triaged_as_no_verdict_and_as_diis_rebound():
    by_id = _findings("fblock_dft_gd_crash.out")
    verdict = by_id[RULE_ABORTED_NO_VERDICT]
    assert verdict.severity == "error"
    assert "TRAH" in verdict.message  # AutoTRAH had taken over before the crash
    assert "Segmentation fault" in verdict.message
    assert verdict.refusals
    rebound = by_id[RULE_DIIS_REBOUND]
    assert rebound.severity == "warn"
    assert "3.300e-01" in rebound.message and "1.300e+00" in rebound.message
    assert "2.370e+01" in rebound.message  # largest |Delta-E| in the DIIS phase


def test_casscf_only_run_has_no_scf_triage():
    """No SCF table in the output -> nothing for this layer to say (not a silent pass)."""
    assert triage(parse_auto(FIXTURES / "n2_casscf_nevpt2.out")) == ()


def test_triage_findings_render_in_the_report():
    """The triage must be usable through the existing report renderer (no new plumbing)."""
    markdown = render_markdown(
        triage(parse_auto(FIXTURES / "scf_noconv.out")), subject="scf_noconv.out"
    )
    assert RULE_NOT_CONVERGED in markdown
    assert SEVERITY_LABELS["error"] in markdown
    assert "scf_noconv.out" in markdown


def test_every_finding_carries_severity_evidence_and_id():
    for name in (
        "scf_noconv.out",
        "fblock_dft_la_complex.out",
        "fblock_dft_gd_crash.out",
    ):
        for finding in triage(parse_auto(FIXTURES / name)):
            assert finding.rule_id in TRIAGE_IDS
            assert finding.severity in ("info", "warn", "error", "refuse")
            assert finding.evidence
            for item in finding.evidence:
                assert item.text and item.ref
            assert finding.suggested_fix


# --- triage on synthetic outputs --------------------------------------------


def test_clean_short_run_stays_silent(tmp_path):
    text = _output(cycles=12, diis_errors=[9.0e-02, 5.0e-02, 1.0e-04, 1.0e-06, 1.0e-08], criteria=_MET)
    assert triage(parse_auto(_write(tmp_path, "clean.out", text))) == ()


def test_pseudo_convergence_is_triaged(tmp_path):
    text = _output(cycles=2, diis_errors=[9.0e-02, 5.0e-02], criteria=_MET)
    findings = triage(parse_auto(_write(tmp_path, "pseudo.out", text)))
    by_id = {finding.rule_id: finding for finding in findings}
    assert set(by_id) == {RULE_PSEUDO_CONVERGENCE}
    assert by_id[RULE_PSEUDO_CONVERGENCE].severity == "warn"
    assert "-108.9" in by_id[RULE_PSEUDO_CONVERGENCE].message  # last energy read


def test_informational_rows_do_not_fire_under_the_default_mode(tmp_path):
    """Three rows far above tolerance, but the enforced energy change is met -> silent."""
    criteria = list(_MET)
    criteria[1] = ("MAX-Density change", 1.0e-03, 1.0e-05)
    criteria[2] = ("RMS-Density change", 1.0e-03, 1.0e-06)
    criteria[3] = ("DIIS Error", 1.0e-01, 1.0e-06)
    text = _output(cycles=12, diis_errors=_CALM, criteria=criteria)
    assert triage(parse_auto(_write(tmp_path, "informational.out", text))) == ()


def test_unmet_energy_criterion_fires_with_its_numbers(tmp_path):
    criteria = list(_MET)
    criteria[0] = ("Energy change", 2.0e-06, 1.0e-08)  # 200 times above tolerance
    text = _output(cycles=12, diis_errors=_CALM, criteria=criteria)
    findings = triage(parse_auto(_write(tmp_path, "unmet.out", text)))
    assert {finding.rule_id for finding in findings} == {RULE_CRITERIA_UNMET}
    message = findings[0].message
    assert "Energy change 2.000e-06 vs tolerance 1.000e-08" in message
    assert "mode 2" in message and "printed mode label" in message
    assert "informational" in message  # the five other rows are named as not triaged


def test_mode_unknown_in_the_output_falls_back_to_the_default_and_says_so(tmp_path):
    criteria = list(_MET)
    criteria[0] = ("Energy change", 2.0e-06, 1.0e-08)
    text = _output(
        cycles=12, diis_errors=_CALM, criteria=criteria, mode_label=None
    )
    findings = triage(parse_auto(_write(tmp_path, "nomode.out", text)))
    assert {finding.rule_id for finding in findings} == {RULE_CRITERIA_UNMET}
    assert "assumed default" in findings[0].message


def test_unknown_mode_label_falls_back_to_the_default_and_says_so(tmp_path):
    criteria = list(_MET)
    criteria[0] = ("Energy change", 2.0e-06, 1.0e-08)
    text = _output(
        cycles=12, diis_errors=_CALM, criteria=criteria, mode_label="Some-New-Mode"
    )
    findings = triage(parse_auto(_write(tmp_path, "newmode.out", text)))
    assert {finding.rule_id for finding in findings} == {RULE_CRITERIA_UNMET}
    assert "not one we know" in findings[0].message


def test_criteria_margin_ignores_marginal_excursions(tmp_path):
    """A criterion only a hair above its tolerance is last-digit noise, not a finding."""
    criteria = list(_MET)
    criteria[0] = ("Energy change", 1.2e-08, 1.0e-08)  # 1.2 times above
    text = _output(cycles=12, diis_errors=_CALM, criteria=criteria)
    assert triage(parse_auto(_write(tmp_path, "marginal.out", text))) == ()


def test_mode_zero_from_the_echoed_input_checks_every_row(tmp_path):
    criteria = list(_MET)
    criteria[1] = ("MAX-Density change", 1.0e-03, 1.0e-05)  # informational in mode 2
    text = _output(
        cycles=12,
        diis_errors=_CALM,
        criteria=criteria,
        echo=["! BP86 def2-SVP", "%scf", "  ConvCheckMode 0", "end"],
    )
    findings = triage(parse_auto(_write(tmp_path, "mode0.out", text)))
    assert {finding.rule_id for finding in findings} == {RULE_CRITERIA_UNMET}
    message = findings[0].message
    assert "mode 0" in message and "ConvCheckMode 0" in message
    assert "MAX-Density change 1.000e-03 vs tolerance 1.000e-05" in message
    assert "1 of 6" in message


def test_extreme_scf_keyword_also_means_mode_zero(tmp_path):
    """Manual Table 2.10: !ExtremeSCF sets ConvCheckMode 0, i.e. every criterion counts."""
    criteria = list(_MET)
    criteria[1] = ("MAX-Density change", 1.0e-03, 1.0e-05)
    text = _output(
        cycles=12,
        diis_errors=_CALM,
        criteria=criteria,
        echo=["! ExtremeSCF def2-SVP"],
    )
    findings = triage(parse_auto(_write(tmp_path, "extreme.out", text)))
    assert {finding.rule_id for finding in findings} == {RULE_CRITERIA_UNMET}
    assert "ExtremeSCF" in findings[0].message


def test_mode_one_reads_no_row(tmp_path):
    """Mode 1 stops as soon as one criterion is met, so no printed row can be a failure."""
    criteria = list(_MET)
    criteria[1] = ("MAX-Density change", 1.0e-03, 1.0e-05)
    text = _output(
        cycles=12,
        diis_errors=_CALM,
        criteria=criteria,
        echo=["%scf", "  ConvCheckMode 1", "end"],
    )
    assert triage(parse_auto(_write(tmp_path, "mode1.out", text))) == ()


def test_diis_rebound_needs_the_rise_to_come_after_the_switch(tmp_path):
    """A large error in the start-up cycles (before AO-DIIS) is not a rebound."""
    text = _output(cycles=12, diis_errors=[9.0, 3.0, 1.0e-02, 1.0e-04, 1.0e-06, 1.0e-08], criteria=_MET)
    assert triage(parse_auto(_write(tmp_path, "startup.out", text))) == ()


# --- proposals --------------------------------------------------------------


def test_slowconv_proposal_edits_only_the_simple_input_line():
    proposals = _proposals_for("scf_noconv.out")
    assert [proposal.name for proposal in proposals] == [FIX_SLOWCONV, FIX_PRESCF]
    slowconv = proposals[0]
    assert "! HF def2-SVP SlowConv" in slowconv.content
    # The manual's caution travels with the proposal (wrapped over several comment lines).
    assert "SlowConv may converge to a local minimum solution" in _prose(slowconv.content)
    assert slowconv.evidence and any(item.kind == "manual" for item in slowconv.evidence)
    assert len(slowconv.changes) == 2
    original = SCF_NOCONV_INP.read_text(encoding="utf-8")
    original_payload = _active_lines(original)
    assert _active_lines(slowconv.content) == [
        original_payload[0] + " SlowConv",
        *original_payload[1:],
    ]
    assert not any("MaxIter" in line for line in _changed_payload(original, slowconv.content))


def test_proposals_keep_the_geometry_and_the_original_file(tmp_path):
    before = SCF_NOCONV_INP.read_bytes()
    proposals = _proposals_for("scf_noconv.out")
    assert proposals
    for proposal in proposals:
        active = _active_lines(proposal.content)
        for geometry_line in ("* xyz 0 1", "N 0 0 0", "N 0 0 1.094", "*"):
            assert geometry_line in active, proposal.name
    assert SCF_NOCONV_INP.read_bytes() == before


def test_generated_files_stay_ascii_apart_from_the_user_lines():
    """Nothing this layer writes may add typography to an input file (byte-safe files)."""
    original = SCF_NOCONV_INP.read_text(encoding="utf-8")
    user_characters = set(original)
    for proposal in _proposals_for("scf_noconv.out"):
        extra = {
            character for character in proposal.content if ord(character) > 127
        } - user_characters
        assert not extra, (proposal.name, extra)


def test_no_proposal_merely_raises_maxiter():
    original = SCF_NOCONV_INP.read_text(encoding="utf-8")
    proposals = _proposals_for("scf_noconv.out")
    assert proposals
    for proposal in proposals:
        assert "maxiter" not in proposal.name.lower()
        changed = _changed_payload(original, proposal.content)
        assert changed, proposal.name
        assert any("MaxIter" not in line for line in changed), proposal.name


def test_maxiter_refusal_is_attached_to_the_finding():
    finding = _findings("scf_noconv.out")[RULE_NOT_CONVERGED]
    refusals = " ".join(finding.refusals)
    assert "Increasing MaxIter will not help in many cases" in refusals


def test_prescf_proposal_is_the_two_step_route():
    prescf = _proposals_for("scf_noconv.out")[1]
    assert prescf.name == FIX_PRESCF
    assert "! BP86 def2-SVP def2/J RI LooseSCF SlowConv" in prescf.content
    assert '%moinp "prescf.gbw"' in prescf.content
    assert "MaxIter 200" in prescf.content
    # The orbital read-in and the user's own settings are in the commented step 2.
    assert "# ! HF def2-SVP moread" in prescf.content
    assert "# %moinp \"prescf.gbw\"" in prescf.content
    assert "#   MaxIter 3" in prescf.content
    assert prescf.content.count("MaxIter") >= 3
    active_maxiter = [line for line in _active_lines(prescf.content) if "MaxIter" in line]
    assert active_maxiter == ["MaxIter 200"]  # the user's value is only kept as a comment
    assert len(prescf.changes) == 3


def test_prescf_handles_the_one_line_scf_block():
    """The heavy-metal fixtures write '%scf MaxIter 300 end' on one line: keep that form."""
    findings = triage(parse_auto(FIXTURES / "fblock_dft_la_complex.out"))
    text = (FIXTURES / "fblock_dft_la_complex.inp").read_text(encoding="utf-8")
    proposals = propose_fixes(text, findings)
    # A long (but converged) SCF gets the pre-SCF route; the damping keyword is not proposed.
    assert [proposal.name for proposal in proposals] == [FIX_PRESCF]
    prescf = proposals[0]
    assert "# %scf MaxIter 300 end" in prescf.content  # original kept, as a comment
    assert [line for line in _active_lines(prescf.content) if "MaxIter" in line] == ["MaxIter 200"]
    assert "# ! wB97M-V def2-SVPD def2/J TightSCF moread" in prescf.content


def test_pseudo_convergence_gets_the_slowconv_proposal_only(tmp_path):
    text = _output(cycles=2, diis_errors=[9.0e-02, 5.0e-02], criteria=_MET)
    findings = triage(parse_auto(_write(tmp_path, "pseudo.out", text)))
    proposals = propose_fixes(SCF_NOCONV_INP.read_text(encoding="utf-8"), findings)
    assert [proposal.name for proposal in proposals] == [FIX_SLOWCONV]


def test_proposals_are_skipped_when_the_input_already_has_them():
    findings = triage(parse_auto(FIXTURES / "scf_noconv.out"))
    original = SCF_NOCONV_INP.read_text(encoding="utf-8")
    damped = original.replace("! HF def2-SVP", "! HF def2-SVP SlowConv")
    assert [p.name for p in propose_fixes(damped, findings)] == [FIX_PRESCF]
    reading = original.replace("! HF def2-SVP", "! HF def2-SVP moread")
    assert [p.name for p in propose_fixes(reading, findings)] == [FIX_SLOWCONV]
    assert propose_fixes(original, ()) == ()


# --- error paths ------------------------------------------------------------


def test_empty_input_text_is_rejected():
    with pytest.raises(ScfRescueError, match="Next step"):
        propose_fixes("\n\n", ())


def test_input_without_simple_line_is_rejected():
    findings = triage(parse_auto(FIXTURES / "scf_noconv.out"))
    geometry_only = "* xyz 0 1\nN 0 0 0\nN 0 0 1.094\n*\n"  # no '!' line at all
    with pytest.raises(ScfRescueError, match="Next step"):
        propose_fixes(geometry_only, findings)


def test_other_programs_are_refused():
    with pytest.raises(ScfRescueError, match="Next step"):
        triage(ParseResult(program="molcas", path="job.out"))


def test_triage_is_pure_in_memory(tmp_path):
    """v0.2: the DIIS table and the SCF convergence summary are read by the parser, so
    the triage works on the parse result alone -- a result whose file is gone is still
    triaged, with the same findings (it never re-reads the path)."""
    source = FIXTURES / "scf_noconv.out"
    copy = tmp_path / source.name
    shutil.copy(source, copy)
    result = parse_auto(copy)
    findings_before = triage(result)

    copy.unlink()
    findings_after = triage(result)

    assert [f.rule_id for f in findings_after] == [f.rule_id for f in findings_before]
    assert findings_after, "the non-converged fixture must still be triaged"


# --- helpers ----------------------------------------------------------------


def _changed_payload(original: str, proposal: str) -> list[str]:
    """Active lines the proposal adds or replaces, relative to the original input.

    Only the lines ORCA would read are compared, so comment re-wrapping cannot be
    mistaken for an edit of a setting.
    """
    before, after = _active_lines(original), _active_lines(proposal)
    matcher = difflib.SequenceMatcher(a=before, b=after)
    changed: list[str] = []
    for tag, _i1, _i2, j1, j2 in matcher.get_opcodes():
        if tag in ("insert", "replace"):
            changed.extend(after[j1:j2])
    return changed


def _active_lines(content: str) -> list[str]:
    """Lines ORCA would read: everything that is neither blank nor a ``#`` comment."""
    return [
        line.strip()
        for line in content.splitlines()
        if line.strip() and not line.strip().startswith("#")
    ]


def _prose(content: str) -> str:
    """Comment text with the ``#`` markers removed and the wrapping collapsed."""
    text = " ".join(line.strip().lstrip("#").strip() for line in content.splitlines())
    return " ".join(text.split())
