"""D4 diagnosis rules: frequency-analysis and transition-state delayed-onset checks.

Rule logic is checked on hand-built fact tables (a positive and a negative case per
rule, plus the boundary of the 50 cm^-1 threshold), and the end-to-end behaviour is
checked on the real ORCA 6.1.1 fixtures under fixtures/orca/.
"""

from __future__ import annotations

from pathlib import Path

from fblockkit.diagnosis import diagnose, diagnose_facts
from fblockkit.knowledge.loader import load_rules
from fblockkit.parsers import parse_auto

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"

D4_IDS = {
    "D4-FREQ-NO-IMAGINARY-MODE",
    "D4-TS-ONE-IMAGINARY-MODE",
    "D4-TS-HIGHER-ORDER-SADDLE",
    "D4-FREQ-SMALL-IMAGINARY-MODE",
    "D4-TS-NO-FREQUENCY-VERIFICATION",
    "D4-TS-OPTIMISATION-NOT-CONVERGED",
    "D4-BARRIER-REFERENCE-POINT",
}


def _ids(facts: dict) -> set[str]:
    """Rule ids that fire on a hand-built fact table."""
    return {finding.rule_id for finding in diagnose_facts(facts)}


def _fixture_ids(name: str) -> set[str]:
    return {finding.rule_id for finding in diagnose(parse_auto(FIXTURES / name))}


# --- shape of the D4 rule table ----------------------------------------------


def test_d4_rules_are_registered_with_expected_severities():
    rules = {rule.id: rule for rule in load_rules() if rule.id.startswith("D4-")}
    assert set(rules) == D4_IDS
    assert all(rule.kind == "diagnosis" for rule in rules.values())
    assert rules["D4-FREQ-NO-IMAGINARY-MODE"].severity == "info"
    assert rules["D4-TS-ONE-IMAGINARY-MODE"].severity == "info"
    assert rules["D4-TS-HIGHER-ORDER-SADDLE"].severity == "warn"
    assert rules["D4-FREQ-SMALL-IMAGINARY-MODE"].severity == "warn"
    assert rules["D4-TS-NO-FREQUENCY-VERIFICATION"].severity == "warn"
    assert rules["D4-TS-OPTIMISATION-NOT-CONVERGED"].severity == "warn"
    assert rules["D4-BARRIER-REFERENCE-POINT"].severity == "info"
    # the low-frequency rule is our inference, not a manual statement
    assert rules["D4-FREQ-SMALL-IMAGINARY-MODE"].confidence == "provisional"
    for rule in rules.values():
        assert rule.evidence, rule.id
        for item in rule.evidence:
            assert item.text.strip() and item.ref.strip(), rule.id


def test_ts_no_frequency_rule_carries_the_tutorial_source():
    rules = {rule.id: rule for rule in load_rules()}
    evidence = rules["D4-TS-NO-FREQUENCY-VERIFICATION"].evidence
    manual = [item for item in evidence if item.kind == "manual"]
    assert manual, "the tutorial statement must be recorded as manual evidence"
    assert any("faccts.de/docs/orca/6.1/tutorials/react/tsopt.html" in item.url for item in manual)
    assert any("requesting a prior exact Hessian calculation" in item.text for item in manual)


# --- rule logic on hand-built fact tables ------------------------------------


def test_rule1_fires_only_without_imaginary_modes():
    fire = {"frequency_present": True, "frequency_imaginary_count": 0}
    assert "D4-FREQ-NO-IMAGINARY-MODE" in _ids(fire)
    # one imaginary mode is a saddle, not a minimum
    assert "D4-FREQ-NO-IMAGINARY-MODE" not in _ids(
        {"frequency_present": True, "frequency_imaginary_count": 1}
    )
    # no frequency job at all: nothing to conclude about the geometry
    assert "D4-FREQ-NO-IMAGINARY-MODE" not in _ids({"frequency_present": False})


def test_rule2_fires_only_for_exactly_one_imaginary_mode():
    assert "D4-TS-ONE-IMAGINARY-MODE" in _ids(
        {"frequency_imaginary_count": 1, "frequency_min_imaginary": -724.56}
    )
    assert "D4-TS-ONE-IMAGINARY-MODE" not in _ids(
        {"frequency_imaginary_count": 2, "frequency_min_imaginary": -1493.31}
    )
    assert "D4-TS-ONE-IMAGINARY-MODE" not in _ids({"frequency_imaginary_count": 0})


def test_rule3_fires_from_two_imaginary_modes_up():
    assert "D4-TS-HIGHER-ORDER-SADDLE" in _ids(
        {"frequency_imaginary_count": 2, "frequency_min_imaginary": -1493.31}
    )
    assert "D4-TS-HIGHER-ORDER-SADDLE" in _ids(
        {"frequency_imaginary_count": 3, "frequency_min_imaginary": -840.0}
    )
    assert "D4-TS-HIGHER-ORDER-SADDLE" not in _ids(
        {"frequency_imaginary_count": 1, "frequency_min_imaginary": -724.56}
    )


def test_rule4_boundary_is_at_minus_50_cm1():
    """A mode at exactly -50.0 must not fire; -49.9 (a smaller magnitude) must."""
    assert "D4-FREQ-SMALL-IMAGINARY-MODE" not in _ids(
        {"frequency_imaginary_count": 1, "frequency_min_imaginary": -50.0}
    )
    assert "D4-FREQ-SMALL-IMAGINARY-MODE" in _ids(
        {"frequency_imaginary_count": 1, "frequency_min_imaginary": -49.9}
    )
    # well inside the low-frequency region, and with more than one mode counted
    assert "D4-FREQ-SMALL-IMAGINARY-MODE" in _ids(
        {"frequency_imaginary_count": 2, "frequency_min_imaginary": -30.0}
    )
    # a strong imaginary mode is outside the region this rule talks about
    assert "D4-FREQ-SMALL-IMAGINARY-MODE" not in _ids(
        {"frequency_imaginary_count": 1, "frequency_min_imaginary": -724.56}
    )
    # no imaginary mode at all (and no count): the rule must stay silent
    assert "D4-FREQ-SMALL-IMAGINARY-MODE" not in _ids(
        {"frequency_present": True, "frequency_imaginary_count": 0}
    )
    assert "D4-FREQ-SMALL-IMAGINARY-MODE" not in _ids({})


def test_rule5_fires_without_a_verification_of_the_optimised_geometry():
    assert "D4-TS-NO-FREQUENCY-VERIFICATION" in _ids(
        {"ts_optimization": True, "frequency_present": False}
    )
    # a frequency block that does not describe the optimised geometry counts as none
    assert "D4-TS-NO-FREQUENCY-VERIFICATION" in _ids(
        {"ts_optimization": True, "frequency_present": True, "frequency_after_geometry": False}
    )
    # a verified TS job: the frequency block follows the optimisation
    assert "D4-TS-NO-FREQUENCY-VERIFICATION" not in _ids(
        {
            "ts_optimization": True,
            "frequency_present": True,
            "frequency_after_geometry": True,
            "frequency_imaginary_count": 1,
            "frequency_min_imaginary": -90.48,
        }
    )
    # an ordinary optimisation without frequencies is not a TS job
    assert "D4-TS-NO-FREQUENCY-VERIFICATION" not in _ids(
        {"ts_optimization": False, "frequency_present": False}
    )


def test_rule6_fires_only_for_an_unconverged_ts_optimisation():
    assert "D4-TS-OPTIMISATION-NOT-CONVERGED" in _ids(
        {"ts_optimization": True, "optimization_converged": False}
    )
    assert "D4-TS-OPTIMISATION-NOT-CONVERGED" not in _ids(
        {"ts_optimization": True, "optimization_converged": True}
    )
    # an unconverged ordinary optimisation is reported by the SCF / optimisation
    # statements, not by this TS-specific rule
    assert "D4-TS-OPTIMISATION-NOT-CONVERGED" not in _ids(
        {"ts_optimization": False, "optimization_converged": False}
    )


def test_rule7_fires_when_a_ts_like_species_was_found():
    assert "D4-BARRIER-REFERENCE-POINT" in _ids(
        {"frequency_imaginary_count": 1, "frequency_min_imaginary": -724.56}
    )
    assert "D4-BARRIER-REFERENCE-POINT" not in _ids(
        {"frequency_imaginary_count": 2, "frequency_min_imaginary": -1493.31}
    )
    assert "D4-BARRIER-REFERENCE-POINT" not in _ids({"frequency_imaginary_count": 0})


def test_rule7_records_the_measured_barrier_pair():
    rule = next(rule for rule in load_rules() if rule.id == "D4-BARRIER-REFERENCE-POINT")
    text = " ".join(item.text for item in rule.evidence)
    refs = " ".join(item.ref for item in rule.evidence)
    assert "+5.90" in text and "-22.14" in text
    assert "SCINE 能力评估报告.md" in refs and "4.1" in refs


# --- end to end on the real ORCA fixtures ------------------------------------


def test_water_minimum_end_to_end():
    findings = {f.rule_id: f.severity for f in diagnose(parse_auto(FIXTURES / "h2o_freq_min.out"))}
    assert findings["D4-FREQ-NO-IMAGINARY-MODE"] == "info"
    for absent in (
        "D4-TS-ONE-IMAGINARY-MODE",
        "D4-TS-HIGHER-ORDER-SADDLE",
        "D4-FREQ-SMALL-IMAGINARY-MODE",
        "D4-TS-NO-FREQUENCY-VERIFICATION",
        "D4-TS-OPTIMISATION-NOT-CONVERGED",
        "D4-BARRIER-REFERENCE-POINT",
    ):
        assert absent not in findings


def test_planar_ammonia_one_imaginary_mode_end_to_end():
    findings = {f.rule_id: f.severity for f in diagnose(parse_auto(FIXTURES / "nh3_planar_freq.out"))}
    assert findings["D4-TS-ONE-IMAGINARY-MODE"] == "info"
    # the umbrella mode at -724.56 cm^-1 is far outside the low-frequency region
    assert "D4-FREQ-SMALL-IMAGINARY-MODE" not in findings
    assert "D4-TS-HIGHER-ORDER-SADDLE" not in findings
    # the same species is what the barrier-reference note keys on
    assert findings["D4-BARRIER-REFERENCE-POINT"] == "info"


def test_linear_water_two_imaginary_modes_end_to_end():
    findings = {f.rule_id: f.severity for f in diagnose(parse_auto(FIXTURES / "h2o_linear_freq.out"))}
    assert findings["D4-TS-HIGHER-ORDER-SADDLE"] == "warn"
    assert "D4-TS-ONE-IMAGINARY-MODE" not in findings
    assert "D4-FREQ-SMALL-IMAGINARY-MODE" not in findings


def test_fhh_optts_with_freq_end_to_end():
    findings = {f.rule_id: f.severity for f in diagnose(parse_auto(FIXTURES / "fhh_optts_freq.out"))}
    assert findings["D4-TS-ONE-IMAGINARY-MODE"] == "info"
    assert findings["D4-BARRIER-REFERENCE-POINT"] == "info"
    # the frequency block verifies the optimised geometry -> no missing-verification doubt
    assert "D4-TS-NO-FREQUENCY-VERIFICATION" not in findings
    assert "D4-TS-OPTIMISATION-NOT-CONVERGED" not in findings
    # -90.48 cm^-1 is below the 50 cm^-1 magnitude threshold
    assert "D4-FREQ-SMALL-IMAGINARY-MODE" not in findings


def test_fhh_optts_without_freq_end_to_end():
    findings = {f.rule_id: f.severity for f in diagnose(parse_auto(FIXTURES / "fhh_optts_nofreq.out"))}
    assert findings["D4-TS-NO-FREQUENCY-VERIFICATION"] == "warn"
    # the optimisation itself converged, so the non-convergence rule must stay silent
    assert "D4-TS-OPTIMISATION-NOT-CONVERGED" not in findings
    assert "D4-TS-ONE-IMAGINARY-MODE" not in findings


def test_planar_ammonia_optts_not_converged_end_to_end():
    """The optimisation hit the cycle limit, so ORCA skipped the requested Freq step:
    both the missing verification and the non-convergence are reported."""
    result = parse_auto(FIXTURES / "nh3_planar_optts_freq.out")
    findings = {f.rule_id: f.severity for f in diagnose(result)}
    assert findings["D4-TS-NO-FREQUENCY-VERIFICATION"] == "warn"
    assert findings["D4-TS-OPTIMISATION-NOT-CONVERGED"] == "warn"


def test_clean_hf_run_has_no_d4_findings():
    assert not {rule_id for rule_id in _fixture_ids("n2_hf_clean.out") if rule_id.startswith("D4-")}
