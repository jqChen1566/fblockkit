"""Checks of the two-step CAS-CI core-excited XAS input writer
(recipe/casci_xas.py; menu 37's third mode).

The byte anchor is the [FeCl4]2- L-edge probe pair (fixtures/rocis/
fecl4_casci_xas.step1.* / .step2.*): the generated inputs ran on ORCA
6.1.1 to normal termination; the step-2 output carries the valence
transitions and the Fe 2p to 3d L-edge states at 719.358864 eV, and the
window mechanism is the measured positional one (42 = (96-12)/2 for the
probe; the manual's own example lands on 87 for Fe(acac)3).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from fblockkit.recipe import casci_xas

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
ROCI = FIXTURES / "rocis"

PROBE_COORDINATES = (
    ("Fe", -0.009575, 0.000087, 0.011550),
    ("Cl", -1.774578, -1.558533, -0.219791),
    ("Cl", -0.681385, 1.874633, 1.289404),
    ("Cl", 0.666770, 0.847828, -2.091305),
    ("Cl", 1.756927, -1.164015, 1.071196),
)


def _probe_plan(**overrides):
    arguments = dict(
        charge=-2,
        multiplicity=5,
        valence_nel=6,
        valence_norb=5,
        step1_nroots=5,
        core_orbitals=[6, 7, 8],
    )
    arguments.update(overrides)
    return casci_xas.build_inputs(PROBE_COORDINATES, **arguments)


def test_the_probe_plan_writes_the_frozen_bytes():
    """Both inputs freeze as fixtures/rocis/fecl4_casci_xas.step*.inp (the
    walk's own parameters: %pal nprocs 48 end, and a step-2 %moinp naming
    the gbw the step-1 run actually produces)."""
    plan = _probe_plan(nprocs=48, gbw_name="fecl4.casci_xas.step1.gbw")
    assert plan.step1_text == (ROCI / "fecl4_casci_xas.step1.inp").read_text(
        encoding="utf-8"
    )
    assert plan.step2_text == (ROCI / "fecl4_casci_xas.step2.inp").read_text(
        encoding="utf-8"
    )
    assert plan.window_start == 42  # (96 - 12)/2
    assert plan.step2_nel == 12 and plan.step2_norb == 8
    assert "%pal nprocs 48 end" in plan.step1_text
    assert "{6,42,90,0,0}" in plan.step2_text
    assert '%moinp "fecl4.casci_xas.step1.gbw"' in plan.step2_text
    assert "FrozenCore FC_NONE" in plan.step2_text
    assert "maxiter 1" in plan.step2_text


def test_the_new_parameter_refusals():
    with pytest.raises(casci_xas.CasciXasError, match="unknown mode"):
        _probe_plan(mode="rixs")
    with pytest.raises(casci_xas.CasciXasError, match="XES mode's index"):
        _probe_plan(xas_mo=42)
    with pytest.raises(casci_xas.CasciXasError, match="not a positive count"):
        _probe_plan(nprocs=0)


def test_the_refusals():
    with pytest.raises(casci_xas.CasciXasError, match="no atoms"):
        casci_xas.build_inputs(
            [], charge=0, multiplicity=1, valence_nel=2, valence_norb=2,
            step1_nroots=1, core_orbitals=[0],
        )
    with pytest.raises(casci_xas.CasciXasError, match="not physical"):
        _probe_plan(valence_nel=11, valence_norb=5)
    with pytest.raises(casci_xas.CasciXasError, match="no core orbital"):
        _probe_plan(core_orbitals=[])
    with pytest.raises(casci_xas.CasciXasError, match="repeats"):
        _probe_plan(core_orbitals=[6, 6, 8])
    with pytest.raises(casci_xas.CasciXasError, match="not be integer"):
        _probe_plan(charge=-1)  # 95 electrons against a 12-electron space
    with pytest.raises(casci_xas.CasciXasError, match="differ in"):
        _probe_plan(step2_mult="5,3", step2_nroots="20,20,20")
    with pytest.raises(casci_xas.CasciXasError, match="non-integer"):
        _probe_plan(step2_nroots="twenty")


def test_the_plan_render_states_the_protocol_and_boundaries():
    body = casci_xas.render_plan(_probe_plan(), step1_path="s1.inp", step2_path="s2.inp")
    assert "SOCABS" in body and "FrozenCore" not in body  # the mapspc mode is named
    assert "785" in body  # the measured peak count
    assert "719.36" in body  # the measured L-edge anchor
    assert "positional" in body
    assert "SC-NEVPT2" in body


def test_the_mapspc_recipe_lines():
    assert casci_xas.mapspc_recipe("run.out").startswith("orca_mapspc run.out SOCABS")
    assert "ABS" in casci_xas.mapspc_recipe("run.out", mode="ABS")


def test_the_frozen_outputs_carry_the_measured_anchors():
    """The engine-side anchors of the frozen pair: the step-1 orbital table
    holds the Fe 2p at indices 6-8 (-27.05 Eh), and the step-2 output holds
    the positional window and the L-edge transitions."""
    step1 = (ROCI / "fecl4_casci_xas.step1.out").read_text(encoding="utf-8")
    assert "   6   2.0000     -27.053959" in step1
    assert "   7   2.0000     -27.053956" in step1
    assert "   8   2.0000     -27.053711" in step1
    step2 = (ROCI / "fecl4_casci_xas.step2.out").read_text(encoding="utf-8")
    assert "Active        42 -   49 (   8 orbitals)" in step2
    assert "0-5A  ->  5-5A   719.358864" in step2
    assert "SOC CORRECTED ABSORPTION SPECTRUM" in step2


def test_the_xes_probe_plan_writes_the_frozen_bytes():
    """The emission-mode plan freezes as fixtures/rocis/fecl4_casci_xes.step2.inp
    (the K-beta core set {0,26,27,28}; the walk's own parameters: %pal nprocs 48
    end and a %moinp naming the step-1 gbw).  The step-1 half equals the CAS-CI
    XAS step-1 fixture byte-for-byte (same valence space), so only the step-2
    half carries a frozen copy."""
    plan = _probe_plan(
        mode="xes",
        core_orbitals=[0, 26, 27, 28],
        nprocs=48,
        gbw_name="fecl4.casci_xes.step1.gbw",
    )
    assert plan.step1_text == (ROCI / "fecl4_casci_xas.step1.inp").read_text(
        encoding="utf-8"
    )
    assert plan.step2_text == (ROCI / "fecl4_casci_xes.step2.inp").read_text(
        encoding="utf-8"
    )


def test_the_frozen_xes_output_carries_the_measured_anchors():
    """The engine-side anchors of the frozen RAS-CI XES output (menu 37 mode
    4, the K-beta protocol: the Fe 1s + 3p rotated into the 41-49 window):
    the K-beta emission lines around 7 keV, with the main 319 -> 0 line at
    7086.9 eV (experimental Fe K-beta1 is 7058 eV; the def2-SVP level
    accounts for the offset, measured 2026-10-03)."""
    step2 = (ROCI / "fecl4_casci_xes.step2.out").read_text(encoding="utf-8")
    assert "Active        41 -   49 (   9 orbitals)" in step2
    assert "295-5.0A  -> 294-5.0A   7025.483764" in step2
    assert "319-5.0A  ->  0-5.0A   7086.908295" in step2
    assert "SOC CORRECTED   EMISSION SPECTRUM" in step2


def test_evidence_states_the_protocol_and_the_measurements():
    text = " ".join(
        entry.ref + " " + entry.text + " " + entry.url for entry in casci_xas.evidence()
    )
    assert "3.13.18" in text
    assert "719.36" in text and "785" in text
    assert "XAS/XASSOC" in text


def test_the_xes_plan_text_shape():
    """The emission mode's step-2 text (K-beta core set {0,26,27,28}; the
    byte anchor freezes as fixtures/rocis/fecl4_casci_xes.step2.inp): the
    RAS hole restriction, the velocity/fosc requests, the XESSOC/XASMOs
    rel block, the 40,40 default roots (the manual's saturated 1000,1000
    recipe grows superlinearly in the QDPT transition-density stage and
    did not finish within 24 h at 290 states; measured 2026-10-03) -- and
    the deliberate absence of rel/DoVelocity, which the 6.1.1 input
    scanner rejects (measured)."""
    plan = _probe_plan(
        mode="xes",
        core_orbitals=[0, 26, 27, 28],
        nprocs=48,
        gbw_name="fecl4.casci_xes.step1.gbw",
    )
    text = plan.step2_text
    assert "refs ras(14:4 1/5/0 0) end" in text
    assert "DoDipoleVelocity true" in text
    assert "DecomposeFosc true" in text
    assert "XESSOC true" in text
    assert "XASMOs 41" in text
    assert "DoDTensor false" in text
    assert "nroots 40,40" in text
    assert "%pal nprocs 48 end" in text
    assert 'moinp "fecl4.casci_xes.step1.gbw"' in text
    assert "FrozenCore FC_NONE" in text
    assert "DoVelocity" not in text
