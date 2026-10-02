"""Checks of the ROCIS XES input writer (recipe/rocis_xes.py; menu 37's
writing mode).

The byte anchor is the [FeCl4]2- L-edge probe: the built input runs on
ORCA 6.1.1 to normal termination and prints the emission tables with the
ground-state rows at the core-emission energies (744.2/744.3 eV); the
frozen pair is fixtures/rocis/fecl4_xes.{inp,out}.  The measured boundaries
(root-count sufficiency of the plain channel, the SOC channel's storage,
the KHD printout abort) are pinned in the module and the report.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from fblockkit.analysis import xas
from fblockkit.parsers import parse_auto
from fblockkit.recipe import rocis_xes

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
        xas_element=0,
        nroots=30,
        window=(6, 6, 7, 8, 0, 2000),
        rohf_electrons=4,
    )
    arguments.update(overrides)
    return rocis_xes.build_input(PROBE_COORDINATES, **arguments)


def test_the_probe_plan_writes_the_frozen_bytes():
    """The input frozen as fixtures/rocis/fecl4_xes.inp (engine-validated)."""
    plan = _probe_plan()
    frozen = (ROCI / "fecl4_xes.inp").read_text(encoding="utf-8")
    assert plan.text == frozen
    assert "DoRIXS true" in plan.text and "DoRIXSSOC false" in plan.text
    assert "OrbWin 6,6,7,8,0,2000" in plan.text
    assert plan.text.endswith("*\n")


def test_the_soc_flavour_and_the_no_scf_shape():
    plan = _probe_plan(do_rixssoc=True, rohf_electrons=None, do_elastic=False)
    assert "DoRIXSSOC true" in plan.text and "DoElastic false" in plan.text
    assert "%scf" not in plan.text


def test_the_refusals():
    with pytest.raises(rocis_xes.RocisXesError, match="no atoms"):
        rocis_xes.build_input(
            [], charge=0, multiplicity=1, xas_element=0, nroots=30,
            window=(6, 6, 7, 8, 0, 2000),
        )
    with pytest.raises(rocis_xes.RocisXesError, match="not positive"):
        _probe_plan(multiplicity=0)
    with pytest.raises(rocis_xes.RocisXesError, match="outside the structure"):
        _probe_plan(xas_element=5)
    with pytest.raises(rocis_xes.RocisXesError, match="too small for a spectrum"):
        _probe_plan(nroots=1)
    with pytest.raises(rocis_xes.RocisXesError, match="six integers"):
        _probe_plan(window=(6, 6, 7, 8))
    with pytest.raises(rocis_xes.RocisXesError, match="not an increasing"):
        _probe_plan(window=(7, 6, 7, 8, 0, 2000))
    with pytest.raises(rocis_xes.RocisXesError, match="not positive"):
        _probe_plan(rohf_electrons=0)


def test_the_plan_render_states_the_recipe_and_boundaries():
    body = rocis_xes.render_plan(_probe_plan(), path="fecl4.xes.inp")
    assert "orca_mapspc" in body and "XES.stk" in body
    assert "plain RIXS channel" in body
    assert "36 GB" in body  # the measured SOC-channel storage boundary
    assert "KHD variant" in body and "documented termination" in body
    assert "donor 1, donor 2, acceptor" in body
    soc_body = rocis_xes.render_plan(
        _probe_plan(do_rixssoc=True), path="fecl4.xes.inp"
    )
    assert "XESSOC" in soc_body


def test_the_mapspc_recipe_lines():
    assert rocis_xes.mapspc_recipe("run.out") == (
        "orca_mapspc run.out XES -x0<lo> -x1<hi> -w<fwhm> -eV -n<npoints>"
    )
    assert "XESSOC" in rocis_xes.mapspc_recipe("run.out", mode="XESSOC")


def test_the_reader_still_handles_the_xes_run_output():
    """Menu 37's reading mode on the XES run's own output (the emission
    tables do not disturb the absorption-block parsing)."""
    parsed = parse_auto(ROCI / "fecl4_xes.out").sections["rocis"]
    key, rows = xas.primary_spectrum(parsed)
    nonzero = [row for row in rows if row["fosc"] > 1e-6]
    assert key == "soc_dipole_length"
    assert len(rows) == 2235 and len(nonzero) == 934
    body = xas.render(parsed, source="fecl4_xes.out")
    assert "934 with non-zero fosc" in body


def test_evidence_states_the_route_and_the_measurements():
    text = " ".join(
        entry.ref + " " + entry.text + " " + entry.url for entry in rocis_xes.evidence()
    )
    assert "7.31.4" in text
    assert "56 peaks" in text and "2376 peaks" in text
    assert "aborts in the engine's printout" in text
