"""Checks of the property-file reader and the CASSCF state data (menu 31).

The fixtures are the ORCA 6.1.1 property files of two state-averaged CASSCF
runs (N2 and H2O, CAS(6,6), 3 singlet roots) together with their outputs: the
per-state energies the reader extracts must reproduce the output's final
``ROOT n:`` lines to the printed precision, and the transition energies' first
two columns must be the eV/cm-1 pair (cross-checked against the output's own
``eV``/``cm**-1`` print).  The availability survey behind menu 31 is recorded
in ``analysis/state_data.py``.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from fblockkit.analysis import state_data
from fblockkit.parsers import parse_auto
from fblockkit.parsers.orca_property import PropertyError, parse_property

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"
N2_PROPERTY = FIXTURES / "n2_sa.property.txt"
H2O_PROPERTY = FIXTURES / "h2o_sa.property.txt"


def _sections(path=N2_PROPERTY):
    return parse_property(path)


def test_the_sections_and_the_types_are_read():
    sections = _sections()
    names = [section.name for section in sections]
    assert "CAS_SCF_Energies" in names and "CASSCF_Absorption_Spectrum" in names
    assert names.count("CAS_SCF_Energies") == 1
    energies = next(s for s in sections if s.name == "CAS_SCF_Energies")
    final = energies.field("finalEnergy")
    assert final.type == "Double" and final.value == pytest.approx(-108.70707440711183)
    # a bare integer carries no &Type annotation (measured: &GeometryIndex 1)
    geometry = energies.field("GeometryIndex")
    assert geometry.type == "Integer" and geometry.value == 1
    method = energies.field("Method")
    assert method.type == "String" and method.value == "CASSCF"


def test_the_arrays_stitch_their_column_groups():
    """A wide array is printed as column groups (8 + 3 here): the reader must
    put the groups back together, not read the second header as a row."""
    sections = _sections()
    spectrum = next(s for s in sections if s.name == "CASSCF_Absorption_Spectrum")
    rows = spectrum.field("ExcitationEnergies").value
    assert len(rows) == 2 and all(len(row) == 11 for row in rows)
    assert spectrum.field("ExcitationEnergies").dim == (2, 11)
    states = spectrum.field("States").value
    assert states == ((0, 0, 1, 0), (0, 0, 2, 0))
    assert spectrum.field("Multiplicities").value == ((1.0, 1.0), (1.0, 1.0))


def test_the_reader_refusals_carry_next_steps(tmp_path):
    broken = tmp_path / "broken.property.txt"
    broken.write_text(
        "$S\n   &x [&Type \"Complex128\"] 1\n$End\n", encoding="utf-8"
    )
    with pytest.raises(PropertyError, match="does not know"):
        parse_property(broken)
    empty = tmp_path / "empty.property.txt"
    empty.write_text("nothing here\n", encoding="utf-8")
    with pytest.raises(PropertyError, match=r"no '\$SECTION'"):
        parse_property(empty)


def test_the_state_table_matches_the_output_roots():
    table = state_data.state_table(_sections())
    assert table["method"] == "CASSCF"
    assert table["final_energy"] == pytest.approx(-108.70707440711183)
    assert (table["n_active_electrons"], table["n_active_orbitals"]) == (6, 6)
    energies = [state.energy_eh for state in table["states"]]
    assert energies == pytest.approx(
        (-108.9756667043, -108.5872679717, -108.5582885454), abs=1e-9
    )
    # the same numbers as the output's final ROOT lines (cross-check)
    parsed = parse_auto(FIXTURES / "n2_sa.out")
    root_energies = tuple(
        root["energy"] for root in parsed.sections["casscf"]["states"][0]["roots"]
    )
    assert tuple(energies) == pytest.approx(root_energies, abs=1e-9)
    assert [state.root for state in table["states"]] == [0, 1, 2]
    assert all(state.mult == 1 for state in table["states"])


def test_the_transitions_and_their_cross_checks():
    table = state_data.state_table(_sections())
    lines = state_data.transitions(_sections())
    assert len(lines) == 2
    first = lines[0]
    assert (first.initial, first.final) == (0, 1)
    assert first.de_ev == pytest.approx(10.568866823824711)
    assert first.de_cm1 == pytest.approx(85243.66864952973)
    # the eV/cm-1 pair satisfies the conversion (the module's own gate)
    assert first.de_ev * state_data.EV_TO_CM1 == pytest.approx(first.de_cm1, rel=2e-5)
    # relative energies of the state table agree with the transitions
    relative = table["states"][1].energy_eh - table["states"][0].energy_eh
    assert relative * 27.211386 == pytest.approx(first.de_ev, rel=1e-6)
    # the measured column layout: wavelength, fosc, D2 and the dipole pairs
    assert first.wavelength_nm == pytest.approx(117.31076522661095)
    assert first.fosc == 0.0 and first.d2_au == 0.0  # g -> g, symmetry-forbidden
    assert first.dipoles is not None and len(first.dipoles) == 3
    assert first.extra == ()
    text = state_data.render(table, lines, source="n2_sa.property.txt")
    assert "fosc" in text
    assert "not persistable" in text


def test_the_h2o_fixture_reads_too():
    sections = _sections(H2O_PROPERTY)
    table = state_data.state_table(sections)
    assert (table["n_active_electrons"], table["n_active_orbitals"]) == (6, 6)
    assert [round(s.energy_eh, 6) for s in table["states"]] == [
        -75.982540, -75.698462, -75.632204
    ]
    lines = state_data.transitions(sections)
    assert len(lines) == 2 and lines[1].de_ev > lines[0].de_ev


def test_missing_sections_refuse_with_a_next_step(tmp_path):
    partial = tmp_path / "partial.property.txt"
    partial.write_text("$Geometry\n   &GeometryIndex 1\n$End\n", encoding="utf-8")
    with pytest.raises(state_data.StateDataError, match="CAS_SCF_Energies"):
        state_data.state_table(parse_property(partial))


# --- the measured absorption columns (nonzero values) ---------


def test_the_nonzero_probe_pins_the_column_layout():
    sections = parse_property(FIXTURES / "h2o_absp.property.txt")
    lines = state_data.transitions(sections)
    assert len(lines) == 3
    first, strong = lines[0], lines[2]
    # wavelength column
    assert first.wavelength_nm == pytest.approx(165.41808013500193)
    # fosc and D2 columns (the strong line carries real intensity)
    assert first.fosc == pytest.approx(0.0112119101246488)
    assert first.d2_au == pytest.approx(0.0610573809153212)
    assert strong.fosc == pytest.approx(0.0776394781480265)
    assert strong.d2_au == pytest.approx(0.3094891937900124)
    # the complex dipole pairs: (re, im) x (DX, DY, DZ); aligned with the
    # output's ABSORPTION SPECTRUM block (DX of the strong line at column 5,
    # DY at 7, DZ at 9)
    assert first.dipoles[2].real == pytest.approx(-0.24709791766690625)
    assert strong.dipoles[0].real == pytest.approx(0.37712944632790035)
    assert strong.dipoles[1].real == pytest.approx(0.4089774743215373)
    # the module's internal gates re-derive f from D2 and D2 from the components
    expected_f = (2.0 / 3.0) * (strong.de_ev / state_data.EV_PER_HARTREE) * strong.d2_au
    assert expected_f == pytest.approx(strong.fosc, rel=1e-6)
    squared = sum(abs(component) ** 2 for component in strong.dipoles)
    assert squared == pytest.approx(strong.d2_au, rel=1e-6)


def test_the_soc_run_reads_its_highest_rel_correction():
    sections = parse_property(FIXTURES / "h2o_absp_soc.property.txt")
    names = [section.name for section in sections]
    # measured: a SOC run prints the absorption section twice
    assert names.count("CASSCF_Absorption_Spectrum") == 2
    lines = state_data.transitions(sections)
    # the chosen section is the SOC-corrected one: its &States irreps are -1
    # (the measured marker of the post-SOC copy) and its density sidecar name
    # is the QDSOC one
    assert all(line.initial_irrep == -1 for line in lines)
    chosen = state_data._section(sections, "CASSCF_Absorption_Spectrum")
    assert chosen.field("RelCorrection").value == 2
    assert chosen.field("Density_name").value == "Tdens-CASQDSOC"
