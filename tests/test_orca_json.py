"""Tests for the ``orca_2json`` orbital-export reader (parsers/orca_json.py).

The fixtures are two exports of the same N2/def2-SVP CAS(6,6) orbital file:
``n2_fcidump.canonical.json`` (canonical orbitals) and ``n2_fcidump.localized.json``
(the same orbitals after an ``orca_loc`` IAO-IBO localization).  The localized range
is the one printed in ``n2_fcidump_step_c.loc.out`` as "Orbital range for
localization ... 4 to 9", 0-based.  Every asserted number is measured from the
fixtures; the natural occupations are cross-checked against the last ``N(occ)=``
line printed by the CASSCF step in ``n2_fcidump_step_a.out``.

The reader is not part of ``parse_auto`` (an orbital export is auxiliary data, not a
program output); the last test pins that decision down.

The ``n2_apc`` / ``h2_apc`` fixtures (one RHF run each, 2026-09-27) additionally
carry the ``FockMatrix`` family and the windowed ``MO_IAJB`` block the APC
ranking consumes; their section pins the measured shapes (one spin matrix for
RHF, the four-integer inclusive window record, the (i, j, a, b, value) entry
layout) and the refusal paths for unmeasured variants.
"""

from __future__ import annotations

import json
from dataclasses import FrozenInstanceError
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pytest

from fblockkit.parsers import ParserError, available, detect_program, parse_auto
from fblockkit.parsers.orca_json import parse_orca_json

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"
CANONICAL = FIXTURES / "n2_fcidump.canonical.json"
LOCALIZED = FIXTURES / "n2_fcidump.localized.json"
APC_N2 = FIXTURES / "n2_apc.json"
APC_H2 = FIXTURES / "h2_apc.json"
ASS1ST_ST = FIXTURES / "n2_ass1st.json"
ASS1ST_SA = FIXTURES / "n2_ass1st_sa.json"

# orca_loc localized this window only (n2_fcidump_step_c.loc.out: "4 to 9", 0-based)
LOCALIZED_RANGE = range(4, 10)

# Last N(occ)= line of n2_fcidump_step_a.out, in the same order as the MO block
ACTIVE_OCCUPATIONS = (1.99345, 1.93634, 1.93634, 0.06599, 0.06599, 0.00190)


def _variant(tmp_path: Path, mutate: Callable[[dict[str, Any]], None]) -> Path:
    """Write a copy of the canonical export with ``mutate`` applied to its Molecule block."""
    document = json.loads(CANONICAL.read_text(encoding="utf-8"))
    mutate(document["Molecule"])
    target = tmp_path / "variant.json"
    target.write_text(json.dumps(document), encoding="utf-8")
    return target


# --- shape and metadata -----------------------------------------------------


def test_shapes_of_the_export():
    export = parse_orca_json(CANONICAL)
    assert export.n_mo == 28
    assert export.n_ao == 28
    assert len(export.mo_coefficients) == 28
    assert {len(row) for row in export.mo_coefficients} == {28}
    assert len(export.mo_occupations) == 28
    assert len(export.mo_energies) == 28
    assert export.overlap is not None
    assert len(export.overlap) == 28
    assert {len(row) for row in export.overlap} == {28}


def test_metadata_and_atom_labels():
    canonical = parse_orca_json(CANONICAL)
    assert canonical.base_name == "step_a"
    assert canonical.charge == 0
    assert canonical.multiplicity == 1
    assert canonical.hftyp == "CASSCF"
    assert canonical.point_group == "C1"
    assert canonical.atoms == ("N", "N")
    # the localized export comes from the .loc.gbw: same molecule, different base name
    localized = parse_orca_json(LOCALIZED)
    assert localized.base_name == "step_a.loc"
    assert localized.atoms == canonical.atoms


def test_export_is_immutable():
    export = parse_orca_json(CANONICAL)
    with pytest.raises(FrozenInstanceError):
        export.n_mo = 12  # type: ignore[misc]


# --- AO labels --------------------------------------------------------------


def test_ao_labels_are_split_into_center_angular_and_component():
    export = parse_orca_json(CANONICAL)
    assert export.ao_labels is not None
    assert len(export.ao_labels) == export.n_ao
    # def2-SVP on N2: 3 s + 2 p shells per atom, then the d polarization shell
    first = export.ao_labels[:6]
    assert [(label.center, label.element) for label in first] == [(0, "N")] * 6
    assert [(label.shell, label.angular, label.component) for label in first] == [
        (1, "s", ""),
        (2, "s", ""),
        (3, "s", ""),
        (1, "p", "z"),
        (1, "p", "x"),
        (1, "p", "y"),
    ]
    second_atom = [label for label in export.ao_labels if label.center == 1]
    assert len(second_atom) == 14
    assert {label.element for label in second_atom} == {"N"}
    # the raw label is kept verbatim next to the split fields
    assert first[0].raw == "0N   1s"


def test_missing_ao_labels_are_accepted(tmp_path):
    path = _variant(tmp_path, lambda molecule: molecule["MolecularOrbitals"].pop("OrbitalLabels"))
    assert parse_orca_json(path).ao_labels is None


def test_ao_label_count_mismatch_is_rejected(tmp_path):
    def drop_one(molecule: dict[str, Any]) -> None:
        molecule["MolecularOrbitals"]["OrbitalLabels"].pop()

    path = _variant(tmp_path, drop_one)
    with pytest.raises(ParserError, match="do not belong to these coefficients"):
        parse_orca_json(path)


def test_an_unparsable_ao_label_is_rejected(tmp_path):
    def break_one(molecule: dict[str, Any]) -> None:
        molecule["MolecularOrbitals"]["OrbitalLabels"][3] = "N pz"

    path = _variant(tmp_path, break_one)
    with pytest.raises(ParserError, match="does not split into the measured ORCA grammar"):
        parse_orca_json(path)


# --- numbers ----------------------------------------------------------------


def test_natural_occupations_match_the_casscf_printout():
    export = parse_orca_json(CANONICAL)
    assert export.mo_occupations[:4] == (2.0,) * 4
    assert export.mo_occupations[4:10] == pytest.approx(ACTIVE_OCCUPATIONS, abs=1e-4)
    assert export.mo_occupations[10:] == (0.0,) * 18


def test_orbital_energies_are_read_in_hartree():
    """The unit is the export's Molecule.MolecularOrbitals.EnergyUnit ("Eh").  The MO
    order is the ``.gbw`` order and is **not** monotonous in energy (measured: MO 9 is
    0.865795 Eh and MO 10 is 0.723735 Eh), so an orbital index must never be derived
    from the energy ordering."""
    export = parse_orca_json(CANONICAL)
    assert export.mo_energies[0] == pytest.approx(-15.653375, abs=1e-6)
    assert export.mo_energies[1] == pytest.approx(-15.652812, abs=1e-6)
    assert export.mo_energies[9] > export.mo_energies[10]


@pytest.mark.parametrize("path", [CANONICAL, LOCALIZED])
def test_coefficients_are_orthonormal_in_the_ao_metric(path):
    """C S C^T = I (rows are MOs): the strong check that the coefficients and the
    overlap matrix belong to the same basis set (measured 9.6e-15 canonical and
    4.1e-13 localized; the tolerance leaves room for a different BLAS).  It holds after
    the localization as well, which is what makes the localized export usable for a
    basis transformation."""
    export = parse_orca_json(path)
    coefficients = np.array(export.mo_coefficients)
    overlap = np.array(export.overlap)
    assert np.abs(overlap - overlap.T).max() < 1e-12
    deviation = np.abs(coefficients @ overlap @ coefficients.T - np.eye(export.n_ao)).max()
    assert deviation < 1e-8


# --- the canonical / localized pair -----------------------------------------


def test_localization_changes_only_the_localized_range():
    """Two exports of the same orbital file differ exactly inside the localized window;
    outside it the coefficients are identical element by element (measured bit-for-bit,
    the block contains no signed zeros)."""
    canonical = parse_orca_json(CANONICAL)
    localized = parse_orca_json(LOCALIZED)
    for index in range(canonical.n_mo):
        left = canonical.mo_coefficients[index]
        right = localized.mo_coefficients[index]
        if index in LOCALIZED_RANGE:
            assert left != right
            assert max(abs(a - b) for a, b in zip(left, right)) > 0.1
        else:
            assert left == right


def test_localized_export_repeats_the_source_occupations_and_energies():
    """Measured caveat recorded in the module docstring: the Occupancy and OrbitalEnergy
    fields of the localized export are carried over from the source orbitals (all 28
    MOs), so they are not the occupations of the rotated orbitals."""
    canonical = parse_orca_json(CANONICAL)
    localized = parse_orca_json(LOCALIZED)
    assert localized.mo_occupations == canonical.mo_occupations
    assert localized.mo_energies == canonical.mo_energies


# --- optional and malformed exports -----------------------------------------


def test_export_without_overlap_is_accepted(tmp_path):
    """No S-Matrix block is a documented export variant, not a defect: the reader
    returns overlap=None and the rest of the export unchanged."""
    variant = _variant(tmp_path, lambda molecule: molecule.pop("S-Matrix"))
    export = parse_orca_json(variant)
    assert export.overlap is None
    assert export.n_mo == 28 and export.n_ao == 28
    assert export.mo_coefficients == parse_orca_json(CANONICAL).mo_coefficients


@pytest.mark.parametrize("mutate", [
    pytest.param(lambda molecule: molecule["S-Matrix"].pop(), id="one-row-missing"),
    pytest.param(lambda molecule: molecule["S-Matrix"][0].pop(), id="ragged-row"),
], )
def test_overlap_dimension_mismatch_is_rejected(tmp_path, mutate):
    variant = _variant(tmp_path, mutate)
    with pytest.raises(ParserError, match="Next step: ") as excinfo:
        parse_orca_json(variant)
    assert "S-Matrix" in str(excinfo.value)


def test_missing_molecule_block_is_rejected(tmp_path):
    target = tmp_path / "no_molecule.json"
    target.write_text(json.dumps({"ORCA Header": {"Version": "6.1.1"}}), encoding="utf-8")
    with pytest.raises(ParserError, match="Next step: ") as excinfo:
        parse_orca_json(target)
    assert "Molecule" in str(excinfo.value)


def test_missing_orbitals_are_rejected(tmp_path):
    variant = _variant(tmp_path, lambda molecule: molecule.pop("MolecularOrbitals"))
    with pytest.raises(ParserError, match="Next step: ") as excinfo:
        parse_orca_json(variant)
    assert "MolecularOrbitals" in str(excinfo.value)


def test_a_wrong_value_type_is_rejected(tmp_path):
    """The JSON is machine-written: a string where a number belongs means the file is
    not the export the caller thinks it is, so it must not be coerced silently."""
    def mutate(molecule: dict[str, Any]) -> None:
        molecule["MolecularOrbitals"]["MOs"][0]["Occupancy"] = "2.0"

    variant = _variant(tmp_path, mutate)
    with pytest.raises(ParserError, match="Next step: ") as excinfo:
        parse_orca_json(variant)
    assert "Occupancy" in str(excinfo.value)


# --- error paths on the file itself -----------------------------------------


def test_a_byte_order_mark_does_not_break_the_read(tmp_path):
    """Some Windows tools prefix a UTF-8 byte-order mark.  json.loads rejects it, so the
    reader strips it instead of reporting a valid export as malformed JSON."""
    target = tmp_path / "bom.json"
    target.write_bytes(CANONICAL.read_bytes().decode("utf-8").encode("utf-8-sig"))
    export = parse_orca_json(target)
    assert export.n_mo == 28
    assert export.base_name == "step_a"


def test_missing_file_error_has_next_action(tmp_path):
    with pytest.raises(ParserError, match="Next step: "):
        parse_orca_json(tmp_path / "does_not_exist.json")


def test_an_orca_output_is_not_a_json_export():
    with pytest.raises(ParserError, match="Next step: ") as excinfo:
        parse_orca_json(FIXTURES / "n2_fcidump_step_a.out")
    assert "not valid JSON" in str(excinfo.value)


# --- the Fock family and the MO_IAJB window (menu 21) ------------------------


def _apc_variant(tmp_path: Path, mutate: Callable[[dict[str, Any]], None]) -> Path:
    """Write a copy of the n2_apc export with ``mutate`` applied to its Molecule block."""
    document = json.loads(APC_N2.read_text(encoding="utf-8"))
    mutate(document["Molecule"])
    target = tmp_path / "apc_variant.json"
    target.write_text(json.dumps(document), encoding="utf-8")
    return target


def test_the_fock_family_blocks_are_read():
    """The density-derived blocks carry one spin dimension (manual: N_spin x N_AO x
    N_AO); the n2_apc request asked for J and K only, the h2_apc request adds F."""
    n2 = parse_orca_json(APC_N2)
    assert n2.hamiltonian is not None and len(n2.hamiltonian) == n2.n_ao
    assert n2.coulomb is not None and len(n2.coulomb) == 1
    assert n2.exchange is not None and len(n2.exchange) == 1
    assert len(n2.exchange[0]) == n2.n_ao
    assert n2.coulomb_exchange is None
    h2 = parse_orca_json(APC_H2)
    assert h2.coulomb_exchange is not None and len(h2.coulomb_exchange) == 1


def test_the_mo_iajb_entries_and_window_are_read():
    """Entries are (i, j, a, b, value) with the internal indices first and value the
    chemist-notation (ia|jb); the recorded window is the four-integer inclusive form."""
    n2 = parse_orca_json(APC_N2)
    assert n2.mo_iajb_window == (0, 6, 7, 16)
    assert n2.mo_iajb is not None and len(n2.mo_iajb) == 2425
    assert n2.mo_iajb[0] == (0, 0, 7, 7, pytest.approx(0.01160806221516865))
    assert all(0 <= i <= 6 and 0 <= j <= 6 for i, j, _, _, _ in n2.mo_iajb)
    assert all(7 <= a <= 16 and 7 <= b <= 16 for _, _, a, b, _ in n2.mo_iajb)
    h2 = parse_orca_json(APC_H2)
    assert h2.mo_iajb == ((0, 0, 1, 1, 0.18121046221962572),)
    assert h2.mo_iajb_window == (0, 0, 1, 1)


def test_exports_without_the_window_blocks_are_accepted():
    """The older exports carry no Fock family and no 2elIntegrals: documented variants,
    returned as None, not defects.  Whether the core Hamiltonian is present depends on
    the ``1elIntegrals`` request of the run (the mapping fixture asked for it, the
    entropy fixture did not)."""
    canonical = parse_orca_json(CANONICAL)
    assert canonical.hamiltonian is None
    assert canonical.coulomb is None and canonical.exchange is None
    assert canonical.coulomb_exchange is None
    assert canonical.mo_iajb is None and canonical.mo_iajb_window is None
    mapping = parse_orca_json(FIXTURES / "n2_scan_1.600.json")
    assert mapping.hamiltonian is not None
    assert mapping.coulomb is None and mapping.exchange is None
    assert mapping.mo_iajb is None


def test_an_eight_integer_window_record_is_reported(tmp_path):
    """The input syntax is eight integers, but the *record* the measurements found is
    the four-integer inclusive form; a record of any other shape is refused so the
    reader cannot silently mis-read a window."""
    def mutate(molecule: dict[str, Any]) -> None:
        molecule["2elIntegrals"]["OrbWin"] = [0, 6, 7, 16, 0, 0, 0, 0]

    variant = _apc_variant(tmp_path, mutate)
    with pytest.raises(ParserError, match="four inclusive integers"):
        parse_orca_json(variant)


def test_an_unmeasured_spin_key_is_reported(tmp_path):
    def mutate(molecule: dict[str, Any]) -> None:
        block = molecule["2elIntegrals"]["MO_IAJB"]
        block["beta/beta"] = block.pop("alpha/alpha")

    variant = _apc_variant(tmp_path, mutate)
    with pytest.raises(ParserError, match="spin keys") as excinfo:
        parse_orca_json(variant)
    assert "beta/beta" in str(excinfo.value)


def test_a_malformed_iajb_entry_is_rejected(tmp_path):
    def mutate(molecule: dict[str, Any]) -> None:
        molecule["2elIntegrals"]["MO_IAJB"]["alpha/alpha"][0] = [0, 0, 7, 0.5]

    variant = _apc_variant(tmp_path, mutate)
    with pytest.raises(ParserError, match="five-item"):
        parse_orca_json(variant)


def test_an_out_of_range_iajb_index_is_rejected(tmp_path):
    def mutate(molecule: dict[str, Any]) -> None:
        molecule["2elIntegrals"]["MO_IAJB"]["alpha/alpha"][0] = [0, 0, 7, 28, 0.5]

    variant = _apc_variant(tmp_path, mutate)
    with pytest.raises(ParserError, match="molecular orbitals"):
        parse_orca_json(variant)


def test_a_fock_block_of_the_wrong_size_is_rejected(tmp_path):
    def mutate(molecule: dict[str, Any]) -> None:
        molecule["K-Matrix"][0][0].pop()

    variant = _apc_variant(tmp_path, mutate)
    with pytest.raises(ParserError, match="K-Matrix"):
        parse_orca_json(variant)


# --- the density sidecar block (menus 23) -------------------------------------


def _ass1st_variant(tmp_path: Path, mutate: Callable[[dict[str, Any]], None]) -> Path:
    """Write a copy of the n2_ass1st export with ``mutate`` applied to its Molecule block."""
    document = json.loads(ASS1ST_ST.read_text(encoding="utf-8"))
    mutate(document["Molecule"])
    target = tmp_path / "ass1st_variant.json"
    target.write_text(json.dumps(document), encoding="utf-8")
    return target


def test_the_density_blocks_are_read():
    """The Densities request returns the named sidecar matrices: the CASSCF
    reference (Tdens-CAS.*), the NEVPT2 unrelaxed density per state
    (Tdens-CASNEV.*) and scfp; one square AO matrix per title."""
    st = parse_orca_json(ASS1ST_ST)
    assert st.densities is not None
    names = [name for name, _ in st.densities]
    assert names == [
        "Tdens-CAS.mult.1.root.0.p",
        "Tdens-CASNEV.mult.1.root.0.p",
        "scfp",
    ]
    assert all(len(block) == st.n_ao for _, block in st.densities)
    sa = parse_orca_json(ASS1ST_SA)
    assert sa.densities is not None
    assert len([n for n, _ in sa.densities if n.startswith("Tdens-CASNEV")]) == 3


def test_exports_without_densities_are_accepted():
    """The density request is optional, like the other request-dependent blocks."""
    assert parse_orca_json(CANONICAL).densities is None
    assert parse_orca_json(APC_N2).densities is None


def test_a_density_of_the_wrong_size_is_rejected(tmp_path):
    def mutate(molecule: dict[str, Any]) -> None:
        molecule["Densities"]["scfp"].pop()

    variant = _ass1st_variant(tmp_path, mutate)
    with pytest.raises(ParserError, match="scfp"):
        parse_orca_json(variant)


def test_an_unmeasured_density_shape_is_reported(tmp_path):
    def mutate(molecule: dict[str, Any]) -> None:
        molecule["Densities"]["scfp"] = {"alpha": molecule["Densities"]["scfp"]}

    variant = _ass1st_variant(tmp_path, mutate)
    with pytest.raises(ParserError, match="one square matrix per title"):
        parse_orca_json(variant)


# --- registry decision ------------------------------------------------------


def test_the_export_is_not_registered_as_a_program_output():
    """Design (module docstring): an orbital export is auxiliary data, not an engine
    output.  No registered parser may claim it, and parse_auto must refuse the file with
    its own actionable error; the reader is reachable only through its own name."""
    assert detect_program(CANONICAL) is None
    assert all(parser.program != "orca_2json" for parser in available())
    with pytest.raises(ParserError, match="Next step: "):
        parse_auto(CANONICAL)
