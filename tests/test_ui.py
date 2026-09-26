"""UI-layer tests: menu, session (interaction is a script), CLI, tool index.

Acceptance criterion (architecture §6 step 7): an interactive run can be exported as a
script and replays identically -- that is, "replaying the same script twice gives
byte-identical output", and "the recorded interaction script replays to the same output
as the original interaction". Every test writes only inside tmp_path and never touches
fixtures/.
"""

from __future__ import annotations

import io
import shutil
from pathlib import Path

import pytest

from fblockkit.toolindex import (
    ToolIndexError,
    guide,
    load_tools,
    search,
)
from fblockkit.ui import Session, load_menu, main, run_session, strip_comments
from fblockkit.ui.handlers import HANDLERS

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"


def _session_run(lines: list[str]) -> tuple[str, Session]:
    out = io.StringIO()
    session = Session(lines=lines, out=out)
    session.run(HANDLERS)
    return out.getvalue(), session


# --- menu --------------------------------------------------------------------


def test_menu_numbers_unique_and_handlers_exist():
    menu = load_menu()
    numbers = [item.number for item in menu]
    assert len(numbers) == len(set(numbers))
    assert numbers[:6] == ["1", "2", "3", "4", "5", "6"]
    assert "0" in numbers
    for item in menu:
        assert item.handler in HANDLERS, item.handler


def test_strip_comments_keeps_blank_lines():
    lines = ["# comment", "1", "", "# another comment", "0"]
    assert strip_comments(lines) == ["1", "", "0"]


# --- interaction is a script (acceptance) ------------------------------------


def _report_script(fixture_copy: Path) -> list[str]:
    return ["1", str(fixture_copy), "0"]


def test_replay_is_deterministic(tmp_path):
    """Replaying the same script twice gives byte-identical output (once paths are
    normalised), and identical product files."""
    source = FIXTURES / "n2_casscf_nevpt2.out"
    copy_a = tmp_path / "a" / source.name
    copy_b = tmp_path / "b" / source.name
    copy_a.parent.mkdir()
    copy_b.parent.mkdir()
    shutil.copy(source, copy_a)
    shutil.copy(source, copy_b)

    out_a, _ = _session_run(_report_script(copy_a))
    out_b, _ = _session_run(_report_script(copy_b))
    # the two runs differ only in the file path; everything else must be byte-identical
    assert out_a.replace(str(copy_a), "X") == out_b.replace(str(copy_b), "X")
    report_a = Path(str(copy_a) + ".fbk.md").read_text(encoding="utf-8")
    report_b = Path(str(copy_b) + ".fbk.md").read_text(encoding="utf-8")
    assert report_a.replace(str(copy_a), "X") == report_b.replace(str(copy_b), "X")
    assert "## Findings" in report_a


def test_recorded_interaction_replays_identically(tmp_path):
    """The recorded interaction (the complete input sequence) replays to output that is
    byte-identical to the original run -- the core acceptance criterion."""
    source = FIXTURES / "n2_casscf_orbcomp.out"
    copy = tmp_path / source.name
    shutil.copy(source, copy)

    out_original = io.StringIO()
    session = Session(lines=["1", str(copy), "0"], out=out_original)
    session.run(HANDLERS)
    script_lines = session.record  # "record it once and you have a script": the complete input sequence
    assert script_lines == ["1", str(copy), "0"]

    out_replay, _ = _session_run(script_lines)
    assert out_replay == out_original.getvalue()


def test_menu_eight_saves_replayable_script(tmp_path):
    """The script menu 8 writes to disk is replayed: its content is the input sequence
    from before the saving action."""
    copy = tmp_path / "orb.out"
    shutil.copy(FIXTURES / "n2_casscf_orbcomp.out", copy)
    script_path = tmp_path / "session.txt"
    out, _ = _session_run(["1", str(copy), "8", str(script_path), "0"])
    assert "Saved" in out
    recorded = script_path.read_text(encoding="utf-8").splitlines()
    # the saving interaction itself (the path line and the "0") does not enter the script;
    # the script ends with the input sequence of menu 1
    assert recorded == ["1", str(copy), "8"]
    # the script can be replayed (at menu 8 the second path prompt hits EOF -> an explicit
    # cancellation, no broken file is written)
    out_replay, _ = _session_run(recorded)
    assert "Cancelled (input finished)" in out_replay


def test_menu_one_writes_markdown_and_json(tmp_path):
    copy = tmp_path / "gen.out"
    shutil.copy(FIXTURES / "generated_ce3_sarc2.out", copy)
    out, _ = _session_run(["1", str(copy), "0"])
    assert "Report written" in out and "Data written" in out
    md = Path(str(copy) + ".fbk.md").read_text(encoding="utf-8")
    payload = Path(str(copy) + ".fbk.json").read_text(encoding="utf-8")
    assert "## Summary" in md
    assert '"program": "orca"' in payload


def test_menu_two_geometry(tmp_path):
    xyz = tmp_path / "oct.xyz"
    xyz.write_text(
        "7\ncomment\nCe 0 0 0\nO 2.4 0 0\nO -2.4 0 0\nO 0 2.4 0\nO 0 -2.4 0\nO 0 0 2.4\nO 0 0 -2.4\n",
        encoding="utf-8",
    )
    out, _ = _session_run(["2", str(xyz), "Oh", "0"])
    assert "Coordination shell: 6 ligands" in out
    assert "4 B_k^q allowed" in out


def test_menu_three_generates_input(tmp_path):
    xyz = tmp_path / "ce.xyz"
    xyz.write_text("1\ncerium\nCe 0 0 0\n", encoding="utf-8")
    lines = [
        "3",
        "Ce",       # element list
        "3",        # charge
        "2",        # multiplicity
        "3",        # valence
        "energy",   # targets
        str(xyz),   # structure
        "1",        # basis tier
        "1,7,2,1",  # active space
        "default",  # convergence tier
        "2500",     # MaxCore per process (MB)
        "0",
    ]
    out, _ = _session_run(lines)
    assert "Input file written" in out
    inp = tmp_path / "ce.fbk.inp"
    assert inp.is_file()
    text = inp.read_text(encoding="ascii")
    assert "SARC2-DKH-QZVP" in text
    assert "%casscf" in text and "nel 1" in text
    assert "%maxcore 2500" in text  # the chosen value reaches the input
    assert "Run guidance" in out
    assert "%maxcore 2500 MB" in out  # and the guidance quotes it


def test_menu_four_basis_query(tmp_path):
    out, _ = _session_run(["4", "Pu,Cl", "0", "2", "3", "energy", "0"])
    assert "Basis-set / ECP recommendation" in out
    assert "SARC-DKH-TZVPP" in out
    assert "Not applicable: 5f-in-core" in out


def test_menu_five_and_six_tools():
    out, _ = _session_run(["5", "mokit", "6", "openmolcas", "0"])
    assert "MOKIT" in out
    assert "LGPL" in out  # the license line of OpenMolcas


def test_invalid_menu_choice_is_reported():
    out, _ = _session_run(["99", "0"])
    assert "Invalid number: 99" in out


def test_menu_seven_cross_level_reachable():
    """Menu 7 (A7) is reachable and reproduces the literature case end-to-end."""
    records = Path(__file__).resolve().parents[1] / "fixtures" / "literature" / "pucl3_s18.json"
    out, _ = _session_run(["7", str(records), "0"])
    assert "cross-level" in out
    assert "No.10" in out  # the HF-second solution dropping out of the CCSD(T) top-3
    assert "10.1021/acs.jctc.4c01189" in out  # evidence line with the DOI


def test_menu_nine_scf_rescue(tmp_path):
    """Menu 9: triage a failed SCF output and write a corrected input."""
    out_copy = tmp_path / "scf_noconv.out"
    in_copy = tmp_path / "scf_noconv.inp"
    shutil.copy(FIXTURES / "scf_noconv.out", out_copy)
    shutil.copy(FIXTURES / "inputs" / "scf_noconv.inp", in_copy)
    out, _ = _session_run(["9", str(out_copy), str(in_copy), "0"])
    assert "not converged" in out.lower() or "NOT CONVERGED" in out
    assert "Corrected input written" in out
    fixed = list(tmp_path.glob("scf_noconv.fix_*.inp"))
    assert fixed, "expected at least one corrected input file"
    assert all("!SlowConv" in p.read_text(encoding="utf-8") or "moread" in p.read_text(encoding="utf-8") for p in fixed)
    # the original input is untouched
    assert (in_copy.read_text(encoding="utf-8")) == (FIXTURES / "inputs" / "scf_noconv.inp").read_text(encoding="utf-8")


def test_menu_nine_stays_silent_on_clean_run(tmp_path):
    copy = tmp_path / "n2_hf_clean.out"
    shutil.copy(FIXTURES / "n2_hf_clean.out", copy)
    out, _ = _session_run(["9", str(copy), "0"])
    assert "SCF looks healthy: no triage finding." in out


def test_menu_nine_prints_the_convergence_block_verbatim(tmp_path):
    """The informational rows (which the check mode does not enforce) are shown
    exactly as printed, on healthy runs too, and never interpreted here."""
    copy = tmp_path / "n2_hf_clean.out"
    shutil.copy(FIXTURES / "n2_hf_clean.out", copy)
    out, _ = _session_run(["9", str(copy), "0"])
    assert "SCF CONVERGENCE block, verbatim from the output" in out
    assert "Last DIIS Error" in out
    assert "Last Energy change" in out
    # and it is still reported as healthy: the raw DIIS row above tolerance does
    # not turn into a finding (it is informational under the default check mode)
    assert "SCF looks healthy: no triage finding." in out


def test_menu_ten_crystal_field_fit(tmp_path):
    """Menu 10: round-trip a known C3 Hamiltonian through the JSON menu path."""
    import json

    import numpy as np

    from fblockkit.analysis import crystal_field as cf

    true = {(2, 0): 1000.0, (4, 0): 10.0, (4, 3): 50.0, (6, 0): 1.0}
    J = 8.0
    matrix = cf.hamiltonian(true, J)
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    coefficients = [[[float(v.real), float(v.imag)] for v in row] for row in eigenvectors.T]
    payload = {
        "point_group": "C3",
        "J": J,
        "levels": [float(v) for v in eigenvalues],
        "coefficients": coefficients,
    }
    path = tmp_path / "cf.json"
    path.write_text(json.dumps(payload), encoding="utf-8")

    out, _ = _session_run(["10", str(path), "0"])
    assert "Report written" in out
    for (k, q), value in true.items():
        assert f"({k},{q:>2})" in out  # the parameter appears in the printed table
    assert "max |residual|" in out
    assert "Refused" not in out  # nothing refused at this level of the report
    report = Path(str(path) + ".fbk.md").read_text(encoding="utf-8")
    assert "A4 crystal-field fit" in report
    assert "References" in report
    assert "10.1021/acs.jpclett.5c02971" in report  # complete citation with DOI


def test_menu_ten_appends_the_declaration_check(tmp_path):
    """Menu 10 now carries the A5 projection-basis declaration check: the five items
    are reported, and a bundled comparison set is judged per (k, q)."""
    import json

    import numpy as np

    from fblockkit.analysis import crystal_field as cf

    true = {(2, 0): 1000.0, (4, 0): 10.0, (4, 3): 50.0, (6, 0): 1.0}
    J = 8.0
    matrix = cf.hamiltonian(true, J)
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    coefficients = [[[float(v.real), float(v.imag)] for v in row] for row in eigenvectors.T]
    payload = {
        "point_group": "C3",
        "J": J,
        "levels": [float(v) for v in eigenvalues],
        "coefficients": coefficients,
        "declaration": {
            "convention": "Stevens",
            "projection": "J = 8",
            "units": "cm^-1",
            "z_axis": "main symmetry axis",
            "origin": "metal site",
        },
        "compare": {
            "label": "Table S1 (L = 5)",
            "parameters": {"2,0": 1879.0, "2,2": -43.0},
            "declaration": {"projection": "L = 5"},
        },
    }
    path = tmp_path / "cf_decl.json"
    path.write_text(json.dumps(payload), encoding="utf-8")

    out, _ = _session_run(["10", str(path), "0"])
    assert "Projection-basis declaration: 5 of 5 item(s) recorded" in out
    report = Path(str(path) + ".fbk.md").read_text(encoding="utf-8")
    assert "A5 projection-basis declaration check" in report
    assert "comparable-roughly" in report  # the leading axial (2, 0) pair
    assert "10.1039/d5cs00493d" in report  # the A5 citation with its DOI


def test_menu_one_includes_the_mr_panel_and_the_local_spin_section(tmp_path):
    """Menu 1 picks up the two new output analysers: A3 (any CASSCF/CC output) and A6
    (an output whose input divided the molecule into fragments)."""
    copy = tmp_path / "n2.out"
    shutil.copy(FIXTURES / "n2_casscf_nevpt2.out", copy)
    _session_run(["1", str(copy), "0"])
    report = Path(str(copy) + ".fbk.md").read_text(encoding="utf-8")
    assert "## A3 multi-reference character" in report
    assert "Combined verdict" in report

    copy = tmp_path / "n2s.out"
    shutil.copy(FIXTURES / "n2_stretch_local_spin.out", copy)
    _session_run(["1", str(copy), "0"])
    report = Path(str(copy) + ".fbk.md").read_text(encoding="utf-8")
    assert "## A6 local spin analysis" in report
    assert "NOT the environment spin polarisation entropy" in report  # the boundary line


def test_menu_eleven_point_charge_estimate(tmp_path):
    """Menu 11: geometry-only output without radial moments, cm^-1 parameters with them."""
    xyz = tmp_path / "ceo.xyz"
    xyz.write_text("3\ncomment\nCe 0 0 0\nO 2.4 0 0\nO -2.4 0 0\n", encoding="utf-8")
    out, _ = _session_run(["11", str(xyz), "O=-2", "", "0"])
    assert "Centre: Ce" in out
    assert "Crystal-field parameters in cm^-1: not produced" in out
    assert "Report written" in out
    report = Path(str(xyz) + ".fbk.md").read_text(encoding="utf-8")
    assert "S2 point-charge crystal-field estimate" in report
    assert "References" in report
    # with radial moments the cm^-1 parameters are produced
    out, _ = _session_run(["11", str(xyz), "O=-2", "0.9,2.0,6.0", "0"])
    assert "Crystal-field parameters" in out and "(2,0)" in out
    assert "not produced" not in out


def test_menu_eleven_missing_charge_is_refused(tmp_path):
    """A structure element without a charge is refused with the missing symbols named
    (never silently charged as zero)."""
    xyz = tmp_path / "ceo.xyz"
    xyz.write_text("3\ncomment\nCe 0 0 0\nO 2.4 0 0\nO -2.4 0 0\n", encoding="utf-8")
    out, _ = _session_run(["11", str(xyz), "", "0"])
    assert "refused" in out
    assert "'O'" in out  # the missing element is named


def test_eof_exits_cleanly():
    out, session = _session_run(["1"])  # EOF before menu 1 was even given a path
    assert session.eof is True
    assert "(end of input, exiting.)" in out


def _fcidump_chain(tmp_path):
    fixtures = Path(__file__).resolve().parents[1] / "fixtures" / "orca"
    out_file = tmp_path / "n2.out"
    shutil.copy(fixtures / "n2_fcidump_step_a.out", out_file)
    dump_file = tmp_path / "FCIDUMP"
    shutil.copy(fixtures / "n2_fcidump.fcidump", dump_file)
    canonical = tmp_path / "canonical.json"
    shutil.copy(fixtures / "n2_fcidump.canonical.json", canonical)
    localized = tmp_path / "localized.json"
    shutil.copy(fixtures / "n2_fcidump.localized.json", localized)
    return out_file, dump_file, canonical, localized


def test_menu_twelve_exact_entropy_full_chain(tmp_path):
    """Menu 12: the exact four-state entropy with the inferred window and the
    localized step; the report carries the engine cross-checks."""
    out_file, dump_file, canonical, localized = _fcidump_chain(tmp_path)
    out, _ = _session_run(
        ["12", str(out_file), str(dump_file), str(canonical), str(localized), "", "0", "0"]
    )
    assert "Active window inferred: [4, 5, 6, 7, 8, 9]" in out
    assert "Cross-checks against the engine" in out
    assert "-108.950671945" in out  # the printed CASSCF energy, reproduced
    assert "IAO-IBO" in out
    assert "Report written" in out
    report = Path(str(dump_file) + ".fbk.md").read_text(encoding="utf-8")
    assert "A2x exact four-state single-orbital entropy" in report
    assert "References" in report
    # the environment-spin section runs from the localized pair (cluster centre 0;
    # N2 has no f-block element, so the atomic-term section stays absent)
    assert "Environment spin-polarisation entropy" in report
    assert "Delta S_E = 0.000000" in report  # a singlet has no polarisation
    assert "Atomic-term check" not in report


def test_menu_twelve_without_exports_skips_the_localized_step(tmp_path):
    out_file, dump_file, _, _ = _fcidump_chain(tmp_path)
    out, _ = _session_run(["12", str(out_file), str(dump_file), "", "", "4 9", "0"])
    assert "Cross-checks against the engine" in out
    assert "IAO-IBO" not in out


def test_menu_thirteen_orbital_space_identity(tmp_path):
    """Menu 13: the same active space before and after a localisation reads as
    unchanged (sigma_F = 1, smallest singular value 1) -- the check must be blind
    to a rotation inside one space."""
    canonical = tmp_path / "canonical.json"
    localized = tmp_path / "localized.json"
    shutil.copy(FIXTURES / "n2_fcidump.canonical.json", canonical)
    shutil.copy(FIXTURES / "n2_fcidump.localized.json", localized)
    out, _ = _session_run(
        ["13", str(canonical), "4 9", str(localized), "4 9", "0"]
    )
    assert "sigma_F = ||M||_F / sqrt(min(n_A, n_B)) = 1.000000" in out
    assert "essentially the same space" in out
    report = Path(str(canonical) + ".fbk.md").read_text(encoding="utf-8")
    assert "sayfutyarova2017avas" in report or "10.1021/acs.jctc.7b00128" in report
    assert "References" in report


def test_menu_thirteen_flags_a_shifted_window(tmp_path):
    """A window shifted by one orbital shares five of six directions: the report
    must name the deviation instead of reading 'contained'."""
    canonical = tmp_path / "canonical.json"
    shutil.copy(FIXTURES / "n2_fcidump.canonical.json", canonical)
    out, _ = _session_run(["13", str(canonical), "4 9", str(canonical), "5 10", "0"])
    assert "sigma_F = ||M||_F / sqrt(min(n_A, n_B)) = 0." in out
    assert "contained only approximately" in out or "misses a substantial part" in out


def test_menu_thirteen_reports_a_basis_mismatch(tmp_path):
    canonical = tmp_path / "canonical.json"
    localized = tmp_path / "localized.json"
    shutil.copy(FIXTURES / "n2_fcidump.canonical.json", canonical)
    shutil.copy(FIXTURES / "n2_fcidump.localized.json", localized)
    # drop one AO from the overlap of the second export: the two files then no
    # longer belong to the same basis and the reader must refuse the second one
    import json

    document = json.loads(localized.read_text(encoding="utf-8"))
    matrix = document["Molecule"]["S-Matrix"]
    document["Molecule"]["S-Matrix"] = [row[:-1] for row in matrix[:-1]]
    tampered = tmp_path / "tampered.json"
    tampered.write_text(json.dumps(document), encoding="utf-8")
    out, _ = _session_run(["13", str(canonical), "", str(tampered), "", "0"])
    assert "Orbital-space comparison failed" in out


def test_menu_twelve_refuses_a_non_casscf_output(tmp_path):
    fixtures = Path(__file__).resolve().parents[1] / "fixtures" / "orca"
    out_file = tmp_path / "hf.out"
    shutil.copy(fixtures / "n2_hf_clean.out", out_file)
    dump_file = tmp_path / "FCIDUMP"
    shutil.copy(fixtures / "n2_fcidump.fcidump", dump_file)
    out, _ = _session_run(["12", str(out_file), str(dump_file), "", "", "4 9", "0"])
    assert "no CASSCF section" in out


# --- CLI --------------------------------------------------------------------


def test_cli_search_and_guide(capsys):
    assert main(["search", "dmrg"]) == 0
    assert "DMRG" in capsys.readouterr().out
    assert main(["guide", "openmolcas"]) == 0
    assert "magnetic-property chain" in capsys.readouterr().out
    assert main(["guide", "no-such-tool"]) == 2


def test_cli_unknown_argument_shows_usage(capsys):
    assert main(["wat"]) == 2
    assert "Usage:" in capsys.readouterr().out


def test_package_main_entry_importable():
    """The package-level entry point of python -m fblockkit is importable (equivalent to
    the installed command)."""
    import importlib

    module = importlib.import_module("fblockkit.__main__")
    assert callable(module.main)


def test_cli_run_missing_script(capsys):
    assert main(["run", "does-not-exist.txt"]) == 2
    assert "script does not exist" in capsys.readouterr().out


# --- tool index --------------------------------------------------------------


def test_tool_index_loads_with_evidence():
    records = load_tools()
    assert len(records) >= 15
    relations = {record.relation for record in records}
    assert relations == {"absorb", "interface", "index"}
    for record in records:
        assert record.evidence, record.id
        assert record.status in ("active", "planned"), record.id


def test_only_integrated_tools_are_marked_active():
    """status=active must live up to its name: the only thing this version really uses is
    the autoCAS protocol (the A2 entropy spectrum).

    When a tool is genuinely wired in later (MOKIT handing orbitals over, the OpenMolcas
    magnetic-chain template, say), change it to active here -- this step is the mechanical
    gate for "registered != used".
    """
    active = {record.id for record in load_tools() if record.status == "active"}
    assert active == {"autocas"}


def test_tool_entry_without_status_rejected(tmp_path):
    bad = tmp_path / "tools.yaml"
    bad.write_text(
        "- id: x\n  name: X\n  purpose: p\n  relation: index\n  license: MIT\n  source: s\n"
        "  evidence:\n    - {kind: literature, text: t, ref: r}\n",
        encoding="utf-8",
    )
    with pytest.raises(ToolIndexError, match="status"):
        load_tools(bad)


def test_tool_index_requires_evidence(tmp_path):
    bad = tmp_path / "tools.yaml"
    bad.write_text(
        "- id: x\n  name: X\n  purpose: p\n  relation: index\n  status: planned\n"
        "  license: MIT\n  source: s\n  evidence: []\n",
        encoding="utf-8",
    )
    with pytest.raises(ToolIndexError, match="provenance"):
        load_tools(bad)


def test_tool_search_multi_term():
    hits = search("lanthanide DMET")
    assert any(record.id == "liblan" for record in hits)
    assert search("no-such-keyword-zz") == ()
    assert search("   ") == ()


def test_tool_guide_unknown_id():
    with pytest.raises(ToolIndexError, match="Next step: "):
        guide("no-such-id")
