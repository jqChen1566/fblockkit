"""Manual gates: the manual's examples replay, and its structure matches the menu.

Two mechanical checks, in the spirit of tests/test_docs.py:

1. every menu number has a chapter in the manual (and no chapter labels a number
   that is not a menu entry);
2. every example script under docs/manual/examples/scripts/ replays to the
   captured output under examples/expected/ -- nothing in the manual is typed by
   hand. The captures were recorded on Windows and the program prints
   platform-dependent path separators, so both sides are normalised before the
   byte comparison (documented in the manual).

The example preparation mirrors examples/run-examples.sh: the fixture files are
copied into examples/work/ and the deterministic input generator runs there.
"""

from __future__ import annotations

import io
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from fblockkit.ui import load_menu
from fblockkit.ui.cli import main as cli_main

REPO = Path(__file__).resolve().parents[1]
MANUAL = REPO / "docs" / "manual"
EXAMPLES = MANUAL / "examples"
FIXTURES = REPO / "fixtures" / "orca"


def _normalise(text: str) -> str:
    """Platform artefacts in the captures: path separators and line endings.

    The captures were recorded on Windows, where a redirected stdout writes
    CRLF and the program prints backslash paths (inside the JSON companion a
    separator appears JSON-escaped as two backslashes); the comparison is about
    content, so both are normalised (the repository stores LF, and a Linux
    checkout sees LF).
    """
    text = text.replace("\r\n", "\n")
    text = text.replace("\\\\", "/")
    return text.replace("\\", "/")


_NUMBER = re.compile(r"[-+]?\d+(?:\.\d+)?(?:e[-+]?\d+)?")


def _lines_match(expected: str, produced: str) -> bool:
    """Two transcript lines match: identical text, or identical once every number
    is compared within the cross-BLAS reproducibility band.

    Rationale (measured on the crystal-field example): a least-squares solve is
    reproducible to ~1e-11 relative between BLAS builds, which shows up in the
    printed digits of coefficients that are numerically zero (1e-10) and in the
    last digit of a 4.5e-6 residual.  The tolerances below (rtol 1e-6, atol 1e-8)
    are orders of magnitude below any quantity the manual quotes, and the rest of
    the line must match exactly, so a real change in wording or in a meaningful
    digit still fails.
    """
    if expected == produced:
        return True
    expected_numbers = _NUMBER.findall(expected)
    produced_numbers = _NUMBER.findall(produced)
    if len(expected_numbers) != len(produced_numbers):
        return False
    skeleton_expected = _NUMBER.sub("N", expected)
    skeleton_produced = _NUMBER.sub("N", produced)
    if skeleton_expected != skeleton_produced:
        return False
    for left, right in zip(expected_numbers, produced_numbers):
        a, b = float(left), float(right)
        if abs(a - b) > 1e-6 * max(abs(a), abs(b)) + 1e-8:
            return False
    return True


def _matches(expected: str, produced: str) -> bool:
    expected_lines = _normalise(expected).splitlines()
    produced_lines = _normalise(produced).splitlines()
    if len(expected_lines) != len(produced_lines):
        return False
    return all(_lines_match(a, b) for a, b in zip(expected_lines, produced_lines))


def _prepare_work() -> None:
    work = EXAMPLES / "work"
    work.mkdir(exist_ok=True)
    # the menu-4 datasource leg queries the deployed basis library; the same
    # relative path run-examples.sh exports, so the library path printed into
    # the capture matches on replay (the comparison normalises separators,
    # not absolute/relative forms)
    os.environ["FBK_BASISDB"] = "../../../fixtures/basisdb/mini_basis.db"
    subprocess.run(
        [sys.executable, "generate_inputs.py"],
        cwd=EXAMPLES,
        check=True,
        capture_output=True,
    )
    for name in (
        "n2_casscf_nevpt2.out",
        "generated_ce3_sarc2.out",
        "n2_stretch_local_spin.out",
        "fhh_optts_freq.out",
        "n2_hf_clean.out",
        "n2_diffuse.out",
        "scf_noconv.out",
    ):
        shutil.copy(FIXTURES / name, work / name)
    shutil.copy(FIXTURES / "inputs" / "scf_noconv.inp", work / "scf_noconv.inp")
    shutil.copy(FIXTURES / "inputs" / "fhh_reopt.inp", work / "fhh_reopt.inp")
    shutil.copy(REPO / "fixtures" / "literature" / "pucl3_s18.json", work / "pucl3_s18.json")
    shutil.copy(FIXTURES / "n2_fcidump_step_a.out", work / "n2_fcidump_step_a.out")
    shutil.copy(FIXTURES / "n2_fcidump.fcidump", work / "FCIDUMP")
    shutil.copy(FIXTURES / "n2_fcidump.canonical.json", work / "canonical.json")
    shutil.copy(FIXTURES / "n2_fcidump.localized.json", work / "localized.json")
    # the cross-structure mapping chain of menu 17 and the menu-18 template pair
    for distance in ("1.094", "1.600", "1.610", "2.600"):
        shutil.copy(
            FIXTURES / f"n2_scan_{distance}.loc.json",
            work / f"n2_scan_{distance}.loc.json",
        )
    shutil.copy(FIXTURES / "n2_scan_1.600.json", work / "n2_scan_1.600.json")
    shutil.copy(FIXTURES / "n2_scan_1.600.mkl", work / "n2_scan_1.600.mkl")
    # the dipole-moment chain of menus 19/20 (reference, prep and the candidates)
    for name in (
        "h2o_dm_ref_pbe0.out",
        "h2o_dm_prep_mp2.out",
        "h2o_dm_casci_e6o6.out",
        "h2o_dm_casci_e6o7.out",
        "h2o_dm_casci_e6o8.out",
        "h2o_dm_casci_e8o7.out",
        "h2o_dm_casci_e8o8.out",
        "h2o_dm_casci_e10o8.out",
        "h2o_dm_casci_e6o6_sa4.out",
    ):
        shutil.copy(FIXTURES / name, work / name)
    # the APC ranking export of menu 21 and the ASS1ST round export of menu 23
    shutil.copy(FIXTURES / "n2_apc.json", work / "n2_apc.json")
    shutil.copy(FIXTURES / "n2_ass1st.json", work / "n2_ass1st.json")
    # the AEGISS benzene platform of menu 25
    shutil.copy(FIXTURES / "benzene.json", work / "benzene.json")
    shutil.copy(FIXTURES / "benzene.fcidump", work / "benzene.fcidump")
    shutil.copy(FIXTURES / "benzene.out", work / "benzene.out")
    # the menu-29 reference pair (the UKS CH4 wrong-convergence fixture)
    shutil.copy(FIXTURES / "ch4_diss_prop.mkl", work / "ch4_diss_prop.mkl")
    shutil.copy(FIXTURES / "inputs" / "ch4_diss_prop.inp", work / "ch4_diss_prop.inp")
    # the menu-31 property-file fixture (N2 SA-CASSCF)
    shutil.copy(FIXTURES / "n2_sa.property.txt", work / "n2_sa.property.txt")
    # the pysisyphus pair of menus 32/33 (a structure for the generator, and a
    # real run laid out flat so menu 33 can read the directory directly)
    pysisyphus = REPO / "fixtures" / "pysisyphus" / "h2o_opt"
    shutil.copy(pysisyphus / "h2o.xyz", work / "h2o_start.xyz")
    shutil.copy(pysisyphus / "h2o.xyz", work / "h2o_ts.xyz")
    for name in (
        "run_stdout.log",
        "optimization.trj",
        "final_geometry.xyz",
        "RUN.yaml",
        "optimization.h5",
        "h2o_opt.yaml",
    ):
        shutil.copy(pysisyphus / name, work / name)
    if (work / "qm_calcs").exists():
        shutil.rmtree(work / "qm_calcs")
    shutil.copytree(pysisyphus / "qm_calcs", work / "qm_calcs")
    # the menu-34 Judd-Ofelt dataset
    shutil.copy(
        REPO / "fixtures" / "judd_ofelt" / "babu2000_eu3.yaml",
        work / "babu2000_eu3.yaml",
    )
    # the menu-35 pNMR trio (structure, QDPT output, run file)
    pnmr = REPO / "fixtures" / "pnmr" / "co_plus"
    for name in ("co_plus.xyz", "co_plus_qdpt.out", "co_plus_qdpt_g.yaml"):
        shutil.copy(pnmr / name, work / name)
    # the menu-36 relaxation pair (the CO+ KD SINGLE_ANISO output and its magrelax run)
    shutil.copy(
        REPO / "fixtures" / "single_aniso" / "co_aniso2.out", work / "co_aniso2.out"
    )
    shutil.copy(
        REPO / "fixtures" / "magrelax" / "co_magrelax.out", work / "co_magrelax.out"
    )
    # the menu-37 core-excited-spectra fixture (the [FeCl4]2- ROCIS run)
    shutil.copy(
        REPO / "fixtures" / "rocis" / "fecl4_xas.out", work / "fecl4_xas.out"
    )
    # the menu-38 AILFT fixture (the Ni(2+) d8 free-ion run)
    shutil.copy(
        REPO / "fixtures" / "ailft" / "ni_ailft.out", work / "ni_ailft.out"
    )
    # the menu-39 polynuclear-magnetism fixture (the two-center POLY_ANISO probe)
    shutil.copy(
        REPO / "fixtures" / "poly_aniso" / "two_center_probe.out",
        work / "two_center_probe.out",
    )
    # the menu-40 hyperfine fixtures (DFT and CASSCF EPRNMR probes)
    shutil.copy(
        REPO / "fixtures" / "hyperfine" / "cef3_epr_dft.out", work / "cef3_epr_dft.out"
    )
    shutil.copy(
        REPO / "fixtures" / "hyperfine" / "cef3_epr_casscf.out",
        work / "cef3_epr_casscf.out",
    )
    # the menu-41 magnetocaloric fixture (the POLY_ANISO M(H) probe)
    shutil.copy(
        REPO / "fixtures" / "magnetocaloric" / "poly_mh.out", work / "poly_mh.out"
    )
    # the menu-1 Gaussian leg (the G09 TS-optimization probe)
    shutil.copy(
        REPO / "fixtures" / "gaussian" / "g09_h2co_ts.out", work / "g09_h2co_ts.out"
    )


@pytest.fixture(scope="module")
def prepared_examples():
    _prepare_work()
    return EXAMPLES


# --- structure: the manual covers the menu -----------------------------------


def test_manual_covers_every_menu_number():
    """Negative control: deleting a menu chapter (or its \\menulabel) fails here."""
    chapters = {
        path.name: path.read_text(encoding="utf-8")
        for path in (MANUAL / "chapters").glob("*.tex")
    }
    assert chapters
    for item in load_menu():
        if item.number == "0":
            continue  # quitting is documented in the "Running" chapter, not a chapter of its own
        marker = f"\\menulabel{{{item.number}}}"
        assert any(marker in text for text in chapters.values()), (
            f"the manual has no chapter for menu {item.number}"
        )
    numbers = {item.number for item in load_menu()}
    for name, text in chapters.items():
        for found in re.findall(r"\\menulabel\{(\w+)\}", text):
            assert found in numbers, f"{name} labels menu {found}, which is not in menu.yaml"


# --- the examples replay ------------------------------------------------------


def test_manual_examples_replay_byte_identically(prepared_examples, monkeypatch, capsys):
    scripts = sorted((EXAMPLES / "scripts").glob("*.txt"))
    assert len(scripts) >= 10
    monkeypatch.chdir(EXAMPLES)
    for script in scripts:
        capsys.readouterr()
        code = cli_main(["run", f"scripts/{script.name}"])
        captured = capsys.readouterr()
        produced = captured.out + captured.err
        expected = (EXAMPLES / "expected" / f"{script.stem}.txt").read_text(encoding="utf-8")
        assert code == 0, script.name
        assert _matches(expected, produced), (
            f"the replay of {script.name} differs from its capture -- refresh the "
            f"captures with examples/run-examples.sh if the program output changed"
        )


def test_manual_example_products_match_the_captures(prepared_examples, monkeypatch, capsys):
    """The report files quoted in the manual are part of the capture set."""
    monkeypatch.chdir(EXAMPLES)
    expected_dir = EXAMPLES / "expected" / "products"
    for expected in sorted(expected_dir.iterdir()):
        produced = EXAMPLES / "work" / expected.name
        if expected.is_dir():
            # a product directory (the menu-19 batch): its shape is checked by the
            # replay capture and by the acceptance walk; here only its presence
            assert produced.is_dir(), f"{expected.name}/ was not produced by the examples"
            continue
        assert produced.is_file(), f"{expected.name} was not produced by the examples"
        assert _matches(
            expected.read_text(encoding="utf-8"),
            produced.read_text(encoding="utf-8"),
        ), f"{expected.name} differs from its capture"
