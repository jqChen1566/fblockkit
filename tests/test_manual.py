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
    shutil.copy(REPO / "fixtures" / "literature" / "pucl3_s18.json", work / "pucl3_s18.json")
    shutil.copy(FIXTURES / "n2_fcidump_step_a.out", work / "n2_fcidump_step_a.out")
    shutil.copy(FIXTURES / "n2_fcidump.fcidump", work / "FCIDUMP")
    shutil.copy(FIXTURES / "n2_fcidump.canonical.json", work / "canonical.json")
    shutil.copy(FIXTURES / "n2_fcidump.localized.json", work / "localized.json")
    # the cross-structure mapping chain of menu 17 (three localized scan exports)
    for distance in ("1.094", "1.600", "2.600"):
        shutil.copy(
            FIXTURES / f"n2_scan_{distance}.loc.json",
            work / f"n2_scan_{distance}.loc.json",
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
        assert produced.is_file(), f"{expected.name} was not produced by the examples"
        assert _matches(
            expected.read_text(encoding="utf-8"),
            produced.read_text(encoding="utf-8"),
        ), f"{expected.name} differs from its capture"
