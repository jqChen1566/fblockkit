"""End-to-end acceptance: exercise every menu and every CLI form on real inputs.

This is the release-level "walk the whole program" check, complementary to the
unit suite: it drives the actual command-line entry point as a subprocess, on a
temporary copy of the registered fixtures, and verifies the *products* --
report content, generated inputs, corrected inputs, JSON validity -- plus replay
determinism and the CLI exit codes.

Run from the repository root:

    python scripts/acceptance.py            (local checkout, src layout)

or with the package installed anywhere:

    python <repo>/scripts/acceptance.py

Exit code 0 = every check passed; 1 = at least one check failed (the failing
checks are listed at the end).
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
FIXTURES = REPO / "fixtures" / "orca"
EXAMPLES = REPO / "docs" / "manual" / "examples"

PY = sys.executable
ENV = dict(os.environ)
ENV["PYTHONUTF8"] = "1"
if (REPO / "src").is_dir():
    ENV["PYTHONPATH"] = str(REPO / "src")

RESULTS: list[tuple[bool, str]] = []


def check(ok: bool, label: str) -> None:
    RESULTS.append((bool(ok), label))
    print(f"{'PASS' if ok else 'FAIL'}  {label}")


def run_cli(args: list[str], cwd: Path, stdin: str | None = None) -> subprocess.CompletedProcess:
    return subprocess.run(
        [PY, "-m", "fblockkit", *args],
        cwd=cwd,
        env=ENV,
        input=stdin,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )


def norm(text: str) -> str:
    return text.replace("\\", "/").replace("\r\n", "\n")


def main() -> int:
    work = Path(tempfile.mkdtemp(prefix="fbk_acceptance_"))
    try:
        # --- inputs: registered fixtures + the manual's deterministic generator.
        # The session script refers to them as work/<name>, so they live in a
        # work/ subdirectory of the run directory.
        inputs = work / "work"
        inputs.mkdir()
        for name in (
            "n2_casscf_nevpt2.out",
            "n2_caspt2.out",
            "co_plus_soc.out",
            "scf_noconv.out",
            "n2_hf_clean.out",
            "generated_ce3_sarc2.out",
            "n2_stretch_local_spin.out",
            "fhh_optts_freq.out",
            "ce3_orbcomp.out",
        ):
            shutil.copy(FIXTURES / name, inputs / name)
        shutil.copy(FIXTURES / "inputs" / "scf_noconv.inp", inputs / "scf_noconv.inp")
        shutil.copy(REPO / "fixtures" / "literature" / "pucl3_s18.json", inputs / "pucl3_s18.json")
        subprocess.run(
            [PY, str(EXAMPLES / "generate_inputs.py")], cwd=work, env=ENV, check=True,
            capture_output=True,
        )
        for name in ("octahedron.xyz", "ce_atom.xyz", "ceo6.xyz", "cf_c3.json"):
            # the generator writes into its own examples/work directory
            shutil.copy(EXAMPLES / "work" / name, inputs / name)

        # --- one session that walks the whole menu -----------------------------
        script = work / "walk.txt"
        script.write_text(
            "\n".join(
                [
                    "1", "work/n2_casscf_nevpt2.out",
                    "1", "work/n2_caspt2.out",
                    "1", "work/co_plus_soc.out",
                    "1", "work/fhh_optts_freq.out",
                    "1", "work/n2_stretch_local_spin.out",
                    "2", "work/octahedron.xyz", "Oh",
                    "3", "Ce", "3", "2", "3", "energy", "work/ce_atom.xyz", "1", "1,7,2,1", "default",
                    "4", "Pu,Cl", "0", "2", "3", "energy",
                    "5", "dmrg",
                    "6", "openmolcas",
                    "7", "work/pucl3_s18.json",
                    "9", "work/scf_noconv.out", "work/scf_noconv.inp",
                    "10", "work/cf_c3.json",
                    "11", "work/ceo6.xyz", "O=-2", "",
                    "8", "work/saved.txt",
                    "0",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        first = run_cli(["run", str(script)], cwd=work)
        check(first.returncode == 0, f"the full-menu session exits 0 (got {first.returncode})")
        out = first.stdout

        # --- products and their content ----------------------------------------
        def product(name: str) -> Path:
            return inputs / name

        def read(name: str) -> str:
            path = product(name)
            assert path.is_file(), f"missing product {name}"
            return path.read_text(encoding="utf-8")

        try:
            report = read("n2_casscf_nevpt2.out.fbk.md")
            check(
                all(marker in report for marker in
                    ("## A1", "## A2", "## A3", "## Findings", "## Provenance", "## References")),
                "check-up report carries A1/A2/A3, findings, provenance, references",
            )
            check("DG-CASSCF-ENERGY-ONLY-CONVERGENCE" in report,
                  "the CASSCF energy-only-convergence warning is raised")
            payload = json.loads(read("n2_casscf_nevpt2.out.fbk.json"))
            check(payload.get("program") == "orca" and payload["findings"],
                  "the JSON companion parses and carries findings")

            caspt2 = read("n2_caspt2.out.fbk.md")
            check("CASPT2" in caspt2 and ("W0" in caspt2 or "weight" in caspt2.lower()),
                  "the CASPT2 check-up reports the reference-weight findings")

            soc = read("co_plus_soc.out.fbk.md")
            check("## A2" in soc or "## A1" in soc, "the SOC output is characterised")

            ts = read("fhh_optts_freq.out.fbk.md")
            check("D4-TS-ONE-IMAGINARY-MODE" in ts and "D4-BARRIER-REFERENCE-POINT" in ts,
                  "the transition-state findings fire on the OptTS fixture")

            local_spin = read("n2_stretch_local_spin.out.fbk.md")
            check("## A6 local spin analysis" in local_spin
                  and "NOT the environment spin polarisation entropy" in local_spin,
                  "the local-spin section and its Delta S_E boundary are present")

            geometry = read("octahedron.xyz.fbk.md")
            check("Coordination shell: 6 ligands" in geometry and "Oh" in geometry,
                  "the geometry report finds the octahedral shell and the Oh parameter count")

            generated = read("ce_atom.fbk.inp")
            check(
                ("SARC2-DKH-QZVP" in generated and "%casscf" in generated and "nel 1" in generated
                 and generated.isascii()),
                "the generated ORCA input is ASCII with the SARC2 basis and the CASSCF block",
            )

            fixed = sorted(p.name for p in inputs.glob("scf_noconv.fix_*.inp"))
            check(len(fixed) >= 2 and any("slowconv" in n for n in fixed)
                  and any("prescf" in n for n in fixed),
                  f"SCF rescue wrote the corrected inputs {fixed}")

            cf = read("cf_c3.json.fbk.md")
            cf_ok = "A4 crystal-field fit" in cf and "A5 projection-basis declaration check" in cf
            for label, value in (("(2, 0)", 1000.0), ("(4, 0)", 10.0), ("(4, 3)", 50.0), ("(6, 0)", 1.0)):
                match = re.search(rf"\({re.escape(label[1:-1])}\)\s+(-?\d+\.?\d*(?:e[-+]\d+)?)", cf)
                cf_ok = cf_ok and bool(match) and abs(float(match.group(1)) - value) < 1e-3
            check(cf_ok, "the crystal-field fit recovers the known parameters and runs the A5 check")

            pc = read("ceo6.xyz.fbk.md")
            check("S2 point-charge crystal-field estimate" in pc and "(4,0)" in pc and "(4,4)" in pc,
                  "the point-charge estimate reports the Oh non-zero pattern")

            check(product("saved.txt").is_file(), "menu 8 wrote the session script")
        except AssertionError as exc:
            check(False, f"product check: {exc}")

        # --- replay determinism -------------------------------------------------
        second = run_cli(["run", str(script)], cwd=work)
        check(second.returncode == 0 and norm(second.stdout) == norm(out),
              "replaying the same script is byte-identical")

        # --- CLI forms and exit codes ------------------------------------------
        search = run_cli(["search", "dmrg"], cwd=work)
        check(search.returncode == 0 and "DMRG" in search.stdout, "search returns the DMRG tools")
        guide = run_cli(["guide", "openmolcas"], cwd=work)
        check(guide.returncode == 0 and "OpenMolcas" in guide.stdout, "guide returns the OpenMolcas entry")
        bad_guide = run_cli(["guide", "no-such-tool"], cwd=work)
        check(bad_guide.returncode == 2, "an unknown tool id exits 2")
        missing = run_cli(["run", "does-not-exist.txt"], cwd=work)
        check(missing.returncode == 2, "a missing script exits 2")
        unknown = run_cli(["wat"], cwd=work)
        check(unknown.returncode == 2 and "Usage:" in unknown.stdout, "an unknown argument prints usage")

        record = work / "recorded.txt"
        recorded_run = run_cli(["--record", str(record)], cwd=work, stdin="5\ndmrg\n0\n")
        check(recorded_run.returncode == 0 and record.is_file(), "--record writes the session script")
        if record.is_file():
            lines = [line.strip() for line in record.read_text(encoding="utf-8").splitlines()]
            check(lines == ["5", "dmrg", "0"], "the recorded script is the exact input sequence")

        # --- the healthy-run negative control ----------------------------------
        clean_script = work / "clean.txt"
        clean_script.write_text("9\nwork/n2_hf_clean.out\n0\n", encoding="utf-8")
        clean = run_cli(["run", str(clean_script)], cwd=work)
        check(clean.returncode == 0 and "SCF looks healthy: no triage finding." in clean.stdout,
              "the SCF triage stays silent on the clean run")
    finally:
        shutil.rmtree(work, ignore_errors=True)

    failed = [label for ok, label in RESULTS if not ok]
    print(f"\n{len(RESULTS) - len(failed)}/{len(RESULTS)} checks passed")
    if failed:
        print("FAILED:")
        for label in failed:
            print(f"  - {label}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
