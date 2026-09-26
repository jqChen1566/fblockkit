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
            "n2_diffuse.out",
            "fhh_optts_freq.out",
            "ce3_orbcomp.out",
        ):
            shutil.copy(FIXTURES / name, inputs / name)
        shutil.copy(FIXTURES / "inputs" / "scf_noconv.inp", inputs / "scf_noconv.inp")
        shutil.copy(REPO / "fixtures" / "literature" / "pucl3_s18.json", inputs / "pucl3_s18.json")
        # the exact-entropy chain (menu 12): the converged CASSCF output, the
        # FCIDUMP it dumped, and the two orca_2json exports of the gbw
        shutil.copy(FIXTURES / "n2_fcidump_step_a.out", inputs / "n2_fcidump_step_a.out")
        shutil.copy(FIXTURES / "n2_fcidump.fcidump", inputs / "FCIDUMP")
        shutil.copy(FIXTURES / "n2_fcidump.canonical.json", inputs / "canonical.json")
        shutil.copy(FIXTURES / "n2_fcidump.localized.json", inputs / "localized.json")
        subprocess.run(
            [PY, str(EXAMPLES / "generate_inputs.py")], cwd=work, env=ENV, check=True,
            capture_output=True,
        )
        for name in (
            "octahedron.xyz",
            "ce_atom.xyz",
            "ceo6.xyz",
            "cf_c3.json",
            "dy_doublets.json",
            "mapping_manifest.json",
            "guess_manifest.json",
            "n2_scan_1.094.loc.json",
            "n2_scan_1.600.loc.json",
            "n2_scan_1.610.loc.json",
            "n2_scan_2.600.loc.json",
            "n2_scan_1.600.json",
            "n2_scan_1.600.mkl",
        ):
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
                    "1", "work/n2_diffuse.out",
                    "2", "work/octahedron.xyz", "Oh",
                    "3", "Ce", "3", "2", "3", "energy", "work/ce_atom.xyz", "1", "1,7,2,1", "default", "5000",
                    "4", "Pu,Cl", "0", "2", "3", "energy",
                    "5", "dmrg",
                    "6", "openmolcas",
                    "6", "liblan",
                    "7", "work/pucl3_s18.json",
                    "9", "work/scf_noconv.out", "work/scf_noconv.inp",
                    "10", "work/cf_c3.json",
                    "11", "work/ceo6.xyz", "O=-2", "",
                    "12", "work/n2_fcidump_step_a.out", "work/FCIDUMP",
                    "work/canonical.json", "work/localized.json", "4 9", "0",
                    "13", "work/canonical.json", "4 9", "work/localized.json", "4 9",
                    "14", "work/canonical.json", "0", "p", "2", "", "",
                    "15", "work/canonical.json", "",
                    "16", "work/dy_doublets.json",
                    "17", "work/mapping_manifest.json",
                    "18", "work/guess_manifest.json",
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

            diffuse = read("n2_diffuse.out.fbk.md")
            check(
                "A8 diffuse-orbital (Rydberg) check" in diffuse
                and "0.18764592" in diffuse
                and "ranking (most diffuse channel first)" in diffuse,
                "menu 1 report: the rule-G panel ranks the active orbitals by channel "
                "diffuseness with the run's own exponents",
            )

            check(
                "DMET embedding recipe for Dy" in out and "42x" in out,
                "menu 6 prints the DMET workflow recipe with the liblan entry (the "
                "source's gates included)",
            )

            check(
                "Relativistic tier:" in out and "SFX2C-1e" in out,
                "menu 3 adds the relativistic tier line the f-block profile implies",
            )

            portrait = read("canonical.json.portrait.fbk.md")
            check(
                "A9 orbital portrait" in portrait
                and "bonding" in portrait
                and "Jaccard" not in portrait,
                "menu 15 report: the descriptor panel names the N2 bonding labels",
            )

            geometry = read("octahedron.xyz.fbk.md")
            check("Coordination shell: 6 ligands" in geometry and "Oh" in geometry,
                  "the geometry report finds the octahedral shell and the Oh parameter count")

            generated = read("ce_atom.fbk.inp")
            check(
                ("SARC2-DKH-QZVP" in generated and "%casscf" in generated and "nel 1" in generated
                 and "%maxcore 5000" in generated and generated.isascii()),
                "the generated ORCA input is ASCII with the SARC2 basis, the CASSCF block "
                "and the chosen MaxCore",
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

            exact = read("FCIDUMP.fbk.md")
            check(
                all(
                    marker in exact
                    for marker in (
                        "A2x exact four-state single-orbital entropy",
                        "-108.950671945",
                        "IAO-IBO",
                        "Cross-checks against the engine",
                    )
                ),
                "menu 12 report: the exact entropy route reproduces the engine numbers",
            )
            check(
                "Environment spin-polarisation entropy" in exact
                and "Delta S_E = 0.000000" in exact
                and "partitioned at centre" not in exact,
                "menu 12 report: the environment-spin section runs from the localised pair",
            )

            space = read("canonical.json.fbk.md")
            check(
                "sigma_F = ||M||_F / sqrt(min(n_A, n_B)) = 1.000000" in space
                and "essentially the same space" in space,
                "menu 13 report: the same space before and after a localisation reads unchanged",
            )

            avas = read("canonical.json.avas.fbk.md")
            check(
                "AVAS target projection" in avas
                and "recommended active space" in avas
                and "6 orbitals" in avas,
                "menu 14 report: the N 2p target reads back the textbook (6, 6) space",
            )

            magnetic = read("dy_doublets.json.magnetic.fbk.md")
            check(
                "A10 magnetic-doublet criterion" in magnetic
                and "facilitates QTM" in magnetic
                and "supports excitation" in magnetic
                and "DEGREES" in magnetic
                and "15.31" in magnetic and "30.25" in magnetic,
                "menu 16 report: the source's calibration row gets both sides of "
                "the g_T*theta_3 line, with the domain and the gap",
            )

            mapping = read("mapping_manifest.json.mapping.fbk.md")
            check(
                "A11 cross-structure orbital mapping" in mapping
                and "non-matchable" in mapping
                and "r=1.094: [4, 5, 6]" in mapping
                and "tau: occupied 0.5, virtual 0.5" in mapping,
                "menu 17 report: the N2 scan maps its cores and puts the hybrids in "
                "the non-matchable block; the selection implies the bond triad",
            )

            check(
                "Active-space overlap" in mapping
                and "0.995" in mapping
                and "preserved" in mapping,
                "menu 17 report: the active-space overlap block reads the scan's "
                "small step as preserved",
            )

            guess = read("guess_manifest.json.guess.fbk.md")
            # the residual is an eigh-level quantity (3.10e-15 locally, 4.66e-15 on
            # the 101 BLAS): the check confines it to the round-off scale, not a value
            check(
                "G5 WASP initial guess" in guess
                and "0.664" in guess
                and "orthonormalisation residual:" in guess
                and "e-15" in guess,
                "menu 18 report: the 1/d weights and the orthonormalisation residual",
            )
            written = read("n2_scan_1.600.fbk.mkl")
            check(
                "$COEFF_ALPHA" in written and "$BASIS" in written and "$OCC_ALPHA" in written,
                "menu 18 wrote the gbw-ready mkl next to the template",
            )

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
