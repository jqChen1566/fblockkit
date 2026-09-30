<div align="center">

# fBlockKit

**the f-block calculation toolkit** -- input generation and characterisation
analysis for strongly correlated systems and lanthanide/actinide calculations

[![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-blue)](LICENSE)
[![Python](https://img.shields.io/badge/python-%E2%89%A53.11-blue)](https://www.python.org/)
[![Release](https://img.shields.io/badge/release-v0.1.0-brightgreen)](https://github.com/jqChen1566/fblockkit/releases)
[![Status](https://img.shields.io/badge/status-early%20development-yellow)](#status)

</div>

A zero-barrier workbench for strongly correlated systems and lanthanide/actinide
calculations.  It writes the input file and reads the results back; running the
calculation itself stays with you.  Every conclusion it draws carries its
provenance: a program-manual quote (with section and URL), a literature
reference (with DOI and a paste-ready BibTeX entry), or a measured record.  The
decision layer is a deterministic rule engine -- no machine learning, no network
access, no file uploads.

## Functionality

**Input generation** (a system profile in, runnable files out):

- ORCA input from a few questions about your system: method chain, basis set /
  ECP, active space, convergence settings, plus plain run guidance (menus 3, 4);
- focused generators: DeltaSCF / MOM excited-state SCF (27), RAS / ORMAS model
  spaces (28), multistart and imaginary-mode restart inputs (29, 30), multireference
  recipes with their gates (6);
- an OpenMolcas magnetic-chain template (SA-CASSCF / RASSI / SINGLE_ANISO);
- external-program inputs: pysisyphus PES exploration (32) and MOKIT/automr (44).

**Characterisation analysis** (a file in, a provenance-carrying report out):

- the check-up report and SCF rescue: composition, active-space analysis,
  multi-reference diagnostics, oscillation/trajectory checks, wrong-solution
  warnings (1, 9);
- exact active-space analysis from an FCIDUMP (four-state entropy, environment
  spin entropy, atomic-term check; 12) and the selection family: AVAS (14),
  orbital portrait (15), cross-structure mapping (17), WASP guess transfer (18),
  dipole-moment selection (19, 20), APC ranking (21), ASS1ST (22, 23), QICAS
  (24), AEGISS (25), TNASS (26);
- crystal-field fitting (10) and point-charge estimates (11); cross-level
  solution consistency (7); CASSCF state data and transitions (31);
- magnetic and spectroscopic side: Judd-Ofelt intensities (34), pNMR shifts
  (35), relaxation / QTM (36), hyperfine / EFG (40), magnetocaloric effect
  (41), tunnelling-relaxation prediction (46), XAS/RIXS (37), AILFT (38),
  polynuclear magnetism (39);
- external-program readings: xTB campaigns (42), CREST conformer ensembles
  (43), MOKIT/automr runs (45), pysisyphus runs (33).

## Status

v0.1 -- **early development**. The interfaces and the report formats may still
change between releases.

## Install

Pick the distribution form that suits you (all are attached to the repository's
**Releases** page):

- **Windows green package** (no Python needed): download
  `fblockkit-<version>-windows-x64.zip`, unpack it anywhere and run
  `fblockkit.exe` (or just double-click it).  The folder is self-contained and
  portable.
- **pip** (Python >= 3.11; dependencies PyYAML, Jinja2, NumPy):
  `pip install fblockkit` -- then the `fblockkit` command is on your PATH from
  any directory.
- **from source** (for development): `pip install .` inside a checkout; the
  build/test scripts are described under Development below.

The compiled manual (one chapter per menu entry, the criteria and format
references) is attached as a PDF to the **manual release** on the Releases page,
and `docs/USER_GUIDE.md` ships inside the program itself.

## Use

### First five minutes

1. Start it: run `fblockkit` (or double-click `fblockkit.exe`).  You get a
   numbered menu; type a number and answer the prompts -- every prompt shows its
   default, and Enter accepts it.  No programming is involved.
2. **Already have an engine output?**  Choose menu 1 and give the path: a
   check-up report appears next to the file (`<file>.fbk.md`, plus a
   machine-readable `.fbk.json`) with the orbital and active-space analysis,
   multi-reference diagnostics and findings -- each conclusion carrying its
   provenance.
3. **Need an input file?**  Choose menu 3 (input generation) or menu 4 (basis
   set / ECP advice): a few questions about your system yield a runnable ORCA
   input, a run-guidance block and an optional OpenMolcas magnetic-chain
   template.  Run the input on your machine or cluster, then feed the output
   back to menu 1.
4. Reports are plain Markdown (read them in any editor); the manual's menus and
   the guide's section numbers match one-to-one, so a report's section can be
   looked up directly.

### A worked example

The whole session is a replayable script -- the same construction the manual's
worked examples use:

```text
# session.txt  --  replay with:  fblockkit run session.txt
1
n2_casscf_nevpt2.out
0
```

```console
$ fblockkit run session.txt
Choose a number (Enter = leave empty)> 1
Output file path (ORCA or Gaussian) (Enter = leave empty)> n2_casscf_nevpt2.out
Report written: n2_casscf_nevpt2.out.fbk.md
Data written: n2_casscf_nevpt2.out.fbk.json
Analysis sections: 3.
Diagnosis: Warning 3 (details in the report file).
```

The report it writes (excerpt):

```markdown
# fBlockKit report

## Summary

Subject: n2_casscf_nevpt2.out
Counts: Warning 3

## A1 orbital composition (per MO, Loewdin)

Dominant shells of the active orbitals (0.02 < occ < 1.98):
  MO    5  occ 1.7092  N1 p 48.8%; next: N2 p 48.8%
  MO    6  occ 1.7092  N1 p 48.8%; next: N2 p 48.8%
  MO    7  occ 0.2936  N1 p 48.2%; next: N2 p 48.2%
  MO    8  occ 0.2936  N1 p 48.2%; next: N2 p 48.2%
```

### The whole program as a script

```text
fblockkit                     interactive menu (numbered, append-only)
fblockkit --record FILE       interactive, recording your input as a replayable script
fblockkit run SCRIPT.txt      replay a script (same script -> byte-identical output)
fblockkit search KEYWORDS     search the external-tool index
fblockkit guide TOOL_ID       tool onboarding notes
```

Anything the interactive session does is scriptable: record a session
(`--record`), replay it later (`run`), and get byte-identical reports -- this is
also how the manual's worked examples are produced.  `search` and `guide` give
the onboarding notes for the external programs the toolkit pairs with.

## Documentation

- `docs/USER_GUIDE.md` -- the quick-start guide (section numbers are the menu
  numbers); it ships inside the program;
- the full manual (`docs/manual/`, built to `main.pdf` with
  `bash docs/manual/build.sh`; one chapter per menu entry, the criterion
  reference, the format reference).  Every session transcript and report
  excerpt in it is a replayable script under `docs/manual/examples/`, and the
  test suite replays them (`tests/test_manual.py`).  The built PDF ships with
  the release archives.

## Scope and guarantees

- calculations run outside the program: it writes the input file and the run
  guidance, and you submit the job;
- your files are read-only: every product is a new file (`*.fbk.md` /
  `*.fbk.json` / `*.fbk.inp`);
- the decision layer is a deterministic rule engine -- no machine learning, no
  online service, and every conclusion carries its provenance.

## Author

fBlockKit is developed by Jianqi Chen, Wen Zeng and Guanghui Song, College
of Chemistry and Materials Science, Jinan University, Guangzhou, China.

## Funding

This project is supported by the Jinan University Provincial College Students'
Innovation and Entrepreneurship Training Program (Project No. S202610559049).

本项目受暨南大学省级大学生创新创业训练计划资助（项目编号：S202610559049）。

## Citing

If you use this software, cite it as described in `CITATION.cff`:

```bibtex
@software{fblockkit,
  author  = {Jianqi Chen and Wen Zeng and Guanghui Song},
  title   = {{fBlockKit}: the f-block calculation toolkit},
  year    = {2026},
  version = {0.1.0},
  url     = {https://github.com/jqChen1566/fblockkit},
  license = {Apache-2.0},
}
```

The reports list the complete citation and BibTeX entry of every reference they
use.

## License

Apache-2.0 (see `LICENSE` and `NOTICE`). The internal design and research notes of
this project live outside this repository; the user guide shipped here is the
authoritative user-facing documentation.

## Development

```bash
# unified check: the layering contract, the 3.11 syntax gate, then the tests
bash scripts/verify.sh
```
