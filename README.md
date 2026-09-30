# fBlockKit (the f-block calculation toolkit)

A zero-barrier workbench for strongly correlated systems and lanthanide/actinide
calculations. It does two things, and it never runs a calculation for you:

1. **Generates input files** -- system profile in, runnable ORCA input (method chain,
   basis set / ECP, active space, convergence settings) plus plain run guidance out;
2. **Reads results** -- an output file (or a structure file) in, a characterisation
   report out: orbital composition, occupation and entropy analysis, multi-reference
   character, local spin analysis, crystal-field fitting and point-charge estimates,
   cross-level consistency, SCF rescue, transition-state checks.

Every conclusion carries its provenance: a program-manual quote (with section and
URL), a literature reference (with DOI and a paste-ready BibTeX entry), or a measured
record. The decision layer is a deterministic rule engine -- no machine learning, no
network access, no file uploads.

## Status

v0.1 -- **early development**. The interfaces and the report formats may still change
between releases.

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

fBlockKit is developed by Jianqi Chen (陈建棋), College of Chemistry and
Materials Science, Jinan University, Guangzhou, China.

## Citing

If you use this software, cite it as described in `CITATION.cff`. The reports list
the complete citation and BibTeX entry of every reference they use.

## License

Apache-2.0 (see `LICENSE` and `NOTICE`). The internal design and research notes of
this project live outside this repository; the user guide shipped here is the
authoritative user-facing documentation.

## Development

```bash
# unified check: the layering contract, the 3.11 syntax gate, then the tests
bash scripts/verify.sh
```
