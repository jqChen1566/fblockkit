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

- **pip**: `pip install fblockkit` (Python >= 3.11; needs PyYAML, Jinja2, NumPy);
- **green package**: unpack the archive for your platform and run `fblockkit.exe`
  (Windows) or `fblockkit` (Linux) -- no Python installation needed.

## Use

```text
fblockkit                     interactive menu (numbered, append-only)
fblockkit --record FILE       interactive, recording your input as a replayable script
fblockkit run SCRIPT.txt      replay a script (same script -> byte-identical output)
fblockkit search KEYWORDS     search the external-tool index
fblockkit guide TOOL_ID       tool onboarding notes
```

The guide (`docs/USER_GUIDE.md`) is isomorphic to the menu: its section numbers are
the menu numbers. You need no programming to use the menu; every prompt shows its
default and pressing Enter accepts it.

## Documentation

- `docs/USER_GUIDE.md` -- the quick-start guide (section numbers are the menu
  numbers); it ships inside the program;
- the full manual (`docs/manual/`, built to `main.pdf` with
  `bash docs/manual/build.sh`; one chapter per menu entry, the criterion
  reference, the format reference).  Every session transcript and report
  excerpt in it is a replayable script under `docs/manual/examples/`, and the
  test suite replays them (`tests/test_manual.py`).  The built PDF ships with
  the release archives.

## What it does not do

- it never runs a calculation engine (you submit the generated input yourself);
- it never modifies your files: every product is a new file (`*.fbk.md` /
  `*.fbk.json` / `*.fbk.inp`);
- it uses no machine learning and no online service: the decision layer is a
  deterministic rule engine.

## Provenance and verification

- every rule carries non-empty evidence, and literature evidence must point into
  `src/fblockkit/knowledge/sources.bib` (DOI-verified entries) or loading fails;
- the test suite is the acceptance contract: a layering contract (import-linter)
  plus the tests, run with `bash scripts/verify.sh`;
- the parser regular expressions are fixed against real ORCA 6.1.1 output
  (`fixtures/orca/`, each fixture registered with its source in the fixture README).

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
