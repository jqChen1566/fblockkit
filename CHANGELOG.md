# Changelog

All notable changes to this project are documented in this file.
The format follows Keep a Changelog; versions are date-agnostic semantic versions.

## [0.2.0] - Unreleased

*(release date set at the cut; the entry conditions are in
`packaging/RELEASE_CHECKLIST.md`)*

The second development snapshot: the external-program interface layer, the
f-block post-processing family, the X-ray-spectroscopy routes and the
reproducibility hardening on top of the 0.1.0 core.

### Input generation

- External-program interfaces: MOKIT automr input generation with its run
  reader, the CREST conformer-ensemble report and the CREST-to-ORCA upgrade
  generator, the xTB pre-screening run report, and the pysisyphus run report;
- Core-excited spectra (menu 37): the ROCIS XES writer (the six-element
  donor/acceptor window; the plain RIXS channel as the emission carrier) and
  the two-step CAS-CI/RAS-CI family -- the XAS pair (manual 3.13.18) and its
  RAS-CI XES counterpart (3.13.19; the `refs ras` saturation and the
  XESSOC/XASMOs request).  Both two-step writers ask for the MPI process
  count (`%pal`), and the step-2 `%moinp` names the gbw the step-1 run
  produces, so a generated pair runs as written;
- RAS/ORMAS library (menu 28's writer; not yet exposed as menu questions):
  the manual's arbitrary-CFG references (`{2 2 2 0 0 0}`), `irrep` lists and
  the `%rasci` `douv` / `rel dosoc` coupling lines (module level only; engine
  acceptance pending its probe round);
- POLY_ANISO input writer: the SYMM block for equivalent centres, the
  isotropic PAIR and axis-diagonal LIN3 pair models, plus the couplings
  checklist; the AOP rotation guess and the PiOS pi-orbital space write
  gbw-ready orbital files.

### Analysis

- f-block post-processing: Judd-Ofelt intensity parameters (menu 34), pNMR
  pseudocontact shifts with the optional Bleaney comparator (menu 35),
  magnetic relaxation / QTM metrics (menu 36), AILFT ligand-field analysis
  with a built-in free-ion reference table (menu 38), polynuclear magnetism
  (menu 39), hyperfine/EFG parameters and magnetic entropy/magnetocaloric
  quantities (menus 40/41), tunnelling-relaxation prediction (menu 46),
  cross-run state tracking (menu 49) and the TNASS subset selection with its
  energy-ranked best-k sweep (menu 26);
- MOKIT automr side products: the `.fch` formatted-checkpoint reader
  (column-major MO blocks and row-major density triangles pinned against the
  file's own density; stage notes derived from the file names);
- curve companions: plot-ready tidy CSVs beside the reports of menus 34, 36,
  37 and 41 (one shared formatting module);
- the two-step spectrum outputs are rendered with `orca_mapspc`
  (SOCABS/ABS for the XAS tables, XESSOC for the emission tables).

### Fixed

- the two-step writers' step-2 `%moinp` now names the gbw the generated
  step-1 input actually produces (previously a rename was needed);
- negative-zero normalisation in the orbital write-back routes;
- degenerate-subspace bases in the written orbital files are rebuilt
  deterministically (cross-machine byte reproducibility: AOP, PiOS, WASP and
  the quasi-natural-orbital write-back);
- the manual's own XES example carries `rel / DoVelocity`, which the 6.1.1
  input scanner rejects (measured); the generated inputs omit it.

### Documentation

- the manual grows to 263 pages (menus 34-49 chapters, the criterion table,
  the format chapter) with every example transcript replayable against the
  frozen fixtures; the user guide mirrors all added routes and their measured
  boundaries;
- fixtures: real engine outputs for every added route, each with its README
  recording the run conditions and the measured boundaries.

## [0.1.0] - 2026-09-26

Initial release (early development).

### Input generation (ORCA)

- Basis-set / ECP recommendation for f-block elements (element, valence and target
  aware; hard refusals for structurally impossible combinations such as 5f-in-core
  pseudopotentials for f-f transitions), with auxiliary-basis matching;
- f^n active-space template with a validation protocol (occupation windows, entropy
  caveats, the CASCI-as-identity trap) and starting-space suggestions;
- convergence plan (defaults first, TRAH only when asked for and only with its /C
  auxiliary basis, otherwise explicitly downgraded) and initial-guess guidance;
- ORCA input rendering (ASCII-checked, directly runnable) plus run guidance, and an
  OpenMolcas magnetic-property chain template (SA-CASSCF / RASSI / SINGLE_ANISO).

### Analysis (reads ORCA output or a structure)

- A1 orbital composition by shell (partitioned by occupation);
- A2 occupation spectrum and the single-orbital entropy **bound** (the literature
  four-state entropy is not determined by spin-summed occupations; the bound is used
  one-sidedly as an exclusion test) with plateau/cliff ranking;
- A3 multi-reference character panel (T1 diagnostic, fractional occupations, the
  entropy bound as an exclusion test; combined verdict with a completeness line);
- A4 crystal-field fit from sampled levels and projection coefficients (extended
  Stevens operators; rank/condition diagnostics);
- A5 projection-basis declaration check (the five items a B_k^q set must record, and
  the per-(k,q) comparison verdict: only the leading axial (2,0) may be compared
  roughly across schemes);
- A6 local spin analysis (ORCA fragment tables, total <S^2>, spin contamination, the
  metal-fragment share) with the explicit Delta S_E boundary and the two-step
  route to the real environment spin-polarisation entropy;
- A7 cross-level solution consistency (including "a cheap-level top solution drops
  out of the expensive-level top N");
- S1 coordination geometry and symmetry hints (symmetry -> number of CF parameters);
- S2 point-charge crystal-field estimate (exact multipole expansion, validated by
  independent numerical quadrature; radial moments user-supplied).

### Diagnosis (rule engine over parser facts)

- CASSCF/NEVPT2 check-up (convergence via gradient, active-space sanity, NEVPT2
  intruder signs), CASPT2 mandatory checks (reference weights, smallest energy
  denominators), SCF rescue (triage plus corrected input files; refuses
  "just raise MaxIter"), and seven frequency / transition-state rules (imaginary-mode
  reading, the small-imaginary-mode caveat, missing frequency verification, an
  unconverged TS optimisation, the barrier reference point).

### Documentation

- the full manual (`docs/manual/`, LaTeX source + built PDF): one chapter per
  menu entry following a fixed walkthrough template, the command-line and
  format references, the complete criterion/threshold reference with sources,
  and a development chapter. Every session transcript and report excerpt in it
  is a replayable example (`docs/manual/examples/`, gated by
  `tests/test_manual.py`).

### Infrastructure

- Parser layer for ORCA 6.1.1 output (SCF, orbitals, CASSCF, NEVPT2/CASPT2, SOC,
  T1 diagnostic, vibrational frequencies, local spin analysis, geometry-optimisation
  context), every pattern fixed against registered real fixtures;
- knowledge layer with a strict rule loader (every rule carries evidence; literature
  evidence must resolve in `sources.bib`);
- interactive menu with "the interaction is the script" replay, a searchable external
  tool index (relation x status), and a Markdown/JSON report with a References block
  (complete citations and paste-ready BibTeX);
- packaging: pip wheel and a PyInstaller green package, both exercised end to end.
