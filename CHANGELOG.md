# Changelog

All notable changes to this project are documented in this file.
The format follows Keep a Changelog; versions are date-agnostic semantic versions.

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
