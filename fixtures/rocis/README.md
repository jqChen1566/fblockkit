# ROCIS fixtures (fixtures/rocis/)

Real ORCA 6.1.1 ROCIS outputs (the core-excited-spectra protocol of the
manual's section 5.7), read by `src/fblockkit/parsers/rocis_spectra.py` and
`src/fblockkit/analysis/xas.py` (menu 37).  All runs on 101, 2026-09-29.

## The samples

- **fecl4_xas.{inp,out}** -- [FeCl4]2- (x2c, x2c-SVPall, ROHF high-spin d6,
  `OrbWin 6,8,0,2000`, NRoots 30, `DoPNO`/`TCutPNO 1e-11`, `XASelems 0`,
  `rel DoSOC`, `DoHigherMult`/`DoLowerMult`/`DoRI`/`DecomposeFosc`),
  terminated normally.  The RIXS flags are off in this sample.  The output
  holds **fourteen absorption blocks**: seven plain (`dipole_length`,
  `dipole_velocity` and the five `combined_*` D2/M2/Q2 variants) and the
  same seven again SOC-corrected after the spin-orbit step; the SOC tables
  weight fosc by the initial-state population (column header
  `fosc(D2) (*population)`) and list the **full state-pair matrix** (2235
  rows, of which 935 carry fosc > 1e-6).  Ninety excitation-table rows
  (three spin blocks of 30).  The primary-block branching (largest-gap
  split) lands at 735.602 eV with clusters [663, 272] and ratio 56.355.
- **fecl4_xas_rixs.{inp,out}** -- the same system with the RIXS flags on
  (`DoRIXS`/`DoRIXSSOC`/`DoElastic`) and the 4-element `OrbWin` (NRoots 10,
  for a short run; the SOC tables hold the correspondingly smaller
  state-pair matrix): the
  measured **refusal mode**, printed by the engine as
  `Making the RIXS files ... WARNING!: Flag for RIXS property calculation
  was identified but there is zero number of intermediate and/or final
  states: No Cross-Section properties will be evaluated ...Skipping this
  part / TIP: Increase the number of NRoots and/or decrease or increase the
  acceptor orbital space`.

## Measured window semantics (the wave's probing record)

- The 4-element `OrbWin` (`donor_start, donor_end, acceptor_start,
  acceptor_end`) selects the 2p donor for XAS; every absorption block comes
  out with the core excitations from ~715 eV up.
- The **6-element** window (two donor spaces plus the acceptor; the RIXS
  variant's requirement) could not be pinned down this round: the tested
  variants were refused (`zero intermediate/final states`) or crashed --
  an acceptor window covering doubly occupied orbitals (e.g. `45,49`)
  segfaults, and a window pair like `46,49` turns the spectrum into
  valence transitions.  Generating RIXS inputs is a registered candidate
  increment; the engine's refusal text is carried verbatim by the report.
- Absolute excitation energies carry the method's own shift (the fixture's
  first core excitation prints 718.86 eV); the menu reports energies as
  printed.

## Regenerating

Inputs are stored verbatim; the runs need only the ORCA 6.1.1 binary (no
external files).
