# ROCIS fixtures (fixtures/rocis/)

Real ORCA 6.1.1 ROCIS outputs (the core-excited-spectra protocol of the
manual's section 5.7), read by `src/fblockkit/parsers/rocis_spectra.py` and
`src/fblockkit/analysis/xas.py` (menu 37).  All runs on 101: the XAS/RIXS
pair 2026-09-29, the XES pair 2026-10-02; the CAS-CI two-step pairs on the
84 server.

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
- **fecl4_xes.{inp,out}** (2026-10-02) -- the off-resonance XES recipe, the
  menu-37 writer's own bytes (the same system; `OrbWin 6,6,7,8,0,2000` --
  the two spin-orbit-split core ranges of the Fe 2p shell plus a wide
  acceptor; `DoRIXS true`, `DoRIXSSOC false`, `DoElastic true`; NRoots 30),
  terminated normally on 101 in ~56 s.  The output carries the emission
  tables with the ground-state rows (`<root>-<mult>A -> 0-<mult>A`) at
  744.2/744.3 eV; `orca_mapspc <out> XES -x0720 -x1760 -w1.5 -eV -n800`
  renders 56 sticks to `<out>.XES.stk/.dat` (measured).  The SOC-corrected
  flavour (a run with `DoRIXSSOC true`) renders with the `XESSOC` mode
  (measured: 2376 sticks at NRoots 10).
- **fecl4.xyz** -- the probe geometry (the structure input of the XES leg).
- **fecl4_casci_xas.step1.{inp,out} / .step2.{inp,out}** (2026-10-02;
  refreshed the same day with the parallel/`%moinp` fix) -- the
  two-step CAS-CI core-excited XAS protocol (ORCA manual section 3.13.18) on
  the same system: step 1 is the valence SA-CASSCF(6,5) (122 s with
  `%pal nprocs 48 end` on the 84 server; its orbital table carries the Fe 2p
  at indices 6-8, -27.05 Eh); step 2 (`MOREAD` + `%scf rotate {6,42,90,0,0}
  {7,43,90,0,0} {8,44,90,0,0}` + `FrozenCore FC_NONE` + nel 12 norb 8 +
  `maxiter 1` + DoSOC; 74 s) moves the core into the positional window 42-49
  -- (96-12)/2 = 42, the manual's own example lands on 87 = (185-11)/2 for
  Fe(acac)3 -- and prints the L-edge states at 719.358864 eV (0.5 eV from
  the ROCIS result of the same system, 718.87 eV; the refresh reproduced it
  digit-for-digit).  `orca_mapspc <out> SOCABS` renders 785 peaks
  (the SOC-corrected table) and `ABS` 19 (the plain one); `XAS`/`XASSOC`
  refuse these table titles (measured).  Both inputs are the menu-37
  writer's own bytes, carrying the writer's `%pal nprocs 48 end` choice and
  a step-2 `%moinp` naming the gbw the step-1 run produces
  (`fecl4.casci_xas.step1.gbw`; the stored copies rename the walk's files
  per the fixtures/ convention).
- **fecl4_casci_xes.step2.{inp,out}** (2026-10-03) -- the two-step RAS-CI
  core-excited XES protocol (ORCA manual section 3.13.19) on the same
  system, the K-beta shape: the Fe 1s and 3p (indices 0/26/27/28, -7102 /
  -65 eV in the step-1 orbital table) rotated into the positional window
  41-49 ((96-14)/2 = 41) with `refs ras(14:4 1/5/0 0)` (at most one hole in
  the rotated core set), the `XESSOC`/`XASMOs 41` rel block and
  `nroots 40,40`.  Measured (84 server, `%pal nprocs 48 end`): 10m29s to
  `****ORCA TERMINATED NORMALLY****`, with the SOC-corrected emission
  blocks at the Fe K-beta scale (the main `319-5.0A -> 0-5.0A` line at
  7086.908295 eV; the experimental K-beta1 is 7058 eV; the def2-SVP level
  accounts for the offset).  `orca_mapspc <out> XESSOC -w5.0 -eV -n800`
  detects the XESSOC spectrum (7375 peaks) and writes
  `<out>.XESSOC.stk/.dat`; the mode's own window (5000-7200 eV) covers the
  K-beta region and `-x` overrides did not shift it (measured); the
  `ABS`/`SOCABS` modes refuse these table titles (measured).  The step-1
  half equals `fecl4_casci_xas.step1.inp` byte-for-byte (same valence
  space), so only the step-2 half is stored.  **Performance boundary
  (measured 2026-10-03)**: 20 roots miss the K-beta core hole entirely;
  40 roots complete in 8m44s and print the emission blocks; 60 roots in
  18m35s; the manual's saturated recipe (`1000,1000`, raised by the engine
  to the 290-state restricted space) grows superlinearly in the QDPT
  transition-density stage and **did not finish within 24 h** (killed by
  its own clock).  **Directory discipline (measured)**: rerunning the same
  basename over the killed run's stale `transition.densities*.tmp` residue
  aborts immediately in `TDensityContainer::RetrieveMetaData`
  (`qcdenstore.cpp:910`) -- run each step-2 job in a fresh directory.

## Measured window semantics (the probing record)

- The 4-element `OrbWin` (`donor_start, donor_end, acceptor_start,
  acceptor_end`) selects the 2p donor for XAS; every absorption block comes
  out with the core excitations from ~715 eV up.
- The **6-element** window is the RIXS form (manual section 7.31.4): the
  first four elements are the ranges of the TWO donor spaces, the last two
  the acceptor range.  Measured (2026-10-01/02): `6,6,7,8,0,2000` -- the
  spin-orbit-split Fe 2p donors (2p1/2, 2p3/2) with a wide acceptor --
  opens both RIXS channels; a valence donor range first
  (`33,42,6,8,0,2000`) gives the engine's zero-states refusal; wide
  valence donors (`6,8,33,42,0,2000`, `6,8,29,43,0,2000`) drive the SOC
  machinery into impractical territory (the SOC reduced-matrix and QDPT
  stages ran tens of minutes without finishing, and the
  transition-density storage grew past 36 GB at ~1 GB/min); an acceptor
  window over doubly occupied orbitals (`45,49`) segfaults (an earlier
  probe).
- The **KHD XES variant** (`DoKHDXESSOC` with `NIStates`/`IStates`,
  `NFStates`/`FStates`, `NXESMultiplets`/`XESMultiplets`) computes the
  emission intensities for every probed input and then aborts in its own
  printout (`PrintNXESKHDSpectrum`, a null stream, signal 11; independent
  of the state values, the list layout and the print level) -- a
  documented termination.
- Absolute excitation energies carry the method's own shift (the fixture's
  first core excitation prints 718.86 eV); the menu reports energies as
  printed.

## Regenerating

Inputs are stored verbatim; the runs need only the ORCA 6.1.1 binary (no
external files).
