# Orca_Magrelax fixtures (fixtures/magrelax/)

A real ORCA 6.1.1 Orca_Magrelax run (CO+, CAS(9,8), `def2-SVP`,
terminated normally, 6 min 41 s on 101, 2026-09-28), read by
`src/fblockkit/analysis/relaxation.py` (menu 36).

## The file set

- **co_hess.inp / co_magrelax.hess** -- the CASSCF numerical-frequency run
  that produced the `.hess` file magrelax reads (6 x 6 Hessian, one
  vibration at 2299.92 cm-1, the C-O stretch; the five zero modes are
  translations/rotation).  This is also the file the pysisyphus boundary
  probes used: ORCA 6.1.1's `.hess` carries the `$multiplicity`
  block.
- **co_magrelax.inp / co_magrelax.out** -- the magrelax run: CASSCF with
  `projectHSOC true`, `projectedstates 4`, `rel DoSOC/DoMagrelax true`,
  `%qgprop` (MODE DERIV, COORDSTYPE ENWNORMAL, STEPSIZE 0.00005, DERHSOC
  true) and `%magrelax` (magfld 0 0 0.01, temperature 2-30 step 1,
  doRaman false).
- **co_magrelax.magrelaxinp** -- the editable standalone input the module
  writes (`orca_magrelax <file>.magrelaxinp` reruns from it); its format
  matches the manual's example line for line.

## Measured facts

- The QGPROP stage generates **two displaced geometries per vibrational
  mode** (CO+ has one mode; the run computed `geom-0`/`geom-1` singles
  points with the SOC data sidecars `*.soc{inp,lx,ly,lz,sx,sy}.magrelax`);
  the main program then reads `co_magrelax.socder.magrelax` (one
  derivative, 4 x 4 matrices).
- The `Relaxation rates` table prints `T (K) | Rate (s^-1) | tau (s)` for
  the temperature list; **the rate is identically zero here** (tau =
  -inf at every temperature).  Measured cause: the one-phonon (Orbach)
  channel needs a phonon energy matching a Zeeman gap -- the only mode
  (2299.9 cm-1) against the Zeeman eigenvalues (0.0, 14680.6, 15089.1,
  0.0 cm-1) has no resonance, and the module's frequency list prints
  empty.  A structure fact of a single-mode toy, kept as the boundary
  fixture; a real multi-mode SMM run is a registered candidate increment.
- `Eigenvalues of zeeman in cm-1` are printed in the module's own
  projected (4-state) basis.

## Regenerating

The inputs are stored verbatim; the run needs the converged CO+ CASSCF
gbw of `fixtures/single_aniso/` (or the pNMR chain's `co_plus_qdpt.gbw`)
and the `.hess` file above.
