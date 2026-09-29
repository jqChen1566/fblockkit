# xTB fixtures (plan item 6.1; the reader behind menu 42)

All files are real xTB 6.7.1 output (101, 2026-09-30). The engine is
installed at `/home/chen_jianqi/xtb-dist/bin/xtb` on 101 (user-installed);
the runs below took under a second each. Commands (run in this directory;
`PATH` carries the xtb distribution):

```
xtb h2o.xyz --gfn2 > h2o_sp.out            # single point
xtb h2o.xyz --gfn2 --ohess > h2o_ohess.out # opt + numerical Hessian
xtb nh3_planar.xyz --gfn2 --hess > nh3_planar_hess.out  # Hessian at the input geometry
xtb h2o.xyz --gfn2 --opt --cycles 2 > h2o_noconv.out    # optimisation cut at two cycles
xtb butane.xyz --gfn2 --opt                # the optimisation whose xtbopt.xyz is shipped
```

Files:

- `h2o.xyz`, `nh3_planar.xyz`, `butane.xyz` — the input geometries;
- `h2o_sp.out` — single point on the distorted water; carries the result box
  (TOTAL ENERGY / GRADIENT NORM / HOMO-LUMO GAP), the dipole/quadrupole and
  WBO sections. Note: `normal termination of xtb` goes to **stderr**, so a
  stdout-only capture lacks it (the two `--ohess`/`--hess`/`--cycles`
  captures above were taken with `2>&1` and carry it inline);
- `h2o_ohess.out` — converged optimisation (5 cycles) plus the Hessian: nine
  frequency values in a doubly printed `eigval :` block (six
  translation/rotation zeros and the three water modes 1538.60 / 3642.93 /
  3657.89 cm-1), the thermochemistry summary and the closing marker;
- `nh3_planar_hess.out` — the imaginary-mode positive: planar ammonia at the
  input geometry carries one imaginary mode (-1337.66 cm-1; the engine
  reports "found 1 significant imaginary frequency", imag cut-off 5 cm-1);
- `h2o_noconv.out` — the failed-optimisation positive:
  `*** FAILED TO CONVERGE GEOMETRY OPTIMIZATION IN 2 ITERATIONS ***` after
  `--cycles 2`; the geometry is a partial step, not a minimum;
- `xtbopt.xyz` — the optimised butane geometry; its comment line carries
  `energy: ... gnorm: ... xtb: ...` (the engine's own annotated XYZ);
- `vibspectrum` — the `$vibrational spectrum` side file of the last Hessian
  run (planar ammonia): mode / symmetry / wavenumber / IR intensity / IR
  selection columns, `$end`-terminated.

A recorded boundary (not shipped as a fixture): linear water gives five zero
modes and a **positive** degenerate bend (707.64 cm-1 twice) rather than an
imaginary one -- the harmonic picture of this engine does not show the linear
bend instability; the report only relays what the file carries.
