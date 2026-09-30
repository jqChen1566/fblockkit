# fixtures/poly_aniso -- the ORCA POLY_ANISO driver probe (menu 39)

Real two-center probe of `otool_poly_aniso` (POLY_ANISO v1.0.0, compiled
2025-12-04, shipped with ORCA 6.1.1), run on 101 on 2026-09-29.

## Contents

- `two_center_probe.polyinp` -- the driver input (`&POLY_ANISO` head, NNEQ
  block `2 T` / `1 1` / `2 2`, one PAIR with J = 0.1 cm^-1, COOR at z = 0
  and 3.7 A, TINT 0..300 K / 101 points), fed by redirection;
- `aniso_1.input`, `aniso_2.input` -- the per-center data files, byte copies
  of `co_aniso2.CASSCF.anisofile` / `co_aniso.CASSCF.anisofile` from the
  menu-36 SINGLE_ANISO probe (wave53 run on 101).  The names are mandatory:
  the program looks them up as aniso_N.input (a renamed file aborts in
  inquire_key_presence -- measured 2026-09-29);
- `two_center_probe.out` -- the driver output (750 lines), reproduced with

      otool_poly_aniso < two_center_probe.polyinp > two_center_probe.out

## What the probe demonstrates

- the full input format: the `&POLY_ANISO` head is mandatory (without it the
  read aborts in fetch_init with "Unexpected end of input file", STOP 128);
- the per-center echo cross-checks the menu-36 analysis: the g values
  (2.00000 / 2.00000 / 1.99952) and the spin-orbit spectrum
  (0 / 0 / 29361.197 / 29361.197 cm^-1) match the SINGLE_ANISO results of
  the same data;
- the interaction-matrix decomposition (LINES-1: Isotropic 32.708 % /
  Symmetric 93.839 % / Anti-Symmetric 11.155 %), the coupled-state table
  (relatives 0 / 8.224e-06 / 0.051287 / 0.117104 cm^-1), the chiT table
  (101 points; 0.50154932 -> 0.75042950 cm3 K mol-1) and the Van Vleck
  tensor sequence (Z main value 1.500582 at 0.0001 K).

## Boundary

The two centers are not the two halves of one physical cluster: this is a
format-and-flow probe (two independent Co(II) single-ion runs placed at
invented relative coordinates with an invented J).  It validates the file
formats, the mandatory data-file naming and the report parsing; it does not
represent a physical exchange calculation.  No such ab initio calculation
exists in the ORCA ecosystem anyway -- the exchange constants are user
input, and the LDF-CAHF / many-state PNO-CASPT2 route is outside this
tool's scope (registered as a documented termination in the wave 5.6
records).

## The generated-input acceptance run (2026-09-30, server 101, the same ORCA 6.1.1)

The menu-39 input writer (mode 2) was validated end to end: the plan
reproducing this fixture's probe (2 types x 1 centre, spin-orbit basis 2+2,
one J = 0.1 cm-1 pair, coordinates (0,0,0)/(0,0,3.7), TINT 0..300 K / 101
points) wrote `poly_aniso.input`; feeding it to `otool_poly_aniso` with this
directory's `aniso_1.input`/`aniso_2.input` returns rc = 0, `POLY_ANISO
finished sucessfully!`, and reproduces `two_center_probe.out` byte for byte
(0 differing lines).  The input format itself is the ORCA manual section
7.18 (the NNEQ flag, the per-type count lines, PAIR under the Lines-type
Hamiltonian, optional COOR/TINT).
