## Polynuclear magnetism input (POLY_ANISO exchange cluster)

POLY_ANISO exchange-cluster input (Lines-type exchange + exact dipole-dipole):
  centre types: 2 (equivalent centres per type: 1, 1)
  total centres: 2
  spin-orbit functions per centre type: 2, 2  (exchange basis size 4)
  coupled pairs: 1 (J in cm-1: 1-2 0.1)
  dipole-dipole (COOR): included
  temperature grid: 0..300 K, 101 points

Input written: work\poly_aniso.input

Next steps:
  1. place one SINGLE_ANISO data file per magnetic centre next to this input, named aniso_1.input, aniso_2.input, ... (each is the <job>.CASSCF.anisofile of a menu-36 chain run; symmetry-equivalent centres of one type share a file)
  2. run:  otool_poly_aniso < work\poly_aniso.input > poly_aniso.output
  3. read the result back with this menu's first mode (the poly_aniso.output report)

Boundaries:
  - the J values and the coordinates are the caller's (the writer does not fit or judge them); the exchange model is the Lines-type Hamiltonian sum -J s_i.s_j over the declared pairs
  - the COOR block computes the dipolar coupling only for the declared pairs; coordinates are the symmetrised per-type positions in Angstrom
  - the plan covers the manual's NNEQ/PAIR/COOR/TINT blocks; the symmetry (SYMM) and anisotropic-coupling (LIN3/LIN9) variants are registered
