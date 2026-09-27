## CASSCF state data (property file)

CASSCF state data from n2_sa.property.txt
  method: CASSCF; final (SA) energy: -108.707074407 Eh; active space: (6e, 6o)

  state  block  root  mult  irrep        energy (Eh)      (relative, eV)
      0      0     0     1      -    -108.975666704          0.0000
      1      0     1     1      -    -108.587267972         10.5689
      2      0     2     1      -    -108.558288545         11.3574

  absorption transitions (electric-dipole route, density: see the section):
  initial -> final (irrep)   mult (i->f)      dE (eV)      dE (cm**-1)
  0 -> 1 (0->0)   1.0->1.0           10.5689        85243.7
  0 -> 2 (0->0)   1.0->1.0           11.3574        91603.9
  note: this transition carries further columns the manual's schema leaves unnamed (117.311, 0, 0, 0, 0, -1.33643e-12, 0, 0, 0); they are printed in the report file, not interpreted.

Boundaries (the Wave-4.2 availability survey):
  - the CI vectors are not persistable (run-time temporaries only; the .cis file belongs to the CIS/STEOM modules), so state identity along a series cannot be read from one output;
  - the .out side prints the per-root dominant CSF occupations twice (initial state check and final states block), and this property file adds the structured per-state energies and transitions;
  - per-state observables for identity tracking: single-root runs carry the per-state dipole (the menu-19 route), the FIC-NEVPT2 sidecar carries per-root densities (the menu-23 chain) -- a tracker across a series is registered as a candidate increment.
