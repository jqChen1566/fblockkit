# OpenMolcas fixtures

The first non-ORCA engine fixtures.  Kept as produced (no trimming), one
directory per chain.

## `dy_smoke` — Dy(III) 4f9 single ion, RASSCF + RASSI-SO + SINGLE_ANISO

| File | Content |
|---|---|
| `dy_smoke.in` | The full chain input, verbatim (Gateway `Group = C1`, `&SEWARD`, `&RASSCF` 9-in-7, `&RASSI` with `SpinOrbit`/`MEES`/`PROP`/`EJOB`, `&SINGLE_ANISO` MLTP 1 2). |
| `dy_smoke.out` | The complete driver output, verbatim (6,226 lines; ends with `Happy landing!`, Wall 34 s). |

**Provenance.** OpenMolcas v26.06 (serial + OpenMP build, thread count
pinned to 1), Dy at the origin, ANO-RCC-VDZP, charge 3, multiplicity 6
(4f9 sextet).  Produced on the project's 101 server, 2026-10-01, and
retrieved with the chain that regenerates it recorded here.

**Why C1.**  Two measured reasons.  (1) In the full D2h symmetry the
Lucia CI kernel returned `Number of CSFs 0` whenever the state-symmetry
irrep carried no active orbital (measured: adding one ag active orbital
restores the count), which aborts RASSCF with "You can't ask for more
roots than there are configurations"; C1 avoids the boundary entirely.
(2) The g_T criterion's calibration set (menu 16, Chilton S2) uses
9-in-7 SA-CASSCF-SO (Open)Molcas calculations, so a C1 single-ion chain
is the same-methodology shape.

**Physics notes.**  The 66 spin-orbit states reproduce the free-ion 6H
J-multiplets (J = 15/2 at 0, 13/2 near 3013 cm^-1, then 5626, 7836,
9644, 11050 cm^-1 -- 16+14+12+10+8+6 = 66).  The bare ion has no ligand
field, so the J = 15/2 manifold is fully degenerate (states 1-16 at
0.000 cm^-1) and the pseudospin construction of the first two states is
representational, not physical: the tunnelling splitting is 0 and the
g-tensor (0.9132 / 1.9906 / 13.1443) carries the degeneracy artifact.
The fixture's purpose is the *parse chain* (menu 16's OpenMolcas input
route); the ligand-field-bearing multi-doublet case now
lives in fixtures/single_aniso/dy_acac.* (the real SMM Dy(acac)3(H2O)2,
ORCA on the 84 server).

**MCSCF quality checks recorded with the fixture.**  `Number of CSFs 21`
(= C(7,5), the sextet count of f9 in seven orbitals), convergence after
9 macro-iterations; the state-averaged RASSCF energy is printed
identically for all 11 roots (-12145.07209494 Eh) as expected for a
state average.

**Not shipped here (regenerable from the input):** `dy_smoke.aniso`,
`dy_smoke.rassi.h5` (binary products; the parser reads the text output
only).
