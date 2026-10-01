# AILFT fixtures (fixtures/ailft/)

Real ORCA 6.1.1 AILFT outputs (the CASSCF module's ab initio ligand-field
driver, manual section 3.13.16), read by `src/fblockkit/parsers/ailft.py`
and `src/fblockkit/analysis/ailft.py` (menu 38).  Both runs on 101,
2026-09-29, terminated normally.

## The samples

- **ni_ailft.{inp,out}** -- the manual's Ni(2+) d8 free-ion example
  (`!NEVPT2 def2-SVP`, `%casscf nel 8 norb 5 ActOrbs dOrbs`, `mult 3,1`,
  `nroots 10,15`, `rel dosoc`).  Two parameter levels: CASSCF-level Racah
  B = 1328.1, C = 4865.5, C/B = 3.663 (F0dd from the two-electron
  integrals, "(fixed)"); NEVPT2-level B = 1218.0, C/B = 3.706 (F0 and
  Racah A are not printed at this level).  The CASSCF-level fit RMS is
  0.0 cm^-1 -- **intrinsic**: the LFT parametrization is exact for that
  level (the manual's own note); the NEVPT2-level total RMS is 457.4 cm^-1
  (blocks 315.9 / 517.9), Pearson 1.000, and the SOC section fits
  ZETA_D = 664.14 cm^-1.
- **cef3_ailft.{inp,out}** -- the f-shell probe: CeF3 (C3v, DKH-free
  def2-SVP), `%casscf nel 1 norb 7 LFTCase 4f mult 2 nroots 7`.  The
  **f-shell VLFT table goes straight to its header** (no "Ligand
  field one-electron matrix VLFT" title line -- that line is d-shell
  only); a single f electron has no electron repulsion, so only F0ff
  prints (fixed) and Racah B = 0.0; the ligand-field eigenfunction spread
  is 2333.0 cm^-1 (CASSCF) and the SOC section gives ZETA_F = 638.98 cm^-1.
- **dy3_freeion_unconv.{inp,out}** -- the Dy(3+) f9 free-ion probe, first
  attempt (`!NEVPT2 def2-SVP def2-SVP/C TightSCF SlowConv`, `LFTCase 4f`,
  `mult 6`, `nroots 21`): the CASSCF ran 243 macro-iterations without meeting
  the criterion; the AILFT module still completed and printed
  F0ff = 269559.5 (fixed) / F2ff = 115890.8 / Racah B = 515.1 / ZETA_F =
  2054.43 cm^-1, and ORCA then aborted with "the wavefunction IS NOT FULLY
  CONVERGED" (the file ends at the abort banner; the process was stopped
  after it).  This is the measured positive for the run-status flag: the
  parser records `converged = False, via = "wavefunction not fully
  converged"` from that text, and menu 38 marks the parameters diagnostic.
  The TRAH retry is the converged successor of the same probe.
- **nd3_freeion.{inp,out}** -- the Nd(3+) f3 free-ion probe, no `!TRAH`
  (measured engine boundary: TRAH runs abort inside the AILFT driver with
  "failed to retrieve the FAO matrix", reproduced on the Dy and Nd probes;
  without TRAH the driver runs normally).  CASSCF(3,7), mult 4, nroots 35;
  the CASSCF converged via the **energy marker** (gradient 3.4e-4 against
  the 2.5e-4 threshold) -- the measured positive for the report's
  provisional caution (`converged_via = "energy"`).  AILFT: casscf F0ff =
  219463.8 (fixed) / F2ff = 94914.9 / Racah B = 421.8; nevpt2 F2ff =
  75681.6 / B = 336.4; ZETA_F = 907.41.  This run is a source of the
  built-in free-ion reference table (`knowledge/ailft_references.py`).

## Measured format points

- The confidence-interval rows of the fit section reuse the
  Slater-Condon line shape (`F2dd = 0.000000000 a.u. = ...`); they are not
  parameters.
- Unrelated log lines also start with "Orbital" ("Orbital improvement
  steps"); the table headers are recognised by their second token
  ("Energy" for the eigenfunction table, an orbital label for VLFT).
- The SOC fit uses the CASSCF orbitals ("SPIN ORBIT COUPLING (based on
  CASSCF orbitals)"); the `-----SOC-CONSTANTS-----` summary carries
  ZETA_D or ZETA_F.
- The run-ending "the wavefunction IS NOT FULLY CONVERGED" / "Wavefunction
  not fully converged" text is a third kind of casscf verdict (alongside the
  ENERGY/GRADIENT markers and the NOT CONVERGED AFTER n CYCLES banner): it can
  appear with no convergence marker at all, so the section's converged flag
  would otherwise stay unknown.

## Regenerating

Inputs are stored verbatim; both runs need only the ORCA 6.1.1 binary.
