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

## Regenerating

Inputs are stored verbatim; both runs need only the ORCA 6.1.1 binary.
