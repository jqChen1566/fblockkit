# fBlockKit report

## Summary

Subject: work\n2_stretch_local_spin.out
Counts: no findings

## A6 local spin analysis

Local spin source: the "LOCAL SPIN ANALYSIS (Loewdin* projector)" block printed by ORCA's SCF / CASSCF module (the projector is named in the block header, and the values depend on it).
Blocks found in this output: 3.
Fragments declared in the input: 2
  1: N
  2: N
Block selection: this output carries no state-resolved block, so it is read as an SCF job. ORCA prints the analysis for the initial guess, after SCF and once more for the final density; only the LAST block describes the converged wavefunction, so block 3 of 3 is used and the earlier ones are left out.
Unused earlier block(s), for comparison only: block 1 (Seff 1.5231, 1.5231); block 2 (Seff 1.4699, 1.4699). Their numbers are not the converged wavefunction's -- reading them instead of the last block changes the values.

--- SCF final wavefunction (block 3 of 3) ---
Block scope: 2 atoms, 28 basis functions.
Per-fragment local spin (fragment numbers as printed by ORCA):
  fragment  elements                    <SzA>   Seff(A)
         1  N                         +1.4371    1.4699
         2  N                         -1.4371    1.4699
<SA*SB> matrix (row = fragment A, column = fragment B; the ORCA table is printed lower-triangular, this is the full symmetric matrix):
     A\B       B=1       B=2
     A=1    3.6304   -2.1835
     A=2   -2.1835    3.6304
Interpretation:
- f-block fragment(s): none -- no declared fragment contains a lanthanide or actinide element, so this is not a metal/ligand decomposition. The tables above are printed as they stand and no metal spin-share is reported.
- Local-spin magnitude: Seff = 1.4699 (fragment 1, N), 1.4699 (fragment 2, N).
- Total <S^2> = sum over all fragment pairs (A, B) of <SA*SB> = 2.8938.
- Spin purity: total <S^2> = 2.8938; this block declares no multiplicity (an SCF local spin block prints none), so no S(S+1) is declared to compare against. 2.8938 is not S(S+1) of any spin eigenstate (the values run 0, 0.75, 2, 3.75, 6, ...); the printed <SzA> values give M_S = +0.0000, so this solution is spin-contaminated.
  Against the smallest spin the printed M_S = +0.0000 admits -- a singlet, S = 0, S(S+1) = 0.0000 -- the deviation is +2.8938. A spin-unrestricted single determinant is not in general a spin eigenstate, so a deviation here is a property of the method, not a defect of the file.
  This output carries a broken-symmetry warning, and a broken-symmetry solution is deliberately not spin-pure: the contamination is the intended non-eigenstate behaviour of that method, not a defect to repair.
- Spin-share (this module's own construction, not from the cited sources -- provisional):
    share(M) = (sum over all fragments B of <S_M*S_B>) / <S^2>_total
    share(fragment 1, N) = 1.4469 / 2.8938 = 0.5000
    share(fragment 2, N) = 1.4469 / 2.8938 = 0.5000
    The shares sum to 1.0000 by construction (the rows of <SA*SB> partition the total <S^2>).
- Sum of <SzA> = +0.0000, the values being +1.4371 / -1.4371: the local projections add up to the total M_S of the state, so a zero sum is exact whenever the state has as many alpha as beta electrons -- the fragments carry equal and opposite local spin even though M_S = 0 (a spin-polarised density in a state with no net spin).
- Printed-number check (this module's own, provisional): the diagonal <SA*SA> equals Seff(A)(Seff(A)+1) in every fragment, the largest deviation being 1.06e-04 at fragment 1 (limit 3e-03), so the printed columns are mutually consistent.

Boundary of this analysis (travels with the numbers above):
- Local spin analysis is NOT the environment spin polarisation entropy Delta S_E of Eq. (9). It partitions the spin over fragments; Delta S_E measures how far the environment density is from the spin-unpolarised point D_E^alpha = D_E^beta = D_E/2, and it needs the two-step off-line route orca_2json (RDM2_aa / RDM2_ab / RDM2_bb; the 2-RDMs are not stored inside the .densities file) plus orca_loc (IAO-IBO localised orbitals, .loc.gbw). No printed ORCA table carries the spin-resolved environment density.
    Eq. (9): Delta S_E = -2 Tr[(D_E/2) ln(D_E/2)] + Tr[D_E^alpha ln D_E^alpha] + Tr[D_E^beta ln D_E^beta], with D_E = D_E^alpha + D_E^beta.
- Because this block is a surrogate, its discriminating power is limited: an unpolarised environment and a strongly polarised one both give finite fragment Seff values here, and nothing in this report separates them on the scale below.
- The scale this surrogate stands in for (measured in the reference paper on the Dy single-ion magnet 1Dy, 6H, 11 states): the wrong SCF solution gave Delta S_E = 2.766 against 0.007 for the correct one, and the wrong solution degraded the downstream pre-SOC magnetic MAE from 8.6 to 357.7 cm^-1. Treat this report as a first look at the spin partition, not as that criterion.
- The numbers depend on two user choices: the fragment definition (which atoms were assigned to which fragment) and the projector named in the block header (this project's fixtures use 'Loewdin*'). Values from different choices are not comparable.

## Findings

No diagnostic rule matched.

## Provenance

- [Manual] ORCA 6.1 manual, section 5.1.10 'Local Spin Analysis' (read 2026-09-26): https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/population.html
- [Literature] Ai Y., Li Z.-W., Guan Z.-B., Jiang H., J. Chem. Theory Comput., 2025, 21(19), 9631-9640, DOI 10.1021/acs.jctc.5c01336
- [Literature] Herrmann C., Reiher M., Hess B. A., J. Chem. Phys., 2005, 122(3), 034102, DOI 10.1063/1.1829050
- [Measured] reproduced by tests/test_local_spin.py (2026-09-26 fixtures, ORCA 6.1.1)

## References

Complete citations:
- [ai2025density] Ai, Y.; Li, Z.; Guan, Z.; Jiang, H. (2025). Density Matrix Embedding Theory-Based Multiconfigurational Quantum Chemistry Approach to Lanthanide Single-Ion Magnets. Journal of Chemical Theory and Computation, 21(19), 9631-9640. DOI: 10.1021/acs.jctc.5c01336
- [herrmann2005comparative] Herrmann, C.; Reiher, M.; Hess, B. A. (2005). Comparative analysis of local spin definitions. The Journal of Chemical Physics, 122(3), 034102. DOI: 10.1063/1.1829050

BibTeX (paste-ready):

```bibtex
@article{ai2025density,
  author  = {Ai, Yuhang and Li, Ze-Wei and Guan, Zhe-Bin and Jiang, Hong},
  title   = {Density Matrix Embedding Theory-Based Multiconfigurational Quantum Chemistry Approach to Lanthanide Single-Ion Magnets},
  journal = {Journal of Chemical Theory and Computation},
  year    = {2025},
  volume  = {21},
  number  = {19},
  pages   = {9631--9640},
  doi     = {10.1021/acs.jctc.5c01336},
}
```

```bibtex
@article{herrmann2005comparative,
  author  = {Herrmann, Carmen and Reiher, Markus and Hess, Bernd Artur},
  title   = {Comparative analysis of local spin definitions},
  journal = {The Journal of Chemical Physics},
  year    = {2005},
  volume  = {122},
  number  = {3},
  pages   = {034102},
  doi     = {10.1063/1.1829050},
}
```
