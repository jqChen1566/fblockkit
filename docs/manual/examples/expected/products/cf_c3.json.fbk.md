## A4 crystal-field fit

Point group C3, J = 8.0
Sampled states: 24; fitted parameters: 9
Condition number: 4.2e+05; rank 10
max |residual| = 4.466e-06 (caller's energy unit)

B_k^q:
  (2, 0)   1000
  (4,-3)   4.742679582e-10
  (4, 0)   10
  (4, 3)   50
  (6,-6)  -6.592698609e-12
  (6,-3)  -3.87506971e-11
  (6, 0)   1
  (6, 3)  -5.557741767e-13
  (6, 6)   5.339957643e-12
  const    -1.896225259e-07
note: well posed: full column rank and cond below the ill-conditioning limit.

## A5 projection-basis declaration check

Set: this set
B_k^q given: 9 (k, q) entries (axial q = 0: 3; transverse |q| >= 1: 6)

Declaration (all five items are required; a missing one is why the numbers stop
being comparable with another program's output):
  1. convention present   Stevens (Rudowicz/Ryabov lineage, cosine/sine)
  2. projection present   J = 8 manifold
  3. units      present   cm^-1
  4. z_axis     present   the three-fold axis of the model structure
  5. origin     present   Ce site
  stated: 5 of 5 items

Comparison: this set (A) vs Chilton Table S1, SINGLE_ANISO (L = 5) column (B)
  (2, 0)  comparable-roughly    A = 1000, B = 1879 (87.9% apart)
      why: (k, q) = (2, 0) is the leading axial parameter, the one class this measurement licenses for a rough comparison across projection schemes: Table S1 moved it from 1938 (SINGLE_ANISO, J = 15/2) to 1879 (SINGLE_ANISO, L = 5), 3.04%, for the same wavefunction. The two sides do not establish one scheme (not stated on both sides: convention, units, z_axis, origin; stated differently: projection). The 3% is the scale of the effect on one compound, not a bound, and it does not carry over to any other (k, q) -- in particular not to the high-order axial terms.
  (4,-3)  missing-on-one-side   A = 4.742679582e-10, B = not given
      why: (k, q) = (4, -3) is given by the first set (A) only (value 4.742679582e-10); the second set (B) does not list it, so there is nothing to compare. Next step: supply this (k, q) on both sides, or confirm that the side without it sets it to zero (a symmetry zero) before reading the two sets as one expansion.
  (4, 0)  refused               A = 10, B = 125 (1.15e+03% apart)
      why: (k, q) = (4, 0) is a high-order axial parameter (k = 4 >= 4), and high-order parameters are not comparable across projection schemes either -- only the leading axial (2, 0) is. Table S1 measured this (k, q) directly: 66 (SINGLE_ANISO, J = 15/2) vs 125 (SINGLE_ANISO, L = 5) (a factor 1.9). The two sides do not establish one declared scheme (not stated on both sides: convention, units, z_axis, origin; stated differently: projection). Next step: state all five declaration items on both sides -- if they then agree the verdict becomes 'comparable'; if they differ, compare within one scheme only, or restrict the cross-scheme reading to (2, 0) and treat it as rough.
  (4, 3)  refused               A = 50, B = -1 (102% apart, signs differ)
      why: (k, q) = (4, 3) is a transverse component (q != 0), and transverse parameters are frame and manifold dependent: two such numbers from different projection schemes can differ in sign and magnitude with neither being wrong. Table S1 measured, for the (k, q) it lists as non-zero, moves between two manifolds of one wavefunction: B_2^2 = +198 vs -43 (sign flip), B_2^-2 = 124 vs 239, B_4^0 = 66 vs 125 (a factor ~1.9), B_6^2 = 18 vs -4 (sign flip). The two sides do not establish one declared scheme (not stated on both sides: convention, units, z_axis, origin; stated differently: projection). Next step: state all five declaration items on both sides -- if they then agree the verdict becomes 'comparable'; if they differ, compare within one scheme only, or restrict the cross-scheme reading to (2, 0) and treat it as rough.
  (6,-6)  missing-on-one-side   A = -6.592698609e-12, B = not given
      why: (k, q) = (6, -6) is given by the first set (A) only (value -6.592698609e-12); the second set (B) does not list it, so there is nothing to compare. Next step: supply this (k, q) on both sides, or confirm that the side without it sets it to zero (a symmetry zero) before reading the two sets as one expansion.
  (6,-3)  missing-on-one-side   A = -3.87506971e-11, B = not given
      why: (k, q) = (6, -3) is given by the first set (A) only (value -3.87506971e-11); the second set (B) does not list it, so there is nothing to compare. Next step: supply this (k, q) on both sides, or confirm that the side without it sets it to zero (a symmetry zero) before reading the two sets as one expansion.
  (6, 0)  refused               A = 1, B = 21 (2e+03% apart)
      why: (k, q) = (6, 0) is a high-order axial parameter (k = 6 >= 4), and high-order parameters are not comparable across projection schemes either -- only the leading axial (2, 0) is. Table S1 measured, for the (k, q) it lists as non-zero, moves between two manifolds of one wavefunction: B_2^2 = +198 vs -43 (sign flip), B_2^-2 = 124 vs 239, B_4^0 = 66 vs 125 (a factor ~1.9), B_6^2 = 18 vs -4 (sign flip). The two sides do not establish one declared scheme (not stated on both sides: convention, units, z_axis, origin; stated differently: projection). Next step: state all five declaration items on both sides -- if they then agree the verdict becomes 'comparable'; if they differ, compare within one scheme only, or restrict the cross-scheme reading to (2, 0) and treat it as rough.
  (6, 3)  missing-on-one-side   A = -5.557741767e-13, B = not given
      why: (k, q) = (6, 3) is given by the first set (A) only (value -5.557741767e-13); the second set (B) does not list it, so there is nothing to compare. Next step: supply this (k, q) on both sides, or confirm that the side without it sets it to zero (a symmetry zero) before reading the two sets as one expansion.
  (6, 6)  missing-on-one-side   A = 5.339957643e-12, B = not given
      why: (k, q) = (6, 6) is given by the first set (A) only (value 5.339957643e-12); the second set (B) does not list it, so there is nothing to compare. Next step: supply this (k, q) on both sides, or confirm that the side without it sets it to zero (a symmetry zero) before reading the two sets as one expansion.
  verdicts: comparable 0, comparable-roughly 1, refused 3, missing-on-one-side 5

Rule of thumb (measured, one compound): Chilton, Chem. Soc. Rev. 2025, 54(24), 11468-11487, Table S1 -- one SA-CASSCF-SO wavefunction of one Dy(III) compound projected five ways gives different B_k^q (cm^-1):
  axial (2, 0):  1938 (SINGLE_ANISO, J = 15/2) vs 1879 (SINGLE_ANISO, L = 5), 3.04% apart -> comparable across projection schemes, roughly only;
  transverse:    (2, 2) 198 vs -43 (sign flip) and (2, -2) 124 vs 239 -> not comparable across schemes;
  high order:    (4, 0) 66 vs 125 and (6, 2) 18 vs -4 -> not comparable across schemes;
  same manifold, different program (SINGLE_ANISO, J = 15/2 vs angmom_suite, J = 15/2): (2, 0) 1938/1938, (2, 2) 198/201, (4, 0) 66/66 -> comparable.
  Rule applied here: only (2, 0) may be compared roughly when the two sides do not establish the same declared scheme; every other (k, q) needs both sides to declare the same projection scheme. Provisional: the 3% comes from one compound and one wavefunction, and it is a scale rather than a bound.

## References

Complete citations:
- [peng2025accurate] Peng, L.; Liu, S.; Zhang, X.; Chen, X.; Li, C.; Ung, S. F.; Cheng, H.; Chan, G. K. (2025). Accurate Crystal Field Hamiltonians of Single-Ion Magnets at Mean-Field Cost. The Journal of Physical Chemistry Letters, 16(47), 12312-12320. DOI: 10.1021/acs.jpclett.5c02971
- [chilton2025abinitio] Chilton, N. F. (2025). Ab initio electronic structure calculations of lanthanide single-molecule magnets; a practical guide. Chemical Society Reviews, 54(24), 11468-11487. DOI: 10.1039/d5cs00493d

BibTeX (paste-ready):

```bibtex
@article{peng2025accurate,
  author  = {Peng, Linqing and Liu, Shuanglong and Zhang, Xing and Chen, Xiao and Li, Chenghan and Ung, Shu Fay and Cheng, Hai-Ping and Chan, Garnet Kin-Lic},
  title   = {Accurate Crystal Field Hamiltonians of Single-Ion Magnets at Mean-Field Cost},
  journal = {The Journal of Physical Chemistry Letters},
  year    = {2025},
  volume  = {16},
  number  = {47},
  pages   = {12312--12320},
  doi     = {10.1021/acs.jpclett.5c02971},
}
```

```bibtex
@article{chilton2025abinitio,
  author  = {Chilton, Nicholas F.},
  title   = {Ab initio electronic structure calculations of lanthanide single-molecule magnets; a practical guide},
  journal = {Chemical Society Reviews},
  year    = {2025},
  volume  = {54},
  number  = {24},
  pages   = {11468--11487},
  doi     = {10.1039/d5cs00493d},
}
```
