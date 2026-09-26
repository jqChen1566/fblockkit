# fBlockKit report

## Summary

Subject: work\generated_ce3_sarc2.out
Counts: Warning 2

## A1 orbital composition (per MO, Loewdin)

Composition source: LOEWDIN ORBITAL-COMPOSITIONS table (176 orbitals in total)
Occupation partition: occ > 0.02 counts as occupied (28), otherwise as a virtual orbital (148)

f composition ranking (occupied partition, top 8):
  (no f-composition orbital in the occupied partition: the highest f weight is only 0.0% -- if this system should have electrons in that shell, check the initial guess / active-orbital window, it may be the wrong branch of a multiple-solution problem)

Dominant shells of the active orbitals (0.02 < occ < 1.98):
  MO   27  occ 1.0000  Ce1 d 100.0%

Criterion and convention: the assignment is the (atom, shell) with the largest weight; every ranking is partitioned by occupation first (measured on this group's EuF: ranking the whole space by weight picks up virtual orbitals).

## A2 occupation spectrum and entropy bound spectrum

Active space: 7 orbitals; occupations (input order): 1.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000

Single-orbital entropy BOUND spectrum (descending; s_bound = -2[x ln x + (1-x) ln(1-x)], x = n/2 -- the maximum-entropy completion of the four occupation weights for a spin-summed occupation):
  1.3863 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
  max(s_bound) = 1.3863
  Exclusion test: max(s_bound) = 1.3863 > 0.14 -> NOT decidable from the occupations alone. The bound is loose (it is attained only when the two spin channels of the orbital are uncorrelated), so it is not evidence of multireference character; the true four-state entropy needs the 2-RDM route (orca_2json RDM2_aa/ab/bb plus orca_loc orbitals; see the A6 section).
  Plateau/cliff (ranking proxy): the largest relative gap (100.0%) sits at item 1 of the descending order -- the first 1 items form the plateau candidate (active-orbital indices, input order from 0: [0]). The proxy ranks orbitals by the bound; relative order is not guaranteed to match the true four-state entropies.

Threshold-based candidates (s_bound > 0.14, active-orbital indices): [0]

Preconditions (must travel with the conclusion):
- The entropy spectrum must be read in a localized-orbital basis (autoCAS); a "plateau" in a non-localized or truncated space can be numerical noise.
- A closed-shell single determinant gives bound 0 for every orbital; for an open-shell occupation n = 1 the bound is ln 4 = 1.386294 by construction (a spin-restricted description of one unpaired electron), so an f1 active space always reads as "not decidable", never as multireference.
- The 0.14 line applies to the true four-state entropy; the threshold 0.14 and Zs(1) > 0.2 mean different things and must not be mixed.

## A3 multi-reference character

One panel over the multi-reference (MR) character indicators this output carries; every line is labelled with its source. The verdict below is drawn from the available indicators only.

Indicator table (3 indicators; line = the published criterion, n/a = not available in this output):
  indicator                         value                   line                verdict
  T1 diagnostic                     n/a                     0.02                unavailable
  fractional active orbitals        0 of 7 in 0.02-1.98     window 0.02-1.98    none
  max single-orbital entropy bound  1.386294                0.14 (exclusion)    above: indecisive

Per-indicator verdicts:
- T1 diagnostic: not available: T1 needs a coupled-cluster (CCSD) output: ORCA prints it inside the "COUPLED CLUSTER ENERGY" block. Next step: run a CCSD single point at the same geometry and rerun this panel on that output.
- fractional active orbitals: no active orbital has a fractional occupation (all 7 occupations are integer-valued: 0, 1 or 2) -> this indicator gives no evidence of multi-reference character. The window is ORCA's active-space convention (auto-ICE); the count itself has no published threshold. An occupation of exactly 1 (a single unpaired electron) is not fractional.
- max single-orbital entropy bound: max(s_bound) = 1.386294 > 0.14: indecisive -- the bound is above the line, which is not evidence of multireference character. The true four-state entropy needs the 2-RDM route (orca_2json RDM2_aa/ab/bb plus orca_loc orbitals; see the A6 section). Convention: s_bound is the maximum-entropy completion of the four occupation weights for a spin-summed occupation (A2 module), so it never under-estimates the true four-state entropy; it is tight only for uncorrelated spin channels. An orbital with occupation 1 gives s_bound = ln 4 = 1.386294 by construction. The 0.14 line (Stein & Reiher 2016) applies to the true four-state entropy; Wardzala et al. 2026 quotes the same number as the M-diagnostic line. The Zs(1) > 0.2 line means something else again and is not used here.

Combined verdict: MR character not indicated
  fired: (none)
  not fired: fractional active orbitals (0 of 7 fractional); max single-orbital entropy bound (1.386294 > 0.14 (indecisive))
  unavailable: T1 diagnostic (needs a CCSD output)
Panel completeness: 2 of 3 indicators available; the verdict rests on the available indicator(s) only and the missing one(s) are named here, never replaced silently.

Caveats (carried with the verdict):
- Threshold calibration range: the T1 0.02 line and the max(s1) 0.14 line were calibrated on main-group and 3d systems; the review that quotes the 0.14 M-diagnostic line contains no lanthanide or actinide example at all (Wardzala et al. 2026). For an f-block system these are a default starting value that must be verified, not a verdict.
- Two diagnostics share the number 0.14: max(s1) > 0.14 (Stein & Reiher 2016) and the M diagnostic > 0.14 (Wardzala et al. 2026). The Zs(1) > 0.2 line means something else again and is not used by this panel.
- Entropy bound convention: s_bound is the maximum-entropy completion of the four occupation weights for a spin-summed occupation (A2 module); it never under-estimates the true four-state entropy and is used here as an exclusion test only. An orbital with occupation 1 gives s_bound = ln 4 = 1.386294 by construction (measured on the Ce3+ fixture: N(occ) = (1, 0, 0, 0, 0, 0, 0), 0 fractional orbitals, max(s_bound) = 1.386294 -> indecisive). The true four-state entropy needs the 2-RDM route (orca_2json RDM2 + orca_loc).
- The verdict rests on the available indicators only; the ones missing from this output are named above and are not replaced by another indicator.

## Findings

### 1. [Warning] Active space contains a nearly empty orbital (occupation < 0.02)

- Suggested action: Check the active-orbital window: an active orbital with an occupation below 0.02 does not really take part in the correlation, so the space may be too large or the window poorly chosen; consider narrowing the active space by occupation and rerunning, and confirm whether that orbital was meant to be included at all.
- Rule: DG-ACTIVE-SPACE-NEAR-EMPTY
- Evidence:
  - [Manual] "All orbitals between occupation number say 1.98 down to 0.02 will be included in the active space." (auto-ICE's automatic window convention; this rule borrows the 0.02/1.98 interval as a check-up threshold)
    - Source: ORCA 6.1 manual §3.14 (ICE-CI and auto-ICE) (https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/iceci.html)

### 2. [Warning] No active orbital carries f character on an f-block system (a d-type solution branch?)

- Suggested action: Check the initial guess and the active-orbital window before using anything from this f-block calculation: the active orbitals carry essentially no f character, which is the signature of an SCF/CASSCF that found a d-type solution branch instead of the f^n branch. Measured examples: a Ce3+ CASSCF(1,7) run whose singly occupied active orbital came out 100% Ce-d with all seven 4f orbitals virtual (largest f weight over the active orbitals 0.0%), and the multi-solution SCF/CASSCF records of this group (EuF branches differ by 33 kcal/mol in De). Next steps, in the manual's order: give the SCF better starting orbitals (a fragment guess, or the target orbitals by hand), revise the active-orbital window, or restart from several starting points and keep the branch that matches the target configuration -- then rerun this check-up.
- Rule: DG-ACTIVE-ORBITALS-WITHOUT-F-CHARACTER
- Evidence:
  - [Measured] Both Ce3+ fixtures (ce3_orbcomp.out and the recipe-generated generated_ce3_sarc2.out, ORCA 6.1.1) landed on the d1 solution: the singly occupied active orbital (MO 27, occupation 1.0) is 100% Ce-d, every 4f orbital is virtual, and the largest f weight over the active orbitals is 0.0%.
    - Source: Fixtures fixtures/orca/ce3_orbcomp.out and generated_ce3_sarc2.out (the A1 note in fixtures/orca/README.md)
  - [Measured] This group's multi-solution records: f-block SCF/CASSCF runs converge to different solution branches from different starting guesses (EuF: two branches 33 kcal/mol apart in De; EuO: the branch judgement has to be made at the method level), so the branch actually found must be checked, not assumed.
    - Source: Group's records (EuF/EuO multi-solution notes, B1_可微分DFT_cjq6)
  - [Measured] Threshold choice (ours, provisional): the flagged case sits at 0.0% while a genuine f window would show ~100% (the Ce fixture's 4f orbitals are 100% Ce-f, merely virtual), so 10% separates the two sides with a wide margin; f/d-mixed windows are handled correctly because the maximum is taken over all active orbitals (one f-carrying orbital is enough to keep the check silent).
    - Source: measured on the fixtures above (the threshold is this module's own choice)

## Provenance

- [Manual] ORCA 6.1 manual §3.14 (ICE-CI and auto-ICE)
- [Measured] Fixtures fixtures/orca/ce3_orbcomp.out and generated_ce3_sarc2.out (the A1 note in fixtures/orca/README.md)
- [Measured] Group's records (EuF/EuO multi-solution notes, B1_可微分DFT_cjq6)
- [Measured] measured on the fixtures above (the threshold is this module's own choice)
- [Measured] This group's differentiable multireference stack records: composition-split-needs-two-shells / rank-within-partition-not-across (2026-09)
- [Literature] Wardzala J. J. et al., Chem. Rev., 2026, 126(8), 4592-4618, DOI 10.1021/acs.chemrev.5c00866 (Eq. (6))
- [Literature] Stein C. J., Reiher M., J. Chem. Theory Comput., 2016, 12(4), 1760-1771, DOI 10.1021/acs.jctc.6b00156
- [Literature] Stein C. J., Reiher M., Mol. Phys., 2017, 115(17-18), 2110-2119, DOI 10.1080/00268976.2017.1288934
- [Literature] Stein C. J., Reiher M., J. Comput. Chem., 2019, 40(25), 2216-2226, DOI 10.1002/jcc.25869
- [Measured] B1_可微分DFT_cjq6/t177熵谱失效_文献定案_20260922.md; bound property: tests/test_analysis.py
- [Literature] Lee T. J., Taylor P. R., Int. J. Quantum Chem., 1989, 36(S23), 199-207, DOI 10.1002/qua.560360824
- [Literature] Stein C. J., Reiher M., J. Chem. Theory Comput., 2016, 12(4), 1760-1771, DOI 10.1021/acs.jctc.6b00156
- [Literature] Wardzala J. J. et al., Chem. Rev., 2026, 126(8), 4592-4618, DOI 10.1021/acs.chemrev.5c00866 (Sec. 2.3.2)
- [Literature] Wardzala J. J. et al., Chem. Rev., 2026, 126(8), 4592-4618, DOI 10.1021/acs.chemrev.5c00866; scope noted in this group's reading record 文献细读/细读_自动活性空间综述_ChemRev2026.md (2026-09-25, fact 14)
- [Manual] ORCA 6.1 manual Sec. 3.14 (ICE-CI and auto-ICE)
- [Measured] fixtures/orca/ (ORCA 6.1.1 outputs; the D1 measurement is also recorded in src/fblockkit/parsers/orca.py)

## References

Complete citations:
- [wardzala2026multireference] Wardzala, J. J.; Hennefarth, M. R.; Agarawal, V.; Jangid, B.; Seal, A.; Hermes, M. R.; King, D. S.; Gagliardi, L. (2026). Multireference Methods for Chemistry and Materials Science: Automated Active Spaces, Efficient Dynamic Correlation, and Extended Systems. Chemical Reviews, 126(8), 4592-4618. DOI: 10.1021/acs.chemrev.5c00866
- [stein2016automated] Stein, C. J.; Reiher, M. (2016). Automated Selection of Active Orbital Spaces. Journal of Chemical Theory and Computation, 12(4), 1760-1771. DOI: 10.1021/acs.jctc.6b00156
- [stein2017measuring] Stein, C. J.; Reiher, M. (2017). Measuring multi-configurational character by orbital entanglement. Molecular Physics, 115(17-18), 2110-2119. DOI: 10.1080/00268976.2017.1288934
- [stein2019autocas] Stein, C. J.; Reiher, M. (2019). autoCAS: A Program for Fully Automated Multiconfigurational Calculations. Journal of Computational Chemistry, 40(25), 2216-2226. DOI: 10.1002/jcc.25869
- [lee1989diagnostic] Lee, T. J.; Taylor, P. R. (1989). A diagnostic for determining the quality of single-reference electron correlation methods. International Journal of Quantum Chemistry, 36(S23), 199-207. DOI: 10.1002/qua.560360824

BibTeX (paste-ready):

```bibtex
@article{wardzala2026multireference,
  author  = {Wardzala, Jacob J. and Hennefarth, Matthew R. and Agarawal, Valay and Jangid, Bhavnesh and Seal, Aniruddha and Hermes, Matthew R. and King, Daniel S. and Gagliardi, Laura},
  title   = {Multireference Methods for Chemistry and Materials Science: Automated Active Spaces, Efficient Dynamic Correlation, and Extended Systems},
  journal = {Chemical Reviews},
  year    = {2026},
  volume  = {126},
  number  = {8},
  pages   = {4592--4618},
  doi     = {10.1021/acs.chemrev.5c00866},
}
```

```bibtex
@article{stein2016automated,
  author  = {Stein, Christopher J. and Reiher, Markus},
  title   = {Automated Selection of Active Orbital Spaces},
  journal = {Journal of Chemical Theory and Computation},
  year    = {2016},
  volume  = {12},
  number  = {4},
  pages   = {1760--1771},
  doi     = {10.1021/acs.jctc.6b00156},
}
```

```bibtex
@article{stein2017measuring,
  author  = {Stein, Christopher J. and Reiher, Markus},
  title   = {Measuring multi-configurational character by orbital entanglement},
  journal = {Molecular Physics},
  year    = {2017},
  volume  = {115},
  number  = {17-18},
  pages   = {2110--2119},
  doi     = {10.1080/00268976.2017.1288934},
}
```

```bibtex
@article{stein2019autocas,
  author  = {Stein, Christopher J. and Reiher, Markus},
  title   = {autoCAS: A Program for Fully Automated Multiconfigurational Calculations},
  journal = {Journal of Computational Chemistry},
  year    = {2019},
  volume  = {40},
  number  = {25},
  pages   = {2216--2226},
  doi     = {10.1002/jcc.25869},
}
```

```bibtex
@article{lee1989diagnostic,
  author  = {Lee, Timothy J. and Taylor, Peter R.},
  title   = {A diagnostic for determining the quality of single-reference electron correlation methods},
  journal = {International Journal of Quantum Chemistry},
  year    = {1989},
  volume  = {36},
  number  = {S23},
  pages   = {199--207},
  doi     = {10.1002/qua.560360824},
  note    = {Year of record 1989 (Symposium issue); the Crossref record carries a digitisation date of 2009.},
}
```
