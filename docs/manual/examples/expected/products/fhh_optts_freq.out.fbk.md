# fBlockKit report

## Summary

Subject: work\fhh_optts_freq.out
Counts: Warning 1; Info 2

## Findings

### 1. [Warning] Too few SCF cycles is treated as a suspected pseudo-convergence

- Suggested action: Counter-check with SlowConv, and check the N of "SCF CONVERGED AFTER N CYCLES".
- Rule: DG-SCF-PSEUDO-CONVERGENCE
- Evidence:
  - [Measured] For a transition-metal complex ORCA can report "convergence" after only 2 SCF cycles (a pseudo-convergence signal).
    - Source: SCINE-stack project notes pit table row 35 ("ORCA may report SCF convergence after only 2 cycles for Pd complexes"); restated as item 11, section 8, of SCINE-stack/SCINE 能力评估报告.md

### 2. [Info] Exactly one imaginary mode -- a first-order saddle, i.e. a transition-state candidate

- Suggested action: One imaginary mode is the signature of a first-order saddle point, so this structure is a transition-state candidate. Before it is used as a transition state, confirm that the imaginary mode is the mode of interest and that it connects the intended minima -- inspect the normal mode or follow the intrinsic reaction coordinate and check which minima it links, which is what the tutorial asks for; confirm that the whole study stayed at one level of theory, because a saddle point of one surface is not a stationary point of another; and confirm that the count came from the last frequency block of the job rather than from the initial Hessian block an OptTS run prints first. If the mode is small, treat it as unverified and check it by the independent route of the small-imaginary-mode rule (D4-FREQ-SMALL-IMAGINARY-MODE).
- Rule: D4-TS-ONE-IMAGINARY-MODE
- Evidence:
  - [Manual] "However, it is further important to check if the optimized structure is a minimum or a saddle point (TS) on the potential energy surface. The latter is characterized by exactly one imaginary frequency in the frequency calculations on the finally optimized TS geometry. This imaginary frequency must further represent the structural transition of interest ..."
    - Source: ORCA 6.1 tutorial, Reaction path / TS optimization, TS Verification (https://www.faccts.de/docs/orca/6.1/tutorials/react/tsopt.html)
  - [Measured] Single-level rule, measured on the Pd(PH3)2 + CH3I oxidative addition: giving the GFN2 transition-state geometry to r2SCAN-3c gave 7 imaginary modes on the DFT surface (it is not a stationary point there at all), while the same-level r2SCAN-3c verification run gave exactly 1 imaginary mode (at -240.56 cm^-1) together with three minima reporting 0 imaginary modes.
    - Source: SCINE-stack/SCINE capability assessment report, section 4.1 and pitfall list item 9 (single level of theory throughout) of D:/project/SCINE-stack/SCINE 能力评估报告.md; the same record is restated as pit table row 32 of D:/project/SCINE-stack (project notes)

### 3. [Info] A transition-state-like species was found -- a barrier derived from it must state its reference point

- Suggested action: When a barrier derived from this species is reported, state explicitly what the energies are measured from (separated reactants, pre-complex, the endpoint of a path profile, ...), and check that the transition state and the reference structures were computed at the same level of theory. Measured: the same barrier data can be read as +5.90 kcal/mol (a barrier) or -22.14 kcal/mol (no barrier) depending on the reference point, because the routine that produced the number measured from a path endpoint that was a pre-complex rather than the separated reactants.
- Rule: D4-BARRIER-REFERENCE-POINT
- Evidence:
  - [Measured] The reference point of a barrier: get_barrier_from_spline() returns the barrier relative to the two spline endpoints, and an endpoint can be a pre-complex rather than the separated reactants, so the same data set can give the two opposite readings "+5.90 kcal/mol, there is a barrier" and "-22.14 kcal/mol, there is no barrier"; whenever this number is quoted, the reference point must be stated. The same source records the effect at the DFT level as well: the r2SCAN-3c barrier of that step is +12.30 kcal/mol relative to the pre-complex but +6.54 kcal/mol relative to the separated reactants. (Translated from the source, which is written in Chinese.)
    - Source: SCINE-stack/SCINE capability assessment report, section 4.1 and pitfall list item 8 ("a barrier must state its reference point"): D:/project/SCINE-stack/SCINE 能力评估报告.md

## Provenance

- [Measured] SCINE-stack project notes pit table row 35 ("ORCA may report SCF convergence after only 2 cycles for Pd complexes"); restated as item 11, section 8, of SCINE-stack/SCINE 能力评估报告.md
- [Manual] ORCA 6.1 tutorial, Reaction path / TS optimization, TS Verification
- [Measured] SCINE-stack/SCINE capability assessment report, section 4.1 and pitfall list item 9 (single level of theory throughout) of D:/project/SCINE-stack/SCINE 能力评估报告.md; the same record is restated as pit table row 32 of D:/project/SCINE-stack (project notes)
- [Measured] SCINE-stack/SCINE capability assessment report, section 4.1 and pitfall list item 8 ("a barrier must state its reference point"): D:/project/SCINE-stack/SCINE 能力评估报告.md
