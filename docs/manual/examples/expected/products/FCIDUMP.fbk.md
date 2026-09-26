## A2x exact four-state single-orbital entropy (FCIDUMP route)

Active space: 6 orbitals / 6 electrons (M_s sector with n_alpha=3, n_beta=3); determinants solved: 400; selected root 0 with <S^2> = 0.0000

Cross-checks against the engine (the reconstruction is trusted only with them):
  - CI energy + core energy = -108.950671945279 Eh vs the run's -108.950671945000 Eh (deviation 2.79e-10; the solver refuses a match beyond 1e-6)
  - natural occupations vs the printed N(occ): max |deviation| = 4.56e-06 (the print precision is 1e-5)
  - localization rotation orthogonality: max |u^T u - I| = 7.41e-14 (two exports of the same gbw, one localized window)

Four-state single-orbital entropy (canonical orbital basis of the dump; informational -- the 0.14 line and the plateau reading are defined in a localized basis):
  0.0223 0.2258 0.2258 0.2302 0.2302 0.0078

Four-state single-orbital entropy (localized basis: IAO-IBO rotation of the active window [4, 5, 6, 7, 8, 9] from the user's orca_loc):
  1.3547 1.3547 1.3547 1.3547 0.3732 0.3885
  max s1 = 1.3547 vs the Stein & Reiher line 0.14: candidate orbitals (original numbering): [4, 5, 6, 7, 8, 9]
  Plateau/cliff: the largest relative gap (71.3%) sits after item 4 of the descending order.

Four weights per localized orbital (empty, up, down, double):
  orbital 4: 0.1869 0.3125 0.3125 0.1881
  orbital 5: 0.1869 0.3125 0.3125 0.1881
  orbital 6: 0.1869 0.3125 0.3125 0.1881
  orbital 7: 0.1869 0.3125 0.3125 0.1881
  orbital 8: 0.9078 0.0459 0.0459 0.0003
  orbital 9: 0.0026 0.0459 0.0459 0.9055

Preconditions carried with these numbers:
- The four-state density is diagonal for a fixed-particle-number state (an exact statement, not an approximation), so the three numbers n_up, n_down and p2 are complete -- and they exist only in a reconstruction like this one; a printed N(occ) line cannot supply them.
- Read the threshold against the localized-basis spectrum only; the canonical spectrum is reported for orientation.
- The solver is exact within the dumped active space; the usual active-space caveats (size, composition, solution branch) transfer unchanged.

## Environment spin-polarisation entropy (Delta S_E)

Partition of the active orbitals (largest Löwdin population per centre; the centres treated as the cluster: 0):
  orbital 0: centre 0 (population 0.9517) -> cluster
  orbital 1: centre 0 (population 0.9517) -> cluster
  orbital 2: centre 1 (population 0.9517) -> environment
  orbital 3: centre 1 (population 0.9517) -> environment
  orbital 4: centre 1 (population 0.5000) -> environment
  orbital 5: centre 1 (population 0.5000) -> environment

Environment: 4 active orbital(s) [2, 3, 4, 5]; electrons in the environment block: 1.998836 alpha + 1.998836 beta
  eigenvalues of D_E (alpha+beta): 0.001897 1.001164 1.001164 1.993449
  Delta S_E = -2 Tr[(D/2) ln(D/2)] + Tr[D_a ln D_a] + Tr[D_b ln D_b]
            = 1.405324 - 0.702662 - 0.702662 = 0.000000

Reading (provisional line; the source's anchors are 0.007 for its correct solution and 2.766 for its wrong one, and the criterion's zero is exact):
  - the environment carries no appreciable spin polarisation (Delta S_E = 0.000000)

Boundary of this route: the inactive orbitals of the CASSCF wave function are doubly occupied and carry no spin by construction, so this is the *active* environment. A mean-field solution's ligand spin polarisation is diagnosed by the local-spin analysis of the SCF output (menu 1) instead.

## References

Complete citations:
- [stein2016automated] Stein, C. J.; Reiher, M. (2016). Automated Selection of Active Orbital Spaces. Journal of Chemical Theory and Computation, 12(4), 1760-1771. DOI: 10.1021/acs.jctc.6b00156
- [ai2025density] Ai, Y.; Li, Z.; Guan, Z.; Jiang, H. (2025). Density Matrix Embedding Theory-Based Multiconfigurational Quantum Chemistry Approach to Lanthanide Single-Ion Magnets. Journal of Chemical Theory and Computation, 21(19), 9631-9640. DOI: 10.1021/acs.jctc.5c01336

BibTeX (paste-ready):

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
