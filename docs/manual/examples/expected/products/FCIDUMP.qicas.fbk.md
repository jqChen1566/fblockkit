## 3.4 QICAS active-space optimization

QICAS orbital optimization (out-of-CAS correlation over orbital rotations):
  window: 6 orbitals, 6 electrons; target (2e, 4o) with 2 closed, 4 active, 0 virtual
  requested partition (by natural occupation): closed [0, 2], active [1, 3, 4, 5], virtual []

  orbital  occupation  S initial  S optimized
        0     1.99345     0.0223       0.0078
        1     1.93634     0.2258       0.2258
        2     1.93634     0.2258       0.0223
        3     0.06599     0.2302       0.2302
        4     0.06599     0.2302       0.2302
        5     0.00190     0.0078       0.2258

  F_QI (non-active entropies): 0.2481 -> 0.0302  (2 rotations accepted in 2 cycle(s), set 'touch', eps1 = 1e-08)

Energy checks (window: exact CI; CASCI in the partition):
  E(window FCI)          = -108.950671945 Eh
  E(CASCI, initial)      = -108.883698113 Eh  (+6.697383e-02 vs FCI)
  E(CASCI, optimized)    = -108.945340635 Eh  (+5.331310e-03 vs FCI)
  Theorem-1 check: gap 5.331e-03 <= dE_max/ln4 * F_QI = 1.034e-01 (holds)
  final partition (occupancy reading of the optimized basis): closed [2], active [1, 3, 4, 5], virtual [0]  [the optimizer re-shuffled the non-active slots: the final space reads as (4e, 4o)]

Boundaries and checks:
  - the RDMs are the *exact* ground state of the FCIDUMP window (the four-state entropy route); the source drives QICAS with a low-bond-dimension DMRG ground state over the full space -- this is the source's own subset application, so F_QI here sums over the window's non-active orbitals only and is not comparable with the source's full-space values
  - the partition is this tool's documented ordering choice (the dumped orbitals re-sorted by occupation); the source fixes closed/active/virtual as index ranges of an ordered basis
  - the optimized rotation is reported as a matrix and used for the CASCI check; writing it back into a .gbw (the mkl route of menu 18) is not implemented
  - the source's size-selection variant (minimize the total orbital entropy and read the plateau of the threshold diagram; their Appendix C) is not implemented
  - the run's printed CASSCF energy cross-check: applied (the CI root was matched to it)

## References

Complete citations:
- [ding2023qicas] Ding, L.; Knecht, S.; Schilling, C. (2023). Quantum Information-Assisted Complete Active Space Optimization (QICAS). The Journal of Physical Chemistry Letters, 14(49), 11022-11029. DOI: 10.1021/acs.jpclett.3c02536

BibTeX (paste-ready):

```bibtex
@article{ding2023qicas,
  author  = {Ding, Lexin and Knecht, Stefan and Schilling, Christian},
  title   = {Quantum Information-Assisted Complete Active Space Optimization (QICAS)},
  journal = {The Journal of Physical Chemistry Letters},
  year    = {2023},
  volume  = {14},
  number  = {49},
  pages   = {11022--11029},
  doi     = {10.1021/acs.jpclett.3c02536},
}
```
