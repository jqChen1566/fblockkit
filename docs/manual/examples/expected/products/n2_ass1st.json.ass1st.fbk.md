## ASS1ST selection round (NEVPT2-density quasi-NOONs)

ASS1ST selection round (NEVPT2 quasi-natural occupation numbers):
  system: n2_ass1st; 6 active orbitals (6e) against 4 internal and 18 external; single state
  band: external line 0.03, internal line 1.97 (the source's quasi-NOON window; 2 - T for a single threshold T)

  internal block    occupation  in band
                       2.00000  
                       2.00000  
                       1.98196  
                       1.97680  
  external block    occupation  in band
                       0.01566  
                       0.00818  
                       0.00531  
                       0.00492  
                       0.00492  
                       0.00377  
                       0.00377  
                       0.00344  
                       0.00344  
                       0.00153  
                       0.00153  
                       0.00148  
                       0.00148  
                       0.00136  
                       0.00087  
                       0.00070  
                       0.00070  
                       0.00012  

Active orbitals (CASSCF natural occupations of this round):
  MO   4  occupation 1.99345
  MO   5  occupation 1.93634
  MO   6  occupation 1.93634
  MO   7  occupation 0.06599
  MO   8  occupation 0.06599
  MO   9  occupation 0.00190

Suggestion:
  reassign active to internal (occupation at the top of the band): 4
  reassign active to external (occupation at the bottom of the band): 9
  next space: (4e, 4o) -- from (6e, 6o)

Boundaries and checks:
  - the density is ORCA's FIC-NEVPT2 unrelaxed density; the source works with the SC-NEVPT2 first-order density -- same object family (the 0+1-order density), different contraction, so counts near the thresholds may shift (carried as an honest substitution; no literature comparison exists)
  - the source's own caveats: the outcome depends on the initial space (different sensible initials converge together only with conservative thresholds), and growing plus shrinking can cycle between two spaces -- break ties by chemical judgment
  - every block table is this tool's own block diagonalization of the exported density (the source's construction); ORCA's printed NOON list is the whole-space naturalization and is not used
  - multiplicity 1; the next round's input is generated next to this report by the menu

## References

Complete citations:
- [khedkar2019ass1st] Khedkar, A.; Roemelt, M. (2019). Active Space Selection Based on Natural Orbital Occupation Numbers from n-Electron Valence Perturbation Theory. Journal of Chemical Theory and Computation, 15(6), 3522-3536. DOI: 10.1021/acs.jctc.8b01293
- [khedkar2020sa] Khedkar, A.; Roemelt, M. (2020). Extending the ASS1ST Active Space Selection Scheme to Large Molecules and Excited States. Journal of Chemical Theory and Computation, 16(8), 4993-5005. DOI: 10.1021/acs.jctc.0c00332

BibTeX (paste-ready):

```bibtex
@article{khedkar2019ass1st,
  author  = {Khedkar, Abhishek and Roemelt, Michael},
  title   = {Active Space Selection Based on Natural Orbital Occupation Numbers from n-Electron Valence Perturbation Theory},
  journal = {Journal of Chemical Theory and Computation},
  year    = {2019},
  volume  = {15},
  number  = {6},
  pages   = {3522--3536},
  doi     = {10.1021/acs.jctc.8b01293},
}
```

```bibtex
@article{khedkar2020sa,
  author  = {Khedkar, Abhishek and Roemelt, Michael},
  title   = {Extending the ASS1ST Active Space Selection Scheme to Large Molecules and Excited States},
  journal = {Journal of Chemical Theory and Computation},
  year    = {2020},
  volume  = {16},
  number  = {8},
  pages   = {4993--5005},
  doi     = {10.1021/acs.jctc.0c00332},
}
```
