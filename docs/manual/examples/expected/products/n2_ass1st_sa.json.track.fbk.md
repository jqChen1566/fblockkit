## Cross-run state tracking

Cross-run state tracking (source: n2_ass1st_sa.json)

Runs (in tracking order):
  1. n2_ass1st_sa.json
  2. n2_ass1st_sa4.json
Target: root 0 (mult 1) of run 1; -108.940993 Eh

Step 1 -> 2 (n2_ass1st_sa4.json): candidates against the tracked density, rotated into this run's basis
  mult  root  energy (Eh)     W0          D           Q
  1     0     -108.939702     1.6665e-06  7.2159e-03  7.2176e-03  [chosen]
  1     1     -108.557060     1.4740e-01  8.9991e-01  1.0473e+00
  1     2     -108.527455     1.7101e-01  8.8692e-01  1.0579e+00
  1     3     -108.527455     1.7101e-01  8.8692e-01  1.0579e+00
  margin over the runner-up: 1.0401e+00

Boundaries (declared substitutions against the source's criterion
Q = W0 + W1 + D, J. Chem. Theory Comput. 2019, 15, 4790):
  - the W1 term (the active-to-virtual stationarity measure) is not
    evaluated: no ORCA export carries the coupling it needs. For a
    converged candidate run W1 tends to zero (the source's equivalence:
    W1 = 0 iff the energy is stationary), which is the intended input;
  - the source's 1/n_CAS scaling of D is not applied; the ranking within
    one sequence is unaffected and the values here are unscaled
    Frobenius norms of the full-space 1-RDM difference;
  - omega (the source's energy target) is taken as the tracked state's
    own energy, updated each step;
  - the geometry gate: the runs must share atoms, coordinates (1e-6
    Angstrom) and the AO dimension -- cross-geometry tracking is outside
    the current scope.

## References

Complete citations:
- [tran2019tracking] Tran, L. N.; Shea, J. A. R.; Neuscamman, E. (2019). Tracking Excited States in Wave Function Optimization Using Density Matrices and Variational Principles. Journal of Chemical Theory and Computation, 15(9), 4790-4803. DOI: 10.1021/acs.jctc.9b00351

BibTeX (paste-ready):

```bibtex
@article{tran2019tracking,
  author  = {Tran, Lan Nguyen and Shea, Jacqueline A. R. and Neuscamman, Eric},
  title   = {Tracking Excited States in Wave Function Optimization Using Density Matrices and Variational Principles},
  journal = {Journal of Chemical Theory and Computation},
  year    = {2019},
  volume  = {15},
  number  = {9},
  pages   = {4790--4803},
  doi     = {10.1021/acs.jctc.9b00351},
}
```
