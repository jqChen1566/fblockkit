## pysisyphus run report (B-layer reading)

pysisyphus run: work
  program: pysisyphus 1.0.0 (Python 3.13.9); executed at Mon Sep 28 11:41:22 2026 on 'localhost.localdomain'
  job: OPTIMIZATION; system: H₂O (3 atoms); coordinate system: redund (3)
  calculator: ORCA (calculator_000); charge 0, multiplicity 1; optimizer: rfo
  convergence thresholds: max(|force|) <= 0.00045 (overachieved 9e-05), rms(force) <= 0.0003 (overachieved 6e-05), max(|step|) <= 0.0018, rms(step) <= 0.0012
  run record (RUN.yaml): orca5 'hf def2-svp', pal 2, mem 1500 MB, thresh gau, max_cycles 150 (record version 1.0.0)

  cycle   d(energy)  max|force|  rms|force|   max|step|   rms|step|  s/cycle
      0           nan*   0.018610    0.012830    0.035705    0.028671     5.250
      1     -0.000535*   0.003988    0.003372    0.014404    0.008848     4.663
      2     -0.000056*   0.001703    0.001289    0.010590    0.006336     5.864
      3     -0.000011*   0.000060*   0.000042*   0.000167*   0.000117*    5.105

outcome: converged (the program's own marker)
  final: energy -75.96133849 Eh; max|force| 0.000060 (internal) / 0.000060 (cartesian)
  artifacts: trajectory (4 frames); closing geometry; run record RUN.yaml; 2 calculator output(s), last calculator_000.003.orca.out; structured history optimization.h5 (not read: HDF5)

cross-checks:
  - trajectory length vs cycle count: ok (4 frame(s) vs 4 cycle row(s))
  - closing frame energy vs its calculator call: ok (frame -75.96133849 Eh vs calculator_000.003.orca.out -75.961338488 Eh)
  - closing state (final summary) vs the last calculator call: ok (summary -75.96133849 Eh vs calculator_000.003.orca.out -75.961338488 Eh)
  - closing geometry file vs the closing state: ok (geometry -75.96133849 vs summary -75.96133849 Eh)
  - closing geometry coordinates vs the closing frame: ok (atoms and coordinates within 1e-8 Angstrom)

boundaries (the Wave-4.6 survey of the ORCA 6.1.1 pairing):
  - gradient-driven work is exact through ORCA's .engrad sidecar (the checks above); a TS search is reliable with a model Hessian;
  - frequency-dependent steps inside pysisyphus (hessian_init: calc, do_hess) fail on the '$multiplicity' block -- use a model Hessian and check frequencies outside pysisyphus;
  - the structured history (optimization.h5) is HDF5 and is not read here (h5py is not a dependency); the text artifacts carry the same numbers;
  - qm_calcs/cur_out is a symlink into the scratch directory, which pysisyphus cleans after every call, so it dangles after the run; the per-call copies (calculator_*.orca.*) are the durable ones.

## References

Complete citations:
- [steinmetzer2021pysisyphus] Steinmetzer, J.; Kupfer, S.; Gr{\"a}fe, S. (2021). pysisyphus: Exploring potential energy surfaces in ground and excited states. International Journal of Quantum Chemistry, 121(3), e26390. DOI: 10.1002/qua.26390

BibTeX (paste-ready):

```bibtex
@article{steinmetzer2021pysisyphus,
  author  = {Steinmetzer, Johannes and Kupfer, Stephan and Gr{\"a}fe, Stefanie},
  title   = {pysisyphus: Exploring potential energy surfaces in ground and excited states},
  journal = {International Journal of Quantum Chemistry},
  year    = {2021},
  volume  = {121},
  number  = {3},
  pages   = {e26390},
  doi     = {10.1002/qua.26390},
}
```
