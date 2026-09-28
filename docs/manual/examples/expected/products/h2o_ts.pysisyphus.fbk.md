## pysisyphus input (B-layer generation)

pysisyphus input for 3 atom(s)
  job: transition-state search (tsopt/rsprfo, model Hessian)
  calculator: orca5, keywords 'HF 3-21G', charge 0, multiplicity 1, pal 1, mem 1500 MB per core
  threshold: baker
  hessian: model (fischer)
  written: h2o_ts.pysisyphus.xyz (structure, XYZ) and h2o_ts.pysisyphus.yaml (the input)

- run it on a machine with pysisyphus installed and the engine configured: the program reads the engine command from ~/.pysisyphusrc (an INI file: [orca5] cmd=/path/to/orca), or from $PATH when the section is absent
- give pysisyphus a real scratch disk before starting (it stages every calculator call in a temporary directory under $TMPDIR and copies the results back into qm_calcs/ at the end of each call)
- keep the console:  pysis <file> > run.out  -- that capture is the run's record (the cycle table, the outcome marker, the closing summary), and it is what menu 33 reads back
- frequency verification stays outside pysisyphus in this pairing: the quantum-Hessian route (hessian_init: calc / do_hess) crashes on ORCA 6's .hess format (measured), so run a plain ORCA Freq job on the closing geometry -- menu 30 is the toolkit's own route back from an imaginary mode
- bring the whole run directory back to menu 33: it cross-checks the console, the trajectory, the closing geometry and the per-cycle ORCA outputs against each other

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
