## CREST conformer ensemble (B-layer reading)

CREST conformer ensemble report
Source: crest_run

Conformers: 2 (ensemble order; the energies are the frame comments)
Best conformer (crest_best.xyz): index 1

 idx   E (Eh)               rel (kcal/mol)   weight (298.15 K)
    1  -13.665127580000         0.000         0.7322
    2  -13.664177370000         0.596         0.2678

Cross-check: the frame energies reproduce crest.energies within 0.0003 kcal/mol (print precision).

Reading notes: the ensemble energies are GFN2-xTB level (Tier 1) -- the ordering and the weights select candidates, they do not rank final stabilities; for a quantitative comparison re-optimise the leading conformers at the target level (ORCA). The ensemble was deduplicated by CREST's own RMSD/energy criteria (CREGEN), so indices are post-filter positions, not sampling order. The Boltzmann weights use R = 1.9872042586e-3 kcal mol-1 K-1 over the listed set only; a run with several distinct basins at a close energy has its weights split across the mirror images if the engine did not merge them.

## References

Complete citations:
- [pracht2020crest] Pracht, P.; Bohle, F.; Grimme, S. (2020). Automated Exploration of the Low-Energy Chemical Space with Fast Quantum Chemical Methods. Physical Chemistry Chemical Physics, 22(14), 7169-7192. DOI: 10.1039/C9CP06869D

BibTeX (paste-ready):

```bibtex
@article{pracht2020crest,
  author  = {Pracht, Philipp and Bohle, Fabian and Grimme, Stefan},
  title   = {Automated Exploration of the Low-Energy Chemical Space with Fast Quantum Chemical Methods},
  journal = {Physical Chemistry Chemical Physics},
  year    = {2020},
  volume  = {22},
  number  = {14},
  pages   = {7169--7192},
  doi     = {10.1039/C9CP06869D},
}
```
