## TNASS best-k sweep (energy-ranked)

TNASS best-k sweep (the selection path's prefixes ranked by the CASCI energy):
  method: greedy; k = 1..6 of the window's spatial orbitals

    k  subset                    n_elec         E(k) [Eh]   dE [mEh]
    4  0, 3, 4, 5                     2    -230.655257977    138.561
    5  0, 1, 3, 4, 5                  4    -230.716007717     77.811
    6  0, 1, 2, 3, 4, 5               6    -230.793818898      0.000
    1  (skipped: the prefix carries 0 electrons)
    2  (skipped: the prefix carries 0 electrons)
    3  (skipped: the prefix carries 0 electrons)

  best: k = 6, orbitals [0, 1, 2, 3, 4, 5] (6e, 6o), E = -230.793818898 Eh
  compact pick (this tool): k = 6, orbitals [0, 1, 2, 3, 4, 5] -- the smallest prefix within 1 mEh of the sweep minimum

Boundaries and checks:
  - the environment convention is this tool's own (the source's main text leaves it to its CASCI implementation): the complement carries its natural occupations as a fractional closed-shell environment (h'_pq = h_pq + 1/2 sum w_i [2(pq|ii) - (pi|qi)]; E_core' = E_core + sum w_i h_ii + 1/4 sum w_i w_j [2(ii|jj) - (ij|ji)]), and each prefix's CI runs with round(sum of its occupations) electrons
  - the source selects k manually off the energy curve (in its words, 'while we are selecting the best k manually with respect to energy'); the table and the minimum marker serve that choice, and the compact pick is this tool's convenience
  - at k = the full window the environment is empty and E(k) reproduces the window's own FCI energy -- the identity pinned in the test suite
  - the prefix path is greedy: each k's subset is the k-prefix of the single greedy run (prefix-consistent); the source's best-k is defined over its greedy family the same way
  - delivering the winning space to a CASSCF needs the orbital-order machinery (menu 22's boundary note)

## References

Complete citations:
- [mingare2026tnass] Mingare, A.; Heuz{\'e}, I.; Coveney, P. V. (2026). {TNASS}: Tensor Network Active Space Selection with the Entanglement Feature. arXiv preprint arXiv:2608.03645. DOI: 10.48550/arXiv.2608.03645

BibTeX (paste-ready):

```bibtex
@article{mingare2026tnass,
  author  = {Mingare, Angus and Heuz{\'e}, Isabelle and Coveney, Peter V.},
  title   = {{TNASS}: Tensor Network Active Space Selection with the Entanglement Feature},
  journal = {arXiv preprint arXiv:2608.03645},
  year    = {2026},
  doi     = {10.48550/arXiv.2608.03645},
}
```
