## xTB run report (B-layer reading)

xTB run report
Source: nh3_planar_hess.out

Program: xTB 6.7.1 (build edcfbbe, compiled on 2024-07-22)
Task: frequencies

Total energy: -4.407989776233 Eh
Gradient norm: 0.105528858018 Eh/alpha
HOMO-LUMO gap: 12.788786 eV
SCF: converged in 8 iteration(s) (1 convergence marker(s): 8)

Frequencies: 12 mode(s) in the harmonic set; imaginary: 1
  Engine counters: 5 vibration(s), 1 imaginary
  -0.00  -0.00  0.00  0.00  0.00  0.00  -1337.66  1227.22  1227.47  4400.04  4525.52  4526.04
  (the block is printed 2 times; the identical copies are kept as measured)
  Most negative mode: -1337.66 cm-1; engine: found 1 significant imaginary frequency
  An imaginary mode marks a saddle or a non-stationary geometry at this level; a transition-state candidate must be re-verified at the target level.

Thermochemistry (from the engine's printed table):
  total free energy: -4.388893464062 Eh
  total energy: -4.407989776233 Eh
  zero point energy: 0.036237202996 Eh
Run finished: 2026/09/29 23:35:20.495

Reading notes: xTB (GFN2-xTB here) is the Tier-1 pre-screening level of this group's protocol chain -- its structures and energies rank and seed, they are not carried to a higher level as results (a measured example: a GFN2 transition state re-optimised at r2SCAN-3c showed seven imaginary modes, so a pre-screening structure is a starting point only). GFN2 covers all fifteen lanthanides; the actinides have no semi-empirical coverage in this stack (GFN-FF fails silently for them), so pre-screening actinide systems is out of scope for this engine. The 'normal termination of xtb' line goes to stderr -- keep it by redirecting 2>&1 when capturing the run.

## References

Complete citations:
- [bannwarth2021xtb] Bannwarth, C.; Caldeweyher, E.; Ehlert, S.; Hansen, A.; Pracht, P.; Seibert, J.; Spicher, S.; Grimme, S. (2021). Extended Tight-Binding Quantum Chemistry Methods. WIREs Computational Molecular Science, 11(2), e01493. DOI: 10.1002/wcms.1493
- [bannwarth2019gfn2] Bannwarth, C.; Ehlert, S.; Grimme, S. (2019). {GFN2-xTB}---An Accurate and Broadly Parametrized Self-Consistent Tight-Binding Quantum Chemical Method with Multipole Electrostatics and Density-Dependent Dispersion Contributions. Journal of Chemical Theory and Computation, 15(3), 1652-1671. DOI: 10.1021/acs.jctc.8b01176

BibTeX (paste-ready):

```bibtex
@article{bannwarth2021xtb,
  author  = {Bannwarth, Christoph and Caldeweyher, Eike and Ehlert, Sebastian and Hansen, Andreas and Pracht, Philipp and Seibert, Jakob and Spicher, Steffen and Grimme, Stefan},
  title   = {Extended Tight-Binding Quantum Chemistry Methods},
  journal = {WIREs Computational Molecular Science},
  year    = {2021},
  volume  = {11},
  number  = {2},
  pages   = {e01493},
  doi     = {10.1002/wcms.1493},
}
```

```bibtex
@article{bannwarth2019gfn2,
  author  = {Bannwarth, Christoph and Ehlert, Sebastian and Grimme, Stefan},
  title   = {{GFN2-xTB}---An Accurate and Broadly Parametrized Self-Consistent Tight-Binding Quantum Chemical Method with Multipole Electrostatics and Density-Dependent Dispersion Contributions},
  journal = {Journal of Chemical Theory and Computation},
  year    = {2019},
  volume  = {15},
  number  = {3},
  pages   = {1652--1671},
  doi     = {10.1021/acs.jctc.8b01176},
}
```
