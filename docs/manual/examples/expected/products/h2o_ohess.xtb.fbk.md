## xTB run report (B-layer reading)

xTB run report
Source: h2o_ohess.out

Program: xTB 6.7.1 (build edcfbbe, compiled on 2024-07-22)
Task: geometry optimisation + frequencies

Total energy: -5.070544013094 Eh
Gradient norm: 0.000932444846 Eh/alpha
HOMO-LUMO gap: 14.397676 eV
SCF: converged in 3 iteration(s) (2 convergence marker(s): 3, 3)
Geometry optimisation: converged after 5 iteration(s)

Frequencies: 9 mode(s) in the harmonic set; imaginary: 0
  Engine counters: 3 vibration(s), 0 imaginary
  -0.00  -0.00  -0.00  -0.00  0.00  0.00  1538.60  3642.93  3657.89
  (the block is printed 2 times; the identical copies are kept as measured)

Thermochemistry (from the engine's printed table):
  total free energy: -5.068023573707 Eh
  total energy: -5.070544013094 Eh
  zero point energy: 0.020137684918 Eh
Run finished: 2026/09/29 23:33:55.899

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
