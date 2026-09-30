## 3.1 ranked-orbital active-space selection (APC-2)

Ranked-orbital active-space selection (APC-2), orbital energies as the model gap:
  system: n2_apc: doubly occupied 7, virtual window 7-27 (21 of 21 virtuals, lowest in energy).
  scheme: the balanced variant set 2 high-entropy virtual(s) aside during the entropy evaluation (MO 8, 7); they stay candidates at the top of the entropy scale (King et al. 2022).

  rank   MO  role                    S
     1    7  virtual(aside)     0.1914
     2    8  virtual(aside)     0.1914
     3    5  doubly occupied    0.1814
     4    6  doubly occupied    0.1814
     5    4  doubly occupied    0.1790
     6    3  doubly occupied    0.1618
     7   17  virtual            0.1170
     8    9  virtual            0.1169
     9    2  doubly occupied    0.1022
    10   10  virtual            0.0918
    11   11  virtual            0.0752
    12   13  virtual            0.0725
    13   14  virtual            0.0725
    14   12  virtual            0.0721
    15   15  virtual            0.0632
    16   16  virtual            0.0632
    17   18  virtual            0.0562
    18   19  virtual            0.0562
    19   20  virtual            0.0419
    20   21  virtual            0.0419
    21   24  virtual            0.0245
    22   27  virtual            0.0221
    23   22  virtual            0.0191
    24   23  virtual            0.0191
    25   25  virtual            0.0139
    26   26  virtual            0.0139
    27    1  doubly occupied    0.0054
    28    0  doubly occupied    0.0054

Selection:
  cap max(10,10) = 19404 CSFs; selected (10e, 10o), N_CSF = 19404 <= 19404
  active orbitals: 2, 3, 4, 5, 6, 7, 8, 9, 10, 17
  dropped (in drop order): 0, 1, 25, 26, 22, 23, 27, 24, 20, 21, 18, 19, 15, 16, 12, 13, 14, 11

Boundaries and checks:
  - closed-shell RHF exports only; the source's singly-occupied-orbital rule (assign the maximum approximated entropy) and its UNO variant are not covered
  - APC systematically overestimates the doubly-occupied orbital entropies (R^2 = 0.64, MAE 0.0240 against DMRG; virtuals 0.83/0.0064; ranking precision about 88% on the source's set)
  - the source expects the scheme to perform worse in much larger systems and where the HF determinant is a poor approximation; the ranking is a screening device, not a converged answer
  - the entropies and the ranking depend on the candidate window; the window is printed above, and N_CSF is checked with eq. (2) at every drop

## References

Complete citations:
- [king2021ranked] King, D. S.; Gagliardi, L. (2021). A Ranked-Orbital Approach to Select Active Spaces for High-Throughput Multireference Computation. Journal of Chemical Theory and Computation, 17(5), 2817-2831. DOI: 10.1021/acs.jctc.1c00037
- [king2022benchmark] King, D. S.; Hermes, M. R.; Truhlar, D. G.; Gagliardi, L. (2022). Large-Scale Benchmarking of Multireference Vertical-Excitation Calculations via Automated Active-Space Selection. Journal of Chemical Theory and Computation, 18(10), 6065-6076. DOI: 10.1021/acs.jctc.2c00630

BibTeX (paste-ready):

```bibtex
@article{king2021ranked,
  author  = {King, Daniel S. and Gagliardi, Laura},
  title   = {A Ranked-Orbital Approach to Select Active Spaces for High-Throughput Multireference Computation},
  journal = {Journal of Chemical Theory and Computation},
  year    = {2021},
  volume  = {17},
  number  = {5},
  pages   = {2817--2831},
  doi     = {10.1021/acs.jctc.1c00037},
}
```

```bibtex
@article{king2022benchmark,
  author  = {King, Daniel S. and Hermes, Matthew R. and Truhlar, Donald G. and Gagliardi, Laura},
  title   = {Large-Scale Benchmarking of Multireference Vertical-Excitation Calculations via Automated Active-Space Selection},
  journal = {Journal of Chemical Theory and Computation},
  year    = {2022},
  volume  = {18},
  number  = {10},
  pages   = {6065--6076},
  doi     = {10.1021/acs.jctc.2c00630},
}
```
