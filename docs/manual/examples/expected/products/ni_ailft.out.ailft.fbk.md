## Ab initio ligand-field analysis

Ab initio ligand-field analysis (ni_ailft.out)
  d8 configuration, 2 CI blocks, MOs 9 to 13, metal center atom 0

  == CASSCF level ==
    ligand-field eigenfunctions (cm-1): 0.0 0.0 0.0 0.0 0.0   (spread 0.0 cm-1)
    Slater-Condon (cm-1): F0dd(from 2el Ints) = 215295.8 (fixed)  F2dd = 99138.0  F4dd = 61305.7
    Racah (cm-1): B = 1328.1  C = 4865.5  C/B = 3.663
    LFT orbitals stored in ni_ailft.casscf.lft.gbw

  == NEVPT2 level ==
    ligand-field eigenfunctions (cm-1): 0.0 0.3 1.4 2.6 2.8   (spread 2.8 cm-1)
    Slater-Condon (cm-1): F2dd = 91278.1  F4dd = 56869.6
    Racah (cm-1): B = 1218.0  C = 4513.5  C/B = 3.706
    LFT orbitals stored in ni_ailft.nevpt2.lft.gbw

  fit quality (CASSCF): total RMS = 0.0 cm-1  [block 1 0.0  block 2 0.0]  Pearson = 1.000

  fit quality (NEVPT2): total RMS = 457.4 cm-1  [block 1 315.9  block 2 517.9]  Pearson = 1.000
  SOC constant: ZETA_D = 664.14 cm-1 (SOC based on casscf orbitals)

  Reading notes: the near-zero CASSCF-level RMS is intrinsic -- the LFT parametrization is exact for that level; the correlation-level RMS reflects (among other things) the neglected anisotropy of electron-electron repulsion in covalent complexes (Lang, Atanasov & Neese, J. Phys. Chem. A 2020). The parameters are model quantities of the ligand-field Hamiltonian; this menu reads the engine's fit, it never refits. Nephelauxetic comparisons (Jung, Atanasov & Neese, Inorg. Chem. 2017) need the caller's free-ion references -- pass them to get the ratios.

## References

Complete citations:
- [lang2020ailft] Lang, L.; Atanasov, M.; Neese, F. (2020). Improvement of Ab Initio Ligand Field Theory by Means of Multistate Perturbation Theory. The Journal of Physical Chemistry A, 124(6), 1025-1037. DOI: 10.1021/acs.jpca.9b11227

BibTeX (paste-ready):

```bibtex
@article{lang2020ailft,
  author  = {Lang, Lucas and Atanasov, Mihail and Neese, Frank},
  title   = {Improvement of Ab Initio Ligand Field Theory by Means of Multistate Perturbation Theory},
  journal = {The Journal of Physical Chemistry A},
  year    = {2020},
  volume  = {124},
  number  = {6},
  pages   = {1025--1037},
  doi     = {10.1021/acs.jpca.9b11227},
}
```
