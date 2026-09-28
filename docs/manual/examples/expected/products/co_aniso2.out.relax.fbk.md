## Magnetic relaxation and QTM

Single_aniso segment 1 of 1 (co_aniso2.out)
  spin-orbit spectrum (cm-1): 0.000 0.000 29361.197 29361.197

  KD / group   S     E (cm-1)      g1        g2        g3      theta(deg)  g_T   g_T*theta
     1         0.5       0.000    1.99952   2.00000   2.00000       0.00  1.333       n/a  *
     2         0.5   29361.197    2.00000   2.00000   2.00048      90.00  2.000       n/a  *
  * the g anisotropy of this group is below 1% of its largest value: its main-axis direction is set by the print's digits, so the axis-angle criteria do not apply to it (CO+ fixture: Delta g/g ~ 2e-4)

  U_eff estimate: no computed group is flagged (1 of them are skipped: main axes not defined); U_eff is at least 29361.2 cm-1 as far as the computed groups go
  The flags are the guide's empirical thresholds (15 deg non-collinearity; g_T*theta > 20 from 20 Dy(III) SMMs) -- they bound a plausible barrier, not a measured one.

  UBAR (ab initio blocking barrier, qualitative per the engine):
    group 1: m+ = -0.999761  m- = +0.999761  E = 0.0000 cm-1
    group 2: m+ = -1.000239  m- = +1.000239  E = 29361.1969 cm-1
    intra-group |<i+|mu|i->| averages (mu_B): g1=6.667e-01  g2=6.667e-01
    largest inter-group element: |<1.1+|mu_Y|2.1+>| = 7.2538e-02 mu_B
    the 'check the tunnelling splitting instead' note is the engine's fixed template sentence (printed for odd-electron systems too)

  the segment also carries the static susceptibility / magnetisation tables

## References

Complete citations:
- [chilton2025abinitio] Chilton, N. F. (2025). Ab initio electronic structure calculations of lanthanide single-molecule magnets; a practical guide. Chemical Society Reviews, 54(24), 11468-11487. DOI: 10.1039/d5cs00493d

BibTeX (paste-ready):

```bibtex
@article{chilton2025abinitio,
  author  = {Chilton, Nicholas F.},
  title   = {Ab initio electronic structure calculations of lanthanide single-molecule magnets; a practical guide},
  journal = {Chemical Society Reviews},
  year    = {2025},
  volume  = {54},
  number  = {24},
  pages   = {11468--11487},
  doi     = {10.1039/d5cs00493d},
}
```
