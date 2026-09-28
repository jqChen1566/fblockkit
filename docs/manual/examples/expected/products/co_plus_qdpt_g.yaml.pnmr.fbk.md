## pNMR pseudocontact shifts

pNMR pseudocontact shifts (point-dipole approximation)
  structure: co_plus.xyz (2 atoms); centre: atom 1 (C)
  susceptibility from g-tensor (spin 0.5); route: nonsymmetric; T = 300 K
  chi tensor (1e-32 m^3):
        2.6150     -0.0003      0.0000
       -0.0003      2.6152      0.0000
        0.0000      0.0000      2.6153
  Delta chi_ax = 0.0004e-32 m^3; Delta chi_rh = -0.0005e-32 m^3; |antisymmetric|_max = 0.000e+00e-32 m^3
  note: the g-matrix was read from 'co_plus_qdpt.out' (effective-Hamiltonian block)

  nucleus        r (A)     theta (deg)   phi (deg)    delta_PCS (ppm)
  O2             1.130          90.0        90.0            0.085

boundaries (measured and from the sources):
  - the PDA is the long-range limit: the source's benchmark shows <10 % deviation of the dipolar part beyond ~7 A from the centre, while the contact (through-bond) contribution dominates below ~4-5 A -- a PCS value predicted for a nucleus close to the metal is a lower bound on the total shift error, not a full prediction;
  - only the traceless symmetric part of chi enters the PCS (the antisymmetric part contributes nothing; tested exactly) -- the route's difference (nonsymmetric vs symmetric construction) changes the symmetric part itself and hence the values;
  - the ZFS route uses the effective spin dyadic <SS> at the given temperature; near-axial systems can cross the node condition <SS>|| g|| = <SS>perp gperp where the axiality (and the PCS) vanish;
  - magnetic-dipole-free, isotropic-shift contexts only: the contact shift is not computed here (it needs hyperfine coupling).

## References

Complete citations:
- [mares2018pcs] Mare{\v{s}}, J.; Vaara, J. (2018). Ab initio paramagnetic {NMR} shifts via point-dipole approximation in a large magnetic-anisotropy {Co}({II}) complex. Physical Chemistry Chemical Physics, 20(35), 22547-22555. DOI: 10.1039/c8cp04123g
- [parker2020ligandfield] Parker, D.; Suturina, E. A.; Kuprov, I.; Chilton, N. F. (2020). How the Ligand Field in Lanthanide Coordination Complexes Determines Magnetic Susceptibility Anisotropy, Paramagnetic {NMR} Shift, and Relaxation Behavior. Accounts of Chemical Research, 53(8), 1520-1534. DOI: 10.1021/acs.accounts.0c00275

BibTeX (paste-ready):

```bibtex
@article{mares2018pcs,
  author  = {Mare{\v{s}}, Ji{\v{r}}{\'i} and Vaara, Juha},
  title   = {Ab initio paramagnetic {NMR} shifts via point-dipole approximation in a large magnetic-anisotropy {Co}({II}) complex},
  journal = {Physical Chemistry Chemical Physics},
  year    = {2018},
  volume  = {20},
  number  = {35},
  pages   = {22547--22555},
  doi     = {10.1039/c8cp04123g},
}
```

```bibtex
@article{parker2020ligandfield,
  author  = {Parker, David and Suturina, Elizaveta A. and Kuprov, Ilya and Chilton, Nicholas F.},
  title   = {How the Ligand Field in Lanthanide Coordination Complexes Determines Magnetic Susceptibility Anisotropy, Paramagnetic {NMR} Shift, and Relaxation Behavior},
  journal = {Accounts of Chemical Research},
  year    = {2020},
  volume  = {53},
  number  = {8},
  pages   = {1520--1534},
  doi     = {10.1021/acs.accounts.0c00275},
}
```
