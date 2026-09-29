## Polynuclear magnetism (POLY_ANISO)

Polynuclear magnetism (POLY_ANISO) report (two_center_probe.out)
  centers: 2 independent, 2 in total
    center 1: data file aniso_1.input; at (0.000, 0.000, 0.000) A; 4 spin-orbit states; g = 2.0000 / 2.0000 / 1.9995
      spin-orbit spectrum (cm-1): 0.00 0.00 29361.20 29361.20
    center 1: data file aniso_2.input; at (0.000, 0.000, 3.700) A; 4 spin-orbit states; g = 2.0000 / 2.0000 / 1.9995
      spin-orbit spectrum (cm-1): 0.00 0.00 29361.20 29361.20

  exchange: 4 coupled states; 1 pairs; Lines-1 INCLUDED; dipole-dipole INCLUDED; ITO decomposition INCLUDED
    pair 1 (centers 1-2): J = 0.10000 cm-1

  first-order anisotropic coupling (weights of the decomposition):
    pair 1-2 [LINES-1]: Isotropic 32.708% / Symmetric 93.839% / Anti-Symmetric 11.155%
       7.447400e-02  -5.560449e-03  -6.650324e-02
       4.224047e-03  -9.905969e-02   1.301287e-02
      -6.660148e-02  -1.250034e-02  -7.353884e-02
    pair 1-2 [DIPOLE-DIPOLE]: Isotropic 93.448% / Symmetric 33.524% / Anti-Symmetric 11.982%
       6.367713e-02   1.900942e-03  -6.197318e-03
      -1.443735e-03   3.385331e-02  -4.447663e-03
      -1.688008e-02   4.273470e-03   3.801499e-02

  coupled states (cm-1):
    state 1: lines -0.025000  dipole -0.017100  total -0.042100  relative 0.000000
    state 2: lines -0.025000  dipole -0.017091  total -0.042091  relative 0.000008
    state 3: lines -0.025000  dipole 0.000004  total 0.009187  relative 0.051287
    state 4: lines 0.075000  dipole 0.034187  total 0.075004  relative 0.117104

  chiT(T): 101 points, T = 0.0001 to 300 K; chiT = 0.501549 -> 0.750429 cm3 K mol-1
  Van Vleck susceptibility tensors at 101 temperatures; main values at 0.0001 K: 0.001914 / 0.002153 / 1.500582; at 300 K: 0.750350 / 0.750442 / 0.750496
  note: Computation of the magnetization torque ... skipped by the user
  note: Computation of the molar magnetization ... skipped by the user
  note: finished ok

  Reading notes: the exchange constants are the caller's input (measured elsewhere or fitted) -- this workflow never computes or fits them; a joint ab initio computation of exchange splittings (the LDF-CAHF / many-state PNO-CASPT2 route) is outside the ORCA ecosystem and is registered as a documented termination, not implemented here. The Lines model is exact only for two isotropic spins, one Ising plus one isotropic spin, or two Ising spins, and approximate otherwise; the dipole-dipole coupling is evaluated exactly from the ab initio moments and usually dominates in strongly anisotropic lanthanides (ORCA manual section 7.18; method reference: Ungur & Chibotaru, Chem. Eur. J. 2017).

## References

Complete citations:
- [ungur2017abinitio] Ungur, L.; Chibotaru, L. F. (2017). Ab Initio Crystal Field for Lanthanides. Chemistry -- A European Journal, 23(15), 3708-3718. DOI: 10.1002/chem.201605102

BibTeX (paste-ready):

```bibtex
@article{ungur2017abinitio,
  author  = {Ungur, Liviu and Chibotaru, Liviu F.},
  title   = {Ab Initio Crystal Field for Lanthanides},
  journal = {Chemistry -- A European Journal},
  year    = {2017},
  volume  = {23},
  number  = {15},
  pages   = {3708--3718},
  doi     = {10.1002/chem.201605102},
}
```
