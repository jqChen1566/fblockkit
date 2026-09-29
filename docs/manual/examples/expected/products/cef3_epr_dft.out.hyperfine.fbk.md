## Hyperfine / EFG report

Hyperfine / EFG report (cef3_epr_dft.out)
  method: SCF; multiplicity 2; nuclei computed: A(iso) 1  A(dip) 1  A(orb) 0  A(dia) 0  EFG 1  Rho(0) 1

  nucleus Ce (engine index 0):
    nuclear parameters: I = 0.0, P = 0.0 MHz/au^3, Q = 0.0 barn
    A tensor (MHz): all zeros -- nuclear parameters (I, P) were not supplied in the input (measured engine behaviour)
    EFG: V(Tot) = (0.677198, 0.733772, -1.410970) a.u.; |Vzz| = 1.410970, eta = 0.0401
      V(El) = (0.489266, 0.545878, -1.035145)  V(Nuc) = (0.187932, 0.187894, -0.375826)  (sum vs V(Tot): max deviation 1.00e-07)
      with Q = 0.5 barn: eQVzz/h = -165.7642 MHz; first-order DeltaE_Q = 82.9043 MHz (MHz == NQR nu_Q)
    Rho(0) = 727217 a.u.**-3  [basis-domain dependent: all-electron carries the core, small-core ECP reports the valence density only]

  Reading notes: eta = (Vxx - Vyy)/Vzz with |Vxx| <= |Vyy| <= |Vzz|; the quadrupole splitting is first order (I >= 3/2, higher-order terms neglected). The A tensor values are the engine's 'SAI' spin-Hamiltonian convention; their MHz scale needs the per-nucleus parameters (I, P) in the input -- zeros mean they were not supplied (measured). Rho(0) must not be compared across basis domains (all-electron vs ECP).

## References

Complete citations:
- [aerts2019nqcc] Aerts, A.; Brown, A. (2019). A Revised Nuclear Quadrupole Moment for Aluminum: Theoretical Nuclear Quadrupole Coupling Constants of Aluminum Compounds. The Journal of Chemical Physics, 150, 224302. DOI: 10.1063/1.5097151

BibTeX (paste-ready):

```bibtex
@article{aerts2019nqcc,
  author  = {Aerts, Antoine and Brown, Alex},
  title   = {A Revised Nuclear Quadrupole Moment for Aluminum: Theoretical Nuclear Quadrupole Coupling Constants of Aluminum Compounds},
  journal = {The Journal of Chemical Physics},
  year    = {2019},
  volume  = {150},
  pages   = {224302},
  doi     = {10.1063/1.5097151},
}
```
