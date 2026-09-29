## Hyperfine / EFG report

Hyperfine / EFG report (cef3_epr_casscf.out)
  method: CASSCF/ALL STATES AVERAGE; nuclei computed: A(iso) 1  A(dip) 1  A(orb) 0  A(dia) 0  EFG 1  Rho(0) 1

  nucleus Ce (engine index 0):
    nuclear parameters: I = 0.0, P = 0.0 MHz/au^3, Q = 0.0 barn
    A tensor: not computed by the engine for this reference (measured: CASSCF switches the A components off; DFT with nuclear parameters P/I prints them)
    EFG: V(Tot) = (-0.086135, -0.086148, 0.172283) a.u.; |Vzz| = 0.172283, eta = 0.0001
      V(El) = (-0.274029, -0.274080, 0.548109)  V(Nuc) = (0.187894, 0.187932, -0.375826)  (sum vs V(Tot): max deviation 2.78e-17)
      with Q = 0.5 barn: eQVzz/h = 20.2402 MHz; first-order DeltaE_Q = 10.1201 MHz (MHz == NQR nu_Q)
    Rho(0) = 0.0756188 a.u.**-3  [basis-domain dependent: all-electron carries the core, small-core ECP reports the valence density only]

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
