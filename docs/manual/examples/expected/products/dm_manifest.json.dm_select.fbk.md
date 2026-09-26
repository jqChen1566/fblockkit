## A12 dipole-moment active-space selection

Dipole-moment active-space selection (DM-AS):
  reference: 2.0801 D  (h2o_dm_ref_pbe0.out (SCF block))
  protocol:  GDM-AS

  rank   active space     mu (D)   deviation (D)
     1   ( 6e,  8o)     2.0720        0.0081  <- selected
     2   ( 8e,  8o)     2.0960        0.0159
     3   (10e,  8o)     2.0974        0.0173
     4   ( 8e,  7o)     2.0974        0.0174
     5   ( 6e,  7o)     2.1017        0.0217
     6   ( 6e,  6o)     2.1046        0.0245

Selected active space: (6e, 8o)  [deviation 0.0081 D]

Criteria and boundaries:
  - the protocol is GDM-AS: the absolute difference of the total dipole moments against the reference (2.0801 D from h2o_dm_ref_pbe0.out (SCF block)); the smallest deviation wins and ties go to the smaller space
  - each candidate is read from its single-root dipole block (State: 0, relaxed density); state-averaged outputs are refused (they print only the average)
  - The source runs the candidates as state-averaged calculations over six states and takes the S0 dipole of that ensemble; the ORCA route here runs each candidate with nroots 1 (the only form whose S0 dipole ORCA prints). The two differ by the state-averaging relaxation of the density -- the ranking between candidates is what the protocol uses.
  - The reference value drives the choice: the source tested experimental (NIST) dipole moments and KS-DFT values with several functionals, and the 2026 recommendations name R2SCAN-family references as unsuitable for the EDM variant. State the reference and its origin in the report when publishing.
  - The protocol needs a nonzero dipole moment, and the source excludes charged systems (the dipole of a charged molecule depends on the coordinate origin).
  - For a potential-energy scan the protocol offers no guarantee that every geometry selects the same active space; select at one geometry and carry the space with the cross-structure mapping (menu 17) when consistency matters.
  - The per-state variants (EDM-AS, D2DM-AS) are not implemented: their data is not printed by ORCA (see the module docstring). The source's own recommendation names CASCI-GDM-AS the default protocol.

## References

Complete citations:
- [kaufold2023dipole] Kaufold, B. W.; Chintala, N.; Pandeya, P.; Dong, S. S. (2023). Automated Active Space Selection with Dipole Moments. Journal of Chemical Theory and Computation, 19(8), 2469-2483. DOI: 10.1021/acs.jctc.2c01128
- [kaufold2026casci] Kaufold, B. W.; Dong, S. S. (2026). Automated Active Space Selection with CASCI Dipole Moments. Journal of Chemical Theory and Computation, 22(15), 7624-7641. DOI: 10.1021/acs.jctc.6c00473

BibTeX (paste-ready):

```bibtex
@article{kaufold2023dipole,
  author  = {Kaufold, Benjamin W. and Chintala, Nithin and Pandeya, Pratima and Dong, Sijia S.},
  title   = {Automated Active Space Selection with Dipole Moments},
  journal = {Journal of Chemical Theory and Computation},
  year    = {2023},
  volume  = {19},
  number  = {8},
  pages   = {2469--2483},
  doi     = {10.1021/acs.jctc.2c01128},
}
```

```bibtex
@article{kaufold2026casci,
  author  = {Kaufold, Benjamin W. and Dong, Sijia S.},
  title   = {Automated Active Space Selection with CASCI Dipole Moments},
  journal = {Journal of Chemical Theory and Computation},
  year    = {2026},
  volume  = {22},
  number  = {15},
  pages   = {7624--7641},
  doi     = {10.1021/acs.jctc.6c00473},
}
```
