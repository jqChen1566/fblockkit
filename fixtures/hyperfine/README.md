# fixtures/hyperfine -- the ORCA EPRNMR probes (menu 40)

Real ORCA 6.1.1 `%eprnmr` runs on CeF3 (Ce 3+ f1; the geometry of the
ailft fixtures), 2026-09-29, on 101.  Four probes cover the interface and
its boundaries:

| file | route | shows |
|---|---|---|
| `cef3_epr_dft.out` (+`.inp`) | PBE0/x2c-SVPall (all-electron), `Nuclei = all Ce { aiso, adip, fgrad, rho }` | the full section without nuclear parameters: EFG (V(Tot), V(El)/V(Nuc), orientation), RHO(0) = 727216.913536801 a.u.**-3, A matrix all zeros |
| `cef3_epr_nuc.out` (+`.inp`) | same, with `PPP=1.5, III=2.5` (values arbitrary -- a mechanism probe) | the A tensor printed: A(FC) = -11.8696 MHz isotropic, A(SD) anisotropic; the engine "Quadrupole tensor eigenvalues" block with eta = 0.040223 |
| `cef3_epr_nuc_q.out` (+`.inp`) | same, with `PPP=1.5, QQQ=0.5, III=2.5` | engine `e**2qQ = -165.764347 MHz` for Q = 0.5 barn -- the cross-check anchor for the CODATA-derived conversion constant (234.9648 MHz per barn per a.u.; the derived value on the printed Vzz agrees to ~1e-6 relative) |
| `cef3_epr_casscf.out` (+`.inp`) | CASSCF (def2-SVP = small-core ECP28) + `rel DoSOC` | the CASSCF route: method "CASSCF/ALL STATES AVERAGE", EFG and RHO(0) computed (0.075618828 -- the valence density of the ECP basis), A components switched OFF |

## Measured interface facts

- **the coordinates block must precede `%eprnmr` when `Nuclei` is given**
  (otherwise the engine stops with "[EPRNMR] block: nuclear properties are
  requested but no coordinates have been read!"); `Nuclei` counts atoms
  from 1 (the printed nucleus labels use the engine's 0-based index);
- the A tensor in MHz needs the nuclear parameters (I, P; `PPP`/`III`);
  without them the engine prints zeros -- not an error;
- RHO(0) is basis-domain dependent: the all-electron probe reports the
  core density (order 10^5), the ECP probe the valence density (order
  10^-1) -- never compare across domains;
- **the EFG is sensitive to SCF micro-solutions**: probes 1 and 3 (same
  input energy, -8527.705850861650 Eh, reproduced digit-for-digit) share
  V(Tot) = 0.6771983 / 0.7337722 / -1.4109705, while probe 2 converged to
  a slightly lower micro-solution (-8527.705853214089 Eh) whose EFG
  differs by ~1e-4 relative -- converged energy does not imply a
  reproducible EFG (recorded, not a defect);
- the nuclear EFG triplet V(Nuc) is identical across the DFT and CASSCF
  probes of the same geometry (sorted), and V(El) + V(Nuc) reproduces
  V(Tot) to the printed digits -- both are used as report cross-checks.
