## AVAS target projection (active-space construction)

Target and blocks:
  - target: 3 l=1 AO(s) of N at centre 0, shell(s) [2]
  - occupied block: 7 orbital(s) (occupation > 0.5); virtual block: 21
  - the export carries 3 fractional occupation(s) (0.0660, 0.0660, 0.0019): this is a correlated export (CASSCF), so the occupied/virtual split is a convention, not the SCF aufbau. For a pure AVAS construction use an SCF export

Occupied-side target overlaps (descending):
  0.6963 0.5595 0.5595 0.0000 0.0000 0.0000 0.0000
Virtual-side target overlaps (descending):
  0.4405 0.4405 0.3037 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000

Recommendation (0.1 truncation, option 3): (6.13388 electrons, 6 orbitals) -- occupied [0, 1, 2], virtual [7, 8, 9] (0-based orbital indices of the export)

Reading:
  - occupied side: 3 orbital(s) above the 0.1 threshold; the largest overlap(s) 0.6963 0.5595 0.5595 0.0000 ... -- values near 1 are essentially the target AO itself, and low values among the kept ones name orbitals strongly mixed with their environment
  - virtual side: 3 orbital(s) above the threshold (largest 0.4405) -- these are the target-centred antibonding combinations that belong in the active space
  - recommended active space: (6.13388 electrons, 6 orbitals) = 3 occupied + 3 virtual above the threshold (the electron count uses each active orbital's own occupation from the export, so a correlated export contributes its fractional occupations)
  - both AVAS criteria are of the falsify-only kind (the source states it): a reading here that shows nothing wrong does not establish that the space is the best one -- the space-change SVD and the 'CASCI below the variational HF energy' check can still falsify it

For the same target, ORCA's built-in AVAS block (when its own minimal basis covers the element):
  %scf
    avas
      system
        shell  2, 2, 2
        l      1, 1, 1
        m_l    z, x, y
        center 0, 0, 0
      end
    end
  end

Boundary of this route: the target set here is a subset of the calculation's own AOs (the source's non-minimal-ANO variant), which is what makes the projection computable from an export alone and is the only route open to the f block -- ORCA's AVAS minimal basis was measured to have no f-block entries (Eu). And ORCA's CASSCF takes its active space by orbital order, not by an index list: the (n_el, n_orb) above sizes the input and the window, it does not name the orbitals to a restart.

## References

Complete citations:
- [sayfutyarova2017automated] Sayfutyarova, E. R.; Sun, Q.; Chan, G. K.; Knizia, G. (2017). Automated Construction of Molecular Active Spaces from Atomic Valence Orbitals. Journal of Chemical Theory and Computation, 13(9), 4063-4078. DOI: 10.1021/acs.jctc.7b00128

BibTeX (paste-ready):

```bibtex
@article{sayfutyarova2017automated,
  author  = {Sayfutyarova, Elvira R. and Sun, Qiming and Chan, Garnet Kin-Lic and Knizia, Gerald},
  title   = {Automated Construction of Molecular Active Spaces from Atomic Valence Orbitals},
  journal = {Journal of Chemical Theory and Computation},
  year    = {2017},
  volume  = {13},
  number  = {9},
  pages   = {4063--4078},
  doi     = {10.1021/acs.jctc.7b00128},
}
```
