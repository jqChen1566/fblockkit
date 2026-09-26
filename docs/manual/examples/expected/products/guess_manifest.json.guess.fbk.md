## G5 WASP initial guess (interpolated orbitals for a geometry series)

WASP initial guess for a geometry series (weighted orbital interpolation, 1/d over the neighbourhood):
  neighbour            d (Angstrom)      weight
  n2_scan_1.094.loc        0.357796     0.664011
  n2_scan_2.600.loc        0.707107     0.335989

  orthonormalisation residual: 3.10e-15
  occupations/energies from:   n2_scan_1.094.loc

Criteria and boundaries:
  - the mixture C = sum_beta w_beta C_beta / sum w_beta with w_beta = 1/d, over the neighbourhood delta = all given neighbours (the source's Eqs. (11)-(13)); d is the RMSD without alignment (the coefficients live in each structure's own AO frame)
  - the mixture is orthonormalised in the target geometry's overlap metric (Loewdin); the residual max |C^T S C - I| is 3.10e-15
  - the occupations and orbital energies written with the guess come from the nearest neighbour ('n2_scan_1.094.loc'); the orbitals themselves are the mixture
  - The written file is a Molekel mkl: run ``orca_2mkl <name>.fbk -gbw`` to turn it into a ``.gbw`` and read it with ``!moread`` + ``%moinp`` (measured: ORCA accepts the converted file as INITIAL GUESS: MOREAD).
  - This module builds the guess; ORCA re-optimises the orbitals. The mixture is not a converged orbital set and is not claimed to be one.

## References

Complete citations:
- [wardzala2026multireference] Wardzala, J. J.; Hennefarth, M. R.; Agarawal, V.; Jangid, B.; Seal, A.; Hermes, M. R.; King, D. S.; Gagliardi, L. (2026). Multireference Methods for Chemistry and Materials Science: Automated Active Spaces, Efficient Dynamic Correlation, and Extended Systems. Chemical Reviews, 126(8), 4592-4618. DOI: 10.1021/acs.chemrev.5c00866

BibTeX (paste-ready):

```bibtex
@article{wardzala2026multireference,
  author  = {Wardzala, Jacob J. and Hennefarth, Matthew R. and Agarawal, Valay and Jangid, Bhavnesh and Seal, Aniruddha and Hermes, Matthew R. and King, Daniel S. and Gagliardi, Laura},
  title   = {Multireference Methods for Chemistry and Materials Science: Automated Active Spaces, Efficient Dynamic Correlation, and Extended Systems},
  journal = {Chemical Reviews},
  year    = {2026},
  volume  = {126},
  number  = {8},
  pages   = {4592--4618},
  doi     = {10.1021/acs.chemrev.5c00866},
}
```
