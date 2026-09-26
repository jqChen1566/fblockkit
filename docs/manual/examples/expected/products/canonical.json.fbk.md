## Orbital-space comparison (sigma_F / space-change SVD)

Space step_a [4 9]: 6 orbitals; space step_a.loc [4 9]: 6 orbitals.

Singular values of C_step_a [4 9]^T S C_step_a.loc [4 9] (descending, 6 values):
  1.0000 1.0000 1.0000 1.0000 1.0000 1.0000
  sigma_F = ||M||_F / sqrt(min(n_A, n_B)) = 1.000000
  deficit 1 - sigma_F = 7.327e-15; smallest singular value = 1.000000

Reading (provisional bands; the sigma_F source -- Guan & Jiang, Eq. (10), arXiv:2607.08178 as recorded in the project notes, not re-verified -- and the S_change source report 1 - sigma_F = 2.9e-3/6.8e-5 and a 0.65-0.99 SVD range respectively, and the stated use is ranking candidate spaces, not an absolute pass mark):
  - containment: the smaller space is contained in the larger to numerical precision (deficit vs the 1e-4 / 1e-2 provisional lines)
  - space change (equal sizes; the source's measured range for its own runs is 0.65-0.99): the two spaces are essentially the same space (smallest singular value vs the 0.9 / 0.5 provisional lines)

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
