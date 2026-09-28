## Judd-Ofelt intensity parameters

Judd-Ofelt intensity parameters (standard JO fit, S-space, weighted (1/S))
  dataset: babu2000_eu3.yaml
  host: Sellmeier n0 = 1 (A=1.24278 B=0.023833 um^2); chi_ED = (n^2+2)^2/9 (virtual cavity)
  Omega_2 = 16.409e-20 cm^2; Omega_4 = 9.658e-20 cm^2; Omega_6 = 1.698e-20 cm^2
  fit: 9 transitions, sigma = 4.2717e-05 a.u., sigma/S_max = 13.25 %

  transition            energy (cm-1)      f_exp      S_exp (a.u.)    S_ED/S_exp
  7F6 <- 7F1                   4541.3       1.9830e-06     3.2241e-04         0.709
  7F6 <- 7F0                   4791.6       1.2320e-06     6.3271e-05         1.389
  5D1 <- 7F1                  18726.6       4.5000e-07     1.7278e-05         0.882
  5D2 <- 7F0                  21459.2       3.3300e-07     3.6828e-06         1.273
  5D3 <- 7F1                  24154.6       3.0200e-07     8.7984e-06         0.776
  5L6 <- 7F1                  25000.0       1.3830e-06     3.8769e-05         0.150
  5L6 <- 7F0                  25380.7       3.3380e-06     3.0664e-05         0.291
  5G2 <- 7F0                  26178.0       5.2300e-07     4.6385e-06         0.758
  5D4 <- 7F0                  27700.8       4.8900e-07     4.0631e-06         1.103

boundaries (measured and documented in the close-reading record):
  - the local field is the squared virtual-cavity form (n^2+2)^2/9; the source papers print it without the square -- the square is what their reference code and their published numbers use (regression-tested here);
  - the U^(lambda) values are host-independent table values for the ion and must come from the same source as the fit (transcribe them together with the dataset and check them against the published fit);
  - magnetic-dipole contributions are not computed: give f_md per transition where it is significant (the fit subtracts it);
  - the extended (perturbative X_k) JO models of the 2022/2024 papers need free-ion atomic-structure wave functions and are outside this toolkit's post-processing scope; this menu implements the standard JO fit.

## References

Complete citations:
- [hovhannesyan2024extension] Hovhannesyan, G.; Boudon, V.; Lepers, M. (2024). Extension of {J}udd-{O}felt theory: Application on {Eu}$^{3+}$, {Nd}$^{3+}$ and {Er}$^{3+}$. Journal of Luminescence, 266, 120234. DOI: 10.1016/j.jlumin.2023.120234

BibTeX (paste-ready):

```bibtex
@article{hovhannesyan2024extension,
  author  = {Hovhannesyan, Gohar and Boudon, Vincent and Lepers, Maxence},
  title   = {Extension of {J}udd-{O}felt theory: Application on {Eu}$^{3+}$, {Nd}$^{3+}$ and {Er}$^{3+}$},
  journal = {Journal of Luminescence},
  year    = {2024},
  volume  = {266},
  pages   = {120234},
  doi     = {10.1016/j.jlumin.2023.120234},
}
```
