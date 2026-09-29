## Magnetic entropy / magnetocaloric report

Magnetic entropy / magnetocaloric report (poly_mh.out)
  route: magnetization table (6 temperatures x 71 fields, powder-averaged, molar M in Bohr magnetons)

  -DeltaS(T, H) [J mol-1 K-1] (positive = direct MCE):
    T_mid [K] |   1.0 T |   2.0 T |   3.0 T |   5.0 T |   7.0 T
        1.900 |   1.012 |   3.392 |   5.963 |   9.431 |  10.840
        2.250 |   0.738 |   2.587 |   4.805 |   8.381 |  10.242
        2.750 |   0.497 |   1.818 |   3.572 |   6.961 |   9.225
        3.500 |   0.313 |   1.181 |   2.428 |   5.254 |   7.642
        4.500 |   0.188 |   0.728 |   1.549 |   3.659 |   5.821
  maximum: -DeltaS = 10.8404 J mol-1 K-1 at T = 1.90 K, H = 7.000 T

  Reading notes: DeltaS(T, H) = N_A mu_B int (dM/dT)_H dH' (the Maxwell relation; trapezoid over the printed field grid, adjacent-temperature differences; the printed sign convention is DeltaS = S(T,H) - S(T,0), so the table shows -DeltaS). The data must be converged (a fine temperature grid); the field grid starts at 0.0001 T.

## References

Complete citations:
- [szalowski2020mce] Sza{\l}owski, K.; Kowalewska, P. (2020). Magnetocaloric Effect in {Cu5-NIPA} Molecular Magnet: A Theoretical Study. Materials, 13(2), 485. DOI: 10.3390/ma13020485

BibTeX (paste-ready):

```bibtex
@article{szalowski2020mce,
  author  = {Sza{\l}owski, Karol and Kowalewska, Pamela},
  title   = {Magnetocaloric Effect in {Cu5-NIPA} Molecular Magnet: A Theoretical Study},
  journal = {Materials},
  year    = {2020},
  volume  = {13},
  number  = {2},
  pages   = {485},
  doi     = {10.3390/ma13020485},
}
```
