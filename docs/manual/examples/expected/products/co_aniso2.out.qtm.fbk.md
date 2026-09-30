## Quantum-tunnelling relaxation prediction

Quantum-tunnelling relaxation prediction
Source: co_aniso2.out (segment 1 of 1)
Model: equivalent-Zeeman (Yin & Li 2020), B_ave = 20.0 mT

Kramers doublets (axial g = the largest principal value):
   KD      E (cm-1)          g_xy          g_z      tau_QTM (s)
    1       0.000     2.82808773    1.99999879   1.547e-09 s (= 10^-8.81)
    2   29361.197     2.82842541    2.00047880   1.547e-09 s (= 10^-8.81)

Ground-state tau_QTM: 1.547e-09 s (= 10^-8.81)

U_eff(T) (thermal-activation weighting over the doublets):
     T (K)   U_eff (cm-1)   contributions (KD: weight)
      2.0         0.000   1: 1.000, 2: 0.000
      4.0         0.000   1: 1.000, 2: 0.000
      6.0         0.000   1: 1.000, 2: 0.000
      8.0         0.000   1: 1.000, 2: 0.000
     10.0         0.000   1: 1.000, 2: 0.000
     15.0         0.000   1: 1.000, 2: 0.000
     20.0         0.000   1: 1.000, 2: 0.000
     30.0         0.000   1: 1.000, 2: 0.000
     40.0         0.000   1: 1.000, 2: 0.000
     60.0         0.000   1: 1.000, 2: 0.000
     80.0         0.000   1: 1.000, 2: 0.000
    100.0         0.000   1: 1.000, 2: 0.000
    150.0         0.000   1: 1.000, 2: 0.000
    200.0         0.000   1: 1.000, 2: 0.000
    300.0         0.000   1: 1.000, 2: 0.000
  (the first excited doublet lies at 29361.2 cm-1 = 42244 K, far above this temperature grid: the thermally activated regime of the model is not reached in the printed window, and U_eff stays at the ground-state value)

Spin-dipolar model (Aravena 2018; non-collinear form 2026)
Neighbour table: qtm_neighbours.txt (12 centres)
  sigma_r = 2.519437e-25  sigma_i = 2.518835e-25  (J)
  <|E_sf|> = 2.842540e-25 J
  tau_QT = 1.434e-09 s (= 10^-8.84)
  Dilution variant (median over repeats, seeded):
        x     log10(tau/s)   active centres (median)
     1.000       -8.843   12
     0.300       -8.366   4
     0.050       -7.764   1

Reading notes: the models are zero-field, Kramers-ion, single-centre descriptions of the tunnelling (QTM) relaxation; B_ave is an empirical field scale (20 mT in the source, adjustable -- the paper states it can be treated as an empirical parameter). The low-temperature U_eff tends to zero (the ground doublet does not contribute to the barrier) and rises to the Orbach plateau. These predictions are absolute-value estimates from ab initio parameters; menu 36's reading of an experimental-style plateau is a separate, data-side quantity -- keep the two apart when quoting them side by side. The dipolar model needs the neighbour geometry as a table (crystal-structure parsing sits outside this menu's scope).

## References

Complete citations:
- [yin2020qtm] Yin, B.; Li, C. (2020). A Method to Predict Both the Relaxation Time of Quantum Tunneling of Magnetization and the Effective Barrier of Magnetic Reversal for a Kramers Single-Ion Magnet. Physical Chemistry Chemical Physics, 22(18), 9923-9933. DOI: 10.1039/D0CP00933D
- [aravena2018tunneling] Aravena, D. (2018). Ab Initio Prediction of Tunneling Relaxation Times and Effective Demagnetization Barriers in Kramers Lanthanide Single-Molecule Magnets. The Journal of Physical Chemistry Letters, 9(18), 5327-5333. DOI: 10.1021/acs.jpclett.8b02359
- [aravena2026packing] Aravena, D. (2026). Lanthanide Single Molecule Magnets: Relation Between Crystal Packing and Tunnelling Relaxation Time. Dalton Transactions, 55, 7848-7857. DOI: 10.1039/d6dt00521g
- [llanos2019dilution] Llanos, L.; Aravena, D. (2019). Relaxation Time Enhancement by Magnetic Dilution in Single-Molecule Magnets: An Ab Initio Study. Journal of Magnetism and Magnetic Materials, 489, 165456. DOI: 10.1016/j.jmmm.2019.165456

BibTeX (paste-ready):

```bibtex
@article{yin2020qtm,
  author  = {Yin, Bing and Li, Chao-Chao},
  title   = {A Method to Predict Both the Relaxation Time of Quantum Tunneling of Magnetization and the Effective Barrier of Magnetic Reversal for a Kramers Single-Ion Magnet},
  journal = {Physical Chemistry Chemical Physics},
  year    = {2020},
  volume  = {22},
  number  = {18},
  pages   = {9923--9933},
  doi     = {10.1039/D0CP00933D},
}
```

```bibtex
@article{aravena2018tunneling,
  author  = {Aravena, Daniel},
  title   = {Ab Initio Prediction of Tunneling Relaxation Times and Effective Demagnetization Barriers in Kramers Lanthanide Single-Molecule Magnets},
  journal = {The Journal of Physical Chemistry Letters},
  year    = {2018},
  volume  = {9},
  number  = {18},
  pages   = {5327--5333},
  doi     = {10.1021/acs.jpclett.8b02359},
}
```

```bibtex
@article{aravena2026packing,
  author  = {Aravena, Daniel},
  title   = {Lanthanide Single Molecule Magnets: Relation Between Crystal Packing and Tunnelling Relaxation Time},
  journal = {Dalton Transactions},
  year    = {2026},
  volume  = {55},
  pages   = {7848--7857},
  doi     = {10.1039/d6dt00521g},
}
```

```bibtex
@article{llanos2019dilution,
  author  = {Llanos, Leonel and Aravena, Daniel},
  title   = {Relaxation Time Enhancement by Magnetic Dilution in Single-Molecule Magnets: An Ab Initio Study},
  journal = {Journal of Magnetism and Magnetic Materials},
  year    = {2019},
  volume  = {489},
  pages   = {165456},
  doi     = {10.1016/j.jmmm.2019.165456},
}
```
