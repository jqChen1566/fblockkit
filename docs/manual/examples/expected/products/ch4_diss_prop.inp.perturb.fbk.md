## 4.1 perturbed multistart batch (randomized occupied-virtual mixing)

perturbed start 1 (seed 20260927):
  occupied window (mkl columns): 0, 1, 2, 3, 4
  virtual window  (mkl columns): 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19
  alpha mixing (10 pair(s)): 1<-12 28.19 deg; 0<-6 6.98 deg; 3<-10 61.35 deg; 2<-6 16.65 deg; 3<-9 22.86 deg; 1<-10 74.46 deg; 4<-8 57.31 deg; 4<-7 46.87 deg; 0<-9 18.41 deg; 3<-6 38.15 deg
  beta  mixing (10 pair(s)): 0<-19 55.09 deg; 3<-5 11.50 deg; 2<-14 28.49 deg; 1<-12 12.42 deg; 1<-7 9.50 deg; 0<-5 14.36 deg; 4<-6 4.22 deg; 3<-6 40.13 deg; 0<-15 27.69 deg; 1<-18 41.64 deg
perturbed start 2 (seed 20260928):
  occupied window (mkl columns): 0, 1, 2, 3, 4
  virtual window  (mkl columns): 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19
  alpha mixing (10 pair(s)): 4<-19 25.22 deg; 0<-8 85.25 deg; 3<-10 82.27 deg; 1<-16 30.19 deg; 0<-12 2.62 deg; 3<-15 51.84 deg; 4<-18 80.06 deg; 0<-5 82.34 deg; 2<-13 2.25 deg; 3<-15 52.94 deg
  beta  mixing (10 pair(s)): 2<-10 15.56 deg; 1<-18 86.00 deg; 1<-8 87.72 deg; 0<-8 15.05 deg; 2<-9 38.27 deg; 0<-19 70.75 deg; 2<-13 53.19 deg; 0<-8 50.94 deg; 1<-19 0.06 deg; 3<-19 65.43 deg
perturbed start 3 (seed 20260929):
  occupied window (mkl columns): 0, 1, 2, 3, 4
  virtual window  (mkl columns): 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19
  alpha mixing (10 pair(s)): 3<-13 66.22 deg; 3<-7 39.25 deg; 4<-10 43.41 deg; 0<-13 0.82 deg; 4<-14 53.31 deg; 3<-12 89.42 deg; 4<-12 60.53 deg; 0<-19 16.73 deg; 3<-7 16.87 deg; 2<-12 81.54 deg
  beta  mixing (10 pair(s)): 0<-7 45.24 deg; 4<-12 75.22 deg; 2<-16 26.53 deg; 3<-12 17.02 deg; 0<-14 33.72 deg; 3<-9 68.28 deg; 4<-8 30.01 deg; 3<-16 87.32 deg; 1<-5 87.25 deg; 1<-7 63.97 deg

Next steps:
  - convert with ``orca_2mkl <name> -gbw`` and run: ch4_diss_prop.p1.fbk, ch4_diss_prop.p2.fbk, ch4_diss_prop.p3.fbk
  - a lower converged energy than the reference's identifies wrong convergence
  - the perturbation cannot guarantee detection (the source's boundary)

## References

Complete citations:
- [vaucher2017steering] Vaucher, A. C.; Reiher, M. (2017). Steering Orbital Optimization out of Local Minima and Saddle Points Toward Lower Energy. Journal of Chemical Theory and Computation, 13(3), 1219-1228. DOI: 10.1021/acs.jctc.7b00011

BibTeX (paste-ready):

```bibtex
@article{vaucher2017steering,
  author  = {Vaucher, Alain C. and Reiher, Markus},
  title   = {Steering Orbital Optimization out of Local Minima and Saddle Points Toward Lower Energy},
  journal = {Journal of Chemical Theory and Computation},
  year    = {2017},
  volume  = {13},
  number  = {3},
  pages   = {1219--1228},
  doi     = {10.1021/acs.jctc.7b00011},
}
```
