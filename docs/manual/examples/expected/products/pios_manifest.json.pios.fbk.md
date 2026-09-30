## A13 PiOS pi-orbital active space (oriented-p projection)

PiOS pi-orbital active space (Huckel-style construction on the oriented valence p orbitals):
  system: benzene_rhf
  atom    element   sigma-bonds   pi electrons
     0    C              3            1
     1    C              3            1
     2    C              3            1
     3    C              3            1
     4    C              3            1
     5    C              3            1

  pi plane normal:       (+0.000000, +0.000000, +1.000000)
  centroid:              (-0.000000, +0.000000, +0.000000)
  max out-of-plane:      0.000000 Angstrom
  pi electrons:          6 (charge +0; 6 from the atoms)
  active space:          CAS(6e, 6o) = 3 occupied + 3 virtual pi orbitals

Occupied projection spectrum (selected first):
   orbital      sigma
         0   0.778934  <- selected
         1   0.764871  <- selected
         2   0.764871  <- selected
         3   0.000000
         4   0.000000
         5   0.000000
         6   0.000000
         7   0.000000
         8   0.000000
         9   0.000000
        10   0.000000
        11   0.000000
        12   0.000000
        13   0.000000
        14   0.000000
        15   0.000000
        16   0.000000
        17   0.000000
        18   0.000000
        19   0.000000
        20   0.000000

Virtual projection spectrum (selected first):
   orbital      sigma
         0   1.000000  <- selected
         1   1.000000  <- selected
         2   1.000000  <- selected
         3   0.235129
         4   0.235129
         5   0.221066
         6   0.000000
         7   0.000000
         8   0.000000
         9   0.000000
        10   0.000000
        11   0.000000
        12   0.000000
        13   0.000000
        14   0.000000
        15   0.000000
        16   0.000000
        17   0.000000
        18   0.000000
        19   0.000000
        20   0.000000
        21   0.000000
        22   0.000000
        23   0.000000
        24   0.000000
        25   0.000000
        26   0.000000
        27   0.000000
        28   0.000000
        29   0.000000
        30   0.000000
        31   0.000000
        32   0.000000
        33   0.000000
        34   0.000000
        35   0.000000
        36   0.000000
        37   0.000000
        38   0.000000
        39   0.000000
        40   0.000000
        41   0.000000
        42   0.000000
        43   0.000000
        44   0.000000
        45   0.000000
        46   0.000000
        47   0.000000
        48   0.000000
        49   0.000000
        50   0.000000
        51   0.000000
        52   0.000000
        53   0.000000
        54   0.000000
        55   0.000000
        56   0.000000
        57   0.000000
        58   0.000000
        59   0.000000
        60   0.000000
        61   0.000000
        62   0.000000
        63   0.000000
        64   0.000000
        65   0.000000
        66   0.000000
        67   0.000000
        68   0.000000
        69   0.000000
        70   0.000000
        71   0.000000
        72   0.000000
        73   0.000000
        74   0.000000
        75   0.000000
        76   0.000000
        77   0.000000
        78   0.000000
        79   0.000000
        80   0.000000
        81   0.000000
        82   0.000000
        83   0.000000
        84   0.000000
        85   0.000000
        86   0.000000
        87   0.000000
        88   0.000000
        89   0.000000
        90   0.000000
        91   0.000000
        92   0.000000

Selected pi orbitals (semicanonical energies, parent SCF orbitals):
   occ 0:  E = -0.498453 Eh   parents: 16 (1.000), 17 (0.000)
   occ 1:  E = -0.333508 Eh   parents: 19 (1.000), 17 (0.000)
   occ 2:  E = -0.333508 Eh   parents: 20 (1.000), 13 (0.000)
   vir 0:  E = +0.281786 Eh   parents: 1 (0.790), 26 (0.208)
   vir 1:  E = +0.281786 Eh   parents: 0 (0.790), 25 (0.208)
   vir 2:  E = +0.470139 Eh   parents: 8 (0.789), 31 (0.209)

Criteria and boundaries:
  - the construction follows the source's Eqs. (1)-(16): the plane from the inertia tensor, the oriented valence p orbitals, the occupied/virtual projection spectra, the selection of N_pi,occ = N_pi,e / 2 and its virtual partner, and the Fock semicanonicalisation of the selected blocks
  - the oriented orbitals live on each atom's valence p shell (the outermost p shell of the export; used shells: 2, 2, 2, 2, 2, 2); the source builds them in an auxiliary MINAO basis through IAOs, which needs the cross-basis overlap ORCA does not export -- this is the calculation's own basis, the same substitution menu 14 carries
  - the Fock matrix is assembled as H + J + K from the export's blocks; the written set is orthonormal to 1.1e-13 in the export's overlap
  - the written basis of every block is the Gram-Schmidt of the parent SCF orbitals projected onto it (pi character decides the partition; index order within each pool), so the files do not inherit the diagonaliser's BLAS-dependent choice of a degenerate cluster's basis (a determinism convention, not part of the source)
  - The written file is a Molekel mkl in ORCA's partition order -- inactive occupied | pi occupied | pi virtual | inactive virtual -- so a %casscf run with nel 6, norb 6 takes the pi space as its active window: run ``orca_2mkl <name>.fbk -gbw`` and start with ``!moread`` + ``%moinp`` (the menu 47 route)
  - Scope: one pi system per call and closed-shell exports; the projection eigenvalues in the report are the source's own validity measure for the chosen atoms and electron count (a selected eigenvalue below an excluded one means the choices deserve a second look)

## References

Complete citations:
- [sayfutyarova2019constructing] Sayfutyarova, E. R.; Hammes-Schiffer, S. (2019). Constructing Molecular pi-Orbital Active Spaces for Multireference Calculations of Conjugated Systems. Journal of Chemical Theory and Computation, 15(3), 1679-1689. DOI: 10.1021/acs.jctc.8b01196

BibTeX (paste-ready):

```bibtex
@article{sayfutyarova2019constructing,
  author  = {Sayfutyarova, Elvira R. and Hammes-Schiffer, Sharon},
  title   = {Constructing Molecular pi-Orbital Active Spaces for Multireference Calculations of Conjugated Systems},
  journal = {Journal of Chemical Theory and Computation},
  year    = {2019},
  volume  = {15},
  number  = {3},
  pages   = {1679--1689},
  doi     = {10.1021/acs.jctc.8b01196},
}
```
