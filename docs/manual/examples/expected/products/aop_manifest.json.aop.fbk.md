## G7 AOP rotation guess (a target orbital set aligned onto a reference active space)

AOP rotation guess (two-step SVD rotation of the target's orbitals onto a reference active space):
  reference: n2_cas666_1.600
  target:    n2_scan_1.600
  partition at the target: closed 4, active 6, virtual 18

  reference containment in the window: 0.958044
  O_min of the built active block:    1.000000
  orthonormality residual:            6.44e-15

Criteria and boundaries:
  - the two-step rotation follows the source's Eqs. (4)-(11): the first SVD rotates the lowest closed+active target orbitals (Eq. (6)), the second concatenates the rotated active block with the untouched virtuals (Eq. (7)) and rotates again; the final set takes the closed block from the first step and the active and virtual blocks from the second (Eqs. (8)-(11))
  - the overlap block uses the *target's* S matrix -- the source's own setting (chromophore orbitals inside a same-geometry condensed-phase basis); between distant geometries it is the small-step approximation the source uses along its interpolations, so the alignment there is a demonstration, not a calibrated reading
  - the reference's representation in the target's closed+active window (containment) is 0.958044: the smallest singular value of the reference-to-window overlap; the gate requires the source's own alignment line (0.85), below which the closed block the construction would build is not the closed manifold (measured: a non-corresponding cross-geometry reference reads 0.022 and sent the CASSCF to a wrong solution)
  - O_min = 1.000000 is the smallest singular value of the built active block's overlap with the reference -- the source's alignment diagnostic (its Eq. (3)); its working criterion O_min >= 0.85 is calibrated on its own condensed-phase dataset and read there off *converged* orbitals, while this value belongs to the guess (measured: with a same-geometry reference the built block reproduces the reference exactly, O_min = 1.000000)
  - the guess is a unitary rotation of the target's own orbital set; the orthonormality residual max |C' S C'^T - I| in the target's metric is 6.44e-15, and the occupation/energy tags written with it are the target's
  - Within a degenerate reference window the individual orbital assignment is arbitrary (any rotation among degenerate reference orbitals aligns equally well); the SVD fixes the subspace, which is what the O_min diagnostic reads
  - The written file is a Molekel mkl: run ``orca_2mkl <name>.fbk -gbw`` to turn it into a ``.gbw`` and read it with ``!moread`` + ``%moinp`` (the route measured for menu 18).
  - This module builds the guess; ORCA re-optimises the orbitals. The rotated set is not a converged orbital set and is not claimed to be one.

## References

Complete citations:
- [paz2021active] Paz, A. S. P.; Baleeva, N. S.; Glover, W. J. (2021). Active orbital preservation for multiconfigurational self-consistent field. The Journal of Chemical Physics, 155(7), 071103. DOI: 10.1063/5.0058673

BibTeX (paste-ready):

```bibtex
@article{paz2021active,
  author  = {Paz, Amiel S. P. and Baleeva, Nadezhda S. and Glover, William J.},
  title   = {Active orbital preservation for multiconfigurational self-consistent field},
  journal = {The Journal of Chemical Physics},
  year    = {2021},
  volume  = {155},
  number  = {7},
  pages   = {071103},
  doi     = {10.1063/5.0058673},
}
```
