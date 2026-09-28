# SINGLE_ANISO fixtures (fixtures/single_aniso/)

Real ORCA 6.1.1 outputs carrying the SINGLE_ANISO embedded section (the
`ANISO` sub-block of `%casscf`), read by `src/fblockkit/parsers/single_aniso.py`
and `src/fblockkit/analysis/relaxation.py` (menu 36).  All runs on 101,
2026-09-28, `def2-SVP`, `TightSCF`, terminated normally.

## The three samples

- **o2_aniso.{inp,out}** -- O2, CAS(12,8), `mult 3,1`, three roots each,
  `DoSOC` + `DoSSC`, `ANISO doaniso true MLTP 3,3,3,3 UBAR true TINT/HINT/
  TMAG`.  An **integer pseudospin** sample (S = 1, four groups) and the
  **two-segment** case: with `DoSSC` the output holds two complete
  SINGLE_ANISO runs, the first on the SOC-only QDPT data (segment D =
  2.175287) and the second after the spin-spin pass (D = 3.185964); the
  intermediate `QDPT WITH CASSCF DIAGONAL ENERGIES` /
  `Calculating Spin-Spin Coupling Integrals` sections separate them.  The
  Zeeman-eigenstate table carries the `m = 0` column for each group.
- **co_aniso2.{inp,out}** -- CO+, CAS(9,8), `mult 2`, two roots,
  `MLTP 2,2`: two Kramers doublets; a Kramers doublet prints **no ZFS/D
  block** (the group block goes g table -> ITO decomposition -> angular
  moments).  UBAR prints `m+, m-, E` columns.
- **co_aniso.{inp,out}** -- the same CO+ input **without MLTP**: the header
  reads `Computation of the EPR g-tensor ... NOT INCLUDED` and
  `Estimation of blocking barrier for SMMs ... NOT INCLUDED`, and no
  pseudospin group block is printed at all.  Measured lesson: give `MLTP`
  explicitly (the manual's "default = ground term multiplicity" wording
  notwithstanding).

## Measured format points

- Each segment: program banner -> `Information about this calculation`
  (data file name, number of spin-orbit states, the **per-state SOC energy
  spectrum**, spin-free spectrum, the g-tensor/UBAR inclusion echoes) ->
  `EXPECTATION VALUES` (M, S, L, J moments per SOC state) -> one block per
  pseudospin group -> `AB INITIO BLOCKING BARRIER` (UBAR: the qualitative
  statement, Zeeman eigenstates, intra-group "opposite magnetization" and
  inter-group "neighbouring multiplets" magnetic-moment matrix elements,
  averaged as printed `(|mu_X| + |mu_Y| + |mu_Z|)/3`) ->
  `CALCULATION OF THE MOLAR MAGNETIZATION` -> `HAPPY LANDING`.
- The `in cases with even number of electrons ... check the tunnelling
  splitting instead` sentence is a **fixed template** warning: the CO+
  sample is odd-electron and prints it too -- record, never interpret.
- The engine's matrix-element row labels repeat (`mu_Y` printed for the
  third component); the component names are therefore taken in row order
  and the printed average is the value to use.
- Cross-check anchors: segment-1 D = 2.175287 and segment-2 D = 3.185964
  equal the QDPT ZFS variants that ORCA's own QDPT driver prints in
  `fixtures/orca/o2_qdpt2.out` (`EFFECTIVE HAMILTONIAN SOC CONTRIBUTION`
  and `... SOC and SSC CONTRIBUTION`).  The CO+ groups are nearly isotropic
  (Delta g / g ~ 2.4e-4): their main-axis directions are set by the
  printed digits -- the measured case behind the analysis module's 1 %
  axis-definition guard.
- For O2's group 2 the printed g values really are (0, 0, 2): the transverse
  components are zero to the print (recorded as printed).

## Regenerating

The inputs are stored verbatim; rerunning on ORCA 6.1.1 with the same
versions reproduces the outputs (the CO+ samples reuse the converged
CASSCF gbw via `%moinp`; the `co_aniso2` run reads `co_aniso.gbw`).
