## MOKIT automr run report (B-layer reading)

MOKIT automr run report
Source: h2o_generated_automr.out

Version: 1.2.8 (2026-Aug-18)
Program paths: bdf=NOT FOUND, dalton=NOT FOUND, gau=/home/chen_jianqi/gaussian/g16/g16, gms=NOT FOUND, molcas=NOT FOUND, molpro=NOT FOUND, orca=/home/chen_jianqi/orca/orca_6_1_1/orca, psi4=NOT FOUND
Run settings: memory 4GB, nproc 2, method/basis casscf/cc-pvdz
MOKIT options (merged): gvb_prog=gaussian
Strategy (No. 1) active flags: UNO, GVB, CASSCF
Stages: HF -> UNO rotation -> GVB -> CASSCF/CASCI

Energy chain:
  RHF     -75.78429273 Eh  <S**2> = 0.000
  UHF     -75.82292560 Eh  <S**2> = 1.066
  UHF2    -75.82292560 Eh  <S**2> = 1.066
  GVB     -75.91806514 Eh
  CASCI   -75.90802996 Eh
  CASSCF  -75.90823052 Eh
GVB: 4 pair(s), program gaussian
Active space (automatically determined): CAS(4e,4o), program pyscf

Radical index (table 1 of 3):
  biradical character y0: 0.109
  tetraradical character y1: 0.041
  Yamaguchi unpaired electrons: 2.132
  Head-Gordon unpaired (min): 1.279
  Head-Gordon unpaired (squared): 1.165

Radical index (table 2 of 3):
  biradical character y0: 0.154
  tetraradical character y1: 0.154
  Yamaguchi unpaired electrons: 1.181
  Head-Gordon unpaired (min): 0.638
  Head-Gordon unpaired (squared): 0.325

Radical index (table 3 of 3):
  biradical character y0: 0.165
  tetraradical character y1: 0.139
  Yamaguchi unpaired electrons: 1.121
  Head-Gordon unpaired (min): 0.607
  Head-Gordon unpaired (squared): 0.315

Side products echoed by the run (audit trail, not parsed here):
  mokit_gen2_uhf_uno_asrot2gvb4.gjf
  mokit_gen2_uhf_uno_asrot2gvb4.dat
  mokit_gen2_uhf_uno_asrot2gvb4_s.dat
  mokit_gen2_uhf_uno_asrot2gvb4_s.fch

Termination: Normal termination of AutoMR -- the workflow ran to its end.

Reading notes: the active space on the CASSCF line is the automatically determined selection (GVB natural-orbital occupations above the 0.02 threshold); writing CASSCF(n,m) in the route line pins the size instead. The energy chain shows the flow RHF/UHF (the lower one is kept) -> GVB -> CASCI/CASSCF, so consecutive entries are different wave-function levels, not an error. GVB runs need a backend program (GAMESS by default; Gaussian and QChem are the alternates -- PySCF is not a GVB backend); the CASSCF stage defaults to PySCF. The natural-orbital .fch side products are listed for the audit trail but this report does not parse them (stated boundary); dynamic-correlation stages (CASPT2/NEVPT2/DMRG) appear in the strategy table only when requested.
