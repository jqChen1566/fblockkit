## RAS-CI core-excited XES input (two-step protocol)

RAS-CI core-excited XES input (two-step protocol; ORCA manual 3.13.19):
  step 1: valence SA-CASSCF (6e, 5o), multiplicity 5
  step 2: rotate the core orbital(s) [0, 26, 27, 28] into the window head; active (14e, 9o); one CAS-CI iteration over the single-hole-saturated space (refs ras, XESSOC rel block)
  step-2 window: orbitals 41..49 of the orbital list ((N_electrons - nel)/2 = 41)

Step-1 input written: work\fecl4.casci_xes.step1.inp
Step-2 input written: work\fecl4.casci_xes.step2.inp

Next steps:
  1. run step 1 (it writes fecl4.casci_xes.step1.gbw next to itself), then step 2 (its %moinp reads that gbw)
  2. the step-2 output carries the valence transitions and, in the SOC-corrected blocks, the emission spectra -- for the K-beta protocol (1s + 3p rotated in) the main line sits at 7086.9 eV in the probe's output (the experimental Fe K-beta1 is 7058 eV; the def2-SVP level accounts for the offset)
  3. render the emission spectrum:  orca_mapspc <step2-name>.out XESSOC -w<fwhm> -eV -n<npoints>
     (measured on the probe's 40-root output: 7375 peaks; the window is the mode's own, 5000-7200 eV, covering the K-beta region, and -x overrides did not shift it in this mode; the ABS/SOCABS modes do not read these tables -- measured refusals)

Boundaries:
  - the core set selects the edge: 1s + 3p for K-beta emission, or the 2p group for L-edge XES -- the same protocol covers both
  - the core-orbital indices are SYSTEM-SPECIFIC: read them from the step-1 output's orbital table (an L-edge 2p sits near -700 eV; a K-edge 1s near -7000 eV, a 3p near -65 eV)
  - the step-2 root count gates whether the core hole enters the listed roots (measured on the K-beta probe: 20 roots miss it, the 40,40 default reaches it in 8m44s); the saturated recipe (1000,1000 -> the restricted CSF count, 290 states here) grows superlinearly in the QDPT transition-density stage and did not finish within 24 h (measured 2026-10-03)
  - run the step-2 job in a fresh directory (measured: stale transition-density residue of a previous run under the same basename aborts the new run in TDensityContainer)
  - the manual extends the protocol with SC-NEVPT2 on the same active space (PTMethod SC_NEVPT2); the generated input stops at the CAS-CI level -- add the PT line when dynamical correlation is wanted
