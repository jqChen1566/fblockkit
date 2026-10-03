## CAS-CI core-excited XAS input (two-step protocol)

CAS-CI core-excited XAS input (two-step protocol; ORCA manual 3.13.18):
  step 1: valence SA-CASSCF (6e, 5o), multiplicity 5
  step 2: rotate the core orbital(s) [6, 7, 8] into the window head; active (12e, 8o); one CAS-CI iteration over the core-saturated space
  step-2 window: orbitals 42..49 of the orbital list ((N_electrons - nel)/2 = 42)

Step-1 input written: work\fecl4.casci_xas.step1.inp
Step-2 input written: work\fecl4.casci_xas.step2.inp

Next steps:
  1. run step 1 (it writes fecl4.casci_xas.step1.gbw next to itself), then step 2 (its %moinp reads that gbw)
  2. the step-2 output carries the valence transitions and the core transitions -- the core-excited CAS-CI roots sit far above the ground state (the probe's Fe 2p to 3d L-edge states at 719.36 eV)
  3. render the spectrum:  orca_mapspc <step2-name>.out SOCABS -x0<lo> -x1<hi> -w<fwhm> -eV -n<npoints>
     for the SOC-corrected table (measured: 785 peaks on the probe), or the ABS mode for the plain one (19); the XAS/XASSOC modes do not read these tables (measured refusals)

Boundaries:
  - the window selection is positional (measured): with the unfrozen core the active window is the norb consecutive orbitals starting at (N_electrons - nel)/2 -- the manual's own example lands on 87 for Fe(acac)3, the probe on 42 for [FeCl4]2-; the rotation targets are the leading slots of that window
  - the core-orbital indices are SYSTEM-SPECIFIC: read them from the step-1 output's orbital table (an L-edge 2p sits near -700 eV; a K-edge 1s near -7000 eV in light elements)
  - enough step-2 roots must be requested for the core manifolds (the probe: nroots 20,20 reached the L3/L2 region); the manual's Fe(acac)3 example uses 16,173
  - the manual extends the protocol with SC-NEVPT2 on the same active space (PTMethod SC_NEVPT2); the generated input stops at the CAS-CI level -- add the PT line when dynamical correlation is wanted
  - cross-engine anchor (measured): the probe's first core excitation (719.36 eV) sits 0.5 eV from the ROCIS result of the same system (718.87 eV) -- two methods, one edge
