## ROCIS XES input (off-resonance X-ray emission)

ROCIS XES input (off-resonance X-ray emission):
  charge -2, multiplicity 5 (also ReferenceMult)
  XASelems 0 (the 0-based position of the core element among the atoms)
  NRoots 30; OrbWin 6,6, 7,8, 0,2000 (donor 1, donor 2, acceptor)
  plain RIXS channel: requested (it carries the off-resonance XES table)
  SOC-corrected RIXS channel: off
  elastic line: included
  SCF preparation: ROHF high-spin, ROHF_NEL[1] = 4

Input written: work\fecl4.xes.inp

Next steps:
  1. run:  orca <name>.inp > <name>.out
  2. the run's emission tables carry the ground-state rows ('<root>-<mult>A -> 0-<mult>A') at the core-emission energies (they appear at the end of the emission tables, after the between-excited-state rows)
  3. render the spectrum, e.g.:  orca_mapspc <name>.out XES -x0<lo> -x1<hi> -w<fwhm> -eV -n<npoints>
     -> <name>.out.XES.stk (the sticks) and .XES.dat (the broadened curve)

Boundaries:
  - the off-resonance XES is automatic in a RIXS-requested run (the manual's sentence, measured on the probe); no dedicated keyword is needed
  - the plain channel needs enough roots: too few and the engine skips it with its own warning, printing no XES table (measured: 10 skipped, 30 covered, on the probe)
  - the SOC-corrected channel stores per-pair transition densities (measured: the transient storage grew past 36 GB at about 1 GB/min on the probe at NRoots 30) -- requested only for small systems or small root sets
  - the KHD variant (DoKHDXESSOC and its state lists) is a documented termination: every probed input computes and then aborts in the engine's own printout (measured)
  - the window indices are the SYSTEM's own orbital numbers (the probe's 6,6,7,8,0,2000 are not transferable); read the two spin-orbit-split core ranges and a wide acceptor from the target run's own orbital table
