# MOKIT fixtures (plan item 6.1; the reader behind menu 45, the generator behind menu 44)

All files are real MOKIT 1.2.8 output (101, 2026-09-30). Deployment and
environment (the full chain, no sudo):

- MOKIT: the official prebuilt `py310_gcc10` artifact from GitLab
  (`mokit-master_linux_py310_gcc10`), unpacked on the 101 server
  (a MOKIT v1.28 directory under /data1/cjq12/). The build's extension modules are
  CPython-3.10 and compiled against NumPy 1.x.
- Python: an isolated venv `mokitenv` (built from a Python 3.10 that exists
  on the machine) with **numpy 1.26.4** + pyscf 2.6.2 + scipy 1.13.1 +
  h5py 3.11.0 installed from local wheels. Two measured traps: an
  anaconda 3.13 python cannot import the extension modules
  (`No module named 'mokit.lib.py2fch'`), and a Python 3.10 with NumPy 2
  fails at import (`numpy.core.multiarray failed to import`).
- Backends: GAMESS is not installed, so GVB ran with the Gaussian 16
  backend (`mokit{GVB_prog=Gaussian}`; Gaussian 16 at
  `/home/chen_jianqi/gaussian/g16`), and the CASSCF stage with its PySCF
  default. Environment used for every run:

```
export MOKIT_ROOT=<the unpacked MOKIT>
export PATH=$MOKIT_ROOT/bin:$MOKIT_ROOT/../mokitenv/bin:$PATH
export PYTHONPATH=$MOKIT_ROOT
export g16root=/home/chen_jianqi/gaussian
export GAUSS_SCRDIR=<workdir>/scratch
automr <name>.gjf > <name>.out
```

Files:

- `h2o_gvb.gjf` / `h2o_gvb_automr.out` — the hand-written probe
  (CASSCF/cc-pVDZ, `mokit{GVB_prog=Gaussian}`) and its run: RHF/UHF →
  UNO rotation → GVB(4) by Gaussian → CASSCF(4e,4o) by PySCF, three
  Radical-index tables, `Normal termination of AutoMR`;
- `h2o_generated.gjf` — the input produced by **menu 44** for the same
  geometry (byte-reproducible from the generator; the test asserts this);
- `h2o_generated_automr.out` — the acceptance run of that generated input
  (same chain, CAS(4e,4o) automatically determined, Normal termination):
  this is the generator's verification record.
- `h2o_gvb_rhf.fch`, `h2o_gvb_uhf_gvb4_CASSCF_NO.fch`,
  `h2o_gvb_uhf_uno_asrot2gvb4.fch` — three `.fch` side products of the same
  probe chain (the RHF stage, the CASSCF natural orbitals, and the UNO
  active-space rotation towards GVB(4)), read by `parsers/mokit_fch.py`
  (menu 45's side-product section).  Measured format facts: the MO
  coefficient arrays are column-major over (basis, MO); the density
  triangles are row-major lower (the RHF file's five closed-shell MOs
  reconstruct its own `Total SCF Density` to 9e-9 — the reader's regression
  anchor); coordinates are in Bohr; the transformed-stage files carry alpha
  sections only (no beta block); the energy scalars are the stage's
  SCF/total values (-75.78429 Eh RHF, -75.91807 Eh GVB stages).

Measured boundary (kept out of the fixture set): an N2 probe with the
active-space size unpinned made automr choose GVB(5), and Gaussian 16
`l506.exe` crashed on that step on this machine (stack dump, no
`Normal termination`); the record stays in the project scratch notes. The
shipped probes are the two water runs that ran to completion.
