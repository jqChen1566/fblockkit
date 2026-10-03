# fixtures/gaussian -- the Gaussian minimal-parser probes (menu 1's Gaussian leg)

Two real Gaussian 09 Rev D.01 outputs, generated 2018-05-30 and taken from
the gau_orca example set (a public worked example of driving Gaussian
optimizations/frequencies with an attached energy program).  Recorded here
2026-09-29 as the fixture anchor for `parsers/gaussian.py` (plan item 6.4,
the minimal fact set).  The molecules are formaldehyde-sized and carry no
sensitive content; the sessions are ordinary named output files.  Two real
**G16 Rev C.01** outputs (run 2026-10-02 on the 101 server, Gaussian 16 at
`/home/chen_jianqi/gaussian/g16`; H2O/STO-3G) supply the conventional-SCF
and CASSCF anchors.

| file | job | anchors |
|---|---|---|
| `g09_h2co_freq.out` | frequency run (external driver) | one imaginary mode at -1089.0060 cm-1 (the engine's "1 imaginary frequencies" mark agrees); `Normal termination`; the closing optimization activity follows the frequency block (after_geometry False) |
| `g09_h2co_ts.out` | TS optimization, `# opt(nomicro,calcfc,ts,noeigen)` | `ts_optimization` True from the input echo; 7 steps; `Optimization completed.` + `Stationary point found.`; no frequency block |
| `g16_h2o_rhf.log` | conventional RHF/STO-3G single point | `SCF Done:  E(RHF) =  -74.9630631539     A.U. after    7 cycles` -- the `SCF Done` regex's fixture anchor; `Normal termination` |
| `g16_h2o_cas22.log` | CAS(2,2)/STO-3G single point (`nosymm`) | the measured CASSCF facts: `no. active orbitals (n) 2` / `no. active ELECTRONS (N)= 2`, the `NBasis= 7 NCore= 4 NVal= 2 NVirt= 1` split, seven `ITN=` rows (MaxIt 64, final E -74.9643165190), `MCSCF converged.`, and the final symbolic-density diagonal 1.99789 / 0.00211 (the smaller one below the 0.02 documentation line -- a positive case for that diagnostic); `Normal termination` |

The two G09 probes are **external-driver runs**: the energy came from an
attached program, so there is **no `SCF Done` line** anywhere in them; they
anchor the termination / frequency / optimization facts *and* the
absent-field semantics (the SCF facts stay absent, and ORCA-only rules do
not fire).  The `SCF Done` regex is now **fixture-verified** by the G16
conventional-SCF probe, and the Gaussian CASSCF facts (active space,
core/valence/virtual split, per-iteration energies, `MCSCF converged.`,
active natural occupations) are emitted from the measured G16 sample --
the non-convergence spelling is unmeasured and stays unparsed (nothing is
guessed).  The G16 `.log` files are the engines' own outputs (Gaussian
writes `<name>.log` for `g16 <name>.gjf`); the routes are echoed in the
files themselves.

The parser is deliberately minimal: menu 1 (check-up and characterisation)
consumes these facts through the shared diagnosis vocabulary; the deeper
ORCA-only menus reject Gaussian files naturally (their sections are absent).

A third G16 probe (2026-10-03: `n2_cas66.log`, N2 CASSCF(6,6)/STO-3G, normal
termination) covers the **multi-block column pages** of the symbolic density
matrix: G16 pages the lower triangle in five-column blocks, and the sixth
row's first-page last value is the (6,5) element rather than its diagonal --
the fixture pins the diagonal-by-column-position read (N2 occupations
1.98270 / 0.0176416 at the extremes).  The ORCA-only energy-only-convergence
rule does not fire on this file (G16's marker states no criterion; nothing
is guessed).
