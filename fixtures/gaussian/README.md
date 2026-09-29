# fixtures/gaussian -- the Gaussian minimal-parser probes (menu 1's Gaussian leg)

Two real Gaussian 09 Rev D.01 outputs, generated 2018-05-30 and taken from
the gau_orca example set (a public worked example of driving Gaussian
optimizations/frequencies with an attached energy program).  Recorded here
2026-09-29 as the fixture anchor for `parsers/gaussian.py` (plan item 6.4,
the minimal fact set).  The molecules are formaldehyde-sized and carry no
sensitive content; the sessions are ordinary named output files.

| file | job | anchors |
|---|---|---|
| `g09_h2co_freq.out` | frequency run (external driver) | one imaginary mode at -1089.0060 cm-1 (the engine's "1 imaginary frequencies" mark agrees); `Normal termination`; the closing optimization activity follows the frequency block (after_geometry False) |
| `g09_h2co_ts.out` | TS optimization, `# opt(nomicro,calcfc,ts,noeigen)` | `ts_optimization` True from the input echo; 7 steps; `Optimization completed.` + `Stationary point found.`; no frequency block |

Both are **external-driver runs**: the energy came from an attached program,
so there is **no `SCF Done` line** anywhere in them.  They therefore anchor
the termination / frequency / optimization facts *and* the absent-field
semantics (the SCF facts stay absent, and ORCA-only rules do not fire).
The `SCF Done` regex in the parser follows the published line shape; a
conventional-SCF fixture anchor is **registered as a follow-up** -- do not
treat that regex as fixture-verified until such a sample is added.  Gaussian
CASSCF facts are registered likewise (no sample in hand; nothing guessed).

The parser is deliberately minimal: menu 1 (check-up and characterisation)
consumes these facts through the shared diagnosis vocabulary; the deeper
ORCA-only menus reject Gaussian files naturally (their sections are absent).
