# CREST fixtures (plan item 6.1; the reader behind menu 43)

The files are a real CREST 3.0.2 run (101, 2026-09-30, n-butane on GFN2;
CREST deployed from the official `crest-gnu-12-ubuntu-latest.tar.xz`
release, LGPL-3.0). Sequence (run in this directory; the geometric input
was pre-optimised because CREST aborts on a topology change of an
unpreoptimised start -- see below):

```
xtb butane.xyz --gfn2 --opt          # pre-optimisation (butane_clean.xyz)
crest butane_clean.xyz --gfn2 --quick -T 4 > crest_butane_q4.out 2>&1
```

Files:

- `butane_clean.xyz` — the pre-optimised input geometry;
- `crest_butane_q4.out` — the console capture of the quick iMTD-GC search
  (sampling banner, MTD blocks, CREGEN de-duplication lines with the lowest
  energies, `CREST terminated normally.` and the engine's call counter);
- `crest_conformers.xyz` — the final ensemble, two conformers (the
  anti/gauche pair); atom-count lines carry leading spaces, the comment
  line is the absolute energy in hartree in a wide right-aligned field;
- `crest.energies` — the relative-energy table
  (`<index>  <kcal/mol>` with three decimals: 0.000, 0.596);
- `crest_best.xyz` — the single lowest conformer.

Measured boundary (kept out of the fixture set): an **unpreoptimised**
input whose initial optimisation changes the bonding topology aborts with
`*WARNING* Change in topology detected!` and offers options A/B/C (read
`crestopt.log`, `--noreftopo`, or a fixed-topology pre-optimisation); with
stdin at /dev/null the run stops there and leaves no ensemble.
