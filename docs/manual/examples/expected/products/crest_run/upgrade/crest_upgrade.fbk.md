## CREST conformer upgrade (upgrade inputs for the ensemble)

CREST conformer upgrade (one ORCA optimisation input per conformer):
  ensemble conformers upgraded: 2
  method line: ! Opt r2SCAN-3c   charge 0, multiplicity 1
  directory: work\crest_run\upgrade

  input                ensemble E (Eh)      relative (kcal/mol)
  conf_01.opt.inp            -13.66512758        0.000
  conf_02.opt.inp            -13.66417737        0.596

Next steps:
  1. run one input after the other (or submit the directory as a batch): orca conf_01.opt.inp > conf_01.opt.out
  2. read each output back with menu 1 (the check-up report carries the final geometry and energy)

Boundaries:
  - the geometries are the pre-screen ensemble's, used as-is; the selection order is the ensemble's own (energy-ascending), so a count takes the lowest conformers first
  - the upgrade re-ranks nothing by itself: comparing the upgraded energies and re-weighting the ensemble is registered as the next increment
