# Judd-Ofelt fixtures (fixtures/judd_ofelt/)

The regression dataset of the Judd-Ofelt fit (menu 34;
`src/fblockkit/analysis/judd_ofelt.py`).

## babu2000_eu3.yaml

Eu3+ in lithium fluoroborate glass, the standard benchmark of the modern
Judd-Ofelt literature:

- **oscillator strengths and wavelengths**: Babu et al., Physica B 279, 262
  (2000), Table 3, **SET B** (thermal-corrected) -- transcribed from the file
  the authors of the extended JO theory distribute with their reference code
  (jo_so, GPL-3.0, `gitlab.com/labicb/joso`, `Eu3+/Eu3+_spec_babu2000.txt`;
  data only, no code; read 2026-09-28);
- **[U^(lambda)]^2 table**: Hovhannesyan, Boudon, Lepers, J. Lumin. 266
  (2024) 120234, Table 1 (where each value is cross-checked against Carnall
  et al. 1968);
- **host**: Sellmeier n0 = 1, A = 1.242781, B = 0.023833 um^2 (Adamiv et al.
  2011, via the same paper).

Published values this file regresses against (same paper / its predecessor
J. Lumin. 241 (2022) 118456):

| quantity | published | this toolkit |
|---|---|---|
| relative sigma (standard JO, 9 transitions) | 8.52 % | 8.52 % |
| Omega_6 | 2.253e-20 cm^2 | 2.253e-20 cm^2 |
| Omega_2 | 18.73e-20 cm^2 | 18.43e-20 cm^2 (1.6 %) |
| Omega_4 | 12.58e-20 cm^2 | 11.02e-20 cm^2 (12 %, documented band) |
| r = S_theory/S_exp, six strong lines | Paper I Table 1 | within 0.05 |

The Omega_4 band is wide on purpose: Omega_4 is carried by the two
U^(4)-only weak lines, and the reference implementation additionally
subtracts per-transition magnetic-dipole oscillator strengths computed from
free-ion eigenvectors (an atomic-structure kernel outside this toolkit's
scope; the dataset schema accepts an optional `f_md` column for callers who
have those numbers).
