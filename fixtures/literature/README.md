# Literature data fixtures (fixtures/literature/)

This directory holds regression data transcribed from literature tables. Unlike the raw
output in `fixtures/orca/`, this data has been transcribed once, so its source and the
way it was cross-checked must be recorded.

| File | Content | Source | Cross-check |
|---|---|---|---|
| `pucl3_s18.json` | Total energies (Hartree) of the 21 5f-occupation solutions of PuCl3 at the HF and CCSD(T) levels | Lu J.-B. et al., *J. Chem. Theory Comput.*, 2025, **21**, 170–182, DOI 10.1021/acs.jctc.4c01189 (SI Table S18) | Compared item by item with the numbers quoted in the literature review §4.6 (the three points No.1, No.3 and No.10, plus the full HF/CCSD(T) span and the 1.76 kcal/mol difference) |
| `peng_s1_cf.json` | The 9 $B_k^q$ of Er-trensal (C3) (columns HF@HF / PBE0@HF / CASPT2) and 8 Kramers doublet levels (two columns) | Peng L. et al., *J. Phys. Chem. Lett.*, 2025, **16**, 12312–12320, DOI 10.1021/acs.jpclett.5c02971 (SI Table S1, S6) | Checked value by value against the converted SI text (`jz5c02971_si_001.md`); used as the known answer of the CF fitter (`analysis/crystal_field.py`) |
| `peng_s2_s7_oh_cf.json` | The 4 $B_k^q$ of Cs2NaDyCl6 (Oh) (two columns) and 8 Kramers doublet levels (two columns) | As above (SI Table S2 first block, S7) | As above; two further independent pieces of evidence confirm that this block belongs to compound 2 -- cubic symmetry allows only 4 parameters, and the caption of Figure S7 says "4 $B_k^q$ parameters" |

Use:
- `pucl3_s18.json` -- the regression fixture of `diagnosis.cross_level_check`, checking the
  risk of "picking a solution by the cheap level's ordering alone";
- `peng_s1_cf.json`, `peng_s2_s7_oh_cf.json` -- the known answers of
  `analysis.crystal_field` (A4): the basic self-check is to build and diagonalise
  $\hat H_{\mathrm{CF}}$ from the tabulated $B_k^q$, which should reproduce the levels in
  the same table.

Transcription discipline: the table data comes from the MinerU conversion of the SI text;
any number cited in a test must agree with a number appearing in a close-reading file
(`文献细读/细读_5f入芯ECP_Lu2025.md`, `文献细读/细读_CF哈密顿量_Chan2025.md`). A disagreement
counts as a transcription doubt and must be settled against the original PDF before the
number is changed.

Known internal inconsistency in the source data (measured 2026-09-25, not a transcription
error): the two HF@HF columns of `peng_s1_cf.json` do not agree with each other -- the
spectrum obtained by diagonalising the nine $B_k^q$ of S1 is always 0.76336 times the
levels listed in S6 (the ratio is the same to 4e-5 across the seven excited doublets).
The two PBE0@HF columns of the same paper and both columns of `peng_s2_s7_oh_cf.json`
agree exactly (deviations below the tabulated precision of 0.005 cm$^{-1}$), so this is a
one-off in the source table; the details and how it is handled are in that fixture's
`cross_check` block.
