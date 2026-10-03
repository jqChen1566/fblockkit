# Literature data fixtures (fixtures/literature/)

This directory holds regression data transcribed from literature tables. Unlike the raw
output in `fixtures/orca/`, this data has been transcribed once, so its source and the
way it was cross-checked must be recorded.

| File | Content | Source | Cross-check |
|---|---|---|---|
| `pucl3_s18.json` | Total energies (Hartree) of the 21 5f-occupation solutions of PuCl3 at the HF and CCSD(T) levels | Lu J.-B. et al., *J. Chem. Theory Comput.*, 2025, **21**, 170–182, DOI 10.1021/acs.jctc.4c01189 (SI Table S18) | Compared item by item with the numbers quoted in the literature review §4.6 (the three points No.1, No.3 and No.10, plus the full HF/CCSD(T) span and the 1.76 kcal/mol difference) |
| `peng_s1_cf.json` | The 9 $B_k^q$ of Er-trensal (C3) (columns HF@HF / PBE0@HF / CASPT2) and 8 Kramers doublet levels (two columns) | Peng L. et al., *J. Phys. Chem. Lett.*, 2025, **16**, 12312–12320, DOI 10.1021/acs.jpclett.5c02971 (SI Table S1, S6) | Checked value by value against the converted SI text (`jz5c02971_si_001.md`); used as the known answer of the CF fitter (`analysis/crystal_field.py`) |
| `peng_s2_s7_oh_cf.json` | The 4 $B_k^q$ of Cs2NaDyCl6 (Oh) (two columns) and 8 Kramers doublet levels (two columns) | As above (SI Table S2 first block, S7) | As above; two further independent pieces of evidence confirm that this block belongs to compound 2 -- cubic symmetry allows only 4 parameters, and the caption of Figure S7 says "4 $B_k^q$ parameters" |
| `chilton_s2_dy19.json` | The 19 mononuclear Dy(III) SMM rows (a "safe" and a "definite QTM" doublet each: $g_1,g_2,g_3,\theta_3$ plus the tabulated $g_T$, and the experimental $U_{eff}$) | Chilton, *Chem. Soc. Rev.*, 2025, DOI 10.1039/d5cs00493d (SI Table S2 and the seven methodology items above it) | All 38 tabulated $g_T$ values recomputed from the raw columns (largest disagreement $8\times10^{-4}$); the claimed separation reproduced (largest safe product 15.31, smallest definite-QTM 30.25) with exactly one doublet on the wrong side of the line -- the source's recorded outlier. Regression fixture of `analysis.magnetic_doublets` (menu 16) |
| `lnsim_precision.json` | The DMET-vs-all-electron errors of Ai et al. (Tables 1--3): pre-SOC state MAE by low-level solution, and the CASSCF-SO / NEVPT2 final accuracy (with the expanded-cluster variant) | Ai Y. et al., *J. Chem. Theory Comput.*, 2025, **21**, 9631--9640, DOI 10.1021/acs.jctc.5c01336 | Transcribed from the close reading (`文献细读/细读_LnSIM_DMET_CASSCF-SO_JCTC2025.md` section 5.3); the recipe's accuracy expectations (`recipe/dmet.py`) are checked against it in `tests/test_dmet.py` |

Use:
- `pucl3_s18.json` -- the regression fixture of `diagnosis.cross_level_check`, checking the
  risk of "picking a solution by the cheap level's ordering alone";
- `peng_s1_cf.json`, `peng_s2_s7_oh_cf.json` -- the known answers of
  `analysis.crystal_field` (A4): the basic self-check is to build and diagonalise
  $\hat H_{\mathrm{CF}}$ from the tabulated $B_k^q$, which should reproduce the levels in
  the same table;
- `chilton_s2_dy19.json` -- the regression fixture of `analysis.magnetic_doublets`
  (menu 16): the criterion's $g_T$ arithmetic and its separation are checked against the
  source's own table; a slice of it is also the menu-16 manual example;
- `lnsim_precision.json` -- the numbers the DMET recipe's "accuracy expectations" section
  prints, so that the product and the transcription cannot drift apart.

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

## Spectroscopic anchor tables (digitized 2026-10-03)

Thirteen tables from the classic f-block spectroscopic literature were digitized cell-by-cell
into CSV (the spectroscopy-extension round; the full process record, in Chinese, is in
`调研/侦察2026-10-03_光谱扩展/E0_锚获取现状.md`).

| File | Source | Table | Rows | Units |
|---|---|---|---|---|
| `carnall1968_eu3_aquo_levels_U.csv` | Carnall, Fields, Rajnak, *J. Chem. Phys.* **49**, 4450–4455 (1968), DOI 10.1063/1.1669896 | TABLE II | 56 | cm$^{-1}$ |
| `carnall1968_eu3_aquo_parameters.csv` | as above | TABLE III | 16 | cm$^{-1}$ |
| `carnall1988_laf3_eu3_cf_levels.csv` | Carnall, Goodman, Rajnak, Rana, ANL-88-8 (1988), DOI 10.2172/6995771 | Appendix V | 78 | cm$^{-1}$ |
| `carnall1988_laf3_ln_energy_parameters.csv` | as above | Table 4 | 297 | cm$^{-1}$ |
| `gorller1991_eu3_md_strengths.csv` | Görller-Walrand, Fluyt, Ceulemans, Carnall, *J. Chem. Phys.* **95**, 3099–3106 (1991), DOI 10.1063/1.460867 | TABLE IV | 73 | $10^{-7}$ D$^2$ |
| `gorller1991_eu3_md_expt_vs_theor.csv` | as above | TABLE V | 14 | $10^{-7}$ D$^2$ |
| `carnall1978_laf3_omega_lambda.csv` | Carnall, Crosswhite & Crosswhite, ANL-78-XX-95 (1978), DOI 10.2172/6417825 | App. XIV Table 1 | 11 | $10^{-20}$ cm$^2$ |
| `carnall1978_laf3_eu3_U2.csv` | as above | App. VII Table 2 | 187 | dimensionless |
| `carnall1989_an3_lacl3_parameters.csv` | Carnall, ANL-89/39 (1989), DOI 10.2172/7017566 | Table 3 | 226 | cm$^{-1}$ |
| `carnall1968_ln3_aquo_levels_U.csv` | Carnall, Fields, Rajnak, *J. Chem. Phys.* **49**, 4424–4442 (1968) part I, DOI 10.1063/1.1669893 | TABLE IV–XII | 381 | cm$^{-1}$ |
| `carnall1988_laf3_ho_tb_gd_cf_levels.csv` | Carnall, Goodman, Rajnak, Rana, ANL-88-8 (1988), DOI 10.2172/6995771 | App. VIII/X/XI | 752 | cm$^{-1}$ |
| `carnall1988_laf3_remaining_cf_levels.csv` | as above | App. I/II/III/IV/VI/VII/IX (Pr, Nd, Pm, Sm, Er, Tm, Dy) | 1448 | cm$^{-1}$ |
| `carnall1978_laf3_lifetimes.csv` | Carnall, Crosswhite & Crosswhite, ANL-78-XX-95 (1978), DOI 10.2172/6417825 | App. XIV Table 2 | 17 | microseconds |

Coverage:

- **Carnall 1968 (Eu$^{3+}$ aquo)** -- $^7$F and excited J levels (observed / calculated /
  deviations) and the $U^{(2)}, U^{(4)}, U^{(6)}$ reduced matrix elements of the
  $^7$F$_0$ and $^7$F$_1$ transition groups (the standard input of Judd–Ofelt intensity
  analysis), plus the free-ion/CI least-squares parameters.
- **Carnall 1988 (Ln$^{3+}$:LaF$_3$)** -- all 78 crystal-field components of the 12 J
  states of Eu$^{3+}$:LaF$_3$, and the complete energy-parameter table of 13 Ln$^{3+}$
  ions ($F^2, F^4, F^6, \zeta, \alpha, \beta, \gamma, T^2$–$T^8$, $M^0$, $P^2$, nine
  $B_q^k$, $n$, $\sigma$). The separate file `carnall1988_laf3_ho_tb_gd_cf_levels.csv`
  adds the observed/calculated levels and O–C of Ho$^{3+}$ (App. VIII), Tb$^{3+}$
  (App. X) and Gd$^{3+}$ (App. XI), reconstructed from the two independent side-by-side
  streams with empty-label-continues-previous grouping and per-state expansion of the
  printed leading-digit omissions. Some states legitimately appear as two separate
  printed runs (e.g. Ho $^5$G$_4$ and $^5$F$_2$; each run individually satisfies
  $2J+1$). Five Tb blocks from the source's high-energy shorthand section carry a
  `grouping caveat` note: the values are as printed, but the state assignment within
  those narrow packed blocks is uncertain.
- **Görller-Walrand 1991 (Eu$^{3+}$ magnetic-dipole standards)** -- the MD strengths of
  $^5$D$_1\leftarrow{}^7$F$_0$, $^5$D$_0\leftarrow{}^7$F$_1$ and
  $^5$D$_2\leftarrow{}^7$F$_1$ in five symmetries (IC and J-mixing schemes), plus the
  experiment-vs-theory comparison; the source's means 18, 94, 9 $\times10^{-7}$ D$^2$
  (transitions in order) are the widely quoted MD standards.
- **Carnall 1968 part I (Ln$^{3+}$ aquo ions)** -- level assignments (observed / calculated /
  $\Delta E$) and the $U^{(2)}/U^{(4)}/U^{(6)}$ squared reduced matrix elements for the Pr,
  Nd, Pm, Sm, Dy, Ho, Er and Tm aquo ions (TABLE IV–XII; 381 entries). Where the source
  table prints a second environment it is carried in the `env2_*` columns (LaCl$_3$ for
  Nd/Sm/Dy/Ho; Tm(C$_2$H$_5$SO$_4$)$_3\cdot$9H$_2$O for Tm). Twenty-one J-label misreads of
  the conversion (a systematic 3/5/9$\to$8 pattern) were adjudicated on page renders and
  corrected; the list is in the CSV header.
- **Carnall 1988, remaining seven appendices (Pr, Nd, Pm, Sm, Er, Tm, Dy)** -- the same
  "Experimental and Computed Energy Level Structure" tables as the two files above,
  digitized with a header-driven generic parser (each SLJ group in a table header becomes
  an independent stream; empty label continues the previous state). Pm has computed
  levels only; Nd/Sm/Dy print reference columns from other papers (Rast, Dieke, Fry,
  Caspers, Wong, Voron'ko), of which only the report's own columns are taken.
  Verification scope note: this file's blocks received structural validation plus targeted
  page spot-checks, not the full page-by-page adjudication applied to the other files.
- **Carnall 1978 (Ln$^{3+}$:LaF$_3$ Judd–Ofelt parameters)** -- the fitted $\Omega_2$,
  $\Omega_4$, $\Omega_6$ of 11 trivalent lanthanides in LaF$_3$ (ANL-78-XX-95,
  Appendix XIV Table 1), with the per-ion literature source flags (a=Krupke 1966;
  b=approximate values from that report; c–h = Weber/Krupke/Pappalardo). Checked
  double-source: text-layer extraction asserted equal to 170-dpi page-render reading
  (33/33 values).
- **Carnall 1978 (Eu$^{3+}$:LaF$_3$ unit-tensor matrix elements)** -- the
  $[U^{(2)}]^2, [U^{(4)}]^2, [U^{(6)}]^2$ between pairs of J levels of Eu$^{3+}$:LaF$_3$
  (ANL-78-XX-95, Appendix VII Table 2; 187 printed entries). By the report's stated
  convention, absent rows are zero matrix elements. Cross-checked against 25 rows of a
  160-dpi page render (all equal, including the printed 0.0 vs 0.0000 distinction) and
  against three rows of the Carnall 1968 aquo U-table (differences only in the last
  printed digit, e.g. $0\to2$: 0.1374 here vs 0.1375 aquo).
- **Carnall 1989 (An$^{3+}$:LaCl$_3$ energy-level parameters)** -- the full
  effective-operator parameter set ($F^2/F^4/F^6$, $\zeta$, $\alpha$, $\beta$, $\gamma$,
  $T^2$–$T^8$, $M^0/M^2/M^4$, $P^2$, four D$_{3h}$ crystal-field parameters
  $B_0^2/B_0^4/B_0^6/B_6^6$, $\sigma$, $n$) for 11 trivalent actinides U–No
  (ANL-89/39, Table 3; 226 entries). Eight conversion/OCR artefacts were adjudicated
  against 150–1600 dpi page renders and fixed in a recorded whitelist: e.g. Pu$^{3+}$
  $M^4$ $[.388]\to[.338]$ and Fm$^{3+}$ $M^4$ $.312\to.612$, both proven by the invariant
  $M^4/M^0=0.385$ ratio across nine ions; the four B rows (superscripts lost in
  conversion) restored as $B_0^2/B_0^4/B_0^6/B_6^6$ per D$_{3h}$ symmetry; ion label
  Ee$\to$Es; the "Na" in the footnote read as No.

Method: MinerU 4.0.10 (hybrid-OCR) PDF-to-Markdown, then the parse scripts in
`scratch_anchor/parse_*.py` (rowspan expansion, superscript recovery, footnote
stripping, missing-value normalisation); every key table cross-checked cell-by-cell
against 180–1200 dpi page renders; the Eu$^{3+}$:LaF$_3$ table checked against $2J+1$
component counts for all 12 J states; the Görller-Walrand TABLE IV blocks recomputed as
$2\Sigma D(\pi)+\Sigma D(\sigma)$ against the printed SUM rows (28/30 agree to ±0.01).

Data discipline: numbers are the printed values -- no interpolation, no correction, no
unit conversion; ellipsis/dash cells are empty with the original mark kept in a note
column. Source-level notes: two SUM rows of Görller-Walrand TABLE IV are internally
inconsistent (0.36 and 0.02; stored as printed after 700–900 dpi magnification -- source
inconsistencies, not transcription errors); the Pm$^{3+}$ column of the ANL-88-8 Table 4
is interpolated in the source (footnote a) and flagged. One cross-source warning: the
U-table values quoted by a later paper differ from the Carnall 1968 original in the last
digit (e.g. 0.1449 vs 0.1450); Judd–Ofelt regressions must use one source consistently
and state which.

Use: these CSVs are the regression anchors of the planned f–f absorption chain (menu 34
Judd–Ofelt route) and of the core-excited-spectroscopy work; they are also the known
answer against which the planned crystal-field-level JO mapping is checked.

Not yet digitized (source data on hand; the scripts are reusable as-is): the remaining
per-ion $[U^{(\lambda)}]^2$ tables and radiative lifetimes of ANL-78-XX-95; and the
energy-level appendices / (MD)$^2$ table of ANL-89/39. Both ANL reports (195 pp. and
289 pp.) were secured as OSTI open full texts, with MinerU conversions archived in
`scratch_anchor/anl_md/`.
