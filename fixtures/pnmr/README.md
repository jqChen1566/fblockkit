# pNMR fixtures (fixtures/pnmr/)

The regression data and one end-to-end run of the pseudocontact-shift menu
(35; `src/fblockkit/analysis/pnmr.py`).

## mares2018_co2.yaml

Transcribed values of Mares & Vaara, Phys. Chem. Chem. Phys. 2018, 20,
22547-22555 (CC-BY, DOI 10.1039/c8cp04123g) for their Co(II) clathrochelate
benchmark: the **two printed susceptibility tensors** (the non-symmetric
eq. (4) and the symmetric variant of eq. (5), in 1e-32 m^3 in the file's
comments) and the magnetic parameter sets behind them (NEVPT2/TZVP-DKH:
D = -85.5 cm^-1, E/D = 0.00281, g_iso = 2.325, Delta g_ax = 0.979; the
experimental set: D = -64 cm^-1, Delta g_ax = 0.2).  Published anchors the
tests regress against: Delta chi_ax = 14.8e-32 m^3 (non-symmetric) and
27.3e-32 m^3 (symmetric) for the parameter set, 7.1e-32 m^3 for the
experimental set.

## co_plus/ (the end-to-end run)

CO+ (the QDPT probe geometry, an f-block-like doublet is not needed for the
format; a real f-ion example is registered as the next fixture) with its own
QDPT output: `co_plus_qdpt_g.yaml` reads the effective-Hamiltonian g from
`co_plus_qdpt.out` (copy kept here so the run directory is self-contained),
builds the non-symmetric susceptibility at 300 K and prints the PCS of the
oxygen.  The tiny PCS (0.085 ppm) and axiality (0.0004e-32 m^3) are the
expected result: CO+ is almost isotropic (Delta g ~ 5e-4), which is exactly
why this run doubles as the invariant: our chi from the parsed g reproduces
ORCA's own printed susceptibility to the printed precision (the 4 pi
cgs-emu conversion, tested).
