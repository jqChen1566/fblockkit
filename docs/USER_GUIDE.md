# fBlockKit User Guide

> For chemists who do not program. Section numbers in this guide are the menu
> numbers of the program (the number is the path); menu numbers are append-only
> and this guide follows them.
>
> Version: v0.1.0 (early development). The provenance of every methodological
> criterion is in the program's output ("provenance" section and the References
> block) and in the data files under `src/fblockkit/knowledge/`.
>
> The full manual is built from `docs/manual/` (`bash docs/manual/build.sh`):
> the same structure, one chapter per menu entry, with worked examples that are
> replayable scripts; the built PDF ships with the release archives.

## 0 Before you start

fBlockKit does two things, and it never runs calculations for you:

1. **Generates input files** - it turns "which method, basis set, active space and
   convergence settings" into an input file you can submit to ORCA, plus plain
   guidance on how to run it;
2. **Characterises results** - it reads an output file you already have (or a
   structure file) and reports composition, entropy spectra, geometry and
   diagnostic conclusions. Every conclusion carries its provenance.

It does **not** run engines, does not use the network, and does not upload your
files anywhere.

Three ways to start:

- Green package: unpack and run `fblockkit.exe` (Windows) or `fblockkit` (Linux);
- Installed package: `pip install fblockkit`, then the `fblockkit` command;
- From source: `python -m fblockkit`.

You get a numbered menu. Type your answer at each prompt: **pressing Enter alone
takes the default**; a text file holding your inputs, one per line, is a replayable
script (Section 8).

## 1 Check-up and characterisation (read an ORCA output)

**What it is for**: turn one ORCA output file into an "analysis + diagnosis" report.

**What you need**: an ORCA 6.x output file (`.out`/`.log`). The per-MO composition
analysis needs the `LOEWDIN ORBITAL-COMPOSITIONS` table in the output (add
`%output Print[P_ReducedOrbPopMO_L] 1` to the input; it is printed at the normal
print level); when the table is absent, that analysis section is skipped with an
explanation.

**How**: choose menu 1, then give the file path.

**Optional prints**: two analyses read keywords of the run itself - the A1
composition table needs `%output Print[P_ReducedOrbPopMO_L] 1 end`, and the A8
diffuse-orbital (Rydberg) check needs `!PrintBasis` in addition (it reads the
shell exponents the run prints). Without them those sections are skipped, never
guessed.

**What you get** (written next to the input file; nothing you already have is
overwritten):

- `name.fbk.md` - the report: summary -> analysis sections (A1 orbital
  composition, A2 occupations and the entropy bound, A3 the multi-reference
  character panel; A6 the local spin analysis when the input divided the
  molecule into fragments; A8 the diffuse-orbital check when the basis is
  printed) -> diagnostic findings (graded refuse / error /
  warning / info) -> provenance summary -> References;
- `name.fbk.json` - the same content, machine-readable;
- The References block lists, for every cited work, the **complete citation**
  (authors, year, title, journal, volume(issue), pages, DOI) and a
  **paste-ready BibTeX entry**.

**How to read the report**:

- Every finding carries a rule id, a suggested action, the refusals that apply
  ("do not use this"), and its evidence (a manual quote, a literature DOI, or a
  measured record of this group);
- An "info" finding is a mandatory check (for example the CASPT2 reference
  weights), not a detected problem;
- Entropy thresholds and their precondition (a localized-orbital basis) travel
  with the report - never read the numbers without the caveats. The A2 spectrum
  is the maximum-entropy *bound* on the literature single-orbital entropy (an
  ORCA output prints only spin-summed occupations): a bound at or below 0.14
  excludes strong multireference character on that metric, a bound above it
  decides nothing and the section says so;
- The A5 projection-basis declaration check is attached to every crystal-field
  fit (Section 10); the A6 section states explicitly that the local spin
  analysis is not the environment spin-polarisation entropy Delta S_E.

## 2 Coordination geometry and symmetry hints (read an XYZ structure)

**What it is for**: see what the system looks like before computing anything.

**What you need**: a standard XYZ structure file.

**How**: menu 2 -> file path -> point group (optional; Enter to skip).

**What you get**: the coordination shell (metal centre, coordination number,
per-ligand distances) and symmetry hints (direction distribution, inversion
pairing, planarity), plus the "symmetry -> number of crystal-field parameters"
table (C3 -> 9, Oh -> 4, no real symmetry -> 27).

**Note**: these are geometric-distribution hints, not a rigorous point-group
determination; confirm the point group yourself and type it at the prompt - the
program then reports the parameter count for it.

## 3 Generate an ORCA input (system profile)

**What it is for**: answer a few questions and get a directly runnable ORCA input
plus plain run guidance.

**How**: menu 3, then: element list -> charge -> spin multiplicity -> f-block
metal valence -> targets (energy/geometry/excited/magnetic/spectra) -> structure
file (XYZ) -> basis tier number -> active space `nel,norb,mult,nroots`
(Enter = no `%casscf` block) -> convergence tier -> MaxCore per process in MB
(Enter = 2000; a measured f-block TRAH-CASSCF needed 9345 MB, so raise it for
large-basis hard cases).

Starting active spaces are suggested first (G2): the f-only minimal space with
its documented precedent, and an f+d double-shell option marked *provisional*,
each with its rationale.

**What you get**: `structure.fbk.inp` (pure ASCII, directly runnable) and, on
screen, the run guidance (method chain, basis, convergence discipline, refusals,
and what to check afterwards). When the profile names an f-block element or a
spin-orbit target (`magnetic`, `spectra`), the guidance adds the **relativistic
tier** the need implies: the one-electron variants for scalar work, the
two-electron correction for light main-group spin-orbit work, the
atomic-mean-field tier for the f block (where the two-electron correction is
worth only 1.6-4.9 cm^-1 on Nd3+ splittings -- so a higher tier must be
justified by the target accuracy), and the amfX2C variant is ruled out
everywhere for its unstable scalar approximation.

**Built-in discipline** (all from the ORCA manual or measured by this group; the
sources are in the generated guidance):

- Geometry optimisation / frequencies may only use a CASSCF reference layer -
  correlated layers are never put into an optimisation loop;
- Default convergence settings first; `TRAH` only for difficult cases, and it
  requires a matching /C auxiliary basis (otherwise the request is explicitly
  downgraded, never silently kept);
- Basis tiers are matched by element, valence and target (the SARC2 tier is the
  documented default for lanthanide wavefunction work); unsuitable tiers (for
  example 5f-in-core pseudopotentials for f-f transitions) are listed as "not
  applicable" with the reason.

## 4 Basis-set / ECP recommendation

**What it is for**: look up recommendations without generating a file.

**How**: menu 4, then element list, charge, multiplicity, valence, targets.

**What you get**: recommendation tiers (each with its matching auxiliary basis
and rationale), boundary cautions (for example: the standard def2 family stops
at Rn and does not cover the actinides), and not-applicable tiers with their hard
refusal conditions. Everything carries its provenance.

## 5 Search the tool index

**What it is for**: look up external tools. Each record has two dimensions that
must be read separately:

- **Tier (relation)** - how we intend to use it: A absorb (file/algorithm level),
  B interface (generate its input, read its output), C index (calculation kernels,
  including algorithmic borrowing);
- **Status** - whether it is actually used *now*: `active` / `registered (not
  enabled)`. "Registered" is a dot on the map, not a feature of this program.

**How**: menu 5 -> keywords (several words, space separated; all must match).

**What you get**: tier, status, licence and purpose of each hit; licences we have
not verified are marked "to verify".

## 6 Tool onboarding notes (by index id)

**What it is for**: the tier, status, where to obtain it, licence and integration
notes of a single tool.

**How**: menu 6 -> tool id (for example `openmolcas`; search first with Section 5).
Some entries carry a **workflow recipe** as well; the DMET code (`liblan`) is the
first one, and menu 6 prints the recipe after the guide -- the source's settings
(the loose double convergence criterion, the Loewdin cluster/environment split at
the f centre, subspace R-DIIS with the residual Delta S_E, the cluster SA-CASSCF,
the SOMF spin-orbit step, the ANO-RCC basis tiers) plus the gates that decide
whether a low-level wave function may be used. Two of those gates are this
toolkit's own sections: Delta S_E (Section 12) and the local-spin table
(Section 1).

## 7 Cross-level solution consistency

**What it is for**: when the same molecule has several candidate solutions each
computed at a cheap level (say HF) and an expensive level (say CCSD(T)), check
whether "pick the best solution by the cheap level" is safe.

**What you need**: a records JSON such as

```json
[
  {"label": "No.1", "occupation": "...", "energies": {"HF": -1930.87692, "CCSD(T)": -1932.97951}},
  {"label": "No.2", "occupation": "...", "energies": {"HF": -1930.87392, "CCSD(T)": -1932.97745}}
]
```

(The top level may also be an object with a `records` list; extra fields such as
`occupation` are kept for readability and otherwise ignored.)

**How**: menu 7 -> records JSON path.

**What you get**: warnings when the best solution differs between levels, and
when solutions in the cheap-level top N drop out of the expensive-level top N -
each with its evidence. The bundled `fixtures/literature/pucl3_s18.json`
(literature Table S18: 21 PuCl3 5f occupations x HF/CCSD(T)) is an instance of
this format; run menu 7 on it to see the effect:

- the HF-best solution (No.3) is not the CCSD(T)-best solution (No.1);
- the HF second-best (No.10) drops out of the CCSD(T) top three (about
  1.76 kcal/mol behind).

**How to use the conclusion**: re-rank the candidates at the expensive level
before choosing; do not pick a solution from the cheap-level ordering alone.

## 8 Save this session as a script

**What it is for**: keep every input you typed as a script, for one-command
replay ("record once, replay forever").

**How**: menu 8 -> script path (Enter = `fbk_session.txt`).

**Replay**: `fblockkit run script.txt`. Script format: one line per input, an
empty line means "just press Enter", `#` starts a comment. Replaying the same
script twice gives byte-identical output - this is also the project's regression
acceptance criterion.

## 9 SCF rescue (triage + corrected input files)

**What it is for**: read one output file whose SCF struggled, name what happened,
and - if you give the matching input file - write corrected inputs for you.

**How**: menu 9 -> ORCA output path -> input file path (Enter = skip the
corrected inputs).

**What you get**: the `SCF CONVERGENCE` block printed verbatim first (so the
informational rows - density, DIIS error - can be read exactly as ORCA printed
them; the check mode decides which rows are actually enforced, and nothing is
interpreted there), then the findings from the SCF triage (each with the measured
numbers - cycle counts, DIIS errors, achieved-vs-tolerance comparisons - and a
suggested action), then one corrected input per applicable fix, written next to
your input as `<name>.fix_<variant>.inp`. Your original files are never modified.

**What it will and will not propose**: damping (`SlowConv`) with the manual's
warning that it can converge closer to the initial guess; a two-step pre-SCF
route that reads the orbitals into your original input. It will *not* propose
merely raising `MaxIter` - the manual states that will not help in many cases -
that refusal is reported instead. The criteria used and every self-set threshold
are listed in the module and its tests.

## 10 Crystal-field fit (levels + coefficients JSON)

**What it is for**: fit the crystal-field parameters `B_k^q` from sampled state
energies and their projection coefficients (the mean-field-cost route of the
reference paper: cheap states plus exact linear algebra).

**What you need**: a JSON file such as

```json
{
  "point_group": "C3",
  "J": 7.5,
  "levels": [0.0, 89.85, 124.25],
  "coefficients": [[[0.0, 0.0], 1.0], ...]
}
```

`coefficients[i][m]` is the amplitude of state `i` on `|J M>` with `M = m - J`,
given as a number or a `[re, im]` pair. Pass the amplitudes themselves, not their
complex conjugates (the fitter warns about this trap).

**How**: menu 10 -> JSON path.

**What you get**: the fitted `B_k^q` (respecting the point group's allowed set:
C3-like 9, Oh 4, none 27), the constant, residual per state, the design-matrix
condition number, and rank warnings; plus a report file with the complete
citation and a paste-ready BibTeX entry. Sequence-to-level agreement is the
sanity check: diagonalising the fitted Hamiltonian must reproduce your input
levels (the bundled literature fixtures do this to 0.005 cm^-1).

**Projection-basis declaration check (A5)**: the same JSON may carry a
`declaration` object and a `compare` block, and the menu then appends the A5
check to the report:

```json
{
  "point_group": "C3", "J": 7.5, "levels": [...], "coefficients": [...],
  "declaration": {"convention": "Stevens", "projection": "J = 15/2",
                  "units": "cm^-1", "z_axis": "KD ground-state easy axis",
                  "origin": "metal site"},
  "compare": {"label": "Table S1 (L = 5)",
              "parameters": {"2,0": 1879.0, "2,2": -43.0, "...": 0.0},
              "declaration": {"projection": "L = 5"}}
}
```

The five declaration items are the ones a `B_k^q` set must record before its
numbers can be compared with another scheme (operator convention, projection
manifold, units, z-axis definition, origin/chirality). The comparison verdict
is per `(k, q)`: only the leading axial `(2, 0)` may be compared roughly
across schemes (the source's one measured example moved it by 3%), everything
else is refused without a matching declaration - the same wavefunction's
`B_2^2` reads +198 or -43 cm^-1 depending on the projection scheme (Chilton's
Table S1).

## 11 Point-charge crystal-field estimate (read an XYZ structure + charges)

**What it is for**: a first, model-level estimate of the crystal-field potential
created by the surrounding point charges - the classic point-charge model,
computed from the geometry alone (no quantum-chemical calculation). Use it to
see the sign pattern, the size and the symmetry of the field the ligands
produce, and as a sanity check on a fitted parameter set.

**What you need**: an XYZ structure, a point charge per element (in units of the
elementary charge, e.g. `O=-2,H=0.4`), and - for values in cm^-1 - the radial
expectation values `r2,r4,r6` in Angstrom^k. The radial moments are
ion- and method-specific and are deliberately **not** shipped with this tool:
without them the menu still reports the geometry-only lattice sums and says
explicitly that no cm^-1 parameters were produced.

**How**: menu 11 -> XYZ path -> charges -> radial moments (Enter = geometry
only).

**What you get**: the charges used, the lattice sums per `(k, q)`, the
parameters in cm^-1 (when the radial moments are given), the seven f-orbital
energies (the potential matrix's eigenvalues) and the overall splitting, plus
the symmetry note (which `(k, q)` survive) and a report file with the complete
citations.

**What it is not**: the parameters are the exact one-electron point-charge
matrix elements in the f-orbital (`l = 3`) basis - they are *not* the
`J`-manifold Stevens parameters of Section 10 (the operator-equivalent
projection is not carried out here), and the point-charge model itself is
qualitative: measured deviations make it unusable for quantitative parameter
extraction (the same caveat the section prints, with its citation). The
`l`-basis numbers and Section 10's `B_k^q` are comparable only through their
spectra, not parameter by parameter.

## 12 Exact active-space analysis (CASSCF output + FCIDUMP; optional orca_2json exports)

**What it is for**: Exact four-state entropy of the active space, rebuilt from
your own calculation - the autoCAS single-orbital entropy itself, the
four-state quantity `s_i = -sum w ln w` over the (empty / up / down / double)
occupations of each active orbital. Section 1
reports a rigorous upper bound for the same number (an ORCA output alone cannot
supply it, and ORCA's own 2-RDM export was measured to be unreachable for
CAS-type methods on 6.1.1); this section reports the number, with its
cross-checks.

**What you need**: (a) the converged CASSCF output; (b) the FCIDUMP that ORCA
writes when the same job is rerun with the converged orbitals plus the
`!FCIDUMP` keyword (`!moread` and `%moinp "previous.gbw"`). The dump run then
reports "IS NOT FULLY CONVERGED" and stops - that is how the keyword behaves;
the cross-checks below verify the dumped data against the original run. For the
localized-basis step (required before reading the 0.14 line): run `orca_loc` on
the active window, then `orca_2json` on both gbw files (canonical and
localized) with a configuration requesting the MO coefficients and the
`S-Matrix`. On lanthanides give the IAO basis explicitly - the default minimal
basis is not defined for them; ANO-RCC-MB (option 4) works, and the exact
working `orca_loc` line is recorded in the fixtures README (Eu3+ chain).

**How**: menu 12 -> output path -> FCIDUMP path -> canonical export (Enter =
skip the localized step) -> localized export -> active window (Enter = infer
from the occupation table). The inference grows the window away from the core
when an active orbital prints as 0.0000 and validates it against the `N(occ)=`
line; give the window explicitly whenever you prefer.

**What you get**: the four occupation weights and entropy per orbital (canonical
basis, informational), the same spectrum in the localized basis when the two
exports are given, the engine cross-checks (reconstructed energy vs the printed
CASSCF energy - the gate is 1e-6 Eh; the agreement is ~1e-10 against the
summary line and reaches 1e-13 with a fuller-precision reference - and natural
occupations vs `N(occ)`), the 0.14-line candidates and the plateau readout, and
a report file next to the FCIDUMP with the complete citations.

**What it is not**: the active space is taken as given (the usual window and
solution-branch caveats apply; Section 1 checks them); the solver is a dense
determinant CI capped at a few thousand determinants per M_s sector - far above
the f-block windows it targets, but not a DMRG; and the numbers do not replace
Section 1's diagnostics.

**Two further analyses of the same state** run in the same pass when the
exports are given:

- *environment spin-polarisation entropy* (Delta S_E of Ai et al., Eq. 9):
  with the localized export, every active orbital is assigned to the centre of
  its largest Löwdin population, the orbitals not on the cluster centre form
  the environment, and Delta S_E = -2 Tr[(D/2) ln(D/2)] + Tr[D_a ln D_a] +
  Tr[D_b ln D_b] measures the spin polarisation that sits there instead of on
  the metal. The criterion's zero is exact; the source's own anchors are 0.007
  for its correct solution and 2.766 for its wrong one (a factor 42 in the
  fitted crystal-field MAE). Two boundaries are reported with the number: the
  inactive orbitals of a CASSCF wave function are doubly occupied by
  construction, so this is the *active* environment (a mean-field solution's
  ligand spin polarisation is diagnosed by Section 1's local-spin table
  instead), and a closed-shell state has Delta S_E = 0 in every basis. The
  cluster centre is asked for explicitly (Enter = the f-block centre);
- *atomic-term check*: the effective L, S and J of the computed state from
  <L^2>, <S^2>, <L.S> - the operators are built by projecting the active
  orbitals onto the f AOs of the f-block centre (the export's AO labels carry
  the centre and the component), and the expectations are exact for the CI
  state. The verdict compares (L_eff, S_eff) with Hund's rules for the shell
  occupation, the check that decides whether a mean-field density may be used
  for a crystal-field fit; the J value is classified as the maximal-J
  component (what a scalar stretched density carries), the minimal-J one, or a
  non-stretched component, and only a mismatch in L or S is a defect. The
  f character of each active orbital is reported as the quality measure, and
  the Hund verdict is withheld below the provisional 0.9 line - which is also
  how the f-block wrong-solution mode shows up (the electrons sit in d/s/p
  orbitals and the f shell is empty).

Both readings need the canonical export (and the localized one for Delta S_E)
and the active window; they are appended to the same report file.

## 13 Orbital-space comparison (two orca_2json exports)

**What it is for**: two zero-external-reference checks of an orbital space, both
read from the singular values of `M = C_A^T S C_B` between two orbital sets
(each set is orthonormalized internally first; the source states the result is
independent of the orthonormalization scheme):

- the *subspace fraction* sigma_F = ||M||_F / sqrt(min(|A|, |B|)) answers "how
  much of space B is contained in space A" - use it to rank a recommended
  active space against a fuller reference space (the deficit 1 - sigma_F
  orders like the energy error in the source; a saturated sigma_F does not
  exclude that A is too large);
- the *space-change SVD* compares two sets of the same size (an initial and a
  final active space of one optimization): singular values near 1 mean the
  space did not move, and a value near 0 means one initial active orbital was
  replaced by an unrelated one, i.e. the initial space lacked an element.

**What you need**: two `orca_2json` exports of the same system in the same
basis (for the identity case: the same gbw before and after `orca_loc`, which
exercises the check's blind spot - a rotation inside one space must come out as
"unchanged").

**How**: menu 13 -> first export path -> its window `first last` in ORCA's
0-based orbital numbering (Enter = all orbitals) -> second export path -> its
window. A window that is not given covers every orbital of the export.

**What you get**: the singular values, sigma_F and the deficit, the smallest
singular value, and the two readings with their provisional bands (anchored on
the DMET source's own numbers: 1 - sigma_F = 2.9e-3 for a deliberately loose
start and 6.8e-5 after a localization; the AVAS source's SVD range 0.65-0.99).
The report file is written next to the first export, citations included.

**What it is not**: an energy - the checks rank spaces and flag replaced
orbitals; they do not evaluate the spaces. The bands are provisional and the
intended use is comparison across candidate spaces for one system.

## 14 AVAS target projection (an orbital export + a target AO shell)

**What it is for**: the AVAS construction of an active space from target atomic
orbitals (Sayfutyarova et al. 2017). The target shell is projected onto the
occupied and the virtual orbital blocks of an SCF/CASSCF export; the overlap
eigenvalues say which orbitals carry target character, and the report gives the
resulting `(n_electrons, n_orbitals)` at the chosen truncation threshold.

**What you need**: one `orca_2json` export (the MO coefficients, the
occupations and the S matrix), plus the target: a centre, an angular momentum
and the shell numbers as the export's AO labels spell them.

**How**: menu 14 -> export path -> target centre (Enter = the f-block element)
-> angular-momentum letter (Enter = `f`) -> shell number(s) (Enter = every shell
of that momentum) -> truncation threshold (Enter = 0.1, the source's range is
0.05-0.1) -> open-shell option (Enter = 3).

**What you get**: the occupied-side and virtual-side overlap spectra, the
orbitals kept at the threshold with their indices, the `(n_el, n_orb)`
recommendation, the ready `%scf avas` block for ORCA, and a report file next to
the export with the citations.

**How to read it**:

- values near 1 on the occupied side are essentially the target AO itself;
  low values among the kept ones name orbitals strongly mixed with their
  environment (in the N2 example the three kept occupied overlaps are 0.70,
  0.56, 0.56 - the sigma and the pi pair - and the recommendation comes out as
  the textbook (6 electrons, 6 orbitals));
- **the virtual side decides whether AVAS suits the system**: if no virtual
  orbital passes the threshold, the antibonding target character is
  ligand-centred (the source's [CuCl4]2- case) and the remedy is to add ligand
  AOs to the target, or to read it as "this target is not suited";
- the open-shell **option 3** (the default here) keeps every singly occupied
  orbital in the active space regardless of its target overlap. The source's
  option 2 can force an unoccupied beta orbital into the core, after which
  CASCI may lie above the variational HF energy - the reason the f block is
  advised to use option 3;
- both AVAS quality criteria are falsify-only (the source states it): a clean
  reading does not establish that the space is best. Section 13's space-change
  SVD and the "CASCI below the variational HF energy" check can still falsify
  it.

**Boundaries**: the target set here is a subset of the calculation's own AOs -
the source's non-minimal-ANO variant, which is what makes the projection
computable from an export alone and is the only route open to the f block
(ORCA's own AVAS minimal basis was measured to have no f-block entries, Eu).
ORCA's CASSCF takes its active space by orbital order, not by an index list, so
the deliverable is the size and the verdict on a window, not a restart file;
the emitted `%scf avas` block covers the systems where ORCA's built-in AVAS
works.

## 15 Orbital portrait (an export -> descriptor table)

**What it is for**: the deterministic part of RLEASE's orbital descriptor -- a
per-orbital panel of the things a chemist reads off one by one: occupation,
orbital energy, the dominant centre and its Löwdin share, the angular-momentum
composition, a **bonding label**, and an AO-centre estimate of the orbital's
extent and charge-centroid displacement. Nothing here needs a model; the
source's neural predictor and threshold policy are explicitly not adopted.

**What you need**: one `orca_2json` export.

**How**: menu 15 -> export path -> window `first last` in ORCA's 0-based
numbering (Enter = the orbitals with fractional occupations, i.e. the active
space of a correlated export).

**What you get**: the table on screen and a report file
(`<export>.portrait.fbk.md`) with the citations.

**How to read it**:

- the **bonding label** is the source's own criterion: accumulate the atom-pair
  population `S_AB = sum_{mu in A} sum_{nu in B} c_mu S_mu,nu c_nu` over pairs
  within 6 Angstrom; a single atom carrying more than 95% of the Löwdin
  population makes the orbital non-bonding, otherwise the sign of the cross
  terms decides. Two measured additions are documented with the criterion: a
  cross term within 0.01 counts as non-bonding (a symmetric 1s core measures
  3e-4, and a bare sign test would call it "bonding"), and a row marked `!` is
  cancellation-heavy (the absolute sum of its atom-pair blocks exceeds 5 -- the
  N2 orbitals run up to 2.2 while one node-heavy virtual measures 158), so its
  cross magnitude is not a bond order even though its sign still reads;
- the **shares** are Löwdin atomic populations: non-negative and summing to one
  for every orbital, unlike the raw block sums.
- the extent and centroid displacement are **AO-centre estimates** (the exact
  `<r^2>` and the diagonal integrals need integrals an export does not carry;
  an FCIDUMP has them when needed).

**What it is not**: a selection criterion by itself -- the source's own
benchmark found that a larger active space is not automatically a better one
(autoCAS picked 17 orbitals for CH4 where a smaller space gave half the error),
so read this panel as evidence for a choice, not as the choice.

## 16 Magnetic-doublet criterion (a Kramers-doublet table JSON)

**What it is for**: the empirical metric a Dy(III) single-molecule-magnet study
reads off the Kramers-doublet ladder -- which doublet can still act as a step of
the relaxation barrier and which one opens quantum tunnelling (QTM):

```text
g_T = (g1 + g2 + g3 * sin(theta3)) / 3        line:  g_T * theta3 = 20
```

with `g1 <= g2 <= g3` the doublet's principal g values and `theta3` (in
**degrees**) the angle between its `g3` axis and the ground doublet's `g3` axis.
Below the line the doublet supports excitation above it; at or above it, QTM is
opened. Nothing here is computed by an engine: this is pure post-processing of
the g tensors, so it runs from a small table you prepare.

**What you need**: a JSON table of the doublets. Per doublet, the three
principal g values and **either** `theta3` (as published tables give it)
**or** `axis3` (the unit axis of the largest value, when you have the tensor
directions), plus an optional `energy`:

```json
{
  "system": "label, optional",
  "reference": 0,
  "doublets": [
    {"label": "ground",       "g": [0.41, 0.44, 9.04], "energy": 0.0,
     "axis3": [0.0, 0.0, 1.0]},
    {"label": "1st excited",  "g": [0.41, 0.44, 9.04], "energy": 101.0,
     "theta3": 9.52}
  ]
}
```

`reference` (an index or a label, default 0) names the doublet whose `g3` axis
is the quantisation axis the computed angles are measured against; its own
`theta3` is 0 by definition. If every doublet carries `theta3` given, no
reference axis is needed.

**How**: menu 16 -> table JSON path.

**What you get**: the doublet table on screen with `theta3`, `g_T`,
`g_T * theta3` and the verdict per doublet, plus a report file
(`<table>.magnetic.fbk.md`) with the criteria, the applicability domain and the
citations.

**How to read it**:

- the angle is in **degrees** -- that is a measured fact of the calibration
  table (recomputing its values in radians disagrees in the second digit, and
  the regression check pins all 38 tabulated doublets down);
- an angle computed from `axis3` is folded to [0, 90] degrees, because a
  principal axis's sign is arbitrary (the source's own values run to 89.9);
- the reference doublet is marked `*`; for doublets that give `theta3`
  directly, the number is used as given (the axis is not needed);
- **the criterion is empirical, not derived**: the source states that the
  `g1`/`g2` weights should in principle be reduced as `g3` rotates into the
  plane, but those historic calculations do not carry that information, and the
  line's physical background (the material's internal dipolar field) varies
  with packing and dilution -- so the line is not a universal constant;
- **applicability domain**: the line was calibrated on 19 mononuclear Dy(III)
  complexes (2016-2025, one consistent methodology); other ions or nuclearities
  need their own calibration;
- **a known failure mode the criterion cannot see**: in the source's own set
  exactly one doublet is misread (`[Dy(Cp^ttt)2]+`'s QTM doublet falls deep in
  the safe zone, product 2.52). A "safe" verdict is the criterion's reading,
  not a proof of barrier behaviour.

**Boundaries**: this menu reads a table, it does not compute g tensors. The
tensors come from your own magnetic-property calculation (for the OpenMolcas
chain that menu 3's guidance generates, the doublet g values are read from the
SINGLE_ANISO output; the toolkit's `g_T` check is the post-processing step of
that chain). The line is printed as source evidence, and the report always
carries the domain and the outlier note.

## 17 Cross-structure orbital mapping (a manifest of exports)

**What it is for**: keeping one active space *consistent* along a path -- a
reaction coordinate, a bond-dissociation curve, a scan. Active spaces chosen
per structure drift apart, the total correlation energy is then computed
inconsistently, and the relative energies become erratic. The protocol maps
the localized orbitals of the structures onto each other (no interpolation of
structures needed) and transfers a selection through the map: any orbital that
changes along the path (a bond being broken or formed) ends up in a
non-matchable block, and the whole block joins the active space together
whenever one of its orbitals is selected anywhere.

**What you need**: one `orca_2json` export per structure -- *localized*, since
localized orbitals are what transfers between structures. The exports must
carry the `T-Matrix` (the kinetic criterion) and the `S-Matrix`/labels (the
shell-population criterion): export with a `<basename>.json.conf` file of
`{"MOCoefficients": true, "1elIntegrals": ["H", "S", "T", "V"]}`. Localize
with `orca_loc` first (the frozen example localizes the occupied block with
`LocMet 4`, the IBO scheme; ORCA's IAO-based methods refuse virtual orbitals,
so the virtual block uses Foster-Boys, `LocMet 2`). All structures must use
the same basis set, atom order and orientation.

**How**: menu 17 -> a manifest JSON:

```json
{
  "structures": [
    {"name": "r=1.094", "export": "n2_scan_1.094.loc.json"},
    {"name": "r=2.600", "export": "n2_scan_2.600.loc.json"}
  ],
  "tau": 0.5,
  "selections": {"r=1.094": [4]}
}
```

Export paths are relative to the manifest. `tau` (optional) is the threshold
of both criteria; without it the source's rule is applied per orbital set: the
threshold minimising the non-matchable count, with the whole curve printed.
`selections` (optional) names the orbitals to make consistent, per structure
(0-based indices; e.g. the output of an entropy-based selection for one
structure).

**What you get**: per structure the descriptor table (occupation, kinetic
energy, largest shell-wise populations), the tau curve, the classes the map
built, the non-matchable sets, and -- with `selections` -- each structure's
consistent active space; plus a report file (`<manifest>.mapping.fbk.md`) with
the citations.

**How to read it**:

- the criteria are the source's Eq. (1): same kinetic energy
  (`|t_i - t_j| < tau`, in Eh) and the same shell-wise populations (summed
  absolute difference `< tau`), with one `tau` for both;
- a **class** is a set of orbitals that map onto each other across *all*
  structures, in both directions and as a whole set -- degenerate orbitals
  (e.g. a symmetric trio) ride together. The listed non-matchable block is the
  changing set: in the N2 scan example the two 1s cores and the three bond
  orbitals map, while the two one-centre hybrids (whose s/p ratio follows the
  bond length) are the non-matchable set at `tau = 0.5`;
- the **tau curve** is the feature, not one number: its plateaus are where the
  non-matchable count is stable (the source reads its own example at
  `tau <= 0.5`, and so does the example here);
- **populations here are shell-wise Löwdin populations** computed from the
  export; the source uses IAO ones (its own text licenses any orbital-wise
  population analysis). The population *scale* therefore differs from the
  source's, so `tau` is a data-calibrated parameter here, never copied from
  the source's absolute values -- read the curve, not a recalled threshold.

**Boundaries**: the mapping compares structures in their own frames -- same
basis, same atom order, same orientation; two structures that differ in the AO
shell set or the orbital count are refused. The virtual manifold of a tiny
system can be genuinely unmatchable (diffuse orbitals reorganize with
geometry): the report shows it rather than hiding it. The selection itself is
not made here: menu 12's entropy protocol (or any other) provides it, and this
menu reports the consistent space it implies.

## Appendix A Command line

```text
fblockkit                     interactive menu
fblockkit --record FILE       interactive, recording your input as a script
fblockkit run SCRIPT.txt      replay a script
fblockkit search KEYWORDS     search the tool index
fblockkit guide TOOL_ID       tool onboarding notes
```

## Appendix B Installation and distribution

- **Green package**: unpack and run; no Python needed (built for Windows and
  Linux from `packaging/fblockkit.spec`);
- **pip package**: `pip install fblockkit` (Python >= 3.11). Data files (rules,
  templates, menu, indexes, BibTeX database) ship with the package - see the
  package-data declaration in `pyproject.toml`.

## Appendix C Citations and references

Every criterion in the output carries provenance in three classes: **manual**
(section number + URL), **literature** (complete citation with DOI, plus a
paste-ready BibTeX entry in the References block), and **measured** (a record of
this group). The reference database is `src/fblockkit/knowledge/sources.bib`;
literature evidence is rejected at load time unless it points to an entry there.
When you cite an external tool, use its own official citation format (given in
its onboarding notes).
