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

## 1 Check-up and characterisation (read an ORCA or Gaussian output)

**What it is for**: turn one ORCA output file (deep analysis) -- or a
Gaussian 09/16 output (the minimal fact set: termination, SCF, frequency and
optimization facts) -- into an "analysis + diagnosis" report.

**What you need**: an ORCA 6.x output file (`.out`/`.log`), or a Gaussian
09/16 output; the diagnosis rules run off the shared fact vocabulary, so
Gaussian files get the engine-independent checks (termination, SCF,
frequencies, optimization) while ORCA-only analyses are skipped. The per-MO
composition analysis needs the `LOEWDIN ORBITAL-COMPOSITIONS` table in the
ORCA output (add `%output Print[P_ReducedOrbPopMO_L] 1` to the input; it is
printed at the normal print level); when the table is absent, that analysis
section is skipped with an explanation.

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

**What it is for**: look up recommendations without generating a file --
plus an optional datasource leg that queries the deployed basis library.

**How**: menu 4, then element list, charge, multiplicity, valence, targets;
finally (optionally) an element symbol for the deployed basis library
(Enter skips).  The library is auto-detected: `$FBK_BASISDB`, then
`~/projects/orca_basis_sets` (the basisdb deployment, a Basis Set Exchange
snapshot with a SQLite index).

**What you get**: recommendation tiers (each with its matching auxiliary basis
and rationale), boundary cautions (for example: the standard def2 family stops
at Rn and does not cover the actinides), and not-applicable tiers with their hard
refusal conditions. Everything carries its provenance.  The datasource leg
lists, for the queried element, the sets that cover it -- name, family,
relativistic method, whether ORCA reaches them through a built-in keyword or
through GTOName, ECP pairing, contracted-function size and reference count.
No library content ships with this program: the query runs against the
deployment you point at.

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
suggested action): no convergence, a crash without a verdict, pseudo-convergence,
a long run, a DIIS error that rebounds after the AO-DIIS switch, an energy
trajectory that oscillates without collapsing, and unmet enforced criteria. Then
one corrected input per applicable fix, written next to your input as
`<name>.fix_<variant>.inp`. Your original files are never modified.

**What it will and will not propose**: damping (`SlowConv`) with the manual's
warning that it can converge closer to the initial guess; a two-step pre-SCF
route that reads the orbitals into your original input; and, when the DIIS
trajectory rebounds or oscillates, a `TRAH` variant (the manual's robust
second-order SCF - ORCA has no keyword spelled ARH, TRAH is its
trust-region augmented-Hessian route with AutoTRAH on by default; an input
that already carries TRAH gets no such proposal, and a CASSCF input carries
the /C auxiliary-basis caveat in the proposal). It will *not* propose merely
raising `MaxIter` - the manual states that will not help in many cases - that
refusal is reported instead. The criteria used and every self-set threshold
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

**Optional: the active-space overlap (`active` block)**: give per-structure
active-orbital lists (`"active": {"r=1.094": [4, 5, 6], ...}`) and the report
adds the *overlap-preservation* scalar of the second source (its Supporting
Information): between each adjacent pair, |det S_act| ~ 1 means the active
space survived, ~ 0 means at least one active orbital exchanged with the
inactive space (measured on the frozen scan: the cores give 1.0000 for every
pair, the bond triad 0.9955 across a 0.01-Angstrom step and 0.72 across
0.5-Angstrom steps). The number uses the *second* structure's overlap -- the
source's own small-step approximation -- so between distant structures it
demonstrates degradation rather than certifying preservation, and the printed
bands are this project's reading of the source's qualitative scale. Note the
two checks answer different questions: the determinant tracks the active
*subspace* (a rotation inside a degenerate window leaves it near 1), while the
mapping follows individual orbitals.

## 18 WASP guess transfer (neighbours -> a gbw-ready mkl)

**What it is for**: starting every structure of a series from an orbital set
that already lies in the right basin. The WASP protocol (the review's
Eqs. (11)-(13)) builds the guess for a target geometry as the 1/d-weighted
interpolation of the *neighbouring* structures' orbital coefficients, d the
RMSD between structures, and orthonormalises the mixture in the target's
overlap metric. This menu writes that guess into a **gbw-ready Molekel mkl**,
which `orca_2mkl` converts into a `.gbw` that ORCA reads through `!moread` +
`%moinp` -- the toolkit's write-back route to ORCA's initial guess (measured
on 6.1.1: the converted file is accepted as `INITIAL GUESS: MOREAD`; the
round trip preserves coefficients to the print precision).

**What you need**: the neighbour exports (localized orbitals of the
already-computed structures, same basis and atom order as the target), the
target geometry's own export (**it supplies the overlap metric**) and its mkl
(`orca_2mkl <target> -mkl` after any cheap run at that geometry -- the
template supplies geometry, basis and the file layout).

**How**: menu 18 -> a manifest JSON:

```json
{
  "structures": [
    {"name": "r=1.094", "export": "n2_scan_1.094.loc.json"},
    {"name": "r=2.600", "export": "n2_scan_2.600.loc.json"}
  ],
  "template": {"export": "n2_scan_1.600.loc.json", "mkl": "n2_scan_1.600.mkl"},
  "delta": 1.0
}
```

`delta` (optional, Angstrom) restricts the neighbourhood; without it every
given neighbour counts. Paths are relative to the manifest.

**What you get**: the weights table (distance and 1/d weight per neighbour),
the orthonormalisation residual, the file
`<template stem>.fbk.mkl` (run `orca_2mkl <template stem>.fbk -gbw` next to
it), and a report (`<manifest>.guess.fbk.md`) with the citations. The
occupations and orbital energies written along come from the nearest
neighbour; the orbitals are the mixture.

**How to read it**:

- the residual should sit at round-off (3e-15 in the example): the mixture is
  orthonormal in the target's metric before it is written;
- a neighbour *at* the target geometry takes the whole weight (the 1/d weight
  would diverge): the guess then is that structure's orbital set as it stands,
  and the report says so;
- **a guess selects the basin, not just the starting point.** In the frozen
  example the interpolated guess converged the target geometry's RHF to a
  *second* solution 19 mEh *below* the branch the scan itself used -- a real
  multiple-solution instance (RHF/def2-SVP at 1.6 Angstrom). After any guess
  transfer, verify that the solution still belongs to the series (energy,
  orbital character, the mapping): basis continuity is not guaranteed by a
  good guess.

**Boundaries**: same basis, atom order and orientation across the structures
(the coefficients live in each structure's own AO frame; nothing is aligned
silently -- the RMSD is taken as given). The written mkl is a *guess*: ORCA
re-optimises; nothing here is a converged result. The machine-readable form of
the mkl (and the reader/writer) is documented in the format appendix.

## 19 Dipole-moment candidate batch (a structure -> the DM-AS scan inputs)

**What it is for**: preparing the scan behind the dipole-moment active-space
selection (Kaufold et al. 2023; the CASCI sequel 2026). The protocol runs the
ground state of every candidate space, compares each space's dipole moment
with a reliable reference, and keeps the closest -- a selection whose numbers
are, in principle, experimentally checkable. Menu 19 writes the batch; menu 20
reads it and makes the choice.

**What you need**: the structure (XYZ) and a few settings. The protocol is
built for neutral singlet ground states with a nonzero dipole moment (both
restrictions are the source's, and the menu refuses the rest with the reason).

**How**: menu 19 -> structure file -> charge (Enter = 0) -> multiplicity
(Enter = 1) -> basis keyword (Enter = def2-TZVP) -> orbital preparation
(Enter = MP2 natural orbitals, the source's choice; `hf` for plain HF
orbitals) -> largest active space in orbitals (Enter = 14) -> include the
PASS+ extras (Enter = no).

**What you get**: a directory next to the structure
(`<structure>.fbk.dm/`) with the preparation input (`prep.inp`: RHF + MP2
with `NatOrbs true` -- the natural orbitals stay in the gbw), one CASCI input
per candidate (`cand_e{ne}o{no}.inp`: `!NoIter moread` + `%moinp` the prep
gbw + `%casscf nel/norb/nroots 1`), a run script (`run_scan.sh`) and a
manifest (`manifest.json`) whose candidate outputs are already named; fill in
the reference after the runs and hand it to menu 20.

**How to read it**:

- the candidate family is the source's PASS set: 6-14 active electrons (even),
  with `n_e/2 + 3 <= n_o <= 14`; the PASS+ extras (the 4-electron and
  two-virtual rows) can be included -- the source kept them in its scanning
  set but flagged them as generally poor (35 candidates by default, 51 with
  the extras);
- the candidates share one orbital source: the MP2 natural orbitals of the
  preparation run. The source's own benchmark found MP2 the best starting
  point for the CASCI dipole moments; `hf` is the documented fallback;
- each candidate is a CAS-CI (no orbital optimisation -- that is what makes
  the scan cheap) and prints its own dipole moment, which is what the
  selection reads.

**Boundaries**: the batch runs nothing (this toolkit never does); the script is
yours to run. The source excludes charged systems (the dipole of a charged
molecule depends on the coordinate origin), and the protocol needs a nonzero
dipole moment.

## 20 Dipole-moment selection (the batch outputs + a reference)

**What it is for**: the second half -- rank the candidate spaces by
`|mu(candidate) - mu(reference)|` and report the chosen space.

**What you need**: the filled manifest. The reference is either another run's
output (`{"output": "dft.out"}` -- its SCF block is used, a KS-DFT value is
what the source used) or a number you supply (`{"magnitude_debye": 1.85,
"source": "NIST (gas phase)"}`), optionally with its `vector` for the
directional protocol.

**How**: menu 20 -> manifest path. `protocol` (optional) is `gdm` (default;
the source's own recommendation names CASCI-GDM-AS the sensible default) or
`vgdm` (the directional variant, which needs the reference vector).

**What you get**: the ranking table on screen and a report
(`<manifest>.dm_select.fbk.md`) with the citations.

**How to read it**:

- the smallest deviation wins; ties go to the smaller space -- the protocol's
  purpose is a suitable space, not an unnecessarily large one;
- the ranking, not the winner alone, is the result: in the manual's water
  example the six candidates run from 0.0081 D ((6e, 8o), selected) to
  0.0245 D against the PBE0 reference (2.0801 D), and the spread is the
  protocol's resolution;
- the reference value drives the choice -- the source tested experimental
  (NIST) and several functionals' values, and to publish a selection you
  should state the reference and its origin;
- each candidate is read from its **single-root** dipole block (State: 0). A
  state-averaged output prints only the average (State: -1) and is refused
  with the rerun instruction; the source averaged six states, and the
  difference from the single-root densities is documented in the report;
- for a potential-energy scan, the protocol does not guarantee the same space
  at every geometry: select at one geometry and carry the space with menu 17
  when consistency matters.

**Not implemented, with the reason**: the per-state variants (EDM-AS,
D2DM-AS) compare per-state dipole moments with per-state TD-DFT values;
measured on ORCA 6.1.1, a state-averaged CASSCF/CASCI run prints only the
state-average dipole and the TDDFT module prints transition dipoles but no
state dipoles, so their data is not obtainable from ORCA outputs. The source's
own recommendation puts CASCI-GDM-AS first; the per-state route (per-state
densities crossed with the exported dipole integrals) is recorded as future
work.

## 21 APC orbital ranking (an export -> a ranked active space)

**What it is for**: the ranked-orbital scheme of King & Gagliardi (JCTC 2021):
score every candidate orbital (all doubly occupied ones plus a window of the
lowest virtuals) by the approximate pair coefficient, then keep the top of the
ranking up to a cap on the number of configuration state functions (CSFs).

**What you need**: one `orca_2json` export of a converged RHF run (closed
shell), with a `<base>.json.conf` requesting the Fock blocks:

```json
{ "MOCoefficients": true, "1elIntegrals": ["H"], "FockMatrix": ["J", "K"] }
```

Produce the export where the run wrote it: the Fock terms are read from the
run's `<base>.densities` / `<base>.densitiesinfo` pair (measured on ORCA
6.1.1), so a bare `.gbw` copy cannot be exported for them. The default `apc`
variant needs only the `K` block; the `apcx` variant (the source's
exact-integral comparison) additionally needs the windowed two-electron
block, requested as a *second window*:

```json
{ "MOCoefficients": true, "2elIntegrals": ["MO_IAJB"],
  "OrbWin": [0, 6, 7, 16, 0, 0, 0, 0] }
```

Write the window as eight integers -- the second window all zeros for a single
window, because the four-integer form is rejected by `orca_2json`; the numbers
are inclusive 0-based indices, first/last internal then first/last external.

**How**: menu 21 -> export path; the ranking variant (`apc` default, `apcx`
for the exact-integral comparison); the candidate window (default 23 lowest
virtuals in energy, the source's general-scheme choice); the CSF cap
(`max(7,6)` = 490, `max(8,8)` = 1764, `max(10,10)` = 19404, `max(12,12)` =
226512, or any integer); the model gap source (`energies` default, or `fock`
for the Fock diagonal, which is the correct choice for localized orbitals).
For the default there is no extra ORCA run and no integral transformation: the
ranking needs only the export.

**What you get**: the ranked table on screen and a report
(`<export>.apc.fbk.md`) with the citations.

**How to read it**:

- the cap is the only binding constraint: the selection drops the
  lowest-ranked orbital until the CSF count (equation 2 of the source) fits,
  never leaving fewer than one occupied and two unoccupied orbitals in the
  space (the source's CASSCF-stability floor; anything kept by the floor is
  listed);
- the ranking is the result, not the winner alone; the entropies depend on
  the candidate window, so the window is printed with the table;
- APC is a cheap screening scheme calibrated on small molecules: the source
  measures it to overestimate doubly-occupied orbital entropies (R^2 0.64,
  MAE 0.0240 against DMRG) and expects degradation in much larger systems or
  where the HF determinant is a poor reference. The report states this with
  every run; treat a selection as a screening answer, not a converged one;
- APC and APCX on the same window usually agree on the top set: in the
  manual's N2 example the pi/pi* quartet ranks in the top four of both while
  the magnitudes differ by an order of magnitude (the diagonal-sum
  approximation overestimates; the source finds APC still ranks better
  because the errors cancel);
- the source's `max(a,b)` labels are ambiguous in its own text ((7e, 6o)
  evaluates to 210 CSFs; its listed 490 is a seven-orbital space). This tool
  carries the labels and numbers verbatim and treats the CSF count as the
  constraint.

**Not implemented, with the reason**: open-shell references (ROHF/UHF). The
source assigns singly occupied orbitals the maximum approximated entropy, and
the UNO variant needs unrestricted exports; the per-spin Fock-block
conventions are not measured yet, so the menu refuses non-RHF exports with
the instruction to run RHF.

## 22 ASS1ST round-1 input (a structure + an initial active space)

**What it is for**: first round of the ASS1ST active-space construction
(Khedkar & Roemelt 2019/2020): a CASSCF calculation whose
perturbation-theory density is kept, so that menu 23 can read the
quasi-natural occupation numbers and propose the next space.

**What you need**: a structure (XYZ) and a *small but chemically reasonable*
initial active space -- the source's advice: singly occupied orbitals for open
shells, the metal d shell for first-row transition metals, a few pi/pi*
orbitals for conjugated systems. The outcome depends on this choice, but with
conservative thresholds different sensible initials converge together.

**How**: menu 22 -> structure path; charge and multiplicity (Enter = 0 and 1);
the initial space `nel,norb`; the number of states (Enter = 1; more than one
makes it a state-averaged round for menu 23); the method/basis keyword line
(Enter = `RHF def2-SVP TightSCF`); MaxCore in MB (Enter = 2000). The menu
writes `<stem>.r1.inp` and, next to it, `<stem>.r1.json.conf`.

The generated input carries the blocks the round needs (measured on ORCA
6.1.1): `FIC-NEVPT2` (the unrelaxed density exists only for the FIC ansatz),
`KeepDens` (the density is read from the run's sidecar), and the CASSCF
`PTSettings` block with `Density Unrelaxed` + `NatOrbs true`.

**Then**: run the input (`orca <stem>.r1.inp`), export it
(`orca_2json <stem>.r1.gbw`; the request file is already beside the `.gbw`),
and feed `<stem>.r1.json` to menu 23.

**Boundary to know**: the round's active orbitals are selected by ORCA's own
default window from its starting guess -- the source's near-ideal
quasi-natural restart is not written yet. Check the active orbitals in menu
23's output (it prints their CASSCF occupations) and, if the window missed the
intended set, adjust the input before continuing the chain.

## 23 ASS1ST selection round (an export -> the next round's input)

**What it is for**: one selection step of ASS1ST: read the round's NEVPT2
unrelaxed density, diagonalize its internal/internal and external/external
blocks separately (the source's construction), and propose the next active
space from the threshold band.

**What you need**: the `orca_2json` export of a finished round, with the
request `{"MOCoefficients": true, "1elIntegrals": ["S"], "Densities": ["all"]}`
(menu 22 writes it); the run must have produced its density sidecar
(`KeepDens` + the `PTSettings` block, both in the generated input).

**How**: menu 23 -> export path; the threshold band (`0.05` means the window
[0.05, 1.95] -- the source's conservative value, `0.03` is its relaxed one;
two independent lines `T_ext,T_int` are allowed, e.g. `0.03,1.96`); state
weights (Enter = equal; used when the round averaged several states); the
spaces you have visited so far (`ne,no`, space-separated -- this is what
enables the cycle warning); the method/basis keywords for the next round
(Enter = as in the example).

**What you get**: the block quasi-occupation tables with the band marked, the
active orbitals with their CASSCF occupations, the next-space suggestion, and
a report (`<export>.ass1st.fbk.md`) with the citations. Unless the round is
self-consistent, the menu also writes the next round's input
`<stem>.r{N+1}.inp` (and its `.json.conf`) with the suggested space.

**How to read it**:

- every quasi-occupation is a *block* eigenvalue of the exported NEVPT2
  unrelaxed density -- the source's construction. ORCA's own printed
  "Natural Orbital Occupation Numbers" block is the whole-space
  naturalization and is *not* what the scheme reads (measured); the external
  block in particular is far from diagonal, so the block step matters;
- band bookkeeping (the source's): an internal quasi-orbital inside the band
  adds +2 electrons and +1 orbital to the active space, an external one adds
  +0 electrons and +1 orbital; an active orbital whose CASSCF occupation has
  drifted to the top of the band is reassigned to internal (-2 e, -1 o), one
  at the bottom to external (-0 e, -1 o).  The source allows the space to
  shrink as well as grow;
- `self-consistent` means no orbital crossed the band in either block and no
  active orbital drifted: use that space for the production calculation. If
  the same space was visited before, the source warns of cycling between
  spaces -- break the tie by chemical judgment;
- **the density is ORCA's FIC-NEVPT2 unrelaxed density**: the source uses the
  SC-NEVPT2 first-order density, and ORCA answers the unrelaxed density only
  for the FIC ansatz (measured: SC + `Density Unrelaxed` is refused).  Same
  object family (the density of the zeroth-plus-first-order wavefunction),
  different contraction -- counts near the thresholds may shift, and no
  literature comparison exists for the exchange.  The report states this;
- the source's own caveats apply: results depend on the initial space (stated
  above) and on the potential-energy point (run the scheme at the geometries
  you care about, not just one).

**Measured chain behaviour (in the fixtures)**: on N2/def2-SVP, round 1
CAS(6,6) suggests shrinking the sigma pair out -- (6e, 6o) -> (4e, 4o), the
pi/pi* quartet -- and a rerun of that space from the default guess converged
to a *different-shaped* (4,4) solution, which the next analysis then reports
as a further shrink.  It is the multi-solution behaviour this project records
elsewhere, and the reason for the boundary note in menu 22: check the round's
active orbitals each time.

## 24 QICAS active-space optimization (an FCIDUMP -> F_QI-minimized orbitals)

**What it is for**: QICAS (Ding, Knecht & Schilling 2023) optimizes the orbital
basis of a CAS(n,m) scheme *without touching the Hamiltonian*: the cost function
is the out-of-CAS correlation `F_QI` -- the sum of the four-state entropies of
the non-active orbitals -- which depends only on the 1- and 2-RDM, and the
source's Theorem 1 bounds the CASCI error by it
(`E_CASCI - E_FCI <= dE_max/ln4 * F_QI`).  Minimizing F_QI over orbital
rotations therefore realigns the active space so that the CASCI approaches the
FCI.

**What you need**: the FCIDUMP of a converged CASSCF run (the same file menu 12
reads), the target active space `nel,norb` within the window, and -- optionally
-- the run's output so the CI root is matched to the printed CASSCF energy (the
menu-12 cross-check).

**How**: menu 24 -> FCIDUMP path; output path (Enter = skip the cross-check);
the target `nel,norb`; the rotation set (`touch` = every pair touching a
non-active orbital, the source's chemical-accuracy choice; `exclusive` =
active/non-active pairs only, its economical variant for large spaces).

**What you get**: the per-orbital entropy profile before/after the optimization,
the F_QI values, the CASCI energies in the initial and the optimized basis
against the window FCI, the Theorem-1 check with the run's own numbers, and a
report (`<FCIDUMP>.qicas.fbk.md`) with the citations.

**How to read it**:

- the partition is read two ways: the *requested* one (this tool orders the
  dumped orbitals by natural occupation and cuts the highest as closed, the
  lowest as virtual) and the *final* one -- the occupancy reading of the
  optimized basis (an orbital counts as closed if its occupancy ends above 1).
  The optimizer may rotate a strongly correlated orbital out of a closed slot
  and an uncorrelated one in, so the final space can read as a different
  (n,m) -- that is the mechanism, and the manual example shows it: an awkward
  (2,4) request on N2 is repaired to the (4,4) pi space, and the CASCI gap
  drops from 6.7e-2 Eh to 5.3e-3 Eh;
- the Theorem-1 line is a genuine check, not a decoration: it is evaluated with
  `dE_max` from the CI spectrum and F_QI; if a target sits so far from the
  window ground state that the reference overlap deficit exceeds 1/2 (the
  source's proof assumption), the inequality is reported as not informative
  rather than silently skipped;
- **scope**: this is QICAS applied *on the FCIDUMP window* (the source's own
  subset application) with exact RDMs from the four-state-entropy route.  The
  source drives it from a full-space DMRG ground state, which this toolbox does
  not have -- so F_QI here is not comparable with the source's full-space
  numbers;
- the optimized rotation is reported as a matrix and used for the CASCI check;
  writing it back into a `.gbw` (the mkl route of menu 18) is not implemented,
  and neither is the source's size-selection variant (its Appendix C: minimize
  the total orbital entropy, read the plateau of the threshold diagram).

## 25 AEGISS selection (entropy screening + an AO projection)

**What it is for**: the AEGISS workflow (Tarocco et al., arXiv 2026) joins
the two complementary selection ideas: *how correlated* an orbital is (the
single-orbital entropy of the AutoCAS family) and *whether it belongs to the
chemistry of interest* (the atomic-orbital projection of the AVAS family).
The result is a compact, chemically anchored active space.

**What you need**: the `orca_2json` export of the CASSCF gbw (needs S and the
orbital labels), the same run's FCIDUMP (for the exact entropies -- the
four-state-entropy route), optionally the run's output for the energy
cross-check, and the AO label of the chemistry of interest (`C pz`, `Fe d`,
`N p`).

**How**: menu 25 -> export path; FCIDUMP path; output path (Enter = skip);
the AO label; the entropy fraction (`tau_E = tau * S_max`; Enter = 0.1, the
source's AutoCAS-style default -- its benzene example used 0.2); the
projection threshold (Enter = 0.5, the source's value).

**What you get**: the five-step trace -- every window orbital's occupation,
exact entropy, screen verdict, projection weight and final verdict -- the
resulting `(n_electrons, n_orbitals)`, and a report
(`<export>.aegiss.fbk.md`) with the citations.

**How to read it**:

- the entropy screen keeps every orbital strictly above `tau * S_max` (the
  relative line, so it transfers between systems); the projection keeps the
  entropy-screened orbitals whose weight exceeds the threshold;
- **the projection weight is the projection norm** `sum |O[eta, p]|^2` (the
  standard AVAS reading).  The source writes it as a signed row sum; measured
  on the benzene fixture that sum cancels by symmetry for every nodal pi
  orbital (only the nodeless a2u survives), so this tool uses the norm and
  says so in every report.  On the same fixture the norm separates cleanly:
  sigma orbitals measure exactly 0.0000, the pi manifold 0.05-3.09;
- the AO label's shell index is optional: `C pz` matches every pz shell of
  every carbon (12 functions here), `C 2pz` only the second shell (6).
  Measured on benzene: the whole-family label keeps all six pi orbitals at
  the 0.5 threshold -- the textbook (6e, 6o) -- while the shell-resolved
  label narrows the space, which is the resolution knob the label offers;
- **scope, stated in the report**: the entropies are exact over the FCIDUMP
  window (the source estimates them from a DMRG over a larger window), and
  the target is the calculation's own AO subset (the non-minimal AVAS route
  of menu 14 -- the only one open to the f block).  Multi-group unions are run
  one group per pass; the selection is delivered as the window orbitals to
  include, and feeding it to a CASSCF needs the orbital-order machinery (see
  menu 22's boundary note).

**Measured anchor**: on the benzene/cc-pVDZ pi platform the exact window FCI
energy reproduces the engine's printed CASSCF energy to eleven digits
(-230.793818898, and the source's own benzene CASSCF -230.793770 agrees to
5e-5), and the whole-family label recovers the (6e, 6o) pi space.

## 26 TNASS subset selection (an FCIDUMP -> the S2-maximizing subset)

**What it is for**: TNASS (Mingare, Heuzé & Coveney, arXiv 2026) picks the
active space as the spatial-orbital subset whose bipartition with its
complement has the **largest Rényi-2 entropy** `S2(A) = -log Tr(rho_A^2)` --
a multi-orbital entanglement measure beyond the single-orbital entropy of the
AutoCAS family.

**What you need**: the FCIDUMP of a converged CASSCF run (the same dump the
entropy and QICAS menus read), optionally the run's output for the energy
cross-check, the target size (number of active spatial orbitals) and the
method.

**How**: menu 26 -> FCIDUMP path; output path (Enter = skip); target size;
method (`greedy` default, `block K` for the block-greedy variant, `brute` for
the exact search, capped at C(N,n) <= 20000).

**What you get**: the one-orbital S2 seed ranking, the step-by-step trace
(each addition and the subset's S2), the resulting space with its electron
count, and a report (`<FCIDUMP>.tnass.fbk.md`) with the citations.

**How to read it**:

- **the S2 oracle is exact over the FCIDUMP window** (the four-state-entropy
  route's CI, grouped by the subset's occupation pattern; the one-orbital S2
  reproduces `-log(sum w^2)` from the four-state weights to 1e-12).  The
  source builds its entanglement feature as a tensor network from a
  low-bond-dimension DMRG state over the full space; there is no
  bond-dimension truncation here, but the window boundary is the restriction;
- the bipartition symmetry `S2(A) = S2(complement)` holds for every subset --
  useful as a self-consistency check of the trace;
- greedy cannot backtrack: on both example windows the trace shows an S2
  *drop* on the final addition, and the brute-force optimum beats it
  (benzene 0.23280 vs 0.23000; N2 0.13134 vs 0.12663).  Use `brute` when the
  window is small enough, and treat a greedy result as the source's practical
  approximation;
- the subset is a spatial-orbital set (both spins together), matching the
  source's selection domain; the source's best-k variant (choose k by the
  CASCI energy) is not implemented -- run the menu at several target sizes
  and compare;
- delivering the space to a CASSCF needs the orbital-order machinery (see
  menu 22's boundary note).

## 27 DeltaSCF / MOM excited-state SCF input (a structure + an excitation spec)

**What it is for**: ORCA's DeltaSCF route converges the SCF to a chosen
excited-state solution (a higher stationary point of the SCF energy surface)
by constraining the frontier occupations and following them with a
maximum-overlap metric.  The menu generates the input; the engine run stays
with you.

**What you need**: a structure (XYZ) and the excitation spec -- an `ALPHACONF`
occupation list (`0,1` = HOMO->LUMO, `0,0,1` = HOMO->LUMO+1, `0,1,1` =
HOMO-1->LUMO; `0,2` = double excitation for RHF), optionally a `BETACONF`
list, or `ionize N` for a core ionization (`IONIZEALPHA N`; charge and
multiplicity stay those of the reference system).  A converged ground-state
`.gbw` to start from is recommended (the source's own advice).

**How**: menu 27 -> structure path; charge and multiplicity of the reference
system; the occupation spec; optional BETACONF; the method/basis keywords
(Enter = `PBE0 def2-TZVP UHF`; every source example uses UHF); the MOM metric
(`mom` regular, `pmom` the Hratchian projection-operator variant, `imom` =
`KeepInitialRef`); hard-case tactics (`none`, `freeze` = `FreezeAndRelease`,
`gmf`); the ground-state gbw path (Enter = skip the MORead lines).  The input
is written as `<stem>.dscf.inp`.

**How to read it** (the report prints these with every input):

- run the input and **check the converged state** (occupations and character):
  a DeltaSCF run that silently fell back to the ground state carries no
  warning in the output;
- single-determinant states only: open-shell singly excited states inherently
  break spin symmetry and need spin purification; multi-determinant cases are
  out of scope (the manual's boundary);
- not apt for most pi->pi* states -- the manual names benzene's HOMO->LUMO
  explicitly -- but reasonable for particle-hole states, spatially separated
  occupied/virtual pairs, and closed-shell doubly excited states;
- for a core ionization, localize the core orbital first if it is not atomic;
- do not feed the wavefunction to single-reference correlation methods as if
  it were a ground state.

**Measured anchor**: the generated input (formaldehyde, `ALPHACONF 0,1`,
UHF) converged to the n->pi* saddle at -114.294975 Eh against the clean
ground state's -114.418617 Eh -- an excitation of 3.364 eV, the textbook
vertical n->pi* value.

## 28 RAS / ORMAS model-space input (a structure + a partition mask)

**What it is for**: the generalized active space (GAS) concept in ORCA's
spelling -- an incomplete model space defined by a partition mask in the
`refs` sub-block.  Two masks: `RAS(Nel: NRAS1 MaxHoles / NRAS2 / NRAS3
MaxParticles)` and `ORMAS(nel: m1 min1 max1, m2 min2 max2, ...)` (up to 25
sub-spaces, commas or slashes).  Two routes: `%casscf` (orbital-optimized
RASSCF / ORMAS-SCF) and `%rasci` (a standalone CI after the frozen core;
the manual's MRCISD-style usage).  The menu generates the input; the engine
run stays with you.

**What you need**: a structure (XYZ) and the mask.  The mask is the single
source of the active-space numbers -- the RAS orbital counts must sum to
`norb` exactly (the engine refuses otherwise), and both masks override the
`nel`/`norb` lines -- so the generated input always describes the mask's
own space, with the numbers written out consistently.

**How**: menu 28 -> structure path; charge and multiplicity; the route
(`casscf` or `rasci`); the mask type (`ras` or `ormas`); the mask (e.g.
`6:2 2/2/2 2`, or `6: 2 0 4, 2 0 4, 2 0 4`); multiplicities and nroots
(Enter = the structure's multiplicity, 1 root); on the CI-only route a
CIStep (`accci`/`csfci`/`detci`/`treecsf`) and an ExcLevel (Enter = the
module's defaults); the method/basis keywords (Enter = `RHF def2-SVP`).
The input is written as `<stem>.rasormas.inp` (MCSCF route) or
`<stem>.rasci.inp` (CI-only route).

**How to read it** (the report prints these with every input):

- check the partition in the output before trusting the number: the `%casscf`
  route prints its configuration counts ("Building the RAS space ... done (N
  configurations)", the ORMAS sub-space table), the `%rasci` route echoes
  "Number of active orbitals" and "First active orbital";
- the orbital optimization omits the active-active rotation for incomplete
  model spaces (the manual's own warning), so the energy is sensitive to the
  orbital canonicalization (`actorbs` / `actconstraints`); compare spaces,
  not only energies;
- for unconstrained sub-spaces ORMAS reproduces the full CAS to all printed
  digits (measured 12-digit identity on N2) -- a deviation from the CAS
  energy is exactly the restriction you asked for;
- ORCA silently overwrites `nel`/`norb` from either mask (measured); the
  generator refuses inconsistent numbers instead, so the input you get
  always carries the mask's own parameters;
- the `%rasci` module runs no orbital optimization: its energy sits above
  the corresponding RASSCF unless you supply converged orbitals.

**Measured anchors** (N2/def2-SVP at 1.10 Angstrom): the plain CASSCF(6,6)
reference gives -108.989034756374 Eh; `ORMAS(6: 2 0 4, 2 0 4, 2 0 4)`
reproduces it to all 12 printed digits; `RAS(6:2 2/2/2 2)` raises the MCSCF
energy to -108.985623236851 Eh; the CI-only route gives -108.921051085808 Eh
(RAS) and -108.921178474942 Eh (the superset ORMAS space).

**Boundary**: arbitrary-CFG references (`{2 2 2 0 0 0}`), `irrep` lists and
the `%rasci` module's QDPT/OPA couplings (`rel`/`douv`) are documented in
the same manual sections but not generated here.

## 29 Perturbed multistart batch (a converged reference mkl + its input)

**What it is for**: the inexpensive test that an SCF solution is not a
non-global minimum (the source: Vaucher & Reiher 2017, section IV).  The
converged reference's orbitals are perturbed by random occupied-virtual pair
mixings (10 random pairs out of the 15 highest occupied and the 15 lowest
unoccupied orbitals, random angles in [0, 90) degrees; their Eqs. (2)-(3)),
and each perturbed set becomes a restart.  A start either returns to the same
solution or finds another stationary point: a lower energy identifies the
reference as wrongly converged.  Stability analysis detects unstable
solutions (saddle points) but cannot distinguish local from global minima, so
this is complementary.  The menu generates the starts; the engine run stays
with you.

**What you need**: the reference's Molekel mkl (export it with
`orca_2mkl <base> -mkl`) and the ORCA input that produced it (the menu
rewrites that input's guess; its inline geometry is cross-checked against the
reference's coordinates).

**How**: menu 29 -> the reference mkl path; the base input path; the number
of starts (Enter = 3); the random seed (Enter = 20260927; it fixes every pair
and angle, so the batch replays byte-identically); the pairs per start
(Enter = 10, the source's value); the window per side (Enter = 15, the
source's value).  For each start you get
`<reference>.p<k>.fbk.mkl` (the perturbed orbitals) and
`<input>.p<k>.inp` (the base input with `!MORead` and
`%moinp "<reference>.p<k>.fbk.gbw"`), plus a report
`<input>.perturb.fbk.md` with every window, pair and angle.

**How to read it**:

- convert each orbital file next to itself (`orca_2mkl <name>.fbk -gbw`),
  run the inputs, and compare the converged energies with the reference's: a
  lower energy heals a wrongly converged reference; the same (or a higher)
  energy is no information -- the source is explicit that the test cannot
  guarantee detection;
- the perturbed columns are mixtures, not eigenfunctions: the occupations and
  orbital energies in the file are the reference's, and ORCA re-determines
  everything after reading the guess;
- an unrestricted reference is perturbed in both spins independently.

**Measured anchor** (CH4 dissociation, UKS PBE0/def2-SVP, one C-H at
2.6 Angstrom): the guess propagated from the equilibrium orbitals keeps the
restricted solution (-40.111230225721 Eh, <S**2> 0.000000); ORCA's default
guess also lands on a restricted solution (-40.183584827774 Eh, <S**2>
0.000000); all three generated starts reach the broken-symmetry solution
(-40.2416 Eh, <S**2> 0.971943) -- a healing of 0.130 Eh and a further
0.058 Eh below the default guess.

**Boundary**: the source validated the mechanism at the SCF level (it notes
that MC-SCF solutions can in principle be caught in local minima as well);
this menu edits any input's guess, so the same procedure applies to a CASSCF
reference, but no MC-SCF anchor is claimed here.

## 30 Imaginary-mode displacement (a frequency output + its input)

**What it is for**: the displacement stage of the automated
no-imaginary-frequency workflow (the source's NIFREC, ChemRxiv 2026; the
software is MIT-licensed and its protocol was read from the source).  A
converged optimisation that shows an imaginary frequency is not the
stationary point it claims to be; displacing the geometry along that mode
and re-running `opt freq` moves the search off the saddle, and the source's
success criterion is exactly the rerun: every frequency real.

**What you need**: the frequency output that carries the imaginary mode
(any Opt+Freq output -- the menu reads the last frequency block, the
normal-mode vectors and the run's own final geometry and masses) and the
base input of that run (its job keywords are kept; its coordinate block is
replaced).

**How**: menu 30 -> the frequency output path; the base input path; the
vector choice (`sum` = all imaginary modes, the source's default, or
`lowest` = the most negative mode alone); the displacement amplitude
(Enter = 0.1 Angstrom, the source's `base_disp`).  The menu writes
`<input>.disp_p.inp` and `<input>.disp_m.inp` -- the two signs of the
displacement, so both branches of the mode are available in one round --
plus a report `<output>.imagdisp.fbk.md` with the mode record.

**How to read it**:

- the two signs are an addition to the source (which grows one direction
  0.1, 0.2, ... 0.5 Angstrom); if the rerun still shows imaginary
  frequencies, escalate the amplitude along that schedule or switch the
  vector choice;
- the job keywords are the base input's own, so a plain `Opt` heals toward
  the minimum on each side of the mode while `OptTS` follows the mode; the
  `Freq` token is forced in, because the rerun's frequency check is the
  criterion;
- the displaced structure is a starting point, not a stationary point.

**Measured anchor** (the F + H2 fixture): the TS's imaginary mode
(-90.48 cm**-1) displaced +/-0.1 Angstrom (largest single-atom move
0.055 Angstrom, masses de-weighted with the run's own values) relaxes, both
signs, to geometries with no imaginary frequencies at -100.895757 Eh.  This
mode is the near-linear bend of the shallow saddle, so the two signs land on
the two mirror images of the same minimum -- the two-sided reaction
separation needs an asymmetric-stretch mode, which this fixture does not
carry.

## 31 CASSCF state data (a property file -> per-state energies and transitions)

**What it is for**: read the per-state table and the absorption spectrum out
of an ORCA property file (`<base>.property.txt`, written automatically by any
CASSCF job; the manual's property-file appendix lists the schema, and
`orca_2json <base> -property` gives the same data as JSON).

**How**: menu 31 -> the property file path.  You get the per-state table
(block, root, multiplicity, irrep, energy and the relative energy in eV) and
the absorption transitions (state pairs with irreps and multiplicities,
dE in eV and cm-1), plus a report `<property>.states.fbk.md`.

**How to read it**:

- the per-state energies reproduce the output's final `ROOT n: E=` lines to
  the printed precision (measured cross-check), and the transitions' first
  two columns satisfy the eV/cm-1 conversion (8065.544) and the state table's
  own energy differences;
- columns the manual's schema leaves unnamed are printed under an explicit
  note, never interpreted;
- a run without a spectrum request carries no `CASSCF_Absorption_Spectrum`
  section; the menu then reports the state table alone.

**State tracking boundary (the Wave-4.2 survey)**: CI vectors are not
persistable (run-time temporaries only; the `.cis` file belongs to the
CIS/STEOM modules), so state identity along a geometry series cannot be read
from one output.  What exists: the output's initial/final dominant-CSF
snapshots, this menu's structured per-state energies and transitions, and --
for identity fingerprints -- per-state dipoles from single-root runs (the
menu-19 route) or per-root densities from the FIC-NEVPT2 sidecar (the menu-23
chain); a tracker across a series built on those is a registered candidate
increment.  The state-averaged dipole in the property file is one x/y/z
vector (`State -1`), not a per-root table.

## 32 pysisyphus input (a structure + a method -> a PES-exploration YAML)

**What it is for**: write the input of a pysisyphus job -- minimum
optimisation or transition-state search -- that drives ORCA as its
calculator.  pysisyphus (GPL-3.0) is an external PES-exploration program
(RFO/GDIIS optimisers, growing-string / NEB / dimer routes); this menu covers
the two jobs the rescue chains of this toolkit use, and nothing is linked in:
the program stays external.

**How**: menu 32 -> a structure source (an XYZ file, an ORCA input with
inline coordinates, or an ORCA output whose final geometry is taken) ->
keywords, charge, multiplicity, job (min or ts), threshold, pal/mem, an
optional ORCA block string.  You get `<stem>.pysisyphus.xyz` (the structure),
`<stem>.pysisyphus.yaml` (the input) and a report `<stem>.pysisyphus.fbk.md`;
guidance covers the environment (`~/.pysisyphusrc` with `[orca5] cmd=...`),
the scratch disk (`$TMPDIR`) and the command line (`pysis <file>`, keep the
console).

**How to read it**:

- the keywords string is the `!` line of every ORCA call; the engine adds
  `engrad` itself for gradients, so any method ORCA supports can be driven;
- the TS job is written with a **model Hessian** (`hessian_init: fischer`);
  `hessian_init: calc` and `do_hess` are **refused/not offered**, because
  the quantum-Hessian route crashes with ORCA 6 (measured: ORCA 6 writes a
  `$multiplicity` block into its `.hess` file and the interface's grammar
  stops at it -- identical in pysisyphus 1.0.0 and on master).  Verify
  frequencies outside pysisyphus: a plain ORCA `Freq` job on the closing
  geometry (menu 30 is the toolkit's own route back from an imaginary mode);
- the thresholds are pysisyphus's own (`gau_loose` ... `baker`); the dummy
  `never` is not offered (it sets a billion cycles and disables the dump
  files).

## 33 pysisyphus run report (a run directory or console capture -> cross-checked record)

**What it is for**: read a finished pysisyphus run (**pysis** writes no
single result file; the record is the console plus a small set of artifacts).

**How**: menu 33 -> the run directory (or the console capture itself, e.g.
`pysis x.yaml > run.out`).  You get the run's own facts (program version,
job kind, system, calculator, charge/multiplicity, thresholds, the cycle
table), the outcome from the program's own marker (`Converged!` /
`Number of cycles exceeded!`), the closing energy and forces, the artifacts
found, and a report `<capture>.pysisyphus.fbk.md` next to the run.

**How to read it**:

- the outcome is the program's marker, never the force values: a run stopped
  by its cycle limit can close with forces already below the thresholds and
  is still "not converged" (measured on the stopping fixture);
- the cross-checks tie the artifacts together: the trajectory has one frame
  per cycle, the closing frame's energy equals its calculator call (ORCA's
  own outputs under `qm_calcs/`, parsed by the ORCA parser of this toolkit),
  the `Final summary` equals the run's last call -- which in a stopped run is
  the **extra closing evaluation** after the last cycle -- and the closing
  geometry file reproduces the closing state;
- a crashed run is classified: the measured ORCA-6 Hessian parse stop is
  recognised by its `$multiplicity` signature (with the `crashed_*` backup
  searched for the offending `.hess`), together with the way forward; other
  crashes point at the backup directory;
- boundaries: `optimization.h5` holds the full history but is HDF5 and is not
  read here (the text artifacts carry the same numbers), and
  `qm_calcs/cur_out` dangles after the run (a symlink into the cleaned
  scratch directory).

## 34 Judd-Ofelt intensity parameters (a transition dataset)

**What it is for**: the standard Judd-Ofelt (JO) fit of f-f intensities:
from a table of transitions (oscillator strengths + the host-independent
U^(lambda) matrix elements) obtain Omega_2/Omega_4/Omega_6 and the fit
quality; optionally, from the fitted parameters, the emission-side radiative
rates, branching ratios and radiative lifetime.

**How**: menu 34 -> a dataset file (YAML; schema in the formats chapter) ->
weighting (unweighted or the normalized 1/S variant).  You get the fitted
parameters (a.u. and the conventional 10^-20 cm^2), sigma and sigma/S_max,
the per-transition S_exp and r = S_theory/S_exp table, and a report
`<dataset>.jo.fbk.md`.

**How to read it**:

- the local field is the **squared** virtual-cavity form
  `chi_ED = (n_r^2+2)^2/9`; the source papers print it without the square —
  the square is what their reference code and their published numbers use
  (this menu's regression test).  Where a paper quotes parameters under a
  different local-field convention, expect differences of the order of that
  factor — state the convention before comparing numbers;
- the U^(lambda) values are table values for the ion (they do not depend on
  the host); transcribe them together with the intensities, from the same
  source, and check the fit against that source's published parameters;
- a host may be a constant refractive index or Sellmeier terms (B in um^2,
  wavelength cut at 6000 nm as the reference implementation); the measured
  wavelength drives both the Sellmeier lookup and the energy when no energy
  column is given;
- a parameter whose U^(lambda) column is zero over all transitions is
  reported as *not determined* (never as a silent zero);
- magnetic-dipole parts are not computed (they need free-ion wave
  functions): give `f_md` per transition where they are significant — the
  fit subtracts them, as the source implementations do.

**Boundaries**: the extended (perturbative X_k / configuration-interaction)
JO models of the 2022/2024 literature need free-ion atomic-structure wave
functions and are outside this toolkit's post-processing scope; the ORCA-side
data for *ab initio* intensities are the per-transition fosc/D2/dipole
columns of the property file's absorption section (read by menu 31; layout
measured).

## 35 pNMR pseudocontact shifts (a susceptibility run file + a structure)

**What it is for**: the pseudocontact shift (PCS) of every nucleus around a
paramagnetic centre, in the point-dipole approximation: from a magnetic
susceptibility tensor (given directly, built from a g-tensor, or built from
g + zero-field splitting at the run temperature) and a structure, the menu
prints the tensor, its axiality and rhombicity, and the per-nucleus PCS
(ppm) with the nucleus's r, theta, phi in the tensor's eigenframe.

**How**: menu 35 -> a run file (YAML; schema in the formats chapter) that
names the structure, the 1-based centre atom, the temperature and the
susceptibility source.  The source may be a 3x3 tensor (m^3, or the
literature's 1e-32 m^3 unit), a g-matrix, or g + D/E/D -- each of these may
also be read straight from an ORCA QDPT output (`orca_output:`); the menu's
ORCA reader takes the effective-Hamiltonian g-matrix and the preferred ZFS
variant (effective Hamiltonian with the spin-spin contribution when
present).  The report is `<runfile>.pnmr.fbk.md`.

**How to read it**:

- the two susceptibility constructions are both available by name: the
  non-symmetric `chi' = mu_B^2 mu_0 g_e/(k_B T) g.<SS>` (the default; the
  one consistent with the modern pNMR shielding theory) and the symmetric
  `g.<SS>.g^T` alternative -- for the source's Co(II) benchmark they give
  Delta chi_ax = 14.8 vs 27.3 (1e-32 m^3).  State which one produced a
  number before comparing with literature;
- only the traceless symmetric part of chi enters the PCS (adding an
  isotropic part, or the antisymmetric part of a non-symmetric tensor,
  changes nothing -- tested exactly); the axial limit is exactly the
  classical `(1/12 pi r^3) Delta chi_ax (3 cos^2 theta - 1)` form;
- the point-dipole form is the long-range limit: the source's benchmark
  shows below-10 % deviation of the dipolar part beyond ~7 A, while within
  ~4-5 A the contact (through-bond) term dominates -- a PCS for a nucleus
  close to the metal is not a prediction of the total shift;
- the ZFS route uses the effective-spin dyadic <SS> at the given
  temperature; near-axial systems can hit the node condition
  `<SS>|| g|| = <SS>perp gperp` where the axiality (and the PCS) vanish;
- ORCA cross-check: for an ORCA QDPT output with `DoSusceptibility`, the
  menu's chi (converted: cgs-emu molar chi*T = NA chi_SI/(4 pi) x T x 1e6)
  reproduces ORCA's own printed susceptibility to the printed precision.

**Boundaries**: the contact shift is not computed (it needs hyperfine
coupling); Bleaney's analytic anisotropy theory is not implemented as a
construction (the ligand-field review's limits of it frame the menu's
documentation; a Bleaney comparator is a registered candidate increment);
the OpenMolcas-side g/chi data belong to the OpenMolcas chain (menus 1.4/
0.5) when that deployment lands.

## 36 Magnetic relaxation and QTM (an ORCA output with SINGLE_ANISO or MAGRELAX)

**What it is for**: the static QTM metrics and the relaxation-rate table
that an ORCA calculation already contains, summarised against the practical
guide's criteria -- per-Kramers-doublet g-tensors and their mutual
orientation (the empirical barrier estimate), the UBAR magnetic-moment
matrix elements, the spin-orbit level structure, and, when the
Orca_Magrelax section is present, the tau(T) table with an Arrhenius fit.

**How**: menu 36 -> one file name: an ORCA output carrying either the
SINGLE_ANISO embedded section (the `ANISO` sub-block of `%casscf`; manual
section 5.32) or the `* ORCA MAGRELAX *` section (manual section 5.30).  A
file with both gets both parts of the report, written as
`<output>.relax.fbk.md`.

**How to read it**:

- the per-KD table lists each pseudospin group's energy, g principal
  values, the angle theta between its largest-g axis and the ground
  group's, and the guide's empirical `g_T = (g1 + g2 + g3 sin(theta_3))/3`
  with its product `g_T * theta`.  The barrier estimate is the first
  excited group with theta > 15 deg (non-collinear) or `g_T theta` > 20
  (the guide's rule from a survey of 20 Dy(III) SMMs); both thresholds
  are **empirical** and are restated with every report;
- groups whose relative g spread is below 1 % are marked `*`: their
  main-axis direction is decided by the print's digits, so the axis-angle
  criteria do not apply to them (measured on the CO+ fixture, where
  Delta g/g is about 2e-4);
- the UBAR block gives the engine's Zeeman eigenstates and magnetic-moment
  matrix elements `(|mu_X| + |mu_Y| + |mu_Z|)/3`; the engine itself calls
  the printed barrier "only a qualitative relaxation path", and its
  "even number of electrons" sentence is a fixed template warning (it is
  printed for odd-electron systems too);
- the magrelax part tabulates tau(T) in seconds and fits
  `tau = tau_0 exp(U_eff/k_B T)` over the finite positive points (at least
  three are required).  The CO+ fixture has an all-zero rate table: with a
  single vibrational mode (2299.9 cm-1) no phonon matches its Zeeman gaps,
  so there is no one-phonon channel -- an all-infinity column is a
  structure fact, not a fit target;
- with `DoSSC true` the output contains two complete SINGLE_ANISO segments
  (SOC-only and SOC+SSC; the O2 fixture prints D = 2.175287 vs 3.185964) --
  compare the spectra before mixing numbers across segments.  `MLTP` must
  be given explicitly: without it the engine prints no g/D analysis at all
  (both the g-tensor and the barrier estimation read NOT INCLUDED).

**Boundaries**: nothing here is a dynamical simulation -- the guide's own
section 7.3 calls the UBAR-permutation picture "not realistic", and the
empirical thresholds bound a plausible barrier rather than a measured one;
the magnetic-dilution tau_QT model of the Aravena group (dipolar-field
statistics; the source of the dilution design rules) is a registered
candidate increment pending its closed-access papers; polynuclear exchange
and the POLY_ANISO route belong to the multinuclear item (5.6).

## 37 Core-excited spectra XAS/RIXS (a ROCIS output)

**What it is for**: the core-excited absorption spectra an ORCA ROCIS
calculation already contains (the transition-metal L-edge protocol of the
manual's section 5.7, and the CAS-CI/RAS-CI protocol of section 3.13.18):
the transition table of the best available spectrum block, the edge
branching ratio of the spin-orbit-split pair, and the run's RIXS
bookkeeping.

**How**: menu 37 -> the ORCA output path of a ROCIS run (`%rocis` with
`DoGenROCIS`), then optionally a statistical branching ratio for the
comparison (Enter to report the ratio only).  The report is
`<output>.xas.fbk.md`.

**How to read it**:

- the primary block is the SOC-corrected electric-dipole table when the run
  has `DoSOC` (its fosc column is **weighted by the initial-state
  population**, as the engine's column header says), otherwise the plain
  electric-dipole table; with `DecomposeFosc` up to seven blocks exist
  (electric and velocity dipoles plus five combined D2/M2/Q2 variants);
- the transition table lists the non-zero-fosc transitions (up to 80; the
  SOC-corrected tables hold thousands of state pairs, most of them zero --
  the report says how many);
- the **branching ratio** splits the transitions at their largest energy
  gap (the natural two-cluster default for a spin-orbit-split edge) and
  reports each cluster's centroid and summed oscillator strength and the
  low/high ratio.  The ratio is a data fact; the deviation from the
  statistical value (2:1 for a 2p edge, 3:2 for a 3d edge -- from the
  (2j+1) degeneracies) indicates electrostatic and spin-orbit effects
  (Thole & van der Laan 1988, whose rules this menu cites rather than
  reimplements).  Pass your reference value as the menu's second answer to
  get the ratio/stat column;
- the RIXS bookkeeping states which of the three situations the run is in:
  not requested; the engine's **refusal mode** (the RIXS flags with zero
  intermediate/final states -- measured on a 4-element OrbWin, which lacks
  the second donor space the RIXS variant needs); or cross sections
  evaluated, with the transition counts and the `orca_mapspc` recipe
  (manual section 5.7.4.2) for the 2D data files.

**Boundaries**: ROCIS applies several approximations (the manual says the
results are qualitatively correct); the menu does not name the edge
clusters (their assignment depends on the orbital windows), it reports the
low/high clusters and the ratio.  The six-element `OrbWin` semantics of a
successful RIXS run were not established in this round's probing (the
measured variants were refused or crashed); generating RIXS inputs is a
registered candidate increment.  The CAS-CI/RAS-CI XAS protocol
(section 3.13.18: rotate the core orbitals in, `FrozenCore FC_NONE`,
`maxiter 1`) is documented in the manual chapter but this menu reads ROCIS
outputs only.

## 38 Ab initio ligand-field analysis (an AILFT output)

**What it is for**: the ligand-field parameters ORCA's AILFT module has
fitted to an ab initio (SA-CASSCF/NEVPT2) effective Hamiltonian (manual
section 3.13.16): the ligand-field one-electron eigenfunctions (the LF
splitting), the Slater-Condon and Racah parameters at each theoretical
level, the fit quality and the SOC constant.

**How**: menu 38 -> the ORCA output path of a CASSCF run with the AILFT
driver (`%casscf` with `ActOrbs dOrbs`/`fOrbs` or `LFTCase 3d`/`4f` ...),
then optionally the free-ion Racah B and/or SOC constant zeta0 for the
nephelauxetic ratios (Enter skips).  The report is `<output>.ailft.fbk.md`.

**How to read it**:

- one block per level (CASSCF and NEVPT2): the ligand-field eigenfunction
  energies (cm-1) with their spread, the Slater-Condon parameters (F0 is
  marked "(fixed)" when taken from the raw two-electron integrals; f shells
  carry F6), the Racah parameters B, C and C/B, and the `*.lft.gbw` file
  name (openable with `orca_plot`);
- the fit quality: total and per-block RMS errors and Pearson's
  correlation.  The near-zero CASSCF-level RMS is intrinsic -- the LFT
  parametrization is exact for that level -- while the correlation-level
  RMS reflects, among other things, the neglected anisotropy of
  electron-electron repulsion in covalent complexes (Lang, Atanasov &
  Neese 2020);
- the SOC constant (ZETA_D for d shells, ZETA_F for f shells; the engine
  fits it against the CASSCF-orbital SOC matrix elements);
- when you pass free-ion references, the nephelauxetic ratio
  beta = B/B0 and the relativistic ratio zeta/zeta0 -- reductions are the
  classic covalency indicators (Jung, Atanasov & Neese 2017, the
  actinide/lanthanide AILFT reference).

**Boundaries**: this menu reads the engine's fit and never refits; the
parameters are model quantities of the ligand-field Hamiltonian (their
interpretation is the ligand-field model's business).  LFDFT (ligand-field
DFT) lives inside ADF, a commercial package that is not available here --
registered as a termination: the AILFT route covers the same analysis
needs from open programs.

## 39 Polynuclear magnetism (a POLY_ANISO output)

**What it is for**: the cluster magnetic analysis ORCA's POLY_ANISO driver
produces from single-ion ab initio data (manual section 7.18): the
exchange-coupled states of a polynuclear complex, the interaction-matrix
decomposition, chiT(T) and the Van Vleck susceptibility tensors.

**How**: menu 39 -> the `poly_aniso.output` path of a run of
`$ORCA/otool_poly_aniso < poly_aniso.input > poly_aniso.output` (the driver
is called independently of ORCA).  Each magnetic center is first treated as
an isolated fragment (`CASSCF/NEVPT2 + SOC + SINGLE_ANISO`, the menu-36
data); the per-center data files (`<job>.CASSCF.anisofile`) must be placed
as `aniso_1.input`, `aniso_2.input`, ... -- the names are mandatory
(measured).  The report is `<output>.polyaniso.fbk.md`.

**How to read it**:

- the per-center echo (data file, coordinates, spin-orbit states and
  spectrum, g values) -- a cross-check against the menu-36 analysis of the
  same data;
- the exchange block: the coupled-state count, the pair list with the J
  values you supplied, and what models were included (Lines-1,
  dipole-dipole, the ITO decomposition);
- the first-order anisotropic coupling: per pair and per model the full 3x3
  interaction matrix and its decomposition into isotropic / symmetric /
  anti-symmetric terms with weights;
- the coupled-state table (Lines / dipole-dipole / total, absolute and
  relative cm-1), the population analysis and the expectation-value tables
  (moments per exchange state and center);
- chiT(T) (101 points by default) and the Van Vleck susceptibility tensor
  sequence (one 3x3 with main values and main axes per printed temperature).

**Boundaries**: the exchange constants J are your input (measured elsewhere
or fitted) -- this workflow never computes or fits them; a joint ab initio
computation of exchange splittings (the LDF-CAHF / many-state PNO-CASPT2
route) is outside the ORCA ecosystem and is registered as a documented
termination.  The Lines model is exact only for two isotropic spins, one
Ising plus one isotropic spin, or two Ising spins, and approximate
otherwise; the dipole-dipole coupling is evaluated exactly from the ab
initio moments and usually dominates in strongly anisotropic lanthanides.

## 40 Hyperfine and EFG parameters (an EPRNMR output)

**What it is for**: the electric and magnetic hyperfine structure ORCA's
`%eprnmr` prints (manual section 7.51.3): per nucleus the A-tensor
components, the electric field gradient with its electron/nuclear
decomposition, and the density at the nucleus -- the quantities behind
Mossbauer quadrupole splittings and hyperfine parameters.

**How**: menu 40 -> the ORCA output path (a run with `%eprnmr` and a
`Nuclei` list; the coordinates block must come before `%eprnmr` --
measured), then optionally the nuclear quadrupole moment Q (barn) for the
quadrupole-splitting conversion (Enter skips).  The report is
`<output>.hyperfine.fbk.md`.

**How to read it**:

- per nucleus the nuclear parameters (I, P, Q as printed) and the A tensor
  (A(iso), A(Tot), the FC/SD split) -- the MHz values need the nuclear
  parameters in the input; zeros mean they were not supplied (measured);
- the EFG: the V(Tot) principal values, |Vzz|, the asymmetry parameter
  eta = (Vxx - Vyy)/Vzz, the principal axes, and the V(El)/V(Nuc)
  decomposition with the sum cross-check against V(Tot);
- Rho(0), the density at the nucleus (the isomer-shift core quantity) --
  basis-domain dependent: all-electron bases carry the core density, ECP
  bases report the valence density only; never compare across domains;
- with Q supplied: eQVzz/h and the first-order quadrupole splitting
  DeltaE_Q = eQVzz/2 sqrt(1 + eta^2/3).  The conversion constant is
  CODATA-derived (234.9648 MHz per barn per a.u.; the engine's own
  quadrupole block agrees to ~1e-6).

**Boundaries**: the engine computes the A components on the DFT routes;
on CASSCF it switches them off while keeping EFG/Rho(0) (measured).  The
EFG is sensitive to SCF micro-solutions (measured: two runs of the same
input differing by 2e-6 Eh in energy moved V(Tot) by ~1e-4 relative) --
converged energy does not imply a reproducible EFG.  Core-level
photoemission (XPS) has no ORCA module (binary search: zero hits) --
registered as a termination; the XES/SOC-XES modes seen in the
relativistic-CASSCF strings belong to the ROCIS family (menu 37's domain).

## 41 Magnetic entropy and magnetocaloric effect

**What it is for**: the magnetic entropy S(T) and the isothermal entropy
change -DeltaS(T, H) -- the thermodynamic figures of merit of the
magnetocaloric effect (molecular coolers).

**How**: menu 41 -> the ORCA output path.  The route is chosen by content:
a POLY_ANISO output with the HINT/TMAG magnetization table goes through
the Maxwell relation; a SINGLE_ANISO output (or the per-center spectra of
a POLY_ANISO run) goes through the partition function of the spin-orbit
levels.  The report is `<output>.mce.fbk.md`.

**How to read it**:

- levels route: S(T) = R (ln Z + <E>/kT) at B = 0 from the printed
  spin-orbit levels (each printed line counts as one state, so the
  engine's level list handles degeneracies); the high-temperature check
  R ln(N) is printed -- the level list is truncated, so it is an upper
  bound: do not trust S above the saturation temperature;
- magnetization route: -DeltaS(T, H) from the Maxwell relation
  (dS/dH)_T = (dM/dT)_H, integrated over the printed field grid with
  adjacent-temperature differences; positive values are the direct
  magnetocaloric effect (the cooling capacity quoted in the literature);
  the probe maximum is reported.  Converged data (a fine temperature grid)
  are required.

**Boundaries**: both routes are pure post-processing of printed values;
the levels route is B = 0 (the levels are field-free) and the
magnetization route uses the powder-averaged molar table (Bohr magnetons).
Method reference: Szalowski & Kowalewska 2020 (conventions; the V6
companion study 2020 uses the same scheme).

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
