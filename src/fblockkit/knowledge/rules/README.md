# Rule data format (knowledge/rules/)

This directory holds fBlockKit's rule data (YAML), loaded and validated by
`fblockkit.knowledge.loader`. Every rule must carry its provenance (a non-empty
`evidence`), otherwise loading rejects it outright -- that is how "rules are
auditable" is implemented.

**Writing convention**: `evidence.text` always uses the block scalar (`>-`). A
YAML double-quoted scalar ends at the closing quote, and nothing else may follow
on that line (including non-ASCII punctuation) -- merge several quotations into
one block scalar instead, which avoids the parsing trap.

## Rule structure

```yaml
- id: A2-f-block-energy-nevpt2        # unique id (a duplicate is an error)
  kind: recipe                        # recipe | diagnosis
  title: a one-line title
  condition:                          # see "condition syntax" below
    all:
      - { field: f_block, op: eq, value: true }
      - { field: targets, op: contains, value: energy }
  action: the action / recommendation when the rule matches (written for the user)
  refusal:                            # optional; statements for a match that must be refused
    - A correlation-layer method must not be put into an optimisation loop (no analytic gradient).
  evidence:                           # required, non-empty; each entry has kind/text/ref (url optional)
    - kind: manual                    # manual | literature | measured
      text: "...(a manual sentence, a literature conclusion, or a measured record)"
      ref: ORCA 6.1 manual §3.17
      url: https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/caspt2.html
  confidence: confirmed               # confirmed | provisional
  severity: warn                      # required for diagnosis rules only (info | warn | error | refuse);
                                      # recipe rules must not carry it; enforced at load time
```

**severity semantics** (diagnosis rules): `info` = a check that must be carried
out (no concrete problem found yet); `warn` = a concrete doubt was found;
`error` = the conclusion cannot be relied on; `refuse` = refusing to give a
conclusion (insufficient data, and so on).

## Condition syntax

- Single test: `{ field: <fact field>, op: <operator>, value: <value>}`
- Conjunction: `{ all: [<condition>, ...] }` (an empty `all: []` is always true --
  used for an unconditional preliminary note)
- Disjunction: `{ any: [<condition>, ...] }` (an empty `any: []` is always false)
- Operators: `eq ne in not_in contains gt ge lt le`; a missing field evaluates to
  None (no error)

## Fact field vocabulary

**Recipe rules** (derived by the recipe layer from `SystemProfile`):

| Field | Type | Meaning |
|---|---|---|
| `f_block` | bool | the system has an open 4f/5f shell |
| `targets` | list | target list: energy / geometry / excited / magnetic / spectra |
| `geometry_task` | bool | the targets include a geometry optimisation or frequencies |
| `soc_task` | bool | the targets require spin-orbit coupling |
| `active_space_orbitals` | int or null | number of active-space orbitals planned |
| `open_shells` | int | number of open shells |

**Diagnosis rules** (extracted by the parser layer from the output file; the
producing function is `parsers.facts_from`):

| Field | Type | Meaning |
|---|---|---|
| `terminated_normally` | bool | did the program terminate normally (False = aborted/crashed, numbers unusable) |
| `caspt2_present` | bool | the output contains a CASPT2 calculation |
| `soc_present` | bool | the output contains an SOC calculation |
| `scf_cycles` | int | number of cycles the SCF needed to converge |
| `scf_converged` | bool | did the SCF converge (False when unconverged and not aborted) |
| `casscf_present` | bool | the output contains a CASSCF calculation |
| `casscf_converged` | bool | did CASSCF report convergence (energy or gradient criterion; the two are recorded separately in the parse result) |
| `caspt2_min_reference_weight` | float | smallest CASPT2 Reference Weight over the states (manual recommends >0.9) |
| `caspt2_min_denominator` | float | smallest CASPT2 smallest-energy-denominator (a small value points to intruder-state risk) |
| `casscf_converged_via` | str | source of the CASSCF convergence marker: `energy` (energy criterion) or `gradient` (gradient criterion) |
| `casscf_active_occ_min` | float | smallest natural occupation in the active space (<0.02 counts as a nearly empty orbital) |
| `casscf_active_occ_max` | float | largest natural occupation in the active space (>1.98 counts as a nearly full orbital) |
| `nevpt2_max_hole_particle` | float | largest V1_i (ITUV) and Vm1_a (TUVA) class contribution over the NEVPT2 states (>0 is a sign of an intruder state) |
| `ts_optimization` | bool | the job is a transition-state optimisation (an OptTS marker appears in the output, or OptTS in the echoed input keywords) |
| `optimization_converged` | bool | outcome of the (last) geometry optimisation; absent when the output contains no geometry optimisation |
| `frequency_present` | bool | the output contains a vibrational-frequency analysis (a VIBRATIONAL FREQUENCIES block) |
| `frequency_after_geometry` | bool | the last frequency block comes after the last geometry-optimisation cycle (i.e. it describes the optimised structure); absent without a frequency block |
| `frequency_imaginary_count` | int | number of imaginary modes in the last frequency block |
| `frequency_min_imaginary` | float | the most negative wavenumber (the strongest imaginary mode, negative by definition); only present when at least one imaginary mode exists |

A field that is missing from the file does not enter the fact table; False is a
valid fact and is kept. When adding a field: register it in this table first,
then add an evaluation case to `tests/test_knowledge.py`.
