## AEGISS selection (entropy + AO projection)

AEGISS active-space selection (entropy screening + atomic-orbital projection):
  system: benzene; window 6 orbitals, AO label 'C pz' (12 target functions)
  entropy line: S > 0.03417 = 0.1 * S_max; projection threshold: w > 0.5 (projection norm)

    MO  occupation         S   entropy    weight  selected
    18     1.96074    0.1813      keep    3.0874       yes
    19     1.90108    0.3393      keep    1.8672       yes
    20     1.90108    0.3393      keep    1.8672       yes
    21     0.10026    0.3417      keep    0.9920       yes
    22     0.10026    0.3417      keep    0.9920       yes
    23     0.03659    0.1743      keep    0.8096       yes

  final active space: (6e, 6o) -- 6 orbital(s) passed the entropy screen, 6 of those passed the projection
  exact window FCI energy: -230.793818898 Eh  (engine CASSCF print: -230.793818898; difference 1.72e-10)

Boundaries and checks:
  - the entropies are exact (the four-state-entropy route over the FCIDUMP window); the source estimates them from a DMRG over a larger window
  - the projection weight is the projection norm sum_eta |O[eta, p]|^2 -- a measured deviation from the source's signed row sum, which cancels by symmetry for nodal pi orbitals (on the benzene fixture only a2u survives the signed sum, while the norm separates sigma (0.0000) from pi cleanly)
  - the target is the calculation's own AO subset (the non-minimal AVAS route); omitting the shell index from the label matches every shell of that (element, angular, component) family -- the analogue of the source's minimal-basis label
  - the cluster labeling of the source's step 1 is folded into the AO label; multi-group unions are run one group per menu pass
  - the deliverable is the selection over the dumped window; feeding it to a CASSCF needs the orbital-order machinery (menu 22's boundary note)

## References

Complete citations:
- [tarocco2026aegiss] Tarocco, F.; Haase, P. A. B.; Pavo{\v{s}}evi{\'c}, F.; Krishna, V.; Guidoni, L.; Knecht, S.; Stella, M. (2026). {AEGISS} -- Atomic Orbital and Entropy-based Guided Inference for Space Selection. arXiv preprint arXiv:2508.10671. DOI: 10.48550/arXiv.2508.10671

BibTeX (paste-ready):

```bibtex
@article{tarocco2026aegiss,
  author  = {Tarocco, Fabio and Haase, Pi A. B. and Pavo{\v{s}}evi{\'c}, Fabijan and Krishna, Vijay and Guidoni, Leonardo and Knecht, Stefan and Stella, Martina},
  title   = {{AEGISS} -- Atomic Orbital and Entropy-based Guided Inference for Space Selection},
  journal = {arXiv preprint arXiv:2508.10671},
  year    = {2026},
  doi     = {10.48550/arXiv.2508.10671},
}
```
