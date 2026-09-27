## 3.5 TNASS selection (Renyi-2 bipartition)

TNASS active-space selection (Renyi-2 entropy of the bipartition):
  method: greedy; target size 4 of the window's spatial orbitals

  orbital  S2(single)  (the seed ranking)
        0     0.06769
        1     0.14837
        2     0.14837
        3     0.14977
        4     0.14977
        5     0.06493

  step  added  subset                       S2(A)
     1      3  3                          0.14977
     2      4  3, 4                       0.23280
     3      5  3, 4, 5                    0.24044
     4      0  0, 3, 4, 5                 0.23000

  selected: orbitals [0, 3, 4, 5] (2e, 4o)
  exact window FCI energy: -230.793818898 Eh

Boundaries and checks:
  - the S2 oracle is exact over the FCIDUMP window (the determinant CI of the four-state-entropy route); the source builds it as a tensor-network entanglement feature over the full orbital space with an approximate MPS (its useful bond dimensions are 4-6) -- there is no bond-dimension truncation here, but the window boundary is the restriction instead
  - the selection maximizes the bipartition entanglement with the window complement; the source's best-k variant (choose k by the CASCI energy) is not implemented -- run the menu at several target sizes and compare
  - the subset is a spatial-orbital set (both spins travel together), matching the source's selection domain
  - delivering the space to a CASSCF needs the orbital-order machinery (menu 22's boundary note)

## References

Complete citations:
- [mingare2026tnass] Mingare, A.; Heuz{\'e}, I.; Coveney, P. V. (2026). {TNASS}: Tensor Network Active Space Selection with the Entanglement Feature. arXiv preprint arXiv:2608.03645. DOI: 10.48550/arXiv.2608.03645

BibTeX (paste-ready):

```bibtex
@article{mingare2026tnass,
  author  = {Mingare, Angus and Heuz{\'e}, Isabelle and Coveney, Peter V.},
  title   = {{TNASS}: Tensor Network Active Space Selection with the Entanglement Feature},
  journal = {arXiv preprint arXiv:2608.03645},
  year    = {2026},
  doi     = {10.48550/arXiv.2608.03645},
}
```
