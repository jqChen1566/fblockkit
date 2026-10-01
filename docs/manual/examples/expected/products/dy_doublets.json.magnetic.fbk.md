## A10 magnetic-doublet criterion (g_T * theta_3)

Magnetic-doublet criterion (g_T = (g1 + g2 + g3 sin(theta3)) / 3, degrees; line g_T * theta3 = 20):
  #   doublet                             E(cm-1)      g1      g2      g3  theta3(deg)      g_T  g_T*theta3  reading
  0 * doublet that supports excitation          -    0.41    0.44    9.04         9.52    0.782        7.44  supports excitation
  1   doublet that facilitates QTM              -    3.45    3.73    5.96        17.58    2.993       52.62  facilitates QTM
  (* = the reference row; no angle needed measuring from an axis)

Criteria and boundaries:
  - the criterion: g_T = (g1 + g2 + g3 sin(theta3)) / 3 with theta3 in DEGREES, compared with g_T * theta3 = 20: below the line the doublet can still act as a step of the barrier (supports excitation above it), at or above it the doublet facilitates quantum tunnelling
  - the three principal values are sorted internally, the largest taking the sin term (the source's convention); a 'axis3' entry is the axis of that largest value
  - every theta3 is used as given (as a published table measures it, against the ground doublet it refers to); entry 0 ('doublet that supports excitation') is marked as the table's reference row -- set 'reference' to mark another entry
  - an angle computed from axes is folded to [0, 90] degrees (a principal axis's sign is arbitrary); theta3 given in the table is used as given
  - calibration context: on the source's 19-system table the line sits in the gap between the largest 'safe' product (15.31) and the smallest 'definite QTM' product (30.25)
  - The metric is empirical, not derived: the source states that the g1/g2 weights should in principle be reduced as g3 rotates into the plane, but those historic calculations do not carry that information. The line's physical background is the material's internal dipolar field, which varies with packing and dilution, so the line is not a universal constant.
  - Applicability domain: the line was calibrated on 19 mononuclear Dy(III) complexes (2016-2025, one consistent methodology -- 9-in-7 active space, ANO-RCC three-tier basis, SA-CASSCF-SO in (Open)Molcas). Other ions or nuclearities would need their own calibration.
  - Known failure mode the criterion cannot see: in the source's set exactly one doublet is misread -- [Dy(Cp^ttt)2]+'s QTM doublet falls deep in the safe zone (product 2.52), and the source offers no explanation. A 'safe' verdict here is therefore the criterion's reading, not a proof of barrier behaviour.

## References

Complete citations:
- [chilton2025abinitio] Chilton, N. F. (2025). Ab initio electronic structure calculations of lanthanide single-molecule magnets; a practical guide. Chemical Society Reviews, 54(24), 11468-11487. DOI: 10.1039/d5cs00493d

BibTeX (paste-ready):

```bibtex
@article{chilton2025abinitio,
  author  = {Chilton, Nicholas F.},
  title   = {Ab initio electronic structure calculations of lanthanide single-molecule magnets; a practical guide},
  journal = {Chemical Society Reviews},
  year    = {2025},
  volume  = {54},
  number  = {24},
  pages   = {11468--11487},
  doi     = {10.1039/d5cs00493d},
}
```
