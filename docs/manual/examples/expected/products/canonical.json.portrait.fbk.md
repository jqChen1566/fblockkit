## A9 orbital portrait (descriptor panel)

Orbital portrait (occupation, energy, dominant centre, angular-momentum composition, bonding label, AO-centre extent):
  orbital   occ        E(Eh)   centre        share   bonding        shell composition            extent(A)  shift(A)
        4   1.9934    -0.70408   N1      50.0%   antibonding    d0% p36% s64%                   0.547    0.547
        5   1.9363    -0.59891   N0      50.0%   bonding        d2% p98% s0%                    0.547    0.547
        6   1.9363    -0.59891   N0      50.0%   bonding        d2% p98% s0%                    0.547    0.547
        7   0.0660     0.26659   N1      50.0%   antibonding    d0% p100% s0%                   0.547    0.547
        8   0.0660     0.26659   N1      50.0%   antibonding    d0% p100% s0%                   0.547    0.547
        9   0.0019     0.86579   N1      50.0%   bonding        d4% p30% s66%                   0.547    0.547

Notes:
  - atom pairs within 6 Angstrom enter the bonding label; a single atom carrying more than 95% of the Loewdin population makes it non-bonding, and a cross term within 0.01 counts as non-bonding too (the source's criterion plus a band that its own examples never needed: a symmetric core orbital here measures 3e-4)
  - the population sums check: the Loewdin atomic populations of each orbital add up to one within round-off
  - a row marked '!' is cancellation-heavy (the absolute sum of its atom-pair blocks exceeds 5; the N2 export's orbitals run up to 2.2 while its one node-heavy virtual measures 158): the sign of its cross term still says what it is, but the magnitude is not a bond order
  - Two of the source's descriptors need integrals an orbital export does not carry (the diagonal one- and two-electron integrals and the exact <r^2>): the extent and the centroid displacement here are AO-centre estimates, named as such, and the diagonal integrals are available from an FCIDUMP (menu 12's route) when they are needed.
  - The angular-momentum composition is the weighted form of the source's binary shell encoding -- the reading notes record that the binary version is deliberately coarse and that the weighted one is what the shell-classification checks need.
  - The neural predictor and the threshold policy of the source are not adopted here (the project runs a deterministic rule engine).
