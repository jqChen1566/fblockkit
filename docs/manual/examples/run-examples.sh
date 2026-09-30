#!/usr/bin/env bash
# Run the manual's example sessions and refresh their captured outputs.
#
# Discipline: every session transcript and report excerpt printed in the manual
# comes from these runs -- nothing is typed by hand.  The captures live in
# expected/ and are replayed by tests/test_manual.py.
set -euo pipefail
cd "$(dirname "$0")"
REPO=../../..
export PYTHONUTF8=1
export PYTHONPATH="$REPO/src"
# the menu-4 datasource leg queries the deployed basis library; point it at
# the shipped miniature so the capture is byte-identical on both ends
export FBK_BASISDB="$REPO/fixtures/basisdb/mini_basis.db"
PY=python

rm -rf work
mkdir -p work expected expected/products

"$PY" generate_inputs.py

cp "$REPO/fixtures/orca/n2_casscf_nevpt2.out" work/
cp "$REPO/fixtures/orca/generated_ce3_sarc2.out" work/
cp "$REPO/fixtures/orca/n2_stretch_local_spin.out" work/
cp "$REPO/fixtures/orca/fhh_optts_freq.out" work/
cp "$REPO/fixtures/orca/inputs/fhh_reopt.inp" work/
cp "$REPO/fixtures/orca/n2_hf_clean.out" work/
cp "$REPO/fixtures/orca/n2_diffuse.out" work/
cp "$REPO/fixtures/orca/scf_noconv.out" work/
cp "$REPO/fixtures/orca/inputs/scf_noconv.inp" work/
cp "$REPO/fixtures/literature/pucl3_s18.json" work/
# the exact-entropy chain of menu 12 (converged CASSCF output, its FCIDUMP and
# the two orca_2json exports of the gbw)
cp "$REPO/fixtures/orca/n2_fcidump_step_a.out" work/
cp "$REPO/fixtures/orca/n2_fcidump.fcidump" work/FCIDUMP
cp "$REPO/fixtures/orca/n2_fcidump.canonical.json" work/canonical.json
cp "$REPO/fixtures/orca/n2_fcidump.localized.json" work/localized.json
# the cross-structure mapping chain of menu 17 (four localized scan exports,
# the canonical export of the template geometry, and the template mkl of menu 18)
cp "$REPO/fixtures/orca/n2_scan_1.094.loc.json" work/
cp "$REPO/fixtures/orca/n2_scan_1.600.loc.json" work/
cp "$REPO/fixtures/orca/n2_scan_1.610.loc.json" work/
cp "$REPO/fixtures/orca/n2_scan_2.600.loc.json" work/
cp "$REPO/fixtures/orca/n2_scan_1.600.json" work/
cp "$REPO/fixtures/orca/n2_scan_1.600.mkl" work/
# the dipole-moment chain of menus 19/20 (six single-root candidates + reference)
cp "$REPO/fixtures/orca/h2o_dm_ref_pbe0.out" work/
cp "$REPO/fixtures/orca/h2o_dm_prep_mp2.out" work/
cp "$REPO/fixtures/orca/h2o_dm_casci_e6o6.out" work/
cp "$REPO/fixtures/orca/h2o_dm_casci_e6o7.out" work/
cp "$REPO/fixtures/orca/h2o_dm_casci_e6o8.out" work/
cp "$REPO/fixtures/orca/h2o_dm_casci_e8o7.out" work/
cp "$REPO/fixtures/orca/h2o_dm_casci_e8o8.out" work/
cp "$REPO/fixtures/orca/h2o_dm_casci_e10o8.out" work/
cp "$REPO/fixtures/orca/h2o_dm_casci_e6o6_sa4.out" work/
# the APC ranking export of menu 21 (N2/def2-SVP RHF with the Fock family and
# the windowed MO_IAJB block)
cp "$REPO/fixtures/orca/n2_apc.json" work/
# the ASS1ST round export of menu 23 (N2/def2-SVP CAS(6,6) + FIC-NEVPT2 with
# the unrelaxed NEVPT2 density in the sidecar)
cp "$REPO/fixtures/orca/n2_ass1st.json" work/
# the AEGISS benzene platform of menu 25 (the pi-window CASSCF export, its
# FCIDUMP and the run output for the energy cross-check)
cp "$REPO/fixtures/orca/benzene.json" work/
cp "$REPO/fixtures/orca/benzene.fcidump" work/
cp "$REPO/fixtures/orca/benzene.out" work/
# the menu-29 reference pair: the UKS CH4 wrong-convergence fixture (the mkl
# exported from the propagated solution and the input that produced it)
cp "$REPO/fixtures/orca/ch4_diss_prop.mkl" work/
cp "$REPO/fixtures/orca/inputs/ch4_diss_prop.inp" work/
# the menu-31 property-file fixture (N2 SA-CASSCF, 3 roots)
cp "$REPO/fixtures/orca/n2_sa.property.txt" work/
# the pysisyphus pair of menus 32/33: a structure for the input generator, and a
# real run laid out flat (capture, trajectory, closing geometry, record, h5 and
# the per-call ORCA outputs) so that menu 33 can read the directory directly
cp "$REPO/fixtures/pysisyphus/h2o_opt/h2o.xyz" work/h2o_start.xyz
cp "$REPO/fixtures/pysisyphus/h2o_opt/h2o.xyz" work/h2o_ts.xyz
for name in run_stdout.log optimization.trj final_geometry.xyz RUN.yaml \
            optimization.h5 h2o_opt.yaml; do
  cp "$REPO/fixtures/pysisyphus/h2o_opt/$name" work/$name
done
cp -r "$REPO/fixtures/pysisyphus/h2o_opt/qm_calcs" work/qm_calcs
# the menu-34 Judd-Ofelt dataset (the Eu3+ Babu-2000 regression fixture)
cp "$REPO/fixtures/judd_ofelt/babu2000_eu3.yaml" work/
# the menu-35 pNMR pair (the CO+ QDPT run: structure, output, run file)
cp "$REPO/fixtures/pnmr/co_plus/co_plus.xyz" work/
cp "$REPO/fixtures/pnmr/co_plus/co_plus_qdpt.out" work/
cp "$REPO/fixtures/pnmr/co_plus/co_plus_qdpt_g.yaml" work/
# the menu-36 relaxation pair (the CO+ Kramers-doublet SINGLE_ANISO output and
# its single-mode Orca_Magrelax run)
cp "$REPO/fixtures/single_aniso/co_aniso2.out" work/
cp "$REPO/fixtures/magrelax/co_magrelax.out" work/
# the menu-37 core-excited-spectra fixture (the [FeCl4]2- ROCIS run)
cp "$REPO/fixtures/rocis/fecl4_xas.out" work/
# the menu-38 AILFT fixture (the Ni(2+) d8 free-ion run)
cp "$REPO/fixtures/ailft/ni_ailft.out" work/
# the menu-39 polynuclear-magnetism fixture (the two-center POLY_ANISO probe)
cp "$REPO/fixtures/poly_aniso/two_center_probe.out" work/
# the menu-40 hyperfine fixtures (the DFT and CASSCF EPRNMR probes) and the
# menu-41 magnetocaloric fixture (the POLY_ANISO M(H) probe)
cp "$REPO/fixtures/hyperfine/cef3_epr_dft.out" work/
cp "$REPO/fixtures/hyperfine/cef3_epr_casscf.out" work/
cp "$REPO/fixtures/magnetocaloric/poly_mh.out" work/
# the menu-1 Gaussian leg (the G09 TS-optimization probe)
cp "$REPO/fixtures/gaussian/g09_h2co_ts.out" work/
# the menu-42..45 external-program fixtures: the xTB captures and the water
# geometry (menus 42/44), the CREST ensemble (menu 43) and the automr run of
# the generator's own input (menu 45)
cp "$REPO/fixtures/xtb/h2o_ohess.out" work/
cp "$REPO/fixtures/xtb/nh3_planar_hess.out" work/
cp "$REPO/fixtures/xtb/h2o.xyz" work/h2o_probe.xyz
mkdir -p work/crest_run
cp "$REPO/fixtures/crest/crest_conformers.xyz" work/crest_run/
cp "$REPO/fixtures/crest/crest.energies" work/crest_run/
cp "$REPO/fixtures/crest/crest_best.xyz" work/crest_run/
cp "$REPO/fixtures/mokit/h2o_generated_automr.out" work/
# the menu-46 tunnelling fixture: the constructed neighbour-table probe
cp "$REPO/fixtures/qtm/neighbours_example.txt" work/qtm_neighbours.txt

for script in scripts/*.txt; do
  name=$(basename "$script" .txt)
  "$PY" -m fblockkit run "$script" > "expected/${name}.txt" 2>&1
  echo "captured expected/${name}.txt"
done

# the menu-8 example saves a script; replaying that script is part of the example
"$PY" -m fblockkit run work/saved_session.txt > expected/m8_replay.txt 2>&1 || true
echo "captured expected/m8_replay.txt"

# the script saved by the menu-8 example (its content is quoted in the manual)
cp work/saved_session.txt expected/saved_session.txt

# products written next to the inputs (reports, generated inputs, corrected inputs,
# the gbw-ready mkl of menu 18, and the menu-19 batch directory)
cp work/*.fbk.md expected/products/
cp work/*.fbk.json expected/products/
cp work/*.fbk.inp expected/products/
cp work/*.fbk.mkl expected/products/
cp work/*.fix_*.inp expected/products/
cp work/*.pysisyphus.xyz work/*.pysisyphus.yaml expected/products/
cp -r work/*.fbk.dm expected/products/
# the menu-44 generated input and the menu-43 ensemble report (the latter lives
# inside its run directory; the products copy keeps just the report)
cp work/*_automr.gjf expected/products/
mkdir -p expected/products/crest_run
cp work/crest_run/crest_ensemble.fbk.md expected/products/crest_run/

echo "done; products in expected/products/"
