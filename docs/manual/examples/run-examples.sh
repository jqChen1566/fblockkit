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
PY=python

rm -rf work
mkdir -p work expected expected/products

"$PY" generate_inputs.py

cp "$REPO/fixtures/orca/n2_casscf_nevpt2.out" work/
cp "$REPO/fixtures/orca/generated_ce3_sarc2.out" work/
cp "$REPO/fixtures/orca/n2_stretch_local_spin.out" work/
cp "$REPO/fixtures/orca/fhh_optts_freq.out" work/
cp "$REPO/fixtures/orca/n2_hf_clean.out" work/
cp "$REPO/fixtures/orca/scf_noconv.out" work/
cp "$REPO/fixtures/orca/inputs/scf_noconv.inp" work/
cp "$REPO/fixtures/literature/pucl3_s18.json" work/
# the exact-entropy chain of menu 12 (converged CASSCF output, its FCIDUMP and
# the two orca_2json exports of the gbw)
cp "$REPO/fixtures/orca/n2_fcidump_step_a.out" work/
cp "$REPO/fixtures/orca/n2_fcidump.fcidump" work/FCIDUMP
cp "$REPO/fixtures/orca/n2_fcidump.canonical.json" work/canonical.json
cp "$REPO/fixtures/orca/n2_fcidump.localized.json" work/localized.json

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

# products written next to the inputs (reports, generated inputs, corrected inputs)
cp work/*.fbk.md expected/products/
cp work/*.fbk.json expected/products/
cp work/*.fbk.inp expected/products/
cp work/*.fix_*.inp expected/products/

echo "done; products in expected/products/"
