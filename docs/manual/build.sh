#!/usr/bin/env bash
# Build the fBlockKit manual.
# Pass sequence (project rule for documents with a bibliography + index):
#   xelatex -> biber -> makeindex -> xelatex -> xelatex
# The bibliography is the program's own sources.bib (single source of truth):
# it is copied into generated/ on every build.
set -euo pipefail
cd "$(dirname "$0")"

mkdir -p generated
cp ../../src/fblockkit/knowledge/sources.bib generated/sources.bib

step() { printf '\n== %s ==\n' "$*"; "$@"; }

step xelatex -interaction=nonstopmode -halt-on-error main.tex
step biber main
if [ -f main.idx ]; then
  step makeindex -q main.idx
fi
step xelatex -interaction=nonstopmode -halt-on-error main.tex
step xelatex -interaction=nonstopmode -halt-on-error main.tex

printf '\n== acceptance summary ==\n'
fail=0
if grep -E "LaTeX Warning: (There were undefined references|Label\(s\) may have changed|Citation .* undefined|Reference .* undefined)" main.log; then
  echo "!! unresolved references or a pending rerun" >&2
  fail=1
else
  echo "references: resolved (no undefined/citation/rerun warnings)"
fi
# Overfull boxes: the acceptance line is "none above 20pt"; the count is printed
# for both bands so regressions are visible.
python - main.log <<'PYEOF'
import re, sys
log = open(sys.argv[1], encoding="utf-8", errors="replace").read()
vals = [float(m) for m in re.findall(r"Overfull \\hbox \(([0-9.]+)pt", log)]
big = [v for v in vals if v > 20.0]
print(f"overfull boxes: {len(vals)} total, {len(big)} above 20pt (max {max(vals, default=0):.1f}pt)")
sys.exit(1 if big else 0)
PYEOF
rc=$?
[ $rc -ne 0 ] && fail=1
ls -l main.pdf
exit $fail
