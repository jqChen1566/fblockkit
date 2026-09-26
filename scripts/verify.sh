#!/usr/bin/env bash
# fBlockKit unified check: the layering contract (import-linter), a byte-compile of
# every source and test file under THIS interpreter, then the tests.
set -euo pipefail
cd "$(dirname "$0")/.."
export PYTHONPATH=src
# On Windows (GBK console default) Python must be told to use UTF-8, otherwise
# lint-imports crashes while reading the sources; harmless on Linux/macOS.
export PYTHONUTF8=1
lint-imports
# The package supports Python >= 3.11 while development machines often run newer
# interpreters: compileall here is the syntax gate for the supported minimum (it
# caught a PEP 701 f-string that Python 3.12+ accepts and 3.11 rejects).
python -m compileall -q src tests
pytest -q
