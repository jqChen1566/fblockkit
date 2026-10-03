# fBlockKit release checklist (the v0.1.0 process, distilled; used for v0.2.0)

## Entry conditions (all green before the cut)

1. `scripts/verify.sh` passes on **both ends** (local and the 101 server):
   import-linter contract, the 3.11 syntax gate, the unit tests, and the
   full acceptance walk (`scripts/acceptance.py`);
2. the example captures are regenerated
   (`docs/manual/examples/run-examples.sh`) and every diff reviewed;
   `expected/products/` in sync;
3. the manual PDF rebuilt (XeLaTeX, twice) with 0 Overfull boxes above
   20 pt and all cross-references resolved; the page count recorded;
4. paper-side counts refreshed when the release accompanies a submission
   (menu count, tests, acceptance walks, manual pages);
5. `git status` clean; the public-facing history per the release discipline
   (push `release:main` only, never `main` directly).

## The cut

6. `CHANGELOG.md`: turn the `## [0.2.0] - Unreleased` heading into the
   release date;
7. `python -m build` -> wheel + sdist into `dist/`; smoke-install the wheel
   into a scratch venv and run `fblockkit --help` plus one menu smoke;
8. Windows green package: `cd packaging && python -m PyInstaller
   fblockkit.spec --noconfirm`; verify the entry runs from a clean path;
   zip the folder;
9. commit the release state and tag `v0.2.0`; push `release:main`;
10. GitHub release `v0.2.0`: attach wheel + sdist + the green-package zip;
    replace the manual Release attachment with the freshly built PDF
    (record its byte size);
11. PyPI: stays user-gated (the standing `一不管` instruction); when
    cleared, `twine upload dist/fblockkit-0.2.0*`;
12. post-release: leave the next CHANGELOG section as Unreleased; update the
    project index and memory with the release SHA and the attachment sizes.
