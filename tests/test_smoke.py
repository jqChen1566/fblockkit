"""Step-1 smoke tests: the package imports and the CLI runs.

The layering contract is checked separately by `lint-imports` (see scripts/verify.sh).
"""

from __future__ import annotations


def test_version_importable():
    import fblockkit

    assert fblockkit.__version__


def test_cli_main_runs(tmp_path, capsys):
    """CLI menu entry: replay a script that quits immediately (the main()
    argument-parsing path of step 7)."""
    from fblockkit.ui.cli import main

    script = tmp_path / "quit.txt"
    script.write_text("0\n", encoding="utf-8")
    assert main(["run", str(script)]) == 0
    out = capsys.readouterr().out
    assert "fBlockKit" in out
