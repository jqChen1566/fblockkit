"""PyInstaller entry script (used by the green package).

Equivalent to the installed ``fblockkit`` command and to ``python -m fblockkit``.
"""

from fblockkit.ui.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
