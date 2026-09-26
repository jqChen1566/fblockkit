"""``python -m fblockkit`` entry point (equivalent to the installed ``fblockkit``
command).

It lives at package level rather than in ``ui.cli`` because
``python -m fblockkit.ui.cli`` triggers runpy's double-import warning (that module
is already imported by the package ``__init__``); the package-level entry avoids it.
"""

from .ui.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
