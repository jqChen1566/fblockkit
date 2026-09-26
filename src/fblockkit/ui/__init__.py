"""User-interface layer: declarative menu, CLI, script replay
(architecture design v0.1 §3 L3).

- ``menu.yaml``: the single source of number -> title -> handler (numbers are
  append-only, never reordered);
- ``session.Session``: interaction and script replay share one code path
  (interaction is a script);
- ``handlers``: the handler functions behind the menu items;
- ``cli.main``: the command-line entry point (interactive / run replay / search / guide).
"""

from .cli import main, run_session
from .menu import MenuError, MenuItem, find_item, load_menu, render_menu
from .session import Session, strip_comments

__all__ = [
    "MenuError",
    "MenuItem",
    "Session",
    "find_item",
    "load_menu",
    "main",
    "render_menu",
    "run_session",
    "strip_comments",
]
