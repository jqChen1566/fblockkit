"""Interactive session engine: interaction and script replay share one code path.

"Interaction is a script" (the Multiwfn pattern):

- all input comes from a **sequence of lines** -- standard input in interactive mode,
  a script file in replay mode;
- every prompt consumes exactly one line; **an empty line means pressing Enter (take
  the default)**; a line starting with ``#`` is a comment and is removed during
  preprocessing (it does not count towards the line count); input lines are echoed to
  the output so the two modes produce byte-identical output;
- the session records the lines it consumes in ``record`` -- "record it once and you
  have a script", which menu 8 can write to disk.

EOF (the script is exhausted or the input stream ends) means exit at the menu; at a
prompt it returns None, which the handler treats as "cancelled", so an incomplete
script behaves predictably.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Callable, Iterable, TextIO

from .menu import MenuItem, find_item, load_menu, render_menu

PROMPT_SUFFIX = "> "


def strip_comments(lines: Iterable[str]) -> list[str]:
    """Remove comment lines (starting with #); blank lines are kept -- a blank line is
    itself the "press Enter" input."""
    kept = []
    for line in lines:
        if line.strip().startswith("#"):
            continue
        kept.append(line.rstrip("\n"))
    return kept


class Session:
    """One session: menu loop + prompt reading + recording."""

    def __init__(
        self,
        *,
        lines: Iterable[str],
        out: TextIO | None = None,
        root: Path | None = None,
        menu: tuple[MenuItem, ...] | None = None,
    ) -> None:
        self.out = out if out is not None else sys.stdout
        self._source = iter(strip_comments(lines))
        self.root = Path(root) if root is not None else Path(".")
        self.menu = menu if menu is not None else load_menu()
        self.record: list[str] = []
        self.eof = False
        self.stop = False  # a handler sets this to True to end the main loop (e.g. "quit")

    # --- input and output ---------------------------------------------------

    def say(self, text: str = "") -> None:
        self.out.write(text + "\n")

    def ask(self, prompt: str, default: str = "") -> str | None:
        """Print a prompt and read one line. A blank line -> the default; EOF -> None
        (the handler treats it as cancelled)."""
        hint = f"(Enter = {default})" if default else "(Enter = leave empty)"
        self.out.write(f"{prompt} {hint}{PROMPT_SUFFIX}")
        try:
            line = next(self._source)
        except StopIteration:
            self.eof = True
            self.out.write("\n")
            return None
        self.record.append(line)
        self.out.write(line + "\n")  # echo: interactive and replay output have the same shape
        value = line.strip()
        return default if value == "" else value

    def snapshot(self) -> list[str]:
        """The lines recorded so far (for "save this session as a script"; it excludes
        the saving action itself)."""
        return list(self.record)

    def save_record(self, path: Path | str) -> Path:
        target = Path(path)
        target.write_text("\n".join(self.record) + "\n", encoding="utf-8")
        return target

    # --- main loop ----------------------------------------------------------

    def run(self, handlers: dict[str, Callable[["Session"], None]]) -> int:
        self.say("fBlockKit -- the f-block calculation toolkit (input generation + characterisation analysis; it runs no calculations)")
        while True:
            self.say("")
            self.say(render_menu(self.menu))
            choice = self.ask("Choose a number")
            if choice is None:
                self.say("(end of input, exiting.)")
                return 0
            if choice == "":
                continue
            item = find_item(self.menu, choice)
            if item is None:
                self.say(f"Invalid number: {choice} (nothing was done)")
                continue
            handler = handlers.get(item.handler)
            if handler is None:
                self.say(
                    f"menu item {choice} has no implementation for its handler "
                    f"{item.handler!r} -- this is a program defect, please report it."
                )
                return 2
            handler(self)
            if self.stop:
                return 0
