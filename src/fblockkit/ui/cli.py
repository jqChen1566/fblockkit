"""fBlockKit command-line entry point.

Usage:

- ``fblockkit``                 interactive menu (input is a script; ``--record FILE``
                                records this session's input)
- ``fblockkit run script.txt``  replay a script (identical to typing the same lines by hand)
- ``fblockkit search KEYWORD``  search the tool index (T1)
- ``fblockkit guide TOOL_ID``   tool onboarding notes (T2)

Everything is written to standard output and contains no unstable information such as
timestamps -- replaying the same script twice is byte-identical.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Sequence, TextIO

from ..toolindex import ToolIndexError
from ..toolindex import guide as tool_guide_text
from ..toolindex import search as tool_search_fn
from .handlers import HANDLERS
from .session import Session

USAGE = """fBlockKit -- the f-block calculation toolkit (input generation + characterisation analysis; it runs no calculations)

Usage:
  fblockkit                      interactive menu (type at the prompts; Enter = default)
  fblockkit --record FILE        the same, and record this session's input as a replayable script
  fblockkit run script.txt       replay a script (format: one input per line, # for comments)
  fblockkit search KEYWORD       search the tool index
  fblockkit guide TOOL_ID        tool onboarding notes
"""


def run_session(lines, out: TextIO, record_path: Path | None = None) -> int:
    session = Session(lines=lines, out=out)
    code = session.run(HANDLERS)
    if record_path is not None:
        session.save_record(record_path)
        out.write(f"(this session's input was recorded to: {record_path})\n")
    return code


def main(argv: Sequence[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args:
        return run_session(sys.stdin, sys.stdout)
    head, rest = args[0], args[1:]
    if head == "--record":
        if not rest:
            sys.stdout.write("--record needs a script path.\n")
            return 2
        return run_session(sys.stdin, sys.stdout, record_path=Path(rest[0]))
    if head == "run":
        if not rest:
            sys.stdout.write("run needs a script path.\n")
            return 2
        script = Path(rest[0])
        if not script.is_file():
            sys.stdout.write(f"script does not exist: {script}\n")
            return 2
        return run_session(script.read_text(encoding="utf-8").splitlines(), sys.stdout)
    if head == "search":
        query = " ".join(rest).strip()
        if not query:
            sys.stdout.write("search needs a keyword.\n")
            return 2
        try:
            hits = tool_search_fn(query)
        except ToolIndexError as exc:
            sys.stdout.write(f"index loading failed: {exc}\n")
            return 2
        if not hits:
            sys.stdout.write(f"no hits (keywords: {query}).\n")
            return 0
        for record in hits:
            sys.stdout.write(f"- {record.line()}\n  {record.purpose}\n")
        return 0
    if head == "guide":
        if not rest:
            sys.stdout.write("guide needs a tool index id (search for one with 'search').\n")
            return 2
        try:
            sys.stdout.write(tool_guide_text(rest[0]) + "\n")
        except ToolIndexError as exc:
            sys.stdout.write(f"{exc}\n")
            return 2
        return 0
    sys.stdout.write(USAGE)
    return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
