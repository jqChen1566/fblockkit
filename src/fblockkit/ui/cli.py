"""fBlockKit command-line entry point.

Usage:

- ``fblockkit``                 interactive menu (input is a script; ``--record FILE``
                                records this session's input)
- ``fblockkit run script.txt``  replay a script (identical to typing the same lines by hand)
- ``fblockkit submit input.inp``  run an input with ORCA (the optional execution bridge;
                                ``--orca PATH`` full pathname, ``--nprocs N``,
                                ``--maxcore MB``, ``--timeout S``)
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

USAGE = """fBlockKit -- the f-block calculation toolkit (input generation + characterisation analysis; the core runs no calculations -- an optional execution bridge, off by default, can run a generated input when explicitly asked)

Usage:
  fblockkit                      interactive menu (type at the prompts; Enter = default)
  fblockkit --record FILE        the same, and record this session's input as a replayable script
  fblockkit run script.txt       replay a script (format: one input per line, # for comments)
  fblockkit submit input.inp     run an input with ORCA (--orca PATH, --nprocs N,
                                 --maxcore MB, --timeout S; the optional execution bridge)
  fblockkit search KEYWORD       search the tool index
  fblockkit guide TOOL_ID        tool onboarding notes
"""


def run_session(lines, out: TextIO, record_path: Path | None = None, is_replay: bool = False) -> int:
    session = Session(lines=lines, out=out, is_replay=is_replay)
    code = session.run(HANDLERS)
    if record_path is not None:
        session.save_record(record_path)
        out.write(f"(this session's input was recorded to: {record_path})\n")
    return code


def _run_submit(args: list[str]) -> int:
    """The execution bridge's CLI door: run one input with ORCA and judge it.

    Kept deliberately thin: parse the flags, build the explicit configuration,
    hand off to the bridge and print the verdict -- the bridge itself owns the
    process handling and the three-source classification.
    """
    from ..bridge import BridgeConfig, BridgeError, run_bridge

    input_path: Path | None = None
    orca_path: str | None = None
    nprocs, maxcore, timeout = 1, 4000, 3600
    i = 0
    while i < len(args):
        token = args[i]
        if token == "--orca" and i + 1 < len(args):
            orca_path = args[i + 1]
            i += 2
            continue
        if token == "--nprocs" and i + 1 < len(args):
            nprocs = int(args[i + 1])
            i += 2
            continue
        if token == "--maxcore" and i + 1 < len(args):
            maxcore = int(args[i + 1])
            i += 2
            continue
        if token == "--timeout" and i + 1 < len(args):
            timeout = int(args[i + 1])
            i += 2
            continue
        if input_path is None and not token.startswith("--"):
            input_path = Path(token)
            i += 1
            continue
        sys.stdout.write(USAGE)
        return 2
    if input_path is None:
        sys.stdout.write("submit needs an input path.\n")
        return 2
    if orca_path is None:
        sys.stdout.write("submit needs --orca with the full pathname of the ORCA executable.\n")
        return 2
    try:
        config = BridgeConfig(
            orca_path=orca_path,
            workdir_root=Path("runs"),
            nprocs=nprocs,
            maxcore_mb=maxcore,
            timeout_s=timeout,
        )
        result = run_bridge(config, input_path=input_path, is_replay=False)
    except BridgeError as exc:
        sys.stdout.write(f"{exc}\n")
        return 2
    sys.stdout.write(f"verdict: {result.verdict}\n")
    if result.output_path is not None:
        sys.stdout.write(f"output: {result.output_path}\n")
    if result.run_json_path is not None:
        sys.stdout.write(f"ledger: {result.run_json_path} (a side product; not part of any byte-compared report)\n")
    return 0


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
        return run_session(
            script.read_text(encoding="utf-8").splitlines(), sys.stdout, is_replay=True
        )
    if head == "submit":
        return _run_submit(rest)
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
