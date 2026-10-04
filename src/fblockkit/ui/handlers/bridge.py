"""Menu 50: the optional execution bridge -- run a generated input with ORCA.

The bridge is deliberately outside the toolkit's core: only the UI layer imports
it (the layered contract in ``.importlinter`` keeps the direction core <-
bridge <- ui, so dropping the subpackage strips the bridge without touching the
generation and analysis layers), it acts only when the user explicitly asks,
and a script replay never re-runs an engine -- the run point prints one fixed
line instead.

The run's own report surface is small: the verdict, the output path and the
run.json ledger (a side product carrying the configuration snapshot and hashes;
it never enters a byte-compared report).  The generated input itself is left
untouched -- the bridge works on a copy inside ``runs/<case>/``.
"""

from __future__ import annotations

from pathlib import Path

from ..session import Session


def _ask_int(session: Session, prompt: str, default: int) -> int | None:
    text = session.ask(f"{prompt} (Enter = {default})").strip()
    if text == "":
        return default
    try:
        value = int(text)
    except ValueError:
        session.say(f"{text!r} is not a whole number. Next step: give an integer, or Enter for {default}.")
        return None
    return value


def bridge_run(session: Session) -> None:
    """Menu 50: run one generated input with the local ORCA and judge it."""
    from ...bridge import REPLAY_SKIP_TEXT, BridgeConfig, BridgeError, run_bridge

    input_text = session.ask("Input file path (an .inp this toolkit generated)")
    if not input_text:
        session.say("Cancelled (no input given).")
        return
    input_path = Path(input_text)
    if not input_path.is_file():
        session.say(f"input file does not exist: {input_path}")
        return
    orca_text = (
        session.ask(
            "Full pathname of the ORCA executable (the manual requires the complete path)"
        )
        or ""
    ).strip()
    if not orca_text:
        session.say("Cancelled (no ORCA path given).")
        return
    nprocs = _ask_int(session, "Parallel processes", 1)
    maxcore = _ask_int(session, "MaxCore in MB (per core)", 4000)
    timeout = _ask_int(session, "Wall-clock timeout in seconds", 3600)
    if nprocs is None or maxcore is None or timeout is None:
        session.say("Cancelled (a non-integer was given).")
        return
    if session.is_replay:
        # A replayed script never re-runs an engine (and does not probe the
        # machine for one): the fixed line replaces the run.
        session.say(REPLAY_SKIP_TEXT)
        return
    try:
        config = BridgeConfig(
            orca_path=orca_text,
            workdir_root=Path("runs"),
            nprocs=nprocs,
            maxcore_mb=maxcore,
            timeout_s=timeout,
        )
    except BridgeError as exc:
        session.say(f"bridge configuration rejected: {exc}")
        return
    try:
        result = run_bridge(config, input_path=input_path, is_replay=False)
    except BridgeError as exc:
        session.say(f"the run failed before launch: {exc}")
        return
    session.say(f"Verdict: {result.verdict}")
    if result.returncode is not None:
        session.say(f"Return code: {result.returncode}")
    if result.output_path is not None:
        session.say(f"Output written: {result.output_path}")
    if result.run_json_path is not None:
        session.say(f"Ledger written: {result.run_json_path}")
    session.say(
        "Analyse the output with menu 1 (the check-up report); the ledger is a "
        "side product and stays out of any byte-compared report. Ions the tool "
        "does not cover are refused with a measurement, never guessed."
    )
