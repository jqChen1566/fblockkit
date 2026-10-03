"""The local execution bridge (E1, A-tier): prepare, launch, monitor, judge.

The bridge is optional and off by default: the core of the toolkit neither
imports it nor depends on it (the layered import contract keeps the direction
core <- bridge <- ui; see .importlinter), and this module is stdlib only
(``subprocess``/``pathlib``/``json``/``hashlib``/``time``).

The judgement follows the measured engine semantics (ORCA 6.1.1 on the probe
machine, 2026-10-03, E0 probe record):

- the exit code carries exactly ONE reliable negative signal: rc = 11 is the
  input scanner's rejection (``Unknown identifier`` in a ``%`` block); every
  other failure measured -- an SCF non-convergence (``ORCA finished by error
  termination in LEANSCF``), a missing gbw, a failing orca_mapspc -- returned
  rc = 0;
- therefore the return code never decides alone: ``normal`` is given only
  when the text carries ``****ORCA TERMINATED NORMALLY****``, and a silent
  rc = 0 exit is ``suspicious``, never success.

The judgement table (three sources: return code, output text, heartbeat):

- rc 0 + the termination banner -> ``normal``;
- rc 0 + ``ORCA finished by error termination in <module>`` -> ``orca_error``;
- rc 11 -> ``input_rejected`` (the input-scanner stage);
- the text ``Cannot open GBW file`` -> ``input_error`` (a missing
  prerequisite, measured at rc 0);
- a signal-style return code (negative, or 128 + n as shells report it) with
  no banner -> ``killed``;
- no banner + the bridge itself terminated the process (wall clock or a
  heartbeat stall) -> ``timeout``;
- anything else that ends without a banner -> ``suspicious``.

Every run follows prepare -> launch -> monitor -> classify -> collect:
``<case>/<case>.inp`` is a fresh copy (injected with the missing ``%pal`` /
``%maxcore`` / ``%base``; the generated original is never edited), the engine
is started by full pathname without a shell and without mpirun, the monitor
polls the ``.tmp``/``.gbw`` mtimes and the ``.out`` size as the heartbeat, and
``<case>/run.json`` records the configuration snapshot, the times, the
verdict, the sha256 sums and the injection record.

Determinism rules: a script replay skips the run with one fixed line
(:data:`REPLAY_SKIP_TEXT`) and starts no process; ``run.json`` is the only
artifact carrying times and is a side product -- it never enters the
byte-comparison surface.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import time
from collections.abc import Callable
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from .config import BridgeConfig, BridgeError, config_to_dict

__all__ = [
    "BridgeResult",
    "REPLAY_SKIP_TEXT",
    "TERMINATED_MARK",
    "VERDICTS",
    "classify",
    "run_bridge",
]

#: the termination banner of a successful run (the parser layer carries the
#: same marker; the bridge keeps its own copy so the package stays stdlib-only)
TERMINATED_MARK = "****ORCA TERMINATED NORMALLY****"

#: the banner of an engine-side error termination (a module name follows)
ERROR_TERMINATION_MARK = "ORCA finished by error termination in"

#: the measured sign of the input scanner: the one reliable non-zero exit code
INPUT_REJECTED_RC = 11

#: a missing prerequisite, reported in the text while rc stays 0
MISSING_FILE_MARK = "Cannot open GBW file"

#: the fixed line a replayed run point prints (a replay never re-runs)
REPLAY_SKIP_TEXT = "(skipped on replay: the execution bridge never re-runs)"

#: after SIGTERM the bridge waits this grace window before SIGKILL
TERMINATION_GRACE_S = 30.0

#: the judgement tokens of the three-source table
VERDICTS = (
    "normal",
    "orca_error",
    "input_rejected",
    "input_error",
    "killed",
    "timeout",
    "suspicious",
)

#: the replay marker (a run state, outside the judgement table)
SKIPPED = "skipped"


@dataclass(frozen=True)
class BridgeResult:
    """The outcome of one bridged run (the ``run.json`` ledger as an object).

    ``verdict`` is one of :data:`VERDICTS`, or ``"skipped"`` for a replayed
    run point (where every other field stays empty).  The time fields are the
    run.json ledger's own content -- the bridge prints nothing with
    timestamps, and the ledger is a side product.
    """

    case: str
    verdict: str
    returncode: int | None = None
    case_dir: Path | None = None
    input_path: Path | None = None
    output_path: Path | None = None
    run_json_path: Path | None = None
    injections: tuple[str, ...] = ()
    injections_skipped: tuple[str, ...] = ()
    input_sha256: str | None = None
    output_sha256: str | None = None
    started_utc: str | None = None
    finished_utc: str | None = None
    terminated_by: str | None = None


@dataclass(frozen=True)
class _Prepared:
    """The case layout prepare() wrote (input copy included)."""

    case_dir: Path
    input_path: Path
    output_path: Path
    injections: tuple[str, ...]
    injections_skipped: tuple[str, ...]
    fresh_dir: bool


@dataclass(frozen=True)
class _Monitored:
    """What the monitor observed: the reaped exit status and the kill cause."""

    returncode: int
    terminated_by: str | None
    wall_s: float


# --- the steps of run_bridge -------------------------------------------------


def _utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _check_case(case) -> None:
    if not isinstance(case, str) or not case:
        raise BridgeError(
            f"the case name {case!r} is empty. Next step: give a plain name "
            "(it becomes <run root>/<case>/<case>.inp)."
        )
    bad = any(char in case for char in '\\/:*?"<>|') or any(
        char < " " for char in case
    )
    if case in (".", "..") or bad or case != case.strip() or case.endswith("."):
        raise BridgeError(
            f"the case name {case!r} cannot be used as a directory name. Next "
            "step: give a plain name without path separators or reserved "
            "characters."
        )


def _injection_point(lines: list[str]) -> int:
    """The index the injected blocks go to: after the leading run of
    simple-input lines (``!``), or the top of the file when there is none."""
    index = 0
    for position, line in enumerate(lines):
        if line.strip().startswith("!"):
            index = position + 1
        elif index:
            break
    return index


def _prepare(config: BridgeConfig, *, case: str, input_text: str) -> _Prepared:
    """Create ``<root>/<case>/`` and write the injected copy of the input.

    The blocks the input already carries are kept verbatim and recorded in
    ``injections_skipped``; the missing ones are injected in the shapes of the
    recipe writers (``%pal nprocs N end`` / ``%maxcore M`` / ``%base "name"``).
    """
    case_dir = config.workdir_root / case
    fresh = not case_dir.exists()
    try:
        case_dir.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise BridgeError(
            f"the case directory {case_dir} cannot be created ({exc}). Next "
            "step: check the run root of the configuration."
        ) from exc
    input_path = case_dir / f"{case}.inp"
    if input_path.exists():
        raise BridgeError(
            f"{input_path} already exists: a previous run left its input "
            "behind. Next step: choose another case name or clear that "
            "directory (the engine manual asks for a fresh directory per "
            "run, so stale scratch files cannot be mistaken for this run's)."
        )
    lines = input_text.splitlines()
    lowered = [line.strip().lower() for line in lines]
    injections_skipped = [
        block
        for block in ("%pal", "%maxcore", "%base")
        if any(line.startswith(block) for line in lowered)
    ]
    injections = []
    if "%pal" not in injections_skipped:
        injections.append(f"%pal nprocs {config.nprocs} end")
    if "%maxcore" not in injections_skipped:
        injections.append(f"%maxcore {config.maxcore_mb}")
    if config.base_name is not None and "%base" not in injections_skipped:
        injections.append(f'%base "{config.base_name}"')
    if injections:
        at = _injection_point(lines)
        lines = lines[:at] + injections + lines[at:]
    text = "\n".join(lines)
    if not text.endswith("\n"):
        text += "\n"
    input_path.write_text(text, encoding="utf-8")
    return _Prepared(
        case_dir=case_dir,
        input_path=input_path,
        output_path=case_dir / f"{case}.out",
        injections=tuple(injections),
        injections_skipped=tuple(injections_skipped),
        fresh_dir=fresh,
    )


def _launch(config: BridgeConfig, prepared: _Prepared, handle) -> subprocess.Popen:
    """Start the engine by full pathname, no shell, no mpirun (the manual's
    calling rules); stdin is diverted so a menu session's input is never
    stolen by the engine."""
    try:
        return subprocess.Popen(
            [config.orca_path, prepared.input_path.name],
            cwd=prepared.case_dir,
            stdin=subprocess.DEVNULL,
            stdout=handle,
            stderr=subprocess.STDOUT,
        )
    except OSError as exc:
        raise BridgeError(
            f"the engine could not be started ({exc}). Next step: check "
            "orca_path -- it must be the executable itself, given by its "
            "full pathname (the bridge never goes through a shell and never "
            "adds mpirun)."
        ) from exc


def _heartbeat(case_dir: Path, output_path: Path) -> tuple[int, int]:
    """The heartbeat signature: the newest ``.tmp``/``.gbw`` mtime plus the
    ``.out`` size (a change means the engine is still working)."""
    newest = 0
    for pattern in ("*.tmp", "*.gbw"):
        for path in case_dir.glob(pattern):
            try:
                newest = max(newest, path.stat().st_mtime_ns)
            except OSError:
                continue
    try:
        size = output_path.stat().st_size
    except OSError:
        size = -1
    return newest, size


def _terminate(proc: subprocess.Popen) -> None:
    """SIGTERM, wait the grace window, then SIGKILL."""
    try:
        proc.terminate()
    except OSError:
        pass
    deadline = time.monotonic() + TERMINATION_GRACE_S
    while proc.poll() is None and time.monotonic() < deadline:
        time.sleep(0.2)
    if proc.poll() is None:
        try:
            proc.kill()
        except OSError:
            pass
        proc.wait()


def _monitor(
    config: BridgeConfig, prepared: _Prepared, proc: subprocess.Popen
) -> _Monitored:
    """Wait for the process while polling the wall clock and the heartbeat at
    the configured interval; terminate on a limit (the wall clock first, then
    a heartbeat stall).

    The wait returns the moment the process exits, so a run shorter than the
    polling interval is not delayed by it; the wait step is the smaller of the
    polling interval and the time left to the nearer deadline, so a limit
    fires at the deadline rather than up to one interval late.
    """
    started = time.monotonic()
    last_activity = started
    signature = _heartbeat(prepared.case_dir, prepared.output_path)
    terminated_by = None
    returncode = proc.poll()
    while returncode is None:
        now = time.monotonic()
        until_wall = config.timeout_s - (now - started)
        until_stall = config.stall_timeout_s - (now - last_activity)
        wait_s = max(min(config.poll_interval_s, until_wall, until_stall), 0.05)
        try:
            returncode = proc.wait(timeout=wait_s)
            break
        except subprocess.TimeoutExpired:
            pass
        now = time.monotonic()
        if now - started > config.timeout_s:
            terminated_by = "wall"
            _terminate(proc)
            returncode = proc.wait()
            break
        current = _heartbeat(prepared.case_dir, prepared.output_path)
        if current != signature:
            signature = current
            last_activity = now
        elif now - last_activity > config.stall_timeout_s:
            terminated_by = "stall"
            _terminate(proc)
            returncode = proc.wait()
            break
    return _Monitored(
        returncode=int(returncode),
        terminated_by=terminated_by,
        wall_s=time.monotonic() - started,
    )


def classify(
    *, text: str, returncode: int | None, terminated_by: str | None = None
) -> str:
    """Judge one finished run from the three sources (see the module table).

    ``normal`` is given only when the termination banner is present -- never
    from the return code alone; ``terminated_by`` ("wall"/"stall") marks a
    run the bridge itself terminated.
    """
    if TERMINATED_MARK in text:
        return "normal"
    if terminated_by is not None:
        return "timeout"
    if returncode == INPUT_REJECTED_RC:
        return "input_rejected"
    if MISSING_FILE_MARK in text:
        return "input_error"
    if ERROR_TERMINATION_MARK in text:
        return "orca_error"
    if returncode is not None and (returncode < 0 or returncode >= 128):
        return "killed"
    return "suspicious"


def _cleanup_failed_launch(prepared: _Prepared) -> None:
    """Remove the artifacts of a run that never started (only when this call
    created the case directory; a pre-existing directory is left alone)."""
    if not prepared.fresh_dir:
        return
    for path in (prepared.input_path, prepared.output_path):
        try:
            path.unlink()
        except OSError:
            continue
    try:
        prepared.case_dir.rmdir()
    except OSError:
        pass


def _collect(
    config: BridgeConfig,
    prepared: _Prepared,
    monitored: _Monitored,
    verdict: str,
    *,
    started_utc: str,
    finished_utc: str,
    raw_output: bytes,
) -> BridgeResult:
    """Write the ``run.json`` ledger and apply the scratch policy."""
    raw_input = prepared.input_path.read_bytes()
    result = BridgeResult(
        case=prepared.case_dir.name,
        verdict=verdict,
        returncode=monitored.returncode,
        case_dir=prepared.case_dir,
        input_path=prepared.input_path,
        output_path=prepared.output_path,
        run_json_path=prepared.case_dir / "run.json",
        injections=prepared.injections,
        injections_skipped=prepared.injections_skipped,
        input_sha256=hashlib.sha256(raw_input).hexdigest(),
        output_sha256=hashlib.sha256(raw_output).hexdigest(),
        started_utc=started_utc,
        finished_utc=finished_utc,
        terminated_by=monitored.terminated_by,
    )
    record = {
        "case": result.case,
        "config": config_to_dict(config),
        "started_utc": started_utc,
        "finished_utc": finished_utc,
        "wall_s": round(monitored.wall_s, 3),
        "returncode": monitored.returncode,
        "terminated_by": monitored.terminated_by,
        "verdict": verdict,
        "injections": list(prepared.injections),
        "injections_skipped": list(prepared.injections_skipped),
        "input_file": prepared.input_path.name,
        "input_sha256": result.input_sha256,
        "output_file": prepared.output_path.name,
        "output_sha256": result.output_sha256,
    }
    result.run_json_path.write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if not config.keep_scratch:
        for pattern in ("*.tmp", "*.gbw"):
            for path in prepared.case_dir.glob(pattern):
                try:
                    path.unlink()
                except OSError:
                    continue
    return result


def _source_text(input_text, input_path) -> str:
    """Exactly one of the two input forms; the file is read, never edited."""
    if (input_text is None) == (input_path is None):
        raise BridgeError(
            "give exactly one input: the generated text (input_text) or a "
            "file (input_path). Next step: pass one of the two."
        )
    if input_path is not None:
        source = Path(input_path)
        if not source.is_file():
            raise BridgeError(
                f"the input file {str(source)!r} does not exist. Next step: "
                "check the path -- the bridge copies the file into "
                "<case>/<case>.inp and edits only the copy."
            )
        try:
            return source.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError) as exc:
            raise BridgeError(
                f"the input file {str(source)!r} cannot be read ({exc}). "
                "Next step: give a UTF-8 engine input file."
            ) from exc
    if not isinstance(input_text, str):
        raise BridgeError(
            f"input_text is a {type(input_text).__name__}, not text. Next "
            "step: give the input as a string, or use input_path."
        )
    return input_text


def run_bridge(
    config: BridgeConfig,
    *,
    input_text: str | None = None,
    input_path: str | Path | None = None,
    case: str | None = None,
    is_replay: bool = False,
    say: Callable[[str], None] | None = None,
) -> BridgeResult:
    """Run one generated input through the local engine and judge it.

    ``config``: a :class:`BridgeConfig`; exactly one of ``input_text`` (the
    text just generated) or ``input_path`` (a file to copy); ``case``: the
    case name -- defaults to the input file's stem (give it explicitly for
    ``input_text``); ``is_replay``: a script replay skips the run entirely
    (no process, no files) and sends the fixed line
    :data:`REPLAY_SKIP_TEXT` to ``say``; ``say``: an optional line sink (the
    session's ``say``) -- the bridge emits nothing else through it, so no
    progress output with times can enter a replay-comparable stream.

    Returns a :class:`BridgeResult` (verdict, return code, paths, the
    injection record, the artifact sha256 sums); ``<case>/run.json`` carries
    the same facts plus the configuration snapshot and the times.
    """
    if not isinstance(config, BridgeConfig):
        raise BridgeError(
            f"the configuration is a {type(config).__name__}, not a "
            "BridgeConfig. Next step: build one (see fblockkit.bridge)."
        )
    source_text = _source_text(input_text, input_path)
    if case is None:
        if input_path is None:
            raise BridgeError(
                "no case name was given for a text input. Next step: pass a "
                "case name (it becomes <run root>/<case>/<case>.inp)."
            )
        case = Path(input_path).stem
    _check_case(case)
    if is_replay:
        if say is not None:
            say(REPLAY_SKIP_TEXT)
        return BridgeResult(case=case, verdict=SKIPPED)
    prepared = _prepare(config, case=case, input_text=source_text)
    started_utc = _utc_now()
    with prepared.output_path.open("wb") as handle:
        try:
            proc = _launch(config, prepared, handle)
        except BridgeError:
            _cleanup_failed_launch(prepared)
            raise
        try:
            monitored = _monitor(config, prepared, proc)
        except BaseException:
            # never leave an engine running on an unexpected exit (a keyboard
            # interrupt included): terminate the process before re-raising
            _terminate(proc)
            raise
    finished_utc = _utc_now()
    raw_output = prepared.output_path.read_bytes()
    verdict = classify(
        text=raw_output.decode("utf-8", errors="replace"),
        returncode=monitored.returncode,
        terminated_by=monitored.terminated_by,
    )
    return _collect(
        config,
        prepared,
        monitored,
        verdict,
        started_utc=started_utc,
        finished_utc=finished_utc,
        raw_output=raw_output,
    )
