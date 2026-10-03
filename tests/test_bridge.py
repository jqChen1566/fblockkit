"""Checks of the optional local execution bridge (bridge/config.py,
bridge/runner.py).

No real engine is needed: every scenario runs a small inline Python snippet as
the "engine".  On Windows the snippet is reached through a .bat shim (measured:
CreateProcess cannot start a .py directly, but a .bat works); elsewhere the
snippet is a shebang script.  The snippets replay the measured engine exit
semantics (ORCA 6.1.1, E0 probe record, 2026-10-03): rc = 0 both on a normal
termination AND on an error termination, rc = 11 on an input-scanner
rejection -- so the bridge must judge from the text and the heartbeat as well.

The monitor knobs (``poll_interval_s``, ``stall_timeout_s``, ``timeout_s``)
are turned down so the timeout and stall paths are exercised in seconds.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
import textwrap
from pathlib import Path

import pytest

from fblockkit.bridge import BridgeConfig, BridgeError, run_bridge
from fblockkit.bridge import runner as bridge_runner
from fblockkit.bridge.config import read_config, write_config


def _fake_engine(directory: Path, body: str) -> Path:
    """An executable that runs the scenario body like the engine would: the
    working directory is the case directory and argv[1] is the case input."""
    directory.mkdir(parents=True, exist_ok=True)
    snippet = textwrap.dedent(body).lstrip("\n")
    if os.name == "nt":
        (directory / "fake_orca.py").write_text(snippet, encoding="utf-8")
        shim = directory / "fake_orca.bat"
        shim.write_text(
            '@echo off\r\n"' + sys.executable + '" "%~dp0fake_orca.py" %*\r\n',
            encoding="ascii",
        )
        return shim
    launcher = directory / "fake_orca"
    launcher.write_text(f"#!{sys.executable}\n" + snippet, encoding="utf-8")
    launcher.chmod(0o755)
    return launcher


def _config(tmp_path: Path, engine: Path, **overrides) -> BridgeConfig:
    values = dict(orca_path=str(engine), workdir_root=tmp_path / "runs")
    values.update(overrides)
    return BridgeConfig(**values)


# 1 --------------------------------------------------------------------------


def test_config_validation_refuses_with_a_next_step(tmp_path):
    good = dict(orca_path=sys.executable, workdir_root=tmp_path / "runs")
    with pytest.raises(BridgeError, match="Next step"):
        BridgeConfig(**{**good, "orca_path": str(tmp_path / "missing_engine")})
    for broken in (
        {"nprocs": 0},
        {"maxcore_mb": 0},
        {"timeout_s": 0},
        {"stall_timeout_s": 0},
        {"poll_interval_s": 0},
        {"base_name": 'a"b'},
    ):
        with pytest.raises(BridgeError, match="Next step"):
            BridgeConfig(**{**good, **broken})
    (tmp_path / "a_file").write_text("x", encoding="utf-8")
    with pytest.raises(BridgeError, match="Next step"):
        BridgeConfig(**{**good, "workdir_root": tmp_path / "a_file"})
    config = BridgeConfig(**good)
    assert (config.nprocs, config.maxcore_mb, config.timeout_s) == (1, 4000, 3600)
    assert config.stall_timeout_s == 3600 and config.poll_interval_s == 25.0
    assert config.keep_scratch is True and config.base_name is None
    # the YAML round trip keeps the configuration a script references
    stored = write_config(config, tmp_path / "bridge.yaml")
    assert read_config(stored) == config
    unknown = tmp_path / "unknown.yaml"
    unknown.write_text(
        "orca_path: x\nworkdir_root: y\nsurprise: 1\n", encoding="utf-8"
    )
    with pytest.raises(BridgeError, match="Next step"):
        read_config(unknown)


# 2 --------------------------------------------------------------------------


def test_injection_keeps_existing_blocks_and_the_source_untouched(tmp_path):
    engine = _fake_engine(
        tmp_path / "normal", 'print("****ORCA TERMINATED NORMALLY****")'
    )
    source = tmp_path / "generated.inp"
    source.write_text(
        "!HF def2-SVP TightSCF\n%maxcore 2000\n*xyz 0 1\nH 0 0 0\n*\n",
        encoding="utf-8",
    )
    before = source.read_bytes()
    config = _config(tmp_path, engine, nprocs=3, base_name="probe")
    result = run_bridge(config, input_path=source, case="inject")
    assert result.verdict == "normal"
    copied = (result.case_dir / "inject.inp").read_text(encoding="utf-8")
    lines = copied.splitlines()
    assert lines[0] == "!HF def2-SVP TightSCF"
    assert lines[1:3] == ["%pal nprocs 3 end", '%base "probe"']
    assert "%maxcore 2000" in copied and "%maxcore 4000" not in copied
    assert result.injections == ("%pal nprocs 3 end", '%base "probe"')
    assert result.injections_skipped == ("%maxcore",)
    assert source.read_bytes() == before  # the generated artifact is a delivery
    # a second run into the same case refuses (the stale-scratch guard)
    with pytest.raises(BridgeError, match="Next step"):
        run_bridge(config, input_path=source, case="inject")
    # an input that already carries %pal keeps it and records the skip
    already = run_bridge(config, input_text="%pal nprocs 48 end\n!HF\n", case="present")
    assert already.injections_skipped == ("%pal",)
    assert already.injections == ("%maxcore 4000", '%base "probe"')
    kept = (already.case_dir / "present.inp").read_text(encoding="utf-8")
    assert "%pal nprocs 48 end" in kept and "%pal nprocs 3 end" not in kept


# 3 --------------------------------------------------------------------------


def test_the_classification_trio_replays_the_measured_semantics(tmp_path):
    scenarios = {
        "normal": ('print("****ORCA TERMINATED NORMALLY****")', 0),
        "orca_error": ('print("ORCA finished by error termination in LEANSCF")', 0),
        "input_rejected": (
            'print("Unknown identifier in REL block line 4 : Last token : '
            'DOVELOCITY")\nimport sys\nsys.exit(11)',
            11,
        ),
    }
    for name, (body, expected_rc) in scenarios.items():
        engine = _fake_engine(tmp_path / name, body)
        result = run_bridge(_config(tmp_path, engine), input_text="!HF\n", case=name)
        assert result.verdict == name, (name, result.verdict)
        assert result.returncode == expected_rc
        if name != "normal":
            assert result.verdict != "normal"  # rc = 0 is not success


# 4 --------------------------------------------------------------------------


def test_a_silent_zero_exit_is_suspicious_never_normal(tmp_path):
    engine = _fake_engine(tmp_path / "quiet", 'print("SCF done")\n')
    result = run_bridge(_config(tmp_path, engine), input_text="!HF\n", case="quiet")
    assert result.returncode == 0
    assert result.verdict == "suspicious"


# 5 --------------------------------------------------------------------------


SLEEPER = """
    import pathlib, time
    here = pathlib.Path.cwd()
    (here / "started.marker").write_text("x", encoding="utf-8")
    time.sleep(10)
    (here / "finished.marker").write_text("x", encoding="utf-8")
"""


def test_the_timeout_and_stall_paths_terminate_the_process(tmp_path):
    engine = _fake_engine(tmp_path / "sleeper", SLEEPER)
    wall = run_bridge(
        _config(tmp_path, engine, timeout_s=3, poll_interval_s=0.05),
        input_text="!HF\n",
        case="wall",
    )
    assert wall.verdict == "timeout" and wall.terminated_by == "wall"
    assert wall.returncode is not None
    assert (wall.case_dir / "started.marker").exists()
    assert not (wall.case_dir / "finished.marker").exists()  # killed mid-sleep
    stalled = run_bridge(
        _config(tmp_path, engine, poll_interval_s=0.05, stall_timeout_s=1),
        input_text="!HF\n",
        case="stall",
    )
    assert stalled.verdict == "timeout" and stalled.terminated_by == "stall"
    assert not (stalled.case_dir / "finished.marker").exists()


# 6 --------------------------------------------------------------------------


def test_an_externally_killed_run_is_killed(tmp_path):
    engine = _fake_engine(
        tmp_path / "crash",
        """
        import os, signal
        if os.name == "nt":
            os._exit(137)  # the 128 + SIGKILL convention shells report
        os.kill(os.getpid(), signal.SIGKILL)
        """,
    )
    result = run_bridge(_config(tmp_path, engine), input_text="!HF\n", case="crash")
    assert result.verdict == "killed"
    assert result.returncode in (-9, 137)


# 7 --------------------------------------------------------------------------


def test_replay_skips_the_run_with_the_fixed_line_and_no_process(tmp_path, monkeypatch):
    def _forbidden(*args, **kwargs):
        raise AssertionError("the bridge started a process on replay")

    monkeypatch.setattr(bridge_runner.subprocess, "Popen", _forbidden)
    lines = []
    config = _config(tmp_path, Path(sys.executable))
    result = run_bridge(
        config, input_text="!HF\n", case="replayed", is_replay=True, say=lines.append
    )
    assert lines == [bridge_runner.REPLAY_SKIP_TEXT]
    assert result.verdict == "skipped"
    assert result.returncode is None and result.case_dir is None
    assert not (tmp_path / "runs").exists()  # a replay touches nothing


# 8 --------------------------------------------------------------------------


def test_run_json_is_the_side_ledger_and_scratch_obeys_keep_scratch(tmp_path):
    engine = _fake_engine(
        tmp_path / "ledger",
        """
        import pathlib
        here = pathlib.Path.cwd()
        (here / "ledger.gbw").write_text("gbw", encoding="utf-8")
        (here / "ledger.tmp").write_text("scratch", encoding="utf-8")
        print("****ORCA TERMINATED NORMALLY****")
        """,
    )
    said = []
    result = run_bridge(
        _config(tmp_path, engine),
        input_text="!HF def2-SVP\n",
        case="ledger",
        say=said.append,
    )
    record = json.loads((result.case_dir / "run.json").read_text(encoding="utf-8"))
    assert record["verdict"] == "normal" and record["returncode"] == 0
    assert record["injections"] == ["%pal nprocs 1 end", "%maxcore 4000"]
    assert record["injections_skipped"] == []
    assert record["input_sha256"] == hashlib.sha256(
        (result.case_dir / "ledger.inp").read_bytes()
    ).hexdigest()
    assert record["output_sha256"] == hashlib.sha256(
        (result.case_dir / "ledger.out").read_bytes()
    ).hexdigest()
    assert record["config"]["nprocs"] == 1
    assert record["config"]["orca_path"] == str(engine)
    assert record["started_utc"] and record["finished_utc"]
    assert said == []  # the bridge writes nothing itself: the ledger is the
    # only carrier of times (a side product, outside the compared stream)
    trimmed = run_bridge(
        _config(tmp_path, engine, keep_scratch=False),
        input_text="!HF\n",
        case="trimmed",
    )
    assert sorted(path.name for path in trimmed.case_dir.iterdir()) == [
        "run.json",
        "trimmed.inp",
        "trimmed.out",
    ]
