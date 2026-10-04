"""Explicit configuration of the optional local execution bridge (E1).

The bridge runs a generated input file on a local engine (ORCA) and judges the
outcome from three sources (return code, output text, heartbeat); this module
carries its configuration object, its validation and a YAML round trip.

The defaults follow the recipe layer: ``nprocs`` 1 and ``maxcore_mb`` 4000
match the ``%pal nprocs N end`` / ``%maxcore M`` shapes the input generators
inject (see recipe/casci_xas.py).  ``stall_timeout_s`` is the no-heartbeat
window, independent of the total wall clock ``timeout_s``; ``poll_interval_s``
is the monitor's polling period (for real runs the design band is 20-30 s).

Every validation failure raises :class:`BridgeError` whose message ends with a
"Next step:" sentence (the library's error convention).  The YAML helpers
import pyyaml lazily inside the functions, so the bridge package imports
stdlib only; the layered import contract keeps the direction core <- bridge <-
ui, so dropping this subpackage removes the execution bridge without touching
the generation and analysis layers.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

__all__ = [
    "BridgeConfig",
    "BridgeError",
    "FIELDS",
    "config_from_dict",
    "config_to_dict",
    "read_config",
    "write_config",
]

#: the configuration fields: the YAML keys, and the order they are written in
FIELDS = (
    "orca_path",
    "workdir_root",
    "nprocs",
    "maxcore_mb",
    "timeout_s",
    "stall_timeout_s",
    "poll_interval_s",
    "keep_scratch",
    "base_name",
)


class BridgeError(ValueError):
    """The bridge cannot run as asked (the message carries a next step)."""


@dataclass(frozen=True)
class BridgeConfig:
    """One bridge configuration (the file a script references).

    ``orca_path``: the engine executable, given as the full pathname (the
    engine manual requires the complete path) and validated to exist and be
    executable; ``workdir_root``: the run root, validated to be a usable
    directory -- every run gets a fresh ``<workdir_root>/<case>/``;
    ``nprocs`` / ``maxcore_mb``: the values of the injected ``%pal`` /
    ``%maxcore`` blocks (both skipped when the input already carries them);
    ``timeout_s``: the wall-clock limit (SIGTERM, then SIGKILL after a grace
    window); ``stall_timeout_s``: the no-heartbeat limit, independent of the
    wall clock (the ``.tmp``/``.gbw`` mtimes and the ``.out`` size are the
    heartbeat); ``poll_interval_s``: the monitor's polling period;
    ``keep_scratch``: keep ``.tmp``/``.gbw`` next to the ``.inp``/``.out``/
    ``run.json``; ``base_name``: when given, the injected ``%base "name"``
    that pins the engine's file prefix.
    """

    orca_path: str
    workdir_root: Path
    nprocs: int = 1
    maxcore_mb: int = 4000
    timeout_s: int = 3600
    stall_timeout_s: int = 3600
    poll_interval_s: float = 25.0
    keep_scratch: bool = True
    base_name: str | None = None

    def __post_init__(self) -> None:
        try:
            # Resolve to an absolute pathname now: validation runs against this
            # process's cwd, but the engine is launched with cwd = the case
            # directory, where a relative pathname would quietly point elsewhere
            # (measured: the launch then dies with the shell's "path not found").
            object.__setattr__(
                self, "orca_path", str(Path(os.fspath(self.orca_path)).resolve())
            )
        except TypeError as exc:
            raise BridgeError(
                f"orca_path {self.orca_path!r} is not a path. Next step: give "
                "the full pathname of the engine executable."
            ) from exc
        try:
            object.__setattr__(self, "workdir_root", Path(self.workdir_root))
        except TypeError as exc:
            raise BridgeError(
                f"workdir_root {self.workdir_root!r} is not a path. Next "
                "step: give the directory the runs go under."
            ) from exc
        for name in ("nprocs", "maxcore_mb", "timeout_s", "stall_timeout_s"):
            value = getattr(self, name)
            try:
                object.__setattr__(self, name, int(value))
            except (TypeError, ValueError) as exc:
                raise BridgeError(
                    f"{name} {value!r} is not a whole number ({exc}). Next "
                    "step: give an integer count."
                ) from exc
        try:
            object.__setattr__(self, "poll_interval_s", float(self.poll_interval_s))
        except (TypeError, ValueError) as exc:
            raise BridgeError(
                f"poll_interval_s {self.poll_interval_s!r} is not a number "
                f"({exc}). Next step: give the polling period in seconds."
            ) from exc
        object.__setattr__(self, "keep_scratch", bool(self.keep_scratch))
        if self.base_name is not None:
            name = str(self.base_name)
            if not name or any(char in name for char in '"\r\n'):
                raise BridgeError(
                    f"base_name {self.base_name!r} cannot be written as a "
                    '%base "..." line. Next step: give a plain name without '
                    "quotes or line breaks, or leave it empty."
                )
            object.__setattr__(self, "base_name", name)
        self._validate()

    def _validate(self) -> None:
        engine = Path(self.orca_path)
        if not engine.is_file():
            raise BridgeError(
                f"the engine executable {self.orca_path!r} does not exist (or "
                "is not a file). Next step: give the full pathname of the "
                "engine executable (the manual requires the complete path)."
            )
        if not os.access(engine, os.X_OK):
            raise BridgeError(
                f"the engine executable {self.orca_path!r} is not executable. "
                "Next step: check its permissions (chmod +x) or name the "
                "launcher instead."
            )
        root = self.workdir_root
        if root.exists():
            if not root.is_dir():
                raise BridgeError(
                    f"the run root {str(root)!r} is not a directory. Next "
                    "step: give a directory (the bridge creates one fresh "
                    "subdirectory per case under it)."
                )
            if not os.access(root, os.W_OK):
                raise BridgeError(
                    f"the run root {str(root)!r} is not writable. Next step: "
                    "pick a writable directory as the run root."
                )
        else:
            probe = root
            while not probe.exists():
                parent = probe.parent
                if parent == probe:
                    break
                probe = parent
            if not probe.is_dir() or not os.access(probe, os.W_OK):
                raise BridgeError(
                    f"the run root {str(root)!r} does not exist and cannot be "
                    "created (no writable parent directory). Next step: "
                    "create the directory (or a writable parent) first."
                )
        if self.nprocs < 1:
            raise BridgeError(
                f"nprocs {self.nprocs} is not a positive count. Next step: "
                "give the number of MPI processes the engine should use "
                "(1 = a serial run)."
            )
        if self.maxcore_mb < 1:
            raise BridgeError(
                f"maxcore_mb {self.maxcore_mb} is not a positive count. Next "
                "step: give the per-core memory in MB (the input generators "
                "default to 4000)."
            )
        if self.timeout_s < 1:
            raise BridgeError(
                f"timeout_s {self.timeout_s} is not a positive count. Next "
                "step: give the wall-clock limit in seconds."
            )
        if self.stall_timeout_s < 1:
            raise BridgeError(
                f"stall_timeout_s {self.stall_timeout_s} is not a positive "
                "count. Next step: give the no-heartbeat limit in seconds."
            )
        if self.poll_interval_s <= 0:
            raise BridgeError(
                f"poll_interval_s {self.poll_interval_s} is not positive. "
                "Next step: give the monitor's polling period in seconds."
            )


def config_to_dict(config: BridgeConfig) -> dict:
    """The YAML-ready primitive mapping of a configuration."""
    data = {name: getattr(config, name) for name in FIELDS}
    data["workdir_root"] = str(config.workdir_root)
    return data


def config_from_dict(data) -> BridgeConfig:
    """Build a configuration from a mapping (unknown keys are refused)."""
    if not isinstance(data, dict):
        raise BridgeError(
            f"the configuration is a {type(data).__name__}, not a mapping. "
            "Next step: give a YAML mapping with the documented keys."
        )
    unknown = sorted(set(data) - set(FIELDS))
    if unknown:
        raise BridgeError(
            f"the configuration carries unknown key(s) {unknown}. Next step: "
            f"use only the documented keys {list(FIELDS)}."
        )
    missing = [name for name in ("orca_path", "workdir_root") if name not in data]
    if missing:
        raise BridgeError(
            f"the configuration misses {missing}. Next step: give orca_path "
            "(the full engine pathname) and workdir_root (the run root)."
        )
    return BridgeConfig(**data)


def write_config(config: BridgeConfig, path) -> Path:
    """Write the configuration as YAML (the file a script references)."""
    import yaml  # lazy: the bridge package stays stdlib-only at import time

    target = Path(path)
    target.write_text(
        yaml.safe_dump(config_to_dict(config), sort_keys=True), encoding="utf-8"
    )
    return target


def read_config(path) -> BridgeConfig:
    """Read a configuration written by :func:`write_config` (or by hand)."""
    import yaml  # lazy: see write_config

    source = Path(path)
    if not source.is_file():
        raise BridgeError(
            f"the configuration file {str(source)!r} does not exist. Next "
            "step: write one (write_config) or check the path."
        )
    try:
        data = yaml.safe_load(source.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise BridgeError(
            f"the configuration file {str(source)!r} is not valid YAML "
            f"({exc}). Next step: fix the file syntax."
        ) from exc
    if data is None:
        raise BridgeError(
            f"the configuration file {str(source)!r} is empty. Next step: "
            "give at least orca_path and workdir_root."
        )
    return config_from_dict(data)
