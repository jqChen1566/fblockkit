"""The optional local execution bridge (E1): run a generated input on a local
engine and judge the outcome from the return code, the output text and the
heartbeat.

The subpackage is deliberately small and stdlib-only, and the core never
imports it: the layered import contract keeps the direction core <- bridge <-
ui, so dropping this package removes the execution bridge without touching the
generation and analysis layers (the rollback anchor of the E1 specification).

Public surface:

- :class:`BridgeConfig` and :class:`BridgeError` (bridge.config): the
  explicit configuration, its validation and its YAML round trip;
- :func:`run_bridge` and :class:`BridgeResult` (bridge.runner): one run --
  ``prepare`` (a fresh ``<case>/`` with the injected copy of the input),
  ``launch`` (full pathname, no shell, no mpirun), ``monitor`` (heartbeat
  polling, timeout and stall termination), ``classify`` (the three-source
  judgement table) and ``collect`` (the ``run.json`` side ledger).

A script replay never re-runs: the run point prints one fixed line and starts
no process.
"""

from .config import BridgeConfig, BridgeError
from .runner import REPLAY_SKIP_TEXT, BridgeResult, run_bridge

__all__ = ["BridgeConfig", "BridgeError", "BridgeResult", "REPLAY_SKIP_TEXT", "run_bridge"]
