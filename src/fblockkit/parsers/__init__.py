"""Parsing layer: ORCA (deep) / OpenMolcas magnetic-chain output parsing
(architecture design v0.1 §3 L1).

Public interface:

- ``parse_auto(path)``: identify the program and parse, returning a ``ParseResult``;
- ``detect_program(path)``: identification only, returning a parser or None;
- ``available()``: the registered parsers;
- ``facts_from(result)``: map a parse result onto the diagnosis rules' fact fields.

Parsers self-register at module import (for ORCA see the end of ``orca.py``).
"""

from __future__ import annotations

from typing import Any

from ..knowledge.models import ParseResult
from .base import (
    Parser,
    ParserError,
    available,
    detect_program,
    parse_auto,
    read_text,
    register,
)
from . import orca  # noqa: F401  (import registers the ORCA parser and exposes orca.facts)

__all__ = [
    "Parser",
    "ParserError",
    "available",
    "detect_program",
    "parse_auto",
    "read_text",
    "register",
    "facts_from",
]


def facts_from(result: ParseResult) -> dict[str, Any]:
    """Dispatch fact extraction by program (the diagnosis layer's input interface).

    The field vocabulary is documented in ``knowledge/rules/README.md``; a program
    with no fact extractor raises an explicit error rather than silently returning
    an empty table.
    """
    if result.program == orca.PROGRAM:
        return orca.facts(result)
    raise ParserError(
        f"no fact extractor for {result.program} yet. Next step: implement that "
        f"program's facts mapping in the parsers layer and register its field "
        f"vocabulary."
    )
