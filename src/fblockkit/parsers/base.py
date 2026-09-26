"""Parser-layer public interface and registry (architecture design v0.1 §3 L1).

Every parser implements two methods:

- ``detect(path) -> bool``: is this file the output of that program (format
  recognition only, no full parse);
- ``parse(path) -> ParseResult``: parse the file into a structured result.

Error discipline: an unrecognised format must raise ParserError with a
"Next step: " hint, and must never silently return an empty result. The callers
(diagnosis layer / user interface) translate that error for the user.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Protocol, runtime_checkable

from ..knowledge.models import ParseResult

# Encoding policy when reading a file: ORCA output is ASCII/UTF-8; undecodable
# bytes are swallowed as replacement characters so that parsing never stops on a
# real production file (which often carries mangled paths and parallel-library
# output).
_ENCODING = "utf-8"
_ENCODING_ERRORS = "replace"


class ParserError(Exception):
    """Parse failure: file unreadable, format unrecognised, or a required block missing."""


@runtime_checkable
class Parser(Protocol):
    """One program-output parser."""

    program: str  # program name (lower case), e.g. "orca"

    def detect(self, path: str | Path) -> bool:  # pragma: no cover - protocol declaration
        ...

    def parse(self, path: str | Path) -> ParseResult:  # pragma: no cover - protocol declaration
        ...


_REGISTRY: list[Parser] = []


def register(parser: Parser) -> Parser:
    """Register a parser (called at module import; one parser per program name)."""
    _REGISTRY[:] = [p for p in _REGISTRY if p.program != parser.program]
    _REGISTRY.append(parser)
    return parser


def available() -> tuple[Parser, ...]:
    """The registered parsers."""
    return tuple(_REGISTRY)


def read_text(path: str | Path) -> str:
    """Read an output file as text; raise ParserError when missing or unreadable."""
    p = Path(path)
    if not p.exists():
        raise ParserError(
            f"file does not exist: {p}. Next step: check the spelling of the path; "
            "if the job has not finished yet, wait for it to complete before parsing."
        )
    if p.is_dir():
        raise ParserError(
            f"the path is a directory, not a file: {p}. Next step: give the concrete "
            "output file name (e.g. job.out)."
        )
    try:
        return p.read_text(encoding=_ENCODING, errors=_ENCODING_ERRORS)
    except OSError as exc:
        raise ParserError(
            f"file is not readable: {p} ({exc}). Next step: check the read permission "
            "or whether the file is locked by another process."
        ) from exc


def detect_program(path: str | Path) -> Parser | None:
    """Identify the program automatically; return None when no parser matches (no guessing)."""
    for parser in _REGISTRY:
        try:
            if parser.detect(path):
                return parser
        except ParserError:
            continue
    return None


def parse_auto(path: str | Path) -> ParseResult:
    """Identify and parse; an unrecognised file gives an actionable error."""
    parser = detect_program(path)
    if parser is None:
        known = ", ".join(p.program for p in _REGISTRY) or "(none)"
        raise ParserError(
            f"cannot recognise the output format of {path} (registered parsers: {known}). "
            "Next step: confirm the file comes from a supported program; if it is a "
            "version this tool does not support yet, submit a sample file so the parser "
            "can be extended."
        )
    return parser.parse(path)


def float_or_none(token: str) -> float | None:
    """Lenient float conversion (Fortran-style D exponents are accepted)."""
    try:
        return float(token.replace("D", "E").replace("d", "e"))
    except (TypeError, ValueError):
        return None


# Numeric field splitting shared by the parsers (runs of whitespace)
_WS = re.compile(r"\s+")


def split_fields(line: str) -> list[str]:
    return [t for t in _WS.split(line.strip()) if t]
