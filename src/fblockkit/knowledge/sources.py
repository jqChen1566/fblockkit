"""Reference database (BibTeX) for literature evidence.

Single source of truth for citations: every evidence item of kind ``literature``
(rules, basis tables, tool index, code constants) must carry a ``bibkey`` that
exists in ``sources.bib``; loaders enforce this and tests cross-check it.

The report layer uses this module to print *complete* citations and paste-ready
BibTeX blocks, so users never have to reconstruct a reference by hand.

Deliberately dependency-free: a small BibTeX subset parser covers our own file
(@article entries with brace- or quote-delimited values; ``%`` line comments).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field as dataclass_field
from pathlib import Path
from typing import Mapping

DEFAULT_BIB_FILE = Path(__file__).resolve().parent / "sources.bib"

REQUIRED_FIELDS = ("author", "title", "journal", "year", "doi")


class BibDataError(ValueError):
    """The BibTeX database is missing or malformed."""


@dataclass(frozen=True)
class BibEntry:
    key: str
    entry_type: str
    fields: Mapping[str, str]

    def field(self, name: str) -> str:
        return self.fields.get(name.lower(), "")


def _strip_comments(text: str) -> str:
    return "\n".join(line for line in text.splitlines() if not line.lstrip().startswith("%"))


def _parse_fields(body: str) -> dict[str, str]:
    fields: dict[str, str] = {}
    index, length = 0, len(body)
    while index < length:
        while index < length and body[index] in " \t\r\n,":
            index += 1
        match = re.match(r"([A-Za-z]+)\s*=\s*", body[index:])
        if match is None:
            break
        name = match.group(1).lower()
        index += match.end()
        if index >= length:
            break
        if body[index] == "{":
            depth, cursor = 1, index + 1
            while cursor < length and depth:
                if body[cursor] == "{":
                    depth += 1
                elif body[cursor] == "}":
                    depth -= 1
                cursor += 1
            value = body[index + 1 : cursor - 1]
            index = cursor
        elif body[index] == '"':
            cursor = body.find('"', index + 1)
            value = body[index + 1 : cursor]
            index = cursor + 1
        else:
            cursor = index
            while cursor < length and body[cursor] not in ",\n":
                cursor += 1
            value = body[index:cursor]
            index = cursor
        fields[name] = " ".join(value.split())
    return fields


def _iter_entries(text: str):
    for match in re.finditer(r"@(\w+)\s*\{", text):
        if match.group(1).lower() in ("comment", "string", "preamble"):
            continue
        start = match.end() - 1
        depth, cursor = 0, start
        while cursor < len(text):
            if text[cursor] == "{":
                depth += 1
            elif text[cursor] == "}":
                depth -= 1
                if depth == 0:
                    break
            cursor += 1
        body = text[start + 1 : cursor]
        key, _, rest = body.partition(",")
        yield match.group(1).lower(), key.strip(), _parse_fields(rest)


def load_bib(path: Path | str = DEFAULT_BIB_FILE) -> tuple[BibEntry, ...]:
    """Load and validate the BibTeX database; all problems reported at once."""
    p = Path(path)
    if not p.is_file():
        raise BibDataError(f"BibTeX database not found: {p}")
    text = _strip_comments(p.read_text(encoding="utf-8"))
    entries: list[BibEntry] = []
    errors: list[str] = []
    seen: set[str] = set()
    for entry_type, key, fields in _iter_entries(text):
        where = f"{p.name}[{key or '?'}]"
        if not key:
            errors.append(f"{where}: entry without a citation key")
            continue
        if key in seen:
            errors.append(f"{where}: duplicate citation key")
            continue
        seen.add(key)
        if entry_type == "article":
            missing = [name for name in REQUIRED_FIELDS if not fields.get(name)]
            if missing:
                errors.append(f"{where}: missing required field(s) {missing}")
                continue
        entries.append(BibEntry(key=key, entry_type=entry_type, fields=fields))
    if errors:
        raise BibDataError(
            "BibTeX database failed validation (%d problem(s)):\n%s"
            % (len(errors), "\n".join(f"- {e}" for e in errors))
        )
    return tuple(entries)


def index_by_key(entries: tuple[BibEntry, ...] | None = None) -> dict[str, BibEntry]:
    return {entry.key: entry for entry in (entries if entries is not None else load_bib())}


def get(key: str, entries: tuple[BibEntry, ...] | None = None) -> BibEntry:
    table = index_by_key(entries)
    if key not in table:
        known = ", ".join(sorted(table)) or "(none)"
        raise BibDataError(
            f"Unknown bibkey {key!r}. Next step: add the entry to sources.bib "
            f"(known keys: {known})."
        )
    return table[key]


def _short_authors(author_field: str) -> str:
    out = []
    for chunk in author_field.split(" and "):
        chunk = chunk.strip()
        if not chunk:
            continue
        family, _, given = chunk.partition(",")
        initials = " ".join(f"{part[0]}." for part in given.split() if part)
        out.append(f"{family.strip()}, {initials}".strip().rstrip(","))
    return "; ".join(out)


def format_citation(entry: BibEntry) -> str:
    """One-line human-readable citation (complete: authors, year, title, venue, DOI)."""
    authors = _short_authors(entry.field("author"))
    pages = entry.field("pages").replace("--", "-")
    location = entry.field("journal")
    volume, number = entry.field("volume"), entry.field("number")
    if volume:
        location += f", {volume}"
        if number:
            location += f"({number})"
    if pages:
        location += f", {pages}"
    return (
        f"{authors} ({entry.field('year')}). {entry.field('title')}. "
        f"{location}. DOI: {entry.field('doi')}"
    )


def to_bibtex(entry: BibEntry) -> str:
    """Normalized, paste-ready BibTeX block."""
    order = ("author", "title", "journal", "year", "volume", "number", "pages", "doi")
    lines = [f"@{entry.entry_type}{{{entry.key},"]
    written = set()
    for name in order:
        value = entry.field(name)
        if value:
            lines.append(f"  {name:<8}= {{{value}}},")
            written.add(name)
    for name, value in sorted(entry.fields.items()):
        if name not in written:
            lines.append(f"  {name:<8}= {{{value}}},")
    lines.append("}")
    return "\n".join(lines)


def citation_for(bibkey: str) -> str:
    return format_citation(get(bibkey))


def bibtex_for(bibkey: str) -> str:
    return to_bibtex(get(bibkey))
