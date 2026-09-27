"""Reader for ORCA's property file (``<base>.property.txt``).

Measured format (fixtures: ``n2_sa.property.txt``, ``h2o_sa.property.txt``,
ORCA 6.1.1; the schema of every section is listed in the manual's property-file
appendix, ``contents/utilitiesvisualization/property_file_list``):

- ``$SECTION_NAME`` opens a section, ``$End`` closes it;
- inside, one field per line::

      &finalEnergy [&Type "Double"]      -1.0870707440711183e+02  "Final GS or SA energy"
      &totalEnergy [&Type "ArrayOfDoubles", &Dim (3,1)] "Total energy of each state"

  the bracketed block carries ``&Type`` (Double, Integer, String, Boolean,
  ArrayOfDoubles, ArrayOfIntegers, Coordinates, ...) and optionally ``&Units``
  and ``&Dim``; the value follows the bracket -- a bare number, a quoted
  string, or, for arrays, a block on the following lines: a column-index
  header (``0`` or ``0 1 2 ...``), a blank line, then one row per element
  (``<index> <values...>``) until the next ``&field`` or ``$End``;
- comments in double quotes may follow a value (they are attached to the
  field, not parsed as data).

The reader returns every section with its fields as typed Python values
(``PropertyField.value``); arrays keep their row structure (a tuple of
row tuples).  Unknown ``&Type`` names raise -- silently dropping a value the
manual defines would hide data, and a new type is a format fact to record.
``orca_2json <base> -property`` converts the same file to JSON (the exported
key names match the field names; measured).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .base import ParserError, read_text

__all__ = ["PropertyError", "PropertyField", "PropertySection", "parse_property"]

_HEADER_RE = re.compile(r"^\$(?!End\s*$)([A-Za-z_][A-Za-z0-9_]*)\s*$")
_FIELD_RE = re.compile(r"^\s*&([A-Za-z0-9_]+)\s*(\[.*?\])?\s*(.*)$")
_TYPE_RE = re.compile(r'&Type\s+"([^"]+)"')
_UNITS_RE = re.compile(r'&Units\s+"([^"]+)"')
_DIM_RE = re.compile(r"&Dim\s*\(\s*(\d+)\s*(?:,\s*(\d+)\s*)?\)")
_QUOTED_RE = re.compile(r'"([^"]*)"')

#: The value types the manual's property-file appendix uses; measured on the
#: two CASSCF fixtures.  Any other type is a format fact to record, not to drop.
_SCALAR_TYPES = ("Double", "Integer", "String", "Boolean")
_ARRAY_TYPES = ("ArrayOfDoubles", "ArrayOfIntegers")


class PropertyError(ParserError):
    """The property file cannot be read (with a next step)."""


@dataclass(frozen=True)
class PropertyField:
    """One ``&field`` with its type annotation and parsed value."""

    name: str
    type: str
    units: str | None
    dim: tuple[int, ...] | None
    comment: str
    value: Any


@dataclass(frozen=True)
class PropertySection:
    """One ``$SECTION`` with its fields in file order."""

    name: str
    fields: tuple[PropertyField, ...]

    def field(self, name: str) -> PropertyField | None:
        for item in self.fields:
            if item.name == name:
                return item
        return None


def _fail(detail: str) -> PropertyError:
    return PropertyError(
        f"{detail} Next step: give the ORCA property file (``<base>.property.txt``, "
        "written automatically by the job)."
    )


def _parse_array_rows(
    lines: list[str], position: int, *, as_double: bool
) -> tuple[tuple[Any, ...], ...]:
    """Rows of one array block; ``position`` points at the line after the field.

    Measured layout: an array is a sequence of **column groups** (up to eight
    columns each, like the normal-modes block); each group is a header line of
    column indices, a blank line, then one row per element (``<index>
    <values...>``).  A header is distinguishable from a row because a row
    carries exactly one more token than the group's column count.  The rows
    are returned row-major with the column index deciding the position, so a
    wide array's groups are stitched back together.
    """
    raw: list[list[str]] = []
    while position < len(lines):
        line = lines[position]
        stripped = line.strip()
        if stripped.startswith("&") or _HEADER_RE.match(line) or stripped == "$End":
            break
        position += 1
        if stripped:
            raw.append(stripped.split())
    collected: dict[int, dict[int, Any]] = {}
    index = 0
    convert = float if as_double else int
    while index < len(raw):
        header = raw[index]
        if not all(token.lstrip("+-").isdigit() for token in header):
            break
        columns = [int(token) for token in header]
        index += 1
        rows_here = 0
        while index < len(raw) and len(raw[index]) == 1 + len(columns):
            row = raw[index]
            if not row[0].lstrip("+-").isdigit():
                break
            row_index = int(row[0])
            slot = collected.setdefault(row_index, {})
            for column, value in zip(columns, row[1:]):
                slot[column] = convert(value)
            index += 1
            rows_here += 1
        if rows_here == 0:
            break
    if not collected:
        return ()
    width = max(max(slot) for slot in collected.values()) + 1
    return tuple(
        tuple(collected[row_index].get(column, 0.0 if as_double else 0) for column in range(width))
        for row_index in sorted(collected)
    )


def parse_property(path: str | Path) -> tuple[PropertySection, ...]:
    """Read an ORCA property file into typed sections (file order kept)."""
    text = read_text(path)
    lines = text.splitlines()
    sections: list[PropertySection] = []
    current_name: str | None = None
    current_fields: list[PropertyField] = []
    position = 0
    while position < len(lines):
        line = lines[position]
        position += 1
        match = _HEADER_RE.match(line)
        if match is not None:
            if current_name is not None:
                sections.append(PropertySection(current_name, tuple(current_fields)))
            current_name = match.group(1)
            current_fields = []
            continue
        stripped = line.strip()
        if stripped == "$End":
            if current_name is not None:
                sections.append(PropertySection(current_name, tuple(current_fields)))
                current_name = None
                current_fields = []
            continue
        if current_name is None:
            continue  # the banner before the first section
        field_match = _FIELD_RE.match(line)
        if field_match is None:
            continue
        name = field_match.group(1)
        bracket = field_match.group(2) or ""
        rest = field_match.group(3).strip()
        type_match = _TYPE_RE.search(bracket)
        if type_match is None:
            # measured: a bare integer carries no annotation (&GeometryIndex 1);
            # anything else without a type is a format fact to record
            tokens = rest.split()
            if tokens and tokens[0].lstrip("+-").isdigit():
                field_type = "Integer"
            else:
                raise _fail(
                    f"the field '{name}' in ${current_name} carries no &Type annotation "
                    "and its value is not a bare integer."
                )
        else:
            field_type = type_match.group(1)
        units_match = _UNITS_RE.search(bracket)
        dim_match = _DIM_RE.search(bracket)
        dim = (
            tuple(int(value) for value in dim_match.groups() if value is not None)
            if dim_match is not None
            else None
        )
        comment_match = _QUOTED_RE.search(rest)
        if field_type in _SCALAR_TYPES:
            if field_type == "String":
                value: Any = comment_match.group(1) if comment_match is not None else ""
                comment = ""
            else:
                tokens = rest.split()
                bare = tokens[0] if tokens else ""
                try:
                    if field_type == "Double":
                        value = float(bare)
                    elif field_type == "Integer":
                        value = int(bare)
                    else:  # Boolean
                        if bare.lower() not in ("true", "false"):
                            raise ValueError(bare)
                        value = bare.lower() == "true"
                except ValueError as exc:
                    raise _fail(
                        f"the field '{name}' in ${current_name} is {rest!r}, not a "
                        f"{field_type} value."
                    ) from exc
                comment = comment_match.group(1) if comment_match is not None else ""
        elif field_type in _ARRAY_TYPES:
            value = _parse_array_rows(
                lines, position, as_double=(field_type == "ArrayOfDoubles")
            )
            position = _advance_past_array(lines, position)
            comment = comment_match.group(1) if comment_match is not None else ""
        elif field_type == "Coordinates":
            value = _parse_array_rows(lines, position, as_double=True)
            position = _advance_past_array(lines, position)
            comment = comment_match.group(1) if comment_match is not None else ""
        else:
            raise _fail(
                f"the field '{name}' in ${current_name} has type {field_type!r}, which "
                "this reader does not know; record the format before extending it."
            )
        current_fields.append(
            PropertyField(
                name=name,
                type=field_type,
                units=units_match.group(1) if units_match is not None else None,
                dim=dim,
                comment=comment,
                value=value,
            )
        )
    if current_name is not None:
        sections.append(PropertySection(current_name, tuple(current_fields)))
    if not sections:
        raise _fail(f"{path} carries no '$SECTION' blocks.")
    return tuple(sections)


def _advance_past_array(lines: list[str], position: int) -> int:
    """The first line at/after ``position`` that starts a new field, section or $End."""
    while position < len(lines):
        line = lines[position]
        stripped = line.strip()
        if stripped.startswith("&") or stripped == "$End" or _HEADER_RE.match(line):
            return position
        position += 1
    return position
