"""The local ORCA basis-set library query (plan item 6.2: the U1 datasource).

The deployed basis library (the basisdb project: a Basis Set Exchange v0.12
snapshot with 776 sets, deployed on 101 under ``~/projects/orca_basis_sets``)
carries a SQLite metadata index, ``basis.db``.  This module opens that index
read-only with the standard library and answers "which sets cover element X"
from its metadata -- names, role (orbital / fit / guess / ...), family,
relativistic method, the ORCA built-in flag, contracted-function size and
reference count.  No library content is bundled with the package: the query
layer runs against the deployment the caller points at.

Lookup order: an explicit path, then ``$FBK_BASISDB``, then
``~/projects/orca_basis_sets``.  Measured schema facts (2026-09-29, library
v0.12): ``elements_json`` is a JSON array of element symbols (matched with a
quoted LIKE so substrings cannot collide); ``orca_builtin = 1`` means the set
is reachable through ORCA's own ``!name`` keyword, otherwise the file must
load through ``GTOName`` with an absolute path (the library's own
BASIS_USAGE.md documents both paths); sizes count contracted functions over
the set's full element range.
"""

from __future__ import annotations

import json
import os
import sqlite3
from pathlib import Path
from typing import Any

__all__ = ["BasisdbError", "find_library", "query_element", "render_query"]

_ROLES = ("admmfit", "dftjfit", "dftxfit", "guess", "jfit", "jkfit", "optri", "rifit")
_QUERY = """
SELECT bse_name, role, family, relativistic_method, orca_builtin,
       has_ecp, total_contracted_functions, ref_count, elements_json
FROM basis
WHERE elements_json LIKE ? AND role = ?
ORDER BY family, bse_name
"""


class BasisdbError(Exception):
    """A refusal with the next step spelled out."""


def find_library(explicit: str | Path | None = None) -> Path | None:
    """The library's basis.db: explicit path, $FBK_BASISDB, then ~/projects/..."""
    candidates: list[Path] = []
    for base in (explicit, os.environ.get("FBK_BASISDB")):
        if not base:
            continue
        path = Path(base)
        # an existing file is taken as the index itself (the shipped
        # miniature is named mini_basis.db); anything else is a directory
        candidates.append(path if path.is_file() else path / "basis.db")
    candidates.append(Path.home() / "projects" / "orca_basis_sets" / "basis.db")
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return None


def _rows(db_path: Path, symbol: str, role: str) -> list[dict[str, Any]]:
    try:
        connection = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    except sqlite3.Error as exc:
        raise BasisdbError(f"cannot open {db_path}: {exc}") from exc
    try:
        cursor = connection.execute(_QUERY, (f'%"{symbol}"%', role))
        entries: list[dict[str, Any]] = []
        for row in cursor.fetchall():
            entries.append(
                {
                    "name": row[0],
                    "role": row[1],
                    "family": row[2],
                    "relativistic": row[3],
                    "orca_builtin": bool(row[4]),
                    "has_ecp": bool(row[5]),
                    "size": row[6],
                    "ref_count": row[7],
                    "elements": json.loads(row[8]) if row[8] else [],
                }
            )
        return entries
    except sqlite3.Error as exc:
        raise BasisdbError(
            f"{db_path} does not look like a basis library ({exc})"
        ) from exc
    finally:
        connection.close()


def query_element(db_path: Path, symbol: str) -> dict[str, Any]:
    """All library entries covering an element, split into orbital/auxiliary."""
    symbol = symbol.strip()
    if not (1 <= len(symbol) <= 2) or not symbol[0].isupper() or not symbol.isalpha():
        raise BasisdbError(f"not an element symbol: {symbol!r}")
    orbital = _rows(db_path, symbol, "orbital")
    auxiliary: dict[str, int] = {}
    for role in _ROLES:
        count = len(_rows(db_path, symbol, role))
        if count:
            auxiliary[role] = count
    return {"symbol": symbol, "orbital": orbital, "auxiliary_counts": auxiliary}


def render_query(db_path: Path, symbol: str, *, max_rows: int = 60) -> str:
    """The menu-4 library-query text."""
    data = query_element(db_path, symbol)
    orbital = data["orbital"]
    auxiliary = data["auxiliary_counts"]
    total = len(orbital) + sum(auxiliary.values())
    if total == 0:
        raise BasisdbError(
            f"the library has no entry covering {data['symbol']}. Next step: "
            "check the element symbol (an ORCA-supported element), or the "
            "library version -- the deployment carries a BSE snapshot."
        )
    lines = [f"Basis-library query: {data['symbol']} at {db_path}"]
    aux_text = (
        " + ".join(f"{count} {role}" for role, count in sorted(auxiliary.items()))
        or "none"
    )
    lines.append(f"  {total} entries: {len(orbital)} orbital + {aux_text} auxiliary")
    lines.append("")
    lines.append(
        "  -- orbital basis sets ([built-in] = ORCA !name; [external] = GTOName) --"
    )
    for entry in orbital[:max_rows]:
        tag = "built-in" if entry["orca_builtin"] else "external"
        ecp = "  ECP" if entry["has_ecp"] else ""
        lines.append(
            f"    [{tag}] {entry['name']:<26} {entry['family'] or '?':<10}"
            f" rel={entry['relativistic'] or '?':<8} {entry['size'] or 0:>6} fn"
            f"  refs {entry['ref_count'] or 0}{ecp}"
        )
    if len(orbital) > max_rows:
        lines.append(
            f"    ... ({len(orbital) - max_rows} more; the full listing is INDEX.md"
            " in the library)"
        )
    if auxiliary:
        lines.append("")
        lines.append(
            "  -- auxiliary roles: "
            + ", ".join(f"{role} {count}" for role, count in sorted(auxiliary.items()))
        )
    lines.append("")
    lines.append(
        "  Reading notes: ECP marks sets that pair with an effective-core "
        "potential for this element; external sets load via GTOName with an "
        "absolute path (the library's BASIS_USAGE.md documents both paths); "
        "sizes count contracted functions over the set's full element range."
    )
    return "\n".join(lines)
