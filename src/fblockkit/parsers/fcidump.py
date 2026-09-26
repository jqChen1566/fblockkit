"""FCIDUMP reader (an auxiliary data file, not an engine output).

An FCIDUMP carries one electronic-structure Hamiltonian in a small orbital
space: the effective one-electron integrals, the two-electron integrals in
chemist notation and the core energy.  The analysis layer uses it to solve the
CAS-CI of a converged ORCA CASSCF reference itself and to rebuild the
spin-resolved density matrices -- that is how fBlockKit reaches the literal
four-state entropy (the spin-summed natural occupations ORCA prints are not
enough, and ORCA's own 2-RDM export is unreachable for CAS-type methods in
6.1.1; see the module docstring of :mod:`fblockkit.analysis.entropy_rdm`).

How to make ORCA write one: run the CASSCF to convergence first, then rerun it
with the converged orbitals and the ``!FCIDUMP`` keyword
(``!moread`` + ``%moinp``); the dump run saves the integrals of the current
orbitals after its first macro-iteration and then exits with the
"IS NOT FULLY CONVERGED" warning -- that message is inherent to the dump mode
and does not affect the file.

The reader is deliberately **not** registered in the ``parse_auto`` registry:
that registry identifies engine *output*; an FCIDUMP is a data file the user
produces on purpose and handles explicitly.

Format (the Molpro/ORCA convention, all indices 1-based in the file):

- header: ``&FCI NORB=...,NELEC=...,MS2=...`` lines up to a lone ``/``;
  ``ORBSYM`` and ``ISYM`` are read when present;
- body: five fields per line, ``value p q r s``; ``r = s = 0`` marks a
  one-electron integral, ``p = q = r = s = 0`` the core energy, everything
  else a two-electron integral ``(pq|rs)`` listed once and expanded here with
  the full 8-fold permutation symmetry.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .base import ParserError, read_text, split_fields


@dataclass(frozen=True)
class Fcidump:
    """One active-space Hamiltonian (all indices 0-based from here on)."""

    norb: int
    nelec: int
    ms2: int
    ecore: float
    h: tuple[tuple[float, ...], ...]
    #: chemist-notation (pq|rs), expanded to the full 8-fold symmetry
    g: dict[tuple[int, int, int, int], float]
    orbsym: tuple[int, ...]
    isym: int


_HEADER_KEYS = ("NORB", "NELEC", "MS2")


def _header_value(header: str, key: str) -> int:
    for token in header.replace(",", " ").split():
        if token.startswith(key + "="):
            value = token.split("=", 1)[1]
            try:
                return int(value)
            except ValueError:
                break
    raise ParserError(
        f"the FCIDUMP header has no usable {key}= entry (header: {header!r}). "
        "Next step: check that the file really comes from ORCA's !FCIDUMP "
        "keyword and was not truncated."
    )


def parse_fcidump(path: str | Path) -> Fcidump:
    """Read an FCIDUMP; raise ParserError with a next-step hint on any defect."""
    text = read_text(path)
    lines = text.splitlines()
    header_lines: list[str] = []
    body_at = None
    for i, line in enumerate(lines):
        if line.strip() == "/":
            body_at = i + 1
            break
        header_lines.append(line)
    if body_at is None or not header_lines:
        raise ParserError(
            f"no FCIDUMP header found in {path} (the file must start with "
            "&FCI ... / ). Next step: check that this is an FCIDUMP and not an "
            "ORCA output file."
        )
    header = " ".join(header_lines)
    if "&FCI" not in header.upper():
        raise ParserError(
            f"the file does not start with an &FCI header: {header_lines[0]!r}. "
            "Next step: pass the FCIDUMP file itself (the one ORCA names FCIDUMP), "
            "not the output of the run that wrote it."
        )
    norb = _header_value(header, "NORB")
    nelec = _header_value(header, "NELEC")
    ms2 = _header_value(header, "MS2")
    if norb <= 0 or nelec < 0:
        raise ParserError(
            f"implausible FCIDUMP header values (NORB={norb}, NELEC={nelec}). "
            "Next step: inspect the header; a valid dump always carries "
            "positive NORB."
        )
    orbsym: tuple[int, ...] = ()
    isym = 1
    for line in header_lines:
        tokens = line.replace(",", " ").split()
        if tokens and tokens[0].upper().startswith("ORBSYM"):
            values = []
            for token in tokens:
                for piece in token.replace("=", " ").split():
                    if piece.isdigit():
                        values.append(int(piece))
            orbsym = tuple(values)
        if tokens and tokens[0].upper().startswith("ISYM"):
            digits = []
            for token in tokens:
                for piece in token.replace("=", " ").split():
                    if piece.isdigit():
                        digits.append(int(piece))
            if digits:
                isym = digits[0]

    h = [[0.0] * norb for _ in range(norb)]
    g: dict[tuple[int, int, int, int], float] = {}
    ecore: float | None = None
    for line in lines[body_at:]:
        fields = split_fields(line)
        if not fields:
            continue
        if len(fields) != 5:
            raise ParserError(
                f"malformed FCIDUMP body line (expected 5 fields): {line!r}. "
                "Next step: check the file for truncation or a concatenation of "
                "two dumps."
            )
        try:
            value = float(fields[0])
            p, q, r, s = (int(t) for t in fields[1:])
        except ValueError as exc:
            raise ParserError(
                f"non-numeric FCIDUMP body line: {line!r}. Next step: check the "
                "file for truncation."
            ) from exc
        if p == 0 and q == 0 and r == 0 and s == 0:
            ecore = value
            continue
        if r == 0 and s == 0:
            if not (1 <= p <= norb and 1 <= q <= norb):
                raise ParserError(
                    f"one-electron index out of range in line {line!r} "
                    f"(NORB={norb}). Next step: the header and the body disagree; "
                    "check the file for truncation."
                )
            h[p - 1][q - 1] = h[q - 1][p - 1] = value
            continue
        if not (1 <= p <= norb and 1 <= q <= norb and 1 <= r <= norb and 1 <= s <= norb):
            raise ParserError(
                f"two-electron index out of range in line {line!r} (NORB={norb}). "
                "Next step: the header and the body disagree; check the file for "
                "truncation."
            )
        p, q, r, s = p - 1, q - 1, r - 1, s - 1
        for perm in (
            (p, q, r, s), (q, p, r, s), (p, q, s, r), (q, p, s, r),
            (r, s, p, q), (s, r, p, q), (r, s, q, p), (s, r, q, p),
        ):
            g[perm] = value
    if ecore is None:
        raise ParserError(
            f"the FCIDUMP has no core-energy line (value 0 0 0 0). Next step: "
            "check the file for truncation; the last line of a complete dump is "
            "the core energy."
        )
    return Fcidump(
        norb=norb,
        nelec=nelec,
        ms2=ms2,
        ecore=ecore,
        h=tuple(tuple(row) for row in h),
        g=g,
        orbsym=orbsym,
        isym=isym,
    )


__all__ = ["Fcidump", "parse_fcidump"]
