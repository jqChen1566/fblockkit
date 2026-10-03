"""Gaussian formatted-checkpoint (.fch) reader -- the MOKIT automr side products.

MOKIT's ``automr`` drops one ``.fch`` per stage next to its input (the file
names carry the stage: ``_rhf``, ``_uhf``, ``_uno``, ``_asrot...``,
``_CASSCF_NO``); they are genuine Gaussian formatted-checkpoint text files
(written by the Gaussian backend or by MOKIT's own writer).  This reader
extracts the facts the automr report shows.

Format, **every convention measured** on MOKIT 1.2.8 products
(fixtures ``fixtures/mokit/*.fch``, 2026-10-02):

- the first two lines are bare: the title string and the level line
  (``Stability RHF ... CC-pVDZ``); every later section is a label line --
  ``<Label>  <T>  N= <count>`` for arrays, ``<Label>  <T>  <value>`` for
  scalars -- with the data on the following lines (whitespace-separated,
  five values per line for arrays).  String sections (``Full Title``,
  ``Route``) carry their text on the line(s) after the label;
- the **MO coefficient arrays are column-major over (basis, MO)**: the flat
  index is ``m * nbf + p``; the **density arrays are the lower triangle in
  row-major order** (row outer).  Both layouts were pinned against the
  file's own ``Total SCF Density``: the RHF fixture's five closed-shell
  molecular orbitals, occupied with factor two, reconstruct the stored
  triangle to 9e-9 -- the regression anchor in the tests;
- **coordinates are in Bohr** (the format's own unit) and the dipole-moment
  section is in atomic units; the reader reports both as stored;
- GVB/UNO-transformed files carry **alpha sections only** (no beta MO
  block); occupancy information is not part of the format (it stays in the
  automr text output).

The reader is deliberately **not** registered in ``parse_auto``: a .fch is a
data file the user produced on purpose (like an FCIDUMP), not an engine
output.  It refuses: a file that does not open with the title/level pair,
missing ``Number of basis functions`` / ``Alpha MO coefficients`` /
coordinates, arrays whose declared count was not reached, and coefficient
arrays whose length is not ``nbf**2``.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from .base import ParserError, read_text

__all__ = ["BOHR_TO_ANGSTROM", "FchFile", "read_fch", "stage_label"]

_ARRAY = re.compile(r"^(?P<label>\S.*?)\s{2,}(?P<kind>[ICRF])\s+N=\s*(?P<count>\d+)\s*$")
_SCALAR = re.compile(
    r"^(?P<label>\S.*?)\s{2,}(?P<kind>[IR])\s+"
    r"(?P<value>[-+]?\d+(?:\.\d*)?(?:[DEde][-+]?\d+)?)\s*$"
)

#: Bohr radius (CODATA 2018), for callers that want Angstrom.
BOHR_TO_ANGSTROM = 0.529177210903

_STAGE_WORDS = (
    ("uno_asrot2gvb", "the UNO active-space rotation towards the GVB stage"),
    ("uno_asrot", "the UNO active-space rotation"),
    ("casscf_no", "the CASSCF natural orbitals"),
    ("gvb", "the GVB stage"),
    ("uno", "the UNO rotation"),
    ("rhf", "the RHF reference"),
    ("uhf", "the UHF reference"),
)


def stage_label(name: str) -> str:
    """A name-derived stage note for a side-product file name (best effort).

    The semantics live in MOKIT's own (auto-generated) naming; this helper
    composes a short note from the observed tokens and returns the bare name
    for unknown patterns.  The report marks these notes as name-derived.
    """
    stem = Path(name).stem.lower()
    for token, text in _STAGE_WORDS:
        if token in stem:
            return text
    return Path(name).stem


@dataclass(frozen=True)
class FchFile:
    """The measured content of one .fch side product."""

    name: str
    title: str
    level: str
    route: str | None
    charge: int
    multiplicity: int
    n_electrons: int
    n_alpha: int | None
    n_beta: int | None
    nbf: int
    atomic_numbers: tuple[int, ...]
    coordinates_bohr: tuple[tuple[float, float, float], ...]
    alpha_energies: tuple[float, ...]
    #: coefficient matrix, indexed ``[basis][MO]`` (column-major as stored)
    alpha_coefficients: tuple[tuple[float, ...], ...]
    beta_energies: tuple[float, ...] | None
    beta_coefficients: tuple[tuple[float, ...], ...] | None
    scf_energy: float | None
    total_energy: float | None
    dipole_au: tuple[float, float, float] | None
    #: the stored lower triangle (row-major), kept for the layout regression check
    total_scf_density: tuple[float, ...] | None


def _scan(lines: list[str]):
    """Yield ``(label, kind, count_or_scalar, data, text)`` per section.

    ``data`` carries the numeric arrays (I/R), ``text`` the joined string
    sections (C); scalars land in ``count_or_scalar`` with both empty.
    """
    index = 2  # the two bare header lines
    while index < len(lines):
        line = lines[index]
        match = _ARRAY.match(line)
        if match is not None:
            label = match.group("label").strip()
            kind = match.group("kind")
            count = int(match.group("count"))
            if kind == "C":
                collected: list[str] = []
                cursor = index + 1
                while (
                    cursor < len(lines)
                    and _ARRAY.match(lines[cursor]) is None
                    and _SCALAR.match(lines[cursor]) is None
                ):
                    stripped = lines[cursor].strip()
                    if stripped:
                        collected.append(stripped)
                    cursor += 1
                yield label, kind, count, [], " ".join(collected)
                index = cursor
                continue
            data: list[float] = []
            cursor = index + 1
            while len(data) < count and cursor < len(lines):
                if (
                    _ARRAY.match(lines[cursor]) is not None
                    or _SCALAR.match(lines[cursor]) is not None
                ):
                    # the array ended early: stop here instead of swallowing
                    # the next section's numbers as data; the caller checks
                    # the declared count and refuses the file
                    break
                for token in lines[cursor].split():
                    try:
                        data.append(float(token.replace("D", "E").replace("d", "e")))
                    except ValueError:
                        pass
                cursor += 1
            yield label, kind, count, data, ""
            index = cursor
            continue
        match = _SCALAR.match(line)
        if match is not None:
            value = float(match.group("value").replace("D", "E").replace("d", "e"))
            yield match.group("label").strip(), match.group("kind"), value, [], ""
            index += 1
            continue
        index += 1


def read_fch(path: Path | str) -> FchFile:
    """Read one .fch file into :class:`FchFile` (see the module docstring)."""
    path = Path(path)
    try:
        text = read_text(path)
    except OSError as exc:
        raise ParserError(
            f"the .fch file cannot be read ({exc}). Next step: check the path; "
            "automr writes its side products next to the input."
        ) from exc
    lines = text.splitlines()
    if len(lines) < 3 or not lines[0].strip() or not lines[1].strip():
        raise ParserError(
            "the file does not open with the .fch title/level line pair. Next "
            "step: check that it is a Gaussian formatted checkpoint (automr's "
            "side products are text files beginning with a title line)."
        )
    arrays: dict[str, list[float]] = {}
    scalars: dict[str, float] = {}
    strings: dict[str, str] = {}
    for label, kind, value, data, text in _scan(lines):
        if kind == "C":
            strings[label] = text
        elif data:
            if len(data) != value:
                raise ParserError(
                    f"the section {label!r} declares {value} values but "
                    f"{len(data)} were read. Next step: check the file for "
                    "truncation."
                )
            arrays[label] = data
        else:
            scalars[label] = value

    def _scalar_int(label: str) -> int | None:
        value = scalars.get(label)
        return None if value is None else int(value)

    nbf = _scalar_int("Number of basis functions")
    if nbf is None:
        raise ParserError(
            "the file carries no 'Number of basis functions' section. Next "
            "step: check that the file is a complete .fch."
        )
    charges = arrays.get("Atomic numbers")
    coordinates = arrays.get("Current cartesian coordinates")
    if not charges or coordinates is None:
        raise ParserError(
            "the file carries no atomic numbers or coordinates. Next step: "
            "check that the file is a complete .fch."
        )
    if len(coordinates) != 3 * len(charges):
        raise ParserError(
            f"the coordinate count ({len(coordinates)}) does not match three per "
            f"atom ({len(charges)}). Next step: check the file for truncation."
        )
    coefficients = arrays.get("Alpha MO coefficients")
    if coefficients is None:
        raise ParserError(
            "the file carries no 'Alpha MO coefficients' section. Next step: "
            "check that the file is a complete .fch."
        )
    if len(coefficients) != nbf * nbf:
        raise ParserError(
            f"the coefficient array has {len(coefficients)} values, not nbf**2 = "
            f"{nbf * nbf}. Next step: check the file for truncation."
        )

    def _matrix(flat: list[float]) -> tuple[tuple[float, ...], ...]:
        # column-major over (basis, MO): flat index m * nbf + p
        return tuple(
            tuple(flat[m * nbf + p] for m in range(nbf)) for p in range(nbf)
        )

    beta_flat = arrays.get("Beta MO coefficients")
    dipole = arrays.get("Dipole Moment")
    charge = _scalar_int("Charge")
    multiplicity = _scalar_int("Multiplicity")
    n_electrons = _scalar_int("Number of electrons")
    for label, value in (
        ("Charge", charge),
        ("Multiplicity", multiplicity),
        ("Number of electrons", n_electrons),
    ):
        if value is None:
            raise ParserError(
                f"the file carries no {label!r} field. Next step: check that "
                "the file is a complete Gaussian formatted checkpoint."
            )
    return FchFile(
        name=path.name,
        title=lines[0].strip(),
        level=lines[1].strip(),
        route=strings.get("Route"),
        charge=charge,
        multiplicity=multiplicity,
        n_electrons=n_electrons,
        n_alpha=_scalar_int("Number of alpha electrons"),
        n_beta=_scalar_int("Number of beta electrons"),
        nbf=nbf,
        atomic_numbers=tuple(int(value) for value in charges),
        coordinates_bohr=tuple(
            (coordinates[3 * i], coordinates[3 * i + 1], coordinates[3 * i + 2])
            for i in range(len(charges))
        ),
        alpha_energies=tuple(arrays.get("Alpha Orbital Energies", [])),
        alpha_coefficients=_matrix(coefficients),
        beta_energies=(
            tuple(arrays["Beta Orbital Energies"])
            if "Beta Orbital Energies" in arrays
            else None
        ),
        beta_coefficients=(
            _matrix(beta_flat)
            if beta_flat is not None and len(beta_flat) == nbf * nbf
            else None
        ),
        scf_energy=scalars.get("SCF Energy"),
        total_energy=scalars.get("Total Energy"),
        dipole_au=tuple(dipole) if dipole is not None and len(dipole) == 3 else None,  # type: ignore[arg-type]
        total_scf_density=(
            tuple(arrays["Total SCF Density"]) if "Total SCF Density" in arrays else None
        ),
    )
