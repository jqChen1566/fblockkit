"""Reader for the JSON orbital export written by ORCA's ``orca_2json`` utility.

``orca_2json`` is a standalone ORCA utility rather than a calculation (manual
§7.58): it rewrites the contents of a ``.gbw`` orbital file into JSON so that
other programs can read the orbitals, their occupations and energies, the AO
overlap matrix and the molecule metadata. The local-basis route of this toolkit
needs exactly that -- an entropy spectrum must be read in a localized orbital
basis, and the two ends of that transformation are two exports of the same
orbital file, one before and one after ``orca_loc``.

What this reader returns
------------------------

``OrcaJson``, with the conventions

- ``mo_coefficients[mo][ao]``: one row per molecular orbital, one column per
  atomic orbital, both in ORCA's own order (the AO order is the one listed in the
  export's ``OrbitalLabels`` block). The AO columns are kept verbatim: no
  reordering and no sign fixing is applied, because two exports of different
  geometries may permute or phase-flip their columns and that alignment is the
  consumer's decision, not the reader's;
- ``atoms``: the element labels of the ``Atoms`` block, in file order;
- ``ao_labels``: the positional ``OrbitalLabels`` block split into
  :class:`AoLabel` records (or ``None`` when the export omits it). The labels
  are the only atom/angular-momentum information the export carries, so the
  analyses that project onto one centre's l shell (atomic terms) or partition
  by centre (environment spin) read them from here;
- ``n_ao`` is the number of AO coefficients per MO; an ``S-Matrix`` present in the
  file must be square with that same dimension, otherwise the export is refused
  (an overlap matrix from another system or basis set cannot be combined with
  these coefficients);
- ``overlap`` is ``None`` when the export has no ``S-Matrix`` block, which is not
  an error: the blocks written into the JSON depend on the ``orca.json.conf``
  configuration and on the requested options. The same holds for ``kinetic``
  (the ``T-Matrix`` block): an export is asked for the one-electron matrices
  with a configuration such as ``{"MOCoefficients": true, "1elIntegrals":
  ["H", "S", "T", "V"]}`` written next to the ``.gbw`` as ``<basename>.json.conf``
  (measured on ORCA 6.1.1: the run picks the file up and writes ``S-Matrix`` /
  ``H-Matrix`` / ``T-Matrix`` / ``V-Matrix`` into ``Molecule``). The kinetic
  matrix is what the cross-structure orbital mapping (``analysis/
  orbital_mapping``) needs: per-orbital kinetic energies are ``c_i^T T c_i``.

The overlap matrix is read from ``Molecule.S-Matrix``, the location measured in
the fixtures below (ORCA 6.1.1). The manual page for ``orca_2json`` additionally
documents the overlap as an optional entry of the requested-integrals block; an
export carrying it only there is not covered yet.

The Fock blocks and the windowed MO integrals (menus 21)
--------------------------------------------------------

Two further request families are read here because the APC orbital ranking
(``analysis/apc``) consumes them:

- ``"FockMatrix"`` entries ``"J"``, ``"K"`` and ``"F"`` become ``coulomb`` /
  ``exchange`` / ``coulomb_exchange``.  Unlike the one-electron blocks, these
  carry an outer spin dimension (manual: ``N_spin x N_AO x N_AO``; one block for
  the RHF exports measured here).  Measured conventions on ORCA 6.1.1
  (``fixtures/orca/n2_apc.json``, ``h2_apc.json``): the exported ``F`` block is
  exactly ``J + K`` and **excludes the core Hamiltonian**, so the full Fock in
  the MO basis is ``diag(C (H + J + K) C^T) = eps`` (matched to 6e-10); and the
  exported ``K`` block is the exchange contribution *as it enters the Fock*
  (prefactors included), i.e. in the MO basis
  ``-diag(C K C^T)[a] = sum_i (a i | a i)`` over the doubly occupied ``i``
  (matched to 4e-12 against the ``MO_IAJB`` entries).  The paper the ranking
  follows writes its ``0.5 K_aa`` for exactly that sum, so the two conventions
  agree up to the sign of this block -- see ``analysis/apc``.
- ``"2elIntegrals": ["MO_IAJB"]`` with an ``"OrbWin"`` becomes ``mo_iajb``
  (entries ``(i, j, a, b, value)`` with ``value = (i a | j b)`` in chemist
  notation, internal indices first) and the recorded ``mo_iajb_window``.
  Measured: the *input* ``OrbWin`` must be written as eight integers (the
  second window all zeros for a single window) -- the four-integer form is
  rejected by ``orca_2json`` with "Something is wrong with the orbital
  windows!", while the *recorded* window is the four-integer, inclusive form;
  entries below ``Thresh`` are simply absent.  The ``MO_IAJB`` block is keyed
  by spin pair; the RHF exports measured here carry ``alpha/alpha`` only, which
  is what this reader takes.
- ``"Densities": ["all"]`` becomes ``densities``: the named density matrices of
  the run's ``<base>.densities`` sidecar, one square AO matrix per entry.  The
  measured names include ``scfp`` (the SCF/reference density) and, for a
  ``FIC-NEVPT2`` run with ``%casscf ... PTSettings: Density Unrelaxed``,
  ``Tdens-CASNEV.mult.<M>.root.<R>.p`` per state -- the NEVPT2 unrelaxed
  density (``Tdens-CAS.mult.<M>.root.<R>.p`` is its CASSCF reference).  The
  measured convention is ``D_AO = S D_json S`` (validated against the known
  CASSCF density: recovered exactly diagonal in its natural-orbital basis with
  the printed occupations, block couplings 1e-15; and the recovered NEVPT2
  density's eigenvalues reproduce ORCA's printed "Natural Orbital Occupation
  Numbers" to the print precision, 5e-9).  Measured companion fact: the printed
  NOON list is the *whole-space* naturalization of the unrelaxed density, not
  a per-orbital diagonal; consumers that need the source's block-wise
  construction diagonalize the blocks themselves (``analysis/ass1st``).

Measured caveat for the localized export
----------------------------------------

On the fixture pair (``fixtures/orca/n2_fcidump.canonical.json`` and
``...localized.json``: the same N2/def2-SVP CAS(6,6) orbitals before and after an
``orca_loc`` IAO-IBO localization of orbitals 4-9, 0-based) the ``Occupancy`` and
``OrbitalEnergy`` fields are identical in all 28 MOs, while the coefficients
differ only inside the localized range. Those two fields are therefore carried
over from the source orbitals rather than recomputed for the rotated ones, and
must not be read as the occupations of the localized orbitals.

Why this reader is not registered with ``parse_auto``
-----------------------------------------------------

``parsers/base.py`` registers one parser per *program output*: ``detect``
recognises a format and ``parse`` returns the ``ParseResult`` whose sections feed
the diagnosis rules. An ``orca_2json`` file is neither. It is auxiliary data
produced by a post-processing step (``orca_loc``, a restart, an integral export)
from a ``.gbw`` that may belong to any calculation; it reports no outcome that a
rule could judge; and its JSON mentions ORCA without being an ORCA output, so
registering it would let a data file be accepted where the user asked for a run.
Callers import ``parse_orca_json`` from this module explicitly, exactly as the
file it reads is produced explicitly.

Error discipline is the parser layer's: every failure raises ``ParserError`` with
a "Next step: " hint, and a partly-filled object is never returned.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

from .base import ParserError, read_text

__all__ = ["ANGULAR_LETTERS", "AoLabel", "OrcaJson", "parse_orca_json"]

#: Angular-momentum letter per l, as it appears in ORCA's AO labels (measured:
#: s, p, d for the light shells and f for the lanthanide exports; g and h follow
#: the same series and are listed for completeness).
ANGULAR_LETTERS = {0: "s", 1: "p", 2: "d", 3: "f", 4: "g", 5: "h"}

# Every failure names a concrete next action; the call form is the manual's
# (``orca_2json <basename>.gbw``, the default output being ``<basename>.json``).
_EXPORT_HINT = (
    "the file must be an orbital export: run ``orca_2json`` on the ``.gbw`` the "
    "orbitals came from (``orca_2json <basename>.gbw``, default output "
    "``<basename>.json``) and pass that file here."
)

# U+FEFF, a byte-order mark as written by some Windows tools; json.loads rejects
# it, so it is removed before parsing rather than reported as malformed JSON.
_BOM = chr(0xFEFF)

# Next action for a dimension mismatch: the two blocks must come from one export.
_SAME_SOURCE_HINT = (
    "re-export both blocks from the same ``.gbw`` (``orca_2json <basename>.gbw``): an "
    "overlap matrix from a different system, basis set or export cannot be combined "
    "with these coefficients."
)


@dataclass(frozen=True)
class AoLabel:
    """One AO column of the export, split from ORCA's own label (``raw``).

    ORCA writes the ``OrbitalLabels`` entries as ``<atom index><element>
    <shell><letter><component>`` with the atom index 0-based and the shell
    counter 1-based within each (centre, angular momentum) pair -- measured on
    the fixtures, e.g. ``"0N   1pz"`` (first p shell of N0, component ``z``),
    ``"0Eu  1f+1"`` (first f shell of Eu0, component ``+1``).  The component
    spelling is the reader's business to interpret: p/d labels use the real
    solid-harmonic names (``z``, ``x``, ``y``; ``z2``, ``xz``, ``yz``,
    ``x2y2``, ``xy``), while f and higher labels use the m_l-like spelling
    (``0``, ``+1``, ``-1``, ... -- measured on the Eu export).
    """

    raw: str
    center: int
    element: str
    shell: int
    angular: str
    component: str


@dataclass(frozen=True)
class OrcaJson:
    """One ``orca_2json`` orbital export (auxiliary data, not a program output).

    ``mo_coefficients[mo][ao]`` and the matrices are stored as tuples of tuples
    so that the object is immutable and hashable like every other model of this
    package.  ``overlap``, ``kinetic``, ``hamiltonian`` are ``None`` when the
    export has no such block; ``coulomb`` / ``exchange`` / ``coulomb_exchange``
    are the ``FockMatrix`` spin blocks (outer index = spin, one for the RHF
    exports measured here); ``mo_iajb`` carries the requested ``MO_IAJB``
    entries and ``mo_iajb_window`` the window record from the export itself.
    """

    base_name: str
    charge: int
    multiplicity: int
    hftyp: str
    point_group: str
    atoms: tuple[str, ...]
    coordinates: tuple[tuple[float, float, float], ...]
    n_mo: int
    n_ao: int
    mo_coefficients: tuple[tuple[float, ...], ...]
    mo_occupations: tuple[float, ...]
    mo_energies: tuple[float, ...]
    overlap: tuple[tuple[float, ...], ...] | None = None
    kinetic: tuple[tuple[float, ...], ...] | None = None
    ao_labels: tuple[AoLabel, ...] | None = None
    hamiltonian: tuple[tuple[float, ...], ...] | None = None
    coulomb: tuple[tuple[tuple[float, ...], ...], ...] | None = None
    exchange: tuple[tuple[tuple[float, ...], ...], ...] | None = None
    coulomb_exchange: tuple[tuple[tuple[float, ...], ...], ...] | None = None
    mo_iajb: tuple[tuple[int, int, int, int, float], ...] | None = None
    mo_iajb_window: tuple[int, int, int, int] | None = None
    densities: tuple[tuple[str, tuple[tuple[float, ...], ...]], ...] | None = None


def parse_orca_json(path: str | Path) -> OrcaJson:
    """Read an ``orca_2json`` export; every defect raises ParserError with a next step."""
    text = read_text(path)
    try:
        document = json.loads(text.removeprefix(_BOM))
    except json.JSONDecodeError as exc:
        raise _fail(f"{path} is not valid JSON: {exc}.", _EXPORT_HINT) from exc
    if not isinstance(document, dict):
        raise _fail(
            f"{path} has a JSON {type(document).__name__} at the top level, not an object.",
            _EXPORT_HINT,
        )
    molecule = _mapping(_entry(document, "Molecule", "the export"), "the Molecule block")

    orbital_block = _mapping(
        _entry(molecule, "MolecularOrbitals", "the Molecule block"),
        "Molecule.MolecularOrbitals",
    )
    entries = _entry(orbital_block, "MOs", "Molecule.MolecularOrbitals")
    if not isinstance(entries, list):
        raise _fail(
            f"Molecule.MolecularOrbitals.MOs is a JSON {type(entries).__name__}, not a list.",
            _EXPORT_HINT,
        )
    if not entries:
        raise _fail(
            "Molecule.MolecularOrbitals.MOs is empty: the export carries no orbitals.",
            _EXPORT_HINT,
        )
    coefficients: list[tuple[float, ...]] = []
    occupations: list[float] = []
    energies: list[float] = []
    for index, entry in enumerate(entries):
        where = f"Molecule.MolecularOrbitals.MOs[{index}]"
        mo = _mapping(entry, where)
        coefficients.append(
            _number_row(_entry(mo, "MOCoefficients", where), f"{where}.MOCoefficients")
        )
        occupations.append(_number(_entry(mo, "Occupancy", where), f"{where}.Occupancy"))
        energies.append(_number(_entry(mo, "OrbitalEnergy", where), f"{where}.OrbitalEnergy"))
    n_ao = len(coefficients[0])
    if n_ao == 0:
        raise _fail(
            "Molecule.MolecularOrbitals.MOs[0].MOCoefficients is empty: the export carries "
            "no AO coefficients.",
            _EXPORT_HINT,
        )
    for index, row in enumerate(coefficients):
        if len(row) != n_ao:
            raise _fail(
                f"Molecule.MolecularOrbitals.MOs[{index}].MOCoefficients has {len(row)} AO "
                f"coefficients while the first MO has {n_ao}: the MO block is not a "
                "rectangular matrix.",
                _EXPORT_HINT,
            )

    atoms, coordinates = _atom_labels(molecule)
    # Absent S-Matrix / T-Matrix / H-Matrix: documented export variants (the
    # requested blocks depend on the export configuration), not defects.
    overlap = None
    if "S-Matrix" in molecule:
        overlap = _square_matrix(molecule["S-Matrix"], "Molecule.S-Matrix", n_ao)
    kinetic = None
    if "T-Matrix" in molecule:
        kinetic = _square_matrix(molecule["T-Matrix"], "Molecule.T-Matrix", n_ao)
    hamiltonian = None
    if "H-Matrix" in molecule:
        hamiltonian = _square_matrix(molecule["H-Matrix"], "Molecule.H-Matrix", n_ao)
    # The density-derived Fock blocks carry an outer spin dimension.
    coulomb = None
    if "J-Matrix" in molecule:
        coulomb = _spin_matrices(molecule["J-Matrix"], "Molecule.J-Matrix", n_ao)
    exchange = None
    if "K-Matrix" in molecule:
        exchange = _spin_matrices(molecule["K-Matrix"], "Molecule.K-Matrix", n_ao)
    coulomb_exchange = None
    if "F-Matrix" in molecule:
        coulomb_exchange = _spin_matrices(molecule["F-Matrix"], "Molecule.F-Matrix", n_ao)
    mo_iajb, mo_iajb_window = _mo_integrals(molecule, len(coefficients))
    densities = _densities(molecule, n_ao)
    # Absent OrbitalLabels: likewise a documented export variant; when present
    # the list is positional (AO column k <-> label k) and must fit n_ao.
    ao_labels = None
    if "OrbitalLabels" in orbital_block:
        ao_labels = _ao_labels(orbital_block["OrbitalLabels"], n_ao)

    return OrcaJson(
        base_name=_text(_entry(molecule, "BaseName", "the Molecule block"), "Molecule.BaseName"),
        charge=_integer(_entry(molecule, "Charge", "the Molecule block"), "Molecule.Charge"),
        multiplicity=_integer(
            _entry(molecule, "Multiplicity", "the Molecule block"), "Molecule.Multiplicity"
        ),
        hftyp=_text(_entry(molecule, "HFTyp", "the Molecule block"), "Molecule.HFTyp"),
        point_group=_text(
            _entry(molecule, "PointGroup", "the Molecule block"), "Molecule.PointGroup"
        ),
        atoms=atoms,
        coordinates=coordinates,
        n_mo=len(coefficients),
        n_ao=n_ao,
        mo_coefficients=tuple(coefficients),
        mo_occupations=tuple(occupations),
        mo_energies=tuple(energies),
        overlap=overlap,
        kinetic=kinetic,
        ao_labels=ao_labels,
        hamiltonian=hamiltonian,
        coulomb=coulomb,
        exchange=exchange,
        coulomb_exchange=coulomb_exchange,
        mo_iajb=mo_iajb,
        mo_iajb_window=mo_iajb_window,
        densities=densities,
    )


# --- value extraction -------------------------------------------------------
# The JSON is machine-written, so a wrong type means the file is not what the
# caller thinks it is; guessing a value here would hide that.


def _fail(detail: str, next_step: str) -> ParserError:
    """Build a parser error whose message carries one concrete next action."""
    return ParserError(f"{detail} Next step: {next_step}")


def _entry(block: Mapping[str, Any], key: str, where: str) -> Any:
    if key not in block:
        raise _fail(f"{where} has no {key!r} entry.", _EXPORT_HINT)
    return block[key]


def _mapping(value: Any, where: str) -> Mapping[str, Any]:
    if not isinstance(value, dict):
        raise _fail(f"{where} is a JSON {type(value).__name__}, not an object.", _EXPORT_HINT)
    return value


def _text(value: Any, where: str) -> str:
    if not isinstance(value, str):
        raise _fail(f"{where} is {value!r}, not a string.", _EXPORT_HINT)
    return value


def _integer(value: Any, where: str) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise _fail(f"{where} is {value!r}, not an integer.", _EXPORT_HINT)
    if isinstance(value, float) and not value.is_integer():
        raise _fail(
            f"{where} is {value!r}; truncating it to an integer would silently change it.",
            _EXPORT_HINT,
        )
    return int(value)


def _number(value: Any, where: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise _fail(f"{where} is {value!r}, not a number.", _EXPORT_HINT)
    return float(value)


def _number_row(value: Any, where: str) -> tuple[float, ...]:
    if not isinstance(value, list):
        raise _fail(
            f"{where} is a JSON {type(value).__name__}, not a list of numbers.", _EXPORT_HINT
        )
    return tuple(_number(item, f"{where}[{index}]") for index, item in enumerate(value))


def _atom_labels(
    molecule: Mapping[str, Any],
) -> tuple[tuple[str, ...], tuple[tuple[float, float, float], ...]]:
    """Element labels and coordinates of the ``Atoms`` block, in file order.

    The coordinates are in the export's own units (``CoordinateUnits``, Angstrom
    in every measured fixture); consumers that cut distances read them together
    with that field.
    """
    entries = _entry(molecule, "Atoms", "the Molecule block")
    if not isinstance(entries, list):
        raise _fail(
            f"the 'Atoms' entry is a JSON {type(entries).__name__}, not a list.", _EXPORT_HINT
        )
    if not entries:
        raise _fail("the 'Atoms' entry is empty: the export carries no atoms.", _EXPORT_HINT)
    labels: list[str] = []
    coordinates: list[tuple[float, float, float]] = []
    for index, entry in enumerate(entries):
        where = f"Molecule.Atoms[{index}]"
        atom = _mapping(entry, where)
        labels.append(_text(_entry(atom, "ElementLabel", where), f"{where}.ElementLabel"))
        coords = _entry(atom, "Coords", where)
        if not isinstance(coords, list) or len(coords) != 3:
            raise _fail(
                f"{where}.Coords is {coords!r}, not a list of three numbers.",
                _EXPORT_HINT,
            )
        coordinates.append(
            tuple(_number(value, f"{where}.Coords[{k}]") for k, value in enumerate(coords))
        )
    return tuple(labels), tuple(coordinates)  # type: ignore[return-value]


# ORCA's AO label grammar, measured on the fixtures (N2 and Eu exports of ORCA
# 6.1.1): "<atom index><element>  <shell><letter><component>".  The element is
# one or two letters directly after the digits of the atom index (``0N``,
# ``0Eu``), which is why the two are not separated by whitespace.
_AO_LABEL_RE = re.compile(
    r"^(?P<center>\d+)(?P<element>[A-Za-z]{1,2})\s+"
    r"(?P<shell>\d+)(?P<angular>[spdfghi])(?P<component>[A-Za-z0-9+-]*)$"
)


def _ao_labels(value: Any, dimension: int) -> tuple[AoLabel, ...]:
    """Parse the positional ``OrbitalLabels`` block against the AO dimension."""
    if not isinstance(value, list):
        raise _fail(
            f"Molecule.MolecularOrbitals.OrbitalLabels is a JSON {type(value).__name__}, "
            "not a list.",
            _EXPORT_HINT,
        )
    if len(value) != dimension:
        raise _fail(
            f"Molecule.MolecularOrbitals.OrbitalLabels has {len(value)} entries while the "
            f"MO coefficients have {dimension} AO columns: the labels do not belong to "
            "these coefficients.",
            _SAME_SOURCE_HINT,
        )
    labels: list[AoLabel] = []
    for index, entry in enumerate(value):
        where = f"Molecule.MolecularOrbitals.OrbitalLabels[{index}]"
        raw = _text(entry, where)
        match = _AO_LABEL_RE.match(raw.strip())
        if match is None:
            raise _fail(
                f"{where} is {raw!r}, which does not split into the measured ORCA grammar "
                "'<atom index><element> <shell><angular letter><component>'.",
                "the AO-labelling consumers (atomic-term and environment-partition "
                "analyses) need the label grammar; report the label so the reader can be "
                "extended.",
            )
        labels.append(
            AoLabel(
                raw=raw,
                center=int(match.group("center")),
                element=match.group("element"),
                shell=int(match.group("shell")),
                angular=match.group("angular"),
                component=match.group("component"),
            )
        )
    return tuple(labels)


def _square_matrix(value: Any, where: str, dimension: int) -> tuple[tuple[float, ...], ...]:
    """Read a square matrix that must match the AO dimension of the MO coefficients."""
    if not isinstance(value, list):
        raise _fail(f"{where} is a JSON {type(value).__name__}, not a matrix.", _EXPORT_HINT)
    if len(value) != dimension:
        raise _fail(
            f"{where} has {len(value)} rows while the MO coefficients have {dimension} AO "
            "columns: the two blocks do not belong to the same basis set.",
            _SAME_SOURCE_HINT,
        )
    rows = tuple(_number_row(row, f"{where}[{index}]") for index, row in enumerate(value))
    for index, row in enumerate(rows):
        if len(row) != dimension:
            raise _fail(
                f"{where}[{index}] has {len(row)} entries while {dimension} are expected: the "
                "matrix is not square.",
                _SAME_SOURCE_HINT,
            )
    return rows


def _spin_matrices(
    value: Any, where: str, dimension: int
) -> tuple[tuple[tuple[float, ...], ...], ...]:
    """Read a Fock-family block: a list of square matrices, one per spin."""
    if not isinstance(value, list):
        raise _fail(f"{where} is a JSON {type(value).__name__}, not a list of matrices.", _EXPORT_HINT)
    if not value:
        raise _fail(
            f"{where} is empty: a Fock-family block carries at least one spin matrix.",
            _EXPORT_HINT,
        )
    return tuple(
        _square_matrix(block, f"{where}[{index}]", dimension)
        for index, block in enumerate(value)
    )


def _mo_integrals(
    molecule: Mapping[str, Any], n_mo: int
) -> tuple[tuple[tuple[int, int, int, int, float], ...] | None, tuple[int, int, int, int] | None]:
    """Read the requested ``MO_IAJB`` window entries and the recorded window.

    Absent ``2elIntegrals`` block: a documented export variant (the request
    depends on the configuration), like the absent one-electron matrices.  A
    block present without the measured ``MO_IAJB`` member or without the
    measured ``alpha/alpha`` spin key is not guessed at: the reader reports the
    key so it can be extended once such an export is measured.
    """
    block = molecule.get("2elIntegrals")
    if block is None:
        return None, None
    block = _mapping(block, "Molecule.2elIntegrals")
    if "MO_IAJB" not in block:
        return None, None
    window_raw = _entry(block, "OrbWin", "Molecule.2elIntegrals")
    window = _int_row(window_raw, "Molecule.2elIntegrals.OrbWin")
    if len(window) != 4:
        raise _fail(
            f"Molecule.2elIntegrals.OrbWin is {window_raw!r}; the exports measured here "
            "record the window as four inclusive integers (first/last internal, first/last "
            "external).",
            "report the window record so the reader can be extended.",
        )
    iajb = _mapping(block["MO_IAJB"], "Molecule.2elIntegrals.MO_IAJB")
    if "alpha/alpha" not in iajb:
        raise _fail(
            f"Molecule.2elIntegrals.MO_IAJB carries the spin keys {sorted(iajb)!r}; the "
            "RHF exports measured here carry 'alpha/alpha' only.",
            "report the spin key so the reader can be extended.",
        )
    entries = iajb["alpha/alpha"]
    if not isinstance(entries, list):
        raise _fail(
            f"Molecule.2elIntegrals.MO_IAJB['alpha/alpha'] is a JSON "
            f"{type(entries).__name__}, not a list of entries.",
            _EXPORT_HINT,
        )
    parsed: list[tuple[int, int, int, int, float]] = []
    for index, entry in enumerate(entries):
        where = f"Molecule.2elIntegrals.MO_IAJB['alpha/alpha'][{index}]"
        if not isinstance(entry, list) or len(entry) != 5:
            raise _fail(
                f"{where} is {entry!r}, not a five-item [i, j, a, b, value] entry.",
                "report the entry shape so the reader can be extended.",
            )
        i, j, a, b = (_integer(item, f"{where}[{k}]") for k, item in enumerate(entry[:4]))
        for label, orbital in (("i", i), ("j", j), ("a", a), ("b", b)):
            if not 0 <= orbital < n_mo:
                raise _fail(
                    f"{where} names orbital index {orbital} for {label} while the export "
                    f"carries {n_mo} molecular orbitals.",
                    _SAME_SOURCE_HINT,
                )
        parsed.append((i, j, a, b, _number(entry[4], f"{where}[4]")))
    return tuple(parsed), window


def _int_row(value: Any, where: str) -> tuple[int, ...]:
    """Read a list of integers (the window record)."""
    if not isinstance(value, list):
        raise _fail(f"{where} is a JSON {type(value).__name__}, not a list.", _EXPORT_HINT)
    return tuple(_integer(item, f"{where}[{index}]") for index, item in enumerate(value))


def _densities(
    molecule: Mapping[str, Any], dimension: int
) -> tuple[tuple[str, tuple[tuple[float, ...], ...]], ...] | None:
    """Read the named AO density matrices of a ``Densities`` request.

    Absent block: a documented export variant, like the other request-dependent
    blocks.  Present: a JSON object title -> square matrix; a matrix that is not
    square with the AO dimension is refused (it cannot pair with these
    orbitals), and a title whose value is not a plain matrix (an unmeasured
    per-spin shape) is reported rather than guessed at.
    """
    block = molecule.get("Densities")
    if block is None:
        return None
    block = _mapping(block, "Molecule.Densities")
    entries: list[tuple[str, tuple[tuple[float, ...], ...]]] = []
    for title, value in block.items():
        where = f"Molecule.Densities[{title!r}]"
        if isinstance(value, dict):
            raise _fail(
                f"{where} is a JSON object with keys {sorted(value)!r}, not a plain "
                "matrix; the exports measured here carry one square matrix per title.",
                "report the density shape so the reader can be extended.",
            )
        entries.append((title, _square_matrix(value, where, dimension)))
    return tuple(entries)
