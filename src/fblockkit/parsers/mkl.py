"""Reader and writer for ORCA's Molekel MKL file -- the gbw write-back route.

``orca_2mkl <base> -mkl`` rewrites a ``.gbw`` orbital file into the Molekel MKL
format, which is **plain text**: the geometry (in Angstrom), the basis-set
block, the MO coefficients with their energies, and the occupations.  Editing
that text and converting back with ``orca_2mkl <new> -gbw`` produces a ``.gbw``
that ORCA accepts -- measured on 6.1.1 (2026-09-27, the N2 scan fixtures): a
gbw -> mkl -> gbw round trip preserves the coefficients to 5e-8 (the mkl's own
print precision), editing the ``$COORD`` block or swapping orbitals carries
through, and a cross-geometry file is accepted as ``INITIAL GUESS: MOREAD``.
That is what makes orbital write-back possible at all: the JSON route
(``orca_2json``) has no reverse direction, this one does.

What this reader returns
------------------------

``MklFile``: the charge/multiplicity pair, the atoms (nuclear charge + position
in Angstrom), the per-atom charges, the ``$BASIS`` block **verbatim** (its
normalisation belongs to ``orca_2mkl`` and is not this module's business), the
coefficient groups, and the occupations.  The coefficient block is stored the
way the file lays it out -- groups of up to five orbitals, each group one
label line, one line of orbital energies and then one row per AO -- with the
:meth:`MklFile.coefficients` view turning it into an ``n_ao x n_mo`` matrix.
An unrestricted file carries the beta pair (``$COEFF_BETA`` + ``$OCC_BETA``,
measured on a UKS reference where orca_2mkl writes both); it is optional and
mirrors the alpha block's group structure.

The writer rebuilds the whole file from those pieces; whitespace is
normalised (the format is token-based -- ``orca_2mkl`` reads numbers, not
columns) while every number keeps the file's seven decimals.  ``$BASIS`` and
``$CHAR_MULT`` are carried over untouched; the charges come along as stored.

Measured row-order convention (why ``coefficients`` permutes)
-------------------------------------------------------------

The mkl's AO *rows* order the components of a p shell differently from the
``.gbw`` (and therefore from an ``orca_2json`` export): on the N2 fixtures the
mkl rows of each p shell are (py, pz, px) where the export lists (pz, px, py)
-- measured by matching the same canonical orbitals element by element
(2026-09-27).  s and d shells were measured to agree (identity, no sign
flips); f and higher are not covered by the fixtures and are carried as
identity, marked unverified.  :meth:`MklFile.coefficients` therefore returns
rows in the **export/``.gbw`` convention** (so that ``c^T S c = I`` against an
export's overlap matrix) and :meth:`MklFile.with_orbitals` accepts that same
convention and maps back into the file's own layout.  The shell layout is read
from the ``$BASIS`` headers (``<components> <angular letter> <scale>``, one
primitive per following line, ``$$`` separating the atoms' blocks).

Error discipline is the parser layer's: every failure raises ``ParserError``
with a "Next step: " hint.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path

from .base import ParserError, read_text

__all__ = ["MklError", "MklGroup", "MklFile", "parse_mkl"]

_MKL_HINT = (
    "the file must be a Molekel MKL export of an orbital file: run "
    "``orca_2mkl <basename> -mkl`` on the ``.gbw`` and pass the produced "
    "``<basename>.mkl`` here."
)

#: Measured component reordering inside a p shell: the mkl row for the export's
#: component ``k`` (0-based within the shell) is ``base + _P_ROW_ORDER[k]``.
_P_ROW_ORDER = (2, 0, 1)


class MklError(ParserError):
    """The MKL file cannot be read or written (with a next step)."""


def _fnum(value: float, decimals: int) -> str:
    """Fixed-point formatting with a print-resolution snap against signed zeros.

    A coefficient that is numerically zero but carries BLAS noise of either
    sign prints as ``-0.0000000`` on one build and ``0.0000000`` on another --
    the minus consumes a column and the line differs even though the value is
    zero.  Snapping below half of the last printed digit removes the signed
    zero without touching any printed digit (the same snap is what makes the
    examples' captures comparable across machines).
    """
    if abs(value) < 0.5 * 10.0 ** (-decimals):
        value = 0.0
    return f"{value:.{decimals}f}"


@dataclass(frozen=True)
class MklGroup:
    """One group of up to five orbitals, laid out the way the file stores them."""

    labels: tuple[str, ...]
    energies: tuple[float, ...]
    rows: tuple[tuple[float, ...], ...]  # one row per AO, one value per label


def _shell_components(basis_lines) -> tuple[tuple[int, str], ...]:
    """The shells in file order, as ``(components, angular letter)`` pairs.

    Reads the ``$BASIS`` headers -- ``<components> <angular letter> <scale>``
    followed by one primitive line per primitive -- and skips the ``$$`` atom
    separators (they carry no component information; the AO order is the shell
    order regardless of atom boundaries).
    """
    shells: list[tuple[int, str]] = []
    index = 0
    lines = [line for line in basis_lines]
    while index < len(lines):
        stripped = lines[index].strip()
        index += 1
        if not stripped or stripped == "$$":
            continue
        tokens = stripped.split()
        if len(tokens) != 3:
            raise _fail(
                f"a $BASIS header line is {stripped!r}; the measured grammar is "
                "'<components> <angular letter> <scale>'."
            )
        try:
            components = int(tokens[0])
        except ValueError as exc:
            raise _fail(f"a $BASIS header line is {stripped!r}, not '<n> <L> <scale>'.") from exc
        letter = tokens[1].upper()
        if components != {"S": 1, "P": 3, "D": 5, "F": 7, "G": 9}.get(letter, components):
            raise _fail(
                f"the $BASIS shell {stripped!r} claims {components} component(s) for "
                f"angular momentum {letter!r}, which disagrees with the 2l+1 rule."
            )
        shells.append((components, letter))
        # skip the shell's primitive lines: one per primitive
        # (measured: the header's first field is NOT the primitive count, so count
        # by consuming lines that do not look like a header)
        while index < len(lines):
            candidate = lines[index].strip()
            if not candidate:
                index += 1
                continue
            if candidate == "$$":
                break
            tokens_next = candidate.split()
            if len(tokens_next) == 3 and tokens_next[1].isalpha():
                break
            index += 1
        if index < len(lines) and lines[index].strip() == "$$":
            index += 1
    return tuple(shells)


def _row_order(basis_lines) -> tuple[int, ...]:
    """For each export-convention AO row, the mkl row that stores it."""
    order: list[int] = []
    base = 0
    for components, letter in _shell_components(basis_lines):
        if letter == "P":
            for component in range(3):
                order.append(base + _P_ROW_ORDER[component])
        else:
            for component in range(components):
                order.append(base + component)
        base += components
    return tuple(order)


@dataclass(frozen=True)
class MklFile:
    """A read Molekel MKL file: geometry, basis (verbatim), orbitals, occupations.

    Unrestricted runs carry a second coefficient/occupation pair; the beta
    block is optional (measured: ORCA writes ``$COEFF_BETA``/``$OCC_BETA``
    only for unrestricted references) and is carried through the reader and
    the writer when present.
    """

    charge_mult: tuple[int, int]
    atoms: tuple[tuple[int, float, float, float], ...]  # Z, x, y, z in Angstrom
    charges: tuple[float, ...]
    basis_lines: tuple[str, ...]
    groups: tuple[MklGroup, ...]
    occupations: tuple[float, ...]
    header_comments: tuple[str, ...] = ()
    beta_groups: tuple[MklGroup, ...] = ()
    beta_occupations: tuple[float, ...] = ()

    @property
    def unrestricted(self) -> bool:
        """Whether the file carries a beta orbital block (an unrestricted run)."""
        return bool(self.beta_groups)

    @property
    def n_ao(self) -> int:
        return len(self.groups[0].rows) if self.groups else 0

    @property
    def n_mo(self) -> int:
        return sum(len(group.labels) for group in self.groups)

    def _stored_matrix(self, groups: tuple[MklGroup, ...] | None = None) -> list[list[float]]:
        """The coefficient matrix exactly as stored (mkl row order)."""
        matrix = [[] for _ in range(self.n_ao)]
        for group in groups if groups is not None else self.groups:
            for row_index, row in enumerate(group.rows):
                matrix[row_index].extend(row)
        return matrix

    def _to_export_order(self, stored: list[list[float]]) -> tuple[tuple[float, ...], ...]:
        order = _row_order(self.basis_lines)
        if len(order) != self.n_ao:
            raise _fail(
                f"the $BASIS shells account for {len(order)} AO rows while the "
                f"coefficient block holds {self.n_ao}."
            )
        return tuple(tuple(stored[mkl_row]) for mkl_row in order)

    def coefficients(self) -> tuple[tuple[float, ...], ...]:
        """The alpha ``n_ao x n_mo`` coefficient matrix in export/``.gbw`` row order."""
        return self._to_export_order(self._stored_matrix())

    def beta_coefficients(self) -> tuple[tuple[float, ...], ...]:
        """The beta coefficient matrix (empty tuple for a restricted file)."""
        if not self.beta_groups:
            return ()
        return self._to_export_order(self._stored_matrix(self.beta_groups))

    def _regroup(self, groups: tuple[MklGroup, ...], matrix, energies, occupations):
        if len(matrix) != self.n_ao or any(len(row) != self.n_mo for row in matrix):
            raise MklError(
                f"the new coefficient matrix is {len(matrix)} x "
                f"{len(matrix[0]) if matrix else 0} while the template holds "
                f"{self.n_ao} AOs x {self.n_mo} MOs. Next step: interpolate within "
                "one and the same basis set (the template's)."
            )
        if energies is not None and len(energies) != self.n_mo:
            raise MklError(
                f"the energy list has {len(energies)} entries while the template holds "
                f"{self.n_mo} orbitals. Next step: give one energy per orbital."
            )
        if occupations is not None and len(occupations) != self.n_mo:
            raise MklError(
                f"the occupation list has {len(occupations)} entries while the template "
                f"holds {self.n_mo} orbitals. Next step: give one occupation per orbital."
            )
        stored_energies = (
            tuple(energies) if energies is not None else
            tuple(value for group in groups for value in group.energies)
        )
        order = _row_order(self.basis_lines)
        if len(order) != self.n_ao:
            raise MklError(
                f"the $BASIS shells account for {len(order)} AO rows while the "
                f"template holds {self.n_ao}. Next step: check the template's basis block."
            )
        generated: list[list[float] | None] = [None] * self.n_ao
        for export_row, mkl_row in enumerate(order):
            generated[mkl_row] = list(matrix[export_row])
        rebuilt: list[MklGroup] = []
        offset = 0
        for group in groups:
            width = len(group.labels)
            block = [
                tuple(row[offset + column] for column in range(width))
                for row in generated  # type: ignore[union-attr]
            ]
            rebuilt.append(
                MklGroup(
                    labels=group.labels,
                    energies=tuple(stored_energies[offset : offset + width]),
                    rows=tuple(block),
                )
            )
            offset += width
        return tuple(rebuilt), (
            tuple(occupations) if occupations is not None else None
        )

    def with_orbitals(
        self,
        coefficients,
        occupations=None,
        energies=None,
        *,
        beta_coefficients=None,
        beta_occupations=None,
        beta_energies=None,
    ) -> "MklFile":
        """A copy carrying new coefficients (``n_ao x n_mo``), grouped the same way.

        ``energies`` (per orbital) and ``occupations`` default to the stored
        ones; the label lines are kept as they are (the labels are symmetry
        tags of the template).  The beta-prefixed arguments must be absent
        together with a beta block (a restricted template refuses a beta
        replacement).
        """
        groups, new_occupations = self._regroup(
            self.groups, coefficients, energies, occupations
        )
        if beta_coefficients is not None or beta_occupations is not None or beta_energies is not None:
            if not self.beta_groups:
                raise MklError(
                    "a beta orbital replacement was given for a template without a "
                    "$COEFF_BETA block (a restricted file). Next step: perturb the "
                    "alpha block alone, or use an unrestricted reference mkl."
                )
            beta_groups, new_beta_occupations = self._regroup(
                self.beta_groups,
                beta_coefficients if beta_coefficients is not None else self.beta_coefficients(),
                beta_energies,
                beta_occupations,
            )
        else:
            beta_groups, new_beta_occupations = self.beta_groups, None
        return replace(
            self,
            groups=groups,
            occupations=new_occupations if new_occupations is not None else self.occupations,
            beta_groups=beta_groups,
            beta_occupations=(
                new_beta_occupations if new_beta_occupations is not None
                else self.beta_occupations
            ),
        )

    def render(self) -> str:
        """Rebuild the MKL text (the form ``orca_2mkl <new> -mkl``-style reads)."""
        lines: list[str] = ["$MKL", "#", "# MKL format file produced by fBlockKit", "#"]
        lines.append("$CHAR_MULT")
        lines.append(f"  {self.charge_mult[0]} {self.charge_mult[1]}")
        lines.append("$END")
        lines.append("")
        lines.append("$COORD")
        for z, x, y, zz in self.atoms:
            lines.append(f"  {z}  {_fnum(x, 6)}  {_fnum(y, 6)}  {_fnum(zz, 6)}")
        lines.append("$END")
        lines.append("")
        lines.append("$CHARGES")
        for value in self.charges:
            lines.append(f"  {_fnum(value, 6)}")
        lines.append("$END")
        lines.append("")
        lines.append("$BASIS")
        lines.extend(self.basis_lines)
        lines.append("$END")
        lines.append("")
        lines.append("$COEFF_ALPHA")
        for group in self.groups:
            lines.append("  " + "  ".join(f"{label:>3}" for label in group.labels))
            lines.append("  " + "  ".join(f"{_fnum(value, 7):>13}" for value in group.energies))
            for row in group.rows:
                lines.append("  " + "  ".join(f"{_fnum(value, 7):>13}" for value in row))
        # the measured file terminates the coefficient block with its own END
        # before the occupations; dropping it made orca_2mkl hang on the read
        lines.append(" $END")
        lines.append("")
        lines.append("$OCC_ALPHA")
        values = [f"{_fnum(value, 7):>13}" for value in self.occupations]
        for start in range(0, len(values), 5):
            lines.append("  " + "  ".join(values[start : start + 5]))
        lines.append(" $END")
        if self.beta_groups:
            # the measured unrestricted layout: the beta pair mirrors the alpha
            # pair (COEFF then OCC), exactly as orca_2mkl writes it
            lines.append("")
            lines.append("$COEFF_BETA")
            for group in self.beta_groups:
                lines.append("  " + "  ".join(f"{label:>3}" for label in group.labels))
                lines.append("  " + "  ".join(f"{_fnum(value, 7):>13}" for value in group.energies))
                for row in group.rows:
                    lines.append("  " + "  ".join(f"{_fnum(value, 7):>13}" for value in row))
            lines.append(" $END")
            lines.append("")
            lines.append("$OCC_BETA")
            beta_values = [f"{_fnum(value, 7):>13}" for value in self.beta_occupations]
            for start in range(0, len(beta_values), 5):
                lines.append("  " + "  ".join(beta_values[start : start + 5]))
            lines.append(" $END")
        return "\n".join(lines) + "\n"


def _is_number(line: str) -> bool:
    tokens = line.split()
    if not tokens:
        return False
    try:
        for token in tokens:
            float(token)
    except ValueError:
        return False
    return True


def _fail(detail: str, next_step: str = _MKL_HINT) -> MklError:
    return MklError(f"{detail} Next step: {next_step}")


def parse_mkl(path: str | Path) -> MklFile:
    """Read a Molekel MKL file; every defect raises ``MklError`` with a next step."""
    text = read_text(path)
    lines = text.splitlines()
    sections: dict[str, list[str]] = {}
    order: list[str] = []
    current: str | None = None
    for line in lines:
        stripped = line.strip()
        # Measured grammar: a section starts with '$' + an uppercase letter
        # ('$END', '$COORD', ...); the '$$' line inside $BASIS is content.
        if len(stripped) >= 2 and stripped[0] == "$" and stripped[1].isalpha() and stripped[1].isupper():
            name = stripped[1:].split()[0]
            if name == "END":
                current = None
                continue
            current = name
            order.append(name)
            sections[name] = []
            continue
        if current is not None:
            sections[current].append(line)
    if "COORD" not in sections or "COEFF_ALPHA" not in sections:
        raise _fail(
            f"{path} does not carry the measured MKL sections (found "
            f"{', '.join(order) or 'none'})."
        )
    charge_mult = _charge_mult(sections.get("CHAR_MULT", []))
    atoms = _atoms(sections["COORD"])
    charges = _floats(sections.get("CHARGES", []), "the $CHARGES block")
    groups = _coefficient_groups(sections["COEFF_ALPHA"])
    occupations = _floats(sections.get("OCC_ALPHA", []), "the $OCC_ALPHA block")
    n_mo = sum(len(group.labels) for group in groups)
    if len(occupations) != n_mo:
        raise _fail(
            f"$OCC_ALPHA carries {len(occupations)} values while the coefficient block "
            f"holds {n_mo} orbitals."
        )
    n_ao = len(groups[0].rows)
    for index, group in enumerate(groups):
        if len(group.rows) != n_ao:
            raise _fail(
                f"coefficient group {index} has {len(group.rows)} AO rows while group 0 "
                f"has {n_ao}: the block is not rectangular."
            )
        if len(group.energies) != len(group.labels):
            raise _fail(
                f"coefficient group {index} has {len(group.energies)} energies for "
                f"{len(group.labels)} labels."
            )
    beta_groups: tuple[MklGroup, ...] = ()
    beta_occupations: tuple[float, ...] = ()
    if "COEFF_BETA" in sections or "OCC_BETA" in sections:
        # measured: orca_2mkl writes both beta sections together for unrestricted
        # references; one without the other is a broken file
        if "COEFF_BETA" not in sections or "OCC_BETA" not in sections:
            raise _fail(
                "the file carries only one of $COEFF_BETA / $OCC_BETA; orca_2mkl "
                "writes the beta pair together for unrestricted references."
            )
        beta_groups = _coefficient_groups(sections["COEFF_BETA"])
        if len(beta_groups) != len(groups) or any(
            len(a.rows) != len(b.rows) or len(a.labels) != len(b.labels)
            for a, b in zip(groups, beta_groups)
        ):
            raise _fail(
                "the $COEFF_BETA block does not mirror the $COEFF_ALPHA block's "
                "group structure; the file is inconsistent."
            )
        beta_occupations = _floats(sections["OCC_BETA"], "the $OCC_BETA block")
        if len(beta_occupations) != n_mo:
            raise _fail(
                f"$OCC_BETA carries {len(beta_occupations)} values while the "
                f"coefficient block holds {n_mo} orbitals."
            )
    return MklFile(
        charge_mult=charge_mult,
        atoms=atoms,
        charges=charges,
        basis_lines=tuple(sections.get("BASIS", [])),
        groups=groups,
        occupations=occupations,
        beta_groups=beta_groups,
        beta_occupations=beta_occupations,
    )


def _charge_mult(lines: list[str]) -> tuple[int, int]:
    tokens = [token for line in lines for token in line.split()]
    if len(tokens) != 2:
        raise _fail(f"the $CHAR_MULT block is {lines!r}, not two integers.")
    try:
        return int(tokens[0]), int(tokens[1])
    except ValueError as exc:
        raise _fail(f"the $CHAR_MULT block is {lines!r}, not two integers.") from exc


def _atoms(lines: list[str]) -> tuple[tuple[int, float, float, float], ...]:
    atoms = []
    for index, line in enumerate(lines):
        tokens = line.split()
        if not tokens:
            continue
        if len(tokens) != 4:
            raise _fail(
                f"$COORD line {index} is {line!r}; four values (charge x y z) are expected."
            )
        try:
            atoms.append(
                (int(tokens[0]), float(tokens[1]), float(tokens[2]), float(tokens[3]))
            )
        except ValueError as exc:
            raise _fail(f"$COORD line {index} is {line!r}, not four numbers.") from exc
    if not atoms:
        raise _fail("the $COORD block is empty.")
    return tuple(atoms)


def _floats(lines: list[str], where: str) -> tuple[float, ...]:
    values = []
    for index, line in enumerate(lines):
        for token in line.split():
            try:
                values.append(float(token))
            except ValueError as exc:
                raise _fail(f"{where} line {index} carries {token!r}, not a number.") from exc
    return tuple(values)


def _coefficient_groups(lines: list[str]) -> tuple[MklGroup, ...]:
    """Parse the group layout: label line, energy line, then one row per AO."""
    rows_of = [(kind, line.split()) for line in lines if (kind := ("n" if _is_number(line) else "l")) and line.split()]
    groups: list[MklGroup] = []
    position = 0
    while position < len(rows_of):
        kind, tokens = rows_of[position]
        if kind != "l":
            raise _fail(
                f"the $COEFF_ALPHA block starts a group with numbers ({tokens!r}) "
                "instead of an orbital-label line."
            )
        labels = tuple(tokens)
        position += 1
        if position >= len(rows_of) or rows_of[position][0] != "n":
            raise _fail("a $COEFF_ALPHA group carries no orbital-energy line.")
        energies = tuple(float(token) for token in rows_of[position][1])
        position += 1
        if len(energies) != len(labels):
            raise _fail(
                f"a $COEFF_ALPHA group has {len(labels)} labels but "
                f"{len(energies)} energies on its energy line."
            )
        rows: list[tuple[float, ...]] = []
        while position < len(rows_of) and rows_of[position][0] == "n":
            values = tuple(float(token) for token in rows_of[position][1])
            if len(values) != len(labels):
                raise _fail(
                    f"a $COEFF_ALPHA row carries {len(values)} values for "
                    f"{len(labels)} orbitals."
                )
            rows.append(values)
            position += 1
        if not rows:
            raise _fail("a $COEFF_ALPHA group carries no AO rows.")
        groups.append(MklGroup(labels=labels, energies=energies, rows=tuple(rows)))
    if len({len(group.rows) for group in groups}) != 1:
        raise _fail("the $COEFF_ALPHA groups disagree about the AO count.")
    return tuple(groups)
