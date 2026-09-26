"""A6: environment spin polarisation -- ORCA's local spin analysis read as a surrogate
for the environment spin polarisation entropy ``Delta S_E`` of the reference paper.

Target quantity (reference paper, Eq. (9)): with the spin-resolved truncated environment
densities ``D_E^alpha`` / ``D_E^beta`` and their spin sum ``D_E = D_E^alpha + D_E^beta``,

    Delta S_E = -2 Tr[(D_E/2) ln(D_E/2)]
                + Tr[D_E^alpha ln D_E^alpha] + Tr[D_E^beta ln D_E^beta],

which is exactly zero when ``D_E^alpha = D_E^beta = D_E/2`` -- that is, when the
environment carries no spin polarisation and the whole spin polarisation sits on the
metal.  Measured scale in that paper (1Dy, a 6H Dy(III) single-ion magnet, 11 states):
the wrong SCF solution gives Delta S_E = 2.766, the correct one 0.007, and the wrong
solution degrades the downstream pre-SOC magnetic MAE from 8.6 to 357.7 cm^-1.

Why this module does not evaluate Eq. (9) itself: no printed ORCA table carries the
spin-resolved environment density.  The real quantity needs an off-line two-step route
(this group's manual check, 2026-09-26): ``orca_2json`` exports the spin-resolved
two-particle density matrices RDM2_aa / RDM2_ab / RDM2_bb (the manual states that the
2-RDMs are not stored inside the ``.densities`` file), and ``orca_loc`` produces the
IAO-IBO localised orbitals (``.loc.gbw``); with those one builds ``D_E^alpha`` and
``D_E^beta`` by hand.

What this module reads instead: the ready-made *Local Spin Analysis* block that ORCA
prints itself (header "LOCAL SPIN ANALYSIS (Loewdin* projector)"; ORCA 6.1 manual section
5.1.10; implemented in the SCF and CASSCF modules).  It answers the same kind of question
-- how the spin is distributed over fragments -- but it is a surrogate, and every report
body carries a boundary paragraph that says so.

Input (parsers layer, ``fblockkit.parsers.orca._parse_local_spin``)::

    sections["local_spin"] = {
        "present": bool,
        "blocks": tuple of {
            "state": int | None,        # None for the state-average block and for SCF jobs
            "state_label": str,         # "State average densities" / "0" / "" (SCF)
            "block_label": str,         # the "State belongs to block" value
            "multiplicity": int | None,
            "n_fragments": int, "n_atoms": int | None, "n_basis_functions": int | None,
            "sab": tuple of rows,       # symmetric <SA*SB> matrix
            "sz": tuple,                # <SzA>; None entries for the printed "n.a."
            "sz_na": bool,              # every <SzA> is "n.a." (a singlet state)
            "seff": tuple,              # Seff(A)
        },
        "fragment_elements": tuple of tuples,  # element symbol(s) per fragment, from the echo
    }

Block selection (measured on this project's fixtures; consistent with the manual's "if
nroots > 1, the printing will contain the state-specific analysis of all roots"):

- an SCF job prints several blocks (initial guess / after SCF / a final duplicate) -- only
  the LAST one belongs to the converged wavefunction, and only that one is used;
- a CASSCF job prints a state-average block (``state`` None, ``state_label`` starting with
  "State average") plus one block per root (``state`` 0, 1, ...); all of them are used,
  the state-average block first, the roots in ascending state order;
- the report says which blocks it used and lists the Seff values of the unused SCF blocks,
  so a reader can see that the choice of block changes the numbers.

Conventions that travel with the numbers: the partition depends on the fragment definition
(which atoms were assigned to which fragment) and on the projector named in the block
header (this project's fixtures use "Loewdin*").  Numbers obtained with a different
fragment definition or a different projector are not comparable.

Checks this module adds to the printed columns (each marked as its own, provisional):
the sum over all fragment pairs of ``<SA*SB>`` is compared with ``S(S+1)`` of the state --
``S`` from the block's multiplicity, ``S = (multiplicity - 1)/2`` -- and a deviation is
reported as spin contamination, which is a property of the method (a spin-unrestricted or
broken-symmetry single determinant is not a spin eigenstate), not an error in the file; the
diagonal ``<SA*SA>`` is compared with ``Seff(A)(Seff(A)+1)``; and the ``<SzA>`` values are
summed, since the local projections add up to the total ``M_S`` of the state.

The spin-share printed by this module (the metal row of <SA*SB>, summed over the other
fragments, divided by the total <S^2>) is **this module's own construction** -- it is not
taken from the cited sources, it is provisional, and the report body marks it as such.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Mapping, Sequence

from ..knowledge.elements import ElementError, is_f_element
from ..knowledge.models import (
    EVIDENCE_LITERATURE,
    EVIDENCE_MANUAL,
    EVIDENCE_MEASURED,
    Evidence,
    ParseResult,
    ReportSection,
)

#: Below this |sum over all fragment pairs of <SA*SB>| the state counts as a total singlet
#: and the spin-share is 0/0.  The printed values carry four decimals, so an exact zero
#: comes out exact; the tolerance only guards rounding.
SINGLET_TOLERANCE = 1e-6

#: Tolerance of the printed-number check ``<SA*SA> == Seff(A) * (Seff(A) + 1)``.  Seff is
#: printed to four decimals, which moves ``Seff(Seff+1)`` by about ``(2 Seff + 1) * 5e-5``
#: -- 1.8e-4 at Seff = 1.3, 7.5e-4 at Seff = 7 -- so the limit sits above the metal range.
SEFF_DIAGONAL_TOLERANCE = 3e-3

#: Tolerance on ``|<S^2> - S(S+1)|`` below which a block counts as a spin eigenstate
#: (spin-pure).  Both sides are built from numbers printed to four decimals, and an
#: off-diagonal term enters the sum twice, so the limit sits above the rounding.
SPIN_PURITY_TOLERANCE = 1e-3

#: Fragment numbering in the report follows the printed table (1-based); the helpers that
#: take an index use the 0-based row index of the ``<SA*SB>`` matrix.
SPIN_SHARE_FORMULA = "share(M) = (sum over all fragments B of <S_M*S_B>) / <S^2>_total"

#: Spins scanned by the "is this <S^2> an eigenvalue S(S+1)?" test when a block declares no
#: multiplicity; covers every spin up to S = 50, far beyond any f-block state.
_SPIN_SCAN = tuple(value / 2.0 for value in range(0, 101))


class LocalSpinError(ValueError):
    """Invalid or missing local spin analysis input."""


_EVIDENCE_MANUAL = Evidence(
    kind=EVIDENCE_MANUAL,
    text=(
        "ORCA 6.1 manual section 5.1.10 'Local Spin Analysis': Clark and Davidson's "
        "equations are \"implemented in the SCF and CASSCF modules of Orca\"; \"All that is "
        "required is to divide the molecule into fragments. The rest happens "
        "automatically.\" -- fragments are declared by tagging an atom's element symbol "
        "with a parenthesised fragment index in the coordinate block (the examples use "
        "N(1), N(2), S(3)). The printed tables are headed \"<SA*SB>\" and "
        "\"<SzA>  Seff(A)\"; \"if nroots > 1, the printing will contain the state-specific "
        "analysis of all roots\"; and a footnote reads \"for a singlet state all <SzA> "
        "values are zero by definition\". The block names its own methods (Clark & Davison "
        "2001; Herrmann, Reiher, Hess 2005), and the manual remarks that the total spin is "
        "not local but a property of the whole system."
    ),
    ref=(
        "ORCA 6.1 manual, section 5.1.10 'Local Spin Analysis' (read 2026-09-26): "
        "https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/population.html"
    ),
    url=(
        "https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/"
        "population.html"
    ),
)

_EVIDENCE_DSE = Evidence(
    kind=EVIDENCE_LITERATURE,
    text=(
        "Environment spin polarisation entropy Delta S_E, Eq. (9): Delta S_E = "
        "-2 Tr[(D_E/2) ln(D_E/2)] + Tr[D_E^alpha ln D_E^alpha] + Tr[D_E^beta ln "
        "D_E^beta], with D_E = D_E^alpha + D_E^beta the spin-summed environment density; it "
        "is exactly zero when the environment carries no spin polarisation (D_E^alpha = "
        "D_E^beta = D_E/2). Measured scale in this paper (1Dy, 6H, 11 states): the wrong "
        "SCF solution gives Delta S_E = 2.766 and the correct one 0.007, and the wrong "
        "solution raises the downstream pre-SOC magnetic MAE from 8.6 to 357.7 cm^-1. No "
        "printed ORCA table carries the spin-resolved environment density, so this module "
        "cannot evaluate Eq. (9) and reports the local spin analysis as a surrogate instead "
        "(see the boundary paragraph of every report body)."
    ),
    ref=(
        "Ai Y., Li Z.-W., Guan Z.-B., Jiang H., J. Chem. Theory Comput., 2025, 21(19), "
        "9631-9640, DOI 10.1021/acs.jctc.5c01336"
    ),
    bibkey="ai2025density",
    url="https://doi.org/10.1021/acs.jctc.5c01336",
)

_EVIDENCE_LOCAL_SPIN = Evidence(
    kind=EVIDENCE_LITERATURE,
    text=(
        "The framework behind the printed tables: the local spin of a molecule is analysed "
        "by partitioning it into fragments and reporting the expectation values of the "
        "fragment spin operators <SA*SB> and of the local projections <SzA> together with "
        "an effective local spin Seff(A). This paper is the comparative analysis of local "
        "spin definitions that the ORCA block itself cites (together with Clark & Davison, "
        "J. Chem. Phys. 2001, 115, 7382); which definition (projector) is used is stated in "
        "the block header."
    ),
    ref=(
        "Herrmann C., Reiher M., Hess B. A., J. Chem. Phys., 2005, 122(3), 034102, "
        "DOI 10.1063/1.1829050"
    ),
    bibkey="herrmann2005comparative",
    url="https://doi.org/10.1063/1.1829050",
)

_EVIDENCE_MEASURED = Evidence(
    kind=EVIDENCE_MEASURED,
    text=(
        "Measured on this project's fixtures (fixtures/orca/n2_stretch_local_spin.out and "
        "n2_stretch_casscf_local_spin.out): the UHF job prints three blocks and only the "
        "last one belongs to the converged wavefunction -- its <SA*SB> off-diagonal is "
        "-2.1835 where the initial-guess block gives +2.1641, so the choice of block changes "
        "the numbers. That final block has <SA*SB> = 3.6304 / -2.1835, <SzA> = +-1.4371 "
        "(summing to 0, the total M_S of the 14-electron state) and Seff = 1.4699 on both "
        "fragments; the sum over all fragment pairs is 2.8938, i.e. not a spin eigenvalue "
        "(S(S+1) would be 0 or 2 here), which is what spin contamination of a broken-"
        "symmetry single determinant looks like in this output -- the job was run with the "
        "broken-symmetry settings and the output carries the corresponding warning. The "
        "CASSCF blocks of the other fixture are spin-pure instead: their sums over all "
        "fragment pairs are 0 to the printed precision, matching S(S+1) = 0 of their "
        "declared singlet multiplicity. The CASSCF job prints a state-average block plus "
        "one block per "
        "root; both roots are singlets, so every <SzA> is the printed \"n.a.\" and the sum "
        "over all pairs is 0 -- the spin-share is then 0/0 and is refused, with the local "
        "spin magnitudes (Seff = 1.4110 and 0.5003) reported instead. In every fixture "
        "block the printed diagonal <SA*SA> equals Seff(A)(Seff(A)+1) to better than 2e-4; "
        "that relation is the check this module applies to the numbers it reads."
    ),
    ref="reproduced by tests/test_local_spin.py (2026-09-26 fixtures, ORCA 6.1.1)",
)


@dataclass(frozen=True)
class FragmentRow:
    """One fragment as printed in a local spin analysis block."""

    number: int                      # 1-based, as printed by ORCA
    elements: tuple[str, ...]        # from fragment_elements; empty when not recovered
    sz: float | None                 # None where the block prints "n.a." (a singlet state)
    seff: float | None
    f_block: bool                    # any known element is a lanthanide/actinide

    @property
    def element_text(self) -> str:
        return " ".join(self.elements) if self.elements else "(not recovered)"


def _matrix(block: Mapping[str, Any]) -> tuple[tuple[float, ...], ...]:
    """The block's ``<SA*SB>`` matrix, validated as square and finite."""
    sab = tuple(tuple(float(value) for value in row) for row in block.get("sab", ()))
    if not sab:
        raise LocalSpinError(
            "this local spin block carries an empty <SA*SB> matrix. Next step: check that "
            "the output was not truncated in the middle of the analysis, or parse the full "
            "file."
        )
    size = len(sab)
    if any(len(row) != size for row in sab):
        raise LocalSpinError(
            f"the <SA*SB> matrix is not square ({size} rows of differing length). "
            f"Next step: check the fragment tables of the output; each fragment needs one "
            f"row."
        )
    if not all(math.isfinite(value) for row in sab for value in row):
        raise LocalSpinError(
            "the <SA*SB> matrix contains non-finite values. Next step: check the printed "
            "table; a broken number suggests a truncated or mixed-up output file."
        )
    return sab


def _classify(symbol: str) -> bool:
    """Whether one element symbol is an f-block element; unknown symbols are not f-block.

    An unknown symbol is not an error: the tables stay readable, only the metal/ligand
    classification of that fragment becomes impossible, which the report body notes.
    """
    try:
        return is_f_element(symbol)
    except ElementError:
        return False


def unknown_symbols(fragment_elements: Sequence[Sequence[str]]) -> tuple[str, ...]:
    """Symbols of the fragment declaration that this toolkit's element table does not know."""
    unknown: list[str] = []
    for fragment in fragment_elements:
        for symbol in fragment:
            try:
                is_f_element(symbol)
            except ElementError:
                if symbol not in unknown:
                    unknown.append(symbol)
    return tuple(unknown)


def fragment_rows(
    block: Mapping[str, Any], fragment_elements: Sequence[Sequence[str]]
) -> tuple[FragmentRow, ...]:
    """The per-fragment rows of one block (``<SzA>`` and ``Seff(A)``, 1-based numbering).

    ``fragment_elements[i]`` supplies the element symbols of fragment ``i`` (0-based, as
    the matrix rows are); a shorter element list leaves the remaining fragments unnamed.
    """
    sab = _matrix(block)
    sz = tuple(block.get("sz", ()))
    seff = tuple(block.get("seff", ()))
    rows: list[FragmentRow] = []
    for index in range(len(sab)):
        elements = tuple(fragment_elements[index]) if index < len(fragment_elements) else ()
        rows.append(
            FragmentRow(
                number=index + 1,
                elements=elements,
                sz=float(sz[index]) if index < len(sz) and sz[index] is not None else None,
                seff=(
                    float(seff[index])
                    if index < len(seff) and seff[index] is not None
                    else None
                ),
                f_block=any(_classify(symbol) for symbol in elements),
            )
        )
    return tuple(rows)


def total_s2(sab: Sequence[Sequence[float]]) -> float:
    """``<S^2>`` of the state as the local spin analysis partitions it.

    It is the sum over **all** fragment pairs ``(A, B)`` of the ``<SA*SB>`` matrix, i.e.
    the total of the full (symmetric) matrix -- not of its lower triangle.
    """
    if not sab:
        raise LocalSpinError(
            "no <SA*SB> matrix given. Next step: pass the parsed local spin block (the "
            "matrix is read from the 'sab' entry of the block)."
        )
    return math.fsum(math.fsum(float(value) for value in row) for row in sab)


def spin_share(sab: Sequence[Sequence[float]], fragment: int) -> float | None:
    """This module's spin-share of one fragment, or ``None`` for a total singlet.

    ``fragment`` is the 0-based row index of ``sab`` (the printed table numbers fragments
    from 1).  The share is ``(sum over B of <S_A*S_B>) / <S^2>`` with ``<S^2>`` the sum over
    all fragment pairs; for a total singlet ``<S^2> = 0`` the quotient is 0/0 and ``None``
    is returned -- the caller must then report the local-spin magnitudes instead of a
    share.

    This construction is provisional: it is not taken from the cited sources.
    """
    matrix = tuple(tuple(float(value) for value in row) for row in sab)
    if not 0 <= fragment < len(matrix):
        raise LocalSpinError(
            f"fragment index {fragment} is outside the <SA*SB> matrix "
            f"({len(matrix)} fragments). Next step: pass a 0-based fragment index between "
            f"0 and {len(matrix) - 1}."
        )
    total = total_s2(matrix)
    if abs(total) <= SINGLET_TOLERANCE:
        return None
    return math.fsum(matrix[fragment]) / total


def f_block_fragments(fragment_elements: Sequence[Sequence[str]]) -> tuple[int, ...]:
    """1-based numbers of the fragments that declare an f-block element (lanthanide/actinide)."""
    return tuple(
        index
        for index, elements in enumerate(fragment_elements, start=1)
        if any(_classify(symbol) for symbol in elements)
    )


def relevant_blocks(result: ParseResult) -> tuple[tuple[str, Mapping[str, Any]], ...]:
    """The blocks to report as ``(label, block)`` pairs, in report order.

    A state-resolved output (one state-average block and/or one block per root) is reported
    in full, state average first.  An output without any state-resolved block is read as an
    SCF job: only its last block is returned, because the earlier ones belong to the
    initial guess and to pre-convergence densities.
    """
    blocks = tuple(result.sections.get("local_spin", {}).get("blocks", ()))
    if not blocks:
        return ()
    roots = [block for block in blocks if block.get("state") is not None]
    averages = [
        block
        for block in blocks
        if block.get("state") is None
        and str(block.get("state_label", "")).strip().lower().startswith("state average")
    ]
    if not roots and not averages:
        return (
            (f"SCF final wavefunction (block {len(blocks)} of {len(blocks)})", blocks[-1]),
        )
    selected: list[tuple[str, Mapping[str, Any]]] = []
    for block in averages:
        label = str(block.get("state_label", "")).strip() or "State average"
        selected.append(
            (f"state average ({label}, multiplicity {_multiplicity(block)})", block)
        )
    for block in sorted(roots, key=lambda item: int(item["state"])):
        state = int(block["state"])
        block_label = str(block.get("block_label", "")).strip()
        selected.append(
            (
                f"root {state} (State to be analyzed = {state}, "
                f"block {block_label or '?'}, multiplicity {_multiplicity(block)})",
                block,
            )
        )
    return tuple(selected)


def accepts(result: ParseResult) -> bool:
    """Whether this analyser applies (a local spin analysis block with at least one block)."""
    section = result.sections.get("local_spin", {})
    return bool(section.get("present")) and bool(section.get("blocks"))


def _multiplicity(block: Mapping[str, Any]) -> str:
    value = block.get("multiplicity")
    return "?" if value is None else str(value)


def _seff_text(block: Mapping[str, Any]) -> str:
    seff = tuple(block.get("seff", ()))
    if not seff:
        return "not printed"
    return ", ".join("?" if value is None else f"{float(value):.4f}" for value in seff)


def _spin_text(spin: float) -> str:
    """A spin quantum number as ``0`` / ``1/2`` / ``1`` (not ``0.5``)."""
    if abs(spin - round(spin)) < 1e-9:
        return str(int(round(spin)))
    if abs(2.0 * spin - round(2.0 * spin)) < 1e-9:
        return f"{int(round(2.0 * spin))}/2"
    return f"{spin:g}"


def _is_spin_eigenvalue(total: float) -> float | None:
    """The ``S`` whose ``S(S+1)`` matches ``total`` (within the purity tolerance), else None."""
    for spin in _SPIN_SCAN:
        if abs(total - spin * (spin + 1.0)) <= SPIN_PURITY_TOLERANCE:
            return spin
    return None


def _spin_value_series(limit: int = 5) -> str:
    """The first ``S(S+1)`` values a spin eigenstate can take (for the eigenvalue test)."""
    values = ", ".join(f"{spin * (spin + 1.0):g}" for spin in _SPIN_SCAN[:limit])
    return values + ", ..."


def _spin_purity_lines(
    block: Mapping[str, Any],
    total: float,
    m_s: float | None,
    *,
    broken_symmetry: bool,
) -> list[str]:
    """The ``<S^2>``-versus-``S(S+1)`` line: spin purity, or spin contamination.

    The total ``<S^2>`` is the sum over all fragment pairs of the ``<SA*SB>`` matrix.  For a
    spin eigenstate that sum equals ``S(S+1)`` of the state, and ``S`` follows from the
    block's ``multiplicity`` (``S = (multiplicity - 1)/2``).  A deviation is reported as
    spin contamination -- a property of the method (a spin-unrestricted or broken-symmetry
    single determinant is not spin-pure), not an error in the file.
    """
    multiplicity = block.get("multiplicity")
    if multiplicity is not None:
        spin = (int(multiplicity) - 1) / 2.0
        expected = spin * (spin + 1.0)
        deviation = total - expected
        head = (
            f"- Spin purity: total <S^2> = {total:.4f} against S(S+1) = {expected:.4f} for "
            f"the declared multiplicity {int(multiplicity)} (S = {_spin_text(spin)})"
        )
        if abs(deviation) <= SPIN_PURITY_TOLERANCE:
            return [
                head + ": the two agree to the printed precision, so this block is spin-pure "
                "(a spin eigenstate of its declared multiplicity)."
            ]
        lines = [
            head + f": deviation {deviation:+.4f}, i.e. spin contamination -- the "
            "wavefunction of this block is not a spin eigenstate of the declared "
            "multiplicity. This is a property of the method, not an error in the file: a "
            "spin-unrestricted or broken-symmetry single determinant mixes multiplicities."
        ]
    else:
        head = (
            f"- Spin purity: total <S^2> = {total:.4f}; this block declares no multiplicity "
            f"(an SCF local spin block prints none), so no S(S+1) is declared to compare "
            f"against."
        )
        eigen_spin = _is_spin_eigenvalue(total)
        m_s_text = "" if m_s is None else f"; the printed <SzA> values give M_S = {m_s:+.4f}"
        if eigen_spin is not None:
            return [
                head + f" The value matches S(S+1) of S = {_spin_text(eigen_spin)}{m_s_text}, "
                "so the sum is at least consistent with a spin eigenstate of that spin."
            ]
        lines = [
            head + f" {total:.4f} is not S(S+1) of any spin eigenstate (the values run "
            f"{_spin_value_series()}){m_s_text}, so this solution is spin-contaminated."
        ]
        if m_s is not None:
            spin_min = abs(m_s)
            reference = (
                "a singlet, S = 0" if spin_min == 0.0 else f"S = {_spin_text(spin_min)}"
            )
            lines.append(
                f"  Against the smallest spin the printed M_S = {m_s:+.4f} admits -- "
                f"{reference}, S(S+1) = {spin_min * (spin_min + 1.0):.4f} -- the deviation is "
                f"{total - spin_min * (spin_min + 1.0):+.4f}. A spin-unrestricted single "
                f"determinant is not in general a spin eigenstate, so a deviation here is a "
                f"property of the method, not a defect of the file."
            )
    if broken_symmetry:
        lines.append(
            "  This output carries a broken-symmetry warning, and a broken-symmetry solution "
            "is deliberately not spin-pure: the contamination is the intended non-eigenstate "
            "behaviour of that method, not a defect to repair."
        )
    return lines


def _value_text(value: float | None, *, unit: str = "") -> str:
    return "?" if value is None else f"{value:.4f}{unit}"


def _fragment_table_lines(rows: Sequence[FragmentRow]) -> list[str]:
    # The "(singlet)" tag belongs to the block, not to one row: ORCA prints "n.a." only
    # when every <SzA> is zero by definition (a singlet state).
    singlet = all(row.sz is None for row in rows)
    lines = [
        "Per-fragment local spin (fragment numbers as printed by ORCA):",
        f"  {'fragment':>8}  {'elements':<16}  {'<SzA>':>15}  {'Seff(A)':>8}",
    ]
    for row in rows:
        if row.sz is None:
            sz_text = "n.a. (singlet)" if singlet else "n.a."
        else:
            sz_text = f"{row.sz:+.4f}"
        lines.append(
            f"  {row.number:>8}  {row.element_text:<16}  {sz_text:>15}  "
            f"{_value_text(row.seff):>8}"
        )
    return lines


def _matrix_lines(sab: Sequence[Sequence[float]]) -> list[str]:
    width = 10
    # NB: the label is hoisted out of the f-string expression on purpose -- a backslash
    # inside one needs Python 3.12 (PEP 701) and this package supports 3.11.
    corner = "A\\B"
    lines = [
        "<SA*SB> matrix (row = fragment A, column = fragment B; the ORCA table is printed "
        "lower-triangular, this is the full symmetric matrix):",
        f"  {corner:>6}" + "".join(f"{f'B={index}':>{width}}" for index in range(1, len(sab) + 1)),
    ]
    for index, row in enumerate(sab, start=1):
        lines.append(
            f"  {f'A={index}':>6}"
            + "".join(f"{float(value):>{width}.4f}" for value in row)
        )
    return lines


def _interpretation_lines(
    rows: Sequence[FragmentRow],
    sab: Sequence[Sequence[float]],
    *,
    elements_known: bool,
    unknown: Sequence[str],
    block: Mapping[str, Any] | None = None,
    broken_symmetry: bool = False,
) -> list[str]:
    """The metal/ligand reading of one block: magnitudes, spin-share, <SzA> sum, checks."""
    metal = [row for row in rows if row.f_block]
    others = [row for row in rows if not row.f_block]
    lines: list[str] = []

    if metal:
        names = ", ".join(f"{row.number} ({row.element_text})" for row in metal)
        lines.append(
            f"- f-block fragment(s): {names} -- 'metal' below means any fragment whose "
            f"declared element list contains a lanthanide or actinide."
        )
    elif elements_known:
        lines.append(
            "- f-block fragment(s): none -- no declared fragment contains a lanthanide or "
            "actinide element, so this is not a metal/ligand decomposition. The tables "
            "above are printed as they stand and no metal spin-share is reported."
        )
    else:
        lines.append(
            "- f-block fragment(s): not determined -- the echoed input of this output "
            "carries no atom line with a fragment index, so the fragments cannot be named "
            "and the metal/ligand reading is skipped (the tables above are printed as they "
            "stand). Next step: check that the input declares fragments by tagging the "
            "atoms (e.g. N(1), N(2))."
        )

    if metal:
        metal_text = ", ".join(
            f"{_value_text(row.seff)} (fragment {row.number}, {row.element_text})"
            for row in metal
        )
        other_text = ", ".join(
            f"{_value_text(row.seff)} (fragment {row.number}, {row.element_text})"
            for row in others
        )
        if others:
            line = (
                f"- Local-spin magnitude: metal fragment(s) Seff = {metal_text}; the "
                f"remaining fragments' Seff = {other_text}."
            )
            metal_values = [row.seff for row in metal]
            other_values = [row.seff for row in others]
            if all(value is not None for value in metal_values + other_values):
                metal_mean = math.fsum(metal_values) / len(metal_values)
                other_mean = math.fsum(other_values) / len(other_values)
                if metal_mean > other_mean:
                    line += (
                        f" The metal carries the larger local spin (mean {metal_mean:.4f} "
                        f"vs {other_mean:.4f} over the remaining fragments)."
                    )
                else:
                    line += (
                        f" The metal does not carry the larger local spin (mean "
                        f"{metal_mean:.4f} vs {other_mean:.4f} over the remaining "
                        f"fragments); check whether the fragments or the solution are the "
                        f"intended ones."
                    )
            lines.append(line)
        else:
            lines.append(
                f"- Local-spin magnitude: all {len(rows)} declared fragments are f-block "
                f"(Seff = {metal_text}); there is no ligand fragment to compare against."
            )
    else:
        magnitudes = ", ".join(
            f"{_value_text(row.seff)} (fragment {row.number}, {row.element_text})"
            for row in rows
        )
        lines.append(f"- Local-spin magnitude: Seff = {magnitudes}.")

    available_sz = [row.sz for row in rows if row.sz is not None]
    m_s = math.fsum(available_sz) if available_sz else None

    total = total_s2(sab)
    lines.append(
        f"- Total <S^2> = sum over all fragment pairs (A, B) of <SA*SB> = {total:.4f}."
    )
    lines += _spin_purity_lines(
        block or {}, total, m_s, broken_symmetry=broken_symmetry
    )
    if abs(total) <= SINGLET_TOLERANCE:
        lines.append(
            "- Spin-share: undefined for a total singlet -- <S^2> = 0, so "
            f"{SPIN_SHARE_FORMULA} is 0/0 in every fragment. The local-spin magnitudes "
            "(Seff above) are reported instead: they stay finite and are the only "
            "fragment-spin information of this block."
        )
    else:
        lines.append(
            "- Spin-share (this module's own construction, not from the cited sources -- "
            "provisional):"
        )
        lines.append(f"    {SPIN_SHARE_FORMULA}")
        shares = []
        for row in rows:
            value = spin_share(sab, row.number - 1)
            assert value is not None  # total is non-zero here
            shares.append(value)
            tag = " (f-block)" if row.f_block else ""
            lines.append(
                f"    share(fragment {row.number}, {row.element_text}) = "
                f"{math.fsum(sab[row.number - 1]):.4f} / {total:.4f} = {value:.4f}{tag}"
            )
        lines.append(
            f"    The shares sum to {math.fsum(shares):.4f} by construction (the rows of "
            f"<SA*SB> partition the total <S^2>)."
        )

    if not available_sz:
        lines.append(
            "- Sum of <SzA>: not available -- every fragment prints \"n.a.\" here (the "
            "manual: \"for a singlet state all <SzA> values are zero by definition\"), so "
            "the sum check is skipped; the block's spin information is the Seff column and "
            "the <SA*SB> matrix."
        )
    else:
        sz_sum = m_s if m_s is not None else math.fsum(available_sz)
        complete = len(available_sz) == len(rows)
        coverage = (
            ""
            if complete
            else f" (over {len(available_sz)} of {len(rows)} fragments; the rest print n.a.)"
        )
        values = " / ".join(f"{value:+.4f}" for value in available_sz)
        if not complete:
            lines.append(
                f"- Sum of <SzA> = {sz_sum:+.4f}{coverage}, the values being {values}: a "
                f"partial sum -- the local projections add up to the total M_S of the state "
                f"only when every fragment contributes, so this number is not M_S."
            )
        elif abs(sz_sum) <= SINGLET_TOLERANCE:
            lines.append(
                f"- Sum of <SzA> = {sz_sum:+.4f}, the values being {values}: the local "
                f"projections add up to the total M_S of the state, so a zero sum is exact "
                f"whenever the state has as many alpha as beta electrons -- the fragments "
                f"carry equal and opposite local spin even though M_S = 0 (a spin-polarised "
                f"density in a state with no net spin)."
            )
        else:
            lines.append(
                f"- Sum of <SzA> = {sz_sum:+.4f}, the values being {values}: the local "
                f"projections add up to the total M_S of the state (the fragment projectors "
                f"sum to the identity), so this sum is the state's net spin projection."
            )

    if unknown:
        lines.append(
            f"- Not classified: the fragment declaration carries {', '.join(unknown)}, "
            f"which this toolkit's element table does not know; those fragments count as "
            f"non-f-block here. Next step: check the element symbols of the input."
        )

    checked = [
        (abs(sab[row.number - 1][row.number - 1] - row.seff * (row.seff + 1.0)), row)
        for row in rows
        if row.seff is not None
    ]
    if not checked:
        lines.append(
            "- Printed-number check: skipped -- this block prints no Seff value to check "
            "the <SA*SA> diagonal against."
        )
    else:
        worst, worst_row = max(checked, key=lambda item: item[0])
        if worst <= SEFF_DIAGONAL_TOLERANCE:
            lines.append(
                f"- Printed-number check (this module's own, provisional): the diagonal "
                f"<SA*SA> equals Seff(A)(Seff(A)+1) in every fragment, the largest deviation "
                f"being {worst:.2e} at fragment {worst_row.number} (limit "
                f"{SEFF_DIAGONAL_TOLERANCE:.0e}), so the printed columns are mutually "
                f"consistent."
            )
        else:
            lines.append(
                f"- Printed-number check (this module's own, provisional): the diagonal "
                f"<SA*SA> deviates from Seff(A)(Seff(A)+1) by {worst:.2e} at fragment "
                f"{worst_row.number}, above the limit {SEFF_DIAGONAL_TOLERANCE:.0e}. Next "
                f"step: check that this block was read correctly (fragment rows in order, "
                f"the analysis block of the intended job)."
            )
    return lines


def _boundary_lines() -> list[str]:
    """The Delta S_E boundary every report body carries (see the module docstring)."""
    return [
        "",
        "Boundary of this analysis (travels with the numbers above):",
        "- Local spin analysis is NOT the environment spin polarisation entropy Delta S_E "
        "of Eq. (9). It partitions the spin over fragments; Delta S_E measures how far the "
        "environment density is from the spin-unpolarised point D_E^alpha = D_E^beta = "
        "D_E/2, and it needs the two-step off-line route orca_2json (RDM2_aa / RDM2_ab / "
        "RDM2_bb; the 2-RDMs are not stored inside the .densities file) plus orca_loc "
        "(IAO-IBO localised orbitals, .loc.gbw). No printed ORCA table carries the "
        "spin-resolved environment density.",
        "    Eq. (9): Delta S_E = -2 Tr[(D_E/2) ln(D_E/2)] + Tr[D_E^alpha ln D_E^alpha] "
        "+ Tr[D_E^beta ln D_E^beta], with D_E = D_E^alpha + D_E^beta.",
        "- Because this block is a surrogate, its discriminating power is limited: an "
        "unpolarised environment and a strongly polarised one both give finite fragment "
        "Seff values here, and nothing in this report separates them on the scale below.",
        "- The scale this surrogate stands in for (measured in the reference paper on the "
        "Dy single-ion magnet 1Dy, 6H, 11 states): the wrong SCF solution gave Delta S_E = "
        "2.766 against 0.007 for the correct one, and the wrong solution degraded the "
        "downstream pre-SOC magnetic MAE from 8.6 to 357.7 cm^-1. Treat this report as a "
        "first look at the spin partition, not as that criterion.",
        "- The numbers depend on two user choices: the fragment definition (which atoms "
        "were assigned to which fragment) and the projector named in the block header "
        "(this project's fixtures use 'Loewdin*'). Values from different choices are not "
        "comparable.",
    ]


def _block_lines(
    label: str,
    block: Mapping[str, Any],
    fragment_elements: Sequence[Sequence[str]],
    *,
    broken_symmetry: bool = False,
) -> list[str]:
    sab = _matrix(block)
    if fragment_elements and len(fragment_elements) != len(sab):
        raise LocalSpinError(
            f"the fragment element list recovered from the echoed input has "
            f"{len(fragment_elements)} entries but this block's <SA*SB> matrix has "
            f"{len(sab)} fragments. Next step: check the fragment declarations of the "
            f"input -- an output holding several jobs with different fragmentations is "
            f"ambiguous, so run one job per output file."
        )
    declared = block.get("n_fragments")
    if declared is not None and int(declared) != len(sab):
        raise LocalSpinError(
            f"this block declares {int(declared)} fragments but its <SA*SB> matrix has "
            f"{len(sab)}. Next step: check the printed block; a truncated table or a mixed "
            f"output file would explain it."
        )
    rows = fragment_rows(block, fragment_elements)
    lines = ["", f"--- {label} ---"]
    extras = []
    if block.get("n_atoms") is not None:
        extras.append(f"{int(block['n_atoms'])} atoms")
    if block.get("n_basis_functions") is not None:
        extras.append(f"{int(block['n_basis_functions'])} basis functions")
    if extras:
        lines.append(f"Block scope: {', '.join(extras)}.")
    lines += _fragment_table_lines(rows)
    lines += _matrix_lines(sab)
    lines += ["Interpretation:"] + _interpretation_lines(
        rows,
        sab,
        elements_known=bool(fragment_elements),
        unknown=unknown_symbols(fragment_elements),
        block=block,
        broken_symmetry=broken_symmetry,
    )
    return lines


def _broken_symmetry_declared(result: ParseResult) -> bool:
    """Whether the output carries a broken-symmetry warning (then contamination is intended)."""
    for warning in result.sections.get("warnings", ()):
        if "broken symmetry" in str(warning).lower():
            return True
    return False


def run(result: ParseResult) -> ReportSection:
    """Build the A6 report section (per-block tables plus the Delta S_E boundary).

    Which blocks are reported follows :func:`relevant_blocks`: for a state-resolved output
    the state-average block first and then every root, for an SCF output only the last
    block (the one belonging to the converged wavefunction).
    """
    section = result.sections.get("local_spin", {})
    blocks = tuple(section.get("blocks", ()))
    if not blocks:
        raise LocalSpinError(
            "this output carries no local spin analysis block. Next step: divide the "
            "molecule into fragments in the ORCA input by tagging each atom's element "
            "symbol with a parenthesised fragment index (e.g. N(1), N(2)) and rerun; the "
            "analysis is then printed automatically (ORCA 6.1 manual section 5.1.10)."
        )
    fragment_elements = tuple(
        tuple(str(symbol) for symbol in fragment)
        for fragment in section.get("fragment_elements", ())
    )
    selected = relevant_blocks(result)
    state_resolved = any(block.get("state") is not None for block in blocks) or any(
        str(block.get("state_label", "")).strip().lower().startswith("state average")
        for block in blocks
    )

    lines = [
        "Local spin source: the \"LOCAL SPIN ANALYSIS (Loewdin* projector)\" block printed "
        "by ORCA's SCF / CASSCF module (the projector is named in the block header, and the "
        "values depend on it).",
        f"Blocks found in this output: {len(blocks)}.",
        f"Fragments declared in the input: {len(fragment_elements)}",
    ]
    for number, elements in enumerate(fragment_elements, start=1):
        lines.append(f"  {number}: {' '.join(elements) if elements else '(no atoms recovered)'}")
    if not fragment_elements:
        lines.append("  (the echoed input carries no atom line with a fragment index)")

    if state_resolved:
        lines.append(
            "Block selection: this output carries state-resolved blocks, so all of them are "
            "reported -- the state-average block first, then each root in ascending state "
            "order (the manual: \"if nroots > 1, the printing will contain the "
            "state-specific analysis of all roots\")."
        )
    else:
        lines.append(
            "Block selection: this output carries no state-resolved block, so it is read as "
            "an SCF job. ORCA prints the analysis for the initial guess, after SCF and once "
            "more for the final density; only the LAST block describes the converged "
            f"wavefunction, so block {len(blocks)} of {len(blocks)} is used and the earlier "
            "ones are left out."
        )
        unused = "; ".join(
            f"block {index} (Seff {_seff_text(block)})"
            for index, block in enumerate(blocks[:-1], start=1)
        )
        if unused:
            lines.append(
                f"Unused earlier block(s), for comparison only: {unused}. Their numbers are "
                f"not the converged wavefunction's -- reading them instead of the last "
                f"block changes the values."
            )

    broken_symmetry = _broken_symmetry_declared(result)
    for label, block in selected:
        lines += _block_lines(
            label, block, fragment_elements, broken_symmetry=broken_symmetry
        )

    lines += _boundary_lines()
    return ReportSection(title="A6 local spin analysis", body="\n".join(lines))


def evidence() -> tuple[Evidence, ...]:
    """Provenance of this module's statements (manual, literature, measured here)."""
    return (_EVIDENCE_MANUAL, _EVIDENCE_DSE, _EVIDENCE_LOCAL_SPIN, _EVIDENCE_MEASURED)
