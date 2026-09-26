"""A5: projection-basis declaration check for crystal-field parameters ``B_k^q``.

A ``B_k^q`` value is a coefficient *in a declared scheme*, not a property of the
molecule: one wavefunction projected two ways gives two different parameter sets
and both may be right.  This module checks that the scheme is declared and
decides, per ``(k, q)``, whether two parameter sets may be compared at all.

The measurement this module encodes (Chilton SI, Table S1)
----------------------------------------------------------
Read 2026-09-25 from this project's conversion of the Supporting Information of
Chilton, *Chem. Soc. Rev.* 2025, **54**(24), 11468-11487 (close reading:
``文献细读/细读_Chilton指南_SI.md`` section 4.4, which gives the table and its line
number).  One SA-CASSCF-SO wavefunction of one Dy(III) compound, projected five
ways, gives different ``B_k^q`` (cm^-1); the three columns used here are::

    (k, q)   SINGLE_ANISO   SINGLE_ANISO   angmom_suite
             J = 15/2       L = 5          J = 15/2
    (2,  0)   1938           1879           1938
    (2,  2)    198            -43            201
    (2, -2)    124            239            125
    (4,  0)     66            125             66
    (6,  2)     18             -4             18

Read off that table:

* same manifold, different program -> comparable: 1938/1938, 66/66, 198 vs 201;
* axial ``(2, 0)`` across manifolds -> 1938 vs 1879, a 3% shift, so the leading
  axial parameter may be compared across schemes *roughly*;
* transverse and high-order parameters may not: ``(2, 2)`` +198 vs -43 (sign
  flip), ``(2, -2)`` 124 vs 239, ``(4, 0)`` 66 vs 125 (factor ~1.9), ``(6, 2)``
  18 vs -4 (sign flip).

The engineering rule implemented (the single licence this measurement supports)::

    both sides state all five items and they agree          -> comparable
    (k, q) = (2, 0), no single scheme established           -> comparable-roughly
    everything else                                          -> refused

The licence is ``(2, 0)`` and **not** "every ``q = 0``": Table S1 moved the
high-order axial term ``(4, 0)`` by a factor ~1.9 (66 -> 125), so a plain
"q == 0" test would wave that comparison through.  The 3% is one compound and one
wavefunction, so it is the *scale* of the cross-scheme effect, not a bound
(provisional inference of this module).

The five items every ``B_k^q`` output must declare
--------------------------------------------------
Distilled in this project's close reading of the CF-Hamiltonian paper (``文献细读/
细读_CF哈密顿量_Chan2025.md`` section 5.3, "projection-basis discipline", citing
Peng et al., *J. Phys. Chem. Lett.* 2025, **16**(47), 12312-12320): (i) operator
convention (Stevens vs Wybourne, normalisation), (ii) the ``|J M>`` projection
manifold (a ``J`` multiplet vs an ``L``-``S`` term vs a pseudospin doublet),
(iii) units, (iv) the ``z``-axis definition (a geometric axis or a magnetic
axis), (v) the frame origin and chirality.  Missing any one of them leaves the
numbers non-comparable with another program's output.  The item list is stated
here with each item's consequence rather than assumed; the paper's own wording
was not re-read for this module (the close reading is the source of the list).

Interface and conventions of this module
----------------------------------------
* No file I/O: the caller reads the JSON and passes a mapping.  Payload::

      {"parameters": {"2,0": 1938.0, "2,2": 198.0, ...},   # or "[2, 2]" / [2, 2]
       "declaration": {"convention": ..., "projection": ..., "units": ...,
                       "z_axis": ..., "origin": ...},        # any subset
       "compare": {"label": ..., "parameters": {...}, "declaration": {...}},
       "label": ...}                                        # optional, names this set

* Parameter keys are accepted exactly as the ``"k, q"`` spelling (``"2,0"``) or
  as a two-element pair (``[2, 2]``, ``(2, 2)``); anything else, a duplicate
  ``(k, q)``, a rank outside ``{2, 4, 6}`` or a non-finite value is refused with
  a ``DeclarationError`` whose message ends in a "Next step:" sentence.
* Verdicts: ``comparable`` / ``comparable-roughly`` / ``refused`` for ``(k, q)``
  present on both sides, ``missing-on-one-side`` when only one side lists it.
* Deterministic: same payload -> same report text; no clocks, no randomness, no
  hidden state.  A missing declaration item is reported as missing, never filled
  in from a default.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass

from ..knowledge.models import (
    EVIDENCE_LITERATURE,
    EVIDENCE_MEASURED,
    Evidence,
    ReportSection,
)

#: Report title (the menu path is wired by the caller; this string is the pinned
#: section title of the analysis layer).
TITLE = "A5 projection-basis declaration check"

#: The five items any ``B_k^q`` output must declare before its numbers can be
#: compared with another program's.  Order is fixed and is part of the output.
DECLARATION_FIELDS: tuple[str, ...] = (
    "convention",
    "projection",
    "units",
    "z_axis",
    "origin",
)

VERDICT_COMPARABLE = "comparable"
VERDICT_COMPARABLE_ROUGHLY = "comparable-roughly"
VERDICT_REFUSED = "refused"
VERDICT_MISSING_ON_ONE_SIDE = "missing-on-one-side"
VERDICTS: tuple[str, ...] = (
    VERDICT_COMPARABLE,
    VERDICT_COMPARABLE_ROUGHLY,
    VERDICT_REFUSED,
    VERDICT_MISSING_ON_ONE_SIDE,
)

#: The one parameter class the Table S1 measurement licenses for a rough
#: comparison across projection schemes: the leading axial parameter.  Not every
#: ``q = 0`` -- ``(4, 0)`` moved by a factor ~1.9 (66 -> 125) between the same two
#: manifolds, so high-order axial terms are excluded too.
ROUGH_COMPARABLE_PARAMETERS: frozenset[tuple[int, int]] = frozenset({(2, 0)})

#: Ranks of the extended Stevens expansion for an f shell (``k <= 2l = 6``, even
#: ``k`` only by time-reversal symmetry).
K_RANKS: tuple[int, ...] = (2, 4, 6)

#: Accepted keys of one parameter set.
_SET_KEYS = ("parameters", "declaration", "label")

#: Accepted keys of the top-level payload.
_PAYLOAD_KEYS = ("parameters", "declaration", "compare", "label")

#: Why a missing declaration item makes the numbers non-comparable.  Carried by
#: every :class:`DeclarationItem` and printed for the missing ones.
_FIELD_CONSEQUENCE: Mapping[str, str] = {
    "convention": (
        "the operator convention (Stevens vs Wybourne, normalisation, cosine/sine phase of "
        "the q <-> -q pair) is not recorded, so the same (k, q) cannot be tied to one "
        "operator: another program's value may differ by a normalisation or phase factor "
        "while looking perfectly well formed"
    ),
    "projection": (
        "the projection manifold (J multiplet vs L-S term vs pseudospin doublet) is not "
        "recorded, so there is no way to know which angular-momentum basis the numbers "
        "belong to; Table S1 measured a sign flip (B_2^2 = +198 vs -43) and a factor ~1.9 "
        "change (B_4^0 = 66 vs 125) between two manifolds of the same wavefunction"
    ),
    "units": (
        "the unit (cm^-1 vs K vs meV) is not recorded, so numerical agreement or "
        "disagreement cannot be interpreted; the cm^-1 <-> K factor 1.4388 would read as a "
        "physics difference"
    ),
    "z_axis": (
        "the z-axis definition (a geometric axis vs the ground-Kramers-doublet magnetic "
        "axis, and which symmetry axis) is not recorded, so a rotation of the frame cannot "
        "be separated from a change of the parameters; transverse (q != 0) parameters are "
        "frame dependent by construction"
    ),
    "origin": (
        "the frame origin and chirality are not recorded, so a translation of the frame or "
        "an improper axis choice cannot be excluded: two sets can differ without any "
        "physical difference"
    ),
}

#: Table S1 values quoted by the report and by the provenance record: three of the
#: five columns, ``(k, q) -> (SINGLE_ANISO J = 15/2, SINGLE_ANISO L = 5,
#: angmom_suite J = 15/2)`` in cm^-1.  Source: Chilton SI Table S1 (see module
#: docstring); the tests transcribe the whole table independently and check these
#: entries against it.
_TABLE_S1_EXAMPLE: Mapping[tuple[int, int], tuple[float, float, float]] = {
    (2, 0): (1938.0, 1879.0, 1938.0),
    (2, 2): (198.0, -43.0, 201.0),
    (2, -2): (124.0, 239.0, 125.0),
    (4, 0): (66.0, 125.0, 66.0),
    (6, 2): (18.0, -4.0, 18.0),
}

_COLUMN_J15 = "SINGLE_ANISO, J = 15/2"
_COLUMN_L5 = "SINGLE_ANISO, L = 5"
_COLUMN_ANGMOM = "angmom_suite, J = 15/2"


class DeclarationError(ValueError):
    """Malformed declaration payload or malformed parameter key.

    Every message ends in a "Next step:" sentence telling the caller what to fix.
    """


@dataclass(frozen=True)
class DeclarationItem:
    """One of the five required declaration items.

    ``field`` is the item name (one of :data:`DECLARATION_FIELDS`), ``present``
    says whether the caller stated it, ``value`` is the stated text (empty when
    absent) and ``consequence`` is why the numbers stop being comparable with
    another program's output when the item is missing.
    """

    field: str
    present: bool
    value: str
    consequence: str


@dataclass(frozen=True)
class ComparisonItem:
    """Verdict for one ``(k, q)`` of two parameter sets.

    ``verdict`` is one of :data:`VERDICTS`; ``reason`` states the rule that
    produced it (and, for a refusal, the measured example behind the rule).
    ``value_a`` / ``value_b`` are the two numbers, ``None`` on the side that does
    not list this ``(k, q)``.
    """

    k: int
    q: int
    verdict: str
    reason: str
    value_a: float | None = None
    value_b: float | None = None


# --- provenance --------------------------------------------------------------

_EVIDENCE_TABLE_S1 = Evidence(
    kind=EVIDENCE_LITERATURE,
    text=(
        "One SA-CASSCF-SO wavefunction of one Dy(III) compound, projected five ways, gives "
        "different B_k^q (cm^-1): axial B_2^0 = 1938 (SINGLE_ANISO, J = 15/2) vs 1879 "
        "(SINGLE_ANISO, L = 5) -- 3%; transverse B_2^2 = +198 vs -43 (sign flip), B_2^-2 = "
        "124 vs 239; high order B_4^0 = 66 vs 125 (factor ~1.9), B_6^2 = 18 vs -4 (sign "
        "flip). Within one manifold two programs agree to the printed precision: "
        "SINGLE_ANISO vs angmom_suite, both J = 15/2 -- B_2^0 1938/1938, B_2^2 198 vs 201, "
        "B_4^0 66/66. Hence the rule enforced here: same-manifold numbers from different "
        "programs are comparable; the leading axial (2, 0) may be compared across schemes "
        "roughly; every other (k, q) needs both sides to declare the same projection "
        "scheme. Read from the Supporting Information of this paper, Table S1 (via this "
        "project's close reading, section 4.4)."
    ),
    ref=(
        "Chilton N. F., Chem. Soc. Rev., 2025, 54(24), 11468-11487, DOI 10.1039/d5cs00493d "
        "(Supporting Information, Table S1)"
    ),
    bibkey="chilton2025abinitio",
    url="https://doi.org/10.1039/d5cs00493d",
)

_EVIDENCE_FIVE_ITEMS = Evidence(
    kind=EVIDENCE_LITERATURE,
    text=(
        "The five items any B_k^q output must declare before its numbers can be compared "
        "with another program's: (i) operator convention (Stevens/Wybourne, normalisation), "
        "(ii) the |J M> projection manifold (J multiplet vs L-S term vs pseudospin "
        "doublet), (iii) units, (iv) z-axis definition (a geometric axis or a magnetic "
        "axis), (v) frame origin and chirality; missing any one of them leaves the numbers "
        "non-comparable, stated there with the same sign-flip example (B_2^2 = +198 vs -43). "
        "Item list taken from this project's close reading of this paper, section 5.3 "
        "('projection-basis discipline'); the paper's own wording was not re-read for this "
        "module."
    ),
    ref=(
        "Peng L., Liu S., Zhang X., Chen X., Li C., Ung S. F., Cheng H.-P., Chan G. K.-L., "
        "J. Phys. Chem. Lett., 2025, 16(47), 12312-12320, DOI 10.1021/acs.jpclett.5c02971"
    ),
    bibkey="peng2025accurate",
    url="https://doi.org/10.1021/acs.jpclett.5c02971",
)

_EVIDENCE_REPRODUCTION = Evidence(
    kind=EVIDENCE_MEASURED,
    text=(
        "Reproduction of the Table S1 numbers by this module's tests: the five columns were "
        "transcribed (Table S1 has 27 (k, q) rows) and the ratios recomputed. Measured: the "
        "two same-manifold columns agree for the large parameters (B_2^0 1938/1938, B_4^0 "
        "66/66; B_2^2 198 vs 201, 1.5% apart), the cross-manifold axial shift is 3.0% (1938 "
        "-> 1879), and the cross-manifold transverse/high-order parameters either flip sign "
        "(B_2^2 +198 -> -43, B_6^2 18 -> -4) or change by a factor ~1.9 (B_4^0 66 -> 125, "
        "B_2^-2 124 -> 239). These ratios are the origin of the module's verdicts: "
        "'comparable' needs one declared scheme, 'comparable-roughly' is licensed for "
        "(2, 0) only, everything else is refused. One compound, one wavefunction -- "
        "provisional, a scale rather than a bound."
    ),
    ref="reproduced by tests/test_cf_declaration.py (Table S1 transcription, 2026-09-26)",
)


# --- declaration items -------------------------------------------------------


def check_declaration(declaration: Mapping[str, str]) -> tuple[DeclarationItem, ...]:
    """The five required declaration items, in :data:`DECLARATION_FIELDS` order.

    One item per field, stating whether the caller declared it and, for a missing
    one, why its absence makes the numbers non-comparable outside this set.  An
    item that is absent, ``None`` or blank text is reported as ``present=False``
    with an empty value -- the caller decides what an undeclared scheme means.
    An item that is present but not text, or a key outside the five fields, is an
    error: a declaration is a short statement, not a number or a structure.
    """
    if not isinstance(declaration, Mapping):
        raise DeclarationError(
            f"declaration must be a mapping of the five items, got "
            f"{type(declaration).__name__}. Next step: pass a mapping with any subset of "
            f"{DECLARATION_FIELDS}, e.g. {{'projection': 'J = 15/2', 'units': 'cm^-1'}}."
        )
    unknown = tuple(sorted(str(key) for key in declaration if key not in DECLARATION_FIELDS))
    if unknown:
        raise DeclarationError(
            f"declaration has unknown item(s) {list(unknown)}; the five required items are "
            f"{list(DECLARATION_FIELDS)}. Next step: rename the unknown item to one of those "
            f"five (fold extra facts into the nearest one, e.g. the normalisation into "
            f"'convention'; the producing program belongs in the set's label) or drop it -- "
            f"a sixth item is not checked and would only look declared."
        )
    items = []
    for field in DECLARATION_FIELDS:
        raw = declaration.get(field)
        if raw is None or (isinstance(raw, str) and not raw.strip()):
            items.append(
                DeclarationItem(
                    field=field,
                    present=False,
                    value="",
                    consequence=_FIELD_CONSEQUENCE[field],
                )
            )
            continue
        if not isinstance(raw, str):
            raise DeclarationError(
                f"declaration item {field!r} must be text, got {raw!r} "
                f"({type(raw).__name__}). Next step: state it in words, e.g. units 'cm^-1' "
                f"or z_axis 'ground Kramers doublet magnetic axis'."
            )
        items.append(
            DeclarationItem(
                field=field,
                present=True,
                value=raw.strip(),
                consequence=_FIELD_CONSEQUENCE[field],
            )
        )
    return tuple(items)


def _scheme_gap(
    items_a: tuple[DeclarationItem, ...], items_b: tuple[DeclarationItem, ...]
) -> tuple[bool, tuple[str, ...], tuple[str, ...]]:
    """``(one scheme established, items missing, items stated differently)``.

    One scheme is established only when both sides state all five items and the
    five values are equal (compared as stripped text, case-sensitive: a
    declaration that differs only in spelling is still a different statement, and
    the conservative reading is to compare less, not more).
    """
    missing = tuple(
        field
        for field, item_a, item_b in zip(DECLARATION_FIELDS, items_a, items_b)
        if not (item_a.present and item_b.present)
    )
    differing = tuple(
        field
        for field, item_a, item_b in zip(DECLARATION_FIELDS, items_a, items_b)
        if item_a.present and item_b.present and item_a.value != item_b.value
    )
    return (not missing and not differing), missing, differing


def _gap_text(missing: tuple[str, ...], differing: tuple[str, ...]) -> str:
    parts = []
    if missing:
        parts.append("not stated on both sides: " + ", ".join(missing))
    if differing:
        parts.append("stated differently: " + ", ".join(differing))
    return "; ".join(parts) if parts else "no gap"


def _scheme_summary(items: tuple[DeclarationItem, ...]) -> str:
    values = {item.field: item.value for item in items}
    return (
        f"manifold '{values['projection']}', operator convention "
        f"'{values['convention']}', units '{values['units']}'"
    )


# --- parameter keys and values ----------------------------------------------


def _as_int(value: object, key: object) -> int:
    """A ``k`` or ``q`` component of a parameter key, as an exact integer."""
    if isinstance(value, bool):
        raise DeclarationError(
            f"parameter key {key!r} has a boolean component {value!r}, which is not a "
            f"(k, q) index. Next step: write the key as 'k,q' (e.g. '2,0')."
        )
    if isinstance(value, int):
        return value
    if isinstance(value, float) and value.is_integer():
        return int(value)
    if isinstance(value, str):
        try:
            return int(value.strip())
        except ValueError:
            raise DeclarationError(
                f"parameter key {key!r} has non-integer component {value!r}. "
                f"Next step: write the key as 'k,q' with integer k and q, e.g. '2,0' or "
                f"'[2, 2]'."
            ) from None
    raise DeclarationError(
        f"parameter key {key!r} has component {value!r} ({type(value).__name__}), which is "
        f"not an integer. Next step: write the key as 'k,q' with integer k and q."
    )


def _parse_key(key: object) -> tuple[int, int]:
    """``"k,q"`` / ``"[k, q]"`` / a two-element pair -> ``(k, q)``.

    Only these spellings are accepted; everything else is refused rather than
    guessed, because a silently mis-read key would compare the wrong parameters.
    """
    if isinstance(key, str):
        text = key.strip()
        if text.startswith("[") and text.endswith("]"):
            text = text[1:-1]
        parts = text.split(",")
        if len(parts) != 2:
            raise DeclarationError(
                f"parameter key {key!r} is not a (k, q) pair. Next step: write it as 'k,q' "
                f"(e.g. '2,0', '-4,-3'), as '[k, q]', or as the pair [k, q] itself."
            )
        k = _as_int(parts[0], key)
        q = _as_int(parts[1], key)
    elif isinstance(key, (tuple, list)) and len(key) == 2:
        k = _as_int(key[0], key)
        q = _as_int(key[1], key)
    else:
        raise DeclarationError(
            f"parameter key {key!r} ({type(key).__name__}) is not a (k, q) pair. "
            f"Next step: write it as 'k,q' (e.g. '2,0'), as '[k, q]', or as the pair [k, q] "
            f"itself."
        )
    if k not in K_RANKS:
        raise DeclarationError(
            f"parameter key {key!r} gives k = {k}; the crystal-field expansion of an f "
            f"shell has even ranks k in {list(K_RANKS)} (k <= 2l = 6, and time reversal "
            f"drops the odd ranks). Next step: use a (k, q) with k in {list(K_RANKS)}."
        )
    if abs(q) > k:
        raise DeclarationError(
            f"parameter key {key!r} gives |q| = {abs(q)} > k = {k}, which is not a "
            f"parameter of this expansion. Next step: use |q| <= k, e.g. '2,0' or '4,-3'."
        )
    return k, q


def _parse_value(value: object, key: object) -> float:
    """A ``B_k^q`` value: a finite real number."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise DeclarationError(
            f"B_k^q value for key {key!r} must be a number, got {value!r} "
            f"({type(value).__name__}). Next step: pass the parameter as a number, e.g. "
            f'"2,0": 1938.0; if the value is unknown, omit the (k, q) entry instead of '
            f"passing null."
        )
    number = float(value)
    if not math.isfinite(number):
        raise DeclarationError(
            f"B_k^q value for key {key!r} is {number}, which is not a comparable number. "
            f"Next step: drop the non-converged (k, q) entry -- a parameter that is NaN or "
            f"infinite cannot be compared with anything."
        )
    return number


def _read_parameters(raw: object, what: str) -> dict[tuple[int, int], float]:
    if not isinstance(raw, Mapping):
        raise DeclarationError(
            f"'parameters' of {what} must be a mapping of (k, q) -> value, got "
            f"{type(raw).__name__}. Next step: pass e.g. {{\"2,0\": 1938.0, \"2,2\": 198.0}}."
        )
    parameters: dict[tuple[int, int], float] = {}
    for key, value in raw.items():
        pair = _parse_key(key)
        if pair in parameters:
            raise DeclarationError(
                f"the (k, q) = {pair} appears twice in {what} (two spellings of one key). "
                f"Next step: give each (k, q) exactly once -- '2,0', '[2, 0]' and the pair "
                f"(2, 0) all denote the same parameter."
            )
        parameters[pair] = _parse_value(value, key)
    return parameters


def _read_set(
    block: object, what: str
) -> tuple[dict[tuple[int, int], float], tuple[DeclarationItem, ...]]:
    """One parameter set: its ``(k, q) -> value`` mapping and its declaration items."""
    if not isinstance(block, Mapping):
        raise DeclarationError(
            f"{what} must be a mapping with 'parameters' (and optionally 'declaration'), "
            f"got {type(block).__name__}. Next step: pass "
            f"{{\"parameters\": {{...}}, \"declaration\": {{...}}}}."
        )
    unknown = tuple(sorted(str(key) for key in block if key not in _SET_KEYS))
    if unknown:
        raise DeclarationError(
            f"{what} has unknown key(s) {list(unknown)}; expected {list(_SET_KEYS)}. "
            f"Next step: rename the key to one of those (the declaration items go inside "
            f"'declaration'), or, if this is a whole payload, pass only its 'parameters' and "
            f"'declaration' here -- the 'compare' block is handled by run(), not by a set."
        )
    if "parameters" not in block:
        raise DeclarationError(
            f"{what} has no 'parameters' mapping, so there is nothing to compare. "
            f"Next step: pass the B_k^q values under 'parameters', keyed 'k,q' (e.g. "
            f"'2,0') or as the pair [k, q]."
        )
    parameters = _read_parameters(block["parameters"], what)
    declaration = block.get("declaration")
    items = check_declaration({} if declaration is None else declaration)
    return parameters, items


# --- comparison --------------------------------------------------------------


def _number(value: float) -> str:
    return f"{value:.10g}"


def _difference_note(value_a: float, value_b: float) -> str:
    """``" (3.04% apart, signs differ)"``; empty when the reference is zero."""
    if value_a == 0.0:
        return ""
    percent = 100.0 * abs(value_a - value_b) / abs(value_a)
    note = f" ({percent:.3g}% apart"
    if value_a * value_b < 0.0:
        note += ", signs differ"
    return note + ")"


def _pair_text(item: ComparisonItem) -> str:
    if item.value_a is None:
        return f"A = not given, B = {_number(item.value_b)}"
    if item.value_b is None:
        return f"A = {_number(item.value_a)}, B = not given"
    return (
        f"A = {_number(item.value_a)}, B = {_number(item.value_b)}"
        f"{_difference_note(item.value_a, item.value_b)}"
    )


def _comparable_reason(summary: str) -> str:
    return (
        f"both sides state all five declaration items and they agree ({summary}), so both "
        f"numbers are coefficients of the same operator in the same manifold and frame: "
        f"they are directly comparable whatever their numerical difference. Same-manifold "
        f"numbers from different programs agree to each program's own precision, not "
        f"exactly (Table S1: B_2^0 1938/1938, B_4^0 66/66, B_2^2 198 vs 201)."
    )


def _rough_reason(gap: str) -> str:
    axial_j15, axial_l5, _ = _TABLE_S1_EXAMPLE[(2, 0)]
    percent = 100.0 * abs(axial_j15 - axial_l5) / abs(axial_j15)
    return (
        f"(k, q) = (2, 0) is the leading axial parameter, the one class this measurement "
        f"licenses for a rough comparison across projection schemes: Table S1 moved it from "
        f"{_number(axial_j15)} ({_COLUMN_J15}) to {_number(axial_l5)} ({_COLUMN_L5}), "
        f"{percent:.3g}%, for the same wavefunction. The two sides do not establish one "
        f"scheme ({gap}). The 3% is the scale of the effect on one compound, not a bound, "
        f"and it does not carry over to any other (k, q) -- in particular not to the "
        f"high-order axial terms."
    )


def _measured_pair(k: int, q: int) -> str:
    """What Table S1 measured for this ``(k, q)``, or the set of measured moves."""
    values = _TABLE_S1_EXAMPLE.get((k, q))
    if values is None:
        return (
            " Table S1 measured, for the (k, q) it lists as non-zero, moves between two "
            "manifolds of one wavefunction: B_2^2 = +198 vs -43 (sign flip), B_2^-2 = 124 vs "
            "239, B_4^0 = 66 vs 125 (a factor ~1.9), B_6^2 = 18 vs -4 (sign flip)."
        )
    j15, l5, _ = values
    if j15 * l5 < 0.0:
        note = " (a sign flip)"
    elif j15 != 0.0 and abs(l5 / j15) >= 1.05:
        note = f" (a factor {abs(l5 / j15):.2g})"
    else:
        note = ""
    return (
        f" Table S1 measured this (k, q) directly: {_number(j15)} ({_COLUMN_J15}) vs "
        f"{_number(l5)} ({_COLUMN_L5}){note}."
    )


def _refusal_reason(k: int, q: int, gap: str) -> str:
    measured = _measured_pair(k, q)
    if q != 0:
        why = (
            f"(k, q) = ({k}, {q}) is a transverse component (q != 0), and transverse "
            f"parameters are frame and manifold dependent: two such numbers from different "
            f"projection schemes can differ in sign and magnitude with neither being "
            f"wrong.{measured}"
        )
    else:
        why = (
            f"(k, q) = ({k}, 0) is a high-order axial parameter (k = {k} >= 4), and "
            f"high-order parameters are not comparable across projection schemes either -- "
            f"only the leading axial (2, 0) is.{measured}"
        )
    return (
        f"{why} The two sides do not establish one declared scheme ({gap}). "
        f"Next step: state all five declaration items on both sides -- if they then agree "
        f"the verdict becomes 'comparable'; if they differ, compare within one scheme only, "
        f"or restrict the cross-scheme reading to (2, 0) and treat it as rough."
    )


def _missing_reason(k: int, q: int, in_first: bool, value: float) -> str:
    holder = "the first set (A)" if in_first else "the second set (B)"
    other = "the second set (B)" if in_first else "the first set (A)"
    return (
        f"(k, q) = ({k}, {q}) is given by {holder} only (value {_number(value)}); {other} "
        f"does not list it, so there is nothing to compare. Next step: supply this (k, q) on "
        f"both sides, or confirm that the side without it sets it to zero (a symmetry zero) "
        f"before reading the two sets as one expansion."
    )


def _compare(
    parameters_a: dict[tuple[int, int], float],
    items_a: tuple[DeclarationItem, ...],
    parameters_b: dict[tuple[int, int], float],
    items_b: tuple[DeclarationItem, ...],
) -> tuple[ComparisonItem, ...]:
    """Per-``(k, q)`` verdicts for two already-parsed parameter sets."""
    established, missing, differing = _scheme_gap(items_a, items_b)
    gap = _gap_text(missing, differing)
    summary = _scheme_summary(items_a) if established else ""
    results = []
    for key in sorted(set(parameters_a) | set(parameters_b)):
        k, q = key
        if key not in parameters_a or key not in parameters_b:
            in_first = key in parameters_a
            value = parameters_a[key] if in_first else parameters_b[key]
            results.append(
                ComparisonItem(
                    k=k,
                    q=q,
                    verdict=VERDICT_MISSING_ON_ONE_SIDE,
                    reason=_missing_reason(k, q, in_first, value),
                    value_a=parameters_a.get(key),
                    value_b=parameters_b.get(key),
                )
            )
            continue
        if established:
            verdict = VERDICT_COMPARABLE
            reason = _comparable_reason(summary)
        elif key in ROUGH_COMPARABLE_PARAMETERS:
            verdict = VERDICT_COMPARABLE_ROUGHLY
            reason = _rough_reason(gap)
        else:
            verdict = VERDICT_REFUSED
            reason = _refusal_reason(k, q, gap)
        results.append(
            ComparisonItem(
                k=k,
                q=q,
                verdict=verdict,
                reason=reason,
                value_a=parameters_a[key],
                value_b=parameters_b[key],
            )
        )
    return tuple(results)


def compare_parameters(set_a: Mapping, set_b: Mapping) -> tuple[ComparisonItem, ...]:
    """Verdict per ``(k, q)`` for two parameter sets, ordered by ``(k, q)``.

    Each set is the ``{"parameters": ..., "declaration": ...}`` mapping form (a
    ``"label"`` is accepted and ignored).  A ``(k, q)`` present on both sides is

    * ``comparable`` when both sides state all five declaration items and they
      agree -- the same declared scheme, whatever program produced the numbers;
    * ``comparable-roughly`` when it is the leading axial ``(2, 0)`` and no single
      scheme is established: the one cross-scheme comparison the Table S1
      measurement supports, to within its measured 3% scale;
    * ``refused`` otherwise (any transverse ``q != 0`` or high-order ``k >= 4``
      parameter), with the measured example behind the refusal;

    and ``missing-on-one-side`` when only one side lists it.
    """
    parameters_a, items_a = _read_set(set_a, "the first set (A)")
    parameters_b, items_b = _read_set(set_b, "the second set (B)")
    return _compare(parameters_a, items_a, parameters_b, items_b)


# --- report ------------------------------------------------------------------


def _read_label(raw: object, default: str, what: str) -> str:
    if raw is None:
        return default
    if not isinstance(raw, str) or not raw.strip():
        raise DeclarationError(
            f"{what} label must be non-empty text, got {raw!r}. Next step: pass a short "
            f"label such as 'SINGLE_ANISO, J = 15/2', or omit the key."
        )
    return raw.strip()


def _rule_of_thumb_lines() -> list[str]:
    (axial_j15, axial_l5, _) = _TABLE_S1_EXAMPLE[(2, 0)]
    (t_j15, t_l5, _) = _TABLE_S1_EXAMPLE[(2, 2)]
    (s_j15, s_l5, _) = _TABLE_S1_EXAMPLE[(2, -2)]
    (h4_j15, h4_l5, _) = _TABLE_S1_EXAMPLE[(4, 0)]
    (h6_j15, h6_l5, _) = _TABLE_S1_EXAMPLE[(6, 2)]
    axial_percent = 100.0 * abs(axial_j15 - axial_l5) / abs(axial_j15)
    return [
        "Rule of thumb (measured, one compound): Chilton, Chem. Soc. Rev. 2025, 54(24), "
        "11468-11487, Table S1 -- one SA-CASSCF-SO wavefunction of one Dy(III) compound "
        "projected five ways gives different B_k^q (cm^-1):",
        f"  axial (2, 0):  {_number(axial_j15)} ({_COLUMN_J15}) vs {_number(axial_l5)} "
        f"({_COLUMN_L5}), {axial_percent:.3g}% apart -> comparable across projection "
        f"schemes, roughly only;",
        f"  transverse:    (2, 2) {_number(t_j15)} vs {_number(t_l5)} (sign flip) and "
        f"(2, -2) {_number(s_j15)} vs {_number(s_l5)} -> not comparable across schemes;",
        f"  high order:    (4, 0) {_number(h4_j15)} vs {_number(h4_l5)} and (6, 2) "
        f"{_number(h6_j15)} vs {_number(h6_l5)} -> not comparable across schemes;",
        f"  same manifold, different program ({_COLUMN_J15} vs {_COLUMN_ANGMOM}): "
        f"(2, 0) 1938/1938, (2, 2) 198/201, (4, 0) 66/66 -> comparable.",
        "  Rule applied here: only (2, 0) may be compared roughly when the two sides do not "
        "establish the same declared scheme; every other (k, q) needs both sides to declare "
        "the same projection scheme. Provisional: the 3% comes from one compound and one "
        "wavefunction, and it is a scale rather than a bound.",
    ]


def run(payload: Mapping) -> ReportSection:
    """Build the A5 report section from a parsed declaration payload.

    The caller reads the JSON file; this function takes the mapping (see the
    module docstring for the shape) and does no I/O.  The body carries the
    declaration table (one line per required item plus the consequence of a
    missing one), the per-``(k, q)`` comparison verdicts when a ``"compare"``
    block is present, and the rule-of-thumb line quoting the Table S1 numbers.
    """
    if not isinstance(payload, Mapping):
        raise DeclarationError(
            f"payload must be a mapping, got {type(payload).__name__}. Next step: pass an "
            f"object with keys {list(_PAYLOAD_KEYS)}, e.g. {{\"parameters\": "
            f"{{\"2,0\": 1938.0}}, \"declaration\": {{\"projection\": \"J = 15/2\"}}}}."
        )
    unknown = tuple(sorted(str(key) for key in payload if key not in _PAYLOAD_KEYS))
    if unknown:
        raise DeclarationError(
            f"payload has unknown key(s) {list(unknown)}; expected {list(_PAYLOAD_KEYS)}. "
            f"Next step: rename the key to one of those (declaration items go inside "
            f"'declaration', a second set goes under 'compare')."
        )
    if "parameters" not in payload:
        raise DeclarationError(
            "payload has no 'parameters' mapping, so there are no B_k^q values to check a "
            "declaration for. Next step: pass the parameters of this set, e.g. "
            "{\"parameters\": {\"2,0\": 1938.0, \"2,2\": 198.0}, \"declaration\": {...}}."
        )
    parameters_a = _read_parameters(payload["parameters"], "the payload")
    declaration_a = payload.get("declaration")
    items_a = check_declaration({} if declaration_a is None else declaration_a)
    label_a = _read_label(payload.get("label"), "this set", "the set")

    comparison: tuple[ComparisonItem, ...] | None = None
    label_b = ""
    compare_block = payload.get("compare")
    if compare_block is not None:
        if not isinstance(compare_block, Mapping):
            raise DeclarationError(
                f"'compare' must be a mapping with 'parameters' (and optionally "
                f"'declaration' and 'label'), got {type(compare_block).__name__}. "
                f"Next step: pass {{\"label\": ..., \"parameters\": {{...}}, "
                f"\"declaration\": {{...}}}}."
            )
        label_b = _read_label(compare_block.get("label"), "the comparison set", "the compare")
        parameters_b, items_b = _read_set(compare_block, "the compare block")
        comparison = _compare(parameters_a, items_a, parameters_b, items_b)

    axial = sum(1 for (_, q) in parameters_a if q == 0)
    lines = [
        f"Set: {label_a}",
        f"B_k^q given: {len(parameters_a)} (k, q) entries (axial q = 0: {axial}; "
        f"transverse |q| >= 1: {len(parameters_a) - axial})",
        "",
        "Declaration (all five items are required; a missing one is why the numbers stop",
        "being comparable with another program's output):",
    ]
    for position, item in enumerate(items_a, start=1):
        if item.present:
            lines.append(f"  {position}. {item.field:<10} present   {item.value}")
        else:
            lines.append(f"  {position}. {item.field:<10} MISSING")
            lines.append(f"      consequence: {item.consequence}")
    stated = sum(1 for item in items_a if item.present)
    lines.append(f"  stated: {stated} of {len(items_a)} items")

    if comparison is None:
        lines += [
            "",
            "Comparison: none requested. Pass a 'compare' block (a second set's parameters "
            "and declaration) to check one set against another.",
        ]
    else:
        lines += ["", f"Comparison: {label_a} (A) vs {label_b} (B)"]
        for item in comparison:
            lines.append(f"  ({item.k},{item.q:>2})  {item.verdict:<21} {_pair_text(item)}")
            lines.append(f"      why: {item.reason}")
        counts = {verdict: 0 for verdict in VERDICTS}
        for item in comparison:
            counts[item.verdict] += 1
        lines.append(
            "  verdicts: "
            + ", ".join(f"{verdict} {counts[verdict]}" for verdict in VERDICTS)
        )

    lines += ["", *_rule_of_thumb_lines()]
    return ReportSection(title=TITLE, body="\n".join(lines))


def evidence() -> tuple[Evidence, ...]:
    """Provenance of this module's rules (architecture design v0.1 section 2)."""
    return (_EVIDENCE_TABLE_S1, _EVIDENCE_FIVE_ITEMS, _EVIDENCE_REPRODUCTION)
