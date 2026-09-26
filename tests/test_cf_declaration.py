"""Projection-basis declaration check (A5) tests.

Evidence chain exercised here:

1. a transcription of the **published table**: Table S1 of the Supporting
   Information of Chilton, *Chem. Soc. Rev.* 2025, 54(24), 11468-11487 ("Ab
   initio electronic structure calculations of lanthanide single-molecule
   magnets; a practical guide"), all five columns and all 27 ``(k, q)`` rows,
   read from this project's close reading of that SI (section 4.4).  The
   transcription below is the external known answer; the module's own example
   entries and the report text are checked against it, so a transcription slip
   in the module cannot survive the suite.
2. the two measured ratios the verdicts rest on, recomputed from that
   transcription: the cross-manifold axial shift (1938 -> 1879, 3.0%) and the
   same-manifold program-to-program agreement (1938/1938, 66/66, 198 vs 201).
3. the refusal behaviour: a malformed parameter key, a duplicate ``(k, q)``, an
   unknown declaration item or a non-finite value must be refused with a
   "Next step:" sentence rather than silently compared.

Verdict rule under test (see the module docstring): one declared scheme ->
``comparable``; ``(2, 0)`` only, when no single scheme is established ->
``comparable-roughly``; anything else -> ``refused``; a ``(k, q)`` listed by one
side only -> ``missing-on-one-side``.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from fblockkit.analysis import cf_declaration
from fblockkit.knowledge import sources
from fblockkit.knowledge.models import EVIDENCE_LITERATURE

# --- Table S1 transcription --------------------------------------------------

#: Column order of Table S1 as published: the same SA-CASSCF-SO wavefunction of one
#: Dy(III) compound, projected five ways.  Values in cm^-1, exactly as printed
#: (the source table has integer entries).
TABLE_S1_COLUMNS = (
    "SINGLE_ANISO, J = 15/2",
    "SINGLE_ANISO, L = 5",
    "angmom_suite, J = 15/2",
    "angmom_suite, 6H_15/2",
    "angmom_suite, 6H",
)

TABLE_S1_ROWS: tuple[tuple[int, int, tuple[int, int, int, int, int]], ...] = (
    (2, -2, (124, 239, 125, 116, 215)),
    (2, -1, (0, -3, -1, -2, -3)),
    (2, 0, (1938, 1879, 1938, 2038, 2035)),
    (2, 1, (-6, -7, -7, -7, -7)),
    (2, 2, (198, -43, 201, 186, -38)),
    (4, -4, (0, 0, -1, 0, 0)),
    (4, -3, (0, 5, 0, 0, 6)),
    (4, -2, (5, 36, 9, 14, 27)),
    (4, -1, (7, 10, 7, 10, 10)),
    (4, 0, (66, 125, 66, 84, 84)),
    (4, 1, (4, -3, 4, 3, -3)),
    (4, 2, (8, -6, 15, 23, -5)),
    (4, 3, (5, -1, 5, 6, -1)),
    (4, 4, (0, 0, 0, 0, 0)),
    (6, -6, (0, 0, 0, 0, 0)),
    (6, -5, (0, 0, 0, 0, 0)),
    (6, -4, (0, 0, 0, 0, 0)),
    (6, -3, (-1, 3, -1, -1, 3)),
    (6, -2, (12, 19, 11, 8, 15)),
    (6, -1, (8, 8, 8, 12, 8)),
    (6, 0, (23, 21, 23, 18, 17)),
    (6, 1, (-1, -7, -1, -2, -8)),
    (6, 2, (18, -4, 18, 13, -3)),
    (6, 3, (3, 0, 3, 3, 0)),
    (6, 4, (0, 0, 0, 0, 0)),
    (6, 5, (0, 0, 0, 0, 0)),
    (6, 6, (0, 0, 0, 0, 0)),
)

#: The five declaration items as the two same-compound columns would state them.
#: ``SI_J15`` and ``SI_L5`` differ in the projection manifold only: two declared
#: schemes for one wavefunction.  ``SI_ANGMOM_J15`` carries the *same* declaration
#: as ``SI_J15`` (the caller asserting both programs used one scheme) -- which is
#: exactly what makes that pair comparable.
SI_J15_DECLARATION = {
    "convention": "Stevens (Rudowicz/Ryabov), cosine/sine tesseral",
    "projection": "J = 15/2 (Dy3+, 6H15/2)",
    "units": "cm^-1",
    "z_axis": "ground Kramers doublet magnetic axis",
    "origin": "Dy centre, right-handed frame",
}
SI_L5_DECLARATION = dict(SI_J15_DECLARATION, projection="L = 5, S = 5/2 (6H term)")
SI_ANGMOM_J15_DECLARATION = dict(SI_J15_DECLARATION)


def _column(name: str) -> dict[str, float]:
    """One Table S1 column as the payload's ``parameters`` mapping."""
    index = TABLE_S1_COLUMNS.index(name)
    return {f"{k},{q}": values[index] for k, q, values in TABLE_S1_ROWS}


def _value(key: tuple[int, int], column: str) -> int:
    index = TABLE_S1_COLUMNS.index(column)
    for k, q, values in TABLE_S1_ROWS:
        if (k, q) == key:
            return values[index]
    raise KeyError(key)


def _set(column: str, declaration: dict[str, str]) -> dict:
    return {"parameters": _column(column), "declaration": declaration}


def _verdicts(items) -> dict[tuple[int, int], str]:
    return {(item.k, item.q): item.verdict for item in items}


def _assert_next_step(excinfo: pytest.ExceptionInfo) -> str:
    """Every refusal must end in a real "Next step:" sentence (house rule)."""
    message = str(excinfo.value)
    assert "Next step:" in message, message
    tail = message[message.rindex("Next step:") :]
    assert len(tail) > len("Next step: ."), message
    assert tail.rstrip().endswith((".", ":")), message
    return message


# --- (a) a complete declaration passes ---------------------------------------


def test_complete_declaration_has_all_five_items_present():
    items = cf_declaration.check_declaration(SI_J15_DECLARATION)
    assert isinstance(items, tuple) and len(items) == 5
    assert tuple(item.field for item in items) == cf_declaration.DECLARATION_FIELDS
    assert all(item.present for item in items)
    assert [item.value for item in items] == [
        SI_J15_DECLARATION[field] for field in cf_declaration.DECLARATION_FIELDS
    ]
    assert all(item.consequence for item in items)


def test_empty_declaration_is_five_missing_items_not_an_error():
    items = cf_declaration.check_declaration({})
    assert all(not item.present and item.value == "" for item in items)
    assert all(item.consequence for item in items)


def test_blank_and_null_declaration_values_count_as_missing():
    items = cf_declaration.check_declaration(
        {"units": "   ", "z_axis": None, "origin": "Dy centre, right-handed frame"}
    )
    by_field = {item.field: item for item in items}
    assert not by_field["units"].present and by_field["units"].value == ""
    assert not by_field["z_axis"].present
    assert by_field["origin"].present and by_field["origin"].value.startswith("Dy centre")
    assert by_field["convention"].present is False


# --- (b) each missing field produces exactly one missing item ----------------


@pytest.mark.parametrize("field", cf_declaration.DECLARATION_FIELDS)
def test_each_missing_field_yields_exactly_one_missing_item(field):
    declaration = {
        key: value for key, value in SI_J15_DECLARATION.items() if key != field
    }
    items = cf_declaration.check_declaration(declaration)
    assert len(items) == 5
    missing = [item for item in items if not item.present]
    assert [item.field for item in missing] == [field]
    assert missing[0].value == ""
    assert missing[0].consequence


def test_missing_field_consequences_name_the_measured_example():
    """The projection item's consequence carries the Table S1 sign flip (why the
    numbers stop being comparable when the manifold is not declared)."""
    by_field = {item.field: item for item in cf_declaration.check_declaration({})}
    assert "+198" in by_field["projection"].consequence
    assert "-43" in by_field["projection"].consequence
    assert "q != 0" in by_field["z_axis"].consequence


# --- the measured ratios the verdicts rest on --------------------------------


def test_measured_ratios_from_the_transcribed_table():
    """Reproduction of Table S1 (this is the module's 'measured' evidence)."""
    # same manifold, different program: the large parameters agree exactly,
    # the mid-sized ones to a couple of percent
    assert _value((2, 0), "SINGLE_ANISO, J = 15/2") == _value((2, 0), "angmom_suite, J = 15/2")
    assert _value((4, 0), "SINGLE_ANISO, J = 15/2") == _value((4, 0), "angmom_suite, J = 15/2")
    transverse_j15 = _value((2, 2), "SINGLE_ANISO, J = 15/2")
    transverse_angmom = _value((2, 2), "angmom_suite, J = 15/2")
    assert 100 * abs(transverse_j15 - transverse_angmom) / transverse_j15 == pytest.approx(
        1.5, abs=0.1
    )

    # cross-manifold: the axial parameter shifts by 3%, ...
    axial_j15 = _value((2, 0), "SINGLE_ANISO, J = 15/2")
    axial_l5 = _value((2, 0), "SINGLE_ANISO, L = 5")
    assert 100 * abs(axial_j15 - axial_l5) / axial_j15 == pytest.approx(3.0, abs=0.1)
    # ... the transverse ones flip sign ...
    assert transverse_j15 == 198 and _value((2, 2), "SINGLE_ANISO, L = 5") == -43
    assert _value((2, -2), "SINGLE_ANISO, J = 15/2") == 124
    assert _value((2, -2), "SINGLE_ANISO, L = 5") == 239
    # ... and the high-order axial term moves by a factor ~1.9, which is why the
    # licence is (2, 0) and not "every q = 0"
    assert _value((4, 0), "SINGLE_ANISO, L = 5") / _value(
        (4, 0), "SINGLE_ANISO, J = 15/2"
    ) == pytest.approx(1.9, abs=0.05)
    assert _value((6, 2), "SINGLE_ANISO, J = 15/2") == 18
    assert _value((6, 2), "SINGLE_ANISO, L = 5") == -4


def test_module_example_matches_the_transcription():
    """Guard against the module's quoted example drifting from Table S1."""
    for (k, q), values in cf_declaration._TABLE_S1_EXAMPLE.items():
        expected = (
            _value((k, q), "SINGLE_ANISO, J = 15/2"),
            _value((k, q), "SINGLE_ANISO, L = 5"),
            _value((k, q), "angmom_suite, J = 15/2"),
        )
        assert tuple(values) == expected, (k, q)


# --- (c) the real Table S1 pair: J = 15/2 vs L = 5 ---------------------------


def test_table_s1_cross_manifold_verdicts():
    items = cf_declaration.compare_parameters(
        _set("SINGLE_ANISO, J = 15/2", SI_J15_DECLARATION),
        _set("SINGLE_ANISO, L = 5", SI_L5_DECLARATION),
    )
    verdicts = _verdicts(items)
    assert len(verdicts) == len(TABLE_S1_ROWS) == 27
    assert verdicts[(2, 0)] == cf_declaration.VERDICT_COMPARABLE_ROUGHLY
    assert verdicts[(2, 2)] == cf_declaration.VERDICT_REFUSED
    assert verdicts[(4, 0)] == cf_declaration.VERDICT_REFUSED
    # (2, 0) is the only licence: every other (k, q) is refused
    assert {key for key, verdict in verdicts.items() if verdict != "refused"} == {(2, 0)}

    by_key = {(item.k, item.q): item for item in items}
    axial = by_key[(2, 0)]
    assert (axial.value_a, axial.value_b) == (1938.0, 1879.0)
    assert "3.04%" in axial.reason and "provisional" not in axial.reason.lower()
    assert "not a bound" in axial.reason

    transverse = by_key[(2, 2)]
    assert (transverse.value_a, transverse.value_b) == (198.0, -43.0)
    assert "sign flip" in transverse.reason
    assert "Next step:" in transverse.reason

    high_order = by_key[(4, 0)]
    assert (high_order.value_a, high_order.value_b) == (66.0, 125.0)
    assert "high-order" in high_order.reason
    assert "a factor 1.9" in high_order.reason  # its own measured move, not a generic line

    # a (k, q) the table lists directly is refused with its own measured pair
    assert "18" in by_key[(6, 2)].reason and "-4" in by_key[(6, 2)].reason
    assert "sign flip" in by_key[(6, 2)].reason


def test_cross_manifold_verdicts_hold_when_both_declarations_are_incomplete():
    """Same verdicts when neither side declares the scheme (the measured case)."""
    items = cf_declaration.compare_parameters(
        _set("SINGLE_ANISO, J = 15/2", {"projection": "J = 15/2"}),
        _set("SINGLE_ANISO, L = 5", {"projection": "L = 5"}),
    )
    verdicts = _verdicts(items)
    assert verdicts[(2, 0)] == cf_declaration.VERDICT_COMPARABLE_ROUGHLY
    assert verdicts[(2, 2)] == cf_declaration.VERDICT_REFUSED
    assert verdicts[(4, 0)] == cf_declaration.VERDICT_REFUSED
    assert set(verdicts.values()) == {"refused", "comparable-roughly"}


def test_identical_but_incomplete_declarations_establish_nothing():
    """Agreeing on two items is not a scheme: only (2, 0) may then be compared."""
    items = cf_declaration.compare_parameters(
        _set("SINGLE_ANISO, J = 15/2", {"units": "cm^-1", "projection": "J = 15/2"}),
        _set("angmom_suite, J = 15/2", {"units": "cm^-1", "projection": "J = 15/2"}),
    )
    verdicts = _verdicts(items)
    assert verdicts[(2, 0)] == cf_declaration.VERDICT_COMPARABLE_ROUGHLY
    assert verdicts[(2, 2)] == cf_declaration.VERDICT_REFUSED
    assert "not stated on both sides: convention, z_axis, origin" in next(
        item.reason for item in items if (item.k, item.q) == (2, 0)
    )


# --- (d) same manifold, different program ------------------------------------


def test_same_manifold_different_program_is_comparable():
    """SINGLE_ANISO vs angmom_suite, both J = 15/2: 1938/1938, 198/201, 66/66."""
    items = cf_declaration.compare_parameters(
        _set("SINGLE_ANISO, J = 15/2", SI_J15_DECLARATION),
        _set("angmom_suite, J = 15/2", SI_ANGMOM_J15_DECLARATION),
    )
    assert len(items) == 27
    assert all(item.verdict == cf_declaration.VERDICT_COMPARABLE for item in items)
    by_key = {(item.k, item.q): item for item in items}
    assert (by_key[(2, 0)].value_a, by_key[(2, 0)].value_b) == (1938.0, 1938.0)
    assert (by_key[(2, 2)].value_a, by_key[(2, 2)].value_b) == (198.0, 201.0)
    assert (by_key[(4, 0)].value_a, by_key[(4, 0)].value_b) == (66.0, 66.0)
    # "comparable" does not require equal numbers: 198 vs 201 is a program
    # difference inside one declared scheme, and the verdict says so
    assert "different programs" in by_key[(2, 2)].reason


def test_declaration_difference_is_what_blocks_the_comparison():
    """Same numbers, one item of the declaration changed -> the verdict changes."""
    parameters = _column("SINGLE_ANISO, J = 15/2")
    same = cf_declaration.compare_parameters(
        {"parameters": parameters, "declaration": SI_J15_DECLARATION},
        {"parameters": parameters, "declaration": dict(SI_J15_DECLARATION)},
    )
    assert all(item.verdict == cf_declaration.VERDICT_COMPARABLE for item in same)
    shifted = cf_declaration.compare_parameters(
        {"parameters": parameters, "declaration": SI_J15_DECLARATION},
        {
            "parameters": parameters,
            "declaration": dict(
                SI_J15_DECLARATION, z_axis="crystallographic c axis"
            ),
        },
    )
    verdicts = _verdicts(shifted)
    assert verdicts[(2, 0)] == cf_declaration.VERDICT_COMPARABLE_ROUGHLY
    assert verdicts[(4, 0)] == cf_declaration.VERDICT_REFUSED
    assert "stated differently: z_axis" in next(
        item.reason for item in shifted if (item.k, item.q) == (2, 0)
    )


# --- missing on one side -----------------------------------------------------


def test_parameter_listed_by_one_side_only_is_reported_as_such():
    items = cf_declaration.compare_parameters(
        {"parameters": {"2,0": 1938.0, "2,2": 198.0}, "declaration": SI_J15_DECLARATION},
        {"parameters": {"2,0": 1879.0}, "declaration": SI_L5_DECLARATION},
    )
    assert [(item.k, item.q) for item in items] == [(2, 0), (2, 2)]  # ordered by (k, q)
    missing = [item for item in items if item.verdict == cf_declaration.VERDICT_MISSING_ON_ONE_SIDE]
    assert len(missing) == 1
    item = missing[0]
    assert (item.k, item.q) == (2, 2)
    assert item.value_a == 198.0 and item.value_b is None
    assert "the first set (A)" in item.reason
    assert "Next step:" in item.reason

    reversed_items = cf_declaration.compare_parameters(
        {"parameters": {"2,0": 1938.0}, "declaration": SI_J15_DECLARATION},
        {"parameters": {"2,0": 1879.0, "6,6": 0.0}, "declaration": SI_L5_DECLARATION},
    )
    item = next(
        entry
        for entry in reversed_items
        if entry.verdict == cf_declaration.VERDICT_MISSING_ON_ONE_SIDE
    )
    assert (item.k, item.q) == (6, 6)
    assert item.value_a is None and item.value_b == 0.0
    assert "the second set (B)" in item.reason


# --- (e) malformed input is refused with a Next step -------------------------

BAD_KEYS = (
    "2",
    "2,0,1",
    "k,q",
    "2.5,0",
    "[2 0]",
    "",
    "two,zero",
    (2, 0, 1),
    (2,),
    (2, "x"),
    3.5,
    None,
    True,
)


@pytest.mark.parametrize("key", BAD_KEYS)
def test_malformed_parameter_key_is_refused(key):
    with pytest.raises(cf_declaration.DeclarationError) as excinfo:
        cf_declaration.run({"parameters": {key: 1.0}})
    _assert_next_step(excinfo)
    with pytest.raises(cf_declaration.DeclarationError) as excinfo:
        cf_declaration.compare_parameters(
            {"parameters": {key: 1.0}}, {"parameters": {"2,0": 1.0}}
        )
    _assert_next_step(excinfo)


@pytest.mark.parametrize("key", ("3,1", "0,0", "8,0", "2,3", "-2,0", "2,-3"))
def test_out_of_range_parameter_index_is_refused(key):
    with pytest.raises(cf_declaration.DeclarationError, match="Next step:"):
        cf_declaration.run({"parameters": {key: 1.0}})


def test_accepted_key_spellings_agree():
    """The pinned spellings ("2,0", "[2, 2]") and the pair form must agree."""
    spellings = (
        {"2,0": 1938.0, "2,2": 198.0},
        {"[2, 0]": 1938.0, "[2, 2]": 198.0},
        {"2, 0": 1938.0, "2, 2": 198.0},
        {(2, 0): 1938.0, (2, 2): 198.0},
        {("[2,0]"): 1938.0, "2,2": 198.0},
    )
    bodies = {
        cf_declaration.run(
            {"parameters": spelling, "declaration": SI_J15_DECLARATION}
        ).body
        for spelling in spellings
    }
    assert len(bodies) == 1


BAD_PAYLOADS = (
    {},
    {"declaration": {"units": "cm^-1"}},
    {"parameters": "2,0"},
    {"parameters": {"2,0": 1.0}, "extra": 1},
    {"parameters": {"2,0": 1.0}, "declaration": {"normalisation": "Rudowicz"}},
    {"parameters": {"2,0": 1.0}, "declaration": {"units": 1}},
    {"parameters": {"2,0": 1.0}, "declaration": ["units"]},
    {"parameters": {"2,0": 1.0}, "label": 7},
    {"parameters": {"2,0": 1.0}, "compare": 5},
    {"parameters": {"2,0": 1.0}, "compare": {"declaration": {}}},
    {"parameters": {"2,0": 1.0}, "compare": {"parameters": {"2,0": 1.0}, "zzz": 1}},
    {"parameters": {"2,0": 1.0}, "compare": {"parameters": {"2,0": 1.0}, "label": ""}},
    {"parameters": {"2,0": 1.0, "[2, 0]": 2.0}},
    {"parameters": {"2,0": "1938"}},
    {"parameters": {"2,0": float("nan")}},
    {"parameters": {"2,0": float("inf")}},
    {"parameters": {"2,0": None}},
)


@pytest.mark.parametrize("payload", BAD_PAYLOADS)
def test_malformed_payload_is_refused(payload):
    with pytest.raises(cf_declaration.DeclarationError) as excinfo:
        cf_declaration.run(payload)
    _assert_next_step(excinfo)


def test_non_mapping_payload_and_set_are_refused():
    with pytest.raises(cf_declaration.DeclarationError, match="Next step:"):
        cf_declaration.run(["parameters"])
    with pytest.raises(cf_declaration.DeclarationError, match="Next step:"):
        cf_declaration.compare_parameters(None, {"parameters": {}})
    with pytest.raises(cf_declaration.DeclarationError, match="Next step:"):
        cf_declaration.check_declaration("units: cm^-1")


def test_declaration_error_is_a_value_error():
    assert issubclass(cf_declaration.DeclarationError, ValueError)


# --- report ------------------------------------------------------------------


def _payload(column: str = "SINGLE_ANISO, J = 15/2") -> dict:
    return {
        "parameters": _column(column),
        "declaration": SI_J15_DECLARATION,
        "label": "this set (SINGLE_ANISO, J = 15/2)",
        "compare": {
            "label": "SINGLE_ANISO, L = 5",
            "parameters": _column("SINGLE_ANISO, L = 5"),
            "declaration": SI_L5_DECLARATION,
        },
    }


def test_run_title_and_declaration_table():
    section = cf_declaration.run(_payload())
    assert section.title == "A5 projection-basis declaration check"
    assert "A5 projection-basis declaration check" in cf_declaration.TITLE
    for field in cf_declaration.DECLARATION_FIELDS:
        assert field in section.body
    assert "stated: 5 of 5 items" in section.body
    assert "convention " in section.body


def test_run_lists_comparison_verdicts():
    section = cf_declaration.run(_payload())
    assert "Comparison: this set (SINGLE_ANISO, J = 15/2) (A) vs SINGLE_ANISO, L = 5 (B)" in (
        section.body
    )
    assert "(2, 0)  comparable-roughly" in section.body
    assert "(2, 2)  refused" in section.body
    assert "verdicts: comparable 0, comparable-roughly 1, refused 26, missing-on-one-side 0" in (
        section.body
    )


def test_run_reports_a_missing_declaration_item():
    payload = _payload()
    payload["declaration"] = {"units": "cm^-1", "projection": "J = 15/2"}
    section = cf_declaration.run(payload)
    assert "projection present" in section.body
    assert "stated: 2 of 5 items" in section.body
    assert "1. convention MISSING" in section.body
    assert "consequence: the operator convention" in section.body


def test_run_without_compare_says_nothing_was_compared():
    section = cf_declaration.run(
        {"parameters": _column("SINGLE_ANISO, J = 15/2"), "declaration": SI_J15_DECLARATION}
    )
    assert "Set: this set" in section.body  # default label
    assert "Comparison: none requested" in section.body
    assert "1938" in section.body  # the rule-of-thumb line is always printed


def test_run_quotes_the_table_s1_numbers_and_rule():
    body = cf_declaration.run(_payload()).body
    for text in ("1938", "1879", "198", "-43", "124", "239", "66", "125", "18", "-4"):
        assert text in body, text
    assert "3.04%" in body
    assert "Table S1" in body
    assert "11468-11487" in body
    assert "same manifold, different program" in body
    assert "not comparable across schemes" in body


def test_run_is_deterministic_and_does_not_mutate_the_payload():
    payload = _payload()
    snapshot = repr(payload)
    first = cf_declaration.run(payload)
    second = cf_declaration.run(payload)
    assert first.body == second.body
    assert repr(payload) == snapshot


def test_compare_parameters_is_deterministic_and_ordered():
    items = cf_declaration.compare_parameters(
        _set("SINGLE_ANISO, J = 15/2", SI_J15_DECLARATION),
        _set("SINGLE_ANISO, L = 5", SI_L5_DECLARATION),
    )
    assert [item.verdict for item in items] == [
        item.verdict
        for item in cf_declaration.compare_parameters(
            _set("SINGLE_ANISO, J = 15/2", SI_J15_DECLARATION),
            _set("SINGLE_ANISO, L = 5", SI_L5_DECLARATION),
        )
    ]
    keys = [(item.k, item.q) for item in items]
    assert keys == sorted(keys)
    assert set(keys) == {(k, q) for k, q, _ in TABLE_S1_ROWS}


def test_module_keeps_the_layer_boundary():
    """analysis may only reach down into knowledge (import-linter contract)."""
    source = Path(cf_declaration.__file__).read_text(encoding="utf-8")
    for forbidden in (
        "fblockkit.recipe",
        "fblockkit.diagnosis",
        "fblockkit.ui",
        "fblockkit.toolindex",
    ):
        assert forbidden not in source


# --- (f) provenance ----------------------------------------------------------


def test_every_evidence_bibkey_exists_in_sources_bib():
    items = cf_declaration.evidence()
    assert items
    kinds = {item.kind for item in items}
    assert kinds == {"literature", "measured"}
    literature = [item for item in items if item.kind == EVIDENCE_LITERATURE]
    assert {item.bibkey for item in literature} == {"chilton2025abinitio", "peng2025accurate"}
    for item in literature:
        entry = sources.get(item.bibkey)  # raises BibDataError for an unknown key
        assert entry.field("doi") in item.ref  # complete citation, as elsewhere
        assert item.url.startswith("https://doi.org/")
    for item in items:
        assert item.text and item.ref
        if item.kind != EVIDENCE_LITERATURE:
            assert not item.bibkey
    measured = [item for item in items if item.kind == "measured"]
    assert len(measured) == 1
    assert "tests/test_cf_declaration.py" in measured[0].ref
