"""G3 tests: OpenMolcas magnetic-chain input generation (SA-CASSCF / RASSI / SINGLE_ANISO).

Discipline under test: the rendered file must be usable as it stands (module chain, values
and basis assignment from the source Supporting Information), deterministic, ASCII, free of
placeholder text, and every line that the source does not show must be tagged
``[fBlockKit:...]`` in the file itself. Every validation error must end with a
"Next step: ..." sentence.
"""

from __future__ import annotations

import math
import re

import pytest

from fblockkit.knowledge.loader import bib_keys
from fblockkit.recipe.openmolcas import (
    SI_BIBKEY,
    OpenMolcasChainSpec,
    OpenMolcasSpecError,
    evidence,
    max_csfs,
    own_deviations,
    render_openmolcas_input,
    run_guidance,
    run_guidance_openmolcas,
    spin_orbit_state_count,
)
from fblockkit.recipe.render import RenderError


def _body(text: str) -> str:
    """Rendered file without the '*' comment lines (the module input proper)."""
    return "\n".join(
        line for line in text.splitlines() if not line.lstrip().startswith("*")
    )


def _spec(**kwargs) -> OpenMolcasChainSpec:
    base = dict(
        element="Dy",
        charge=1,
        nactel=9,
        nactorb=7,
        roots={6: 21, 4: 48, 2: 32},
        coord_file="${CurrDir}/DyBrOH_basis.xyz",
    )
    base.update(kwargs)
    return OpenMolcasChainSpec(**base)


# --- structure of the rendered file --------------------------------------


def test_module_chain_follows_the_source_skeleton():
    text = render_openmolcas_input(_spec())
    modules = re.findall(r"^&([A-Z_]+)", text, re.M)
    assert modules == [
        "GATEWAY",
        "SEWARD",
        "RASSCF",
        "RASSCF",
        "RASSCF",
        "RASSI",
        "SINGLE_ANISO",
    ]


def test_spin_blocks_are_descending_and_chained():
    body = _body(render_openmolcas_input(_spec()))
    assert body.index("Spin= 6") < body.index("Spin= 4") < body.index("Spin= 2")
    # Only the first block is computed from scratch; the others read its orbital file.
    assert body.count("Typeindex") == 2
    assert body.count("FILEORB= ") == 2
    assert body.count("RAS2=") == 1
    assert "FILEORB= 1.RasOrb" in body and "FILEORB= 2.RasOrb" in body
    # CiRoot is written twice with unity weighting, as in the source.
    assert "CiRoot= 21 21 1" in body
    assert "Nactel= 9 0 0" in body
    assert "RAS2= 7" in body


def test_last_block_does_not_copy_an_orbital_file():
    body = _body(render_openmolcas_input(_spec()))
    assert ">> COPY $Project.RasOrb 1.RasOrb" in body
    assert ">> COPY $Project.RasOrb 2.RasOrb" in body
    assert ">> COPY $Project.RasOrb 3.RasOrb" not in body


def test_two_block_spec_still_renders_one_chain():
    text = render_openmolcas_input(_spec(roots={6: 21, 2: 32}))
    body = _body(text)
    assert re.findall(r"^&([A-Z_]+)", body, re.M) == [
        "GATEWAY",
        "SEWARD",
        "RASSCF",
        "RASSCF",
        "RASSI",
        "SINGLE_ANISO",
    ]
    assert "Nr of JobIph= 2 ALL" in body
    assert "'ANGMOM' 1" in body and "'ANGMOM' 2" in body and "'ANGMOM' 3" not in body


def test_basis_assignment_is_present():
    spec = _spec()
    text = render_openmolcas_input(spec)
    for label in (spec.basis_metal, spec.basis_ligand, spec.basis_outer, spec.element):
        assert label in text, label
    assert spec.basis_metal == "ANO-RCC-VTZP"
    assert spec.basis_ligand == "ANO-RCC-VDZP"
    assert spec.basis_outer == "ANO-RCC-VDZ"
    assert "Coord= ${CurrDir}/DyBrOH_basis.xyz" in text


def test_rendered_file_is_ascii_and_has_no_placeholder_text():
    text = render_openmolcas_input(_spec())
    assert text.isascii()
    assert "{{" not in text and "}}" not in text
    assert "None" not in text and "TODO" not in text


def test_render_is_deterministic():
    assert render_openmolcas_input(_spec()) == render_openmolcas_input(_spec())
    # Normalisation sorts the spin blocks, so the order given does not matter.
    assert render_openmolcas_input(_spec(roots={2: 32, 6: 21, 4: 48})) == render_openmolcas_input(
        _spec()
    )
    assert render_openmolcas_input(_spec(roots=((4, 48), (6, 21), (2, 32)))) == (
        render_openmolcas_input(_spec())
    )


def test_spec_is_normalised_and_hashable():
    spec = _spec(element="dy", roots={2: 32, 6: 21, 4: 48})
    assert spec.element == "Dy"
    assert spec.spin_blocks == ((6, 21), (4, 48), (2, 32))
    assert hash(spec) == hash(_spec())
    assert "CRYS= dy" in render_openmolcas_input(spec)


def test_source_choices_are_carried_into_the_file():
    text = _body(render_openmolcas_input(_spec()))
    for token in (
        "AMFI",
        "RICD=acCD",
        "Angmom= 0.0 0.0 0.0 ANGSTROM",
        "EPRG= 7.0D-1",
        "MEES",
        "SPIN",
        "MLTP= 1; 2",
        "TINT= 0.0 330.0 330 0.0001",
        "HINT= 0.0 10.0 201",
        "TMAG= 6 1.8 2 4 5 10 20",
        "QUAX= 1",
    ):
        assert token in text, token


# --- honesty about what is not in the source -----------------------------


def test_every_own_deviation_is_marked_in_the_file():
    text = render_openmolcas_input(_spec())
    deviations = own_deviations()
    assert deviations
    for deviation in deviations:
        assert deviation.text.strip()
        assert f"[fBlockKit:{deviation.tag}]" in text, deviation.tag


def test_evidence_points_at_the_source_supporting_information():
    items = evidence()
    assert items
    keys = bib_keys()
    for item in items:
        assert item.kind == "literature"
        assert item.bibkey == SI_BIBKEY
        assert item.bibkey in keys
        assert "Supporting Information" in item.ref
        assert item.url.startswith("https://doi.org/")


# --- validation -----------------------------------------------------------


def _assert_next_step(excinfo: pytest.ExceptionInfo[OpenMolcasSpecError]) -> str:
    message = str(excinfo.value)
    assert "Next step:" in message, message
    tail = message[message.rindex("Next step:") :]
    assert tail.endswith("."), message
    assert len(tail) > len("Next step: ."), message
    return message


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"element": ""}, "non-empty"),
        ({"element": "Xx"}, "not a known element symbol"),
        ({"element": "Fe"}, "not a lanthanide or actinide"),
        ({"charge": 1.5}, "must be an integer"),
        ({"charge": True}, "must be an integer"),
        ({"charge": 12}, "sanity range"),
        ({"nactel": 0}, "nactel must be an integer >= 1"),
        ({"nactorb": 0}, "nactorb must be an integer >= 1"),
        ({"nactel": 20, "nactorb": 7}, "cannot be hosted"),
        ({"nactorb": 5}, "seven f orbitals"),
        ({"roots": {}}, "roots is empty"),
        ({"roots": "sextets"}, "roots must be a mapping"),
        ({"roots": ((6, 21, 1),)}, "not a (spin multiplicity, nroots) pair"),
        ({"roots": {0: 5}}, "is not a spin multiplicity"),
        ({"roots": {6: 0}}, "positive number of CI roots"),
        ({"roots": ((6, 21), (6, 21))}, "appears twice"),
        ({"roots": {6: 22}}, "contains only 21"),
        ({"roots": {5: 1}}, "not reachable"),
        ({"roots": {8: 1}}, "not reachable"),
        ({"coord_file": "my file.xyz"}, "whitespace"),
        # escaped so that this test file stays ASCII while testing the ASCII guard
        ({"coord_file": "Dyé_basis.xyz"}, "non-ASCII"),
        ({"coord_file": ""}, "coord_file must be a non-empty string"),
        ({"basis_metal": ""}, "basis_metal"),
        ({"basis_ligand": "ANO RCC"}, "whitespace"),
        ({"symmetry": "2"}, "point-group token"),
        ({"title": "two\nlines"}, "whitespace"),
        ({"angmom_origin": (0.0, 0.0)}, "three numbers"),
        ({"angmom_origin": (0.0, 0.0, float("nan"))}, "non-finite"),
    ],
)
def test_validation_errors_carry_a_next_step(kwargs, match):
    with pytest.raises(OpenMolcasSpecError) as excinfo:
        _spec(**kwargs)
    message = _assert_next_step(excinfo)
    assert match in message, message


def test_error_type_is_a_render_error():
    assert issubclass(OpenMolcasSpecError, RenderError)
    assert issubclass(OpenMolcasSpecError, ValueError)


def test_accepted_boundary_values():
    spec = _spec(charge=-1, nactel=7, nactorb=7, roots={8: 1, 6: 7, 4: 20, 2: 35})
    assert spec.nactel == 7 and spec.nactorb == 7
    assert spec.spin_blocks[0] == (8, 1)


# --- CI-space counting (cross-checked against the source numbers) --------


def test_csf_counts_reproduce_the_source_numbers():
    assert max_csfs(9, 7, 6) == 21  # sextets (source section 4.2)
    assert max_csfs(9, 7, 4) == 224  # quartets
    assert max_csfs(9, 7, 2) == 490  # doublets
    assert max_csfs(2, 2, 1) == 3  # two electrons in two orbitals: three singlets
    assert max_csfs(2, 2, 3) == 1
    assert max_csfs(9, 7, 8) == 0  # unreachable in seven orbitals
    assert max_csfs(9, 7, 1) == 0  # parity
    assert max_csfs(-1, 7, 2) == 0
    assert max_csfs(9, 0, 2) == 0


def test_csf_counts_match_an_independent_determinant_count():
    """The binomial count is checked against the branching formula N_S = D(S) - D(S+1).

    ``D(S)`` = number of determinants with M_S = S = C(m, n/2+S) * C(m, n/2-S). This uses
    a different route than the Weyl-Palmer formula, so it catches a wrong implementation.
    """

    def branching(nactel: int, nactorb: int, multiplicity: int) -> int:
        def determinants(two_m_s: int) -> int:
            if (nactel + two_m_s) % 2:
                return 0
            alpha = (nactel + two_m_s) // 2
            beta = (nactel - two_m_s) // 2
            if alpha < 0 or beta < 0 or alpha > nactorb or beta > nactorb:
                return 0
            return math.comb(nactorb, alpha) * math.comb(nactorb, beta)

        two_s = multiplicity - 1
        if two_s < 0:
            return 0
        return determinants(two_s) - determinants(two_s + 2)

    for nactel in range(0, 11):
        for nactorb in range(1, 9):
            for multiplicity in range(1, 14):
                assert max_csfs(nactel, nactorb, multiplicity) == branching(
                    nactel, nactorb, multiplicity
                ), (nactel, nactorb, multiplicity)


def test_spin_orbit_state_count_matches_the_source_example():
    assert spin_orbit_state_count(_spec()) == 382  # 6x21 + 4x48 + 2x32
    assert spin_orbit_state_count(_spec(roots={6: 21})) == 126


# --- run guidance ---------------------------------------------------------


def test_package_exports_do_not_shadow_the_orca_guide():
    """``recipe.run_guidance`` stays the ORCA guide; the G3 guide has its own name."""
    from fblockkit.recipe import OpenMolcasChainSpec as package_spec
    from fblockkit.recipe import render_openmolcas_input as package_render
    from fblockkit.recipe import run_guidance as package_run_guidance
    from fblockkit.recipe import run_guidance_openmolcas
    from fblockkit.recipe.render import run_guidance as orca_guide

    assert package_run_guidance is orca_guide
    assert run_guidance_openmolcas is not orca_guide
    assert package_spec is OpenMolcasChainSpec
    assert package_render is render_openmolcas_input


def test_run_guidance_is_honest_about_parsing_and_lists_hand_checks():
    spec = _spec()
    text = run_guidance_openmolcas(spec)
    assert run_guidance(spec) == text  # the name used by the feature specification
    assert "does not parse OpenMolcas output" in text
    assert "382" in text  # the spin-orbit state arithmetic of this specification
    for topic in ("RASSCF", "RASSI", "SINGLE_ANISO", "g_3", "theta_3"):
        assert topic in text, topic
    assert spec.coord_file in text
    assert spec.basis_metal in text and spec.basis_outer in text
