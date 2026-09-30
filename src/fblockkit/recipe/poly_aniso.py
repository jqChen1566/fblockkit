"""G8: the POLY_ANISO input writer -- the exchange-cluster driver of a magnetic cluster.

The read side (menu 39's first mode) turns an ``otool_poly_aniso`` output
into a report; this module writes the *input* that produces it, so a cluster
workflow is complete inside the toolkit: per-centre SINGLE_ANISO runs first
(menu 36's chain), then this input, then the tool, then the read side.

The input format is the ORCA manual's section 7.18 (measured against the
shipped probe, ``fixtures/poly_aniso/``):

- ``NNEQ``: the number of non-equivalent magnetic centre types, with the
  letter ``T``/``F`` answering whether all of them were computed ab initio;
  the next line carries the number of *equivalent centres* of each type (one
  number per type, one line), the line after the number of low-lying
  spin-orbit functions from each centre that form the local exchange basis
  (the exchange-space size is the product of the two);
- ``PAIR``: the coupled site pairs with their ``J`` (cm-1), under the
  Lines-type Hamiltonian ``sum -J s_i . s_j`` over the pairs;
- ``COOR`` (optional): one symmetrised Cartesian coordinate per centre type
  (Angstrom); activating it adds the exact dipole-dipole coupling for the
  declared pairs;
- ``TINT`` (optional): the susceptibility sampling grid ``t_min t_max
  n_points``.

What the module refuses: a type count outside 1..6, a centre count or
spin-orbit count below its structural floor, a pair index outside the
cluster, a duplicated pair, a coordinate array that does not match the type
count, and a temperature grid that is not increasing.  The *physics* of the
``J`` values and the coordinates is the caller's (the writer does not fit or
judge them; the source's model is the Lines exchange plus exact
dipole-dipole coupling).
"""

from __future__ import annotations

from dataclasses import dataclass

from ..knowledge.models import EVIDENCE_MEASURED, EVIDENCE_MANUAL, Evidence

__all__ = [
    "PolyAnisoPlanError",
    "PolyAnisoPlan",
    "build_input",
    "render_plan",
    "evidence",
]

#: The largest non-equivalent-type count the driver accepts (measured in the
#: wave-56 scouting; the manual's examples use one to three).
MAX_TYPES = 6


class PolyAnisoPlanError(ValueError):
    """The input plan cannot be written (with a next step)."""


@dataclass(frozen=True)
class PolyAnisoPlan:
    """The validated plan and the generated input text."""

    equivalent_centres: tuple[int, ...]
    spin_orbit_states: tuple[int, ...]
    pairs: tuple[tuple[int, int, float], ...]
    coordinates: tuple[tuple[float, float, float], ...] | None
    temperature_grid: tuple[float, float, int] | None
    text: str

    @property
    def total_centres(self) -> int:
        return sum(self.equivalent_centres)

    @property
    def exchange_basis(self) -> int:
        product = 1
        for count in self.spin_orbit_states:
            product *= count
        return product


def build_input(
    *,
    equivalent_centres,
    spin_orbit_states,
    pairs,
    coordinates=None,
    temperature_grid=None,
) -> PolyAnisoPlan:
    """Validate the cluster description and render the ``&POLY_ANISO`` input.

    ``equivalent_centres`` / ``spin_orbit_states``: one number per centre type
    (ints); ``pairs``: ``(i, j, J)`` with 1-based site indices over all
    centres and J in cm-1; ``coordinates``: ``(x, y, z)`` per type in Angstrom
    (``None`` omits the COOR block); ``temperature_grid``: ``(t_min, t_max,
    n_points)`` (``None`` omits the TINT block).
    """
    centres = tuple(int(value) for value in equivalent_centres)
    states = tuple(int(value) for value in spin_orbit_states)
    if not centres or len(centres) != len(states):
        raise PolyAnisoPlanError(
            "the centre-count line and the spin-orbit-count line must be given "
            "together, one number per centre type. Next step: give both lines."
        )
    if len(centres) > MAX_TYPES:
        raise PolyAnisoPlanError(
            f"the plan has {len(centres)} non-equivalent centre types while the "
            f"driver covers up to {MAX_TYPES}. Next step: group symmetry-equivalent "
            "centres under one type (their data file is shared)."
        )
    for count in centres:
        if count < 1:
            raise PolyAnisoPlanError(
                f"a centre type carries {count} equivalent centres. Next step: give "
                "at least one centre per type."
            )
    for count in states:
        if count < 2:
            raise PolyAnisoPlanError(
                f"a centre contributes {count} spin-orbit function(s) to the exchange "
                "basis; an exchange basis needs at least two (a pseudospin). Next "
                "step: give the number of low-lying spin-orbit functions per centre."
            )
    total = sum(centres)
    seen = set()
    checked_pairs = []
    if not pairs:
        raise PolyAnisoPlanError(
            "no coupled pairs were given. Next step: give at least one 'i j J' pair "
            "(the driver requires the PAIR block)."
        )
    for pair in pairs:
        first, second, coupling = pair
        first, second = int(first), int(second)
        if not (1 <= first <= total and 1 <= second <= total):
            raise PolyAnisoPlanError(
                f"the pair {first}-{second} is outside the cluster's {total} centres "
                "(site indices start at 1). Next step: check the pair list."
            )
        if first == second:
            raise PolyAnisoPlanError(
                f"the pair {first}-{second} couples a centre with itself. Next step: "
                "check the pair list."
            )
        key = (min(first, second), max(first, second))
        if key in seen:
            raise PolyAnisoPlanError(
                f"the pair {first}-{second} appears twice. Next step: keep one line "
                "per pair."
            )
        seen.add(key)
        checked_pairs.append((first, second, float(coupling)))
    coordinates_checked = None
    if coordinates is not None:
        coordinates_checked = tuple(
            tuple(float(value) for value in row) for row in coordinates
        )
        if len(coordinates_checked) != len(centres):
            raise PolyAnisoPlanError(
                f"the COOR block carries {len(coordinates_checked)} coordinate line(s) "
                f"for {len(centres)} centre type(s). Next step: give one 'x y z' line "
                "per type."
            )
    grid_checked = None
    if temperature_grid is not None:
        minimum, maximum, points = temperature_grid
        minimum, maximum, points = float(minimum), float(maximum), int(points)
        if not (0.0 <= minimum < maximum) or points < 2:
            raise PolyAnisoPlanError(
                f"the temperature grid {minimum}..{maximum} K with {points} points is "
                "not increasing or too short. Next step: give 't_min t_max n_points' "
                "with 0 <= t_min < t_max and at least two points."
            )
        grid_checked = (minimum, maximum, points)

    def number(value: float) -> str:
        return f"{value:g}"

    lines = ["&POLY_ANISO", "", "NNEQ"]
    lines.append(f"  {len(centres)}  T")
    lines.append("  " + "  ".join(str(count) for count in centres))
    lines.append("  " + "  ".join(str(count) for count in states))
    lines.append("")
    lines.append("PAIR")
    lines.append(f"  {len(checked_pairs)}")
    for first, second, coupling in checked_pairs:
        lines.append(f"  {first} {second} {number(coupling)}")
    if coordinates_checked is not None:
        lines.append("")
        lines.append("COOR")
        for x, y, z in coordinates_checked:
            lines.append(f"  {number(x)} {number(y)} {number(z)}")
    if grid_checked is not None:
        lines.append("")
        lines.append("TINT")
        lines.append(f"  {number(grid_checked[0])} {number(grid_checked[1])} {grid_checked[2]}")
    lines.append("")
    lines.append("End of Input")
    lines.append("")
    return PolyAnisoPlan(
        equivalent_centres=centres,
        spin_orbit_states=states,
        pairs=tuple(checked_pairs),
        coordinates=coordinates_checked,
        temperature_grid=grid_checked,
        text="\n".join(lines),
    )


def render_plan(plan: PolyAnisoPlan, *, path) -> str:
    """The checklist a cluster workflow needs, as the menu prints it."""
    lines = [
        "POLY_ANISO exchange-cluster input (Lines-type exchange + exact "
        "dipole-dipole):",
        f"  centre types: {len(plan.equivalent_centres)} "
        f"(equivalent centres per type: "
        f"{', '.join(str(count) for count in plan.equivalent_centres)})",
        f"  total centres: {plan.total_centres}",
        f"  spin-orbit functions per centre type: "
        f"{', '.join(str(count) for count in plan.spin_orbit_states)}"
        f"  (exchange basis size {plan.exchange_basis})",
        f"  coupled pairs: {len(plan.pairs)} "
        + "(J in cm-1: "
        + ", ".join(f"{i}-{j} {coupling:g}" for i, j, coupling in plan.pairs)
        + ")",
        "  dipole-dipole (COOR): "
        + ("included" if plan.coordinates is not None else "not included"),
        "  temperature grid: "
        + (
            f"{plan.temperature_grid[0]:g}..{plan.temperature_grid[1]:g} K, "
            f"{plan.temperature_grid[2]} points"
            if plan.temperature_grid is not None
            else "not included"
        ),
        "",
        f"Input written: {path}",
        "",
        "Next steps:",
        "  1. place one SINGLE_ANISO data file per magnetic centre next to this "
        "input, named aniso_1.input, aniso_2.input, ... (each is the "
        "<job>.CASSCF.anisofile of a menu-36 chain run; symmetry-equivalent "
        "centres of one type share a file)",
        f"  2. run:  otool_poly_aniso < {path} > poly_aniso.output",
        "  3. read the result back with this menu's first mode (the "
        "poly_aniso.output report)",
        "",
        "Boundaries:",
        "  - the J values and the coordinates are the caller's (the writer does "
        "not fit or judge them); the exchange model is the Lines-type "
        "Hamiltonian sum -J s_i.s_j over the declared pairs",
        "  - the COOR block computes the dipolar coupling only for the declared "
        "pairs; coordinates are the symmetrised per-type positions in Angstrom",
        "  - the plan covers the manual's NNEQ/PAIR/COOR/TINT blocks; the "
        "symmetry (SYMM) and anisotropic-coupling (LIN3/LIN9) variants are "
        "registered",
    ]
    return "\n".join(lines)


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the input format and of its acceptance run."""
    return (
        Evidence(
            kind=EVIDENCE_MANUAL,
            text=(
                "The POLY_ANISO input format (NNEQ with the ab initio flag and the "
                "per-type counts of equivalent centres and low-lying spin-orbit "
                "functions, the PAIR list under the Lines-type Hamiltonian, the "
                "optional COOR block for the exact dipole-dipole coupling, and the "
                "optional TINT susceptibility grid) is section 7.18 of the ORCA "
                "manual; the per-centre data files are the SINGLE_ANISO outputs "
                "placed as aniso_1.input, aniso_2.input, ..."
            ),
            ref="ORCA 6 manual section 7.18 (Interface to POLY_ANISO module)",
            url="https://www.faccts.de/docs/orca/6.0/manual/contents/detailed/poly_aniso.html",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "The generated input is accepted end to end, measured on 6.1.1 "
                "(otool_poly_aniso v1.0.0, 2026-09-30): the plan reproducing the "
                "fixture's two-centre probe (2 types x 1 centre, spin-orbit basis "
                "2+2, one J = 0.1 cm-1 pair, the two coordinates, TINT 0..300 K / "
                "101 points) wrote an input whose run with the fixture's "
                "aniso_1/aniso_2 files returns rc = 0, 'POLY_ANISO finished "
                "sucessfully!', and reproduces the frozen "
                "fixtures/poly_aniso/two_center_probe.out byte for byte."
            ),
            ref="tests/test_poly_aniso_recipe.py; fixtures/poly_aniso/",
        ),
    )
