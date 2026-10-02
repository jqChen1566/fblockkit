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
- ``SYMM`` (required whenever a type carries more than one equivalent
  centre): per type, a count line and that many rotation matrices (three
  rows of three numbers each), mapping the representative centre onto each
  site of the type (the first is the identity);
- ``PAIR`` (alias ``LIN1``) or ``LIN3``: the coupled site pairs.  ``PAIR``
  carries one ``J`` per pair under the Lines-type Hamiltonian
  ``sum -J s_i . s_j``; ``LIN3`` carries three, one per Cartesian axis
  (``sum -J_alpha s_i,alpha s_j,alpha``).  The full anisotropic ``LIN9``
  form (nine values per pair) is *refused with the measured engine defect*:
  the otool_poly_aniso of ORCA 6.1.1 aborts in its own printout for every
  probed LIN9 input (Fortran format/type mismatch at
  ``otool_aniso/poly/input_process.f90`` line 286, format ``(15x,3F9.5)``,
  exit status 2 -- the manual's own example values included; measured
  2026-10-01, fixture ``lin9_probe.out``);
- ``COOR`` (optional): one symmetrised Cartesian coordinate per centre type
  (Angstrom); activating it adds the exact dipole-dipole coupling for the
  declared pairs;
- ``TINT`` (optional): the susceptibility sampling grid ``t_min t_max
  n_points``.

What the module refuses: a type count outside 1..6, a centre count or
spin-orbit count below its structural floor, a pair index outside the
cluster, a duplicated pair, a pair line of the wrong arity, a coordinate
array that does not match the type count, a temperature grid that is not
increasing, symmetry matrices that do not match the per-type centre counts,
**a type with several equivalent centres but no symmetry matrices** (the
driver's own measured check is "SYMM is mandatory for cases when
neq(:) > 1!", but the driver prints a serious-error banner and *still exits
0*, so the guard must be client-side), and the LIN9 form (engine defect
above).  The *physics* of the ``J`` values, the coordinates and the
rotation matrices is the caller's (the writer does not fit or judge them;
the source's model is the Lines exchange plus exact dipole-dipole
coupling).
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
#: driver's format scouting; the manual's examples use one to three).
MAX_TYPES = 6

#: The pair models the writer renders (the manual's Lines family):
#: "lines" = PAIR/LIN1 (one J per pair), "lin3" = one J per Cartesian axis.
#: "lin9" is accepted by :func:`build_input` as a value only to refuse it
#: with the measured engine defect (see the module docstring).
PAIR_MODELS = ("lines", "lin3")


class PolyAnisoPlanError(ValueError):
    """The input plan cannot be written (with a next step)."""


@dataclass(frozen=True)
class PolyAnisoPlan:
    """The validated plan and the generated input text.

    ``pairs`` rows are ``(i, j, J)`` for the ``lines`` model and
    ``(i, j, Jx, Jy, Jz)`` for ``lin3``; ``symmetry`` is one tuple of
    matrices per centre type (each matrix nine numbers, row-major), or
    ``None`` when the cluster needs no SYMM block.
    """

    equivalent_centres: tuple[int, ...]
    spin_orbit_states: tuple[int, ...]
    pairs: tuple[tuple, ...]
    coordinates: tuple[tuple[float, float, float], ...] | None
    temperature_grid: tuple[float, float, int] | None
    text: str
    symmetry: tuple[tuple[tuple[float, ...], ...], ...] | None = None
    pair_model: str = "lines"

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
    symmetry=None,
    pair_model="lines",
) -> PolyAnisoPlan:
    """Validate the cluster description and render the ``&POLY_ANISO`` input.

    ``equivalent_centres`` / ``spin_orbit_states``: one number per centre type
    (ints); ``pairs``: ``(i, j, J)`` (``pair_model="lines"``) or
    ``(i, j, Jx, Jy, Jz)`` (``"lin3"``) with 1-based site indices over all
    centres and the J values in cm-1; ``coordinates``: ``(x, y, z)`` per type
    in Angstrom (``None`` omits the COOR block); ``temperature_grid``:
    ``(t_min, t_max, n_points)`` (``None`` omits the TINT block);
    ``symmetry``: one tuple of rotation matrices per centre type (each matrix
    nine numbers, row-major; the first maps the representative onto itself,
    i.e. the identity) -- required whenever a type carries more than one
    equivalent centre, omitted otherwise.
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
    if pair_model == "lin9":
        raise PolyAnisoPlanError(
            "the LIN9 (full anisotropic) pair model cannot be written: the "
            "otool_poly_aniso of ORCA 6.1.1 aborts in its own printout for every "
            "probed LIN9 input (Fortran format/type mismatch at "
            "otool_aniso/poly/input_process.f90 line 286, format (15x,3F9.5), "
            "exit status 2; the manual's own example values included -- measured "
            "2026-10-01). Next step: use the axis-diagonal LIN3 form or the "
            "isotropic Lines (PAIR) form."
        )
    if pair_model not in PAIR_MODELS:
        raise PolyAnisoPlanError(
            f"unknown pair model {pair_model!r}; the writer renders "
            + " and ".join(repr(model) for model in PAIR_MODELS)
            + ". Next step: pick one of those."
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
        if pair_model == "lines":
            if len(pair) != 3:
                raise PolyAnisoPlanError(
                    f"a pair line under the Lines model reads 'i j J' (three "
                    f"numbers); this one carries {len(pair)}. Next step: check the "
                    "pair list."
                )
            first, second, coupling = pair
            coupling_checked: tuple[float, ...] = (float(coupling),)
        else:  # lin3
            if len(pair) != 5:
                raise PolyAnisoPlanError(
                    f"a pair line under the LIN3 model reads 'i j Jx Jy Jz' (five "
                    f"numbers); this one carries {len(pair)}. Next step: check the "
                    "pair list."
                )
            first, second = pair[0], pair[1]
            coupling_checked = tuple(float(value) for value in pair[2:5])
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
        checked_pairs.append((first, second, *coupling_checked))
    symmetry_checked = None
    if symmetry is not None:
        if len(symmetry) != len(centres):
            raise PolyAnisoPlanError(
                f"the SYMM block carries matrix list(s) for {len(symmetry)} type(s) "
                f"while the cluster has {len(centres)}. Next step: give one matrix "
                "list per centre type (count-1 types get the identity)."
            )
        rows = []
        for type_index, (matrices, count) in enumerate(zip(symmetry, centres), start=1):
            matrices_checked = tuple(
                tuple(float(value) for value in matrix) for matrix in matrices
            )
            if len(matrices_checked) != count:
                raise PolyAnisoPlanError(
                    f"type {type_index} carries {count} equivalent centre(s) but "
                    f"{len(matrices_checked)} rotation matrix(es). Next step: give "
                    "one matrix per equivalent centre of the type (the first maps "
                    "the representative onto itself, i.e. the identity)."
                )
            for matrix in matrices_checked:
                if len(matrix) != 9:
                    raise PolyAnisoPlanError(
                        f"a rotation matrix of type {type_index} carries "
                        f"{len(matrix)} number(s); a matrix needs nine (three rows "
                        "of three, row-major). Next step: check the matrix list."
                    )
            rows.append(matrices_checked)
        symmetry_checked = tuple(rows)
    else:
        for type_index, count in enumerate(centres, start=1):
            if count > 1:
                raise PolyAnisoPlanError(
                    f"type {type_index} carries {count} equivalent centres but no "
                    "symmetry matrices were given. The driver's own check "
                    '(measured): "SYMM is mandatory for cases when: neq(:) > 1!" '
                    "-- it prints a serious-error banner but still exits 0, so this "
                    "guard is client-side. Next step: give one rotation matrix (nine "
                    "numbers, row-major) per equivalent centre of the type, mapping "
                    "the representative onto each of them (the first is the "
                    "identity)."
                )
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
        # a signed zero prints as 0 (the same rule the mkl writer carries:
        # "-0" is a distinct text token and only invites confusion)
        return f"{value if value else 0.0:g}"

    lines = ["&POLY_ANISO", "", "NNEQ"]
    lines.append(f"  {len(centres)}  T")
    lines.append("  " + "  ".join(str(count) for count in centres))
    lines.append("  " + "  ".join(str(count) for count in states))
    if symmetry_checked is not None:
        lines.append("")
        lines.append("SYMM")
        for matrices in symmetry_checked:
            lines.append(f"  {len(matrices)}")
            for matrix in matrices:
                for row in range(3):
                    cells = matrix[3 * row : 3 * row + 3]
                    lines.append("  " + "  ".join(number(value) for value in cells))
    lines.append("")
    lines.append("PAIR" if pair_model == "lines" else "LIN3")
    lines.append(f"  {len(checked_pairs)}")
    if pair_model == "lines":
        for first, second, coupling in checked_pairs:
            lines.append(f"  {first} {second} {number(coupling)}")
    else:
        for first, second, jx, jy, jz in checked_pairs:
            lines.append(
                f"  {first} {second} {number(jx)} {number(jy)} {number(jz)}"
            )
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
        symmetry=symmetry_checked,
        pair_model=pair_model,
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
        + (
            "(Lines isotropic, J in cm-1: "
            + ", ".join(f"{i}-{j} {coupling:g}" for i, j, coupling in plan.pairs)
            + ")"
            if plan.pair_model == "lines"
            else "(LIN3 axis-diagonal, Jx/Jy/Jz in cm-1: "
            + ", ".join(
                f"{i}-{j} [{jx:g} {jy:g} {jz:g}]" for i, j, jx, jy, jz in plan.pairs
            )
            + ")"
        ),
        "  symmetry (SYMM): "
        + (
            "not included (unique centres only)"
            if plan.symmetry is None
            else f"included ({sum(len(matrices) for matrices in plan.symmetry)} "
            f"matrices across {len(plan.symmetry)} type(s))"
        ),
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
        "  - the J values, the coordinates and the SYMM rotation matrices are "
        "the caller's (the writer does not fit or judge them); the exchange "
        "model is the Lines-type Hamiltonian over the declared pairs "
        "(isotropic 'sum -J s_i.s_j', or axis-diagonal under LIN3)",
        "  - the COOR block computes the dipolar coupling only for the declared "
        "pairs; coordinates are the symmetrised per-type positions in Angstrom",
        "  - the SYMM matrices are required whenever a type carries more than "
        "one equivalent centre: the driver's own check (measured) is 'SYMM is "
        "mandatory for cases when: neq(:) > 1!', but the driver still exits 0 "
        "on that error, so the writer refuses the omission client-side",
        "  - the full anisotropic LIN9 form is a documented termination: the "
        "otool_poly_aniso of ORCA 6.1.1 aborts in its own printout for every "
        "probed LIN9 input (Fortran format/type mismatch at input_process.f90 "
        "line 286, exit status 2; measured 2026-10-01)",
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
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "The SYMM and LIN3 variants are accepted end to end (measured "
                "2026-10-01, otool_poly_aniso v1.0.0, ORCA 6.1.1): a one-type "
                "two-centre cluster with the identity+inversion matrices and a "
                "LIN3 pair both return rc = 0 and 'POLY_ANISO finished "
                "sucessfully!' (fixtures poly_aniso/symm_probe.out and "
                "lin3_probe.out). The omitted-SYMM case is a measured silent "
                "failure of the driver: it prints 'SYMM is mandatory for cases "
                "when: neq(:) > 1!' and 'THE CODE WILL STOP NOW', yet exits 0 "
                "-- the writer refuses the omission client-side. The LIN9 form "
                "aborts for every probed input shape (Fortran format/type "
                "mismatch in the tool's own printout, exit status 2; fixture "
                "lin9_probe.out)."
            ),
            ref="tests/test_poly_aniso_recipe.py; fixtures/poly_aniso/",
        ),
    )
