"""A11: cross-structure orbital mapping -- consistent active spaces along a path.

Active spaces must be *consistent* between the structures of a reaction path:
the total correlation energy is computed inconsistently otherwise, and the
relative energies become erratic.  The source (Bensberg & Reiher 2023) maps the
localized orbitals of all structures onto each other -- no interpolation of
structures is needed -- and then transfers the active space selected for one (or
a few) structure(s) to all the others through that map:

1. localize the occupied and the virtual valence orbitals of every structure
   (intrinsic bond orbitals, for transferability);
2. map every orbital onto every other structure's orbitals through the
   similarity criteria of Eq. (1): same kinetic energy, same shell-wise
   populations,

       |t_iL - t_jK| / Eh < tau_kin   and   sum_a |q_iL^a - q_jK^a| < tau_loc

   with ``t`` the orbital kinetic energy and ``q^a`` the orbital-wise
   population of shell ``a``;
3. collect the orbitals that cannot be matched in every structure into a
   non-matchable set -- exactly the orbitals that change along the path (the
   bonds being broken or formed);
4. build the set maps: the self-map ``S_iL`` (orbitals of one structure that
   match each other, e.g. a degenerate pair), the cross map ``m_LK`` between
   self-map sets (defined only when every orbital of one set maps onto exactly
   one whole other set and not to anything else, in both directions), the
   bijective map ``M_LK`` on the unambiguous sets, and the set map ``A_LK``
   covering the non-matchable orbitals;
5. the consistent active space of a structure is the union of all mapped sets
   that contain at least one orbital selected (e.g. by the entropy protocol) in
   any structure -- i.e. once one member of an equivalence class is selected,
   the whole class joins the active space in every structure.

The single parameter is the threshold ``tau`` (``tau = tau_kin = tau_loc``).
When it is not given, the source's Eq. (8) is used per orbital set: the
threshold that minimizes the number of non-matchable orbitals
(``tau_min``); the whole curve is reported because its plateaus are the
feature, not one number.  The source notes that ``tau <= 0.5`` is the
reasonable range for comparing orbital populations and reads its own example
at ``tau_min = 0.5`` for both sets.

Measured boundaries of this ORCA route
--------------------------------------
The source computes both criteria in Serenity (kinetic energies; shell-wise
*IAO* populations).  From an ORCA route this module takes (measured on the
frozen N2 scan fixtures):

- the kinetic energy from the ``T-Matrix`` block of the ``orca_2json`` export
  (requested with a ``<basename>.json.conf`` carrying ``"1elIntegrals":
  ["H", "S", "T", "V"]``); the values were cross-checked against the printed
  ``KINETIC ENERGY MATRIX (AU)`` of the same run (agreement 5e-7, the print's
  own precision);
- the populations as shell-wise **Loewdin** populations computed here from the
  export's coefficients and overlap.  The source's own text licenses any
  orbital-wise population analysis (it chose IAO for basis-set insensitivity);
  a same-basis comparison across structures -- the only comparison this module
  performs -- is not basis-set-insensitive by nature, but the population scale
  differs from the source's IAO one, so ``tau`` is a data-calibrated parameter
  here, never copied from their absolute values;
- the localization: ORCA's ``orca_loc``.  Measured on 6.1.1: the IAO-based
  methods cannot localize virtual orbitals ("impossible to localize so many
  MOs by IAO based Methods (maybe virtual MOs)"), so the frozen fixtures
  localize the occupied block with IAO-BOYS (``LocMet 4``, the IBO scheme) and
  the virtual block with Foster-Boys (``LocMet 2``); both with the random seed
  switched off for reproducibility.

The module is pure post-processing: exports in, maps out.  It does not run
ORCA and it does not select the active space itself -- the selection (menu
12's entropy protocol or any other) enters as an optional annotation, and the
module reports the *consistent* space it implies.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass

import numpy as np

from ..knowledge.models import EVIDENCE_LITERATURE, EVIDENCE_MEASURED, Evidence, ReportSection

__all__ = [
    "MappingError",
    "ShellRef",
    "OrbitalDescriptor",
    "StructureDescription",
    "MappingResult",
    "descriptors_from_export",
    "match_orbitals",
    "active_overlap_determinant",
    "active_overlap_series",
    "analyze",
    "render",
    "render_active_overlap",
    "run",
    "evidence",
]

#: Threshold grid used when ``tau`` is not given (the source reads its examples
#: at ``tau_min = 0.5``; the whole curve is reported).
DEFAULT_TAU_GRID = tuple(round(0.025 * step, 3) for step in range(2, 21))  # 0.05 .. 0.5

#: The source's own reading of the reasonable range: "while it remains reasonable
#: for comparing orbital populations (in this work tau <= 0.5)".
TAU_CEILING = 0.5

#: Occupation splitting: HF-like occupations make a definite occupied/virtual
#: split; anything fractional is the mapping's undefined case (see analyze()).
OCCUPIED_FLOOR = 1.9
VIRTUAL_CEILING = 0.1

OCCUPIED = "occupied"
VIRTUAL = "virtual"


class MappingError(ValueError):
    """The mapping cannot be built from the given structures (with a next step)."""


@dataclass(frozen=True, order=True)
class ShellRef:
    """One AO shell of the basis: ``(centre, element, shell, angular)``.

    The mapping's population criterion is defined per shell, and the same shell
    must exist -- under the same key -- in every structure.  The element is
    carried for the labels, the centre index for the ordering: both structures
    must use the same atom order and basis set, otherwise the criterion has no
    meaning and the analysis refuses.
    """

    centre: int
    element: str
    shell: int
    angular: str

    def render(self) -> str:
        return f"{self.element}{self.centre} {self.shell}{self.angular}"


@dataclass(frozen=True)
class OrbitalDescriptor:
    """One orbital's mapping descriptors: occupation, kinetic energy, shell populations."""

    index: int
    occupation: float
    kinetic: float
    shells: tuple[tuple[ShellRef, float], ...]

    def dominant(self, count: int = 3) -> tuple[tuple[ShellRef, float], ...]:
        """The ``count`` largest shell populations (by magnitude), for the report."""
        ranked = sorted(self.shells, key=lambda item: -abs(item[1]))
        return tuple(ranked[:count])


@dataclass(frozen=True)
class StructureDescription:
    """One structure's orbital set, reduced to the mapping descriptors."""

    name: str
    orbitals: tuple[OrbitalDescriptor, ...]

    @property
    def shell_keys(self) -> tuple[ShellRef, ...]:
        return tuple(shell for shell, _ in self.orbitals[0].shells)


@dataclass(frozen=True)
class MappingResult:
    """The maps, the non-matchable sets, the tau curve and the consistent space."""

    names: tuple[str, ...]
    tau_occupied: float
    tau_virtual: float
    tau_given: bool
    tau_curve: tuple[tuple[float, int, int], ...]  # (tau, no-match occupied, no-match virtual)
    no_match: dict[str, tuple[tuple[int, ...], ...]]  # kind -> per structure
    classes: dict[str, tuple[tuple[tuple[int, int], ...], ...]]  # kind -> classes of (structure, orbital)
    consistent: tuple[tuple[int, ...], ...] | None  # per structure, when selections were given
    occupied_counts: tuple[tuple[int, ...], ...]  # per structure, the indices in each kind
    virtual_counts: tuple[tuple[int, ...], ...]
    checks: tuple[str, ...]
    notes: tuple[str, ...]


# --- descriptors from an export ----------------------------------------------


def descriptors_from_export(export, name: str) -> StructureDescription:
    """Build the mapping descriptors of one ``orca_2json`` export.

    Needs the ``T-Matrix`` block (kinetic), the ``S-Matrix`` and the AO labels,
    and HF-like occupations; each missing piece is refused with the next step.
    """
    if export.kinetic is None:
        raise MappingError(
            f"the export '{name}' carries no T-Matrix, so the kinetic-energy criterion "
            "cannot be computed. Next step: re-export with a <basename>.json.conf "
            'carrying {"MOCoefficients": true, "1elIntegrals": ["H", "S", "T", "V"]} '
            "and run orca_2json on the .gbw again."
        )
    if export.overlap is None or export.ao_labels is None:
        raise MappingError(
            f"the export '{name}' carries no S-Matrix or no AO labels, so the "
            "shell-wise populations cannot be formed. Next step: re-export with the "
            "S-Matrix and the orbital labels (the orca_2json defaults plus the "
            "1elIntegrals configuration)."
        )
    coefficients = np.array(export.mo_coefficients, dtype=float)  # (n_mo, n_ao)
    overlap = np.array(export.overlap, dtype=float)
    kinetic = np.array(export.kinetic, dtype=float)
    occupations = np.array(export.mo_occupations, dtype=float)
    fractional = [
        index
        for index, occupation in enumerate(occupations)
        if VIRTUAL_CEILING <= occupation <= OCCUPIED_FLOOR
    ]
    if fractional:
        raise MappingError(
            f"the export '{name}' has fractional occupations on orbital(s) "
            f"{fractional} (an active space), so the occupied/virtual split the "
            "mapping needs is not defined. Next step: map the Hartree-Fock orbitals "
            "of the path first (localize them with orca_loc), and select the active "
            "space from the map."
        )
    energies = np.einsum("ma,ab,mb->m", coefficients, kinetic, coefficients)
    populations = coefficients.T * (overlap @ coefficients.T)  # (n_ao, n_mo)
    groups: dict[ShellRef, list[int]] = {}
    for column, label in enumerate(export.ao_labels):
        key = ShellRef(centre=label.center, element=label.element, shell=label.shell, angular=label.angular)
        groups.setdefault(key, []).append(column)
    shell_keys = sorted(groups)
    shell_populations = np.array(
        [[populations[groups[key], mo].sum() for mo in range(export.n_mo)] for key in shell_keys]
    )
    orbitals = tuple(
        OrbitalDescriptor(
            index=index,
            occupation=float(occupations[index]),
            kinetic=float(energies[index]),
            shells=tuple(
                (key, float(shell_populations[row, index])) for row, key in enumerate(shell_keys)
            ),
        )
        for index in range(export.n_mo)
    )
    return StructureDescription(name=name, orbitals=orbitals)


# --- the mapping arithmetic ---------------------------------------------------


def match_orbitals(a: OrbitalDescriptor, b: OrbitalDescriptor, tau: float) -> bool:
    """The source's Eq. (1): both criteria, with ``tau = tau_kin = tau_loc``."""
    if abs(a.kinetic - b.kinetic) >= tau:
        return False
    return sum(abs(x - y) for (_, x), (_, y) in zip(a.shells, b.shells)) < tau


# --- the active-space overlap check (AOP, SI-derived) -------------------------


def active_overlap_determinant(export_a, export_b, active_a, active_b=None) -> float:
    """``|det S_act|`` between two structures' active blocks (the AOP scalar).

    The second source's overlap-preservation check, from its Supporting
    Information (Eq. (2) there; the main text is not yet available): the active
    orbitals of the two structures are matched by their overlap,

        S_act[i, j] = <psi_i^(A) | psi_j^(B)>  ~=  c_i^(A)T  S^(B)  c_j^(B),

    using the *current* (second) structure's AO overlap -- justified in the
    source because the geometries change only minimally between adjacent steps
    of its geodesic interpolation.  ``|det S_act| ~ 1`` means the active space
    was preserved between the two structures; a value towards 0 means at least
    one active orbital exchanged with the inactive space.  The determinant is
    taken in absolute value: orbital phases are arbitrary.

    Applicability: the source uses this between *adjacent* small steps.  For
    distant structures the one-metric approximation behind ``S_act`` no longer
    holds and the number is a demonstration of degradation, not a calibrated
    reading; the report says so and prints the value without a threshold.
    """
    coefficients_a = np.array(export_a.mo_coefficients, dtype=float).T
    coefficients_b = np.array(export_b.mo_coefficients, dtype=float).T
    if export_b.overlap is None:
        raise MappingError(
            "the second export carries no S-Matrix, so the active overlap cannot be "
            "formed. Next step: re-export with the S-Matrix."
        )
    overlap = np.array(export_b.overlap, dtype=float)
    indices_a = [int(index) for index in active_a]
    indices_b = indices_a if active_b is None else [int(index) for index in active_b]
    for export, coefficients, indices in (
        (export_a, coefficients_a, indices_a),
        (export_b, coefficients_b, indices_b),
    ):
        for index in indices:
            if not 0 <= index < coefficients.shape[1]:
                raise MappingError(
                    f"the active index {index} is outside the export "
                    f"'{export.base_name}'s {coefficients.shape[1]} orbitals. Next step: "
                    "check the active-space list."
                )
    block = coefficients_a[:, indices_a].T @ overlap @ coefficients_b[:, indices_b]
    return float(abs(np.linalg.det(block)))


def active_overlap_series(exports, actives, names=None) -> tuple[tuple[str, str, float], ...]:
    """The AOP scalar for every adjacent pair of a series (given order).

    ``actives``: mapping structure name -> active-orbital indices; ``names``:
    the display/lookup names (default: the exports' own base names).  A pair
    whose two lists disagree in length is refused (the determinant needs
    square blocks from the same-size spaces).
    """
    labels = [export.base_name for export in exports] if names is None else [str(n) for n in names]
    if len(labels) != len(exports):
        raise MappingError(
            "the name list and the export list differ in length. Next step: check the "
            "manifest."
        )
    for name in labels:
        if name not in actives:
            raise MappingError(
                f"the structure '{name}' has no active-space list. Next step: give the "
                "'active' entry for every structure of the series."
            )
    rows = []
    for position, (first, second) in enumerate(zip(exports, exports[1:])):
        name_a, name_b = labels[position], labels[position + 1]
        active_a = [int(index) for index in actives[name_a]]
        active_b = [int(index) for index in actives[name_b]]
        if len(active_a) != len(active_b):
            raise MappingError(
                f"'{name_a}' has {len(active_a)} active orbitals while '{name_b}' has "
                f"{len(active_b)}: the overlap determinant needs two equally sized "
                "spaces. Next step: check the active-space lists (the mapping's "
                "consistent space gives equal sizes by construction)."
            )
        rows.append(
            (
                name_a,
                name_b,
                active_overlap_determinant(first, second, active_a, active_b),
            )
        )
    return tuple(rows)


def _kind_sets(structure: StructureDescription) -> dict[str, tuple[int, ...]]:
    occupied = tuple(
        orbital.index for orbital in structure.orbitals if orbital.occupation > OCCUPIED_FLOOR
    )
    virtual = tuple(
        orbital.index for orbital in structure.orbitals if orbital.occupation < VIRTUAL_CEILING
    )
    return {OCCUPIED: occupied, VIRTUAL: virtual}


def _maps_for_sets(targets, tau):
    """Build every map of the source's Eqs. (4)-(7) for one orbital set.

    ``targets``: one list of ``OrbitalDescriptor`` per structure (same length).
    Returns ``(classes, no_match_per_structure)`` where a class is a tuple of
    ``(structure, position)`` pairs -- positions within ``targets`` -- and the
    non-matchable orbitals of every structure form one class together, because
    the source's rule includes them as a block ("all of those valence orbitals
    that are varying along a reaction path are assigned to the active orbital
    space as soon as one of them is selected for it").

    The unit of the maps is the self-map *set* (the orbitals of one structure
    that match each other, e.g. a degenerate pair), never the individual
    orbital: a class is accepted only when, for every pair of structures,
    every member of the class in one structure maps onto exactly the class's
    whole membership in the other and onto nothing else -- the source's
    condition "for some structure combination", checked here for all of them.
    Classes are grown from structure 0's self-map sets; a set that fails for
    any combination stays unassigned and lands in the non-matchable block.
    """
    n_struct = len(targets)
    n_orb = len(targets[0])
    matches = lambda l, i, k, j: match_orbitals(targets[l][i], targets[k][j], tau)  # noqa: E731
    cross: dict[tuple[int, int], list[frozenset[int]]] = {
        (l, k): [frozenset(j for j in range(n_orb) if matches(l, i, k, j)) for i in range(n_orb)]
        for l in range(n_struct)
        for k in range(n_struct)
        if k != l
    }
    self_sets = [
        [
            frozenset(j for j in range(n_orb) if matches(l, i, l, j))
            for i in range(n_orb)
        ]
        for l in range(n_struct)
    ]
    taken = [set() for _ in range(n_struct)]
    classes: list[tuple[tuple[int, int], ...]] = []
    seen: set[frozenset[int]] = set()
    for anchor in range(n_orb):
        own = self_sets[0][anchor]
        if own in seen:
            continue
        seen.add(own)
        candidate = [(0, member) for member in own]
        ok = True
        for k in range(1, n_struct):
            target = cross[(0, k)][anchor]
            if not target:
                ok = False
                break
            candidate += [(k, member) for member in target]
        if ok:
            candidate = sorted(set(candidate))
            membership = {
                l: frozenset(orbital for s, orbital in candidate if s == l)
                for l in range(n_struct)
            }
            for l in range(n_struct):
                for k in range(n_struct):
                    if l == k or not ok:
                        continue
                    for orbital in membership[l]:
                        if cross[(l, k)][orbital] != membership[k]:
                            ok = False
                            break
        if ok and not any(orbital in taken[s] for s, orbital in candidate):
            for structure, orbital in candidate:
                taken[structure].add(orbital)
            classes.append(tuple(candidate))
    unmatched = [
        [i for i in range(n_orb) if i not in taken[structure]]
        for structure in range(n_struct)
    ]
    if any(unmatched):
        classes.append(
            tuple((structure, i) for structure in range(n_struct) for i in unmatched[structure])
        )
    return tuple(classes), tuple(tuple(items) for items in unmatched)


def analyze(structures, *, tau: float | None = None, selections=None) -> MappingResult:
    """Map the structures' orbitals and (optionally) make a selection consistent.

    ``structures``: :class:`StructureDescription` objects (same atom order and
    basis set).  ``tau``: an explicit threshold for both criteria, or ``None``
    to minimise the non-matchable set per orbital set (the source's Eq. (8)).
    ``selections``: optional mapping ``structure name -> iterable of orbital
    indices`` (0-based, the whole orbital list); when given, the consistent
    active space is the union of every mapped class containing a selected
    orbital, per structure.
    """
    if len(structures) < 2:
        raise MappingError(
            "the mapping needs at least two structures. Next step: give the exports "
            "of two (or more) structures along the path."
        )
    reference = structures[0]
    for structure in structures[1:]:
        if structure.shell_keys != reference.shell_keys:
            raise MappingError(
                f"the structure '{structure.name}' has a different AO shell set than "
                f"'{reference.name}', so the shell-wise populations are not comparable. "
                "Next step: map structures with the same basis set and atom order."
            )
        if len(structure.orbitals) != len(reference.orbitals):
            raise MappingError(
                f"the structure '{structure.name}' has {len(structure.orbitals)} orbitals "
                f"while '{reference.name}' has {len(reference.orbitals)}. Next step: map "
                "structures of the same system (same basis, same electron count)."
            )
    kinds = {name: [_kind_sets(structure)[name] for structure in structures] for name in (OCCUPIED, VIRTUAL)}
    for name in (OCCUPIED, VIRTUAL):
        counts = {len(indices) for indices in kinds[name]}
        if len(counts) != 1:
            raise MappingError(
                f"the structures carry different numbers of {name} orbitals "
                f"({sorted(counts)}). Next step: check the occupations -- the mapping is "
                "defined for one and the same orbital count per set."
            )
    # the tau curve, or the single given tau
    if tau is None:
        grid = DEFAULT_TAU_GRID
        curve = []
        for value in grid:
            row = []
            for name in (OCCUPIED, VIRTUAL):
                targets = [
                    [structure.orbitals[index] for index in indices]
                    for structure, indices in zip(structures, kinds[name])
                ]
                _, unmatched = _maps_for_sets(targets, value)
                row.append(sum(len(items) for items in unmatched))
            curve.append((value, row[0], row[1]))
        # Eq. (8): the threshold minimising the non-matchable count; the largest
        # such threshold is taken (the flat plateau, capped at the source's 0.5)
        def pick(column: int) -> float:
            best = min(item[column] for item in curve)
            candidates = [item[0] for item in curve if item[column] == best and item[0] <= TAU_CEILING]
            return max(candidates) if candidates else grid[-1]

        tau_occupied = pick(1)
        tau_virtual = pick(2)
        tau_map = {OCCUPIED: tau_occupied, VIRTUAL: tau_virtual}
    else:
        tau_map = {OCCUPIED: float(tau), VIRTUAL: float(tau)}
        curve = []
        row = []
        for name in (OCCUPIED, VIRTUAL):
            targets = [
                [structure.orbitals[index] for index in indices]
                for structure, indices in zip(structures, kinds[name])
            ]
            _, unmatched = _maps_for_sets(targets, tau_map[name])
            row.append(sum(len(items) for items in unmatched))
        curve.append((float(tau), row[0], row[1]))
    # the final maps
    no_match: dict[str, tuple[tuple[int, ...], ...]] = {}
    classes: dict[str, tuple[tuple[tuple[int, int], ...], ...]] = {}
    for name in (OCCUPIED, VIRTUAL):
        targets = [
            [structure.orbitals[index] for index in indices]
            for structure, indices in zip(structures, kinds[name])
        ]
        class_list, unmatched = _maps_for_sets(targets, tau_map[name])
        # translate positions back to the global orbital indices
        translated = []
        for klass in class_list:
            translated.append(
                tuple((structure, kinds[name][structure][position]) for structure, position in klass)
            )
        classes[name] = tuple(translated)
        no_match[name] = tuple(
            tuple(kinds[name][structure][position] for position in unmatched[structure])
            for structure in range(len(structures))
        )
    # the consistent active space, when a selection was given
    consistent = None
    if selections is not None:
        chosen = []
        for structure in range(len(structures)):
            chosen.append(set(selections.get(structures[structure].name, ())))
        consistent_sets = [set() for _ in structures]
        for name in (OCCUPIED, VIRTUAL):
            for klass in classes[name]:
                if any(orbital in chosen[structure] for structure, orbital in klass):
                    for structure, orbital in klass:
                        consistent_sets[structure].add(orbital)
        consistent = tuple(tuple(sorted(items)) for items in consistent_sets)
    checks = (
        "the criteria are the source's Eq. (1): |t_iL - t_jK| < tau (kinetic energy, in "
        "Eh) and the summed absolute shell-wise population difference < tau, with one tau "
        "for both (tau = tau_kin = tau_loc)",
        "the consistent classes are the source's maps M_LK (unambiguous sets, bijective) "
        "and A_LK (the non-matchable block); a class is built only when every orbital of "
        "a self-map set maps onto exactly one whole self-map set of the other structure "
        "and not to anything else, in both directions",
        f"the non-matchable orbitals of all structures form one class together -- the "
        "source's rule: a path-varying orbital joins the active space as soon as one of "
        "its class is selected anywhere",
        f"when tau is not given, the source's Eq. (8) is applied per orbital set over "
        f"the grid {DEFAULT_TAU_GRID[0]:g}-{DEFAULT_TAU_GRID[-1]:g}: the threshold "
        "minimising the non-matchable count; among ties the largest is taken (the flat "
        f"plateau the source itself reads, capped at {TAU_CEILING:g})",
    )
    notes = (
        "The populations here are shell-wise Loewdin populations computed from the "
        "export's coefficients and overlap; the source uses shell-wise IAO populations "
        "("
        "chosen for basis-set insensitivity, and its own text licenses any orbital-wise "
        "population analysis). The scale of the population criterion therefore differs "
        "from the source's, and tau is a data-calibrated parameter here -- the report "
        "shows the whole tau curve instead of copying an absolute threshold.",
        "The orbitals must be localized before mapping (the source: the intrinsic bond "
        "orbital scheme provides the transferable sets, and the maps are built "
        "separately for the occupied and the virtual set). Measured ORCA boundary: the "
        "IAO-based localization refuses the virtual block ('impossible to localize so "
        "many MOs by IAO based Methods'), so the frozen example localizes the occupied "
        "block with IAO-BOYS and the virtual block with Foster-Boys.",
        "The mapping compares the structures in their own frames: same basis set, same "
        "atom order, same orientation of the whole system. Rotated or permuted "
        "structures are refused by the shell-key check or give a meaningless map.",
        "The selection itself is not part of this module: menu 12's entropy protocol (or "
        "any other) provides the annotation, and the module reports the consistent "
        "active space it implies through the classes.",
    )
    return MappingResult(
        names=tuple(structure.name for structure in structures),
        tau_occupied=tau_map[OCCUPIED],
        tau_virtual=tau_map[VIRTUAL],
        tau_given=tau is not None,
        tau_curve=tuple(curve),
        no_match=no_match,
        classes=classes,
        consistent=consistent,
        occupied_counts=tuple(kinds[OCCUPIED]),
        virtual_counts=tuple(kinds[VIRTUAL]),
        checks=checks,
        notes=notes,
    )


# --- output -------------------------------------------------------------------


def render(result: MappingResult, structures) -> str:
    """The descriptor tables, the tau curve, the classes and the consistent space."""
    lines = [
        "Cross-structure orbital mapping (kinetic energy + shell-wise populations, "
        "per the source's Eq. (1)):",
        f"  structures: {', '.join(result.names)}",
        f"  tau: occupied {result.tau_occupied:g}, virtual {result.tau_virtual:g}"
        + ("" if result.tau_given else "  (from the Eq. (8) minimum over the grid)"),
        "",
        "Tau curve (non-matchable orbitals per set):",
        "    tau    occupied   virtual",
    ]
    for value, occ, virt in result.tau_curve:
        mark = ""
        if not result.tau_given and value == result.tau_occupied:
            mark = "  <- occupied tau_min"
        lines.append(f"  {value:>5.3f}   {occ:>7}   {virt:>6}{mark}")
    for position, structure in enumerate(structures):
        lines += ["", f"Structure '{structure.name}':"]
        lines.append(
            "   orbital   occ      t(Eh)   shells (largest Loewdin shell populations)"
        )
        for orbital in structure.orbitals:
            shells = "  ".join(
                f"{shell.render()}:{value:+.3f}" for shell, value in orbital.dominant()
            )
            lines.append(
                f"   {orbital.index:>7}  {orbital.occupation:>5.2f}  {orbital.kinetic:>9.4f}   {shells}"
            )
    for name in (OCCUPIED, VIRTUAL):
        lines += ["", f"{name.capitalize()} classes:"]
        for klass in result.classes[name]:
            members = ", ".join(
                f"{structure}:{orbital}" for structure, orbital in klass
            )
            lines.append(f"  - class ({len(klass)} orbital-structure pairs): {members}")
        unmatched = result.no_match[name]
        lines.append(
            "  non-matchable: "
            + "; ".join(
                f"{result.names[structure]} {list(unmatched[structure])}"
                for structure in range(len(result.names))
            )
        )
    if result.consistent is not None:
        lines += ["", "Consistent active space implied by the selection:"]
        for structure, indices in enumerate(result.consistent):
            lines.append(f"  {result.names[structure]}: {list(indices)}")
    lines += ["", "Criteria and boundaries:"]
    for item in result.checks + result.notes:
        lines.append(f"  - {item}")
    return "\n".join(lines)


def render_active_overlap(rows) -> str:
    """The active-space overlap (AOP) block: one |det S_act| per adjacent pair."""
    lines = [
        "Active-space overlap between adjacent structures (|det S_act|, the "
        "overlap-preservation scalar of the second source's Supporting Information):",
        "  pair                                        |det S_act|   reading",
    ]
    for first, second, value in rows:
        if value >= 0.9:
            reading = "preserved"
        elif value <= 0.1:
            reading = "an active/inactive exchange is likely"
        else:
            reading = "partial"
        lines.append(f"  {first} -> {second:<28} {value:>10.6f}   {reading}")
    lines += [
        "",
        "Boundaries of this check:",
        "  - S_act uses the *second* structure's AO overlap, the source's own "
        "approximation for small steps of its geodesic interpolation; between distant "
        "structures the number demonstrates degradation rather than certifying "
        "preservation",
        "  - the source gives the qualitative reading only (near 1 preserved, near 0 an "
        "exchange); the printed bands are this project's reading of that scale, not the "
        "source's thresholds (measured on the frozen N2 scan: the two 1s core orbitals "
        "give 1.0000 for every pair, the bond triad 0.9955 across a 0.01-Angstrom step "
        "and 0.72 across 0.5-Angstrom steps)",
        "  - the determinant is taken in absolute value (orbital phases are arbitrary); "
        "it measures whether the active *subspace* survived, while the mapping above "
        "follows individual orbitals -- the two answers can differ (a rotation inside "
        "the window leaves |det| near 1 while the mapping reports no match)",
    ]
    return "\n".join(lines)


def run(
    structures,
    *,
    tau: float | None = None,
    selections=None,
    active_overlap=None,
) -> ReportSection:
    """The analyser entry point for menu 17 (``active_overlap``: series rows)."""
    result = analyze(structures, tau=tau, selections=selections)
    body = render(result, structures)
    if active_overlap:
        body += "\n\n" + render_active_overlap(active_overlap)
    return ReportSection(
        title="A11 cross-structure orbital mapping (consistent active spaces)",
        body=body,
    )


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the mapping protocol and of its ORCA route."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "The orbital-mapping protocol for consistent active spaces along a "
                "reaction path: localize the occupied and virtual valence orbitals, map "
                "them between structures through Eq. (1) (|delta t| < tau and the summed "
                "absolute shell-wise population difference < tau), collect the "
                "non-matchable orbitals (the bonds changing along the path), build the "
                "set maps M_LK / A_LK, and transfer the selected active space through "
                "the maps (Eqs. (4)-(7); the tau_min rule of Eq. (8)). Demonstrated on "
                "the homolytic C-C dissociation and the CH2 rotation of 1-pentene."
            ),
            ref=(
                "Bensberg M., Reiher M., J. Phys. Chem. Lett., 2023, 14, 2112-2118, "
                "DOI 10.1021/acs.jpclett.2c03905 (eqs. (1)-(8); the tau <= 0.5 range)"
            ),
            url="https://doi.org/10.1021/acs.jpclett.2c03905",
            bibkey="bensberg2023corresponding",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "The ORCA route and its boundary, measured on the frozen N2 scan "
                "(fixtures/orca/n2_scan_1.094/1.600/2.600.*): the kinetic matrix comes "
                "from the export's T-Matrix block (cross-checked against the printed "
                "KINETIC ENERGY MATRIX of the same run, agreement 5e-7) and the "
                "populations are computed here as shell-wise Loewdin populations; the "
                "IAO-based localization of ORCA refuses virtual orbitals, so the frozen "
                "chain localizes them with Foster-Boys. On that scan the cores map "
                "uniquely, the N-N bond orbitals (the changing set) come out "
                "non-matchable, and the tau curve shows the source's plateau behaviour."
            ),
            ref="tests/test_orbital_mapping.py; fixtures/orca/n2_scan_*",
        ),
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "The active-space overlap scalar: between adjacent (small-step) "
                "structures, match the active orbitals by their overlap -- approximated "
                "with the current geometry's AO overlap -- and read |det S_act|: near 1 "
                "the active space was preserved, near 0 at least one active/inactive "
                "exchange happened. From the second source's Supporting Information "
                "(its Eq. (2)); the main text is not yet available, so the protocol's "
                "applicability beyond adjacent small steps is recorded as unverified and "
                "the report prints the value without a threshold."
            ),
            ref=(
                "Supporting Information, section I ('Determination of active space "
                "consistency between interpolated geometries', Eq. (2)); main text "
                "pending (DOI 10.1063/5.0058673)"
            ),
            url="https://doi.org/10.1063/5.0058673",
            bibkey="paz2021active",
        ),
    )
