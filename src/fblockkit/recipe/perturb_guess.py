"""randomized orbital perturbation -- the multi-start cure for SCF
multi-solutions (menu 29).

Source: Vaucher & Reiher, "Steering Orbital Optimization out of Local Minima
and Saddle Points Toward Lower Energy", J. Chem. Theory Comput. 2017, 13,
1219-1228 (arXiv:1701.00128; section IV.C; page numbers Crossref-checked):
after a converged reference calculation, perturb the converged orbitals by
randomly selecting pairs of occupied-unoccupied molecular orbitals out of the
15 highest occupied and the 15 lowest unoccupied ones, and mixing each pair
by a random angle in [0, 90) degrees,

    phi_o,new = cos(alpha) phi_o,old + sin(alpha) phi_v,old    (their Eq. 2)
    phi_v,new = cos(alpha) phi_v,old - sin(alpha) phi_o,old    (their Eq. 3)

re-running from the perturbed guess either returns to the same solution or
finds another stationary point; a found lower energy identifies the reference
as wrongly converged.  The source's own boundary: the test cannot guarantee
detection (it is an inexpensive verification, complementary to a stability
analysis, which detects unstable solutions but cannot distinguish local from
global minima).  Their single-point variant ("it can then automatically
inspect whether a solution with lower energy is found") is what this menu
generates.

The division of labour: this module perturbs the converged reference's
orbitals and writes gbw-ready Molekel mkl files (the measured write-back
route, ``parsers/mkl.py``) plus one ORCA input per start; ORCA then decides
where the perturbed guess converges.  The perturbation is fully recorded (the
seed, every pair and every angle) so a batch replays byte-identically.
"""

from __future__ import annotations

import random
import re
from dataclasses import dataclass

import numpy as np

from ..knowledge.elements import ELEMENT_Z, ElementError, element_z
from ..knowledge.models import EVIDENCE_LITERATURE, EVIDENCE_MEASURED, Evidence
from ..parsers.mkl import MklError, MklFile

__all__ = [
    "PerturbError",
    "PairMix",
    "PerturbationPlan",
    "plan_perturbation",
    "mix_orbitals",
    "perturb_mkl",
    "mo_read_variant",
    "check_geometry_match",
    "render",
    "evidence",
]

#: Occupations above this count as occupied, below as virtual (the mkl stores
#: 1.0/0.0 for the converged aufbau solutions this menu targets).
OCCUPIED_THRESHOLD = 0.5


class PerturbError(ValueError):
    """The perturbation cannot be planned or written (with a next step)."""


@dataclass(frozen=True)
class PairMix:
    """One occupied-virtual pair and the random angle it is mixed with."""

    occupied: int  # 0-based orbital index (mkl column)
    virtual: int
    angle_deg: float


@dataclass(frozen=True)
class PerturbationPlan:
    """One start's full perturbation record (reproducible from the seed)."""

    pairs: tuple[PairMix, ...]
    window_occupied: tuple[int, ...]
    window_virtual: tuple[int, ...]


def plan_perturbation(
    occupations,
    energies,
    *,
    n_pairs: int = 10,
    window: int = 15,
    rng: random.Random,
) -> PerturbationPlan:
    """Draw ``n_pairs`` occupied-virtual pairs (with replacement) from the window.

    The windows are the ``window`` occupied orbitals of highest energy and the
    ``window`` virtual orbitals of lowest energy (the source: "the 15 highest
    occupied and the 15 lowest unoccupied molecular orbitals"); the draws are
    independent, so a pair may repeat -- the source's own description ("ten
    random pairs ... mixed with a random angle") is a draw, not a matching.
    """
    if n_pairs < 1:
        raise PerturbError("the number of pairs must be at least one.")
    if window < 1:
        raise PerturbError("the window must be at least one orbital.")
    occ_values = [float(value) for value in occupations]
    energy_values = [float(value) for value in energies]
    if len(occ_values) != len(energy_values):
        raise PerturbError(
            "the occupation and energy lists disagree about the orbital count."
        )
    occupied = [i for i, value in enumerate(occ_values) if value > OCCUPIED_THRESHOLD]
    virtual = [i for i, value in enumerate(occ_values) if value <= OCCUPIED_THRESHOLD]
    if not occupied or not virtual:
        raise PerturbError(
            "the reference has no occupied-virtual separation to perturb ("
            f"{len(occupied)} occupied, {len(virtual)} virtual orbital(s))."
        )
    window_occupied = tuple(
        sorted(occupied, key=lambda i: energy_values[i])[-window:]
    )
    window_virtual = tuple(sorted(virtual, key=lambda i: energy_values[i])[:window])
    pairs = tuple(
        PairMix(
            occupied=rng.choice(window_occupied),
            virtual=rng.choice(window_virtual),
            angle_deg=rng.uniform(0.0, 90.0),
        )
        for _ in range(n_pairs)
    )
    return PerturbationPlan(
        pairs=pairs,
        window_occupied=window_occupied,
        window_virtual=window_virtual,
    )


def mix_orbitals(coefficients, plan: PerturbationPlan):
    """Apply the plan's mixings to a coefficient matrix (Eqs. (2)-(3), column-wise).

    Each mixing is an orthogonal rotation of two columns, so the result stays
    an orthonormal orbital set; the source's angle convention is kept (0 deg =
    no mixing, 90 deg = the two orbitals swap).
    """
    matrix = np.array(coefficients, dtype=float).copy()
    if matrix.ndim != 2 or matrix.shape[1] == 0:
        raise PerturbError("the coefficient matrix is not a two-dimensional grid.")
    n_mo = matrix.shape[1]
    for pair in plan.pairs:
        if not (0 <= pair.occupied < n_mo) or not (0 <= pair.virtual < n_mo):
            raise PerturbError(
                f"the pair ({pair.occupied}, {pair.virtual}) is outside the "
                f"{n_mo}-orbital coefficient matrix."
            )
        if pair.occupied == pair.virtual:
            raise PerturbError("a mixing pair needs two different orbitals.")
        angle = np.deg2rad(pair.angle_deg)
        old_occ = matrix[:, pair.occupied].copy()
        old_vir = matrix[:, pair.virtual].copy()
        matrix[:, pair.occupied] = np.cos(angle) * old_occ + np.sin(angle) * old_vir
        matrix[:, pair.virtual] = np.cos(angle) * old_vir - np.sin(angle) * old_occ
    return tuple(tuple(float(value) for value in row) for row in matrix)


@dataclass(frozen=True)
class PerturbedStart:
    """One perturbed start: the file payload and its full record."""

    mkl: MklFile
    alpha: PerturbationPlan
    beta: PerturbationPlan | None


def perturb_mkl(
    reference: MklFile,
    *,
    seed: int,
    n_pairs: int = 10,
    window: int = 15,
    perturb_beta: bool = True,
) -> PerturbedStart:
    """Perturb the reference's orbitals once (the seed fixes every draw).

    An unrestricted reference is perturbed in both spin blocks independently
    (``perturb_beta=False`` mixes the alpha block alone); a restricted
    reference has only the alpha block and the flag is ignored.  Occupations
    and the stored orbital energies travel unchanged: the perturbed columns
    are no longer eigenfunctions, and ORCA re-determines everything after
    reading the guess.
    """
    rng = random.Random(seed)
    alpha_plan = plan_perturbation(
        reference.occupations,
        tuple(value for group in reference.groups for value in group.energies),
        n_pairs=n_pairs,
        window=window,
        rng=rng,
    )
    alpha_mixed = mix_orbitals(reference.coefficients(), alpha_plan)
    beta_plan = None
    beta_mixed = None
    if reference.unrestricted and perturb_beta:
        beta_plan = plan_perturbation(
            reference.beta_occupations,
            tuple(value for group in reference.beta_groups for value in group.energies),
            n_pairs=n_pairs,
            window=window,
            rng=rng,
        )
        beta_mixed = mix_orbitals(reference.beta_coefficients(), beta_plan)
    try:
        updated = reference.with_orbitals(
            alpha_mixed,
            beta_coefficients=beta_mixed if beta_mixed is not None else None,
        )
    except MklError as exc:
        raise PerturbError(f"the reference mkl does not accept the perturbation: {exc}") from exc
    return PerturbedStart(mkl=updated, alpha=alpha_plan, beta=beta_plan)


_MOINP_RE = re.compile(r"^\s*%moinp\b(.*)$", re.IGNORECASE)
_MOREAD_RE = re.compile(r"\bmoread\b", re.IGNORECASE)
#: An inline coordinate header -- ``* xyz <charge> <mult>``; ``* xyzfile`` is a
#: different directive and must not match (the alphabetic token right after
#: ``xyz`` is what distinguishes them).
_XYZ_HEADER_RE = re.compile(r"^\s*\*\s*xyz\s+[-\d]", re.IGNORECASE)


def mo_read_variant(input_text: str, gbw_name: str) -> str:
    """The base input with ``!MORead`` and ``%moinp "<gbw_name>"`` in place.

    An existing ``%moinp`` line is replaced; otherwise the line is inserted
    after the last leading ``!`` line (before the first block or coordinate
    line).  ``MORead`` is appended to the first simple line when absent.
    """
    lines = input_text.splitlines()
    if not lines:
        raise PerturbError("the base input is empty.")
    simple_index = next(
        (i for i, line in enumerate(lines) if line.lstrip().startswith("!")), None
    )
    if simple_index is None:
        raise PerturbError(
            "the base input carries no simple '!' line, so MORead cannot be added. "
            "Next step: give the ORCA input that produced the reference run."
        )
    if not any(
        _MOREAD_RE.search(line)
        for line in lines
        if line.lstrip().startswith("!")
    ):
        lines[simple_index] = lines[simple_index].rstrip() + " MORead"
    moinp_line = f'%moinp "{gbw_name}"'
    replaced = False
    for i, line in enumerate(lines):
        if _MOINP_RE.match(line):
            lines[i] = moinp_line
            replaced = True
            break
    if not replaced:
        insertion = simple_index + 1
        while insertion < len(lines):
            stripped = lines[insertion].strip()
            if stripped.startswith("!") or not stripped:
                insertion += 1
                continue
            break
        lines.insert(insertion, moinp_line)
    return "\n".join(lines) + "\n"


def check_geometry_match(input_text: str, reference: MklFile) -> str | None:
    """Compare an inline ``* xyz`` block against the mkl's atoms.

    Returns ``None`` when the coordinates agree (same order, 1e-4 Angstrom),
    a note when the input has no inline xyz block (nothing was checked), and
    raises when they disagree -- a perturbed start must run at the reference
    geometry, and the likeliest mistake is a base input from another point.
    """
    lines = input_text.splitlines()
    start = None
    for i, line in enumerate(lines):
        if _XYZ_HEADER_RE.match(line):
            start = i + 1
            break
    if start is None:
        return (
            "the base input has no inline '* xyz' block (an external coordinate "
            "file?), so the geometry was not cross-checked against the reference mkl"
        )
    atoms = []
    for line in lines[start:]:
        if line.strip().startswith("*"):
            break
        tokens = line.split()
        if not tokens:
            continue
        if len(tokens) != 4:
            raise PerturbError(
                f"the base input's coordinate line {line!r} is not '<element> x y z'."
            )
        atoms.append(tokens)
    if len(atoms) != len(reference.atoms):
        raise PerturbError(
            f"the base input holds {len(atoms)} atom(s) while the reference mkl holds "
            f"{len(reference.atoms)}. Next step: give the reference run's own input."
        )
    for index, (tokens, (z, x, y, zz)) in enumerate(zip(atoms, reference.atoms)):
        try:
            z_in = element_z(tokens[0])
        except ElementError as exc:
            raise PerturbError(f"the base input's atom {tokens[0]!r} is unknown.") from exc
        if z_in != z or max(
            abs(float(tokens[1]) - x), abs(float(tokens[2]) - y), abs(float(tokens[3]) - zz)
        ) > 1e-4:
            symbol = next(
                (sym for sym, value in ELEMENT_Z.items() if value == z), str(z)
            )
            raise PerturbError(
                f"the base input's atom {index + 1} ({tokens[0]} {tokens[1]} "
                f"{tokens[2]} {tokens[3]}) disagrees with the reference mkl "
                f"({symbol} {x:.6f} {y:.6f} {zz:.6f}). Next step: a perturbed start "
                "must run at the reference geometry -- give the reference run's input."
            )
    return None


def render(start: PerturbedStart, *, seed: int, index: int) -> str:
    """One start's record: the windows, every pair and every angle."""
    lines = [
        f"perturbed start {index} (seed {seed}):",
        f"  occupied window (mkl columns): {', '.join(str(i) for i in start.alpha.window_occupied)}",
        f"  virtual window  (mkl columns): {', '.join(str(i) for i in start.alpha.window_virtual)}",
        f"  alpha mixing ({len(start.alpha.pairs)} pair(s)): "
        + "; ".join(
            f"{pair.occupied}<-{pair.virtual} {pair.angle_deg:.2f} deg"
            for pair in start.alpha.pairs
        ),
    ]
    if start.beta is not None:
        lines.append(
            f"  beta  mixing ({len(start.beta.pairs)} pair(s)): "
            + "; ".join(
                f"{pair.occupied}<-{pair.virtual} {pair.angle_deg:.2f} deg"
                for pair in start.beta.pairs
            )
        )
    return "\n".join(lines)


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the perturbation algorithm and of the write-back route."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "The randomized occupied-virtual mixing of the converged orbitals "
                "(10 random pairs out of the 15 highest occupied and 15 lowest "
                "unoccupied orbitals, random angles in [0, 90) degrees; their "
                "Eqs. (2)-(3)) as an inexpensive verification that an SCF solution "
                "is not a non-global minimum: the perturbed restart either returns "
                "to the same solution or finds another one, and a lower energy "
                "identifies the reference as wrongly converged.  The source notes "
                "that stability analysis detects unstable solutions but cannot "
                "distinguish local from global minima, and that the test cannot "
                "guarantee detection."
            ),
            ref=(
                "Vaucher & Reiher, J. Chem. Theory Comput. 2017, 13, 1219-1228, "
                "section IV (Eqs. (2)-(3)); the project's close reading "
                "scratch_wave41/vaucher2017.md"
            ),
            url="https://doi.org/10.1021/acs.jctc.7b00011",
            bibkey="vaucher2017steering",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "The write-back route this menu needs, measured on ORCA 6.1.1: "
                "``orca_2mkl <base> -mkl`` exports the converged orbitals as plain "
                "text, the perturbed file converts back with ``orca_2mkl <name> "
                "-gbw``, and ORCA accepts it as INITIAL GUESS: MOREAD.  Measured "
                "on the CH4 dissociation fixture: the guess propagated from the "
                "equilibrium orbitals keeps the restricted solution "
                "(-40.111230225721 Eh, <S**2> 0.000000) and the default ORCA guess "
                "also lands on a restricted solution (-40.183584827774 Eh, "
                "<S**2> 0.000000), while the perturbed starts find the "
                "broken-symmetry solution (-40.2416 Eh, <S**2> 0.971943) -- a "
                "healing of 0.130 Eh and a further 0.058 Eh below the default "
                "guess."
            ),
            ref="fixtures/orca/ch4_*; tests/test_perturb_guess.py",
        ),
    )
