"""Cross-run state tracking (menu 49) -- the post-processing form of the
density-matrix tracking measure of Tran, Shea and Neuscamman.

The source method (J. Chem. Theory Comput. 2019, 15, 4790) follows one
excited state through the macro iterations of a state-specific CASSCF
optimization by scoring every newly generated Davidson root with

    Q(c) = W0(c) + W1(c) + D(c),
    D(c) = ||Gamma_t - Gamma(c)||_F / n_CAS,

where ``Gamma_t`` is the target state's one-body reduced density matrix
(rotated into the current orbital basis before the difference, to reduce the
sensitivity to orbital changes) and ``W0 = (omega - E)^2`` aims the search at
the wanted energy.  The source's central finding is that the density-matrix
difference is a far more robust similarity measure than CI-vector dot
products, and that it is essential when two roots have similarly structured
densities.

This module is the post-processing form of that measure.  fBlockKit does not
take over the engine's optimization loop; instead, given a sequence of
already-computed runs of one geometry (each an ``orca_2json`` export with
per-root densities), it walks a chosen target state down the sequence and
picks, in every run, the root that matches the tracked state.  Three
substitutions against the source, all declared in the report:

- **W1 is not evaluated.**  It needs the active-to-virtual single-excitation
  coupling (the source notes its cost equals that of the CASPT2 first-order
  wave-function right-hand side), which no ORCA export carries.  For a
  converged candidate run W1 tends to zero -- the source's equivalence, W1 = 0
  iff the energy is stationary -- and converged runs are the intended input.
  The report does not fabricate the term.
- **the 1/n_CAS scaling of D is not applied.**  It serves the source's
  portability across active-space sizes; within one sequence the ranking is
  unaffected, and the report quotes the unscaled Frobenius difference.
- **omega is the tracked state's own energy** (the previous run's chosen
  root), updated each step -- the tracking semantics of the post-processing
  use.  The source's fixed-omega variant remains available by construction:
  the quantity enters only through W0.

Data contract (measured on ORCA 6.1.1): per-root CASSCF densities come from
``orca_2json`` with ``Densities: ["all"]`` (keys
``Tdens-CAS.mult.<m>.root.<k>.p``; a state-averaged export ships one per
root -- three for the N2 fixture); the 1-RDM convention is
``Gamma = C (S D_json S) C^T`` (pinned by the ASS1ST round suite); the
rotation into a candidate run's basis is ``M = C_cand S C_target^T``
(orthogonal to 2e-14 on the fixtures); per-root energies come from the
sibling ``.out`` ``CAS-SCF STATES`` block.  The geometry gate: the runs must
share their atoms and coordinates (tolerance 1e-6 Angstrom) and their AO
dimension -- the shared AO space the rotation presupposes.  The N2 probe
values the tests pin: within one run D(root 0, root 1) = 0.90035, D(root 0,
root 2) = 0.88714 and D(root 1, root 2) = 0.01325 (the near-degenerate pi
pair); across runs, the same state measures D = 0.05212 while the wrong
roots measure 0.90 and 0.89 -- a factor-17 separation.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, replace
from typing import Mapping, Sequence

import numpy as np

from ..knowledge.models import EVIDENCE_LITERATURE, EVIDENCE_MEASURED, Evidence
from ..parsers.orca_json import OrcaJson

__all__ = [
    "StateTrackingError",
    "TrackRoot",
    "TrackRun",
    "CandidateScore",
    "TrackStep",
    "TrackResult",
    "build_run",
    "track",
    "render",
    "evidence",
]

#: Atoms and coordinates must agree to this tolerance (Angstrom) for two runs
#: to count as one AO space (the coordinates are compared elementwise).
GEOMETRY_TOL = 1e-6

#: The per-root CASSCF density key of the orca_2json sidecar.
_ROOT_DENSITY_RE = re.compile(r"^Tdens-CAS\.mult\.(\d+)\.root\.(\d+)\.p$")

_JSON_URL = (
    "https://www.faccts.de/docs/orca/6.1/manual/contents/utilitiesvisualization/"
    "json.html"
)


class StateTrackingError(ValueError):
    """The sequence cannot be tracked as given (each message carries a next step)."""


@dataclass(frozen=True)
class TrackRoot:
    """One root of one run: its labels, its energy (if read) and its 1-RDM.

    ``root`` is ORCA's own 0-based root index (the export key's ``root.<k>``
    and the output's ``ROOT k`` line agree on it); ``gamma`` is the 1-RDM in
    that run's MO basis, ``Gamma = C (S D_json S) C^T``.
    """

    mult: int
    root: int
    energy: float | None
    gamma: tuple[tuple[float, ...], ...]


@dataclass(frozen=True)
class TrackRun:
    """One run of the sequence: its geometry signature and its roots."""

    name: str
    atoms: tuple[str, ...]
    coordinates: tuple[tuple[float, float, float], ...]
    n_mo: int
    n_ao: int
    coefficients: tuple[tuple[float, ...], ...]
    overlap: tuple[tuple[float, ...], ...]
    roots: tuple[TrackRoot, ...]

    def root_labels(self) -> tuple[tuple[int, int], ...]:
        """The (mult, root) labels, in the order the roots are held."""
        return tuple((root.mult, root.root) for root in self.roots)


@dataclass(frozen=True)
class CandidateScore:
    """One candidate root's score against the step's target."""

    mult: int
    root: int
    energy: float | None
    w0: float | None  # (E_target - E_root)^2 when both energies are known
    d: float  # Frobenius difference of the 1-RDMs (unscaled)
    q: float  # W0 + D (the same ordering as D when energies are absent)
    chosen: bool


@dataclass(frozen=True)
class TrackStep:
    """One step of the walk: every candidate of one run and the choice."""

    from_run: str
    run_name: str
    candidates: tuple[CandidateScore, ...]
    margin: float | None  # second-best minus best Q; None with a single candidate
    energy_used: bool  # whether W0 entered Q in this step


@dataclass(frozen=True)
class TrackResult:
    """The tracked lineage across the sequence."""

    runs: tuple[str, ...]
    target_mult: int
    target_root: int
    target_energy: float | None
    steps: tuple[TrackStep, ...]
    notes: tuple[str, ...]


def build_run(
    export: OrcaJson,
    *,
    energies: Mapping[tuple[int, int], float] | None = None,
    name: str | None = None,
) -> TrackRun:
    """Assemble one run of the sequence from its ``orca_2json`` export.

    ``energies`` maps ``(mult, root)`` to the per-root energy (the handler
    reads them from the sibling ``.out``); missing entries leave the energy
    unknown for that root, which the walk reports.
    """
    if export.overlap is None:
        raise StateTrackingError(
            f"{name or 'the export'} carries no S matrix, so the density "
            "rotation between runs cannot be built. Next step: regenerate the "
            'export with 1elIntegrals: ["S"] (the toolkit\'s export recipe).'
        )
    entries: list[tuple[int, int, tuple[tuple[float, ...], ...]]] = []
    for key, matrix in dict(export.densities or ()).items():
        match = _ROOT_DENSITY_RE.match(key)
        if match is not None:
            entries.append((int(match.group(1)), int(match.group(2)), matrix))
    if not entries:
        raise StateTrackingError(
            f"{name or 'the export'} has no Tdens-CAS.mult.<m>.root.<k>.p "
            "densities to track. Next step: regenerate the export with "
            'Densities: ["all"] on a CASSCF run that kept its densities '
            "(KeepDens)."
        )
    coefficients = np.asarray(export.mo_coefficients, dtype=float)
    overlap = np.asarray(export.overlap, dtype=float)
    roots: list[TrackRoot] = []
    for mult, root, matrix in sorted(entries):
        d_json = np.asarray(matrix, dtype=float)
        gamma = coefficients @ (overlap @ d_json @ overlap) @ coefficients.T
        energy = None if energies is None else energies.get((mult, root))
        roots.append(
            TrackRoot(
                mult=mult,
                root=root,
                energy=energy,
                gamma=tuple(tuple(row) for row in gamma),
            )
        )
    return TrackRun(
        name=name or "run",
        atoms=tuple(export.atoms),
        coordinates=tuple(tuple(float(v) for v in xyz) for xyz in export.coordinates),
        n_mo=export.n_mo,
        n_ao=export.n_ao,
        coefficients=tuple(tuple(row) for row in coefficients),
        overlap=tuple(tuple(row) for row in overlap),
        roots=tuple(roots),
    )


def _check_shared_space(base: TrackRun, run: TrackRun) -> None:
    """Refuse runs that do not share the base run's AO space."""
    if run.n_ao != base.n_ao or run.n_mo != base.n_mo:
        raise StateTrackingError(
            f"run {run.name!r} has {run.n_ao} AOs / {run.n_mo} MOs while "
            f"{base.name!r} has {base.n_ao} / {base.n_mo}; the density rotation "
            "needs one AO space. Next step: track runs with the same basis and "
            "beyond that the same geometry."
        )
    if tuple(run.atoms) != tuple(base.atoms):
        raise StateTrackingError(
            f"run {run.name!r} has atoms {tuple(run.atoms)} while {base.name!r} "
            f"has {tuple(base.atoms)}; the density rotation presupposes one AO "
            "space. Next step: track runs of one geometry (cross-geometry "
            "tracking is outside the current scope)."
        )
    delta = np.asarray(run.coordinates, dtype=float) - np.asarray(
        base.coordinates, dtype=float
    )
    worst = float(np.abs(delta).max()) if delta.size else 0.0
    if worst > GEOMETRY_TOL:
        raise StateTrackingError(
            f"run {run.name!r} sits at coordinates up to {worst:.3e} Angstrom "
            f"away from {base.name!r}; the density rotation presupposes one AO "
            "space. Next step: track runs of one geometry (cross-geometry "
            "tracking is outside the current scope)."
        )


def track(runs: Sequence[TrackRun], target_root: int) -> TrackResult:
    """Walk one target state down the sequence.

    ``target_root`` is the 0-based root index of the first run (the number
    ORCA prints on its ``ROOT`` lines and the export's ``root.<k>`` keys use).
    Each following run is scored candidate by candidate against the tracked
    state's rotated density; the best Q continues as the new target.
    """
    if len(runs) < 2:
        raise StateTrackingError(
            "a tracking sequence needs at least two runs, got "
            f"{len(runs)}. Next step: give the runs in tracking order (comma "
            "separated), first the one that holds the target state."
        )
    base = runs[0]
    notes: list[str] = []
    matches = [root for root in base.roots if root.root == target_root]
    if not matches:
        available = ", ".join(f"(mult {m}, root {r})" for m, r in base.root_labels())
        raise StateTrackingError(
            f"the first run {base.name!r} has no root {target_root}; its roots "
            f"are {available}. Next step: pick one of the listed root numbers "
            "(0-based, as ORCA prints them)."
        )
    target = min(matches, key=lambda root: root.mult)
    if len(matches) > 1:
        notes.append(
            f"root {target_root} appears in multiplicities "
            f"{[root.mult for root in matches]}; the lowest was tracked."
        )
    for run in runs[1:]:
        _check_shared_space(base, run)

    target_gamma = np.asarray(target.gamma, dtype=float)
    target_coeff = np.asarray(base.coefficients, dtype=float)
    target_energy = target.energy
    steps: list[TrackStep] = []
    previous_name = base.name
    for run in runs[1:]:
        overlap = np.asarray(run.overlap, dtype=float)
        coefficients = np.asarray(run.coefficients, dtype=float)
        rotation = coefficients @ overlap @ target_coeff.T
        rotated_target = rotation @ target_gamma @ rotation.T
        energies_used = target_energy is not None and all(
            root.energy is not None for root in run.roots
        )
        scores: list[CandidateScore] = []
        for root in run.roots:
            difference = float(
                np.linalg.norm(rotated_target - np.asarray(root.gamma, dtype=float))
            )
            w0 = (
                float((target_energy - root.energy) ** 2)
                if energies_used and root.energy is not None and target_energy is not None
                else None
            )
            q = difference + w0 if w0 is not None else difference
            scores.append(
                CandidateScore(
                    mult=root.mult,
                    root=root.root,
                    energy=root.energy,
                    w0=w0,
                    d=difference,
                    q=q,
                    chosen=False,
                )
            )
        best = min(range(len(scores)), key=lambda index: scores[index].q)
        scores[best] = replace(scores[best], chosen=True)
        ordered = sorted(score.q for score in scores)
        margin = ordered[1] - ordered[0] if len(ordered) > 1 else None
        steps.append(
            TrackStep(
                from_run=previous_name,
                run_name=run.name,
                candidates=tuple(scores),
                margin=margin,
                energy_used=energies_used,
            )
        )
        chosen = run.roots[best]
        target_gamma = np.asarray(chosen.gamma, dtype=float)
        target_coeff = coefficients
        target_energy = chosen.energy
        previous_name = run.name

    return TrackResult(
        runs=tuple(run.name for run in runs),
        target_mult=target.mult,
        target_root=target.root,
        target_energy=target.energy,
        steps=tuple(steps),
        notes=tuple(notes),
    )


def render(result: TrackResult, *, source: str) -> str:
    """The menu-49 report: the lineage, the per-step tables and the declarations."""
    out = [f"Cross-run state tracking (source: {source})", ""]
    out.append("Runs (in tracking order):")
    for index, name in enumerate(result.runs, start=1):
        out.append(f"  {index}. {name}")
    start_energy = (
        f"{result.target_energy:.6f} Eh"
        if result.target_energy is not None
        else "energy not read"
    )
    out += [
        f"Target: root {result.target_root} (mult {result.target_mult}) of run 1; "
        f"{start_energy}",
        "",
    ]
    for index, step in enumerate(result.steps, start=1):
        out.append(
            f"Step {index} -> {index + 1} ({step.run_name}): candidates against the "
            "tracked density, rotated into this run's basis"
        )
        out.append(
            "  mult  root  energy (Eh)     W0          D           Q"
        )
        for score in step.candidates:
            energy = f"{score.energy:.6f}" if score.energy is not None else "n/a"
            w0 = f"{score.w0:.4e}" if score.w0 is not None else "n/a"
            tail = "  [chosen]" if score.chosen else ""
            out.append(
                f"  {score.mult:<4d}  {score.root:<4d}  {energy:<14s}  "
                f"{w0:<10s}  {score.d:<10.4e}  {score.q:<10.4e}{tail}"
            )
        if not step.energy_used:
            out.append(
                "  (per-root energies were incomplete for this step; Q = D, the "
                "density-difference ordering)"
            )
        if step.margin is not None:
            margin_note = ""
            if step.margin < 1e-3:
                margin_note = " -- the candidates are nearly degenerate in this metric"
            out.append(
                f"  margin over the runner-up: {step.margin:.4e}{margin_note}"
            )
        else:
            out.append("  a single candidate: the margin is undefined")
        out.append("")
    if result.notes:
        out.append("Notes:")
        for note in result.notes:
            out.append(f"  - {note}")
        out.append("")
    out += [
        "Boundaries (declared substitutions against the source's criterion",
        "Q = W0 + W1 + D, J. Chem. Theory Comput. 2019, 15, 4790):",
        "  - the W1 term (the active-to-virtual stationarity measure) is not",
        "    evaluated: no ORCA export carries the coupling it needs. For a",
        "    converged candidate run W1 tends to zero (the source's equivalence:",
        "    W1 = 0 iff the energy is stationary), which is the intended input;",
        "  - the source's 1/n_CAS scaling of D is not applied; the ranking within",
        "    one sequence is unaffected and the values here are unscaled",
        "    Frobenius norms of the full-space 1-RDM difference;",
        "  - omega (the source's energy target) is taken as the tracked state's",
        "    own energy, updated each step;",
        "  - the geometry gate: the runs must share atoms, coordinates (1e-6",
        "    Angstrom) and the AO dimension -- cross-geometry tracking is outside",
        "    the current scope.",
    ]
    return "\n".join(out)


def evidence() -> tuple[Evidence, ...]:
    """Provenance: the source paper, the export schema, the measured anchors."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "The tracking measure Q = W0 + W1 + D with "
                "D = ||Gamma_t - Gamma(c)||_F / n_CAS, the rotation of the "
                "target density into the current orbital basis, and the finding "
                "that the density-matrix difference beats CI-vector overlaps "
                "(essential for roots with similarly structured densities) are "
                "from Tran, Shea and Neuscamman."
            ),
            ref=(
                "Tran, L. N.; Shea, J. A. R.; Neuscamman, E. J. Chem. Theory "
                "Comput. 2019, 15 (9), 4790-4803"
            ),
            bibkey="tran2019tracking",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "Measured on ORCA 6.1.1 (fixtures n2_ass1st_sa.*, "
                "n2_ass1st.*): the state-averaged orca_2json export carries one "
                "Tdens-CAS density per root (three for the SA(6,6) run); the "
                "1-RDM convention Gamma = C (S D_json S) C^T reproduces the "
                "printed natural occupations (the ASS1ST suite's pin); the "
                "cross-run rotation M = C_cand S C_target^T comes out "
                "orthogonal to 2e-14 with S identical between runs of one "
                "geometry; the same state measures D = 0.05212 across runs "
                "while the wrong roots measure 0.90 / 0.89, and the N2 pi pair "
                "sits at D(root 1, root 2) = 0.01325."
            ),
            ref="tests/test_state_tracking.py; fixtures/orca/n2_ass1st_sa.*, n2_ass1st.*",
        ),
    )
