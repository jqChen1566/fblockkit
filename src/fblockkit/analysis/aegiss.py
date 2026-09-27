"""3.9: AEGISS -- entropy screening joined to an atomic-orbital projection.

AEGISS (Tarocco, Haase, Pavošević, Krishna, Guidoni, Knecht & Stella, arXiv
2026) combines the two complementary selection ideas: *how correlated* an
orbital is (the single-orbital entropy of AutoCAS) and *whether it belongs to
the chemistry of interest* (the atomic-orbital projection of AVAS).  The
five-step workflow:

1. (optional) label the molecule's atoms into chemical clusters;
2. a mean-field calculation defines the orbital basis;
3. **entropy screening**: from a correlated wavefunction, keep every orbital
   whose four-state entropy exceeds ``tau_E = tau * S_max`` (the source's
   default ``tau = 0.1`` -- its AutoCAS-style 10 % rule; the benzene example
   used 0.2);
4. (optional) define one or more AO groups from chemical knowledge;
5. **AO projection**: for each AO group D, the weight of an entropy-screened
   orbital is its projection onto the AO family of D; keep the orbitals above
   ``epsilon_D = 0.5`` (the source's benzene and ferrocene values).

The final active space is the union of the per-group selections.

Measured deviation on the projection weight
-------------------------------------------

The source writes the weight as the *signed sum* of the projected-overlap
rows, ``w_p = sum_eta [O]_{eta p}``.  Measured on the benzene fixture, that
signed sum cancels by symmetry for every pi orbital with a nodal plane (only
the nodeless a2u survives: weights -3.54, 0, 0, 0, 0, 0 over the six C 2pz
targets) -- an artifact of summing signed rows, not a chemistry statement.
This tool therefore uses the projection *norm* ``w_p = sum_eta |[O]_{eta p}|^2``
(the standard AVAS reading of the same overlap matrix), on which the benzene
pi manifold separates cleanly: the sigma orbitals measure exactly 0.0000 and
the six pi orbitals 0.047-2.087 (with the label's shell resolution deciding
how much of each).  The threshold semantics are kept (epsilon_D = 0.5 against
these weights); the carrier of the deviation is recorded with every report.

This tool's variant, stated honestly
------------------------------------

The entropy comes from the *exact* determinant CI of the FCIDUMP route (the
four-state entropy machinery, menu 12), over the FCIDUMP's orbital window --
the source runs a DMRG (or another correlated method) over a larger window and
estimates the entropy there.  The AO projection follows the source's
Equations (5)-(7), but targets the *calculation's own AO subset* (the
non-minimal route of the AVAS implementation, menu 14 -- the only route open
to the f block); omitting the shell index from the AO label matches every
shell of that (element, angular momentum, component) family, which is the
analogue of the source's minimal-basis ``2pz`` label.  The cluster labeling
(step 1) is folded into the AO-label choice, and the multi-group union is
supported by running the menu once per group.

The regression anchors are the benzene/cc-pVDZ pi platform (fixtures
``benzene.json`` + ``benzene.fcidump``, produced from an AVAS-prepped CASSCF
run): the exact window FCI energy reproduces the engine's printed CASSCF
energy to eleven digits (-230.793818898, the source's own value -230.793770 to
5e-5); the exact entropies [0.174-0.342] all clear both the 10 % and the 20 %
screen lines; and the projection on the label ``C pz`` (the whole family)
keeps all six pi orbitals above 0.5, recovering the textbook (6e, 6o) space.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np

from ..knowledge.models import EVIDENCE_LITERATURE, EVIDENCE_MEASURED, Evidence, ReportSection
from ..parsers.fcidump import Fcidump
from ..parsers.orca_json import OrcaJson
from . import entropy_rdm
from .entropy_rdm import EntropyRdmError

__all__ = [
    "AEGISS_DEFAULT_TAU",
    "AEGISS_DEFAULT_EPSILON",
    "AegissError",
    "AoGroup",
    "AegissResult",
    "parse_ao_group",
    "matching_aos",
    "entropy_screen",
    "projection_weights",
    "analyze",
    "render",
    "run",
    "evidence",
]

#: The source's AutoCAS-style entropy line: tau_E = 0.1 * S_max (its benzene
#: example used 0.2).
AEGISS_DEFAULT_TAU = 0.1

#: The source's projection threshold (its benzene and ferrocene values).
AEGISS_DEFAULT_EPSILON = 0.5

#: The AO-label grammar: "<element> [shell]<angular>[component]", e.g. "C 2pz",
#: "Fe 1d", "N p" (no shell index = every shell of the family).
_GROUP_RE = re.compile(r"^(?P<shell>\d*)(?P<angular>[spdfghi])(?P<component>[A-Za-z0-9+-]*)$")


class AegissError(ValueError):
    """The selection cannot proceed on the given data (with a next step)."""


@dataclass(frozen=True)
class AoGroup:
    """One AO-family label: element plus shell/angular/component spelling."""

    raw: str
    element: str
    shell: int | None
    angular: str
    component: str | None


@dataclass(frozen=True)
class AegissResult:
    """The workflow trace, ready to render."""

    base_name: str
    label: AoGroup
    n_group_aos: int
    tau: float
    tau_line: float
    epsilon: float
    window: tuple[int, ...]
    occupations: tuple[float, ...]
    entropies: tuple[float, ...]
    kept: tuple[int, ...]  # window positions passing the entropy screen
    weights: tuple[float, ...]  # projection norm per kept orbital
    selected: tuple[int, ...]  # window positions passing both screens
    n_electrons: int
    n_orbitals: int
    energy_fci: float
    engine_energy: float | None


# --- the AO label -------------------------------------------------------------


def parse_ao_group(text: str) -> AoGroup:
    """Parse the AO-family label (``C 2pz``, ``Fe d``, ``N p``)."""
    tokens = text.split()
    if len(tokens) != 2:
        raise AegissError(
            f"the AO label {text!r} does not split into '<element> [shell]<angular>"
            "[component]' (for example 'C 2pz' or 'Fe d'). Next step: give the element "
            "and the angular-momentum spelling."
        )
    element, second = tokens
    match = _GROUP_RE.match(second)
    if match is None:
        raise AegissError(
            f"the AO label {text!r}: {second!r} does not parse as [shell]<angular>"
            "[component]. Next step: examples are '2pz', 'pz', 'd'."
        )
    shell = int(match.group("shell")) if match.group("shell") else None
    component = match.group("component") or None
    return AoGroup(
        raw=text,
        element=element,
        shell=shell,
        angular=match.group("angular"),
        component=component,
    )


def matching_aos(labels, group: AoGroup) -> tuple[int, ...]:
    """The AO columns matching the label family (shell and component optional)."""
    hits = []
    for index, label in enumerate(labels):
        if label.element.lower() != group.element.lower():
            continue
        if label.angular != group.angular:
            continue
        if group.component is not None and label.component != group.component:
            continue
        if group.shell is not None and label.shell != group.shell:
            continue
        hits.append(index)
    return tuple(hits)


# --- the two screens ----------------------------------------------------------


def entropy_screen(entropies: tuple[float, ...], tau: float) -> tuple[tuple[int, ...], float]:
    """Keep the orbitals with S > tau * S_max (the source's relative line)."""
    if not 0.0 < tau < 1.0:
        raise AegissError(f"the entropy fraction {tau} is not in (0, 1).")
    if not entropies:
        raise AegissError("the entropy list is empty; nothing to screen.")
    peak = max(entropies)
    if peak <= 0.0:
        raise AegissError(
            "every single-orbital entropy is zero, so the relative line is undefined. "
            "Next step: check that the FCIDUMP describes a correlated window (a "
            "closed-shell single determinant has a zero entropy spectrum)."
        )
    line = tau * peak
    kept = tuple(i for i, value in enumerate(entropies) if value > line)
    if not kept:
        raise AegissError(
            f"no orbital exceeds the entropy line {line:.6g} = {tau:g} * {peak:.6g}; "
            "Next step: lower the entropy fraction (the source's default is 0.1)."
        )
    return kept, line


def projection_weights(
    export: OrcaJson, window: tuple[int, ...], aos: tuple[int, ...]
) -> tuple[float, ...]:
    """The projection norm of each window orbital onto the AO family."""
    if export.overlap is None:
        raise AegissError(
            "the export carries no S-Matrix, which the projection needs. Next step: "
            're-export with \'{"MOCoefficients": true, "1elIntegrals": ["S"]}\'.'
        )
    coefficients = np.asarray(export.mo_coefficients)
    overlap = np.asarray(export.overlap)
    projected = overlap[np.ix_(aos, range(export.n_ao))] @ coefficients[list(window)].T
    return tuple(float(value) for value in (projected**2).sum(axis=0))


# --- the workflow -------------------------------------------------------------


def analyze(
    export: OrcaJson,
    dump: Fcidump,
    *,
    label: str,
    tau: float = AEGISS_DEFAULT_TAU,
    epsilon: float = AEGISS_DEFAULT_EPSILON,
    reference_energy: float | None = None,
) -> AegissResult:
    """One AEGISS selection: exact window entropies, then the AO projection."""
    if export.ao_labels is None:
        raise AegissError(
            "the export carries no OrbitalLabels block, which the AO projection needs. "
            'Next step: re-export with \'{"MOCoefficients": true, "1elIntegrals": '
            '["S"]}\' (the labels come with the coefficients).'
        )
    if not 0.0 < epsilon < 1.0:
        raise AegissError(f"the projection threshold {epsilon} is not in (0, 1).")
    window = tuple(
        i for i, occ in enumerate(export.mo_occupations) if 1e-6 < occ < 2.0 - 1e-6
    )
    if len(window) != dump.norb:
        raise AegissError(
            f"the export has {len(window)} fractionally occupied orbitals while the "
            f"FCIDUMP carries {dump.norb}: the two files do not belong to the same run "
            "and window. Next step: export the gbw whose orbitals were dumped."
        )
    try:
        state = entropy_rdm.solve_fci(dump, reference_energy=reference_energy)
        densities = entropy_rdm.spin_densities(state)
    except EntropyRdmError as exc:
        raise AegissError(str(exc)) from exc
    entropies = entropy_rdm.entropy_spectrum(densities)
    if len(entropies) != len(window):
        raise AegissError(
            "internal bookkeeping defect: the entropy list does not match the window."
        )

    group = parse_ao_group(label)
    aos = matching_aos(export.ao_labels, group)
    if not aos:
        raise AegissError(
            f"the AO label {label!r} matches no orbital label of the export. Next step: "
            "check the element and the angular spelling against the export's "
            "OrbitalLabels (menu 15 prints them)."
        )

    kept, line = entropy_screen(entropies, tau)
    weights = projection_weights(export, window, aos)
    selected = tuple(i for i in kept if weights[i] > epsilon)
    if not selected:
        raise AegissError(
            f"no entropy-screened orbital exceeds the projection threshold {epsilon:g} "
            f"for the label {label!r}. Next step: lower the threshold, or widen the AO "
            "label (for example drop the shell index)."
        )
    n_electrons = int(round(sum(export.mo_occupations[window[i]] for i in selected)))
    return AegissResult(
        base_name=export.base_name,
        label=group,
        n_group_aos=len(aos),
        tau=tau,
        tau_line=line,
        epsilon=epsilon,
        window=window,
        occupations=tuple(float(export.mo_occupations[i]) for i in window),
        entropies=entropies,
        kept=kept,
        weights=weights,
        selected=selected,
        n_electrons=n_electrons,
        n_orbitals=len(selected),
        energy_fci=float(state.energy_total),
        engine_energy=reference_energy,
    )


# --- output -------------------------------------------------------------------


def render(result: AegissResult) -> str:
    """The five-step trace: entropy screen, projection, and the final space."""
    lines = [
        "AEGISS active-space selection (entropy screening + atomic-orbital projection):",
        f"  system: {result.base_name}; window {len(result.window)} orbitals, "
        f"AO label {result.label.raw!r} ({result.n_group_aos} target functions)",
        f"  entropy line: S > {result.tau_line:.5f} = {result.tau:g} * S_max; "
        f"projection threshold: w > {result.epsilon:g} (projection norm)",
        "",
        f"  {'MO':>4}  {'occupation':>10}  {'S':>8}  {'entropy':>8}  {'weight':>8}  {'selected':>8}",
    ]
    for position, mo in enumerate(result.window):
        entropy_mark = "keep" if position in result.kept else "drop"
        if position in result.kept:
            weight_text = f"{result.weights[position]:>8.4f}"
            selected_mark = "yes" if position in result.selected else ""
        else:
            weight_text = f"{'-':>8}"
            selected_mark = ""
        lines.append(
            f"  {mo:>4}  {result.occupations[position]:>10.5f}  "
            f"{result.entropies[position]:>8.4f}  {entropy_mark:>8}  {weight_text}  "
            f"{selected_mark:>8}"
        )
    lines += [
        "",
        f"  final active space: ({result.n_electrons}e, {result.n_orbitals}o) -- "
        f"{len(result.kept)} orbital(s) passed the entropy screen, "
        f"{result.n_orbitals} of those passed the projection",
        f"  exact window FCI energy: {result.energy_fci:.9f} Eh"
        + (
            ""
            if result.engine_energy is None
            else f"  (engine CASSCF print: {result.engine_energy:.9f}; difference "
            f"{abs(result.energy_fci - result.engine_energy):.2e})"
        ),
        "",
        "Boundaries and checks:",
    ]
    for item in (
        "the entropies are exact (the four-state-entropy route over the FCIDUMP "
        "window); the source estimates them from a DMRG over a larger window",
        "the projection weight is the projection norm sum_eta |O[eta, p]|^2 -- "
        "a measured deviation from the source's signed row sum, which cancels by "
        "symmetry for nodal pi orbitals (on the benzene fixture only a2u survives "
        "the signed sum, while the norm separates sigma (0.0000) from pi cleanly)",
        "the target is the calculation's own AO subset (the non-minimal AVAS route); "
        "omitting the shell index from the label matches every shell of that "
        "(element, angular, component) family -- the analogue of the source's "
        "minimal-basis label",
        "the cluster labeling of the source's step 1 is folded into the AO label; "
        "multi-group unions are run one group per menu pass",
        "the deliverable is the selection over the dumped window; feeding it to a "
        "CASSCF needs the orbital-order machinery (menu 22's boundary note)",
    ):
        lines.append(f"  - {item}")
    return "\n".join(lines)


def run(export: OrcaJson, dump: Fcidump, **kwargs) -> ReportSection:
    """The analyser entry point for menu 25."""
    return ReportSection(
        title="3.9 AEGISS selection (entropy + AO projection)",
        body=render(analyze(export, dump, **kwargs)),
    )


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the workflow, the thresholds and the measured deviation."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "AEGISS: entropy screening (keep S > tau_E = tau * S_max, the AutoCAS "
                "relative line; tau = 0.1 default, the benzene example used 0.2) joined "
                "to an AVAS-style atomic-orbital projection (keep the projection weight "
                "above epsilon_D = 0.5, the source's benzene and ferrocene value); the "
                "final space is the union of the AO-group selections. Source anchors: "
                "benzene -> textbook (6e, 6o) pi space with CASSCF -230.793770 Eh; "
                "ferrocene -> (10e, 7o) matching the reference AVAS paper."
            ),
            ref=(
                "arXiv:2508.10671v3, section 3 (the five-step workflow, eqs. (4)-(7)) "
                "and section 4.1-4.2 (the benzene and ferrocene settings)"
            ),
            url="https://doi.org/10.48550/arXiv.2508.10671",
            bibkey="tarocco2026aegiss",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "Regression anchors on the benzene/cc-pVDZ pi platform: the exact window "
                "FCI energy reproduces the engine's printed CASSCF energy to eleven "
                "digits (-230.793818898; the source's own value -230.793770 agrees to "
                "5e-5); the exact entropies [0.174, 0.342] clear both the 10 % and the "
                "20 % screen lines; the projection on the label 'C pz' (the whole pz "
                "family) keeps all six pi orbitals above 0.5, recovering (6e, 6o). "
                "Measured deviation from the source's weight definition: the signed row "
                "sum cancels for nodal pi orbitals (only a2u survives on benzene), so "
                "the projection norm is used instead."
            ),
            ref="tests/test_aegiss.py; fixtures/orca/benzene.json, benzene.fcidump",
        ),
    )
