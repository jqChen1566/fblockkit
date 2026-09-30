"""G9: the CREST-ensemble upgrade inputs -- re-ranking the conformers at a higher level.

Menu 43's first mode reads a CREST ensemble (the GFN2 pre-screen); this
module writes the *second stage* of the workflow: one ORCA optimisation
input per selected conformer, so the pre-screened ensemble can be re-ranked
at a production level (the toolkit's Tier-2 line by default, r2SCAN-3c).

The inputs carry each conformer's ensemble geometry verbatim (the double
format of the house generators), the user's method line under ``! Opt``,
and the charge/multiplicity; they are written into an ``upgrade/``
subdirectory of the run, one file per conformer in the ensemble's own
(energy-ascending) order, so a batch submission is a plain loop over the
directory.  The report lists each input with its ensemble energy and
relative energy and states the run guidance.

Boundaries, stated once: the pre-screen's force-field ranking is what the
selection inherits (this module re-ranks nothing -- it prepares the inputs
that will); the re-ranking step itself (reading the upgraded outputs back
and comparing) is registered, and the conformer geometries are used as-is
(no symmetry or constraint handling is added).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..knowledge.models import EVIDENCE_MEASURED, Evidence
from ..parsers.crest import CrestEnsemble

__all__ = [
    "UpgradeError",
    "UpgradeInput",
    "UpgradePlan",
    "write_upgrade_inputs",
    "render_plan",
    "evidence",
]

#: The default upgrade level: the project's Tier-2 line for organic
#: conformer re-ranking.
DEFAULT_METHOD = "r2SCAN-3c"


class UpgradeError(ValueError):
    """The upgrade inputs cannot be written (with a next step)."""


@dataclass(frozen=True)
class UpgradeInput:
    """One generated optimisation input and its ensemble provenance."""

    name: str  # file name inside the upgrade directory
    index: int  # 1-based ensemble position
    energy_Eh: float
    relative_kcal: float
    text: str


@dataclass(frozen=True)
class UpgradePlan:
    """The written batch and everything the report needs to say about it."""

    directory: Path
    method: str
    charge: int
    multiplicity: int
    inputs: tuple[UpgradeInput, ...]


def write_upgrade_inputs(
    ensemble: CrestEnsemble,
    directory,
    *,
    method: str = DEFAULT_METHOD,
    charge: int = 0,
    multiplicity: int = 1,
    count: int | None = None,
) -> UpgradePlan:
    """Write one ``! Opt`` input per selected conformer into ``directory``.

    ``count``: how many conformers (lowest first) to upgrade; ``None`` takes
    the whole ensemble.
    """
    method = method.strip()
    if not method:
        raise UpgradeError(
            "the method line is empty. Next step: give the upgrade level, e.g. "
            f"'{DEFAULT_METHOD}' or 'PBE0 D4 def2-TZVP'."
        )
    frames = ensemble.frames
    if not frames:
        raise UpgradeError(
            "the ensemble carries no conformers. Next step: check the crest_conformers.xyz "
            "file."
        )
    if count is not None and not 1 <= count <= len(frames):
        raise UpgradeError(
            f"the upgrade count {count} is outside the ensemble's {len(frames)} "
            "conformer(s). Next step: give a count between 1 and the ensemble size "
            "(or leave it empty for the whole ensemble)."
        )
    selected = frames if count is None else frames[:count]
    target = Path(directory)
    target.mkdir(parents=True, exist_ok=True)
    inputs = []
    for frame in selected:
        name = f"conf_{frame.index:02d}.opt.inp"
        lines = [
            f"! Opt {method}",
            "%maxcore 2000",
            f"* xyz {int(charge)} {int(multiplicity)}",
        ]
        for element, x, y, z in frame.atoms:
            lines.append(f"{element:<2} {x:>16.10f} {y:>16.10f} {z:>16.10f}")
        lines.append("*")
        text = "\n".join(lines) + "\n"
        (target / name).write_text(text, encoding="utf-8")
        relative = (
            ensemble.relative_kcal[frame.index - 1]
            if frame.index - 1 < len(ensemble.relative_kcal)
            else float("nan")
        )
        inputs.append(
            UpgradeInput(
                name=name,
                index=frame.index,
                energy_Eh=frame.energy_Eh,
                relative_kcal=relative,
                text=text,
            )
        )
    return UpgradePlan(
        directory=target,
        method=method,
        charge=int(charge),
        multiplicity=int(multiplicity),
        inputs=tuple(inputs),
    )


def render_plan(plan: UpgradePlan) -> str:
    """The batch listing and the run guidance, as the menu prints them."""
    lines = [
        "CREST conformer upgrade (one ORCA optimisation input per conformer):",
        f"  ensemble conformers upgraded: {len(plan.inputs)}",
        f"  method line: ! Opt {plan.method}   charge {plan.charge}, multiplicity "
        f"{plan.multiplicity}",
        f"  directory: {plan.directory}",
        "",
        "  input                ensemble E (Eh)      relative (kcal/mol)",
    ]
    for item in plan.inputs:
        lines.append(
            f"  {item.name:<20} {item.energy_Eh:>18.8f} {item.relative_kcal:>12.3f}"
        )
    lines += [
        "",
        "Next steps:",
        f"  1. run one input after the other (or submit the directory as a batch): "
        f"orca {plan.inputs[0].name} > {plan.inputs[0].name[:-4]}.out",
        "  2. read each output back with menu 1 (the check-up report carries the final "
        "geometry and energy)",
        "",
        "Boundaries:",
        "  - the geometries are the pre-screen ensemble's, used as-is; the selection "
        "order is the ensemble's own (energy-ascending), so a count takes the lowest "
        "conformers first",
        "  - the upgrade re-ranks nothing by itself: comparing the upgraded energies "
        "and re-weighting the ensemble is registered as the next increment",
    ]
    return "\n".join(lines)


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the upgrade route."""
    return (
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "The generated input is accepted end to end, measured on 6.1.1 "
                "(2026-09-30): the first conformer of the n-butane fixture ensemble "
                "(the anti form, GFN2 geometry from the CREST quick run) written as "
                "``! Opt r2SCAN-3c`` converged in five geometry cycles "
                "(the HURRAY / optimization-converged marker) to "
                "-158.392455872472 Eh; the run is frozen as "
                "fixtures/crest/upgrade_acceptance.out."
            ),
            ref="tests/test_crest_upgrade.py; fixtures/crest/",
        ),
    )
