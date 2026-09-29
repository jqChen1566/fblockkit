"""6.1: the CREST conformer-ensemble report (menu 43).

Renders the measured ensemble (frames, relative energies, the best
conformer) into the sorted table with Boltzmann weights at a chosen
temperature, plus the reading notes.  Weight convention: ``w_i ~
exp(-dE_i / (R T))`` with dE in kcal/mol and ``R = 1.9872042586e-3
kcal mol-1 K-1`` (CODATA-derived, the same constant the manual's worked
conversions use); the weights are normalised over the ensemble.
"""

from __future__ import annotations

import math

from ..knowledge.models import EVIDENCE_LITERATURE, EVIDENCE_MEASURED, Evidence
from ..parsers.crest import CrestEnsemble

__all__ = ["boltzmann_weights", "render", "evidence"]

_R_KCAL_PER_MOL_K = 1.98720425864083e-3  # CODATA molar gas constant in kcal/mol/K


def boltzmann_weights(
    ensemble: CrestEnsemble, *, temperature_K: float = 298.15
) -> tuple[float, ...]:
    """Normalised Boltzmann weights over the ensemble at the given temperature."""
    if temperature_K <= 0:
        raise ValueError("temperature must be positive")
    factors = [
        math.exp(-energy / (_R_KCAL_PER_MOL_K * temperature_K))
        for energy in ensemble.relative_kcal
    ]
    total = sum(factors)
    return tuple(factor / total for factor in factors)


def render(
    ensemble: CrestEnsemble,
    *,
    source: str,
    temperature_K: float = 298.15,
) -> str:
    """The report body (plain text, English, no timestamps of the day)."""
    weights = boltzmann_weights(ensemble, temperature_K=temperature_K)
    count = len(ensemble.frames)
    lines: list[str] = [
        "CREST conformer ensemble report",
        f"Source: {source}",
        "",
        f"Conformers: {count} (ensemble order; the energies are the frame comments)",
    ]
    if ensemble.best_index is not None:
        lines.append(f"Best conformer (crest_best.xyz): index {ensemble.best_index}")
    lines.append("")
    lines.append(f" idx   E (Eh)               rel (kcal/mol)   weight ({temperature_K:.2f} K)")
    for frame, relative, weight in zip(ensemble.frames, ensemble.relative_kcal, weights):
        lines.append(
            f" {frame.index:>4}  {frame.energy_Eh:.12f}      {relative:>8.3f}       {weight:>8.4f}"
        )
    lines.append("")
    lines.append(
        f"Cross-check: the frame energies reproduce crest.energies within "
        f"{ensemble.max_rel_deviation_kcal:.4f} kcal/mol (print precision)."
    )
    lines.append("")
    lines.append(
        "Reading notes: the ensemble energies are GFN2-xTB level (Tier 1) -- the "
        "ordering and the weights select candidates, they do not rank final "
        "stabilities; for a quantitative comparison re-optimise the leading "
        "conformers at the target level (ORCA). The ensemble was deduplicated by "
        "CREST's own RMSD/energy criteria (CREGEN), so indices are post-filter "
        "positions, not sampling order. The Boltzmann weights use R = 1.9872042586e-3 "
        "kcal mol-1 K-1 over the listed set only; a run with several distinct basins "
        "at a close energy has its weights split across the mirror images if the "
        "engine did not merge them."
    )
    return "\n".join(lines)


def evidence() -> tuple[Evidence, ...]:
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "The CREST program paper: the sampling/deduplication workflow whose "
                "ensemble this menu reads."
            ),
            ref="Pracht P., Bohle F., Grimme S., Phys. Chem. Chem. Phys., 2020, 22, 7169-7192, DOI 10.1039/C9CP06869D",
            bibkey="pracht2020crest",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "Fixture anchors measured on CREST 3.0.2 (101, 2026-09-30, n-butane "
                "GFN2): crest_conformers.xyz with leading-space atom-count lines and "
                "hartree frame comments, crest.energies (index + kcal/mol, three "
                "decimals), crest_best.xyz, and the topology-change abort of an "
                "unpreoptimised input (the capture offers options A/B/C)."
            ),
            ref="fixtures/crest/ (README records the run commands)",
        ),
    )
