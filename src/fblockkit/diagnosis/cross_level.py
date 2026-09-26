"""Cross-level solution consistency check (multiple-solution settings: one molecule,
one batch of solutions, two theoretical levels).

Use: for the case where several solutions of the same molecule (different 5f
occupations / spin states / starting orbitals, say) have each been given an
energy at a cheap level (HF, say) and an expensive level (CCSD(T), say), check
whether the two levels agree on which solution is best. Taking the "lowest
solution" by the cheap level alone picks a solution that is not the best one
whenever the ordering changes at the expensive level.

Literature basis (the fixture `fixtures/literature/pucl3_s18.json` is that table):
Lu, J.-B. et al., "Norm-Conserving 5f-in-Core Pseudopotentials and Gaussian Basis Sets
Optimized for Tri- and Tetra-Valent Actinides", J. Chem. Theory Comput. 2025, 21,
170-182, DOI 10.1021/acs.jctc.4c01189 (SI Table S18). Among the 21 PuCl3 5f
occupations in that table, the lowest HF solution and the lowest CCSD(T) solution
are not the same one; the second-lowest HF solution falls 1.76 kcal/mol behind the
best at the CCSD(T) level.

Input record format (a JSON file or an in-memory list):

    [{"label": "No.1", "energies": {"HF": -1930.87692, "CCSD(T)": -1932.97951}}, ...]
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from ..knowledge.models import EVIDENCE_LITERATURE, Evidence, Finding

HARTREE_TO_KCAL = 627.509474

_PUCL3_EVIDENCE = Evidence(
    kind=EVIDENCE_LITERATURE,
    text=(
        "The 21 5f-occupation solutions of PuCl3: the lowest HF solution (No.3) and the "
        "lowest CCSD(T) solution (No.1) are not the same one; the second-lowest HF "
        "solution (No.10) falls 0.00281 Ha (about 1.76 kcal/mol) behind the best at the "
        "CCSD(T) level. Conclusion: a solution must not be picked by the low-level "
        "ordering alone."
    ),
    ref="Lu J.-B., Zhang Y.-Y., Liu J.-B., Li J., J. Chem. Theory Comput., 2025, 21(1), 170-182, DOI 10.1021/acs.jctc.4c01189 (SI Table S18)",
    bibkey="lu2025normconserving",
    url="https://doi.org/10.1021/acs.jctc.4c01189",
)


class CrossLevelError(ValueError):
    """Invalid input for the cross-level check."""


def load_records(path: str | Path) -> list[dict[str, Any]]:
    """Read a records JSON file; an invalid structure is an error (never skipped silently).

    Two top-level forms are accepted: a list of records; or an object with a
    ``records`` key (which may carry a ``meta`` bibliography alongside).
    """
    p = Path(path)
    if not p.exists():
        raise CrossLevelError(f"record file does not exist: {p}. Next step: check the path.")
    data = json.loads(p.read_text(encoding="utf-8"))
    records = data.get("records") if isinstance(data, Mapping) else data
    if not isinstance(records, list) or not records:
        raise CrossLevelError(
            f"{p}: the top level must be a non-empty list of records, or an object with "
            f"a records list (each item carrying label and energies)."
        )
    return records


def _ranking(records: Sequence[Mapping[str, Any]], level: str) -> list[tuple[str, float]]:
    ranked = []
    for index, record in enumerate(records):
        label = record.get("label")
        energies = record.get("energies", {})
        if not isinstance(energies, Mapping) or level not in energies:
            raise CrossLevelError(
                f"record [{index}] has no {level!r} energy. Next step: check that this "
                f"level was computed for every record."
            )
        energy = energies[level]
        if not isinstance(energy, (int, float)):
            raise CrossLevelError(
                f"the {level!r} energy of record [{index}] is not a number: {energy!r}"
            )
        ranked.append((str(label), float(energy)))
    ranked.sort(key=lambda item: item[1])
    return ranked


def cross_level_check(
    records: Sequence[Mapping[str, Any]],
    levels: Sequence[str] = ("HF", "CCSD(T)"),
    top_n: int = 3,
    penalty_tol_kcal: float = 0.5,
) -> tuple[Finding, ...]:
    """Check that the two levels order the solutions the same way; return the
    conclusions (possibly empty).

    Trigger conditions (any one):
    - the two levels do not have the same lowest record;
    - a record in the low level's top ``top_n`` drops out of the high level's top ``top_n``;
    - the low level's best solution is penalised by more than ``penalty_tol_kcal``
      relative to the high level's best solution.
    """
    if len(levels) < 2:
        raise CrossLevelError("at least two levels are needed (HF and CCSD(T), say).")
    low_level, high_level = levels[0], levels[-1]
    low_rank = _ranking(records, low_level)
    high_rank = _ranking(records, high_level)
    high_position = {label: index for index, (label, _) in enumerate(high_rank)}
    high_energy = dict(high_rank)
    best_low_label = low_rank[0][0]
    best_high_label = high_rank[0][0]

    findings: list[Finding] = []
    penalty = high_energy[best_low_label] - high_energy[best_high_label]
    penalty_kcal = penalty * HARTREE_TO_KCAL
    if best_low_label != best_high_label or penalty_kcal > penalty_tol_kcal:
        findings.append(
            Finding(
                severity="warn",
                message=(
                    f"the two levels disagree on the best solution: {low_level} favours "
                    f"{best_low_label} while {high_level} favours {best_high_label}; "
                    f"picking the best by {low_level} leaves it {penalty_kcal:.2f} kcal/mol "
                    f"behind at the {high_level} level."
                ),
                evidence=(_PUCL3_EVIDENCE,),
                suggested_fix=(
                    f"Re-rank the candidate solutions at the {high_level} level before "
                    f"choosing one; do not use the {low_level} ordering directly as the "
                    "basis for solution selection."
                ),
                rule_id="XL-CROSS-LEVEL-BEST",
            )
        )
    dropped = [
        label
        for label, _ in low_rank[:top_n]
        if high_position.get(label, len(high_rank)) >= top_n
    ]
    if dropped:
        details = "; ".join(
            f"{label} ({low_level} position {index + 1} -> {high_level} position {high_position[label] + 1})"
            for index, (label, _) in enumerate(low_rank[:top_n])
            if label in dropped
        )
        findings.append(
            Finding(
                severity="warn",
                message=(
                    f"a solution in the low level's top {top_n} drops out of the top "
                    f"{top_n} at the higher level: {details}."
                ),
                evidence=(_PUCL3_EVIDENCE,),
                suggested_fix=(
                    "Re-check the energy, occupation and configuration assignment of the "
                    "dropped solution at the higher level (it may be a different "
                    "electronic state)."
                ),
                rule_id="XL-CROSS-LEVEL-TOP-N",
            )
        )
    return tuple(findings)
