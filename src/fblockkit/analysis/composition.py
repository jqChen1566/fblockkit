"""A1: orbital composition analysis -- per-MO shell assignment (classified by the
largest weight, and partitioned by occupation first).

Criterion (this group's measured EuF lesson; both intuitive formulations fail):

- an absolute weight threshold is unusable: Eu's ANO basis set carries 49 f-type
  AOs, so 48 of the 76 orbitals in the full space have an f weight > 0.5 -- the
  threshold does not discriminate;
- the ratio of two shells is unusable: it produces 0/0-style false signals (one
  orbital gave 0.0039/0.0039 = 1.0);
- the right answer: assign each orbital to the **(atom, shell)** with the largest
  weight;
- and ranking/selection must be **partitioned by occupation** first: ranking the
  whole space by weight to pick the "most f-like" orbitals puts three virtual
  orbitals in the top three (measured on EuF) -- every ranking in this module is
  partitioned.

Data source: the ``LOEWDIN ORBITAL-COMPOSITIONS`` table of the ORCA output
(requires ``%output Print[P_ReducedOrbPopMO_L] 1``; printed by default at the
normal print level).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence

from ..knowledge.models import EVIDENCE_MEASURED, Evidence, ParseResult, ReportSection

OCCUPIED_THRESHOLD = 0.02  # same convention as the active space: occupation <= 0.02 counts as empty

_EVIDENCE_CRITERION = Evidence(
    kind=EVIDENCE_MEASURED,
    text=(
        "Measured on the molecular end of EuF: the absolute-weight threshold (49 f-type "
        "AOs give 48 of 76 orbitals an f weight > 0.5) and the ratio of two shells (the "
        "0/0-style false signal 0.0039/0.0039 = 1.0) both fail; the right answer is to "
        "take the shell with the largest weight. Ranking the whole space by weight puts "
        "three virtual orbitals in the top three -- orbital selection must be partitioned "
        "by occupation first."
    ),
    ref="This group's differentiable multireference stack records: composition-split-needs-two-shells / rank-within-partition-not-across (2026-09)",
)


@dataclass(frozen=True)
class OrbitalRow:
    """The composition assignment of one orbital."""

    index: int
    occupation: float
    energy: float | None
    shells: tuple[Mapping[str, Any], ...]  # already sorted by weight, descending

    @property
    def dominant(self) -> Mapping[str, Any] | None:
        return self.shells[0] if self.shells else None

    @property
    def partition(self) -> str:
        return "occupied" if self.occupation > OCCUPIED_THRESHOLD else "virtual"

    def label(self) -> str:
        dom = self.dominant
        if dom is None:
            return "(no composition data)"
        return f"{dom['element']}{dom['atom'] + 1} {dom['shell']}"

    def shell_weight(self, shell: str) -> float:
        return sum(item["weight"] for item in self.shells if item["shell"] == shell)


def accepts(result: ParseResult) -> bool:
    """Whether this analyser applies (the output must have a per-MO composition table)."""
    return bool(result.sections.get("orbital_composition", {}).get("present"))


def orbital_rows(result: ParseResult) -> tuple[OrbitalRow, ...]:
    section = result.sections.get("orbital_composition", {})
    return tuple(
        OrbitalRow(
            index=int(entry["index"]),
            occupation=float(entry["occupation"]),
            energy=entry.get("energy"),
            shells=tuple(entry.get("shells", ())),
        )
        for entry in section.get("orbitals", ())
    )


def shell_ranking(
    rows: Sequence[OrbitalRow],
    shell: str,
    *,
    partition: str | None = None,
    top_n: int | None = None,
) -> tuple[OrbitalRow, ...]:
    """Rank by the weight of the given shell; ``partition`` in {None, "occupied", "virtual"}."""
    selected = [row for row in rows if partition is None or row.partition == partition]
    ranked = sorted(selected, key=lambda row: row.shell_weight(shell), reverse=True)
    return tuple(ranked[:top_n] if top_n is not None else ranked)


def _row_line(row: OrbitalRow) -> str:
    dom = row.dominant
    if dom is None:
        return f"MO {row.index:>4}  occ {row.occupation:.4f}  (no composition data)"
    secondary = next(
        (item for item in row.shells[1:] if item["weight"] >= 1.0), None
    )
    tail = (
        f"; next: {secondary['element']}{secondary['atom'] + 1} {secondary['shell']} "
        f"{secondary['weight']:.1f}%"
        if secondary
        else ""
    )
    return (
        f"MO {row.index:>4}  occ {row.occupation:.4f}  "
        f"{dom['element']}{dom['atom'] + 1} {dom['shell']} {dom['weight']:.1f}%{tail}"
    )


def run(result: ParseResult, *, shell: str = "f", top_n: int = 8) -> ReportSection:
    """Build the A1 report section."""
    rows = orbital_rows(result)
    if not rows:
        raise ValueError(
            "this output has no per-MO composition table. Next step: add "
            "%output Print[P_ReducedOrbPopMO_L] 1 to the input (printed by default at "
            "the normal print level) and rerun."
        )
    occupied = [row for row in rows if row.partition == "occupied"]
    virtual = [row for row in rows if row.partition == "virtual"]
    active = [row for row in rows if OCCUPIED_THRESHOLD < row.occupation < 2.0 - OCCUPIED_THRESHOLD]

    lines = [
        f"Composition source: LOEWDIN ORBITAL-COMPOSITIONS table ({len(rows)} orbitals in total)",
        f"Occupation partition: occ > {OCCUPIED_THRESHOLD} counts as occupied ({len(occupied)}), "
        f"otherwise as a virtual orbital ({len(virtual)})",
        "",
        f"{shell} composition ranking (occupied partition, top {top_n}):",
    ]
    ranking = shell_ranking(rows, shell, partition="occupied", top_n=top_n)
    top_weight = ranking[0].shell_weight(shell) if ranking else 0.0
    if ranking and top_weight >= 0.5:
        for position, row in enumerate(ranking, start=1):
            lines.append(f"  {position:>2}. {_row_line(row)}")
    else:
        lines.append(
            f"  (no {shell}-composition orbital in the occupied partition: the highest "
            f"{shell} weight is only {top_weight:.1f}% -- if this system should have "
            "electrons in that shell, check the initial guess / active-orbital window, "
            "it may be the wrong branch of a multiple-solution problem)"
        )

    lines += ["", "Dominant shells of the active orbitals (0.02 < occ < 1.98):"]
    if active:
        for row in active:
            lines.append(f"  {_row_line(row)}")
    else:
        lines.append("  (no active orbital among the occupations of this output -- a pure single-reference calculation)")

    lines += [
        "",
        "Criterion and convention: the assignment is the (atom, shell) with the largest "
        "weight; every ranking is partitioned by occupation first (measured on this "
        "group's EuF: ranking the whole space by weight picks up virtual orbitals).",
    ]
    return ReportSection(title="A1 orbital composition (per MO, Loewdin)", body="\n".join(lines))


def evidence() -> tuple[Evidence, ...]:
    """This analyser's provenance (for the report's provenance summary)."""
    return (_EVIDENCE_CRITERION,)
