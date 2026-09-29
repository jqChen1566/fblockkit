"""5.4: XAS report from ROCIS core-excited spectra (menu 37).

Reads the absorption-spectrum blocks the ROCIS module prints (manual section
5.7; measured format notes in ``parsers/rocis_spectra.py`` and the fixtures)
and reports:

- the excitation-energy table summary;
- the transition-level table of one primary block (the SOC-corrected
  electric-dipole block when present -- its fosc column is population
  weighted as printed -- otherwise the plain electric-dipole block);
- the **branching ratio** of the spin-orbit-split edge: the transitions are
  split at their largest energy gap (the natural two-cluster default for an
  L3/L2 or M5/M4 pair) and the two clusters' centroids, integrated
  oscillator strengths and their ratio are reported.  The ratio is a data
  fact; comparing it with the statistical value (2:1 for a 2p edge, 3:2 for
  a 3d edge, from the (2j+1) degeneracies) is left to the caller's ``r_stat``
  parameter -- deviations indicate electrostatic and spin-orbit effects
  (Thole & van der Laan, Phys. Rev. B 38, 3158 (1988), whose rules this
  menu cites rather than reimplements);
- the RIXS bookkeeping: whether the run produced RIXS cross-section data
  files or hit the engine's measured refusal mode (zero intermediate/final
  states), with the engine's own TIP text and the ``orca_mapspc`` recipe
  from the manual (section 5.7.4.2).

Nothing here is a simulation: ROCIS's own accuracy statement (section 5.7.3:
several approximations, qualitatively correct results) is carried in every
report's boundaries.
"""

from __future__ import annotations

from typing import Any

from ..knowledge.models import EVIDENCE_LITERATURE, EVIDENCE_MEASURED, Evidence

__all__ = [
    "XasError",
    "primary_spectrum",
    "branching",
    "render",
    "evidence",
]

#: preference order of the spectrum blocks (population-weighted SOC first)
PRIMARY_PREFERENCE = (
    "soc_dipole_length",
    "soc_combined_length",
    "soc_combined_origin_independent",
    "dipole_length",
    "combined_length",
    "combined_origin_independent",
)

#: the default minimal gap (eV) that separates two edge clusters; the measured
#: L3/L2 separation is ~10-15 eV, M4/M5 ~4-7 eV (instrumental default, not a
#: literature constant)
DEFAULT_GAP_EV = 3.0

#: the RIXS processing recipe from manual section 5.7.4.2 (verbatim flags)
MAPSPC_RECIPE = (
    "orca_mapspc <file>.rixssoc RIXS -x0<lo> -x1<hi> -x2<lo> -x3<hi> "
    "-w<gamma_x> -g<gamma_y> -l -n<nx> -m<ny> -dx<shift> -eaxis1"
)


class XasError(Exception):
    """A refusal with the next step spelled out."""


def primary_spectrum(data: dict[str, Any]) -> tuple[str, list[dict[str, Any]]]:
    """The best available spectrum block (key, rows)."""
    spectra = data.get("spectra") or {}
    if not spectra:
        raise XasError(
            "no absorption-spectrum block found. Next step: this menu reads the "
            "output of an ORCA ROCIS calculation (%rocis with DoGenROCIS); check "
            "that the run reached the 'ROCIS-EXCITATION SPECTRA' section."
        )
    for key in PRIMARY_PREFERENCE:
        if key in spectra and spectra[key]:
            return key, spectra[key]
    key = next(iter(spectra))
    return key, spectra[key]


def branching(
    rows: list[dict[str, Any]], *, gap_min_ev: float = DEFAULT_GAP_EV
) -> dict[str, Any]:
    """Split the transitions at their largest energy gap and compare clusters.

    Returns ``{"split_ev": float | None, "clusters": [...]}`` with one cluster
    when no gap reaches ``gap_min_ev``.  Cluster entries carry n, the energy
    range, the fosc-weighted centroid (and the plain mean), and the summed
    oscillator strength; the ratio (low/high) is added when there are two.
    """
    ordered = sorted(rows, key=lambda row: row["ev"])
    gaps = [
        (ordered[index + 1]["ev"] - ordered[index]["ev"], index)
        for index in range(len(ordered) - 1)
    ]
    if not gaps:
        clusters = [ordered]
        split_ev = None
    else:
        largest, index = max(gaps, key=lambda item: item[0])
        if largest >= gap_min_ev:
            clusters = [ordered[: index + 1], ordered[index + 1 :]]
            split_ev = (ordered[index]["ev"] + ordered[index + 1]["ev"]) / 2.0
        else:
            clusters = [ordered]
            split_ev = None
    summary: list[dict[str, Any]] = []
    for cluster in clusters:
        strengths = [row["fosc"] for row in cluster]
        total = sum(strengths)
        if total > 0:
            centroid = sum(row["ev"] * row["fosc"] for row in cluster) / total
        else:
            centroid = sum(row["ev"] for row in cluster) / len(cluster)
        summary.append(
            {
                "n": len(cluster),
                "ev_min": cluster[0]["ev"],
                "ev_max": cluster[-1]["ev"],
                "centroid_ev": centroid,
                "sum_fosc": total,
            }
        )
    result: dict[str, Any] = {"split_ev": split_ev, "clusters": summary}
    if len(summary) == 2 and summary[1]["sum_fosc"] > 0:
        result["ratio_low_high"] = summary[0]["sum_fosc"] / summary[1]["sum_fosc"]
    return result


def render(
    data: dict[str, Any],
    *,
    source: str,
    r_stat: float | None = None,
    gap_min_ev: float = DEFAULT_GAP_EV,
) -> str:
    """The menu-37 report body."""
    lines: list[str] = []
    lines.append(f"Core-excited spectra report ({source})")
    table = data.get("excitation_table") or []
    if table:
        mults = sorted({row["mult"] for row in table})
        lines.append(
            f"  excitation table: {len(table)} roots, multiplicities {mults}; "
            f"lowest excitation {min(row['ev'] for row in table if row['ev'] > 0):.3f} eV"
        )
    key, rows = primary_spectrum(data)
    # the SOC-corrected tables carry the full state-pair matrix; all-zero rows
    # (and the print-precision dregs) are noise for the table and the split
    nonzero = [row for row in rows if row["fosc"] > 1e-6]
    lines.append("")
    lines.append(
        f"  primary block: {key} ({len(rows)} transitions, {len(nonzero)} with non-zero fosc)"
    )
    listed = nonzero if len(nonzero) <= 80 else nonzero[:80]
    lines.append("  i_root -> j_root   energy (eV)   fosc")
    for row in listed:
        lines.append(
            f"  {row['i_root']:>4} ({row['i_label']}) -> {row['j_root']:>3} ({row['j_label']})"
            f"   {row['ev']:>11.4f}   {row['fosc']:.6f}"
        )
    if len(nonzero) > len(listed):
        lines.append(f"  ... {len(nonzero) - len(listed)} further non-zero transitions omitted")
    lines.append("")
    edge = branching(nonzero if nonzero else rows, gap_min_ev=gap_min_ev)
    if edge["split_ev"] is None:
        lines.append(
            f"  branching: no edge split detected (largest gap below {gap_min_ev:.1f} eV);"
            f" all {edge['clusters'][0]['n']} transitions treated as one cluster"
            f" (sum fosc = {edge['clusters'][0]['sum_fosc']:.4f})"
        )
    else:
        low, high = edge["clusters"]
        lines.append(
            f"  branching at {edge['split_ev']:.2f} eV (largest gap):"
        )
        lines.append(
            f"    low-energy cluster : n = {low['n']:>3}  "
            f"{low['ev_min']:.2f}-{low['ev_max']:.2f} eV  centroid {low['centroid_ev']:.2f} eV"
            f"  sum fosc = {low['sum_fosc']:.4f}"
        )
        lines.append(
            f"    high-energy cluster: n = {high['n']:>3}  "
            f"{high['ev_min']:.2f}-{high['ev_max']:.2f} eV  centroid {high['centroid_ev']:.2f} eV"
            f"  sum fosc = {high['sum_fosc']:.4f}"
        )
        ratio = edge.get("ratio_low_high")
        if ratio is not None:
            lines.append(f"    ratio (low/high) = {ratio:.3f}")
            if r_stat is not None and r_stat > 0:
                lines.append(
                    f"    against the statistical value {r_stat:.3f}: ratio/stat = {ratio / r_stat:.3f}"
                    " (deviations indicate electrostatic and spin-orbit effects;"
                    " Thole & van der Laan, PRB 38, 3158 (1988))"
                )
    lines.append("")
    if data.get("riqs_refused"):
        lines.append(
            "  RIXS: the run identified the RIXS flag but produced zero intermediate/final"
            " states, so no cross sections were evaluated (the engine's refusal mode)."
            " The measured cause here is the orbital-window (OrbWin) content: the RIXS"
            " variant needs the 6-element window with the two donor spaces; widen the"
            " windows or increase NRoots as the engine's own TIP says."
        )
    elif data.get("riqs_intermediate") is not None or data.get("riqs_final") is not None:
        lines.append(
            f"  RIXS: cross sections evaluated (intermediate states"
            f" {data.get('riqs_intermediate')}, final states {data.get('riqs_final')});"
            " the 2D data files are prepared for processing with"
        )
        lines.append(f"    {MAPSPC_RECIPE}")
    else:
        lines.append(
            "  RIXS: not requested in this run (no RIXS bookkeeping found)."
            " To request it, add DoRIXS/DoRIXSSOC/DoElastic in the ROCIS rel block"
            " with a 6-element OrbWin (manual section 5.7.4)."
        )
    lines.append("")
    lines.append(
        "  Boundaries: ROCIS applies several approximations (manual section 5.7.3 says the"
        " results are qualitatively correct); the branching ratio is a data fact, and its"
        " statistical-value comparison is the caller's stated reference."
    )
    return "\n".join(lines)


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the report's criteria and the measured format facts."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "The branching-ratio semantics (nonstatistical ratios from initial-state"
                " spin-orbit splitting and core-hole/valence electrostatics; the ratio"
                " as a spin-state indicator) follow Thole & van der Laan; this menu"
                " reports the ratio and cites the reference rather than reimplementing"
                " its atomic rules."
            ),
            ref="Thole B. T., van der Laan G., Phys. Rev. B 1988, 38, 3158, DOI 10.1103/PhysRevB.38.3158",
            url="https://doi.org/10.1103/PhysRevB.38.3158",
            bibkey="thole1988branching",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "ORCA 6.1.1 fixtures rocis/ (2026-09-29): the ROCIS module prints up to"
                " seven absorption blocks (electric dipole, velocity dipole and five"
                " combined D2/M2/Q2 variants from DecomposeFosc); with DoSOC the"
                " SOC-corrected blocks weight fosc by the initial-state population; the"
                " RIXS variant needs the 6-element OrbWin (a 4-element window gives the"
                " engine's zero-intermediate-states refusal)."
            ),
            ref="ORCA 6.1.1, measured 2026-09-29 (fixtures rocis/)",
        ),
    )
