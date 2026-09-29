"""5.6: the polynuclear-magnetism workflow from a POLY_ANISO output (menu 39).

POLY_ANISO (the ``otool_poly_aniso`` program shipped with ORCA, v1.0.0,
Ungur & Chibotaru) combines the single-ion ab initio data of each magnetic
center (the SINGLE_ANISO ``aniso_N.input`` files) with the inter-site
exchange and dipole-dipole coupling into cluster magnetic properties.  Two
boundaries are stated in the report and must stay attached to it:

- the exchange constants J are the caller's input (measured elsewhere or
  fitted); this workflow never computes or fits them.  A joint ab initio
  computation of exchange splittings (e.g. the LDF-CAHF / many-state
  PNO-CASPT2 route of Hallmen et al., PCCP 2019) is outside the ORCA
  ecosystem and is registered as a documented termination, not implemented;
- the Lines model (H = -sum J_p s_i.s_j) is exact only for two isotropic
  spins, one Ising plus one isotropic spin, or two Ising spins, and
  approximate otherwise; the dipole-dipole coupling is evaluated exactly
  from the ab initio moments and usually dominates in strongly anisotropic
  lanthanides (ORCA manual section 7.18).

The method reference is Ungur & Chibotaru, Chem. Eur. J. 2017 (the ab
initio crystal-field framework behind the SINGLE_ANISO/POLY_ANISO tools).
"""

from __future__ import annotations

from typing import Any

from ..knowledge.models import EVIDENCE_LITERATURE, EVIDENCE_MEASURED, Evidence

__all__ = ["PolyAnisoError", "render", "evidence"]


class PolyAnisoError(Exception):
    """A refusal with the next step spelled out."""


def render(data: dict[str, Any], *, source: str) -> str:
    """The menu-39 report body."""
    if not data.get("present"):
        raise PolyAnisoError(
            "no POLY_ANISO output found. Next step: this menu reads the output "
            "of the ORCA polynuclear-magnetism driver ($ORCA/otool_poly_aniso < "
            "poly_aniso.input > poly_aniso.output); the input must start with "
            "&POLY_ANISO and the per-center SINGLE_ANISO data files must be "
            "named aniso_1.input, aniso_2.input, ..."
        )
    lines: list[str] = [f"Polynuclear magnetism (POLY_ANISO) report ({source})"]
    counts = data.get("counts") or {}
    if counts:
        lines.append(
            f"  centers: {counts.get('independent', '?')} independent, "
            f"{counts.get('all', '?')} in total"
        )
    for center in data.get("centers") or []:
        parts = [f"data file {center.get('data_file', '?')}"]
        if center.get("coords"):
            parts.append(
                "at (" + ", ".join(f"{value:.3f}" for value in center["coords"]) + ") A"
            )
        if center.get("n_so_states") is not None:
            parts.append(f"{center['n_so_states']} spin-orbit states")
        g = center.get("g_tensor") or {}
        if g:
            parts.append(
                "g = " + " / ".join(f"{g[axis]['value']:.4f}" for axis in "XYZ" if axis in g)
            )
        lines.append(f"    center {center.get('type', '?')}: " + "; ".join(parts))
        spectrum = center.get("so_spectrum_cm1") or []
        if spectrum:
            if len(spectrum) > 8:
                shown = " ".join(f"{value:.2f}" for value in spectrum[:8]) + f" ... ({len(spectrum)} values)"
            else:
                shown = " ".join(f"{value:.2f}" for value in spectrum)
            lines.append(f"      spin-orbit spectrum (cm-1): {shown}")
    exchange = data.get("exchange") or {}
    if exchange:
        pieces = []
        if "coupled_states" in exchange:
            pieces.append(f"{exchange['coupled_states']} coupled states")
        if "pairs_total" in exchange:
            pieces.append(f"{exchange['pairs_total']} pairs")
        if "lines1" in exchange:
            pieces.append(f"Lines-1 {exchange['lines1']}")
        if "dipole_dipole" in exchange:
            pieces.append(f"dipole-dipole {exchange['dipole_dipole']}")
        if "ito" in exchange:
            pieces.append(f"ITO decomposition {exchange['ito']}")
        lines.append("")
        lines.append("  exchange: " + "; ".join(pieces))
        for pair in exchange.get("pairs_j") or []:
            first, second = pair["centers"]
            lines.append(
                f"    pair {pair['pair']} (centers {first}-{second}): "
                f"J = {pair['J_cm1']:.5f} cm-1"
            )
    interaction = data.get("interaction") or []
    if interaction:
        lines.append("")
        lines.append("  first-order anisotropic coupling (weights of the decomposition):")
        grouped: dict[Any, dict[str, Any]] = {}
        order: list[Any] = []
        for entry in interaction:
            key = (entry.get("pair"), entry.get("model"))
            if key not in grouped:
                grouped[key] = {"terms": [], "full": None}
                order.append(key)
            if (entry.get("term") or "").startswith("Full"):
                grouped[key]["full"] = entry.get("matrix")
            else:
                grouped[key]["terms"].append((entry.get("term"), entry.get("weight_percent")))
        for key in order:
            pair, model = key
            group = grouped[key]
            weight_text = " / ".join(
                f"{term.split()[0]} {weight:.3f}%"
                for term, weight in group["terms"]
                if weight is not None
            )
            lines.append(f"    pair {pair[0]}-{pair[1]} [{model}]: {weight_text}")
            for row in group["full"] or []:
                lines.append("      " + "  ".join(f"{value: .6e}" for value in row))
    coupled = data.get("coupled_states") or []
    if coupled:
        lines.append("")
        lines.append("  coupled states (cm-1):")
        for state in coupled[:16]:
            lines.append(
                f"    state {state['state']}: lines {state['lines_cm1']:.6f}  "
                f"dipole {state['dipole_cm1']:.6f}  total {state['total_cm1']:.6f}  "
                f"relative {state['relative_cm1']:.6f}"
            )
        if len(coupled) > 16:
            lines.append(f"    ... ({len(coupled)} states in total)")
    chit = data.get("chiT") or []
    if chit:
        lines.append("")
        lines.append(
            f"  chiT(T): {len(chit)} points, T = {chit[0]['T']:.4g} to "
            f"{chit[-1]['T']:.4g} K; chiT = {chit[0]['chiT']:.6f} -> "
            f"{chit[-1]['chiT']:.6f} cm3 K mol-1"
        )
    vv = data.get("van_vleck") or []
    if vv:
        first, last = vv[0], vv[-1]
        lines.append(
            f"  Van Vleck susceptibility tensors at {len(vv)} temperatures; "
            f"main values at {first['T']:.4g} K: "
            + " / ".join(f"{value:.6f}" for value in first["main_values"])
            + f"; at {last['T']:.4g} K: "
            + " / ".join(f"{value:.6f}" for value in last["main_values"])
        )
    for note in data.get("notes") or []:
        lines.append(f"  note: {note}")
    lines.append("")
    lines.append(
        "  Reading notes: the exchange constants are the caller's input "
        "(measured elsewhere or fitted) -- this workflow never computes or "
        "fits them; a joint ab initio computation of exchange splittings "
        "(the LDF-CAHF / many-state PNO-CASPT2 route) is outside the ORCA "
        "ecosystem and is registered as a documented termination, not "
        "implemented here. The Lines model is exact only for two isotropic "
        "spins, one Ising plus one isotropic spin, or two Ising spins, and "
        "approximate otherwise; the dipole-dipole coupling is evaluated "
        "exactly from the ab initio moments and usually dominates in strongly "
        "anisotropic lanthanides (ORCA manual section 7.18; method reference: "
        "Ungur & Chibotaru, Chem. Eur. J. 2017)."
    )
    return "\n".join(lines)


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the method reference and the measured format facts."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "The polynuclear magnetic framework (non-perturbative spin-orbit"
                " single-ion data combined with Lines-model exchange and exact"
                " dipole-dipole coupling) is the ab initio crystal-field approach"
                " of Ungur & Chibotaru 2017; in the POLY_ANISO workflow the"
                " exchange constants are user input or fitted parameters, never"
                " computed by the program (ORCA manual section 7.18)."
            ),
            ref=(
                "Ungur L., Chibotaru L. F., Chem. Eur. J. 2017, 23, 3708-3718, "
                "DOI 10.1002/chem.201605102"
            ),
            url="https://doi.org/10.1002/chem.201605102",
            bibkey="ungur2017abinitio",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "ORCA 6.1.1 fixture poly_aniso/ (2026-09-29): otool_poly_aniso"
                " v1.0.0 (compiled 2025-12-04) runs independently of ORCA; the"
                " input must start with &POLY_ANISO (otherwise the read aborts"
                " in fetch_init); the per-center data files must be named"
                " aniso_1.input, aniso_2.input, ... (a renamed file aborts in"
                " inquire_key_presence); the two-center probe returns rc=0 with"
                " the per-center echo (g values matching the menu-36 SINGLE_ANISO"
                " run of the same data), the interaction-matrix decomposition"
                " (Isotropic / Symmetric / Anti-Symmetric weights), the coupled"
                " states, the chiT table (101 points) and the Van Vleck tensor"
                " sequence (101 temperatures)."
            ),
            ref="ORCA 6.1.1, otool_poly_aniso v1.0.0, measured 2026-09-29 (fixtures poly_aniso/)",
        ),
    )
