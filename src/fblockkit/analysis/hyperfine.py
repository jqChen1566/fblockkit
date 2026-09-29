"""5.7: the hyperfine / EFG report from an ORCA EPRNMR run (menu 40).

Reads the "ELECTRIC AND MAGNETIC HYPERFINE STRUCTURE" section (manual
section 7.51.3): per nucleus the A-tensor components (Fermi contact, spin
dipole; second-order SOC where requested), the electric field gradient with
its electron/nuclear decomposition, and the density at the nucleus.  The
report computes two derived quantities from the printed EFG principal
values:

- the largest-|.| principal value ``Vzz`` and the asymmetry parameter
  ``eta = (Vxx - Vyy) / Vzz`` with |Vxx| <= |Vyy| <= |Vzz| (the standard NQR
  ordering);
- when the caller passes the nuclear quadrupole moment Q (barn), the
  nuclear quadrupole coupling constant ``eQVzz/h`` (MHz) and the
  first-order quadrupole splitting ``DeltaE_Q = eQVzz/2 * sqrt(1 + eta^2/3)``
  (the classic Mossbauer relation; first order only, I >= 3/2).  The
  conversion constant is derived from CODATA constants, not hard-coded.

Boundaries kept in the report: the A tensor in MHz requires the per-nucleus
nuclear parameters (I, P) in the input -- without them the engine prints
zeros (measured); the density at the nucleus RHO(0) is basis-domain
dependent (all-electron vs small-core ECP: measured 727216.9 vs 0.0756 on
the same geometry) and must not be compared across domains; on CASSCF
reference the engine computes the EFG / RHO(0) but switches the A-tensor
components off (measured: "CASSCF/ALL STATES AVERAGE", HFC iso/dip = NO).
"""

from __future__ import annotations

import math
from typing import Any

from ..knowledge.models import EVIDENCE_LITERATURE, EVIDENCE_MEASURED, Evidence

__all__ = ["HyperfineError", "efg_metrics", "nqcc_mhz", "delta_e_q_mhz", "render", "evidence"]

# CODATA-derived conversion: e Q Vzz / h with Q in barn and Vzz in a.u.
# 1 a.u. of EFG = E_h / (e * a0^2) = 9.717362e21 V/m^2;  1 barn = 1e-28 m^2.
_E = 1.602176634e-19  # C
_H = 6.62607015e-34  # J s
_EFG_AU_SI = 9.717362e21  # V m^-2 per a.u.
_NQCC_PER_BARN_AU = _E * 1e-28 * _EFG_AU_SI / _H  # Hz per (barn * a.u.)


def nqcc_mhz(Q_barn: float, Vzz_au: float) -> float:
    """The nuclear quadrupole coupling constant eQVzz/h in MHz."""
    return _NQCC_PER_BARN_AU * Q_barn * Vzz_au / 1e6


def delta_e_q_mhz(Q_barn: float, Vzz_au: float, eta: float) -> float:
    """First-order quadrupole splitting |eQVzz|/2 * sqrt(1 + eta^2/3), MHz."""
    return abs(nqcc_mhz(Q_barn, Vzz_au)) / 2.0 * math.sqrt(1.0 + (eta * eta) / 3.0)


def efg_metrics(nucleus: dict[str, Any]) -> dict[str, Any] | None:
    """Vzz, eta and the ordered principal values from the V(Tot) row."""
    values = nucleus.get("V_Tot")
    if not values or len(values) != 3:
        return None
    by_abs = sorted(values, key=abs, reverse=True)
    Vzz = by_abs[0]
    rest = sorted(by_abs[1:], key=abs)
    Vxx, Vyy = rest[0], rest[1]
    eta = (Vxx - Vyy) / Vzz
    return {"Vzz": Vzz, "eta": eta, "ordered": [Vxx, Vyy, Vzz]}


class HyperfineError(Exception):
    """A refusal with the next step spelled out."""


def _fmt_vec(values: list[float] | None) -> str:
    if not values:
        return "n/a"
    return "(" + ", ".join(f"{value:.6f}" for value in values) + ")"


def render(data: dict[str, Any], *, source: str, nuclear_Q_barn: float | None = None) -> str:
    """The menu-40 report body."""
    if not data.get("present"):
        raise HyperfineError(
            "no hyperfine / EFG section found. Next step: this menu reads an "
            "ORCA output with the EPRNMR properties requested (%eprnmr with a "
            "Nuclei list, e.g. 'Nuclei = all Ce { aiso, adip, fgrad, rho }'); "
            "note the coordinates block must precede %eprnmr when Nuclei is "
            "given (measured)."
        )
    lines: list[str] = [f"Hyperfine / EFG report ({source})"]
    header = data.get("header") or {}
    if header:
        counts = "  ".join(
            f"{key} {value}"
            for key, value in header.items()
            if key not in ("method", "multiplicity")
        )
        lines.append(
            f"  method: {header.get('method', '?')}"
            + (f"; multiplicity {header['multiplicity']}" if "multiplicity" in header else "")
            + (f"; nuclei computed: {counts}" if counts else "")
        )
    for nucleus in data.get("nuclei") or []:
        lines.append("")
        lines.append(f"  nucleus {nucleus['element']} (engine index {nucleus['index']}):")
        a = nucleus.get("a") or {}
        q = nucleus.get("q") or {}
        lines.append(
            f"    nuclear parameters: I = {a.get('I', '?')}, "
            f"P = {a.get('P_MHz_au3', '?')} MHz/au^3, Q = {q.get('Q_barn', '?')} barn"
        )
        a_zero = (
            nucleus.get("A_Tot") is not None
            and all(abs(value) < 1e-12 for value in nucleus["A_Tot"])
            and (a.get("P_MHz_au3") or 0.0) == 0.0
        )
        if a_zero:
            lines.append(
                "    A tensor (MHz): all zeros -- nuclear parameters (I, P) were"
                " not supplied in the input (measured engine behaviour)"
            )
        elif nucleus.get("A_Tot") is not None:
            lines.append(
                f"    A tensor (MHz): A(iso) = {nucleus['A_iso']:.4f}  "
                f"A(Tot) = {_fmt_vec(nucleus['A_Tot'])}"
            )
            if nucleus.get("A_FC") is not None:
                lines.append(f"      A(FC) = {_fmt_vec(nucleus['A_FC'])}")
            if nucleus.get("A_SD") is not None:
                lines.append(f"      A(SD) = {_fmt_vec(nucleus['A_SD'])}")
        else:
            flags = nucleus.get("hfc_flags") or {}
            if flags and not any(flags.values()):
                lines.append(
                    "    A tensor: not computed by the engine for this reference"
                    " (measured: CASSCF switches the A components off; DFT with"
                    " nuclear parameters P/I prints them)"
                )
        metrics = efg_metrics(nucleus)
        if metrics is not None:
            eta_text = f"{metrics['eta']:.4f}"
            engine_eta = nucleus.get("engine_eta")
            if engine_eta is not None:
                eta_text += f" (engine {engine_eta:.4f})"
            lines.append(
                f"    EFG: V(Tot) = {_fmt_vec(nucleus['V_Tot'])} a.u.; "
                f"|Vzz| = {abs(metrics['Vzz']):.6f}, eta = {eta_text}"
            )
            if nucleus.get("V_el") is not None and nucleus.get("V_nuc") is not None:
                summed = [
                    el + nuc
                    for el, nuc in zip(nucleus["V_el"], nucleus["V_nuc"])
                ]
                deviation = max(
                    abs(a - b) for a, b in zip(summed, nucleus["V_Tot"])
                )
                lines.append(
                    f"      V(El) = {_fmt_vec(nucleus['V_el'])}  "
                    f"V(Nuc) = {_fmt_vec(nucleus['V_nuc'])}  "
                    f"(sum vs V(Tot): max deviation {deviation:.2e})"
                )
            engine_nqcc = nucleus.get("e2qQ_MHz")
            if engine_nqcc is not None and abs(engine_nqcc) > 1e-9:
                lines.append(
                    f"      engine quadrupole block: e**2qQ = {engine_nqcc:.4f} MHz "
                    f"(Q = {(nucleus.get('quad_tensor') or {}).get('Q_barn', '?')} barn, "
                    f"I = {(nucleus.get('quad_tensor') or {}).get('I', '?')})"
                )
            if nuclear_Q_barn:
                nqcc = nqcc_mhz(nuclear_Q_barn, metrics["Vzz"])
                deq = delta_e_q_mhz(nuclear_Q_barn, metrics["Vzz"], metrics["eta"])
                lines.append(
                    f"      with Q = {nuclear_Q_barn} barn: eQVzz/h = {nqcc:.4f} MHz; "
                    f"first-order DeltaE_Q = {deq:.4f} MHz (MHz == NQR nu_Q)"
                )
        if nucleus.get("rho0") is not None:
            lines.append(
                f"    Rho(0) = {nucleus['rho0']:.6g} a.u.**-3"
                "  [basis-domain dependent: all-electron carries the core,"
                " small-core ECP reports the valence density only]"
            )
    lines.append("")
    lines.append(
        "  Reading notes: eta = (Vxx - Vyy)/Vzz with |Vxx| <= |Vyy| <= |Vzz|; "
        "the quadrupole splitting is first order (I >= 3/2, higher-order terms "
        "neglected). The A tensor values are the engine's 'SAI' spin-Hamiltonian "
        "convention; their MHz scale needs the per-nucleus parameters (I, P) in "
        "the input -- zeros mean they were not supplied (measured). Rho(0) must "
        "not be compared across basis domains (all-electron vs ECP)."
    )
    return "\n".join(lines)


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the definitions and the measured format facts."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "The electric field gradient enters hyperfine observables through"
                " the nuclear quadrupole coupling constant eQVzz/h (the NQR"
                " definition and the nuclear quadrupole moment it references are"
                " the Aerts & Brown 2019 framework); the ionic-crystal EFG /"
                " quadrupole-coupling relation is Bersohn 1958. The first-order"
                " Mossbauer quadrupole splitting follows as eQVzz/2 *"
                " sqrt(1 + eta^2/3) with the standard NQR ordering of the"
                " principal values."
            ),
            ref=(
                "Aerts A., Brown A., J. Chem. Phys. 2019, 150, 224302, "
                "DOI 10.1063/1.5097151; Bersohn R., J. Chem. Phys. 1958, 29, "
                "326-333, DOI 10.1063/1.1744480"
            ),
            url="https://doi.org/10.1063/1.5097151",
            bibkey="aerts2019nqcc",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "ORCA 6.1.1 fixture hyperfine/ (2026-09-29): the EPRNMR section"
                " prints per nucleus the A components, the raw EFG matrix, the"
                " V(Tot) principal-value row with orientation, the V(El)/V(Nuc)"
                " decomposition and RHO(0); the coordinates block must precede"
                " %eprnmr when Nuclei is given; without nuclear parameters the A"
                " matrix is zero; on CASSCF the A components are off while"
                " EFG/RHO(0) are computed; V(Nuc) is identical for the DFT and"
                " CASSCF probes of the same geometry (cross-check), and"
                " V(El)+V(Nuc) reproduces V(Tot) to the printed digits."
            ),
            ref="ORCA 6.1.1, measured 2026-09-29 (fixtures hyperfine/)",
        ),
    )
