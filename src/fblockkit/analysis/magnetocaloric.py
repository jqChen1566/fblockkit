r"""5.7: magnetic entropy and the magnetocaloric effect (menu 41).

Two data routes, both pure post-processing:

- **levels route** (a SINGLE_ANISO output, or the per-center echo of a
  POLY_ANISO run): the spin-orbit level spectrum (cm^-1) gives the magnetic
  entropy through the canonical ensemble,
  ``S(T) = R (ln Z + <E> / kT)`` with ``Z = sum_k g_k exp(-E_k/kT)`` (the
  Gibbs-function route of Szalowski & Kowalewska 2020, Eq. (4)-(5)).  The
  high-temperature check ``S(T -> inf) = R ln(sum g)`` is printed -- but the
  engine's level list is truncated (the lowest n_soc states only), so the
  limit bounds the temperature range rather than validating it: the report
  says so and the printed S must not be trusted above the temperature where
  the truncated sum saturates;
- **magnetization route** (a POLY_ANISO run with HINT/TMAG): the
  powder-averaged molar magnetization table M(H, T) (Bohr magnetons) gives
  the isothermal entropy change through the Maxwell relation
  ``(dS/dH)_T = (dM/dT)_H``:
  ``DeltaS(T, H) = N_A mu_B \int_0^H (dM/dT)_{H'} dH'``, evaluated with
  adjacent-temperature differences and a trapezoid in H (converged data
  required; the table starts at H = 0.0001 T).

Conventions: ``DeltaS = S(T, H) - S(T, 0)``; the report prints ``-DeltaS``
(direct magnetocaloric effect when positive -- the cooling capacity quoted
in the literature; inverse MCE appears as negative values).  The method
reference is Szalowski & Kowalewska, Materials 2020 (conventions follow
their Eq. (11) ``DeltaS_T = S(T, Bi) - S(T, Bf)``).
"""

from __future__ import annotations

import math
from typing import Any

from ..knowledge.models import EVIDENCE_LITERATURE, EVIDENCE_MEASURED, Evidence
from .plot_csv import plot_number

__all__ = [
    "MagnetocaloricError",
    "entropy_from_levels",
    "maxwell_delta_s",
    "levels_plot_csv",
    "maxwell_plot_csv",
    "render_levels",
    "render_maxwell",
    "evidence",
]

_R_J_MOL_K = 8.314462618  # J mol^-1 K^-1
_N_A_MU_B = 5.584937  # N_A * mu_B, J T^-1 mol^-1 (molar moment unit)
_K_B_CM1_PER_K = 0.695034800  # cm^-1 per K; 1/1.438776877 (as in analysis/relaxation.py)
#: the report's default temperature grid (K) of the levels route; the CSV
#: companion carries the same grid
_DEFAULT_TEMPERATURES = (1, 2, 5, 10, 20, 50, 100, 200, 300)


class MagnetocaloricError(Exception):
    """A refusal with the next step spelled out."""


def entropy_from_levels(
    energies_cm1: list[float],
    temperatures_K: list[float],
    degeneracies: list[float] | None = None,
) -> list[float]:
    """Magnetic entropy S(T) in J mol^-1 K^-1 from a level list (B = 0)."""
    if not energies_cm1:
        raise MagnetocaloricError("no level energies given")
    weights = degeneracies or [1.0] * len(energies_cm1)
    z_lin = [energy * _K_B_CM1_PER_K for energy in energies_cm1]
    result: list[float] = []
    for temperature in temperatures_K:
        z = 0.0
        e_avg = 0.0
        for theta, weight in zip(z_lin, weights):
            boltzmann = weight * math.exp(-theta / temperature)
            z += boltzmann
            e_avg += theta * boltzmann
        e_avg /= z
        result.append(_R_J_MOL_K * (math.log(z) + e_avg / temperature))
    return result


def maxwell_delta_s(
    temperatures_K: list[float],
    fields_T: list[float],
    M_muB: list[list[float]],
) -> list[dict[str, Any]]:
    """DeltaS(T, H) rows via the Maxwell relation, one row per T pair.

    ``M_muB`` is the aligned magnetisation matrix (rows: fields, columns:
    temperatures).  Returns rows ``{"T_K": mid temperature, "delta_S": [per
    field]}`` in J mol^-1 K^-1 (Sigma S(T, H) - S(T, 0)); the trapezoid adds
    the physical anchor (H = 0, dM/dT = 0) in front of the table.
    """
    if len(temperatures_K) < 2:
        raise MagnetocaloricError("the Maxwell route needs at least two temperatures")
    rows: list[dict[str, Any]] = []
    for i in range(len(temperatures_K) - 1):
        t_low, t_high = temperatures_K[i], temperatures_K[i + 1]
        dT = t_high - t_low
        if dT == 0:
            continue
        dmdT = [
            (M_muB[j][i + 1] - M_muB[j][i]) / dT for j in range(len(fields_T))
        ]
        # cumulative trapezoid over H, with the anchor (0, 0)
        cumulative: list[float] = []
        acc = 0.0
        h_prev = 0.0
        v_prev = 0.0
        for h, v in zip(fields_T, dmdT):
            acc += 0.5 * (v + v_prev) * (h - h_prev)
            cumulative.append(_N_A_MU_B * acc)
            h_prev, v_prev = h, v
        rows.append({"T_K": 0.5 * (t_low + t_high), "delta_S": cumulative})
    return rows


def _require_levels(energies_cm1: list[float]) -> None:
    if not energies_cm1:
        raise MagnetocaloricError(
            "no spin-orbit level spectrum found. Next step: use a SINGLE_ANISO "
            "output (or a POLY_ANISO output with per-center spectra) whose "
            "'Spin-orbit energy spectra' block is present."
        )


def _magnetization_table(
    magnetization: dict[str, Any],
) -> tuple[list[float], list[float], list[list[float]]]:
    temperatures = magnetization.get("temperatures_K") or []
    fields = magnetization.get("fields_T") or []
    matrix = magnetization.get("M_muB") or []
    if len(temperatures) < 2 or not fields or not matrix:
        raise MagnetocaloricError(
            "no magnetization table found. Next step: this route needs a "
            "POLY_ANISO output with the HINT/TMAG requests (the 'HIGH-FIELD"
            " POWDER MAGNETIZATION' table)."
        )
    return temperatures, fields, matrix


def levels_plot_csv(
    energies_cm1: list[float], temperatures_K: list[float] | None = None
) -> str:
    """The S_mag(T) table as a plot-ready CSV companion (the levels route).

    One row per grid temperature (the report's default grid, or the caller's
    when given) with the magnetic entropy in J mol^-1 K^-1, in the same
    canonical-ensemble convention the report prints (the formats chapter,
    "Plot-ready CSV companions").
    """
    _require_levels(energies_cm1)
    grid = list(temperatures_K) if temperatures_K else list(_DEFAULT_TEMPERATURES)
    values = entropy_from_levels(energies_cm1, grid)
    lines = ["temperature_K,entropy_J_per_mol_per_K"]
    lines += [
        f"{plot_number(temperature)},{plot_number(value)}"
        for temperature, value in zip(grid, values)
    ]
    return "\n".join(lines) + "\n"


def maxwell_plot_csv(magnetization: dict[str, Any]) -> str:
    """The -DeltaS(T, H) table as a plot-ready CSV companion (the Maxwell route).

    Tidy long form -- one row per (temperature midpoint, field) point over
    the full printed field grid (the report shows five probe fields; the
    companion carries them all), values in J mol^-1 K^-1 with the report's
    sign convention: the column is -DeltaS, positive = direct MCE.
    """
    temperatures, fields, matrix = _magnetization_table(magnetization)
    rows = maxwell_delta_s(temperatures, fields, matrix)
    lines = ["T_mid_K,field_T,minus_delta_S_J_per_mol_per_K"]
    for row in rows:
        for field, value in zip(fields, row["delta_S"]):
            lines.append(
                f"{plot_number(row['T_K'])},{plot_number(field)},{plot_number(-value)}"
            )
    return "\n".join(lines) + "\n"


def _summary_points(temperatures: list[float]) -> list[int]:
    wanted = (1, 2, 5, 10, 20, 50, 100, 200, 300)
    return [
        index
        for index, temperature in enumerate(temperatures)
        if any(abs(temperature - value) < 1e-6 for value in wanted)
    ]


def render_levels(
    energies_cm1: list[float],
    *,
    source: str,
    temperatures_K: list[float] | None = None,
) -> str:
    """The menu-41 report for the levels route."""
    _require_levels(energies_cm1)
    grid = temperatures_K or list(_DEFAULT_TEMPERATURES)
    values = entropy_from_levels(energies_cm1, grid)
    high_limit = _R_J_MOL_K * math.log(len(energies_cm1))
    lines = [f"Magnetic entropy / magnetocaloric report ({source})"]
    lines.append(
        f"  route: spin-orbit levels ({len(energies_cm1)} levels,"
        " degeneracy 1 assumed); B = 0"
    )
    lines.append("")
    lines.append("  S_mag(T) [J mol-1 K-1]:")
    for temperature, value in zip(grid, values):
        lines.append(f"    T = {temperature:>6.1f} K   S = {value:9.4f}")
    lines.append(
        f"  high-temperature check: R ln(N) = {high_limit:.4f} J mol-1 K-1;"
        " the engine's level list is truncated, so this is an upper bound --"
        " do not trust S above the saturation temperature"
    )
    lines.append("")
    lines.append(
        "  Reading notes: S(T) = R (ln Z + <E>/kT) from the printed level"
        " spectrum (Szalowski & Kowalewska 2020, Eq. (4)-(5)); only the"
        " truncated lowest levels enter, so the high-T limit is an upper"
        " bound; the entropy is per mole of magnetic units."
    )
    return "\n".join(lines)


def render_maxwell(
    magnetization: dict[str, Any],
    *,
    source: str,
    field_probe_T: tuple[float, ...] = (1.0, 2.0, 3.0, 5.0, 7.0),
) -> str:
    """The menu-41 report for the magnetization (Maxwell) route."""
    temperatures, fields, matrix = _magnetization_table(magnetization)
    rows = maxwell_delta_s(temperatures, fields, matrix)
    lines = [f"Magnetic entropy / magnetocaloric report ({source})"]
    lines.append(
        f"  route: magnetization table ({len(temperatures)} temperatures x"
        f" {len(fields)} fields, powder-averaged, molar M in Bohr magnetons)"
    )
    lines.append("")
    lines.append("  -DeltaS(T, H) [J mol-1 K-1] (positive = direct MCE):")
    header = "    T_mid [K] | " + " | ".join(f"{h:>5.1f} T" for h in field_probe_T)
    lines.append(header)
    for row in rows:
        cells = []
        for probe in field_probe_T:
            index = min(range(len(fields)), key=lambda j: abs(fields[j] - probe))
            cells.append(f"{-row['delta_S'][index]:>7.3f}")
        lines.append(f"    {row['T_K']:>9.3f} | " + " | ".join(cells))
    best_row, best_index, best_value = None, None, -math.inf
    for row in rows:
        for index, value in enumerate(row["delta_S"]):
            if -value > best_value:
                best_row, best_index, best_value = row, index, -value
    if best_row is not None:
        lines.append(
            f"  maximum: -DeltaS = {best_value:.4f} J mol-1 K-1 at"
            f" T = {best_row['T_K']:.2f} K, H = {fields[best_index]:.3f} T"
        )
    lines.append("")
    lines.append(
        "  Reading notes: DeltaS(T, H) = N_A mu_B int (dM/dT)_H dH' (the"
        " Maxwell relation; trapezoid over the printed field grid, adjacent-"
        "temperature differences; the printed sign convention is"
        " DeltaS = S(T,H) - S(T,0), so the table shows -DeltaS). The"
        " data must be converged (a fine temperature grid); the field grid"
        " starts at 0.0001 T."
    )
    return "\n".join(lines)


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the method reference and the measured data chain."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "The magnetocaloric quantities follow the canonical-ensemble"
                " route of Szalowski & Kowalewska 2020: the Gibbs function"
                " G = -kT ln Z, the magnetic entropy S = (U - G)/T, and the"
                " isothermal entropy change DeltaS_T = S(T, Bi) - S(T, Bf);"
                " the same conventions are used in the V6 companion study."
            ),
            ref=(
                "Szalowski K., Kowalewska P., Materials 2020, 13(2), 485, "
                "DOI 10.3390/ma13020485; Kowalewska P., Szalowski K., "
                "J. Magn. Magn. Mater. 2020, 496, 165933, "
                "DOI 10.1016/j.jmmm.2019.165933"
            ),
            url="https://doi.org/10.3390/ma13020485",
            bibkey="szalowski2020mce",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "ORCA 6.1.1 fixture magnetocaloric/ (2026-09-29): the POLY_ANISO"
                " driver with HINT/TMAG prints the 'HIGH-FIELD POWDER"
                " MAGNETIZATION' table in chunks of temperature columns (5 per"
                " chunk; 71 field points over 0-7 T in the probe), which the"
                " parser merges into one aligned matrix; the SINGLE_ANISO"
                " 'Spin-orbit energy spectra' block provides the levels route."
            ),
            ref="ORCA 6.1.1, measured 2026-09-29 (fixtures magnetocaloric/)",
        ),
    )
