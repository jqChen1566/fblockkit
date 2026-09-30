"""The literature batch: quantum-tunnelling relaxation prediction (menu 46).

Two published models over one data chain -- this toolkit's SINGLE_ANISO
reading (per-doublet principal g values and energies, menu 36's
``group_metrics``) plus, for the dipolar model, a neighbour table:

- **equivalent-Zeeman model** (Yin & Li, Phys. Chem. Chem. Phys. 2020, 22,
  9923): the dipolar field is taken as one empirical isotropic scale
  B_ave (their 20 mT); the ground-state tunnelling time is

      omega_QTM = (beta B_ave / h) * gXY^2 / (2 sqrt(gXY^2 + gZ^2)),
      tau_QTM = 1 / (2 omega_QTM),

  with the thermally activated weighting over the doublets
  ``omega_i(T) ~ exp(-E_i/kT)/Z * omega_eff,i`` and
  ``U_eff(T) = sum_i (omega_i/sum omega) E_i``.

- **spin-dipolar model** (Aravena, J. Phys. Chem. Lett. 2018, 9, 5327; the
  non-collinear update: Aravena, Dalton Trans. 2026, 55, 7848): each
  neighbour moment mu_b (pseudospin +/-) drives a central spin flip with
  the matrix element ``<down|H|up> = (beta/2)(g_x B_x + i g_y B_y)`` where
  ``B = (mu_0/4 pi r^3)(3 (r_hat . mu_b) r_hat - mu_b)``.  Independent +/-
  environment spins make the real/imaginary parts sums of zero-mean
  independent variables, so the variances are closed-form:

      sigma_r^2 = sum_j Re(c_j)^2,  sigma_i^2 = sum_j Im(c_j)^2,
      sigma_t^2 = sigma_r^2 + sigma_i^2,  <|E_sf|> = sqrt(2/pi) sigma_t,
      k = (2 pi / hbar) <|E_sf|>^2,  tau_QT = 1 / (2 k).

  The magnetic-dilution variant (Llanos & Aravena, J. Magn. Magn. Mater.
  2019, 489, 165456) keeps each neighbour with probability x and takes the
  median tau over repeats.

Axial convention: the model formulas assume an Ising-type ground doublet;
this module takes the **largest principal g** as the axial component and
the remaining two as the transverse pair (the two enter gXY^2
symmetrically).  Regression: the 18-complex literature table
(``fixtures/qtm/tau_zeeman_benchmark.yaml``; g from the JPCL SI, log tau
from the PCCP Table 1) reproduces to a maximum deviation of 0.006 in
log10(tau).

Boundaries kept in the report: the models are zero-field, Kramers-ion and
single-centre; B_ave is an empirical scale (the paper offers it as
adjustable); the dipolar model needs the neighbour geometry, which the
caller supplies as a table (r vectors and moment vectors in the central
g-frame; crystal-structure parsing sits outside this menu).
"""

from __future__ import annotations

import math
import random
import re
from dataclasses import dataclass

from ..knowledge.models import EVIDENCE_LITERATURE, EVIDENCE_MEASURED, Evidence

__all__ = [
    "TunnellingError",
    "KdLevel",
    "Neighbour",
    "axial_parts",
    "zeeman_rate",
    "zeeman_tau",
    "ueff_curve",
    "parse_neighbour_table",
    "dipolar_sigma",
    "dipolar_tau",
    "dilution_medians",
    "render",
    "evidence",
]

# CODATA-derived constants (SI base units).
_MU_B = 9.2740100783e-24       # J/T (= A m^2)
_HBAR = 1.054571817e-34        # J s
_H = 6.62607015e-34            # J s
_MU_0 = 1.25663706212e-6       # N/A^2
_K_B = 1.380649e-23            # J/K
_K_B_CM1_PER_K = 0.695034800   # cm^-1 per K (as in analysis/relaxation.py)

_BOHR_TO_ANGSTROM = 0.529177210903  # not used for input (neighbours are in Angstrom)

DEFAULT_B_AVE_MT = 20.0

# the summary temperature grid (log-spaced like the magnetocaloric menu's use)
DEFAULT_TEMPERATURES_K = (2.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 150.0, 200.0, 300.0)


class TunnellingError(ValueError):
    """The requested prediction cannot be computed (with a next step)."""


@dataclass(frozen=True)
class KdLevel:
    """One Kramers doublet: index, energy (cm^-1, relative) and principal g."""

    index: int
    energy_cm1: float
    g: tuple[float, float, float]


@dataclass(frozen=True)
class Neighbour:
    """One environment centre: position (Angstrom) and moment (Bohr magnetons)."""

    r: tuple[float, float, float]
    moment: tuple[float, float, float]


# --- the equivalent-Zeeman model (Yin & Li 2020) ------------------------------


def axial_parts(g: tuple[float, float, float]) -> tuple[float, float]:
    """(gXY, gZ) with the largest principal value taken as axial."""
    values = sorted((abs(value) for value in g), reverse=True)
    gz = values[0]
    gxy = math.hypot(values[1], values[2])
    return gxy, gz


def zeeman_rate(g: tuple[float, float, float], *, B_ave_mT: float = DEFAULT_B_AVE_MT) -> float:
    """The ground-doublet QTM rate omega (s^-1), Yin & Li eq (9b)."""
    if B_ave_mT <= 0:
        raise TunnellingError("B_ave must be positive.")
    gxy, gz = axial_parts(g)
    denominator = 2.0 * math.hypot(gxy, gz)
    if gxy == 0.0 or denominator == 0.0:
        raise TunnellingError(
            "the transverse g of this doublet is zero: the equivalent-Zeeman model "
            "has no tunnelling channel (perfectly axial, no relaxation by this route)."
        )
    return (_MU_B * (B_ave_mT * 1e-3) / _H) * (gxy * gxy) / denominator


def zeeman_tau(g: tuple[float, float, float], *, B_ave_mT: float = DEFAULT_B_AVE_MT) -> float:
    """The ground-doublet tunnelling time tau_QTM (s), Yin & Li eq (9d)."""
    return 1.0 / (2.0 * zeeman_rate(g, B_ave_mT=B_ave_mT))


def kd_levels(rows: list[dict]) -> list[KdLevel]:
    """The per-doublet levels from menu 36's ``group_metrics`` rows.

    Only Kramers doublets (half-integer effective spin) enter; the energies
    are taken relative to the lowest doublet (the printed spin-orbit spectra
    start at the ground state, but the subtraction is defensive).  Raises
    when no Kramers doublet is present (the models are restricted).
    """
    kramers = [row for row in rows if row.get("kramers")]
    if not kramers:
        raise TunnellingError(
            "no Kramers doublet found in this SINGLE_ANISO segment. The models are "
            "restricted to Kramers (half-integer spin) ions. Next step: give an "
            "output of a Kramers ion, or use menu 36 for the relaxation reading."
        )
    ground = min(row["energy_cm1"] for row in kramers if row["energy_cm1"] is not None)
    levels: list[KdLevel] = []
    for row in kramers:
        if row["energy_cm1"] is None:
            continue
        levels.append(
            KdLevel(
                index=row["index"],
                energy_cm1=float(row["energy_cm1"]) - ground,
                g=tuple(float(v) for v in row["g_principal"]),
            )
        )
    levels.sort(key=lambda level: level.energy_cm1)
    return levels


def ueff_curve(
    levels: list[KdLevel],
    *,
    B_ave_mT: float = DEFAULT_B_AVE_MT,
    temperatures_K: tuple[float, ...] = DEFAULT_TEMPERATURES_K,
) -> list[dict]:
    """U_eff(T) and the per-doublet contributions (Yin & Li eq (10)-(11)).

    Each row is ``{"T_K", "ueff_cm1", "contributions": [(index, weight), ...]}``
    where the weights are ``omega_i(T) / sum omega``.  At low temperature the
    ground doublet dominates and U_eff tends to zero (the ground state does
    not contribute to the barrier), rising to the Orbach plateau.
    """
    weights_per_kd: list[float] = []
    for level in levels:
        gxy, gz = axial_parts(level.g)
        denominator = 2.0 * math.hypot(gxy, gz)
        weights_per_kd.append(0.0 if denominator == 0.0 else (gxy * gxy) / denominator)
    rows: list[dict] = []
    for temperature in temperatures_K:
        if temperature <= 0:
            raise TunnellingError("temperatures must be positive.")
        activated: list[float] = []
        for level, weight in zip(levels, weights_per_kd):
            exponent = -level.energy_cm1 / (_K_B_CM1_PER_K * temperature)
            # a strongly negative exponent underflows to exactly 0.0, which is
            # the physically wanted value here (no clamping: a clamp would
            # leave a spurious non-zero floor)
            activated.append(weight * math.exp(exponent))
        total = sum(activated)
        if total == 0.0:
            rows.append({"T_K": temperature, "ueff_cm1": 0.0, "contributions": []})
            continue
        contributions = [
            (level.index, value / total) for level, value in zip(levels, activated)
        ]
        ueff = sum(
            (value / total) * level.energy_cm1
            for level, value in zip(levels, activated)
        )
        rows.append({"T_K": temperature, "ueff_cm1": ueff, "contributions": contributions})
    return rows


# --- the spin-dipolar model (Aravena 2018/2026; dilution 2019) ----------------


def parse_neighbour_table(text: str) -> tuple[Neighbour, ...]:
    """Parse the neighbour table: ``dx dy dz mx my mz`` per line.

    Positions in Angstrom, moment components in Bohr magnetons -- both in the
    central ion's principal-g frame (x, y across, z along the main axis).
    Lines starting with ``#`` and blank lines are skipped.  A line with a
    different column count raises with the next step.
    """
    neighbours: list[Neighbour] = []
    for number, line in enumerate(text.splitlines(), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = stripped.split()
        if len(fields) != 6:
            raise TunnellingError(
                f"neighbour table line {number} carries {len(fields)} columns; the "
                "format is 'dx dy dz mx my mz' (Angstrom; Bohr magnetons). Next "
                "step: give the six numbers per line, one neighbour per line."
            )
        try:
            values = [float(field) for field in fields]
        except ValueError:
            raise TunnellingError(
                f"neighbour table line {number} has a non-numeric field: {stripped!r}."
            ) from None
        r = tuple(values[0:3])
        if math.dist((0.0, 0.0, 0.0), r) == 0.0:
            raise TunnellingError(
                f"neighbour table line {number} sits at the central ion (zero "
                "position vector); drop it."
            )
        neighbours.append(Neighbour(r=r, moment=tuple(values[3:6])))
    if not neighbours:
        raise TunnellingError(
            "the neighbour table is empty. Next step: one line per environment "
            "centre, 'dx dy dz mx my mz'."
        )
    return tuple(neighbours)


def _dipolar_coefficient(
    neighbour: Neighbour, g_center: tuple[float, float, float]
) -> complex:
    """The complex spin-flip coefficient c for one neighbour (one sign choice).

    ``c = (beta/2)(g_x B_x + i g_y B_y)`` with the dipolar field
    ``B = (mu_0 / 4 pi r^3)(3 (r_hat . mu_b) r_hat - mu_b)``; the moment is
    taken in Bohr magnetons, the position in Angstrom.  The opposite
    environment-pseudospin sign flips c's sign, so both signs are covered by
    the variance sum.
    """
    rx, ry, rz = neighbour.r
    r_m = [value * 1e-10 for value in (rx, ry, rz)]
    distance = math.sqrt(sum(value * value for value in r_m))
    if distance == 0.0:
        raise TunnellingError("a zero-length neighbour vector reached the field model.")
    if distance < 0.1e-10:
        raise TunnellingError(
            f"a neighbour sits {distance * 1e10:.4f} Angstrom from the central ion: "
            "closer than any physical contact and numerically explosive in 1/r^3. "
            "Next step: check the table's units (Angstrom, in the central g-frame)."
        )
    moment = [value * _MU_B for value in neighbour.moment]
    rhat = [value / distance for value in r_m]
    r_dot_mu = sum(a * b for a, b in zip(rhat, moment))
    field = [
        (_MU_0 / (4.0 * math.pi * distance**3)) * (3.0 * r_dot_mu * rh - mu)
        for rh, mu in zip(rhat, moment)
    ]
    gxy, _ = axial_parts(g_center)
    # the transverse pair, in the table's own x/y assignment after the axial
    # convention: keep the two non-largest principal values
    values = sorted((abs(value) for value in g_center), reverse=True)
    gx, gy = values[1], values[2]
    return (_MU_B / 2.0) * complex(gx * field[0], gy * field[1])


def dipolar_sigma(
    neighbours: tuple[Neighbour, ...], g_center: tuple[float, float, float]
) -> tuple[float, float]:
    """(sigma_r^2, sigma_i^2) over the independent-sign environment sums."""
    sigma_r2 = 0.0
    sigma_i2 = 0.0
    for neighbour in neighbours:
        coefficient = _dipolar_coefficient(neighbour, g_center)
        sigma_r2 += coefficient.real**2
        sigma_i2 += coefficient.imag**2
    return sigma_r2, sigma_i2


def dipolar_tau(
    neighbours: tuple[Neighbour, ...], g_center: tuple[float, float, float]
) -> dict:
    """The spin-dipolar tunnelling time, from the closed-form variance sum.

    Returns ``{"tau_s", "sigma_r", "sigma_i", "e_sf_av_J", "e_sf_av_K"}``
    (+/- environment signs enter as independent zero-mean variables; the
    average absolute matrix element is sqrt(2/pi) sigma_t).

    Fermi-golden-rule constant (calibrated against the source tables): the
    published rates reproduce with the density factor 1/k_B, i.e.
    ``k = (2 pi / hbar) (k_B <|E_sf|>_K)^2`` -- the source's <|E_sf|> is
    quoted in kelvin (JPCL SI Table S1) and 1/k_B is the per-kelvin state
    factor of their one-unit window.  The check: NAFMIT's published
    <|E_sf|> = 8.671e-9 K gives log10(tau) = 3.906 against the published
    3.91, and BAJSIQ's 1.530e-4 K gives -4.58 against -4.59.
    """
    sigma_r2, sigma_i2 = dipolar_sigma(neighbours, g_center)
    sigma_t = math.sqrt(sigma_r2 + sigma_i2)
    e_sf_av = math.sqrt(2.0 / math.pi) * sigma_t
    rate = (2.0 * math.pi / _HBAR) * e_sf_av**2 / _K_B
    tau = 1.0 / (2.0 * rate)
    return {
        "tau_s": tau,
        "sigma_r": math.sqrt(sigma_r2),
        "sigma_i": math.sqrt(sigma_i2),
        "e_sf_av_J": e_sf_av,
        "e_sf_av_K": e_sf_av / _K_B,
    }


def dilution_medians(
    neighbours: tuple[Neighbour, ...],
    g_center: tuple[float, float, float],
    concentrations: tuple[float, ...],
    *,
    repeats: int = 60,
    seed: int = 20180930,
) -> list[dict]:
    """The dilution variant: keep each neighbour with probability x, median tau.

    One row per concentration: ``{"x", "log_tau_median", "n_active_median"}``.
    The random stream is seeded (default fixed) so the table replays
    byte-identically.
    """
    if repeats < 1:
        raise TunnellingError("repeats must be at least 1.")
    rng = random.Random(seed)
    rows: list[dict] = []
    for x in concentrations:
        if not 0.0 < x <= 1.0:
            raise TunnellingError("dilution concentrations must lie in (0, 1].")
        taus: list[float] = []
        actives: list[int] = []
        for _ in range(repeats):
            kept = tuple(neighbour for neighbour in neighbours if rng.random() < x)
            if not kept:
                continue
            taus.append(dipolar_tau(kept, g_center)["tau_s"])
            actives.append(len(kept))
        if not taus:
            rows.append({"x": x, "log_tau_median": None, "n_active_median": 0})
            continue
        taus.sort()
        middle = len(taus) // 2
        median = (
            taus[middle]
            if len(taus) % 2 == 1
            else math.sqrt(taus[middle - 1] * taus[middle])
        )
        actives.sort()
        rows.append(
            {
                "x": x,
                "log_tau_median": math.log10(median),
                "n_active_median": actives[len(actives) // 2],
            }
        )
    return rows


# --- the report ---------------------------------------------------------------


def _format_seconds(tau: float) -> str:
    if tau is None:
        return "n/a"
    if not math.isfinite(tau) or tau <= 0.0:
        return "effectively infinite (no active tunnelling channel in this window)"
    if math.log10(tau) > 300:
        return "> 10^300 s (effectively infinite)"
    if -3 <= math.floor(math.log10(tau)) <= 6:
        return f"{tau:.3e} s"
    return f"{tau:.3e} s (= 10^{math.log10(tau):.2f})"


def render(
    *,
    source: str,
    levels: list[KdLevel],
    B_ave_mT: float,
    ueff_rows: list[dict],
    neighbour_report: dict | None = None,
) -> str:
    """The report body (plain text, English, no timestamps of the day)."""
    lines: list[str] = [
        "Quantum-tunnelling relaxation prediction",
        f"Source: {source}",
        f"Model: equivalent-Zeeman (Yin & Li 2020), B_ave = {B_ave_mT:.1f} mT",
        "",
        "Kramers doublets (axial g = the largest principal value):",
        "   KD      E (cm-1)          g_xy          g_z      tau_QTM (s)",
    ]
    ground_tau = None
    for level in levels:
        gxy, gz = axial_parts(level.g)
        try:
            tau = zeeman_tau(level.g, B_ave_mT=B_ave_mT)
        except TunnellingError:
            tau = None
        if ground_tau is None:
            ground_tau = tau
        tau_text = _format_seconds(tau) if tau is not None else "n/a (no transverse g)"
        lines.append(
            f"  {level.index:>3}  {level.energy_cm1:>10.3f}  {gxy:>13.8f}  {gz:>12.8f}   {tau_text}"
        )
    if ground_tau is not None:
        lines.append("")
        lines.append(f"Ground-state tau_QTM: {_format_seconds(ground_tau)}")
    lines.append("")
    lines.append("U_eff(T) (thermal-activation weighting over the doublets):")
    lines.append("     T (K)   U_eff (cm-1)   contributions (KD: weight)")
    for row in ueff_rows:
        contributions = ", ".join(
            f"{index}: {weight:.3f}" for index, weight in row["contributions"][:4]
        )
        lines.append(f"  {row['T_K']:>7.1f}   {row['ueff_cm1']:>11.3f}   {contributions}")
    if len(levels) > 1 and ueff_rows and max(r["ueff_cm1"] for r in ueff_rows) < 1e-9:
        first = levels[1].energy_cm1
        lines.append(
            f"  (the first excited doublet lies at {first:.1f} cm-1 "
            f"= {first / _K_B_CM1_PER_K:.0f} K, far above this temperature grid: "
            "the thermally activated regime of the model is not reached in the "
            "printed window, and U_eff stays at the ground-state value)"
        )
    if neighbour_report is not None:
        lines.extend(_render_neighbours(neighbour_report))
    lines.append("")
    lines.append(
        "Reading notes: the models are zero-field, Kramers-ion, single-centre "
        "descriptions of the tunnelling (QTM) relaxation; B_ave is an empirical "
        "field scale (20 mT in the source, adjustable -- the paper states it can "
        "be treated as an empirical parameter). The low-temperature U_eff tends "
        "to zero (the ground doublet does not contribute to the barrier) and "
        "rises to the Orbach plateau. These predictions are absolute-value "
        "estimates from ab initio parameters; menu 36's reading of an "
        "experimental-style plateau is a separate, data-side quantity -- keep "
        "the two apart when quoting them side by side. The dipolar model needs "
        "the neighbour geometry as a table (crystal-structure parsing sits "
        "outside this menu's scope)."
    )
    return "\n".join(lines)


def _render_neighbours(report: dict) -> list[str]:
    lines: list[str] = [
        "",
        "Spin-dipolar model (Aravena 2018; non-collinear form 2026)",
        f"Neighbour table: {report['table']} ({report['n_neighbours']} centres)",
        f"  sigma_r = {report['sigma_r']:.6e}  sigma_i = {report['sigma_i']:.6e}  (J)",
        f"  <|E_sf|> = {report['e_sf_av_J']:.6e} J",
        f"  tau_QT = {_format_seconds(report['tau_s'])}",
    ]
    if report.get("dilution"):
        lines.append("  Dilution variant (median over repeats, seeded):")
        lines.append("        x     log10(tau/s)   active centres (median)")
        for row in report["dilution"]:
            value = "n/a" if row["log_tau_median"] is None else f"{row['log_tau_median']:+.3f}"
            lines.append(f"    {row['x']:>6.3f}   {value:>10}   {row['n_active_median']}")
    return lines


def evidence() -> tuple[Evidence, ...]:
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "The equivalent-Zeeman tunnelling model (B_ave = 20 mT; eqs (9b), "
                "(10)-(11)) whose regression table this menu reproduces."
            ),
            ref="Yin B., Li C.-C., Phys. Chem. Chem. Phys., 2020, 22, 9923-9933, DOI 10.1039/D0CP00933D",
            bibkey="yin2020qtm",
        ),
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "The spin-dipolar tunnelling model (spin-flip matrix elements, "
                "Fermi golden rule) and its validation set."
            ),
            ref="Aravena D., J. Phys. Chem. Lett., 2018, 9, 5327-5333, DOI 10.1021/acs.jpclett.8b02359",
            bibkey="aravena2018tunneling",
        ),
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "The non-collinear update of the dipolar model (effective "
                "transverse field, eqs (11)-(12)) used for the neighbour-sum "
                "implementation."
            ),
            ref="Aravena D., Dalton Trans., 2026, 55, 7848-7857, DOI 10.1039/d6dt00521g",
            bibkey="aravena2026packing",
        ),
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text="The magnetic-dilution variant (probabilistic neighbours, median over repeats).",
            ref="Llanos L., Aravena D., J. Magn. Magn. Mater., 2019, 489, 165456, DOI 10.1016/j.jmmm.2019.165456",
            bibkey="llanos2019dilution",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "Regression measured 2026-09-30: the 18-complex literature table "
                "(g from the JPCL SI, log tau from the PCCP Table 1) reproduces to "
                "a maximum deviation of 0.006 in log10(tau) with B_ave = 20 mT "
                "and the CODATA constants of this module."
            ),
            ref="fixtures/qtm/tau_zeeman_benchmark.yaml",
        ),
    )
