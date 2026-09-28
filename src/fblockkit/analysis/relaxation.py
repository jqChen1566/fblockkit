"""5.3: magnetic relaxation and QTM metrics (menu 36).

Both data sources are ORCA 6.1.1 outputs (fixtures ``single_aniso/`` and
``magrelax/``; the close reading is ``文献细读/细读_磁弛豫-QTM_Chilton2025.md``):

- the **SINGLE_ANISO** embedded output (the ``ANISO`` sub-block of
  ``%casscf``): per-KD g-tensors, D/E, the SOC spectrum and the UBAR
  magnetic-moment matrix elements;
- the **Orca_Magrelax** output: the relaxation-rate table tau(T).

The metrics follow the practical guide (N. F. Chilton, Chem. Soc. Rev.
2025, 54, 11468).  For Kramers ions the guide bounds the effective barrier
by the first excited Kramers doublet that is either non-collinear with the
ground state (the guide's working threshold: more than 15 degrees) or has a
significant transverse g -- quantified by its empirical rule from a survey
of 20 high-performance Dy(III) SMMs, ``g_T * theta_3 > 20`` with
``g_T = (g1 + g2 + g3 sin(theta_3))/3`` (g1 <= g2 <= g3 the excited KD's
principal values, theta_3 the angle between the two largest-g axes).
**Both thresholds are empirical** and are recorded with every report;
nothing here is a dynamical simulation.  The UBAR matrix elements
(|mu_X| + |mu_Y| + |mu_Z|)/3 (the engine's printed average) are the
qualitative QTM intensity of the guide's section 7.3 (whose own text calls
the underlying magnetic-perturbation model "not realistic").

The Orbach analysis is the textbook Arrhenius fit tau(T) =
tau_0 exp(U_eff / k_B T) over the finite positive points (optionally a
temperature window), reporting U_eff in K and cm^-1, tau_0 in seconds and
the fit's R^2 -- a data reduction, not a prediction.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..knowledge.models import EVIDENCE_LITERATURE, EVIDENCE_MEASURED, Evidence

__all__ = [
    "RelaxationError",
    "MagrelaxData",
    "read_magrelax",
    "group_metrics",
    "ueff_estimate",
    "orbach_fit",
    "render_single_aniso",
    "render_magrelax",
    "evidence",
]

#: h c / k_B in cm K (the same constant as the pnmr module)
_CM1_PER_K = 1.438776877

#: the guide's working non-collinearity threshold (degrees; Chilton section 7.1)
COLLINEARITY_THRESHOLD_DEG = 15.0

#: the guide's empirical transverse-g rule (Chilton section 7.2, 20 Dy(III) SMMs)
G_T_THETA_THRESHOLD = 20.0

#: instrumental threshold (not from the source): below this relative spread of
#: the g principal values the main-axis direction is decided by the print's own
#: digits, so the axis-angle criteria are not applicable to that group.  The CO+
#: fixture (Delta g / g ~ 2e-4) is the measured case that motivated it.
AXIS_DEFINITION_MIN = 1e-2

_ZEEMAN_ROW_RE = re.compile(r"^(\d+)\s+([+-]?[\d.]+(?:E[+-]?\d+)?)\s*$")
_RATE_ROW_RE = re.compile(
    r"^([\d.]+)\s+([+-][\d.]+[eE][+-]\d+)\s+([+-]?(?:inf|[\d.]+[eE]?[+-]?\d*))\s*$"
)


class RelaxationError(Exception):
    """A refusal or a malformed input, with the next step spelled out."""


# --- the SINGLE_ANISO metrics -------------------------------------------------


def _angle_deg(source: np.ndarray, reference: np.ndarray) -> float:
    cosine = abs(float(np.dot(source, reference))) / (
        float(np.linalg.norm(source)) * float(np.linalg.norm(reference))
    )
    return float(np.degrees(np.arccos(min(1.0, cosine))))


def _is_kramers(spin: float) -> bool:
    """Half-integer effective spin (Kramers theorem applies)."""
    twice = round(spin * 2)
    return abs(spin * 2 - twice) < 1e-9 and twice % 2 == 1


def group_metrics(segment: dict) -> list[dict]:
    """Per-pseudospin-group QTM metrics, all angles relative to the ground group.

    A group is Kramers when its effective spin is half-integer.  The g-tensor
    axes come from the printed main-magnetic-axes table; the angle is measured
    between the two largest-g axes (an axis direction is only defined up to
    sign, so the absolute cosine is taken).
    """
    groups = segment.get("groups") or []
    rows: list[dict] = []
    for group in groups:
        values = group.get("g_values")
        axes = group.get("g_axes")
        if not values or not axes or values[0] is None:
            continue
        principal = np.array(values, dtype=float)
        axis_vectors = np.array(axes, dtype=float)
        rows.append(
            {
                "index": group["index"],
                "spin": group["spin"],
                "kramers": _is_kramers(group["spin"]),
                "energy_cm1": min(group["soc_state_energies_cm1"])
                if group["soc_state_energies_cm1"]
                else None,
                "state_energies_cm1": list(group["soc_state_energies_cm1"]),
                "g_principal": principal,
                "g_max_axis": axis_vectors[int(np.argmax(principal))],
                "g_axes": axis_vectors,
                "d_cm1": group.get("d_cm1"),
                "e_cm1": group.get("e_cm1"),
            }
        )
    if not rows:
        return rows
    for row in rows:
        principal = row["g_principal"]
        row["anisotropy_ratio"] = float(
            (principal.max() - principal.min()) / principal.max()
        ) if principal.max() > 0 else 0.0
        row["axis_defined"] = row["anisotropy_ratio"] >= AXIS_DEFINITION_MIN
    ground = rows[0]
    ground_axis = ground["g_max_axis"]
    for row in rows:
        principal = row["g_principal"]
        axis = row["g_max_axis"]
        row["theta_deg"] = _angle_deg(axis, ground_axis)
        ordered = np.sort(principal)
        # the guide's empirical average transverse g (Chilton section 7.2):
        # g_T = (g1 + g2 + g3 sin(theta_3))/3, the smallest two principal
        # values plus the largest one weighted by the tilt of its axis; the
        # naive min-two-norm is reported alongside as a second reading
        row["g_t"] = float(
            (ordered[0] + ordered[1] + ordered[2] * np.sin(np.radians(row["theta_deg"]))) / 3.0
        )
        row["theta_reliable"] = bool(row["axis_defined"] and ground["axis_defined"])
        row["g_t_theta"] = (
            float(row["g_t"] * row["theta_deg"]) if row["theta_reliable"] else None
        )
        row["g_perp_naive"] = float(np.hypot(ordered[0], ordered[1]))
        energies = row["state_energies_cm1"]
        row["intra_splitting_cm1"] = (
            float(max(energies) - min(energies)) if len(energies) > 1 else None
        )
    return rows


def ueff_estimate(rows: list[dict]) -> dict | None:
    """The first excited group flagged by the guide's empirical criteria.

    Flags: theta > 15 deg (non-collinear) or g_T * theta > 20 (the guide's
    transverse-g rule), applied only to groups whose main axes are defined
    (relative g spread at least ``AXIS_DEFINITION_MIN``); for a nearly
    isotropic group the axis direction is decided by the print's digits.
    Returns the flagged row, or None when no computed group qualifies.
    """
    for row in rows[1:]:
        if not row["theta_reliable"]:
            continue
        if row["theta_deg"] > COLLINEARITY_THRESHOLD_DEG or (
            row["g_t_theta"] is not None and row["g_t_theta"] > G_T_THETA_THRESHOLD
        ):
            return row
    return None


# --- the Orca_Magrelax output -------------------------------------------------


@dataclass
class MagrelaxData:
    """The parsed Orca_Magrelax section of an ORCA output."""

    magrelax_input: str = ""
    gbw_file: str | None = None
    derivatives_file: str | None = None
    n_derivatives: int | None = None
    derivative_size: int | None = None
    hessian_file: str | None = None
    zeeman_eigenvalues_cm1: list[float] = field(default_factory=list)
    frequencies_cm1: list[float] = field(default_factory=list)
    temperatures_k: list[float] = field(default_factory=list)
    rates_s1: list[float] = field(default_factory=list)
    tau_s: list[float] = field(default_factory=list)


def _tail_value(line: str) -> str:
    if "=" in line:
        tail = line.rsplit("=", 1)[1]
    else:  # e.g. "Reading Derivatives from ......file" (no equals sign)
        tail = line.split()[-1]
    return tail.strip().lstrip(".").strip()


def read_magrelax(path: str | Path) -> MagrelaxData:
    """Read the Orca_Magrelax block (the relaxation-rate table and its header).

    Raises ``RelaxationError`` when the file carries no MAGRELAX section, so a
    wrong file cannot silently produce an empty report.
    """
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    if "ORCA MAGRELAX" not in text:
        raise RelaxationError(
            f"{path} has no Orca_Magrelax section (no 'ORCA MAGRELAX' banner). "
            "Next step: this menu reads the output of an ORCA run with "
            "%magrelax settings, or with DoMagrelax in the CASSCF rel block; "
            "for the static QTM metrics give a SINGLE_ANISO output instead."
        )
    data = MagrelaxData()
    lines = text.splitlines()
    in_zeeman = False
    in_freq = False
    in_rates = False
    for index, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("Magrelax Input Name"):
            data.magrelax_input = _tail_value(stripped)
        elif stripped.startswith("GBWFile Name"):
            data.gbw_file = _tail_value(stripped)
        elif stripped.startswith("Reading Derivatives from"):
            data.derivatives_file = _tail_value(stripped)
        elif stripped.startswith("Number of derivatives"):
            data.n_derivatives = int(_tail_value(stripped))
        elif stripped.startswith("Size of derivative matrices"):
            data.derivative_size = int(_tail_value(stripped).split("X")[0])
        elif stripped.startswith("Hessianfile"):
            data.hessian_file = _tail_value(stripped)
        elif stripped.startswith("Eigenvalues of zeeman"):
            in_zeeman = True
            in_freq = in_rates = False
        elif stripped.startswith("List of frequencies"):
            in_zeeman = False
            in_freq = True
            in_rates = False
        elif stripped.startswith("Relaxation rates"):
            in_freq = False
            in_rates = True
        elif in_zeeman and (match := _ZEEMAN_ROW_RE.match(stripped)) is not None:
            data.zeeman_eigenvalues_cm1.append(float(match.group(2)))
        elif in_freq and stripped and not stripped.startswith("*"):
            tokens = stripped.split()
            if len(tokens) == 1:
                try:
                    data.frequencies_cm1.append(float(tokens[0]))
                except ValueError:
                    in_freq = False
            elif not stripped.startswith("T (K)"):
                in_freq = False
        elif in_rates and (match := _RATE_ROW_RE.match(stripped)) is not None:
            data.temperatures_k.append(float(match.group(1)))
            data.rates_s1.append(float(match.group(2)))
            data.tau_s.append(float(match.group(3)))
    return data


def orbach_fit(
    temperatures_k: list[float],
    tau_s: list[float],
    window_k: tuple[float, float] | None = None,
) -> dict:
    """The Arrhenius/Orbach fit ln tau = ln tau_0 + U_eff/(k_B T).

    Only finite positive tau points inside the window are used; at least three
    are required (with fewer, a straight line is not a test of anything).
    """
    points = [
        (temperature, tau)
        for temperature, tau in zip(temperatures_k, tau_s)
        if np.isfinite(tau) and tau > 0.0 and (window_k is None or window_k[0] <= temperature <= window_k[1])
    ]
    if len(points) < 3:
        raise RelaxationError(
            f"the Orbach fit needs at least three finite positive tau points "
            f"(got {len(points)}). Next step: widen the temperature window, or check "
            "whether this run has any active one-phonon channel at all (a single-"
            "mode toy may give an all-zero rate table)."
        )
    x = np.array([1.0 / temperature for temperature, _ in points])
    y = np.array([np.log(tau) for _, tau in points])
    slope, intercept = np.polyfit(x, y, 1)
    residual = y - (slope * x + intercept)
    total = y - y.mean()
    r_squared = 1.0 - float(residual @ residual) / float(total @ total) if len(points) > 1 else 1.0
    return {
        "n_points": len(points),
        "ueff_k": float(slope),
        "ueff_cm1": float(slope) / _CM1_PER_K,
        "tau0_s": float(np.exp(intercept)),
        "r_squared": r_squared,
    }


# --- the report -----------------------------------------------------------------


def render_single_aniso(segment: dict, *, segment_index: int, n_segments: int, source: str) -> str:
    """The QTM part of the menu-36 report for one SINGLE_ANISO segment."""
    lines: list[str] = []
    lines.append(f"Single_aniso segment {segment_index + 1} of {n_segments} ({source})")
    if n_segments > 1:
        lines.append(
            "  Note: multiple segments are printed when the QDPT pass has more than one"
        )
        lines.append(
            "  variant (e.g. DoSSC on: the later segment includes the spin-spin part);"
        )
        lines.append("  compare their spectra before reading numbers across segments.")
    spectrum = segment.get("soc_spectrum_cm1") or []
    if spectrum:
        lines.append(
            "  spin-orbit spectrum (cm-1): "
            + " ".join(f"{value:.3f}" for value in spectrum[:12])
            + (" ..." if len(spectrum) > 12 else "")
        )
    rows = group_metrics(segment)
    if not rows:
        lines.append("  no pseudospin groups were printed (MLTP is required for the g/D analysis)")
    else:
        lines.append("")
        lines.append(
            "  KD / group   S     E (cm-1)      g1        g2        g3      theta(deg)  g_T   g_T*theta"
        )
        for row in rows:
            ordered = np.sort(row["g_principal"])
            energy = f"{row['energy_cm1']:.3f}" if row["energy_cm1"] is not None else "n/a"
            gt_theta = (
                f"{row['g_t_theta']:8.2f}" if row["g_t_theta"] is not None else "     n/a"
            )
            flag = "" if row["axis_defined"] else "  *"
            lines.append(
                f"  {row['index']:>4}        {row['spin']:>4}  {energy:>10}  "
                f"{ordered[0]:9.5f} {ordered[1]:9.5f} {ordered[2]:9.5f}  "
                f"{row['theta_deg']:9.2f}  {row['g_t']:5.3f}  {gt_theta}{flag}"
            )
        if any(not row["axis_defined"] for row in rows):
            lines.append(
                f"  * the g anisotropy of this group is below {AXIS_DEFINITION_MIN:.0%} of its"
                " largest value: its main-axis direction is set by the print's digits, so"
                " the axis-angle criteria do not apply to it (CO+ fixture: Delta g/g ~ 2e-4)"
            )
        estimate = ueff_estimate(rows)
        lines.append("")
        if estimate is not None:
            energy_text = (
                f"{estimate['energy_cm1']:.3f} cm-1"
                if estimate["energy_cm1"] is not None
                else "a group without a printed energy"
            )
            lines.append(
                f"  U_eff estimate (Chilton 2025 criteria): first flagged excited group is"
                f" {estimate['index']} at {energy_text}"
                f" (theta = {estimate['theta_deg']:.2f} deg, g_T*theta = {estimate['g_t_theta']:.2f})"
            )
        else:
            highest = max(row["energy_cm1"] or 0.0 for row in rows)
            skipped = sum(1 for row in rows[1:] if not row["theta_reliable"])
            lines.append(
                "  U_eff estimate: no computed group is flagged"
                + (
                    f" ({skipped} of them are skipped: main axes not defined)"
                    if skipped
                    else ""
                )
                + f"; U_eff is at least {highest:.1f} cm-1 as far as the computed groups go"
            )
        lines.append(
            "  The flags are the guide's empirical thresholds (15 deg non-collinearity;"
            " g_T*theta > 20 from 20 Dy(III) SMMs) -- they bound a plausible barrier,"
            " not a measured one."
        )
    ubar = segment.get("ubar")
    if ubar and ubar.get("present"):
        lines.append("")
        lines.append("  UBAR (ab initio blocking barrier, qualitative per the engine):")
        for state in ubar["zeeman_states"]:
            lines.append(
                f"    group {state['mult']}: m+ = {state['m_plus']:+.6f}"
                + (f"  m0 = {state['m_zero']:+.6f}" if state["m_zero"] is not None else "")
                + f"  m- = {state['m_minus']:+.6f}  E = {state['energy_cm1']:.4f} cm-1"
            )
        intra = [
            element
            for element in ubar["matrix_elements"]
            if element["delta_mult"] == 0 and element["average"] is not None
        ]
        if intra:
            lines.append("    intra-group |<i+|mu|i->| averages (mu_B): "
                         + "  ".join(
                             f"g{element['source_mult']}={element['average']:.3e}" for element in intra
                         ))
        inter = [
            element
            for element in ubar["matrix_elements"]
            if element["delta_mult"] != 0 and element["average"] is not None
        ]
        if inter:
            largest = max(inter, key=lambda element: element["average"])
            lines.append(
                f"    largest inter-group element: |<{largest['source_mult']}.{largest['source_state']}"
                f"|mu_{largest['component']}|{largest['target_mult']}.{largest['target_state']}>|"
                f" = {largest['average']:.4e} mu_B"
            )
        if ubar.get("even_electron_note"):
            lines.append(
                "    the 'check the tunnelling splitting instead' note is the engine's fixed"
                " template sentence (printed for odd-electron systems too)"
            )
    if segment.get("chi_present") or segment.get("magnetization_present"):
        lines.append("")
        lines.append(
            "  the segment also carries the static susceptibility / magnetisation tables"
        )
    return "\n".join(lines)


def render_magrelax(data: MagrelaxData, source: str, *, window_k: tuple[float, float] | None = None) -> str:
    """The tau(T) part of the menu-36 report."""
    lines: list[str] = []
    lines.append(f"Orca_Magrelax ({source})")
    lines.append(f"  input {data.magrelax_input}  derivatives {data.derivatives_file}"
                 f" ({data.n_derivatives} of size {data.derivative_size})"
                 f"  hessian {data.hessian_file}")
    if data.zeeman_eigenvalues_cm1:
        lines.append(
            "  Zeeman eigenvalues (cm-1): "
            + " ".join(f"{value:.1f}" for value in data.zeeman_eigenvalues_cm1)
        )
    lines.append(
        f"  vibrational modes listed by magrelax: {len(data.frequencies_cm1)}"
        + (" (none -- no one-phonon channel can be resonant)" if not data.frequencies_cm1 else "")
    )
    lines.append("")
    lines.append("  T (K)        tau (s)")
    for temperature, tau in zip(data.temperatures_k, data.tau_s):
        tau_text = f"{tau:.6e}" if np.isfinite(tau) else f"{tau}"
        lines.append(f"  {temperature:8.3f}   {tau_text}")
    positive = [(t, tau) for t, tau in zip(data.temperatures_k, data.tau_s) if np.isfinite(tau) and tau > 0]
    lines.append("")
    if len(positive) >= 3:
        try:
            fit = orbach_fit(data.temperatures_k, data.tau_s, window_k=window_k)
        except RelaxationError as error:
            lines.append(f"  Orbach fit: {error}")
        else:
            lines.append(
                f"  Orbach fit over {fit['n_points']} points: U_eff = {fit['ueff_k']:.1f} K"
                f" = {fit['ueff_cm1']:.1f} cm-1, tau_0 = {fit['tau0_s']:.3e} s,"
                f" R^2 = {fit['r_squared']:.4f}"
            )
    else:
        lines.append(
            "  Orbach fit: not applicable -- fewer than three finite positive tau points."
            " A single-mode toy system can give an all-zero rate table: the one-phonon"
            " channel needs a phonon energy matching a Zeeman gap (this run's modes: "
            + (", ".join(f"{f:.1f}" for f in data.frequencies_cm1) if data.frequencies_cm1 else "none")
            + " cm-1)."
        )
    return "\n".join(lines)


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the relaxation metrics and the measured format facts."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "The tunnel-gap/QTM framing (tau_QTM^-1 = 4 (Delta/hbar)^2 tau_m for "
                "non-Kramers pseudo-doublets; transverse dipolar fields as the Kramers "
                "QTM mechanism), the 15 deg non-collinearity working threshold and the "
                "empirical g_T*theta > 20 rule (g_T = (g1 + g2 + g3 sin theta_3)/3, "
                "from a survey of 20 Dy(III) SMMs) follow the guide's sections 7.1-7.2; "
                "its section 7.3 supplies the UBAR-based qualitative propagation and "
                "calls the underlying magnetic-perturbation model 'not realistic'."
            ),
            ref=(
                "Chilton N. F., Chem. Soc. Rev. 2025, 54, 11468-11487, "
                "DOI 10.1039/d5cs00493d"
            ),
            url="https://doi.org/10.1039/d5cs00493d",
            bibkey="chilton2025abinitio",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "ORCA 6.1.1 fixtures single_aniso/ and magrelax/: the DoSSC run prints "
                "two complete SINGLE_ANISO segments (SOC-only and SOC+SSC, D = 2.175287 "
                "vs 3.185964); MLTP is required for the g/D analysis (without it both "
                "the g-tensor and the barrier estimation read NOT INCLUDED); a Kramers "
                "doublet prints no ZFS/D block; the single-mode CO+ magrelax run has an "
                "all-zero rate table (no phonon matches its Zeeman gaps)."
            ),
            ref="ORCA 6.1.1, measured 2026-09-28 (fixtures single_aniso/, magrelax/)",
        ),
    )
