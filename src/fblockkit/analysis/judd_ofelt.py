"""Judd-Ofelt intensity parameters from a transition dataset (menu 34).

The standard Judd-Ofelt (JO) theory expresses an f-f electric-dipole line
strength as a linear combination of three phenomenological parameters,

    S_ED = sum_{lambda=2,4,6} Omega_lambda * |<Psi_1 || U^(lambda) || Psi_2>|^2,

with the doubly-reduced unit-tensor matrix elements U^(lambda) supplied per
transition (they are host-independent, ion-specific table values).  A
least-squares fit of experimental line strengths gives Omega_2, Omega_4,
Omega_6 -- the quantities behind Ln3+ absorption/emission intensities,
radiative rates, branching ratios and lifetimes.

Every convention below was pinned on the sources of this entry (the close reading is ``文献细读/细读_Judd-Ofelt_两篇_20260928.md``):

- **Conversion of a measured oscillator strength to a line strength**
  (Hovhannesyan/Boudon/Lepers, J. Lumin. 2022/2024):

      S_exp = 3 (2 J_1 + 1) hbar^2 / (2 m_e a_0^2 (E_2 - E_1))
              * (n_r / chi_ED) * f_exp,

  with S in atomic units (e^2 a_0^2), E in Hartree, and the virtual-cavity
  local-field correction **chi_ED = (n_r^2 + 2)^2 / 9**.  Both papers *print*
  chi_ED = (n_r^2 + 2)/9 -- a missing square; the square is what their
  reference implementation uses (``jo_so.f90``, GPL-3.0, read-only reference:
  ``calc_chi_ed(n,1) = (n**2+2)**2/9``) and what reproduces their published
  numbers (this module's regression test).  The mismatch is recorded here
  because it is exactly the kind of trap a fitter must not silently inherit.
- The same reference code evaluates the conversion through the equivalent
  form S = f * (2J_1+1) * lambda_cm * n_r / (chi_ED * S2gf) with
  S2gf = 1/(1.5 * 219474.63 cm^-1) (the module uses the explicit a.u. form).
- **Omega in cm^2 vs S in a.u.**:  Omega[cm^2] = Omega[a.u.] * a_0^2
  (a_0^2 = 2.8002851e-17 cm^2); equivalently S[a.u.] = e^2[esu^2] *
  Omega[cm^2] * U^2 / (e^2 a_0^2)[esu^2 cm^2].  The literature prints Omega
  in 10^-20 cm^2.
- **Refractive index**: a constant n_r, or the Sellmeier form
  n_r^2(lambda) = n0^2 + sum_i A_i lambda^2 / (lambda^2 - B_i) with B_i in
  micrometres squared (the reference code additionally cuts lambda at
  6000 nm; this module reproduces that cut).
- **Fit quality**: sigma = sqrt(sum (S_exp - S_ED)^2 / (N_tr - 3)) and the
  conventional relative value sigma / S_max (the papers print percentages).
- **Emission side**: A_ED = [e^2 a_0^2 (E_2-E_1)^3 / (3 pi eps0 hbar^4 c^3
  (2 J_2 + 1))] * n_r * chi_ED * S_ED; branching ratio beta = A_i / sum(A);
  radiative lifetime tau = 1 / sum(A).  The magnetic-dipole part of a
  transition is *not* computed here (it needs free-ion wave functions): a
  measured f_md may be given per transition and is subtracted before the
  fit, as the reference implementation does.

Measured anchors (regression, fixtures/judd_ofelt/): the Eu3+ dataset of
Babu et al. (Physica B 279, 262 (2000), SET B -- the very file the reference
code ships) with the U table of the 2024 paper: the fit reproduces the
published relative standard deviation (8.52 %) and Omega_6 (2.253e-20 cm^2)
digit for digit, and the published standard-JO ratios of the strong
transitions to better than 0.05.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..knowledge.models import EVIDENCE_LITERATURE, EVIDENCE_MEASURED, Evidence
from .plot_csv import csv_field, plot_number

__all__ = [
    "JudOError",
    "Host",
    "Transition",
    "Dataset",
    "FitResult",
    "read_dataset",
    "fit",
    "emission_rates",
    "jo_plot_csv",
    "render",
    "evidence",
]

#: cm^-1 per Hartree (CODATA; computed constants use the explicit a.u. formula)
_CM1_PER_HARTREE = 219474.6313705
#: Bohr radius in cm (CODATA 2018): 0.529177210903e-10 m
_A0_CM = 0.529177210903e-8
#: the reference implementation's cutoff on Sellmeier wavelengths (nm)
_SELLMEIER_CUT_NM = 6000.0
#: fine-structure constant and the atomic unit of frequency (s^-1); the a.u.
#: emission rate of an electric-dipole line is (4/3) alpha^3 (dE)^3 S/(2J+1),
#: numerically identical to the papers' SI prefactor
#: e^2 a0^2 / (3 pi eps0 hbar^4 c^3) to 0.02 %
_ALPHA = 7.2973525693e-3
_ATOMIC_FREQUENCY_S1 = 4.1341373335e16

_U_COLUMNS = ("u2", "u4", "u6")


class JudOError(ValueError):
    """The dataset cannot be fitted (with a next step)."""


@dataclass(frozen=True)
class Host:
    """The host material's refractive index: a constant or Sellmeier terms."""

    refractive_index: float | None = None
    sellmeier: tuple[tuple[float, float], ...] = ()  # (A_i, B_i[um^2]) plus n0
    sellmeier_n0: float = 1.0

    def n_at(self, wavelength_nm: float) -> float:
        """The refractive index at a wavelength (Sellmeier; cut at 6000 nm)."""
        if self.refractive_index is not None:
            return float(self.refractive_index)
        wl = min(float(wavelength_nm), _SELLMEIER_CUT_NM) / 1000.0  # micrometres
        n2 = self.sellmeier_n0**2
        for a_i, b_i in self.sellmeier:
            n2 += a_i * wl**2 / (wl**2 - b_i)
        if n2 <= 0.0:
            raise JudOError(
                f"the Sellmeier terms give n^2 = {n2:g} at {wavelength_nm:g} nm; "
                "next step: check the A/B coefficients and their units (B in um^2)."
            )
        return float(np.sqrt(n2))

    @staticmethod
    def chi_ed(n_refractive: float) -> float:
        """The virtual-cavity electric-dipole local field, (n^2+2)^2/9."""
        return (n_refractive**2 + 2.0) ** 2 / 9.0


@dataclass(frozen=True)
class Transition:
    """One f-f transition of the dataset (absorption unless labelled otherwise)."""

    label: str
    energy_cm1: float
    f_exp: float
    j_low: int
    u2: float = 0.0
    u4: float = 0.0
    u6: float = 0.0
    wavelength_nm: float | None = None
    f_md: float = 0.0
    j_up: int | None = None  # needed for the emission side only

    @property
    def wavelength_for_n(self) -> float:
        """The wavelength used by the Sellmeier lookup (measured one if given)."""
        return self.wavelength_nm if self.wavelength_nm is not None else 1e7 / self.energy_cm1


@dataclass(frozen=True)
class Dataset:
    """A JO dataset: the host, the fitted transitions, the emission block."""

    host: Host
    transitions: tuple[Transition, ...]
    emission: tuple[Transition, ...] = ()
    reference: str = ""


@dataclass(frozen=True)
class FitResult:
    """The outcome of one least-squares fit."""

    omega_au: tuple[float | None, float | None, float | None]
    omega_cm2: tuple[float | None, float | None, float | None]
    determined: tuple[bool, bool, bool]
    s_exp: tuple[float, ...]
    s_ed: tuple[float, ...]
    ratios: tuple[float, ...]
    sigma: float
    sigma_relative: float
    n_transitions: int
    weights: str
    host: Host
    transitions: tuple[Transition, ...]
    emission: tuple[Transition, ...] = field(default=())


def _require(value, label: str):
    if value is None:
        raise JudOError(f"the dataset is missing '{label}'.")
    return value


def read_dataset(path: str | Path) -> Dataset:
    """Read a JO dataset file (YAML; the schema is documented in the formats chapter).

    Required: a ``host`` (``refractive_index`` or ``sellmeier``) and a
    ``transitions`` list; each transition needs ``label``, ``energy_cm1`` (or
    ``wavelength_nm``), ``f``, ``j_low`` and the three ``u2``/``u4``/``u6``
    values.  Optional: ``f_md`` (subtracted before the ED fit), ``j_up``, an
    ``emission`` list (same fields; ``j_up`` required there).
    """
    import yaml

    source = Path(path)
    try:
        raw = yaml.safe_load(source.read_text(encoding="utf-8"))
    except OSError as exc:
        raise JudOError(f"the dataset file cannot be read: {exc}") from exc
    except yaml.YAMLError as exc:
        raise JudOError(f"the dataset file does not parse as YAML: {exc}") from exc
    if not isinstance(raw, dict):
        raise JudOError(
            "the dataset's top level is not a mapping. Next step: see the formats "
            "chapter for the Judd-Ofelt dataset schema."
        )
    host_raw = raw.get("host")
    if not isinstance(host_raw, dict):
        raise JudOError(
            "the dataset carries no 'host' mapping (refractive_index or sellmeier)."
        )
    sellmeier = host_raw.get("sellmeier")
    if sellmeier is not None:
        terms = sellmeier.get("terms") if isinstance(sellmeier, dict) else sellmeier
        try:
            host = Host(
                sellmeier=tuple((float(a), float(b)) for a, b in terms),
                sellmeier_n0=float(sellmeier.get("n0", 1.0)) if isinstance(sellmeier, dict) else 1.0,
            )
        except (TypeError, ValueError) as exc:
            raise JudOError(
                "the sellmeier block must carry 'n0' and 'terms' = [[A, B], ...] "
                f"with B in um^2 ({exc})."
            ) from exc
    elif "refractive_index" in host_raw:
        host = Host(refractive_index=float(host_raw["refractive_index"]))
    else:
        raise JudOError(
            "the host must give either 'refractive_index' or 'sellmeier'. Next step: "
            "state the host form and its values (a constant n, or Sellmeier A/B terms)."
        )

    def parse_transitions(entries, *, emission: bool) -> tuple[Transition, ...]:
        if entries is None:
            return ()
        if not isinstance(entries, list) or not entries:
            raise JudOError(
                f"the '{'emission' if emission else 'transitions'}' block must be a "
                "non-empty list."
            )
        parsed = []
        for index, entry in enumerate(entries):
            if not isinstance(entry, dict):
                raise JudOError(f"transition {index} is not a mapping.")
            label = str(entry.get("label", f"#{index}"))
            if "energy_cm1" in entry:
                energy = float(entry["energy_cm1"])
            elif "wavelength_nm" in entry:
                energy = 1e7 / float(entry["wavelength_nm"])
            else:
                raise JudOError(
                    f"transition '{label}' carries neither energy_cm1 nor "
                    "wavelength_nm. Next step: give the transition energy (cm^-1) "
                    "or the measured wavelength (nm)."
                )
            if energy <= 0.0:
                raise JudOError(f"transition '{label}' has a non-positive energy.")
            j_up = entry.get("j_up")
            if emission and j_up is None:
                raise JudOError(
                    f"emission transition '{label}' needs j_up (the emitting level's "
                    "J). Next step: add j_up."
                )
            parsed.append(
                Transition(
                    label=label,
                    energy_cm1=energy,
                    f_exp=float(entry.get("f", 0.0)),
                    j_low=int(_require(entry.get("j_low"), f"j_low of '{label}'")),
                    u2=float(entry.get("u2", 0.0)),
                    u4=float(entry.get("u4", 0.0)),
                    u6=float(entry.get("u6", 0.0)),
                    wavelength_nm=(
                        float(entry["wavelength_nm"]) if "wavelength_nm" in entry else None
                    ),
                    f_md=float(entry.get("f_md", 0.0)),
                    j_up=int(j_up) if j_up is not None else None,
                )
            )
        return tuple(parsed)

    transitions = parse_transitions(raw.get("transitions"), emission=False)
    emission = parse_transitions(raw.get("emission"), emission=True)
    if not transitions:
        raise JudOError(
            "the dataset carries no 'transitions' list. Next step: give the fitted "
            "transitions (label, energy, f, j_low and the U^(lambda) values)."
        )
    return Dataset(
        host=host,
        transitions=transitions,
        emission=emission,
        reference=str(raw.get("reference", "")),
    )


def _s_exp(transition: Transition, host: Host) -> float:
    """The measured oscillator strength as a line strength (atomic units)."""
    n = host.n_at(transition.wavelength_for_n)
    chi = Host.chi_ed(n)
    delta_e_au = transition.energy_cm1 / _CM1_PER_HARTREE
    return (
        3.0
        * (2 * transition.j_low + 1)
        * (0.5)  # hbar^2/(2 m_e a_0^2) = 1/2 Hartree
        / delta_e_au
        * (n / chi)
        * (transition.f_exp - transition.f_md)
    )


def fit(dataset: Dataset, *, weights: str = "none") -> FitResult:
    """Least-squares Omega_2/Omega_4/Omega_6 from the dataset's transitions.

    Unweighted by default (the papers' procedure); ``weights='1/S'`` is the
    so-called normalized variant some authors use -- recorded as a choice,
    because the two give different parameters from the same data.

    A parameter whose U^(lambda) column is zero over every fitted transition
    is *not determined* by the dataset; it is reported as such (never as a
    silent zero).
    """
    if weights not in ("none", "1/S"):
        raise JudOError(f"unknown weighting {weights!r}; use 'none' or '1/S'.")
    transitions = dataset.transitions
    n_tr = len(transitions)
    if n_tr < 4:
        raise JudOError(
            f"the dataset holds {n_tr} transition(s); a three-parameter fit needs "
            "at least four (the standard deviation divides by N-3). Next step: add "
            "transitions."
        )
    s_exp = np.array([_s_exp(t, dataset.host) for t in transitions])
    matrix = np.array([[t.u2, t.u4, t.u6] for t in transitions])
    column_scale = np.abs(matrix).max(axis=0)
    determined = tuple(bool(scale > 0.0) for scale in column_scale)
    if not any(determined):
        raise JudOError(
            "every transition carries U^(lambda) = 0 for all lambda, so no Omega "
            "can be fitted. Next step: give the U^(lambda) values (they are "
            "host-independent table values for the ion and transition)."
        )
    if weights == "1/S":
        with np.errstate(divide="ignore"):
            w = 1.0 / s_exp
        if not np.all(np.isfinite(w)):
            raise JudOError(
                "a transition carries zero experimental line strength; the 1/S "
                "weighting cannot use it. Next step: drop it or use the unweighted fit."
            )
        matrix_w = matrix * w[:, None]
        target_w = s_exp * w
        omega, *_ = np.linalg.lstsq(matrix_w, target_w, rcond=None)
    else:
        omega, *_ = np.linalg.lstsq(matrix, s_exp, rcond=None)
    omega = np.where(np.array(determined), omega, np.nan)
    s_ed = matrix @ np.nan_to_num(omega)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratios = np.where(s_exp != 0.0, s_ed / s_exp, np.nan)
    residual = s_exp - s_ed
    sigma = float(np.sqrt(np.sum(residual**2) / (n_tr - 3)))
    s_max = float(np.max(np.abs(s_exp)))
    omega_au = tuple(None if not d else float(v) for d, v in zip(determined, omega))
    omega_cm2 = tuple(
        None if v is None else float(v) * _A0_CM**2 for v in omega_au
    )
    return FitResult(
        weights=weights,
        omega_au=omega_au,
        omega_cm2=omega_cm2,
        determined=determined,
        s_exp=tuple(float(v) for v in s_exp),
        s_ed=tuple(float(v) for v in s_ed),
        ratios=tuple(float(v) for v in ratios),
        sigma=sigma,
        sigma_relative=sigma / s_max,
        n_transitions=n_tr,
        host=dataset.host,
        transitions=transitions,
        emission=dataset.emission,
    )


def emission_rates(result: FitResult) -> tuple[dict, ...]:
    """Per-transition A_ED, branching ratios and the radiative lifetime.

    Uses the fitted Omega over the dataset's ``emission`` block (each entry
    needs ``j_up``).  Magnetic-dipole rates are not computed (free-ion
    wave functions are needed); where an emission transition is known to be
    MD-dominated the caller should keep it out of this sum, as the sources do.
    """
    if not result.emission:
        return ()
    rows = []
    total = 0.0
    for transition in result.emission:
        s_ed = (
            transition.u2 * (result.omega_au[0] or 0.0)
            + transition.u4 * (result.omega_au[1] or 0.0)
            + transition.u6 * (result.omega_au[2] or 0.0)
        )
        n = result.host.n_at(transition.wavelength_for_n)
        chi = Host.chi_ed(n)
        delta_e_au = transition.energy_cm1 / _CM1_PER_HARTREE
        j_up = transition.j_up or 0
        rate = (
            (4.0 / 3.0)
            * _ALPHA**3
            * delta_e_au**3
            * n
            * chi
            * s_ed
            / (2 * j_up + 1)
        ) * _ATOMIC_FREQUENCY_S1
        rows.append(
            {
                "label": transition.label,
                "energy_cm1": transition.energy_cm1,
                "s_ed_au": s_ed,
                "a_ed_s1": rate,
            }
        )
        total += rate
    if total <= 0.0:
        raise JudOError(
            "the fitted parameters give a non-positive total emission rate; next "
            "step: check the emission block's transitions and the fit."
        )
    for row in rows:
        row["branching"] = row["a_ed_s1"] / total
    lifetime_us = 1e6 / total
    return tuple(rows) + ({"total_a_s1": total, "lifetime_us": lifetime_us},)


def render(result: FitResult, *, source: str = "") -> str:
    """The menu's report: the fit, the per-transition table and the boundaries."""
    weighting = (
        "unweighted" if result.weights == "none" else f"weighted ({result.weights})"
    )
    out = [f"Judd-Ofelt intensity parameters (standard JO fit, S-space, {weighting})"]
    if source:
        out.append(f"  dataset: {source}")
    host = result.host
    if host.refractive_index is not None:
        n_text = f"constant n = {host.refractive_index:g}"
    else:
        terms = ", ".join(f"A={a:g} B={b:g} um^2" for a, b in host.sellmeier)
        n_text = f"Sellmeier n0 = {host.sellmeier_n0:g} ({terms})"
    out.append(f"  host: {n_text}; chi_ED = (n^2+2)^2/9 (virtual cavity)")
    omega_parts = []
    for name, value, ok in zip(("Omega_2", "Omega_4", "Omega_6"), result.omega_cm2, result.determined):
        if ok and value is not None:
            omega_parts.append(f"{name} = {value / 1e-20:.3f}e-20 cm^2")
        else:
            omega_parts.append(f"{name} = not determined by this dataset (all U = 0)")
    out.append("  " + "; ".join(omega_parts))
    out.append(
        f"  fit: {result.n_transitions} transitions, sigma = {result.sigma:.4e} a.u., "
        f"sigma/S_max = {result.sigma_relative * 100:.2f} %"
    )
    out.append("")
    out.append("  transition            energy (cm-1)      f_exp      S_exp (a.u.)    S_ED/S_exp")
    for transition, s_exp, ratio in zip(result.transitions, result.s_exp, result.ratios):
        f_note = f"{transition.f_exp:>12.4e}"
        if transition.f_md:
            f_note += f" (-md {transition.f_md:.1e})"
        out.append(
            f"  {transition.label:<20s} {transition.energy_cm1:>14.1f} {f_note:>16s} "
            f"{s_exp:>14.4e} {ratio:>13.3f}"
        )
    boundaries = [
        "",
        "boundaries (measured and documented in the close-reading record):",
        "  - the local field is the squared virtual-cavity form (n^2+2)^2/9; the "
        "source papers print it without the square -- the square is what their "
        "reference code and their published numbers use (regression-tested here);",
        "  - the U^(lambda) values are host-independent table values for the ion "
        "and must come from the same source as the fit (transcribe them together "
        "with the dataset and check them against the published fit);",
        "  - magnetic-dipole contributions are not computed: give f_md per "
        "transition where it is significant (the fit subtracts it);",
        "  - the extended (perturbative X_k) JO models of the 2022/2024 papers "
        "need free-ion atomic-structure wave functions and are outside this "
        "toolkit's post-processing scope; this menu implements the standard JO fit.",
    ]
    out += boundaries
    return "\n".join(out)


def jo_plot_csv(result: FitResult) -> str:
    """The fitted transition table as a plot-ready CSV companion.

    One row per fitted transition (``side`` = absorption: the label, the
    energy, the measured oscillator strength and the measured/fitted line
    strengths) followed by the emission block's per-transition rows when the
    dataset carries one (``side`` = emission: label, energy, A_ED and the
    branching ratio); columns that do not apply to a side are left empty.
    The derived summaries (Omega, sigma, the radiative lifetime) stay in the
    report -- the companion carries the per-transition numbers the report's
    tables print (the formats chapter, "Plot-ready CSV companions").
    """
    lines = ["side,label,energy_cm1,f_exp,s_exp,s_ed,a_ed_s1,branching"]
    for transition, s_exp, s_ed in zip(result.transitions, result.s_exp, result.s_ed):
        lines.append(
            f"absorption,{csv_field(transition.label)},"
            f"{plot_number(transition.energy_cm1)},{plot_number(transition.f_exp)},"
            f"{plot_number(s_exp)},{plot_number(s_ed)},,"
        )
    # emission_rates appends a trailing total row (the report's summary line);
    # the companion carries only the per-transition rows
    for row in emission_rates(result)[:-1]:
        lines.append(
            f"emission,{csv_field(row['label'])},{plot_number(row['energy_cm1'])},"
            f",,,{plot_number(row['a_ed_s1'])},{plot_number(row['branching'])}"
        )
    return "\n".join(lines) + "\n"


def render_emission(rows: tuple[dict, ...]) -> str:
    """The optional emission-side table (A_ED, branching ratios, lifetime)."""
    if not rows:
        return ""
    out = ["", "emission side (from the fitted Omega over the emission block):"]
    out.append("  transition            energy (cm-1)    A_ED (s^-1)   branching")
    for row in rows[:-1]:
        out.append(
            f"  {row['label']:<20s} {row['energy_cm1']:>14.1f} "
            f"{row['a_ed_s1']:>13.4g} {row['branching']:>11.4f}"
        )
    total = rows[-1]
    out.append(
        f"  total A_ED = {total['total_a_s1']:.4g} s^-1; radiative lifetime = "
        f"{total['lifetime_us']:.4g} us (magnetic-dipole channels not included)"
    )
    return "\n".join(out)


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the formalism, the conventions and the regression."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "The standard Judd-Ofelt equations used here (ED line strength as "
                "sum of Omega_lambda * |<U^(lambda)>|^2; the measured-to-S "
                "conversion with the virtual-cavity local field; sigma; the "
                "Einstein A_ED expression) are those of Hovhannesyan, Boudon and "
                "Lepers, including their wavelength-dependent Sellmeier host; the "
                "U^(lambda) table of the regression dataset is that paper's table, "
                "cross-checked there against Carnall et al. (1968)."
            ),
            ref=(
                "Hovhannesyan G., Boudon V., Lepers M., J. Lumin. 2024, 266, 120234, "
                "DOI 10.1016/j.jlumin.2023.120234 (extension; U table, sigma, "
                "Sellmeier) and J. Lumin. 2022, 241, 118456, "
                "DOI 10.1016/j.jlumin.2021.118456 (paper I, r_0 table)"
            ),
            url="https://doi.org/10.1016/j.jlumin.2023.120234",
            bibkey="hovhannesyan2024extension",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "Measured on the authors' own reference implementation and dataset "
                "(jo_so, GPL-3.0, read-only; gitlab.com/labicb/joso): its "
                "calc_chi_ed is (n^2+2)^2/9 (the papers print the un-squared form); "
                "its Eu3+ Babu-2000 SET B dataset and U table reproduce the "
                "published relative sigma (8.52 %) and Omega_6 (2.253e-20 cm^2) "
                "digit for digit under this module's equations."
            ),
            ref="tests/test_judd_ofelt.py; fixtures/judd_ofelt/",
        ),
    )
