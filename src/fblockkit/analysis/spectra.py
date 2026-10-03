"""Spectra rendering: stick tables in, broadened curves and tidy CSVs out (W0).

A spectral calculation ends in a stick table -- transition energies and
intensities (menu 31's per-state table, menu 34's Judd-Ofelt table, menu 37's
ROCIS blocks).  Comparing that table with experiment, or with another
calculation, needs one shared step: broaden the sticks with a line shape and
hand the curve to a plotting tool.  This module is that step's one home in
code.  It renders no pictures (Multiwfn, orca_plot or the caller's plotter do
that); it produces the curve's numbers:

- unit-area line shapes parameterised by FWHM: :func:`gaussian`,
  :func:`lorentzian`, :func:`voigt` (the Voigt profile is evaluated from the
  Faddeeva function, not from an empirical fit -- see :func:`voigt`);
- :func:`broaden` -- a stick table summed into a curve on an energy grid;
- :func:`boltzmann_weights` -- thermal populations at a temperature;
- :func:`weighted_spectrum` -- weighting and broadening in one call (thermal
  populations and/or an energy power such as the emission nu^3 factor);
- :func:`write_spectrum_csv` -- the curve as a tidy CSV file, in the
  plot-ready companion style of :mod:`fblockkit.analysis.plot_csv`.

Energy units are the caller's; pass cm^-1 (the f-block habit) or eV, as long
as the grid, the stick positions and the FWHM share one unit.  Every profile
is normalised to unit area, so a broadened curve's area equals the sum of the
stick intensities; broadening never renormalises by itself
(:func:`weighted_spectrum` renormalises only with ``normalize=True``).
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .plot_csv import plot_number

__all__ = [
    "boltzmann_weights",
    "broaden",
    "gaussian",
    "lorentzian",
    "voigt",
    "weighted_spectrum",
    "write_spectrum_csv",
]

# Exact SI constants (2019 redefinition); k_B / (h c) converts cm^-1 to kelvin.
_K_B_J_PER_K = 1.380649e-23
_H_J_S = 6.62607015e-34
_C_M_PER_S = 299792458.0
_CM1_PER_K = _K_B_J_PER_K / (_H_J_S * _C_M_PER_S) / 100.0  # 0.6950348 cm^-1 / K

# The Faddeeva evaluation switches from the power series (which covers the
# small-|z| region with mild cancellation) to the asymptotic series at
# |z| = 4.2; that crossover is where the two measured error curves cross
# (both sit below 3e-8 of the profile's peak height there).
_FADDEEVA_SERIES_LIMIT = 4.2
_FADDEEVA_SERIES_TERMS = 160
_FADDEEVA_ASYMPTOTIC_MAX_TERMS = 80

# Power-series coefficients 1 / Gamma(1 + n/2), n = 1 .. _FADDEEVA_SERIES_TERMS.
_FADDEEVA_COEFFICIENTS = tuple(
    1.0 / math.gamma(1.0 + 0.5 * n) for n in range(1, _FADDEEVA_SERIES_TERMS + 1)
)

_SHAPES = ("gaussian", "lorentzian", "voigt")


# --- line shapes -----------------------------------------------------------------


def gaussian(x, center: float = 0.0, fwhm: float = 1.0) -> np.ndarray:
    """The unit-area Gaussian (normal) profile of the given FWHM.

    ``exp(-(x-center)^2 / (2 sigma^2)) / (sigma sqrt(2 pi))`` with
    ``sigma = fwhm / (2 sqrt(2 ln 2))``, so the full width at half maximum is
    exactly ``fwhm`` and the integral over the real line is 1.
    """
    center = _finite_center(center)
    fwhm = _positive_fwhm(fwhm, shape="gaussian")
    sigma = fwhm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    values = np.asarray(x, dtype=float)
    scaled = (values - center) / sigma
    return np.exp(-0.5 * scaled * scaled) / (sigma * math.sqrt(2.0 * math.pi))


def lorentzian(x, center: float = 0.0, fwhm: float = 1.0) -> np.ndarray:
    """The unit-area Lorentzian (Cauchy) profile of the given FWHM.

    ``(gamma / pi) / ((x-center)^2 + gamma^2)`` with ``gamma = fwhm / 2``, so
    the full width at half maximum is exactly ``fwhm`` and the integral over
    the real line is 1.
    """
    center = _finite_center(center)
    fwhm = _positive_fwhm(fwhm, shape="lorentzian")
    half = fwhm / 2.0
    values = np.asarray(x, dtype=float)
    return (half / math.pi) / ((values - center) ** 2 + half * half)


def voigt(x, center: float = 0.0, fwhm_g: float = 1.0, fwhm_l: float = 0.0) -> np.ndarray:
    """The unit-area Voigt profile: a Gaussian of FWHM ``fwhm_g`` convoluted
    with a Lorentzian of FWHM ``fwhm_l`` (``fwhm_l = 0`` gives the pure
    Gaussian, ``fwhm_g = 0`` the pure Lorentzian, both returned analytically).

    The profile is evaluated from the Faddeeva function, not from an empirical
    pseudo-Voigt fit: ``V(x) = Re[w(z)] / (sqrt(2 pi) sigma)`` with ``z =
    ((x - center) + i gamma) / (sqrt(2) sigma)``, where ``sigma`` and
    ``gamma`` are the two components' half-width parameters and ``w(z) =
    exp(-z^2) erfc(-i z)`` (Abramowitz & Stegun, *Handbook of Mathematical
    Functions*, Ch. 7).  ``w`` is evaluated with two coefficient-free standard
    expansions, the power series ``sum_n (i z)^n / Gamma(1 + n/2)`` for
    ``|z| <= 4.2`` and the asymptotic series ``(i/sqrt(pi)) sum_m
    (2m-1)!!/(2^m z^(2m+1))`` truncated at its smallest term above it; the
    profile departs from the true Voigt function by less than 3e-8 of its own
    peak height everywhere (the crossover is where the two measured error
    curves meet; verified against ``scipy.special.wofz`` over a wide (x, y)
    sweep during development -- scipy is not a dependency of this toolkit, and
    no fitted line-shape constants are used).  Accuracy is irrelevant to the
    pure limits, which are returned analytically.
    """
    center = _finite_center(center)
    fwhm_g = _nonnegative_fwhm(fwhm_g, shape="voigt (Gaussian component)")
    fwhm_l = _nonnegative_fwhm(fwhm_l, shape="voigt (Lorentzian component)")
    if fwhm_g == 0.0 and fwhm_l == 0.0:
        raise ValueError(
            "a Voigt profile needs at least one non-zero width, got "
            "fwhm_g = 0 and fwhm_l = 0. Next step: pass the two component "
            "FWHMs, e.g. voigt(x, center, fwhm_g=300.0, fwhm_l=50.0)."
        )
    if fwhm_l == 0.0:
        return gaussian(x, center, fwhm_g)
    if fwhm_g == 0.0:
        return lorentzian(x, center, fwhm_l)
    sigma = fwhm_g / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    gamma = fwhm_l / 2.0
    values = np.asarray(x, dtype=float)
    z = ((values - center) + 1j * gamma) / (math.sqrt(2.0) * sigma)
    return _faddeeva_w(z).real / (math.sqrt(2.0 * math.pi) * sigma)


def _faddeeva_w(z: np.ndarray) -> np.ndarray:
    """The Faddeeva function ``w(z) = exp(-z^2) erfc(-i z)`` on ``Im z >= 0``.

    Two coefficient-free standard expansions (Abramowitz & Stegun, Ch. 7):
    the power series below ``|z| = 4.2`` and the asymptotic series above it.
    Measured against ``scipy.special.wofz`` over a wide sweep, both stay
    within 3e-8 of the profile's peak scale (``|w(i y)|``), the crossover
    sitting where their error curves meet.
    """
    z = np.atleast_1d(np.asarray(z, dtype=complex))
    w = np.empty(z.shape, dtype=complex)
    small = np.abs(z) <= _FADDEEVA_SERIES_LIMIT
    if small.any():
        w[small] = _faddeeva_series(z[small])
    if (~small).any():
        w[~small] = _faddeeva_asymptotic(z[~small])
    return w


def _faddeeva_series(z: np.ndarray) -> np.ndarray:
    """``sum_n (i z)^n / Gamma(1 + n/2)`` to ``_FADDEEVA_SERIES_TERMS`` terms.

    The series converges everywhere; the truncation point is chosen so that
    the discarded tail is negligible up to the region's ``|z| = 4.2`` edge,
    and the residual error is the mild cancellation of that size (its terms
    peak near ``e^(|z|^2)`` while the sum itself is ``O(1)``).
    """
    iz = 1j * z
    total = np.ones_like(z)
    power = np.ones_like(z)
    for coefficient in _FADDEEVA_COEFFICIENTS:
        power = power * iz
        total = total + coefficient * power
    return total


def _faddeeva_asymptotic(z: np.ndarray) -> np.ndarray:
    """``(i/sqrt(pi)) sum_m (2m-1)!! / (2^m z^(2m+1))``, truncated at its
    smallest term (per element, so a mixed batch stops each argument where
    that argument's series stops improving).

    The asymptotic series diverges once ``m`` passes ``|z|^2``; stopping at
    the smallest term leaves a remainder of order ``exp(-|z|^2)``, which is
    below 3e-8 for every ``|z|`` at and above the crossover (and negligible
    for the far arguments a spectrum grid actually carries).
    """
    z = np.asarray(z, dtype=complex)
    total = np.zeros_like(z)
    term = 1.0 / z
    previous = np.full(z.shape, np.inf)
    active = np.ones(z.shape, dtype=bool)
    for m in range(_FADDEEVA_ASYMPTOTIC_MAX_TERMS):
        magnitude = np.abs(term)
        active = active & (magnitude <= previous)
        total = total + np.where(active, term, 0.0)
        previous = np.where(active, magnitude, previous)
        term = term * (2.0 * m + 1.0) / (2.0 * z * z)
    return (1j / math.sqrt(math.pi)) * total


# --- stick tables -> curves ------------------------------------------------------


def broaden(
    sticks,
    grid,
    *,
    shape: str = "gaussian",
    fwhm,
) -> np.ndarray:
    """Sum a stick table into a broadened curve on an energy grid.

    ``sticks`` is a sequence of ``(energy, intensity)`` pairs (or an ``(n, 2)``
    array); each stick contributes ``intensity`` times a unit-area profile of
    the given ``shape`` centred at its energy, so the curve's area equals the
    sum of the intensities.  ``shape`` is ``"gaussian"``, ``"lorentzian"`` or
    ``"voigt"``; ``fwhm`` is one positive width for the first two and the pair
    ``(fwhm_g, fwhm_l)`` for ``"voigt"``, in the grid's energy unit.  The grid
    must be strictly increasing.  The cost is one profile evaluation per stick
    and per grid point.
    """
    shape_name, profile, widths = _line_shape(shape, fwhm)
    grid_values = _monotonic_grid(grid)
    out = np.zeros(grid_values.shape, dtype=float)
    for energy, intensity in _stick_rows(sticks):
        out = out + intensity * profile(grid_values, energy, *widths)
    return out


def weighted_spectrum(
    sticks,
    grid,
    *,
    shape: str = "gaussian",
    fwhm,
    weights=None,
    temperature_k: float | None = None,
    degeneracies=None,
    energy_exponent: float = 0.0,
    normalize: bool = False,
) -> np.ndarray:
    """Weight a stick table, then broaden it -- the W1/W2 one-call combiner.

    The intensity of each stick is scaled by

    - a population factor: either explicit ``weights`` (one per stick, signed
      if the caller wants a difference spectrum) or the Boltzmann populations
      of the stick energies at ``temperature_k`` with the given
      ``degeneracies`` (the sticks are then read as transitions from a
      thermally populated manifold to one common lower level -- the emission
      case); and
    - ``energy ** energy_exponent``, the energy weighting of a transition rate
      (``energy_exponent=3`` is the emission nu^3 factor); this needs positive
      stick energies.

    ``normalize=True`` divides the curve by its integral over the grid (the
    area outside the grid is not recovered), for band-shape comparison with a
    normalised experimental spectrum.  Every input is validated here; the
    error messages say what to pass instead.
    """
    shape_name, profile, widths = _line_shape(shape, fwhm)
    grid_values = _monotonic_grid(grid)
    rows = _stick_rows(sticks)
    energies = np.array([energy for energy, _intensity in rows], dtype=float)
    intensities = np.array([intensity for _energy, intensity in rows], dtype=float)

    if weights is not None and temperature_k is not None:
        raise ValueError(
            "weighted_spectrum takes either explicit weights or a temperature, "
            "not both. Next step: pass weights=... for measured/derived "
            "populations, or temperature_k=... to compute Boltzmann "
            "populations from the stick energies."
        )
    if degeneracies is not None and temperature_k is None:
        raise ValueError(
            "degeneracies only enter through the Boltzmann populations. "
            "Next step: pass temperature_k=... together with degeneracies, or "
            "fold the degeneracies into weights=... yourself."
        )
    if weights is not None:
        population = _numbers(
            weights,
            what="weights",
            next_step="give one weight per stick (the same order as sticks).",
        )
        if population.shape != intensities.shape:
            raise ValueError(
                f"weights has {population.size} entries for {intensities.size} "
                f"sticks. Next step: give one weight per stick (the same order "
                f"as sticks)."
            )
        if not np.all(np.isfinite(population)):
            raise ValueError(
                "weights contains a non-finite entry. Next step: drop the "
                "non-finite sticks before weighting, or fix the population "
                "table."
            )
    elif temperature_k is not None:
        population = boltzmann_weights(energies, temperature_k, degeneracies)
    else:
        population = np.ones_like(intensities)

    try:
        exponent = float(energy_exponent)
    except (TypeError, ValueError):
        raise ValueError(
            f"energy_exponent must be a number, got {energy_exponent!r}. Next "
            f"step: pass 0 for no energy weighting or 3 for the emission nu^3 "
            f"factor."
        ) from None
    if not math.isfinite(exponent):
        raise ValueError(
            "energy_exponent must be a finite number. Next step: pass 0 for "
            "no energy weighting or 3 for the emission nu^3 factor."
        )
    if exponent != 0.0:
        if np.any(energies <= 0.0):
            raise ValueError(
                "an energy power needs positive stick energies, and at least "
                "one stick is at or below zero. Next step: reference the "
                "energies to a common lower level, or drop energy_exponent."
            )
        population = population * energies**exponent

    out = np.zeros(grid_values.shape, dtype=float)
    for energy, intensity, factor in zip(energies, intensities, population):
        out = out + (intensity * factor) * profile(grid_values, float(energy), *widths)
    if normalize:
        area = _trapezoid(grid_values, out)
        if area > 0.0:
            out = out / area
    return out


def boltzmann_weights(
    energies_cm1,
    temperature_k: float,
    degeneracies=None,
) -> np.ndarray:
    """Normalised Boltzmann populations of levels given in cm^-1.

    ``w_i = g_i exp(-(E_i - E_min) / k_B T)`` normalised to sum 1, with
    ``k_B = 0.6950348 cm^-1 / K`` (the exact SI constants are combined in this
    module).  Only energy differences matter, so the energies may sit on any
    zero; the minimum is subtracted before exponentiating so that a table far
    from zero cannot underflow to 0/0.  ``degeneracies`` defaults to one per
    level; pass the multiplicities to weight a degenerate level by its
    degeneracy (at high temperature the weights then tend to
    ``g_i / sum(g)``).  ``temperature_k`` must be positive: at T = 0 every
    level but the ground one is empty and the population is not a limit of
    this formula.
    """
    energies = _numbers(
        energies_cm1,
        what="energies_cm1",
        next_step="pass the levels' energies in cm^-1, one per level.",
    )
    if energies.ndim != 1 or energies.size == 0:
        raise ValueError(
            "energies_cm1 must be a non-empty one-dimensional sequence of "
            "level energies. Next step: pass the levels' energies in cm^-1, "
            "one per level."
        )
    if not np.all(np.isfinite(energies)):
        raise ValueError(
            "energies_cm1 contains a non-finite entry. Next step: drop the "
            "non-finite levels before weighting, or fix the level table."
        )
    temperature = float(temperature_k)
    if not math.isfinite(temperature) or temperature <= 0.0:
        raise ValueError(
            f"temperature_k must be positive and finite, got {temperature_k!r}. "
            f"Next step: pass the temperature in kelvin (for example 300.0); "
            f"at T <= 0 the ground level alone is populated and no Boltzmann "
            f"population exists."
        )
    if degeneracies is None:
        g = np.ones_like(energies)
    else:
        g = _numbers(
            degeneracies,
            what="degeneracies",
            next_step="give one degeneracy per level (the same order as energies_cm1).",
        )
        if g.shape != energies.shape:
            raise ValueError(
                f"degeneracies has {g.size} entries for {energies.size} "
                f"levels. Next step: give one degeneracy per level (the same "
                f"order as energies_cm1)."
            )
        if not np.all(np.isfinite(g)) or np.any(g <= 0.0):
            raise ValueError(
                "degeneracies must be positive and finite. Next step: give "
                "each level its multiplicity (2J + 1 or 2S + 1 as appropriate)."
            )
    shifted = (energies - energies.min()) * (_CM1_PER_K / temperature)
    weights = g * np.exp(-shifted)
    total = weights.sum()
    if total <= 0.0:  # unreachable with g > 0 and a zero (ground) shift
        raise ValueError(  # pragma: no cover
            "the Boltzmann weights underflowed to zero. Next step: check the "
            "temperature and the level energies."
        )
    return weights / total


# --- the tidy CSV ------------------------------------------------------------------


def write_spectrum_csv(
    path,
    grid,
    intensity,
    meta: Mapping[str, Any] | None = None,
) -> Path:
    """Write a curve as a tidy CSV file (the plot-ready companion style).

    One header line with underscore-separated column names, one data row per
    grid point at ten significant digits (:func:`plot_csv.plot_number`), and
    nothing else unless ``meta`` is given; a non-finite intensity is written
    as its token (``inf``, ``-inf``, ``nan``), as in the menu companions.
    ``meta`` carries the context the menu companions leave to their report
    (the line shape, the FWHM, the temperature, the source file): two keys
    are reserved -- ``x_column`` and ``y_column`` name the columns (defaults
    ``energy_cm1`` and ``intensity``) -- and every other entry is written as a
    ``# key: value`` comment line above the header, in the mapping's order.
    There is no timestamp and no version stamp, so the same input gives a
    byte-identical file; comment lines are skipped by pandas, numpy.loadtxt
    and gnuplot, so the file stays readable by every plotting tool.  The
    parent directory must already exist (the file lands beside the report or
    under the case directory that cites it).
    """
    target = Path(path)
    target.write_text(_spectrum_csv_text(grid, intensity, meta), encoding="utf-8")
    return target


def _spectrum_csv_text(
    grid,
    intensity,
    meta: Mapping[str, Any] | None,
) -> str:
    grid_values = _monotonic_grid(grid)
    intensities = _numbers(
        intensity,
        what="intensity",
        next_step="pass one intensity per grid point (the output of "
        "broaden/weighted_spectrum fits directly).",
    )
    if intensities.ndim != 1 or intensities.shape != grid_values.shape:
        raise ValueError(
            f"intensity has shape {intensities.shape} for a grid of "
            f"{grid_values.size} points. Next step: pass one intensity per "
            f"grid point (the output of broaden/weighted_spectrum fits "
            f"directly)."
        )
    x_column, y_column, comments = _split_meta(meta)
    lines = [f"# {comment}" for comment in comments]
    lines.append(f"{x_column},{y_column}")
    for energy, value in zip(grid_values, intensities):
        lines.append(f"{plot_number(energy)},{plot_number(value)}")
    return "\n".join(lines) + "\n"


def _split_meta(meta: Mapping[str, Any] | None) -> tuple[str, str, list[str]]:
    """Reserved column names and the rendered comment lines of a meta mapping."""
    x_column, y_column = "energy_cm1", "intensity"
    comments: list[str] = []
    if meta is None:
        return x_column, y_column, comments
    if not isinstance(meta, Mapping):
        raise ValueError(
            f"meta must be a mapping of context entries, got {type(meta).__name__}. "
            f"Next step: pass for example "
            f'meta={{"shape": "gaussian", "fwhm_cm1": 100.0}}.'
        )
    for key, value in meta.items():
        name = str(key)
        if "\n" in name or "\r" in name:
            raise ValueError(
                "a meta key contains a line break. Next step: use a one-line "
                "key (it becomes the comment's label)."
            )
        if name == "x_column":
            x_column = _column_name(value)
            continue
        if name == "y_column":
            y_column = _column_name(value)
            continue
        if isinstance(value, bool):
            rendered = str(value)
        elif isinstance(value, (int, float)):
            rendered = plot_number(value)
        elif isinstance(value, str):
            rendered = value
        elif value is None:
            raise ValueError(
                f"meta[{name!r}] is None. Next step: drop the key instead of "
                f"passing an empty value (the comment lines are the file's "
                f"only context)."
            )
        else:
            raise ValueError(
                f"meta[{name!r}] has type {type(value).__name__}; only "
                f"numbers, booleans and strings can be written as comments. "
                f"Next step: format the value yourself (for example a list as "
                f"a comma-joined string)."
            )
        if "\n" in rendered or "\r" in rendered:
            raise ValueError(
                f"meta[{name!r}] contains a line break. Next step: keep each "
                f"comment to one line."
            )
        comments.append(f"{name}: {rendered}")
    return x_column, y_column, comments


def _column_name(value: Any) -> str:
    if not isinstance(value, str) or not value or "," in value:
        raise ValueError(
            f"a column name must be a non-empty string without commas, got "
            f"{value!r}. Next step: pass a name like 'energy_eV' or "
            f"'intensity'."
        )
    if "\n" in value or "\r" in value:
        raise ValueError(
            "a column name contains a line break. Next step: keep the name to "
            "one line."
        )
    return value


# --- validation helpers ------------------------------------------------------------


def _numbers(values, *, what: str, next_step: str) -> np.ndarray:
    """One-dimensional float array from a user table, with a Next-step message."""
    try:
        return np.atleast_1d(np.asarray(values, dtype=float))
    except (TypeError, ValueError):
        raise ValueError(
            f"{what} must be a sequence of numbers. Next step: {next_step}"
        ) from None


def _finite_center(center: Any) -> float:
    try:
        value = float(center)
    except (TypeError, ValueError):
        raise ValueError(
            f"center must be a number, got {center!r}. Next step: pass the "
            f"profile centre in the grid's energy unit."
        ) from None
    if not math.isfinite(value):
        raise ValueError(
            f"center must be finite, got {value!r}. Next step: pass the "
            f"profile centre in the grid's energy unit."
        )
    return value


def _positive_fwhm(fwhm: Any, *, shape: str) -> float:
    value = _finite_fwhm(fwhm, shape=shape)
    if value <= 0.0:
        raise ValueError(
            f"the {shape} FWHM must be positive, got {fwhm!r}; a zero-width "
            f"profile is a delta function and has no curve on a grid. Next "
            f"step: pass the experimental or chosen line width in the grid's "
            f"energy unit."
        )
    return value


def _nonnegative_fwhm(fwhm: Any, *, shape: str) -> float:
    value = _finite_fwhm(fwhm, shape=shape)
    if value < 0.0:
        raise ValueError(
            f"the {shape} FWHM must be non-negative, got {fwhm!r}. Next step: "
            f"pass the width in the grid's energy unit (0 is allowed here and "
            f"degenerates the profile)."
        )
    return value


def _finite_fwhm(fwhm: Any, *, shape: str) -> float:
    try:
        value = float(fwhm)
    except (TypeError, ValueError):
        raise ValueError(
            f"the {shape} FWHM must be a number, got {fwhm!r}. Next step: "
            f"pass the width in the grid's energy unit."
        ) from None
    if not math.isfinite(value):
        raise ValueError(
            f"the {shape} FWHM must be finite, got {value!r}. Next step: pass "
            f"the width in the grid's energy unit."
        )
    return value


def _line_shape(shape: Any, fwhm: Any):
    """The profile function and its widths for one (shape, fwhm) pair."""
    shape_name = str(shape).strip().lower()
    if shape_name not in _SHAPES:
        raise ValueError(
            f"unknown line shape {shape!r}; the shapes are "
            f"{', '.join(repr(name) for name in _SHAPES)}. Next step: pass "
            f"shape='gaussian', 'lorentzian' or 'voigt'."
        )
    if shape_name == "gaussian":
        return shape_name, gaussian, (_positive_fwhm(fwhm, shape="gaussian"),)
    if shape_name == "lorentzian":
        return shape_name, lorentzian, (_positive_fwhm(fwhm, shape="lorentzian"),)
    fwhm_g, fwhm_l = _voigt_widths(fwhm)
    return shape_name, voigt, (fwhm_g, fwhm_l)


def _voigt_widths(fwhm: Any) -> tuple[float, float]:
    if isinstance(fwhm, (str, bytes)) or np.ndim(fwhm) != 1:
        raise ValueError(
            f"shape='voigt' needs the two component widths, got {fwhm!r}. "
            f"Next step: pass fwhm=(fwhm_g, fwhm_l), for example (300.0, 50.0) "
            f"in cm^-1."
        )
    widths = list(fwhm)
    if len(widths) != 2:
        raise ValueError(
            f"shape='voigt' needs exactly two widths (fwhm_g, fwhm_l), got "
            f"{len(widths)}. Next step: pass the Gaussian and the Lorentzian "
            f"component FWHMs as a pair."
        )
    fwhm_g = _nonnegative_fwhm(widths[0], shape="voigt (Gaussian component)")
    fwhm_l = _nonnegative_fwhm(widths[1], shape="voigt (Lorentzian component)")
    if fwhm_g == 0.0 and fwhm_l == 0.0:
        raise ValueError(
            "shape='voigt' needs at least one non-zero component width, got "
            "(0, 0). Next step: pass the Gaussian and the Lorentzian FWHMs, "
            "e.g. fwhm=(300.0, 50.0)."
        )
    return fwhm_g, fwhm_l


def _monotonic_grid(grid) -> np.ndarray:
    try:
        values = np.atleast_1d(np.asarray(list(grid), dtype=float))
    except (TypeError, ValueError):
        raise ValueError(
            "grid must be a sequence of energies. Next step: pass the energy "
            "grid the curve is evaluated on, for example "
            "numpy.linspace(lo, hi, n)."
        ) from None
    if values.ndim != 1 or values.size == 0:
        raise ValueError(
            "grid must be a non-empty one-dimensional sequence of energies. "
            "Next step: pass the energy grid the curve is evaluated on."
        )
    if not np.all(np.isfinite(values)):
        raise ValueError(
            "grid contains a non-finite entry. Next step: build the grid from "
            "finite energies (for example numpy.linspace(lo, hi, n))."
        )
    if values.size > 1 and not np.all(np.diff(values) > 0.0):
        raise ValueError(
            "grid must be strictly increasing. Next step: sort the grid and "
            "drop duplicate energies before broadening."
        )
    return values


def _stick_rows(sticks) -> list[tuple[float, float]]:
    rows: list[tuple[float, float]] = []
    try:
        arr = np.asarray(list(sticks), dtype=float)
    except (TypeError, ValueError):
        raise ValueError(
            "sticks must be a sequence of (energy, intensity) pairs. Next "
            "step: pass for example [(17200.0, 0.8), (17350.0, 0.2)]."
        ) from None
    if arr.size == 0:
        raise ValueError(
            "the stick table is empty, so the curve would be identically "
            "zero. Next step: check the transition filter that produced the "
            "table (an empty table usually means the filter, not the "
            "spectrum)."
        )
    if arr.ndim != 2 or arr.shape[1] != 2:
        raise ValueError(
            f"sticks must be a sequence of (energy, intensity) pairs, got "
            f"shape {arr.shape}. Next step: pass for example "
            f"[(17200.0, 0.8), (17350.0, 0.2)]."
        )
    if not np.all(np.isfinite(arr)):
        raise ValueError(
            "the stick table contains a non-finite energy or intensity. Next "
            "step: drop the non-finite transitions before broadening."
        )
    for energy, intensity in arr:
        rows.append((float(energy), float(intensity)))
    return rows


def _trapezoid(x: np.ndarray, y: np.ndarray) -> float:
    """The trapezoidal integral of ``y`` over ``x`` (no numpy API-version risk)."""
    if x.size < 2:
        return 0.0
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))
