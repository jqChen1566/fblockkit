"""Checks of the shared spectra rendering module (analysis/spectra.py; W0 of the
spectroscopy extension).

Every numeric assertion has an analytic reference or an independent
recomputation; none is a snapshot of the implementation's own output:

- the Gaussian and the Lorentzian are compared with their closed forms, their
  FWHM is measured from the curve by half-height interpolation, and their
  areas with the exact finite-range integrals (the Lorentzian's over
  ``[-L, L]`` is ``(2/pi) arctan(2L/fwhm)``);
- the Voigt profile is checked against its line-centre closed form
  ``V(0) = exp(y^2) erfc(y) / (sqrt(2 pi) sigma)`` (``math.erfc``), against an
  independent quadrature of the defining convolution integral, against its two
  analytic limits, and against the Lorentzian far-tail asymptote;
- the Boltzmann weights are recomputed from the closed form with the CODATA
  2018 value ``k_B / (h c) = 0.6950348 cm^-1 / K``;
- the CSV writer is checked by reading the file back, by comparing with the
  exact text it must produce (which pins the absence of a timestamp), and by
  byte-comparing two runs.

The one tolerance not set by an analytic bound is the Voigt evaluation's own
accuracy; the module documents 3e-8 of the peak height (two coefficient-free
standard expansions of the Faddeeva function) and the assertions there use
1e-6.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from pytest import approx

from fblockkit.analysis import plot_csv, spectra

# CODATA 2018: k_B / (h c) in cm^-1 / K (the module combines the exact SI
# constants; this is the literature value used to recompute the weights).
KB_CM1_PER_K = 0.6950348


def _measured_fwhm(x, y):
    """The full width at half maximum read off a sampled symmetric curve."""
    half = y.max() / 2.0
    above = np.where(y >= half)[0]
    left, right = int(above[0]), int(above[-1])
    x_left = x[left - 1] + (half - y[left - 1]) / (y[left] - y[left - 1]) * (
        x[left] - x[left - 1]
    )
    x_right = x[right] + (half - y[right]) / (y[right + 1] - y[right]) * (
        x[right + 1] - x[right]
    )
    return x_right - x_left


def _voigt_quadrature(x, center, fwhm_g, fwhm_l, points=200001):
    """The definition V = G (*) L, evaluated by an independent quadrature.

    With ``s = tan(theta)`` the convolution integral becomes
    ``V(x) = (1/pi) integral exp(-(xt + yt tan(theta))^2) dtheta`` over
    ``(-pi/2, pi/2)`` in the Gaussian's reduced units ``xt = (x - center) /
    (sqrt(2) sigma)``, ``yt = gamma / (sqrt(2) sigma)``; the midpoint rule on
    that smooth integrand is a route to the profile that shares no code with
    the module.
    """
    sigma = fwhm_g / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    gamma = fwhm_l / 2.0
    xt = (np.asarray(x, dtype=float) - center) / (math.sqrt(2.0) * sigma)
    yt = gamma / (math.sqrt(2.0) * sigma)
    theta = (np.arange(points) + 0.5) * (math.pi / points) - math.pi / 2.0
    values = np.empty_like(xt)
    for index, reduced in enumerate(xt):
        values[index] = np.mean(np.exp(-((reduced + yt * np.tan(theta)) ** 2)))
    return values / (math.sqrt(2.0 * math.pi) * sigma)


# --- the line shapes ---------------------------------------------------------------


def test_gaussian_matches_its_closed_form():
    x = np.linspace(-500.0, 500.0, 4001)
    fwhm = 120.0
    sigma = fwhm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    closed_form = np.exp(-(x - 40.0) ** 2 / (2.0 * sigma**2)) / (
        sigma * math.sqrt(2.0 * math.pi)
    )
    assert spectra.gaussian(x, 40.0, fwhm) == approx(closed_form, rel=1e-13)


def test_gaussian_half_height_width_is_the_fwhm():
    fwhm = 120.0
    x = np.linspace(-400.0, 400.0, 16001)
    y = spectra.gaussian(x, 0.0, fwhm)
    assert _measured_fwhm(x, y) == approx(fwhm, rel=1e-4)
    # the peak height of a unit-area Gaussian is 2 sqrt(ln 2 / pi) / fwhm
    assert y.max() == approx(2.0 * math.sqrt(math.log(2.0) / math.pi) / fwhm, rel=1e-12)


def test_gaussian_area_matches_the_error_function():
    fwhm = 120.0
    half_span = 600.0
    x = np.linspace(-half_span, half_span, 200001)
    sigma = fwhm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    area = float(np.trapezoid(spectra.gaussian(x, 0.0, fwhm), x))
    assert area == approx(math.erf(half_span / (math.sqrt(2.0) * sigma)), rel=1e-9)


def test_lorentzian_matches_its_closed_form():
    x = np.linspace(-1000.0, 1000.0, 4001)
    fwhm = 80.0
    half = fwhm / 2.0
    closed_form = (half / math.pi) / ((x + 25.0) ** 2 + half**2)
    assert spectra.lorentzian(x, -25.0, fwhm) == approx(closed_form, rel=1e-13)


def test_lorentzian_half_height_width_is_the_fwhm():
    fwhm = 80.0
    x = np.linspace(-4000.0, 4000.0, 400001)
    y = spectra.lorentzian(x, 0.0, fwhm)
    assert _measured_fwhm(x, y) == approx(fwhm, rel=1e-4)
    # the peak height of a unit-area Lorentzian is 2 / (pi fwhm)
    assert y.max() == approx(2.0 / (math.pi * fwhm), rel=1e-12)


def test_lorentzian_area_matches_the_arctangent_integral():
    fwhm = 80.0
    half_span = 5000.0
    x = np.linspace(-half_span, half_span, 400001)
    area = float(np.trapezoid(spectra.lorentzian(x, 0.0, fwhm), x))
    expected = (2.0 / math.pi) * math.atan(2.0 * half_span / fwhm)
    # the trapezoid rule truncates the last cells of a 1/x^2 tail
    assert area == approx(expected, rel=1e-7)


def test_voigt_degenerates_exactly_to_its_components():
    x = np.linspace(-300.0, 300.0, 1201)
    assert np.array_equal(spectra.voigt(x, 0.0, 120.0, 0.0), spectra.gaussian(x, 0.0, 120.0))
    assert np.array_equal(spectra.voigt(x, 0.0, 0.0, 120.0), spectra.lorentzian(x, 0.0, 120.0))


def test_voigt_approaches_the_pure_shapes_for_vanishing_widths():
    x = np.linspace(-400.0, 400.0, 4001)
    gaussian = spectra.gaussian(x, 0.0, 120.0)
    lorentzian = spectra.lorentzian(x, 0.0, 120.0)
    nearly_gaussian = spectra.voigt(x, 0.0, 120.0, 1e-6)
    nearly_lorentzian = spectra.voigt(x, 0.0, 1e-6, 120.0)
    assert np.max(np.abs(nearly_gaussian - gaussian)) / gaussian.max() < 1e-5
    assert np.max(np.abs(nearly_lorentzian - lorentzian)) / lorentzian.max() < 1e-3


def test_voigt_centre_matches_the_closed_form():
    # V(0) = exp(y^2) erfc(y) / (sqrt(2 pi) sigma) with y = gamma / (sqrt(2) sigma)
    for fwhm_g, fwhm_l in [(300.0, 50.0), (250.0, 250.0), (100.0, 400.0), (5.0, 0.05)]:
        sigma = fwhm_g / (2.0 * math.sqrt(2.0 * math.log(2.0)))
        gamma = fwhm_l / 2.0
        y = gamma / (math.sqrt(2.0) * sigma)
        expected = math.exp(y * y) * math.erfc(y) / (math.sqrt(2.0 * math.pi) * sigma)
        got = float(spectra.voigt(np.array([0.0]), 0.0, fwhm_g, fwhm_l)[0])
        assert got == approx(expected, rel=1e-6)


def test_voigt_matches_an_independent_quadrature():
    x = np.array([-500.0, -120.0, 0.0, 130.0, 800.0])
    for fwhm_g, fwhm_l in [(300.0, 300.0), (300.0, 50.0)]:
        assert spectra.voigt(x, 0.0, fwhm_g, fwhm_l) == approx(
            _voigt_quadrature(x, 0.0, fwhm_g, fwhm_l), rel=1e-4
        )


def test_voigt_is_symmetric_with_the_peak_at_the_centre():
    x = np.linspace(-60.0, 60.0, 12001)
    y = spectra.voigt(x, 30.0, 10.0, 3.0)
    mirrored = spectra.voigt(-x, -30.0, 10.0, 3.0)
    assert y == approx(mirrored, rel=1e-13)
    assert x[int(np.argmax(y))] == approx(30.0, abs=1e-9)


def test_voigt_fwhm_brackets_the_component_widths():
    x = np.linspace(-2000.0, 2000.0, 400001)
    for fwhm_g, fwhm_l in [(300.0, 50.0), (300.0, 300.0), (10.0, 500.0)]:
        measured = _measured_fwhm(x, spectra.voigt(x, 0.0, fwhm_g, fwhm_l))
        assert max(fwhm_g, fwhm_l) <= measured <= fwhm_g + fwhm_l
    # a pure component recovers its own width through the same machinery
    assert _measured_fwhm(x, spectra.voigt(x, 0.0, 300.0, 0.0)) == approx(300.0, rel=1e-4)


def test_voigt_area_is_one():
    fwhm_g, fwhm_l = 300.0, 300.0
    half_span = 40000.0
    x = np.linspace(-half_span, half_span, 800001)
    area = float(np.trapezoid(spectra.voigt(x, 0.0, fwhm_g, fwhm_l), x))
    tail = fwhm_l / (2.0 * math.pi) * 2.0 / half_span  # the 1/x^2 wing beyond the grid
    assert area + tail == approx(1.0, rel=1e-4)


def test_voigt_far_wings_follow_the_lorentzian_asymptote():
    x = np.array([3.0e4, 1.0e5])
    gamma = 25.0
    values = spectra.voigt(x, 0.0, 300.0, 2.0 * gamma)
    assert values == approx(gamma / (math.pi * x**2), rel=1e-2)


# --- broadening the stick table ------------------------------------------------------


def test_broaden_reproduces_a_single_stick_and_superposes_two():
    grid = np.linspace(16800.0, 17800.0, 20001)
    one = spectra.broaden([(17200.0, 1.0)], grid, shape="gaussian", fwhm=120.0)
    assert one == approx(spectra.gaussian(grid, 17200.0, 120.0), rel=1e-13)
    two = spectra.broaden(
        [(17200.0, 1.0), (17400.0, 0.5)], grid, shape="gaussian", fwhm=120.0
    )
    assert two == approx(one + spectra.gaussian(grid, 17400.0, 120.0) * 0.5, rel=1e-12)


def test_broaden_is_linear_in_the_intensity_scale():
    grid = np.linspace(0.0, 1000.0, 4001)
    sticks = [(200.0, 0.4), (500.0, 1.3)]
    doubled = spectra.broaden(
        [(200.0, 0.8), (500.0, 2.6)], grid, shape="lorentzian", fwhm=40.0
    )
    single = spectra.broaden(sticks, grid, shape="lorentzian", fwhm=40.0)
    assert doubled == approx(2.0 * single, rel=1e-13)


def test_broaden_preserves_the_total_intensity_as_area():
    # a grid wide enough that the truncated Lorentzian wings (which fall off
    # only as gamma / (pi x^2)) sit below the tolerance: the leading missing
    # fraction is (fwhm_l / 2 pi) (1/L_low + 1/L_high) ~ 3e-4 here.
    grid = np.linspace(-30000.0, 50000.0, 400001)
    sticks = [(17200.0, 0.8), (17350.0, 0.2)]
    gaussian = spectra.broaden(sticks, grid, shape="gaussian", fwhm=100.0)
    assert float(np.trapezoid(gaussian, grid)) == approx(1.0, rel=1e-9)
    for shape, fwhm in [("lorentzian", 100.0), ("voigt", (100.0, 60.0))]:
        curve = spectra.broaden(sticks, grid, shape=shape, fwhm=fwhm)
        assert float(np.trapezoid(curve, grid)) == approx(1.0, rel=1e-3)


def test_broaden_keeps_the_peaks_at_the_stick_energies():
    grid = np.linspace(0.0, 100.0, 1001)
    curve = spectra.broaden([(30.0, 1.0), (70.0, 1.0)], grid, shape="gaussian", fwhm=2.0)
    spacing = grid[1] - grid[0]
    left = np.where(grid < 50.0, curve, -np.inf)
    right = np.where(grid > 50.0, curve, -np.inf)
    assert grid[np.argmax(left)] == approx(30.0, abs=spacing)
    assert grid[np.argmax(right)] == approx(70.0, abs=spacing)


def test_broaden_voigt_with_a_zero_component_matches_the_pure_shape():
    grid = np.linspace(-500.0, 500.0, 10001)
    assert spectra.broaden([(0.0, 1.0)], grid, shape="voigt", fwhm=(120.0, 0.0)) == approx(
        spectra.gaussian(grid, 0.0, 120.0), rel=1e-13
    )
    assert spectra.broaden([(0.0, 1.0)], grid, shape="voigt", fwhm=(0.0, 120.0)) == approx(
        spectra.lorentzian(grid, 0.0, 120.0), rel=1e-13
    )


BAD_BROADEN_CALLS = [
    # (sticks, grid, kwargs, fragment)
    ([(1.0, 1.0)], [0.0, 1.0, 0.5], dict(shape="gaussian", fwhm=1.0), "strictly increasing"),
    ([(1.0, 1.0)], [0.0, 0.0, 1.0], dict(shape="gaussian", fwhm=1.0), "strictly increasing"),
    ([(1.0, 1.0)], [0.0, 1.0], dict(shape="gaussian", fwhm=-1.0), "must be positive"),
    ([(1.0, 1.0)], [0.0, 1.0], dict(shape="gaussian", fwhm=0.0), "must be positive"),
    ([(1.0, 1.0)], [0.0, 1.0], dict(shape="gaussian", fwhm=float("nan")), "must be finite"),
    ([(1.0, 1.0)], [0.0, 1.0], dict(shape="lorentzian", fwhm="wide"), "must be a number"),
    ([(1.0, 1.0)], [0.0, 1.0], dict(shape="triangle", fwhm=1.0), "unknown line shape"),
    ([], [0.0, 1.0], dict(shape="gaussian", fwhm=1.0), "empty"),
    ([(1.0, 1.0, 1.0)], [0.0, 1.0], dict(shape="gaussian", fwhm=1.0), "pairs"),
    ([(1.0, float("inf"))], [0.0, 1.0], dict(shape="gaussian", fwhm=1.0), "non-finite"),
    ([(1.0, 1.0)], [0.0, 1.0], dict(shape="voigt", fwhm=120.0), "two component widths"),
    ([(1.0, 1.0)], [0.0, 1.0], dict(shape="voigt", fwhm=(120.0,)), "exactly two"),
    ([(1.0, 1.0)], [0.0, 1.0], dict(shape="voigt", fwhm=(0.0, 0.0)), "at least one non-zero"),
    ([(1.0, 1.0)], [0.0, 1.0], dict(shape="voigt", fwhm=(120.0, -1.0)), "non-negative"),
    ([(1.0, 1.0)], [float("nan"), 1.0], dict(shape="gaussian", fwhm=1.0), "non-finite"),
]


@pytest.mark.parametrize("sticks, grid, kwargs, fragment", BAD_BROADEN_CALLS)
def test_broaden_rejects_bad_inputs_with_a_next_step(sticks, grid, kwargs, fragment):
    with pytest.raises(ValueError) as excinfo:
        spectra.broaden(sticks, grid, **kwargs)
    message = str(excinfo.value)
    assert fragment in message
    assert "Next step:" in message


# --- Boltzmann weights ---------------------------------------------------------------


def test_boltzmann_weights_follow_the_closed_form():
    energies = np.array([0.0, 1000.0, 5000.0])
    degeneracies = np.array([2.0, 4.0, 6.0])
    temperature = 350.0
    boltzmann = degeneracies * np.exp(-(energies - energies.min()) * KB_CM1_PER_K / temperature)
    expected = boltzmann / boltzmann.sum()
    got = spectra.boltzmann_weights(energies, temperature, degeneracies)
    assert got == approx(expected, rel=1e-7)
    assert float(got.sum()) == approx(1.0, rel=1e-14)
    assert np.all(np.diff(got) < 0.0)


def test_boltzmann_high_temperature_limit_is_degeneracy_weighted():
    energies = np.array([0.0, 1000.0, 5000.0])
    degeneracies = np.array([2.0, 4.0, 6.0])
    # the approach to the limit is O(max(E) / kT) = 1.4e-7 at 1e10 K
    weights = spectra.boltzmann_weights(energies, 1.0e10, degeneracies)
    assert weights == approx(degeneracies / degeneracies.sum(), rel=1e-6)


def test_boltzmann_low_temperature_leaves_the_ground_level_alone():
    weights = spectra.boltzmann_weights([0.0, 1000.0, 2000.0], 1.0)
    assert weights[0] == 1.0
    # 1000 cm^-1 at 1 K is exp(-695) of the ground population
    assert weights[1] < 1e-100 and weights[2] < 1e-100


def test_boltzmann_weights_ignore_a_common_energy_offset():
    energies = np.array([0.0, 1000.0, 5000.0])
    plain = spectra.boltzmann_weights(energies, 300.0)
    shifted = spectra.boltzmann_weights(energies + 1.0e6, 300.0)
    assert shifted == approx(plain, rel=1e-13)
    assert float(shifted.sum()) == approx(1.0, rel=1e-14)


def test_boltzmann_degeneracy_is_a_multiplicity_factor():
    energies = np.array([0.0, 1000.0])
    weights = spectra.boltzmann_weights(energies, 300.0, [3.0, 1.0])
    ratio = weights[0] / weights[1]
    expected = 3.0 * math.exp(1000.0 * KB_CM1_PER_K / 300.0)
    assert float(ratio) == approx(expected, rel=1e-7)


@pytest.mark.parametrize("temperature", [0.0, -300.0, float("nan"), float("inf")])
def test_boltzmann_rejects_a_non_positive_temperature(temperature):
    with pytest.raises(ValueError) as excinfo:
        spectra.boltzmann_weights([0.0, 100.0], temperature)
    assert "temperature_k must be positive" in str(excinfo.value)
    assert "Next step:" in str(excinfo.value)


def test_boltzmann_rejects_bad_energy_and_degeneracy_tables():
    with pytest.raises(ValueError) as excinfo:
        spectra.boltzmann_weights([], 300.0)
    assert "non-empty" in str(excinfo.value) and "Next step:" in str(excinfo.value)
    with pytest.raises(ValueError) as excinfo:
        spectra.boltzmann_weights([0.0, float("nan")], 300.0)
    assert "non-finite" in str(excinfo.value) and "Next step:" in str(excinfo.value)
    with pytest.raises(ValueError) as excinfo:
        spectra.boltzmann_weights(["ground", "excited"], 300.0)
    assert "energies_cm1 must be a sequence of numbers" in str(excinfo.value)
    with pytest.raises(ValueError) as excinfo:
        spectra.boltzmann_weights([0.0, 100.0], 300.0, ["two", "four"])
    assert "degeneracies must be a sequence of numbers" in str(excinfo.value)
    with pytest.raises(ValueError) as excinfo:
        spectra.boltzmann_weights([0.0, 100.0], 300.0, [1.0])
    assert "one degeneracy per level" in str(excinfo.value)
    for bad in ([1.0, 0.0], [1.0, -2.0], [1.0, float("inf")]):
        with pytest.raises(ValueError) as excinfo:
            spectra.boltzmann_weights([0.0, 100.0], 300.0, bad)
        assert "positive and finite" in str(excinfo.value)


# --- the one-call combiner ------------------------------------------------------------


def test_weighted_spectrum_selects_sticks_through_explicit_weights():
    grid = np.linspace(0.0, 1000.0, 5001)
    both = spectra.weighted_spectrum(
        [(200.0, 1.0), (600.0, 1.0)], grid, fwhm=50.0, weights=[1.0, 0.0]
    )
    lone = spectra.broaden([(200.0, 1.0)], grid, fwhm=50.0)
    assert both == approx(lone, rel=1e-13)


def test_weighted_spectrum_temperature_route_recomputes_boltzmann_and_broadening():
    grid = np.linspace(0.0, 1000.0, 5001)
    energies = [(200.0, 1.0), (600.0, 1.0)]
    curve = spectra.weighted_spectrum(energies, grid, shape="lorentzian", fwhm=50.0, temperature_k=300.0)
    weights = spectra.boltzmann_weights([200.0, 600.0], 300.0)
    expected = weights[0] * spectra.lorentzian(grid, 200.0, 50.0) + weights[1] * spectra.lorentzian(
        grid, 600.0, 50.0
    )
    assert curve == approx(expected, rel=1e-13)


def test_weighted_spectrum_energy_exponent_scales_by_the_energy_power():
    grid = np.linspace(0.0, 1000.0, 5001)
    sticks = [(200.0, 1.0), (600.0, 1.0)]
    curve = spectra.weighted_spectrum(sticks, grid, fwhm=50.0, energy_exponent=3.0)
    expected = spectra.broaden(
        [(200.0, 200.0**3), (600.0, 600.0**3)], grid, fwhm=50.0
    )
    assert curve == approx(expected, rel=1e-12)


def test_weighted_spectrum_normalize_rescales_to_unit_area():
    grid = np.linspace(0.0, 1000.0, 200001)
    curve = spectra.weighted_spectrum(
        [(200.0, 1.0), (600.0, 3.0)], grid, fwhm=50.0, normalize=True
    )
    assert float(np.trapezoid(curve, grid)) == approx(1.0, rel=1e-6)


@pytest.mark.parametrize(
    "kwargs, fragment",
    [
        (dict(weights=[1.0, 1.0], temperature_k=300.0), "either explicit weights or a temperature"),
        (dict(degeneracies=[1.0, 1.0]), "degeneracies only enter through the Boltzmann populations"),
        (dict(weights=[1.0]), "one weight per stick"),
        (dict(weights=[1.0, float("nan")]), "non-finite entry"),
        (dict(weights=["strong", "weak"]), "weights must be a sequence of numbers"),
        (dict(energy_exponent=3.0), "positive stick energies"),
        (dict(energy_exponent=float("inf")), "finite number"),
    ],
)
def test_weighted_spectrum_rejects_conflicting_or_impossible_requests(kwargs, fragment):
    grid = np.linspace(0.0, 1000.0, 101)
    sticks = [(200.0, 1.0), (600.0, 1.0)] if kwargs.get("energy_exponent") is None else [(0.0, 1.0), (600.0, 1.0)]
    with pytest.raises(ValueError) as excinfo:
        spectra.weighted_spectrum(sticks, grid, fwhm=50.0, **kwargs)
    assert fragment in str(excinfo.value)
    assert "Next step:" in str(excinfo.value)


# --- the tidy CSV ----------------------------------------------------------------------


def test_csv_text_is_exact_and_carries_no_timestamp(tmp_path):
    grid = np.array([0.0, 1.0, 2.5])
    intensity = np.array([0.0, 1.0, 0.5])
    path = tmp_path / "curve.csv"
    returned = spectra.write_spectrum_csv(
        path, grid, intensity, {"shape": "gaussian", "fwhm_cm1": 100.0}
    )
    assert returned == path
    assert path.read_text(encoding="utf-8") == (
        "# shape: gaussian\n"
        "# fwhm_cm1: 100\n"
        "energy_cm1,intensity\n"
        "0,0\n"
        "1,1\n"
        "2.5,0.5\n"
    )


def test_csv_is_byte_identical_between_runs(tmp_path):
    grid = np.linspace(17000.0, 17600.0, 601)
    curve = spectra.broaden([(17200.0, 1.0)], grid, shape="voigt", fwhm=(120.0, 40.0))
    meta = {"shape": "voigt", "fwhm_g_cm1": 120.0, "fwhm_l_cm1": 40.0, "source": "menu-31 probe"}
    first = spectra.write_spectrum_csv(tmp_path / "a.csv", grid, curve, meta)
    second = spectra.write_spectrum_csv(tmp_path / "b.csv", grid, curve, meta)
    assert first.read_bytes() == second.read_bytes()


def test_csv_round_trip_reproduces_the_written_numbers(tmp_path):
    grid = np.linspace(0.0, 7.0, 15)
    intensity = np.linspace(-1.0, 1.0, 15) ** 3
    path = spectra.write_spectrum_csv(tmp_path / "round.csv", grid, intensity, None)
    lines = path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "energy_cm1,intensity"
    assert len(lines) == 1 + grid.size
    for line, energy, value in zip(lines[1:], grid, intensity):
        written = line.split(",")
        assert float(written[0]) == float(plot_csv.plot_number(energy))
        assert float(written[1]) == float(plot_csv.plot_number(value))


def test_csv_meta_columns_and_comment_lines(tmp_path):
    path = tmp_path / "meta.csv"
    spectra.write_spectrum_csv(
        path,
        [1.0, 2.0],
        [0.25, 0.75],
        {"x_column": "energy_eV", "y_column": "cross_section_Mb", "shape": "lorentzian", "temperature_k": 300},
    )
    lines = path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "# shape: lorentzian"
    assert lines[1] == "# temperature_k: 300"
    assert lines[2] == "energy_eV,cross_section_Mb"
    assert lines[3] == "1,0.25" and lines[4] == "2,0.75"


def test_csv_writes_non_finite_intensities_as_their_tokens(tmp_path):
    path = tmp_path / "edge.csv"
    spectra.write_spectrum_csv(path, [1.0, 2.0, 3.0], [float("inf"), float("-inf"), float("nan")], None)
    assert path.read_text(encoding="utf-8").splitlines()[1:] == [
        "1,inf",
        "2,-inf",
        "3,nan",
    ]


@pytest.mark.parametrize(
    "grid, intensity, meta, fragment",
    [
        ([0.0, 1.0], [1.0], None, "one intensity per grid point"),
        ([0.0, 1.0], ["a", "b"], None, "intensity must be a sequence of numbers"),
        ([1.0, 0.0], [1.0, 1.0], None, "strictly increasing"),
        ([0.0, 1.0], [1.0, 1.0], {"rows": [1, 2]}, "only numbers, booleans and strings"),
        ([0.0, 1.0], [1.0, 1.0], {"source": None}, "is None"),
        ([0.0, 1.0], [1.0, 1.0], {"x_column": "a,b"}, "without commas"),
        ([0.0, 1.0], [1.0, 1.0], {"note": "two\nlines"}, "line break"),
        ([0.0, 1.0], [1.0, 1.0], 7, "meta must be a mapping"),
    ],
)
def test_csv_rejects_bad_inputs_with_a_next_step(tmp_path, grid, intensity, meta, fragment):
    with pytest.raises(ValueError) as excinfo:
        spectra.write_spectrum_csv(tmp_path / "x.csv", grid, intensity, meta)
    assert fragment in str(excinfo.value)
    assert "Next step:" in str(excinfo.value)


# --- the chain end to end ----------------------------------------------------------------


def test_sticks_to_weighted_curve_to_csv(tmp_path):
    """The W1-shaped chain: a stick table, thermal weighting, broadening, CSV."""
    energies = np.array([17200.0, 17350.0, 17900.0])
    intensities = np.array([0.62, 0.28, 0.10])
    grid = np.linspace(16800.0, 18300.0, 150001)
    temperature = 300.0
    weights = spectra.boltzmann_weights(energies, temperature)
    curve = spectra.weighted_spectrum(
        list(zip(energies, intensities)), grid, shape="voigt", fwhm=(90.0, 40.0),
        weights=weights,
    )
    # the ground stick dominates the band (its Boltzmann factor is e^-0.72
    # against e^-5.04 for the 17900 stick, and its intensity is the largest);
    # the overlapping 17350 line pulls the maximum slightly towards itself
    assert grid[np.argmax(curve)] == approx(energies[0], abs=5.0)
    # the curve is the weighted superposition of the three Voigt profiles
    manual = sum(
        intensity * weight * spectra.voigt(grid, energy, 90.0, 40.0)
        for energy, intensity, weight in zip(energies, intensities, weights)
    )
    assert curve == approx(manual, rel=1e-12)
    # the thermal weights match the closed form with the literature k_B
    closed = np.exp(-(energies - energies[0]) * KB_CM1_PER_K / temperature)
    assert weights == approx(closed / closed.sum(), rel=1e-7)
    path = spectra.write_spectrum_csv(
        tmp_path / "chain.csv", grid, curve, {"shape": "voigt", "temperature_k": temperature}
    )
    lines = path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "# shape: voigt" and lines[1] == "# temperature_k: 300"
    assert lines[2] == "energy_cm1,intensity"
    assert len(lines) == 3 + grid.size
