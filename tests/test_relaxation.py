"""Checks of the magnetic relaxation machinery (menu 36).

Two fixture layers, both real ORCA 6.1.1 output:

- ``fixtures/single_aniso/``: the SINGLE_ANISO-embedded output -- O2 (integer
  pseudospin S=1, four groups, two segments because DoSSC was on), CO+ with
  MLTP 2,2 (two Kramers doublets, UBAR present) and CO+ with the default MLTP
  (no groups, g/D analysis "NOT INCLUDED");
- ``fixtures/magrelax/``: the Orca_Magrelax run of CO+ (a single vibrational
  mode at 2299.9 cm^-1 against a 14680.6 cm^-1 Zeeman gap -- the measured
  all-zero rate table).

Anchors: the segment-1 D of O2 equals the 5.2 QDPT "effective hamiltonian
soc" value 2.175287 and the segment-2 D the SOC+SSC value 3.185964 (the same
numbers the ZFS blocks print); the UBAR averages are the engine's printed
``(|mu_X| + |mu_Y| + |mu_Z|)/3``; the Orbach fit round-trips a synthetic
Arrhenius table; the near-isotropic CO+ g-tensors are rejected by the
axis-definition guard (their Delta g/g ~ 2e-4 is set by the printed digits).
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest
from pytest import approx

from fblockkit.analysis import relaxation
from fblockkit.diagnosis import references_section
from fblockkit.parsers import parse_auto

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
SA = FIXTURES / "single_aniso"
MR = FIXTURES / "magrelax"


def _segment(fixture: str, index: int = 0) -> dict:
    parsed = parse_auto(SA / fixture)
    return parsed.sections["single_aniso"]["segments"][index]


# --- the SINGLE_ANISO parse ----------------------------------------------------


def test_the_orca_doublets_payload_feeds_menu_16():
    """The ORCA SINGLE_ANISO segment converts into the menu-16 table -- the
    engine printout enters the criterion by the same door as the manual
    table (g ascending, the axis of the largest g, the group's lowest
    spin-orbit energy)."""
    from fblockkit.analysis import magnetic_doublets
    from fblockkit.parsers import single_aniso as orca_single_aniso

    parsed = parse_auto(SA / "co_aniso2.out")
    segments = parsed.sections["single_aniso"]["segments"]
    payload = orca_single_aniso.doublets_payload(segments, system="co_aniso2.out")
    assert payload["system"] == "co_aniso2.out" and payload["reference"] == 0
    assert len(payload["doublets"]) == 2
    row = payload["doublets"][0]
    assert row["g"][0] <= row["g"][1] <= row["g"][2]
    assert row["energy"] == 0.0
    assert payload["doublets"][1]["energy"] == approx(29361.197)
    table = magnetic_doublets.parse_doublets(payload)
    report = magnetic_doublets.analyze(table)
    assert len(report.verdicts) == 2


def test_the_orca_doublets_payload_refuses_a_groupless_segment():
    from fblockkit.parsers import single_aniso as orca_single_aniso

    with pytest.raises(ValueError, match="MLTP"):
        orca_single_aniso.doublets_payload([{"groups": []}])
    with pytest.raises(ValueError, match="no SINGLE_ANISO segment"):
        orca_single_aniso.doublets_payload([])


def test_the_o2_output_has_two_segments_with_the_expected_spectra():
    parsed = parse_auto(SA / "o2_aniso.out")
    data = parsed.sections["single_aniso"]
    assert data["present"] and len(data["segments"]) == 2
    first, second = data["segments"]
    assert first["n_soc_states"] == 12 and second["n_soc_states"] == 12
    assert first["soc_spectrum_cm1"][:3] == approx([0.0, 2.1753, 2.1753])
    assert second["soc_spectrum_cm1"][:3] == approx([0.0, 3.186, 3.186])
    assert first["n_spin_free_states"] == 6
    assert first["ubar_included"] is True


def test_the_o2_groups_carry_g_d_and_the_ubar_table():
    segment = _segment("o2_aniso.out")
    assert [group["spin"] for group in segment["groups"]] == [1.0] * 4
    ground = segment["groups"][0]
    assert ground["soc_state_energies_cm1"] == approx([0.0, 2.175, 2.175])
    assert ground["g_values"] == approx((1.99981604, 1.99981604, 2.0))
    assert ground["g_sign_product"] == -1
    assert ground["d_cm1"] == approx(2.175287)
    assert ground["e_cm1"] == approx(0.0)
    assert ground["angular_moments_z"] == approx((0.0, 0.0, -2.0))
    ubar = segment["ubar"]
    assert ubar["present"] and len(ubar["zeeman_states"]) == 4
    assert ubar["zeeman_states"][0]["columns"] == 4  # integer pseudospin keeps m = 0
    assert ubar["zeeman_states"][-1]["energy_cm1"] == approx(49025.6424797)
    assert ubar["even_electron_note"] is True
    intra = [
        element
        for element in ubar["matrix_elements"]
        if element["delta_mult"] == 0 and element["average"] is not None
    ]
    assert len(intra) == 4
    assert intra[3]["average"] == approx(0.667963206844)


def test_the_d_values_match_the_qdpt_zfs_variants():
    parsed = parse_auto(SA / "o2_aniso.out")
    first, second = parsed.sections["single_aniso"]["segments"]
    # the same numbers the 5.2 QDPT ZFS blocks print (fixtures/orca/o2_qdpt2.out):
    # the SOC-only variant and the SOC+SSC effective-Hamiltonian variant
    assert first["groups"][0]["d_cm1"] == approx(2.175287)
    assert second["groups"][0]["d_cm1"] == approx(3.185964)
    zfs = parsed.sections["zfs"]
    variants = {block["variant"]: block["D_cm1"] for block in zfs["blocks"]}
    assert variants["EFFECTIVE HAMILTONIAN SOC CONTRIBUTION"] == approx(2.175287)
    assert variants["EFFECTIVE HAMILTONIAN SOC and SSC CONTRIBUTION"] == approx(3.185964)


def test_the_co_plus_kd_sample_and_its_missing_zfs_block():
    segment = _segment("co_aniso2.out")
    assert [group["spin"] for group in segment["groups"]] == [0.5, 0.5]
    ground = segment["groups"][0]
    assert ground["g_values"] == approx((1.99999879, 1.99999879, 1.9995212))
    assert ground["d_cm1"] is None  # a Kramers doublet prints no ZFS/D block
    assert ground["angular_moments_z"] == approx((0.0, 0.0, -0.999760598))
    ubar = segment["ubar"]
    assert ubar["zeeman_states"][0]["columns"] == 3  # Kramers: m+, m-, E
    assert ubar["zeeman_states"][0]["m_plus"] == approx(-0.9997606)
    intra = [
        element
        for element in ubar["matrix_elements"]
        if element["delta_mult"] == 0 and element["average"] is not None
    ]
    assert intra[0]["average"] == approx(0.666666263248)


def test_the_default_mltp_sample_prints_no_groups():
    segment = _segment("co_aniso.out")
    assert segment["g_tensor_included"] is False
    assert segment["ubar_included"] is False
    assert segment["groups"] == []
    assert segment["soc_spectrum_cm1"][:3] == approx([0.0, 0.0, 29361.197])


# --- the QTM metrics -----------------------------------------------------------


def test_the_nearly_isotropic_co_plus_groups_are_not_flagged():
    rows = relaxation.group_metrics(_segment("co_aniso2.out"))
    assert all(not row["axis_defined"] for row in rows)
    assert all(row["g_t_theta"] is None for row in rows)
    assert relaxation.ueff_estimate(rows) is None


def _group(index: int, spin: float, energies: list[float], g, axes) -> dict:
    return {
        "index": index,
        "spin": spin,
        "soc_state_energies_cm1": energies,
        "g_values": g,
        "g_axes": axes,
        "g_sign_product": -1,
        "d_values": None,
        "d_axes_main": None,
        "d_axes_cartesian": None,
        "d_cm1": None,
        "e_cm1": None,
        "angular_moments_z": None,
    }


def test_the_guide_criteria_flag_the_first_noncollinear_group():
    # a synthetic Dy-like ladder: ground easy-axis (0.2, 0.2, 19) along z;
    # the first excited doublet is tilted by 45 degrees with a large transverse g
    import numpy as np

    tilt = [float(np.sin(np.pi / 4)), 0.0, float(np.cos(np.pi / 4))]
    axes = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], tilt]
    rows = relaxation.group_metrics(
        {
            "groups": [
                _group(1, 1.5, [0.0, 0.0], (0.2, 0.2, 19.0), [[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
                _group(2, 1.5, [500.0, 500.0], (3.0, 3.0, 12.0), axes),
                _group(3, 1.5, [900.0, 900.0], (1.0, 1.0, 15.0), [[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
            ]
        }
    )
    assert rows[1]["axis_defined"] and rows[1]["theta_deg"] == approx(45.0)
    # g_T = (3 + 3 + 12 sin 45)/3; the flag fires through both criteria
    assert rows[1]["g_t"] == approx((3.0 + 3.0 + 12.0 * math.sqrt(0.5)) / 3.0)
    assert rows[1]["g_t_theta"] > relaxation.G_T_THETA_THRESHOLD
    estimate = relaxation.ueff_estimate(rows)
    assert estimate is not None and estimate["index"] == 2
    assert estimate["energy_cm1"] == approx(500.0)


def test_a_collinear_ladder_has_no_flag():
    rows = relaxation.group_metrics(
        {
            "groups": [
                _group(1, 1.5, [0.0, 0.0], (0.2, 0.2, 19.0), [[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
                _group(2, 1.5, [300.0, 300.0], (0.1, 0.1, 17.0), [[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
            ]
        }
    )
    assert rows[1]["theta_deg"] == approx(0.0)
    assert relaxation.ueff_estimate(rows) is None


# --- the Orbach fit ------------------------------------------------------------


def test_the_orbach_fit_round_trips_the_arrhenius_table():
    temperatures = [3.0, 4.0, 5.0, 6.0, 7.0]
    tau0, ueff_k = 1e-7, 200.0
    taus = [tau0 * math.exp(ueff_k / temperature) for temperature in temperatures]
    fit = relaxation.orbach_fit(temperatures, taus)
    assert fit["n_points"] == 5
    assert fit["ueff_k"] == approx(ueff_k, rel=1e-9)
    assert fit["ueff_cm1"] == approx(ueff_k / 1.438776877, rel=1e-9)
    assert fit["tau0_s"] == approx(tau0, rel=1e-6)
    assert fit["r_squared"] == approx(1.0, abs=1e-10)


def test_the_orbach_fit_window_and_refusals():
    temperatures = [2.0, 3.0, 4.0, 5.0, 6.0]
    taus = [1e-7 * math.exp(150.0 / temperature) for temperature in temperatures]
    narrowed = relaxation.orbach_fit(temperatures, taus, window_k=(4.0, 6.0))
    assert narrowed["n_points"] == 3
    assert narrowed["ueff_k"] == approx(150.0, rel=1e-6)
    with pytest.raises(relaxation.RelaxationError, match="at least three"):
        relaxation.orbach_fit([2.0, 3.0], taus[:2])


# --- the Magrelax reader -------------------------------------------------------


def test_the_magrelax_output_reads_with_the_all_zero_rate_table():
    data = relaxation.read_magrelax(MR / "co_magrelax.out")
    assert data.magrelax_input == "co_magrelax.magrelaxinp"
    assert data.derivatives_file == "co_magrelax.socder.magrelax"
    assert data.n_derivatives == 1 and data.derivative_size == 4
    assert data.hessian_file == "co_hess.hess"
    assert data.zeeman_eigenvalues_cm1 == approx([0.0, 14680.6, 15089.1, 0.0])
    assert data.frequencies_cm1 == []
    assert len(data.temperatures_k) == 28
    assert data.temperatures_k[0] == approx(2.0) and data.temperatures_k[-1] == approx(29.0)
    assert all(math.isinf(tau) and tau < 0 for tau in data.tau_s)


def test_the_magrelax_render_reports_the_missing_channel():
    data = relaxation.read_magrelax(MR / "co_magrelax.out")
    body = relaxation.render_magrelax(data, "co_magrelax.out")
    assert "Orbach fit: not applicable" in body
    assert "no one-phonon channel can be resonant" in body


def test_the_magrelax_reader_refuses_a_plain_casscf_output():
    with pytest.raises(relaxation.RelaxationError, match="MAGRELAX"):
        relaxation.read_magrelax(SA / "co_aniso2.out")


# --- the report and evidence ---------------------------------------------------


def test_the_single_aniso_render_carries_the_guard_note_and_the_flags():
    segment = _segment("co_aniso2.out")
    body = relaxation.render_single_aniso(
        segment, segment_index=0, n_segments=1, source="co_aniso2.out"
    )
    assert "do not apply" in body  # the axis-definition guard note
    assert "U_eff estimate" in body and "no computed group is flagged" in body
    assert "fixed template sentence" in body


def test_the_evidence_is_citable():
    evidence = relaxation.evidence()
    assert references_section(evidence) is not None
    literature = [item for item in evidence if item.kind == "literature"]
    assert literature and "10.1039/d5cs00493d" in literature[0].ref


def test_the_magrelax_plot_csv_carries_the_tau_table():
    """The plot-ready companion: plain rows mirroring the engine's table,
    including the probe's uniformly-zero-rate boundary (tau = -inf, written
    as the token the number actually is)."""
    data = relaxation.read_magrelax(MR / "co_magrelax.out")
    csv_text = relaxation.magrelax_plot_csv(data)
    lines = csv_text.splitlines()
    assert lines[0] == "temperature_K,rate_per_s,tau_s"
    assert len(lines) == 1 + len(data.tau_s) == 29
    first = lines[1].split(",")
    assert float(first[0]) == approx(2.0)
    assert float(first[1]) == approx(0.0)
    assert float(first[2]) == float("-inf")
    assert lines[-1].split(",")[0] == "29"
    assert csv_text.endswith("\n") and "\r" not in csv_text
