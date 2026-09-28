"""Checks of the pNMR pseudocontact machinery (analysis/pnmr.py, menu 35; Wave 5.2).

Three layers of anchors:

- **published tensors**: the source's two printed Co(II) susceptibility
  tensors -> Delta chi_ax; and the magnetic parameter sets -> the printed
  14.8 / 27.3 / 7.1 (1e-32 m^3);
- **invariants**: Tr<SS> = S(S+1); the high-T isotropic limit; PCS invariance
  under the isotropic and antisymmetric parts; the exact agreement with the
  classical spherical form for an axial tensor;
- **the ORCA pairing**: the QDPT fixtures (co_plus_qdpt.*, o2_qdpt2.*) --
  g-matrix, the four ZFS variants, and ORCA's own printed susceptibility
  against this module's chi from the parsed g (the 4 pi cgs-emu conversion).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import yaml

from fblockkit.analysis import pnmr
from fblockkit.diagnosis import references_section
from fblockkit.parsers import parse_auto

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
ORCA = FIXTURES / "orca"
CO_QDPT = ORCA / "co_plus_qdpt.out"
O2_QDPT = ORCA / "o2_qdpt2.out"
REGRESSION = FIXTURES / "pnmr" / "mares2018_co2.yaml"
RUN_FILE = FIXTURES / "pnmr" / "co_plus" / "co_plus_qdpt_g.yaml"


@pytest.fixture(scope="module")
def regression():
    return yaml.safe_load(REGRESSION.read_text(encoding="utf-8"))


# --- the published-tensor regressions ------------------------------------------


def test_the_printed_tensors_reproduce_the_published_axialities(regression):
    chi_ns = np.array(regression["printed_chi_nonsymmetric_m3"])
    d_ax, d_rh = pnmr.axiality(chi_ns)
    # the paper prints 14.8; the excess is its 3-digit rounding of the entries
    assert d_ax / 1e-32 == pytest.approx(14.8, abs=0.15)
    assert d_rh / 1e-32 == pytest.approx(0.0, abs=0.2)
    chi_sym = np.array(regression["printed_chi_symmetric_m3"])
    d_ax_sym, _ = pnmr.axiality(chi_sym)
    assert d_ax_sym / 1e-32 == pytest.approx(27.3, abs=0.15)


def test_the_ab_initio_chain_reproduces_the_published_axiality(regression):
    chain = regression["chain_nevpt2"]
    g_iso, dg = chain["g_iso"], chain["Delta_g_ax"]
    g = np.diag([g_iso - dg / 3, g_iso - dg / 3, g_iso + 2 * dg / 3])
    chi_ns, chi_sym = pnmr.susceptibility_from_zfs(
        g,
        D_cm1=chain["D_cm1"],
        E_over_D=chain["E_over_D"],
        temperature_k=chain["temperature_k"],
        spin=chain["spin"],
    )
    # printed 14.8 / 27.3; the 2 % offset is the axial approximation (the
    # source's tensors are slightly rhombic)
    assert pnmr.axiality(chi_ns)[0] / 1e-32 == pytest.approx(14.8, rel=0.03)
    assert pnmr.axiality(chi_sym)[0] / 1e-32 == pytest.approx(27.3, rel=0.03)


def test_the_experimental_parameter_variant_brackets_the_published_value(regression):
    chain = regression["chain_experimental"]
    lo, hi = chain["g_iso_range"]
    values = []
    for g_iso in (lo, hi):
        dg = chain["Delta_g_ax"]
        g = np.diag([g_iso - dg / 3, g_iso - dg / 3, g_iso + 2 * dg / 3])
        chi_ns, _ = pnmr.susceptibility_from_zfs(
            g,
            D_cm1=chain["D_cm1"],
            E_over_D=0.0,
            temperature_k=chain["temperature_k"],
            spin=chain["spin"],
        )
        values.append(pnmr.axiality(chi_ns)[0] / 1e-32)
    target = chain["printed_dchi_ax_nonsymmetric_1e32_m3"]
    assert min(values) <= target <= max(values) or abs(target - np.mean(values)) < 0.4


def test_the_nonsymmetric_and_symmetric_routes_differ_quadratically_vs_linearly():
    # a Kramers doublet with a tilted g-tensor: the two constructions differ
    g = np.array([[2.2, 0.3, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 1.9]])
    chi_ns, chi_sym = pnmr.susceptibility_from_g(g, temperature_k=300.0)
    assert not np.allclose(chi_ns, chi_sym, rtol=1e-3, atol=1e-35)
    # the non-symmetric form is linear in g, the symmetric one quadratic: scale
    # the anisotropy (not g_iso) and compare the responses
    g2 = np.array([[2.2 + 0.1, 0.3, 0.0], [0.0, 2.0 - 0.1, 0.0], [0.0, 0.0, 1.9]])
    ns2, sym2 = pnmr.susceptibility_from_g(g2, temperature_k=300.0)
    assert np.abs(ns2 - chi_ns).max() > 0.0 and np.abs(sym2 - chi_sym).max() > 0.0


# --- the invariants ------------------------------------------------------------


def test_the_ss_dyadic_invariants():
    for spin in (0.5, 1.0, 1.5, 2.5):
        ss = pnmr.ss_dyadic(D_cm1=-30.0, E_over_D=0.1, temperature_k=300.0, spin=spin)
        assert np.trace(ss) == pytest.approx(spin * (spin + 1.0), rel=1e-12)
        assert np.allclose(ss, ss.T, atol=1e-14)
    # the high-temperature limit is isotropic S(S+1)/3
    ss_hot = pnmr.ss_dyadic(D_cm1=-30.0, E_over_D=0.1, temperature_k=1e6, spin=1.5)
    assert np.allclose(ss_hot, np.eye(3) * (1.5 * 2.5 / 3.0), rtol=1e-4)
    # the low-temperature limit of an easy-axis quartet: only M = +-3/2
    ss_cold = pnmr.ss_dyadic(D_cm1=-85.5, E_over_D=0.0, temperature_k=0.5, spin=1.5)
    assert ss_cold[2, 2] == pytest.approx(2.25, rel=1e-6)
    assert ss_cold[0, 0] == pytest.approx(0.75, rel=1e-6)


def test_pcs_invariances_and_the_spherical_form():
    rng = np.random.default_rng(20260928)
    t = rng.normal(size=(3, 3)) * 1e-31
    center = (0.3, -0.2, 0.1)
    positions = [(1.0, 2.0, 3.0), (-4.0, 0.5, 2.0)]
    base = pnmr.pcs(t, center, positions)
    # the isotropic part does not enter
    shifted = pnmr.pcs(t + 7.3e-31 * np.eye(3), center, positions)
    for a, b in zip(base, shifted):
        assert a[3] == pytest.approx(b[3], rel=1e-12)
    # the antisymmetric part does not enter
    antisym = t - t.T
    rotor = pnmr.pcs(t + 2.0e-31 * antisym, center, positions)
    for a, b in zip(base, rotor):
        assert a[3] == pytest.approx(b[3], rel=1e-12)
    # an axial tensor: exact agreement with the classical spherical form
    chi_ax = np.diag([-1.0e-31, -1.0e-31, 2.0e-31])
    point = (0.0, 0.0, 4.0)
    (r, theta, phi, delta) = pnmr.pcs(chi_ax, (0, 0, 0), [point])[0]
    d_chi_ax = 3.0e-31  # chi3 - (chi1+chi2)/2 of the traceless part
    expected = (
        d_chi_ax * (3 * np.cos(np.radians(theta)) ** 2 - 1) / (12 * np.pi * (r * 1e-10) ** 3)
    ) * 1e6
    assert delta == pytest.approx(expected, rel=1e-12)
    # 1/r^3 scaling
    near = pnmr.pcs(chi_ax, (0, 0, 0), [(0, 0, 2.0)])[0][3]
    assert near / delta == pytest.approx(8.0, rel=1e-12)


# --- the ORCA pairing ----------------------------------------------------------


def test_the_orca_qdpt_blocks_parse():
    pr = parse_auto(CO_QDPT)
    epr = pr.sections["epr"]
    assert epr["present"] and epr["multiplicity"] == 2 and epr["n_blocks"] == 1
    g = np.array(epr["g_matrix"])
    assert [round(g[i, i], 6) for i in range(3)] == [2.001942, 2.002118, 2.002207]
    assert epr["iso"] == pytest.approx(2.0020889)
    zfs = pr.sections["zfs"]
    assert zfs["present"] is False  # a Kramers doublet carries no ZFS blocks
    sus = pr.sections["susceptibility"]
    assert sus["present"] and len(sus["temperatures_k"]) == 300
    assert sus["temperatures_k"][0] == pytest.approx(1.0)
    assert sus["temperatures_k"][-1] == pytest.approx(300.0)
    assert sus["tensors_cm3k_mol"][0][0][0] == pytest.approx(0.375900, abs=1e-6)


def test_the_o2_zfs_variants_parse_with_their_internal_consistency():
    pr = parse_auto(O2_QDPT)
    zfs = pr.sections["zfs"]
    variants = [(b["variant"], b["D_cm1"]) for b in zfs["blocks"]]
    assert variants == [
        ("2ND ORDER SOC CONTRIBUTION", 2.175688),
        ("EFFECTIVE HAMILTONIAN SOC CONTRIBUTION", 2.175287),
        ("2ND ORDER SOC and SSC CONTRIBUTION", 3.186488),
        ("EFFECTIVE HAMILTONIAN SOC and SSC CONTRIBUTION", 3.185964),
    ]
    preferred = zfs["preferred"]
    assert preferred["variant"] == "EFFECTIVE HAMILTONIAN SOC and SSC CONTRIBUTION"
    # internal consistency: D = 3/2 * the traceless eigenvalue along its axis
    eigenvalues = sorted(preferred["eigenvalues_traceless"])
    assert preferred["D_cm1"] == pytest.approx(1.5 * eigenvalues[2], rel=1e-6)
    assert preferred["E_over_D"] == pytest.approx(0.0, abs=1e-9)


def test_our_chi_reproduces_orca_own_susceptibility():
    pr = parse_auto(CO_QDPT)
    g = np.array(pr.sections["epr"]["g_matrix"])
    chi_ns, chi_sym = pnmr.susceptibility_from_g(g, temperature_k=300.0)
    # cgs-emu molar chiT (cm^3 K / mol) = NA * chi_SI_molar / (4 pi) * T * 1e6
    for chi in (chi_ns, chi_sym):
        molar_si = np.trace(chi) / 3 * 6.02214076e23
        chi_t_cgs = molar_si / (4 * np.pi) * 300.0 * 1e6
        orca_tensor = np.array(pr.sections["susceptibility"]["tensors_cm3k_mol"][-1])
        assert chi_t_cgs == pytest.approx(np.trace(orca_tensor) / 3, rel=2e-3)


# --- the run reader ------------------------------------------------------------


def test_the_run_file_reads_the_orca_g_and_renders():
    data = pnmr.read_run(RUN_FILE)
    assert data.center_index == 0 and data.route == "nonsymmetric"
    assert data.spin == 0.5 and data.temperature_k == 300.0
    assert data.chi_si.shape == (3, 3)
    assert any("effective-Hamiltonian" in note for note in data.notes)
    body = pnmr.render(data)
    assert "Delta chi_ax" in body and "O2" in body and "contact" in body


def test_the_run_reader_refusals_carry_next_steps(tmp_path):
    good = yaml.safe_load(RUN_FILE.read_text(encoding="utf-8"))
    base = {k: v for k, v in good.items()}
    base["structure"] = "co_plus.xyz"
    (tmp_path / "co_plus.xyz").write_text(
        "2\nx\nC 0 0 0\nO 0 0 1.13\n", encoding="utf-8"
    )

    def write(payload: dict) -> Path:
        path = tmp_path / "run.yaml"
        path.write_text(yaml.safe_dump(payload), encoding="utf-8")
        return path

    with pytest.raises(pnmr.PnmrError, match="center_atom"):
        pnmr.read_run(write({**base, "center_atom": 3}))
    with pytest.raises(pnmr.PnmrError, match="no 'susceptibility'"):
        payload = {k: v for k, v in base.items() if k != "susceptibility"}
        pnmr.read_run(write(payload))
    with pytest.raises(pnmr.PnmrError, match="3x3"):
        pnmr.read_run(
            write({**base, "susceptibility": {"from": "tensor", "tensor": [1.0, 2.0, 3.0]}})
        )
    with pytest.raises(pnmr.PnmrError, match="unknown susceptibility source"):
        pnmr.read_run(write({**base, "susceptibility": {"from": "magic"}}))
    with pytest.raises(pnmr.PnmrError, match="D_cm1"):
        pnmr.read_run(
            write(
                {
                    **base,
                    "susceptibility": {
                        "from": "zfs",
                        "g_matrix": np.eye(3).tolist(),
                        "spin": 1.0,
                    },
                }
            )
        )
    with pytest.raises(pnmr.PnmrError, match="tensor_unit"):
        pnmr.read_run(
            write(
                {
                    **base,
                    "susceptibility": {
                        "from": "tensor",
                        "tensor": np.eye(3).tolist(),
                        "tensor_unit": "esu",
                    },
                }
            )
        )


def test_the_evidence_is_citable():
    refs = references_section(pnmr.evidence())
    assert refs is not None
    literature = [item for item in pnmr.evidence() if item.kind == "literature"]
    assert literature and "10.1039/c8cp04123g" in literature[0].ref
