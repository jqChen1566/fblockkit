"""pNMR pseudocontact shifts from a magnetic susceptibility tensor (menu 35).

The pseudocontact shift (PCS) of a nucleus in the point-dipole approximation
(PDA) is the dipole field of a point magnetic moment at the paramagnetic
centre.  The modern, quantum-chemistry-consistent form (Mares & Vaara, PCCP
2018, eqs. (1) and (2); the close reading is
``文献细读/细读_pNMR_两篇_20260928.md``) is

    sigma_dip = - chi' . d / (4 pi r_IS^3),   d = 3 n_IS n_IS - 1,
    sigma_PCS = Tr(sigma_dip)/3,   delta_PCS = -sigma_PCS = Tr(chi' . d)/(12 pi r^3),

with chi' the (generally **non-symmetric**) magnetic susceptibility tensor in
SI volume units (m^3) and n_IS the unit vector from the paramagnetic centre to
the nucleus.  Because d is symmetric and traceless, only the traceless
symmetric part of chi' enters delta_PCS -- which is exactly why the
*non-symmetric* construction of chi' matters: its symmetric part differs from
the textbook symmetric formula.

The two constructions (both measured against the source's printed tensors):

- **non-symmetric** (eq. (3) of the source, the PDA-consistent one):

      chi' = (mu_B^2 mu_0 g_e / (k_B T)) g . <SS>

  with <SS> the dyadic of the effective spin built from the thermal average
  over the zero-field-splitting eigenstates (the ZFS Hamiltonian
  H = S.D.S with D = D_zz - (D_xx + D_yy)/2 and E = (D_yy - D_xx)/2) and g
  the g-tensor (full matrix; the free-electron factor g_e is scalar);
- **symmetric** (eq. (5), "the alternative approach"):

      chi = (mu_B^2 mu_0 / (k_B T)) g . <SS> . g^T

  which makes the anisotropy quadratic in the g-anisotropy where the
  non-symmetric form is linear.  For the source's Co(II) benchmark the two
  give 14.8 and 27.3 x 1e-32 m^3 -- this module's regression reproduces both
  from the printed tensors and the printed parameter set (see
  ``fixtures/pnmr/``).

Regression anchors (all digits as printed in the source; see the fixture):

- the two printed chi tensors -> Delta chi_ax = 14.9065e-32 / 27.3120e-32 m^3
  (printed 14.8 / 27.3; the small excess is the paper's 3-digit rounding of
  the matrix entries);
- the ab initio chain {D = -85.5 cm^-1, E/D = 0.00281, g_iso = 2.325,
  Delta g_ax = 0.979} at 300 K -> Delta chi'_ax = 14.47e-32 m^3 (printed
  14.8; the 2 % offset is the axial approximation -- the source's g and D
  tensors are slightly rhombic);
- the experimental-parameter variant {D = -64 cm^-1, Delta g_ax = 0.2,
  g_iso ~ 2.3-2.4} -> 6.8-7.1e-32 m^3 (printed 7.1).

ORCA side (measured on ORCA 6.1.1, fixtures co_plus_qdpt.*, o2_qdpt2.*):
the QDPT driver prints the effective-Hamiltonian g-matrix
("ELECTRONIC G-MATRIX FROM EFFECTIVE HAMILTONIAN"), the zero-field-splitting
blocks (2nd-order and effective-Hamiltonian variants, with and without the
spin-spin contribution) and -- with ``DoSusceptibility`` -- the molar
susceptibility tensor per temperature in cm^3 K / mol (the cgs-emu molar
value: chi_SI_molar = 4 pi chi_cgs).  This module's chi from the parsed g of
the CO+ fixture reproduces ORCA's own printed susceptibility to the printed
precision (cross-check in the tests).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..knowledge.models import EVIDENCE_LITERATURE, EVIDENCE_MEASURED, Evidence

__all__ = [
    "PnmrError",
    "RunData",
    "spin_matrices",
    "ss_dyadic",
    "susceptibility_nonsymmetric",
    "susceptibility_symmetric",
    "susceptibility_from_g",
    "susceptibility_from_zfs",
    "chi_tensors",
    "axiality",
    "pcs",
    "read_run",
    "render",
    "evidence",
]

#: physical constants (CODATA 2018)
_MU_B = 9.2740100783e-24  # J/T
_MU_0 = 1.25663706212e-6  # N/A^2
_K_B = 1.380649e-23  # J/K
_G_E = 2.00231930436256  # free-electron g value
_CM1_PER_K = 1.438776877  # 1 cm^-1 in K (for the ZFS Hamiltonian)
_A0_CM = 0.529177210903e-8


class PnmrError(ValueError):
    """The pNMR input is not usable (with a next step)."""


@dataclass(frozen=True)
class RunData:
    """One pNMR run: structure, centre, susceptibility source and settings."""

    atoms: tuple[tuple[str, float, float, float], ...]
    center_index: int  # 0-based
    temperature_k: float
    route: str  # "nonsymmetric" | "symmetric"
    spin: float
    chi_si: np.ndarray  # 3x3, m^3
    chi_source: str  # a human-readable description of where chi came from
    nuclei: tuple[int, ...] = ()  # 0-based indices to report; () = all others
    structure_name: str = ""
    notes: tuple[str, ...] = field(default=())


def spin_matrices(spin: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The S_x, S_y, S_z matrices for the given spin (ascending m convention)."""
    dim = int(round(2 * spin + 1))
    if dim < 2 or abs(dim - (2 * spin + 1)) > 1e-9:
        raise PnmrError(f"invalid spin {spin!r}.")
    m = np.array([-spin + k for k in range(dim)])
    sp = np.zeros((dim, dim), dtype=complex)
    for k in range(dim - 1):
        sp[k, k + 1] = np.sqrt(spin * (spin + 1) - m[k] * (m[k] + 1))
    sm = sp.conj().T
    return (sp + sm) / 2, (sp - sm) / (2j), np.diag(m).astype(complex)


def ss_dyadic(*, D_cm1: float, E_over_D: float, temperature_k: float, spin: float):
    """The <SS> dyadic: thermal average of S_i S_j over the ZFS eigenstates.

    The ZFS Hamiltonian is S.D.S with the principal values D_zz = 2D/3,
    D_xx = -D/3 + E, D_yy = -D/3 - E (E = E_over_D * D).  Invariant checks
    (tested): Tr<SS> = S(S+1) at every temperature, and the high-temperature
    limit is the isotropic S(S+1)/3.
    """
    if temperature_k <= 0.0:
        raise PnmrError("the temperature must be positive.")
    sx, sy, sz = spin_matrices(spin)
    d_k = D_cm1 * _CM1_PER_K
    e_k = E_over_D * d_k
    h = (-d_k / 3 + e_k) * (sx @ sx) + (-d_k / 3 - e_k) * (sy @ sy) + (2 * d_k / 3) * (sz @ sz)
    w, v = np.linalg.eigh(h)
    w = w - w.min()
    weights = np.exp(-w / temperature_k)
    ss = np.zeros((3, 3), dtype=complex)
    for weight, vec in zip(weights, v.T):
        ket = vec.reshape(-1, 1)
        for i, si in enumerate((sx, sy, sz)):
            for j, sj in enumerate((sx, sy, sz)):
                ss[i, j] += weight * (ket.conj().T @ (si @ sj) @ ket)[0, 0]
    ss = (ss / weights.sum()).real
    return ss


def _chi_prefactor(temperature_k: float) -> float:
    return _MU_B**2 * _MU_0 * _G_E / (_K_B * temperature_k)


def susceptibility_nonsymmetric(g_matrix, ss, *, temperature_k: float) -> np.ndarray:
    """chi' = (mu_B^2 mu_0 g_e / k_B T) g . <SS>  (the PDA-consistent form)."""
    g = np.asarray(g_matrix, dtype=float)
    if g.shape != (3, 3):
        raise PnmrError("the g-matrix must be a 3x3 matrix.")
    return _chi_prefactor(temperature_k) * (g @ np.asarray(ss, dtype=float))


def susceptibility_symmetric(g_matrix, ss, *, temperature_k: float) -> np.ndarray:
    """chi = (mu_B^2 mu_0 / k_B T) g . <SS> . g^T  (the textbook alternative)."""
    g = np.asarray(g_matrix, dtype=float)
    if g.shape != (3, 3):
        raise PnmrError("the g-matrix must be a 3x3 matrix.")
    ss_arr = np.asarray(ss, dtype=float)
    return (_MU_B**2 * _MU_0 / (_K_B * temperature_k)) * (g @ ss_arr @ g.T)


def susceptibility_from_g(g_matrix, *, temperature_k: float, spin: float = 0.5) -> tuple[np.ndarray, np.ndarray]:
    """The (chi', chi) pair for a spin system with the given g-tensor.

    ``chi'`` uses the non-symmetric construction (eq. (3)) and ``chi`` the
    symmetric one (eq. (5)); for a Kramers doublet (spin = 1/2) this is the
    classic chi' = (mu_B^2 mu_0 g_e / 4 k_B T) g against
    chi = (mu_B^2 mu_0 / 4 k_B T) g g^T.
    """
    g = np.asarray(g_matrix, dtype=float)
    if g.shape != (3, 3):
        raise PnmrError("the g-matrix must be a 3x3 matrix.")
    ss = np.eye(3) * (spin * (spin + 1.0) / 3.0)  # isotropic without ZFS
    return (
        susceptibility_nonsymmetric(g, ss, temperature_k=temperature_k),
        susceptibility_symmetric(g, ss, temperature_k=temperature_k),
    )


def susceptibility_from_zfs(
    g_matrix, *, D_cm1: float, E_over_D: float, temperature_k: float, spin: float
) -> tuple[np.ndarray, np.ndarray]:
    """The (chi', chi) pair with <SS> taken over the ZFS eigenstates."""
    ss = ss_dyadic(
        D_cm1=D_cm1, E_over_D=E_over_D, temperature_k=temperature_k, spin=spin
    )
    return (
        susceptibility_nonsymmetric(g_matrix, ss, temperature_k=temperature_k),
        susceptibility_symmetric(g_matrix, ss, temperature_k=temperature_k),
    )


def chi_tensors(g_matrix, *, temperature_k: float, route: str = "nonsymmetric", spin: float = 0.5,
                D_cm1: float | None = None, E_over_D: float = 0.0) -> np.ndarray:
    """The susceptibility tensor for the chosen route (no ZFS -> doublet form)."""
    if route not in ("nonsymmetric", "symmetric"):
        raise PnmrError(f"unknown route {route!r}; use nonsymmetric or symmetric.")
    if D_cm1 is None:
        chi_ns, chi_sym = susceptibility_from_g(
            g_matrix, temperature_k=temperature_k, spin=spin
        )
    else:
        chi_ns, chi_sym = susceptibility_from_zfs(
            g_matrix,
            D_cm1=D_cm1,
            E_over_D=E_over_D,
            temperature_k=temperature_k,
            spin=spin,
        )
    return chi_ns if route == "nonsymmetric" else chi_sym


def axiality(tensor) -> tuple[float, float]:
    """(Delta chi_ax, Delta chi_rh) from the symmetrized tensor's eigenvalues.

    Eigenvalues are sorted ascending (chi1 <= chi2 <= chi3):
    Delta chi_ax = chi3 - (chi1 + chi2)/2 and Delta chi_rh = chi1 - chi2 (so
    the rhombicity is <= 0; the source's magnitude convention |chi_x| <
    |chi_y| < |chi_z| differs only in labelling).  Measured convention: it
    reproduces the source's printed 14.8e-32 / 27.3e-32 from their printed
    tensors.
    """
    t = np.asarray(tensor, dtype=float)
    sym = (t + t.T) / 2
    w = np.linalg.eigvalsh(sym)
    return float(w[2] - (w[0] + w[1]) / 2), float(w[0] - w[1])


def pcs(chi_si, center, positions) -> tuple[tuple[float, float, float, float], ...]:
    """Per-position PCS in ppm: (r_angstrom, theta_deg, phi_deg, delta_ppm).

    ``chi_si`` in m^3 (SI volume), positions in Angstrom in the same frame as
    the tensor; theta/phi are reported in the eigenframe of the symmetrized
    tensor (the PCS is invariant under adding an isotropic part and under the
    antisymmetric part -- both tested).
    """
    t = np.asarray(chi_si, dtype=float)
    sym = (t + t.T) / 2
    w, v = np.linalg.eigh(sym)
    chi_traceless = sym - np.trace(sym) / 3.0 * np.eye(3)
    center = np.asarray(center, dtype=float)
    out = []
    for position in positions:
        vec = (np.asarray(position, dtype=float) - center) * 1e-10  # m
        r = float(np.linalg.norm(vec))
        if r <= 0.0:
            out.append((0.0, 0.0, 0.0, float("nan")))
            continue
        n = vec / r
        d = 3.0 * np.outer(n, n) - np.eye(3)
        delta = float(np.trace(chi_traceless @ d) / (12.0 * np.pi * r**3)) * 1e6  # ppm
        # angles in the eigenframe of the symmetrized tensor
        local = v.T @ n
        theta = float(np.degrees(np.arccos(np.clip(local[2], -1.0, 1.0))))
        phi = float(np.degrees(np.arctan2(local[1], local[0])) % 360.0)
        out.append((r * 1e10, theta, phi, delta))
    return tuple(out)


def read_run(path: str | Path) -> RunData:
    """Read a pNMR run file (YAML; schema in the formats chapter).

    The susceptibility block takes one of three forms:

    - ``from: tensor`` with ``tensor`` (3x3, m^3 -- or ``tensor_unit:
      "1e-32 m3"`` for the literature unit);
    - ``from: g`` with ``g_matrix`` (3x3) and ``spin`` (default 0.5, the
      Kramers doublet);
    - ``from: zfs`` with ``g_matrix``, ``D_cm1``, ``E_over_D`` and ``spin``.

    ``g_matrix``/``D_cm1``/``E_over_D`` may be given directly, or through
    ``orca_output: <path>`` whose QDPT blocks are parsed (the g-matrix from
    the effective-Hamiltonian block; for ``from: zfs`` the preferred ZFS
    variant: effective Hamiltonian with the spin-spin contribution when
    present).
    """
    import yaml

    from ..parsers import ParserError, parse_auto

    source = Path(path)
    try:
        raw = yaml.safe_load(source.read_text(encoding="utf-8"))
    except OSError as exc:
        raise PnmrError(f"the run file cannot be read: {exc}") from exc
    except yaml.YAMLError as exc:
        raise PnmrError(f"the run file does not parse as YAML: {exc}") from exc
    if not isinstance(raw, dict):
        raise PnmrError("the run file's top level is not a mapping.")

    from . import geometry as geometry_analysis

    structure_text = raw.get("structure")
    if structure_text is None:
        raise PnmrError(
            "the run file carries no 'structure'. Next step: give the path of an "
            "XYZ file with the molecule (the centre's coordinates come from it)."
        )
    try:
        atoms = geometry_analysis.parse_xyz(source.parent / str(structure_text))
    except Exception as exc:  # StructureError is a ValueError subclass
        raise PnmrError(f"the structure cannot be read: {exc}") from exc
    tuples = tuple((atom.element, atom.x, atom.y, atom.z) for atom in atoms)

    center_atom = raw.get("center_atom")
    if not isinstance(center_atom, int) or not (1 <= center_atom <= len(tuples)):
        raise PnmrError(
            "the run file needs 'center_atom', a 1-based index into the structure "
            f"(the structure holds {len(tuples)} atoms)."
        )
    center_index = center_atom - 1

    temperature = float(raw.get("temperature", 300.0))
    suscept = raw.get("susceptibility")
    if not isinstance(suscept, dict):
        raise PnmrError("the run file carries no 'susceptibility' mapping.")
    route = str(suscept.get("route", "nonsymmetric")).lower()
    source_kind = str(suscept.get("from", "tensor")).lower()
    notes: list[str] = []

    g_matrix = suscept.get("g_matrix")
    D_cm1 = suscept.get("D_cm1")
    E_over_D = suscept.get("E_over_D", 0.0)
    spin = float(suscept.get("spin", 0.5))
    orca_output = suscept.get("orca_output")
    if orca_output is not None:
        try:
            result = parse_auto(source.parent / str(orca_output))
        except ParserError as exc:
            raise PnmrError(f"the ORCA output cannot be parsed: {exc}") from exc
        epr = result.sections.get("epr", {})
        if g_matrix is None and epr.get("g_matrix") is not None:
            g_matrix = epr["g_matrix"]
            notes.append(
                f"the g-matrix was read from '{orca_output}' (effective-Hamiltonian "
                "block)"
            )
        if source_kind == "zfs" and D_cm1 is None:
            zfs = result.sections.get("zfs", {})
            preferred = zfs.get("preferred")
            if preferred is not None:
                D_cm1 = preferred["D_cm1"]
                E_over_D = preferred["E_over_D"] or 0.0
                notes.append(
                    f"D/E were read from '{orca_output}' (ZFS variant "
                    f"'{preferred['variant']}')"
                )

    if source_kind == "tensor":
        tensor = suscept.get("tensor")
        unit = str(suscept.get("tensor_unit", "m3")).lower()
        try:
            chi = np.asarray(tensor, dtype=float)
        except (TypeError, ValueError) as exc:
            raise PnmrError(f"the tensor is not a 3x3 matrix of numbers: {exc}") from exc
        if chi.shape != (3, 3):
            raise PnmrError("the tensor must be a 3x3 matrix.")
        if unit in ("1e-32 m3", "1e-32m3", "10^-32 m3"):
            chi = chi * 1e-32
        elif unit not in ("m3", "m^3", "si"):
            raise PnmrError(
                f"unknown tensor_unit {unit!r}; use 'm3' or '1e-32 m3'."
            )
        chi_source = "a given 3x3 tensor"
    elif source_kind in ("g", "zfs"):
        if g_matrix is None:
            raise PnmrError(
                "the susceptibility block needs a g_matrix (3x3) or an "
                "orca_output to read it from."
            )
        try:
            g_arr = np.asarray(g_matrix, dtype=float)
        except (TypeError, ValueError) as exc:
            raise PnmrError(f"the g_matrix is not a 3x3 matrix of numbers: {exc}") from exc
        if g_arr.shape != (3, 3):
            raise PnmrError("the g_matrix must be a 3x3 matrix.")
        if source_kind == "zfs":
            if D_cm1 is None:
                raise PnmrError(
                    "from: zfs needs D_cm1 (and E_over_D) -- give them, or point "
                    "orca_output at a QDPT output whose ZFS blocks carry them."
                )
            chi = chi_tensors(
                g_arr,
                temperature_k=temperature,
                route=route,
                spin=spin,
                D_cm1=float(D_cm1),
                E_over_D=float(E_over_D),
            )
            chi_source = (
                f"g + ZFS (D = {float(D_cm1):g} cm^-1, E/D = {float(E_over_D):g}, "
                f"spin {spin:g})"
            )
        else:
            chi = chi_tensors(
                g_arr, temperature_k=temperature, route=route, spin=spin
            )
            chi_source = f"g-tensor (spin {spin:g})"
    else:
        raise PnmrError(
            f"unknown susceptibility source {source_kind!r}; use tensor, g or zfs."
        )

    nuclei_raw = raw.get("nuclei", "all")
    if nuclei_raw in ("all", None):
        nuclei: tuple[int, ...] = ()
    elif isinstance(nuclei_raw, list):
        nuclei_list = []
        for value in nuclei_raw:
            if not isinstance(value, int) or not (1 <= value <= len(tuples)):
                raise PnmrError(
                    f"nuclei entry {value!r} is not a 1-based index into the structure."
                )
            nuclei_list.append(value - 1)
        nuclei = tuple(nuclei_list)
    else:
        raise PnmrError("'nuclei' must be 'all' or a list of 1-based indices.")

    return RunData(
        atoms=tuples,
        center_index=center_index,
        temperature_k=temperature,
        route=route,
        spin=spin,
        chi_si=chi,
        chi_source=chi_source,
        nuclei=nuclei,
        structure_name=str(structure_text),
        notes=tuple(notes),
    )


def render(data: RunData) -> str:
    """The menu's report: the susceptibility, the axiality and the PCS table."""
    chi = data.chi_si
    d_ax, d_rh = axiality(chi)
    sym_part = (chi + chi.T) / 2
    asym = float(np.abs(chi - chi.T).max())
    out = [
        "pNMR pseudocontact shifts (point-dipole approximation)",
        f"  structure: {data.structure_name} ({len(data.atoms)} atoms); centre: "
        f"atom {data.center_index + 1} ({data.atoms[data.center_index][0]})",
        f"  susceptibility from {data.chi_source}; route: {data.route}; "
        f"T = {data.temperature_k:g} K",
        f"  chi tensor (1e-32 m^3):",
    ]
    for row in chi:
        out.append(
            "    " + "  ".join(f"{value / 1e-32:>10.4f}" for value in row)
        )
    out.append(
        f"  Delta chi_ax = {d_ax / 1e-32:.4f}e-32 m^3; "
        f"Delta chi_rh = {d_rh / 1e-32:.4f}e-32 m^3; "
        f"|antisymmetric|_max = {asym / 1e-32:.3e}e-32 m^3"
    )
    for note in data.notes:
        out.append(f"  note: {note}")
    center = data.atoms[data.center_index][1:]
    indices = [
        index
        for index in range(len(data.atoms))
        if index != data.center_index and (not data.nuclei or index in data.nuclei)
    ]
    positions = [data.atoms[index][1:] for index in indices]
    values = pcs(chi, center, positions)
    out.append("")
    out.append("  nucleus        r (A)     theta (deg)   phi (deg)    delta_PCS (ppm)")
    for index, (r, theta, phi, delta) in zip(indices, values):
        label = f"{data.atoms[index][0]}{index + 1}"
        out.append(
            f"  {label:<10s} {r:>9.3f} {theta:>13.1f} {phi:>11.1f} {delta:>16.3f}"
        )
    out += [
        "",
        "boundaries (measured and from the sources):",
        "  - the PDA is the long-range limit: the source's benchmark shows <10 % "
        "deviation of the dipolar part beyond ~7 A from the centre, while the "
        "contact (through-bond) contribution dominates below ~4-5 A -- a PCS "
        "value predicted for a nucleus close to the metal is a lower bound on "
        "the total shift error, not a full prediction;",
        "  - only the traceless symmetric part of chi enters the PCS (the "
        "antisymmetric part contributes nothing; tested exactly) -- the "
        "route's difference (nonsymmetric vs symmetric construction) changes "
        "the symmetric part itself and hence the values;",
        "  - the ZFS route uses the effective spin dyadic <SS> at the given "
        "temperature; near-axial systems can cross the node condition "
        "<SS>|| g|| = <SS>perp gperp where the axiality (and the PCS) vanish;",
        "  - magnetic-dipole-free, isotropic-shift contexts only: the contact "
        "shift is not computed here (it needs hyperfine coupling).",
    ]
    return "\n".join(out)


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the PDA form, the susceptibility constructions and the survey."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "The point-dipole pseudocontact-shift form and the two "
                "susceptibility constructions (non-symmetric chi' = "
                "mu_B^2 mu_0 g_e/(k_B T) g.<SS>, consistent with the modern "
                "pNMR shielding theory, versus the symmetric g.<SS>.g^T "
                "alternative) follow the source exactly; its printed Co(II) "
                "benchmark tensors (Delta chi_ax = 14.8 / 27.3 x 1e-32 m^3) "
                "are this module's regression anchors."
            ),
            ref=(
                "Mares J., Vaara J., Phys. Chem. Chem. Phys. 2018, 20, 22547-22555, "
                "DOI 10.1039/c8cp04123g"
            ),
            url="https://doi.org/10.1039/c8cp04123g",
            bibkey="mares2018pcs",
        ),
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "The ligand-field context and the limits of Bleaney's "
                "anisotropy theory (PCS treated as a point susceptibility, "
                "validity broken when the ligand-field splitting exceeds kT; "
                "spin delocalization and tag mobility as the two reasons for "
                "distributed 4f density) frame the boundaries this menu reports."
            ),
            ref=(
                "Parker D., Suturina E. A., Kuprov I., Chilton N. F., Acc. Chem. "
                "Res. 2020, 53, 1522-1534, DOI 10.1021/acs.accounts.0c00275"
            ),
            url="https://doi.org/10.1021/acs.accounts.0c00275",
            bibkey="parker2020ligandfield",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "Measured on ORCA 6.1.1: the QDPT driver prints the "
                "effective-Hamiltonian g-matrix, the ZFS blocks (2nd-order / "
                "effective-Hamiltonian, with and without the SSC contribution) "
                "and, with DoSusceptibility, the molar susceptibility tensor per "
                "temperature in cm^3 K/mol; the module's chi from the parsed "
                "CO+ g reproduces ORCA's own printed susceptibility to the "
                "printed precision (4 pi conversion pinned in the tests)."
            ),
            ref="tests/test_pnmr.py; fixtures/orca/co_plus_qdpt.*, o2_qdpt2.*",
        ),
    )
