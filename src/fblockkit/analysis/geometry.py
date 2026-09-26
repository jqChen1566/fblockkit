"""S1: coordination geometry and symmetry hints (read a structure, no calculation first).

Scope statement: this module gives **geometric-distribution hints** (coordination
number, direction distribution, inversion pairing, planarity), not a rigorous
point-group determination -- the point group is confirmed by the user, and the
"symmetry -> number of CF parameters" rule is applied to that.

CF parameter-count rule (from the internal close-reading record
"文献细读/细读_CF哈密顿量_Chan2025.md" and the source paper behind it): keep only even
k (time-reversal symmetry) and k <= 6 (4f: l = 3); a true C3 fit has 9
parameters and an Oh fit 4; a fit that needs 27 means the frame carries no real
symmetry (the z-axis convention must then be recorded).

Dependencies: the standard library only (3x3 symmetric matrices use analytic
eigenvalues, so numpy is not introduced).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from ..knowledge.elements import ElementError, element_z, is_f_element
from ..knowledge.models import EVIDENCE_LITERATURE, Evidence, ReportSection

_EVIDENCE = Evidence(
    kind=EVIDENCE_LITERATURE,
    text=(
        "The number of non-zero parameters in a CF fit is fixed by the point group: a true "
        "C3 system fits 9 B_k^q, Oh fits 4, and with no true symmetry the fit takes 27 "
        "(the z-axis convention must then be declared). Only even k is kept (time reversal) "
        "and k <= 6 (4f, l = 3)."
    ),
    ref="Peng L., Liu S., Zhang X., Chen X., Li C., Ung S. F., Cheng H.-P., Chan G. K.-L., J. Phys. Chem. Lett., 2025, 16(47), 12312-12320, DOI 10.1021/acs.jpclett.5c02971",
    bibkey="peng2025accurate",
    url="https://doi.org/10.1021/acs.jpclett.5c02971",
)

CF_PARAM_COUNTS = {"C3": 9, "Oh": 4, "none": 27}


class StructureError(ValueError):
    """Invalid structure input."""


@dataclass(frozen=True)
class Atom:
    element: str
    x: float
    y: float
    z: float

    def distance_to(self, other: "Atom") -> float:
        return math.dist((self.x, self.y, self.z), (other.x, other.y, other.z))


def parse_xyz(source: str | Path) -> tuple[Atom, ...]:
    """Parse XYZ text or a file (coordinates only; the comment line is ignored)."""
    text = source if isinstance(source, str) and "\n" in source else None
    if text is None:
        path = Path(source)
        if not path.exists():
            raise StructureError(f"file does not exist: {path}. Next step: check the path.")
        text = path.read_text(encoding="utf-8", errors="replace")
    lines = [line for line in text.splitlines() if line.strip()]
    if not lines:
        raise StructureError(
            "the structure is empty. Next step: give coordinates in XYZ format."
        )
    try:
        declared = int(lines[0].split()[0])
    except ValueError:
        raise StructureError(
            f"the first line is not an atom count: {lines[0]!r}. Next step: use standard "
            f"XYZ (atom count on the first line, a comment on the second)."
        ) from None
    atoms: list[Atom] = []
    for line in lines[2:]:
        fields = line.split()
        if len(fields) < 4:
            continue
        try:
            element_z(fields[0])
            x, y, z = (float(value) for value in fields[1:4])
        except (ElementError, ValueError) as exc:
            raise StructureError(f"invalid coordinate line: {line!r} ({exc})") from exc
        atoms.append(Atom(fields[0].capitalize(), x, y, z))
    if len(atoms) != declared:
        raise StructureError(
            f"{declared} atoms were declared but {len(atoms)} were parsed. "
            f"Next step: check the structure file."
        )
    return tuple(atoms)


def dominant_center(atoms: Sequence[Atom]) -> int:
    """Pick the central atom: the first Ln/An, otherwise the largest Z."""
    for index, atom in enumerate(atoms):
        if is_f_element(atom.element):
            return index
    heaviest = max(range(len(atoms)), key=lambda i: element_z(atoms[i].element))
    return heaviest


def coordination_shell(
    atoms: Sequence[Atom], center: int, max_cn: int = 12, gap_ratio: float = 1.25
) -> tuple[list[tuple[int, float]], float]:
    """Determine the coordination shell by distance sorting plus the largest relative
    gap; returns ([(atom index, distance)], cutoff distance)."""
    pairs = sorted(
        ((index, atoms[center].distance_to(atom)) for index, atom in enumerate(atoms) if index != center),
        key=lambda item: item[1],
    )
    window = pairs[:max_cn]
    cutoff = window[-1][1] if window else 0.0
    best_ratio, best_at = 0.0, len(window)
    for i in range(len(window) - 1):
        ratio = window[i + 1][1] / window[i][1] if window[i][1] > 0 else 0.0
        if ratio > best_ratio:
            best_ratio, best_at = ratio, i + 1
    if best_ratio >= gap_ratio:
        cutoff = (window[best_at - 1][1] + window[best_at][1]) / 2.0
        shell = [pair for pair in pairs if pair[1] <= cutoff]
    else:
        shell = window
    return shell, cutoff


def _symmetric_eigenvalues(matrix: Sequence[Sequence[float]]) -> tuple[float, float, float]:
    """Analytic eigenvalues of a 3x3 symmetric matrix (descending)."""
    (a, b, c), (_, d, e), (_, _, f) = matrix
    trace = a + d + f
    q = trace / 3.0
    p2 = (a - q) ** 2 + (d - q) ** 2 + (f - q) ** 2 + 2.0 * (b * b + c * c + e * e)
    if p2 <= 0.0:
        return (q, q, q)
    p = math.sqrt(p2 / 6.0)
    # B = (M - qI)/p; det(B)/2
    b11, b12, b13 = (a - q) / p, b / p, c / p
    b22, b23 = (d - q) / p, e / p
    b33 = (f - q) / p
    det_b = (
        b11 * (b22 * b33 - b23 * b23)
        - b12 * (b12 * b33 - b23 * b13)
        + b13 * (b12 * b23 - b22 * b13)
    )
    r = max(-1.0, min(1.0, det_b / 2.0))
    phi = math.acos(r) / 3.0
    eig1 = q + 2.0 * p * math.cos(phi)
    eig3 = q + 2.0 * p * math.cos(phi + 2.0 * math.pi / 3.0)
    eig2 = 3.0 * q - eig1 - eig3
    return (eig1, eig2, eig3)


def _shape_label(eigenvalues: Sequence[float]) -> str:
    largest = eigenvalues[0] if eigenvalues[0] > 0 else 1.0
    r2, r3 = eigenvalues[1] / largest, eigenvalues[2] / largest
    if r2 < 0.05 and r3 < 0.05:
        return "nearly linear distribution"
    if r3 < 0.05:
        return "nearly planar distribution (the ligands are nearly coplanar)"
    if r2 > 0.8 and r3 > 0.8:
        return "direction distribution nearly isotropic (octahedral-style high-symmetry coordination, or nearly spherical)"
    if r2 < 0.6 and r3 < 0.3:
        return "nearly axial distribution (weak transverse component)"
    return "low-symmetry distribution"


def run(
    source: str | Path, *, center: int | None = None, cf_point_group: str | None = None
) -> ReportSection:
    """Build the S1 report section. ``cf_point_group`` is passed in once the user has
    confirmed it, and is used to give the number of CF parameters."""
    atoms = parse_xyz(source)
    center_index = dominant_center(atoms) if center is None else center
    center_atom = atoms[center_index]
    shell, cutoff = coordination_shell(atoms, center_index)

    vectors = []
    for index, distance in shell:
        atom = atoms[index]
        vectors.append(
            (
                (atom.x - center_atom.x) / distance,
                (atom.y - center_atom.y) / distance,
                (atom.z - center_atom.z) / distance,
            )
        )
    count = len(vectors) or 1
    moment = [[0.0] * 3 for _ in range(3)]
    for vector in vectors:
        for i in range(3):
            for j in range(3):
                moment[i][j] += vector[i] * vector[j] / count
    eigenvalues = _symmetric_eigenvalues(moment)
    shape = _shape_label(eigenvalues)

    inversion_pairs = 0
    for i, (index_i, distance_i) in enumerate(shell):
        vector_i = vectors[i]
        for j, (index_j, distance_j) in enumerate(shell):
            if j <= i:
                continue
            vector_j = vectors[j]
            cosine = sum(a * b for a, b in zip(vector_i, vector_j))
            if cosine < -math.cos(math.radians(30)) and abs(distance_i - distance_j) / distance_i < 0.15:
                inversion_pairs += 1
                break

    lines = [
        f"Structure: {len(atoms)} atoms; central atom {center_atom.element}{center_index + 1}",
        f"Coordination shell: {len(shell)} ligands (cutoff distance {cutoff:.2f} Å)",
        "  Ligands and distances (Å): "
        + ", ".join(
            f"{atoms[index].element}{index + 1} {distance:.2f}" for index, distance in shell
        ),
        "",
        "Symmetry hints (geometric distribution, not a rigorous point-group determination):",
        f"  - Direction distribution: {shape} (second-moment eigenvalue ratios "
        f"λ2/λ1 = {eigenvalues[1] / eigenvalues[0]:.2f}, "
        f"λ3/λ1 = {eigenvalues[2] / eigenvalues[0]:.2f})",
        f"  - Inversion pairs: {inversion_pairs}/{len(shell)} (pairs of ligands in an approximately centrosymmetric arrangement)",
        "",
        "Number of CF parameters (applies once the point group is confirmed): "
        + "; ".join(f"{group} → {count_value} B_k^q" for group, count_value in CF_PARAM_COUNTS.items()),
        "  Rule: only even k is kept (time-reversal symmetry), k <= 6 (4f); a fit that needs 27 means the frame used carries no true symmetry -- the z-axis convention must then be declared.",
    ]
    if cf_point_group:
        known = CF_PARAM_COUNTS.get(cf_point_group)
        lines.append(
            f"  Point group confirmed by the user, {cf_point_group}: "
            + (
                f"{known} B_k^q allowed."
                if known
                else "this table has no parameter count for that point group -- work it out from the rule above."
            )
        )
    return ReportSection(title="S1 coordination geometry and symmetry hints", body="\n".join(lines))


def evidence() -> tuple[Evidence, ...]:
    return (_EVIDENCE,)
