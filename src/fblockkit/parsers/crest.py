"""The CREST conformer ensemble (plan item 6.1; the reader behind the ensemble menu).

CREST (Pracht, Bohle, Grimme, Phys. Chem. Chem. Phys. 2020, 22, 7169;
LGPL-3.0) samples the low-energy conformer-rotamer space on the GFN2-xTB
surface.  This module reads the final ensemble of a run back; every
convention below is **measured** on CREST 3.0.2 (fixtures ``fixtures/crest/``,
2026-09-30, n-butane GFN2):

- ``crest_conformers.xyz``: one XYZ frame per conformer, in ensemble order;
  the atom-count line carries **leading spaces** (``  14``); the comment
  line is the **absolute energy in hartree**, right-aligned in a wide field
  (``        -13.66512773``); atom lines carry leading spaces too;
- ``crest.energies``: ``<index>  <relative energy in kcal/mol>`` rows, index
  one-based, energies three decimals (``  1         0.000``);
- ``crest_best.xyz``: the single lowest conformer, same frame convention.

The reader cross-checks the two sources: the relative energies implied by
the frame energies must reproduce ``crest.energies`` within the print
precision (kcal/mol against hartree at the CODATA conversion).  A run whose
input geometry was not pre-optimised may abort after a topology-change
warning (the capture offers options A/B/C; measured on the first probe,
2026-09-30) -- such a run leaves no ensemble and is handled upstream.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

__all__ = [
    "CrestError",
    "CrestFrame",
    "CrestEnsemble",
    "read_crest_ensemble",
]

_HARTREE_TO_KCAL_PER_MOL = 627.509474  # CODATA-derived; matches crest's print precision


class CrestError(ValueError):
    """A CREST artifact is missing or of an unrecognised layout (with a next step)."""


@dataclass(frozen=True)
class CrestFrame:
    """One conformer frame: absolute energy plus the atom block."""

    index: int  # 1-based ensemble position
    energy_Eh: float
    atoms: tuple[tuple[str, float, float, float], ...]


@dataclass(frozen=True)
class CrestEnsemble:
    """The measured ensemble: frames, the relative-energy table, the best frame."""

    frames: tuple[CrestFrame, ...]
    relative_kcal: tuple[float, ...]  # aligned with frames (from crest.energies)
    best_index: int | None           # 1-based frame index of the best conformer
    max_rel_deviation_kcal: float    # cross-check: |relative - (E-E_min)*conv|


def read_crest_ensemble(
    conformers_text: str,
    energies_text: str | None = None,
    best_text: str | None = None,
) -> CrestEnsemble:
    """Parse the CREST ensemble files (conformers required; energies optional)."""
    frames = _xyz_frames(conformers_text)
    if not frames:
        raise CrestError(
            "no conformer frames found: the first frame of a CREST ensemble must "
            "start with the atom count line. Next step: give the run directory that "
            "contains crest_conformers.xyz."
        )
    relative: tuple[float, ...]
    if energies_text is not None:
        relative = _relative_energies(energies_text, expected=len(frames))
    else:
        relative = tuple(
            (frame.energy_Eh - frames[0].energy_Eh) * _HARTREE_TO_KCAL_PER_MOL for frame in frames
        )
    minimum = min(frame.energy_Eh for frame in frames)
    deviation = max(
        abs(rel - (frame.energy_Eh - minimum) * _HARTREE_TO_KCAL_PER_MOL)
        for rel, frame in zip(relative, frames)
    )
    best_index: int | None = None
    if best_text is not None:
        best_frames = _xyz_frames(best_text)
        if best_frames:
            target = best_frames[0].energy_Eh
            closest = min(
                range(len(frames)), key=lambda position: abs(frames[position].energy_Eh - target)
            )
            if abs(frames[closest].energy_Eh - target) < 1e-6:
                best_index = closest + 1
    return CrestEnsemble(
        frames=frames,
        relative_kcal=relative,
        best_index=best_index,
        max_rel_deviation_kcal=deviation,
    )


def read_crest_directory(directory: str | Path) -> CrestEnsemble:
    """Read a run directory: crest_conformers.xyz (required), crest.energies, crest_best.xyz."""
    root = Path(directory)
    if root.is_file():
        root = root.parent
    conformers = root / "crest_conformers.xyz"
    if not conformers.is_file():
        raise CrestError(
            f"{conformers} not found. Next step: give the CREST run directory (the one "
            "holding crest_conformers.xyz), or the crest_conformers.xyz path itself."
        )
    energies = root / "crest.energies"
    best = root / "crest_best.xyz"
    return read_crest_ensemble(
        conformers.read_text(encoding="utf-8", errors="replace"),
        energies.read_text(encoding="utf-8", errors="replace") if energies.is_file() else None,
        best.read_text(encoding="utf-8", errors="replace") if best.is_file() else None,
    )


def _xyz_frames(text: str) -> tuple[CrestFrame, ...]:
    lines = text.splitlines()
    frames: list[CrestFrame] = []
    position = 0
    while position < len(lines):
        head = lines[position].strip()
        if not head:
            position += 1
            continue
        try:
            atom_count = int(head)
        except ValueError:
            raise CrestError(
                f"line {position + 1} should be the atom count but reads {head!r}; "
                "this is not a CREST ensemble frame."
            ) from None
        if position + 1 >= len(lines):
            break
        comment = lines[position + 1].strip()
        try:
            energy = float(comment.split()[0])
        except (ValueError, IndexError):
            raise CrestError(
                f"frame comment on line {position + 2} should carry the energy in "
                f"hartree but reads {comment!r}."
            ) from None
        atoms: list[tuple[str, float, float, float]] = []
        for offset in range(atom_count):
            row = lines[position + 2 + offset].split()
            atoms.append((row[0], float(row[1]), float(row[2]), float(row[3])))
        frames.append(
            CrestFrame(index=len(frames) + 1, energy_Eh=energy, atoms=tuple(atoms))
        )
        position += 2 + atom_count
    return tuple(frames)


def _relative_energies(text: str, expected: int) -> tuple[float, ...]:
    rows: dict[int, float] = {}
    for line in text.splitlines():
        fields = line.split()
        if len(fields) != 2:
            continue
        try:
            rows[int(fields[0])] = float(fields[1])
        except ValueError:
            continue
    missing = [index for index in range(1, expected + 1) if index not in rows]
    if missing:
        raise CrestError(
            f"crest.energies covers {len(rows)} of the {expected} conformers "
            f"(missing index {missing[0]}); the ensemble and the table disagree."
        )
    return tuple(rows[index] for index in range(1, expected + 1))
