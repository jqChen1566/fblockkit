"""OpenMolcas SINGLE_ANISO output: the g tensors of the pseudospin multiplets.

The ``SINGLE_ANISO`` module reads the RASSI spin-orbit data and prints, for every
requested pseudospin multiplet, the effective spin, the constituting spin-orbit
states with their relative energies, the tunnelling splitting, the g tensor
(three principal values with their main magnetic axes in the initial Cartesian
frame), and the sign check of the g-value product.  The format below was
measured on OpenMolcas v26.06 (fixture ``fixtures/openmolcas/dy_smoke.out``, the
Dy(III) 4f9 single-ion chain; the module prints its own banner as
``SINGLE_ANISO (OPEN)``).

Parsed per multiplet:

- the multiplet index and the effective spin (as printed, e.g. ``1/2``);
- the spin-orbit states that span it, with energies in cm^-1;
- the tunnelling splitting in cm^-1;
- the three principal g values and the three main magnetic axes (unit vectors
  in the initial Cartesian frame; ``Xm``/``Ym``/``Zm`` pair with
  ``gX``/``gY``/``gZ``);
- the ``CHECK-SIGN`` parameter and the sign of the product gX*gY*gZ.

The irreducible-tensor (ITO/ESO) decomposition tables that follow the g tensor,
and the module's other property blocks (susceptibility, magnetisation), are not
parsed by this reader.

``doublets_payload`` converts the parsed multiplets into the JSON table menu 16
consumes (``{"system", "reference", "doublets": [{"label", "g", "energy",
"axis3"}]}``): the g values are re-ordered ascending (menu 16's g1 <= g2 <= g3
convention, with each axis following its value), and every doublet carries the
axis of its largest g value, so the criterion measures each multiplet's angle
against the reference one.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

__all__ = [
    "MolcasAnisoError",
    "Multiplet",
    "doublets_payload",
    "parse_single_aniso",
]


class MolcasAnisoError(ValueError):
    """Refusal raised when the text carries no parseable SINGLE_ANISO block."""


_FLOAT = r"-?\d+\.\d+(?:[eE][+-]?\d+)?"

_HEADER_RE = re.compile(
    r"CALCULATION OF PSEUDOSPIN HAMILTONIAN TENSORS FOR THE MULTIPLET\s+(\d+)\s*"
    r"\(\s*effective S\s*=\s*([0-9./]+)\s*\)"
)
_SO_STATE_RE = re.compile(
    rf"spin-orbit state\s+(\d+);\s+energy\((\d+)\)\s*=\s*({_FLOAT})\s*cm-1"
)
_TUNNEL_RE = re.compile(rf"Tunnelling splitting:\s*({_FLOAT})\s*cm-1")
_G_LINE_RE = re.compile(
    rf"\bg([XYZ])\s*=\s*({_FLOAT})\s*\|\s*([XYZ])m\s*\|\s*"
    rf"({_FLOAT})\s+({_FLOAT})\s+({_FLOAT})\s*\|"
)
_CHECK_SIGN_RE = re.compile(rf"CHECK-SIGN parameter\s*=\s*({_FLOAT})")
_SIGN_RE = re.compile(
    r"The sign of the product gX \* gY \* gZ for multiplet\s+(\d+):\s*([<>])\s*0"
)


@dataclass(frozen=True)
class Multiplet:
    """One pseudospin multiplet block of a SINGLE_ANISO output."""

    index: int
    effective_spin: str
    so_states: tuple[tuple[int, float], ...]
    tunnelling_splitting_cm1: float | None
    g_values: tuple[float, float, float]
    axes: tuple[
        tuple[float, float, float],
        tuple[float, float, float],
        tuple[float, float, float],
    ]
    check_sign: float | None
    sign_product: int | None


def parse_single_aniso(text: str) -> tuple[Multiplet, ...]:
    """Parse every pseudospin multiplet block in a SINGLE_ANISO output.

    Blocks are located by their ``CALCULATION OF PSEUDOSPIN ...`` headers and read
    up to the next header (or the end of the text).  Each block must carry all
    three g-tensor lines; a partial table is refused rather than partially read.
    """
    headers = list(_HEADER_RE.finditer(text))
    if not headers:
        raise MolcasAnisoError(
            "no 'CALCULATION OF PSEUDOSPIN HAMILTONIAN TENSORS' block was found. "
            "Next step: give the output of an OpenMolcas SINGLE_ANISO run -- the "
            "printed 'g TENSOR' table comes from this block."
        )
    multiplets: list[Multiplet] = []
    for position, header in enumerate(headers):
        end = headers[position + 1].start() if position + 1 < len(headers) else len(text)
        block = text[header.start():end]
        index = int(header.group(1))
        effective_spin = header.group(2)
        states = tuple(
            (int(m.group(1)), float(m.group(3))) for m in _SO_STATE_RE.finditer(block)
        )
        tunnel = _TUNNEL_RE.search(block)
        g_by_component: dict[str, float] = {}
        axis_by_component: dict[str, tuple[float, float, float]] = {}
        for m in _G_LINE_RE.finditer(block):
            g_by_component[m.group(1)] = float(m.group(2))
            axis_by_component[m.group(1)] = (
                float(m.group(4)),
                float(m.group(5)),
                float(m.group(6)),
            )
        if sorted(g_by_component) != ["X", "Y", "Z"]:
            raise MolcasAnisoError(
                f"multiplet {index} carries an incomplete g-tensor table "
                f"(found the {sorted(g_by_component)} lines). Next step: check that "
                "the pseudo-spin format was requested (MLTP) for this multiplet -- "
                "the module prints no g tensor for blocks it did not process."
            )
        check_sign = _CHECK_SIGN_RE.search(block)
        sign = _SIGN_RE.search(block)
        multiplets.append(
            Multiplet(
                index=index,
                effective_spin=effective_spin,
                so_states=states,
                tunnelling_splitting_cm1=(
                    float(tunnel.group(1)) if tunnel is not None else None
                ),
                g_values=(
                    g_by_component["X"],
                    g_by_component["Y"],
                    g_by_component["Z"],
                ),
                axes=(
                    axis_by_component["X"],
                    axis_by_component["Y"],
                    axis_by_component["Z"],
                ),
                check_sign=float(check_sign.group(1)) if check_sign is not None else None,
                sign_product=(
                    (1 if sign.group(2) == ">" else -1) if sign is not None else None
                ),
            )
        )
    return tuple(multiplets)


def doublets_payload(
    multiplets: tuple[Multiplet, ...], system: str | None = None
) -> dict:
    """Build the menu-16 JSON table from parsed multiplets (see module docstring).

    The first multiplet is the reference (its theta3 is 0 by definition); every
    row carries the axis of its largest g value so the criterion can measure the
    remaining multiplets against it.
    """
    if not multiplets:
        raise MolcasAnisoError(
            "no multiplets to convert. Next step: parse a SINGLE_ANISO output first."
        )
    rows = []
    for multiplet in multiplets:
        triples = sorted(zip(multiplet.g_values, multiplet.axes), key=lambda t: t[0])
        lowest_energy = (
            min(energy for _, energy in multiplet.so_states)
            if multiplet.so_states
            else None
        )
        row: dict = {
            "label": f"multiplet {multiplet.index} (effective S = {multiplet.effective_spin})",
            "g": [value for value, _ in triples],
            "axis3": list(triples[-1][1]),
        }
        if lowest_energy is not None:
            row["energy"] = lowest_energy
        rows.append(row)
    return {
        "system": system or "OpenMolcas SINGLE_ANISO output",
        "reference": 0,
        "doublets": rows,
    }
