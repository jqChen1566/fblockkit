"""4.5: displacing a geometry along its imaginary mode (menu 30).

Source: NIFREC (Tanaka & Miyao), the automated no-imaginary-frequency
workflow (ChemRxiv preprint 10.26434/chemrxiv.15004859/v1; the software is
MIT-licensed, ``github.com/tanaka-hideya/NIFREC``; the algorithm was read
from its source, 2026-09-27 -- no code is copied).  Its displacement stage:
build a vector from the imaginary modes -- the **sum of all imaginary
modes** by default, or the **most negative mode alone** -- normalize it, and
apply ``coords + vec * disp`` with successively larger scalar displacements
``base_disp * (i + 1)`` (defaults: ``base_disp`` 0.1 Angstrom, up to five
rounds, i.e. 0.1 ... 0.5), rerunning ``opt freq`` until every frequency is
real.  Their success criterion is exactly that: no imaginary frequency in
the rerun.

The ORCA side of the same idea: ORCA prints, after the frequency list, the
``NORMAL MODES`` block -- "the Cartesian displacements weighted by the
diagonal matrix M(i,i)=1/sqrt(m[i])", unit Euclidean norm (measured on the
F + H2 fixture: the printed columns have norm 1.000000).  The Cartesian
displacement pattern is therefore ``v_i / sqrt(m_i)``, and the masses are
taken from the run's own ``CARTESIAN COORDINATES (A.U.)`` block (so the
de-weighting uses the values ORCA's Hessian was built with).  The geometry
displaced is the run's **final** structure -- the one the modes belong to.

Two deviations from the source, both recorded: the menu writes **both signs**
of the displacement (the source only grows one direction) so that both
branches of a double-well mode are available in one round; and the menu is a
generator, not a driver -- the escalation schedule is printed as guidance
and the user runs the jobs.

The job type of the given base input is kept verbatim (an ``OptTS`` run stays
mode-following, an ``Opt`` run heals toward a minimum); the menu only
replaces the coordinate block with the displaced structure and makes sure the
``Freq`` verification is in the simple line.
"""

from __future__ import annotations

import re

import numpy as np

from ..knowledge.elements import ElementError, element_z
from ..knowledge.models import EVIDENCE_LITERATURE, EVIDENCE_MEASURED, Evidence

__all__ = [
    "ImagDispError",
    "imaginary_modes",
    "displacement_vector",
    "displaced_atoms",
    "disp_input",
    "render",
    "evidence",
]


class ImagDispError(ValueError):
    """The displacement cannot be built (with a next step)."""


def imaginary_modes(result) -> dict:
    """The final frequency list, its normal-mode vectors, geometry and masses.

    ``result`` is an :func:`fblockkit.parsers.parse_auto` result.  The last
    frequency block and the last normal-mode block must agree about the
    structure (the parser's measured convention); the geometry and masses
    come from the final ``(ANGSTROEM)``/``(A.U.)`` pair.
    """
    frequencies = result.sections.get("frequencies", {})
    normal_modes = result.sections.get("normal_modes", {})
    geometry = result.sections.get("final_geometry", {})
    if not frequencies.get("present"):
        raise ImagDispError(
            "the output carries no VIBRATIONAL FREQUENCIES block. Next step: run "
            "the frequency job first (or give the output that has it)."
        )
    if not geometry.get("present"):
        raise ImagDispError(
            "the output carries no final CARTESIAN COORDINATES block, so the "
            "geometry the modes belong to is unknown."
        )
    if not normal_modes.get("present"):
        raise ImagDispError(
            "the output carries no NORMAL MODES block (needed for the mode vectors)."
        )
    block = frequencies["blocks"][-1]
    imaginary = block["imaginary"]
    if not imaginary:
        raise ImagDispError(
            "the last frequency block has no imaginary mode (nothing to displace "
            "along); the geometry is already a minimum in this metric."
        )
    if not geometry.get("masses"):
        raise ImagDispError(
            "the output's CARTESIAN COORDINATES (A.U.) block does not mirror the "
            "(ANGSTROEM) one, so the atomic masses ORCA used are unavailable."
        )
    vectors = {mode["index"]: mode["vector"] for mode in normal_modes["modes"]}
    missing = [mode["index"] for mode in imaginary if mode["index"] not in vectors]
    if missing:
        raise ImagDispError(
            f"the imaginary mode(s) {missing} have no vector in the NORMAL MODES "
            "block (the last frequency block and the last mode block disagree)."
        )
    atoms = geometry["atoms"]
    if normal_modes["n_coord"] != 3 * len(atoms):
        raise ImagDispError(
            f"the normal-mode vectors span {normal_modes['n_coord']} coordinates "
            f"while the final geometry holds {len(atoms)} atom(s) (3N coordinates)."
        )
    return {
        "imaginary": imaginary,
        "vectors": vectors,
        "atoms": atoms,
        "masses": geometry["masses"],
        "scaling": block["scaling"],
    }


def displacement_vector(data, *, selection: str = "sum") -> tuple[float, ...]:
    """The unit-norm Cartesian displacement pattern of the chosen imaginary modes.

    ``selection='sum'``: the sum of every imaginary mode's vector (the source's
    default); ``'lowest'``: the most negative mode alone.  Each printed vector
    is de-weighted to Cartesian (``v_i / sqrt(m_i)``) and the result is
    normalized to unit 2-norm.
    """
    if selection not in ("sum", "lowest"):
        raise ImagDispError(f"unknown mode selection {selection!r}; use sum or lowest.")
    imaginary = data["imaginary"]
    if not imaginary:
        raise ImagDispError(
            "no imaginary mode was selected (the frequency data carry none). Next "
            "step: run the frequency job and pass its output."
        )
    if selection == "sum":
        chosen = [mode["index"] for mode in imaginary]
    else:
        chosen = [min(imaginary, key=lambda mode: mode["wavenumber"])["index"]]
    vector = np.zeros(len(data["atoms"]) * 3, dtype=float)
    for index in chosen:
        vector += np.asarray(data["vectors"][index], dtype=float)
    masses = np.repeat(np.asarray(data["masses"], dtype=float), 3)
    if np.any(masses <= 0.0):
        raise ImagDispError("a non-positive atomic mass appeared in the (A.U.) block.")
    cartesian = vector / np.sqrt(masses)
    norm = float(np.linalg.norm(cartesian))
    if norm <= 1e-12:
        raise ImagDispError(
            "the chosen imaginary mode vector(s) de-weight to a vanishing vector."
        )
    return tuple(float(value) for value in cartesian / norm)


def displaced_atoms(atoms, vector, amplitude: float):
    """``coords + vector * amplitude`` (both signs are the caller's business).

    Returns ``(atoms, max_atomic_displacement)``; the unit vector's largest
    component times the amplitude is the largest single-atom move.
    """
    if amplitude <= 0.0:
        raise ImagDispError("the displacement amplitude must be positive.")
    values = np.asarray(vector, dtype=float)
    if values.shape != (3 * len(atoms),):
        raise ImagDispError(
            f"the vector spans {values.shape[0]} coordinates while the geometry "
            f"holds {3 * len(atoms)}."
        )
    moved = []
    for index, (symbol, x, y, z) in enumerate(atoms):
        dx, dy, dz = values[3 * index: 3 * index + 3] * amplitude
        moved.append((symbol, x + dx, y + dy, z + dz))
    return tuple(moved), float(np.abs(values).max() * amplitude)


_XYZ_HEADER_RE = re.compile(r"^\s*\*\s*xyz\s+[-\d]", re.IGNORECASE)
_FREQ_TOKEN_RE = re.compile(r"\bfreq\b", re.IGNORECASE)


def disp_input(base_text: str, atoms, *, ensure_freq: bool = True) -> str:
    """The base input with its inline ``* xyz`` block replaced by ``atoms``.

    The charge/multiplicity of the base block are kept; the ``!`` lines are
    kept verbatim (the job type decides whether the rerun follows the mode or
    heals to a minimum).  With ``ensure_freq`` the first simple line gains the
    ``Freq`` token when it has none -- the source's success criterion is the
    frequency check of the rerun, so a job without it cannot verify anything.
    """
    lines = base_text.splitlines()
    start = next((i for i, line in enumerate(lines) if _XYZ_HEADER_RE.match(line)), None)
    if start is None:
        raise ImagDispError(
            "the base input has no inline '* xyz' block, so its coordinate block "
            "cannot be replaced. Next step: give the input of the frequency run."
        )
    header = lines[start].split()
    # "* xyz <charge> <mult>" (or "*xyz ...") -- the header line itself is kept
    # verbatim; only its presence with the two numbers is required
    star = next((i for i, token in enumerate(header) if token.startswith("*")), None)
    if star is None or len(header) < star + 3:
        raise ImagDispError(
            f"the base input's coordinate header {lines[start]!r} carries no "
            "'* xyz <charge> <mult>'."
        )
    end = next(
        (i for i in range(start + 1, len(lines)) if lines[i].strip().startswith("*")),
        None,
    )
    if end is None:
        raise ImagDispError("the base input's '* xyz' block has no closing '*'.")
    block = [f"{symbol:<2} {x:>16.10f} {y:>16.10f} {z:>16.10f}" for symbol, x, y, z in atoms]
    if ensure_freq:
        simple = next(
            (i for i, line in enumerate(lines) if line.lstrip().startswith("!")), None
        )
        if simple is not None and not any(
            _FREQ_TOKEN_RE.search(line)
            for line in lines
            if line.lstrip().startswith("!")
        ):
            lines[simple] = lines[simple].rstrip() + " Freq"
    lines = lines[: start + 1] + block + lines[end:]
    text = "\n".join(lines) + "\n"
    try:
        text.encode("ascii")
    except UnicodeEncodeError as exc:  # pragma: no cover - ASCII inputs by construction
        raise ImagDispError(
            f"the generated input is not ASCII ({exc}); ORCA inputs must be pure ASCII."
        ) from exc
    return text


def render(data, *, selection: str, amplitude: float, report) -> str:
    """The menu's record: the modes, the choice, and the two written sides."""
    lines = [
        f"imaginary mode(s) in the last frequency block: "
        + ", ".join(
            f"#{mode['index']} {mode['wavenumber']:.2f} cm**-1"
            for mode in data["imaginary"]
        ),
        f"vector: {'sum of all imaginary modes' if selection == 'sum' else 'the most negative mode'} "
        "(de-weighted to Cartesian with the run's own masses, unit norm)",
        f"amplitude: {amplitude:g} Angstrom",
    ]
    for label, (path, max_move) in report.items():
        lines.append(f"  {label}: {path} (largest single-atom move {max_move:.3f} Angstrom)")
    return "\n".join(lines)


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the displacement protocol and of the ORCA-side convention."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "The imaginary-frequency cure by geometry displacement: a vector "
                "built from the imaginary modes (sum of all, or the most negative "
                "alone), normalized, applied as coords + vec * disp with "
                "successively larger displacements base_disp * (i + 1) -- defaults "
                "0.1 Angstrom up to five rounds -- rerunning opt+freq until every "
                "frequency is real; the source's success criterion is the "
                "all-real frequency check of the rerun."
            ),
            ref=(
                "Tanaka & Miyao, NIFREC (ChemRxiv preprint, 2026; MIT-licensed "
                "software github.com/tanaka-hideya/NIFREC, displacement stage read "
                "2026-09-27)"
            ),
            url="https://doi.org/10.26434/chemrxiv.15004859/v1",
            bibkey="tanaka2026nifrec",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "The ORCA-side conventions, measured on ORCA 6.1.1 (fixtures "
                "fhh_optts_freq.*): the NORMAL MODES block prints the mass-weighted "
                "vectors ('Cartesian displacements weighted by ... 1/sqrt(m[i])', "
                "the block's own header) with unit Euclidean norm -- the Cartesian "
                "pattern is v_i/sqrt(m_i) -- and the mass ORCA used appears in the "
                "CARTESIAN COORDINATES (A.U.) block (18.998 for F, 1.008 for H)."
            ),
            ref="tests/test_imag_disp.py; fixtures/orca/fhh_optts_freq.out",
        ),
    )
