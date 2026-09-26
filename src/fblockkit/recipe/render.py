"""Input rendering: Recommendation + G1 entry + G4 convergence plan -> ORCA input text.

Design trade-offs:

- the input file itself is ASCII only (comments included) -- ORCA's input parsing is
  aimed at ASCII, so prose explanations go into the run guidance returned by
  ``run_guidance`` and never into the .inp;
- the ``%casscf`` block is generated only when the caller gives an active space
  (nel/norb/mult/nroots): every generated file is then directly runnable, with no
  placeholders;
- the TRAH route needs a /C-type auxiliary basis (manual: "add an auxiliary basis"), so
  when ``auxiliary`` is not given explicitly the renderer raises an error instead of
  producing an unrunnable input.
"""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

from jinja2 import Environment, FileSystemLoader, StrictUndefined

from ..knowledge.models import Recommendation
from .basis_ecp import BasisEntry
from .convergence import ConvergencePlan

TEMPLATES_DIR = Path(__file__).resolve().parent.parent / "knowledge" / "templates"
TEMPLATE_NAME = "orca_mr.inp.j2"


class RenderError(ValueError):
    """Invalid input-rendering parameters."""


def _environment() -> Environment:
    return Environment(
        loader=FileSystemLoader(str(TEMPLATES_DIR)),
        undefined=StrictUndefined,
        trim_blocks=True,
        lstrip_blocks=True,
        keep_trailing_newline=True,
    )


def _simple_keywords(
    basis_entry: BasisEntry | None,
    method_keywords: Sequence[str],
    convergence: ConvergencePlan | None,
    auxiliary: str,
) -> list[str]:
    keywords: list[str] = []
    if basis_entry is not None:
        if basis_entry.hamiltonian:
            keywords.append(basis_entry.hamiltonian)
        if basis_entry.basis:
            keywords.append(basis_entry.basis)
        elif basis_entry.ecp:
            keywords.append(basis_entry.ecp)
    if auxiliary:
        keywords.append(auxiliary)
    keywords.extend(method_keywords)
    if convergence is not None:
        keywords.extend(convergence.simple_keywords)
    return list(dict.fromkeys(keywords))


def render_orca_input(
    recommendation: Recommendation,
    geometry: str,
    charge: int,
    mult: int,
    *,
    basis_entry: BasisEntry | None = None,
    method_keywords: Sequence[str] = (),
    casscf: Mapping[str, int] | None = None,
    convergence: ConvergencePlan | None = None,
    auxiliary: str = "",
    nprocs: int = 8,
    maxcore: int = 2000,
) -> str:
    """Render the ORCA input text.

    ``geometry``: the XYZ coordinate block (one ``element x y z`` per line, without the
    ``* xyz`` line itself). ``casscf``: ``{nel, norb, mult, nroots}``; when omitted no
    %casscf block is generated.
    """
    if not geometry.strip():
        raise RenderError(
            "the geometry is empty. Next step: give an XYZ coordinate block (one "
            "'element x y z' per line)."
        )
    if convergence is not None and convergence.needs_auxiliary and not auxiliary:
        raise RenderError(
            "this convergence plan needs a /C-type auxiliary basis (TRAH: the manual "
            "requires 'add an auxiliary basis'), but no auxiliary was given. Next step: "
            "give the /C auxiliary basis matching the orbital basis, or switch to the "
            "default convergence plan."
        )
    if casscf is not None:
        required = ("nel", "norb", "mult", "nroots")
        missing = [key for key in required if key not in casscf]
        if missing:
            raise RenderError(
                f"casscf is missing the field(s) {missing}. Next step: give all of "
                f"nel/norb/mult/nroots."
            )
    for line in geometry.splitlines():
        if line.strip() and not line.strip()[0].isalpha():
            raise RenderError(
                f"a coordinate line must start with an element symbol: {line!r}. "
                f"Next step: rewrite it as 'element x y z'."
            )
    context = {
        "input_notes": _ascii_notes(basis_entry, convergence),
        "simple_keywords": _simple_keywords(basis_entry, method_keywords, convergence, auxiliary),
        "scf_block": [line for line in (convergence.scf_block if convergence else ()) if not line.startswith("#")],
        "casscf": dict(casscf) if casscf else None,
        "charge": int(charge),
        "mult": int(mult),
        "geometry": geometry.strip("\n"),
        "nprocs": int(nprocs),
        "maxcore": int(maxcore),
    }
    text = _environment().get_template(TEMPLATE_NAME).render(**context)
    try:
        text.encode("ascii")
    except UnicodeEncodeError as exc:
        raise RenderError(
            f"the generated input contains a non-ASCII character "
            f"({exc.object[exc.start:exc.end]!r}) -- ORCA inputs must be ASCII. "
            "Next step: check whether keywords/geometry/notes contain non-ASCII text."
        ) from exc
    return text


def _ascii_notes(
    basis_entry: BasisEntry | None, convergence: ConvergencePlan | None
) -> list[str]:
    """The short comments inside the input file (ASCII): only the points directly
    relevant to running it."""
    notes: list[str] = []
    if basis_entry is not None and basis_entry.ecp:
        notes.append(f"ECP: {basis_entry.ecp}")
    if convergence is not None and convergence.needs_auxiliary:
        notes.append("TRAH requested; auxiliary (type /C) is included as given")
    notes.append("run the diagnostics afterwards (see the run guide)")
    return notes


def run_guidance(
    recommendation: Recommendation,
    convergence: ConvergencePlan | None = None,
    basis_entry: BasisEntry | None = None,
) -> str:
    """Build the run guidance (which never enters the .inp): which analysers and
    diagnostics to use once the job has finished."""
    lines = ["Run guidance:", f"- Method chain: {' / '.join(recommendation.method_chain) or '(undecided)'}"]
    if basis_entry is not None:
        lines.append(f"- Basis set / ECP: {basis_entry.line()}")
    lines.append(
        "- Once it has finished, hand the output to the diagnosis layer (diagnosis): "
        "SCF/CASSCF convergence, active-space occupations (outside 0.02-1.98 points to a "
        "space that is too large), CASPT2 reference weights and smallest denominator, and "
        "the sign of the NEVPT2 class contributions; with SOC, verify the state order and "
        "identity."
    )
    if convergence is not None:
        for note in convergence.notes:
            lines.append(f"- {note}")
    if recommendation.refusals:
        lines.append("- This route includes refusals:")
        for refusal in recommendation.refusals:
            lines.append(f"  - {refusal}")
    if recommendation.warnings:
        lines.append("- Reminders:")
        for warning in recommendation.warnings:
            lines.append(f"  - {warning}")
    return "\n".join(lines)
