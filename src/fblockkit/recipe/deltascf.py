"""4.3: DeltaSCF / MOM excited-state SCF inputs (menu 27).

ORCA's DeltaSCF route converges the SCF to a chosen excited-state solution (a
higher-energy stationary point of the SCF energy surface) by constraining the
frontier occupations and following them with a maximum-overlap metric.  The
measured grammar (ORCA 6.1.1 manual section 5.34, its own examples):

- the simple line carries the base method keywords plus ``DeltaSCF``;
  ``UHF`` appears in every manual example (open-shell excited solutions);
- ``%scf ALPHACONF <list>`` sets the alpha frontier occupations -- "any
  orbital below the first zero is assumed occupied, any above the last
  occupied empty": ``0,1`` = HOMO->LUMO, ``0,0,1`` = HOMO->LUMO+1,
  ``0,1,1`` = HOMO-1->LUMO.  ``BETACONF`` mirrors it for beta;
- ``IONIZEALPHA <n>`` / ``IONIZEBETA <n>`` remove an electron from MO ``n``
  (the core-ionization route; charge and multiplicity stay those of the
  *reference* system);
- ``PMOM true`` swaps the overlap metric for the Hratchian group's
  projection-operator variant; ``KeepInitialRef true`` is the literature's
  IMOM; ``DoMOM`` defaults to true;
- hard cases use the tactic keywords ``FreezeAndRelease`` (implicit
  SOSCF/NoDIIS/NoTRAH) or ``GMF`` (needs a target saddle order, e.g.
  ``SOSCFSPO``), and the manual recommends starting from a converged
  ground-state's orbitals (``!MORead`` + ``%moinp``).

The generated file is ASCII and carries the source's caveats in the run
guidance: single-determinant states only; open-shell singly excited states
break spin symmetry and need spin purification; not apt for most pi-pi*
states (the manual names benzene's HOMO->LUMO explicitly); orbitals should
come from a converged ground state.
"""

from __future__ import annotations

import re

__all__ = [
    "DeltaScfError",
    "parse_conf_list",
    "deltascf_input",
    "run_guidance_lines",
]

_CONF_RE = re.compile(r"^\d+(,\d+)*$")


class DeltaScfError(ValueError):
    """Invalid DeltaSCF parameters; the message carries the next step."""


def parse_conf_list(text: str) -> str:
    """Normalise an occupation list like ``0,1`` (validation, spacing)."""
    tokens = [token for token in text.replace(" ", "").split(",") if token != ""]
    if not tokens:
        raise DeltaScfError(
            "the occupation list is empty. Next step: give it like the manual's "
            "examples -- '0,1' (HOMO->LUMO), '0,0,1' (HOMO->LUMO+1), '0,1,1' "
            "(HOMO-1->LUMO), '0,2' (double HOMO->LUMO, only for RHF with "
            "ALPHACONF)."
        )
    joined = ",".join(tokens)
    if _CONF_RE.match(joined) is None:
        raise DeltaScfError(
            f"the occupation list {text!r} is not a comma-separated list of small "
            "integers."
        )
    return joined


def deltascf_input(
    coordinates: tuple[tuple[str, float, float, float], ...],
    *,
    charge: int,
    multiplicity: int,
    alpha_conf: str | None = None,
    beta_conf: str | None = None,
    ionize_alpha: int | None = None,
    keywords: str = "PBE0 def2-TZVP UHF",
    mom: str = "mom",
    tactics: str = "none",
    gs_gbw: str | None = None,
    maxcore: int = 2000,
) -> str:
    """Render the DeltaSCF input (one excitation spec: conf lists or ionize)."""
    if (alpha_conf is None) == (ionize_alpha is None):
        raise DeltaScfError(
            "give exactly one occupation spec: ALPHACONF (a list like '0,1') or an "
            "IONIZEALPHA orbital index."
        )
    if beta_conf is not None and alpha_conf is None:
        raise DeltaScfError("BETACONF needs an ALPHACONF list beside it.")
    if mom not in ("mom", "pmom", "imom"):
        raise DeltaScfError(f"unknown MOM metric {mom!r}; use mom, pmom or imom.")
    if tactics not in ("none", "freeze", "gmf"):
        raise DeltaScfError(f"unknown tactic {tactics!r}; use none, freeze or gmf.")
    if multiplicity < 1:
        raise DeltaScfError(f"multiplicity {multiplicity} is not positive.")
    if ionize_alpha is not None and ionize_alpha < 0:
        raise DeltaScfError("the ionization orbital index cannot be negative.")

    simple = f"! {keywords} DeltaSCF"
    if tactics == "freeze":
        simple += " FreezeAndRelease"
    elif tactics == "gmf":
        simple += " GMF"
    lines = [simple, f"%maxcore {maxcore}"]
    if gs_gbw:
        stem = gs_gbw[:-4] if gs_gbw.endswith(".gbw") else gs_gbw
        lines += ["! MORead", f'%moinp "{stem}.gbw"']
    lines.append("%scf")
    if alpha_conf is not None:
        lines.append(f"  ALPHACONF {parse_conf_list(alpha_conf)}")
        if beta_conf is not None:
            lines.append(f"  BETACONF {parse_conf_list(beta_conf)}")
    else:
        lines.append(f"  IONIZEALPHA {ionize_alpha}")
    if mom != "mom":
        lines.append("  DoMOM true")
        lines.append("  PMOM true" if mom == "pmom" else "  KeepInitialRef true")
    if tactics == "freeze":
        lines.append("  SOSCFMaxStep 0.1")
    lines.append("end")
    lines.append(f"* xyz {charge} {multiplicity}")
    for element, x, y, z in coordinates:
        lines.append(f"{element:<2} {x:>16.10f} {y:>16.10f} {z:>16.10f}")
    lines.append("*")
    text = "\n".join(lines) + "\n"
    try:
        text.encode("ascii")
    except UnicodeEncodeError as exc:  # pragma: no cover - callers keep ASCII input
        raise DeltaScfError(
            f"the generated input is not ASCII ({exc}); ORCA inputs must be pure ASCII."
        ) from exc
    return text


def run_guidance_lines() -> tuple[str, ...]:
    """The source's caveats, printed with every generated DeltaSCF input."""
    return (
        "run it with ORCA and check the converged state (occupations and character): "
        "DeltaSCF targets a higher SCF stationary point, and a run that fell back to "
        "the ground state shows no warning",
        "single-determinant states only: open-shell singly excited states inherently "
        "break spin symmetry and need spin purification, and multi-determinant cases "
        "are out of scope (the manual's own boundary)",
        "not apt for most pi->pi* states -- the manual names benzene's HOMO->LUMO "
        "explicitly -- but reasonable for particle-hole states, spatially separated "
        "occupied/virtual pairs, and closed-shell doubly excited states",
        "start from a converged ground-state calculation's orbitals (the MORead + "
        "%moinp lines) whenever one exists; for a core ionization, localize the core "
        "orbital first if it is not atomic",
        "the wavefunction is an SCF solution: do not feed it to single-reference "
        "correlation methods as if it were a ground state",
    )
