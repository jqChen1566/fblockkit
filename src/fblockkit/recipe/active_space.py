"""G2: starting active-space templates for f-block systems + verification protocol.

This module only *suggests* starting spaces - choosing active orbitals requires
user insight (ORCA manual Sec. 3.13), and every suggestion is anchored to a
documented case rather than invented:

- f-only minimal space = the seven 4f/5f orbitals. Worked precedent: Dy(III)
  9-in-7, i.e. CAS(9,7), in the Chilton 19-Dy(III) workflow (chilton2025abinitio).
- tetravalent actinides: 5f and 6d must both sit in the valence space - the
  documented failure mode of 5f-in-core pseudopotentials (lu2025normconserving).
- f + d double-shell windows for f-block are documented in our records only via
  AVAS (sayfutyarova2017automated) and are marked *provisional*.

Electron counts assume regular f^n configurations
    f = Z - 57 + (3 - valence)   (lanthanides)
    f = Z - 89 + (3 - valence)   (actinides)
Documented exceptions (5f^(n-1) 6d^1 ground states) are flagged, never guessed.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..knowledge.elements import LANTHANIDES, element_z, is_f_element
from ..knowledge.models import (
    CONFIDENCE_CONFIRMED,
    CONFIDENCE_PROVISIONAL,
    EVIDENCE_LITERATURE,
    EVIDENCE_MANUAL,
    EVIDENCE_MEASURED,
    Evidence,
)


class ActiveSpaceError(ValueError):
    """Active-space suggestion input is invalid."""


_VALENCES = (2, 3, 4)

_EV_LN_PRECEDENT = Evidence(
    kind=EVIDENCE_LITERATURE,
    text=(
        "Worked precedent for the minimal f-only space: the 19 Dy(III) single-molecule-magnet "
        "workflow uses 9 electrons in 7 orbitals, CAS(9,7) - the seven 4f orbitals, matching "
        "4f^9 for Dy(III)."
    ),
    ref="Chilton N. F., Chem. Soc. Rev., 2025, 54(24), 11468-11487, DOI 10.1039/d5cs00493d",
    bibkey="chilton2025abinitio",
    url="https://doi.org/10.1039/d5cs00493d",
)

_EV_USER_INSIGHT = Evidence(
    kind=EVIDENCE_MANUAL,
    text=(
        "Often, the default !PModel guess does not provide good starting orbitals or the "
        "correct active space. Which orbitals enter as active, depends on the system and "
        "requires the users insight."
    ),
    ref="ORCA 6.1 manual Sec. 3.13",
    url="https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/CASSCF.html",
)

_EV_AN_IV_VALENCE = Evidence(
    kind=EVIDENCE_LITERATURE,
    text=(
        "For tetravalent actinides the 5f and 6d shells must both be kept in the valence space: "
        "5f-in-core pseudopotentials fail for An(IV) (AnCl4 MAD 11.8 kcal/mol vs 2.9 for the "
        "trivalent series), and the authors attribute this to the stronger 5f/6d radial overlap "
        "at higher oxidation states."
    ),
    ref="Lu J.-B., Zhang Y.-Y., Liu J.-B., Li J., J. Chem. Theory Comput., 2025, 21(1), 170-182, DOI 10.1021/acs.jctc.4c01189",
    bibkey="lu2025normconserving",
    url="https://doi.org/10.1021/acs.jctc.4c01189",
)

_EV_AVAS_DOUBLE_SHELL = Evidence(
    kind=EVIDENCE_LITERATURE,
    text=(
        "AVAS constructs active spaces from selected atomic valence orbitals; the double-shell "
        "option is the interface the paper offers for f-block problems (our records mark the "
        "f-block double-shell route as not yet verified)."
    ),
    ref="Sayfutyarova E. R., Sun Q., Chan G. K.-L., Knizia G., J. Chem. Theory Comput., 2017, 13(9), 4063-4078, DOI 10.1021/acs.jctc.7b00128",
    bibkey="sayfutyarova2017automated",
    url="https://doi.org/10.1021/acs.jctc.7b00128",
)

_EV_OCC_WINDOW = Evidence(
    kind=EVIDENCE_MANUAL,
    text=(
        "All orbitals between occupation number say 1.98 down to 0.02 will be included in the "
        "active space. (auto-ICE selection convention; used here as the verification window)"
    ),
    ref="ORCA 6.1 manual Sec. 3.14",
    url="https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/iceci.html",
)

_EV_ENTROPY = Evidence(
    kind=EVIDENCE_LITERATURE,
    text=(
        "Multireference criterion: at least one single-orbital entropy above 0.14; weak-"
        "correlation cut-off 1-2% of max(s1); read the spectrum in a localized-orbital basis."
    ),
    ref="Stein C. J., Reiher M., J. Chem. Theory Comput., 2016, 12(4), 1760-1771, DOI 10.1021/acs.jctc.6b00156; J. Comput. Chem., 2019, 40(25), 2216-2226, DOI 10.1002/jcc.25869",
    bibkey="stein2016automated",
)

_EV_SCHANGE = Evidence(
    kind=EVIDENCE_LITERATURE,
    text=(
        "S_change = SVD of (C_act^final)^T S C_act^initial: singular values near zero mean the "
        "initial guess space was missing components."
    ),
    ref="Sayfutyarova E. R., Sun Q., Chan G. K.-L., Knizia G., J. Chem. Theory Comput., 2017, 13(9), 4063-4078, DOI 10.1021/acs.jctc.7b00128",
    bibkey="sayfutyarova2017automated",
)

_EV_CASCI_IDENTITY = Evidence(
    kind=EVIDENCE_MEASURED,
    text=(
        "CASCI <= reference SCF energy is an identity; a violation means the CASSCF was handed "
        "the wrong orbital window (measured on Eu systems: wrong window gave CASCI 18.67 Ha "
        "above ROHF)."
    ),
    ref="B1_可微分DFT_cjq6 record: mcscf-default-orbital-slice-trap (2026-09)",
)

_EV_ORBITAL_CONVERGENCE = Evidence(
    kind=EVIDENCE_MEASURED,
    text=(
        "NEVPT2 depends on orbital (not just energy) convergence: an energy-flat but "
        "orbital-drifting CASSCF can shift NEVPT2 by ~8 kcal/mol."
    ),
    ref="B1_可微分DFT_cjq6 record: 经验教训汇编_20260915.md",
)


@dataclass(frozen=True)
class ActiveSpaceSuggestion:
    """One starting-space suggestion (suggestion only - the user decides)."""

    label: str
    n_electrons: int
    n_orbitals: int
    rationale: str
    confidence: str
    evidence: tuple[Evidence, ...]
    caveats: tuple[str, ...] = ()

    def casscf_line(self) -> str:
        return f"CAS({self.n_electrons}e,{self.n_orbitals}o)"


def f_electron_count(element: str, valence: int) -> int:
    """f-electron count for a regular f^n configuration (documented exceptions flagged
    by the caller, never guessed here)."""
    z = element_z(element)
    if not is_f_element(element):
        raise ActiveSpaceError(
            f"{element} is not a lanthanide or actinide. Next step: G2 covers f-block "
            "systems only; for d-block systems choose the active space manually."
        )
    if valence not in _VALENCES:
        raise ActiveSpaceError(
            f"valence {valence} is not supported. Next step: give the metal valence as "
            f"one of {_VALENCES} (the standard-filling formula assumes these)."
        )
    base = 57 if z in LANTHANIDES else 89
    count = z - base + (3 - valence)
    if count < 0:
        raise ActiveSpaceError(
            f"{element}({valence}) has a negative f count from the standard-filling formula "
            f"({count}). Next step: check the element/valence pair (e.g. Ce(IV) is 4f^0, "
            "Th(IV) is 5f^0 - these are the boundary cases)."
        )
    return count


def suggest_active_space(
    element: str, valence: int, *, f_count: int | None = None
) -> tuple[ActiveSpaceSuggestion, ...]:
    """Starting-space suggestions for one f-block ion, most conservative first."""
    z = element_z(element)
    electrons = f_electron_count(element, valence) if f_count is None else f_count
    is_actinide = z not in LANTHANIDES
    shell = "5f" if is_actinide else "4f"

    minimal = ActiveSpaceSuggestion(
        label=f"{shell}-only (minimal)",
        n_electrons=electrons,
        n_orbitals=7,
        rationale=(
            f"The seven {shell} orbitals with the {electrons} f electrons: the minimal space "
            "that captures the f manifold; confirmed by the Dy(III) CAS(9,7) precedent. "
            "Choosing the orbitals still requires user insight."
        ),
        confidence=CONFIDENCE_CONFIRMED,
        evidence=(_EV_LN_PRECEDENT, _EV_USER_INSIGHT),
        caveats=(
            "Assumes a regular f^n ground configuration; for 5f^(n-1) 6d^1 ground states "
            "(documented for some actinides) the d shell must be treated explicitly.",
        ),
    )

    double_shell_caveats = [
        f"Same electron count, extended orbital window ({shell} + d): the double-shell route "
        "for f-block is documented in our records only via AVAS and is not yet verified.",
        "Provisional: confirm with the verification protocol (occupations 0.02-1.98, "
        "entropy spectrum) before trusting a double-shell result.",
    ]
    if is_actinide and valence == 4:
        double_shell_caveats.insert(
            0,
            "Tetravalent actinides: 5f and 6d must both be in the valence space - this is the "
            "documented boundary of 5f-in-core pseudopotentials.",
        )
    double_shell = ActiveSpaceSuggestion(
        label=f"{shell} + d double shell",
        n_electrons=electrons,
        n_orbitals=12,
        rationale=(
            f"{shell} (7 orbitals) plus the d shell (5 orbitals) at the same electron count: "
            "use when the d shell participates (bonding, An(IV) valence requirements, or "
            "when occupation analysis shows d weight in the correlated orbitals)."
        ),
        confidence=CONFIDENCE_PROVISIONAL,
        evidence=(_EV_AVAS_DOUBLE_SHELL, _EV_AN_IV_VALENCE) if is_actinide and valence == 4
        else (_EV_AVAS_DOUBLE_SHELL,),
        caveats=tuple(double_shell_caveats),
    )
    return (minimal, double_shell)


def verification_protocol() -> tuple[str, ...]:
    """Checklist to run after the calculation (each line names what to check and why)."""
    return (
        "Run the diagnosis layer on the output: CASSCF must signal convergence (gradient "
        "criterion preferred) before any NEVPT2/CASPT2 on top.",
        "Check the active natural occupations: all should lie inside 0.02-1.98 "
        "(auto-ICE convention); orbitals outside the window mean the space is too wide "
        "or misplaced.",
        "Run A2 (entropy spectrum): max(s1) > 0.14 indicates multireference character; a "
        "gap in the descending spectrum suggests the correlated subset; read it only in a "
        "localized-orbital basis.",
        "Identity check: CASCI must be at or below the reference SCF energy - a violation "
        "means the wrong orbital window was handed to CASSCF.",
        "Orbital-convergence check before NEVPT2/CASPT2: an energy-flat but orbital-"
        "drifting CASSCF can shift the perturbation layer by kcal/mol.",
        "Initial-guess consistency (S_change: SVD of (C_act^final)^T S C_act^initial): "
        "singular values near zero mean the initial space missed components. NOT implemented "
        "in v0.1 - it needs both orbital coefficient matrices and the overlap matrix; "
        "export them from your engine and check by hand.",
    )


def verification_evidence() -> tuple[Evidence, ...]:
    return (
        _EV_ORBITAL_CONVERGENCE,
        _EV_OCC_WINDOW,
        _EV_ENTROPY,
        _EV_CASCI_IDENTITY,
        _EV_SCHANGE,
    )
