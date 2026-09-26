"""Relativistic Hamiltonian tiers for f-block work, and the two-component SO-CASSCF template.

Two things are encoded here, both from Wang, Wang, Wu & Cheng (2026, the
X2Ccorr two-component CASSCF paper):

1. the tier decision tree -- which relativistic Hamiltonian a calculation
   actually needs.  Its point is avoiding over-configuration: for lanthanide
   spin-orbit work the two-electron correction (X2Ccorr) moves the Nd3+ 4I
   splittings by only 1.6-4.9 cm^-1, while for light main-group molecules the
   same correction can reach 25% of the spin-orbit splitting (O2).  The tree
   also rules out one option outright: the atomic-mean-field variant (amfX2C)
   has a numerically unstable single-centre approximation in its scalar part;
2. the two-component SO-CASSCF template defaults -- the Ln3+ ground-term
   manifold state-averaged in full (Nd3+: the 4I term, 52 states), which is the
   precondition the source names for reliable spin-orbit splittings, together
   with the basis tiers, the Cholesky threshold and the quantities to report
   (Kramers doublet levels and the centre of gravity).

A self-check comes with it: reporting a one-electron-only (X2C-1e) splitting
where the two-electron term roughly doubles it (Nd3+ 4I11/2: 4349.4 vs
1850.3 cm^-1) overstates the splitting systematically, so the Hamiltonian tier
must be stated with any spin-orbit number.

Boundaries: this toolkit generates ORCA input and reads ORCA output.  ORCA
offers the one-electron variants (X2C, DKH2) -- the two-component SO-CASSCF
route itself is the source's own code, so the template here is guidance for
that route, while the tier tree applies directly to the ORCA-side choice.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..knowledge.atomic_terms import TermError, hund_term
from ..knowledge.elements import ACTINIDES, LANTHANIDES, ElementError, element_z
from ..knowledge.models import EVIDENCE_LITERATURE, EVIDENCE_MEASURED, Evidence

__all__ = [
    "RelativisticError",
    "RelativisticPlan",
    "SoCasscfTemplate",
    "plan_relativistic",
    "plan_so_casscf",
]


class RelativisticError(ValueError):
    """The advice cannot be built for the request (with a next step)."""


@dataclass(frozen=True)
class Tier:
    """One row of the decision tree."""

    need: str
    tier: str
    evidence: str


@dataclass(frozen=True)
class RelativisticPlan:
    """The recommended tier plus the tree it was chosen from."""

    need: str
    tier: str
    tree: tuple[Tier, ...]
    notes: tuple[str, ...]
    evidence: tuple[Evidence, ...]


@dataclass(frozen=True)
class SoCasscfTemplate:
    """Defaults for the two-component SO-CASSCF route on one Ln/An centre."""

    element: str
    n_f_electrons: int
    s: float
    l_value: float
    n_states: int
    active_note: str
    basis_notes: tuple[str, ...]
    numerical_notes: tuple[str, ...]
    outputs: tuple[str, ...]
    checks: tuple[str, ...]
    evidence: tuple[Evidence, ...]


_TREE = (
    Tier(
        "only scalar relativity (geometry, energies)",
        "SFX2C-1e / X2C-1e",
        "the standard scalar treatment; the two-electron term is not needed for "
        "geometries and scalar energetics",
    ),
    Tier(
        "spin-orbit coupling, light and medium main-group elements",
        "X2Ccorr",
        "the two-electron spin-same-orbit correction can carry a quarter of the "
        "splitting in light molecules (O2: ~25%)",
    ),
    Tier(
        "spin-orbit coupling, lanthanides and actinides",
        "X2CAMF (Dirac-Coulomb-Breit)",
        "already sufficient: the X2Ccorr increment on Nd3+ splittings is only "
        "1.6-4.9 cm^-1",
    ),
    Tier(
        "highest accuracy, budget permitting",
        "X2CMP + QED",
        "in the f block the QED contribution reaches tens of cm^-1 and exceeds "
        "the gauge and X2CMP differences",
    ),
    Tier(
        "core spectroscopy / core excitations",
        "X2CMP + QED",
        "there the QED and X2CMP increments are comparable in size and often of "
        "opposite sign",
    ),
)


def plan_relativistic(*, needs_soc: bool, f_block: bool, core_spectroscopy: bool = False) -> RelativisticPlan:
    """Pick the tier from the decision tree (no f-block default of X2Ccorr)."""
    if not isinstance(needs_soc, bool) or not isinstance(f_block, bool):
        raise RelativisticError(
            "needs_soc and f_block must be booleans. Next step: state whether the "
            "calculation needs spin-orbit coupling and whether the system is f-block."
        )
    if core_spectroscopy:
        need = "core spectroscopy / core excitations"
    elif not needs_soc:
        need = "only scalar relativity (geometry, energies)"
    elif f_block:
        need = "spin-orbit coupling, lanthanides and actinides"
    else:
        need = "spin-orbit coupling, light and medium main-group elements"
    tier = next(item for item in _TREE if item.need == need)
    notes = list(tier.evidence.split("; "))
    notes.append(
        "never use amfX2C: its scalar part rests on a single-centre approximation "
        "that is numerically unstable"
    )
    if needs_soc and f_block:
        notes.append(
            "this is the avoid-over-configuration case: paying for X2Ccorr in an "
            "f-block spin-orbit calculation buys only 1.6-4.9 cm^-1 (Nd3+), so the "
            "one-electron-plus-AMF tier is the default and the extra tier must be "
            "justified by the target accuracy"
        )
    return RelativisticPlan(
        need=need,
        tier=tier.tier,
        tree=_TREE,
        notes=tuple(notes),
        evidence=(
            Evidence(
                kind=EVIDENCE_LITERATURE,
                text=(
                    "The X2C Hamiltonian hierarchy for two-component CASSCF, with the "
                    "tier tree: scalar-only work uses the one-electron variants; the "
                    "two-electron spin-same-orbit correction (X2Ccorr) matters for light "
                    "main-group spin-orbit splittings (about a quarter in O2) but only "
                    "1.6-4.9 cm^-1 on Nd3+; the atomic-mean-field variant's scalar part "
                    "is numerically unstable; QED reaches tens of cm^-1 in the f block "
                    "and matters for core spectroscopy."
                ),
                ref=(
                    "Wang X., Wang S., Wu Y., Cheng L., 'Relativistic Complete Active "
                    "Space Self-consistent-Field Method with a Hierarchy of Exact "
                    "Two-Component Hamiltonians', arXiv:2602.24236v1 (publication status "
                    "unchecked in the project's reading notes)"
                ),
            ),
            Evidence(
                kind=EVIDENCE_MEASURED,
                text=(
                    "The template's state count is computed from the Hund ground term "
                    "rather than quoted: (2S+1)(2L+1) over the whole manifold gives 52 "
                    "for Nd3+ (4I), which is the source's own state-average window; the "
                    "same rule sizes any Ln3+/An3+ template."
                ),
                ref="tests/test_relativistic.py; analysis/atomic_terms.py (hund_term)",
            ),
        ),
    )


def _f_count(element: str) -> int:
    z = element_z(element)
    if z in LANTHANIDES:
        return z - 57
    if z in ACTINIDES:
        return z - 89
    raise RelativisticError(
        f"{element} is not an f-block element (a lanthanide or actinide). Next step: "
        "give the f-block element whose spin-orbit template is wanted."
    )


def plan_so_casscf(element: str = "Nd") -> SoCasscfTemplate:
    """The two-component SO-CASSCF template defaults for one Ln/An centre.

    The state-average window is the whole Hund ground-term manifold,
    ``(2S+1)(2L+1)`` states -- the source's precondition for reliable
    spin-orbit splittings.
    """
    symbol = element.strip().capitalize()
    try:
        n_f = _f_count(symbol)
    except ElementError as exc:
        raise RelativisticError(
            f"{element!r} is not an element symbol. Next step: give the f-block element "
            "by its symbol, e.g. Nd or U."
        ) from exc
    if n_f == 0:
        raise RelativisticError(
            f"{symbol}3+ has an empty f shell, so there is no ground term to average "
            "over. Next step: use a lanthanide with f electrons."
        )
    try:
        term = hund_term(3, n_f)
    except TermError as exc:  # pragma: no cover - hund_term covers 1..13
        raise RelativisticError(str(exc)) from exc
    n_states = int((2 * term.s + 1) * (2 * term.l + 1))
    return SoCasscfTemplate(
        element=symbol,
        n_f_electrons=n_f,
        s=term.s,
        l_value=term.l,
        n_states=n_states,
        active_note=(
            f"the whole 4f shell: {n_f} electrons in 7 orbitals (the two-component "
            "calculation carries 14 spinors); the state average runs over the "
            f"{n_states} states of the ground term (2S+1)(2L+1) = "
            f"({2 * int(term.s) if term.s == int(term.s) else 2 * term.s:g})"
            f"({2 * int(term.l) if term.l == int(term.l) else 2 * term.l:g})"
        ),
        basis_notes=(
            "spin-orbit contracted basis for the heavy elements (for example Dyall "
            "VTZ-SO)",
            "the light atoms (O, H, ...) use a scalar-recontracted basis (SFX2C-1e)",
        ),
        numerical_notes=(
            "Cholesky threshold 1x10^-7 Hartree for the integrals",
            "the spin-orbit treatment is two-component from the start: this is a "
            "different route from 'scalar CASSCF -> state-interaction RASSI-SO', and "
            "the two must not be mixed up when quoting numbers",
        ),
        outputs=(
            "the Kramers doublet levels",
            "each multiplet's centre of gravity (the source reports both)",
        ),
        checks=(
            "state the Hamiltonian tier with every spin-orbit number: on Nd3+ the "
            "X2C-1e splitting of the 4I11/2 multiplet is 4349.4 cm^-1 where the "
            "two-electron treatment gives 1850.3 -- reporting the one-electron value "
            "alone overstates the splitting by more than a factor of two",
            "group the levels by J multiplet (4I9/2, 4I11/2, 4I13/2, 4I15/2 for Nd3+) "
            "and compute the centre of gravity before comparing with experiment or "
            "literature",
        ),
        evidence=(
            Evidence(
                kind=EVIDENCE_LITERATURE,
                text=(
                    "Two-component SO-CASSCF defaults: the full 4f shell with the "
                    "ground-term manifold state-averaged in full (Nd3+: 52 states of "
                    "4I) as the precondition for reliable spin-orbit splittings; "
                    "spin-orbit contracted bases for the heavy elements and "
                    "scalar-recontracted ones for the light atoms; Cholesky threshold "
                    "1x10^-7 Hartree; report Kramers doublet levels and centres of "
                    "gravity."
                ),
                ref=(
                    "Wang X., Wang S., Wu Y., Cheng L., arXiv:2602.24236v1 (the "
                    "state-average window, the basis tiers and the numerical settings; "
                    "the 4349.4/1850.3 cm^-1 sensitivity example)"
                ),
            ),
        ),
    )
