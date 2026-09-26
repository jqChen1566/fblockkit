"""The DMET embedding recipe for f-block single-ion magnets (read from Ai et al. 2025).

Density matrix embedding theory in the form this group's literature uses: a
low-level determinant (ROHF) is partitioned into an impurity (the lanthanide's
own orbitals) and an environment (everything else) by Loewdin-orthogonalised
localised orbitals, bath orbitals come from the Schmidt decomposition of the
environment block, the static correlation is solved inside the cluster space,
and the magnetic parameters are read from the resulting multiplets.  The source
validates this against full all-electron CASSCF-SO on three real 4f SIMs
(MAE 0.6-7.8 cm^-1 at the CASSCF-SO level; with SC-NEVPT2 the error grows to
8.8-62.7 cm^-1 on the metal-only cluster and comes back to 13.2 cm^-1 when the
cluster is expanded to the metal plus its nearest neighbours -- so that
expansion is the recipe's default whenever NEVPT2 is used).

The recipe exists here because two of its pieces are already implemented in the
toolkit and they gate each other: the Delta S_E criterion (menu 12's
environment-spin section is the same quantity) decides whether the low-level
wave function is physically the right one -- in the source's 1Dy case the
default DIIS solution puts four of five singly occupied orbitals on the
*ligands*, Delta S_E = 2.766 instead of 0.007, and the resulting crystal-field
parameters carry a 42x worse MAE (357.7 against 8.6 cm^-1).  The other piece is
the accounting of the core contribution, which the project settled in a
separate note: the two expressions that appear across the sources are the same
quantity written with two spin accountings (spatial orbitals with an explicit
spin sum, or spin orbitals counted individually), so an implementation must pair
each expression with its own convention -- or, simplest, use the explicit Fock
form that is convention-free.

Scope: geometry-free guidance.  This module only describes the workflow, its
gates and its citations; the external code (liblan) has not been licence-checked
and is not called or bundled.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..knowledge.elements import ACTINIDES, LANTHANIDES, ElementError, element_z
from ..knowledge.models import EVIDENCE_LITERATURE, EVIDENCE_MEASURED, Evidence

__all__ = ["DmetError", "DmetPlan", "plan_dmet"]


class DmetError(ValueError):
    """The recipe cannot be built for the requested element (with a next step)."""


@dataclass(frozen=True)
class DmetStep:
    """One workflow step: what to set, and why it is set that way."""

    number: int
    title: str
    setting: str
    note: str


@dataclass(frozen=True)
class DmetPlan:
    """The DMET workflow: prerequisites, steps, gates, expectations and boundaries."""

    element: str
    n_f_electrons: int
    steps: tuple[DmetStep, ...]
    basis_notes: tuple[str, ...]
    gates: tuple[str, ...]
    expectations: tuple[str, ...]
    boundaries: tuple[str, ...]
    evidence: tuple[Evidence, ...]


def _f_count(element: str) -> int:
    z = element_z(element)
    if z in LANTHANIDES:
        return z - 57
    if z in ACTINIDES:
        return z - 89
    raise DmetError(
        f"{element} is not an f-block element (a lanthanide or actinide). Next step: "
        "give the f-block element whose impurity the embedding should use."
    )


def plan_dmet(element: str = "Dy") -> DmetPlan:
    """The DMET recipe for one Ln/An centre (the source's own settings plus its gates)."""
    symbol = element.strip().capitalize()
    try:
        n_f = _f_count(symbol)
    except ElementError as exc:
        raise DmetError(
            f"{element!r} is not an element symbol. Next step: give the f-block element "
            "by its symbol, e.g. Dy or U."
        ) from exc
    if n_f == 0:
        raise DmetError(
            f"{symbol}3+ has an empty f shell, so there is no f impurity to embed. Next "
            "step: use a lanthanide with f electrons."
        )
    steps = (
        DmetStep(
            0,
            "Scalar relativity",
            "SFX2C-1e (the one-electron spin-free X2C variant)",
            "The source uses PySCF's built-in one-electron variant; the project's basis "
            "and ECP notes apply unchanged (ORCA 6.1 has no DKH3).",
        ),
        DmetStep(
            1,
            "Low-level SCF",
            "ROHF from the default minimal-atomic-orbital guess",
            "ROHF over UHF, so the orbitals stay spin-pure for the later spin-orbit "
            "step; the guess may be wrong -- that is what the gates below catch.",
        ),
        DmetStep(
            2,
            "Loose convergence (deliberate)",
            "Default DIIS to delta E < 1e-4 AND |g| < 0.01",
            "Both bounds: the project's standing rule is that an energy alone does not "
            "show convergence, and the source uses the same pair here.",
        ),
        DmetStep(
            3,
            "Embedding construction",
            "Loewdin-orthogonalised localised orbitals; cluster A = the lanthanide's "
            "orbitals, environment B = the rest; bath threshold 1e-13",
            "One-shot DMET (no correlation-potential self-consistency). The Loewdin "
            "orthogonality between cluster and core orbitals is a precondition of the "
            "mean-field potential, not a detail; the impurity integrals scale as |A|^4 "
            "rather than N_AO^4 -- the source of the cost advantage.",
        ),
        DmetStep(
            4,
            "Subspace R-DIIS",
            "Apply the regularised residual R = Delta S_E inside the cluster subspace",
            "R-DIIS forbids convergence to a solution with R != 0 instead of trying to "
            "start nearer the right one; inside the subspace each iteration is 12-64x "
            "cheaper than in the full space (and the subspace integrals fit in memory).",
        ),
        DmetStep(
            5,
            "Full-space finish",
            "Return to the full space for a few tight R-DIIS rounds",
            "By then Delta S_E is near zero and the iteration behaves like plain DIIS.",
        ),
        DmetStep(
            6,
            "Static correlation",
            f"SA-CASSCF inside the cluster space: CAS({n_f}e,7o) for a {symbol}3+ centre",
            "For f^n the seven-orbital window is the f shell itself; the source averages "
            "over the whole ground-term manifold (Dy: 21 sextets, Er: 35 quartets -- "
            "their settings; the multiplet decomposition itself is marked unchecked in "
            "the project's reading notes).",
        ),
        DmetStep(
            7,
            "Spin-orbit coupling",
            "Diagonalise the cluster Hamiltonian (scalar + SOMF spin-orbit) in the CSF "
            "basis",
            "The state-interaction route, not a variational SOC; the CSFs come from "
            "step 6.",
        ),
        DmetStep(
            8,
            "Dynamic correlation (optional)",
            "Strongly contracted NEVPT2 inside the cluster space; when NEVPT2 is used, "
            "expand the cluster of step 3 to the metal plus its nearest coordinating "
            "atoms (the CAS window of step 6 is unchanged)",
            "Optional; with it the embedding error grows systematically (dynamic "
            "correlation is less local than the static part), and the metal-only "
            "cluster that serves CASSCF-SO is not enough -- expanding 3Dy's cluster "
            "to Dy + nearest C brings its error back by 4.75x (62.7 -> 13.2 cm^-1).",
        ),
        DmetStep(
            9,
            "Post-processing",
            "Extract the model spin-Hamiltonian parameters from the SI multiplet "
            "energies and wave functions",
            "The source reports levels only; relaxation and fit-quality diagnostics are "
            "the crystal-field workflow's business (menus 10/11 and A4/A5 here).",
        ),
    )
    basis_notes = (
        "ANO-RCC-VTZP for the magnetic centre and its nearest coordinating atoms (the "
        "source: Dy with O; Er with B; Dy with the nearest C).",
        "ANO-RCC-VDZP for every other atom.",
        "Cholesky decomposition for the integral storage.",
    )
    gates = (
        "Delta S_E is the selection gate, and this toolkit computes the same quantity: "
        "the environment-spin section of menu 12 (from an exact active-space state) and "
        "menu 1's local-spin table (from an SCF output) both answer 'is the spin "
        "polarisation on the metal?'. The source's 1Dy numbers: 0.007 for the correct "
        "ROHF solution against 2.766 for the default one.",
        "Using a wrong low-level solution costs 42x in the fitted crystal-field MAE "
        "(357.7 against 8.6 cm^-1) -- so the gate runs before any parameter is fitted, "
        "not as a sanity check afterwards.",
        "Convergence needs both the energy and the orbital gradient; the energy alone "
        "can look flat on a wrong solution.",
        "The correct solution was lower in energy in the source's case, so a "
        "multi-start 'take the lowest' works there -- but that is this system's "
        "conclusion, not a rule (the project has counterexamples: the lowest stationary "
        "point is not always the target one).",
    )
    expectations = (
        "CASSCF-SO inside the cluster, embedded against the source's all-electron "
        "values (its Table 2): MAE 7.8 cm^-1 (1Dy), 0.6 (2Er), 6.7 (3Dy); relative "
        "error within 2.3% -- for SIM modelling the source calls this negligible.",
        "SC-NEVPT2, metal-only cluster (its Table 3): 10.5 / 8.8 / 62.7 cm^-1 -- the "
        "embedding error grows systematically, because dynamic correlation is less "
        "local than the static part the CASSCF handles.",
        "SC-NEVPT2 with the cluster expanded to the metal plus its nearest "
        "coordinating atoms: 3Dy falls 4.75x (62.7 -> 13.2 cm^-1, relative error "
        "1.2%). That expansion is this recipe's default whenever NEVPT2 is used.",
    )
    boundaries = (
        "Core-contribution accounting (settled in the project's reading notes): the two "
        "expressions that appear across the DMET sources (sum over spatial orbitals with "
        "an explicit spin sum versus a per-spin-orbital sum) are the same quantity under "
        "two accountings; pair each expression with its own convention, or use the "
        "explicit Fock form, which is convention-free.",
        "One-shot DMET: the source does not make the correlation potential "
        "self-consistent, and neither does this recipe.",
        "The external code (liblan) has not been licence-checked and is neither called "
        "nor bundled; this recipe is method guidance only.",
    )
    evidence = (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "The DMET + SA-CASSCF-SO workflow for 4f single-ion magnets, with R-DIIS "
                "and sR-DIIS to steer the low-level SCF onto the physically correct "
                "solution: ROHF at the loose pair of criteria, Loewdin-localised "
                "cluster/environment split, subspace R-DIIS with R = Delta S_E, "
                "SA-CASSCF inside the cluster, SOMF state-interaction spin-orbit, "
                "optional SC-NEVPT2; validated against all-electron CASSCF-SO on three "
                "real Dy/Er SIMs (MAE 0.6-7.8 cm^-1 at the CASSCF-SO level; with "
                "SC-NEVPT2 8.8-62.7 cm^-1 on the metal-only cluster, back to 13.2 "
                "cm^-1 when the cluster includes the nearest neighbours)."
            ),
            ref=(
                "Ai Y., Li Z.-W., Guan Z.-B., Jiang H., J. Chem. Theory Comput., 2025, "
                "21(19), 9631-9640, DOI 10.1021/acs.jctc.5c01336"
            ),
            bibkey="ai2025density",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "The gate is implemented here, not just quoted: menu 12's "
                "environment-spin section computes Delta S_E from the exact active-space "
                "state (zero exactly when the environment is spin-balanced; ln 2 for one "
                "unpaired electron confined to the environment), and menu 1's local-spin "
                "table gives the SCF-side fragment spin distribution that the source's "
                "criterion is applied to."
            ),
            ref="tests/test_environment_spin.py; analysis/environment_spin.py",
        ),
    )
    return DmetPlan(
        element=symbol,
        n_f_electrons=n_f,
        steps=steps,
        basis_notes=basis_notes,
        gates=gates,
        expectations=expectations,
        boundaries=boundaries,
        evidence=evidence,
    )


def render(plan: DmetPlan) -> str:
    """Render the plan as the text the tool guide prints."""
    lines = [
        f"DMET embedding recipe for {plan.element} (f^{plan.n_f_electrons} centre; from "
        "Ai et al., JCTC 2025, the fields above) -- method guidance, the external code "
        "is not bundled:",
        "",
    ]
    for step in plan.steps:
        lines.append(f"  {step.number}. {step.title}: {step.setting}")
        lines.append(f"     why: {step.note}")
    lines += ["", "Basis rules:"]
    lines += [f"  - {note}" for note in plan.basis_notes]
    lines += ["", "Gates (run these before trusting any parameter):"]
    lines += [f"  - {item}" for item in plan.gates]
    lines += ["", "Accuracy expectations (the source's own tables):"]
    lines += [f"  - {item}" for item in plan.expectations]
    lines += ["", "Boundaries and settled conventions:"]
    lines += [f"  - {item}" for item in plan.boundaries]
    return "\n".join(lines)
