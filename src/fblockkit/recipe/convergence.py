"""G4: SCF/CASSCF convergence and initial-guess setting templates (fragment generation).

Discipline (manual sentences, kept in step with group A of the decision table):

- run the default settings first; for a difficult case use ``!TRAH`` (which needs an
  auxiliary basis); raising MaxIter does not help in most cases;
- initial guess: the default ``!PModel`` often does not give the right starting
  orbitals (the manual says so itself); in the f block prefer a fragment guess, or
  AVAS; when all of that fails, read existing orbitals with ``!moread`` + ``%moinp``;
- ``!SlowConv`` works but pulls the result closer to the initial guess (the manual
  warns about this) -- getting the initial guess right matters more than convergence
  tricks;
- routes that include NEVPT2: this group measured that NEVPT2 depends on orbital
  convergence (a flat energy with drifting orbitals can cause errors of the order of
  kcal/mol), so the convergence criteria must be tightened (A5).

This module only generates "setting fragments and written guidance"; it never
modifies the user's files.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..knowledge.models import EVIDENCE_MANUAL, EVIDENCE_MEASURED, Evidence, SystemProfile

DIFFICULTY_DEFAULT = "default"
DIFFICULTY_DIFFICULT = "difficult"
_DIFFICULTIES = (DIFFICULTY_DEFAULT, DIFFICULTY_DIFFICULT)

_DEFAULT_FIRST = Evidence(
    kind=EVIDENCE_MANUAL,
    text=(
        "Please try the program with default settings before playing with the more advanced "
        "options. If you encounter convergence problems, have a look into your output, read the "
        "warning and see how the gradient and energy evolves. Try !TRAH. Increasing MaxIter will "
        "not help in many cases."
    ),
    ref="ORCA 6.1 manual §3.13 (CASSCF convergence advice)",
    url="https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/CASSCF.html",
)

_TRAH_NEEDS_AUX = Evidence(
    kind=EVIDENCE_MANUAL,
    text=(
        "To activate TRAH for your CASSCF calculation, you just need to add !TRAH in one of the "
        "simple input lines and add an auxiliary basis."
    ),
    ref="ORCA 6.1 manual §3.13",
    url="https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/CASSCF.html",
)


_EVIDENCE_MEMORY = Evidence(
    kind=EVIDENCE_MEASURED,
    text=(
        "A large-basis f-block CASSCF can need several GB per process while the generated "
        "input sets %maxcore 2000: a Yb3+ CASSCF(13,7) with SARC2-DKH-QZVP, run under TRAH, "
        "aborted after several hours of macro iterations with 'MINIMUM REQUIRED: 9345.2 MB "
        "/ MAXCORE: 2000.0 MB' (ORCA 6.1.1, measured 2026-09-26). Raise %maxcore when ORCA "
        "reports an out-of-memory abort."
    ),
    ref="measured on server 101 (fixture generated_yb3_sarc2_trah.out)",
)

_GUESS_STRATEGY = Evidence(
    kind=EVIDENCE_MANUAL,
    text=(
        "Often, the default !PModel guess does not provide good starting orbitals or the correct "
        "active space. ... in the case of lanthanides/actinides, the fragment guess, is "
        "particularly successful and leads to a rapid convergence. ... the atomic valence space "
        "selection (AVAS) offers a user-friendly guess alternative."
    ),
    ref="ORCA 6.1 manual §3.13",
    url="https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/CASSCF.html",
)

_SLOWCONV_CAUTION = Evidence(
    kind=EVIDENCE_MANUAL,
    text=(
        "While ORCA offers keywords like !SlowConv, this might not be the best option. "
        "Specifically, !SlowConv may converge to a local minimum solution that is closer to that "
        "of the initial guess."
    ),
    ref="ORCA 6.1 manual §1 quick start (Converging SCF Calculations)",
    url="https://www.faccts.de/docs/orca/6.1/tutorials/",
)

_NEVPT2_ORBITAL = Evidence(
    kind=EVIDENCE_MEASURED,
    text=(
        "NEVPT2 depends on converged orbitals, not on a converged energy: a flat energy with "
        "still-drifting orbitals can cause errors of the order of 8 kcal/mol; the gradient "
        "criterion must be tightened and the macro-iteration limit enlarged."
    ),
    ref="Group's lessons-learned compilation (B1_可微分DFT_cjq6/经验教训汇编_20260915.md, 2026-09-15)",
)


@dataclass(frozen=True)
class ConvergencePlan:
    """Convergence and initial-guess setting fragments + written guidance + provenance."""

    difficulty: str
    simple_keywords: tuple[str, ...]  # keywords appended to the ! line
    scf_block: tuple[str, ...]        # %scf block content lines
    needs_auxiliary: bool             # whether a matching auxiliary basis is needed (TRAH)
    notes: tuple[str, ...]
    evidence: tuple[Evidence, ...]


def plan_convergence(
    profile: SystemProfile | None = None,
    *,
    difficulty: str = DIFFICULTY_DEFAULT,
    pt2: bool = False,
    maxiter: int = 300,
    maxcore: int = 2000,
) -> ConvergencePlan:
    """Generate the G4 convergence/initial-guess settings.

    ``difficulty="difficult"`` turns on ``!TRAH`` (and requires a matching auxiliary
    basis); ``pt2=True`` (the route includes NEVPT2/CASPT2) attaches the orbital
    convergence reminder.
    """
    if difficulty not in _DIFFICULTIES:
        raise ValueError(
            f"difficulty={difficulty!r} is not one of {_DIFFICULTIES}. "
            f"Next step: pick default (run the defaults first) or difficult (use TRAH "
            f"for a hard case)."
        )
    simple: list[str] = []
    scf_block: list[str] = [f"MaxIter {maxiter}"]
    notes: list[str] = [
        "Run the default convergence settings first (the manual states plainly that the defaults come before hand-tuning).",
        "Initial guess: the default !PModel often does not give the right starting orbitals; "
        "in the f block prefer a fragment guess (compute the fragments first, then combine), or "
        "AVAS to select the atomic valence space; when all of that fails, read existing orbitals "
        "with !moread + %moinp.",
    ]
    evidence: list[Evidence] = [_DEFAULT_FIRST, _GUESS_STRATEGY]
    needs_aux = False
    if difficulty == DIFFICULTY_DIFFICULT:
        simple.append("TRAH")
        needs_aux = True
        notes.append("For a hard case switch to !TRAH as the manual recommends (an auxiliary basis must be given in the input).")
        evidence.append(_TRAH_NEEDS_AUX)
    else:
        notes.append("If the default settings struggle to converge, switch to !TRAH then (do not add MaxIter as the first step).")
    if pt2:
        scf_block.append("# NEVPT2/CASPT2 route: orbital convergence matters more than energy convergence (measured by this group)")
        notes.append("This route includes a perturbation layer: afterwards use the diagnosis layer to check orbital convergence and the reference weights/denominators.")
        evidence.append(_NEVPT2_ORBITAL)
    notes.append(
        f"Memory: the generated input sets %maxcore {maxcore} MB. A large-basis f-block "
        "CASSCF can need several GB per process (measured: a Yb3+ CASSCF(13,7) with "
        "SARC2-DKH-QZVP under TRAH asked for 9345 MB per process and aborted); raise "
        "%maxcore if ORCA reports out of memory."
    )
    evidence.append(_EVIDENCE_MEMORY)
    notes.append(
        "!SlowConv is available but needs caution: the manual warns that it can converge to a "
        "solution closer to the initial guess -- a correct initial guess comes first."
    )
    evidence.append(_SLOWCONV_CAUTION)
    return ConvergencePlan(
        difficulty=difficulty,
        simple_keywords=tuple(simple),
        scf_block=tuple(scf_block),
        needs_auxiliary=needs_aux,
        notes=tuple(notes),
        evidence=tuple(evidence),
    )
