"""SCF rescue triage and corrected-input proposals (feature D2).

What this module is for
-----------------------
Two entry points for a job whose SCF cannot be trusted:

- :func:`triage` reads one parsed ORCA output and reports what is wrong with its SCF:
  no convergence, a crash inside the SCF module, a suspiciously short or long
  convergence, a DIIS error that rises again after AO-DIIS was switched on, an energy
  trajectory that oscillates without collapsing, or a run declared converged while the
  criterion its convergence check enforces is still above the printed tolerance.
- :func:`propose_fixes` turns those findings into *new* input files: a ``SlowConv``
  variant, a two-step route that runs a cheap pre-SCF and reads its orbitals back,
  and (for a rising or oscillating trajectory) a ``TRAH`` variant -- the manual's
  robust second-order SCF (the capability map: ORCA has no keyword spelled
  ARH; TRAH is its trust-region augmented-Hessian route, AutoTRAH default on, and
  SOSCF the approximate second-order one).

Discipline
----------
- Proposals are new strings for new files; the caller's input file is never touched.
- Editing is textual and conservative: every edit is listed in ``FixProposal.changes``
  and nothing else is reordered or dropped (a ``%scf`` block that is not carried over
  is kept in the file as comments).
- No proposal claims the fix will work, and no proposal may consist of a MaxIter
  change alone -- the manual states "Increasing MaxIter will not help in many cases."
- Every finding carries its provenance: manual quotes (section + URL) or this group's
  measured records.

Why these checks are Python and not rows of ``knowledge/rules/*.yaml``
--------------------------------------------------------------------
A rule row can only compare one fact against a threshold under a fixed title. These
checks must report the *measured* numbers (cycle counts, DIIS errors, the
achieved-vs-tolerance table) inside the message, exactly like
:mod:`fblockkit.diagnosis.cross_level`. Keeping them out of the rule table also keeps
``diagnose()`` unchanged: a healthy 64-cycle run still yields no finding there, while
the triage reports it (see ``tests/test_scf_rescue.py``).

Convergence-check scope (why the other printed rows are not triaged)
--------------------------------------------------------------------
ORCA's ``ConvCheckMode`` decides which printed criteria are pass/fail: mode 2 (the default,
and what every standard preset sets) checks the energy change only, mode 0 checks all of
them, and mode 1 stops as soon as one is met. The remaining rows of the ``SCF CONVERGENCE``
block are informational, and on a healthy run they sit above their tolerances routinely
(measured on ``n2_hf_clean.out``: 8 cycles, terminated normally, three rows above
tolerance). The triage therefore tests only the rows the mode enforces; the mode is read
from the echoed input setting, else from the printed label, else assumed to be the
documented default -- and the finding always names which of the three it used, so an
unrecognized future label degrades to the default instead of to a warning.

Reader location (v0.2)
----------------------
The DIIS/SOSCF iteration tables and the ``SCF CONVERGENCE`` summary are read by the
parser layer (``parsers/orca.py``, the ``scf`` section: ``diis_rows``, ``criteria``,
``check_mode`` and the rest) and consumed here from the parse result, so ``triage()`` is
a pure in-memory function and never re-reads the file. Both readers are header-driven
and were developed against the fixtures in ``fixtures/orca/``.
"""

from __future__ import annotations

import textwrap
from dataclasses import dataclass
from typing import Iterable, Sequence

from ..knowledge.models import (
    EVIDENCE_MANUAL,
    EVIDENCE_MEASURED,
    Evidence,
    Finding,
    ParseResult,
)
# The SCF tables and the ConvCheckMode semantics are read by the parser layer; the mode
# constants are re-exported here because the findings below reason about them.
from ..parsers.orca import (
    CONVCHECK_ALL,
    CONVCHECK_ENERGY,
    CONVCHECK_ONE_IS_ENOUGH,
    DEFAULT_CONVCHECK_MODE,
)

# --- public constants -------------------------------------------------------

PROGRAM_ORCA = "orca"

#: Convergence after at most this many cycles counts as a pseudo-convergence signal
#: (group's measured record, see ``EV_PSEUDO_CONVERGENCE_RECORD``).
PSEUDO_CYCLE_LIMIT = 3

#: Triage threshold (ours, not a manual value): from this cycle count on, the run was
#: expensive enough that the manual's pre-SCF route is worth proposing.
LONG_CYCLE_LIMIT = 50

#: Triage threshold (ours): the DIIS error must rise to at least this factor of its
#: value when AO-DIIS was switched on before it counts as "rising again".
DIIS_REBOUND_FACTOR = 2.0

#: Triage margin (ours): a printed criterion counts as unmet only when it exceeds the
#: tolerance ORCA printed next to it by at least this factor, so that last-digit
#: excursions do not produce noise.
CRITERIA_MARGIN = 5.0

#: Oscillation triage (ours, calibrated on the fixture set): the energy tail analyzed
#: for sign alternation (in DIFFERENCES), the minimum number of direction changes,
#: and the persistence floor -- the last step must still carry at least this fraction
#: of the window's largest step, so an alternating-but-collapsing tail (the normal
#: endgame of a converging run) is not an oscillation.
OSCILLATION_WINDOW = 12
OSCILLATION_MIN_FLIPS = 3
OSCILLATION_PERSISTENCE = 0.1
#: Steps below this are numerical noise, not oscillation (Eh).
OSCILLATION_FLOOR = 1e-6

# ConvCheckMode values and the printed-label mapping live in the parser layer
# (parsers/orca.py) and are imported above; the findings here read them.

RULE_NOT_CONVERGED = "SCF-NOT-CONVERGED"
RULE_ABORTED_NO_VERDICT = "SCF-ABORTED-NO-VERDICT"
RULE_PSEUDO_CONVERGENCE = "SCF-PSEUDO-CONVERGENCE"
RULE_LONG_CONVERGENCE = "SCF-LONG-CONVERGENCE"
RULE_DIIS_REBOUND = "SCF-DIIS-REBOUND"
RULE_ENERGY_OSCILLATION = "SCF-ENERGY-OSCILLATION"
RULE_CRITERIA_UNMET = "SCF-CONVERGED-CRITERIA-UNMET"

FIX_SLOWCONV = "slowconv"
FIX_PRESCF = "prescf"
FIX_TRAH = "trah"

# --- manual quotes (verbatim; section numbers refer to the ORCA 6.1 manual) --
# One quote per Evidence entry: quotes from different sentences are never glued
# together, so every ``text`` below is a contiguous piece of the manual.
#
# Rendering notes (the local manual is a Markdown conversion of the PDF):
# - ``_Q_LOOK_AT_ORBITALS`` writes ``Print[P_GuessOrb]`` where the Markdown escapes the
#   underscore (``P\_GuessOrb``); the quote is the rendered text.
# - ``_Q_NOTRAH`` joins two code blocks ("! NOTRAH" / "%scf AutoTRAH false end") that the
#   PDF prints on separate lines; the quotation marks only delimit the code snippets.

_Q_PRESCF_ROUTE = (
    "Perform a small basis set (SV) calculation in using the LSD or BP functional and RI "
    "approximation with a cheap auxiliary basis set. Set Convergence=Loose and MaxIter=200 "
    "or so. The key point is to use a large damping factor and damp until the DIIS comes "
    "into a domain of convergence. This is accomplished by SlowConv or even VerySlowConv."
)
_Q_BETTER_GUESS = (
    "Despite all efforts you may still find molecules where SCF convergence is poor. These "
    "are almost invariably related to open-shell situations and the answer is almost always "
    "to provide “better” starting orbitals."
)
_Q_LOOK_AT_ORBITALS = (
    "Carefully look at the starting orbitals (Print[P_GuessOrb]=1) and see if they make "
    "sense for your molecule."
)
_Q_SLOWCONV = (
    "While ORCA offers keywords like !SlowConv, this might not be the best option. "
    "Specifically, !SlowConv may converge to a local minimum solution that is closer to "
    "that of the initial guess."
)
_Q_DIIS_STUCK = (
    "If the DIIS gets stuck at some error 0.001 or so the SOSCF (or even better TRAH) could "
    "be put in operation from this point on."
)
_Q_TRAH = (
    "Note that for troublesome or lacking SCF convergence the TRAH algorithm should be used "
    "(see Sec. Trust-Region Augmented Hessian (TRAH) SCF). If not turned off explicitly, "
    "TRAH is switched on automatically whenever convergence problems are present by means "
    "of the AutoTRAH feature (see Sec. Trust-Region Augmented Hessian (TRAH) SCF)."
)
_Q_NOTRAH = 'To disable automatic activation: "! NOTRAH" or "%scf AutoTRAH false end".'
_Q_LEVELSHIFT = (
    "Use large level shifts. This increases the number of iterations but stabilizes the "
    "converger. (shift shift 0.5 erroff 0 end)"
)
_Q_CONVCHECK_MODES = (
    "ConvCheckMode  2   # = 0: check all convergence criteria\n"
    "# = 1: stop if one of criterion is met, this is sloppy!\n"
    "# = 2: check change in total energy and in one-electron energy\n"
    "#       Converged if delta(Etot)<TolE and delta(E1)<1e3*TolE"
)
_Q_CONVCHECK_SEMANTICS = (
    "If ConvCheckMode=0, all convergence criteria have to be satisfied for the program to "
    "accept the calculation as converged, which is a quite rigorous criterion. In this mode, "
    "the program also has mechanisms to decide that a calculation is converged even if one "
    "convergence criterion is not fulfilled but the others are overachieved. ConvCheckMode=1 "
    "means that one criterion is enough. This is quite dangerous, so ensure that none of the "
    "criteria are too weak, otherwise the result will be unreliable. The default "
    "ConvCheckMode=2 is a check of medium rigor — the program checks for the change in total "
    "energy and for the change in the one-electron energy."
)
_Q_DENSITY_TOL = (
    "If you have small eigenvalues of the overlap matrix, the density may not be converged "
    "to the number of significant figures requested by TolMaxP and TolRMSP."
)
_Q_CONVFORCED = (
    "Irrespective of the ConvForced value that has been chosen, properties or numerical "
    "calculations (NumGrad, NumFreq) will not be performed on non-converged wavefunctions!"
)
_Q_MAXITER = (
    "Please try the program with default settings before playing with the more advanced "
    "options. If you encounter convergence problems, have a look into your output, read the "
    "warning and see how the gradient and energy evolves. Try !TRAH. Increasing MaxIter will "
    "not help in many cases."
)
_Q_CMATRIX = (
    "Use the orbitals of this calculation and GuessMode=CMatrix to start a calculation with "
    "the target basis set."
)
_Q_CMATRIX_ANION = (
    "This is always required when the orbital energies of the small basis set calculation "
    "are positive, as will be the case for anions."
)

_MANUAL_URL_ROOT = "https://www.faccts.de/docs/orca/6.1/manual/"
_MANUAL_URL_TUTORIALS = "https://www.faccts.de/docs/orca/6.1/tutorials/"
_MANUAL_URL_CASSCF = (
    "https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/CASSCF.html"
)

EV_PRESCF_ROUTE = Evidence(
    kind=EVIDENCE_MANUAL,
    text=_Q_PRESCF_ROUTE,
    ref="ORCA 6.1 manual §2.6.9 (Tips and Tricks: Converging SCF Calculations)",
    url=_MANUAL_URL_ROOT,
)
EV_BETTER_GUESS = Evidence(
    kind=EVIDENCE_MANUAL,
    text=_Q_BETTER_GUESS,
    ref="ORCA 6.1 manual §2.6.9",
    url=_MANUAL_URL_ROOT,
)
EV_LOOK_AT_ORBITALS = Evidence(
    kind=EVIDENCE_MANUAL,
    text=_Q_LOOK_AT_ORBITALS,
    ref="ORCA 6.1 manual §2.6.9",
    url=_MANUAL_URL_ROOT,
)
EV_SLOWCONV_CAUTION = Evidence(
    kind=EVIDENCE_MANUAL,
    text=_Q_SLOWCONV,
    ref="ORCA 6.1 manual §1.7.13 (SCF Convergence Problems)",
    url=_MANUAL_URL_TUTORIALS,
)
EV_DIIS_STUCK = Evidence(
    kind=EVIDENCE_MANUAL,
    text=_Q_DIIS_STUCK,
    ref="ORCA 6.1 manual §2.6.9",
    url=_MANUAL_URL_ROOT,
)
EV_TRAH = Evidence(
    kind=EVIDENCE_MANUAL,
    text=_Q_TRAH,
    ref="ORCA 6.1 manual §2.6.4 (Direct Inversion in Iterative Subspace)",
    url=_MANUAL_URL_ROOT,
)
EV_NOTRAH = Evidence(
    kind=EVIDENCE_MANUAL,
    text=_Q_NOTRAH,
    ref="ORCA 6.1 manual §2.6.7 (Trust-Region Augmented Hessian (TRAH) SCF)",
    url=_MANUAL_URL_ROOT,
)
_Q_TRAH_WHEN = (
    "for troublesome or lacking SCF convergence the TRAH algorithm should be used ... "
    "If not turned off explicitly, TRAH is switched on automatically whenever convergence "
    "problems are present by means of the AutoTRAH feature"
)
_Q_TRAH_SOSCF = (
    "On the other hand, SOSCF is useful when DIIS gets stuck at some error around ~0.001 "
    "or 0.0001. Such cases were the primary motive for the implementation of SOSCF into "
    "ORCA."
)
EV_TRAH_WHEN = Evidence(
    kind=EVIDENCE_MANUAL,
    text=_Q_TRAH_WHEN,
    ref="ORCA 6.1 manual §2.6.4 (Direct Inversion in Iterative Subspace)",
    url=_MANUAL_URL_ROOT,
)
EV_TRAH_WHEN_SOSCF = Evidence(
    kind=EVIDENCE_MANUAL,
    text=(
        _Q_TRAH_SOSCF
        + " TRAH (the trust-region augmented-Hessian second-order SCF, the section after "
        "SOSCF) is the manual's robust answer for difficult cases; the capability "
        "map: ORCA has no keyword spelled ARH -- TRAH is its second-order augmented-Hessian "
        "route (AutoTRAH default on), and SOSCF the approximate one."
    ),
    ref="ORCA 6.1 manual §2.6.6 (Approximate Second Order SCF)",
    url=_MANUAL_URL_ROOT,
)
EV_LEVELSHIFT = Evidence(
    kind=EVIDENCE_MANUAL,
    text=_Q_LEVELSHIFT,
    ref="ORCA 6.1 manual §2.6.9",
    url=_MANUAL_URL_ROOT,
)
EV_CONVFORCED = Evidence(
    kind=EVIDENCE_MANUAL,
    text=_Q_CONVFORCED,
    ref="ORCA 6.1 manual §2.6.1 (Convergence Tolerances)",
    url=_MANUAL_URL_ROOT,
)
EV_CONVCHECK_MODES = Evidence(
    kind=EVIDENCE_MANUAL,
    text=_Q_CONVCHECK_MODES,
    ref="ORCA 6.1 manual §2.6.1 (Convergence Tolerances, %scf example block)",
    url=_MANUAL_URL_ROOT,
)
EV_CONVCHECK_SEMANTICS = Evidence(
    kind=EVIDENCE_MANUAL,
    text=_Q_CONVCHECK_SEMANTICS,
    ref="ORCA 6.1 manual §2.6.1 (Convergence Tolerances)",
    url=_MANUAL_URL_ROOT,
)
EV_DENSITY_TOL = Evidence(
    kind=EVIDENCE_MANUAL,
    text=_Q_DENSITY_TOL,
    ref="ORCA 6.1 manual §2.6.1 (Convergence Tolerances)",
    url=_MANUAL_URL_ROOT,
)
EV_MAXITER_REFUSAL = Evidence(
    kind=EVIDENCE_MANUAL,
    text=_Q_MAXITER,
    ref="ORCA 6.1 manual §3.13 (convergence-problems note, CASSCF chapter)",
    url=_MANUAL_URL_CASSCF,
)
EV_CMATRIX = Evidence(
    kind=EVIDENCE_MANUAL,
    text=_Q_CMATRIX,
    ref="ORCA 6.1 manual §2.6.9",
    url=_MANUAL_URL_ROOT,
)
EV_CMATRIX_ANION = Evidence(
    kind=EVIDENCE_MANUAL,
    text=_Q_CMATRIX_ANION,
    ref="ORCA 6.1 manual §1.7.13 (Converging DFT for Open-Shell Transition Metals)",
    url=_MANUAL_URL_TUTORIALS,
)

# --- measured records of this group ----------------------------------------

EV_PSEUDO_CONVERGENCE_RECORD = Evidence(
    kind=EVIDENCE_MEASURED,
    text=(
        "ORCA can report SCF convergence after only 2 cycles for a transition-metal "
        "complex; a very small cycle count is a pseudo-convergence signal. Counter-check "
        "'SCF CONVERGED AFTER N CYCLES', re-run with SlowConv, or use a different guess."
    ),
    ref=(
        "SCINE-stack pit table row 35 (SCINE-stack project notes); full text in "
        "SCINE 能力评估报告.md §8 item 11"
    ),
)
EV_NEVPT2_NEEDS_ORBITALS = Evidence(
    kind=EVIDENCE_MEASURED,
    text=(
        "NEVPT2 depends on converged orbitals, not on converged energies: a flat energy "
        "with drifting orbitals was measured to shift results by about 8 kcal/mol, so "
        "tighten the convergence criteria (gradient criterion, larger macro-iteration "
        "limit) instead of trusting an energy-only verdict."
    ),
    ref=(
        "Group's lessons-learned compilation "
        "(B1_可微分DFT_cjq6/经验教训汇编_20260915.md, 2026-09-15)"
    ),
)
EV_NOT_CONVERGED_MEASURED = Evidence(
    kind=EVIDENCE_MEASURED,
    text=(
        "Measured on fixture scf_noconv.out (ORCA 6.1.1, N2/def2-SVP with '%scf MaxIter 3'): "
        "the output prints 'SCF NOT CONVERGED AFTER 2 CYCLES' and the run aborts in LEANSCF; "
        "the DIIS table holds 3 rows and its last DIIS error is 4.71e-02."
    ),
    ref="Fixture scf_noconv.out (fixtures/orca/README.md)",
)
EV_LONG_CONVERGENCE_MEASURED = Evidence(
    kind=EVIDENCE_MEASURED,
    text=(
        "Measured on fixture fblock_dft_la_complex.out (ORCA 6.1.1, wB97M-V/def2-SVPD): the "
        "SCF needed 64 cycles, the DIIS history was reset once and the converger was "
        "switched to SOSCF."
    ),
    ref="Fixture fblock_dft_la_complex.out (fixtures/orca/README.md)",
)
EV_HARD_RUN_CRITERIA_MEASURED = Evidence(
    kind=EVIDENCE_MEASURED,
    text=(
        "Measured on fixture fblock_dft_la_complex.out (64 cycles, mode label "
        "'Total+1el-Energy'): the SCF CONVERGENCE block prints MAX-Density 1.2283e-04 "
        "(tolerance 1.0000e-07), RMS-Density 6.2238e-06 (5.0000e-09), DIIS Error 2.6225e-03 "
        "(5.0000e-07) and Orbital Rotation 9.8997e-05 (1.0000e-05) above their tolerances, "
        "while the enforced energy criterion is met (7.4660e-09 vs 1.0000e-08). Those rows "
        "are informational under the default check, so the triage does not report them."
    ),
    ref="Fixture fblock_dft_la_complex.out (fixtures/orca/README.md)",
)
EV_CLEAN_RUN_CRITERIA_MEASURED = Evidence(
    kind=EVIDENCE_MEASURED,
    text=(
        "Measured on fixture n2_hf_clean.out (ORCA 6.1.1, N2/def2-SVP HF, 8 cycles, "
        "terminated normally): the SCF CONVERGENCE block prints MAX-Density 5.8573e-05 "
        "(tolerance 1.0000e-05), RMS-Density 8.9188e-06 (1.0000e-06) and DIIS Error "
        "1.8058e-03 (1.0000e-06) above their tolerances while the enforced energy criterion "
        "is met (2.0488e-07 vs 1.0000e-06). A healthy short run must stay silent, which is "
        "why only the criteria the check enforces are triaged."
    ),
    ref="Fixture n2_hf_clean.out (fixtures/orca/README.md)",
)
EV_DIIS_REBOUND_MEASURED = Evidence(
    kind=EVIDENCE_MEASURED,
    text=(
        "Measured on fixture fblock_dft_gd_crash.out (ORCA 6.1.1, Gd complex): the DIIS "
        "error was 3.30e-01 when AO-DIIS was switched on (cycle 12) and rose to 1.30e+00 "
        "(cycle 50); the energy jumped by up to 2.37e+01 Eh between cycles; the DIIS "
        "history was reset twice; AutoTRAH then took over and the run died with SIGSEGV in "
        "TRAHIterator::Solve, without printing any SCF verdict."
    ),
    ref="Fixture fblock_dft_gd_crash.out (fixtures/orca/README.md)",
)
EV_ENERGY_OSCILLATION_MEASURED = Evidence(
    kind=EVIDENCE_MEASURED,
    text=(
        "Measured on fixture generated_yb3_sarc2_trah.out (ORCA 6.1.1, Yb3+/SARC2 "
        "TRAH-CASSCF, the group's own f-block record): the module ended with the "
        "banner 'CASSCF NOT CONVERGED AFTER 783 CYCLES' (OOM inside the orbital "
        "optimisation) and its TRAH iteration table -- the energy rows this rule "
        "reads -- alternates between two values in the tail (steps of +-4.30e-02 Eh, "
        "several direction changes). An alternating tail with non-vanishing steps is "
        "the signature this rule reports; a converging run's alternating endgame "
        "collapses toward its tolerance instead. (Reading that banner exposed a "
        "parser defect, fixed here: the SCF verdict patterns had no left anchor and "
        "matched the 'SCF' inside 'CASSCF' -- the verdicts are now one letter-"
        "boundary apart, and the CASSCF verdict is captured in the casscf section.)"
    ),
    ref="Fixture generated_yb3_sarc2_trah.out (fixtures/orca/README.md)",
)


class ScfRescueError(ValueError):
    """Invalid input for the SCF rescue layer."""


# --- data models ------------------------------------------------------------

@dataclass(frozen=True)
class FixProposal:
    """One corrected input file, as text. The caller decides whether to write it.

    ``changes`` lists the edits applied to the original text (the auditable form of
    "conservative editing"); ``content`` is a complete input file for a *new* file and
    never a patch of the caller's file.
    """

    name: str
    content: str
    rationale: str
    evidence: tuple[Evidence, ...]
    changes: tuple[str, ...] = ()


@dataclass(frozen=True)
class _ScfScan:
    """Measured SCF signals, as a view over the parser's ``scf`` section.

    Everything is optional; the parse-result fields carry the same names (see
    ``parsers/orca.py``). The view exists so the findings below keep a typed interface.
    """

    solver_seen: bool = False
    diis_rows: tuple[tuple[int, float], ...] = ()
    diis_blocks: tuple[tuple[tuple[int, float], ...], ...] = ()
    diis_block_switches: tuple[tuple[int | None, float | None], ...] = ()
    diis_error_at_switch: float | None = None
    diis_switch_cycle: int | None = None
    diis_resets: int = 0
    converger_switches: tuple[str, ...] = ()
    max_abs_energy_step: float | None = None
    energies: tuple[float, ...] = ()
    criteria: tuple[tuple[str, float, float], ...] = ()
    check_mode: int = DEFAULT_CONVCHECK_MODE
    check_mode_source: str = "assumed default (the output does not print the mode)"
    converged: bool | None = None


def _scan_from_result(result: ParseResult) -> _ScfScan:
    """Build the scan view from one parse result (no file access)."""
    scf = result.sections.get("scf") or {}
    return _ScfScan(
        solver_seen=bool(scf.get("solver_seen", False)),
        diis_rows=tuple((int(row[0]), float(row[1])) for row in scf.get("diis_rows") or ()),
        diis_blocks=tuple(
            tuple((int(row[0]), float(row[1])) for row in block)
            for block in scf.get("diis_blocks") or ()
        ),
        diis_block_switches=tuple(
            (
                None if row[0] is None else int(row[0]),
                None if row[1] is None else float(row[1]),
            )
            for row in scf.get("diis_block_switches") or ()
        ),
        diis_error_at_switch=scf.get("diis_error_at_switch"),
        diis_switch_cycle=scf.get("diis_switch_cycle"),
        diis_resets=int(scf.get("diis_resets") or 0),
        converger_switches=tuple(scf.get("converger_switches") or ()),
        max_abs_energy_step=scf.get("max_abs_energy_step"),
        energies=tuple(float(value) for value in scf.get("energies") or ()),
        criteria=tuple(tuple(row) for row in scf.get("criteria") or ()),
        check_mode=int(scf.get("check_mode", DEFAULT_CONVCHECK_MODE)),
        check_mode_source=str(
            scf.get("check_mode_source", "assumed default (the output does not print the mode)")
        ),
        converged=scf.get("converged"),
    )

# --- findings ---------------------------------------------------------------

def _sci(value: float) -> str:
    return f"{value:.3e}"


def _not_converged_finding(
    cycles: int | None, terminated: bool | None, scan: _ScfScan
) -> Finding:
    printed = (
        f"SCF NOT CONVERGED AFTER {cycles} CYCLES"
        if cycles is not None
        else "SCF NOT CONVERGED"
    )
    # Only claim an abort when the parser says so; "unknown" must not become "aborted".
    tail = " and the run aborted inside the SCF module" if terminated is False else ""
    measured = ""
    if scan.diis_rows:
        measured = (
            f" Measured from the DIIS table: {len(scan.diis_rows)} cycle(s), last DIIS "
            f"error {_sci(scan.diis_rows[-1][1])}."
        )
    return Finding(
        severity="error",
        message=f"SCF reported no convergence: ORCA printed '{printed}'{tail}.{measured}",
        evidence=(
            EV_NOT_CONVERGED_MEASURED,
            EV_BETTER_GUESS,
            EV_PRESCF_ROUTE,
            EV_CONVFORCED,
        ),
        suggested_fix=(
            "Fix path, in the manual's order: (1) give the SCF better starting orbitals -- "
            "the manual says the answer is almost always better starting orbitals, see the "
            "'prescf' proposal; (2) if you want the damping keyword, take its caution with "
            "it, see the 'slowconv' proposal; (3) for a stuck DIIS the manual points to "
            "SOSCF or TRAH and to large level shifts."
        ),
        refusals=(
            'Raising MaxIter is not the remedy: the manual states "Increasing MaxIter will '
            'not help in many cases." Choose a better starting guess instead.',
            "No energy, density or property may be taken from this output as an SCF "
            "result; the wavefunction is not converged.",
        ),
        rule_id=RULE_NOT_CONVERGED,
    )


def _aborted_no_verdict_finding(scan: _ScfScan, errors: Sequence[str]) -> Finding:
    detail = f" First error line: {errors[0]!r}." if errors else ""
    switches = (
        f" Converger switch(es) seen before the end: {', '.join(scan.converger_switches)}."
        if scan.converger_switches
        else ""
    )
    return Finding(
        severity="error",
        message=(
            "The run ended inside the SCF module without any convergence verdict: the SCF "
            f"solver started ({len(scan.diis_rows)} iteration row(s) read) and the job did "
            "not terminate normally, but neither 'SCF CONVERGED' nor 'SCF NOT CONVERGED' "
            f"was printed.{switches}{detail}"
        ),
        evidence=(EV_DIIS_REBOUND_MEASURED, EV_TRAH, EV_NOTRAH),
        suggested_fix=(
            "Read the error lines first: a crash inside the SCF module is a program-level "
            "failure, not a convergence setting -- if it is reproducible, report it with "
            "the input and this output. Since ORCA switches TRAH on automatically when "
            "convergence problems appear, the manual documents how to disable that "
            "(section 2.6.7), which is the first thing to try if the crash is in the TRAH "
            "solver. If the job was killed from outside, simply resubmit it."
        ),
        refusals=(
            "No SCF verdict means no SCF numbers: do not read energy, density or orbitals "
            "from this output.",
        ),
        rule_id=RULE_ABORTED_NO_VERDICT,
    )


def _pseudo_convergence_finding(cycles: int, energies: Sequence[float]) -> Finding:
    last = f" Last energy read: {energies[-1]:.10f} Eh." if energies else ""
    return Finding(
        severity="warn",
        message=(
            f"SCF reported convergence after only {cycles} cycle(s) "
            f"({PSEUDO_CYCLE_LIMIT} or fewer), which our records treat as a "
            f"pseudo-convergence signal.{last}"
        ),
        evidence=(EV_PSEUDO_CONVERGENCE_RECORD, EV_SLOWCONV_CAUTION),
        suggested_fix=(
            "Counter-check 'SCF CONVERGED AFTER N CYCLES' and the orbital energies, then "
            "re-run once with the damping keyword from the 'slowconv' proposal (or from a "
            "different guess) and compare the two solutions before using either."
        ),
        refusals=(
            "Do not treat a 1-3 cycle convergence as evidence that the SCF is right; it "
            "only says the convergence criterion was met.",
        ),
        rule_id=RULE_PSEUDO_CONVERGENCE,
    )


def _long_convergence_finding(cycles: int, scan: _ScfScan) -> Finding:
    details = []
    if scan.diis_rows:
        details.append(
            f"{len(scan.diis_rows)} DIIS cycles (last DIIS error "
            f"{_sci(scan.diis_rows[-1][1])})"
        )
    if scan.diis_resets:
        details.append(f"{scan.diis_resets} DIIS history reset(s)")
    if scan.converger_switches:
        details.append("converger switch(es): " + ", ".join(scan.converger_switches))
    tail = f" Measured: {'; '.join(details)}." if details else ""
    return Finding(
        severity="warn",
        message=(
            f"SCF needed {cycles} cycles to converge ({LONG_CYCLE_LIMIT} or more is our "
            f"triage threshold).{tail}"
        ),
        evidence=(
            EV_LONG_CONVERGENCE_MEASURED,
            EV_BETTER_GUESS,
            EV_LOOK_AT_ORBITALS,
            EV_PRESCF_ROUTE,
        ),
        suggested_fix=(
            "Before the next run, give the SCF better starting orbitals: the manual's "
            "two-step route (cheap GGA small-basis pre-SCF, then read those orbitals) is "
            "spelled out in the 'prescf' proposal; also look at the guess orbitals "
            "(Print[P_GuessOrb]=1) once."
        ),
        rule_id=RULE_LONG_CONVERGENCE,
    )


def _diis_rebound_finding(
    peak_cycle: int, peak: float, reference: float, scan: _ScfScan
) -> Finding:
    steps = (
        f" The largest |Delta-E| in the DIIS phase was {_sci(scan.max_abs_energy_step)} Eh."
        if scan.max_abs_energy_step is not None
        else ""
    )
    resets = (
        f" The DIIS history was reset {scan.diis_resets} time(s)."
        if scan.diis_resets
        else ""
    )
    return Finding(
        severity="warn",
        message=(
            f"DIIS error rose again while AO-DIIS was active: it was {_sci(reference)} when "
            f"AO-DIIS was switched on and reached {_sci(peak)} later (cycle {peak_cycle}), "
            f"a factor of {peak / reference:.1f} (our threshold: "
            f"{DIIS_REBOUND_FACTOR:.1f}).{steps}{resets}"
        ),
        evidence=(EV_DIIS_REBOUND_MEASURED, EV_DIIS_STUCK, EV_TRAH, EV_LEVELSHIFT),
        suggested_fix=(
            "This is the situation the manual describes as DIIS getting stuck: put SOSCF or "
            "TRAH in operation from here, and/or use a large level shift (sections 2.6.9 and "
            "2.6.4). Watch energy and error together; if the error keeps rising, fix the "
            "starting orbitals instead of the converger -- see the 'prescf' proposal."
        ),
        refusals=(
            "An SCF whose DIIS error is rising must not be left running as if it were "
            "converging: judge it on the error, not on the energy alone.",
        ),
        rule_id=RULE_DIIS_REBOUND,
    )


def _oscillation_finding(scan: _ScfScan) -> Finding | None:
    """The energy tail alternates without collapsing while the run has no verdict.

    Deterministic surrogate of "the iteration is bouncing": over the last
    ``OSCILLATION_WINDOW`` reported steps, count the direction changes of the
    successive differences; report when at least ``OSCILLATION_MIN_FLIPS`` of them
    appear, the largest step is above the noise floor, and the final step still
    carries at least ``OSCILLATION_PERSISTENCE`` of that largest step.  A converging
    run's alternating endgame collapses toward the tolerance and fails the last
    condition; a converged run is excluded outright (its verdict is trusted).
    """
    energies = scan.energies
    if scan.converged is True or len(energies) < 5:
        return None
    window = energies[-OSCILLATION_WINDOW:]
    diffs = [b - a for a, b in zip(window, window[1:])]
    flips = sum(1 for a, b in zip(diffs, diffs[1:]) if a * b < 0)
    peak = max((abs(value) for value in diffs), default=0.0)
    last = abs(diffs[-1]) if diffs else 0.0
    if flips < OSCILLATION_MIN_FLIPS or peak < OSCILLATION_FLOOR:
        return None
    if last < OSCILLATION_PERSISTENCE * peak:
        return None
    return Finding(
        severity="warn",
        message=(
            f"The iteration energy is oscillating: over the last {len(window)} reported steps "
            f"the energy changed direction {flips} time(s), and the final step still "
            f"carries {_sci(last)} Eh against a largest step of {_sci(peak)} Eh -- the "
            "iteration is bouncing instead of converging."
        ),
        evidence=(
            EV_ENERGY_OSCILLATION_MEASURED,
            EV_BETTER_GUESS,
            EV_LEVELSHIFT,
            EV_DIIS_STUCK,
        ),
        suggested_fix=(
            "Do not leave an oscillating SCF running: the manual's first answer for poor "
            "SCF convergence is better starting orbitals (see the 'prescf' proposal); on "
            "the converger side the manual offers large level shifts (shift 0.5, erroff "
            "0) and SOSCF/TRAH in place of plain DIIS; damping with a deliberately small "
            "DampErr is the third lever."
        ),
        refusals=(
            "An oscillating SCF must not be accepted because a single cycle finally "
            "printed a small energy change: judge it on the trajectory.",
        ),
        rule_id=RULE_ENERGY_OSCILLATION,
    )


def _enforced_criteria(scan: _ScfScan) -> tuple[tuple[str, float, float], ...]:
    """The printed criteria the active ConvCheckMode actually enforces.

    Mode 2 (ORCA's default) enforces the energy change; mode 0 enforces every criterion;
    mode 1 stops as soon as one criterion is met, so no single printed row can be read as a
    failure there. The other rows of the block are informational either way.
    """
    if scan.check_mode == CONVCHECK_ALL:
        return scan.criteria
    if scan.check_mode == CONVCHECK_ENERGY:
        return tuple(row for row in scan.criteria if "energy" in row[0].lower())
    return ()


def _unmet_criteria(scan: _ScfScan) -> tuple[tuple[str, float, float], ...]:
    return tuple(
        (name, value, tolerance)
        for name, value, tolerance in _enforced_criteria(scan)
        if tolerance > 0 and value > tolerance * CRITERIA_MARGIN
    )


def _criteria_unmet_finding(
    unmet: Sequence[tuple[str, float, float]], scan: _ScfScan
) -> Finding:
    listed = "; ".join(
        f"{name} {_sci(value)} vs tolerance {_sci(tolerance)}"
        for name, value, tolerance in unmet
    )
    if scan.check_mode == CONVCHECK_ALL:
        subject = (
            f"{len(unmet)} of {len(scan.criteria)} printed criteria are above their printed "
            f"tolerance"
        )
        scope = f" (mode 0: {scan.check_mode_source})"
        extra = ""
    else:
        subject = (
            "the energy change -- the only criterion the default check enforces -- is above "
            "its printed tolerance"
        )
        scope = f" (mode {scan.check_mode}: {scan.check_mode_source})"
        informational = len(scan.criteria) - len(_enforced_criteria(scan))
        extra = (
            f" The other {informational} printed row(s) are informational in this mode and "
            "are not triaged."
            if informational
            else ""
        )
    return Finding(
        severity="warn",
        message=(
            f"SCF is reported converged, but {subject}, by at least a factor of "
            f"{CRITERIA_MARGIN:.0f}{scope}: {listed}.{extra}"
        ),
        evidence=(
            EV_CONVCHECK_MODES,
            EV_CONVCHECK_SEMANTICS,
            EV_CLEAN_RUN_CRITERIA_MEASURED,
            EV_HARD_RUN_CRITERIA_MEASURED,
            EV_DENSITY_TOL,
            EV_NEVPT2_NEEDS_ORBITALS,
        ),
        suggested_fix=(
            "Re-converge before building anything orbital-dependent on this wavefunction "
            "(CASSCF guess, NEVPT2, properties): a better starting guess first, tighter "
            "criteria second. If you only need the energy, the printed non-energy rows are "
            "informational under the default check -- but note that a mode-0 run may also be "
            "accepted when one criterion is missed while the others are overachieved, so "
            "mode-0 hits need a look at the block by hand."
        ),
        refusals=(
            "Do not report this run as meeting the convergence criteria its check "
            "enforces: the block shows at least one of them above its printed tolerance.",
        ),
        rule_id=RULE_CRITERIA_UNMET,
    )


def triage(result: ParseResult) -> tuple[Finding, ...]:
    """Report what is wrong with the SCF of one parsed output (may be empty).

    Only ORCA outputs are supported. Findings come in a fixed order: no convergence,
    crash without a verdict, pseudo-convergence, long convergence, rising DIIS error,
    energy oscillation, unmet convergence criteria; ``build_report`` re-sorts them by
    severity for display.
    """
    if result.program != PROGRAM_ORCA:
        raise ScfRescueError(
            f"SCF triage is implemented for {PROGRAM_ORCA!r} only, got {result.program!r}. "
            "Next step: add a reader for that program's SCF tables in this module, or run "
            "the triage on the ORCA output of the job."
        )
    sections = result.sections
    scf = sections.get("scf") or {}
    converged = scf.get("converged")
    cycles = scf.get("cycles")
    energies = tuple(scf.get("energies") or ())
    terminated = sections.get("terminated_normally")
    errors = tuple(sections.get("errors") or ())
    scan = _scan_from_result(result)

    findings: list[Finding] = []
    if converged is False:
        findings.append(_not_converged_finding(cycles, terminated, scan))
    elif converged is None and scan.solver_seen and terminated is False:
        findings.append(_aborted_no_verdict_finding(scan, errors))

    if converged is True and isinstance(cycles, int):
        if cycles <= PSEUDO_CYCLE_LIMIT:
            findings.append(_pseudo_convergence_finding(cycles, energies))
        elif cycles >= LONG_CYCLE_LIMIT:
            findings.append(_long_convergence_finding(cycles, scan))

    if scan.diis_rows:
        # The analysis runs on the FINAL SCF's table (measured: a geometry optimisation
        # prints one iteration table per SCF cycle and the cycle count restarts at 1,
        # so the last block describes the final SCF).  The reference is the DIIS error
        # when AO-DIIS was switched on in that block (or its first row if the marker is
        # absent); the peak is taken *after* that point only, so the large error of the
        # pre-DIIS start-up cycles cannot be mistaken for a rebound.
        block = scan.diis_blocks[-1] if scan.diis_blocks else scan.diis_rows
        if len(block) < 2:
            block = scan.diis_rows
        block_switch = (
            scan.diis_block_switches[-1]
            if len(scan.diis_block_switches) == len(scan.diis_blocks) and scan.diis_blocks
            else (scan.diis_switch_cycle, scan.diis_error_at_switch)
        )
        switch_cycle, reference = block_switch
        if switch_cycle is None:
            switch_cycle, reference = block[0]
        elif reference is None:
            reference = block[0][1]
        later = [row for row in block if row[0] > switch_cycle]
        if later and reference > 0:
            peak_cycle, peak = max(later, key=lambda row: row[1])
            if peak >= reference * DIIS_REBOUND_FACTOR:
                findings.append(_diis_rebound_finding(peak_cycle, peak, reference, scan))

    oscillation = _oscillation_finding(scan)
    if oscillation is not None:
        findings.append(oscillation)

    unmet = _unmet_criteria(scan)
    if unmet:
        findings.append(_criteria_unmet_finding(unmet, scan))

    return tuple(findings)


# --- fix proposals ----------------------------------------------------------

_COMPOUND_MARKER = "%compound"
_PRESCF_SIMPLE_LINE = "! BP86 def2-SVP def2/J RI LooseSCF SlowConv"
_PRESCF_SCF_BLOCK = (
    "%scf",
    "  # manual section 2.6.9 for this step: 'Set Convergence=Loose and MaxIter=200 or so.'",
    "  # \"Convergence=Loose\" is carried by !LooseSCF on the simple-input line above.",
    "  MaxIter 200",
    "end",
)
_PRESCF_MOINP = '%moinp "prescf.gbw"'

#: Findings that motivate each proposal. A proposal is skipped when the input already
#: contains what it would add (an existing damping keyword, or an orbital read-in).
_SLOWCONV_TRIGGERS = frozenset({RULE_NOT_CONVERGED, RULE_PSEUDO_CONVERGENCE})
_TRAH_TRIGGERS = frozenset({RULE_DIIS_REBOUND, RULE_ENERGY_OSCILLATION})
_PRESCF_TRIGGERS = frozenset(
    {
        RULE_NOT_CONVERGED,
        RULE_ABORTED_NO_VERDICT,
        RULE_LONG_CONVERGENCE,
        RULE_DIIS_REBOUND,
        RULE_ENERGY_OSCILLATION,
        RULE_CRITERIA_UNMET,
    }
)


def _comment(items: Iterable[str], width: int = 92) -> list[str]:
    """Wrap text into ``#`` comment lines (ORCA ignores ``#`` comments)."""
    out: list[str] = []
    for item in items:
        out.extend(textwrap.wrap(item, width=width) or [""])
    return [f"# {line}".rstrip() for line in out]


def _quote(text: str) -> str:
    """Quote a manual sentence for a *generated input file*.

    ASCII double quotes only: a generated file must stay byte-safe for parsers that are
    stricter than ORCA's, and the manual's own typography lives in the Evidence text.
    """
    return f'"{text}"'


def _lines_of(input_text: str) -> list[str]:
    if not input_text.strip():
        raise ScfRescueError(
            "the input text is empty. Next step: pass the content of the original ORCA "
            "input file (the .inp you ran), not the output file."
        )
    return input_text.splitlines()


def _simple_input_index(lines: Sequence[str]) -> int | None:
    """Index of the simple-input line (``!``) before any ``%compound`` block, if any."""
    for index, line in enumerate(lines):
        stripped = line.strip()
        if _COMPOUND_MARKER in stripped:
            return None
        if not stripped or stripped.startswith("#") or stripped.startswith("%"):
            continue
        if stripped.startswith("!"):
            return index
    return None


def _scf_block_range(lines: Sequence[str]) -> tuple[int, int] | None:
    """Inclusive line range of the ``%scf`` block, one-line form included."""
    for index, line in enumerate(lines):
        stripped = line.strip().lower()
        if not stripped.startswith("%scf"):
            continue
        if stripped[len("%scf"):].strip().endswith("end"):
            return index, index
        for end in range(index + 1, len(lines)):
            if lines[end].strip().lower() == "end":
                return index, end
        return index, len(lines) - 1
    return None


def _dedupe_evidence(items: Iterable[Evidence]) -> tuple[Evidence, ...]:
    seen: set[tuple[str, str, str]] = set()
    out: list[Evidence] = []
    for item in items:
        key = (item.kind, item.text, item.ref)
        if key not in seen:
            seen.add(key)
            out.append(item)
    return tuple(out)


def _evidence_from(
    findings: Sequence[Finding], extra: Sequence[Evidence]
) -> tuple[Evidence, ...]:
    return _dedupe_evidence(
        [item for finding in findings for item in finding.evidence] + list(extra)
    )


def _slowconv_proposal(lines: Sequence[str], findings: Sequence[Finding]) -> FixProposal:
    index = _simple_input_index(lines)
    if index is None:
        raise ScfRescueError(
            "no simple-input line (a line starting with '!') was found before any "
            "%compound block. Next step: put the method and basis keywords on a "
            "simple-input line and append 'SlowConv' by hand -- editing a line this "
            "reader cannot see could silently change the wrong job."
        )
    original = lines[index].strip()
    header = _comment(
        [
            "fBlockKit SCF rescue: SlowConv variant. This is a new file; the original "
            "input is unchanged.",
            "'SlowConv' was appended to the simple-input line below. Every other keyword "
            "and block is exactly as you wrote it, and MaxIter is untouched.",
            f"Caution (ORCA 6.1 manual, section 1.7.13): {_quote(_Q_SLOWCONV)}",
        ]
    )
    content = header + [lines[index].rstrip() + " SlowConv"] + list(lines[index + 1:])
    return FixProposal(
        name=FIX_SLOWCONV,
        content="\n".join(content) + "\n",
        rationale=(
            "Appends the manual's damping keyword without touching anything else, so the "
            "result can be compared with the original run. The caution is part of the "
            "proposal: SlowConv is documented to be able to converge to a local minimum "
            "solution that is closer to the initial guess, so if the starting orbitals are "
            "the real problem this proposal is the weaker of the two."
        ),
        evidence=_evidence_from(findings, (EV_SLOWCONV_CAUTION, EV_MAXITER_REFUSAL)),
        changes=(
            f"appended 'SlowConv' to the simple-input line: {original!r} -> "
            f"{original + ' SlowConv'!r}",
            "prepended comment lines (purpose + the manual's SlowConv caution)",
        ),
    )


def _trah_proposal(lines: Sequence[str], findings: Sequence[Finding]) -> FixProposal:
    index = _simple_input_index(lines)
    if index is None:
        raise ScfRescueError(
            "no simple-input line (a line starting with '!') was found before any "
            "%compound block. Next step: add a simple-input line carrying your method "
            "and basis keywords plus 'TRAH' by hand."
        )
    original = lines[index].strip()
    has_casscf = any(
        line.lstrip().lower().startswith("%casscf") for line in lines
    )
    header = _comment(
        [
            "fBlockKit SCF rescue: TRAH variant. This is a new file; the original input "
            "is unchanged.",
            "'TRAH' was appended to the simple-input line below. Every other keyword and "
            "block is exactly as you wrote it, and MaxIter is untouched.",
            f"Source (ORCA 6.1 manual, section 2.6.4): {_quote(_Q_TRAH_WHEN)}",
            "Note: in an SCF job TRAH needs nothing else; AutoTRAH is on by default, so "
            "this keyword mainly matters when NOTRAH was set. A CASSCF job's TRAH route "
            "additionally needs a matching /C auxiliary basis (the recipe layer's rule).",
        ]
    )
    content = header + [lines[index].rstrip() + " TRAH"] + list(lines[index + 1:])
    return FixProposal(
        name=FIX_TRAH,
        content="\n".join(content) + "\n",
        rationale=(
            "The manual's robust second-order SCF: TRAH uses the electronic-Hessian "
            "information and is described for exactly the situation this run shows "
            "(troublesome convergence, a rising or oscillating DIIS trajectory). Appended "
            "to the simple line so the result can be compared with the original run."
            + (" This input carries a %casscf block: the CASSCF TRAH route additionally "
               "needs a matching /C auxiliary basis set." if has_casscf else "")
        ),
        evidence=_evidence_from(findings, (EV_TRAH_WHEN, EV_TRAH_WHEN_SOSCF, EV_TRAH)),
        changes=(
            f"appended 'TRAH' to the simple-input line: {original!r} -> "
            f"{original + ' TRAH'!r}",
            "prepended comment lines (purpose + the manual's TRAH/AutoTRAH statement)"
            + (" and the CASSCF /C auxiliary-basis note" if has_casscf else ""),
        ),
    )


def _prescf_proposal(lines: Sequence[str], findings: Sequence[Finding]) -> FixProposal:
    index = _simple_input_index(lines)
    if index is None:
        raise ScfRescueError(
            "no simple-input line (a line starting with '!') was found before any "
            "%compound block. Next step: build the pre-SCF by hand -- copy your input, "
            f"replace the method line with {_PRESCF_SIMPLE_LINE!r}, and read the result "
            f"back with 'moread' plus {_PRESCF_MOINP!r}."
        )
    body = list(lines)
    body[index] = _PRESCF_SIMPLE_LINE
    changes = [
        "replaced the simple-input line with the manual's pre-SCF method: "
        f"{_PRESCF_SIMPLE_LINE!r}",
    ]
    block = _scf_block_range(body)
    if block is None:
        body[index + 1:index + 1] = list(_PRESCF_SCF_BLOCK)
        changes.append(
            "inserted a %scf block with the manual's pre-SCF setting (MaxIter 200); the "
            "input had none"
        )
    else:
        start, end = block
        kept = ["", "# --- your %scf block, kept for reference (step 2 uses it unchanged) ---"]
        kept += [f"# {line}".rstrip() for line in body[start : end + 1]]
        body[start : end + 1] = kept + list(_PRESCF_SCF_BLOCK)
        changes.append(
            "replaced the %scf block with the manual's pre-SCF setting (MaxIter 200); your "
            "%scf lines are kept as comments and are carried unchanged into step 2"
        )

    step2 = [f"# {line}".rstrip() if line.strip() else "#" for line in lines]
    step2[index] = "# " + lines[index].rstrip() + " moread"
    step2.insert(index + 1, "# " + _PRESCF_MOINP)
    changes.append(
        "appended step 2 as a commented block: your input with 'moread' on the simple-input "
        f"line and {_PRESCF_MOINP!r} after it"
    )

    header = _comment(
        [
            "fBlockKit SCF rescue: two-step pre-SCF route (manual section 2.6.9). This is "
            "a new file; the original input is unchanged.",
            "Step 1 (this file): save it as prescf.inp and run it. It is deliberately "
            "crude -- small basis, pure GGA, RI, loose convergence, strong damping -- and "
            "its only job is to produce prescf.gbw, i.e. starting orbitals good enough for "
            "the real job.",
            f"Manual: {_quote(_Q_PRESCF_ROUTE)}",
            "Step 2 (below, commented out): your own input with 'moread' and the %moinp "
            "line added, saved under a different name and run after step 1. Keep the "
            "basenames different from step 1's -- the manual's section 1.7.12 covers the "
            "AutoStart pitfalls.",
            "Step 2 keeps your %scf settings exactly as written. MaxIter is not a remedy "
            f"in itself (manual section 3.13: {_quote('Increasing MaxIter will not help in many cases.')})"
            " -- the only MaxIter in this file is the pre-SCF value section 2.6.9 "
            "prescribes for step 1.",
        ]
    )
    step2_header = _comment(
        [
            "Step 2 of 2 -- the real job, started from the pre-SCF orbitals. Save these "
            "lines as a new file and run them after step 1 has finished.",
            f"Guess guidance (manual section 2.6.9): {_quote(_Q_CMATRIX)} Section 1.7.13 "
            f"adds: {_quote(_Q_CMATRIX_ANION)} Uncomment a GuessMode line in your %scf "
            "block if that applies to your system.",
        ]
    )
    content = header + body + [""] + step2_header + step2
    return FixProposal(
        name=FIX_PRESCF,
        content="\n".join(content) + "\n",
        rationale=(
            "Implements the manual's two-step route: the cheap pre-SCF exists only to "
            "supply better starting orbitals, which the manual calls the answer almost "
            "always when SCF convergence is poor. The target job is preserved verbatim "
            "(including your %scf settings) and only gains 'moread' plus the %moinp line; "
            "the damping keyword used in step 1 is the manual's own pre-SCF ingredient, "
            "with its caution in view. Nothing here promises convergence: compare the "
            "result with the original run before using it."
        ),
        evidence=_evidence_from(
            findings, (EV_PRESCF_ROUTE, EV_BETTER_GUESS, EV_CMATRIX, EV_MAXITER_REFUSAL)
        ),
        changes=tuple(changes),
    )


def propose_fixes(
    input_text: str, findings: Sequence[Finding]
) -> tuple[FixProposal, ...]:
    """Build corrected inputs for the SCF findings (may be empty).

    ``input_text`` is the content of the original ``.inp``; the returned proposals are
    complete files for *new* names, so the original file on disk is never touched. A
    proposal is skipped when the input already contains what it would add (an existing
    damping keyword, or an existing ``moread``/``%moinp`` orbital read).
    """
    lines = _lines_of(input_text)
    selected = {finding.rule_id for finding in findings}
    damping_present = any(
        "SlowConv" in line for line in lines if line.strip().startswith("!")
    )
    orbital_read_present = any(
        "moread" in line.lower() or line.strip().lower().startswith("%moinp")
        for line in lines
    )
    proposals: list[FixProposal] = []
    if selected & _SLOWCONV_TRIGGERS and not damping_present:
        proposals.append(_slowconv_proposal(lines, findings))
    if selected & _PRESCF_TRIGGERS and not orbital_read_present:
        proposals.append(_prescf_proposal(lines, findings))
    trah_present = any(
        "trah" in line.lower() for line in lines if line.strip().startswith("!")
    )
    if selected & _TRAH_TRIGGERS and not trah_present:
        proposals.append(_trah_proposal(lines, findings))
    return tuple(proposals)
