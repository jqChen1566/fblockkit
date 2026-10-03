"""Two-step CAS-CI core-excited XAS input generation (menu 37's third mode).

The measured protocol (ORCA 6.1 manual section 3.13.18; every point below
re-measured on the 6.1.1 engine with the [FeCl4]2- L-edge probe,
2026-10-02):

- step 1 optimizes the valence active-space orbitals with a state-averaged
  CASSCF (the manual: valence excited states in the 6-15 eV range);
- step 2 rotates the core orbital(s) into the active window, unfreezes the
  core (``%method FrozenCore FC_NONE``) and solves ONE CAS-CI iteration
  (``maxiter 1``) over an active space that now also contains the core, so
  the saturated space carries the singly core-excited configurations and
  the higher CAS-CI roots are the core-excited states (measured: the Fe 2p
  to 3d L-edge states at 719.36 eV above the ground state of the probe);
- the rotation targets are positional: with ``nel``/``norb`` and the
  unfrozen core, the active window is the ``norb`` consecutive orbitals of
  the orbital list starting at index ``(N_electrons - nel)/2`` -- the
  manual's own example lands on 87 = (185 - 11)/2 for Fe(acac)3, and the
  probe confirms 42 = (96 - 12)/2 for [FeCl4]2-.  Each core orbital moves
  to the leading slots of that window through the five-field
  ``{i, j, 90, 0, 0}`` rotation of the manual's example;
- the run prints the valence transitions and the core transitions (plain
  and SOC-corrected absorption tables).  Measured processing:
  ``orca_mapspc <out> SOCABS`` renders the SOC-corrected table (785 peaks
  on the probe) and ``ABS`` the plain one (19); the ``XAS``/``XASSOC``
  modes do NOT read these tables (measured refusals: "Could not detect an
  XAS Spectrum").
- cross-engine anchor: the probe's first core excitation (719.36 eV) sits
  0.5 eV from the ROCIS first excitation of the same system (718.87 eV,
  ``fixtures/rocis/fecl4_xas.out``) -- two methods, one edge.
- the manual's examples also write ``rel / DoVelocity`` in the step-2
  ``rel`` block; the 6.1.1 input scanner rejects that line ("Unknown
  identifier in CASSCF::REL block ... DOVELOCITY", measured 2026-10-02),
  so the generator omits it -- do not restore it from the manual.

What the module refuses: an element outside its Z table, a physical
valence space, an empty or repeated core-orbital list, an odd electron
bookkeeping (charge/space mismatch), differently sized multiplicity and
root lists, and a window that would start below the first orbital.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..knowledge.elements import ElementError, element_z
from ..knowledge.models import EVIDENCE_MANUAL, EVIDENCE_MEASURED, Evidence

__all__ = [
    "CasciXasError",
    "CasciXasPlan",
    "DEFAULT_KEYWORDS",
    "build_inputs",
    "render_plan",
    "mapspc_recipe",
    "evidence",
]

#: the probe's method line (a light double-zeta pair; the manual's example
#: uses similar defaults)
DEFAULT_KEYWORDS = "def2-SVP def2-SVP/C TightSCF"


class CasciXasError(ValueError):
    """The two-step XAS input cannot be written (with a next step)."""


@dataclass(frozen=True)
class CasciXasPlan:
    """The validated plan and the two generated input texts."""

    charge: int
    multiplicity: int
    valence_nel: int
    valence_norb: int
    core_orbitals: tuple[int, ...]
    window_start: int
    step2_nel: int
    step2_norb: int
    step1_text: str
    step2_text: str
    mode: str = "xas"
    xas_mo: int | None = None
    gbw_name: str = "step1.gbw"


def _parse_state_list(text, label: str) -> tuple[int, ...]:
    try:
        values = tuple(int(token) for token in str(text).replace(",", " ").split())
    except ValueError as exc:
        raise CasciXasError(
            f"the {label} list carries a non-integer ({exc}). Next step: give a "
            "comma-separated list of positive counts."
        ) from exc
    if not values or any(value < 1 for value in values):
        raise CasciXasError(
            f"the {label} list {values} is empty or carries non-positive counts. "
            "Next step: give a comma-separated list of positive counts."
        )
    return values


def _render_coordinates(rows, charge: int, multiplicity: int) -> list[str]:
    lines = [f"*xyz {charge} {multiplicity}"]
    for element, x, y, z in rows:
        lines.append(f"{element} {x:.6f} {y:.6f} {z:.6f}")
    lines.append("*")
    return lines


def build_inputs(
    coordinates,
    *,
    charge,
    multiplicity,
    valence_nel,
    valence_norb,
    step1_nroots,
    core_orbitals,
    step2_mult=None,
    step2_nroots=None,
    keywords=DEFAULT_KEYWORDS,
    step2_keywords=None,
    maxcore=4000,
    gbw_name="step1.gbw",
    mode="xas",
    xas_mo=None,
    nprocs=None,
) -> CasciXasPlan:
    """Validate the two-step request and render both ORCA inputs.

    ``coordinates``: ``(element, x, y, z)`` rows (Angstrom); ``valence_*``:
    the valence active space of step 1 (its electron and orbital counts);
    ``core_orbitals``: the 0-based indices of the core orbital(s) to rotate
    into the window (the step-1 output's orbital table lists them);
    ``step2_mult`` / ``step2_nroots``: the step-2 state-averaging lists
    (multiplicities and roots of the core-excited solve).  Defaults
    ``20,20`` for ``xas`` and ``40,40`` for ``xes``: the manual's saturated
    recipe ``1000,1000`` makes the engine raise the count to the restricted
    CSF space, whose QDPT transition-density stage grows superlinearly and
    did not finish within 24 h for the 290-state [FeCl4]2- case (measured
    2026-10-03: 40 states 8m44s and 60 states 18m35s complete with the
    emission blocks; 20 states are too few for the K-beta core hole to
    enter the listed roots).  Pass a larger list explicitly for small
    systems; ``mode``: ``"xas"`` (section 3.13.18;
    the plain core-saturated CAS-CI) or ``"xes"`` (section 3.13.19; the
    emission side -- the ``refs ras(...)`` saturation with at most one hole
    in the rotated core set, the ``DoDipoleVelocity``/``DecomposeFosc``
    requests and the ``XESSOC``/``XASMOs`` rel block); ``xas_mo`` (XES
    only): the global index of the rotated 1s MO -- defaults to the window
    head (the lowest-index core is rotated into the leading slot);
    ``nprocs``: when given, both inputs carry ``%pal nprocs N end``;
    ``gbw_name``: the step-1 gbw the step-2 ``%moinp`` reads -- give the
    step-1 input's own basename (the menu passes
    ``<stem>.casci_xas.step1.gbw``) so the generated pair runs as written.
    """
    if mode not in ("xas", "xes"):
        raise CasciXasError(
            f"unknown mode {mode!r}. Next step: use 'xas' (manual section "
            "3.13.18) or 'xes' (section 3.13.19)."
        )
    if mode != "xes" and xas_mo is not None:
        raise CasciXasError(
            "XASMOs is the XES mode's index. Next step: drop it, or ask for "
            "mode='xes'."
        )
    rows = tuple(
        (str(element), float(x), float(y), float(z))
        for element, x, y, z in coordinates
    )
    if not rows:
        raise CasciXasError(
            "no atoms were given. Next step: give the structure (an XYZ file)."
        )
    charge = int(charge)
    multiplicity = int(multiplicity)
    valence_nel = int(valence_nel)
    valence_norb = int(valence_norb)
    maxcore = int(maxcore)
    if valence_nel < 1 or valence_norb < 1 or valence_nel > 2 * valence_norb:
        raise CasciXasError(
            f"the valence active space ({valence_nel}e, {valence_norb}o) is not "
            "physical. Next step: check nel/norb of the step-1 valence CASSCF."
        )
    cores = tuple(sorted(int(value) for value in core_orbitals))
    if not cores:
        raise CasciXasError(
            "no core orbital was given. Next step: list the 0-based indices of "
            "the core orbital(s) from the step-1 output's orbital table (for an "
            "L-edge, the metal 2p orbitals sit near -700 eV)."
        )
    if len(set(cores)) != len(cores):
        raise CasciXasError(
            f"the core-orbital list {cores} repeats an index. Next step: give "
            "each core orbital once."
        )
    step2_nel = valence_nel + 2 * len(cores)
    step2_norb = valence_norb + len(cores)
    try:
        n_electrons = sum(element_z(element) for element, _, _, _ in rows) - charge
    except ElementError as exc:
        raise CasciXasError(
            f"the electron count cannot be derived ({exc}). Next step: check the "
            "element symbols of the structure."
        ) from exc
    if (n_electrons - step2_nel) % 2 != 0 or n_electrons < step2_nel:
        raise CasciXasError(
            f"the structure carries {n_electrons} electrons and the step-2 active "
            f"space wants {step2_nel}; the inactive count would not be integer. "
            "Next step: check the charge, the valence space and the core count."
        )
    window_start = (n_electrons - step2_nel) // 2
    if window_start < 1:
        raise CasciXasError(
            f"the step-2 window would start at orbital {window_start}; there is "
            "no room for the inactive space. Next step: check the charge and "
            "the active-space counts."
        )
    xas_mo_checked = None
    if mode == "xes":
        if xas_mo is None:
            xas_mo_checked = window_start
        else:
            if isinstance(xas_mo, float) and not float(xas_mo).is_integer():
                raise CasciXasError(
                    f"XASMOs {xas_mo!r} is not a whole number. Next step: give "
                    "the global index of the rotated 1s MO as an integer."
                )
            try:
                xas_mo_checked = int(xas_mo)
            except (TypeError, ValueError) as exc:
                raise CasciXasError(
                    f"XASMOs {xas_mo!r} is not an integer ({exc}). Next step: "
                    "give the global index of the rotated 1s MO."
                ) from exc
        if not window_start <= xas_mo_checked < window_start + step2_norb:
            raise CasciXasError(
                f"the XASMOs index {xas_mo_checked} is outside the step-2 "
                f"active window {window_start}..{window_start + step2_norb - 1}. "
                "Next step: give the global index of the rotated 1s MO (the "
                "window head when the lowest-index core is rotated first)."
            )
    if nprocs is not None:
        nprocs = int(nprocs)
        if nprocs < 1:
            raise CasciXasError(
                f"nprocs {nprocs} is not a positive count. Next step: give the "
                "number of MPI processes, or leave it out for a serial run."
            )
    if step2_nroots is None:
        step2_nroots = "20,20" if mode == "xas" else "40,40"
    mults = _parse_state_list(step2_mult or f"{multiplicity},{max(multiplicity - 2, 1)}", "multiplicity")
    nroots = _parse_state_list(step2_nroots, "root")
    if len(mults) != len(nroots):
        raise CasciXasError(
            f"the multiplicity list {mults} and the root list {nroots} differ in "
            "length. Next step: give one root count per multiplicity."
        )

    step1_keywords = keywords
    step2_keywords = step2_keywords or keywords
    step1 = [f"!{step1_keywords}"]
    if nprocs is not None:
        step1.append(f"%pal nprocs {nprocs} end")
    step1 += [f"%maxcore {maxcore}", "%casscf"]
    step1 += [
        f"  nel {valence_nel}",
        f"  norb {valence_norb}",
        f"  mult {multiplicity}",
        f"  nroots {int(step1_nroots)}",
        "end",
    ]
    step1 += _render_coordinates(rows, charge, multiplicity)
    step1.append("")

    step2 = [f"!{step2_keywords} MOREAD"]
    if nprocs is not None:
        step2.append(f"%pal nprocs {nprocs} end")
    step2 += [f"%maxcore {maxcore}", f'%moinp "{gbw_name}"']
    step2.append("%scf")
    step2.append("  rotate")
    for offset, core in enumerate(cores):
        step2.append(f"  {{{core},{window_start + offset},90,0,0}}")
    step2 += ["  end", "end", "%method FrozenCore FC_NONE", "end", "%casscf"]
    step2 += [
        f"  nel {step2_nel}",
        f"  norb {step2_norb}",
        "  mult " + ",".join(str(value) for value in mults),
        "  nroots " + ",".join(str(value) for value in nroots),
        "  maxiter 1",
    ]
    if mode == "xes":
        step2 += [
            f"  refs ras({step2_nel}:{len(cores)} 1/{valence_norb}/0 0) end",
            "  DoDipoleVelocity true",
            "  DecomposeFosc true",
            "  rel",
            "    DoSOC true",
            "    XESSOC true",
            f"    XASMOs {xas_mo_checked}",
            "    DoDTensor false",
            "  end",
        ]
    else:
        step2 += [
            "  rel",
            "    DoSOC true",
            "  end",
        ]
    step2.append("end")
    step2 += _render_coordinates(rows, charge, multiplicity)
    step2.append("")

    return CasciXasPlan(
        charge=charge,
        multiplicity=multiplicity,
        valence_nel=valence_nel,
        valence_norb=valence_norb,
        core_orbitals=cores,
        window_start=window_start,
        step2_nel=step2_nel,
        step2_norb=step2_norb,
        step1_text="\n".join(step1),
        step2_text="\n".join(step2),
        mode=mode,
        xas_mo=xas_mo_checked,
        gbw_name=gbw_name,
    )


def mapspc_recipe(out_name: str, *, mode: str = "SOCABS") -> str:
    """The measured rendering line (SOCABS for the SOC table, ABS for the plain)."""
    return f"orca_mapspc {out_name} {mode} -x0<lo> -x1<hi> -w<fwhm> -eV -n<npoints>"


def render_plan(plan: CasciXasPlan, *, step1_path, step2_path) -> str:
    """The two-step checklist, as the menu prints it (XAS and XES variants)."""
    if plan.mode == "xes":
        headline = (
            "RAS-CI core-excited XES input (two-step protocol; ORCA manual 3.13.19):"
        )
        step2_line = (
            f"  step 2: rotate the core orbital(s) {list(plan.core_orbitals)} into "
            f"the window head; active ({plan.step2_nel}e, {plan.step2_norb}o); one "
            "CAS-CI iteration over the single-hole-saturated space (refs ras, "
            "XESSOC rel block)"
        )
    else:
        headline = (
            "CAS-CI core-excited XAS input (two-step protocol; ORCA manual 3.13.18):"
        )
        step2_line = (
            f"  step 2: rotate the core orbital(s) {list(plan.core_orbitals)} into "
            f"the window head; active ({plan.step2_nel}e, {plan.step2_norb}o); one "
            "CAS-CI iteration over the core-saturated space"
        )
    lines = [
        headline,
        f"  step 1: valence SA-CASSCF ({plan.valence_nel}e, {plan.valence_norb}o), "
        f"multiplicity {plan.multiplicity}",
        step2_line,
        f"  step-2 window: orbitals {plan.window_start}.."
        f"{plan.window_start + plan.step2_norb - 1} of the orbital list "
        f"((N_electrons - nel)/2 = {plan.window_start})",
        "",
        f"Step-1 input written: {step1_path}",
        f"Step-2 input written: {step2_path}",
        "",
        "Next steps:",
        f"  1. run step 1 (it writes {plan.gbw_name} next to itself), then step 2 "
        "(its %moinp reads that gbw)",
    ]
    if plan.mode == "xes":
        lines += [
            "  2. the step-2 output carries the valence transitions and, in the "
            "SOC-corrected blocks, the emission spectra -- for the K-beta "
            "protocol (1s + 3p rotated in) the main line sits at 7086.9 eV in "
            "the probe's output (the experimental Fe K-beta1 is 7058 eV; the "
            "def2-SVP level accounts for the offset)",
            "  3. render the emission spectrum:  orca_mapspc <step2-name>.out "
            "XESSOC -w<fwhm> -eV -n<npoints>",
            "     (measured on the probe's 40-root output: 7375 peaks; the "
            "window is the mode's own, 5000-7200 eV, covering the K-beta "
            "region, and -x overrides did not shift it in this mode; the "
            "ABS/SOCABS modes do not read these tables -- measured refusals)",
        ]
        boundaries = [
            "  - the core set selects the edge: 1s + 3p for K-beta emission, or "
            "the 2p group for L-edge XES -- the same protocol covers both",
            "  - the core-orbital indices are SYSTEM-SPECIFIC: read them from "
            "the step-1 output's orbital table (an L-edge 2p sits near -700 eV; "
            "a K-edge 1s near -7000 eV, a 3p near -65 eV)",
            "  - the step-2 root count gates whether the core hole enters the "
            "listed roots (measured on the K-beta probe: 20 roots miss it, the "
            "40,40 default reaches it in 8m44s); the saturated recipe "
            "(1000,1000 -> the restricted CSF count, 290 states here) grows "
            "superlinearly in the QDPT transition-density stage and did not "
            "finish within 24 h (measured 2026-10-03)",
            "  - run the step-2 job in a fresh directory (measured: stale "
            "transition-density residue of a previous run under the same "
            "basename aborts the new run in TDensityContainer)",
            "  - the manual extends the protocol with SC-NEVPT2 on the same "
            "active space (PTMethod SC_NEVPT2); the generated input stops at "
            "the CAS-CI level -- add the PT line when dynamical correlation is "
            "wanted",
        ]
    else:
        lines += [
            "  2. the step-2 output carries the valence transitions and the core "
            "transitions -- the core-excited CAS-CI roots sit far above the ground "
            "state (the probe's Fe 2p to 3d L-edge states at 719.36 eV)",
            "  3. render the spectrum:  "
            + mapspc_recipe("<step2-name>.out", mode="SOCABS"),
            "     for the SOC-corrected table (measured: 785 peaks on the probe), or "
            "the ABS mode for the plain one (19); the XAS/XASSOC modes do not read "
            "these tables (measured refusals)",
        ]
        boundaries = [
            "  - the window selection is positional (measured): with the unfrozen "
            "core the active window is the norb consecutive orbitals starting at "
            "(N_electrons - nel)/2 -- the manual's own example lands on 87 for "
            "Fe(acac)3, the probe on 42 for [FeCl4]2-; the rotation targets are the "
            "leading slots of that window",
            "  - the core-orbital indices are SYSTEM-SPECIFIC: read them from the "
            "step-1 output's orbital table (an L-edge 2p sits near -700 eV; a K-edge "
            "1s near -7000 eV in light elements)",
            "  - enough step-2 roots must be requested for the core manifolds (the "
            "probe: nroots 20,20 reached the L3/L2 region); the manual's Fe(acac)3 "
            "example uses 16,173",
            "  - the manual extends the protocol with SC-NEVPT2 on the same active "
            "space (PTMethod SC_NEVPT2); the generated input stops at the CAS-CI "
            "level -- add the PT line when dynamical correlation is wanted",
            "  - cross-engine anchor (measured): the probe's first core excitation "
            "(719.36 eV) sits 0.5 eV from the ROCIS result of the same system "
            "(718.87 eV) -- two methods, one edge",
        ]
    lines += ["", "Boundaries:"] + boundaries
    return "\n".join(lines)


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the protocol and of its measured mechanism."""
    return (
        Evidence(
            kind=EVIDENCE_MANUAL,
            text=(
                "The two-step CAS-CI/RAS-CI protocol (valence SA-CASSCF first, "
                "then the core orbitals rotated into the active space with "
                "FrozenCore FC_NONE and maxiter 1), the five-field rotate form "
                "and the RAS(1 hole) variant are ORCA manual section 3.13.18 "
                "(XAS/RIXS) with its Fe(acac)3 L-edge example, and section "
                "3.13.19 (XES) for the emission side."
            ),
            ref="ORCA 6.1 manual section 3.13.18 (Core Excited Spectra: CAS-CI/RAS-CI XAS/RIXS)",
            url="https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/CASSCF.html",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "Measured on the [FeCl4]2- L-edge probe (ORCA 6.1.1, 2026-10-02, "
                "fixtures/rocis/fecl4_casci_xas.*): the window mechanism is "
                "positional ((96-12)/2 = 42 and the rotations to 42/43/44 put "
                "the Fe 2p into the active space); the CAS-CI roots carry the "
                "L-edge states at 719.36 eV; orca_mapspc renders the "
                "SOC-corrected table with the SOCABS mode (785 peaks) and the "
                "plain one with ABS (19), while XAS/XASSOC refuse these table "
                "titles."
            ),
            ref="tests/test_casci_xas.py; fixtures/rocis/",
        ),
    )
