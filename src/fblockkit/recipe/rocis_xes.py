"""ROCIS XES input generation (menu 37's writing mode) -- the off-resonance
X-ray emission route of the ROCIS module.

The measured recipe (ORCA 6.1.1, the [FeCl4]2- L-edge probe, 2026-10-01):
a ROCIS run that carries the plain RIXS request prints the off-resonance
XES automatically -- the manual's sentence (section 7.31.4: in every RIXS
requested calculation the off-resonant XES spectrum is automatically
generated), and measured on the probe: the emission tables then carry the
ground-state rows (``<root>-<mult>A -> 0-<mult>A`` at the core-emission
energies) and ``orca_mapspc <out> XES -x0<lo> -x1<hi> -w<fwhm> -eV
-n<npoints>`` renders the spectrum to ``<out>.XES.stk`` and
``<out>.XES.dat`` (measured: 56 peaks from the probe).  A run whose
SOC-corrected RIXS channel is on prints the "XESSOC" flavour as well; the
same tool renders it with the ``XESSOC`` mode (``.XESSOC.stk/.dat``).

The six-element ``OrbWin`` is the RIXS form of the manual (section 7.31.4):
the first four elements are the ranges of the TWO donor spaces, the last
two the acceptor space.  For an L-edge the donors are the spin-orbit-split
core (the probe: ``6,6,7,8`` = the two 2p1/2 and two 2p3/2 orbitals of the
Fe 2p shell in that run's own numbering) and the acceptor is kept wide
(``0,2000``).  Measured boundaries:

- the plain (non-SOC) RIXS channel carries the XES table; a root count too
  small for it is skipped with the engine's own warning and the XES table
  does not print (measured: NRoots 10 skipped the channel, NRoots 30
  covered it, on the probe);
- the SOC-corrected channel (``DoRIXSSOC``) additionally stores per-pair
  transition densities on disk -- measured on the probe at NRoots 30 the
  transient storage grew past 36 GB at about 1 GB/min and the run was not
  practical -- so it defaults off;
- the KHD variant (``DoKHDXESSOC`` with the state lists) is a documented
  termination: the engine computes the emission intensities of every probed
  input and then aborts in its own printout (a NULL stream inside
  ``PrintNXESKHDSpectrum``, signal 11; independent of the state values, the
  list layout and the print level) -- the writer refuses it.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..knowledge.models import EVIDENCE_MANUAL, EVIDENCE_MEASURED, Evidence

__all__ = [
    "RocisXesError",
    "RocisXesPlan",
    "DEFAULT_KEYWORDS",
    "build_input",
    "render_plan",
    "mapspc_recipe",
    "evidence",
]

#: the probe's method line (x2c, the all-electron SVP basis set, the
#: auto-generated auxiliary basis and a tight SCF)
DEFAULT_KEYWORDS = "x2c x2c-SVPall AutoAux TightSCF"


class RocisXesError(ValueError):
    """The XES input cannot be written (with a next step)."""


@dataclass(frozen=True)
class RocisXesPlan:
    """The validated plan and the generated input text."""

    charge: int
    multiplicity: int
    xas_element: int
    nroots: int
    window: tuple[int, int, int, int, int, int]
    do_rixssoc: bool
    do_elastic: bool
    rohf_electrons: int | None
    text: str


def build_input(
    coordinates,
    *,
    charge,
    multiplicity,
    xas_element,
    nroots,
    window,
    do_rixssoc=False,
    do_elastic=True,
    rohf_electrons=None,
    keywords=DEFAULT_KEYWORDS,
    nprocs=8,
    maxcore=4000,
) -> RocisXesPlan:
    """Validate the XES request and render the complete ORCA input.

    ``coordinates``: ``(element, x, y, z)`` rows (Angstrom); ``xas_element``:
    the 0-based position of the core electron's element among the structure's
    atoms (the ``XASelems`` index of the probe, where the Fe of [FeCl4]2- is
    atom 0); ``window``: the six ``OrbWin`` integers (two donor ranges, then
    the acceptor range); ``rohf_electrons``: when given, the input carries
    the probe's ROHF high-spin preparation with that many unpaired electrons
    (``ROHF_NEL[1]``); ``None`` leaves the SCF to the keyword line.
    """
    rows = tuple(
        (str(element), float(x), float(y), float(z))
        for element, x, y, z in coordinates
    )
    if not rows:
        raise RocisXesError(
            "no atoms were given. Next step: give the structure (an XYZ file)."
        )
    try:
        charge = int(charge)
        multiplicity = int(multiplicity)
        xas_element = int(xas_element)
        nroots = int(nroots)
        nprocs = int(nprocs)
        maxcore = int(maxcore)
    except (TypeError, ValueError) as exc:
        raise RocisXesError(f"a non-integer parameter ({exc}).") from exc
    if multiplicity < 1:
        raise RocisXesError(
            f"the multiplicity {multiplicity} is not positive. Next step: give "
            "the structure's spin multiplicity (it also becomes ReferenceMult)."
        )
    if not (0 <= xas_element < len(rows)):
        raise RocisXesError(
            f"the XASelems index {xas_element} is outside the structure's "
            f"{len(rows)} atoms (0-based). Next step: give the position of the "
            "core electron's element among the atoms, counting from zero."
        )
    if nroots < 2:
        raise RocisXesError(
            f"NRoots {nroots} is too small for a spectrum (the ground plus at "
            "least one excited root). Next step: give a positive root count; "
            "the engine skips the plain RIXS channel (the carrier of the XES "
            "table) when the count is too small for it -- measured on the "
            "probe: 10 skipped the channel, 30 covered it."
        )
    window_checked = tuple(int(value) for value in window)
    if len(window_checked) != 6:
        raise RocisXesError(
            f"the OrbWin for a RIXS-requested run carries six integers (two "
            f"donor ranges, then the acceptor range); this one carries "
            f"{len(window_checked)}. Next step: give the two spin-orbit-split "
            "core ranges and a wide acceptor, e.g. 6,6,7,8,0,2000."
        )
    for start, stop, label in (
        (window_checked[0], window_checked[1], "the first donor range"),
        (window_checked[2], window_checked[3], "the second donor range"),
        (window_checked[4], window_checked[5], "the acceptor range"),
    ):
        if min(start, stop) < 0 or start > stop:
            raise RocisXesError(
                f"{label} {start},{stop} is not an increasing non-negative "
                "range. Next step: check the six OrbWin integers."
            )
    rohf_checked = None
    if rohf_electrons is not None:
        rohf_checked = int(rohf_electrons)
        if rohf_checked < 1:
            raise RocisXesError(
                f"the ROHF high-spin electron count {rohf_checked} is not "
                "positive. Next step: leave it empty (no SCF block) or give "
                "the unpaired-electron count."
            )
    if nprocs < 1 or maxcore < 1:
        raise RocisXesError("nprocs and maxcore must be positive.")

    lines = ["# ROCIS XES input (off-resonance X-ray emission; menu 37's writing mode)"]
    lines.append(f"!{keywords}")
    lines.append(f"%pal nprocs {nprocs} end")
    lines.append(f"%maxcore {maxcore}")
    if rohf_checked is not None:
        lines.append("%scf")
        lines.append("  HFTyp ROHF")
        lines.append("  ROHF_CASE HIGHSPIN")
        lines.append(f"  ROHF_NEL[1] {rohf_checked}")
        lines.append("end")
    lines.append("%rocis")
    lines.append("  DoGenROCIS true")
    lines.append(f"  ReferenceMult {multiplicity}")
    lines.append(f"  NRoots {nroots}")
    lines.append("  OrbWin " + ",".join(str(value) for value in window_checked))
    lines.append("  DoPNO true")
    lines.append("  TCutPNO 1e-11")
    lines.append(f"  XASelems {xas_element}")
    lines.append("  rel")
    lines.append("    DoSOC true")
    lines.append("    DoRIXS true")
    lines.append(f"    DoRIXSSOC {'true' if do_rixssoc else 'false'}")
    lines.append(f"    DoElastic {'true' if do_elastic else 'false'}")
    lines.append("  end")
    lines.append("  DoHigherMult true")
    lines.append("  DoLowerMult true")
    lines.append("  DoRI true")
    lines.append("  Decomposefosc true")
    lines.append("end")
    lines.append(f"*xyz {charge} {multiplicity}")
    for element, x, y, z in rows:
        lines.append(f"{element} {x:.6f} {y:.6f} {z:.6f}")
    lines.append("*")
    lines.append("")
    return RocisXesPlan(
        charge=charge,
        multiplicity=multiplicity,
        xas_element=xas_element,
        nroots=nroots,
        window=window_checked,
        do_rixssoc=do_rixssoc,
        do_elastic=do_elastic,
        rohf_electrons=rohf_checked,
        text="\n".join(lines),
    )


def mapspc_recipe(out_name: str, *, mode: str = "XES") -> str:
    """The measured rendering line (``XES``, or ``XESSOC`` for SOC runs)."""
    return (
        f"orca_mapspc {out_name} {mode} -x0<lo> -x1<hi> -w<fwhm> -eV -n<npoints>"
    )


def render_plan(plan: RocisXesPlan, *, path) -> str:
    """The checklist the XES workflow needs, as the menu prints it."""
    d1 = f"{plan.window[0]},{plan.window[1]}"
    d2 = f"{plan.window[2]},{plan.window[3]}"
    acceptor = f"{plan.window[4]},{plan.window[5]}"
    lines = [
        "ROCIS XES input (off-resonance X-ray emission):",
        f"  charge {plan.charge}, multiplicity {plan.multiplicity} "
        f"(also ReferenceMult)",
        f"  XASelems {plan.xas_element} (the 0-based position of the core "
        "element among the atoms)",
        f"  NRoots {plan.nroots}; OrbWin {d1}, {d2}, {acceptor} "
        "(donor 1, donor 2, acceptor)",
        "  plain RIXS channel: requested (it carries the off-resonance XES "
        "table)",
        "  SOC-corrected RIXS channel: "
        + ("requested (the XESSOC flavour)" if plan.do_rixssoc else "off"),
        "  elastic line: " + ("included" if plan.do_elastic else "excluded"),
        "  SCF preparation: "
        + (
            f"ROHF high-spin, ROHF_NEL[1] = {plan.rohf_electrons}"
            if plan.rohf_electrons is not None
            else "none (the keyword line decides)"
        ),
        "",
        f"Input written: {path}",
        "",
        "Next steps:",
        "  1. run:  orca <name>.inp > <name>.out",
        "  2. the run's emission tables carry the ground-state rows "
        "('<root>-<mult>A -> 0-<mult>A') at the core-emission energies "
        "(they appear at the end of the emission tables, after the "
        "between-excited-state rows)",
        "  3. render the spectrum, e.g.:  "
        + mapspc_recipe("<name>.out", mode="XESSOC" if plan.do_rixssoc else "XES"),
        "     -> <name>.out.XES.stk (the sticks) and .XES.dat (the broadened "
        "curve)" + (
            "; a run with the SOC-corrected channel also renders with the "
            "XESSOC mode" if plan.do_rixssoc else ""
        ),
        "",
        "Boundaries:",
        "  - the off-resonance XES is automatic in a RIXS-requested run (the "
        "manual's sentence, measured on the probe); no dedicated keyword is "
        "needed",
        "  - the plain channel needs enough roots: too few and the engine "
        "skips it with its own warning, printing no XES table (measured: 10 "
        "skipped, 30 covered, on the probe)",
        "  - the SOC-corrected channel stores per-pair transition densities "
        "(measured: the transient storage grew past 36 GB at about 1 GB/min "
        "on the probe at NRoots 30) -- requested only for small systems or "
        "small root sets",
        "  - the KHD variant (DoKHDXESSOC and its state lists) is a "
        "documented termination: every probed input computes and then aborts "
        "in the engine's own printout (measured)",
        "  - the window indices are the SYSTEM's own orbital numbers (the "
        "probe's 6,6,7,8,0,2000 are not transferable); read the two "
        "spin-orbit-split core ranges and a wide acceptor from the target "
        "run's own orbital table",
    ]
    return "\n".join(lines)


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the route and of its measured boundaries."""
    return (
        Evidence(
            kind=EVIDENCE_MANUAL,
            text=(
                "The off-resonance XES route and the six-element OrbWin are "
                "section 7.31.4 of the ORCA manual: the RIXS form of OrbWin "
                "carries the two donor ranges first and the acceptor range "
                "last, and every RIXS-requested calculation generates the "
                "off-resonance XES automatically; the spectrum is rendered "
                "from the output with orca_mapspc and the XES/XESSOC "
                "spectrum-type argument."
            ),
            ref="ORCA 6 manual section 7.31.4 (RIXS and the off-resonance XES)",
            url="https://www.faccts.de/docs/orca/6.0/manual/contents/detailed/rocis.html",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "Measured on the [FeCl4]2- L-edge probe (ORCA 6.1.1, "
                "2026-10-01, fixtures/rocis/fecl4_xes.out): the generated "
                "input runs to normal termination and prints the emission "
                "tables with the ground-state rows at the core-emission "
                "energies (744.2/744.3 eV for the two intermediate roots); "
                "orca_mapspc renders 56 peaks (XES mode); the SOC-corrected "
                "variant (DoRIXSSOC true, NRoots 10) renders 2376 peaks with "
                "the XESSOC mode; the plain channel was skipped at NRoots 10 "
                "(the engine's own zero-states warning) and covered at 30; "
                "the KHD variant aborts in the engine's printout for every "
                "probed shape."
            ),
            ref="tests/test_rocis_xes.py; fixtures/rocis/",
        ),
    )
