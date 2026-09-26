"""Rule G: the diffuse-orbital (Rydberg) check for the active space.

The rule (Lee & Rondinelli's decision ladder, rule G): atomic-orbital selection
alone cannot see Rydberg contamination, because a Rydberg state differs from a
valence one by the *spatial extent* of its orbitals rather than by an
occupation pattern.  The check is therefore to evaluate, for every candidate
active orbital, the smallest Gaussian exponent associated with its centre and
angular-momentum channel -- the channel's diffusest function -- and to rank the
orbitals by it; the orbitals whose target channel is built from highly diffuse
functions are the Rydberg candidates, and the remedy is to extend the active
space with them and recompute the target states with SA-CASSCF + NEVPT2.

What the source asks for is identification *and sorting*, and no threshold: the
exponent's absolute value is not transferable between elements (measured: the
def2-SVP diffusest s exponent on N is 0.19, the SARC2 diffusest on Eu is 0.022,
and Eu's *valence* s functions are the diffuse end of that basis).  This module
therefore reports the numbers and the ranking -- the channel's diffusest and
tightest exponent next to the orbital's weight on the channel -- and leaves the
line to the reader, with the source's remedy quoted.

Data: the run's own printed basis (``!PrintBasis``; the block reader lives in
the ORCA parser) and the per-orbital composition table (the same print the
d-solution rule reads).  The composition table labels its atomic orbitals by
angular-momentum letter only, so the weight is the orbital's total weight on
the (element, angular momentum) channel, which is exactly the quantity the
source's wording names; a per-shell split would need the numbered labels of an
``orca_2json`` export (menu 14's data), not the output.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..knowledge.models import EVIDENCE_LITERATURE, EVIDENCE_MEASURED, Evidence, ReportSection

__all__ = ["DiffuseError", "DiffuseAnalysis", "accepts", "analyze", "report", "run", "evidence"]


class DiffuseError(ValueError):
    """The diffuse-orbital check cannot run on the given output (with a next step)."""

#: An orbital counts as active when its printed occupation is not within this
#: distance of 0 or 2 (the same "fractional" convention the window inference uses).
_INTEGER_TOLERANCE = 1e-4

#: The angular-momentum letters the composition table uses, in shell order.
_LETTERS = "spdfgh"


@dataclass(frozen=True)
class DiffuseAnalysis:
    """The per-orbital channel diffuseness, ranked."""

    entries: tuple[dict, ...]
    active_indices: tuple[int, ...]
    basis_elements: tuple[str, ...]
    checks: tuple[str, ...]
    verdicts: tuple[str, ...]


def _channel_of(entry: dict) -> tuple[str, str, float] | None:
    """The (element, angular letter) channel with the largest weight of an orbital."""
    totals: dict[tuple[str, str], float] = {}
    for shell in entry.get("shells", ()):
        key = (shell["element"], shell["shell"])
        totals[key] = totals.get(key, 0.0) + shell["weight"]
    if not totals:
        return None
    (element, letter), weight = max(totals.items(), key=lambda item: item[1])
    return element, letter, weight


def analyze(sections) -> DiffuseAnalysis:
    """Run rule G on one parsed ORCA output (needs !PrintBasis and the composition print)."""
    basis = sections.get("basis") or {}
    composition = sections.get("orbital_composition") or {}
    if not basis.get("present"):
        raise DiffuseError(
            "the output carries no printed basis, so no Gaussian exponent is available. "
            "Next step: rerun (or reprint) the job with the !PrintBasis keyword; the "
            "check needs the shell exponents of the basis actually used."
        )
    if not composition.get("present"):
        raise DiffuseError(
            "the output carries no per-orbital composition table. Next step: add "
            "'%output Print[P_ReducedOrbPopMO_L] 1 end' to the input; the check needs "
            "the orbital weights per centre and angular momentum."
        )
    orbitals = composition.get("orbitals") or ()
    occupations = (sections.get("orbitals") or {}).get("occupations") or ()
    if not occupations or len(occupations) != len(orbitals):
        raise DiffuseError(
            f"the orbital table has {len(occupations)} occupation(s) while the "
            f"composition table has {len(orbitals)} orbital(s). Next step: check that "
            "both prints come from the same run."
        )
    # the exponents per (element, letter) channel: the diffusest and the tightest
    channel_exponents: dict[tuple[str, str], list[float]] = {}
    elements: list[str] = []
    for entry in basis.get("elements", ()):
        elements.append(entry["element"])
        for shell in entry["shells"]:
            key = (entry["element"], shell["angular"])
            channel_exponents.setdefault(key, []).extend(shell["exponents"])

    active = [
        index
        for index, occupation in enumerate(occupations)
        if min(abs(occupation), abs(occupation - 2.0)) > _INTEGER_TOLERANCE
    ]
    entries: list[dict] = []
    for index in active:
        channel = _channel_of(orbitals[index])
        if channel is None:
            continue
        element, letter, weight = channel
        exponents = channel_exponents.get((element, letter))
        if not exponents:
            entries.append(
                {
                    "index": index,
                    "element": element,
                    "shell": letter,
                    "weight": weight,
                    "diffuse": None,
                    "tight": None,
                    "occupation": occupations[index],
                }
            )
            continue
        entries.append(
            {
                "index": index,
                "element": element,
                "shell": letter,
                "weight": weight,
                "diffuse": min(exponents),
                "tight": max(exponents),
                "occupation": occupations[index],
            }
        )
    entries.sort(key=lambda item: (item["diffuse"] is None, item["diffuse"]))
    checks = [
        f"basis: {len(basis.get('elements', ()))} element block(s) "
        f"({', '.join(elements)}); the printed channels are "
        + ", ".join(
            f"{element} {letter} ({len(values)} primitive(s))"
            for (element, letter), values in sorted(channel_exponents.items())
        ),
        f"active orbitals (printed occupation not within {_INTEGER_TOLERANCE:g} of 0 or "
        f"2): {active}",
    ]
    if basis.get("ecp_blocks"):
        checks.append(
            f"{basis['ecp_blocks']} ECP block(s) were skipped: an effective core "
            "potential carries no Gaussian exponents, so its core has no diffuseness to "
            "report"
        )
    verdicts = [
        "ranking (most diffuse channel first) -- the source asks for identification and "
        "sorting, not for a line: the exponent's absolute value does not transfer "
        "between elements, so compare the diffusest exponent with the tightest one of "
        "the *same* channel (both printed) and with the other active orbitals here",
        "for a candidate whose channel is the diffuse end of the basis, the source's "
        "remedy is to extend the active space with it and recompute the target states "
        "with SA-CASSCF + NEVPT2 (its rule G, which exists because an occupation-based "
        "selection cannot see spatial extent)",
    ]
    return DiffuseAnalysis(
        entries=tuple(entries),
        active_indices=tuple(active),
        basis_elements=tuple(elements),
        checks=tuple(checks),
        verdicts=tuple(verdicts),
    )



def accepts(result) -> bool:
    """Whether rule G applies: both the printed basis and the composition table are there."""
    sections = result.sections
    return bool((sections.get("basis") or {}).get("present")) and bool(
        (sections.get("orbital_composition") or {}).get("present")
    )


def run(result) -> ReportSection:
    """The analyser entry point (the contract of the analysis layer)."""
    return report(analyze(result.sections))

def report(analysis: DiffuseAnalysis) -> ReportSection:
    """Render the diffuse-orbital (rule G) section."""
    lines = ["Inputs:"]
    for check in analysis.checks:
        lines.append(f"  - {check}")
    lines += [
        "",
        "Active orbitals by channel diffuseness (diffusest exponent of the orbital's "
        "dominant (element, angular momentum) channel):",
        "  orbital  element shell   weight%   diffusest exp   tightest exp   occupation",
    ]
    for entry in analysis.entries:
        if entry["diffuse"] is None:
            lines.append(
                f"  {entry['index']:>7}  {entry['element']:<7} {entry['shell']:<6} "
                f"{entry['weight']:>7.1f}   (no exponent in the printed basis: an ECP "
                "core?)"
            )
            continue
        lines.append(
            f"  {entry['index']:>7}  {entry['element']:<7} {entry['shell']:<6} "
            f"{entry['weight']:>7.1f}   {entry['diffuse']:>13.8f}   {entry['tight']:>13.3f} "
            f"   {entry['occupation']:>10.6f}"
        )
    if not analysis.entries:
        lines.append(
            "  (no active orbital: the output's occupations are all integer-valued, so "
            "there is no active space to check -- rule G applies to the orbitals of an "
            "active space)"
        )
    lines += ["", "Reading:"]
    for verdict in analysis.verdicts:
        lines.append(f"  - {verdict}")
    lines += [
        "",
        "Boundary: the composition table labels atomic orbitals by angular-momentum "
        "letter only, so the weight above is the orbital's total weight on the channel "
        "-- the quantity the source names -- and the per-shell split would need the "
        "numbered labels of an orca_2json export rather than the output.",
    ]
    return ReportSection(title="A8 diffuse-orbital (Rydberg) check", body="\n".join(lines))


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the rule-G check."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "Rule G of the automated multi-reference decision ladder: Rydberg "
                "contamination is a matter of spatial extent, not of occupations, so the "
                "check evaluates the smallest Gaussian exponent associated with the "
                "relevant atomic centre and angular-momentum channel, identifies the "
                "orbitals with significant contributions from highly diffuse basis "
                "functions as Rydberg candidates and sorts them; the remedy is to extend "
                "the active space with them and recompute the target states with "
                "SA-CASSCF + NEVPT2."
            ),
            ref=(
                "Lee & Rondinelli, arXiv:2609.13357 (the A-H decision ladder, rule G; "
                "the identification is stated there, the threshold is not)"
            ),
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "The exponents come from ORCA's own printed basis (!PrintBasis, the "
                "'BASIS SET IN INPUT FORMAT' block) and are checked on two measured "
                "fixtures: N2/def2-SVP prints 7s4p1d with the diffusest s exponent "
                "0.1876, and Eu/SARC2-DKH-QZVP prints the four f contractions "
                "(4.870, 1.956, 0.786, 0.316) matching the four f shells the export "
                "labels carry. The same measurement is why the check has no default "
                "line: Eu's valence s functions reach 0.0222, below any fixed "
                "'Rydberg' line that would work for N."
            ),
            ref="tests/test_diffuse.py; fixtures/orca/n2_diffuse.out, eu_basis.out",
        ),
    )
