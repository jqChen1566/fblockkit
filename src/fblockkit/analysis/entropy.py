"""A2: occupation-number spectrum + single-orbital entropy *bound* + plateau detection.

What the literature quantity is
-------------------------------
The single-orbital entropy of the autoCAS protocol (Stein & Reiher) is defined from the
one-orbital reduced density matrix of an orbital, in the basis of its four occupation
states {|0>, |up>, |down>, |up-down>}:

    s_i = -sum_{j in {0,up,down,2}} w_j ln w_j ,      w_j >= 0,  sum w_j = 1,

bounded above by ln 4 (Wardzala et al. 2026, Eq. (6); Stein & Reiher 2016).  The four
weights are fixed by three numbers: the alpha and beta occupations n_a, n_b and the
double-occupancy probability p2,

    w = [1 - n_a - n_b + p2,  n_a - p2,  n_b - p2,  p2] .

**An ORCA output does not print them.**  The CASSCF ``N(occ)=`` line gives only the
spin-summed natural occupations n = n_a + n_b; n_a/n_b/p2 would come from the 2-RDM
(the two-step route documented in the A6 section: ``orca_2json`` RDM2_aa/ab/bb plus
``orca_loc`` localised orbitals).  Computing the four-state entropy from n alone is
therefore impossible, and this module does not pretend otherwise.

What this module computes instead: an upper bound
-------------------------------------------------
This module computes the **maximum-entropy completion** of the four weights consistent
with the spin-summed occupation n,

    s_bound(n) = -2 [x ln x + (1-x) ln(1-x)] ,        x = n/2,

which is the four-state entropy of the spin-independent (product) distribution obtained
at n_a = n_b = n/2 and p2 = n_a n_b = n**2/4.  The maximum-entropy principle makes the
product distribution the entropic maximum among all admissible (n_a, n_b, p2) with
n_a + n_b = n (the split n/2 is itself imposed by Jensen's inequality on the concave
binary entropy), so

    s_bound(n) >= s_i   for the true four-state entropy of any admissible completion,

with equality when the two spin channels are uncorrelated.  Consequences:

- ``max(s_bound) < 0.14`` -> the true entropy is **also** below the Stein & Reiher line:
  a valid one-sided *exclusion* of the multireference reading on this metric;
- ``max(s_bound) >= 0.14`` -> **not** evidence of multireference character: the bound is
  loose (it is attained only by the uncorrelated completion).  The true entropy then
  needs the 2-RDM route, which this file-based tool does not do;
- a single determinant gives bound 0 for n in {0, 2} (closed-shell case) but bound
  ln 4 for n = 1 -- a spin-restricted description of a single unpaired electron
  maximises the occupation uncertainty even though the true entropy is 0 (polarised
  single occupancy) or ln 2 (equal alpha/beta weights).  The bound is tight only for
  closed-shell occupations; do not read it as the literature quantity.

The *relative* structure of the bound spectrum (plateau detection, weak-correlation
cut-off) is offered as a **ranking proxy** for active-space selection, marked as this
tool's approximation: relative order is not guaranteed to agree with the relative order
of the true four-state entropies.

Thresholds and preconditions (each labelled with its source; mixing them is forbidden):

- multireference line: max s_i > 0.14 (Stein & Reiher, JCTC 2016, 12, 1760) -- applies
  to the true four-state entropy and is used here only as the exclusion line of the
  bound;
- weak-correlation cut-off: below 1-2% of max (same source) -- with no plateau,
  selection follows this line;
- the other threshold set, Zs(1) > 0.2 (Stein & Reiher, Mol. Phys. 2017, 115, 2110),
  means something different and is not used by this module;

- an entropy spectrum must be read in a **localized-orbital basis** (autoCAS: in a
  localized basis an orbital usually has strong mutual information with only one or
  two others) -- in a non-localized or truncated space a "plateau" can be numerical
  noise;
- this group's measurement (t177): after truncating a large space to a small one the
  spectrum degrades into noise (all s1 <~ 2e-6).
"""

from __future__ import annotations

import math
from typing import Sequence

from ..knowledge.models import (
    EVIDENCE_LITERATURE,
    EVIDENCE_MEASURED,
    Evidence,
    ParseResult,
    ReportSection,
)

#: Literature line for the TRUE four-state entropy (max s_i > 0.14, Stein & Reiher
#: 2016).  This module uses it one-sided: only the exclusion direction transfers to the
#: occupation-derived bound.
MULTIREFERENCE_THRESHOLD = 0.14
#: Weak-correlation cut-off of the ranking proxy: 1-2% of max (Stein & Reiher 2016).
WEAK_FRACTION = 0.02

_EVIDENCE_FOUR_STATE = Evidence(
    kind=EVIDENCE_LITERATURE,
    text=(
        "The single-orbital entropy is defined from the one-orbital reduced density matrix "
        "in the basis of the four occupation states {0, up, down, up-down}, "
        "S_i = -sum_j rho_jj ln rho_jj (Eq. (6) of this review); its four weights need the "
        "spin-resolved occupations and the double occupancy, which a spin-summed natural "
        "occupation does not determine."
    ),
    ref=(
        "Wardzala J. J. et al., Chem. Rev., 2026, 126(8), 4592-4618, "
        "DOI 10.1021/acs.chemrev.5c00866 (Eq. (6))"
    ),
    url="https://doi.org/10.1021/acs.chemrev.5c00866",
    bibkey="wardzala2026multireference",
)

_EVIDENCE_THRESHOLD = Evidence(
    kind=EVIDENCE_LITERATURE,
    text=(
        "Multireference criterion: at least one single-orbital entropy above 0.14; "
        "weak-correlation cut-off: orbitals below 1-2% of max(s1). The line applies to the "
        "four-state entropy; this module uses it only as the exclusion line of the "
        "occupation-derived upper bound."
    ),
    ref="Stein C. J., Reiher M., J. Chem. Theory Comput., 2016, 12(4), 1760-1771, DOI 10.1021/acs.jctc.6b00156",
    bibkey="stein2016automated",
)

_EVIDENCE_ZS1 = Evidence(
    kind=EVIDENCE_LITERATURE,
    text=(
        "The Zs(1) thresholds (0.1/0.2) come from a different diagnostic; they must not be mixed "
        "with the 0.14 max(s1) criterion."
    ),
    ref="Stein C. J., Reiher M., Mol. Phys., 2017, 115(17-18), 2110-2119, DOI 10.1080/00268976.2017.1288934",
    bibkey="stein2017measuring",
)

_EVIDENCE_PRECONDITION = Evidence(
    kind=EVIDENCE_LITERATURE,
    text=(
        "\"In a localized basis, a given orbital is often connected to only one or at most to a "
        "few other orbitals by a strong mutual information element.\" - the localized-orbital "
        "precondition for reading an entropy spectrum."
    ),
    ref="Stein C. J., Reiher M., J. Comput. Chem., 2019, 40(25), 2216-2226, DOI 10.1002/jcc.25869",
    bibkey="stein2019autocas",
)

_EVIDENCE_MEASURED = Evidence(
    kind=EVIDENCE_MEASURED,
    text=(
        "Measured by this group (t177): after a large space is truncated to a small one the "
        "entropy spectrum degrades into numerical noise (all s1 <~ 2e-6, and the plateau "
        "recommendation set contains virtual orbitals with occupation 0). The upper-bound "
        "property s_bound(n) >= s_true(n_a, n_b, p2) is verified numerically over the "
        "admissible weight grid in tests/test_analysis.py (2026-09-26)."
    ),
    ref=(
        "B1_可微分DFT_cjq6/t177熵谱失效_文献定案_20260922.md; bound property: "
        "tests/test_analysis.py"
    ),
)


def single_orbital_entropy_bound(occupation: float) -> float:
    """Spin-summed occupation n in [0, 2] -> the upper bound on the four-state
    single-orbital entropy (0 at n = 0 and n = 2, ln 4 at n = 1).

    This is NOT the literature single-orbital entropy: it is the maximum-entropy
    completion of the four occupation weights for a spin-summed occupation (see the
    module docstring); it never under-estimates the true value.
    """
    if occupation < 0.0 or occupation > 2.0:
        raise ValueError(
            f"occupation {occupation} is out of range (it must be in [0, 2]). "
            f"Next step: check where the data came from."
        )
    x = occupation / 2.0
    if x <= 0.0 or x >= 1.0:
        return 0.0
    return -2.0 * (x * math.log(x) + (1.0 - x) * math.log(1.0 - x))


def entropy_bound_spectrum(occupations: Sequence[float]) -> tuple[float, ...]:
    """The bound spectrum in input order (one value per occupation)."""
    return tuple(single_orbital_entropy_bound(n) for n in occupations)


def _largest_gap(values: Sequence[float]) -> tuple[int, float]:
    """Largest relative gap of a descending sequence: (break position, relative gap)."""
    if len(values) < 2:
        return len(values), 0.0
    gaps = [
        (values[i] - values[i + 1]) / values[i] if values[i] > 0 else 0.0
        for i in range(len(values) - 1)
    ]
    index = max(range(len(gaps)), key=lambda i: gaps[i])
    return index + 1, gaps[index]


def accepts(result: ParseResult) -> bool:
    return bool(result.sections.get("casscf", {}).get("active_occupations"))


def run(
    result: ParseResult,
    *,
    threshold: float = MULTIREFERENCE_THRESHOLD,
    weak_fraction: float = WEAK_FRACTION,
) -> ReportSection:
    """Build the A2 report section (occupation spectrum + entropy bound spectrum +
    exclusion test + plateau proxy + precondition reminders)."""
    occupations = tuple(result.sections.get("casscf", {}).get("active_occupations", ()))
    if not occupations:
        raise ValueError(
            "this output has no CASSCF active-space occupations (no N(occ)= line). "
            "Next step: confirm this is a CASSCF output, or use an output file that "
            "contains an active space."
        )
    spectrum = entropy_bound_spectrum(occupations)
    order = sorted(range(len(spectrum)), key=lambda i: spectrum[i], reverse=True)
    sorted_s = [spectrum[i] for i in order]
    max_s = sorted_s[0] if sorted_s else 0.0
    split, gap = _largest_gap(sorted_s)
    selected = [i for i in order if spectrum[i] > threshold]
    weak_line = max_s * weak_fraction

    lines = [
        f"Active space: {len(occupations)} orbitals; occupations (input order): "
        + " ".join(f"{n:.4f}" for n in occupations),
        "",
        "Single-orbital entropy BOUND spectrum (descending; s_bound = -2[x ln x + "
        "(1-x) ln(1-x)], x = n/2 -- the maximum-entropy completion of the four occupation "
        "weights for a spin-summed occupation):",
        "  " + " ".join(f"{value:.4f}" for value in sorted_s),
        f"  max(s_bound) = {max_s:.4f}",
    ]
    if max_s <= threshold:
        lines.append(
            f"  Exclusion test: max(s_bound) = {max_s:.4f} <= {threshold} -> the true "
            "four-state entropy (which this bound never under-estimates) is also below "
            "the Stein & Reiher line: no orbital shows significant multireference "
            "character in this metric."
        )
    else:
        lines.append(
            f"  Exclusion test: max(s_bound) = {max_s:.4f} > {threshold} -> NOT decidable "
            "from the occupations alone. The bound is loose (it is attained only when the "
            "two spin channels of the orbital are uncorrelated), so it is not evidence of "
            "multireference character; the true four-state entropy needs the 2-RDM route "
            "(orca_2json RDM2_aa/ab/bb plus orca_loc orbitals; see the A6 section)."
        )
    if len(sorted_s) >= 2:
        if gap >= 0.25:
            plateau = order[:split]
            lines.append(
                f"  Plateau/cliff (ranking proxy): the largest relative gap ({gap:.1%}) "
                f"sits at item {split} of the descending order -- the first {split} items "
                f"form the plateau candidate (active-orbital indices, input order from 0: "
                f"{sorted(plateau)}). The proxy ranks orbitals by the bound; relative "
                "order is not guaranteed to match the true four-state entropies."
            )
        else:
            lines.append(
                f"  No significant plateau (the largest relative gap is only {gap:.1%}): "
                f"filtering by the weak-correlation cut-off "
                f"({weak_fraction:.0%} of max = {weak_line:.4f}) leaves "
                f"{len(selected)} candidate(s) (indices: {sorted(selected)}). This is the "
                "ranking proxy, not the literature quantity."
            )
    lines += [
        "",
        f"Threshold-based candidates (s_bound > {threshold}, active-orbital indices): {sorted(selected)}" if selected else "Threshold-based candidates: none",
        "",
        "Preconditions (must travel with the conclusion):",
        "- The entropy spectrum must be read in a localized-orbital basis (autoCAS); a \"plateau\" in a non-localized or truncated space can be numerical noise.",
        "- A closed-shell single determinant gives bound 0 for every orbital; for an open-shell occupation n = 1 the bound is ln 4 = 1.386294 by construction (a spin-restricted description of one unpaired electron), so an f1 active space always reads as \"not decidable\", never as multireference.",
        "- The 0.14 line applies to the true four-state entropy; the threshold 0.14 and Zs(1) > 0.2 mean different things and must not be mixed.",
    ]
    return ReportSection(title="A2 occupation spectrum and entropy bound spectrum", body="\n".join(lines))


def evidence() -> tuple[Evidence, ...]:
    return (
        _EVIDENCE_FOUR_STATE,
        _EVIDENCE_THRESHOLD,
        _EVIDENCE_ZS1,
        _EVIDENCE_PRECONDITION,
        _EVIDENCE_MEASURED,
    )
