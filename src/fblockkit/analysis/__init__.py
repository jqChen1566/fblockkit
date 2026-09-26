"""Analysis layer: characterisation analysis from an output file / a structure
(architecture design v0.1 §3 L2').

- A1 ``composition``: orbital composition (reads ORCA's per-MO Loewdin composition table);
- A2 ``entropy``: occupation spectrum + single-orbital entropy *bound* spectrum + plateau
  detection (reads the CASSCF active occupations; the bound and its one-sided use are
  defined in the module);
- A3 ``mr_diagnostics``: the multi-reference character panel (T1, fractional occupations,
  the entropy bound as an exclusion test);
- A4 ``crystal_field``: linear fit of the crystal-field parameters B_k^q (Eq. (3)) +
  extended Stevens operators;
- A5 ``cf_declaration``: the projection-basis declaration check for crystal-field
  parameter sets (pure mapping input, called by the menus);
- A6 ``local_spin``: ORCA's local spin analysis (a surrogate for the environment spin
  polarisation entropy Delta S_E);
- S1 ``geometry``: coordination geometry and symmetry hints (reads an XYZ structure);
- A8 ``diffuse``: the diffuse-orbital (Rydberg) check of rule G, ranked per active
  orbital (reads the printed basis and the composition table);
- S2 ``point_charge``: the point-charge crystal-field estimate from a structure.

Common convention: every output-reading module provides ``accepts(input) -> bool`` and
``run(input, ...) -> ReportSection``; each module's ``evidence()`` supplies the
provenance for its conclusions. The analysers only read files and never change the
user's data. A4/A5 read no file (their input is sampled state energies + projection
coefficients + symmetry, or a mapping), S1/S2 read a structure, and they are called
directly by the caller.
"""

from . import (
    cf_declaration,
    composition,
    crystal_field,
    diffuse,
    entropy,
    geometry,
    local_spin,
    mr_diagnostics,
    point_charge,
)
from .composition import OrbitalRow, orbital_rows, shell_ranking
from .crystal_field import (
    CFResult,
    CrystalFieldError,
    allowed_parameters,
    design_matrix,
    fit_crystal_field,
    hamiltonian,
    stevens_matrices,
)
from .entropy import entropy_bound_spectrum, single_orbital_entropy_bound
from .geometry import Atom, coordination_shell, dominant_center, parse_xyz

# Analysers that read a "program output" (only those whose accepts() passes are run);
# the structure-based S1/S2 and the mapping-based A5 are called directly by the caller
_OUTPUT_ANALYZERS = (composition, entropy, mr_diagnostics, local_spin, diffuse)


def run_all(result) -> tuple:
    """Run every applicable analyser on one parse result and return the report sections."""
    sections = []
    for module in _OUTPUT_ANALYZERS:
        if module.accepts(result):
            sections.append(module.run(result))
    return tuple(sections)


def evidence_for(result) -> tuple:
    """Provenance of every analyser that accepts this parse result (merged into the
    report's References block by the caller)."""
    return tuple(
        item
        for module in _OUTPUT_ANALYZERS
        if module.accepts(result)
        for item in module.evidence()
    )


__all__ = [
    "run_all",
    "evidence_for",
    "Atom",
    "CFResult",
    "CrystalFieldError",
    "OrbitalRow",
    "allowed_parameters",
    "cf_declaration",
    "composition",
    "diffuse",
    "coordination_shell",
    "crystal_field",
    "design_matrix",
    "dominant_center",
    "entropy",
    "entropy_bound_spectrum",
    "fit_crystal_field",
    "geometry",
    "hamiltonian",
    "local_spin",
    "mr_diagnostics",
    "orbital_rows",
    "parse_xyz",
    "point_charge",
    "shell_ranking",
    "single_orbital_entropy_bound",
    "stevens_matrices",
]
