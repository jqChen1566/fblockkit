"""Recipe layer: system profile -> recommendation -> input rendering
(architecture design v0.1 §3 L2).

Public interface:

- ``recommend(profile)`` / ``recommend_with_advice(profile)``: rule matching + G1 basis advice;
- ``recommend_basis_ecp(profile)``: basis-set/ECP matching only (structured output);
- ``plan_convergence(...)``: G4 convergence/initial-guess setting fragments;
- ``render_orca_input(...)`` / ``run_guidance(...)``: generate an ORCA input and the
  written run guidance;
- ``render_openmolcas_input(...)`` / ``run_guidance_openmolcas(...)``: the G3
  magnetic-property chain (SA-CASSCF / RASSI / SINGLE_ANISO) input and run guidance.

This layer never changes the user's files: its products are always new strings / new
file contents, and the caller decides whether to write them to disk.
"""

from .active_space import (
    ActiveSpaceError,
    ActiveSpaceSuggestion,
    f_electron_count,
    suggest_active_space,
    verification_evidence,
    verification_protocol,
)
from .basis_ecp import (
    BasisAdvice,
    BasisDataError,
    BasisEntry,
    element_z,
    is_f_element,
    load_basis_entries,
    recommend_basis_ecp,
)
from .convergence import ConvergencePlan, plan_convergence
from .dmet import DmetError, DmetPlan, plan_dmet, render as render_dmet
from .openmolcas import (
    Deviation,
    OpenMolcasChainSpec,
    OpenMolcasSpecError,
    max_csfs,
    own_deviations,
    render_openmolcas_input,
    run_guidance_openmolcas,
    spin_orbit_state_count,
)
from .openmolcas import evidence as openmolcas_evidence
from .recommend import facts_from_profile, recommend, recommend_with_advice
from .render import RenderError, render_orca_input, run_guidance

__all__ = [
    "ActiveSpaceError",
    "ActiveSpaceSuggestion",
    "BasisAdvice",
    "BasisDataError",
    "BasisEntry",
    "ConvergencePlan",
    "DmetError",
    "DmetPlan",
    "plan_dmet",
    "render_dmet",
    "Deviation",
    "OpenMolcasChainSpec",
    "OpenMolcasSpecError",
    "RenderError",
    "element_z",
    "f_electron_count",
    "facts_from_profile",
    "is_f_element",
    "load_basis_entries",
    "max_csfs",
    "openmolcas_evidence",
    "own_deviations",
    "plan_convergence",
    "recommend",
    "recommend_basis_ecp",
    "recommend_with_advice",
    "render_openmolcas_input",
    "render_orca_input",
    "run_guidance",
    "run_guidance_openmolcas",
    "spin_orbit_state_count",
    "suggest_active_space",
    "verification_evidence",
    "verification_protocol",
]
