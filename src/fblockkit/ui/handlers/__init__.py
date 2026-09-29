"""Menu handler functions: every handler takes the session and produces text and files.

Discipline: the user's input files are only ever read; every product is a new file
(``*.fbk.md`` / ``*.fbk.json`` / ``*.fbk.inp``) and never overwrites the original; the
output contains no unstable content such as timestamps (so a script replays
byte-identically).

The package is grouped by menu family, one module per theme; the registry at the
bottom of this file is the only thing the CLI imports:

- ``core``: the check-up report, the geometry report, script saving, quitting;
- ``inputs``: the input-generation side (menu 3 and the basis/tool/cross-level
  helpers, the SCF rescue);
- ``crystal``: the crystal-field fit and the point-charge estimate;
- ``orbitals``: the export-driven analyses (exact entropy, orbital space, AVAS,
  portrait, magnetic doublets, cross-structure mapping, WASP);
- ``selection``: the active-space selection family (DM-AS batch/select, APC,
  ASS1ST, QICAS, AEGISS, TNASS);
- ``external``: the other-program pair (write a pysisyphus run input, read the
  run back);
- ``common``: the report-file plumbing and the shared CASSCF-reference reader.
"""


from __future__ import annotations

from .core import (
    geometry_report,
    quit_session,
    report_output,
    save_script,
)
from .crystal import (
    crystal_field_fit,
    point_charge_estimate,
)
from .external import (
    pysisyphus_generate,
    pysisyphus_report,
)
from .inputs import (
    basis_query,
    cross_level,
    deltascf_generate,
    generate_input,
    imag_disp_generate,
    ras_ormas_generate,
    scf_rescue,
    tool_guide,
    tool_search,
)
from .orbitals import (
    ailft_report,
    avas_target,
    exact_entropy,
    judd_ofelt_fit,
    magnetic_doublets_report,
    orbital_mapping_report,
    orbital_portrait_report,
    orbital_space,
    perturb_batch,
    pnmr_report,
    relaxation_report,
    state_data_report,
    wasp_guess,
    xas_report,
)
from .selection import (
    aegiss_select,
    apc_ranking,
    ass1st_round,
    ass1st_start,
    dm_batch_generate,
    dm_select,
    qicas_optimize,
    tnass_select,
)

__all__ = ["HANDLERS"]

HANDLERS = {
    "report_output": report_output,
    "geometry_report": geometry_report,
    "generate_input": generate_input,
    "deltascf_generate": deltascf_generate,
    "ras_ormas_generate": ras_ormas_generate,
    "imag_disp_generate": imag_disp_generate,
    "basis_query": basis_query,
    "tool_search": tool_search,
    "tool_guide": tool_guide,
    "cross_level": cross_level,
    "save_script": save_script,
    "scf_rescue": scf_rescue,
    "crystal_field_fit": crystal_field_fit,
    "point_charge_estimate": point_charge_estimate,
    "exact_entropy": exact_entropy,
    "orbital_space": orbital_space,
    "avas_target": avas_target,
    "orbital_portrait": orbital_portrait_report,
    "magnetic_doublets": magnetic_doublets_report,
    "orbital_mapping": orbital_mapping_report,
    "wasp_guess": wasp_guess,
    "perturb_batch": perturb_batch,
    "state_data_report": state_data_report,
    "judd_ofelt_fit": judd_ofelt_fit,
    "pnmr_report": pnmr_report,
    "relaxation_report": relaxation_report,
    "ailft_report": ailft_report,
    "xas_report": xas_report,
    "pysisyphus_generate": pysisyphus_generate,
    "pysisyphus_report": pysisyphus_report,
    "dm_batch": dm_batch_generate,
    "dm_select": dm_select,
    "apc_ranking": apc_ranking,
    "ass1st_start": ass1st_start,
    "ass1st_round": ass1st_round,
    "qicas_optimize": qicas_optimize,
    "aegiss_select": aegiss_select,
    "tnass_select": tnass_select,
    "quit": quit_session,
}
