"""G6: the dipole-moment (DM-AS) candidate batch -- inputs for the space scan.

The source protocol (Kaufold et al., 2023; Kaufold & Dong, 2026) selects an
active space by how well it reproduces a reliable dipole moment: run the
ground state of every candidate space ``(n_e, n_o)`` from one cheap orbital
source, compare each space's dipole moment with a reference (experiment or
KS-DFT), and keep the space with the smallest absolute deviation.  The 2026
sequel makes the scan cheap by using **CASCI** (no orbital optimisation) on
**MP2 orbitals**; its own recommendation names CASCI-GDM-AS the sensible
default protocol.

This module builds the input batch for that scan on the ORCA route (measured
on 6.1.1, 2026-09-27):

- one **preparation** run (:file:`prep.inp`) produces the shared orbital source
  -- RHF followed by MP2 with ``NatOrbs`` (the natural orbitals stay in the
  gbw; the CASCI runs read it);
- one CASCI run per candidate (:file:`cand_e{ne}o{no}.inp`):
  ``! <scf> <basis> NoIter moread`` with ``%moinp`` pointing at the prep gbw
  and ``%casscf nel .. norb .. nroots 1 mult .. end``.  Measured: ``NoIter``
  makes ORCA take the occupations from the input (a CAS-CI, "MaxMacroIter 1
  detected"), and a single-root run prints its own dipole moment ("State: 0",
  relaxed density) -- the quantity the protocol needs;
- a run script and a manifest skeleton for menu 20 (the selection).

The candidate family is the source's PASS set: ``n_e`` even, ``6 <= n_e <=
14``, ``n_e/2 + 3 <= n_o <= max_norb``; the PASS+ extras (the 4-electron rows
and the ``n_o = n_e/2 + 2`` rows) can be added -- the source found them
"generally poor performers" but kept them in its scanning set.

The batch runs no calculation itself (this toolkit never does): it writes the
inputs, the script and the manifest; the user runs them and comes back to
menu 20.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field

from ..knowledge.models import EVIDENCE_LITERATURE, EVIDENCE_MEASURED, Evidence

__all__ = [
    "DmBatchError",
    "DmCandidate",
    "DmBatch",
    "pass_candidates",
    "plan_batch",
    "evidence",
]


class DmBatchError(ValueError):
    """The batch cannot be built for the given settings (with a next step)."""


@dataclass(frozen=True)
class DmCandidate:
    """One candidate active space of the scan."""

    nel: int
    norb: int

    def render(self) -> str:
        return f"e{self.nel}o{self.norb}"


#: The source's PASS constraints.
PASS_MIN_ELECTRONS = 6
PASS_MAX_ELECTRONS = 14
PASS_ORBITAL_MARGIN = 3  # n_o >= n_e/2 + 3


def pass_candidates(
    *,
    include_pass_plus: bool = False,
    max_norb: int = 14,
) -> tuple[DmCandidate, ...]:
    """The PASS (or PASS+) candidate set of the source.

    PASS: ``n_e`` even in [6, 14], ``n_e/2 + 3 <= n_o <= max_norb``.  PASS+
    adds the rows the source kept in its scanning set but flagged as generally
    poor: the ``n_e = 4`` rows and the ``n_o = n_e/2 + 2`` rows.
    """
    if max_norb < PASS_MIN_ELECTRONS // 2 + PASS_ORBITAL_MARGIN:
        raise DmBatchError(
            f"max_norb = {max_norb} leaves no PASS candidate (the smallest needs "
            f"n_o = {PASS_MIN_ELECTRONS // 2 + PASS_ORBITAL_MARGIN}). Next step: raise "
            "max_norb."
        )
    candidates = []
    for nel in range(PASS_MIN_ELECTRONS, PASS_MAX_ELECTRONS + 1, 2):
        for norb in range(nel // 2 + PASS_ORBITAL_MARGIN, max_norb + 1):
            candidates.append(DmCandidate(nel, norb))
    if include_pass_plus:
        for norb in range(4, max_norb + 1):
            candidates.append(DmCandidate(4, norb))
        for nel in range(PASS_MIN_ELECTRONS, PASS_MAX_ELECTRONS + 1, 2):
            norb = nel // 2 + 2
            if norb <= max_norb:
                candidates.append(DmCandidate(nel, norb))
    unique = sorted(set(candidates), key=lambda item: (item.nel, item.norb))
    return tuple(unique)


@dataclass(frozen=True)
class DmBatch:
    """The generated batch: file names and their texts, plus the manifest payload."""

    xyz_text: str
    charge: int
    multiplicity: int
    basis: str
    scf: str
    prep: str  # "mp2" | "hf"
    candidates: tuple[DmCandidate, ...]
    prep_name: str
    manifest_name: str
    script_name: str
    manifest: dict = field(default_factory=dict)

    def prep_input(self) -> str:
        lines = [f"! {self.scf} {self.basis}", "%maxcore 1000"]
        if self.prep == "mp2":
            lines += ["%mp2", " NatOrbs true", " Density relaxed", "end"]
        lines += ["* xyz %d %d" % (self.charge, self.multiplicity)]
        lines += self.xyz_text.strip().splitlines()
        lines += ["*"]
        return "\n".join(lines) + "\n"

    def prep_stem(self) -> str:
        return self.prep_name[:-4] if self.prep_name.endswith(".inp") else self.prep_name

    def candidate_input(self, candidate: DmCandidate) -> str:
        lines = [
            f"! {self.scf} {self.basis} NoIter moread",
            f'%moinp "{self.prep_stem()}.gbw"',
            "%casscf",
            f" nel {candidate.nel}",
            f" norb {candidate.norb}",
            " nroots 1",
            f" mult {self.multiplicity}",
            "end",
            "%maxcore 1000",
            "* xyz %d %d" % (self.charge, self.multiplicity),
            *self.xyz_text.strip().splitlines(),
            "*",
        ]
        return "\n".join(lines) + "\n"

    def candidate_name(self, candidate: DmCandidate) -> str:
        return f"cand_{candidate.render()}.inp"

    def render_script(self, orca: str = "orca") -> str:
        lines = [
            "#!/bin/bash",
            "# The DM-AS candidate scan: the preparation first, then every candidate.",
            "set -u",
            f'ORCA="${{{orca.upper()}:={orca}}}"',
            f'"$ORCA" {self.prep_name} > {self.prep_stem()}.out 2>&1 || {{ echo "prep failed"; exit 1; }}',
        ]
        for candidate in self.candidates:
            name = self.candidate_name(candidate)
            lines.append(
                f'"$ORCA" {name} > {name[:-4]}.out 2>&1 || echo "candidate '
                f'{candidate.render()} failed"'
            )
        lines.append('echo "done; run the selection with menu 20 on ' + self.manifest_name + '"')
        return "\n".join(lines) + "\n"


def plan_batch(
    *,
    xyz_text: str,
    charge: int,
    multiplicity: int,
    basis: str,
    scf: str = "RHF",
    prep: str = "mp2",
    include_pass_plus: bool = False,
    max_norb: int = 14,
) -> DmBatch:
    """Plan the candidate batch (texts and manifest; nothing is written here)."""
    if prep not in ("mp2", "hf"):
        raise DmBatchError(
            f"the preparation level {prep!r} is not supported. Next step: use 'mp2' (the "
            "source's choice) or 'hf'."
        )
    if multiplicity != 1:
        raise DmBatchError(
            f"the source protocol is built for singlet ground states (multiplicity 1); "
            f"got {multiplicity}. Next step: use the protocol on closed-shell singlets, "
            "or select the space another way (menus 12-17)."
        )
    if charge != 0:
        raise DmBatchError(
            "the source excludes charged systems: the dipole moment of a charged molecule "
            "depends on the coordinate origin. Next step: use a neutral system, or select "
            "the space another way (menus 12-17)."
        )
    candidates = pass_candidates(include_pass_plus=include_pass_plus, max_norb=max_norb)
    manifest = {
        "candidates": [
            {"nel": candidate.nel, "norb": candidate.norb, "output": f"cand_{candidate.render()}.out"}
            for candidate in candidates
        ],
        "reference": {"output": "REFERENCE_OUTPUT.out"},
        "protocol": "gdm",
    }
    return DmBatch(
        xyz_text=xyz_text,
        charge=charge,
        multiplicity=multiplicity,
        basis=basis,
        scf=scf,
        prep=prep,
        candidates=candidates,
        prep_name="prep.inp",
        manifest_name="manifest.json",
        script_name="run_scan.sh",
        manifest=manifest,
    )


def manifest_text(batch: DmBatch) -> str:
    """The manifest skeleton (menu 20 reads the filled-in version)."""
    return json.dumps(batch.manifest, indent=1) + "\n"


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the protocol and of the ORCA batch route."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "The dipole-moment active-space selection (DM-AS): run the ground state "
                "of every candidate space, compare its dipole moment with a reference "
                "(experiment or KS-DFT), keep the smallest deviation.  The 2026 sequel "
                "replaces CASSCF by CASCI on MP2 orbitals (the scan becomes cheap) and "
                "its recommendations name CASCI-GDM-AS the sensible default, with the "
                "CASCI-EDM-AS (per-state) and D2DM (charge-transfer) variants for the "
                "cases they were designed for.  The candidate family is the source's "
                "PASS set."
            ),
            ref=(
                "Kaufold B. W., Chintala N., Pandeya P., Dong S. S., J. Chem. Theory "
                "Comput., 2023, 19, 2469-2483, DOI 10.1021/acs.jctc.2c01128; Kaufold "
                "B. W., Dong S. S., J. Chem. Theory Comput., 2026, 22, 7624-7641, "
                "DOI 10.1021/acs.jctc.6c00473"
            ),
            url="https://doi.org/10.1021/acs.jctc.6c00473",
            bibkey="kaufold2026casci",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "The ORCA batch route, measured on 6.1.1 (2026-09-27, H2O fixtures): "
                "MP2 with NatOrbs keeps the natural orbitals in the gbw; `!NoIter` with "
                "`moread` runs a CAS-CI (ORCA prints 'MaxMacroIter 1 detected >>> CAS-CI') "
                "taking the input occupations; and every single-root CASSCF/CASCI run "
                "prints its own dipole moment (block 'DIPOLE MOMENT', State: 0, relaxed "
                "density).  A state-averaged run prints only the state-average dipole "
                "(State: -1) -- so the scan uses nroots 1, and the selection refuses an "
                "average-only output."
            ),
            ref="tests/test_dm_batch.py, tests/test_dm_selection.py; fixtures/orca/h2o_dm_*",
        ),
    )
