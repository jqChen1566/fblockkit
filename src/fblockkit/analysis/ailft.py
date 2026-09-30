"""ab initio ligand-field analysis from an ORCA AILFT output (menu 38).

The AILFT module fits the parameters of a ligand-field model (a 5x5 or 7x7
one-electron matrix plus Slater-Condon/electron-repulsion parameters and the
SOC constant) to an ab initio effective Hamiltonian built from SA-CASSCF
(and NEVPT2) wave functions (manual section 3.13.16; the original articles
are Jung, Atanasov & Neese, Inorg. Chem. 2017 -- the actinide/lanthanide
reference -- and Lang, Atanasov & Neese, J. Phys. Chem. A 2020).  This menu
reads what the engine already fitted; it never refits:

- per theoretical level (CASSCF, NEVPT2): the ligand-field one-electron
  eigenfunctions (LF splitting), the Slater-Condon parameters, the Racah
  parameters (B, C, C/B) and the ``*.lft.gbw`` file name;
- the fit quality: per-block and total RMS errors, Pearson's correlation;
  the report repeats the source's statement that the near-zero CASSCF-level
  RMS is intrinsic (the model parametrizes that level exactly), while the
  correlation-level RMS reflects, among other things, the neglected
  anisotropy of the electron-electron repulsion in covalent complexes
  (Lang et al. 2020, the rmsd discussion);
- the SOC constant (ZETA_D/ZETA_F) from the SOC section (note: SOC is
  treated with CASSCF orbitals).

When the caller passes free-ion references, the report adds the
nephelauxetic ratios beta = B(cplx)/B(free ion) and the relativistic
nephelauxetic zeta/zeta0 -- reductions are the classic covalency
indicators (Jung et al. 2017).
"""

from __future__ import annotations

from typing import Any

from ..knowledge.models import EVIDENCE_LITERATURE, EVIDENCE_MEASURED, Evidence

__all__ = ["AilftError", "lf_splitting", "render", "evidence"]


class AilftError(Exception):
    """A refusal with the next step spelled out."""


def lf_splitting(level: dict[str, Any]) -> dict[str, Any] | None:
    """The ligand-field eigenfunction energies (cm^-1) and their spread."""
    eigen = level.get("eigenfunctions")
    if not eigen:
        return None
    energies = [entry["cm1"] for entry in eigen]
    return {
        "energies_cm1": energies,
        "spread_cm1": max(energies) - min(energies),
    }


def _racah_value(level: dict[str, Any], label: str) -> float | None:
    entry = level.get("racah", {}).get(label)
    return entry["cm1"] if isinstance(entry, dict) else None


def render(
    data: dict[str, Any],
    *,
    source: str,
    free_ion_B_cm1: float | None = None,
    free_ion_zeta_cm1: float | None = None,
) -> str:
    """The menu-38 report body."""
    if not data.get("present"):
        raise AilftError(
            "no AILFT section found. Next step: this menu reads an ORCA CASSCF "
            "output with the AILFT driver requested (%casscf with ActOrbs "
            "dOrbs/fOrbs or LFTCase); check the run reached the 'AB INITIO "
            "LIGAND FIELD THEORY' section."
        )
    lines: list[str] = []
    lines.append(f"Ab initio ligand-field analysis ({source})")
    header = data.get("header", {})
    if header:
        mo = header.get("mo_range")
        lines.append(
            f"  {header.get('configuration', '?')} configuration, "
            f"{header.get('ci_blocks', '?')} CI blocks, "
            + (f"MOs {mo[0]} to {mo[1]}, " if mo else "")
            + f"metal center atom {header.get('center_atom', '?')}"
        )
    for name, level in data.get("levels", {}).items():
        lines.append("")
        lines.append(f"  == {name.upper()} level ==")
        split = lf_splitting(level)
        if split is not None:
            energies = " ".join(f"{value:.1f}" for value in split["energies_cm1"])
            lines.append(
                f"    ligand-field eigenfunctions (cm-1): {energies}"
                f"   (spread {split['spread_cm1']:.1f} cm-1)"
            )
        scp = level.get("slater_condon") or []
        if scp:
            lines.append(
                "    Slater-Condon (cm-1): "
                + "  ".join(
                    f"{entry['label']} = {entry['cm1']:.1f}"
                    + (" (fixed)" if entry["fixed"] else "")
                    for entry in scp
                )
            )
        b = _racah_value(level, "B")
        c = _racah_value(level, "C")
        cb = level.get("racah", {}).get("C_over_B")
        parts = []
        if b is not None:
            parts.append(f"B = {b:.1f}")
        if c is not None:
            parts.append(f"C = {c:.1f}")
        if cb is not None:
            parts.append(f"C/B = {cb:.3f}")
        if parts:
            lines.append("    Racah (cm-1): " + "  ".join(parts))
            if b is not None and free_ion_B_cm1:
                lines.append(
                    f"    nephelauxetic ratio beta = B/B0 = {b / free_ion_B_cm1:.3f}"
                    f" (B0 = {free_ion_B_cm1:.1f} cm-1, caller's reference)"
                )
        if level.get("lft_gbw"):
            lines.append(f"    LFT orbitals stored in {level['lft_gbw']}")
    for name, stats in data.get("fit", {}).items():
        if stats.get("total_rms_cm1") is None:
            continue
        blocks = stats.get("blocks") or []
        block_text = "  ".join(
            f"block {block['index']} {block['rms_cm1']:.1f}" for block in blocks if block["rms_cm1"] is not None
        )
        lines.append("")
        lines.append(
            f"  fit quality ({name.upper()}): total RMS = {stats['total_rms_cm1']:.1f} cm-1"
            + (f"  [{block_text}]" if block_text else "")
            + (f"  Pearson = {stats['pearson']:.3f}" if stats.get("pearson") is not None else "")
        )
    soc = data.get("soc")
    if soc:
        summary = soc.get("summary") or {}
        for label, value in summary.items():
            lines.append(f"  SOC constant: {label} = {value:.2f} cm-1 (SOC based on {'/'.join(soc['bases'])} orbitals)")
            if free_ion_zeta_cm1:
                lines.append(
                    f"    relativistic nephelauxetic ratio = {value / free_ion_zeta_cm1:.3f}"
                    f" (zeta0 = {free_ion_zeta_cm1:.1f} cm-1, caller's reference)"
                )
    lines.append("")
    lines.append(
        "  Reading notes: the near-zero CASSCF-level RMS is intrinsic -- the LFT"
        " parametrization is exact for that level; the correlation-level RMS reflects"
        " (among other things) the neglected anisotropy of electron-electron repulsion"
        " in covalent complexes (Lang, Atanasov & Neese, J. Phys. Chem. A 2020)."
        " The parameters are model quantities of the ligand-field Hamiltonian; this"
        " menu reads the engine's fit, it never refits. Nephelauxetic comparisons"
        " (Jung, Atanasov & Neese, Inorg. Chem. 2017) need the caller's free-ion"
        " references -- pass them to get the ratios."
    )
    return "\n".join(lines)


def evidence() -> tuple[Evidence, ...]:
    """Provenance of the analysis notes and the measured format facts."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "The actinide/lanthanide AILFT reference (Slater-Condon and SOC"
                " parameter analysis; nephelauxetic and relativistic-nephelauxetic"
                " reductions as covalency indicators) is Jung, Atanasov & Neese 2017;"
                " the rmsd semantics (exact CASSCF-level parametrization; the"
                " electron-repulsion-anisotropy origin of the correlation-level"
                " rmsd) follows Lang, Atanasov & Neese 2020 (CC-BY)."
            ),
            ref=(
                "Jung J., Atanasov M., Neese F., Inorg. Chem. 2017, 56, 8802-8816, "
                "DOI 10.1021/acs.inorgchem.7b00642; Lang L., Atanasov M., Neese F., "
                "J. Phys. Chem. A 2020, 124, 1025-1037, DOI 10.1021/acs.jpca.9b11227"
            ),
            url="https://doi.org/10.1021/acs.jpca.9b11227",
            bibkey="lang2020ailft",
        ),
        Evidence(
            kind=EVIDENCE_MEASURED,
            text=(
                "ORCA 6.1.1 fixtures ailft/ (2026-09-29): the AILFT section prints one"
                " parameter block per level (VLFT matrix, Slater-Condon, Racah with"
                " C/B, ligand-field eigenfunctions, *.lft.gbw), the AI-vs-LF fit"
                " comparison per level (block tables, RMS, Pearson; the confidence-"
                " interval rows reuse the SCF line shape) and the SOC section with"
                " the ZETA_D/ZETA_F summary; for the d8 free ion the CASSCF-level RMS"
                " is 0.0 and the NEVPT2-level 457.4 cm-1."
            ),
            ref="ORCA 6.1.1, measured 2026-09-29 (fixtures ailft/)",
        ),
    )
