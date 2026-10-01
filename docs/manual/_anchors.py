"""List the anchors of every session capture: menu-table headers and 'Choose > N' lines.

Used to re-anchor the historically mis-shifted references (batch 2/3).
"""
from pathlib import Path

BASE = Path(__file__).resolve().parent / "examples" / "expected"

FILES = [
    "m2_geometry.txt", "m3_generate.txt", "m4_basis.txt", "m5_tools.txt",
    "m7_crosslevel.txt", "m9_scf.txt", "m10_crystal_field.txt",
    "m11_point_charge.txt", "m12_exact_entropy.txt", "m13_orbital_space.txt",
    "m14_avas.txt", "m15_portrait.txt", "m17_mapping.txt", "m18_guess.txt",
    "m19_dm_batch.txt", "m20_dm_select.txt", "m21_apc.txt", "m22_ass1st_start.txt",
    "m23_ass1st_round.txt", "m24_qicas.txt", "m25_aegiss.txt", "m26_tnass.txt",
    "m27_deltascf.txt", "m28_rasormas.txt", "m29_perturb.txt", "m30_imagdisp.txt",
    "m31_statedata.txt", "m32_pysisyphus.txt", "m33_pysisyphus.txt",
    "m34_juddofelt.txt", "m35_pnmr.txt", "m36_relaxation.txt", "m37_xas.txt",
    "m38_ailft.txt", "m40_hyperfine.txt", "m41_magnetocaloric.txt",
]
for name in FILES:
    lines = (BASE / name).read_text(encoding="utf-8").splitlines()
    marks = []
    for i, line in enumerate(lines, start=1):
        if line.startswith("Available functions:"):
            marks.append(f"TBL@{i}")
        elif line.startswith("Choose a number"):
            marks.append(f"CHOOSE@{i}:{line.split('> ')[-1]}")
    print(f"{name:28s} {' '.join(marks)}")
