"""Prepare the manual's example inputs under work/ (deterministic; no fixtures copied here).

Run through run-examples.sh, which also copies the fixtures the sessions read.
"""

from __future__ import annotations

import json
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)

sys.path.insert(0, str(HERE.parents[2] / "src"))

from fblockkit.analysis import crystal_field as cf  # noqa: E402


def write(name: str, text: str) -> None:
    (WORK / name).write_text(text, encoding="utf-8")
    print(f"wrote work/{name}")


# --- structures --------------------------------------------------------------

write(
    "octahedron.xyz",
    "7\nCeO6 octahedron (manual example)\n"
    "Ce  0.000  0.000  0.000\n"
    "O   2.400  0.000  0.000\n"
    "O  -2.400  0.000  0.000\n"
    "O   0.000  2.400  0.000\n"
    "O   0.000 -2.400  0.000\n"
    "O   0.000  0.000  2.400\n"
    "O   0.000  0.000 -2.400\n",
)
write("ce_atom.xyz", "1\nCe atom (manual example)\nCe 0.0 0.0 0.0\n")
write(
    "ceo6.xyz",
    "7\nCeO6 octahedral toy structure (manual example)\n"
    "Ce  0.000  0.000  0.000\n"
    "O   2.400  0.000  0.000\n"
    "O  -2.400  0.000  0.000\n"
    "O   0.000  2.400  0.000\n"
    "O   0.000 -2.400  0.000\n"
    "O   0.000  0.000  2.400\n"
    "O   0.000  0.000 -2.400\n",
)

# --- crystal-field fitting input (menu 10), with A5 declaration + compare ----

TRUE = {(2, 0): 1000.0, (4, 0): 10.0, (4, 3): 50.0, (6, 0): 1.0}
J = 8.0
import numpy as np  # noqa: E402

hamiltonian = cf.hamiltonian(TRUE, J)

# Sample oriented (non-eigen-) states, the way the reference method does: random
# states in the |J M> basis, no eigendecomposition involved (the example input is
# then identical on every platform: eigenvector phases and BLAS rounding would
# otherwise leak into the fit).  Sampling eigenstates alone would be rank
# deficient anyway -- the members of a time-reversal pair share every expectation
# value -- so oriented states are required.
rng = random.Random(20260926)
levels, coefficients = [], []
for _ in range(24):
    weights = np.array([rng.gauss(0.0, 1.0) + 1j * rng.gauss(0.0, 1.0) for _ in range(17)])
    state = weights / np.linalg.norm(weights)
    # Quantise the stored coefficients to 1e-7 (well inside the fitter's 1e-6 norm
    # tolerance) so the committed example input is byte-stable across platforms.
    stored = np.array(
        [complex(round(v.real, 7), round(v.imag, 7)) for v in state]
    )
    # The fitter's model is Eq. (3) as written for the coefficients it is given
    # (norms within its 1e-6 tolerance of 1 are accepted without rescaling), so
    # the level recorded here is the expectation value of the STORED vector --
    # same normalisation convention, no artificial residual.
    levels.append(round(float((stored.conj() @ hamiltonian @ stored).real), 5))
    coefficients.append([[float(v.real), float(v.imag)] for v in stored])

payload = {
    "point_group": "C3",
    "J": J,
    "levels": [round(v, 10) for v in levels],
    "coefficients": coefficients,
    "declaration": {
        "convention": "Stevens (Rudowicz/Ryabov lineage, cosine/sine)",
        "projection": "J = 8 manifold",
        "units": "cm^-1",
        "z_axis": "the three-fold axis of the model structure",
        "origin": "Ce site",
    },
    "compare": {
        "label": "Chilton Table S1, SINGLE_ANISO (L = 5) column",
        "parameters": {"2,0": 1879.0, "4,0": 125.0, "4,3": -1.0, "6,0": 21.0},
        "declaration": {"projection": "L = 5 manifold"},
    },
}
(WORK / "cf_c3.json").write_text(json.dumps(payload, indent=1) + "\n", encoding="utf-8")
print("wrote work/cf_c3.json")

# --- the magnetic-doublet table (menu 16): a slice of the source's calibration
# table (fixtures/literature/chilton_s2_dy19.json), prepared the way a user would
# from a published SINGLE_ANISO table -- the two doublets of [Dy(BC4Ph5)2]- with
# theta3 as the source gives it (measured against the ground doublet, which the
# published table does not list).

SOURCE = HERE.parents[2] / "fixtures" / "literature" / "chilton_s2_dy19.json"
calibration = json.loads(SOURCE.read_text(encoding="utf-8"))
row = next(item for item in calibration["rows"] if item["compound"] == "[Dy(BC4Ph5)2]-")
doublets = []
for role, label in (
    ("safe", "doublet that supports excitation"),
    ("qtm", "doublet that facilitates QTM"),
):
    entry = row[role]
    doublets.append(
        {
            "label": label,
            "g": [entry["g1"], entry["g2"], entry["g3"]],
            "theta3": entry["theta3_deg"],
        }
    )
doublet_payload = {
    "system": row["compound"] + " (the calibration table's row 7)",
    "doublets": doublets,
}
(WORK / "dy_doublets.json").write_text(
    json.dumps(doublet_payload, indent=1) + "\n", encoding="utf-8"
)
print("wrote work/dy_doublets.json")

# --- the cross-structure mapping manifest (menu 17): the frozen N2 scan's three
# localized exports, with a selection on one bond orbital of the first structure
# -- the run reads back the consistent space (the degenerate bond triad) that
# selection implies in every structure.

mapping_payload = {
    "structures": [
        {"name": f"r={r}", "export": f"n2_scan_{r}.loc.json"}
        for r in ("1.094", "1.600", "2.600")
    ],
    "selections": {"r=1.094": [4]},
}
(WORK / "mapping_manifest.json").write_text(
    json.dumps(mapping_payload, indent=1) + "\n", encoding="utf-8"
)
print("wrote work/mapping_manifest.json")
