"""G3: OpenMolcas magnetic-property chain -- input generation (SA-CASSCF / RASSI / SINGLE_ANISO).

This module writes text only: fBlockKit never runs an engine. The rendered file is a
complete OpenMolcas job for a mononuclear f-block magnetic centre, following the annotated
input skeleton printed in the Supporting Information of the practical guide:

    N. F. Chilton, "Ab initio electronic structure calculations of lanthanide
    single-molecule magnets; a practical guide", Chem. Soc. Rev. 2025, 54, 11468-11487,
    DOI 10.1039/d5cs00493d -- Supporting Information, section 4.1 (annotated OpenMolcas
    input for [Dy(OH)Br]+, 28 line annotations).

The internal close-reading of that SI (record: ``文献细读/细读_Chilton指南_SI.md``,
2026-09-25) transcribes the input, its annotations and the protocol below; the DOI itself
comes from this project's survey metadata, not from the SI, which does not print it.

What comes from the source (bibkey ``chilton2025abinitio``)
-----------------------------------------------------------
- Module chain ``&GATEWAY`` -> ``&SEWARD`` -> 3 x ``&RASSCF`` -> ``&RASSI`` -> ``&SINGLE_ANISO``,
  with the per-line annotations of the source placed on the lines they describe
  (that placement is an inference from the source's printed order -- see ``own_deviations``).
- ``CiRoot= N N 1`` (root count given twice, unity weighting) and the ``Typeindex`` /
  ``FILEORB`` / ``>> COPY`` chain that hands each spin block the previous block's orbitals.
- ``Nactel= <n> 0 0`` with ``RAS2= <nactorb>``, i.e. a CAS with no RAS1/RAS3 split; the
  source's Dy(III) case is 9 electrons in 7 orbitals (the 4f shell).
- Three-tier ANO-RCC basis assignment: magnetic centre at ANO-RCC-VTZP, first coordination
  sphere at ANO-RCC-VDZP, all remaining atoms at ANO-RCC-VDZ (source section 4.2, item 5).
- ``&RASSI`` keywords ``SPIN`` / ``MEES`` / ``EPRG= 7.0D-1`` / ``Nr of JobIph= n ALL`` /
  ``IPHN`` / ``PROP`` with ``'ANGMOM' i``, and the ``&SINGLE_ANISO`` block with ``CRYS`` /
  ``MLTP`` / ``TINT`` / ``HINT`` / ``TMAG`` / ``QUAX`` (source section 4.1).

Decisions taken here that the source does not show (each is tagged ``[fBlockKit:...]``
inside the rendered file and listed by :func:`own_deviations`): the comment character
(verified against the OpenMolcas documentation on 2026-09-25 -- ``users.guide/emil.html``:
"a comment line identified by an asterisk (*) in the first position on the line", and
comments may not sit between a keyword and its continuation lines, which this template
respects), the ``>> COPY`` of the RASSCF orbital file, the ``group=`` default, the
ordering of the spin blocks, the fixed ``&SINGLE_ANISO`` sampling parameters, the
generalisation of ``'ANGMOM'`` to N blocks, the annotation placement, and the verbatim
(not cross-checked) spelling of the OpenMolcas keywords. This module deliberately does not
invent OpenMolcas syntax beyond the source: where a choice had to be made it is marked in
the file and listed by :func:`own_deviations` instead of passing silently. Apart from the
comment marker, nothing here is verified against the OpenMolcas documentation, which was
not consulted for this version.

Validation
----------
:class:`OpenMolcasChainSpec` validates on construction (element is an f-block element,
active space can host the electrons, each requested spin block is allowed by spin parity
and fits in the CI space of CAS(nactel, nactorb), every free-text field is single-line
ASCII because it is pasted into the input file). The CI-space bound uses the Weyl-Palmer
binomial count, cross-checked against the source's own numbers: 21 sextets, 224 quartets
and 490 doublets for 4f^9 in 7 orbitals.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, cast

from jinja2 import Environment, FileSystemLoader, StrictUndefined

from ..knowledge.elements import ElementError, element_z, is_f_element
from ..knowledge.models import EVIDENCE_LITERATURE, Evidence
from .render import RenderError

TEMPLATES_DIR = Path(__file__).resolve().parent.parent / "knowledge" / "templates"
TEMPLATE_NAME = "openmolcas_magnetic.j2"

SI_BIBKEY = "chilton2025abinitio"
SI_DOI_URL = "https://doi.org/10.1039/d5cs00493d"
SI_SKELETON_REF = (
    "Supporting Information, section 4.1 (annotated OpenMolcas input for [Dy(OH)Br]+, "
    "28 line annotations)"
)
SI_PROTOCOL_REF = (
    "Supporting Information, section 4.2 (seven-point protocol for 19 mononuclear "
    "Dy(III) single-molecule magnets)"
)
SI_CRITERION_REF = (
    "Supporting Information, sections 3.1 and 4.3 (average transverse g value, "
    "g_T * theta_3 = 20 criterion)"
)
SI_SOURCE_NOTE = (
    "Chilton, Chem. Soc. Rev. 2025, 54, 11468-11487, Supporting Information section 4.1 "
    "([Dy(OH)Br]+ example, 28 line annotations)"
)
# Default basis tiers (source section 4.2, item 5).
DEFAULT_BASIS_METAL = "ANO-RCC-VTZP"
DEFAULT_BASIS_LIGAND = "ANO-RCC-VDZP"
DEFAULT_BASIS_OUTER = "ANO-RCC-VDZ"

# &SINGLE_ANISO sampling parameters, reproduced from the source example: pseudospin
# multiplets to treat, temperature grid for the susceptibility, field grid and
# temperatures for the magnetisation, quantisation axis from the ground doublet.
SINGLE_ANISO_MLTP = "1; 2"
SINGLE_ANISO_TINT = "0.0 330.0 330 0.0001"
SINGLE_ANISO_HINT = "0.0 10.0 201"
SINGLE_ANISO_TMAG = "6 1.8 2 4 5 10 20"
SINGLE_ANISO_QUAX = "1"
RASSI_EPRG = "7.0D-1"

# Source example numbers, used in messages and as the documented cross-check.
SI_EXAMPLE_SPIN_BLOCKS = ((6, 21), (4, 48), (2, 32))
SI_DATASET_CSFS = "21 sextets, 224 quartets and 490 doublets"

_POINT_GROUP = re.compile(r"^[A-Za-z][A-Za-z0-9*]*$")


class OpenMolcasSpecError(RenderError):
    """Invalid OpenMolcas chain specification (raised while building or rendering it).

    Subclasses :class:`fblockkit.recipe.render.RenderError` so callers that already catch
    render errors cover this renderer too; every message ends with a "Next step: ..."
    sentence.
    """


@dataclass(frozen=True)
class Deviation:
    """A choice made by this toolkit that the source input does not show."""

    tag: str  # appears in the rendered file as "[fBlockKit:<tag>]"
    text: str  # what was chosen, and why the source does not cover it


@dataclass(frozen=True)
class OpenMolcasChainSpec:
    """Everything needed to render one OpenMolcas magnetic-chain job.

    ``element``: f-block symbol of the magnetic centre (La-Lu, Ac-Lr).
    ``charge``: molecular charge of the isolated ion (the ``&RASSCF`` ``Charge=``).
    ``nactel`` / ``nactorb``: active electrons and active orbitals; the source's Dy(III)
        case is 9 in 7 (the 4f shell), written there as ``Nactel= n 0 0`` and ``RAS2= n``.
    ``roots``: spin blocks as ``{spin multiplicity (2S+1): number of CI roots}``, e.g.
        ``{6: 21, 4: 48, 2: 32}`` for the source example. Stored (and rendered) as a
        tuple of ``(multiplicity, nroots)`` pairs sorted by decreasing multiplicity.
    ``coord_file``: the structure-and-basis file given to ``Coord=`` -- the source keeps
        the geometry and the per-atom basis labels in one file.

    The three basis tiers are labels that must match the ones inside ``coord_file``; they
    are printed into the input as a reminder and are not written into that file.
    """

    element: str
    charge: int
    nactel: int
    nactorb: int
    roots: Mapping[int, int] | tuple[tuple[int, int], ...]
    coord_file: str
    basis_metal: str = DEFAULT_BASIS_METAL
    basis_ligand: str = DEFAULT_BASIS_LIGAND
    basis_outer: str = DEFAULT_BASIS_OUTER
    symmetry: str = "C1"
    title: str = ""
    angmom_origin: tuple[float, float, float] = (0.0, 0.0, 0.0)

    def __post_init__(self) -> None:
        object.__setattr__(self, "element", _check_element(self.element))
        object.__setattr__(self, "charge", _check_charge(self.charge))
        nactel = _check_count(self.nactel, "nactel", "electrons", 9)
        nactorb = _check_count(self.nactorb, "nactorb", "orbitals", 7)
        object.__setattr__(self, "nactel", nactel)
        object.__setattr__(self, "nactorb", nactorb)
        if nactel > 2 * nactorb:
            raise _error(
                f"nactel={nactel} electrons cannot be hosted by nactorb={nactorb} "
                f"orbitals ({2 * nactorb} spin orbitals at most).",
                "enlarge the active space or lower the electron count.",
            )
        if nactorb < 7:
            raise _error(
                f"nactorb={nactorb} is smaller than the seven f orbitals of a "
                f"{self.element} centre, so the f shell would not fit in the active space.",
                "give nactorb >= 7 (the source's 4f protocol uses 9 electrons in 7 "
                "orbitals; add the d shell, e.g. 12 orbitals, if you need it).",
            )
        object.__setattr__(self, "roots", _check_roots(self.roots, nactel, nactorb))
        object.__setattr__(self, "coord_file", _check_one_line(self.coord_file, "coord_file"))
        for name in ("basis_metal", "basis_ligand", "basis_outer"):
            object.__setattr__(self, name, _check_one_line(getattr(self, name), name))
        symmetry = _check_one_line(self.symmetry, "symmetry")
        if not _POINT_GROUP.match(symmetry):
            raise _error(
                f"symmetry={symmetry!r} is not a point-group token (letters, digits and "
                "'*' only, e.g. C1, Ci, C2v, D2h).",
                "give the point group OpenMolcas should use (C1 = no symmetry).",
            )
        object.__setattr__(self, "symmetry", symmetry)
        object.__setattr__(
            self, "title", _check_one_line(self.title, "title") if self.title else ""
        )
        object.__setattr__(self, "angmom_origin", _check_origin(self.angmom_origin))

    @property
    def spin_blocks(self) -> tuple[tuple[int, int], ...]:
        """``((multiplicity, nroots), ...)``, highest spin first."""
        return cast("tuple[tuple[int, int], ...]", self.roots)


# --- validation helpers --------------------------------------------------


def _error(message: str, action: str) -> OpenMolcasSpecError:
    return OpenMolcasSpecError(f"{message} Next step: {action}")


def _check_element(value: object) -> str:
    if not isinstance(value, str) or not value.strip():
        raise _error(
            f"element must be a non-empty element symbol, got {value!r}.",
            "give the f-block element symbol of the magnetic centre, e.g. 'Dy'.",
        )
    try:
        if not is_f_element(value):
            raise _error(
                f"element {value.strip().capitalize()} (Z={element_z(value)}) is not a "
                "lanthanide or actinide, and this chain targets 4f/5f magnetic centres.",
                "use an f-block element of the magnetic centre, or the ORCA recipes of "
                "this package for other elements.",
            )
    except ElementError as exc:
        # The inner message is not repeated here: it is Chinese, and this module's
        # errors are English-only.
        raise _error(
            f"element {value!r} is not a known element symbol.",
            "check the spelling of the element symbol, e.g. 'Dy', 'Ho', 'U'.",
        ) from exc
    return value.strip().capitalize()


def _check_charge(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise _error(
            f"charge must be an integer, got {value!r} of type {type(value).__name__}.",
            "give the molecular charge of the isolated ion as an int, e.g. charge=1.",
        )
    if abs(value) > 9:
        raise _error(
            f"charge={value} is outside the sanity range -9..9 used by this generator.",
            "check the molecular charge of the isolated ion, which the source treats "
            "without counter-ions.",
        )
    return value


def _check_count(value: object, name: str, unit: str, typical: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise _error(
            f"{name} must be an integer >= 1, got {value!r}.",
            f"give the number of active {unit} (the source's Dy(III) 4f case uses "
            f"{typical}).",
        )
    return value


def _check_roots(value: object, nactel: int, nactorb: int) -> tuple[tuple[int, int], ...]:
    if isinstance(value, Mapping):
        items: list[tuple[object, object]] = list(value.items())
    elif isinstance(value, (tuple, list)):
        items = []
        for item in value:
            if not isinstance(item, (tuple, list)) or len(item) != 2:
                raise _error(
                    f"roots entry {item!r} is not a (spin multiplicity, nroots) pair.",
                    "give roots as {multiplicity: nroots}, e.g. {6: 21, 4: 48, 2: 32}.",
                )
            items.append((item[0], item[1]))
    else:
        raise _error(
            "roots must be a mapping {spin multiplicity: number of CI roots}, got "
            f"{type(value).__name__}.",
            "give one entry per spin block, e.g. {6: 21, 4: 48, 2: 32} (sextets, "
            "quartets and doublets as in the source example).",
        )
    if not items:
        raise _error(
            "roots is empty, so the job would contain no &RASSCF block.",
            "give at least one spin block, e.g. {6: 21, 4: 48, 2: 32}.",
        )
    reachable = [
        multiplicity
        for multiplicity in range(1, 2 * nactorb + 2)
        if max_csfs(nactel, nactorb, multiplicity) > 0
    ]
    pairs: list[tuple[int, int]] = []
    seen: set[int] = set()
    for raw_spin, raw_nroots in items:
        if isinstance(raw_spin, bool) or not isinstance(raw_spin, int) or raw_spin < 1:
            raise _error(
                f"roots key {raw_spin!r} is not a spin multiplicity.",
                "use the spin multiplicity 2S+1 (6 for a sextet, 4 for a quartet, "
                "2 for a doublet).",
            )
        if raw_spin in seen:
            raise _error(
                f"spin multiplicity {raw_spin} appears twice in roots.",
                "give each spin multiplicity once; a repeated block would overwrite the "
                "files of the first one.",
            )
        seen.add(raw_spin)
        if isinstance(raw_nroots, bool) or not isinstance(raw_nroots, int) or raw_nroots < 1:
            raise _error(
                f"roots[{raw_spin}]={raw_nroots!r} is not a positive number of CI roots.",
                "give the number of roots as an integer >= 1, e.g. {6: 21}.",
            )
        if raw_spin not in reachable:
            raise _error(
                f"spin multiplicity {raw_spin} is not reachable by {nactel} electrons in "
                f"{nactorb} orbitals (multiplicity and electron count must have opposite "
                "parity, and the spin cannot exceed the unpaired electrons).",
                f"pick one of {reachable} for this active space.",
            )
        available = max_csfs(nactel, nactorb, raw_spin)
        if raw_nroots > available:
            raise _error(
                f"CiRoot requests {raw_nroots} states of multiplicity {raw_spin}, but "
                f"CAS({nactel},{nactorb}) contains only {available} of them.",
                f"lower the root count to <= {available} (the source's 4f^9 9-in-7 space "
                f"has {SI_DATASET_CSFS}).",
            )
        pairs.append((raw_spin, raw_nroots))
    return tuple(sorted(pairs, key=lambda pair: -pair[0]))


def _check_one_line(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise _error(
            f"{name} must be a non-empty string, got {value!r}.",
            f"give {name} as a single line of text (it is pasted into the input file).",
        )
    text = value.strip()
    try:
        text.encode("ascii")
    except UnicodeEncodeError as exc:
        raise _error(
            f"{name}={text!r} contains non-ASCII characters "
            f"({exc.object[exc.start:exc.end]!r}); the OpenMolcas input must be ASCII.",
            f"rewrite {name} in plain ASCII.",
        ) from exc
    if any(character.isspace() for character in text):
        raise _error(
            f"{name}={text!r} contains whitespace, which would end the value early.",
            f"give {name} without spaces or line breaks.",
        )
    return text


def _check_origin(value: object) -> tuple[float, float, float]:
    if (
        not isinstance(value, (tuple, list))
        or len(value) != 3
        or any(isinstance(item, bool) or not isinstance(item, (int, float)) for item in value)
    ):
        raise _error(
            f"angmom_origin must be three numbers (x, y, z), got {value!r}.",
            "give the angular-momentum origin in Angstrom, e.g. (0.0, 0.0, 0.0).",
        )
    origin = tuple(float(item) for item in value)
    if any(not math.isfinite(item) for item in origin):
        raise _error(
            f"angmom_origin={value!r} contains a non-finite number.",
            "give three finite coordinates in Angstrom.",
        )
    return cast("tuple[float, float, float]", origin)


# --- CI-space counting ---------------------------------------------------


def max_csfs(nactel: int, nactorb: int, spin_multiplicity: int) -> int:
    """Number of spin-adapted CSFs of one spin multiplicity in CAS(nactel, nactorb).

    Weyl-Palmer binomial count ``(2S+1)/(m+1) * C(m+1, n/2-S) * C(m+1, n/2+S+1)`` with
    ``m = nactorb``; returns 0 when the combination is not allowed (wrong spin parity,
    spin above the unpaired-electron limit, or electrons exceeding the orbitals).
    Cross-checked against the source's numbers for 4f^9 in 7 orbitals: 21 sextets,
    224 quartets, 490 doublets (source section 4.2, item 6).
    """
    if nactel < 0 or nactorb < 1 or spin_multiplicity < 1:
        return 0
    if nactel > 2 * nactorb:
        return 0
    two_s = spin_multiplicity - 1
    if (nactel - two_s) % 2:
        return 0
    lower = (nactel - two_s) // 2
    upper = (nactel + two_s) // 2 + 1
    if lower < 0 or upper < 0 or lower > nactorb + 1 or upper > nactorb + 1:
        return 0
    return (
        (two_s + 1) * math.comb(nactorb + 1, lower) * math.comb(nactorb + 1, upper)
    ) // (nactorb + 1)


def spin_orbit_state_count(spec: OpenMolcasChainSpec) -> int:
    """Spin-orbit states the ``&RASSI SPIN`` step can produce: sum of (2S+1) x nroots.

    Arithmetic use only -- the source does not state it. The source example gives
    6x21 + 4x48 + 2x32 = 382 states; compare that with the number of states the RASSI and
    SINGLE_ANISO steps report.
    """
    return sum(multiplicity * nroots for multiplicity, nroots in spec.spin_blocks)


# --- evidence ------------------------------------------------------------


def evidence() -> tuple[Evidence, ...]:
    """Provenance of everything this module takes from the source (G3 template)."""
    return (
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "Block chain &GATEWAY -> &SEWARD -> 3 x &RASSCF -> &RASSI -> &SINGLE_ANISO "
                "with 28 line annotations; CiRoot= N N 1 (root count twice, unity "
                "weighting); the Typeindex / FILEORB / >> COPY chain handing each spin "
                "block the previous block's orbitals; Nactel= 9 0 0 with RAS2= 7 for the "
                "4f^9 case; &RASSI SPIN / MEES / EPRG= 7.0D-1 / 'Nr of JobIph= 3 ALL' / "
                "IPHN / PROP with 'ANGMOM' 1..3; &SINGLE_ANISO CRYS / MLTP= 1; 2 / "
                "TINT= 0.0 330.0 330 0.0001 / HINT= 0.0 10.0 201 / "
                "TMAG= 6 1.8 2 4 5 10 20 / QUAX= 1."
            ),
            ref=SI_SKELETON_REF,
            url=SI_DOI_URL,
            bibkey=SI_BIBKEY,
        ),
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "Default recipe behind the generated job: experimental crystal structure "
                "without geometry optimisation, isolated ion without counter-ions, "
                "SA-CASSCF-SO with OpenMolcas, 9-in-7 active space for 4f^9, ANO-RCC "
                "three-tier basis (Dy at VTZP quality, first coordination sphere at VDZP, "
                "remaining atoms at VDZ), spin-free space of " + SI_DATASET_CSFS + "."
            ),
            ref=SI_PROTOCOL_REF,
            url=SI_DOI_URL,
            bibkey=SI_BIBKEY,
        ),
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "Cross-check of the CI-space bound this module enforces: the Weyl-Palmer "
                "count reproduces the source's " + SI_DATASET_CSFS + " for 9 electrons in "
                "7 orbitals."
            ),
            ref=SI_PROTOCOL_REF,
            url=SI_DOI_URL,
            bibkey=SI_BIBKEY,
        ),
        Evidence(
            kind=EVIDENCE_LITERATURE,
            text=(
                "Hand-check recommended after SINGLE_ANISO for mononuclear Dy(III): "
                "g_T = (g_1 + g_2 + g_3 sin(theta_3)) / 3 with theta_3 in degrees, compared "
                "with g_T * theta_3 = 20 (below: the doublet can still act as a step of the "
                "barrier; above: ground-state quantum tunnelling is opened). Calibrated on "
                "19 mononuclear Dy(III) systems only, with one known outlier "
                "([Dy(Cp^ttt)_2]+, whose QTM doublet falls in the safe range)."
            ),
            ref=SI_CRITERION_REF,
            url=SI_DOI_URL,
            bibkey=SI_BIBKEY,
        ),
    )


def own_deviations() -> tuple[Deviation, ...]:
    """Choices this module makes that the source input does not show.

    Each tag appears in the rendered file as ``[fBlockKit:<tag>]``, so the file says which
    of its lines are not from the source. The provenance of these items is this module and
    the internal close-reading record ``文献细读/细读_Chilton指南_SI.md`` (2026-09-25) --
    they are not manual or literature evidence, except G3-1, which was verified against
    the OpenMolcas documentation on 2026-09-25 (see the module docstring).
    """
    return (
        Deviation(
            tag="G3-1",
            text=(
                "Comment character '*': VERIFIED against the OpenMolcas documentation "
                "(users.guide/emil.html, checked 2026-09-25): 'a comment line identified by "
                "an asterisk (*) in the first position on the line'; comments may be inserted "
                "anywhere except between a keyword and its following additional "
                "specifications. This module places comments only on their own lines."
            ),
        ),
        Deviation(
            tag="G3-2",
            text=(
                ">> COPY of the RASSCF orbital file ($Project.RasOrb) to the name the next "
                "block reads via FILEORB: the source shows >> COPY only for the IPH and "
                "rasscf.h5 files, yet its later blocks read 1.RasOrb / 2.RasOrb."
            ),
        ),
        Deviation(
            tag="G3-3",
            text=(
                "group= default 'C1' (no symmetry); the source example writes "
                "'group= NoSym' instead."
            ),
        ),
        Deviation(
            tag="G3-4",
            text=(
                "Spin blocks are emitted in decreasing multiplicity with each block reading "
                "the previous block's orbital file; this follows the source's example "
                "order, which the source does not state as a rule."
            ),
        ),
        Deviation(
            tag="G3-5",
            text=(
                "&SINGLE_ANISO sampling parameters (MLTP / TINT / HINT / TMAG / QUAX) are "
                "fixed to the source's values; they are not options in this version."
            ),
        ),
        Deviation(
            tag="G3-6",
            text=(
                "One 'ANGMOM' property per spin block in &RASSI PROP: the source asks for "
                "three (three blocks); generalised here to N blocks by the same pattern."
            ),
        ),
        Deviation(
            tag="G3-7",
            text=(
                "The source's line annotations are placed on the lines the internal "
                "close-reading aligns them to; that alignment is an inference from the "
                "source's printed order, not a printed fact."
            ),
        ),
        Deviation(
            tag="G3-8",
            text=(
                "Keyword spellings (Nactel=, RAS2=, CiRoot=, Typeindex, EPRG=) are "
                "reproduced verbatim from the source and were not cross-checked against the "
                "OpenMolcas documentation for this version."
            ),
        ),
    )


# --- rendering -----------------------------------------------------------


def _environment() -> Environment:
    return Environment(
        loader=FileSystemLoader(str(TEMPLATES_DIR)),
        undefined=StrictUndefined,
        trim_blocks=True,
        lstrip_blocks=True,
        keep_trailing_newline=True,
    )


def _title(spec: OpenMolcasChainSpec) -> str:
    if spec.title:
        return spec.title
    blocks = ", ".join(f"{spin}x{nroots}" for spin, nroots in spec.spin_blocks)
    return (
        f"{spec.element} magnetic chain, CAS({spec.nactel},{spec.nactorb}), "
        f"spin blocks {blocks}"
    )


def render_openmolcas_input(spec: OpenMolcasChainSpec) -> str:
    """Render the complete OpenMolcas input text for ``spec`` (ASCII, deterministic)."""
    blocks = []
    last = len(spec.spin_blocks)
    for index, (spin, nroots) in enumerate(spec.spin_blocks, start=1):
        blocks.append(
            {
                "index": index,
                "spin": spin,
                "nroots": nroots,
                "iph": f"{index}_IPH",
                "h5": f"{index}.rasscf.h5",
                "prev_rasorb": f"{index - 1}.RasOrb" if index > 1 else "",
                "next_rasorb": f"{index}.RasOrb" if index < last else "",
                "first": index == 1,
            }
        )
    context = {
        "title": _title(spec),
        "source_note": SI_SOURCE_NOTE,
        "spin_block_summary": "; ".join(
            f"Spin= {spin} with CiRoot= {nroots}"
            for spin, nroots in spec.spin_blocks
        ),
        "element": spec.element,
        "element_lower": spec.element.lower(),
        "charge": spec.charge,
        "nactel": spec.nactel,
        "nactorb": spec.nactorb,
        "blocks": blocks,
        "coord_file": spec.coord_file,
        "basis_metal": spec.basis_metal,
        "basis_ligand": spec.basis_ligand,
        "basis_outer": spec.basis_outer,
        "symmetry": spec.symmetry,
        "angmom_origin": " ".join(f"{value:.1f}" for value in spec.angmom_origin),
        "single_aniso_mltp": SINGLE_ANISO_MLTP,
        "single_aniso_tint": SINGLE_ANISO_TINT,
        "single_aniso_hint": SINGLE_ANISO_HINT,
        "single_aniso_tmag": SINGLE_ANISO_TMAG,
        "single_aniso_quax": SINGLE_ANISO_QUAX,
        "rassi_eprg": RASSI_EPRG,
    }
    text = _environment().get_template(TEMPLATE_NAME).render(**context)
    try:
        text.encode("ascii")
    except UnicodeEncodeError as exc:
        # Guard for future template edits: every spec field is ASCII-validated above.
        raise _error(
            "the rendered input contains non-ASCII characters "
            f"({exc.object[exc.start:exc.end]!r}); an OpenMolcas input must be ASCII.",
            "check the template and the spec fields for non-ASCII text.",
        ) from exc
    return text


def run_guidance_openmolcas(spec: OpenMolcasChainSpec) -> str:
    """Plain-text guide: how to run the job and what to check in the output by hand.

    v0.1 of this toolkit does not parse OpenMolcas output -- the diagnosis layer reads ORCA
    output only -- so every check below is manual, and the guide says so.
    """
    states = spin_orbit_state_count(spec)
    arithmetic = " + ".join(
        f"{multiplicity}x{nroots}"
        for multiplicity, nroots in spec.spin_blocks
    )
    source_example_states = sum(
        multiplicity * nroots for multiplicity, nroots in SI_EXAMPLE_SPIN_BLOCKS
    )
    blocks = "\n".join(
        f"     block {index}: Spin= {spin}, CiRoot= {nroots} {nroots} 1"
        for index, (spin, nroots) in enumerate(spec.spin_blocks, start=1)
    )
    return "\n".join(
        [
            f"OpenMolcas magnetic chain for {spec.element} "
            f"(CAS({spec.nactel},{spec.nactorb}), charge {spec.charge})",
            "",
            "1. Before running",
            f"   - Write the structure-and-basis file {spec.coord_file}: one line per atom,",
            "     with the basis label in front of the coordinates, as the source example",
            "     does (geometry and per-atom basis in one file). Assign the three tiers:",
            f"       {spec.element}: {spec.basis_metal}",
            f"       first coordination sphere: {spec.basis_ligand}",
            f"       all remaining atoms: {spec.basis_outer}",
            "     Use the experimental structure without geometry optimisation, and the",
            "     isolated ion without counter-ions, as in the source protocol.",
            "   - Run it with your installation's OpenMolcas driver (normally 'pymolcas')",
            "     from a scratch directory you can write to. The job is one input file with",
            f"     {len(spec.spin_blocks)} RASSCF blocks, then RASSI, then SINGLE_ANISO:",
            blocks,
            "   - The '>> COPY' lines are shell commands run between modules; they must be",
            "     allowed to execute in that directory.",
            "",
            "2. Reading the output (every check below is manual in v0.1)",
            "   - This version does not parse OpenMolcas output: the diagnosis layer of",
            "     fBlockKit reads ORCA output only, so nothing below is done for you.",
            "   - RASSCF: every block must report convergence of the orbital optimisation;",
            "     a block that stops at the iteration limit makes the rest of the chain and",
            "     SINGLE_ANISO unreliable. Check as well that the active orbitals really",
            "     are the seven f orbitals, and that the orbital sets passed between blocks",
            "     look like the same physical solution (the later blocks start from the",
            "     previous block's file).",
            "   - RASSI: check the number of spin-orbit states actually coupled against the",
            "     arithmetic expectation for this specification, sum (2S+1) x nroots =",
            f"     {arithmetic} = {states} states. The source example couples "
            f"{source_example_states} states.",
            "     Check also that all requested roots of every block entered the mixing: a",
            "     block silently contributing fewer roots changes the whole spectrum.",
            "   - SINGLE_ANISO: check the number of Kramers doublets (or pseudospin",
            "     multiplets) against the spin-orbit states you coupled, then the g tensors.",
            "     For a Dy(III) centre the ground doublet should be strongly axial (g_3",
            "     close to 20, g_1 and g_2 small); a g_3 far from that means the ground",
            "     doublet is not the expected |+/-15/2> state, so check the state ordering",
            "     and the pseudospin convention before using the numbers.",
            "   - Optional hand criterion for mononuclear Dy(III) (source sections 3.1 and",
            "     4.3): g_T = (g_1 + g_2 + g_3 sin(theta_3)) / 3 with theta_3 in degrees,",
            "     compared with g_T * theta_3 = 20 -- below that line the excited doublet",
            "     can still act as a step of the barrier, above it ground-state quantum",
            "     tunnelling is opened. The criterion was calibrated on 19 mononuclear",
            "     Dy(III) systems only, so do not transfer it to other ions or to",
            "     polynuclear systems.",
            "",
            "3. What this input does not cover",
            "   - No geometry optimisation and no counter-ions (source protocol:",
            "     experimental structure, isolated ion).",
            "   - No magnetic relaxation rates and no crystal-field parameter projection;",
            "     projected B_k^q values depend on the projection manifold, so they cannot",
            "     be compared across manifolds without saying which one was used.",
            "   - If your OpenMolcas version rejects the comment lines that start with '*'",
            "     (see the [fBlockKit:G3-1] note in the file), delete them: the source's own",
            "     input contains no comment lines and runs without them.",
            "   - The generated file does not contain the structure itself: the geometry and",
            f"     the per-atom basis labels must be in {spec.coord_file}.",
        ]
    )


# The G3 guide under the name used by the feature specification. It is exported from the
# package as ``run_guidance_openmolcas`` because ``fblockkit.recipe.run_guidance`` is
# already the ORCA guide.
run_guidance = run_guidance_openmolcas
