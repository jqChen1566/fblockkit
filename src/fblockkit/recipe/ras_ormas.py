"""4.8: RAS / ORMAS (generalized active space) input generation (menu 28).

ORCA has no GAS keyword; its model-space counterpart of the generalized
active space concept is the partition masks carried by the ``refs``
sub-block.  The measured grammar (ORCA 6.1 manual sections 3.13.3.6-3.13.3.9,
its own examples; every point below re-measured on the 6.1.1 engine):

- ``refs`` is a sub-block closed by its own ``end``: without the inner
  ``end`` the parser reports the next token as an unrecognized symbol of the
  enclosing block (measured both for %casscf and %rasci);
- the RAS mask is ``RAS(Nel: NRAS1 MaxHoles / NRAS2 / NRAS3 MaxParticles)``;
  the three orbital counts must sum to ``norb`` exactly (engine refusal with
  its verbatim message otherwise).  ``NRAS3 0`` (no RAS3) is legal -- the
  manual's RIXS example uses it;
- the ORMAS mask is ``ORMAS(nel: m1 min1 max1, m2 min2 max2, ...)`` with
  commas or slashes between sub-spaces (both measured), up to 25 sub-spaces
  (the manual's limit; the engine prints a maximum-subspace note).  The mask
  *overrides* the ``nel``/``norb`` lines (the manual: "the standard CASSCF
  parameters nel and norb are overwritten with the parameters from the ORMAS
  input mask"; measured: nel 4/norb 6 against a nel-6 mask ran the mask's
  space) -- the generator therefore refuses inconsistent numbers instead of
  letting ORCA silently override them;
- per-sub-space engine checks: minimum <= maximum <= 2*m, m >= 1, minimum >=
  0; across sub-spaces the mask must be able to carry ``nel`` electrons
  (measured refusal: "The sum of maximum occupations (=4) is smaller than
  the number of active electrons");
- ``%casscf`` + ``refs`` runs the orbital-optimized route (RASSCF /
  ORMAS-SCF); the manual warns that for incomplete model spaces the orbital
  optimization omits the active-active rotation, so the energy is sensitive
  to the canonicalization choice;
- ``%rasci`` + ``refs`` runs the standalone CI-only route (the manual's
  MRCISD-style example); the mask then defines a contiguous window of its
  orbital count after the frozen core (measured echo: "Number of active
  orbitals ... 6", "First active orbital ... 4" for N2/def2-SVP).

Measured anchors (N2/def2-SVP at 1.10 Angstrom, all engine runs):
plain CASSCF(6,6) -108.989034756374 Eh; ORMAS(6: 2 0 4, 2 0 4, 2 0 4) --
whose sub-spaces ban nothing, each holding 2 orbitals with 0-4 electrons --
identical to all 12 printed digits; the restricted RAS(6:2 2/2/2 2) raises
the energy to -108.985623236851 Eh; the CI-only %rasci route with the same
RAS mask gives -108.921051085808 Eh (no orbital relaxation).
"""

from __future__ import annotations

import re

__all__ = [
    "RasOrmasError",
    "parse_ras_mask",
    "parse_ormas_mask",
    "format_ras_mask",
    "format_ormas_mask",
    "ras_ormas_input",
    "run_guidance_lines",
]

_MAX_ORMAS_SUBSPACES = 25  # the manual's limit for the ORMAS mask


class RasOrmasError(ValueError):
    """Invalid RAS/ORMAS parameters; the message carries the next step."""


def _ints(text: str) -> list[int]:
    fields = [item for item in re.split(r"[\s,/:]+", text.strip()) if item]
    if not fields:
        raise RasOrmasError("the mask is empty.")
    try:
        return [int(item) for item in fields]
    except ValueError as exc:
        raise RasOrmasError(f"the mask contains a non-integer field ({exc}).") from exc


def parse_ras_mask(text: str) -> tuple[int, int, int, int, int, int]:
    """Parse ``<nel>:<n1> <h1>/<n2>/<n3> <p3>`` into its six integers.

    This is the manual's mask ``RAS(Nel: NRAS1 MaxHoles / NRAS2 / NRAS3
    MaxParticles)`` without the surrounding keyword: the colon and the three
    slashes are the only punctuation that matters, spacing is free.
    """
    stripped = text.strip()
    if not stripped:
        raise RasOrmasError(
            "the RAS mask is empty. Next step: give it like the manual -- "
            "'3 1/5/0 0' or '2 2/2/2 2' (NRAS1 MaxHoles / NRAS2 / NRAS3 "
            "MaxParticles)."
        )
    nel_text, _, rest = stripped.partition(":")
    if not rest:
        raise RasOrmasError(
            f"the RAS mask {text!r} has no ':' after the active-electron count. "
            "Give it as '<nel>:<n1> <h1>/<n2>/<n3> <p3>'."
        )
    parts = rest.split("/")
    if len(parts) != 3:
        raise RasOrmasError(
            f"the RAS mask {text!r} has {len(parts) - 1} slash(es) after the colon; "
            "the manual's form is '<nel>:<n1> <h1>/<n2>/<n3> <p3>' (two slashes: "
            "the RAS1 pair, RAS2, and the RAS3 pair)."
        )
    spans = [[field for field in part.split()] for part in parts]
    if [len(span) for span in spans] != [2, 1, 2]:
        raise RasOrmasError(
            f"the RAS mask {text!r} has the wrong field counts; the manual's form "
            "is '<nel>:<n1> <h1>/<n2>/<n3> <p3>' (two numbers, one, two)."
        )
    try:
        nel = int(nel_text.strip())
        n1, h1 = (int(item) for item in spans[0])
        n2 = int(spans[1][0])
        n3, p3 = (int(item) for item in spans[2])
    except ValueError as exc:
        raise RasOrmasError(f"the RAS mask contains a non-integer field ({exc}).") from exc
    if nel < 0:
        raise RasOrmasError(f"the active-electron count ({nel}) cannot be negative.")
    for label, value in (("NRAS1", n1), ("NRAS2", n2), ("NRAS3", n3)):
        if value < 0:
            raise RasOrmasError(f"{label} ({value}) cannot be negative.")
    if h1 < 0 or p3 < 0:
        raise RasOrmasError(
            f"the hole/particle limits cannot be negative (MaxHoles {h1}, "
            f"MaxParticles {p3})."
        )
    return nel, n1, h1, n2, n3, p3


def format_ras_mask(fields: tuple[int, int, int, int, int, int]) -> str:
    """The canonical mask text (the manual's own spelling)."""
    nel, n1, h1, n2, n3, p3 = fields
    return f"RAS({nel}:{n1} {h1}/{n2}/{n3} {p3})"


def parse_ormas_mask(text: str) -> tuple[int, tuple[tuple[int, int, int], ...]]:
    """Parse ``<nel>: <m1> <min1> <max1>, <m2> <min2> <max2>, ...``.

    Returns the mask's electron count and the sub-space triples.  Commas and
    slashes both separate sub-spaces (measured); each sub-space is exactly
    three integers: orbital count, minimum electrons, maximum electrons.
    """
    stripped = text.strip()
    if not stripped:
        raise RasOrmasError(
            "the ORMAS mask is empty. Next step: give it like the manual -- "
            "'6: 2 0 4, 2 0 4, 2 0 4' (per sub-space: orbitals, min electrons, "
            "max electrons)."
        )
    nel_text, _, rest = stripped.partition(":")
    if not rest:
        raise RasOrmasError(
            f"the ORMAS mask {text!r} has no ':' after the active-electron count. "
            "Give it as '<nel>: <m1> <min1> <max1>, <m2> <min2> <max2>, ...'."
        )
    groups = [group for group in re.split(r"[,/]", rest) if group.strip()]
    if not groups:
        raise RasOrmasError(f"the ORMAS mask {text!r} defines no sub-space.")
    if len(groups) > _MAX_ORMAS_SUBSPACES:
        raise RasOrmasError(
            f"the ORMAS mask defines {len(groups)} sub-spaces; the manual's limit "
            f"is {_MAX_ORMAS_SUBSPACES}."
        )
    subspaces = []
    for group in groups:
        fields = _ints(group)
        if len(fields) != 3:
            raise RasOrmasError(
                f"the ORMAS sub-space {group.strip()!r} has {len(fields)} field(s); "
                "each sub-space is exactly three integers (orbitals, min, max)."
            )
        m, low, high = fields
        if m < 1:
            raise RasOrmasError(f"a sub-space with {m} orbital(s) is impossible.")
        if low < 0 or high < 0:
            raise RasOrmasError(
                f"the sub-space {group.strip()!r} has a negative electron bound."
            )
        if high > 2 * m:
            raise RasOrmasError(
                f"the sub-space {group.strip()!r} asks for up to {high} electrons in "
                f"{m} orbital(s); the engine's own check is 'maximum number of "
                "electrons >2*morb'."
            )
        if low > high:
            raise RasOrmasError(
                f"the sub-space {group.strip()!r} has its minimum above its maximum."
            )
        subspaces.append((m, low, high))
    nel = int(nel_text.strip()) if nel_text.strip().lstrip("+-").isdigit() else None
    if nel is None:
        raise RasOrmasError(
            f"the ORMAS mask {text!r} does not start with an active-electron count."
        )
    if nel < 0:
        raise RasOrmasError(f"the active-electron count ({nel}) cannot be negative.")
    total_min = sum(low for _, low, _ in subspaces)
    total_max = sum(high for _, _, high in subspaces)
    if nel > total_max:
        # the engine's own wording for this refusal
        raise RasOrmasError(
            f"the sum of maximum occupations (={total_max}) is smaller than the "
            f"number of active electrons ({nel}) -- the engine itself rejects this; "
            "widen a sub-space or lower the electron count."
        )
    if nel < total_min:
        raise RasOrmasError(
            f"the sum of minimum occupations (={total_min}) is larger than the "
            f"number of active electrons ({nel}); the mask cannot carry that few "
            "electrons."
        )
    return nel, tuple(subspaces)


def format_ormas_mask(nel: int, subspaces: tuple[tuple[int, int, int], ...]) -> str:
    """The canonical mask text (the manual's own spelling)."""
    body = ", ".join(f"{m} {low} {high}" for m, low, high in subspaces)
    return f"ORMAS({nel}: {body})"


def _int_list(text: str, label: str) -> tuple[int, ...]:
    fields = [item for item in text.replace(" ", "").split(",") if item != ""]
    if not fields:
        raise RasOrmasError(f"the {label} list is empty.")
    try:
        values = tuple(int(item) for item in fields)
    except ValueError as exc:
        raise RasOrmasError(f"the {label} list is not a comma list of integers.") from exc
    if any(value < 1 for value in values):
        raise RasOrmasError(f"the {label} list carries a non-positive entry.")
    return values


def ras_ormas_input(
    coordinates: tuple[tuple[str, float, float, float], ...],
    *,
    charge: int,
    multiplicity: int,
    space: str,
    mask: str,
    nel: int | None = None,
    norb: int | None = None,
    route: str = "casscf",
    mult: str | None = None,
    nroots: str = "1",
    cistep: str | None = None,
    exc_level: int | None = None,
    keywords: str = "RHF def2-SVP",
    maxcore: int = 2000,
) -> str:
    """Render the RAS/ORMAS input (one partition mask, one route).

    The mask is the single source of the active-space numbers (the manual's
    masks carry the electron count, and the orbital count is their sum); the
    ``nel``/``norb`` arguments are optional cross-checks so a caller can
    assert what the user meant -- a mismatch is refused with the engine's own
    wording instead of being silently overridden.  ``route='casscf'``
    optimizes the orbitals inside the partition (RASSCF / ORMAS-SCF);
    ``route='rasci'`` is the standalone CI-only module.
    """
    if space not in ("ras", "ormas"):
        raise RasOrmasError(f"unknown partition type {space!r}; use ras or ormas.")
    if route not in ("casscf", "rasci"):
        raise RasOrmasError(f"unknown route {route!r}; use casscf or rasci.")
    if multiplicity < 1:
        raise RasOrmasError(f"multiplicity {multiplicity} is not positive.")

    if space == "ras":
        nel_m, n1, h1, n2, n3, p3 = parse_ras_mask(mask)
        total = n1 + n2 + n3
        if norb is not None and norb != total:
            raise RasOrmasError(
                f"sum of number of orbitals in RAS1 ({n1}), RAS2 ({n2}), and RAS3 "
                f"({n3}) must sum up to NORB ({norb}) -- the engine's own refusal; "
                "fix norb or the mask."
            )
        if nel is not None and nel_m != nel:
            raise RasOrmasError(
                f"the mask carries {nel_m} active electrons but nel is {nel}; give "
                "one number (the mask's first field is the active-electron count)."
            )
        mask_text = format_ras_mask((nel_m, n1, h1, n2, n3, p3))
        resolved_nel, resolved_norb = nel_m, total
    else:
        nel_m, subspaces = parse_ormas_mask(mask)
        total = sum(m for m, _, _ in subspaces)
        if norb is not None and norb != total:
            raise RasOrmasError(
                f"the ORMAS mask defines {total} orbital(s) but norb is {norb}; ORCA "
                "overwrites the nel/norb lines with the mask parameters, so the "
                "generator refuses the mismatch -- fix norb or the mask."
            )
        if nel is not None and nel != nel_m:
            raise RasOrmasError(
                f"the ORMAS mask carries {nel_m} active electrons but nel is {nel}; "
                "ORCA overwrites the nel line with the mask, so fix the number or "
                "the mask."
            )
        mask_text = format_ormas_mask(nel_m, subspaces)
        resolved_nel, resolved_norb = nel_m, total

    mult_values = _int_list(mult, "mult") if mult is not None else (multiplicity,)
    nroots_values = _int_list(nroots, "nroots")
    if len(mult_values) > 1 and len(nroots_values) not in (1, len(mult_values)):
        raise RasOrmasError(
            f"{len(mult_values)} multiplicities against {len(nroots_values)} nroots "
            "entries; give one nroots per multiplicity (or a single value)."
        )
    if route == "casscf" and len(mult_values) > 1 and len(nroots_values) == 1:
        # ORCA reads one nroots per mult block; broadcasting is unambiguous
        nroots_values = nroots_values * len(mult_values)
    if cistep is not None and cistep not in ("accci", "csfci", "detci", "treecsf"):
        raise RasOrmasError(
            f"unknown CIStep {cistep!r}; the %rasci module documents accci, csfci, "
            "detci and treecsf."
        )
    if exc_level is not None and exc_level < 0:
        raise RasOrmasError("ExcLevel cannot be negative.")

    lines = [f"! {keywords}", f"%maxcore {maxcore}"]
    if route == "casscf":
        lines.append("%casscf")
        lines.append(f"  nel {resolved_nel}")
        lines.append(f"  norb {resolved_norb}")
        lines.append(f"  mult {','.join(str(v) for v in mult_values)}")
        lines.append(f"  nroots {','.join(str(v) for v in nroots_values)}")
        lines.append("  refs")
        lines.append(f"    {mask_text}")
        lines.append("  end")
        lines.append("end")
    else:
        lines.append("%rasci")
        lines.append("  refs")
        lines.append(f"    {mask_text}")
        lines.append("  end")
        lines.append(f"  mult {','.join(str(v) for v in mult_values)}")
        lines.append(f"  nroots {','.join(str(v) for v in nroots_values)}")
        if cistep is not None:
            lines.append(f"  cistep {cistep}")
        if exc_level is not None:
            lines.append(f"  ExcLevel {exc_level}")
        lines.append("end")
    lines.append(f"* xyz {charge} {multiplicity}")
    for element, x, y, z in coordinates:
        lines.append(f"{element:<2} {x:>16.10f} {y:>16.10f} {z:>16.10f}")
    lines.append("*")
    text = "\n".join(lines) + "\n"
    try:
        text.encode("ascii")
    except UnicodeEncodeError as exc:  # pragma: no cover - callers keep ASCII input
        raise RasOrmasError(
            f"the generated input is not ASCII ({exc}); ORCA inputs must be pure ASCII."
        ) from exc
    return text


def run_guidance_lines(route: str = "casscf") -> tuple[str, ...]:
    """The source facts to print with every generated RAS/ORMAS input."""
    common = (
        "check the partition in the output before trusting the number: the %casscf "
        "route prints 'Building the RAS space ... done (N configurations)' or the "
        "ORMAS sub-space configuration counts; a wrong mask shows there",
    )
    if route == "rasci":
        return common + (
            "the %rasci module is a standalone CI: it takes the mask's orbital "
            "count as a contiguous window after the frozen core (the output echoes "
            "'Number of active orbitals' and 'First active orbital') and runs no "
            "orbital optimization -- its energy sits above the corresponding "
            "RASSCF unless you supply converged orbitals",
            "ExcLevel adds excitations on top of the references (the manual's "
            "MRCISD-style example relies on the default); mult/nroots take the "
            "CASSCF list syntax",
            "the module couples to the QDPT (magnetic) and OPA (optical) drivers "
            "via its rel/douv options; those are not generated here",
        )
    return common + (
        "the RASSCF/ORMAS-SCF orbital optimization omits the active-active "
        "rotation for incomplete model spaces (the manual's own warning), so the "
        "energy is sensitive to the orbital canonicalization (actorbs / "
        "actconstraints); compare spaces, not only energies",
        "for unconstrained sub-spaces ORMAS reproduces the full CAS to all printed "
        "digits (measured 12-digit identity on N2); any deviation from the CAS "
        "energy is the restriction you asked for",
        "ORCA silently overwrites nel/norb from an ORMAS mask -- this generator "
        "refuses inconsistent numbers instead, so the input you get always "
        "carries the mask's own parameters",
    )
