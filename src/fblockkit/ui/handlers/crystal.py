"""The crystal handler group."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from ...analysis import cf_declaration, crystal_field, geometry as geometry_analysis, point_charge
from ...diagnosis import references_section
from ..session import Session
# --- 11 point-charge crystal-field estimate ---------------------------------


def point_charge_estimate(session: Session) -> None:
    """S2: point-charge crystal-field estimate from a structure (XYZ + charges)."""
    path_text = session.ask("XYZ structure file path")
    if not path_text:
        session.say("Cancelled (no path given).")
        return
    path = Path(path_text)
    try:
        atoms = geometry_analysis.parse_xyz(path)
    except geometry_analysis.StructureError as exc:
        session.say(f"Structure read failed: {exc}")
        return
    centre = atoms[geometry_analysis.dominant_center(atoms)]
    session.say(
        f"Centre: {centre.element} at ({centre.x}, {centre.y}, {centre.z}) "
        "(the dominant centre of the structure)"
    )
    charges_text = session.ask("Point charges per element (e.g. O=-2,H=0.4)")
    charges: dict[str, float] = {}
    for chunk in (charges_text or "").split(","):
        if not chunk.strip():
            continue
        symbol, separator, value = chunk.partition("=")
        if not separator:
            session.say(f"Charge entry {chunk.strip()!r} is not element=charge; skipped.")
            continue
        try:
            charges[symbol.strip()] = float(value)
        except ValueError:
            session.say(f"Charge value {value.strip()!r} is not a number; skipped.")
    radial_text = session.ask(
        "Radial moments r2,r4,r6 in Angstrom^k (Enter = geometry-only output)", default=""
    )
    radial: dict[int, float] | None = None
    if radial_text:
        values = [item.strip() for item in radial_text.split(",")]
        if len(values) != 3:
            session.say("Expected three values (r2,r4,r6); giving geometry-only output.")
        else:
            try:
                radial = {k: float(v) for k, v in zip((2, 4, 6), values)}
            except ValueError:
                session.say("Radial moments must be numbers; giving geometry-only output.")
                radial = None
    try:
        estimate = point_charge.estimate(
            atoms, charges, (centre.x, centre.y, centre.z), radial
        )
        section = point_charge.report(estimate, charges, radial)
    except point_charge.PointChargeError as exc:
        session.say(f"Point-charge estimate refused: {exc}")
        return
    session.say(section.body)
    body = f"## {section.title}\n\n{section.body}\n"
    refs = references_section(point_charge.evidence())
    if refs is not None:
        body += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = path.with_name(path.name + ".fbk.md")
    md_path.write_text(body, encoding="utf-8")
    session.say(f"Report written: {md_path}")
# --- 10 crystal-field fit ---------------------------------------------------


def _coefficient_value(entry) -> complex:
    """A coefficient may be a real number or a [re, im] pair."""
    if isinstance(entry, (int, float)):
        return complex(entry)
    if isinstance(entry, (list, tuple)) and len(entry) == 2:
        return complex(entry[0], entry[1])
    raise ValueError(f"coefficient {entry!r} must be a number or a [re, im] pair")


def crystal_field_fit(session: Session) -> None:
    """Fit B_k^q from a levels + coefficients JSON file and write a report."""
    path_text = session.ask("CF input JSON path (levels + coefficients; see the user guide)")
    if not path_text:
        session.say("Cancelled (no path given).")
        return
    path = Path(path_text)
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        point_group = str(data["point_group"])
        j_value = float(data["J"])
        levels = [float(value) for value in data["levels"]]
        coefficients = [
            [_coefficient_value(entry) for entry in row] for row in data["coefficients"]
        ]
        result = crystal_field.fit_crystal_field(levels, coefficients, point_group, j_value)
    except (OSError, KeyError, TypeError, ValueError) as exc:
        session.say(
            f"CF fit failed: {exc}. Next step: provide a JSON with keys "
            "point_group, J, levels (list) and coefficients (list of rows of "
            "numbers or [re, im] pairs)."
        )
        return

    lines = [
        f"Point group {result.point_group}, J = {result.j}",
        f"Sampled states: {result.n_states}; fitted parameters: {result.n_parameters}",
        f"Condition number: {result.condition_number:.3g}; rank {result.rank}",
        f"max |residual| = {result.max_abs_residual:.4g} (caller's energy unit)",
        "",
        "B_k^q:",
    ]
    for (k, q), value in sorted(result.parameters.items()):
        lines.append(f"  ({k},{q:>2})  {value: .10g}")
    lines.append(f"  const    {result.const: .10g}")
    for note in result.notes:
        lines.append(f"note: {note}")
    body = "\n".join(lines)
    session.say(body)

    # A5: projection-basis declaration check on the same JSON (the five items any
    # B_k^q set must record before it may be compared with another scheme)
    a5_payload: dict[str, Any] = {
        "parameters": {f"{k},{q}": value for (k, q), value in result.parameters.items()},
        "declaration": data.get("declaration") or {},
    }
    if isinstance(data.get("compare"), dict):
        a5_payload["compare"] = data["compare"]
    a5_section = None
    try:
        a5_section = cf_declaration.run(a5_payload)
    except cf_declaration.DeclarationError as exc:
        session.say(f"Projection-basis declaration check refused: {exc}")

    report_lines = f"## A4 crystal-field fit\n\n{body}\n"
    if a5_section is not None:
        report_lines += f"\n## {a5_section.title}\n\n{a5_section.body}\n"
    refs = references_section(crystal_field.evidence() + cf_declaration.evidence())
    if refs is not None:
        report_lines += f"\n## {refs.title}\n\n{refs.body}\n"
    md_path = path.with_name(path.name + ".fbk.md")
    md_path.write_text(report_lines, encoding="utf-8")
    session.say(f"Report written: {md_path}")
    if a5_section is not None:
        declared = sum(1 for item in cf_declaration.check_declaration(a5_payload["declaration"]) if item.present)
        session.say(
            f"Projection-basis declaration: {declared} of "
            f"{len(cf_declaration.DECLARATION_FIELDS)} item(s) recorded (section A5 in the report)."
        )
        if "compare" in a5_payload:
            refusals = [
                item
                for item in cf_declaration.compare_parameters(
                    {
                        "parameters": a5_payload["parameters"],
                        "declaration": a5_payload["declaration"],
                    },
                    a5_payload["compare"],
                )
                if item.verdict == cf_declaration.VERDICT_REFUSED
            ]
            session.say(
                f"Parameter comparison against the bundled set: "
                f"{len(refusals)} (k, q) refused without a matching declaration."
            )


