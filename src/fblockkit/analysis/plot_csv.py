"""The plot-ready CSV convention in code (the formats chapter's section
"Plot-ready CSV companions").

Menus whose report is a curve or a table of plottable points write a
``<output>.<menu>.fbk.csv`` companion beside the report: one header line
(column names with units, underscored), one data row per point, comma
separated, ten significant digits, no comment lines; the report beside it
carries the context.  Any plotting tool (spreadsheets, Origin, pandas,
gnuplot) reads the file as it stands.

This module holds the convention's two formatting rules so they have one
home in code:

- :func:`plot_number` -- a numeric column value at ten significant digits
  (non-finite values are written as the tokens they are: ``inf``, ``-inf``,
  ``nan``);
- :func:`csv_field` -- a text column value with minimal RFC 4180 quoting,
  applied only when the field would otherwise break the row (it contains a
  comma, a quote or a line break).
"""

from __future__ import annotations

__all__ = ["csv_field", "plot_number"]


def plot_number(value: float) -> str:
    """One numeric column value (ten significant digits, the companion format)."""
    return f"{float(value):.10g}"


def csv_field(value: object) -> str:
    """A text column value, quoted only when its content requires quoting."""
    text = str(value)
    if any(character in text for character in ',"\r\n'):
        return '"' + text.replace('"', '""') + '"'
    return text
