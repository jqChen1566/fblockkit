"""Checks of the plot-ready CSV convention's shared formatting rules
(analysis/plot_csv.py; the formats chapter's section "Plot-ready CSV
companions").

The convention: one header line, one data row per point, ten significant
digits, no comment lines; text fields are quoted only when they would
otherwise break the row.
"""

from __future__ import annotations

from fblockkit.analysis import plot_csv


def test_a_plain_field_is_written_verbatim():
    assert plot_csv.csv_field("5D1 <- 7F1") == "5D1 <- 7F1"
    assert plot_csv.csv_field("7F6") == "7F6"


def test_a_field_with_a_comma_is_quoted():
    assert plot_csv.csv_field("5D0->7F2, electric") == '"5D0->7F2, electric"'


def test_quotes_are_doubled_inside_a_quoted_field():
    assert plot_csv.csv_field('say "hi", ok') == '"say ""hi"", ok"'


def test_a_field_with_a_line_break_is_quoted():
    assert plot_csv.csv_field("two\nlines") == '"two\nlines"'


def test_plot_numbers_carry_ten_significant_digits():
    assert plot_csv.plot_number(0.00032241059450348005) == "0.0003224105945"
    assert plot_csv.plot_number(715.072591) == "715.072591"
    assert plot_csv.plot_number(1.0) == "1"


def test_non_finite_values_are_written_as_their_tokens():
    assert plot_csv.plot_number(float("-inf")) == "-inf"
    assert plot_csv.plot_number(float("nan")) == "nan"
