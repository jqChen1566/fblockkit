"""Tests for the orbital portrait (analysis/orbital_portrait.py).

The anchors are the N2 CAS(6,6) export's own chemistry: the two 1s-derived cores
come out non-bonding (their cross populations are ~1e-3, below the band), the pi
pair is bonding, the pi* pair is antibonding, and the one node-heavy virtual
orbital carries the cancellation flag.  The descriptor arithmetic is checked by
identities: the Loewdin atomic populations of a normalised orbital sum to one,
and the atom-pair blocks sum to the norm.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from fblockkit.analysis import orbital_portrait
from fblockkit.analysis.orbital_portrait import (
    CANCELLATION_LIMIT,
    PortraitError,
    analyze,
    evidence,
    render,
)
from fblockkit.parsers.orca_json import AoLabel, OrcaJson, parse_orca_json

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "orca"
CANONICAL = FIXTURES / "n2_fcidump.canonical.json"


def _export():
    return parse_orca_json(CANONICAL)


# --- the real-data anchors ----------------------------------------------------


def test_the_n2_window_reads_its_own_chemistry():
    portrait = analyze(_export(), window=range(12))
    labels = {row.index: row.bonding_label for row in portrait.rows}
    assert labels[0] == "non-bonding" and labels[1] == "non-bonding"  # the 1s cores
    assert labels[2] == "bonding" and labels[3] == "antibonding"
    assert labels[5] == "bonding" and labels[6] == "bonding"  # the pi pair
    assert labels[7] == "antibonding" and labels[8] == "antibonding"  # the pi* pair


def test_the_cross_terms_are_the_measured_bond_orders():
    portrait = analyze(_export(), window=range(12))
    by_index = {row.index: row for row in portrait.rows}
    assert by_index[5].bonding_total == pytest.approx(0.2572, abs=5e-4)
    assert by_index[7].bonding_total == pytest.approx(-0.3697, abs=5e-4)
    # the cores' cross terms are ~1e-3: below the band, hence non-bonding
    assert abs(by_index[0].bonding_total) < 1e-3


def test_the_cancellation_flag_fires_only_on_the_node_heavy_virtual():
    portrait = analyze(_export(), window=range(12))
    flagged = [row.index for row in portrait.rows if row.cancellation_heavy]
    assert flagged == [10]
    assert portrait.rows[10].bonding_total < -CANCELLATION_LIMIT


def test_the_lowdin_populations_sum_to_one():
    portrait = analyze(_export(), window=range(12))
    for row in portrait.rows:
        assert sum(row.atom_populations) == pytest.approx(1.0, abs=1e-9)
        assert max(row.atom_populations) == pytest.approx(row.dominant_share, abs=1e-12)


def test_the_shell_composition_weights_sum_to_one():
    portrait = analyze(_export(), window=[5])
    row = portrait.rows[0]
    assert sum(value for _, value in row.shell_weights) == pytest.approx(1.0, abs=1e-9)
    letters = dict(row.shell_weights)
    assert letters["p"] > 0.95  # the pi orbital is a p combination


def test_the_default_window_is_the_fractional_occupation_set():
    portrait = analyze(_export())
    assert portrait.window == (4, 5, 6, 7, 8, 9)


def test_the_report_states_the_criteria_and_the_estimates():
    body = render(analyze(_export(), window=[5]))
    assert "bonding label" in body or "bonding" in body
    assert "AO-centre estimates" in body
    assert "not adopted" in body  # the model is explicitly not adopted


# --- errors -------------------------------------------------------------------


def test_an_export_without_an_overlap_is_refused():
    export = parse_orca_json(CANONICAL)
    stripped = OrcaJson(
        base_name=export.base_name,
        charge=export.charge,
        multiplicity=export.multiplicity,
        hftyp=export.hftyp,
        point_group=export.point_group,
        atoms=export.atoms,
        coordinates=export.coordinates,
        n_mo=export.n_mo,
        n_ao=export.n_ao,
        mo_coefficients=export.mo_coefficients,
        mo_occupations=export.mo_occupations,
        mo_energies=export.mo_energies,
        overlap=None,
        ao_labels=export.ao_labels,
    )
    with pytest.raises(PortraitError, match="Next step"):
        analyze(stripped)


def test_an_out_of_range_window_is_refused():
    with pytest.raises(PortraitError, match="outside"):
        analyze(_export(), window=[999])


def test_evidence_carries_the_source():
    text = " ".join(entry.ref + " " + entry.text for entry in evidence())
    assert "2606.07879" in text
    assert "6 Angstrom" in text and "95%" in text
