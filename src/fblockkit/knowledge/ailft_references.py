"""Built-in free-ion AILFT references for the nephelauxetic ratios (menu 38).

The table has two kinds of entries, each marked by its ``source``:

- **Published values** for the trivalent lanthanide and actinide series
  (Ce(3+)--Yb(3+), Th(3+)--No(3+)) transcribed from the supporting information
  of Jung, Atanasov and Neese (``jung2017covalency``; SI sections S3.5.1 and
  S3.1.1): F2/F4/F6 and zeta in cm-1, at the CASSCF and the NEVPT2 level.
  Their method: state-averaged CAS(n,7) over all SO-free states, DKH2 +
  QDPT/SOMF, SARC2 (actinide) / DKH-DEF2-TZPP (ligand) basis sets -- the same
  family as this project's recommended f-block chain.  f1 and f13 ions
  (Ce/Th, Yb/No) carry no electron-repulsion parameters in the source; only
  zeta is quoted.  One transcription note: the SI's Tb(3+) NEVPT2 row repeats
  Gd(3+)'s digits exactly (a copied row in the source); it is kept verbatim
  and flagged here rather than silently corrected.
- **A measured entry** for Ni(2+) d8 (the manual's own free-ion AILFT sample
  run, ``measured:fixtures/ailft/ni_ailft.*``), the only d-shell reference;
  its values are re-derived from the fixture by the test suite.

Cross-check on record: the Nd(3+) NEVPT2 F2 and the zeta of the two kinds of
entry agree to better than 2% -- the published SARC2/DKH2 value (77107 /
922 cm-1) against this toolkit's def2-SVP probe (75681.6 / 907.41 cm-1); the
CASSCF-level F2 differs by 7% (102546 against 94914.9), the expected
basis-set/relativistic spread for that level.  ``tests/test_ailft.py`` pins
both the Ni fixture re-derivation and the Nd cross-check.

Ratio semantics: the report resolves references **per parameter level** (a
complex run's NEVPT2 B against the free-ion NEVPT2 value), and a caller-
entered value overrides the table.  An entry that came from a probe run which
met only the energy convergence criterion carries ``quality: "energy-only"``
and the report marks its ratios provisional; published entries do not carry
the flag.
"""

from __future__ import annotations

__all__ = [
    "FREE_ION_REFERENCES",
    "covered_ions",
    "lookup",
]

#: The bibliography key of the published series (sources.bib).
JUNG2017 = "jung2017covalency"

#: ion label (element symbol + charge, e.g. "Nd3+") -> {"levels": {...}, ...}
#: levels: "casscf" / "nevpt2"; keys F2/F4/F6 (cm-1) and zeta (cm-1, one value
#: per ion, carried on the CASSCF level as in the source tables).
FREE_ION_REFERENCES: dict[str, dict[str, object]] = {
    # --- trivalent lanthanides (Jung et al., SI section S3.5.1) ---------------
    "Ce3+": {"source": JUNG2017, "levels": {"casscf": {"zeta": 678}}},
    "Pr3+": {
        "source": JUNG2017,
        "levels": {
            "casscf": {"F2": 98333, "F4": 61673, "F6": 44361, "zeta": 792},
            "nevpt2": {"F2": 72924, "F4": 56617, "F6": 38118},
        },
    },
    "Nd3+": {
        "source": JUNG2017,
        "levels": {
            "casscf": {"F2": 102546, "F4": 64345, "F6": 46292, "zeta": 922},
            "nevpt2": {"F2": 77107, "F4": 54201, "F6": 40901},
        },
    },
    "Pm3+": {
        "source": JUNG2017,
        "levels": {
            "casscf": {"F2": 106262, "F4": 66683, "F6": 47976, "zeta": 1059},
            "nevpt2": {"F2": 80257, "F4": 56609, "F6": 42627},
        },
    },
    "Sm3+": {
        "source": JUNG2017,
        "levels": {
            "casscf": {"F2": 110028, "F4": 69059, "F6": 49690, "zeta": 1209},
            "nevpt2": {"F2": 83716, "F4": 58508, "F6": 44371},
        },
    },
    "Eu3+": {
        "source": JUNG2017,
        "levels": {
            "casscf": {"F2": 113508, "F4": 71244, "F6": 51262, "zeta": 1369},
            "nevpt2": {"F2": 86964, "F4": 60418, "F6": 45744},
        },
    },
    "Gd3+": {
        "source": JUNG2017,
        "levels": {
            "casscf": {"F2": 117003, "F4": 73441, "F6": 52844, "zeta": 1543},
            "nevpt2": {"F2": 89960, "F4": 63481, "F6": 47648},
        },
    },
    # the source's Tb NEVPT2 row repeats Gd's digits exactly (a copied row);
    # kept verbatim, see the module docstring.
    "Tb3+": {
        "source": JUNG2017,
        "levels": {
            "casscf": {"F2": 120270, "F4": 75480, "F6": 54308, "zeta": 1731},
            "nevpt2": {"F2": 89960, "F4": 63481, "F6": 47648},
        },
    },
    "Dy3+": {
        "source": JUNG2017,
        "levels": {
            "casscf": {"F2": 123549, "F4": 77531, "F6": 55782, "zeta": 1932},
            "nevpt2": {"F2": 96336, "F4": 66843, "F6": 50727},
        },
    },
    "Ho3+": {
        "source": JUNG2017,
        "levels": {
            "casscf": {"F2": 126690, "F4": 79489, "F6": 57187, "zeta": 2149},
            "nevpt2": {"F2": 99273, "F4": 69703, "F6": 51879},
        },
    },
    "Er3+": {
        "source": JUNG2017,
        "levels": {
            "casscf": {"F2": 130949, "F4": 82629, "F6": 59595, "zeta": 2117},
            "nevpt2": {"F2": 103129, "F4": 72592, "F6": 55028},
        },
    },
    "Tm3+": {
        "source": JUNG2017,
        "levels": {
            "casscf": {"F2": 132883, "F4": 83349, "F6": 59958, "zeta": 2630},
            "nevpt2": {"F2": 105124, "F4": 76071, "F6": 55379},
        },
    },
    "Yb3+": {"source": JUNG2017, "levels": {"casscf": {"zeta": 2889}}},
    # --- trivalent actinides (Jung et al., SI section S3.1.1) -----------------
    "Th3+": {"source": JUNG2017, "levels": {"casscf": {"zeta": 1326}}},
    "Pa3+": {
        "source": JUNG2017,
        "levels": {
            "casscf": {"F2": 67473, "F4": 43689, "F6": 31911, "zeta": 1576},
            "nevpt2": {"F2": 48277, "F4": 40375, "F6": 26326},
        },
    },
    "U3+": {
        "source": JUNG2017,
        "levels": {
            "casscf": {"F2": 71379, "F4": 46337, "F6": 33894, "zeta": 1842},
            "nevpt2": {"F2": 52350, "F4": 37903, "F6": 30279},
        },
    },
    "Np3+": {
        "source": JUNG2017,
        "levels": {
            "casscf": {"F2": 74800, "F4": 48644, "F6": 35617, "zeta": 2116},
            "nevpt2": {"F2": 55823, "F4": 40315, "F6": 32249},
        },
    },
    "Pu3+": {
        "source": JUNG2017,
        "levels": {
            "casscf": {"F2": 78142, "F4": 50901, "F6": 37305, "zeta": 2405},
            "nevpt2": {"F2": 59418, "F4": 42389, "F6": 34087},
        },
    },
    "Am3+": {
        "source": JUNG2017,
        "levels": {
            "casscf": {"F2": 81212, "F4": 52965, "F6": 38847, "zeta": 2706},
            "nevpt2": {"F2": 62798, "F4": 44423, "F6": 35447},
        },
    },
    "Cm3+": {
        "source": JUNG2017,
        "levels": {
            "casscf": {"F2": 84227, "F4": 54994, "F6": 40363, "zeta": 3023},
            "nevpt2": {"F2": 66016, "F4": 46998, "F6": 37134},
        },
    },
    "Bk3+": {
        "source": JUNG2017,
        "levels": {
            "casscf": {"F2": 87063, "F4": 56895, "F6": 41781, "zeta": 3300},
            "nevpt2": {"F2": 69419, "F4": 48366, "F6": 38533},
        },
    },
    "Cf3+": {
        "source": JUNG2017,
        "levels": {
            "casscf": {"F2": 89861, "F4": 58772, "F6": 43182, "zeta": 3705},
            "nevpt2": {"F2": 72685, "F4": 50029, "F6": 40142},
        },
    },
    "Es3+": {
        "source": JUNG2017,
        "levels": {
            "casscf": {"F2": 92531, "F4": 60557, "F6": 44513, "zeta": 4072},
            "nevpt2": {"F2": 75685, "F4": 52869, "F6": 40974},
        },
    },
    "Fm3+": {
        "source": JUNG2017,
        "levels": {
            "casscf": {"F2": 95175, "F4": 62326, "F6": 45832, "zeta": 4458},
            "nevpt2": {"F2": 78924, "F4": 54232, "F6": 42572},
        },
    },
    "Md3+": {
        "source": JUNG2017,
        "levels": {
            "casscf": {"F2": 97720, "F4": 64023, "F6": 47096, "zeta": 4862},
            "nevpt2": {"F2": 82868, "F4": 57020, "F6": 41342},
        },
    },
    "No3+": {"source": JUNG2017, "levels": {"casscf": {"zeta": 5285}}},
    # --- the measured d8 probe (the manual's own free-ion sample run) ---------
    "Ni2+": {
        "source": "measured:fixtures/ailft/ni_ailft.*",
        "levels": {
            "casscf": {"B": 1328.1, "C": 4865.5, "zeta": 664.14, "quality": "converged"},
            "nevpt2": {"B": 1218.0, "C": 4513.5, "quality": "converged"},
        },
    },
}


def covered_ions() -> tuple[str, ...]:
    """The ion labels the built-in table covers, in table order."""
    return tuple(FREE_ION_REFERENCES)


def lookup(label: str) -> dict[str, object] | None:
    """The table entry for an ion label; case- and space-insensitive on the label."""
    wanted = label.strip().lower().replace(" ", "")
    for key, entry in FREE_ION_REFERENCES.items():
        if key.lower() == wanted:
            return entry
    return None
