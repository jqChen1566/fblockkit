# fixtures/basisdb -- the miniature basis-library index (menu 4's datasource leg)

`mini_basis.db` is a four-row SQLite miniature baked from the real deployed
library (101: `~/projects/orca_basis_sets`, the basisdb project's BSE v0.12
snapshot with 776 sets; measured 2026-09-29).  It exists so the query layer
(`knowledge/basisdb.py`) and the examples run byte-identically on both
ends without bundling any part of the 776-set library into the package.

Baked rows (verbatim copies of the real index rows, minus the unused
columns): `def2-svp`, `def2-tzvp` (both cover Ce; ahlrichs family, ORCA
built-in, ECP pairing, 1058/1373 contracted functions, 30/31 refs) and
`3-21g` (pople family, does not cover Ce).

Regeneration (on 101, against the deployed library):

    python3 - <<'PYEOF'
    import sqlite3
    src = sqlite3.connect('/home/chen_jianqi/projects/orca_basis_sets/basis.db')
    mini = sqlite3.connect('mini_basis.db')
    mini.execute("DROP TABLE IF EXISTS basis")
    mini.execute("""CREATE TABLE basis (bse_name TEXT, name_upper TEXT, role TEXT,
      family TEXT, relativistic_method TEXT, orca_builtin INTEGER DEFAULT 0,
      has_ecp INTEGER DEFAULT 0, total_contracted_functions INTEGER,
      ref_count INTEGER, elements_json TEXT)""")
    for name in ('DEF2-SVP', 'DEF2-TZVP', 'SARC2-DKH-QZVP', '3-21G'):
        rows = src.execute(
            "select bse_name, name_upper, role, family, relativistic_method,"
            " orca_builtin, has_ecp, total_contracted_functions, ref_count,"
            " elements_json from basis where name_upper=? and role='orbital'",
            (name,)).fetchall()
        for r in rows:
            mini.execute("insert into basis values (?,?,?,?,?,?,?,?,?,?)", r)
    mini.commit()
    PYEOF

(The `SARC2-DKH-QZVP` name matches nothing under `role='orbital'` in the
current library -- recorded, the snapshot simply bakes three rows.)

The query layer never bundles library content: it opens the user's own
deployment read-only (`$FBK_BASISDB`, then `~/projects/orca_basis_sets`).
