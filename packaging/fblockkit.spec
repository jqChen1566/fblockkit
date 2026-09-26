# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller configuration (Windows/Linux; the green package unpacks and runs,
no Python installation needed).

Build (inside an environment with `pip install .` done):

    pyinstaller packaging/fblockkit.spec --noconfirm

Product: ``dist/fblockkit/`` (containing ``fblockkit.exe`` or ``fblockkit``).
Data files (the rule/template/menu/index YAML files plus the bibliography
sources.bib) are collected by collect_data_files through the package-data
declarations in pyproject.toml. Two lessons are baked in here:

* every new data-file suffix must be added to the ``includes`` pattern below --
  .bib once was missing and the green package crashed in menu 1 while the pip
  wheel worked (the wheel's package-data had it), so always rebuild and run the
  green package after adding data files;
* if a packaged build ever reports "rule file not found", first check whether the
  wheel really contains those YAML files.
"""

import os

from PyInstaller.utils.hooks import collect_data_files

# Resolve the script relative to the spec directory (a relative path gets the
# spec directory appended a second time by PyInstaller).
entry_script = os.path.join(SPECPATH, "entry.py")

datas = collect_data_files(
    "fblockkit", includes=["**/*.yaml", "**/*.j2", "**/*.md", "**/*.bib"]
)

a = Analysis(
    [entry_script],
    pathex=[],
    binaries=[],
    datas=datas,
    hiddenimports=[],
    hookspath=[],
    runtime_hooks=[],
    excludes=["tkinter", "unittest", "pydoc", "doctest"],
    noarchive=False,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="fblockkit",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=True,
    disable_windowed_traceback=False,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    name="fblockkit",
)
