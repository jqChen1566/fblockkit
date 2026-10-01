"""Audit every session-capture reference: title vs the content it currently prints.

Prints, per reference: chapter, current window, a title snippet, and the first/last
printed lines of the window.  Used to judge each reference by its title's promise.
"""
import re
from pathlib import Path

BASE = Path(__file__).resolve().parent
PATTERN = re.compile(r"\\lstinputlisting\[(?P<opts>[^\]]*)\]\s*\{(?P<file>[^}]*)\}", re.S)

rows = []
for tex in sorted((BASE / "chapters").glob("*.tex")):
    text = tex.read_text(encoding="utf-8")
    for match in PATTERN.finditer(text):
        file = match.group("file").strip()
        if not file.startswith("examples/expected/m"):
            continue
        first = re.search(r"firstline=(\d+)", match.group("opts"))
        last = re.search(r"lastline=(\d+)", match.group("opts"))
        if not first or not last:
            continue
        # title snippet: from the first 'title=' after the listing start, 90 chars
        tail = text[match.end() : match.end() + 400]
        title = ""
        tm = re.search(r"title=\{(.{0,110})", tail, re.S)
        if tm:
            title = " ".join(tm.group(1).split())
        f, l = int(first.group(1)), int(last.group(1))
        lines = (BASE / file).read_text(encoding="utf-8").splitlines()
        first_line = lines[f - 1].strip()[:88] if f - 1 < len(lines) else "<oob>"
        last_line = lines[l - 1].strip()[:88] if l - 1 < len(lines) else "<oob>"
        rows.append((tex.name, f, l, title, first_line, last_line))

print(f"references: {len(rows)}")
for name, f, l, title, first_line, last_line in rows:
    print(f"\n## {name}  [{f}-{l}]")
    print(f"   title: {title}")
    print(f"   first: {first_line}")
    print(f"   last : {last_line}")
