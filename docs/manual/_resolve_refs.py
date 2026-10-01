"""Resolve every session-capture reference by content, using the chapter's own history.

For each lstinputlisting with a line range into examples/expected/m*.txt:
  1. git blame the listing line -> the commit that wrote it;
  2. take the capture file AS OF that commit (the era when the window was written);
  3. extract the era's window [first,last];
  4. slide that window over the CURRENT capture to find its present location;
  5. report: chapter, era commit, old range, found range(s) or NO MATCH.

Dry run: prints a resolution table.  Pass --apply to rewrite the line numbers when
the match is unique (and report the rest for human handling).
"""
import re
import subprocess
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent
REPO = BASE.parent.parent
PATTERN = re.compile(r"\\lstinputlisting\[(?P<opts>[^\]]*)\]\s*\{(?P<file>[^}]*)\}", re.S)


def git(*args: str) -> str:
    out = subprocess.run(
        ["git", *args], cwd=REPO, capture_output=True, text=True, encoding="utf-8"
    )
    return out.stdout if out.returncode == 0 else ""


apply = "--apply" in sys.argv

for tex in sorted((BASE / "chapters").glob("*.tex")):
    text = tex.read_text(encoding="utf-8")
    fixes: list[tuple[int, int, str, str]] = []  # (match_start, opts_end, old, new)
    out_lines: list[str] = []
    for match in PATTERN.finditer(text):
        file = match.group("file").strip()
        if not file.startswith("examples/expected/m"):
            continue
        first = re.search(r"firstline=(\d+)", match.group("opts"))
        last = re.search(r"lastline=(\d+)", match.group("opts"))
        if not first or not last:
            continue
        f, l = int(first.group(1)), int(last.group(1))
        listing_line = text[: match.start()].count("\n") + 1
        blame = git(
            "blame", "HEAD", "-L", f"{listing_line},{listing_line}", "--porcelain",
            f"docs/manual/chapters/{tex.name}",
        )
        sha = blame.split(" ", 1)[0].lstrip("\ufeff") if blame else "?"
        era_text = git("show", f"{sha}:docs/manual/{file}") if sha != "?" else ""
        era_lines = era_text.splitlines()
        current_lines = (BASE / file).read_text(encoding="utf-8").splitlines()
        if not era_lines:
            out_lines.append(f"?? {tex.name} [{f}-{l}] {file}: no era file at {sha[:8]}")
            continue
        window = [x.rstrip() for x in era_lines[f - 1 : l]]
        hay = [x.rstrip() for x in current_lines]
        hits = [
            start + 1
            for start in range(len(hay) - len(window) + 1)
            if hay[start : start + len(window)] == window
        ]
        if len(hits) == 1:
            nf, nl = hits[0], hits[0] + (l - f)
            status = f"-> [{nf}-{nl}]"
            if (nf, nl) != (f, l):
                old_opts = match.group("opts")
                new_opts = old_opts.replace(
                    f"firstline={f}", f"firstline={nf}"
                ).replace(f"lastline={l}", f"lastline={nl}")
                fixes.append((match.start(), match.end("opts"), old_opts, new_opts))
        elif not hits:
            status = "NO MATCH"
        else:
            status = f"multi {hits[:5]}"
        out_lines.append(f"{tex.name:34s} [{f:>3}-{l:<3}] {status:24s} sha={sha[:8]}")
    if apply and fixes:
        for start, opts_end, old_opts, new_opts in reversed(fixes):
            text = (
                text[:start]
                + "\\lstinputlisting["
                + new_opts
                + "]"
                + text[opts_end + 1 :]
            )
        tex.write_text(text, encoding="utf-8")
    for line in out_lines:
        print(line)
