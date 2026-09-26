"""Menu loading and rendering (declarative; the data lives in ``ui/menu.yaml``).

"Numbers are append-only, never reordered" and "documentation and menu stay
isomorphic" are both carried by menu.yaml alone; this module only loads, validates
and renders, and holds no business logic.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

DEFAULT_MENU_FILE = Path(__file__).resolve().parent / "menu.yaml"


class MenuError(ValueError):
    """The menu file is invalid."""


@dataclass(frozen=True)
class MenuItem:
    number: str
    title: str
    handler: str


def load_menu(path: Path | str = DEFAULT_MENU_FILE) -> tuple[MenuItem, ...]:
    p = Path(path)
    if not p.is_file():
        raise MenuError(f"menu file does not exist: {p}")
    try:
        data = yaml.safe_load(p.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise MenuError(f"{p.name}: YAML parse failed: {exc}") from exc
    items = (data or {}).get("items")
    if not isinstance(items, list) or not items:
        raise MenuError(f"{p.name}: items is missing or empty")
    seen: set[str] = set()
    menu: list[MenuItem] = []
    for index, entry in enumerate(items):
        if not isinstance(entry, dict):
            raise MenuError(f"{p.name}: items[{index}] must be a mapping")
        number = str(entry.get("number", "")).strip()
        title = str(entry.get("title", "")).strip()
        handler = str(entry.get("handler", "")).strip()
        if not number or not title or not handler:
            raise MenuError(f"{p.name}: items[{index}] is missing number/title/handler")
        if number in seen:
            raise MenuError(
                f"{p.name}: duplicate number: {number!r} (numbers are append-only, "
                f"never reordered)"
            )
        seen.add(number)
        menu.append(MenuItem(number=number, title=title, handler=handler))
    return tuple(menu)


def render_menu(items: tuple[MenuItem, ...]) -> str:
    width = max(len(item.number) for item in items)
    lines = ["Available functions:"]
    for item in items:
        lines.append(f"  {item.number:>{width}}  {item.title}")
    return "\n".join(lines)


def find_item(items: tuple[MenuItem, ...], choice: str) -> MenuItem | None:
    for item in items:
        if item.number == choice:
            return item
    return None
