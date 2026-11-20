"""README.md must show the same version as ``darkroom_atlas.ATLAS_VERSION``."""

from pathlib import Path

from darkroom_atlas import ATLAS_VERSION

README = Path(__file__).resolve().parents[1] / "README.md"
TABLE_HEADER = "| 版本 | 日期 | 变更内容 | 类型 |"


def _lines() -> list[str]:
    return README.read_text(encoding="utf-8").splitlines()


def test_readme_badge_matches_version():
    badge = f"https://img.shields.io/badge/version-{ATLAS_VERSION}-blue"
    assert badge in README.read_text(encoding="utf-8")


def test_badge_is_first_line_after_title():
    lines = _lines()
    assert lines[0].startswith("# ")
    first = next(line for line in lines[1:] if line.strip())
    assert first.startswith("![version]")


def test_top_table_row_matches_version():
    lines = _lines()
    start = lines.index(TABLE_HEADER)
    separator = lines[start + 1]
    assert set(separator) <= set("|-: ")
    assert "-" in separator
    top_row = lines[start + 2]
    cells = [cell.strip() for cell in top_row.strip().strip("|").split("|")]
    assert cells[0] == ATLAS_VERSION
