"""docs/physics.md must agree with the code it describes.

The checks look for identifiers and numbers, not for wording, so the prose can change freely.
"""

import re
from pathlib import Path

from darkroom_atlas.manifest import _PATH_FIELDS, _VOLATILE_FIELDS, RunManifest
from darkroom_atlas.spectral import DEFAULT_GRID, cmf_xyz, illuminant

REPO = Path(__file__).resolve().parents[1]
PHYSICS = REPO / "docs" / "physics.md"


def _section(heading: str) -> str:
    """The text of the ``###`` section whose heading starts with ``heading``."""
    text = PHYSICS.read_text(encoding="utf-8")
    match = re.search(rf"^### {re.escape(heading)}.*?(?=^##)", text, flags=re.MULTILINE | re.DOTALL)
    assert match, heading
    return match.group(0)


def _identifiers(text: str) -> set[str]:
    return set(re.findall(r"`([^`\n]+)`", text))


def test_recipe_hash_exclusions_match_manifest():
    section = _section("1.5")
    named = _identifiers(section)
    for field in _VOLATILE_FIELDS:
        assert field in named, field
    for field in _PATH_FIELDS:
        assert f"{field}.path" in named, field
    assert "curve_files" in named
    assert "curve_files" in RunManifest.model_fields


def test_white_point_cause_matches_grid():
    section = _section("3.1")
    white = illuminant("D65", DEFAULT_GRID) @ cmf_xyz(DEFAULT_GRID)
    x, y = white[:2] / white.sum()
    assert f"{x:.6f}" in section and f"{y:.6f}" in section
    assert "0.3127" in section and "0.3290" in section


def test_absolute_photometry_anchors_are_named():
    # The ideal round trip is blind to absolute photometry, so section 3.4 must say which
    # tests anchor it.
    named = _identifiers(_section("3.4"))
    assert {
        "test_photon_energy_555nm",
        "test_photon_energy_across_the_grid",
        "test_luminous_exposure_anchor",
        "test_average_metering_hits_iso",
        "test_roundtrip_is_blind_to_absolute_photometry",
    } <= named


def test_named_tests_exist():
    text = PHYSICS.read_text(encoding="utf-8")
    named = {name for name in _identifiers(text) if name.startswith("test_")}
    assert named
    sources = "\n".join(p.read_text(encoding="utf-8") for p in (REPO / "tests").rglob("test_*.py"))
    for name in named:
        assert re.search(rf"^def {name}\(", sources, flags=re.MULTILINE), name
