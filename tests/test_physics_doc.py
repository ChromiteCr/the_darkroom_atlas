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


# --- D1 chapters (section 4) ---------------------------------------------------------------


def _d1_text() -> str:
    text = PHYSICS.read_text(encoding="utf-8")
    start = text.index("## 4 D1 已实现的物理")
    return text[start : text.index("## 5 已知的不物理之处")]


def _doc_section(heading: str) -> str:
    text = _d1_text()
    match = re.search(
        rf"^### {re.escape(heading)}.*?(?=^### |\Z)", text, flags=re.MULTILINE | re.DOTALL
    )
    assert match, heading
    return match.group(0)


def test_d1_chapters_name_the_public_functions():
    from importlib import import_module

    expected = {
        "4.1": [("emulsion.classes", "size_classes"), ("emulsion.classes", "SizeClasses")],
        "4.2": [
            ("optics.stack", "optical_stack"),
            ("optics.kubelka", "depth_response"),
            ("optics.halation", "halation_kernel"),
            ("optics.halation", "halation_fraction"),
            ("optics.psf", "scatter_kernel"),
            ("optics.expose", "exposure_fields"),
            ("optics.expose", "uniform_exposure"),
            ("models.bw", "expose_bw"),
        ],
        "4.3": [
            ("latent.chain", "ChainParams"),
            ("latent.chain", "cluster_distribution"),
            ("latent.table", "latent_table"),
            ("latent.table", "LatentTable"),
        ],
        "4.4": [
            ("chemistry.species", "component_table"),
            ("chemistry.speciation", "make_speciation"),
            ("chemistry.speciation", "solve_ph"),
            ("chemistry.solution", "working_solution"),
            ("chemistry.solution", "component_charges"),
        ],
        "4.5": [
            ("light", "LowRankLight"),
            ("fieldstore", "FieldStore"),
            ("state", "FrameAssembler"),
            ("blocks", "block_grid"),
            ("blocks", "run_blocks"),
            ("pipeline", "segments"),
        ],
    }
    for heading, items in expected.items():
        named = _identifiers(_doc_section(heading))
        for module, name in items:
            assert hasattr(import_module(f"darkroom_atlas.{module}"), name), (module, name)
            assert any(name in item for item in named), (heading, name)


def test_d1_constants_match_code():
    from darkroom_atlas.chemistry import activity, speciation
    from darkroom_atlas.chemistry.solution import IONIC_STRENGTH_TOL
    from darkroom_atlas.chemistry.species import (
        _MANOV_A,
        _MANOV_B,
        _MANOV_C,
        DEFAULT_GELATIN_PKA,
    )
    from darkroom_atlas.emulsion import classes
    from darkroom_atlas.models.bw import EXPOSURE_TEMPERATURE_C, PROFILE_CLASSES
    from darkroom_atlas.optics.expose import BLOCK_PX, PSF_TAIL
    from darkroom_atlas.optics.stack import PHOTON_ENERGY_EV_NM

    sizes = _doc_section("4.1")
    for value in (
        classes.M_AG,
        classes.M_AGBR,
        classes.M_AGI,
        classes.RHO_AGX,
        classes.RHO_GELATIN,
    ):
        assert f"{value:g}" in sizes, value
    for profile, k in PROFILE_CLASSES.items():
        assert f"`{profile}` {k}" in sizes, profile

    optics = _doc_section("4.2")
    assert f"`BLOCK_PX = {BLOCK_PX}`" in optics
    assert f"`PSF_TAIL = {PSF_TAIL:.0e}`".replace("e-04", "e-4") in optics
    assert f"{PHOTON_ENERGY_EV_NM}" in optics
    assert f"（{EXPOSURE_TEMPERATURE_C}）" in optics

    chemistry = _doc_section("4.4")
    assert f"{activity.PKW_20C}" in chemistry and f"{activity.PKW_25C:.2f}" in chemistry
    for value in (_MANOV_A, _MANOV_B):
        assert f"{value:g}" in chemistry, value
    assert f"{abs(_MANOV_C)}" in chemistry
    assert f"{DEFAULT_GELATIN_PKA}" in chemistry
    assert "1e-6 M" in chemistry and IONIC_STRENGTH_TOL == 1e-6
    assert f"[{speciation.PH_MIN:g}, {speciation.PH_MAX}]" in chemistry


def test_d1_exposure_state_keys_are_the_ones_the_code_writes():
    import inspect

    from darkroom_atlas.models import bw
    from darkroom_atlas.optics import expose

    source = inspect.getsource(bw) + inspect.getsource(expose)
    section = _doc_section("4.2")
    fields = (
        "emulsion/exposure_volume",
        "emulsion/exposure_surface",
        "emulsion/residual_dye",
    )
    for key in fields:
        assert f"`{key}`" in section
        assert f'"{key}"' in source
    meta_keys = (
        "time_s",
        "temperature_c",
        "n_z",
        "profile",
        "film_format",
        "halation_fraction",
        "halation_radius_um",
        "halation_truncated",
        "psf_radius_um",
        "block_px",
        "halo_px",
    )
    for key in meta_keys:
        assert f"`{key}`" in section
        assert f'"{key}"' in source


def test_d1_latent_defaults_match_code():
    import inspect

    from darkroom_atlas.latent.chain import ChainParams, cluster_distribution
    from darkroom_atlas.latent.table import latent_table

    section = _doc_section("4.3")
    table = inspect.signature(latent_table).parameters
    assert table["nodes"].default == 201 and "201" in section
    assert table["log10_range"].default == (-6.0, 4.0) and "[−6, 4]" in section
    assert (
        table["quadrature"].default
        == 9
        == (inspect.signature(cluster_distribution).parameters["quadrature"].default)
    )
    assert "缺省 9" in section
    assert ChainParams.__dataclass_fields__["cluster_cap"].default == 5
    assert "n_cap = 5" in section
