import copy
import json

import pytest
from pydantic import ValidationError

from darkroom_atlas import ATLAS_VERSION
from darkroom_atlas.manifest import DataRef, FileRef, ReadoutRef, RunManifest
from darkroom_atlas.scene.io import EIGHT_BIT_WARNING


def _data() -> dict:
    return {
        "schema_version": 1,
        "atlas_version": ATLAS_VERSION,
        "created_utc": "2026-11-20T08:00:00Z",
        "input": {"path": "in/photo.png", "sha256": "a" * 64},
        "material": {"id": "ideal", "file": {"path": "films/ideal.yaml", "sha256": "b" * 64}},
        "process": {
            "id": "identity",
            "file": {"path": "processes/identity.yaml", "sha256": "c" * 64},
        },
        "baths": [],
        "species_file": {"path": "chemistry/species.yaml", "sha256": "d" * 64},
        "readout": {"device": "ideal_colorimetric", "params": {}},
        "exposure": {
            "illuminant": "D65",
            "iso": 100.0,
            "ev": 0.0,
            "time_s": 0.008,
            "metering": "average",
            "energy_scale_j_per_um2": 1.234e-15,
        },
        "geometry": {"shape_hw": [64, 96], "long_edge_mm": 36.0, "pixel_pitch_um": 375.0},
        "spectral_grid": {"start_nm": 380.0, "stop_nm": 780.0, "step_nm": 5.0},
        "seed": 0,
        "tile_px": None,
        "environment": {
            "python": "3.13.0",
            "numpy": "2.0.0",
            "scipy": "1.0.0",
            "colour": "0.4.7",
            "platform": "test",
        },
        "git_commit": None,
        "output": {"path": "out/photo__ideal__identity.png", "sha256": "e" * 64},
        "warnings": [EIGHT_BIT_WARNING],
    }


def _changed(data: dict, dotted: str, value) -> dict:
    data = copy.deepcopy(data)
    *parents, last = dotted.split(".")
    target = data
    for key in parents:
        target = target[key]
    target[last] = value
    return data


@pytest.fixture
def manifest() -> RunManifest:
    return RunManifest.model_validate(_data())


def test_write_read_roundtrip(tmp_path, manifest):
    path = tmp_path / "photo.run.json"
    manifest.write(path)
    back = RunManifest.read(path)
    assert back == manifest
    assert back.recipe_hash() == manifest.recipe_hash()
    assert json.loads(path.read_text(encoding="utf-8"))["material"]["id"] == "ideal"


def test_roundtrip_normalises_to_json_values(tmp_path):
    # Tuples (as in FrameGeometry.shape_hw) become lists, so that what is written and read
    # back compares equal to what was built.
    data = _changed(_data(), "geometry.shape_hw", (64, 96))
    built = RunManifest.model_validate(data)
    assert built.geometry["shape_hw"] == [64, 96]
    path = tmp_path / "m.run.json"
    built.write(path)
    assert RunManifest.read(path) == built


def test_rejects_non_finite_and_bad_hash():
    with pytest.raises(ValidationError):
        RunManifest.model_validate(_changed(_data(), "exposure.iso", float("nan")))
    with pytest.raises(ValidationError):
        RunManifest.model_validate(_changed(_data(), "input.sha256", "not-a-hash"))
    with pytest.raises(ValidationError):
        RunManifest.model_validate({**_data(), "unexpected": 1})


def test_parts_are_typed(manifest):
    assert isinstance(manifest.input, FileRef)
    assert isinstance(manifest.material, DataRef)
    assert isinstance(manifest.readout, ReadoutRef)
    assert manifest.material.file.path == "films/ideal.yaml"


def test_recipe_hash_is_hex_sha256(manifest):
    digest = manifest.recipe_hash()
    assert len(digest) == 64
    assert int(digest, 16) >= 0


def test_recipe_hash_ignores_volatile(manifest):
    volatile = [
        ("created_utc", "2030-01-01T00:00:00Z"),
        ("environment", {"python": "3.14.0", "platform": "elsewhere"}),
        ("git_commit", "0" * 40),
        ("tile_px", 32),
        ("input.path", "/somewhere/else/photo.png"),
        ("output.path", "/another/place/photo__ideal__identity.png"),
    ]
    data = _data()
    for dotted, value in volatile:
        changed = RunManifest.model_validate(_changed(_data(), dotted, value))
        assert changed != manifest, dotted
        assert changed.recipe_hash() == manifest.recipe_hash(), dotted
        data = _changed(data, dotted, value)
    all_changed = RunManifest.model_validate(data)
    assert all_changed.recipe_hash() == manifest.recipe_hash()

    for dotted, value in [("seed", 1), ("material.file.sha256", "f" * 64)]:
        changed = RunManifest.model_validate(_changed(_data(), dotted, value))
        assert changed.recipe_hash() != manifest.recipe_hash(), dotted


@pytest.mark.parametrize(
    ("dotted", "value"),
    [
        ("process.file.sha256", "f" * 64),
        ("input.sha256", "f" * 64),
        ("output.sha256", "f" * 64),
        ("exposure.ev", 1.0),
        ("readout.device", "other"),
        ("atlas_version", "D1"),
    ],
)
def test_recipe_hash_tracks_recipe(manifest, dotted, value):
    other = RunManifest.model_validate(_changed(_data(), dotted, value))
    assert other.recipe_hash() != manifest.recipe_hash()
