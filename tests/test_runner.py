import shutil
from pathlib import Path

import numpy as np
import pytest

from darkroom_atlas import ATLAS_VERSION
from darkroom_atlas.data.library import ENV_DATA_ROOT, DataNotFound, Library
from darkroom_atlas.manifest import RunManifest
from darkroom_atlas.registry import RegistryError
from darkroom_atlas.runner import manifest_path, run_image
from darkroom_atlas.scene.io import EIGHT_BIT_WARNING, read_image

TESTS = Path(__file__).resolve().parent
SAMPLE = TESTS / "assets" / "sample.png"
FIXTURES = TESTS / "fixtures" / "library"
REPO_DATA = TESTS.parent / "data"
OUTPUT_NAME = "sample__ideal__identity.png"


@pytest.fixture(scope="module")
def lib():
    return Library(REPO_DATA)


def _run(out_dir, lib, **kwargs):
    kwargs.setdefault("material_id", "ideal")
    kwargs.setdefault("process_id", "identity")
    return run_image(SAMPLE, out_dir, library=lib, **kwargs)


def _copy_with_material(tmp_path, source_lib: Path, rel: str, text: str) -> Library:
    """Copy ``source_lib`` to ``tmp_path`` and add the material file ``films/<rel>``."""
    root = tmp_path / "library"
    shutil.copytree(source_lib, root)
    target = root / "films" / rel
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(text, encoding="utf-8")
    return Library(root)


def _ideal_yaml(replacements: dict[str, str]) -> str:
    """``data/films/ideal.yaml`` with each ``old`` replaced by ``new``."""
    text = (REPO_DATA / "films" / "ideal.yaml").read_text(encoding="utf-8")
    for old, new in replacements.items():
        assert old in text
        text = text.replace(old, new)
    return text


def test_run_creates_outputs(tmp_path, monkeypatch):
    monkeypatch.delenv(ENV_DATA_ROOT, raising=False)
    out_dir = tmp_path / "out" / "nested"
    assert not out_dir.exists()
    out, manifest = run_image(SAMPLE, out_dir, material_id="ideal", process_id="identity")
    assert out == out_dir / OUTPUT_NAME
    assert out.is_file()
    assert manifest_path(out) == out_dir / "sample__ideal__identity.run.json"
    assert manifest_path(out).is_file()
    assert isinstance(manifest, RunManifest)


def test_output_roundtrip(tmp_path, lib):
    out, _ = _run(tmp_path, lib)
    written = read_image(out)
    assert written.bit_depth == 16
    error = float(np.max(np.abs(written.linear_rgb - read_image(SAMPLE).linear_rgb)))
    print(f"run_image round trip through 16-bit PNG: max abs error {error:.3e} (limit 5e-3)")
    assert error < 5e-3


def test_rerun_bit_identical(tmp_path, lib):
    out_a, manifest_a = _run(tmp_path / "a", lib, seed=3)
    out_b, manifest_b = _run(tmp_path / "b", lib, seed=3)
    assert out_a.read_bytes() == out_b.read_bytes()
    assert manifest_a.recipe_hash() == manifest_b.recipe_hash()


def test_tiled_bit_identical(tmp_path, lib):
    out_a, manifest_a = _run(tmp_path / "untiled", lib)
    out_b, manifest_b = _run(tmp_path / "tiled", lib, tile_px=32)
    assert manifest_b.tile_px == 32 and manifest_a.tile_px is None
    assert out_a.read_bytes() == out_b.read_bytes()
    assert manifest_a.recipe_hash() == manifest_b.recipe_hash()


def test_seed_changes_recipe_hash(tmp_path, lib):
    _, manifest_a = _run(tmp_path / "a", lib, seed=0)
    _, manifest_b = _run(tmp_path / "b", lib, seed=1)
    assert manifest_a.recipe_hash() != manifest_b.recipe_hash()


def test_manifest_contents(tmp_path, lib):
    out, manifest = _run(tmp_path, lib)
    assert manifest.atlas_version == ATLAS_VERSION
    assert manifest.material.id == "ideal"
    assert manifest.material.file.path == "films/ideal.yaml"
    assert manifest.material.file.sha256 == Library.sha256(lib.path_of("material", "ideal"))
    assert manifest.process.id == "identity"
    assert manifest.process.file.sha256 == Library.sha256(lib.path_of("process", "identity"))
    assert manifest.baths == []
    assert manifest.species_file is not None
    assert manifest.species_file.sha256 == Library.sha256(REPO_DATA / "chemistry/species.yaml")
    assert manifest.readout.device == "ideal_colorimetric"
    assert manifest.exposure["iso"] == 100
    assert manifest.exposure["illuminant"] == "D65"
    assert manifest.geometry["shape_hw"] == [64, 96]
    assert manifest.geometry["pixel_pitch_um"] == pytest.approx(36000.0 / 96)
    assert manifest.spectral_grid == {"start_nm": 380.0, "stop_nm": 780.0, "step_nm": 5.0}
    assert manifest.seed == 0
    assert EIGHT_BIT_WARNING in manifest.warnings
    assert set(manifest.environment) == {"python", "numpy", "scipy", "colour", "platform"}
    assert manifest.input.sha256 == Library.sha256(SAMPLE)
    assert manifest.output.sha256 == Library.sha256(out)
    # What was written is what was returned, and it survives the JSON round trip.
    assert RunManifest.read(manifest_path(out)) == manifest


def test_slash_ids_in_filename(tmp_path):
    material = _ideal_yaml({"id: ideal\n": "id: test/ideal\n"})
    lib = _copy_with_material(tmp_path, FIXTURES, "test/ideal.yaml", material)
    process_dir = lib.root / "processes" / "test"
    process_dir.mkdir(parents=True)
    identity = (REPO_DATA / "processes" / "identity.yaml").read_text(encoding="utf-8")
    (process_dir / "identity.yaml").write_text(
        identity.replace("id: identity\n", "id: test/identity\n"), encoding="utf-8"
    )
    out, manifest = run_image(
        SAMPLE, tmp_path / "out", material_id="test/ideal", process_id="test/identity",
        library=lib,
    )  # fmt: skip
    assert out.name == "sample__test-ideal__test-identity.png"
    assert "test-ideal" in out.name
    assert out.parent == tmp_path / "out"
    assert manifest_path(out).is_file()
    assert manifest.material.id == "test/ideal"
    assert manifest.material.file.path == "films/test/ideal.yaml"


def test_no_default_readout_errors(tmp_path):
    material = _ideal_yaml(
        {"id: ideal\n": "id: no-readout\n", "default_readout: ideal_colorimetric\n": ""}
    )
    lib = _copy_with_material(tmp_path, REPO_DATA, "no-readout.yaml", material)
    out_dir = tmp_path / "out"
    with pytest.raises(ValueError, match="default_readout"):
        run_image(SAMPLE, out_dir, material_id="no-readout", process_id="identity", library=lib)
    assert not out_dir.exists()
    # Naming the device explicitly works.
    out, manifest = run_image(
        SAMPLE, out_dir, material_id="no-readout", process_id="identity",
        readout_device="ideal_colorimetric", library=lib,
    )  # fmt: skip
    assert out.is_file()
    assert manifest.readout.device == "ideal_colorimetric"


def test_iso_from_material(tmp_path):
    iso_line = 'iso_speed: {value: 400, unit: "1", source: assumed}\n'
    material = _ideal_yaml({"id: ideal\n": "id: iso400\n" + iso_line})
    lib = _copy_with_material(tmp_path, REPO_DATA, "iso400.yaml", material)
    _, from_material = run_image(
        SAMPLE, tmp_path / "a", material_id="iso400", process_id="identity", library=lib
    )
    assert from_material.exposure["iso"] == 400
    _, explicit = run_image(
        SAMPLE, tmp_path / "b", material_id="iso400", process_id="identity", iso=50.0,
        library=lib,
    )  # fmt: skip
    assert explicit.exposure["iso"] == 50


@pytest.mark.parametrize(
    ("kwargs", "error", "match"),
    [
        ({"material_id": "nope"}, DataNotFound, "ideal"),
        ({"process_id": "nope"}, DataNotFound, "identity"),
        ({"readout_device": "nope"}, RegistryError, "ideal_colorimetric"),
        ({"input_path": TESTS / "assets" / "missing.png"}, FileNotFoundError, "missing.png"),
        ({"illuminant": "X"}, ValueError, "D65"),
        ({"tile_px": 0}, ValueError, "tile_px"),
    ],
)
def test_errors_leave_no_output(tmp_path, lib, kwargs, error, match):
    args = {"material_id": "ideal", "process_id": "identity", "library": lib, **kwargs}
    input_path = args.pop("input_path", SAMPLE)
    out_dir = tmp_path / "out"
    with pytest.raises(error, match=match):
        run_image(input_path, out_dir, **args)
    assert not out_dir.exists()


def test_output_dir_is_a_file(tmp_path, lib):
    not_a_dir = tmp_path / "taken"
    not_a_dir.write_text("occupied", encoding="utf-8")
    with pytest.raises(NotADirectoryError, match="taken"):
        _run(not_a_dir, lib)
    with pytest.raises(NotADirectoryError, match="taken"):
        _run(not_a_dir / "below", lib)
    assert not_a_dir.read_text(encoding="utf-8") == "occupied"
