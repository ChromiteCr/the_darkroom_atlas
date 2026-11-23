"""End-to-end ``bw-mechanistic`` run: expose, develop, stop, fix, wash, dry and scan (spec D1 10 D1d).

The library is a copy of the fixture library, with the engine tests' agent rates in its
kinetics file, plus the real species the fictional developer
``test-bw-mech-dev`` needs (copied from ``data/chemistry/species.yaml``; the fixture species
file holds only fictional species, and a library reads one root), that developer, the stop,
fix and wash baths and the process ``test-bw-mechanistic``. Frames are small (a 0.75 mm long
edge, so few chemistry cells), development is short and the development blocks and halo are
shrunk as in ``tests/develop/test_step.py``; the spec's 256 x 384 frame and canonical blocks
are left to the benchmark.
"""

import shutil
from pathlib import Path

import imagecodecs
import numpy as np
import pytest
import yaml

from darkroom_atlas.data.library import Library
from darkroom_atlas.develop import step as dstep
from darkroom_atlas.manifest import RunManifest
from darkroom_atlas.pipeline import RunContext, expose, process
from darkroom_atlas.runner import manifest_path, run_image
from darkroom_atlas.scene.camera import scene_light
from darkroom_atlas.scene.io import read_image

TESTS = Path(__file__).resolve().parent
REPO = TESTS.parent
FIXTURES = TESTS / "fixtures" / "library"
SAMPLE = TESTS / "assets" / "sample.png"
MATERIAL = "test/bw-mechanistic"
PROCESS = "test-bw-mechanistic"
LONG_EDGE_MM = 0.75
SMALL_BLOCK_CELLS = 8
SMALL_HALO_UM = 48.0

# The engine tests' agent rates (tests/develop/conftest.py TEST_RATES): inside the prior bounds
# of spec D1 2.4.11; the kinetics.yaml placeholders hardly develop the fixture film.
TEST_RATES = {"k_et": 0.02, "lambda0": 10.0, "K_ads": 1000.0}

REAL_SPECIES = (
    "sulfite",
    "borate",
    "hydroquinone",
    "metol",
    "bromide",
    "iodide",
    "silver",
    "sodium-sulfite",
    "borax-pentahydrate",
    "metol-sulfate",
)

FILES = {
    "developers/test-bw-mech-dev.yaml": """\
# Fictional MQ-borax developer from real species, for the bw-mechanistic end-to-end test.
schema_version: 1
id: test-bw-mech-dev
name: test-bw-mechanistic-developer
kind: developer
formula_confidence: published
components:
  - {species: sodium-sulfite, amount: {value: 100.0, unit: g/L, source: assumed}}
  - {species: hydroquinone, amount: {value: 5.0, unit: g/L, source: assumed}}
  - {species: metol-sulfate, amount: {value: 2.0, unit: g/L, source: assumed}}
  - {species: borax-pentahydrate, amount: {value: 2.0, unit: g/L, source: assumed}}
standard_dilutions: ["stock", "1+1"]
note: Not a real formula.
""",
    "baths/test-mech-stop.yaml": """\
# Fictional stop bath for the bw-mechanistic end-to-end test (spec D1 2.6).
schema_version: 1
id: test-mech-stop
name: test-stop
kind: stop
formula_confidence: published
components: []
constants:
  neutralisation_time: {value: 2.0, unit: s, source: assumed, note: "spec D1 2.6"}
note: Not a real formula.
""",
    "baths/test-mech-fix.yaml": """\
# Fictional fixer for the bw-mechanistic end-to-end test (spec D1 2.6).
schema_version: 1
id: test-mech-fix
name: test-fixer-constants
kind: fixer
formula_confidence: published
components: []
constants:
  clearing_time_constant: {value: 20.0, unit: s, source: assumed, note: "spec D1 2.6"}
  dye_removal_time: {value: 60.0, unit: s, source: assumed, note: "spec D1 2.6"}
note: Not a real formula.
""",
    "baths/test-mech-wash.yaml": """\
# Fictional wash for the bw-mechanistic end-to-end test (spec D1 2.6).
schema_version: 1
id: test-mech-wash
name: test-wash
kind: wash
formula_confidence: published
components: []
constants:
  wash_time_constant: {value: 120.0, unit: s, source: assumed, note: "spec D1 2.6"}
  wash_dye_time_constant: {value: 300.0, unit: s, source: assumed, note: "spec D1 2.6"}
note: Not a real formula.
""",
    "processes/test-bw-mechanistic.yaml": """\
# Fictional bw-mechanistic process with a short development, for the end-to-end test.
schema_version: 1
id: test-bw-mechanistic
name: test-bw-mechanistic-process
family: bw
steps:
  - kind: develop
    bath: test-bw-mech-dev
    dilution: stock
    time: "60s"
    temperature: {value: 20.0, unit: C, source: assumed}
    agitation: {method: inversion, initial_s: 15.0, interval_s: 30.0, duration_s: 5.0}
    tank: {volume_ml: 300.0, film_area_cm2: 80.0}
  - {kind: stop, bath: test-mech-stop, time: "30s"}
  - {kind: fix, bath: test-mech-fix, time: "5min"}
  - {kind: wash, bath: test-mech-wash, time: "10min"}
  - {kind: dry, time: "30min"}
note: Not a real process.
""",
}


def mechanistic_library(root: Path) -> Library:
    """Copy the fixture library to ``root`` and add the bw-mechanistic end-to-end files."""
    shutil.copytree(FIXTURES, root)
    real = yaml.safe_load((REPO / "data" / "chemistry" / "species.yaml").read_text())
    by_id = {s["id"]: s for s in real["species"]}
    species_path = root / "chemistry" / "species.yaml"
    species = yaml.safe_load(species_path.read_text())
    species["species"].extend(by_id[i] for i in REAL_SPECIES)
    species_path.write_text(yaml.safe_dump(species, sort_keys=False, allow_unicode=True))
    kinetics_path = root / "chemistry" / "kinetics.yaml"
    kinetics = yaml.safe_load(kinetics_path.read_text())
    for agent in kinetics["agents"].values():
        for name, value in TEST_RATES.items():
            agent[name]["value"] = value
            agent[name]["note"] = "test rate (tests/develop/conftest.py TEST_RATES)"
    kinetics_path.write_text(yaml.safe_dump(kinetics, sort_keys=False, allow_unicode=True))
    for rel, text in FILES.items():
        (root / rel).write_text(text)
    return Library(root)


def developed_state(library: Library, *, seed: int = 0):
    """The processed state of the sample image (small blocks and halo), for readout tests."""
    ctx = RunContext(library, seed=seed, long_edge_mm=LONG_EDGE_MM)
    image = read_image(SAMPLE)
    light = scene_light(image.linear_rgb, grid=ctx.grid, iso=100.0, low_rank=True)
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(dstep, "BLOCK_CELLS", SMALL_BLOCK_CELLS)
        mp.setattr(dstep, "halo_um", lambda *a: SMALL_HALO_UM)
        material = library.material(MATERIAL)
        state = process(expose(material, light, 1 / 125, ctx), library.process(PROCESS), ctx)
    return state, ctx


@pytest.fixture(scope="module")
def library(tmp_path_factory):
    return mechanistic_library(tmp_path_factory.mktemp("bwlib") / "library")


@pytest.fixture(autouse=True)
def small_blocks(monkeypatch):
    monkeypatch.setattr(dstep, "BLOCK_CELLS", SMALL_BLOCK_CELLS)
    monkeypatch.setattr(dstep, "halo_um", lambda *a: SMALL_HALO_UM)


def _run(library, out_dir: Path, **kw):
    kw.setdefault("workers", 1)
    return run_image(
        SAMPLE,
        out_dir,
        material_id=MATERIAL,
        process_id=PROCESS,
        long_edge_mm=LONG_EDGE_MM,
        library=library,
        **kw,
    )


@pytest.fixture(scope="module")
def baseline(library, tmp_path_factory):
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(dstep, "BLOCK_CELLS", SMALL_BLOCK_CELLS)
        mp.setattr(dstep, "halo_um", lambda *a: SMALL_HALO_UM)
        return _run(library, tmp_path_factory.mktemp("base"))


def test_library_is_valid(library):
    assert library.validate() == []


def test_run_image_bw_end_to_end(library, baseline, tmp_path):
    out_path, manifest = baseline
    pixels = imagecodecs.png_decode(out_path.read_bytes())
    assert pixels.shape == (64, 96, 3)
    assert np.array_equal(pixels[..., 0], pixels[..., 1])
    assert np.array_equal(pixels[..., 0], pixels[..., 2])
    assert pixels.min() < pixels.max()
    assert RunManifest.read(manifest_path(out_path)) == manifest
    assert manifest.readout.device == "scan_standard"
    reference = out_path.read_bytes()
    for name, kw in (
        ("rerun", {}),
        ("tiled", {"tile_px": 64}),
        ("parallel", {"workers": 2}),
    ):
        path, other = _run(library, tmp_path / name, **kw)
        assert path.read_bytes() == reference, name
        assert other.recipe_hash() == manifest.recipe_hash(), name


def test_manifest_records_model_settings(library, baseline):
    _, manifest = baseline
    settings = manifest.model_settings
    assert settings["profile"] == "image"
    assert (settings["n_z"], settings["n_l"], settings["r_l"]) == (2, 4, 2.0)
    assert settings["size_classes"] == 6
    assert settings["time_steps"] == {"h0_s": 0.1, "ratio": 1.6, "h_max_s": 5.0}
    assert settings["develop_block_cells"] == SMALL_BLOCK_CELLS
    assert settings["develop_halo_um"] == SMALL_HALO_UM
    assert manifest.film_format == "135"
    assert manifest.kinetics_file is not None
    assert manifest.kinetics_file.path == "chemistry/kinetics.yaml"
    assert manifest.kinetics_file.sha256 == Library.sha256(
        library.root / "chemistry" / "kinetics.yaml"
    )
    assert manifest.calibration_runs == []
    assert manifest.workers == 1
    changed = manifest.model_copy(update={"workers": 7})
    assert changed.recipe_hash() == manifest.recipe_hash()
