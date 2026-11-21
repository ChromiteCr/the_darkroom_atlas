"""The bw-mechanistic exposure model (spec §2.0.3, §2.1.9, §3.6, D1a done criterion)."""

import copy
import logging
import math
import tracemalloc
from pathlib import Path

import numpy as np
import pytest
import yaml

from darkroom_atlas.data.library import DataError, Library
from darkroom_atlas.data.models import Material
from darkroom_atlas.emulsion import size_classes
from darkroom_atlas.light import LowRankLight
from darkroom_atlas.models import bw
from darkroom_atlas.optics.expose import uniform_exposure
from darkroom_atlas.optics.kubelka import depth_response
from darkroom_atlas.optics.stack import optical_stack
from darkroom_atlas.pipeline import RunContext, expose
from darkroom_atlas.registry import get_exposure
from darkroom_atlas.spectral import DEFAULT_GRID

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "library"
FILM = "test/bw-mechanistic"
TIME_S = 0.01
B = DEFAULT_GRID.n_bands
META_KEYS = {
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
}

log = logging.getLogger(__name__)


@pytest.fixture(scope="module")
def lib():
    return Library(FIXTURES)


@pytest.fixture(scope="module")
def film(lib):
    return lib.material(FILM)


def _light(shape, seed=0, rank=3):
    rng = np.random.default_rng(seed)
    coefficients = rng.uniform(0.1, 1.0, (*shape, rank))
    basis = rng.uniform(0.2, 1.0, (rank, B)) * 50.0
    return LowRankLight(coefficients, basis, DEFAULT_GRID)


def test_registered_as_exposure_model():
    assert bw.MODEL == "bw-mechanistic"
    assert get_exposure("bw-mechanistic") is bw.expose_bw


def test_expose_bw_fields_and_meta(lib, film):
    light = _light((30, 40))
    state = expose(film, light, TIME_S, RunContext(lib, long_edge_mm=0.8))
    assert set(state.fields) == {
        "emulsion/exposure_volume",
        "emulsion/exposure_surface",
        "emulsion/residual_dye",
    }
    for key in ("emulsion/exposure_volume", "emulsion/exposure_surface"):
        field = state.fields[key]
        assert field.shape == (30, 40, 2)
        assert field.dtype == np.float64
        assert np.all(field > 0.0)
    dye = state.fields["emulsion/residual_dye"]
    assert dye.shape == (30, 40)
    assert dye.strides == (0, 0)  # a broadcast view, zero bytes (spec §4.3)
    assert np.all(dye == 1.0)
    meta = state.meta["exposure"]
    assert set(meta) == META_KEYS
    assert meta["time_s"] == TIME_S
    assert meta["temperature_c"] == 23.0
    assert meta["n_z"] == 2
    assert meta["profile"] == "image"
    assert meta["film_format"] == "135"
    assert meta["block_px"] == 512
    assert meta["halation_fraction"] == pytest.approx(0.224, abs=0.002)
    assert state.meta["light"] == light.meta


@pytest.mark.parametrize(("profile", "n_z"), [("image", 2), ("calibration", 4), ("reference", 8)])
def test_profile_sets_depth_cells(lib, film, profile, n_z):
    state = expose(film, _light((8, 9)), TIME_S, RunContext(lib, long_edge_mm=0.5, profile=profile))
    assert state.fields["emulsion/exposure_volume"].shape == (8, 9, n_z)
    assert state.meta["exposure"]["n_z"] == n_z
    assert state.meta["exposure"]["profile"] == profile


def test_film_format_selects_base(lib, film):
    light = _light((8, 9))
    s135 = expose(film, light, TIME_S, RunContext(lib, long_edge_mm=0.5))
    s120 = expose(film, light, TIME_S, RunContext(lib, long_edge_mm=0.5, film_format="120"))
    assert s120.meta["exposure"]["film_format"] == "120"
    # The 120 backing (D = 0.5) absorbs most of the returning light.
    assert (
        s120.meta["exposure"]["halation_fraction"]
        < 0.2 * s135.meta["exposure"]["halation_fraction"]
    )
    with pytest.raises(ValueError, match=r"unknown film format '220'.*135, 120"):
        expose(film, light, TIME_S, RunContext(lib, long_edge_mm=0.5, film_format="220"))


def test_expose_uniform_frame_matches_uniform_exposure(lib, film):
    # A uniform frame through the model: interior equals the one-column path up to the
    # truncated halation tail (<= 1e-3 relative, spec §2.1.5).
    shape = (140, 120)
    basis = _light((1, 1)).basis
    light = LowRankLight(np.ones((*shape, 3)), basis, DEFAULT_GRID)
    ctx = RunContext(lib, long_edge_mm=5.6)  # 40 um pixels
    state = expose(film, light, TIME_S, ctx)
    h = state.meta["exposure"]["halo_px"]
    assert 2 * h + 2 < min(shape)
    classes = size_classes(film.emulsion_layers[0], 6)
    stack = optical_stack(film, "135", classes, DEFAULT_GRID)
    v, s = uniform_exposure(basis.sum(axis=0) * TIME_S, stack, depth_response(stack, 2))
    centre = state.fields["emulsion/exposure_volume"][h + 1 : -h - 1, h + 1 : -h - 1]
    surface = state.fields["emulsion/exposure_surface"][h + 1 : -h - 1, h + 1 : -h - 1]
    assert np.max(np.abs(centre / v - 1.0)) <= 1e-3
    assert np.max(np.abs(surface / s - 1.0)) <= 1e-3


@pytest.mark.parametrize("tile_px", [None, 37, 256])
@pytest.mark.parametrize("workers", [1, 3])
def test_expose_independent_of_tiling_and_workers(lib, film, tile_px, workers, reference_fields):
    ctx = RunContext(lib, long_edge_mm=55.0, tile_px=tile_px, workers=workers)
    state = expose(film, _tiling_light(), TIME_S, ctx)
    for key, ref in reference_fields.items():
        assert np.array_equal(state.fields[key], ref), key


def _tiling_light():
    # 1100 px across at 50 um: three canonical 512-px blocks.
    return _light((40, 1100), seed=11)


@pytest.fixture(scope="module")
def reference_fields(lib, film):
    state = expose(film, _tiling_light(), TIME_S, RunContext(lib, long_edge_mm=55.0))
    assert state.geometry.pixel_pitch_um == pytest.approx(50.0)
    return {
        key: np.array(state.fields[key])
        for key in ("emulsion/exposure_volume", "emulsion/exposure_surface")
    }


def test_expose_bw_reports_missing_fields(lib):
    data = yaml.safe_load((FIXTURES / "films" / "test" / "bw-mechanistic.yaml").read_text())
    data = copy.deepcopy(data)
    del data["layers"][1]["latent"]
    del data["layers"][2]["refractive_index"]
    material = Material.model_validate(data)
    with pytest.raises(DataError) as info:
        expose(material, _light((4, 4)), TIME_S, RunContext(lib))
    message = str(info.value)
    assert "bw-mechanistic.yaml: layers[1].latent: required" in message
    assert "layers[2].refractive_index: required" in message


def test_expose_bw_small_fields_stay_in_ram(lib, film, tmp_path):
    # With a scratch directory the exposure fields live in read-only memory maps once they
    # exceed the store's threshold; small ones stay in RAM.
    state = expose(film, _light((16, 16)), TIME_S, RunContext(lib, scratch_dir=tmp_path))
    assert not isinstance(state.fields["emulsion/exposure_volume"], np.memmap)


@pytest.mark.slow
def test_expose_24mp_within_memory_budget(lib, film, tmp_path, record_property):
    # D1a done criterion: a 24 MP LowRankLight exposed with the fixture film inside 8 GB,
    # peak memory = tracemalloc peak + memory-mapped field files.
    tracemalloc.start()
    try:
        h, w = 4000, 6000
        coefficients = np.empty((h, w, 3))
        ramp = np.linspace(0.05, 1.0, w)
        for r in range(3):
            coefficients[..., r] = (
                ramp[None, :] * (0.5 + 0.5 * np.cos(np.arange(h) / (300 + r)))[:, None]
            )
        coefficients.flags.writeable = False
        basis = np.random.default_rng(0).uniform(0.2, 1.0, (3, B)) * 50.0
        light = LowRankLight(coefficients, basis, DEFAULT_GRID)
        assert light.nbytes <= 0.6e9
        state = expose(film, light, TIME_S, RunContext(lib, scratch_dir=tmp_path))
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    volume = state.fields["emulsion/exposure_volume"]
    surface = state.fields["emulsion/exposure_surface"]
    for field in (volume, surface):
        assert isinstance(field, np.memmap) and field.mode == "r"
        assert field.shape == (h, w, 2)
    mapped = sum(p.stat().st_size for p in tmp_path.rglob("*") if p.is_file())
    assert mapped == 2 * h * w * 2 * 8
    total = peak + mapped
    log.warning(
        "24 MP bw-mechanistic exposure: tracemalloc peak %.2f GB + mapped %.2f GB = %.2f GB",
        peak / 1e9,
        mapped / 1e9,
        total / 1e9,
    )
    record_property("peak_traced_bytes", peak)
    record_property("mapped_bytes", mapped)
    assert total <= 8e9
    assert np.all(np.isfinite(volume[::997, ::997])) and np.all(volume[::997, ::997] > 0.0)
    assert math.isclose(state.geometry.pixel_pitch_um, 6.0)
