"""Densitometer readouts of a developed bw-mechanistic film (spec D1 2.8)."""

import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import pytest

from darkroom_atlas.imaging.silver import silver_density
from darkroom_atlas.pipeline import readout, unexposed_state
from darkroom_atlas.readout import densitometer as dens
from darkroom_atlas.readout.products import visual_product
from darkroom_atlas.scene.camera import scene_light


def _load_e2e():
    name = "d1d_end_to_end_bw"
    if name not in sys.modules:
        path = Path(__file__).resolve().parents[1] / "test_end_to_end_bw.py"
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


_e2e = _load_e2e()


@pytest.fixture(scope="module")
def developed(tmp_path_factory):
    library = _e2e.mechanistic_library(tmp_path_factory.mktemp("lib") / "library")
    return _e2e.developed_state(library)


def test_neutral_film_equal_densities(developed):
    state, ctx = developed
    visual = readout(state, "densitometer_visual", ctx)
    type2 = readout(state, "densitometer_type2", ctx)
    assert visual.kind == type2.kind == "density"
    assert visual.data.shape == state.geometry.shape_hw
    np.testing.assert_allclose(type2.data, visual.data, rtol=0.0, atol=1e-12)
    assert np.ptp(visual.data) > 0.1


def test_silver_density_readout(developed):
    state, ctx = developed
    result = readout(state, "silver_density", ctx)
    assert np.array_equal(result.data, silver_density(state))


def test_visual_density_includes_base(developed):
    state, ctx = developed
    d = readout(state, "densitometer_visual", ctx).data
    assert np.all(d >= silver_density(state) + 0.10 - 1e-12)  # 135 base grey density 0.10


def test_aperture_average_on_transmittance(developed):
    state, ctx = developed
    pitch = state.geometry.pixel_pitch_um
    aperture = 5.0 * pitch
    pixel = readout(state, "densitometer_visual", ctx).data
    result = readout(state, "densitometer_visual", ctx, aperture_um=aperture)
    kernel = dens.aperture_kernel(aperture, pitch)
    assert kernel.shape == (5, 5) and kernel.sum() == 21
    r, c, n = 30, 40, 2
    window = pixel[r - n : r + n + 1, c - n : c + n + 1]
    t_mean = np.sum(10.0**-window * kernel) / kernel.sum()
    assert result.data[r, c] == pytest.approx(-math.log10(t_mean), rel=1e-12)
    d_mean = np.sum(window * kernel) / kernel.sum()
    assert result.data[r, c] < d_mean  # Jensen: averaging T first gives the lower density
    assert result.meta["aperture_um"] == aperture
    # At the edge only the in-frame pixels count.
    t_corner = np.sum(10.0 ** -pixel[:3, :3] * kernel[2:, 2:]) / kernel[2:, 2:].sum()
    assert result.data[0, 0] == pytest.approx(-math.log10(t_corner), rel=1e-12)


def test_small_aperture_is_per_pixel(developed):
    state, ctx = developed
    pixel = readout(state, "densitometer_visual", ctx).data
    small = readout(state, "densitometer_visual", ctx, aperture_um=0.5).data
    np.testing.assert_allclose(small, pixel, rtol=1e-14, atol=0.0)


def test_grain_mode_meta_and_mean(developed):
    state, ctx = developed
    mean = readout(state, "densitometer_visual", ctx)
    sampled = readout(state, "densitometer_visual", ctx, grain=True)
    assert "grain" not in mean.meta and mean.meta["grain_mode"] is False
    grain = sampled.meta["grain"]
    assert grain["pixels"] == mean.data.size
    assert grain["b_ap"] is None
    assert 0.0 <= grain["truncated_fraction"] <= 1.0
    assert not np.array_equal(sampled.data, mean.data)
    t_mean = np.mean(10.0**-mean.data)
    t_sampled = np.mean(10.0**-sampled.data)
    assert t_sampled == pytest.approx(t_mean, rel=0.1)
    with_ap = readout(state, "densitometer_visual", ctx, grain=True, aperture_um=48.0)
    assert with_ap.meta["grain"]["aperture_um"] == 48.0
    assert with_ap.meta["grain"]["b_ap"] is not None


def test_row_blocks_do_not_change_the_result(developed, monkeypatch):
    state, _ = developed
    product = visual_product(state.grid)
    whole, stats = dens.density_frame(state, product, grain=True)
    monkeypatch.setattr(dens, "ROW_BLOCK_PX", 7)
    blocked, blocked_stats = dens.density_frame(state, product, grain=True)
    assert np.array_equal(whole, blocked)
    assert stats.disk_area_parts == blocked_stats.disk_area_parts
    assert stats.truncated_pixels == blocked_stats.truncated_pixels


@pytest.mark.parametrize("aperture", [0.0, -1.0, math.nan, math.inf, True, "48"])
def test_bad_aperture(developed, aperture):
    state, ctx = developed
    with pytest.raises(ValueError, match="aperture_um"):
        readout(state, "densitometer_visual", ctx, aperture_um=aperture)


def test_bad_grain_and_unknown_parameter(developed):
    state, ctx = developed
    with pytest.raises(TypeError, match="grain"):
        readout(state, "densitometer_type2", ctx, grain=1)
    with pytest.raises(TypeError):
        readout(state, "densitometer_visual", ctx, colour=True)
    with pytest.raises(TypeError):
        readout(state, "silver_density", ctx, grain=True)


def test_undeveloped_state_is_rejected(developed):
    state, ctx = developed
    light = scene_light(np.full((4, 6, 3), 0.18), grid=ctx.grid)
    blank = unexposed_state(state.material, light, ctx)
    with pytest.raises((ValueError, KeyError)):
        readout(blank, "densitometer_visual", ctx)
