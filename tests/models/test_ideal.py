from pathlib import Path

import numpy as np
import pytest

from darkroom_atlas.data.library import Library
from darkroom_atlas.light import Light
from darkroom_atlas.pipeline import Readout, RunContext, expose, process, readout
from darkroom_atlas.scene.camera import scene_light
from darkroom_atlas.scene.io import read_image

SAMPLE = Path(__file__).resolve().parents[1] / "assets" / "sample.png"
TIME_S = 1 / 125


@pytest.fixture(scope="module")
def lib():
    return Library()


@pytest.fixture(scope="module")
def ctx(lib):
    return RunContext(lib)


@pytest.fixture(scope="module")
def ideal(lib):
    return lib.material("ideal")


@pytest.fixture(scope="module")
def identity(lib):
    return lib.process("identity")


@pytest.fixture(scope="module")
def img():
    return read_image(SAMPLE).linear_rgb


def _develop_and_read(img, ideal, identity, ctx, **scene):
    time_s = scene.pop("time_s", TIME_S)
    light = scene_light(img, time_s=time_s, **scene)
    state = process(expose(ideal, light, time_s, ctx), identity, ctx)
    return readout(state, "ideal_colorimetric", ctx)


def test_exposure_fields(img, ideal, ctx):
    state = expose(ideal, scene_light(img, time_s=TIME_S), TIME_S, ctx)
    for layer in ("x", "y", "z"):
        exposure = state.get(layer, "exposure")
        assert exposure.shape == (64, 96)
        assert exposure.dtype == np.float64
        assert np.all(exposure >= 0.0)
    assert sorted(state.fields) == ["x/exposure", "y/exposure", "z/exposure"]


def test_y_exposure_matches_luminous(img, ideal, ctx):
    light = scene_light(img, time_s=TIME_S)
    state = expose(ideal, light, TIME_S, ctx)
    np.testing.assert_allclose(
        state.get("y", "exposure") * 683.0 * 1e12,
        light.luminous_exposure_lux_s(TIME_S),
        rtol=1e-9,
        atol=0.0,
    )


def test_end_to_end_roundtrip(img, ideal, identity, ctx):
    out = _develop_and_read(img, ideal, identity, ctx)
    assert isinstance(out, Readout)
    assert out.kind == "image_linear_srgb"
    assert out.data.shape == img.shape
    error = float(np.max(np.abs(out.data - img)))
    print(f"ideal end-to-end round trip: max abs error {error:.3e} (limit 5e-3)")
    assert error < 5e-3


def test_roundtrip_independent_of_iso_ev_time(img, ideal, identity, ctx):
    base = _develop_and_read(img, ideal, identity, ctx).data
    for scene in ({"iso": 400.0}, {"ev": -2.5}, {"time_s": 1 / 4}, {"iso": 25.0, "ev": 1.0}):
        other = _develop_and_read(img, ideal, identity, ctx, **scene).data
        np.testing.assert_allclose(other, base, rtol=0.0, atol=1e-9, err_msg=str(scene))


def test_tungsten_scene_is_warm(img, ideal, identity, ctx):
    rgb = _develop_and_read(img, ideal, identity, ctx, illuminant="A").data
    grey = rgb[32, 24]  # the middle of the grey ramp on the left half
    assert np.ptp(img[32, 24]) == 0.0
    assert grey[0] > grey[2]
    daylight = _develop_and_read(img, ideal, identity, ctx).data[32, 24]
    assert daylight[0] == pytest.approx(daylight[2], abs=5e-3)


def test_readout_requires_signal(img, ideal, ctx):
    state = expose(ideal, scene_light(img, time_s=TIME_S), TIME_S, ctx)
    with pytest.raises(ValueError, match="signal"):
        readout(state, "ideal_colorimetric", ctx)


def test_readout_requires_energy_scale(ideal, identity, ctx):
    bare = Light(np.ones((2, 3, ctx.grid.n_bands)), ctx.grid)  # not from scene_light: no meta
    state = process(expose(ideal, bare, 1.0, ctx), identity, ctx)
    with pytest.raises(ValueError, match="energy_scale_j_per_um2"):
        readout(state, "ideal_colorimetric", ctx)


def test_readout_rejects_params(img, ideal, identity, ctx):
    light = scene_light(img, time_s=TIME_S)
    state = process(expose(ideal, light, TIME_S, ctx), identity, ctx)
    with pytest.raises(TypeError, match="gamma"):
        readout(state, "ideal_colorimetric", ctx, gamma=2.2)
