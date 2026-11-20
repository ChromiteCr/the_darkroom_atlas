import math
from pathlib import Path

import numpy as np
import pytest

from darkroom_atlas.scene.camera import scene_light
from darkroom_atlas.scene.io import read_image
from darkroom_atlas.spectral import DEFAULT_GRID

SAMPLE = Path(__file__).resolve().parents[1] / "assets" / "sample.png"


@pytest.fixture(scope="module")
def image():
    return read_image(SAMPLE).linear_rgb


def _mean_exposure(light):
    return float(np.mean(light.luminous_exposure_lux_s(light.meta["time_s"])))


def test_average_metering_hits_iso(image):
    assert _mean_exposure(scene_light(image, iso=100)) == pytest.approx(0.1, rel=1e-9)
    assert _mean_exposure(scene_light(image, iso=400)) == pytest.approx(0.025, rel=1e-9)
    noise = np.random.default_rng(5).random((7, 9, 3)) * 3.0
    assert _mean_exposure(scene_light(noise, iso=50)) == pytest.approx(0.2, rel=1e-9)


def test_ev_doubles(image):
    base = _mean_exposure(scene_light(image, ev=0.0))
    assert _mean_exposure(scene_light(image, ev=1.0)) == pytest.approx(2.0 * base, rel=1e-9)
    assert _mean_exposure(scene_light(image, ev=-2.0)) == pytest.approx(0.25 * base, rel=1e-9)


def test_fixed_white_exposure():
    rgb = np.array([[[1.0, 1.0, 1.0], [0.18, 0.18, 0.18], [0.0, 0.0, 0.0]]])
    light = scene_light(rgb, white_exposure_lux_s=1.0, time_s=1 / 60)
    exposure = light.luminous_exposure_lux_s(1 / 60)
    assert exposure[0, 0] == pytest.approx(1.0, rel=1e-9)
    assert exposure[0, 1] == pytest.approx(0.18, rel=1e-9)
    assert exposure[0, 2] == 0.0
    assert light.meta["metering"] == "fixed"
    assert light.meta["energy_scale_j_per_um2"] == pytest.approx(1.0 / 683.0 * 1e-12, rel=1e-15)


def test_time_does_not_change_exposure(image):
    fast = scene_light(image, time_s=1 / 500)
    slow = scene_light(image, time_s=1 / 8)
    np.testing.assert_allclose(slow.photon_flux * (1 / 8), fast.photon_flux * (1 / 500), rtol=1e-12)
    np.testing.assert_allclose(
        slow.luminous_exposure_lux_s(1 / 8), fast.luminous_exposure_lux_s(1 / 500), rtol=1e-12
    )


def test_black_frame_metered_as_grey():
    light = scene_light(np.zeros((4, 6, 3)), iso=100)
    np.testing.assert_array_equal(light.photon_flux, 0.0)
    scale = light.meta["energy_scale_j_per_um2"]
    assert math.isfinite(scale)
    assert scale > 0.0
    # An 18% grey frame metered the same way gets the same scale.
    grey = scene_light(np.full((4, 6, 3), 0.18), iso=100)
    assert grey.meta["energy_scale_j_per_um2"] == pytest.approx(scale, rel=1e-9)


def test_meta_recorded(image):
    light = scene_light(image, illuminant="D50", iso=200.0, ev=0.5, time_s=1 / 30)
    meta = light.meta
    assert meta["illuminant"] == "D50"
    assert meta["iso"] == 200.0
    assert meta["ev"] == 0.5
    assert meta["time_s"] == 1 / 30
    assert meta["metering"] == "average"
    assert meta["energy_scale_j_per_um2"] > 0.0
    assert light.grid == DEFAULT_GRID
    assert light.shape_hw == image.shape[:2]


def test_rejects_bad_arguments(image):
    with pytest.raises(ValueError):
        scene_light(image, iso=0.0)
    with pytest.raises(ValueError):
        scene_light(image, time_s=0.0)
    with pytest.raises(ValueError):
        scene_light(image, white_exposure_lux_s=-1.0)
    with pytest.raises(ValueError):
        scene_light(image, illuminant="X")
    with pytest.raises(ValueError):
        scene_light(np.ones((4, 4)))
    with pytest.raises(ValueError):
        scene_light(np.full((2, 2, 3), np.nan))
    with pytest.raises(ValueError):
        scene_light(np.full((2, 2, 3), -0.5))


def test_tungsten_shifts_spectrum():
    white = np.ones((1, 1, 3))
    wl = DEFAULT_GRID.wavelengths
    long, short = wl > 580.0, wl < 580.0

    def long_over_short(name):
        photons = scene_light(white, illuminant=name).photon_flux[0, 0]
        return photons[long].sum() / photons[short].sum()

    assert long_over_short("A") > 1.0
    assert long_over_short("A") > long_over_short("D65")
