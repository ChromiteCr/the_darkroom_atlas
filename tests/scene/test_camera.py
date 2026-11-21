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


@pytest.mark.parametrize("ev", [1010.0, 1023.0, 1024.0, 1100.0, -1060.0, -1100.0])
def test_rejects_ev_outside_float_range(image, ev):
    # 2**ev overflowing (ev >= 1024), the flux overflowing (ev near 1010) and the energy
    # scale underflowing to 0 or a subnormal (ev near -1060) are all one ValueError naming ev.
    with pytest.raises(ValueError, match="ev"):
        scene_light(image, ev=ev)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"iso": 1e-310},
        {"time_s": 1e-320},
        {"white_exposure_lux_s": 1e308},
        {"white_exposure_lux_s": 1e-320},
    ],
)
def test_rejects_settings_outside_float_range(image, kwargs):
    name = next(iter(kwargs))
    with pytest.raises(ValueError, match=name):
        scene_light(image, **kwargs)


@pytest.mark.parametrize("ev", [-100.0, 100.0])
def test_large_but_representable_ev(image, ev):
    light = scene_light(image, ev=ev)
    assert np.all(np.isfinite(light.photon_flux))
    assert _mean_exposure(light) == pytest.approx(0.1 * 2.0**ev, rel=1e-9)


def test_tungsten_shifts_spectrum():
    white = np.ones((1, 1, 3))
    wl = DEFAULT_GRID.wavelengths
    long, short = wl > 580.0, wl < 580.0

    def long_over_short(name):
        photons = scene_light(white, illuminant=name).photon_flux[0, 0]
        return photons[long].sum() / photons[short].sum()

    assert long_over_short("A") > 1.0
    assert long_over_short("A") > long_over_short("D65")


# --------------------------------------------------------------------------------------
# D1a: low-rank scene light (spec D1 section 4.2, section 10)
# --------------------------------------------------------------------------------------

_SETTINGS = [
    {},
    {"illuminant": "A", "iso": 400.0, "ev": -1.5, "time_s": 1 / 30},
    {"illuminant": "D55", "white_exposure_lux_s": 2.0},
]


def _d0_flux(rgb, *, illuminant="D65", iso=100.0, ev=0.0, time_s=1 / 125, white=None):
    """The D0 scene light, written out from the module docs."""
    from darkroom_atlas.scene.upsample import rgb_to_reflectance
    from darkroom_atlas.spectral import cmf_xyz
    from darkroom_atlas.spectral import illuminant as spd
    from darkroom_atlas.units import K_M, photon_energy_j

    energy = rgb_to_reflectance(rgb, DEFAULT_GRID) * spd(illuminant, DEFAULT_GRID)
    if white is not None:
        scale = white / K_M * 1e-12
    else:
        mean = float(np.mean(energy @ cmf_xyz(DEFAULT_GRID)[:, 1]))
        scale = (10.0 / iso) * 2.0**ev / (K_M * mean) * 1e-12
    return energy * (scale / time_s / photon_energy_j(DEFAULT_GRID.wavelengths))


def test_scene_light_low_rank_flag(image):
    from darkroom_atlas.light import Light

    default = scene_light(image)
    explicit = scene_light(image, low_rank=False)
    assert isinstance(explicit, Light)
    np.testing.assert_array_equal(explicit.photon_flux, default.photon_flux)
    np.testing.assert_array_equal(explicit.photon_flux, _d0_flux(image))
    assert explicit.meta == default.meta
    tungsten = scene_light(image, illuminant="A", iso=400.0, ev=-1.5, time_s=1 / 30)
    np.testing.assert_array_equal(
        tungsten.photon_flux,
        _d0_flux(image, illuminant="A", iso=400.0, ev=-1.5, time_s=1 / 30),
    )
    fixed = scene_light(image, illuminant="D55", white_exposure_lux_s=2.0, low_rank=False)
    np.testing.assert_array_equal(fixed.photon_flux, _d0_flux(image, illuminant="D55", white=2.0))


@pytest.mark.parametrize("kwargs", _SETTINGS)
def test_scene_light_low_rank_matches_full(image, kwargs):
    from darkroom_atlas.light import LowRankLight

    full = scene_light(image, **kwargs)
    low = scene_light(image, low_rank=True, **kwargs)
    assert isinstance(low, LowRankLight)
    assert low.rank == 3
    assert low.shape_hw == full.shape_hw
    assert low.meta["energy_scale_j_per_um2"] == pytest.approx(
        full.meta["energy_scale_j_per_um2"], rel=1e-12
    )
    assert {k: v for k, v in low.meta.items() if k != "energy_scale_j_per_um2"} == {
        k: v for k, v in full.meta.items() if k != "energy_scale_j_per_um2"
    }
    expanded = low.to_light().photon_flux
    mask = full.photon_flux > 0.0
    assert np.array_equal(expanded > 0.0, mask)
    rel = np.abs(expanded[mask] - full.photon_flux[mask]) / full.photon_flux[mask]
    assert rel.max() <= 1e-12


def test_scene_light_low_rank_shares_read_only_input():
    rgb = np.random.default_rng(7).random((8, 9, 3))
    rgb.flags.writeable = False
    low = scene_light(rgb, low_rank=True)
    assert low.coefficients is rgb
    writable = np.random.default_rng(8).random((8, 9, 3))
    copied = scene_light(writable, low_rank=True)
    assert not np.shares_memory(copied.coefficients, writable)
    as_float32 = scene_light(writable.astype(np.float32), low_rank=True)
    assert as_float32.coefficients.dtype == np.float64


def test_scene_light_low_rank_black_frame_and_errors(image):
    black = scene_light(np.zeros((4, 6, 3)), iso=100, low_rank=True)
    full_black = scene_light(np.zeros((4, 6, 3)), iso=100)
    assert black.meta["energy_scale_j_per_um2"] == pytest.approx(
        full_black.meta["energy_scale_j_per_um2"], rel=1e-12
    )
    for ev in (1010.0, 1024.0, -1100.0):
        with pytest.raises(ValueError, match="ev"):
            scene_light(image, ev=ev, low_rank=True)
    with pytest.raises(ValueError, match="white_exposure_lux_s"):
        scene_light(image, white_exposure_lux_s=1e308, low_rank=True)
    with pytest.raises(TypeError):
        scene_light(image, low_rank=1)
    with pytest.raises(ValueError):
        scene_light(np.full((2, 2, 3), -0.5), low_rank=True)


def test_scene_light_low_rank_no_full_spectrum_allocation():
    import tracemalloc

    rgb = np.random.default_rng(9).random((600, 400, 3))
    rgb.flags.writeable = False
    tracemalloc.start()
    try:
        scene_light(rgb, low_rank=True)
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    # The full flux would be 600 * 400 * 81 * 8 = 156 MB; the input is 5.8 MB.
    assert peak < 4 * rgb.nbytes
