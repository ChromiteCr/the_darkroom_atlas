import numpy as np
import pytest

from darkroom_atlas.spectral import (
    DEFAULT_GRID,
    LINEAR_SRGB_TO_XYZ,
    XYZ_TO_LINEAR_SRGB,
    SpectralGrid,
    blackbody,
    cmf_xyz,
    illuminant,
)


def _xy(spd, grid=DEFAULT_GRID):
    xyz = spd @ cmf_xyz(grid)
    return xyz[:2] / xyz.sum()


def test_default_grid():
    grid = DEFAULT_GRID
    assert grid == SpectralGrid()
    assert grid.n_bands == 81
    assert grid.wavelengths.shape == (81,)
    assert grid.wavelengths[0] == 380
    assert grid.wavelengths[-1] == 780


def test_grid_rejects_non_integer():
    with pytest.raises(ValueError):
        SpectralGrid(380, 780, 7)
    with pytest.raises(ValueError):
        SpectralGrid(380, 780, 0)
    with pytest.raises(ValueError):
        SpectralGrid(380, 780, -5)
    with pytest.raises(ValueError):
        SpectralGrid(780, 380, 5)
    with pytest.raises(ValueError):
        SpectralGrid(500, 500, 5)


def test_resample_linear_exact():
    x = np.arange(380.0, 781.0, 10.0)
    y = 2.0 * x + 1.0
    out = DEFAULT_GRID.resample(x, y)
    expected = 2.0 * DEFAULT_GRID.wavelengths + 1.0
    assert out.shape == (81,)
    assert np.max(np.abs(out - expected)) < 1e-9


def test_resample_outside_modes():
    x = np.arange(400.0, 701.0, 10.0)
    y = 0.01 * x + 3.0
    with pytest.raises(ValueError):
        DEFAULT_GRID.resample(x, y)
    with pytest.raises(ValueError):
        DEFAULT_GRID.resample(x, y, outside="error")
    wl = DEFAULT_GRID.wavelengths
    i380 = int(np.flatnonzero(wl == 380.0)[0])
    i400 = int(np.flatnonzero(wl == 400.0)[0])
    zero = DEFAULT_GRID.resample(x, y, outside="zero")
    assert zero[i380] == 0.0
    assert zero[i400] == pytest.approx(7.0)
    hold = DEFAULT_GRID.resample(x, y, outside="hold")
    assert hold[i380] == hold[i400]
    assert hold[i400] == pytest.approx(7.0)


def test_resample_axis_and_shape_errors():
    x = np.arange(380.0, 781.0, 10.0)
    y = np.stack([x, 2.0 * x], axis=-1)  # (N, 2): wavelength on axis 0
    with pytest.raises(ValueError):
        DEFAULT_GRID.resample(x, y)
    out = DEFAULT_GRID.resample(x, y, axis=0)
    assert out.shape == (81, 2)
    np.testing.assert_allclose(out[:, 1], 2.0 * DEFAULT_GRID.wavelengths, rtol=0, atol=1e-9)
    with pytest.raises(ValueError):
        DEFAULT_GRID.resample(x[::-1], y[::-1], axis=0)


def test_integrate_sums_bands():
    values = np.ones((2, 3, DEFAULT_GRID.n_bands))
    assert np.array_equal(DEFAULT_GRID.integrate(values), np.full((2, 3), 81.0))
    with pytest.raises(ValueError):
        DEFAULT_GRID.integrate(np.ones((2, 80)))


def test_cmf_shape_and_peak():
    cmf = cmf_xyz(DEFAULT_GRID)
    assert cmf.shape == (81, 3)
    i555 = int(np.flatnonzero(DEFAULT_GRID.wavelengths == 555.0)[0])
    assert cmf[i555, 1] == pytest.approx(1.0, abs=2e-3)


@pytest.mark.parametrize("name", ["D65", "D50", "A", "E"])
def test_illuminant_normalised(name):
    spd = illuminant(name, DEFAULT_GRID)
    ybar = cmf_xyz(DEFAULT_GRID)[:, 1]
    assert spd.shape == (81,)
    assert np.sum(spd * ybar) == pytest.approx(1.0, rel=1e-12)


@pytest.mark.parametrize(
    "name, expected",
    [("D65", (0.3127, 0.3290)), ("A", (0.4476, 0.4074)), ("E", (1 / 3, 1 / 3))],
)
def test_illuminant_chromaticity(name, expected):
    xy = _xy(illuminant(name, DEFAULT_GRID))
    assert xy[0] == pytest.approx(expected[0], abs=1e-3)
    assert xy[1] == pytest.approx(expected[1], abs=1e-3)


def test_blackbody_2856_matches_A():
    bb = blackbody(2856.0, DEFAULT_GRID)
    ybar = cmf_xyz(DEFAULT_GRID)[:, 1]
    assert np.sum(bb * ybar) == pytest.approx(1.0, rel=1e-12)
    diff = _xy(bb) - _xy(illuminant("A", DEFAULT_GRID))
    assert np.max(np.abs(diff)) < 2e-3


def test_blackbody_rejects_bad_temperature():
    for t in (0.0, -10.0, float("nan"), float("inf")):
        with pytest.raises(ValueError):
            blackbody(t, DEFAULT_GRID)


def test_unknown_illuminant():
    with pytest.raises(ValueError) as excinfo:
        illuminant("X", DEFAULT_GRID)
    message = str(excinfo.value)
    for name in ("D65", "D50", "A", "E"):
        assert name in message


def test_cached_data_not_shared():
    cmf = cmf_xyz(DEFAULT_GRID)
    cmf[:] = 0.0
    assert cmf_xyz(DEFAULT_GRID).max() > 0.5
    spd = illuminant("D65", DEFAULT_GRID)
    spd[:] = 0.0
    assert illuminant("D65", DEFAULT_GRID).max() > 0.0


def test_srgb_matrices_inverse():
    np.testing.assert_allclose(
        XYZ_TO_LINEAR_SRGB @ LINEAR_SRGB_TO_XYZ, np.eye(3), rtol=0, atol=1e-6
    )
    assert not XYZ_TO_LINEAR_SRGB.flags.writeable
    assert not LINEAR_SRGB_TO_XYZ.flags.writeable


def test_d65_white_maps_to_srgb_white():
    # Grid-sampled D65 lands within 3e-4 of the matrices' (0.3127, 0.3290) white point.
    xyz = illuminant("D65", DEFAULT_GRID) @ cmf_xyz(DEFAULT_GRID)
    np.testing.assert_allclose(XYZ_TO_LINEAR_SRGB @ xyz, np.ones(3), rtol=0, atol=1e-3)


def test_srgb_matrix_matches_published():
    # IEC 61966-2-1 publishes the RGB -> XYZ matrix rounded to four decimals.
    published = np.array(
        [
            [0.4124, 0.3576, 0.1805],
            [0.2126, 0.7152, 0.0722],
            [0.0193, 0.1192, 0.9505],
        ]
    )
    np.testing.assert_allclose(LINEAR_SRGB_TO_XYZ, published, rtol=0, atol=5e-5)
    np.testing.assert_allclose(LINEAR_SRGB_TO_XYZ[1].sum(), 1.0, rtol=0, atol=1e-12)
