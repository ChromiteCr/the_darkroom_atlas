import numpy as np

from darkroom_atlas.scene.upsample import mallett2019_basis, rgb_to_reflectance
from darkroom_atlas.spectral import DEFAULT_GRID, XYZ_TO_LINEAR_SRGB, cmf_xyz, illuminant

PATCHES = np.array(
    [
        [1.0, 0.0, 0.0],  # red
        [0.0, 1.0, 0.0],  # green
        [0.0, 0.0, 1.0],  # blue
        [1.0, 1.0, 1.0],  # white
        [0.0, 1.0, 1.0],  # cyan
        [1.0, 0.0, 1.0],  # magenta
        [1.0, 1.0, 0.0],  # yellow
        [0.0, 0.0, 0.0],  # black
    ]
)


def test_basis_partition_of_unity():
    basis = mallett2019_basis(DEFAULT_GRID)
    assert basis.shape == (DEFAULT_GRID.n_bands, 3)
    np.testing.assert_allclose(basis.sum(axis=1), 1.0, rtol=0, atol=1e-3)
    assert basis.min() >= 0.0


def test_white_is_flat():
    white = rgb_to_reflectance(np.ones((2, 3, 3)), DEFAULT_GRID)
    assert white.shape == (2, 3, DEFAULT_GRID.n_bands)
    np.testing.assert_allclose(white, 1.0, rtol=0, atol=1e-12)
    black = rgb_to_reflectance(np.zeros((2, 3, 3)), DEFAULT_GRID)
    np.testing.assert_array_equal(black, 0.0)


def test_roundtrip_through_cmf():
    reflectance = rgb_to_reflectance(PATCHES.reshape(2, 4, 3), DEFAULT_GRID)
    energy = reflectance * illuminant("D65", DEFAULT_GRID)
    xyz = energy @ cmf_xyz(DEFAULT_GRID)
    rgb = xyz @ XYZ_TO_LINEAR_SRGB.T
    error = np.max(np.abs(rgb.reshape(8, 3) - PATCHES))
    print(f"max round-trip error: {error:.2e}")
    assert error < 5e-3
