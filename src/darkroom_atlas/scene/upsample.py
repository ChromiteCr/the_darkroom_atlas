"""Spectral upsampling: linear sRGB to reflectance spectra (Mallett & Yuksel 2019).

Mallett & Yuksel, "Spectral Primary Decomposition for Rendering with sRGB Reflectance",
EGSR 2019, give three non-negative basis spectra whose weighted sum, with the linear sRGB
values as weights, reproduces those values under D65 for the CIE 1931 observer. The
mapping is linear, so ``(1, 1, 1)`` gives a flat reflectance of 1.

This is the least physical step of the chain: an RGB triple says nothing about which of
the many spectra with that colour (metamers) the scene really had, and the basis picks one
smooth spectrum per colour. Values above 1 (scene-linear highlights) give "reflectances"
above 1, which is what a light source or a specular highlight needs.

The tabulated basis comes from colour-science and is resampled to the grid (holding the
end values outside 380-780 nm). Its bands sum to 1 only to about 1e-5; each band is divided
by its sum so that white maps to exactly 1 and metering anchors (white, 18% grey) are exact.
"""

from functools import lru_cache

import numpy as np

from darkroom_atlas.spectral import SpectralGrid, _colour

__all__ = ["mallett2019_basis", "rgb_to_reflectance"]


@lru_cache(maxsize=32)
def _basis_cached(grid: SpectralGrid) -> np.ndarray:
    msds = _colour().recovery.MSDS_BASIS_FUNCTIONS_sRGB_MALLETT2019
    basis = grid.resample(msds.wavelengths, msds.values, outside="hold", axis=0)
    basis = basis / basis.sum(axis=1, keepdims=True)
    basis.flags.writeable = False
    return basis


def mallett2019_basis(grid: SpectralGrid) -> np.ndarray:
    """Mallett & Yuksel 2019 sRGB basis spectra on ``grid``, shape ``(B, 3)`` (red, green, blue).

    Each band is renormalised so that the three basis values sum to exactly 1.
    """
    return _basis_cached(grid).copy()


def rgb_to_reflectance(linear_rgb, grid: SpectralGrid) -> np.ndarray:
    """Reflectance spectra for linear sRGB values ``(..., 3)``, shape ``(..., B)``."""
    rgb = np.asarray(linear_rgb, dtype=np.float64)
    if rgb.ndim < 1 or rgb.shape[-1] != 3:
        raise ValueError(f"linear_rgb must have shape (..., 3), got {rgb.shape}")
    return rgb @ _basis_cached(grid).T
