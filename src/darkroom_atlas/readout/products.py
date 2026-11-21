"""Spectral products of the readout devices (spec D1 section 2.8).

A density is ``D = -log10(sum_b P_b T_b / sum_b P_b)``, so a product's scale does not matter;
every product here is normalised to ``sum_b P_b == 1`` and returned read-only, shape ``(B,)``.
"""

import numpy as np
from scipy.interpolate import PchipInterpolator

from darkroom_atlas.spectral import SpectralGrid, cmf_xyz, illuminant

__all__ = ["TYPE2_ANCHORS", "scan_product", "type2_log_product", "type2_product", "visual_product"]

# ISO 5-3:1995 Table 2 anchors of the Type 2 printing product, (nm, log10 P) (spec 2.8, CAL 1.2).
TYPE2_ANCHORS: tuple[tuple[float, float], ...] = (
    (360.0, 4.280),
    (400.0, 4.916),
    (430.0, 5.000),
    (500.0, 3.993),
    (540.0, 1.0),
)
# The Type 2 product is zero outside [380, 540] nm (spec 2.8).
_TYPE2_SUPPORT_NM = (380.0, 540.0)
_TYPE2_INTERP = PchipInterpolator(
    [a[0] for a in TYPE2_ANCHORS], [a[1] for a in TYPE2_ANCHORS], extrapolate=False
)


def _normalised(product: np.ndarray, what: str, grid: SpectralGrid) -> np.ndarray:
    total = float(np.sum(product))
    if not (np.isfinite(total) and total > 0.0):
        raise ValueError(f"{what} product has no weight on {grid}")
    out = product / total
    out.flags.writeable = False
    return out


def visual_product(grid: SpectralGrid) -> np.ndarray:
    """ISO visual density product: CIE illuminant A times CIE 1931 ybar (spec 2.8, CAL 1.4)."""
    return _normalised(illuminant("A", grid) * cmf_xyz(grid)[:, 1], "visual", grid)


def type2_log_product(wavelengths_nm: np.ndarray) -> np.ndarray:
    """log10 of the Type 2 product before the support cut-off (monotone cubic through anchors).

    ``nan`` outside the anchor range 360-540 nm.
    """
    return _TYPE2_INTERP(np.asarray(wavelengths_nm, dtype=np.float64))


def type2_product(grid: SpectralGrid) -> np.ndarray:
    """ISO Type 2 printing density product, approximated from the ISO 5-3 anchors (spec 2.8).

    ``10**P(lambda)`` with ``P`` the monotone cubic (PCHIP) interpolation of the anchors in
    log space, set to 0 outside [380, 540] nm.
    """
    wl = grid.wavelengths
    lo, hi = _TYPE2_SUPPORT_NM
    inside = (wl >= lo) & (wl <= hi)
    product = np.zeros(wl.shape, dtype=np.float64)
    product[inside] = 10.0 ** type2_log_product(wl[inside])
    return _normalised(product, "type2", grid)


def scan_product(grid: SpectralGrid) -> np.ndarray:
    """Standard scan product: CIE D50 times CIE 1931 ybar (spec 2.8, fixed convention)."""
    return _normalised(illuminant("D50", grid) * cmf_xyz(grid)[:, 1], "scan", grid)
