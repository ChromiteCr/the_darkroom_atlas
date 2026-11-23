"""Densitometers of a developed ``bw-mechanistic`` film (spec D1 2.8).

A diffuse transmission density weighted by a spectral product ``P_b``:
``D_P = -log10(sum_b P_b T_b / sum_b P_b)`` (CAL 1.4). The devices

* ``densitometer_visual``: ``P = S_A ybar`` (ISO visual density);
* ``densitometer_type2``: the ISO 5-3 Type 2 printing product (anchor approximation);
* ``silver_density``: the silver density ``D_Ag`` itself (diagnostic).

``aperture_um`` (``None``: per pixel) averages the product-weighted transmittance over a
circular aperture of that diameter before the logarithm is taken. ``grain`` (default
``False``) uses the statistical grain mode (spec 2.7.1) instead of the mean transmittance;
its record is ``Readout.meta["grain"]``.

The frame is read in row blocks; every per-pixel value is independent of the block it is
computed in (spec 3.6). The aperture average is a mean over the frame pixels whose centres
lie within the aperture; at the frame edge only the pixels inside the frame count.
"""

import math
from collections.abc import Callable
from numbers import Real

import numpy as np
from scipy.ndimage import correlate

from darkroom_atlas.grain.statistical import (
    GrainBlockStats,
    grain_meta,
    statistical_transmittance_with_stats,
)
from darkroom_atlas.imaging.silver import product_density, silver_density, transmittance_factors
from darkroom_atlas.pipeline import Readout, RunContext
from darkroom_atlas.readout.products import type2_product, visual_product
from darkroom_atlas.registry import register_readout
from darkroom_atlas.spectral import SpectralGrid
from darkroom_atlas.state import EmulsionState

__all__ = [
    "ROW_BLOCK_PX",
    "aperture_average",
    "aperture_kernel",
    "density_frame",
    "read_silver_density",
    "read_type2",
    "read_visual",
]

ROW_BLOCK_PX = 256
"""Rows read at a time (spec 4.5: the readout keeps row-block temporaries only)."""


def _check_aperture(aperture_um: object) -> float | None:
    if aperture_um is None:
        return None
    if (
        isinstance(aperture_um, bool)
        or not isinstance(aperture_um, Real)
        or not math.isfinite(aperture_um)
        or aperture_um <= 0.0
    ):
        raise ValueError(f"aperture_um must be None or positive and finite, got {aperture_um!r}")
    return float(aperture_um)


def _check_grain(grain: object) -> bool:
    if not isinstance(grain, bool):
        raise TypeError(f"grain must be True or False, got {grain!r}")
    return grain


def aperture_kernel(aperture_um: float, pixel_um: float) -> np.ndarray:
    """0/1 mask of the pixels whose centres lie within a circle of diameter ``aperture_um``."""
    radius_px = 0.5 * aperture_um / pixel_um
    n = math.floor(radius_px)
    offsets = np.arange(-n, n + 1, dtype=np.float64)
    mask = offsets[:, None] ** 2 + offsets[None, :] ** 2 <= radius_px * radius_px
    return mask.astype(np.float64)


def aperture_average(transmittance: np.ndarray, aperture_um: float, pixel_um: float) -> np.ndarray:
    """Mean of ``transmittance`` ``(H, W)`` over the aperture around each pixel (in-frame only)."""
    kernel = aperture_kernel(aperture_um, pixel_um)
    if kernel.shape == (1, 1):
        return np.array(transmittance, dtype=np.float64)
    total = correlate(transmittance, kernel, mode="constant", cval=0.0)
    count = correlate(np.ones_like(transmittance), kernel, mode="constant", cval=0.0)
    return total / count


def density_frame(
    state: EmulsionState,
    product: np.ndarray,
    *,
    aperture_um: float | None = None,
    grain: bool = False,
) -> tuple[np.ndarray, GrainBlockStats | None]:
    """Product density ``(H, W)`` of the whole frame and, with ``grain``, its grain stats."""
    aperture_um = _check_aperture(aperture_um)
    grain = _check_grain(grain)
    h, w = state.geometry.shape_hw
    cols = slice(0, w)
    out = np.empty((h, w), dtype=np.float64)
    stats: GrainBlockStats | None = None
    for r0 in range(0, h, ROW_BLOCK_PX):
        rows = slice(r0, min(h, r0 + ROW_BLOCK_PX))
        factors = transmittance_factors(state, rows, cols)
        if grain:
            t, block = statistical_transmittance_with_stats(state, rows, cols)
            stats = block if stats is None else stats.merge(block)
        else:
            t = 10.0 ** -silver_density(state, rows, cols)
        out[rows] = product_density(t, factors, product)
    if aperture_um is not None:
        t_avg = aperture_average(10.0**-out, aperture_um, state.geometry.pixel_pitch_um)
        with np.errstate(divide="ignore"):
            out = -np.log10(t_avg)
    return out, stats


def _density_readout(
    state: EmulsionState,
    ctx: RunContext,
    device: str,
    product_fn: Callable[[SpectralGrid], np.ndarray],
    aperture_um: float | None,
    grain: bool,
) -> Readout:
    density, stats = density_frame(
        state, product_fn(state.grid), aperture_um=aperture_um, grain=grain
    )
    meta: dict = {
        "device": device,
        "aperture_um": _check_aperture(aperture_um),
        "grain_mode": grain,
    }
    if stats is not None:
        meta["grain"] = grain_meta(stats, meta["aperture_um"])
    return Readout(kind="density", data=density, meta=meta)


@register_readout("densitometer_visual")
def read_visual(
    state: EmulsionState,
    ctx: RunContext,
    *,
    aperture_um: float | None = None,
    grain: bool = False,
) -> Readout:
    """ISO visual diffuse transmission density ``(H, W)`` (spec 2.8)."""
    return _density_readout(state, ctx, "densitometer_visual", visual_product, aperture_um, grain)


@register_readout("densitometer_type2")
def read_type2(
    state: EmulsionState,
    ctx: RunContext,
    *,
    aperture_um: float | None = None,
    grain: bool = False,
) -> Readout:
    """ISO Type 2 printing density ``(H, W)`` (spec 2.8; anchor approximation)."""
    return _density_readout(state, ctx, "densitometer_type2", type2_product, aperture_um, grain)


@register_readout("silver_density")
def read_silver_density(state: EmulsionState, ctx: RunContext) -> Readout:
    """The silver density ``D_Ag`` ``(H, W)`` (spec 2.5.1, 2.8; diagnostic)."""
    return Readout(kind="density", data=silver_density(state), meta={"device": "silver_density"})
