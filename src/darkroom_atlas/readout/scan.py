"""Standard scan of a developed ``bw-mechanistic`` negative (spec D1 2.8).

The scan density ``D_scan`` uses the product ``S_D50 ybar``; the inversion rule is fixed for
every film and developer and does not adapt to the image:

``D_net = D_scan - D_g(format)`` (base only, fog is kept),
``Y_out = 0.18 * 10^((D_net - D_ref) / G_ref)``, ``G_ref = 0.56``,
``D_ref = 0.10 + G_ref log10(10 / 0.8) = 0.714``,

and the output is ``(Y_out, Y_out, Y_out)`` in linear sRGB, which ``write_png16`` clips to
[0, 1]. An opaque pixel (``D_scan = inf``) or one whose ``Y_out`` overflows is given
``Y_out = 1``, its value after clipping.
"""

import math

import numpy as np

from darkroom_atlas.grain.statistical import grain_meta
from darkroom_atlas.pipeline import Readout, RunContext
from darkroom_atlas.readout.densitometer import density_frame
from darkroom_atlas.readout.products import scan_product
from darkroom_atlas.registry import register_readout
from darkroom_atlas.state import EmulsionState

__all__ = [
    "SCAN_D_REF",
    "SCAN_D_SPEED",
    "SCAN_G_REF",
    "SCAN_Y_REF",
    "base_density",
    "read_scan",
    "scan_inversion",
]

SCAN_G_REF: float = 0.56
"""Reference gradient: Kodak's target contrast index (spec 2.8, assumed)."""

SCAN_Y_REF: float = 0.18
"""Output luminance of the metered average (18 % grey, spec 2.8)."""

SCAN_D_SPEED: float = 0.10
"""ISO speed point above base plus fog (spec 2.8)."""

SCAN_D_REF: float = SCAN_D_SPEED + SCAN_G_REF * math.log10(10.0 / 0.8)
"""Net density of the metered average: 0.10 + 0.56 log10(12.5) = 0.714 (spec 2.8)."""


def _scan_raw(d_scan: np.ndarray, base: float) -> np.ndarray:
    """``Y_out`` before non-finite values are replaced."""
    d_net = np.asarray(d_scan, dtype=np.float64) - float(base)
    with np.errstate(over="ignore", invalid="ignore"):
        return SCAN_Y_REF * 10.0 ** ((d_net - SCAN_D_REF) / SCAN_G_REF)


def scan_inversion(d_scan: np.ndarray, base: float) -> np.ndarray:
    """``Y_out = 0.18 * 10^((D_scan - base - D_ref) / G_ref)``; non-finite results become 1."""
    y = _scan_raw(d_scan, base)
    return np.where(np.isfinite(y), y, 1.0)


def base_density(state: EmulsionState) -> float:
    """``D_g`` of the run's film format (``meta["exposure"]["film_format"]``, else default)."""
    exposure = state.meta.get("exposure")
    film_format = None if exposure is None else exposure.get("film_format")
    if film_format is None:
        film_format = state.material.default_format
    formats = state.material.formats
    if film_format not in formats:
        raise ValueError(
            f"material {state.material.id!r}: unknown film format {film_format!r}; available: "
            f"{', '.join(formats) or '(none)'}"
        )
    value = formats[film_format].grey_density.to("1")
    if np.ndim(value) != 0:
        raise ValueError(f"material {state.material.id!r}: grey_density must be one value")
    return float(value)


@register_readout("scan_standard")
def read_scan(state: EmulsionState, ctx: RunContext, *, grain: bool = True) -> Readout:
    """Standard scan as a linear sRGB image ``(H, W, 3)`` (spec 2.8)."""
    d_scan, stats = density_frame(state, scan_product(state.grid), grain=grain)
    base = base_density(state)
    y_raw = _scan_raw(d_scan, base)
    replaced = ~np.isfinite(y_raw)
    y = np.where(replaced, 1.0, y_raw)
    meta: dict = {
        "device": "scan_standard",
        "grain_mode": grain,
        "base_density": base,
        "d_ref": SCAN_D_REF,
        "g_ref": SCAN_G_REF,
        "opaque_pixels": int(np.count_nonzero(replaced)),
    }
    if stats is not None:
        meta["grain"] = grain_meta(stats)
    rgb = np.repeat(y[:, :, None], 3, axis=2)
    return Readout(kind="image_linear_srgb", data=rgb, meta=meta)
