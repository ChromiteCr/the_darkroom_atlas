"""A virtual camera: turn a linear RGB scene into the light on the film plane.

The scene is upsampled to reflectance spectra, lit by a standard illuminant and scaled so
that the film receives a chosen photometric exposure:

* ``e = reflectance * illuminant`` is the relative band-integrated energy per pixel; the
  illuminant is normalised so that an ideal white has ``sum(e * ybar) == 1``.
* ``k`` (J per um^2 for relative energy 1) sets the absolute level. A pixel's photopic
  exposure is then ``H = 683 * 1e12 * k * sum(e * ybar)`` lux s.
* Average metering (the default) chooses ``k`` so that the mean of ``H`` over the frame is
  ``(10 / iso) * 2**ev`` lux s, the mean film-plane exposure ``H = 10 / S`` that the
  design takes from ISO 2721 for a film of speed ``S``. A frame darker than
  ``1e-9`` relative luminance (for example all black) is metered as 18% grey instead.
* A fixed ``white_exposure_lux_s`` instead gives ``k = white_exposure_lux_s / 683 * 1e-12``,
  so an ideal white receives exactly that exposure; ``iso`` and ``ev`` then only appear in
  the metadata.
* ``photon_flux = k * e / photon_energy / time_s``: the exposure time changes the flux,
  not the exposure.
"""

import math
from numbers import Real

import numpy as np

from darkroom_atlas.light import Light
from darkroom_atlas.scene.upsample import rgb_to_reflectance
from darkroom_atlas.spectral import DEFAULT_GRID, SpectralGrid, cmf_xyz
from darkroom_atlas.spectral import illuminant as illuminant_spd
from darkroom_atlas.units import K_M, photon_energy_j

__all__ = ["DARK_FRAME_LUMINANCE", "METERED_GREY", "scene_light"]

DARK_FRAME_LUMINANCE = 1e-9
"""Mean relative luminance below which a frame is metered as :data:`METERED_GREY`."""

METERED_GREY = 0.18
"""Relative luminance assumed for frames too dark to meter."""

_J_PER_M2_TO_J_PER_UM2 = 1e-12


def scene_light(
    linear_rgb,
    *,
    grid: SpectralGrid = DEFAULT_GRID,
    illuminant: str = "D65",
    iso: float = 100.0,
    ev: float = 0.0,
    time_s: float = 1 / 125,
    white_exposure_lux_s: float | None = None,
) -> Light:
    """Light on the film plane for a linear sRGB scene ``(H, W, 3)``; see the module docs.

    ``linear_rgb`` must be finite and non-negative. ``meta`` of the result records
    ``illuminant iso ev time_s metering`` (``"average"`` or ``"fixed"``) and
    ``energy_scale_j_per_um2`` (``k``).
    """
    rgb = np.asarray(linear_rgb, dtype=np.float64)
    if rgb.ndim != 3 or rgb.shape[2] != 3 or rgb.shape[0] < 1 or rgb.shape[1] < 1:
        raise ValueError(f"linear_rgb must have shape (H, W, 3), got {rgb.shape}")
    if not np.all(np.isfinite(rgb)):
        raise ValueError("linear_rgb must be finite (found NaN or inf)")
    if np.any(rgb < 0.0):
        raise ValueError("linear_rgb must be non-negative")
    iso = _positive_finite("iso", iso)
    time_s = _positive_finite("time_s", time_s)
    if isinstance(ev, bool) or not isinstance(ev, Real) or not math.isfinite(ev):
        raise ValueError(f"ev must be a finite number, got {ev!r}")
    ev = float(ev)

    energy = rgb_to_reflectance(rgb, grid) * illuminant_spd(illuminant, grid)
    if white_exposure_lux_s is not None:
        white = _positive_finite("white_exposure_lux_s", white_exposure_lux_s)
        scale = white / K_M * _J_PER_M2_TO_J_PER_UM2
        metering = "fixed"
    else:
        target_lux_s = (10.0 / iso) * 2.0**ev
        mean_luminance = float(np.mean(energy @ cmf_xyz(grid)[:, 1]))
        if not mean_luminance >= DARK_FRAME_LUMINANCE:
            mean_luminance = METERED_GREY
        scale = target_lux_s / (K_M * mean_luminance) * _J_PER_M2_TO_J_PER_UM2
        metering = "average"

    flux = energy * (scale / time_s / photon_energy_j(grid.wavelengths))
    flux.flags.writeable = False  # owned here, so Light can keep it without a copy
    meta = {
        "illuminant": str(illuminant),
        "iso": iso,
        "ev": ev,
        "time_s": time_s,
        "metering": metering,
        "energy_scale_j_per_um2": float(scale),
    }
    return Light(flux, grid, meta)


def _positive_finite(name: str, value: object) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, Real)
        or not math.isfinite(value)
        or value <= 0.0
    ):
        raise ValueError(f"{name} must be a positive, finite number, got {value!r}")
    return float(value)
