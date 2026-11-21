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

With ``low_rank=True`` the same flux is returned factored (spec D1 sections 2.1.1 and 4.2):
``photon_flux = linear_rgb @ Q`` with ``Q[r, b] = k * M[b, r] * S_b / (E_b * time_s)``,
``M`` the Mallett basis, ``S`` the illuminant and ``E_b`` the photon energy, so the
``(H, W, B)`` array is never built.
"""

import math
import sys
from numbers import Real

import numpy as np

from darkroom_atlas._frozen import readonly_array
from darkroom_atlas.light import Light, LowRankLight
from darkroom_atlas.scene.upsample import mallett2019_basis, rgb_to_reflectance
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
    low_rank: bool = False,
) -> Light | LowRankLight:
    """Light on the film plane for a linear sRGB scene ``(H, W, 3)``; see the module docs.

    ``linear_rgb`` must be finite and non-negative. Settings so extreme that the energy scale
    ``k`` is not a normal floating-point number or the photon flux is not finite (for
    example ``|ev|`` near 1000) raise ``ValueError`` naming them. ``meta`` of the result records
    ``illuminant iso ev time_s metering`` (``"average"`` or ``"fixed"``) and
    ``energy_scale_j_per_um2`` (``k``).

    ``low_rank=False`` (the default) returns a :class:`Light`. ``low_rank=True`` returns a
    :class:`LowRankLight` whose coefficients are ``linear_rgb`` itself (shared when it is a
    read-only ``float64`` array, copied otherwise) and whose ``to_light()`` equals the
    ``low_rank=False`` result to rounding.
    """
    if not isinstance(low_rank, bool | np.bool_):
        raise TypeError(f"low_rank must be a bool, got {type(low_rank).__name__}")
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

    low_rank = bool(low_rank)
    if low_rank:
        # Relative energy per unit RGB, (3, B): Mallett basis times the illuminant.
        unit_energy = mallett2019_basis(grid).T * illuminant_spd(illuminant, grid)
        energy = None
    else:
        energy = rgb_to_reflectance(rgb, grid) * illuminant_spd(illuminant, grid)
    if white_exposure_lux_s is not None:
        white = _positive_finite("white_exposure_lux_s", white_exposure_lux_s)
        scale = white / K_M * _J_PER_M2_TO_J_PER_UM2
        metering = "fixed"
        settings = f"white_exposure_lux_s={white!r}"
    else:
        try:
            target_lux_s = (10.0 / iso) * 2.0**ev
        except OverflowError:  # 2.0**ev raises where the other operations give inf
            target_lux_s = math.inf
        ybar = cmf_xyz(grid)[:, 1]
        if energy is None:
            mean_luminance = float(np.mean(rgb, axis=(0, 1)) @ (unit_energy @ ybar))
        else:
            mean_luminance = float(np.mean(energy @ ybar))
        if not mean_luminance >= DARK_FRAME_LUMINANCE:
            mean_luminance = METERED_GREY
        scale = target_lux_s / (K_M * mean_luminance) * _J_PER_M2_TO_J_PER_UM2
        metering = "average"
        settings = f"iso={iso!r} and ev={ev!r}"

    # The energy scale must be a normal float (the readout divides by it) and the flux finite.
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        if energy is None:
            flux = None
            basis = unit_energy * (scale / time_s / photon_energy_j(grid.wavelengths))
            # Per-channel maxima bound every pixel's flux from above (all terms >= 0).
            finite = np.all(np.isfinite(basis)) and np.all(
                np.isfinite(np.max(rgb, axis=(0, 1)) @ basis)
            )
        else:
            flux = energy * (scale / time_s / photon_energy_j(grid.wavelengths))
            finite = np.all(np.isfinite(flux))
    if not (math.isfinite(scale) and scale >= sys.float_info.min and finite):
        raise ValueError(
            f"{settings} with time_s={time_s!r} put the exposure outside the floating-point "
            f"range (energy scale {scale:.3g} J/um^2); use less extreme values"
        )
    meta = {
        "illuminant": str(illuminant),
        "iso": iso,
        "ev": ev,
        "time_s": time_s,
        "metering": metering,
        "energy_scale_j_per_um2": float(scale),
    }
    if flux is None:
        basis.flags.writeable = False
        if rgb is not linear_rgb and rgb.base is None:
            rgb.flags.writeable = False  # a fresh conversion owned here: share it
        return LowRankLight(readonly_array(rgb), basis, grid, meta)
    flux.flags.writeable = False  # owned here, so Light can keep it without a copy
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
