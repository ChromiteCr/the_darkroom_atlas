"""Light arriving at the film plane, as a spectral photon flux per pixel.

A :class:`Light` holds ``photon_flux`` of shape ``(H, W, B)``: photons per um^2 per second in
each band of its :class:`~darkroom_atlas.spectral.SpectralGrid` (band-integrated, so summing
over the last axis integrates over wavelength). Multiplying by an exposure time gives the
exposure in photons per um^2 per band, the unit the exposure models consume.
"""

import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from numbers import Real
from typing import Any

import numpy as np

from darkroom_atlas._frozen import freeze_mapping, readonly_array
from darkroom_atlas.backend import xp_of
from darkroom_atlas.spectral import SpectralGrid, cmf_xyz
from darkroom_atlas.units import K_M, photon_energy_j

__all__ = ["Light", "apply_transmittance"]

_J_PER_UM2_TO_J_PER_M2 = 1e12


@dataclass(frozen=True, eq=False)
class Light:
    """Spectral photon flux on the film plane.

    ``photon_flux`` has shape ``(H, W, B)`` with ``B == grid.n_bands``, in photons per um^2
    per second per band; it must be finite and non-negative, otherwise ``ValueError``. It is
    stored as a read-only ``float64`` array (copied unless it already is one). ``meta``
    describes how the light was made (for scene light: illuminant, ISO, EV, time, metering
    and energy scale) and is stored as a read-only mapping.
    """

    photon_flux: np.ndarray
    grid: SpectralGrid
    meta: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.grid, SpectralGrid):
            raise TypeError(f"grid must be a SpectralGrid, got {type(self.grid).__name__}")
        try:
            flux = readonly_array(self.photon_flux, dtype=np.float64)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"photon_flux must be a real numeric array: {exc}") from exc
        if flux.ndim != 3:
            raise ValueError(f"photon_flux must have shape (H, W, B), got shape {flux.shape}")
        if flux.shape[0] < 1 or flux.shape[1] < 1:
            raise ValueError(f"photon_flux must have at least one pixel, got shape {flux.shape}")
        if flux.shape[2] != self.grid.n_bands:
            raise ValueError(
                f"photon_flux has {flux.shape[2]} bands but the grid has {self.grid.n_bands}"
            )
        if not np.all(np.isfinite(flux)):
            raise ValueError("photon_flux must be finite (found NaN or inf)")
        if np.any(flux < 0.0):
            raise ValueError("photon_flux must be non-negative")
        object.__setattr__(self, "photon_flux", flux)
        object.__setattr__(self, "meta", freeze_mapping(self.meta))

    @property
    def shape_hw(self) -> tuple[int, int]:
        """Frame shape ``(H, W)`` in pixels."""
        h, w, _ = self.photon_flux.shape
        return (h, w)

    def filtered(self, transmittance) -> "Light":
        """Return this light after passing through a filter of the given ``transmittance``.

        ``transmittance`` has shape ``(B,)`` (the same filter for every pixel) or
        ``(H, W, B)`` (per pixel, e.g. a processed negative being printed), with values in
        ``[0, 1]``; anything else raises ``ValueError``. The grid and ``meta`` are kept.
        """
        t = np.asarray(transmittance, dtype=np.float64)
        flux = apply_transmittance(self.photon_flux, t)
        flux.flags.writeable = False
        return Light(flux, self.grid, self.meta)

    def luminous_exposure_lux_s(self, time_s: float) -> np.ndarray:
        """Photopic exposure in lux s accumulated over ``time_s`` seconds, shape ``(H, W)``.

        ``H = 683 * sum_b(photons_b * E_b * 1e12 * ybar_b)``, where ``photons_b`` is the
        photon count per um^2 in band ``b`` (flux times ``time_s``), ``E_b`` the photon energy
        in J at the band centre, ``1e12`` converts per um^2 to per m^2 and ``ybar`` is the CIE
        1931 photopic luminosity function.
        """
        time_s = _positive_finite("time_s", time_s)
        wavelengths = self.grid.wavelengths
        weights = photon_energy_j(wavelengths) * cmf_xyz(self.grid)[:, 1]
        return K_M * _J_PER_UM2_TO_J_PER_M2 * time_s * (self.photon_flux @ weights)


def apply_transmittance(photon_flux, transmittance):
    """Attenuate ``photon_flux`` ``(H, W, B)`` by ``transmittance`` ``(B,)`` or ``(H, W, B)``.

    Written against the Array API standard: both arrays must come from the same namespace,
    and the result does too. ``transmittance`` must lie in ``[0, 1]`` (NaN is rejected);
    otherwise, or if the shapes do not fit, ``ValueError``.
    """
    xp = xp_of(photon_flux, transmittance)
    if photon_flux.ndim != 3:
        raise ValueError(f"photon_flux must have shape (H, W, B), got {photon_flux.shape}")
    if tuple(transmittance.shape) not in (tuple(photon_flux.shape[-1:]), tuple(photon_flux.shape)):
        raise ValueError(
            f"transmittance must have shape (B,) or (H, W, B) = {tuple(photon_flux.shape)}, "
            f"got {tuple(transmittance.shape)}"
        )
    in_range = xp.all((transmittance >= 0.0) & (transmittance <= 1.0))
    if not bool(in_range):
        raise ValueError("transmittance must lie in [0, 1] (NaN is not allowed)")
    return photon_flux * transmittance


def _positive_finite(name: str, value: object) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, Real)
        or not math.isfinite(value)
        or value <= 0.0
    ):
        raise ValueError(f"{name} must be a positive, finite number, got {value!r}")
    return float(value)
