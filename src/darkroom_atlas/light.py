"""Light arriving at the film plane, as a spectral photon flux per pixel.

A :class:`Light` holds ``photon_flux`` of shape ``(H, W, B)``: photons per um^2 per second in
each band of its :class:`~darkroom_atlas.spectral.SpectralGrid` (band-integrated, so summing
over the last axis integrates over wavelength). Multiplying by an exposure time gives the
exposure in photons per um^2 per band, the unit the exposure models consume.

A :class:`LowRankLight` holds the same quantity in factored form (spec D1 section 4.2):
``photon_flux[y, x, b] = sum_r coefficients[y, x, r] * basis[r, b]``. A linear RGB scene
needs only ``R = 3`` coefficients per pixel, so the ``(H, W, B)`` flux is never stored; it
is produced only by :meth:`LowRankLight.window` (one block) and :meth:`LowRankLight.to_light`.
Both classes satisfy the :class:`LightField` protocol the exposure models consume.
"""

import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from numbers import Real
from typing import Any, Protocol, runtime_checkable

import numpy as np

from darkroom_atlas._frozen import freeze_mapping, readonly_array
from darkroom_atlas.backend import xp_of
from darkroom_atlas.spectral import SpectralGrid, cmf_xyz
from darkroom_atlas.units import K_M, photon_energy_j

__all__ = ["Light", "LightField", "LowRankLight", "apply_transmittance"]

_J_PER_UM2_TO_J_PER_M2 = 1e12


@runtime_checkable
class LightField(Protocol):
    """Light on the film plane as the exposure models see it (full or low-rank)."""

    @property
    def grid(self) -> SpectralGrid: ...

    @property
    def meta(self) -> Mapping[str, Any]: ...

    @property
    def shape_hw(self) -> tuple[int, int]: ...

    def window(self, rows: slice, cols: slice) -> "Light":
        """Full-spectrum light of the block ``[rows, cols]``."""
        ...

    def project(
        self, weights: np.ndarray, rows: slice | None = None, cols: slice | None = None
    ) -> np.ndarray:
        """``photon_flux[rows, cols] @ weights``: weights ``(B, M)`` give ``(h, w, M)``."""
        ...

    def luminous_exposure_lux_s(self, time_s: float) -> np.ndarray: ...


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

    def window(self, rows: slice, cols: slice) -> "Light":
        """The block ``[rows, cols]`` as a :class:`Light` sharing this flux (no copy).

        The block must contain at least one pixel, otherwise ``ValueError``.
        """
        _check_slices(rows, cols)
        return Light(self.photon_flux[rows, cols], self.grid, self.meta)

    def project(
        self, weights: np.ndarray, rows: slice | None = None, cols: slice | None = None
    ) -> np.ndarray:
        """``photon_flux[rows, cols] @ weights`` for ``weights`` ``(B, M)``: shape ``(h, w, M)``.

        ``rows`` and ``cols`` default to the whole frame. ``weights`` must be finite.
        """
        w = _check_weights(weights, self.grid.n_bands)
        flux = self.photon_flux[_full(rows), _full(cols)]
        _check_block(flux.shape[:2])
        return flux @ w

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
        weights = _luminous_weights(self.grid)
        return K_M * _J_PER_UM2_TO_J_PER_M2 * time_s * (self.photon_flux @ weights)


@dataclass(frozen=True, eq=False)
class LowRankLight:
    """Spectral photon flux on the film plane in factored form (spec D1 section 4.2).

    ``photon_flux = coefficients @ basis``: ``coefficients`` ``(H, W, R)`` (for scene light,
    the linear RGB values) and ``basis`` ``(R, B)`` in photons per um^2 per second per band
    per unit coefficient, with ``B == grid.n_bands``. Both must be finite and non-negative,
    otherwise ``ValueError``. Both are stored read-only as ``float64``; an input array that is
    already read-only (down its whole ``base`` chain) is shared, not copied.
    """

    coefficients: np.ndarray
    basis: np.ndarray
    grid: SpectralGrid
    meta: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.grid, SpectralGrid):
            raise TypeError(f"grid must be a SpectralGrid, got {type(self.grid).__name__}")
        try:
            coefficients = readonly_array(self.coefficients, dtype=np.float64)
            basis = readonly_array(self.basis, dtype=np.float64)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"coefficients and basis must be real numeric arrays: {exc}") from exc
        if coefficients.ndim != 3 or 0 in coefficients.shape:
            raise ValueError(
                f"coefficients must have shape (H, W, R) with H, W, R >= 1, "
                f"got shape {coefficients.shape}"
            )
        rank = coefficients.shape[2]
        if basis.shape != (rank, self.grid.n_bands):
            raise ValueError(
                f"basis must have shape (R, B) = ({rank}, {self.grid.n_bands}), "
                f"got shape {basis.shape}"
            )
        for name, array in (("coefficients", coefficients), ("basis", basis)):
            if not np.all(np.isfinite(array)):
                raise ValueError(f"{name} must be finite (found NaN or inf)")
            if np.any(array < 0.0):
                raise ValueError(f"{name} must be non-negative")
        object.__setattr__(self, "coefficients", coefficients)
        object.__setattr__(self, "basis", basis)
        object.__setattr__(self, "meta", freeze_mapping(self.meta))

    @property
    def shape_hw(self) -> tuple[int, int]:
        """Frame shape ``(H, W)`` in pixels."""
        h, w, _ = self.coefficients.shape
        return (h, w)

    @property
    def rank(self) -> int:
        """Number of basis spectra ``R``."""
        return self.coefficients.shape[2]

    @property
    def nbytes(self) -> int:
        """Bytes held by ``coefficients`` and ``basis``."""
        return self.coefficients.nbytes + self.basis.nbytes

    def window(self, rows: slice, cols: slice) -> Light:
        """Full-spectrum :class:`Light` of the block ``[rows, cols]``, shape ``(h, w, B)``."""
        _check_slices(rows, cols)
        block = self.coefficients[rows, cols]
        _check_block(block.shape[:2])
        flux = block @ self.basis
        flux.flags.writeable = False
        return Light(flux, self.grid, self.meta)

    def project(
        self, weights: np.ndarray, rows: slice | None = None, cols: slice | None = None
    ) -> np.ndarray:
        """``photon_flux[rows, cols] @ weights`` for ``weights`` ``(B, M)``: shape ``(h, w, M)``.

        Computed as ``coefficients[rows, cols] @ (basis @ weights)``, without the full flux.
        """
        w = _check_weights(weights, self.grid.n_bands)
        block = self.coefficients[_full(rows), _full(cols)]
        _check_block(block.shape[:2])
        return block @ (self.basis @ w)

    def luminous_exposure_lux_s(self, time_s: float) -> np.ndarray:
        """Photopic exposure in lux s over ``time_s`` seconds, shape ``(H, W)``.

        Same definition as :meth:`Light.luminous_exposure_lux_s`.
        """
        time_s = _positive_finite("time_s", time_s)
        weights = _luminous_weights(self.grid)
        return K_M * _J_PER_UM2_TO_J_PER_M2 * time_s * (self.coefficients @ (self.basis @ weights))

    def filtered(self, transmittance) -> "LowRankLight":
        """This light through a uniform filter ``transmittance`` ``(B,)`` in ``[0, 1]``.

        Only the basis changes (the coefficients are shared). A per-pixel transmittance would
        break the factorisation; use :meth:`to_light` first. Bad input raises ``ValueError``.
        """
        t = np.asarray(transmittance, dtype=np.float64)
        if t.shape != (self.grid.n_bands,):
            raise ValueError(
                f"LowRankLight.filtered takes a transmittance of shape (B,) = "
                f"({self.grid.n_bands},), got {t.shape}"
            )
        if not np.all((t >= 0.0) & (t <= 1.0)):
            raise ValueError("transmittance must lie in [0, 1] (NaN is not allowed)")
        basis = self.basis * t
        basis.flags.writeable = False
        return LowRankLight(self.coefficients, basis, self.grid, self.meta)

    def to_light(self) -> Light:
        """The full ``(H, W, B)`` :class:`Light` (materialises the whole flux)."""
        h, w = self.shape_hw
        return self.window(slice(0, h), slice(0, w))


def _luminous_weights(grid: SpectralGrid) -> np.ndarray:
    return photon_energy_j(grid.wavelengths) * cmf_xyz(grid)[:, 1]


def _full(index: slice | None) -> slice:
    if index is None:
        return slice(None)
    if not isinstance(index, slice):
        raise TypeError(f"rows and cols must be slices or None, got {type(index).__name__}")
    return index


def _check_slices(rows: object, cols: object) -> None:
    for name, index in (("rows", rows), ("cols", cols)):
        if not isinstance(index, slice):
            raise TypeError(f"{name} must be a slice, got {type(index).__name__}")


def _check_block(shape_hw: tuple[int, ...]) -> None:
    if shape_hw[0] < 1 or shape_hw[1] < 1:
        raise ValueError(f"the selected block is empty (shape {tuple(shape_hw)})")


def _check_weights(weights: object, n_bands: int) -> np.ndarray:
    try:
        w = np.asarray(weights, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"weights must be a real numeric array: {exc}") from exc
    if w.ndim != 2 or w.shape[0] != n_bands or w.shape[1] < 1:
        raise ValueError(f"weights must have shape (B, M) with B = {n_bands}, got {w.shape}")
    if not np.all(np.isfinite(w)):
        raise ValueError("weights must be finite (found NaN or inf)")
    return w


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
