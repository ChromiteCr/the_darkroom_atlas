"""Spectral grid, CIE colour-matching functions, standard illuminants and sRGB matrices.

Conventions:

* Wavelengths are in nm on a uniform :class:`SpectralGrid` (380-780 nm in 5 nm steps by
  default). The band axis is always the last axis of spectral arrays, except for the
  ``(B, 3)`` colour-matching and basis matrices, which are meant for ``spectrum @ matrix``.
* Spectral quantities (energies, photon counts, illuminant powers) are stored *integrated
  over each band*, so integrating over wavelength is a plain sum over the band axis.
  Colour-matching functions and sensitivities are dimensionless weights sampled at the band
  centres, so ``XYZ = spectrum @ cmf_xyz(grid)``.
* Illuminants and blackbodies are returned as relative band-integrated energy normalised to
  ``sum(spd * ybar) == 1``, i.e. a perfect white reflector under them has ``Y == 1``.

CIE data come from colour-science and are resampled to the grid; no tables are copied here.
"""

import math
import warnings
from dataclasses import dataclass
from functools import lru_cache
from numbers import Real
from typing import Literal

import numpy as np

from darkroom_atlas.backend import xp_of
from darkroom_atlas.units import C_LIGHT, H_PLANCK, K_BOLTZMANN

__all__ = [
    "DEFAULT_GRID",
    "LINEAR_SRGB_TO_XYZ",
    "XYZ_TO_LINEAR_SRGB",
    "SpectralGrid",
    "blackbody",
    "cmf_xyz",
    "illuminant",
]

# Grid points within this distance (nm) of a data endpoint count as inside the data range.
_RANGE_TOLERANCE_NM = 1e-6


@dataclass(frozen=True)
class SpectralGrid:
    """Uniform wavelength grid from ``start_nm`` to ``stop_nm`` inclusive, in ``step_nm`` steps.

    ``stop_nm - start_nm`` must be a whole number of steps.
    """

    start_nm: float = 380.0
    stop_nm: float = 780.0
    step_nm: float = 5.0

    def __post_init__(self) -> None:
        for name in ("start_nm", "stop_nm", "step_nm"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
                raise ValueError(f"SpectralGrid.{name} must be a finite number, got {value!r}")
            object.__setattr__(self, name, float(value))
        if self.step_nm <= 0.0:
            raise ValueError(f"SpectralGrid.step_nm must be positive, got {self.step_nm}")
        if self.start_nm <= 0.0:
            raise ValueError(f"SpectralGrid.start_nm must be positive, got {self.start_nm}")
        if self.stop_nm <= self.start_nm:
            raise ValueError(
                f"SpectralGrid.stop_nm ({self.stop_nm}) must be greater than "
                f"start_nm ({self.start_nm})"
            )
        steps = (self.stop_nm - self.start_nm) / self.step_nm
        if abs(steps - round(steps)) > 1e-9 * max(1.0, steps):
            raise ValueError(
                f"SpectralGrid span {self.start_nm}-{self.stop_nm} nm is not a whole number "
                f"of {self.step_nm} nm steps"
            )

    @property
    def n_bands(self) -> int:
        """Number of bands (grid points)."""
        return round((self.stop_nm - self.start_nm) / self.step_nm) + 1

    @property
    def wavelengths(self) -> np.ndarray:
        """Band-centre wavelengths in nm, shape ``(n_bands,)``."""
        return self.start_nm + self.step_nm * np.arange(self.n_bands, dtype=np.float64)

    def resample(
        self,
        wavelengths_nm,
        values,
        *,
        outside: Literal["zero", "hold", "error"] = "error",
        axis: int = -1,
    ) -> np.ndarray:
        """Linearly interpolate ``values`` sampled at ``wavelengths_nm`` onto this grid.

        ``wavelengths_nm`` must be one-dimensional and strictly increasing; ``values`` holds
        the samples along ``axis`` (the last axis by default) and may have any other axes.
        Grid points outside the data range are handled by ``outside``: ``"error"`` raises
        ``ValueError``, ``"zero"`` sets them to 0 and ``"hold"`` repeats the nearest endpoint
        value. The result has the same shape as ``values`` with ``axis`` resized to
        ``n_bands``.

        This interpolates point samples. To obtain band-integrated values from a spectral
        density (per nm), multiply the result by ``step_nm``.
        """
        if outside not in ("zero", "hold", "error"):
            raise ValueError(f"outside must be 'zero', 'hold' or 'error', got {outside!r}")
        x = np.asarray(wavelengths_nm, dtype=np.float64)
        y = np.asarray(values, dtype=np.float64)
        if x.ndim != 1 or x.size < 2:
            raise ValueError("wavelengths_nm must be one-dimensional with at least two samples")
        if not np.all(np.isfinite(x)):
            raise ValueError("wavelengths_nm must be finite")
        if np.any(np.diff(x) <= 0.0):
            raise ValueError("wavelengths_nm must be strictly increasing")
        if y.ndim == 0:
            raise ValueError("values must have at least one dimension")
        try:
            y = np.moveaxis(y, axis, -1)
        except np.exceptions.AxisError as exc:
            raise ValueError(f"axis {axis} is out of range for values of shape {y.shape}") from exc
        if y.shape[-1] != x.size:
            raise ValueError(
                f"values have {y.shape[-1]} samples along axis {axis} but {x.size} "
                "wavelengths were given"
            )
        if not np.all(np.isfinite(y)):
            raise ValueError("values must be finite")

        wl = self.wavelengths
        out_of_range = (wl < x[0] - _RANGE_TOLERANCE_NM) | (wl > x[-1] + _RANGE_TOLERANCE_NM)
        if outside == "error" and out_of_range.any():
            raise ValueError(
                f"data cover {x[0]:g}-{x[-1]:g} nm but the grid spans "
                f"{self.start_nm:g}-{self.stop_nm:g} nm; pass outside='zero' or 'hold' "
                "to extend them"
            )

        q = np.clip(wl, x[0], x[-1])
        idx = np.clip(np.searchsorted(x, q, side="right") - 1, 0, x.size - 2)
        t = (q - x[idx]) / (x[idx + 1] - x[idx])
        # (1 - t) * y0 + t * y1 reproduces the samples exactly at t == 0 and t == 1.
        out = y[..., idx] * (1.0 - t) + y[..., idx + 1] * t
        if outside == "zero":
            out[..., out_of_range] = 0.0
        return np.moveaxis(out, -1, axis)

    def integrate(self, values, axis: int = -1):
        """Integrate band-integrated spectral ``values`` over wavelength: a sum over ``axis``."""
        try:
            xp = xp_of(values)
        except TypeError:
            values = np.asarray(values, dtype=np.float64)
            xp = xp_of(values)
        try:
            n = values.shape[axis]
        except IndexError as exc:
            raise ValueError(
                f"axis {axis} is out of range for values of shape {values.shape}"
            ) from exc
        if n != self.n_bands:
            raise ValueError(
                f"values have {n} entries along axis {axis}, but the grid has {self.n_bands} bands"
            )
        return xp.sum(values, axis=axis)


DEFAULT_GRID = SpectralGrid()


# --------------------------------------------------------------------------------------
# CIE data (colour-science)
# --------------------------------------------------------------------------------------

_CMF_NAME = "CIE 1931 2 Degree Standard Observer"
# D55 is the ISO sensitometric daylight base (spec D1 section 4.7, section 6.1).
_ILLUMINANT_NAMES = ("D65", "D55", "D50", "A", "E")


@lru_cache(maxsize=1)
def _colour():
    """Import colour-science lazily, silencing its notice about the optional matplotlib."""
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message=r'.*"Matplotlib" related API features are not available'
        )
        import colour
    return colour


def _read_only(array: np.ndarray) -> np.ndarray:
    array.flags.writeable = False
    return array


@lru_cache(maxsize=32)
def _cmf_cached(grid: SpectralGrid) -> np.ndarray:
    cmfs = _colour().MSDS_CMFS[_CMF_NAME]
    # The tabulated CIE 1931 functions span 360-830 nm and are negligible beyond it.
    return _read_only(grid.resample(cmfs.wavelengths, cmfs.values, outside="zero", axis=0))


def cmf_xyz(grid: SpectralGrid) -> np.ndarray:
    """CIE 1931 2 degree colour-matching functions on ``grid``, shape ``(B, 3)``.

    Columns are ``xbar, ybar, zbar`` sampled at the band centres (``ybar(555 nm) == 1``).
    """
    return _cmf_cached(grid).copy()


def _normalise_to_unit_luminance(banded: np.ndarray, grid: SpectralGrid, what: str) -> np.ndarray:
    ybar = _cmf_cached(grid)[:, 1]
    luminance = float(np.sum(banded * ybar))
    if not (math.isfinite(luminance) and luminance > 0.0):
        raise ValueError(f"{what} has no luminous content on {grid}; cannot normalise")
    return banded / luminance


@lru_cache(maxsize=64)
def _illuminant_cached(name: str, grid: SpectralGrid) -> np.ndarray:
    sd = _colour().SDS_ILLUMINANTS[name]
    try:
        spd = grid.resample(sd.wavelengths, sd.values, outside="error")
    except ValueError as exc:
        raise ValueError(f"illuminant {name!r}: {exc}") from exc
    banded = spd * grid.step_nm
    return _read_only(_normalise_to_unit_luminance(banded, grid, f"illuminant {name!r}"))


def illuminant(name: str, grid: SpectralGrid) -> np.ndarray:
    """Relative band-integrated energy of a CIE illuminant on ``grid``, shape ``(B,)``.

    ``name`` is one of ``"D65" "D55" "D50" "A" "E"``. The tabulated spectral power distribution is
    resampled to the band centres, multiplied by ``step_nm`` and normalised so that
    ``sum(spd * ybar) == 1``.
    """
    if not isinstance(name, str) or name not in _ILLUMINANT_NAMES:
        raise ValueError(f"unknown illuminant {name!r}; available: {', '.join(_ILLUMINANT_NAMES)}")
    return _illuminant_cached(name, grid).copy()


def blackbody(temperature_k: float, grid: SpectralGrid) -> np.ndarray:
    """Relative band-integrated energy of a Planckian radiator, shape ``(B,)``.

    Planck's law ``M(lambda) ~ lambda**-5 / (exp(c2 / (lambda T)) - 1)`` with
    ``c2 = h c / k``, multiplied by ``step_nm`` and normalised to ``sum(spd * ybar) == 1``.
    """
    if (
        isinstance(temperature_k, bool)
        or not isinstance(temperature_k, Real)
        or not math.isfinite(temperature_k)
        or temperature_k <= 0.0
    ):
        raise ValueError(f"temperature_k must be positive and finite, got {temperature_k!r}")
    wl = grid.wavelengths
    c2_nm_k = H_PLANCK * C_LIGHT / K_BOLTZMANN * 1e9
    x = c2_nm_k / (wl * float(temperature_k))
    # log(lambda**-5 / expm1(x)), written so that neither large nor small x overflows.
    log_spd = -5.0 * np.log(wl) - (x + np.log(-np.expm1(-x)))
    spd = np.exp(log_spd - log_spd.max())
    return _normalise_to_unit_luminance(spd * grid.step_nm, grid, f"blackbody at {temperature_k} K")


# --------------------------------------------------------------------------------------
# sRGB (IEC 61966-2-1)
# --------------------------------------------------------------------------------------

_SRGB_PRIMARIES_XY = ((0.64, 0.33), (0.30, 0.60), (0.15, 0.06))
_SRGB_WHITE_XY = (0.3127, 0.3290)  # D65


def _xy_to_xyz(xy: tuple[float, float]) -> np.ndarray:
    x, y = xy
    return np.array([x / y, 1.0, (1.0 - x - y) / y])


def _normalised_primary_matrix(primaries_xy, white_xy) -> np.ndarray:
    primaries = np.column_stack([_xy_to_xyz(c) for c in primaries_xy])
    scale = np.linalg.solve(primaries, _xy_to_xyz(white_xy))
    return primaries * scale


LINEAR_SRGB_TO_XYZ: np.ndarray = _read_only(
    _normalised_primary_matrix(_SRGB_PRIMARIES_XY, _SRGB_WHITE_XY)
)
"""Linear sRGB to CIE XYZ (``Y`` of RGB white is 1), shape ``(3, 3)``, for column vectors.

Derived at full precision from the IEC 61966-2-1 primaries and D65 white point. It rounds to
the four-decimal matrix printed in the standard; the printed forward and inverse tables are
not exact inverses of each other, so both matrices here come from the same derivation.
"""

XYZ_TO_LINEAR_SRGB: np.ndarray = _read_only(np.linalg.inv(LINEAR_SRGB_TO_XYZ))
"""CIE XYZ to linear sRGB, the exact inverse of :data:`LINEAR_SRGB_TO_XYZ`."""
