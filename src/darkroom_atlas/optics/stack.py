"""Optical description of a film for the exposure stage (spec §2.1).

:class:`BaseOptics` holds the film base of one format; :class:`OpticalStack` holds everything
the light-transport code needs: the per-band absorption and scattering of the dry emulsion
layer, the interface reflectances of the two-stream model (§2.1.3) and the scatter PSF
parameters (§2.1.6). Spectral arrays have shape ``(B,)`` on the stack's grid, are float64,
finite, non-negative and stored read-only.
"""

import math
from dataclasses import dataclass
from numbers import Real

import numpy as np

from darkroom_atlas._frozen import readonly_array
from darkroom_atlas.spectral import SpectralGrid

__all__ = ["BaseOptics", "OpticalStack"]


def _number(owner: str, name: str, value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
        raise ValueError(f"{owner}.{name} must be a finite number, got {value!r}")
    return float(value)


@dataclass(frozen=True)
class BaseOptics:
    """Film base of one format (spec §2.1.5, §2.1.8).

    ``thickness_um`` > 0; ``refractive_index`` > 1 (the critical angle of the base-air
    interface must exist); ``grey_density`` and ``backing_density`` >= 0 (decadic, for one
    normal pass). Anything else raises ``ValueError``.
    """

    thickness_um: float
    refractive_index: float
    grey_density: float
    backing_density: float

    def __post_init__(self) -> None:
        for name in ("thickness_um", "refractive_index", "grey_density", "backing_density"):
            object.__setattr__(self, name, _number("BaseOptics", name, getattr(self, name)))
        if self.thickness_um <= 0.0:
            raise ValueError(f"BaseOptics.thickness_um must be positive, got {self.thickness_um}")
        if self.refractive_index <= 1.0:
            raise ValueError(
                f"BaseOptics.refractive_index must be greater than 1, got {self.refractive_index}"
            )
        for name in ("grey_density", "backing_density"):
            if getattr(self, name) < 0.0:
                raise ValueError(
                    f"BaseOptics.{name} must be non-negative, got {getattr(self, name)}"
                )

    @property
    def critical_angle_rad(self) -> float:
        """Critical angle of the base-air interface, ``arcsin(1/n_b)`` (spec §2.1.5)."""
        return math.asin(1.0 / self.refractive_index)

    @property
    def critical_radius_um(self) -> float:
        """Return radius of the critical ray, ``r_c = 2 t_b tan(theta_c)`` (spec §2.1.5)."""
        return 2.0 * self.thickness_um * math.tan(self.critical_angle_rad)


_SPECTRAL_FIELDS = (
    "alpha_per_um",
    "dye_cross_section",
    "mu_a_per_um",
    "km_k_per_um",
    "km_s_per_um",
)


@dataclass(frozen=True, eq=False)
class OpticalStack:
    """Per-band optics of the dry emulsion layer plus its base (spec §2.1.2-§2.1.6).

    * ``alpha_per_um``: intrinsic AgX absorption ``alpha_b`` (um^-1).
    * ``dye_cross_section``: dye absorption per unit grain surface ``s_b`` (1).
    * ``mu_a_per_um``: layer absorption coefficient ``mu_a,b`` (um^-1).
    * ``km_k_per_um``: two-stream absorption ``K_b = 2 mu_a,b``; checked against
      ``mu_a_per_um``.
    * ``km_s_per_um``: two-stream scattering ``S_b`` (um^-1).
    * ``t_top``: transmittance of the top surface for the incident light, in (0, 1].
    * ``r_internal_top``, ``r_bottom``: reflectances ``r_i`` and ``R_bot`` for diffuse light
      inside the layer, in [0, 1).
    * ``psf_weight`` in [0, 1] and ``psf_length_um`` > 0: scatter PSF ``w_P`` and ``a_P``.

    Invalid values raise ``ValueError``.
    """

    grid: SpectralGrid
    emulsion_thickness_um: float
    alpha_per_um: np.ndarray
    dye_cross_section: np.ndarray
    mu_a_per_um: np.ndarray
    km_k_per_um: np.ndarray
    km_s_per_um: np.ndarray
    t_top: float
    r_internal_top: float
    r_bottom: float
    base: BaseOptics
    psf_weight: float
    psf_length_um: float

    def __post_init__(self) -> None:
        if not isinstance(self.grid, SpectralGrid):
            raise TypeError(f"grid must be a SpectralGrid, got {type(self.grid).__name__}")
        if not isinstance(self.base, BaseOptics):
            raise TypeError(f"base must be a BaseOptics, got {type(self.base).__name__}")
        n_bands = self.grid.n_bands
        for name in _SPECTRAL_FIELDS:
            try:
                array = readonly_array(getattr(self, name), dtype=np.float64)
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"OpticalStack.{name} must be a real numeric array: {exc}"
                ) from exc
            if array.shape != (n_bands,):
                raise ValueError(
                    f"OpticalStack.{name} must have shape ({n_bands},), got {array.shape}"
                )
            if not np.all(np.isfinite(array)):
                raise ValueError(f"OpticalStack.{name} must be finite (found NaN or inf)")
            if np.any(array < 0.0):
                raise ValueError(f"OpticalStack.{name} must be non-negative")
            object.__setattr__(self, name, array)
        if not np.allclose(self.km_k_per_um, 2.0 * self.mu_a_per_um, rtol=1e-12, atol=0.0):
            raise ValueError("OpticalStack.km_k_per_um must equal 2 * mu_a_per_um (spec §2.1.3)")
        for name in (
            "emulsion_thickness_um",
            "t_top",
            "r_internal_top",
            "r_bottom",
            "psf_weight",
            "psf_length_um",
        ):
            object.__setattr__(self, name, _number("OpticalStack", name, getattr(self, name)))
        if self.emulsion_thickness_um <= 0.0:
            raise ValueError(
                f"OpticalStack.emulsion_thickness_um must be positive, "
                f"got {self.emulsion_thickness_um}"
            )
        if not 0.0 < self.t_top <= 1.0:
            raise ValueError(f"OpticalStack.t_top must be in (0, 1], got {self.t_top}")
        for name in ("r_internal_top", "r_bottom"):
            if not 0.0 <= getattr(self, name) < 1.0:
                raise ValueError(
                    f"OpticalStack.{name} must be in [0, 1), got {getattr(self, name)}"
                )
        if not 0.0 <= self.psf_weight <= 1.0:
            raise ValueError(f"OpticalStack.psf_weight must be in [0, 1], got {self.psf_weight}")
        if self.psf_length_um <= 0.0:
            raise ValueError(
                f"OpticalStack.psf_length_um must be positive, got {self.psf_length_um}"
            )
