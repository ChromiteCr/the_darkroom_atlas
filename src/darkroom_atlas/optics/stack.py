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
from darkroom_atlas.data.models import BaseLayer, EmulsionLayer, Material
from darkroom_atlas.emulsion.classes import SizeClasses
from darkroom_atlas.spectral import SpectralGrid

__all__ = [
    "PHOTON_ENERGY_EV_NM",
    "BaseOptics",
    "OpticalStack",
    "base_optics",
    "interface_reflectances",
    "optical_stack",
]

PHOTON_ENERGY_EV_NM = 1239.84
"""``E_ph = 1239.84 / lambda`` in eV for ``lambda`` in nm (spec §2.1.2)."""

_SCATTER_REFERENCE_NM = 550.0
"""Reference wavelength of ``S_b = S_550 (550/lambda)^q_S`` (spec §2.1.3)."""


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


# --------------------------------------------------------------------------------------
# Builder
# --------------------------------------------------------------------------------------


def interface_reflectances(gelatin_index: float, base_index: float) -> tuple[float, float, float]:
    """Return ``(T_top, r_i, R_bot)`` of the two-stream boundaries (spec §2.1.3).

    ``R_top = ((n_gel - 1)/(n_gel + 1))^2``, ``T_top = 1 - R_top``,
    ``r_i = 1 - T_top / n_gel^2`` and ``R_bot = 1 - (n_b/n_gel)^2``. When ``n_b >= n_gel``
    there is no total reflection into the base and ``R_bot`` is 0 (spec §2.1.3 gives only
    the ``n_b < n_gel`` case). ``n_gel`` must be > 1 and ``n_b`` > 0.
    """
    n_gel = _number("interface_reflectances", "gelatin_index", gelatin_index)
    n_b = _number("interface_reflectances", "base_index", base_index)
    if n_gel <= 1.0:
        raise ValueError(f"gelatin_index must be greater than 1, got {n_gel}")
    if n_b <= 0.0:
        raise ValueError(f"base_index must be positive, got {n_b}")
    r_top = ((n_gel - 1.0) / (n_gel + 1.0)) ** 2
    t_top = 1.0 - r_top
    r_internal = 1.0 - t_top / n_gel**2
    r_bottom = max(0.0, 1.0 - (n_b / n_gel) ** 2)
    return t_top, r_internal, r_bottom


def _layers(material: Material) -> tuple[EmulsionLayer, BaseLayer]:
    if not isinstance(material, Material):
        raise TypeError(f"material must be a Material, got {type(material).__name__}")
    emulsions = material.emulsion_layers
    bases = [layer for layer in material.layers if isinstance(layer, BaseLayer)]
    if len(emulsions) != 1 or len(bases) != 1:
        raise ValueError(
            f"material {material.id!r}: an optical stack needs exactly one emulsion layer and "
            f"one base layer, found {len(emulsions)} and {len(bases)}"
        )
    return emulsions[0], bases[0]


def base_optics(material: Material, film_format: str) -> BaseOptics:
    """Base of ``material`` in ``film_format`` (spec §2.1.5, §2.1.8).

    ``ValueError`` if the format is not one of the base's formats or the base lacks
    ``refractive_index``.
    """
    _, base = _layers(material)
    formats = dict(base.formats)
    if film_format not in formats:
        raise ValueError(
            f"material {material.id!r}: unknown film format {film_format!r}; available: "
            f"{', '.join(formats) or '(none)'}"
        )
    if base.refractive_index is None:
        raise ValueError(f"material {material.id!r}: layer {base.id!r} has no refractive_index")
    fmt = formats[film_format]
    return BaseOptics(
        thickness_um=fmt.thickness.to("um"),
        refractive_index=base.refractive_index.to("1"),
        grey_density=fmt.grey_density.to("1"),
        backing_density=fmt.backing_density.to("1"),
    )


def optical_stack(
    material: Material, film_format: str, classes: SizeClasses, grid: SpectralGrid
) -> OpticalStack:
    """Build the optical stack of ``material`` in ``film_format`` on ``grid`` (spec §2.1).

    * Urbach edge: ``alpha_b = alpha_0`` for ``E_ph,b >= E_0``, otherwise
      ``alpha_0 exp((E_ph,b - E_0)/E_U)``, ``E_ph,b = 1239.84/lambda_b`` eV (§2.1.2).
    * Dye: ``s_b = sum_j A_j exp(-((nu_b - nu_j)/w_j)^2 / 2)``, ``nu_b = 1000/lambda_b``
      um^-1 (§2.1.2).
    * ``mu_a,b = sum_c (nu_c/L_e)(V_c alpha_b + S_c s_b)``, ``K_b = 2 mu_a,b``,
      ``S_b = S_550 (550/lambda_b)^q_S`` (§2.1.2, §2.1.3), ``L_e = classes.dry_thickness_um``.
    * Interfaces from :func:`interface_reflectances`; base from :func:`base_optics`.

    ``classes`` are the size classes of the material's emulsion layer
    (:func:`darkroom_atlas.emulsion.size_classes`). ``ValueError`` if a field the stack needs
    is missing.
    """
    emulsion, _ = _layers(material)
    if not isinstance(classes, SizeClasses):
        raise TypeError(f"classes must be SizeClasses, got {type(classes).__name__}")
    if not isinstance(grid, SpectralGrid):
        raise TypeError(f"grid must be a SpectralGrid, got {type(grid).__name__}")
    optics = emulsion.optics
    if optics is None:
        raise ValueError(f"material {material.id!r}: layer {emulsion.id!r} has no optics")
    base = base_optics(material, film_format)

    wavelengths = grid.wavelengths
    photon_ev = PHOTON_ENERGY_EV_NM / wavelengths
    intrinsic = optics.intrinsic_absorption
    alpha0 = intrinsic.alpha0.to("um-1")
    edge = intrinsic.edge_energy.to("eV")
    urbach = intrinsic.urbach_energy.to("eV")
    if not alpha0 >= 0.0 or not urbach > 0.0:
        raise ValueError(
            f"material {material.id!r}: intrinsic_absorption needs alpha0 >= 0 and "
            f"urbach_energy > 0, got {alpha0} and {urbach}"
        )
    alpha = alpha0 * np.exp(np.minimum(photon_ev - edge, 0.0) / urbach)

    wavenumber = 1000.0 / wavelengths
    dye = np.zeros(grid.n_bands)
    for j, band in enumerate(optics.dye_bands):
        amplitude = band.amplitude.to("1")
        peak = band.peak_wavenumber.to("um-1")
        width = band.width.to("um-1")
        if not amplitude >= 0.0 or not width > 0.0:
            raise ValueError(
                f"material {material.id!r}: dye_bands[{j}] needs amplitude >= 0 and width > 0, "
                f"got {amplitude} and {width}"
            )
        dye += amplitude * np.exp(-0.5 * ((wavenumber - peak) / width) ** 2)

    thickness = classes.dry_thickness_um
    volume_density = float(np.sum(classes.number_per_um2 * classes.volume_um3)) / thickness
    surface_density = float(np.sum(classes.number_per_um2 * classes.surface_um2)) / thickness
    mu_a = volume_density * alpha + surface_density * dye

    s550 = optics.scattering_550.to("um-1")
    q_s = optics.scattering_exponent.to("1")
    scattering = s550 * (_SCATTER_REFERENCE_NM / wavelengths) ** q_s

    t_top, r_internal, r_bottom = interface_reflectances(
        optics.gelatin_index.to("1"), base.refractive_index
    )
    return OpticalStack(
        grid=grid,
        emulsion_thickness_um=thickness,
        alpha_per_um=alpha,
        dye_cross_section=dye,
        mu_a_per_um=mu_a,
        km_k_per_um=2.0 * mu_a,
        km_s_per_um=scattering,
        t_top=t_top,
        r_internal_top=r_internal,
        r_bottom=r_bottom,
        base=base,
        psf_weight=optics.psf_weight.to("1"),
        psf_length_um=optics.psf_length.to("um"),
    )
