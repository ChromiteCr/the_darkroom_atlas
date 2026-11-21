"""Grain size classes of an emulsion layer (spec D1 2.0.2).

The equivalent circular diameter (ECD) is log-normal by number, ``ln d ~ N(ln d_50,
sigma_d^2)``, truncated at +-3 sigma and renormalised. The range is cut into ``K`` classes of
equal width in ``ln d``; class ``c`` is represented by the geometric centre of its bounds.
``conventional`` (and ``cubic``) grains are spheres of diameter ``d_c``; ``tabular`` grains are
flat discs of diameter ``d_c`` and thickness ``t_g``.
"""

from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.special import ndtr

from darkroom_atlas.data.models import EmulsionLayer

__all__ = [
    "M_AG",
    "M_AGBR",
    "M_AGI",
    "RHO_AGX",
    "RHO_GELATIN",
    "TRUNCATION_SIGMA",
    "SizeClasses",
    "size_classes",
]

M_AG = 107.87
"""Molar mass of silver, g/mol (spec D1 2.0.2)."""

M_AGBR = 187.77
"""Molar mass of AgBr, g/mol (spec D1 2.0.2)."""

M_AGI = 234.77
"""Molar mass of AgI, g/mol, from atomic masses (spec D1 2.0.2, assumed)."""

RHO_AGX = 6.473
"""Density of the silver halide, g/cm3: AgBr, iodide ignored (spec D1 2.0.2; LAT 2.3)."""

RHO_GELATIN = 1.35
"""Dry gelatin density, g/cm3, for the dry layer thickness (LAT 2.3, recalled -> assumed).

Spec D1 2.0.2 uses rho_gel without a value; this is the value behind its 4.1 um example.
"""

TRUNCATION_SIGMA = 3.0
"""The size distribution is truncated at +-3 sigma and renormalised (spec D1 2.0.2)."""


@dataclass(frozen=True, eq=False)
class SizeClasses:
    """Per-class grain geometry and numbers of one emulsion layer.

    All arrays are read-only, float64, shape ``(K,)``. Lengths in um, areas in um2,
    volumes in um3, ``number_per_um2`` per um2 of film, ``silver_g_per_grain`` in g of Ag.
    ``thickness_um`` is ``t_g`` for tabular grains and ``d_c`` for spheres.
    ``dry_thickness_um`` is the dry emulsion thickness ``L_e``.
    """

    morphology: Literal["conventional", "tabular", "cubic"]  # cubic is treated as conventional
    ecd_um: np.ndarray
    thickness_um: np.ndarray
    volume_um3: np.ndarray
    surface_um2: np.ndarray
    projected_um2: np.ndarray
    number_fraction: np.ndarray
    number_per_um2: np.ndarray
    silver_g_per_grain: np.ndarray
    iodide_fraction: float
    dry_thickness_um: float
    silver_coverage_g_m2: float

    @property
    def count(self) -> int:
        """Number of classes ``K``."""
        return int(self.ecd_um.shape[0])


def _read_only(array: np.ndarray) -> np.ndarray:
    array = np.ascontiguousarray(array, dtype=np.float64)
    array.flags.writeable = False
    return array


def _scalar(quantity, unit: str, name: str) -> float:
    value = quantity.to(unit)
    if isinstance(value, np.ndarray):
        raise ValueError(f"{name} must be a single value, got a list")  # noqa: TRY004
    return value


def size_classes(layer: EmulsionLayer, n_classes: int | None = None) -> SizeClasses:
    """Size classes of ``layer`` (spec D1 2.0.2).

    ``n_classes`` overrides ``layer.grain.n_classes`` (spec D1 3.1 uses 6 for calibration
    and images, 12 for the reference profile). A missing ``iodide_fraction`` counts as 0
    (pure AgBr). Raises ``ValueError`` if the layer has no grain or no ``gelatin_coverage``,
    or if a parameter is out of range.
    """
    grain = layer.grain
    if grain is None:
        raise ValueError(f"emulsion layer {layer.id!r}: grain is required for size classes")
    if layer.gelatin_coverage is None:
        raise ValueError(
            f"emulsion layer {layer.id!r}: gelatin_coverage is required for size classes"
        )
    k = grain.n_classes if n_classes is None else n_classes
    if isinstance(k, bool) or not isinstance(k, (int, np.integer)) or not 2 <= k <= 64:
        # The same range as Grain.n_classes (spec D1 5.4).
        raise ValueError(f"n_classes must be an integer from 2 to 64, got {k!r}")
    k = int(k)

    d50 = _scalar(grain.mean_diameter, "um", "grain.mean_diameter")
    sigma = _scalar(grain.sigma_log, "1", "grain.sigma_log")
    m_coat = _scalar(grain.silver_coverage, "g/m2", "grain.silver_coverage")
    m_gel = _scalar(layer.gelatin_coverage, "g/m2", "gelatin_coverage")
    x_i = (
        0.0
        if grain.iodide_fraction is None
        else _scalar(grain.iodide_fraction, "1", "grain.iodide_fraction")
    )
    if not d50 > 0.0:
        raise ValueError(f"grain.mean_diameter must be positive, got {d50}")
    if not sigma > 0.0:
        raise ValueError(f"grain.sigma_log must be positive, got {sigma}")
    if not m_coat > 0.0:
        raise ValueError(f"grain.silver_coverage must be positive, got {m_coat}")
    if m_gel < 0.0:
        raise ValueError(f"gelatin_coverage must not be negative, got {m_gel}")

    # Class bounds in standardised log diameter t = (ln d - mu_d) / sigma_d.
    t_edges = -TRUNCATION_SIGMA + 2.0 * TRUNCATION_SIGMA * np.arange(k + 1) / k
    cdf = ndtr(t_edges)
    number_fraction = np.diff(cdf) / (cdf[-1] - cdf[0])
    ln_d50 = np.log(d50)
    ecd = np.exp(ln_d50 + sigma * 0.5 * (t_edges[:-1] + t_edges[1:]))

    projected = np.pi * ecd**2 / 4.0
    if grain.morphology == "tabular":
        t_g = _scalar(grain.thickness, "um", "grain.thickness")
        if not t_g > 0.0:
            raise ValueError(f"grain.thickness must be positive, got {t_g}")
        thickness = np.full(k, t_g)
        volume = np.pi * ecd**2 * t_g / 4.0
        surface = np.pi * ecd**2 / 2.0 + np.pi * ecd * t_g
    else:
        thickness = ecd.copy()
        volume = np.pi * ecd**3 / 6.0
        surface = np.pi * ecd**2

    # Silver bookkeeping (spec D1 2.0.2): 1 g/cm3 = 1e-12 g/um3, 1 g/m2 = 1e-12 g/um2.
    m_agx = (1.0 - x_i) * M_AGBR + x_i * M_AGI
    silver_per_grain = 1e-12 * RHO_AGX * volume * M_AG / m_agx
    n_total = 1e-12 * m_coat / np.sum(number_fraction * silver_per_grain)
    number_per_um2 = number_fraction * n_total
    # g/m2 divided by g/cm3 is um.
    dry_thickness = m_coat * (m_agx / M_AG) / RHO_AGX + m_gel / RHO_GELATIN

    return SizeClasses(
        morphology=grain.morphology,
        ecd_um=_read_only(ecd),
        thickness_um=_read_only(thickness),
        volume_um3=_read_only(volume),
        surface_um2=_read_only(surface),
        projected_um2=_read_only(projected),
        number_fraction=_read_only(number_fraction),
        number_per_um2=_read_only(number_per_um2),
        silver_g_per_grain=_read_only(silver_per_grain),
        iodide_fraction=float(x_i),
        dry_thickness_um=float(dry_thickness),
        silver_coverage_g_m2=float(m_coat),
    )
