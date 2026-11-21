"""Silver image of a developed ``bw-mechanistic`` film (spec D1 2.5).

Nutting's law per size class with a size-dependent covering power:

* ``CP_c = CP_ref (d_c / d_ref)^(-beta_CP)``, ``d_ref = 1 um``;
* ``a_c = ln10 * 1e12 * CP_c * m_Ag,c`` (um^2), the effective projected area of one fully
  developed grain of class ``c``;
* ``D_Ag = sum_c CP_c silver_c + r_CP CP_ref silver_physical``.

The spectral transmittance is ``T_b = 10^-(D_base,b + D_Ag n_Ag,b + CP_hal halide +
D_dye,b residual_dye)`` with neutral silver (``n_Ag,b = 1``). A density for a spectral
product ``P_b`` is ``-log10(sum_b P_b T_b / sum_b P_b)``; because silver and halide are
neutral it factorises into a per-pixel silver and halide term and a spectral term of the
base and the dye, so no ``(h, w, B)`` array is built (:func:`product_density`).

Sums over size classes and bands are written as element-wise multiply-adds (spec 3.6), so a
pixel's value does not depend on the window it is computed in.
"""

import math
from dataclasses import dataclass
from numbers import Integral

import numpy as np

from darkroom_atlas.data.models import EmulsionDevelopment, EmulsionLayer
from darkroom_atlas.emulsion.classes import SizeClasses, size_classes
from darkroom_atlas.state import EmulsionState

__all__ = [
    "D_REF_UM",
    "LN10",
    "ImagingConstants",
    "TransmittanceFactors",
    "chemical_silver_density",
    "class_silver_g_m2",
    "covering_powers",
    "developed_areas_um2",
    "imaging_constants",
    "product_density",
    "resolve_window",
    "silver_density",
    "transmittance_factors",
]

LN10 = math.log(10.0)
"""ln 10: ``T = exp(-N a)`` is ``D = N a / ln 10`` (Nutting, spec 2.5.1)."""

D_REF_UM = 1.0
"""Reference diameter of the covering-power size law, um (spec 2.5.1)."""

# 1 g/um^2 = 1e12 g/m^2 (spec 2.0.2, 2.4.10).
_UM2_PER_M2 = 1e12

_LAYER = "emulsion"


@dataclass(frozen=True, eq=False)
class ImagingConstants:
    """Per-class imaging constants of a developed state (spec 2.5.1, 2.5.3).

    Arrays are read-only, shape ``(K,)``. ``class_silver_g_m2`` is the silver of class ``c``
    when every grain is fully developed, ``1e12 m_Ag,c nu_c`` (g/m2), so the developed
    fraction is ``X_c = silver_c / class_silver_g_m2``.
    """

    classes: SizeClasses
    covering_power: np.ndarray
    developed_areas_um2: np.ndarray
    class_silver_g_m2: np.ndarray
    covering_power_ref: float
    physical_covering_ratio: float
    halide_density: float


@dataclass(frozen=True, eq=False)
class TransmittanceFactors:
    """The non-silver factors of the spectral transmittance of a window (spec 2.5.2).

    ``base_density`` and ``dye_density`` are per band ``(B,)``; ``halide_density`` is the
    neutral density ``CP_hal * halide`` of the remaining silver halide and ``dye_residual``
    the residual-dye fraction, both ``(h, w)``.
    """

    base_density: np.ndarray
    dye_density: np.ndarray
    halide_density: np.ndarray
    dye_residual: np.ndarray


def _read_only(array: np.ndarray) -> np.ndarray:
    array = np.ascontiguousarray(array, dtype=np.float64)
    array.flags.writeable = False
    return array


def _scalar(quantity, unit: str, name: str) -> float:
    value = quantity.to(unit)
    if isinstance(value, np.ndarray):
        raise ValueError(f"{name} must be a single value, got a list")  # noqa: TRY004
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite, got {value}")
    return value


def covering_powers(classes: SizeClasses, dev: EmulsionDevelopment) -> np.ndarray:
    """Covering power ``CP_c = CP_ref (d_c / 1 um)^(-beta_CP)`` per class, 1/(g/m2) (spec 2.5.1).

    Returns a read-only ``(K,)`` array. Raises ``ValueError`` if ``CP_ref`` is not positive.
    """
    cp_ref = _scalar(dev.covering_power, "1/(g/m2)", "development.covering_power")
    beta = _scalar(
        dev.covering_power_size_exponent, "1", "development.covering_power_size_exponent"
    )
    if not cp_ref > 0.0:
        raise ValueError(f"development.covering_power must be positive, got {cp_ref}")
    return _read_only(cp_ref * (np.asarray(classes.ecd_um) / D_REF_UM) ** (-beta))


def developed_areas_um2(classes: SizeClasses, covering_power: np.ndarray) -> np.ndarray:
    """Effective area ``a_c = ln10 * 1e12 * CP_c * m_Ag,c`` of a fully developed grain, um^2.

    ``covering_power`` is ``(K,)`` in 1/(g/m2) (spec 2.5.1). Returns a read-only ``(K,)``
    array. Raises ``ValueError`` on a shape mismatch or negative or non-finite values.
    """
    cp = np.asarray(covering_power, dtype=np.float64)
    if cp.shape != (classes.count,):
        raise ValueError(f"covering_power must have shape ({classes.count},), got {cp.shape}")
    if not np.all(np.isfinite(cp)) or np.any(cp < 0.0):
        raise ValueError("covering_power must be finite and non-negative")
    return _read_only(LN10 * _UM2_PER_M2 * cp * np.asarray(classes.silver_g_per_grain))


def class_silver_g_m2(classes: SizeClasses) -> np.ndarray:
    """Silver of each class when fully developed, ``1e12 m_Ag,c nu_c``, g/m2 (spec 2.4.10)."""
    return _read_only(
        _UM2_PER_M2 * np.asarray(classes.silver_g_per_grain) * np.asarray(classes.number_per_um2)
    )


def _emulsion_layer(state: EmulsionState) -> EmulsionLayer:
    layers = [layer for layer in state.material.emulsion_layers if layer.id == _LAYER]
    if len(layers) != 1:
        raise ValueError(
            f"material {state.material.id!r} has no emulsion layer {_LAYER!r}, which the "
            "bw-mechanistic imaging needs"
        )
    layer = layers[0]
    if layer.development is None or layer.grain is None:
        raise ValueError(
            f"material {state.material.id!r}: emulsion layer {_LAYER!r} needs grain and "
            "development parameters for imaging"
        )
    return layer


def imaging_constants(state: EmulsionState) -> ImagingConstants:
    """Imaging constants of a developed state.

    The number of size classes ``K`` is the trailing dimension of ``emulsion/silver``; the
    classes are rebuilt from the material with that ``K``. Raises ``KeyError`` if the state
    has no ``emulsion/silver`` and ``ValueError`` if its shape is not ``(H, W, K)``.
    """
    silver = state.get(_LAYER, "silver")
    if silver.ndim != 3:
        raise ValueError(f"emulsion/silver must have shape (H, W, K), got {silver.shape}")
    layer = _emulsion_layer(state)
    dev = layer.development
    classes = size_classes(layer, int(silver.shape[2]))
    cp = covering_powers(classes, dev)
    r_cp = _scalar(dev.physical_covering_ratio, "1", "development.physical_covering_ratio")
    cp_hal = _scalar(dev.halide_density, "1/(g/m2)", "development.halide_density")
    if r_cp < 0.0 or cp_hal < 0.0:
        raise ValueError(
            "development.physical_covering_ratio and development.halide_density must not be "
            f"negative, got {r_cp} and {cp_hal}"
        )
    return ImagingConstants(
        classes=classes,
        covering_power=cp,
        developed_areas_um2=developed_areas_um2(classes, cp),
        class_silver_g_m2=class_silver_g_m2(classes),
        covering_power_ref=_scalar(dev.covering_power, "1/(g/m2)", "development.covering_power"),
        physical_covering_ratio=r_cp,
        halide_density=cp_hal,
    )


def resolve_window(
    state: EmulsionState, rows: slice | None, cols: slice | None
) -> tuple[slice, slice]:
    """Resolve ``rows``/``cols`` (``None``: everything) into explicit, non-empty slices."""
    out = []
    for name, s, size in zip(("rows", "cols"), (rows, cols), state.geometry.shape_hw, strict=True):
        if s is None:
            s = slice(None)
        if not isinstance(s, slice):
            raise TypeError(f"{name} must be a slice or None, got {s!r}")
        if s.step not in (None, 1):
            raise ValueError(f"{name} slice must have a step of 1, got {s!r}")
        start = 0 if s.start is None else s.start
        stop = size if s.stop is None else s.stop
        for bound in (start, stop):
            if isinstance(bound, bool) or not isinstance(bound, Integral) or not 0 <= bound <= size:
                raise ValueError(f"{name} slice {s!r} must lie within 0..{size}")
        if stop <= start:
            raise ValueError(f"{name} slice {s!r} is empty")
        out.append(slice(int(start), int(stop)))
    return out[0], out[1]


def _field(state: EmulsionState, name: str, window: tuple[slice, slice]) -> np.ndarray:
    value = np.asarray(state.get(_LAYER, name)[window], dtype=np.float64)
    if not np.all(np.isfinite(value)):
        raise ValueError(f"emulsion/{name} must be finite")
    if np.any(value < 0.0):
        raise ValueError(f"emulsion/{name} must not be negative")
    return value


def chemical_silver_density(silver: np.ndarray, covering_power: np.ndarray) -> np.ndarray:
    """``sum_c CP_c silver_c`` for ``silver`` of shape ``(..., K)``, as element-wise adds."""
    out = covering_power[0] * silver[..., 0]
    for c in range(1, silver.shape[-1]):
        out = out + covering_power[c] * silver[..., c]
    return out


def silver_density(
    state: EmulsionState, rows: slice | None = None, cols: slice | None = None
) -> np.ndarray:
    """Silver density ``D_Ag = sum_c CP_c silver_c + r_CP CP_ref silver_physical`` (spec 2.5.1).

    Over the window ``rows`` x ``cols`` of the state's arrays (``None``: all). Returns a
    read-only ``(h, w)`` array. Needs ``emulsion/silver`` and ``emulsion/silver_physical``.
    """
    window = resolve_window(state, rows, cols)
    const = imaging_constants(state)
    silver = _field(state, "silver", window)
    physical = _field(state, "silver_physical", window)
    if physical.ndim != 2:
        raise ValueError(f"emulsion/silver_physical must have shape (H, W), got {physical.shape}")
    density = chemical_silver_density(silver, const.covering_power)
    density = density + const.physical_covering_ratio * const.covering_power_ref * physical
    return _read_only(density)


def _film_format(state: EmulsionState) -> str:
    exposure = state.meta.get("exposure")
    film_format = None
    if exposure is not None:
        film_format = exposure.get("film_format")
    if film_format is None:
        film_format = state.material.default_format
    formats = state.material.formats
    if film_format not in formats:
        raise ValueError(
            f"material {state.material.id!r}: unknown film format {film_format!r}; available: "
            f"{', '.join(formats) or '(none)'}"
        )
    return film_format


def transmittance_factors(state: EmulsionState, rows: slice, cols: slice) -> TransmittanceFactors:
    """Base, halide and residual-dye factors of the transmittance of a window (spec 2.5.2).

    ``D_base,b`` is the grey density of the run's film format (``meta["exposure"]
    ["film_format"]``, else the material's default), neutral; the 120 backing is removed
    in processing and contributes 0 (spec 2.5.2). ``D_dye,b`` is
    ``development.dye_density``, one value for all bands or one per band.
    """
    window = resolve_window(state, rows, cols)
    const = imaging_constants(state)
    n_bands = state.grid.n_bands
    base_format = state.material.formats[_film_format(state)]
    grey = _scalar(base_format.grey_density, "1", "base grey_density")
    base = np.full(n_bands, grey, dtype=np.float64)

    dev = _emulsion_layer(state).development
    dye = np.asarray(dev.dye_density.to("1"), dtype=np.float64)
    if dye.ndim == 0:
        dye = np.full(n_bands, float(dye))
    elif dye.shape != (n_bands,):
        raise ValueError(
            f"development.dye_density must be one value or {n_bands} values (one per band), "
            f"got {dye.shape[0]}"
        )
    if not np.all(np.isfinite(dye)):
        raise ValueError("development.dye_density must be finite")

    halide = _field(state, "halide", window)
    residual = _field(state, "residual_dye", window)
    if halide.ndim != 2 or residual.ndim != 2:
        raise ValueError("emulsion/halide and emulsion/residual_dye must have shape (H, W)")
    halide_density = const.halide_density * halide
    halide_density.flags.writeable = False
    residual.flags.writeable = False  # a view of the read-only field; may be a broadcast
    return TransmittanceFactors(
        base_density=_read_only(base),
        dye_density=_read_only(dye),
        halide_density=halide_density,
        dye_residual=residual,
    )


def _is_uniform(array: np.ndarray) -> bool:
    if array.strides[:2] == (0, 0):
        return True
    first = array.flat[0]
    return bool(np.all(array == first))


def product_density(
    silver_transmittance: np.ndarray, factors: TransmittanceFactors, product: np.ndarray
) -> np.ndarray:
    """Density ``-log10(sum_b P_b T_b / sum_b P_b)`` of a window for spectral product ``P``.

    ``silver_transmittance`` is the neutral silver transmittance ``(h, w)`` (``10^-D_Ag``, or a
    statistical-mode sample); with neutral silver and halide the density is
    ``-log10(T_Ag) + CP_hal halide + D_spec``, ``D_spec = -log10(sum_b p_b
    10^-(D_base,b + D_dye,b r))`` with ``p = P / sum P`` and ``r`` the residual dye. ``D_spec``
    is accumulated band by band, on one pixel when the dye is neutral or uniform over the
    window, else on ``(h, w)``; no ``(h, w, B)`` array is built. A zero
    transmittance gives ``inf``. Returns a read-only ``(h, w)`` array.
    """
    t = np.asarray(silver_transmittance, dtype=np.float64)
    if not isinstance(factors, TransmittanceFactors):
        raise TypeError(f"factors must be TransmittanceFactors, got {type(factors).__name__}")
    p = np.asarray(product, dtype=np.float64)
    base, dye = factors.base_density, factors.dye_density
    if p.shape != base.shape:
        raise ValueError(f"product must have shape {base.shape}, got {p.shape}")
    if not np.all(np.isfinite(p)) or np.any(p < 0.0):
        raise ValueError("product must be finite and non-negative")
    total = float(np.sum(p))
    if not total > 0.0:
        raise ValueError("product has no weight")
    if t.shape != factors.halide_density.shape:
        raise ValueError(
            f"silver_transmittance must have shape {factors.halide_density.shape}, got {t.shape}"
        )
    if not np.all(np.isfinite(t)) or np.any((t < 0.0) | (t > 1.0)):
        raise ValueError("silver_transmittance must be finite and lie in [0, 1]")
    weights = p / total
    residual = factors.dye_residual
    bands = [b for b in range(weights.shape[0]) if weights[b] > 0.0]

    if not np.any(dye[bands] != 0.0) or _is_uniform(residual):
        # One pixel's worth, by the same element-wise arithmetic as the general case, so a
        # pixel gets the same bits whichever case its window takes.
        r = residual[:1, :1] if np.any(dye[bands] != 0.0) else np.zeros((1, 1))
    else:
        r = residual
    acc = np.zeros(r.shape, dtype=np.float64)
    for b in bands:
        # exp, not 10.0 ** array: ufunc power may round differently by array length.
        acc = acc + weights[b] * np.exp(-LN10 * (base[b] + dye[b] * r))
    d_spec = -np.log10(acc)

    with np.errstate(divide="ignore"):
        density = -np.log10(t)
    density = density + factors.halide_density + d_spec
    return _read_only(density)
