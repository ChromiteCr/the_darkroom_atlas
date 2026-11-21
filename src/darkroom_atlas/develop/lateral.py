"""Exact lateral diffusion step ``L(tau)`` of the development split (spec D1 3.2, 3.5).

Lateral diffusion is linear with constant coefficients per cell level and component, so each
Fourier/cosine mode is advanced exactly::

    C_hat(k, t + tau) = C_hat(k, t) * exp(-D_{l,k} |k|^2 tau)

Boundaries (``boundary``):

- ``"neumann"``: canonical image blocks; DCT-II (orthonormal, ``scipy.fft.dctn``) over the
  block, zero flux at its outer edge. Cell-centred mode ``j`` of ``n`` cells of width ``dx``
  has the continuous wavenumber ``k_j = pi j / (n dx)``.
- ``"periodic"``: MTF and edge strips; real FFT, ``k = 2 pi m / (n dx)``. Exact for sinusoids.
- ``"pair"``: small-signal MTF column pairs (calibration profile, spec D1 3.5); columns
  ``a0, b0, a1, b1, ...`` represent the +1 and -1 phases of ``cos(k x)``; only the difference
  mode decays: ``C_a,b <- Cbar +- (C_a - C_b)/2 * exp(-D k^2 tau)``.

The decay factor is separable, ``exp(-D ky^2 tau) * exp(-D kx^2 tau)``, so only 1-D factors
are exponentiated. Mode 0 has factor 1 exactly, so the block total is conserved to rounding.
Transforms have the block's fixed shape and run single-threaded; the kernels use no BLAS and
no reduction, so the result is bit-identical whatever the memory alignment of the input
(spec D1 3.6). The caller supplies the gel diffusivities at the step midpoint
(``r_k D_k s(t)``, spec D1 3.5); liquid ones are time independent. Kernels and the clocks
``Lambda_z``, ``D_z`` are not passed here (they do not diffuse laterally).

TOTH (spec D1 3.5, revised). The mobile proton excess ``TOTH_mob = TOTH - T_gel nbar_gel(pH)``
(pH solved once at the start of the L step) is diffused instead of TOTH, with one diffusivity
per level::

    D_toth,l = (sum_f D_f,l bbar_f,l + D_H,l bbar_H,l + D_OH,l bbar_OH,l) / bbar_mob,l

(``bbar`` the mean over the block's cells, ``beta_mob`` without gelatin); then ``TOTH <- TOTH
+ (TOTH_mob' - TOTH_mob)``. Gelatin-bound protons stay put, TOTH is conserved and the number of
transforms does not change (TOTH_mob takes TOTH's slot in the same transform).

Public functions::

    lateral_diffuse(conc, diffusion_um2_s, tau_s, *, boundary, spacing_um=None,
                    wavenumber_per_um=None, toth=None) -> ndarray (same shape as conc)
        # grid: conc (ny, nx, *trailing), diffusion broadcast to trailing (e.g. (L, K))
        # pair: conc (2 n_pairs, *trailing), wavenumber_per_um (n_pairs,)
    lateral_factors(n, spacing_um, boundary, diffusion_um2_s, tau_s) -> ndarray (n_modes, ...)
    carrier_capacities(spec, totals, ph, mobile) -> ndarray (n_mobile + 2, *cells)
    toth_lateral_diffusion(beta, diffusion_um2_s, *, inner=None) -> ndarray (L,)
    TothLateral(index, gel_bound, diffusion_um2_s)
"""

import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy import fft

from darkroom_atlas.chemistry.speciation import Speciation, fractions, water_ions

__all__ = [
    "BOUNDARIES",
    "TothLateral",
    "carrier_capacities",
    "lateral_diffuse",
    "lateral_factors",
    "toth_lateral_diffusion",
]

Boundary = Literal["neumann", "periodic", "pair"]
BOUNDARIES = ("neumann", "periodic", "pair")
"""Lateral boundary kinds (spec D1 3.5; ``Lateral.boundary`` of spec D1 10)."""

_LN10 = math.log(10.0)


@dataclass(frozen=True, eq=False)
class TothLateral:
    """TOTH treatment of spec D1 3.5 for :func:`lateral_diffuse`.

    ``index``: position of TOTH on the last (component) axis of ``conc``. ``gel_bound``:
    ``T_gel nbar_gel(pH)`` at the start of the step, shape ``conc.shape[:-1]`` (0 in the
    liquid). ``diffusion_um2_s``: ``D_toth`` per level, shape ``conc.shape[-2:-1]`` (``(L,)``),
    from :func:`toth_lateral_diffusion`.
    """

    index: int
    gel_bound: np.ndarray
    diffusion_um2_s: np.ndarray


def _check_tau(tau_s: float) -> float:
    tau = float(tau_s)
    if not (math.isfinite(tau) and tau >= 0.0):
        raise ValueError(f"tau_s must be finite and >= 0, got {tau_s}")
    return tau


def _check_diffusion(diffusion, trailing: tuple[int, ...]) -> np.ndarray:
    d = np.asarray(diffusion, dtype=np.float64)
    try:
        d = np.broadcast_to(d, trailing)
    except ValueError:
        raise ValueError(
            f"diffusion shape {d.shape} does not broadcast to the trailing shape {trailing}"
        ) from None
    if not (np.all(np.isfinite(d)) and np.all(d >= 0.0)):
        raise ValueError("diffusion must be finite and >= 0 um^2/s")
    return d


def _wavenumbers(n: int, spacing_um: float, boundary: str) -> np.ndarray:
    """Continuous wavenumbers (1/um) of the transform modes along one axis."""
    if boundary == "neumann":
        return math.pi * np.arange(n) / (n * spacing_um)
    if boundary == "periodic":
        return 2.0 * math.pi * fft.fftfreq(n, d=spacing_um)
    raise ValueError(f"boundary must be 'neumann' or 'periodic' here, got {boundary!r}")


def _rfft_wavenumbers(n: int, spacing_um: float) -> np.ndarray:
    return 2.0 * math.pi * fft.rfftfreq(n, d=spacing_um)


def lateral_factors(
    n: int, spacing_um: float, boundary: str, diffusion_um2_s, tau_s: float, *, half: bool = False
) -> np.ndarray:
    """1-D decay factors ``exp(-D k^2 tau)``, shape ``(n_modes, *D.shape)``.

    ``half`` selects the ``rfft`` half spectrum of a periodic axis.
    """
    if isinstance(n, bool) or not isinstance(n, int) or n < 1:
        raise ValueError(f"n must be an int >= 1, got {n!r}")
    spacing = float(spacing_um)
    if not (math.isfinite(spacing) and spacing > 0.0):
        raise ValueError(f"spacing_um must be finite and > 0, got {spacing_um}")
    tau = _check_tau(tau_s)
    d = np.asarray(diffusion_um2_s, dtype=np.float64)
    if half:
        if boundary != "periodic":
            raise ValueError("half spectra exist only for periodic axes")
        k = _rfft_wavenumbers(n, spacing)
    else:
        k = _wavenumbers(n, spacing, boundary)
    k2 = (k * k).reshape((k.size,) + (1,) * d.ndim)
    return np.exp(-(d * tau) * k2)


def lateral_diffuse(
    conc: np.ndarray,
    diffusion_um2_s,
    tau_s: float,
    *,
    boundary: Boundary,
    spacing_um: float | None = None,
    wavenumber_per_um: np.ndarray | None = None,
    toth: TothLateral | None = None,
) -> np.ndarray:
    """Advance ``conc`` by exact lateral diffusion over ``tau_s`` (spec D1 3.5).

    Grid boundaries (``"neumann"``, ``"periodic"``): ``conc`` ``(ny, nx, *trailing)`` with the
    lateral axes first, typically ``(ny, nx, L, K)``; ``diffusion_um2_s`` broadcasts to
    ``trailing`` (``(L, K)``: per cell level and component); ``spacing_um`` is the cell size.
    ``"pair"``: ``conc`` ``(2 n_pairs, *trailing)`` ordered ``a0, b0, a1, b1, ...`` and
    ``wavenumber_per_um`` ``(n_pairs,)``. ``toth`` applies the TOTH treatment (TOTH on the last
    axis). Returns a new array; the input is not modified.
    """
    if boundary not in BOUNDARIES:
        raise ValueError(f"boundary must be one of {BOUNDARIES}, got {boundary!r}")
    c = np.asarray(conc, dtype=np.float64)
    tau = _check_tau(tau_s)
    lateral_ndim = 1 if boundary == "pair" else 2
    if c.ndim < lateral_ndim:
        raise ValueError(f"conc needs at least {lateral_ndim} lateral axes, got shape {c.shape}")
    if not np.all(np.isfinite(c)):
        raise ValueError("conc must be finite")
    trailing = c.shape[lateral_ndim:]
    d = np.array(_check_diffusion(diffusion_um2_s, trailing), copy=True)
    work = c
    mobile = None
    if toth is not None:
        work, d, mobile = _toth_setup(c, d, toth)
    if boundary == "pair":
        out = _diffuse_pair(work, d, tau, wavenumber_per_um)
    else:
        if wavenumber_per_um is not None:
            raise ValueError("wavenumber_per_um is only used with boundary='pair'")
        if spacing_um is None:
            raise ValueError("spacing_um is required for grid boundaries")
        out = _diffuse_grid(work, d, tau, float(spacing_um), boundary)
    if toth is not None:
        i = toth.index
        out[..., i] = c[..., i] + (out[..., i] - mobile)
    return out


def _toth_setup(c: np.ndarray, d: np.ndarray, toth: TothLateral):
    if not isinstance(toth, TothLateral):
        raise TypeError("toth must be a TothLateral")
    if c.ndim < 2:
        raise ValueError("the TOTH treatment needs a component axis last")
    i = toth.index
    if isinstance(i, bool) or not isinstance(i, int) or not (0 <= i < c.shape[-1]):
        raise ValueError(f"toth.index must be an int in [0, {c.shape[-1]}), got {i!r}")
    bound = np.asarray(toth.gel_bound, dtype=np.float64)
    if bound.shape != c.shape[:-1] or not np.all(np.isfinite(bound)):
        raise ValueError(f"toth.gel_bound must be finite with shape {c.shape[:-1]}")
    d_toth = np.asarray(toth.diffusion_um2_s, dtype=np.float64)
    level_shape = d.shape[:-1]
    try:
        d_toth = np.broadcast_to(d_toth, level_shape)
    except ValueError:
        raise ValueError(
            f"toth.diffusion_um2_s shape {d_toth.shape} does not broadcast to {level_shape}"
        ) from None
    if not (np.all(np.isfinite(d_toth)) and np.all(d_toth >= 0.0)):
        raise ValueError("toth.diffusion_um2_s must be finite and >= 0")
    mobile = c[..., i] - bound
    work = np.array(c, copy=True)
    work[..., i] = mobile
    d[..., i] = d_toth
    return work, d, mobile


def _diffuse_grid(c: np.ndarray, d: np.ndarray, tau: float, spacing: float, boundary: str):
    if not (math.isfinite(spacing) and spacing > 0.0):
        raise ValueError(f"spacing_um must be finite and > 0, got {spacing}")
    ny, nx = c.shape[:2]
    if boundary == "neumann":
        fy = lateral_factors(ny, spacing, "neumann", d, tau)
        fx = lateral_factors(nx, spacing, "neumann", d, tau)
        spec = fft.dctn(c, type=2, norm="ortho", axes=(0, 1), workers=1)
        spec *= fy[:, None]
        spec *= fx[None, :]
        return fft.idctn(spec, type=2, norm="ortho", axes=(0, 1), workers=1)
    fy = lateral_factors(ny, spacing, "periodic", d, tau)
    fx = lateral_factors(nx, spacing, "periodic", d, tau, half=True)
    spec = fft.rfftn(c, axes=(0, 1), workers=1)
    spec *= fy[:, None]
    spec *= fx[None, :]
    return fft.irfftn(spec, s=(ny, nx), axes=(0, 1), workers=1)


def _diffuse_pair(c: np.ndarray, d: np.ndarray, tau: float, wavenumber) -> np.ndarray:
    if wavenumber is None:
        raise ValueError("boundary='pair' needs wavenumber_per_um")
    k = np.asarray(wavenumber, dtype=np.float64)
    if k.ndim != 1 or not np.all(np.isfinite(k)):
        raise ValueError("wavenumber_per_um must be a finite 1-D array")
    if c.shape[0] != 2 * k.size:
        raise ValueError(
            f"pair columns must be 2 * n_pairs = {2 * k.size} (a0, b0, a1, b1, ...), "
            f"got {c.shape[0]}"
        )
    a = c[0::2]
    b = c[1::2]
    k2 = (k * k).reshape((k.size,) + (1,) * d.ndim)
    factor = np.exp(-(d * tau) * k2)
    mean = 0.5 * (a + b)
    half_diff = (0.5 * (a - b)) * factor
    out = np.empty_like(c)
    out[0::2] = mean + half_diff
    out[1::2] = mean - half_diff
    return out


def carrier_capacities(
    spec: Speciation, totals: np.ndarray, ph: np.ndarray, mobile: Sequence[bool]
) -> np.ndarray:
    """Buffer capacities of the mobile proton carriers (spec D1 3.3, 3.5), mol/(L pH).

    ``totals`` ``(F, *cells)`` and ``ph`` ``cells`` as for
    :func:`darkroom_atlas.chemistry.speciation.solve_ph`; ``mobile`` ``(F,)`` marks the mobile
    families (gelatin is not). Returns ``(n_mobile + 2, *cells)``: ``beta_f = ln10 T_f
    var_f`` for each mobile family in order (``var_f`` the variance of its proton number, i.e.
    ``-d nbar_f / dpH / ln10``), then ``beta_H = ln10 [H+]`` and ``beta_OH = ln10 [OH-]``.
    """
    t = np.asarray(totals, dtype=np.float64)
    p = np.asarray(ph, dtype=np.float64)
    mob = np.asarray(mobile, dtype=bool)
    if mob.shape != (spec.n_families,):
        raise ValueError(f"mobile must have shape ({spec.n_families},), got {mob.shape}")
    if t.shape != (spec.n_families, *p.shape):
        raise ValueError(f"totals must have shape {(spec.n_families, *p.shape)}, got {t.shape}")
    alpha = fractions(spec, p)  # (F, J + 1, *cells)
    rows = []
    for f, fam in enumerate(spec.families):
        if not mob[f]:
            continue
        n_states = fam.n_states
        nn = np.asarray(fam.protons_rel, dtype=np.float64).reshape((n_states,) + (1,) * p.ndim)
        a = alpha[f, :n_states]
        nbar = (nn * a).sum(axis=0)
        var = (a * (nn - nbar) ** 2).sum(axis=0)
        rows.append(_LN10 * t[f] * var)
    h, oh = water_ions(spec, p)
    rows.append(_LN10 * h)
    rows.append(_LN10 * oh)
    return np.stack(rows)


def toth_lateral_diffusion(
    beta: np.ndarray, diffusion_um2_s: np.ndarray, *, inner: tuple[slice, ...] | None = None
) -> np.ndarray:
    """``D_toth,l = sum_c D_c,l bbar_c,l / sum_c bbar_c,l`` per level (spec D1 3.5).

    ``beta`` ``(C, *cells, L)``: buffer capacity of each mobile carrier (families, then H+ and
    OH-, as :func:`carrier_capacities` gives, moved so that the level axis is last);
    ``diffusion_um2_s`` ``(C, L)``: the carrier's lateral D on each level. ``bbar`` is the mean
    over the cells selected by ``inner`` (slices of the cell axes; default all cells of the
    block). Returns ``(L,)``.
    """
    b = np.asarray(beta, dtype=np.float64)
    d = np.asarray(diffusion_um2_s, dtype=np.float64)
    if b.ndim < 2 or d.shape != (b.shape[0], b.shape[-1]):
        raise ValueError(
            f"beta (C, *cells, L) and diffusion (C, L) disagree: {b.shape} vs {d.shape}"
        )
    if not (np.all(np.isfinite(b)) and np.all(b >= 0.0)):
        raise ValueError("beta must be finite and >= 0")
    if not (np.all(np.isfinite(d)) and np.all(d >= 0.0)):
        raise ValueError("diffusion must be finite and >= 0 um^2/s")
    n_cell_axes = b.ndim - 2
    if inner is not None:
        if len(inner) != n_cell_axes or not all(isinstance(s, slice) for s in inner):
            raise ValueError(f"inner must be {n_cell_axes} slices of the cell axes")
        b = b[(slice(None), *inner, slice(None))]
    if b.size == 0:
        raise ValueError("no cells selected")
    mean = b.reshape(b.shape[0], -1, b.shape[-1]).mean(axis=1)  # (C, L)
    total = mean.sum(axis=0)
    if not np.all(total > 0.0):
        raise ValueError("the mobile buffer capacity of a level is zero")
    return (d * mean).sum(axis=0) / total
