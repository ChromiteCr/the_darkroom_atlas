"""Finite well-mixed tank: bulk balance, agitation mixing and the prepass norms (spec D1 2.4.9, 3.7).

Bulk balance (spec D1 2.4.9), for each transported component and TOTH::

    dC_k,bulk/dt = -1e-15 (A_film / V_tank) Jbar_k(t)        [mol L^-1 s^-1]

``Jbar_k`` is the area-mean flux into the film (mol L^-1 um s^-1) at the current Dirichlet face,
``A_film`` in um^2, ``V_tank`` in L; ``1e-15`` converts um^3 to L, so the coefficient
``1e-15 A / V`` has unit 1/um. It is advanced explicitly once per step with the step's mean flux
(spec D1 3.2): ``c_{n+1} = c_n - coef h_n Jbar_n``.

Agitation mixing (spec D1 2.4.8, 2.4.9). At each agitation start the upper liquid section joins
the tank: its inventory excess ``dI_k = sum_upper Delta_l (C_k,l - C_k,bulk)`` (mol L^-1 um) is a
jump ``dC_bulk = +coef dIbar``. In ``flux_mean`` / ``Jbar_blk`` it is booked as the impulse
``-dIbar / h`` of that step, so re-integrating the tank from the recorded fluxes includes it.

Prepass check (spec D1 3.7, revised): ``delta = max_b |M_b(c') - M_b(c)| / (M_b(c) + 0.3 g/m2)``
on the developed silver of the 65 prepass columns, pass at ``delta <= 1e-3``; diagnostics are
the contraction estimate ``L_hat = delta_{i+1} / delta_i`` (warn above 0.5) and the
depletion-normed trajectory change ``max_{k,n} |c'_k - c_k| / s_k`` with ``s_k = max(|c_k(0)|,
max_n |c_k(t_n) - c_k(0)|, 1e-6 M)``.

Arrays use components on the last axis (``K_c`` components then TOTH, as ``flux_mean`` of spec
D1 10). Public functions::

    tank_coefficient(film_area_um2: float, volume_l: float) -> float          # 1/um
    advance_bulk(bulk (K,), flux (K,), h_s, coefficient) -> (K,)
    integrate_bulk(initial (K,), flux_mean (S, K), times_s (S+1,), coefficient) -> (S+1, K)
    mixing_inventory(conc (..., n_up, K), thickness_um (n_up,), bulk (K,)) -> (..., K)
    mixing_flux(inventory_mean (K,), h_s) -> (K,)
    area_mean(values (M, K), weights (M,)) -> (K,)
    combine_flux(block_flux, zero_flux, image_area_fraction) -> same shape
    bulk_at(times_s (S+1,), values (S+1, K), t_s) -> (K,)
    depletion_scale(trajectory (S+1, K)) -> (K,)
    depletion_change(new (S+1, K), current (S+1, K)) -> float
    prepass_delta(silver_new (B,), silver_current (B,)) -> float
    prepass_passed(delta) -> bool
    contraction(deltas) -> tuple[float, ...]
"""

import itertools
import math
from collections.abc import Sequence

import numpy as np

__all__ = [
    "CONTRACTION_WARN",
    "DEPLETION_FLOOR_MOL_L",
    "MAX_EXTRA_PASSES",
    "SILVER_FLOOR_G_M2",
    "TANK_TOLERANCE",
    "UM3_TO_L",
    "advance_bulk",
    "area_mean",
    "bulk_at",
    "combine_flux",
    "contraction",
    "depletion_change",
    "depletion_scale",
    "integrate_bulk",
    "mixing_flux",
    "mixing_inventory",
    "prepass_delta",
    "prepass_passed",
    "tank_coefficient",
]

UM3_TO_L = 1e-15
"""Litres per cubic micrometre, the factor of the tank balance (spec D1 2.4.9)."""

TANK_TOLERANCE = 1e-3
"""epsilon_tank: the prepass check passes at delta <= 1e-3 (spec D1 3.7)."""

SILVER_FLOOR_G_M2 = 0.3
"""Floor added to M_b in the denominator of delta, g/m2 (spec D1 3.7)."""

DEPLETION_FLOOR_MOL_L = 1e-6
"""Floor of the depletion scale s_k, mol/L (spec D1 3.7)."""

CONTRACTION_WARN = 0.5
"""Warn when the estimated contraction L_hat exceeds this (spec D1 3.7)."""

MAX_EXTRA_PASSES = 2
"""Recomputations of all blocks after a failed check, at most (spec D1 3.7)."""


def _finite(name: str, value, ndim: int | None = None) -> np.ndarray:
    a = np.asarray(value, dtype=np.float64)
    if ndim is not None and a.ndim != ndim:
        raise ValueError(f"{name} must be {ndim}-D, got shape {a.shape}")
    if not np.all(np.isfinite(a)):
        raise ValueError(f"{name} must be finite")
    return a


def tank_coefficient(film_area_um2: float, volume_l: float) -> float:
    """``1e-15 A_film / V_tank`` in 1/um (spec D1 2.4.9); ``ValueError`` unless both > 0."""
    area, volume = float(film_area_um2), float(volume_l)
    if not (math.isfinite(area) and area > 0.0):
        raise ValueError(f"film_area_um2 must be finite and > 0, got {film_area_um2}")
    if not (math.isfinite(volume) and volume > 0.0):
        raise ValueError(f"volume_l must be finite and > 0, got {volume_l}")
    return UM3_TO_L * area / volume


def _coefficient(coefficient: float) -> float:
    c = float(coefficient)
    if not (math.isfinite(c) and c >= 0.0):
        raise ValueError(f"coefficient must be finite and >= 0, got {coefficient}")
    return c


def advance_bulk(bulk: np.ndarray, flux: np.ndarray, h_s: float, coefficient: float) -> np.ndarray:
    """One explicit step ``c - coefficient h Jbar`` (spec D1 2.4.9, 3.2)."""
    c = _finite("bulk", bulk, 1)
    j = _finite("flux", flux, 1)
    if j.shape != c.shape:
        raise ValueError(f"flux shape {j.shape} differs from bulk shape {c.shape}")
    h = float(h_s)
    if not (math.isfinite(h) and h > 0.0):
        raise ValueError(f"h_s must be finite and > 0, got {h_s}")
    return c - (_coefficient(coefficient) * h) * j


def integrate_bulk(
    initial: np.ndarray, flux_mean: np.ndarray, times_s: np.ndarray, coefficient: float
) -> np.ndarray:
    """Bulk trajectory ``(S+1, K)`` from the per-step mean fluxes ``(S, K)`` (spec D1 3.7).

    Uses :func:`advance_bulk` step by step, so a coupled run that advanced its tank with the
    same fluxes is reproduced bit for bit.
    """
    c = _finite("initial", initial, 1)
    j = _finite("flux_mean", flux_mean, 2)
    t = _finite("times_s", times_s, 1)
    if j.shape != (t.size - 1, c.size):
        raise ValueError(f"flux_mean must have shape {(t.size - 1, c.size)}, got {j.shape}")
    h = np.diff(t)
    if not np.all(h > 0.0):
        raise ValueError("times_s must increase strictly")
    coef = _coefficient(coefficient)
    out = np.empty((t.size, c.size))
    out[0] = c
    for n in range(h.size):
        out[n + 1] = out[n] - (coef * h[n]) * j[n]
    return out


def mixing_inventory(conc: np.ndarray, thickness_um: np.ndarray, bulk: np.ndarray) -> np.ndarray:
    """Inventory excess of the upper liquid section over the bath (spec D1 2.4.8).

    ``conc`` ``(..., n_up, K)`` for the cells above ``delta_agit`` just before the agitation
    start, ``thickness_um`` ``(n_up,)``, ``bulk`` ``(K,)``. Returns ``sum_l Delta_l (C_l -
    C_bulk)`` with shape ``(..., K)`` in mol L^-1 um, summed bottom to top.
    """
    c = _finite("conc", conc)
    d = _finite("thickness_um", thickness_um, 1)
    b = _finite("bulk", bulk, 1)
    if c.ndim < 2 or c.shape[-2:] != (d.size, b.size):
        raise ValueError(f"conc must have shape (..., {d.size}, {b.size}), got {c.shape}")
    if not np.all(d > 0.0):
        raise ValueError("thickness_um must be > 0")
    out = np.zeros(c.shape[:-2] + (b.size,))
    for lev in range(d.size):
        out = out + d[lev] * (c[..., lev, :] - b)
    return out


def mixing_flux(inventory_mean: np.ndarray, h_s: float) -> np.ndarray:
    """The mixing impulse booked as a flux of the step that starts the agitation: ``-dIbar/h``
    (spec D1 2.4.9)."""
    i = _finite("inventory_mean", inventory_mean)
    h = float(h_s)
    if not (math.isfinite(h) and h > 0.0):
        raise ValueError(f"h_s must be finite and > 0, got {h_s}")
    return -i / h


def area_mean(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Weighted mean over the first axis: ``sum_m w_m v_m / sum_m w_m`` (columns' cell weights)."""
    v = _finite("values", values)
    w = _finite("weights", weights, 1)
    if v.ndim < 1 or v.shape[0] != w.size:
        raise ValueError(f"values have {v.shape[:1]} rows, weights {w.size}")
    if np.any(w < 0.0) or not np.sum(w) > 0.0:
        raise ValueError("weights must be >= 0 with a positive sum")
    ww = w.reshape((w.size,) + (1,) * (v.ndim - 1))
    return np.sum(ww * v, axis=0) / np.sum(w)


def combine_flux(
    block_flux: np.ndarray, zero_flux: np.ndarray, image_area_fraction: float
) -> np.ndarray:
    """``f_img Jbar_blk + (1 - f_img) Jbar_zero`` (spec D1 3.7 step 4)."""
    jb = _finite("block_flux", block_flux)
    jz = _finite("zero_flux", zero_flux)
    if jb.shape != jz.shape:
        raise ValueError(f"block_flux {jb.shape} and zero_flux {jz.shape} differ in shape")
    f = float(image_area_fraction)
    if not (0.0 <= f <= 1.0):
        raise ValueError(f"image_area_fraction must be in [0, 1], got {image_area_fraction}")
    return f * jb + (1.0 - f) * jz


def bulk_at(times_s: np.ndarray, values: np.ndarray, t_s: float) -> np.ndarray:
    """Bulk composition at ``t_s``, linear within a step and exact at the nodes (spec D1 3.2)."""
    t = _finite("times_s", times_s, 1)
    v = _finite("values", values)
    if v.ndim < 1 or v.shape[0] != t.size:
        raise ValueError(f"values must have {t.size} rows, got shape {v.shape}")
    if t.size < 2 or not np.all(np.diff(t) > 0.0):
        raise ValueError("times_s needs at least 2 strictly increasing times")
    x = float(t_s)
    if not (t[0] <= x <= t[-1]):
        raise ValueError(f"t_s = {t_s} lies outside [{t[0]}, {t[-1]}]")
    n = int(np.searchsorted(t, x, side="right")) - 1
    if n >= t.size - 1:
        return v[-1].copy()
    w = (x - t[n]) / (t[n + 1] - t[n])
    if w == 0.0:
        return v[n].copy()
    return v[n] + w * (v[n + 1] - v[n])


def depletion_scale(trajectory: np.ndarray) -> np.ndarray:
    """``s_k = max(|c_k(0)|, max_n |c_k(t_n) - c_k(0)|, 1e-6 M)`` (spec D1 3.7)."""
    c = _finite("trajectory", trajectory, 2)
    if c.shape[0] < 1:
        raise ValueError("trajectory is empty")
    change = np.max(np.abs(c - c[0]), axis=0)
    return np.maximum(np.maximum(np.abs(c[0]), change), DEPLETION_FLOOR_MOL_L)


def depletion_change(new: np.ndarray, current: np.ndarray) -> float:
    """``max_{k,n} |c'_k(t_n) - c_k(t_n)| / s_k`` with ``s_k`` from the current trajectory
    (spec D1 3.7; a diagnostic)."""
    a = _finite("new", new, 2)
    c = _finite("current", current, 2)
    if a.shape != c.shape:
        raise ValueError(f"trajectories differ in shape: {a.shape} vs {c.shape}")
    return float(np.max(np.abs(a - c) / depletion_scale(c)))


def prepass_delta(silver_new: np.ndarray, silver_current: np.ndarray) -> float:
    """``delta = max_b |M_b(c') - M_b(c)| / (M_b(c) + 0.3 g/m2)`` (spec D1 3.7)."""
    a = _finite("silver_new", silver_new, 1)
    c = _finite("silver_current", silver_current, 1)
    if a.shape != c.shape or a.size == 0:
        raise ValueError(f"silver arrays must share a non-empty shape: {a.shape} vs {c.shape}")
    if np.any(c < 0.0) or np.any(a < 0.0):
        raise ValueError("developed silver must be >= 0 g/m2")
    return float(np.max(np.abs(a - c) / (c + SILVER_FLOOR_G_M2)))


def prepass_passed(delta: float) -> bool:
    """The check of spec D1 3.7: ``delta <= epsilon_tank = 1e-3``."""
    d = float(delta)
    if math.isnan(d) or d < 0.0:
        raise ValueError(f"delta must be >= 0, got {delta}")
    return d <= TANK_TOLERANCE


def contraction(deltas: Sequence[float]) -> tuple[float, ...]:
    """Contraction estimates ``L_hat_i = delta_{i+1} / delta_i`` of successive passes (spec D1
    3.7); ``inf`` where ``delta_i = 0`` and ``delta_{i+1} > 0``, ``0`` where both are 0."""
    d = [float(x) for x in deltas]
    if any(math.isnan(x) or x < 0.0 for x in d):
        raise ValueError("deltas must be >= 0")
    out = []
    for a, b in itertools.pairwise(d):
        if a > 0.0:
            out.append(b / a)
        else:
            out.append(0.0 if b == 0.0 else math.inf)
    return tuple(out)
