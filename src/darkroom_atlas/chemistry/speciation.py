"""Speciation of acid-base families and the pH of a cell (spec D1 2.3.4).

Given the family totals ``T_f`` and the proton excess ``TOTH`` of each cell, pH solves the
proton condition::

    F(pH) = sum_f T_f nbar_f(pH) + [H+] - [OH-] - TOTH = 0

``F`` decreases monotonically in pH and ``-dF/dpH = beta``, the buffer capacity::

    beta = ln 10 (sum_f T_f var_f(pH) + [H+] + [OH-])

where ``var_f`` is the variance of the proton number over the states of family ``f``. Root
finding is a safeguarded Newton iteration on ``pH in [0, 14.5]``, vectorised over cells: it
falls back to bisection of the current bracket when the Newton step leaves the bracket or does
not shrink the residual fast enough. Each cell is frozen once its update is within ``tol``, so
its result does not depend on the other cells in the batch.
"""

import math
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

from darkroom_atlas.chemistry.activity import conditional_pka, water_constants
from darkroom_atlas.chemistry.species import AcidBaseFamily
from darkroom_atlas.data.models import ActivityParams

__all__ = [
    "PH_MAX",
    "PH_MIN",
    "Speciation",
    "buffer_capacity",
    "fractions",
    "make_speciation",
    "mean_protons",
    "solve_ph",
    "water_ions",
]

PH_MIN = 0.0
"""Lower end of the pH bracket (spec D1 2.3.4)."""

PH_MAX = 14.5
"""Upper end of the pH bracket (spec D1 2.3.4)."""

_LN10 = math.log(10.0)


def _readonly(a: np.ndarray) -> np.ndarray:
    a = np.array(a, dtype=np.float64)
    a.flags.writeable = False
    return a


@dataclass(frozen=True, eq=False)
class Speciation:
    """Conditional constants of a set of families at one temperature and ionic strength.

    ``pka`` is ``(F, J_max)``: row ``f`` holds the mixed constants ``pK'`` of family ``f``,
    padded with NaN. ``log10_kw`` is the thermodynamic ``log10 K_w``; ``log10_gamma_h`` and
    ``log10_gamma_oh`` convert activities of the water ions to concentrations.
    """

    families: tuple[AcidBaseFamily, ...]
    pka: np.ndarray  # (F, J_max), conditional, NaN padded
    log10_kw: float
    log10_gamma_h: float
    log10_gamma_oh: float
    # Per family: cumulative pK' sums (n_f + 1,) and proton numbers (n_f + 1,).
    _cum: tuple[np.ndarray, ...] = field(init=False, repr=False)
    _protons: tuple[np.ndarray, ...] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        families = tuple(self.families)
        n_f = len(families)
        j_max = max((f.n_states - 1 for f in families), default=0)
        pka = np.array(self.pka, dtype=np.float64).reshape(n_f, -1) if n_f else np.zeros((0, 0))
        if pka.shape != (n_f, j_max):
            raise ValueError(f"pka must have shape {(n_f, j_max)}, got {pka.shape}")
        cum, protons = [], []
        for f, fam in enumerate(families):
            n = fam.n_states - 1
            row = pka[f, :n]
            if not np.all(np.isfinite(row)) or not np.all(np.isnan(pka[f, n:])):
                raise ValueError(
                    f"pka row {f} ({fam.component}) needs {n} finite values then NaN padding"
                )
            cum.append(_readonly(np.concatenate([[0.0], np.cumsum(row)])))
            protons.append(_readonly(np.asarray(fam.protons_rel, dtype=np.float64)))
        for name in ("log10_kw", "log10_gamma_h", "log10_gamma_oh"):
            if not math.isfinite(getattr(self, name)):
                raise ValueError(f"{name} must be finite")
        object.__setattr__(self, "families", families)
        object.__setattr__(self, "pka", _readonly(pka))
        object.__setattr__(self, "_cum", tuple(cum))
        object.__setattr__(self, "_protons", tuple(protons))

    @property
    def n_families(self) -> int:
        """Number of families ``F``."""
        return len(self.families)

    @property
    def j_max(self) -> int:
        """Largest number of acid steps over the families."""
        return self.pka.shape[1]


def make_speciation(
    families: Sequence[AcidBaseFamily],
    temperature_k: float,
    ionic_strength: float,
    cation_molality: float,
    params: ActivityParams,
) -> Speciation:
    """Conditional constants of ``families`` at ``(T, I, m_cat)`` (spec D1 2.3.3, 2.3.4)."""
    families = tuple(families)
    j_max = max((f.n_states - 1 for f in families), default=0)
    pka = np.full((len(families), j_max), np.nan)
    for i, fam in enumerate(families):
        pka[i, : fam.n_states - 1] = conditional_pka(
            fam, temperature_k, ionic_strength, cation_molality, params
        )
    log10_kw, log10_gh, log10_goh = water_constants(temperature_k, ionic_strength, params)
    return Speciation(
        families=families,
        pka=pka,
        log10_kw=log10_kw,
        log10_gamma_h=log10_gh,
        log10_gamma_oh=log10_goh,
    )


def _family_fractions(cum: np.ndarray, ph: np.ndarray) -> np.ndarray:
    """State fractions of one family, ``(n + 1, N)``, from the log-weights in base 10."""
    j = np.arange(cum.size, dtype=np.float64)[:, None]
    log_w = (j * ph[None, :] - cum[:, None]) * _LN10
    log_w -= log_w.max(axis=0, keepdims=True)
    w = np.exp(log_w)
    return w / w.sum(axis=0, keepdims=True)


def _moments(spec: Speciation, ph: np.ndarray, f: int) -> tuple[np.ndarray, np.ndarray]:
    """Mean and variance of the proton number of family ``f`` at ``ph`` (shape ``(N,)``)."""
    alpha = _family_fractions(spec._cum[f], ph)
    n = spec._protons[f][:, None]
    nbar = (n * alpha).sum(axis=0)
    var = (alpha * (n - nbar[None, :]) ** 2).sum(axis=0)
    return nbar, var


def _water(spec: Speciation, ph: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    h = 10.0 ** (-ph - spec.log10_gamma_h)
    oh = 10.0 ** (spec.log10_kw + ph - spec.log10_gamma_oh)
    return h, oh


def _residual(
    spec: Speciation, totals: np.ndarray, toth: np.ndarray, ph: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """``F(pH)`` and ``beta`` for flat cells: totals ``(F, N)``, toth and ph ``(N,)``."""
    h, oh = _water(spec, ph)
    f_val = h - oh
    var_sum = h + oh
    for f in range(spec.n_families):
        nbar, var = _moments(spec, ph, f)
        f_val = f_val + totals[f] * nbar
        var_sum = var_sum + totals[f] * var
    return f_val - toth, _LN10 * var_sum


def _flat_totals(spec: Speciation, totals) -> tuple[np.ndarray, tuple[int, ...]]:
    t = np.asarray(totals, dtype=np.float64)
    if t.ndim == 0 or t.shape[0] != spec.n_families:
        raise ValueError(
            f"totals must have shape (F, ...) with F = {spec.n_families} families, got {t.shape}"
        )
    if not np.all(np.isfinite(t)) or np.any(t < 0.0):
        raise ValueError("totals must be finite and >= 0 mol/L")
    shape = t.shape[1:]
    return t.reshape(spec.n_families, math.prod(shape)), shape


def _flat_ph(ph, shape: tuple[int, ...] | None = None) -> tuple[np.ndarray, tuple[int, ...]]:
    p = np.asarray(ph, dtype=np.float64)
    if shape is not None:
        p = np.broadcast_to(p, shape)
    if not np.all(np.isfinite(p)):
        raise ValueError("pH must be finite")
    return p.reshape(-1), p.shape


def solve_ph(
    spec: Speciation,
    totals: np.ndarray,
    toth: np.ndarray,
    *,
    ph0: np.ndarray | None = None,
    tol: float = 1e-12,
    max_iter: int = 30,
) -> np.ndarray:
    """pH of each cell (spec D1 2.3.4): totals ``(F, N)`` mol/L, toth ``(N,)`` mol/L.

    Any trailing cell shape works: totals ``(F, *S)``, toth broadcastable to ``S``, result
    ``S``. ``ph0`` is the starting value (the previous pH; default 7), clipped to the
    bracket. Raises ``ValueError`` if a root lies outside ``[0, 14.5]`` and ``RuntimeError``
    if a cell does not converge within ``max_iter`` iterations.
    """
    t, shape = _flat_totals(spec, totals)
    q = np.asarray(toth, dtype=np.float64)
    try:
        q = np.broadcast_to(q, shape)
    except ValueError:
        raise ValueError(f"toth shape {q.shape} does not match the cell shape {shape}") from None
    if not np.all(np.isfinite(q)):
        raise ValueError("toth must be finite")
    q = q.reshape(-1)
    if not (tol > 0.0 and max_iter >= 1):
        raise ValueError("tol must be > 0 and max_iter >= 1")
    n = q.size
    if ph0 is None:
        x = np.full(n, 7.0)
    else:
        x, _ = _flat_ph(ph0, shape)
        x = np.clip(x, PH_MIN, PH_MAX).astype(np.float64, copy=True)
    lo = np.full(n, PH_MIN)
    hi = np.full(n, PH_MAX)
    # Last step and the one before (Numerical Recipes rtsafe); a bisection "step" is the
    # half-width of the bracket, so it measures the remaining uncertainty.
    dx = np.full(n, PH_MAX - PH_MIN)
    dx_old = np.full(n, PH_MAX - PH_MIN)
    active = np.arange(n)
    for _ in range(max_iter):
        if active.size == 0:
            break
        xa = x[active]
        f_val, beta = _residual(spec, t[:, active], q[active], xa)
        lo_a = np.where(f_val > 0.0, xa, lo[active])
        hi_a = np.where(f_val < 0.0, xa, hi[active])
        newton_step = f_val / beta
        newton = xa + newton_step
        bisect = (
            ~((newton >= lo_a) & (newton <= hi_a))  # outside the bracket, or not finite
            | (np.abs(2.0 * f_val) > np.abs(dx_old[active] * beta))  # not shrinking fast enough
        )
        step = np.where(bisect, 0.5 * (hi_a - lo_a), np.abs(newton_step))
        x_new = np.where(bisect, 0.5 * (lo_a + hi_a), newton)
        exact = f_val == 0.0
        x_new = np.where(exact, xa, x_new)
        done = exact | (step <= tol)
        x[active] = x_new
        lo[active] = lo_a
        hi[active] = hi_a
        dx_old[active] = dx[active]
        dx[active] = step
        active = active[~done]
    if active.size:
        _check_bracket(spec, t[:, active], q[active], active, shape)
        raise RuntimeError(
            f"solve_ph: {active.size} cell(s) did not converge to {tol} in {max_iter} "
            f"iterations (first flat index {int(active[0])})"
        )
    edge = np.nonzero((x <= PH_MIN + 1e-9) | (x >= PH_MAX - 1e-9))[0]
    if edge.size:
        _check_bracket(spec, t[:, edge], q[edge], edge, shape)
    return x.reshape(shape)


def _check_bracket(
    spec: Speciation, t: np.ndarray, q: np.ndarray, index: np.ndarray, shape: tuple[int, ...]
) -> None:
    f_lo, _ = _residual(spec, t, q, np.full(q.size, PH_MIN))
    f_hi, _ = _residual(spec, t, q, np.full(q.size, PH_MAX))
    bad = (f_lo < 0.0) | (f_hi > 0.0)
    if np.any(bad):
        k = int(np.argmax(bad))
        where = np.unravel_index(int(index[k]), shape) if shape else ()
        side = "below 0" if f_lo[k] < 0.0 else "above 14.5"
        raise ValueError(
            f"solve_ph: the pH of cell {where} is {side} (TOTH {q[k]:.6g} mol/L lies outside "
            "the range the families can buffer within the pH bracket [0, 14.5])"
        )


def fractions(spec: Speciation, ph: np.ndarray) -> np.ndarray:
    """State fractions ``alpha_f,j``, shape ``(F, J_max + 1, *ph.shape)``.

    States a family does not have are 0, so each family's fractions sum to 1.
    """
    p, shape = _flat_ph(ph)
    out = np.zeros((spec.n_families, spec.j_max + 1, p.size))
    for f in range(spec.n_families):
        alpha = _family_fractions(spec._cum[f], p)
        out[f, : alpha.shape[0]] = alpha
    return out.reshape((spec.n_families, spec.j_max + 1, *shape))


def mean_protons(spec: Speciation, ph: np.ndarray) -> np.ndarray:
    """Mean proton number ``nbar_f`` relative to the reference state, ``(F, *ph.shape)``."""
    p, shape = _flat_ph(ph)
    out = np.zeros((spec.n_families, p.size))
    for f in range(spec.n_families):
        out[f] = _moments(spec, p, f)[0]
    return out.reshape((spec.n_families, *shape))


def water_ions(spec: Speciation, ph: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Concentrations ``[H+]`` and ``[OH-]`` in mol/L at activity pH ``ph``."""
    p, shape = _flat_ph(ph)
    h, oh = _water(spec, p)
    return h.reshape(shape), oh.reshape(shape)


def buffer_capacity(spec: Speciation, totals: np.ndarray, ph: np.ndarray) -> np.ndarray:
    """``beta = -dF/dpH`` in mol/(L pH), shape of the cells (spec D1 2.3.4)."""
    t, shape = _flat_totals(spec, totals)
    p, _ = _flat_ph(ph, shape)
    _, beta = _residual(spec, t, np.zeros(p.size), p)
    return beta.reshape(shape)
