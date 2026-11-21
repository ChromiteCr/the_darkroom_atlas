"""Two-site blocking chain for latent-image formation (spec D1 2.2.1-2.2.4).

Each grain tracks a main site 1 and one lumped competing site 2. The state ``(k1, k2)``,
``k_s in {0, ..., n_cap}`` (``n_cap`` merges ``>= n_cap``), is stored at index
``k1 * (n_cap + 1) + k2``. Generators use the column convention: ``Q[j, i]`` is the rate from
state ``i`` to state ``j`` and every column sums to zero, so ``P(t) = expm(Q t) e_(0,0)``.

Rates (spec D1 2.2.2), with electron rate ``nu`` at the latent sites:

- ``q1 = 1 / (1 + nu tau_i)`` (M/M/1/1 blocking of the main site),
- ``q2 = (1 - q1) / (S_eff - 1)`` (the tracked competing site),
- the rest ``(1 - q1)(S_eff - 2)/(S_eff - 1)`` is lost: no state change,
- ``k_s -> k_s + 1`` at ``nu q_s`` for ``k_s < n_cap``,
- ``k_s = 1 -> 0`` at ``1/tau0 + rho1 nu`` (thermal decay and hole attack).

The single-atom lifetime ``tau0`` is log-normal across grains (spec D1 2.2.3) and is
integrated by probabilists' Gauss-Hermite quadrature.
"""

import math
from dataclasses import dataclass
from functools import cache

import numpy as np
from numpy.polynomial.hermite_e import hermegauss

from darkroom_atlas.data.models import EmulsionLayer

__all__ = [
    "ChainParams",
    "cluster_distribution",
    "generator",
]


def _finite_nonnegative(name: str, value: float) -> float:
    value = float(value)
    if not (math.isfinite(value) and value >= 0.0):
        raise ValueError(f"{name} must be finite and >= 0, got {value!r}")
    return value


@dataclass(frozen=True)
class ChainParams:
    """Parameters of the latent-image chain (spec D1 2.2.7), in SI-like plain floats.

    ``tau0_median_s`` may be ``math.inf`` (no thermal decay). ``cluster_cap`` must be >= 2:
    the output merges ``n = 1`` into ``n = 0`` (spec D1 2.2.4), so a cap of 1 has no
    developable state.
    """

    efficiency: float
    efficiency_size_exponent: float
    tau0_median_s: float  # math.inf allowed (no decay)
    tau0_sigma_ln: float
    ionic_time_s: float
    competing_sites: float  # >= 2
    hole_attack: float = 0.0
    cluster_cap: int = 5

    def __post_init__(self) -> None:
        eta = float(self.efficiency)
        if not (math.isfinite(eta) and eta > 0.0):
            raise ValueError(f"efficiency must be finite and > 0, got {self.efficiency!r}")
        xi = float(self.efficiency_size_exponent)
        if not math.isfinite(xi):
            raise ValueError(
                f"efficiency_size_exponent must be finite, got {self.efficiency_size_exponent!r}"
            )
        tau0 = float(self.tau0_median_s)
        if not (tau0 > 0.0):  # also rejects NaN; +inf is allowed
            raise ValueError(f"tau0_median_s must be > 0 (inf allowed), got {tau0!r}")
        sigma = _finite_nonnegative("tau0_sigma_ln", self.tau0_sigma_ln)
        tau_i = _finite_nonnegative("ionic_time_s", self.ionic_time_s)
        s_eff = float(self.competing_sites)
        if not (math.isfinite(s_eff) and s_eff >= 2.0):
            raise ValueError(f"competing_sites must be finite and >= 2, got {s_eff!r}")
        rho1 = _finite_nonnegative("hole_attack", self.hole_attack)
        cap = self.cluster_cap
        if isinstance(cap, bool) or not isinstance(cap, (int, np.integer)) or cap < 2:
            raise ValueError(f"cluster_cap must be an integer >= 2, got {cap!r}")
        for name, value in (
            ("efficiency", eta),
            ("efficiency_size_exponent", xi),
            ("tau0_median_s", tau0),
            ("tau0_sigma_ln", sigma),
            ("ionic_time_s", tau_i),
            ("competing_sites", s_eff),
            ("hole_attack", rho1),
            ("cluster_cap", int(cap)),
        ):
            object.__setattr__(self, name, value)

    @property
    def states(self) -> int:
        """Number of chain states ``(n_cap + 1)**2``."""
        return (self.cluster_cap + 1) ** 2

    @classmethod
    def from_layer(cls, layer: EmulsionLayer) -> "ChainParams":
        """Chain parameters from ``layer.latent`` (spec D1 2.2.7)."""
        latent = layer.latent
        if latent is None:
            raise ValueError(f"emulsion layer {layer.id!r}: latent parameters are required")

        def scalar(quantity, unit: str, name: str) -> float:
            value = quantity.to(unit)
            if isinstance(value, np.ndarray):
                raise ValueError(f"latent.{name} must be a single value, got a list")  # noqa: TRY004
            return value

        return cls(
            efficiency=scalar(latent.efficiency, "1", "efficiency"),
            efficiency_size_exponent=scalar(
                latent.efficiency_size_exponent, "1", "efficiency_size_exponent"
            ),
            tau0_median_s=scalar(latent.tau0_median, "s", "tau0_median"),
            tau0_sigma_ln=scalar(latent.tau0_sigma_ln, "1", "tau0_sigma_ln"),
            ionic_time_s=scalar(latent.ionic_time, "s", "ionic_time"),
            competing_sites=scalar(latent.competing_sites, "1", "competing_sites"),
            hole_attack=scalar(latent.hole_attack, "1", "hole_attack"),
            cluster_cap=latent.cluster_cap,
        )


@cache
def _structure(cap: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Unit-rate generators for growth at site 1, growth at site 2 and decay of k_s = 1."""
    n = cap + 1
    s = n * n
    grow1 = np.zeros((s, s))
    grow2 = np.zeros((s, s))
    decay = np.zeros((s, s))
    for k1 in range(n):
        for k2 in range(n):
            i = k1 * n + k2
            if k1 < cap:
                grow1[(k1 + 1) * n + k2, i] += 1.0
                grow1[i, i] -= 1.0
            if k2 < cap:
                grow2[k1 * n + k2 + 1, i] += 1.0
                grow2[i, i] -= 1.0
            if k1 == 1:
                decay[k2, i] += 1.0
                decay[i, i] -= 1.0
            if k2 == 1:
                decay[k1 * n, i] += 1.0
                decay[i, i] -= 1.0
    for array in (grow1, grow2, decay):
        array.flags.writeable = False
    return grow1, grow2, decay


def _generators(
    nu_per_s: np.ndarray, tau0_s: np.ndarray, params: ChainParams, scale: float = 1.0
) -> np.ndarray:
    """Stack of generators ``scale * Q`` for broadcast arrays of ``nu`` and ``tau0``."""
    nu, tau0 = np.broadcast_arrays(
        np.asarray(nu_per_s, dtype=np.float64), np.asarray(tau0_s, dtype=np.float64)
    )
    grow1, grow2, decay = _structure(params.cluster_cap)
    q1 = 1.0 / (1.0 + nu * params.ionic_time_s)
    q2 = (1.0 - q1) / (params.competing_sites - 1.0)
    with np.errstate(divide="ignore"):
        d = 1.0 / tau0 + params.hole_attack * nu  # 1/inf = 0: no thermal decay
    r1 = scale * nu * q1
    r2 = scale * nu * q2
    rd = scale * d
    return r1[..., None, None] * grow1 + r2[..., None, None] * grow2 + rd[..., None, None] * decay


def generator(nu_per_s: float, tau0_s: float, params: ChainParams) -> np.ndarray:
    """Generator ``Q`` (S, S), ``S = (n_cap + 1)**2``, of the chain; columns sum to 0.

    ``nu_per_s`` is the electron rate at the latent sites (1/s), ``tau0_s`` the single-atom
    lifetime of this grain (s, ``math.inf`` for no decay).
    """
    nu = float(nu_per_s)
    tau0 = float(tau0_s)
    if not (math.isfinite(nu) and nu >= 0.0):
        raise ValueError(f"nu_per_s must be finite and >= 0, got {nu_per_s!r}")
    if not (tau0 > 0.0):
        raise ValueError(f"tau0_s must be > 0 (inf allowed), got {tau0_s!r}")
    return _generators(np.float64(nu), np.float64(tau0), params)


@cache
def _quadrature(nodes: int) -> tuple[np.ndarray, np.ndarray]:
    """Probabilists' Gauss-Hermite nodes and weights normalised to sum 1 (spec D1 2.2.3)."""
    x, w = hermegauss(nodes)
    w = w / w.sum()
    x.flags.writeable = False
    w.flags.writeable = False
    return x, w


def _check_quadrature(quadrature: int) -> int:
    if (
        isinstance(quadrature, bool)
        or not isinstance(quadrature, (int, np.integer))
        or quadrature < 1
    ):
        raise ValueError(f"quadrature must be a positive integer, got {quadrature!r}")
    return int(quadrature)


def _tau0_nodes(params: ChainParams, quadrature: int) -> tuple[np.ndarray, np.ndarray]:
    """Lifetime nodes and weights; a single node when the spread is absent or moot."""
    if params.tau0_sigma_ln == 0.0 or math.isinf(params.tau0_median_s):
        return np.array([params.tau0_median_s]), np.array([1.0])
    x, w = _quadrature(quadrature)
    return params.tau0_median_s * np.exp(params.tau0_sigma_ln * x), w


def efficiency(params: ChainParams, ecd_um: float | np.ndarray) -> float | np.ndarray:
    """Electron efficiency ``eta_c = eta_0 (d_c / 1 um)**xi_eta`` (spec D1 2.2.1)."""
    return params.efficiency * np.asarray(ecd_um, dtype=np.float64) ** (
        params.efficiency_size_exponent
    )


_TAYLOR_DEGREE = 18
"""Taylor degree for exp of a shifted generator with 1-norm <= 1: remainder <= 1/19! < 1e-17."""


def _site_occupancy(birth: np.ndarray, decay: np.ndarray, cap: int) -> np.ndarray:
    """``expm(A) e_0`` for one site, (..., cap+1), for dimensionless rates ``birth``, ``decay``.

    ``A`` is the single-site generator times ``t_exp``: ``k -> k+1`` at ``birth`` for
    ``k < cap``, ``1 -> 0`` at ``decay``. The exponential is computed without cancellation:
    scale by ``2**-s`` so the largest exit rate ``c`` is <= 1, write ``A/2**s = N - c I`` with
    ``N >= 0``, sum the Taylor series of ``exp(N)`` (all terms >= 0), multiply by ``e**-c``
    and square ``s`` times. Unlike ``scipy.linalg.expm`` this stays accurate for the very
    stiff decay rates of the outer quadrature nodes (``decay`` up to ~1e19).
    """
    n = cap + 1
    birth = np.asarray(birth, dtype=np.float64)
    decay = np.asarray(decay, dtype=np.float64)
    shape = np.broadcast_shapes(birth.shape, decay.shape)
    b = np.broadcast_to(birth, shape).ravel()
    d = np.broadcast_to(decay, shape).ravel()
    exit_max = b + d  # state 1 has the largest exit rate
    if not np.all(np.isfinite(exit_max)):
        # An infinite rate would make the squaring count saturate and the loop never end.
        raise FloatingPointError(
            "latent chain rate overflow: birth or decay rate times t_exp is not finite "
            "(check tau0_median_s, exposure time and photon counts)"
        )
    squarings = np.zeros(b.shape, dtype=np.intp)
    big = exit_max > 1.0
    squarings[big] = np.ceil(np.log2(exit_max[big])).astype(np.intp)
    factor = np.ldexp(1.0, -squarings)
    bs = b * factor
    ds = d * factor
    cs = exit_max * factor  # <= 1
    nonneg = np.zeros((b.shape[0], n, n))
    k = np.arange(n)
    diag = np.empty((b.shape[0], n))
    diag[:] = (cs[:, None] - bs[:, None]) * (k < cap)
    diag[:, cap] = cs
    diag[:, 1] -= ds
    nonneg[:, k, k] = diag
    nonneg[:, k[1:], k[:-1]] = bs[:, None]
    nonneg[:, 0, 1] = ds
    eye = np.eye(n)
    taylor = np.broadcast_to(eye, nonneg.shape).copy()
    for j in range(_TAYLOR_DEGREE, 0, -1):
        taylor = eye + np.matmul(nonneg, taylor) / j
    expo = taylor * np.exp(-cs)[:, None, None]
    expo /= expo.sum(axis=1, keepdims=True)
    for j in range(int(squarings.max(initial=0))):
        idx = np.nonzero(squarings > j)[0]
        sub = expo[idx]
        sub = np.matmul(sub, sub)
        # Squaring doubles a column-sum error each time (2**s ~ c); restoring the exact
        # column sums of a stochastic matrix removes that growing mode.
        expo[idx] = sub / sub.sum(axis=1, keepdims=True)
    return expo[:, :, 0].reshape(*shape, n)


def _distributions(
    electrons: np.ndarray, exposure_time_s: float, params: ChainParams, quadrature: int
) -> np.ndarray:
    """Cluster-size distributions (..., cap+1) for mean electron counts ``eta_c m``.

    The two sites change state independently given ``nu`` and ``tau0`` (``Q`` is the
    Kronecker sum of two single-site generators), so ``expm(Q t) e_(0,0)`` is the outer
    product of the two single-site solutions. This is exact and equals the 36-state
    ``expm`` of spec D1 2.2.2.
    """
    cap = params.cluster_cap
    tau0, weights = _tau0_nodes(params, quadrature)
    t = exposure_time_s
    electrons = np.asarray(electrons, dtype=np.float64)[..., None]  # (..., 1) = nu t
    # Overflow to inf is caught as a FloatingPointError in _site_occupancy.
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        nu = electrons / t
        q1 = 1.0 / (1.0 + nu * params.ionic_time_s)
        q2 = (1.0 - q1) / (params.competing_sites - 1.0)
        decay = t / tau0 + params.hole_attack * electrons  # (..., Q); t/inf = 0: no decay
        birth1 = electrons * q1
        birth2 = electrons * q2
    site1 = _site_occupancy(birth1, decay, cap)  # (..., Q, cap+1)
    site2 = _site_occupancy(birth2, decay, cap)
    # n = max(k1, k2): P(n = j) = p1(j) F2(j) + F1(j-1) p2(j), all terms >= 0.
    cdf1 = np.cumsum(site1, axis=-1)
    cdf2 = np.cumsum(site2, axis=-1)
    below1 = np.concatenate([np.zeros_like(cdf1[..., :1]), cdf1[..., :-1]], axis=-1)
    p_nodes = site1 * cdf2 + below1 * site2
    p = np.einsum("...qn,q->...n", p_nodes, weights)
    # Single atoms decay before development (spec D1 2.2.4): n = 1 joins n = 0.
    p[..., 0] += p[..., 1]
    p[..., 1] = 0.0
    if not np.all(np.isfinite(p)):
        raise FloatingPointError("latent chain produced non-finite probabilities")
    p = np.clip(p, 0.0, None)
    return p / p.sum(axis=-1, keepdims=True)


def _check_time(exposure_time_s: float) -> float:
    t = float(exposure_time_s)
    if not (math.isfinite(t) and t > 0.0):
        raise ValueError(f"exposure_time_s must be finite and > 0, got {exposure_time_s!r}")
    return t


def cluster_distribution(
    mean_photons: float,
    exposure_time_s: float,
    params: ChainParams,
    *,
    ecd_um: float = 1.0,
    quadrature: int = 9,
) -> np.ndarray:
    """Latent cluster-size distribution ``p(n | m)``, shape ``(cap+1,)`` (spec D1 2.2.4).

    ``mean_photons`` is the mean number of photons absorbed per grain during a
    constant-intensity exposure of ``exposure_time_s``; ``ecd_um`` sets the size-dependent
    efficiency. ``n = max(k1, k2)`` at the end of the exposure; entry 1 is 0 (single atoms
    are folded into ``n = 0``) and the entries sum to 1.
    """
    m = float(mean_photons)
    if not (math.isfinite(m) and m >= 0.0):
        raise ValueError(f"mean_photons must be finite and >= 0, got {mean_photons!r}")
    t = _check_time(exposure_time_s)
    d = float(ecd_um)
    if not (math.isfinite(d) and d > 0.0):
        raise ValueError(f"ecd_um must be finite and > 0, got {ecd_um!r}")
    q = _check_quadrature(quadrature)
    electrons = np.array([float(efficiency(params, d)) * m])  # python float: inf, no warning
    return _distributions(electrons, t, params, q)[0]
