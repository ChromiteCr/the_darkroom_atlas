"""Two-stream (Kubelka-Munk) light transport through the depth of the emulsion (spec §2.1.3).

For each band the downward and upward diffuse fluxes ``I+``, ``I-`` obey

    dI+/ds = -(K + S) I+ + S I-,   dI-/ds = (K + S) I- - S I+,   s in [0, L_e]

with ``I+(0) = T_top N + r_i I-(0)`` and ``I-(L_e) = R_bot I+(L_e) + J``. Two exact closed
forms are used:

* ``beta L_e > 1``: the exponential form of §2.1.3 with ``beta = sqrt(K (K + 2S))`` and
  ``R_inf = S / (K + S + beta)``; every exponential is <= 1.
* ``beta L_e <= 1``: the hyperbolic form in ``u = I+ + I-``, ``v = I+ - I-``
  (``du/ds = -(K + 2S) v``, ``dv/ds = -K u``) with ``cosh(beta s)`` and ``sinh(beta s)/beta``.
  It stays well conditioned as ``K/S -> 0`` and at ``K = 0``, where the exponential form's
  two modes merge.

Cell integrals over the ``n_z`` equal depth cells are analytic.
"""

from dataclasses import dataclass
from numbers import Integral

import numpy as np

from darkroom_atlas._frozen import readonly_array
from darkroom_atlas.optics.stack import OpticalStack

__all__ = ["DepthResponse", "depth_response"]

# Branch switch on beta * L_e (spec §2.1.3 gives only the exponential form; see module doc).
_HYPERBOLIC_MAX_BETA_L = 1.0


@dataclass(frozen=True, eq=False)
class DepthResponse:
    """Depth response of the layer per band (spec §2.1.3).

    * ``plus`` ``(n_z, B)``: cell mean of the scalar flux ``psi = 2 (I+ + I-)`` for unit
      incident flux on the top (``N = 1, J = 0``), ``g+``.
    * ``minus`` ``(n_z, B)``: the same for unit diffuse flux entering from the base
      (``N = 0, J = 1``), ``g-``.
    * ``into_base`` ``(B,)``: ``tau = (1 - R_bot) I+(L_e)`` for ``N = 1, J = 0``.
    * ``escaped_top`` ``(B,)``: ``(1 - r_i) I-(0)`` for ``N = 1, J = 0``.
    * ``absorbed`` ``(B,)``: ``integral K (I+ + I-) ds`` for ``N = 1, J = 0``.

    All per unit incident flux ``N``; ``escaped_top + absorbed + into_base == T_top``.
    """

    plus: np.ndarray
    minus: np.ndarray
    into_base: np.ndarray
    escaped_top: np.ndarray
    absorbed: np.ndarray

    def __post_init__(self) -> None:
        for name in ("plus", "minus", "into_base", "escaped_top", "absorbed"):
            object.__setattr__(self, name, readonly_array(getattr(self, name), dtype=np.float64))
        if self.plus.ndim != 2 or self.minus.shape != self.plus.shape:
            raise ValueError(
                f"plus and minus must share one (n_z, B) shape, got {self.plus.shape} "
                f"and {self.minus.shape}"
            )
        n_bands = self.plus.shape[1]
        for name in ("into_base", "escaped_top", "absorbed"):
            if getattr(self, name).shape != (n_bands,):
                raise ValueError(
                    f"{name} must have shape ({n_bands},), got {getattr(self, name).shape}"
                )

    @property
    def n_z(self) -> int:
        """Number of depth cells."""
        return self.plus.shape[0]


def depth_response(stack: OpticalStack, n_z: int) -> DepthResponse:
    """Solve the two-stream model of ``stack`` on ``n_z`` equal depth cells (spec §2.1.3)."""
    if not isinstance(stack, OpticalStack):
        raise TypeError(f"stack must be an OpticalStack, got {type(stack).__name__}")
    if isinstance(n_z, bool) or not isinstance(n_z, Integral) or n_z < 1:
        raise ValueError(f"n_z must be a positive integer, got {n_z!r}")
    common = {
        "k": stack.km_k_per_um,
        "s": stack.km_s_per_um,
        "thickness_um": stack.emulsion_thickness_um,
        "t_top": stack.t_top,
        "r_top": stack.r_internal_top,
        "r_bottom": stack.r_bottom,
        "n_z": int(n_z),
    }
    down = _two_stream(n_in=1.0, j_in=0.0, **common)
    up = _two_stream(n_in=0.0, j_in=1.0, **common)
    return DepthResponse(
        plus=2.0 * down["cell_mean_u"],
        minus=2.0 * up["cell_mean_u"],
        into_base=(1.0 - stack.r_bottom) * down["down_bottom"],
        escaped_top=(1.0 - stack.r_internal_top) * down["up_top"],
        absorbed=down["absorbed"],
    )


def _sinhc(x: np.ndarray) -> np.ndarray:
    """``sinh(x)/x`` with the limit 1 at ``x = 0``."""
    x = np.asarray(x, dtype=np.float64)
    safe = np.where(x == 0.0, 1.0, x)
    return np.where(x == 0.0, 1.0, np.sinh(safe) / safe)


def _two_stream(
    *,
    k: np.ndarray,
    s: np.ndarray,
    thickness_um: float,
    t_top: float,
    r_top: float,
    r_bottom: float,
    n_in: float,
    j_in: float,
    n_z: int,
    branch: str = "auto",
) -> dict[str, np.ndarray]:
    """Two-stream solution per band for incident flux ``n_in`` and base flux ``j_in``.

    ``k``, ``s`` are ``(B,)`` with ``k, s >= 0``. Returns ``cell_mean_u`` ``(n_z, B)`` (cell
    mean of ``I+ + I-``), the boundary fluxes ``down_top`` = I+(0), ``up_top`` = I-(0),
    ``down_bottom`` = I+(L), ``up_bottom`` = I-(L) and ``absorbed`` = integral K (I+ + I-) ds,
    each ``(B,)``. ``branch`` forces ``"exponential"`` or ``"hyperbolic"`` (tests only).
    """
    k = np.asarray(k, dtype=np.float64)
    s = np.asarray(s, dtype=np.float64)
    length = float(thickness_um)
    beta = np.sqrt(k * (k + 2.0 * s))
    if branch == "auto":
        hyper = beta * length <= _HYPERBOLIC_MAX_BETA_L
    elif branch == "hyperbolic":
        hyper = np.ones(k.shape, dtype=bool)
    elif branch == "exponential":
        if np.any(beta == 0.0):
            raise ValueError("the exponential branch needs beta > 0")
        hyper = np.zeros(k.shape, dtype=bool)
    else:
        raise ValueError(f"unknown branch {branch!r}")
    edges = np.linspace(0.0, length, n_z + 1)
    out = {
        "cell_mean_u": np.empty((n_z, k.size)),
        **{key: np.empty(k.size) for key in ("down_top", "up_top", "down_bottom", "up_bottom")},
        "absorbed": np.empty(k.size),
    }
    boundary = (t_top * n_in, j_in, r_top, r_bottom)
    for mask, solver in ((hyper, _hyperbolic), (~hyper, _exponential)):
        if np.any(mask):
            part = solver(k[mask], s[mask], beta[mask], length, edges, *boundary)
            for key, value in part.items():
                out[key][..., mask] = value
    return out


def _exponential(
    k: np.ndarray,
    s: np.ndarray,
    beta: np.ndarray,
    length: float,
    edges: np.ndarray,
    source_top: float,
    j_in: float,
    r_top: float,
    r_bottom: float,
) -> dict[str, np.ndarray]:
    """Exponential closed form of spec §2.1.3 (requires ``beta > 0``)."""
    r_inf = s / (k + s + beta)
    e_l = np.exp(-beta * length)
    a11 = 1.0 - r_top * r_inf
    a12 = e_l * (r_inf - r_top)
    a21 = e_l * (r_inf - r_bottom)
    a22 = 1.0 - r_bottom * r_inf
    det = a11 * a22 - a12 * a21
    p = (source_top * a22 - a12 * j_in) / det
    q = (a11 * j_in - a21 * source_top) / det
    width = edges[1] - edges[0]
    # (1 - e^{-beta width}) / beta, evaluated without cancellation.
    cell_factor = -np.expm1(-beta * width) / beta
    lo = edges[:-1, None]
    hi = edges[1:, None]
    int_decay = np.exp(-beta * lo) * cell_factor  # integral of e^{-beta s} over a cell
    int_grow = np.exp(-beta * (length - hi)) * cell_factor  # integral of e^{-beta (L - s)}
    cell_mean_u = (1.0 + r_inf) * (p * int_decay + q * int_grow) / width
    whole = -np.expm1(-beta * length) / beta
    return {
        "cell_mean_u": cell_mean_u,
        "down_top": p + r_inf * q * e_l,
        "up_top": r_inf * p + q * e_l,
        "down_bottom": p * e_l + r_inf * q,
        "up_bottom": r_inf * p * e_l + q,
        "absorbed": k * (1.0 + r_inf) * (p + q) * whole,
    }


def _hyperbolic(
    k: np.ndarray,
    s: np.ndarray,
    beta: np.ndarray,
    length: float,
    edges: np.ndarray,
    source_top: float,
    j_in: float,
    r_top: float,
    r_bottom: float,
) -> dict[str, np.ndarray]:
    """Closed form in ``u = I+ + I-``, ``v = I+ - I-`` (exact; used for ``beta L <= 1``)."""
    k2s = k + 2.0 * s

    def f1(depth: np.ndarray) -> np.ndarray:  # integral_0^depth cosh(beta x) dx
        return depth * _sinhc(beta * depth)

    def f2(depth: np.ndarray) -> np.ndarray:  # integral_0^depth sinh(beta x)/beta dx
        return 0.5 * depth**2 * _sinhc(0.5 * beta * depth) ** 2

    cosh_l = np.cosh(beta * length)
    sinh_l = length * _sinhc(beta * length)  # sinh(beta L) / beta
    # Top: (u0 + v0)/2 - r_i (u0 - v0)/2 = T_top N.
    a11 = np.full(k.shape, 0.5 * (1.0 - r_top))
    a12 = np.full(k.shape, 0.5 * (1.0 + r_top))
    # Bottom: (uL - vL)/2 - R_bot (uL + vL)/2 = J with
    # uL = cosh u0 - (K + 2S) sinh/beta v0,  vL = cosh v0 - K sinh/beta u0.
    a21 = 0.5 * ((1.0 - r_bottom) * cosh_l + (1.0 + r_bottom) * k * sinh_l)
    a22 = -0.5 * ((1.0 - r_bottom) * k2s * sinh_l + (1.0 + r_bottom) * cosh_l)
    det = a11 * a22 - a12 * a21
    u0 = (source_top * a22 - a12 * j_in) / det
    v0 = (a11 * j_in - a21 * source_top) / det
    u_l = cosh_l * u0 - k2s * sinh_l * v0
    v_l = cosh_l * v0 - k * sinh_l * u0
    lo = edges[:-1, None]
    hi = edges[1:, None]
    cell_int = u0 * (f1(hi) - f1(lo)) - k2s * v0 * (f2(hi) - f2(lo))
    width = edges[1] - edges[0]
    whole = u0 * f1(np.float64(length)) - k2s * v0 * f2(np.float64(length))
    return {
        "cell_mean_u": cell_int / width,
        "down_top": 0.5 * (u0 + v0),
        "up_top": 0.5 * (u0 - v0),
        "down_bottom": 0.5 * (u_l + v_l),
        "up_bottom": 0.5 * (u_l - v_l),
        "absorbed": k * whole,
    }
