"""ROS2 column step and block-Thomas solver for development columns (spec D1 3.3, 3.6, 12.8).

One step of the two-stage W-method ROS2 (Verwer, Spee, Blom & Hundsdorfer 1999) for a batch of
``N`` independent columns, gamma = 1 + 1/sqrt(2)::

    (I - gamma h J) k1 = f(t, y)
    (I - gamma h J) k2 = f(t + h, y + h k1) - 2 k1
    y_new = y + (3/2) h k1 + (1/2) h k2
    err   = y_new - (y + h k1) = (h/2)(k1 + k2)

ROS2 keeps order 2 for any approximation ``J`` (W-method), and for a non-autonomous ``f`` the
dropped ``df/dt`` column is such an approximation, so ``f`` is evaluated at ``t`` and ``t + h``
with no extra term.

Column layout (:class:`ColumnLayout`). Levels are ordered bottom (film base) to top (tank), as
in :class:`darkroom_atlas.develop.grid.ColumnGrid`. The lowest ``n_dense`` levels (the emulsion
cells) carry ``n_comp`` transported components followed by ``n_local`` local unknowns
(``Lambda_z``, ``D_z``, ...), block size ``b = n_comp + n_local``. The ``n_diag`` levels above
(overcoat, liquid) carry the ``n_comp`` components only. A flat state ``y`` of shape
``(N, size)`` holds the dense levels first (level by level, components then local unknowns),
then the diagonal levels.

Jacobian structure (:class:`ColumnJacobian`, spec D1 3.3): a dense ``b x b`` block per dense
level (reaction coupling, plus transport self-terms), a diagonal per diagonal level (each
component's own transport), and diagonal couplings between neighbouring levels acting on the
components only (``lower[l] = d f(level l+1) / d C(level l)``, ``upper[l] = d f(level l) /
d C(level l+1)``, per component). Local unknowns have no neighbour coupling.

Solver. ``I - gamma h J`` is eliminated top-down (from the tank towards the base), so the
reduced blocks of the diagonal levels stay diagonal and dense elimination happens only in the
dense levels (spec D1 3.3). Dense blocks are factored by a hand-written vectorised Gaussian
elimination with partial pivoting over the batch (no ``np.linalg`` batch solves, spec D1 3.6,
12.8). Every operation is element-wise across columns, with no BLAS call and no reduction, so
a column's result is bit-identical whatever the batch size or the memory alignment of the
inputs (spec D1 3.6).

Composition with eliminated unknowns. The cluster kernels G, H, Q are eliminated by a Schur
complement (spec D1 3.3); the caller does that around this solver: ``factor_columns(jac,
gamma_h, dense_update=...)`` adds the Schur term to the dense blocks of ``I - gamma h J``, and
the caller's own ``solve`` closure handed to :func:`ros2_step` reduces and back-substitutes the
eliminated part around :func:`solve_columns`.

Frozen cells. During agitation the upper liquid cells are not solved (spec D1 2.4.8): give them
zero rows in ``f`` and ``J`` (so ``M = I`` there, ``k = 0``) and no coupling across the
``delta_agit`` face, then overwrite them with the bulk value after the step.

Negative concentrations (spec D1 3.3, revised): steps are never rejected or split on a fixed
schedule; :func:`min_scaled_concentration` records ``min C_k / s_k`` with ``s_k = max(c_k,bulk(0),
1 mM)`` (:func:`concentration_scale`). The adaptive reference profile treats ``min C/s_k <
-1e-6`` as a failed error test of the whole batch (:func:`adaptive_accept`). Nothing is
clipped.

Public functions::

    factor_columns(jac: ColumnJacobian, gamma_h: float, *, dense_update=None) -> ColumnFactor
    solve_columns(factor: ColumnFactor, rhs: ndarray (N, size)) -> ndarray (N, size)
    ros2_step(f, t: float, y: ndarray (N, n), h: float, solve) -> (y_new, err)
        # f(t, y) -> (N, n); solve(rhs) solves (I - GAMMA h J) k = rhs for this h
    ros2_step_columns(f, t, y, h, jac: ColumnJacobian) -> (y_new, err)
    toth_transport_jacobian(beta_carriers, beta_total, face_diffusion_um2_s,
                            face_distance_um, thickness_um) -> (lower, upper, self)
    concentration_scale(bulk_initial: ndarray (K,)) -> ndarray (K,)
    min_scaled_concentration(conc: ndarray (..., K), scale: ndarray (K,)) -> float
    error_norm(err, y_old, y_new, *, rtol: float, scale) -> float
    adaptive_accept(err_norm: float, min_scaled: float) -> bool
    next_step_size(h: float, err_norm: float, accepted: bool) -> float
"""

import math
from collections.abc import Callable
from dataclasses import dataclass
from numbers import Integral

import numpy as np

__all__ = [
    "GAMMA",
    "MIN_SCALE_MOL_L",
    "NEGATIVE_REJECT",
    "NEGATIVE_WARN",
    "ColumnFactor",
    "ColumnJacobian",
    "ColumnLayout",
    "adaptive_accept",
    "concentration_scale",
    "error_norm",
    "factor_columns",
    "min_scaled_concentration",
    "next_step_size",
    "ros2_step",
    "ros2_step_columns",
    "solve_columns",
    "toth_transport_jacobian",
]

GAMMA = 1.0 + 1.0 / math.sqrt(2.0)
"""ROS2 parameter gamma = 1 + 1/sqrt(2), L-stable (spec D1 3.3)."""

MIN_SCALE_MOL_L = 1e-3
"""Floor of the concentration scale s_k = max(c_k,bulk(0), 1 mM) (spec D1 3.3)."""

NEGATIVE_REJECT = -1e-6
"""Adaptive ROS2 rejects the batch's step when min C/s_k falls below this (spec D1 3.3)."""

NEGATIVE_WARN = -1e-3
"""min C/s_k over a run below this means the schedule is too coarse; warn (spec D1 3.3)."""

# Step-size controller of the adaptive reference profile. The spec (D1 3.3) fixes only the
# acceptance rule; these are the usual safety factor and growth bounds (Hairer & Wanner IV.8).
_SAFETY = 0.9
_FAC_MIN = 0.2
_FAC_MAX = 5.0
_FAC_REJECT = 0.5  # largest factor after a rejection (also for a negativity-only rejection)


@dataclass(frozen=True)
class ColumnLayout:
    """Unknowns of one column: ``n_dense`` levels of ``n_comp + n_local``, then ``n_diag``
    levels of ``n_comp`` (bottom to top)."""

    n_dense: int
    n_diag: int
    n_comp: int
    n_local: int = 0

    def __post_init__(self) -> None:
        for name in ("n_dense", "n_diag", "n_comp", "n_local"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, Integral) or value < 0:
                raise ValueError(f"{name} must be an int >= 0, got {value!r}")
            object.__setattr__(self, name, int(value))  # numpy integers stored as int
        if self.n_comp < 1:
            raise ValueError("n_comp must be >= 1")
        if self.n_dense + self.n_diag < 1:
            raise ValueError("a column needs at least one level")

    @property
    def block(self) -> int:
        """Unknowns per dense level, ``n_comp + n_local``."""
        return self.n_comp + self.n_local

    @property
    def n_levels(self) -> int:
        """Number of levels ``L = n_dense + n_diag``."""
        return self.n_dense + self.n_diag

    @property
    def size(self) -> int:
        """Unknowns per column."""
        return self.n_dense * self.block + self.n_diag * self.n_comp

    def split(self, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Views of a flat ``(N, size)`` state: dense ``(N, n_dense, b)``, diagonal
        ``(N, n_diag, n_comp)``."""
        y = np.asarray(y)
        if y.ndim != 2 or y.shape[1] != self.size:
            raise ValueError(f"state must have shape (N, {self.size}), got {y.shape}")
        cut = self.n_dense * self.block
        n = y.shape[0]
        return (
            y[:, :cut].reshape(n, self.n_dense, self.block),
            y[:, cut:].reshape(n, self.n_diag, self.n_comp),
        )

    def join(self, dense: np.ndarray, diag: np.ndarray) -> np.ndarray:
        """Flat ``(N, size)`` state from dense ``(N, n_dense, b)`` and diagonal
        ``(N, n_diag, n_comp)`` parts."""
        dense = np.asarray(dense, dtype=np.float64)
        diag = np.asarray(diag, dtype=np.float64)
        n = dense.shape[0]
        if dense.shape != (n, self.n_dense, self.block) or diag.shape != (
            n,
            self.n_diag,
            self.n_comp,
        ):
            raise ValueError(
                f"parts must have shapes (N, {self.n_dense}, {self.block}) and "
                f"(N, {self.n_diag}, {self.n_comp}), got {dense.shape} and {diag.shape}"
            )
        return np.concatenate([dense.reshape(n, -1), diag.reshape(n, -1)], axis=1)


@dataclass(frozen=True, eq=False)
class ColumnJacobian:
    """Approximate Jacobian of a batch of columns in the structure of spec D1 3.3.

    ``dense`` ``(N, n_dense, b, b)``; ``diag`` ``(N, n_diag, n_comp)``; ``lower`` and ``upper``
    ``(N, L - 1, n_comp)`` (or broadcastable to it, e.g. ``(L - 1, n_comp)`` when all columns
    share the transport): ``lower[:, l, k] = d f_k(level l+1) / d C_k(level l)`` and
    ``upper[:, l, k] = d f_k(level l) / d C_k(level l+1)``.
    """

    layout: ColumnLayout
    dense: np.ndarray
    diag: np.ndarray
    lower: np.ndarray
    upper: np.ndarray

    def __post_init__(self) -> None:
        lay = self.layout
        if not isinstance(lay, ColumnLayout):
            raise TypeError("layout must be a ColumnLayout")
        dense = np.asarray(self.dense, dtype=np.float64)
        diag = np.asarray(self.diag, dtype=np.float64)
        if dense.ndim != 4 or dense.shape[1:] != (lay.n_dense, lay.block, lay.block):
            raise ValueError(
                f"dense must have shape (N, {lay.n_dense}, {lay.block}, {lay.block}), "
                f"got {dense.shape}"
            )
        n = dense.shape[0]
        if diag.shape != (n, lay.n_diag, lay.n_comp):
            raise ValueError(
                f"diag must have shape ({n}, {lay.n_diag}, {lay.n_comp}), got {diag.shape}"
            )
        target = (n, lay.n_levels - 1, lay.n_comp)
        couplings = []
        for name in ("lower", "upper"):
            a = np.asarray(getattr(self, name), dtype=np.float64)
            try:
                a = np.broadcast_to(a, target)
            except ValueError:
                raise ValueError(f"{name} must broadcast to {target}, got {a.shape}") from None
            couplings.append(a)
        for name, a in (("dense", dense), ("diag", diag), ("lower", couplings[0])):
            if not np.all(np.isfinite(a)):
                raise ValueError(f"{name} must be finite")
        if not np.all(np.isfinite(couplings[1])):
            raise ValueError("upper must be finite")
        object.__setattr__(self, "dense", dense)
        object.__setattr__(self, "diag", diag)
        object.__setattr__(self, "lower", couplings[0])
        object.__setattr__(self, "upper", couplings[1])

    @property
    def n_columns(self) -> int:
        """Batch size ``N``."""
        return self.dense.shape[0]

    def to_dense(self) -> np.ndarray:
        """The full ``(N, size, size)`` matrix (for tests and debugging; not for hot code)."""
        lay = self.layout
        n, b, nc = self.n_columns, lay.block, lay.n_comp
        out = np.zeros((n, lay.size, lay.size))
        offsets = _level_offsets(lay)
        for lev in range(lay.n_levels):
            o = offsets[lev]
            if lev < lay.n_dense:
                out[:, o : o + b, o : o + b] = self.dense[:, lev]
            else:
                idx = np.arange(nc)
                out[:, o + idx, o + idx] = self.diag[:, lev - lay.n_dense]
        idx = np.arange(nc)
        for lev in range(lay.n_levels - 1):
            lo, hi = offsets[lev], offsets[lev + 1]
            out[:, hi + idx, lo + idx] = self.lower[:, lev]
            out[:, lo + idx, hi + idx] = self.upper[:, lev]
        return out


def _level_offsets(lay: ColumnLayout) -> list[int]:
    out, o = [], 0
    for lev in range(lay.n_levels):
        out.append(o)
        o += lay.block if lev < lay.n_dense else lay.n_comp
    return out


@dataclass(frozen=True, eq=False)
class ColumnFactor:
    """Block-Thomas factorisation of ``I - gamma h J`` (from :func:`factor_columns`).

    Arrays keep the batch axis last, so every update runs over contiguous column vectors:
    ``lu``/``piv``: pivoted LU of the reduced dense blocks ``(n_dense, b, b, N)`` /
    ``(n_dense, b, N)``; ``schur``: the first ``n_comp`` columns of the inverse reduced dense
    blocks ``(n_dense, b, n_comp, N)`` (zero for level 0, which needs none); ``diag``: reduced diagonal blocks ``(n_diag, n_comp,
    N)``; ``lower``/``upper``: couplings of ``I - gamma h J``, ``(L - 1, n_comp, N)``.
    """

    layout: ColumnLayout
    lu: np.ndarray
    piv: np.ndarray
    schur: np.ndarray
    diag: np.ndarray
    lower: np.ndarray
    upper: np.ndarray

    @property
    def n_columns(self) -> int:
        """Batch size ``N``."""
        return self.diag.shape[-1] if self.layout.n_dense == 0 else self.lu.shape[-1]


def _batch_last(a: np.ndarray) -> np.ndarray:
    """Contiguous copy with axis 0 (the columns) moved to the end."""
    return np.ascontiguousarray(np.moveaxis(a, 0, -1))


def _lu_factor(a: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Pivoted LU of a batch ``(b, b, N)`` (batch last), on a copy: returns ``(lu, piv)``.

    ``piv[i]`` is the original row now at row ``i``. Ties in the pivot search go to the first
    maximum (``argmax``), so the result is deterministic.
    """
    lu = np.array(a, dtype=np.float64, copy=True)
    b, _, n = lu.shape
    piv = np.broadcast_to(np.arange(b)[:, None], (b, n)).copy()
    for k in range(b):
        p = k + np.argmax(np.abs(lu[k:, k, :]), axis=0)
        cols = np.nonzero(p != k)[0]
        if cols.size:
            pk = p[cols]
            row_k = lu[k][:, cols].copy()  # (b, s)
            lu[k][:, cols] = lu[pk, :, cols].T
            lu[pk, :, cols] = row_k.T
            piv_k = piv[k, cols].copy()
            piv[k, cols] = piv[pk, cols]
            piv[pk, cols] = piv_k
        pivot = lu[k, k]
        if np.any(pivot == 0.0):
            raise ValueError("singular block in the column system (zero pivot)")
        if k + 1 < b:
            mult = lu[k + 1 :, k] / pivot
            lu[k + 1 :, k] = mult
            lu[k + 1 :, k + 1 :] -= mult[:, None, :] * lu[k, k + 1 :][None, :, :]
    return lu, piv


def _lu_solve(lu: np.ndarray, piv: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    """Solve with a batch LU (batch last): ``rhs`` ``(b, N)`` or ``(b, r, N)``; no reductions."""
    vec = rhs.ndim == 2
    x = rhs[:, None, :] if vec else rhs
    x = np.take_along_axis(x, piv[:, None, :], axis=0).astype(np.float64, copy=True)
    b = lu.shape[0]
    tmp = np.empty_like(x[0])
    for i in range(1, b):  # unit lower triangle: x_i -= L_ij x_j
        xi = x[i]
        for j in range(i):
            np.multiply(lu[i, j], x[j], out=tmp)
            np.subtract(xi, tmp, out=xi)
    for i in range(b - 1, -1, -1):  # upper triangle
        xi = x[i]
        for j in range(i + 1, b):
            np.multiply(lu[i, j], x[j], out=tmp)
            np.subtract(xi, tmp, out=xi)
        np.divide(xi, lu[i, i], out=xi)
    return x[:, 0, :] if vec else x


def factor_columns(
    jac: ColumnJacobian, gamma_h: float, *, dense_update: np.ndarray | None = None
) -> ColumnFactor:
    """Factor ``M = I - gamma_h J`` (plus ``dense_update`` on the dense blocks) top-down.

    ``dense_update`` ``(N, n_dense, b, b)`` is added to the dense blocks of ``M`` before the
    elimination (the Schur complement of eliminated unknowns, spec D1 3.3). ``ValueError`` on a
    non-finite or negative ``gamma_h``, a bad ``dense_update`` shape or a singular block.
    """
    if not isinstance(jac, ColumnJacobian):
        raise TypeError("jac must be a ColumnJacobian")
    gamma_h = float(gamma_h)
    if not (math.isfinite(gamma_h) and gamma_h >= 0.0):
        raise ValueError(f"gamma_h must be finite and >= 0, got {gamma_h}")
    lay = jac.layout
    n, b, nc, nd = jac.n_columns, lay.block, lay.n_comp, lay.n_dense
    m_dense = np.eye(b)[None, None] - gamma_h * jac.dense
    if dense_update is not None:
        upd = np.asarray(dense_update, dtype=np.float64)
        if upd.shape != m_dense.shape:
            raise ValueError(f"dense_update must have shape {m_dense.shape}, got {upd.shape}")
        if not np.all(np.isfinite(upd)):
            raise ValueError("dense_update must be finite")
        m_dense = m_dense + upd
    m_dense = _batch_last(m_dense)  # (nd, b, b, N)
    diag = _batch_last(1.0 - gamma_h * jac.diag)  # (n_diag, nc, N)
    lower = _batch_last(-gamma_h * jac.lower)  # (L - 1, nc, N)
    upper = _batch_last(-gamma_h * jac.upper)
    # Diagonal levels, top-down: d'_l = d_l - upper_l lower_l / d'_{l+1}.
    for i in range(lay.n_diag - 2, -1, -1):
        lev = nd + i
        diag[i] = diag[i] - upper[lev] * lower[lev] / diag[i + 1]
    if np.any(diag == 0.0):
        raise ValueError("singular diagonal level in the column system")
    lu = np.empty((nd, b, b, n))
    piv = np.empty((nd, b, n), dtype=np.intp)
    schur = np.zeros((nd, b, nc, n))  # level 0 is never used
    unit = np.broadcast_to(np.eye(b)[:, :nc, None], (b, nc, n))
    idx = np.arange(nc)
    for lev in range(nd - 1, -1, -1):
        block = m_dense[lev]
        if lev == nd - 1 and lay.n_diag:
            # Level above is diagonal: subtract upper * lower / d' on the component diagonal.
            block[idx, idx] -= upper[lev] * lower[lev] / diag[0]
        elif lev < nd - 1:
            # Level above is dense: subtract diag(upper) inv(D'_{l+1})[:nc, :nc] diag(lower).
            block[:nc, :nc] -= upper[lev][:, None, :] * schur[lev + 1, :nc] * lower[lev][None]
        lu[lev], piv[lev] = _lu_factor(block)
        if lev > 0:  # needed by the level below and by the back substitution
            schur[lev] = _lu_solve(lu[lev], piv[lev], unit)
    return ColumnFactor(
        layout=lay, lu=lu, piv=piv, schur=schur, diag=diag, lower=lower, upper=upper
    )


def solve_columns(factor: ColumnFactor, rhs: np.ndarray) -> np.ndarray:
    """Solve ``M x = rhs`` for a flat ``(N, size)`` right-hand side with a :class:`ColumnFactor`."""
    lay = factor.layout
    r_dense, r_diag = lay.split(np.asarray(rhs, dtype=np.float64))
    n = factor.n_columns
    if r_dense.shape[0] != n:
        raise ValueError(f"rhs has {r_dense.shape[0]} columns, the factor {n}")
    nd, nc = lay.n_dense, lay.n_comp
    rs = _batch_last(r_dense)  # (nd, b, N)
    rd = _batch_last(r_diag)  # (n_diag, nc, N)
    lower, upper, diag = factor.lower, factor.upper, factor.diag
    # Forward (top-down) reduction of the right-hand side.
    for i in range(lay.n_diag - 2, -1, -1):
        lev = nd + i
        rd[i] = rd[i] - upper[lev] * rd[i + 1] / diag[i + 1]
    z = np.empty_like(rs)
    for lev in range(nd - 1, -1, -1):
        if lev == nd - 1 and lay.n_diag:
            rs[lev, :nc] = rs[lev, :nc] - upper[lev] * rd[0] / diag[0]
        elif lev < nd - 1:
            rs[lev, :nc] = rs[lev, :nc] - upper[lev] * z[lev + 1, :nc]
        z[lev] = _lu_solve(factor.lu[lev], factor.piv[lev], rs[lev])
    # Back substitution, bottom-up.
    x_dense = np.empty_like(rs)
    x_diag = np.empty_like(rd)
    for lev in range(nd):
        if lev == 0:
            x_dense[0] = z[0]
            continue
        v = lower[lev - 1] * x_dense[lev - 1, :nc]  # (nc, N)
        acc = z[lev]
        sch = factor.schur[lev]  # (b, nc, N)
        for j in range(nc):
            acc = acc - sch[:, j] * v[j]
        x_dense[lev] = acc
    for i in range(lay.n_diag):
        lev = nd + i
        if lev == 0:
            x_diag[0] = rd[0] / diag[0]
            continue
        below = x_dense[nd - 1, :nc] if i == 0 else x_diag[i - 1]
        x_diag[i] = (rd[i] - lower[lev - 1] * below) / diag[i]
    return lay.join(np.moveaxis(x_dense, -1, 0), np.moveaxis(x_diag, -1, 0))


def ros2_step(
    f: Callable[[float, np.ndarray], np.ndarray],
    t: float,
    y: np.ndarray,
    h: float,
    solve: Callable[[np.ndarray], np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    """One ROS2 step from ``(t, y)`` with step ``h`` (spec D1 3.3).

    ``solve(rhs)`` must solve ``(I - GAMMA h J) k = rhs`` for this ``h``. Returns ``(y_new,
    err)`` with ``err = (h/2)(k1 + k2)``, the embedded first-order error estimate. The step is
    never rejected or split here (spec D1 3.3).
    """
    h = float(h)
    if not (math.isfinite(h) and h > 0.0):
        raise ValueError(f"h must be finite and > 0, got {h}")
    y = np.asarray(y, dtype=np.float64)
    f0 = np.asarray(f(t, y), dtype=np.float64)
    if f0.shape != y.shape:
        raise ValueError(f"f returned shape {f0.shape}, expected {y.shape}")
    k1 = solve(f0)
    f1 = np.asarray(f(t + h, y + h * k1), dtype=np.float64)
    k2 = solve(f1 - 2.0 * k1)
    y_new = y + (1.5 * h) * k1 + (0.5 * h) * k2
    err = (0.5 * h) * (k1 + k2)
    return y_new, err


def ros2_step_columns(
    f: Callable[[float, np.ndarray], np.ndarray],
    t: float,
    y: np.ndarray,
    h: float,
    jac: ColumnJacobian,
) -> tuple[np.ndarray, np.ndarray]:
    """:func:`ros2_step` with the block-Thomas solve of ``jac`` (factored once for both stages)."""
    factor = factor_columns(jac, GAMMA * float(h))
    return ros2_step(f, t, y, h, lambda rhs: solve_columns(factor, rhs))


def toth_transport_jacobian(
    beta_carriers: np.ndarray,
    beta_total: np.ndarray,
    face_diffusion_um2_s: np.ndarray,
    face_distance_um: np.ndarray,
    thickness_um: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """TOTH transport entries of ``J`` (spec D1 3.3), only the TOTH-on-TOTH derivatives.

    The face flux of spec D1 2.4.8 is built from the mobile carriers ``c`` (the mobile families'
    proton amounts ``P_f``, ``[H+]``, ``[OH-]``). Its derivative with respect to the TOTH of the
    cell on either side of face ``l+1/2`` is ``+-G`` with::

        G(face, cell) = sum_c Dbar_c,face beta_c,cell / (beta_cell h_face)

    (``beta_c`` the carrier's buffer capacity, ``beta_cell`` the cell's total including
    gelatin). Inputs: ``beta_carriers`` ``(C, N, L)``, ``beta_total`` ``(N, L)``,
    ``face_diffusion_um2_s`` ``(C, N, L-1)`` (or broadcastable), ``face_distance_um``
    ``(L-1,)`` = ``(Delta_l + Delta_l+1)/2``, ``thickness_um`` ``(L,)``. Returns ``(lower,
    upper, self)`` in the convention of :class:`ColumnJacobian`: ``lower`` ``(N, L-1)`` =
    ``G(l+1/2, l)/Delta_l+1``, ``upper`` ``(N, L-1)`` = ``G(l+1/2, l+1)/Delta_l``, ``self``
    ``(N, L)`` = ``-(G(l-1/2, l) + G(l+1/2, l))/Delta_l``. The flux through the top Dirichlet
    face is the caller's (it adds its own self term to the cell below that face).
    """
    bc = np.asarray(beta_carriers, dtype=np.float64)
    bt = np.asarray(beta_total, dtype=np.float64)
    thick = np.asarray(thickness_um, dtype=np.float64)
    hf = np.asarray(face_distance_um, dtype=np.float64)
    if bc.ndim != 3 or bt.shape != bc.shape[1:]:
        raise ValueError(
            f"beta_carriers (C, N, L) and beta_total (N, L) disagree in shape: {bc.shape} vs "
            f"{bt.shape}"
        )
    n_lev = bc.shape[2]
    if thick.shape != (n_lev,) or hf.shape != (n_lev - 1,):
        raise ValueError(
            f"thickness_um must be ({n_lev},) and face_distance_um ({n_lev - 1},), got "
            f"{thick.shape} and {hf.shape}"
        )
    target = (bc.shape[0], bc.shape[1], n_lev - 1)
    try:
        df = np.broadcast_to(np.asarray(face_diffusion_um2_s, dtype=np.float64), target)
    except ValueError:
        raise ValueError(f"face_diffusion_um2_s must broadcast to {target}") from None
    for name, a in (("beta_carriers", bc), ("face_diffusion_um2_s", df)):
        if not (np.all(np.isfinite(a)) and np.all(a >= 0.0)):
            raise ValueError(f"{name} must be finite and >= 0")
    if not (np.all(np.isfinite(bt)) and np.all(bt > 0.0)):
        raise ValueError("beta_total must be finite and > 0")
    if not (np.all(thick > 0.0) and np.all(hf > 0.0)):
        raise ValueError("thickness_um and face_distance_um must be > 0")
    below = np.zeros(target[1:])  # sum_c Dbar_c beta_c,l over face l+1/2, cell l (below)
    above = np.zeros(target[1:])  # same with the cell l+1 (above)
    for c in range(bc.shape[0]):  # explicit loop: no reduction (spec D1 3.6)
        below = below + df[c] * bc[c, :, :-1]
        above = above + df[c] * bc[c, :, 1:]
    g_below = below / (bt[:, :-1] * hf)
    g_above = above / (bt[:, 1:] * hf)
    lower = g_below / thick[1:]
    upper = g_above / thick[:-1]
    self_ = np.zeros(bt.shape)
    self_[:, :-1] -= g_below / thick[:-1]
    self_[:, 1:] -= g_above / thick[1:]
    return lower, upper, self_


def concentration_scale(bulk_initial: np.ndarray) -> np.ndarray:
    """``s_k = max(c_k,bulk(0), 1 mM)`` per component (spec D1 3.3)."""
    c = np.asarray(bulk_initial, dtype=np.float64)
    if c.ndim != 1 or not np.all(np.isfinite(c)):
        raise ValueError("bulk_initial must be a finite 1-D array (one value per component)")
    return np.maximum(c, MIN_SCALE_MOL_L)


def min_scaled_concentration(conc: np.ndarray, scale: np.ndarray) -> float:
    """``min C_k / s_k`` over all cells (spec D1 3.3); ``conc`` ``(..., K)``, ``scale`` ``(K,)``.

    Pass concentrations only: TOTH (which may be negative), the clocks and the kernels are not
    concentrations.
    """
    c = np.asarray(conc, dtype=np.float64)
    s = np.asarray(scale, dtype=np.float64)
    if s.ndim != 1 or c.ndim < 1 or c.shape[-1] != s.size:
        raise ValueError(f"conc (..., K) and scale (K,) disagree: {c.shape} vs {s.shape}")
    if not np.all(s > 0.0):
        raise ValueError("scale must be > 0")
    if c.size == 0:
        raise ValueError("conc is empty")
    return float(np.min(c / s))


def error_norm(
    err: np.ndarray, y_old: np.ndarray, y_new: np.ndarray, *, rtol: float, scale: np.ndarray
) -> float:
    """Max-norm of ``err / (rtol * max(|y_old|, |y_new|, scale))`` over the whole batch.

    ``scale`` broadcasts against the state (``s_k`` for concentrations, a typical size for
    other unknowns). The batch is accepted or rejected as a whole (spec D1 3.3), hence the
    max over all columns.
    """
    rtol = float(rtol)
    if not (math.isfinite(rtol) and rtol > 0.0):
        raise ValueError(f"rtol must be finite and > 0, got {rtol}")
    s = np.asarray(scale, dtype=np.float64)
    if not np.all(s > 0.0):
        raise ValueError("scale must be > 0")
    w = rtol * np.maximum(np.maximum(np.abs(y_old), np.abs(y_new)), s)
    return float(np.max(np.abs(err) / w))


def adaptive_accept(err_norm: float, min_scaled: float) -> bool:
    """Adaptive acceptance (spec D1 3.3): error test passes and ``min C/s_k >= -1e-6``."""
    if math.isnan(err_norm) or math.isnan(min_scaled):
        return False
    return err_norm <= 1.0 and min_scaled >= NEGATIVE_REJECT


def next_step_size(h: float, err_norm: float, accepted: bool) -> float:
    """Next step of the adaptive profile: ``h * clip(0.9 err^(-1/2), 0.2, 5)``.

    The exponent 1/2 follows from the first-order embedded estimate (``err ~ h^2``). After a
    rejection the factor is at most 0.5, also when only the negativity test failed (spec D1
    3.3: treated as a failed error test).
    """
    h = float(h)
    if not (math.isfinite(h) and h > 0.0):
        raise ValueError(f"h must be finite and > 0, got {h}")
    if not math.isfinite(err_norm):
        return h * _FAC_MIN
    fac = _FAC_MAX if err_norm <= 0.0 else _SAFETY * err_norm**-0.5
    fac = min(_FAC_MAX, max(_FAC_MIN, fac))
    if not accepted:
        fac = min(fac, _FAC_REJECT)
    return h * fac
