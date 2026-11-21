import itertools
import math

import numpy as np
import pytest
from scipy.linalg import expm

from darkroom_atlas.develop.rosenbrock import (
    GAMMA,
    MIN_SCALE_MOL_L,
    NEGATIVE_REJECT,
    NEGATIVE_WARN,
    ColumnJacobian,
    ColumnLayout,
    adaptive_accept,
    concentration_scale,
    error_norm,
    factor_columns,
    min_scaled_concentration,
    next_step_size,
    ros2_step,
    ros2_step_columns,
    solve_columns,
    toth_transport_jacobian,
)

# A small development-like column: 2 dense (emulsion) levels with 2 components and 1 local
# unknown, then 3 diagonal levels (overcoat, liquid), Dirichlet at the top face.
LAYOUT = ColumnLayout(n_dense=2, n_diag=3, n_comp=2, n_local=1)
THICK = np.array([2.0, 2.0, 1.0, 20.0, 60.0])  # um, bottom to top
DIFF = np.array([[50.0, 20.0], [50.0, 20.0], [60.0, 30.0], [800.0, 600.0], [800.0, 600.0]])


def _shifted_copy(array: np.ndarray, offset: int) -> np.ndarray:
    raw = np.empty(array.nbytes + 64, dtype=np.uint8)
    start = (-raw.ctypes.data) % 64 + offset
    out = raw[start : start + array.nbytes].view(array.dtype).reshape(array.shape)
    out[...] = array
    return out


def _transport(thick, diff):
    """Finite-volume transport couplings (spec D1 2.4.8) with a top Dirichlet face."""
    n_lev = thick.size
    h_face = 0.5 * (thick[:-1] + thick[1:])
    d_face = (thick[:-1] + thick[1:])[:, None] / (
        thick[:-1, None] / diff[:-1] + thick[1:, None] / diff[1:]
    )
    g = d_face / h_face[:, None]  # (L-1, K) conductance
    lower = g / thick[1:, None]  # d f(l+1) / d C(l)
    upper = g / thick[:-1, None]  # d f(l) / d C(l+1)
    self_ = np.zeros((n_lev, diff.shape[1]))
    self_[:-1] -= upper
    self_[1:] -= lower
    g_top = diff[-1] / (0.5 * thick[-1])
    self_[-1] -= g_top / thick[-1]
    return lower, upper, self_, g_top / thick[-1]


def _reaction_block():
    # Components A, B and a local product X: A + B -> X at a linearised rate.
    r = np.zeros((3, 3))
    r[0, 0] = -3.0
    r[1, 0] = -1.5
    r[1, 1] = -0.5
    r[2, 0] = 3.0
    r[0, 1] = 0.2
    return r


def _linear_jacobian(n_cols=1, react_scale=1.0):
    lower, upper, self_, _ = _transport(THICK, DIFF)
    dense = np.zeros((n_cols, 2, 3, 3))
    for lev in range(2):
        dense[:, lev] = react_scale * _reaction_block()
        dense[:, lev, 0, 0] += self_[lev, 0]
        dense[:, lev, 1, 1] += self_[lev, 1]
    diag = np.broadcast_to(self_[2:], (n_cols, 3, 2)).copy()
    return ColumnJacobian(LAYOUT, dense, diag, lower, upper)


def _bath(t):
    return np.array([1.0 + 0.5 * math.sin(2.0 * t), 0.4 + 0.2 * math.cos(3.0 * t)])


def _linear_rhs(jac, top):
    full = jac.to_dense()

    def f(t, y):
        out = np.einsum("nij,nj->ni", full, y)
        _, diag = LAYOUT.split(out)
        diag[:, -1] += top * _bath(t)
        return out

    return f


def _initial(n_cols=1):
    y = np.zeros((n_cols, LAYOUT.size))
    _, diag = LAYOUT.split(y)
    diag[:, -1] = _bath(0.0)
    diag[:, -2] = 0.5 * _bath(0.0)
    return y


def _integrate(f, jac, y0, t0, t_end, n_steps):
    """Fixed steps; ``jac`` is a ColumnJacobian or a callable y -> ColumnJacobian."""
    h = (t_end - t0) / n_steps
    y = y0
    for i in range(n_steps):
        j = jac(y) if callable(jac) else jac
        y, _ = ros2_step_columns(f, t0 + i * h, y, h, j)
    return y


def _ratios(solve_n, ns):
    sols = [solve_n(n) for n in ns]
    diffs = [np.max(np.abs(a - b)) for a, b in itertools.pairwise(sols)]
    return [d0 / d1 for d0, d1 in itertools.pairwise(diffs)]


# ---- order -----------------------------------------------------------------------------------

# Order is measured on [1, 3] s from the state reached at t = 1 s (1000 fine steps with the
# exact J), so the initial stiff layer of the dry start has decayed and the solution is smooth.
T0, T1 = 1.0, 3.0
STEPS = [40, 80, 160, 320]


def _offdiag_dropped():
    j = _linear_jacobian()
    d = j.dense.copy()
    d[:, :, 2, 0] = d[:, :, 0, 1] = d[:, :, 1, 0] = 0.0
    return ColumnJacobian(LAYOUT, d, j.diag, j.lower, j.upper)


@pytest.mark.parametrize("jac_kind", ["exact", "half_reaction", "no_cross_terms"])
def test_rosenbrock_second_order(jac_kind):
    # Spec D1 2.4.13: fixed steps halved repeatedly, ratio of successive differences 3.5-4.5.
    # Stiff linear column (transport eigenvalues up to ~1e2 /s, h up to 0.05 s) with a
    # time-dependent Dirichlet value. The inexact Jacobians (half the reaction block; the
    # cross-component reaction terms dropped, as spec D1 3.3 drops couplings) check the
    # W-method property.
    exact = _linear_jacobian()
    _, _, _, top = _transport(THICK, DIFF)
    f = _linear_rhs(exact, top)
    jac = {
        "exact": exact,
        "half_reaction": _linear_jacobian(react_scale=0.5),
        "no_cross_terms": _offdiag_dropped(),
    }[jac_kind]
    assert np.max(np.abs(np.linalg.eigvals(exact.to_dense()[0]))) > 50.0
    y0 = _integrate(f, exact, _initial(), 0.0, T0, 1000)
    ratios = _ratios(lambda n: _integrate(f, jac, y0, T0, T1, n), STEPS)
    for r in ratios:
        assert 3.5 <= r <= 4.5, ratios


def _nonlinear_problem(diagonal_only=False):
    # Prototype of spec D1 12.6: r = k A B / (1 + K_B B) (1 - X) in the dense levels.
    lower, upper, self_, top = _transport(THICK, DIFF)
    k_rate, k_b = 4.0, 2.0

    def f(t, y):
        dense, diag = LAYOUT.split(y)
        c_all = np.concatenate([dense[:, :, :2], diag], axis=1)  # (N, L, 2)
        out_c = self_[None] * c_all
        out_c[:, 1:] += lower[None] * c_all[:, :-1]
        out_c[:, :-1] += upper[None] * c_all[:, 1:]
        out_c[:, -1] += top * _bath(t)
        a, b, x = dense[:, :, 0], dense[:, :, 1], dense[:, :, 2]
        r = k_rate * a * b / (1.0 + k_b * b) * (1.0 - x)
        out_d = np.empty_like(dense)
        out_d[:, :, 0] = out_c[:, :2, 0] - r
        out_d[:, :, 1] = out_c[:, :2, 1] - 0.5 * r
        out_d[:, :, 2] = 0.3 * r
        return LAYOUT.join(out_d, out_c[:, 2:])

    def jacobian(y):
        # Transport exact; reaction Jacobian exact, or only its diagonal entries.
        dense, _ = LAYOUT.split(y)
        a, b, x = dense[:, :, 0], dense[:, :, 1], dense[:, :, 2]
        n = y.shape[0]
        g = 1.0 + k_b * b
        dr = np.stack(  # dr/dA, dr/dB, dr/dX
            [k_rate * b / g * (1.0 - x), k_rate * a / g**2 * (1.0 - x), -k_rate * a * b / g],
            axis=-1,
        )
        jd = np.zeros((n, 2, 3, 3))
        for row, coef in ((0, -1.0), (1, -0.5), (2, 0.3)):
            jd[:, :, row, :] = coef * dr
        if diagonal_only:
            jd = jd * np.eye(3)
        jd[:, :, 0, 0] += self_[:2, 0]
        jd[:, :, 1, 1] += self_[:2, 1]
        diag = np.broadcast_to(self_[2:], (n, 3, 2))
        return ColumnJacobian(LAYOUT, jd, diag, lower, upper)

    return f, jacobian


@pytest.mark.parametrize("diagonal_only", [False, True])
def test_rosenbrock_second_order_nonlinear(diagonal_only):
    f, jacobian = _nonlinear_problem(diagonal_only)
    # Coupled to the stiff transport, the error approaches its h^2 law from below (3.1, 3.3,
    # 3.6 at 20-160 steps), so this problem is measured from 160 steps on.
    y0 = _integrate(f, jacobian, _initial(), 0.0, T0, 1000)
    ratios = _ratios(lambda n: _integrate(f, jacobian, y0, T0, T1, n), [160, 320, 640, 1280])
    for r in ratios:
        assert 3.5 <= r <= 4.5, ratios


# ---- stability -------------------------------------------------------------------------------


def _stability(z):
    """Closed form R(z) = (1 + (1 - 2 gamma) z) / (1 - gamma z)^2 of ROS2."""
    return (1.0 + (1.0 - 2.0 * GAMMA) * z) / (1.0 - GAMMA * z) ** 2


@pytest.mark.parametrize("z", [-1e-3, -0.5, -3.0, -40.0, -1e4, -1e9])
def test_rosenbrock_l_stable(z):
    # One step of y' = lam y with J = lam matches R(z) and |R| <= 1; R -> 0 as z -> -inf.
    lam = -2.0
    h = z / lam
    y0 = np.array([[1.0]])

    def solve(rhs):
        return rhs / (1.0 - GAMMA * h * lam)

    y1, _ = ros2_step(lambda t, y: lam * y, 0.0, y0, h, solve)
    assert y1[0, 0] == pytest.approx(_stability(z), rel=1e-9, abs=1e-15)
    assert abs(y1[0, 0]) <= 1.0
    if z <= -1e9:
        assert abs(y1[0, 0]) < 1e-9


def test_rosenbrock_a_stable_on_left_half_plane():
    re = -np.logspace(-4, 6, 60)
    z = re[:, None] + 1j * np.concatenate([[0.0], np.logspace(-4, 6, 60)])[None, :]
    assert np.max(np.abs(_stability(z))) <= 1.0 + 1e-14
    assert np.max(np.abs(_stability(1j * np.logspace(-4, 6, 200)))) <= 1.0 + 1e-14


def test_rosenbrock_damps_stiff_mode():
    # L-stability in a column: the stiffest transport mode, taken with h = 1e4 / |lambda|, is
    # multiplied by R(h lambda) ~ 1e-4 in one step (exact J, constant bath 0).
    lower, upper, self_, _ = _transport(THICK, DIFF)
    dense = np.zeros((1, 2, 3, 3))
    for lev in range(2):
        dense[:, lev, 0, 0] = self_[lev, 0]
        dense[:, lev, 1, 1] = self_[lev, 1]
    jac = ColumnJacobian(LAYOUT, dense, self_[None, 2:], lower, upper)
    full = jac.to_dense()[0]
    lam, vec = np.linalg.eig(full)
    i = int(np.argmax(np.abs(lam)))
    assert abs(lam[i].imag) == 0.0
    y0 = vec[:, i].real[None]
    h = 1e4 / abs(lam[i].real)
    y1, _ = ros2_step_columns(lambda t, y: np.einsum("ij,nj->ni", full, y), 0.0, y0, h, jac)
    ratio = np.max(np.abs(y1)) / np.max(np.abs(y0))
    assert ratio == pytest.approx(abs(_stability(h * lam[i].real)), rel=1e-6)
    assert ratio < 1e-4


def test_diffusion_integrator():
    # Pure vertical diffusion: ROS2 on the column ODE converges at order 2 to expm(t J) y0.
    lower, upper, self_, top = _transport(THICK, DIFF)
    lay = ColumnLayout(n_dense=0, n_diag=5, n_comp=2)
    n = 1
    jac = ColumnJacobian(
        lay, np.zeros((n, 0, 2, 2)), np.broadcast_to(self_, (n, 5, 2)), lower, upper
    )
    full = jac.to_dense()[0]
    src = np.zeros(lay.size)
    src[-2:] = top * np.array([1.0, 0.5])

    def f(t, y):
        return np.einsum("ij,nj->ni", full, y) + src

    t_end = 1.5
    y0 = np.zeros((n, lay.size))
    # Exact: y(t) = expm(tJ) y0 + J^-1 (expm(tJ) - I) src.
    e = expm(t_end * full)
    exact = e @ y0[0] + np.linalg.solve(full, (e - np.eye(lay.size)) @ src)
    errs = []
    for steps in (60, 120, 240):
        h = t_end / steps
        y = y0
        for i in range(steps):
            y, _ = ros2_step_columns(f, i * h, y, h, jac)
        errs.append(np.max(np.abs(y[0] - exact)))
    assert errs[-1] < 1e-4
    for e0, e1 in itertools.pairwise(errs):
        assert 3.5 <= e0 / e1 <= 4.5, errs


def _slab_series_cell_means(n_z, tau, terms=4000):
    """Cell means of the finite-slab series solution, base zero flux, C0 = 1 at the top."""
    edges = np.arange(n_z + 1) / n_z  # fractional height from the base
    k = (2 * np.arange(terms)[:, None] + 1) * math.pi / 2
    sign = (-1.0) ** np.arange(terms)[:, None]
    s = np.sin(k * edges[None, 1:]) - np.sin(k * edges[None, :-1])
    return 1.0 - np.sum(4.0 * sign / (2.0 * k) * s / (k / n_z) * np.exp(-k * k * tau), axis=0)


def test_diffusion_into_slab():
    # Spec D1 2.4.13 (num16): no reaction, s = 1, constant top value, emulsion n_z = 8/16/32
    # (dense levels, one component); cell means against the finite-slab series. The error
    # ratio between refinements is >= 3 at tD/L^2 in {0.01, 0.1, 0.3}; at n_z = 8 the error is
    # <= 2e-3 C0 for tD/L^2 >= 1. The error is the max over cells of |C - cell mean of the
    # series|, which reproduces the spec's n_z = 8 values (7.3e-2 at 0.01, 3.3e-3 at 0.3).
    length, d = 10.0, 1.0  # um, um^2/s: t = tau L^2 / D
    taus = (0.01, 0.1, 0.3, 1.0)
    errors = {}
    for n_z in (8, 16, 32):
        thick = np.full(n_z, length / n_z)
        lower, upper, self_, top = _transport(thick, np.full((n_z, 1), d))
        lay = ColumnLayout(n_dense=n_z, n_diag=0, n_comp=1)
        jac = ColumnJacobian(
            lay, self_[None, :, :, None], np.zeros((1, 0, 1)), lower, upper
        )  # fmt: skip
        full = jac.to_dense()[0]
        src = np.zeros(lay.size)
        src[-1] = top[0]

        def f(t, y, full=full, src=src):
            return np.einsum("ij,nj->ni", full, y) + src

        y, t, h = np.zeros((1, lay.size)), 0.0, 1e-4
        for tau in taus:
            t_out = tau * length**2 / d
            while t < t_out:
                step = min(h, t_out - t)
                y, _ = ros2_step_columns(f, t, y, step, jac)
                t += step
                h = min(1.03 * h, 0.1)
            errors[n_z, tau] = np.max(np.abs(y[0] - _slab_series_cell_means(n_z, tau)))
    assert errors[8, 0.01] == pytest.approx(7.3e-2, abs=1e-3)
    assert errors[8, 0.3] == pytest.approx(3.3e-3, abs=1e-4)
    for tau in (0.01, 0.1, 0.3):
        assert errors[8, tau] / errors[16, tau] >= 3.0, errors
        assert errors[16, tau] / errors[32, tau] >= 3.0, errors
    assert errors[8, 1.0] <= 2e-3


# ---- TOTH transport Jacobian (spec D1 2.4.8, 3.3) -------------------------------------------


def _toth_column():
    from darkroom_atlas.chemistry.speciation import make_speciation
    from darkroom_atlas.chemistry.species import AcidBaseFamily
    from darkroom_atlas.data.models import ActivityParams

    q = {"value": 0.0, "unit": "1", "source": "assumed", "note": "test"}
    params = ActivityParams(
        debye_huckel_A_20C=q,
        debye_huckel_A_25C=q,
        sit_b={"value": 1.5, "unit": "1", "source": "assumed", "note": "test"},
        delta_eps={},
    )
    sulfite = AcidBaseFamily(
        "sulfite", (0, -1, -2), (2, 1, 0), (1.86, 7.19), "van_t_hoff", (0.0, 0.0), (), "x"
    )
    hq = AcidBaseFamily(
        "hydroquinone", (0, -1, -2), (0, -1, -2), (9.85, 11.45), "van_t_hoff", (0.0, 0.0), (), "x"
    )
    gel = AcidBaseFamily("gelatin", (1, 0), (1, 0), (10.5,), "van_t_hoff", (0.0,), (), "none")
    spec = make_speciation([sulfite, hq, gel], 293.15, 0.0, 0.0, params)
    # 4 levels: 2 gel cells (with gelatin), overcoat, liquid.
    totals = np.array([[0.6, 0.62, 0.7, 0.79], [0.03, 0.035, 0.04, 0.045], [0.02, 0.02, 0.01, 0.0]])
    ph = np.array([8.2, 8.3, 8.4, 8.5])
    thick = np.array([5.0, 5.0, 3.0, 20.0])
    # Face D per carrier (sulfite, HQ, H+, OH-): gel faces smaller.
    face_d = np.array([[500.0, 600.0, 900.0], [80.0, 100.0, 700.0], [4000.0, 5000.0, 8130.0],
                       [2300.0, 3000.0, 4600.0]])  # fmt: skip
    return spec, totals, ph, thick, face_d


def test_toth_transport_jacobian_matches_finite_difference():
    from darkroom_atlas.chemistry.speciation import (
        buffer_capacity,
        mean_protons,
        solve_ph,
        water_ions,
    )
    from darkroom_atlas.develop.lateral import carrier_capacities

    spec, totals, ph, thick, face_d = _toth_column()
    nbar = mean_protons(spec, ph)
    h_ion, oh_ion = water_ions(spec, ph)
    toth = (totals * nbar).sum(axis=0) + h_ion - oh_ion
    h_face = 0.5 * (thick[:-1] + thick[1:])

    def f_toth(toth_vec):
        p = solve_ph(spec, totals, toth_vec, ph0=ph)
        pf = totals[:2] * mean_protons(spec, p)[:2]  # mobile families only
        hh, oo = water_ions(spec, p)
        flux = np.zeros(thick.size - 1)
        for c, x in enumerate((pf[0], pf[1], hh)):
            flux -= face_d[c] * (x[1:] - x[:-1]) / h_face
        flux += face_d[3] * (oo[1:] - oo[:-1]) / h_face
        out = np.zeros(thick.size)  # zero flux at both ends
        out[:-1] -= flux / thick[:-1]
        out[1:] += flux / thick[1:]
        return out

    beta_c = carrier_capacities(spec, totals, ph, [True, True, False])  # (4, L)
    beta = buffer_capacity(spec, totals, ph)
    lower, upper, self_ = toth_transport_jacobian(
        beta_c[:, None, :], beta[None], face_d[:, None, :], h_face, thick
    )
    fd = np.zeros((thick.size, thick.size))
    delta = 1e-7
    for j in range(thick.size):
        e = np.zeros(thick.size)
        e[j] = delta
        fd[:, j] = (f_toth(toth + e) - f_toth(toth - e)) / (2 * delta)
    scale = np.max(np.abs(fd))
    np.testing.assert_allclose(np.diag(fd), self_[0], rtol=1e-6, atol=1e-9 * scale)
    np.testing.assert_allclose(np.diag(fd, -1), lower[0], rtol=1e-6, atol=1e-9 * scale)
    np.testing.assert_allclose(np.diag(fd, 1), upper[0], rtol=1e-6, atol=1e-9 * scale)
    # Only the tridiagonal part is present (no longer-range coupling).
    assert np.max(np.abs(np.triu(fd, 2))) <= 1e-9 * scale


def test_toth_transport_jacobian_rejects_bad_input():
    *_, thick, face_d = _toth_column()
    h_face = 0.5 * (thick[:-1] + thick[1:])
    beta_c = np.ones((4, 1, 4))
    with pytest.raises(ValueError, match="beta_total"):
        toth_transport_jacobian(beta_c, np.zeros((1, 4)), face_d[:, None], h_face, thick)
    with pytest.raises(ValueError, match="shape"):
        toth_transport_jacobian(beta_c, np.ones((1, 3)), face_d[:, None], h_face, thick)
    with pytest.raises(ValueError, match="face"):
        toth_transport_jacobian(beta_c, np.ones((1, 4)), face_d[:, None], h_face[:2], thick)


# ---- block-Thomas ----------------------------------------------------------------------------


def _random_jacobian(rng, layout, n_cols):
    b, nc = layout.block, layout.n_comp
    dense = rng.normal(size=(n_cols, layout.n_dense, b, b)) - 3.0 * np.eye(b)
    diag = -np.abs(rng.normal(size=(n_cols, layout.n_diag, nc))) * 3.0
    lower = np.abs(rng.normal(size=(n_cols, layout.n_levels - 1, nc)))
    upper = np.abs(rng.normal(size=(n_cols, layout.n_levels - 1, nc)))
    return ColumnJacobian(layout, dense, diag, lower, upper)


@pytest.mark.parametrize(
    "layout",
    [
        ColumnLayout(n_dense=4, n_diag=9, n_comp=7, n_local=2),
        ColumnLayout(n_dense=2, n_diag=0, n_comp=3, n_local=1),
        ColumnLayout(n_dense=0, n_diag=6, n_comp=4),
        ColumnLayout(n_dense=1, n_diag=1, n_comp=1),
        ColumnLayout(n_dense=3, n_diag=2, n_comp=5),
    ],
)
def test_block_thomas_matches_dense_solve(layout):
    rng = np.random.default_rng(7)
    n_cols = 25
    jac = _random_jacobian(rng, layout, n_cols)
    gamma_h = 0.8
    update = rng.normal(size=(n_cols, layout.n_dense, layout.block, layout.block)) * 0.1
    factor = factor_columns(jac, gamma_h, dense_update=update)
    rhs = rng.normal(size=(n_cols, layout.size))
    x = solve_columns(factor, rhs)
    m = np.eye(layout.size)[None] - gamma_h * jac.to_dense()
    off = 0
    for lev in range(layout.n_dense):
        b = layout.block
        m[:, off : off + b, off : off + b] += update[:, lev]
        off += b
    ref = np.linalg.solve(m, rhs[:, :, None])[:, :, 0]
    assert np.max(np.abs(x - ref)) <= 1e-12 * np.max(np.abs(ref))
    # Residual check as well.
    res = np.einsum("nij,nj->ni", m, x) - rhs
    assert np.max(np.abs(res)) <= 1e-12 * np.max(np.abs(rhs))


def test_block_thomas_pivots():
    # J = [[1, -1], [-1, 0]] and gamma h = 1 give M = [[0, 1], [1, 1]]: needs a row swap.
    lay = ColumnLayout(n_dense=1, n_diag=0, n_comp=2)
    jac = ColumnJacobian(
        lay,
        np.array([[[[1.0, -1.0], [-1.0, 0.0]]]]),
        np.zeros((1, 0, 2)),
        np.zeros((0, 2)),
        np.zeros((0, 2)),
    )
    x = solve_columns(factor_columns(jac, 1.0), np.array([[1.0, 2.0]]))
    np.testing.assert_allclose(x[0], [1.0, 1.0], rtol=1e-15)


def test_block_thomas_singular_raises():
    lay = ColumnLayout(n_dense=1, n_diag=0, n_comp=2)
    jac = ColumnJacobian(
        lay,
        np.array([[[[1.0, 1.0], [1.0, 1.0]]]]),
        np.zeros((1, 0, 2)),
        np.zeros((0, 2)),
        np.zeros((0, 2)),
    )
    with pytest.raises(ValueError, match="singular"):
        factor_columns(jac, 0.5)  # M = [[0.5, -0.5], [-0.5, 0.5]]


def test_block_thomas_batch_bit_identical():
    # Spec D1 2.4.13 test_ensemble_vectorisation, at the solver level: a batch equals the
    # columns solved one at a time, bit for bit.
    rng = np.random.default_rng(3)
    lay = ColumnLayout(n_dense=3, n_diag=5, n_comp=4, n_local=2)
    n_cols = 37
    jac = _random_jacobian(rng, lay, n_cols)
    rhs = rng.normal(size=(n_cols, lay.size))
    whole = solve_columns(factor_columns(jac, 0.7), rhs)
    for i in (0, 5, 36):
        one = ColumnJacobian(
            lay, jac.dense[i : i + 1], jac.diag[i : i + 1], jac.lower[i : i + 1],
            jac.upper[i : i + 1],
        )  # fmt: skip
        single = solve_columns(factor_columns(one, 0.7), rhs[i : i + 1])
        assert np.array_equal(single[0], whole[i])


def test_block_result_independent_of_alignment():
    # Spec D1 3.6: the same block from copies of its inputs misaligned by 8 bytes is
    # bit-identical (factor, solve and a full ROS2 step).
    rng = np.random.default_rng(11)
    lay = ColumnLayout(n_dense=3, n_diag=6, n_comp=5, n_local=2)
    n_cols = 19
    jac = _random_jacobian(rng, lay, n_cols)
    rhs = rng.normal(size=(n_cols, lay.size))
    y0 = rng.normal(size=(n_cols, lay.size))

    def run(offset):
        j = ColumnJacobian(
            lay,
            _shifted_copy(jac.dense, offset),
            _shifted_copy(jac.diag, offset),
            _shifted_copy(np.ascontiguousarray(jac.lower), offset),
            _shifted_copy(np.ascontiguousarray(jac.upper), offset),
        )
        full = j.to_dense()
        x = solve_columns(factor_columns(j, 0.4), _shifted_copy(rhs, offset))
        y1, err = ros2_step_columns(
            lambda t, y: np.einsum("nij,nj->ni", full, y) + math.sin(t),
            0.0,
            _shifted_copy(y0, offset),
            0.3,
            j,
        )
        return x, y1, err

    for a, b in zip(run(0), run(8), strict=True):
        assert np.array_equal(a, b)


# ---- negative concentrations and step control (spec D1 3.3) ---------------------------------


def test_constants_match_spec():
    assert GAMMA == 1.0 + 1.0 / math.sqrt(2.0)
    assert MIN_SCALE_MOL_L == 1e-3
    assert NEGATIVE_REJECT == -1e-6
    assert NEGATIVE_WARN == -1e-3


def test_concentration_scale_floor():
    s = concentration_scale(np.array([0.0, 5e-4, 0.1]))
    np.testing.assert_array_equal(s, [1e-3, 1e-3, 0.1])


def test_min_scaled_concentration():
    conc = np.array([[[0.05, -2e-9], [0.1, 1e-3]], [[0.2, 0.0], [-1e-5, 2e-3]]])
    scale = np.array([0.1, 1e-3])
    assert min_scaled_concentration(conc, scale) == pytest.approx(-1e-4)


def test_adaptive_accept_rejects_negative_concentration():
    assert adaptive_accept(0.5, -5e-7)
    assert not adaptive_accept(0.5, -2e-6)  # treated as a failed error test
    assert not adaptive_accept(1.5, 0.0)
    assert adaptive_accept(1.0, NEGATIVE_REJECT)
    assert not adaptive_accept(float("nan"), 0.0)


def test_next_step_size():
    assert next_step_size(1.0, 0.25, True) == pytest.approx(1.8)
    assert next_step_size(1.0, 0.0, True) == pytest.approx(5.0)
    assert next_step_size(1.0, 100.0, False) == pytest.approx(0.2)
    # A negativity-only rejection (small error) still shrinks the step.
    assert next_step_size(1.0, 0.01, False) == pytest.approx(0.5)
    assert next_step_size(1.0, float("inf"), False) == pytest.approx(0.2)


def test_error_norm_is_batch_max():
    err = np.array([[1e-5, 0.0], [0.0, 3e-4]])
    y = np.array([[1.0, 0.0], [0.0, 1.0]])
    assert error_norm(err, y, y, rtol=1e-4, scale=np.array([1e-3, 1e-3])) == pytest.approx(3.0)


def test_fixed_step_never_rejects_negative_values():
    # A step that overshoots below zero is returned as is: no rejection, split or clipping.
    y0 = np.array([[0.01]])
    y1, _ = ros2_step(lambda t, y: np.full_like(y, -1.0), 0.0, y0, 1.0, lambda r: r)
    assert y1[0, 0] == pytest.approx(0.01 - 1.5 + 0.5)


# ---- bad input -------------------------------------------------------------------------------


def test_layout_rejects_bad_counts():
    with pytest.raises(ValueError):
        ColumnLayout(n_dense=-1, n_diag=2, n_comp=1)
    with pytest.raises(ValueError):
        ColumnLayout(n_dense=1, n_diag=2, n_comp=0)
    with pytest.raises(ValueError):
        ColumnLayout(n_dense=0, n_diag=0, n_comp=2)
    with pytest.raises(ValueError):
        ColumnLayout(n_dense=1.0, n_diag=2, n_comp=1)
    with pytest.raises(ValueError, match="shape"):
        LAYOUT.split(np.zeros((2, LAYOUT.size + 1)))


def test_jacobian_rejects_bad_shapes_and_values():
    good = _linear_jacobian(2)
    with pytest.raises(ValueError, match="dense"):
        ColumnJacobian(LAYOUT, good.dense[:, :1], good.diag, good.lower, good.upper)
    with pytest.raises(ValueError, match="diag"):
        ColumnJacobian(LAYOUT, good.dense, good.diag[:1], good.lower, good.upper)
    with pytest.raises(ValueError, match="lower"):
        ColumnJacobian(LAYOUT, good.dense, good.diag, good.lower[:, :2], good.upper)
    bad = good.dense.copy()
    bad[0, 0, 0, 0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        ColumnJacobian(LAYOUT, bad, good.diag, good.lower, good.upper)
    with pytest.raises(TypeError):
        ColumnJacobian("layout", good.dense, good.diag, good.lower, good.upper)


def test_step_rejects_bad_input():
    jac = _linear_jacobian()
    f = _linear_rhs(jac, 1.0)
    with pytest.raises(ValueError, match="h must"):
        ros2_step_columns(f, 0.0, _initial(), 0.0, jac)
    with pytest.raises(ValueError, match="h must"):
        ros2_step_columns(f, 0.0, _initial(), float("nan"), jac)
    with pytest.raises(ValueError, match="gamma_h"):
        factor_columns(jac, -1.0)
    with pytest.raises(ValueError, match="shape"):
        ros2_step(lambda t, y: y[:, :1], 0.0, _initial(), 0.1, lambda r: r)
    with pytest.raises(ValueError, match="columns"):
        solve_columns(factor_columns(jac, 0.1), np.zeros((3, LAYOUT.size)))
    with pytest.raises(ValueError):
        min_scaled_concentration(np.zeros((2, 3)), np.ones(2))
    with pytest.raises(ValueError):
        concentration_scale(np.array([[1.0]]))
    with pytest.raises(ValueError):
        error_norm(np.zeros(2), np.zeros(2), np.zeros(2), rtol=0.0, scale=1.0)
