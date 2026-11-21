"""Statistical grain mode (spec D1 2.7.1, 2.7.4)."""

import importlib.util
import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

import darkroom_atlas.grain.statistical as stat
from darkroom_atlas.emulsion.classes import size_classes
from darkroom_atlas.grain.overlap import kappa
from darkroom_atlas.grain.statistical import (
    OVERLAP_STREAM,
    GrainBlockStats,
    aperture_bias,
    class_stream,
    counted_classes,
    epsilon_variance,
    grain_meta,
    grain_statistics,
    statistical_transmittance,
    statistical_transmittance_with_stats,
)
from darkroom_atlas.imaging.silver import covering_powers, imaging_constants, silver_density
from darkroom_atlas.state import EmulsionState


def _load_synthetic():
    name = "d1d_synthetic"
    if name not in sys.modules:
        path = Path(__file__).resolve().parents[1] / "imaging" / "synthetic.py"
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


syn = _load_synthetic()
K = 6
ALL = (slice(None), slice(None))


@pytest.fixture(scope="module")
def film():
    return syn.material()


def _uniform_state(mat, shape, pitch, density, x_bar, x2_ratio, *, k=K, physical=0.0, seed=0):
    """Uniform state with D_Ag (chemical) = density, X = x_bar F, Y = x2_ratio F."""
    f = syn.initiated_for_density(mat, k, density, x_bar)
    assert f <= 1.0
    return syn.state(
        mat,
        shape,
        pitch,
        k=k,
        initiated=f,
        developed=x_bar * f,
        developed_m2=x2_ratio * f,
        physical=physical,
        seed=seed,
    )


def _random_state(mat, shape, pitch=4.0, seed=0, **kw):
    rng = np.random.default_rng(seed)
    f = rng.uniform(0.0, 0.5, (*shape, K))
    x = f * rng.uniform(0.2, 1.0, (*shape, K))
    y = x * rng.uniform(0.5, 1.0, (*shape, K))
    phys = rng.uniform(0.0, 0.1, shape)
    return syn.state(
        mat, shape, pitch, k=K, initiated=f, developed=x, developed_m2=y, physical=phys, **kw
    )


def _nutting(mat, state, k=K):
    """exp(-sum_c nu_c F_c a_c x_c (1 + phi)) from first principles (spec 2.7.4)."""
    layer = mat.emulsion_layers[0]
    classes = size_classes(layer, k)
    cp = covering_powers(classes, layer.development)
    a = math.log(10.0) * 1e12 * cp * classes.silver_g_per_grain
    cap = 1e12 * classes.silver_g_per_grain * classes.number_per_um2
    silver = state.get("emulsion", "silver")[0, 0]
    f = state.get("emulsion", "initiated")[0, 0]
    phys = state.get("emulsion", "silver_physical")[0, 0]
    x_bar = (silver / cap) / f
    chem = float(np.sum(cp * silver))
    phi = 0.5 * 0.6 * phys / chem  # r_CP CP_ref of the fixture film
    return math.exp(-float(np.sum(classes.number_per_um2 * f * a * x_bar)) * (1.0 + phi))


# ---- helpers ---------------------------------------------------------------------------


def test_counted_classes():
    a = np.array([0.1, 1.0, math.pi])  # r = sqrt(a / pi): 0.178, 0.564, 1.0
    np.testing.assert_array_equal(counted_classes(a, 4.0), [True, True, True])
    np.testing.assert_array_equal(counted_classes(a, 3.9), [True, True, False])
    np.testing.assert_array_equal(counted_classes(a, 2.0), [True, False, False])
    assert counted_classes(np.array([math.pi * 0.25]), 2.0)[0]  # p == 4 r is counted
    with pytest.raises(ValueError, match="pixel_um"):
        counted_classes(a, 0.0)
    with pytest.raises(ValueError, match="non-negative"):
        counted_classes(-a, 1.0)
    with pytest.raises(ValueError, match="shape"):
        counted_classes(a[None], 1.0)


def test_epsilon_variance_formula():
    mu = np.array([[3.0, 50.0, 2.0]])
    q = np.array([[0.05, 0.01, 0.3]])
    t = np.exp(-np.sum(mu * q))
    counted = np.array([True, True, False])
    var_count = t * t * math.expm1(3.0 * 0.05**2 + 50.0 * 0.01**2)
    var_bool = np.array([5.0 * var_count])
    got = epsilon_variance(var_bool, np.array([t]), mu, q, counted)
    np.testing.assert_allclose(got, 4.0 * var_count / (t * t + var_count), rtol=1e-14)
    # Var_bool below Var_count: no excess.
    assert epsilon_variance(np.array([0.5 * var_count]), np.array([t]), mu, q, counted)[0] == 0.0
    # T = 0: 0, not nan.
    assert epsilon_variance(np.array([0.0]), np.array([0.0]), mu, q, counted)[0] == 0.0
    with pytest.raises(ValueError, match="counted"):
        epsilon_variance(var_bool, np.array([t]), mu, q, np.array([True, False]))
    with pytest.raises(ValueError, match="same shape"):
        epsilon_variance(var_bool, np.array([t]), mu, q[:, :2], counted)
    with pytest.raises(ValueError, match="non-negative"):
        epsilon_variance(-var_bool, np.array([t]), mu, q, counted)


def test_block_stats_merge_and_meta():
    a = GrainBlockStats((4, 5), 100, 2, 3.0, 6.0)
    b = GrainBlockStats((4, 5), 50, 1, 1.0, 2.0)
    m = a.merge(b)
    assert m.pixels == 150 and m.truncated_pixels == 3
    assert m.truncated_fraction == pytest.approx(0.02)
    assert m.effective_radius_um == pytest.approx(math.sqrt(8.0 / (math.pi * 4.0)))
    with pytest.raises(ValueError, match="mean-only"):
        a.merge(GrainBlockStats((5,), 1, 0, 0.0, 0.0))
    meta = grain_meta(m, aperture_um=48.0)
    assert meta["mean_only_classes"] == [4, 5]
    assert meta["b_ap"] == pytest.approx(1.15 * m.effective_radius_um / 48.0)
    assert meta["truncated_fraction"] == m.truncated_fraction
    assert grain_meta(m)["b_ap"] is None
    with pytest.warns(UserWarning, match="b_ap"):
        grain_meta(GrainBlockStats((), 1, 0, 1.0, 100.0), aperture_um=10.0)
    with pytest.raises(ValueError, match="aperture_um"):
        aperture_bias(1.0, 0.0)
    assert GrainBlockStats((), 0, 0, 0.0, 0.0).effective_radius_um == 0.0


def test_stream_names():
    assert class_stream(0) == "grain/emulsion/0"
    assert class_stream(5) == "grain/emulsion/5"
    assert OVERLAP_STREAM == "grain/emulsion/overlap"


# ---- mean and variance -----------------------------------------------------------------


@pytest.mark.parametrize("physical", [0.0, 0.3])
def test_statistical_mean_is_nutting(film, physical):
    n = 500
    s = _uniform_state(film, (n, n), 4.0, 0.5, 0.6, 0.45, physical=physical)
    t, block = statistical_transmittance_with_stats(s, *ALL)
    expected = _nutting(film, s)
    # The mean is Nutting's transmittance of D_Ag including physical silver (spec 2.5.1).
    d_ag = float(silver_density(s, slice(0, 1), slice(0, 1))[0, 0])
    assert abs(expected - 10.0**-d_ag) <= 1e-12 * expected
    se = float(np.std(t)) / math.sqrt(t.size)
    assert abs(float(np.mean(t)) - expected) <= 3.0 * se
    assert block.truncated_pixels == 0
    assert block.mean_only_classes == ()


def test_physical_silver_without_disks_is_mean_factor(film):
    s = syn.state(film, (4, 5), 4.0, initiated=0.0, developed=0.0, developed_m2=0.0, physical=0.2)
    t = statistical_transmittance(s, *ALL)
    # No chemical silver: phi = 0 and physical silver enters as 10^-(r_CP CP_ref silver_phys).
    np.testing.assert_allclose(t, 10.0 ** -(0.5 * 0.6 * 0.2), rtol=1e-14)
    g = grain_statistics(s, *ALL)
    np.testing.assert_array_equal(g.var_bool, 0.0)


@pytest.mark.parametrize("x_mode", ["broad", "one"])
def test_siedentopf_limit(film, x_mode):
    # D -> 0: Var_bool / (T^2 sum mu a^2 <x^2> / A^2) -> 1 (the ratio is kappa(D));
    # D = 0.001 within 0.1 %. All classes counted at p = 20 um.
    x_bar, x2 = (0.5, 1.0 / 3.0) if x_mode == "broad" else (1.0, 1.0)
    pitch = 20.0
    ratios = []
    for d in (0.1, 0.01, 0.001):
        s = _uniform_state(film, (1, 1), pitch, d, x_bar, x2)
        g = grain_statistics(s, *ALL)
        a = imaging_constants(s).developed_areas_um2
        area = pitch * pitch
        t = float(g.mean_t[0, 0])
        siedentopf = t * t * float(np.sum(g.mu[0, 0] * a * a * x2)) / area**2
        ratios.append(float(g.var_bool[0, 0]) / siedentopf)
        if d == 0.001:
            assert abs(ratios[-1] - 1.0) <= 1e-3
            assert abs(ratios[-1] - kappa(0.001)) <= 2e-4  # the mixture is close to kappa(D)
            if x_mode == "one":
                assert g.counted.all()
                assert abs(float(g.var_bool[0, 0]) / float(g.var_count[0, 0]) - 1.0) <= 1e-3
    assert ratios[0] > ratios[1] > ratios[2] > 1.0


def _single_class_film():
    # Class 1 of K = 2 carries the image; class 0 stays undeveloped (F = 0). CP_ref and the
    # coverage are raised so D = 1 is reachable with X = F / 2.
    return syn.material(grain={"silver_coverage": 6.0}, development={"covering_power": 1.0})


@pytest.mark.parametrize("density", [0.3, 1.0])
def test_statistical_variance_broad_x(density):
    mat = _single_class_film()
    k = 2
    classes = size_classes(mat.emulsion_layers[0], k)
    cp = covering_powers(classes, mat.emulsion_layers[0].development)
    a = math.log(10.0) * 1e12 * cp * classes.silver_g_per_grain
    pitch = math.sqrt(a[1] / 0.02)  # a_c / A_pix = 0.02, p = 12.5 r_c
    cap = syn.class_capacity(mat, k)
    f1 = density / (cp[1] * cap[1] * 0.5)
    assert f1 <= 1.0
    n = 1000
    s = syn.state(
        mat,
        (n, n),
        pitch,
        k=k,
        initiated=[0.0, f1],
        developed=[0.0, 0.5 * f1],
        developed_m2=[0.0, f1 / 3.0],
    )
    t, block = statistical_transmittance_with_stats(s, *ALL)
    g = grain_statistics(s, slice(0, 1), slice(0, 1))
    assert g.counted.all()
    assert abs(float(silver_density(s, slice(0, 1), slice(0, 1))[0, 0]) - density) <= 1e-12
    var_bool = float(g.var_bool[0, 0])
    # sigma_eps from the closed form: one radius r~ = sqrt(a <x^2>^1/2 / pi), so Var_bool is
    # kappa at the disks' own coverage D~ = mu a <x^2>^1/2 / (A ln10) times the Siedentopf
    # term. At D = 1 this is 0.203; spec 2.7.4 quotes "about 0.18", which is kappa(1) instead.
    mu = float(g.mu[0, 0, 1])
    d_tilde = mu * 0.02 * math.sqrt(1.0 / 3.0) / math.log(10.0)
    var_count_rel = math.expm1(mu * 0.01**2)
    var_bool_rel = kappa(d_tilde) * mu * 0.02**2 / 3.0
    eps_expected = (var_bool_rel - var_count_rel) / (1.0 + var_count_rel)
    assert float(g.var_eps[0, 0]) == pytest.approx(eps_expected, rel=1e-9)
    dev = t - float(np.mean(t))
    sample_var = float(np.mean(dev * dev)) * t.size / (t.size - 1)
    se = float(np.std(dev * dev)) / math.sqrt(t.size)
    assert abs(sample_var - var_bool) <= 4.0 * se
    assert block.truncated_fraction < 1e-4
    # The Siedentopf limit of Var_bool equals the random-disk-area value T^2 mu a^2 <x^2> / A^2.
    tiny = syn.state(
        mat,
        (1, 1),
        pitch,
        k=k,
        initiated=[0.0, f1 * 1e-6],
        developed=[0.0, 0.5e-6 * f1],
        developed_m2=[0.0, f1 * 1e-6 / 3.0],
    )
    gt = grain_statistics(tiny, *ALL)
    tt = float(gt.mean_t[0, 0])
    exact = tt * tt * float(gt.mu[0, 0, 1]) * a[1] ** 2 / 3.0 / pitch**4
    assert abs(float(gt.var_bool[0, 0]) / exact - 1.0) <= 1e-4


def test_mean_only_classes(film, monkeypatch):
    pitch = 2.5  # 4 r_c > p for the largest class only
    n = 500
    s = _uniform_state(film, (n, n), pitch, 0.3, 0.6, 0.45)
    a = imaging_constants(s).developed_areas_um2
    counted = counted_classes(a, pitch)
    assert not counted[-1] and counted[:-1].all()

    drawn: list[str] = []
    real_uniform = stat.pixel_uniform
    means: list[np.ndarray] = []
    real_poisson = stat.poisson_inverse

    def spy_uniform(seed, stream, *args, **kw):
        drawn.append(stream)
        return real_uniform(seed, stream, *args, **kw)

    def spy_poisson(u, mean, **kw):
        means.append(np.array(mean))
        return real_poisson(u, mean, **kw)

    monkeypatch.setattr(stat, "pixel_uniform", spy_uniform)
    monkeypatch.setattr(stat, "poisson_inverse", spy_poisson)
    t, block = statistical_transmittance_with_stats(s, *ALL)
    # No Poisson draw for the mean-only class.
    assert drawn == [class_stream(c) for c in range(K - 1)] + [OVERLAP_STREAM]
    assert len(means) == K - 1
    assert block.mean_only_classes == (K - 1,)
    # The mean is still Nutting's.
    expected = _nutting(film, s)
    se = float(np.std(t)) / math.sqrt(t.size)
    assert abs(float(np.mean(t)) - expected) <= 3.0 * se
    # Var_count leaves the mean-only class out.
    g = grain_statistics(s, slice(0, 1), slice(0, 1))
    mu, q = g.mu[0, 0], g.q[0, 0]
    tt = float(g.mean_t[0, 0])
    without = tt * tt * math.expm1(float(np.sum((mu * q * q)[:-1])))
    with_all = tt * tt * math.expm1(float(np.sum(mu * q * q)))
    assert float(g.var_count[0, 0]) == pytest.approx(without, rel=1e-13)
    assert abs(float(g.var_count[0, 0]) - with_all) > 1e-3 * with_all
    # The deterministic factor is exp(-mu q) of the mean-only class.
    assert float(g.mean_factor[0, 0]) == pytest.approx(math.exp(-mu[-1] * q[-1]), rel=1e-14)


@pytest.mark.parametrize("pitch", [6.0, 4.0])
def test_trix_prior_centre_statistical(pitch):
    trix = syn.material(
        grain={
            "morphology": "conventional",
            "mean_diameter": 0.9,
            "sigma_log": 0.45,
            "silver_coverage": 4.0,
            "iodide_fraction": 0.04,
        },
        development={"covering_power": 0.6},
    )
    s = _uniform_state(trix, (64, 64), pitch, 1.0, 0.8, 0.7)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t, block = statistical_transmittance_with_stats(s, *ALL)
        meta = grain_meta(block, aperture_um=48.0)
    assert t.shape == (64, 64) and np.all((t >= 0.0) & (t <= 1.0))
    for key in ("mean_only_classes", "b_ap", "truncated_fraction"):
        assert key in meta
    a = imaging_constants(s).developed_areas_um2
    assert meta["mean_only_classes"] == [int(c) for c in np.flatnonzero(~counted_classes(a, pitch))]
    assert meta["mean_only_classes"]  # the coarse Tri-X classes are mean-only at 4-6 um
    assert 0.0 < meta["b_ap"] and 0.0 <= meta["truncated_fraction"] <= 1.0


# ---- determinism and tiling ------------------------------------------------------------


def test_statistical_deterministic(film):
    s = _random_state(film, (23, 29))
    a = statistical_transmittance(s, *ALL)
    b = statistical_transmittance(s, *ALL)
    assert np.array_equal(a, b)
    assert not a.flags.writeable
    other = EmulsionState(s.material, s.geometry, s.grid, 1, s.fields, s.history, s.meta)
    assert not np.array_equal(a, statistical_transmittance(other, *ALL))


def test_statistical_tiling_invariance(film):
    # 100 x 100 > one variance chunk, so chunk boundaries fall inside windows.
    shape = (100, 100)
    s = _random_state(film, shape)
    full = statistical_transmittance(s, *ALL)
    # Single pixels (a corner and a stretch across the variance-chunk boundary at 8192).
    for r, c in [*((i, j) for i in range(6) for j in range(6)), *((81, j) for j in range(88, 96))]:
        assert statistical_transmittance(s, slice(r, r + 1), slice(c, c + 1))[0, 0] == full[r, c]
    for tile in (7, 32, 64):
        pieced = np.empty(shape)
        for r in range(0, shape[0], tile):
            for c in range(0, shape[1], tile):
                rows, cols = slice(r, min(r + tile, shape[0])), slice(c, min(c + tile, shape[1]))
                if tile == 32:  # through a cropped state: origin offsets
                    pieced[rows, cols] = statistical_transmittance(s.crop((rows, cols)), *ALL)
                else:
                    pieced[rows, cols] = statistical_transmittance(s, rows, cols)
        assert np.array_equal(pieced, full), tile


def test_uniform_window_matches_mixed_window(film):
    # A uniform window takes the one-evaluation path; the same pixels in a non-uniform
    # window take the chunked path. Both give the same bits.
    rng = np.random.default_rng(4)
    shape = (40, 30)
    f = np.full((*shape, K), 0.3)
    f[20:] = rng.uniform(0.0, 0.5, (20, 30, K))
    s = syn.state(film, shape, 4.0, initiated=f, developed=0.6 * f, developed_m2=0.45 * f)
    full = statistical_transmittance(s, *ALL)
    top = statistical_transmittance(s, slice(0, 20), slice(None))
    assert np.array_equal(full[:20], top)
    g_full = grain_statistics(s, *ALL)
    g_top = grain_statistics(s, slice(0, 20), slice(None))
    assert np.array_equal(g_full.var_bool[:20], g_top.var_bool)


def test_statistical_independent_of_alignment(film):
    s = _random_state(film, (17, 19))
    ref = statistical_transmittance(s, *ALL)
    shifted = {}
    for key, value in s.fields.items():
        if value.strides[:2] == (0, 0):
            shifted[key] = value
            continue
        # A copy whose data start 8 bytes past the start of a fresh allocation (spec 3.6).
        raw = np.empty(value.size + 1, dtype=np.float64)
        view = raw[1:].reshape(value.shape)
        assert view.ctypes.data % 16 != value.ctypes.data % 16 or view.ctypes.data % 16 == 8
        view[...] = value
        view.flags.writeable = False
        shifted[key] = view
    moved = s.with_fields(shifted)
    assert np.array_equal(statistical_transmittance(moved, *ALL), ref)


# ---- bad input -------------------------------------------------------------------------


def test_bad_fields(film):
    s = _random_state(film, (5, 6))
    with pytest.raises(ValueError, match="initiated"):
        statistical_transmittance(
            s.with_fields({"emulsion/initiated": np.full((5, 6, K), 1.5)}), *ALL
        )
    with pytest.raises(ValueError, match="non-negative"):
        statistical_transmittance(
            s.with_fields({"emulsion/silver_m2": np.full((5, 6, K), -0.1)}), *ALL
        )
    with pytest.raises(ValueError, match="must match"):
        statistical_transmittance(
            s.with_fields({"emulsion/initiated": np.full((5, 6, K - 1), 0.1)}), *ALL
        )
    with pytest.raises(ValueError, match="finite"):
        statistical_transmittance(
            s.with_fields({"emulsion/silver_physical": np.full((5, 6), np.nan)}), *ALL
        )
    with pytest.raises(KeyError):
        statistical_transmittance(s.drop("emulsion/silver_m2"), *ALL)
    with pytest.raises(ValueError, match="empty"):
        statistical_transmittance(s, slice(2, 2), slice(None))
