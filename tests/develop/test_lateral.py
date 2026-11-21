"""Exact lateral diffusion (spec D1 3.5, 3.6)."""

import math

import numpy as np
import pytest

from darkroom_atlas.chemistry.speciation import buffer_capacity, make_speciation
from darkroom_atlas.chemistry.species import AcidBaseFamily
from darkroom_atlas.data.models import ActivityParams
from darkroom_atlas.develop.lateral import (
    BOUNDARIES,
    TothLateral,
    carrier_capacities,
    lateral_diffuse,
    lateral_factors,
    toth_lateral_diffusion,
)

DX = 24.0  # um, the image chemistry cell (spec D1 3.1)
# Per level (3) and component (2) diffusivities, um^2/s.
DIFF = np.array([[400.0, 50.0], [900.0, 0.0], [1816.0, 1030.0]])


def _shifted_copy(array: np.ndarray, offset: int) -> np.ndarray:
    raw = np.empty(array.nbytes + 64, dtype=np.uint8)
    start = (-raw.ctypes.data) % 64 + offset
    out = raw[start : start + array.nbytes].view(array.dtype).reshape(array.shape)
    out[...] = array
    return out


def _decay(k2, tau):
    return np.exp(-DIFF * k2 * tau)


@pytest.mark.parametrize("modes", [(0, 1), (3, 0), (2, 5), (11, 15)])
def test_lateral_cosine_mode_exact_neumann(modes):
    # A single DCT-II cosine decays by exp(-D |k|^2 tau) with k_j = pi j / (n dx).
    ny, nx, tau = 12, 16, 0.7
    jy, jx = modes
    y = (np.arange(ny) + 0.5)[:, None]
    x = (np.arange(nx) + 0.5)[None, :]
    mode = np.cos(math.pi * jy * y / ny) * np.cos(math.pi * jx * x / nx)
    conc = 0.3 + mode[:, :, None, None] * np.ones((3, 2))
    out = lateral_diffuse(conc, DIFF, tau, boundary="neumann", spacing_um=DX)
    k2 = (math.pi * jy / (ny * DX)) ** 2 + (math.pi * jx / (nx * DX)) ** 2
    expect = 0.3 + mode[:, :, None, None] * _decay(k2, tau)
    assert np.max(np.abs(out - expect)) <= 1e-14


@pytest.mark.parametrize("modes", [(0, 1), (1, 3), (6, 0), (5, 8)])
def test_lateral_fourier_mode_exact_periodic(modes):
    ny, nx, tau = 12, 16, 1.3
    my, mx = modes
    y = np.arange(ny)[:, None]
    x = np.arange(nx)[None, :]
    mode = np.cos(2 * math.pi * (my * y / ny + mx * x / nx) + 0.4)
    conc = 1.0 + mode[:, :, None, None] * np.ones((3, 2))
    out = lateral_diffuse(conc, DIFF, tau, boundary="periodic", spacing_um=DX)
    k2 = (2 * math.pi * my / (ny * DX)) ** 2 + (2 * math.pi * mx / (nx * DX)) ** 2
    expect = 1.0 + mode[:, :, None, None] * _decay(k2, tau)
    assert np.max(np.abs(out - expect)) <= 1e-14


@pytest.mark.parametrize("boundary", ["neumann", "periodic"])
def test_lateral_mass_conserved(boundary):
    rng = np.random.default_rng(1)
    conc = rng.uniform(0.0, 0.2, size=(20, 28, 3, 2))
    out = lateral_diffuse(conc, DIFF, 5.0, boundary=boundary, spacing_um=DX)
    before = conc.sum(axis=(0, 1))
    after = out.sum(axis=(0, 1))
    assert np.max(np.abs(after - before) / before) <= 1e-13
    # Zero-D component (level 1, component 1) does not move.
    np.testing.assert_allclose(out[:, :, 1, 1], conc[:, :, 1, 1], rtol=0, atol=1e-15)
    # Diffusion smooths: the variance never grows.
    assert np.all(out.var(axis=(0, 1)) <= conc.var(axis=(0, 1)) + 1e-18)


def test_lateral_zero_time_is_identity():
    rng = np.random.default_rng(2)
    conc = rng.uniform(size=(8, 6, 3, 2))
    for boundary in ("neumann", "periodic"):
        out = lateral_diffuse(conc, DIFF, 0.0, boundary=boundary, spacing_um=DX)
        np.testing.assert_allclose(out, conc, rtol=0, atol=1e-15)


def test_lateral_long_time_reaches_mean():
    rng = np.random.default_rng(3)
    conc = rng.uniform(size=(10, 10, 3, 2))
    out = lateral_diffuse(conc, 1e3, 1e6, boundary="neumann", spacing_um=DX)
    np.testing.assert_allclose(out, np.broadcast_to(conc.mean(axis=(0, 1)), out.shape), rtol=1e-13)


def test_lateral_semigroup():
    # Exact in time: two half steps equal one full step.
    rng = np.random.default_rng(4)
    conc = rng.uniform(size=(9, 14, 3, 2))
    for boundary in ("neumann", "periodic"):
        one = lateral_diffuse(conc, DIFF, 2.0, boundary=boundary, spacing_um=DX)
        half = lateral_diffuse(conc, DIFF, 1.0, boundary=boundary, spacing_um=DX)
        two = lateral_diffuse(half, DIFF, 1.0, boundary=boundary, spacing_um=DX)
        np.testing.assert_allclose(two, one, rtol=0, atol=1e-14)


def test_lateral_pair_matches_periodic_strip():
    # Spec D1 3.5: the pair (+1, -1 phases of cos(kx)) follows the periodic strip of the same
    # sinusoid, read at the cells where cos = +1 and -1.
    n, tau, eps = 16, 2.0, 1e-3
    x = np.arange(n)
    k = 2 * math.pi * 3 / (n * DX)  # 3 periods across the strip
    strip = 0.1 + eps * np.cos(k * DX * x)[None, :, None, None] * np.ones((1, 1, 3, 2))
    out = lateral_diffuse(strip, DIFF, tau, boundary="periodic", spacing_um=DX)
    pair = np.stack([np.full((3, 2), 0.1 + eps), np.full((3, 2), 0.1 - eps)])
    got = lateral_diffuse(pair, DIFF, tau, boundary="pair", wavenumber_per_um=np.array([k]))
    amp = eps * np.exp(-DIFF * k * k * tau)
    np.testing.assert_allclose(got[0], 0.1 + amp, rtol=0, atol=1e-16)
    np.testing.assert_allclose(got[1], 0.1 - amp, rtol=0, atol=1e-16)
    # The strip's cell 0 sits at cos = +1.
    np.testing.assert_allclose(out[0, 0], got[0], rtol=0, atol=1e-15)


def test_lateral_pair_multiple_and_conserving():
    rng = np.random.default_rng(5)
    k = np.array([0.0, 2 * math.pi / 400.0, 2 * math.pi / 200.0])
    conc = rng.uniform(size=(6, 3, 2))
    out = lateral_diffuse(conc, DIFF, 3.0, boundary="pair", wavenumber_per_um=k)
    np.testing.assert_allclose(out[0::2] + out[1::2], conc[0::2] + conc[1::2], rtol=1e-15)
    np.testing.assert_allclose(out[:2], conc[:2], rtol=0, atol=2e-16)  # k = 0: no change
    diff_ratio = (out[4] - out[5]) / (conc[4] - conc[5])
    np.testing.assert_allclose(diff_ratio, np.exp(-DIFF * k[2] ** 2 * 3.0), rtol=1e-12)


def test_lateral_factors():
    f = lateral_factors(8, DX, "neumann", DIFF, 2.0)
    assert f.shape == (8, 3, 2)
    np.testing.assert_array_equal(f[0], 1.0)
    np.testing.assert_allclose(f[3], np.exp(-DIFF * (3 * math.pi / (8 * DX)) ** 2 * 2.0))
    half = lateral_factors(8, DX, "periodic", DIFF, 2.0, half=True)
    assert half.shape == (5, 3, 2)


# ---- TOTH (spec D1 3.5) ----------------------------------------------------------------------


def _toth_block(rng):
    conc = rng.uniform(0.05, 0.2, size=(10, 12, 3, 3))  # components: A, B, TOTH
    bound = np.zeros((10, 12, 3))
    bound[:, :, :2] = rng.uniform(0.0, 0.05, size=(10, 12, 2))  # gel levels only
    conc[..., 2] += bound
    diff = np.concatenate([DIFF, [[0.0], [0.0], [0.0]]], axis=1)  # TOTH has no own D
    return conc, bound, diff


@pytest.mark.parametrize("boundary", ["neumann", "periodic"])
def test_lateral_toth_moves_only_mobile_protons(boundary):
    rng = np.random.default_rng(6)
    conc, bound, diff = _toth_block(rng)
    d_toth = np.array([300.0, 500.0, 2000.0])
    toth = TothLateral(index=2, gel_bound=bound, diffusion_um2_s=d_toth)
    out = lateral_diffuse(conc, diff, 4.0, boundary=boundary, spacing_um=DX, toth=toth)
    # TOTH conserved per level.
    before = conc[..., 2].sum(axis=(0, 1))
    after = out[..., 2].sum(axis=(0, 1))
    assert np.max(np.abs(after - before) / before) <= 1e-13
    # Equals TOTH + (diffused TOTH_mob - TOTH_mob) with D_toth.
    mob = conc[..., 2] - bound
    mob_new = lateral_diffuse(mob, d_toth, 4.0, boundary=boundary, spacing_um=DX)
    np.testing.assert_allclose(out[..., 2], conc[..., 2] + mob_new - mob, rtol=0, atol=1e-15)
    # The other components are unaffected by the TOTH treatment.
    plain = lateral_diffuse(conc, diff, 4.0, boundary=boundary, spacing_um=DX)
    np.testing.assert_allclose(out[..., :2], plain[..., :2], rtol=0, atol=1e-16)


def test_lateral_toth_gel_bound_protons_stay():
    # Uniform mobile protons over non-uniform gelatin binding: TOTH does not change.
    rng = np.random.default_rng(7)
    _, bound, diff = _toth_block(rng)
    conc = np.full((10, 12, 3, 3), 0.1)
    conc[..., 2] = 0.02 + bound
    toth = TothLateral(index=2, gel_bound=bound, diffusion_um2_s=np.full(3, 1e3))
    out = lateral_diffuse(conc, diff, 50.0, boundary="neumann", spacing_um=DX, toth=toth)
    np.testing.assert_allclose(out[..., 2], conc[..., 2], rtol=0, atol=1e-16)
    # Without the treatment, the gel-bound pattern would diffuse away.
    plain = lateral_diffuse(conc, np.full((3, 3), 1e3), 50.0, boundary="neumann", spacing_um=DX)
    assert np.max(np.abs(plain[..., 2] - conc[..., 2])) > 1e-3


def test_lateral_toth_pair():
    conc = np.array([[[0.1, 0.05]], [[0.2, 0.09]]])  # 1 pair, 1 level, components (A, TOTH)
    bound = np.array([[0.01], [0.03]])
    toth = TothLateral(index=1, gel_bound=bound, diffusion_um2_s=np.array([800.0]))
    k = np.array([2 * math.pi / 100.0])
    out = lateral_diffuse(
        conc, np.array([[400.0, 0.0]]), 1.0, boundary="pair", wavenumber_per_um=k, toth=toth
    )
    mob = conc[:, :, 1] - bound
    f = math.exp(-800.0 * k[0] ** 2)
    mean, half = 0.5 * (mob[0] + mob[1]), 0.5 * (mob[0] - mob[1]) * f
    np.testing.assert_allclose(out[0, :, 1], bound[0] + mean + half, rtol=1e-14)
    np.testing.assert_allclose(out[1, :, 1], bound[1] + mean - half, rtol=1e-14)


def test_toth_lateral_diffusion_formula():
    rng = np.random.default_rng(8)
    beta = rng.uniform(0.001, 0.05, size=(4, 6, 7, 3))  # carriers: 2 families, H+, OH-
    d = np.array([[400.0, 500.0, 900.0], [80.0, 100.0, 700.0], [8130.0] * 3, [4600.0] * 3])
    got = toth_lateral_diffusion(beta, d)
    mean = beta.mean(axis=(1, 2))
    np.testing.assert_allclose(got, (d * mean).sum(axis=0) / mean.sum(axis=0), rtol=1e-14)
    inner = (slice(1, 5), slice(2, 6))
    got_inner = toth_lateral_diffusion(beta, d, inner=inner)
    mean_i = beta[:, 1:5, 2:6].mean(axis=(1, 2))
    np.testing.assert_allclose(got_inner, (d * mean_i).sum(axis=0) / mean_i.sum(axis=0), rtol=1e-14)
    # A single carrier gives its own D.
    np.testing.assert_allclose(toth_lateral_diffusion(beta[:1], d[:1]), d[0], rtol=1e-15)


def _ideal():
    q = {"value": 0.0, "unit": "1", "source": "assumed", "note": "test"}
    return ActivityParams(
        debye_huckel_A_20C=q,
        debye_huckel_A_25C=q,
        sit_b={"value": 1.5, "unit": "1", "source": "assumed", "note": "test"},
        delta_eps={},
    )


def test_carrier_capacities_add_up_to_buffer_capacity():
    sulfite = AcidBaseFamily(
        "sulfite", (0, -1, -2), (2, 1, 0), (1.86, 7.19), "van_t_hoff", (0.0, 0.0), (), "x"
    )
    hq = AcidBaseFamily(
        "hydroquinone", (0, -1, -2), (0, -1, -2), (9.85, 11.45), "van_t_hoff", (0.0, 0.0), (), "x"
    )
    gel = AcidBaseFamily("gelatin", (1, 0), (1, 0), (10.5,), "van_t_hoff", (0.0,), (), "none")
    spec = make_speciation([sulfite, hq, gel], 293.15, 0.0, 0.0, _ideal())
    rng = np.random.default_rng(9)
    ph = rng.uniform(7.5, 10.5, size=(4, 5))
    totals = np.stack(
        [rng.uniform(0.1, 0.8, (4, 5)), rng.uniform(0.0, 0.05, (4, 5)), np.full((4, 5), 0.02)]
    )
    beta = carrier_capacities(spec, totals, ph, [True, True, False])
    assert beta.shape == (4, 4, 5)
    gel_only = totals.copy()
    gel_only[:2] = 0.0
    beta_gel = buffer_capacity(spec, gel_only, ph) - beta[-2] - beta[-1]
    total = buffer_capacity(spec, totals, ph)
    np.testing.assert_allclose(beta.sum(axis=0) + beta_gel, total, rtol=1e-12)
    assert np.all(beta >= 0.0)


# ---- determinism (spec D1 3.6) ---------------------------------------------------------------


@pytest.mark.parametrize("boundary", BOUNDARIES)
def test_block_result_independent_of_alignment(boundary):
    rng = np.random.default_rng(10)
    if boundary == "pair":
        conc = rng.uniform(size=(8, 3, 3))
        bound = rng.uniform(0.0, 0.01, size=(8, 3))
        extra = {"wavenumber_per_um": np.array([0.01, 0.02, 0.05, 0.1])}
    else:
        conc = rng.uniform(size=(30, 22, 3, 3))
        bound = rng.uniform(0.0, 0.01, size=(30, 22, 3))
        extra = {"spacing_um": DX}
    diff = np.concatenate([DIFF, np.zeros((3, 1))], axis=1)

    def run(offset):
        toth = TothLateral(2, _shifted_copy(bound, offset), np.array([300.0, 400.0, 2e3]))
        return lateral_diffuse(
            _shifted_copy(conc, offset),
            _shifted_copy(diff, offset),
            1.7,
            boundary=boundary,
            toth=toth,
            **extra,
        )

    assert np.array_equal(run(0), run(8))


# ---- bad input -------------------------------------------------------------------------------


def test_lateral_rejects_bad_input():
    conc = np.ones((4, 5, 3, 2))
    with pytest.raises(ValueError, match="boundary"):
        lateral_diffuse(conc, DIFF, 1.0, boundary="dirichlet", spacing_um=DX)
    with pytest.raises(ValueError, match="tau"):
        lateral_diffuse(conc, DIFF, -1.0, boundary="neumann", spacing_um=DX)
    with pytest.raises(ValueError, match="spacing"):
        lateral_diffuse(conc, DIFF, 1.0, boundary="neumann")
    with pytest.raises(ValueError, match="spacing"):
        lateral_diffuse(conc, DIFF, 1.0, boundary="neumann", spacing_um=0.0)
    with pytest.raises(ValueError, match="broadcast"):
        lateral_diffuse(conc, np.ones((2, 2)), 1.0, boundary="neumann", spacing_um=DX)
    with pytest.raises(ValueError, match="diffusion"):
        lateral_diffuse(conc, -DIFF, 1.0, boundary="neumann", spacing_um=DX)
    bad = conc.copy()
    bad[0, 0, 0, 0] = np.inf
    with pytest.raises(ValueError, match="finite"):
        lateral_diffuse(bad, DIFF, 1.0, boundary="neumann", spacing_um=DX)
    with pytest.raises(ValueError, match="wavenumber"):
        lateral_diffuse(np.ones((4, 3, 2)), DIFF, 1.0, boundary="pair")
    with pytest.raises(ValueError, match="pair columns"):
        lateral_diffuse(
            np.ones((3, 3, 2)), DIFF, 1.0, boundary="pair", wavenumber_per_um=np.ones(2)
        )
    with pytest.raises(ValueError, match="only used"):
        lateral_diffuse(
            conc, DIFF, 1.0, boundary="neumann", spacing_um=DX, wavenumber_per_um=np.ones(1)
        )
    with pytest.raises(ValueError, match="lateral axes"):
        lateral_diffuse(np.ones(4), 1.0, 1.0, boundary="neumann", spacing_um=DX)


def test_toth_rejects_bad_input():
    conc = np.ones((4, 5, 3, 2))
    with pytest.raises(ValueError, match="index"):
        lateral_diffuse(
            conc, DIFF, 1.0, boundary="neumann", spacing_um=DX,
            toth=TothLateral(2, np.zeros((4, 5, 3)), np.ones(3)),
        )  # fmt: skip
    with pytest.raises(ValueError, match="gel_bound"):
        lateral_diffuse(
            conc, DIFF, 1.0, boundary="neumann", spacing_um=DX,
            toth=TothLateral(1, np.zeros((4, 5)), np.ones(3)),
        )  # fmt: skip
    with pytest.raises(ValueError, match="diffusion_um2_s"):
        lateral_diffuse(
            conc, DIFF, 1.0, boundary="neumann", spacing_um=DX,
            toth=TothLateral(1, np.zeros((4, 5, 3)), np.ones(2)),
        )  # fmt: skip
    with pytest.raises(TypeError):
        lateral_diffuse(conc, DIFF, 1.0, boundary="neumann", spacing_um=DX, toth=(1, 0, 0))
    with pytest.raises(ValueError, match="disagree"):
        toth_lateral_diffusion(np.ones((2, 3, 4)), np.ones((3, 4)))
    with pytest.raises(ValueError, match="zero"):
        toth_lateral_diffusion(np.zeros((2, 3, 4)), np.ones((2, 4)))
    with pytest.raises(ValueError, match="slices"):
        toth_lateral_diffusion(np.ones((2, 3, 3, 4)), np.ones((2, 4)), inner=(slice(0, 2),))
    with pytest.raises(ValueError, match="no cells"):
        toth_lateral_diffusion(
            np.ones((2, 3, 3, 4)), np.ones((2, 4)), inner=(slice(0, 0), slice(None))
        )
    with pytest.raises(ValueError):
        lateral_factors(0, DX, "neumann", 1.0, 1.0)
    with pytest.raises(ValueError, match="half"):
        lateral_factors(4, DX, "neumann", 1.0, 1.0, half=True)
