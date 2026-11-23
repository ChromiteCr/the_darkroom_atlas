"""Development kinetics: resolved constants, surface state and hazards (spec D1 2.4.2-2.4.6,
2.4.13)."""

import dataclasses
import math

import numpy as np
import pytest
from scipy.optimize import brentq

from darkroom_atlas.chemistry.kinetics import (
    ARRHENIUS_REFERENCE_K,
    POSITIVE_PART_EPS,
    arrhenius_factor,
    cluster_hazards,
    positive_part,
    resolve_kinetics,
    silver_saturation,
    surface_state,
)
from darkroom_atlas.chemistry.solution import working_solution
from darkroom_atlas.chemistry.species import GAS_CONSTANT
from darkroom_atlas.data.library import Library
from darkroom_atlas.data.provenance import Quantity


@pytest.fixture(scope="module")
def library():
    return Library()


@pytest.fixture(scope="module")
def species(library):
    return library.species()


@pytest.fixture(scope="module")
def kfile(library):
    return library.kinetics()


def _solution(library, species, kfile, bath_id, dilution="stock", t=20.0):
    return working_solution(
        library.bath(bath_id), species, kfile, dilution=dilution, temperature_c=t
    )


@pytest.fixture(scope="module")
def d76(library, species, kfile):
    return _solution(library, species, kfile, "d-76-published")


@pytest.fixture(scope="module")
def kin_d76(kfile, d76, species):
    return resolve_kinetics(kfile, d76, species)


def _with(kfile, section, name, value, unit):
    """Copy of the kinetics file with one constant replaced."""
    block = getattr(kfile, section)
    q = Quantity(value=value, unit=unit, source="assumed")
    return kfile.model_copy(update={section: block.model_copy(update={name: q})})


def _without_agent_field(kfile, agent, name):
    entry = kfile.agents[agent].model_copy(update={name: None})
    return kfile.model_copy(update={"agents": {**kfile.agents, agent: entry}})


def _col(*values):
    return np.array(values, dtype=np.float64)


def _metol_hq(kin_d76, **changes):
    """The D-76 kinetics with chosen per-agent constants (order: hydroquinone, metol)."""
    assert kin_d76.agent_ids == ("hydroquinone", "metol")
    return dataclasses.replace(kin_d76, **{k: np.asarray(v, float) for k, v in changes.items()})


# ---------------------------------------------------------------------------------------
# resolve_kinetics
# ---------------------------------------------------------------------------------------


def test_resolve_d76_roles_and_reference_values(kin_d76, kfile):
    kin = kin_d76
    assert kin.agent_ids == ("hydroquinone", "metol")
    assert kin.roles == ("regenerator", "eta")
    np.testing.assert_array_equal(kin.charges, [-2.0, -1.0])  # HQ2-, metol phenolate
    np.testing.assert_array_equal(kin.sulfonated, [True, True])
    # 20 C is T_r: no Arrhenius scaling.
    hq, met = kfile.agents["hydroquinone"], kfile.agents["metol"]
    np.testing.assert_array_equal(kin.k_et, [hq.k_et.value, met.k_et.value])
    np.testing.assert_array_equal(kin.K_ads, [hq.K_ads.value, met.K_ads.value])
    np.testing.assert_array_equal(kin.lambda0, [hq.lambda0.value, met.lambda0.value])
    assert math.isnan(kin.k_des[0]) and kin.k_des[1] == met.k_des.value
    assert kin.k_reg[0] == hq.k_reg.value and math.isnan(kin.k_reg[1])
    assert kin.K_br == kfile.halides.K_Br.value
    assert kin.K_i == pytest.approx(30.0 * kin.K_br, rel=1e-15)
    assert kin.n_sat == 5 and isinstance(kin.n_sat, int)
    assert kin.solvent is kfile.solvent and kin.transport is kfile.transport
    assert kin.gelatin is kfile.gelatin
    for name in ("charges", "k_et", "K_ads", "lambda0", "k_des", "k_reg", "sulfonated"):
        assert not getattr(kin, name).flags.writeable


def test_resolve_xtol_and_rodinal_roles(library, species, kfile):
    xtol = resolve_kinetics(kfile, _solution(library, species, kfile, "xtol"), species)
    assert xtol.agent_ids == ("dimezone-s", "ascorbate")
    assert xtol.roles == ("eta", "regenerator")
    np.testing.assert_array_equal(xtol.charges, [-1.0, -1.0])
    np.testing.assert_array_equal(xtol.sulfonated, [False, False])

    rodinal = resolve_kinetics(
        kfile, _solution(library, species, kfile, "rodinal", dilution="1+25"), species
    )
    assert rodinal.agent_ids == ("p-aminophenol",)
    assert rodinal.roles == ("single",)
    # A single agent is an ETA without regenerator (spec D1 2.4.3): it needs k_des.
    assert rodinal.k_des[0] == kfile.agents["p-aminophenol"].k_des.value
    assert math.isnan(rodinal.k_reg[0])
    np.testing.assert_array_equal(rodinal.eta_mask, [True])


@pytest.mark.parametrize(
    ("bath", "dilution", "agents", "roles", "sulfonated"),
    [
        ("d-76", "stock", ("hydroquinone", "metol"), ("regenerator", "eta"), (True, True)),
        ("d-76-published", "1+1", ("hydroquinone", "metol"), ("regenerator", "eta"), (1, 1)),
        ("xtol", "1+1", ("dimezone-s", "ascorbate"), ("eta", "regenerator"), (0, 0)),
        (
            "hc-110-pre2019",
            "1+31",
            ("hydroquinone", "methylphenidone"),
            ("regenerator", "eta"),
            (1, 0),
        ),
        ("hc-110-2019", "1+31", ("hydroquinone", "dimezone-s"), ("regenerator", "eta"), (1, 0)),
        ("t-max-dev", "1+4", ("hydroquinone", "dimezone-s"), ("regenerator", "eta"), (1, 0)),
        ("rodinal", "1+50", ("p-aminophenol",), ("single",), (1,)),
    ],
)
def test_roles_of_every_developer(
    library, species, kfile, bath, dilution, agents, roles, sulfonated
):
    # Spec D1 2.3.5: roles and sulfonated products of every D1 developer.
    kin = resolve_kinetics(kfile, _solution(library, species, kfile, bath, dilution), species)
    assert kin.agent_ids == agents
    assert kin.roles == roles
    np.testing.assert_array_equal(kin.sulfonated, np.array(sulfonated, dtype=bool))
    # Pyrazolidones share one kinetic class (spec D1 2.3.5).
    for i, cid in enumerate(agents):
        if cid in ("dimezone-s", "methylphenidone"):
            assert kin.k_et[i] == kfile.agents["pyrazolidone"].k_et.value


def test_agents_absent_from_the_bath_are_left_out(kin_d76):
    assert "dimezone-s" not in kin_d76.agent_ids and "ascorbate" not in kin_d76.agent_ids


def test_arrhenius_factor():
    assert arrhenius_factor(67.5, ARRHENIUS_REFERENCE_K) == 1.0
    t = 297.15
    expected = math.exp(-(67.5e3 / GAS_CONSTANT) * (1.0 / t - 1.0 / 293.15))
    assert arrhenius_factor(67.5, t) == pytest.approx(expected, rel=1e-15)
    assert arrhenius_factor(67.5, t) > 1.0 > arrhenius_factor(67.5, 291.15)
    with pytest.raises(ValueError, match="kelvin"):
        arrhenius_factor(60.0, 0.0)
    with pytest.raises(ValueError, match="finite"):
        arrhenius_factor(math.nan, 293.15)


def test_resolved_rates_scale_with_each_agents_arrhenius_factor(kfile, d76, species):
    warm = resolve_kinetics(kfile, dataclasses.replace(d76, temperature_c=24.0), species)
    ref = resolve_kinetics(kfile, d76, species)
    for i, agent in enumerate(("hydroquinone", "metol")):
        a = arrhenius_factor(kfile.agents[agent].activation_energy.value, 297.15)
        assert warm.k_et[i] == pytest.approx(ref.k_et[i] * a, rel=1e-14)
        assert warm.lambda0[i] == pytest.approx(ref.lambda0[i] * a, rel=1e-14)
        assert warm.K_ads[i] == ref.K_ads[i]  # adsorption constants do not scale
    a_met = arrhenius_factor(61.8, 297.15)
    a_hq = arrhenius_factor(67.5, 297.15)
    assert warm.k_des[1] == pytest.approx(ref.k_des[1] * a_met, rel=1e-14)
    assert warm.k_reg[0] == pytest.approx(ref.k_reg[0] * a_hq, rel=1e-14)
    assert warm.u_gel == ref.u_gel and warm.u_halide == ref.u_halide


def test_halide_adsorption_follows_van_t_hoff(kfile, d76, species):
    def k_br(t_c):
        return resolve_kinetics(kfile, dataclasses.replace(d76, temperature_c=t_c), species)

    cold, warm = k_br(18.3), k_br(23.9)
    # K(T) = K_20 exp(-(dH/R)(1/T - 1/T_r)), dH = -45 kJ/mol: adsorption weakens when warm.
    expected = math.exp(-(-45e3 / GAS_CONSTANT) * (1.0 / 291.45 - 1.0 / 297.05))
    assert cold.K_br / warm.K_br == pytest.approx(expected, rel=1e-13)
    assert cold.K_br / warm.K_br == pytest.approx(1.4, abs=0.05)  # spec D1 2.4.11
    assert warm.K_i == pytest.approx(30.0 * warm.K_br, rel=1e-15)  # same factor for K_I


def test_u_gel_folds_the_reference_ionic_strength(kfile, d76, species):
    kin = resolve_kinetics(
        _with(kfile, "barriers", "ionic_strength_ref", 4.0, "mol/L"), d76, species
    )
    assert kin.u_gel == pytest.approx(2.0 * kfile.barriers.u_gel.value, rel=1e-15)


def test_resolve_rejects_missing_constants(kfile, d76, species):
    agents = {k: v for k, v in kfile.agents.items() if k != "metol"}
    with pytest.raises(ValueError, match="kinetic class 'metol' is not in kinetics.agents"):
        resolve_kinetics(kfile.model_copy(update={"agents": agents}), d76, species)
    with pytest.raises(ValueError, match=r"needs agents.metol.k_des"):
        resolve_kinetics(_without_agent_field(kfile, "metol", "k_des"), d76, species)
    with pytest.raises(ValueError, match=r"needs agents.hydroquinone.k_reg"):
        resolve_kinetics(_without_agent_field(kfile, "hydroquinone", "k_reg"), d76, species)
    with pytest.raises(ValueError, match="n_sat must be an integer"):
        resolve_kinetics(_with(kfile, "latent", "n_sat", 4.5, "1"), d76, species)
    with pytest.raises(ValueError, match="not in species"):
        resolve_kinetics(kfile, d76, {k: v for k, v in species.items() if k != "sulfite"})


def test_kinetics_validates_its_arrays(kin_d76):
    with pytest.raises(ValueError, match="shape"):
        dataclasses.replace(kin_d76, k_et=np.array([1.0]))
    with pytest.raises(ValueError, match="k_des must be finite and > 0"):
        dataclasses.replace(kin_d76, k_des=np.array([np.nan, 0.0]))
    with pytest.raises(ValueError, match="unknown agent role"):
        dataclasses.replace(kin_d76, roles=("regenerator", "catalyst"))
    with pytest.raises(ValueError, match="n_sat"):
        dataclasses.replace(kin_d76, n_sat=2.5)


# ---------------------------------------------------------------------------------------
# surface_state
# ---------------------------------------------------------------------------------------


def test_no_bromide_limit(kin_d76):
    kin = kin_d76
    a = np.array([[2e-3, 1e-4], [5e-4, 3e-3]])
    zero = np.zeros(2)
    s = surface_state(a, zero, zero, kin, 2.4)
    np.testing.assert_array_equal(s.theta_halide, 0.0)
    supply = kin.k_reg[0] * a[0]
    phi = kin.k_et[1] / (kin.k_des[1] + supply)
    xi = 1.0 + kin.K_ads[0] * a[0] + kin.K_ads[1] * a[1] * (1.0 + phi)
    np.testing.assert_allclose(s.theta[0], kin.K_ads[0] * a[0] / xi, rtol=1e-15)
    np.testing.assert_allclose(s.theta[1], kin.K_ads[1] * a[1] / xi, rtol=1e-15)
    np.testing.assert_allclose(s.regenerated[1], supply / (kin.k_des[1] + supply), rtol=1e-15)
    assert np.all(np.isnan(s.regenerated[0]))
    np.testing.assert_allclose(
        s.kappa_tilde, kin.k_et[0] * s.theta[0] + kin.k_et[1] * s.theta[1], rtol=1e-15
    )
    barrier = kin.u_gel / math.sqrt(2.4)
    lam = kin.lambda0[0] * s.theta[0] * np.exp(-2.0 * barrier) + kin.lambda0[1] * s.theta[
        1
    ] * np.exp(-1.0 * barrier)
    np.testing.assert_allclose(s.lambda_tilde, lam, rtol=1e-15)


def test_bromide_free_start_through_positive_part(kin_d76):
    # A bromide-free bath gives C+ = eps/2, a negligible but finite coverage (spec D1 2.4.6).
    a = np.array([[2e-3], [5e-4]])
    br = positive_part(np.array([0.0]))
    assert br[0] == 0.5 * POSITIVE_PART_EPS
    io = positive_part(np.array([-1e-14]))
    assert 0.0 < io[0] < 0.5 * POSITIVE_PART_EPS
    s = surface_state(a, br, io, kin_d76, 2.4)
    # theta_X <= (K_Br + K_I) eps/2: negligible but finite.
    assert 0.0 < s.theta_halide[0] < (kin_d76.K_br + kin_d76.K_i) * 0.5 * POSITIVE_PART_EPS


def test_saturating_bromide(kin_d76):
    a = np.array([[2e-3], [5e-4]])
    zero = np.zeros(1)
    free = surface_state(a, zero, zero, kin_d76, 2.4)
    sat = surface_state(a, np.array([1e5]), zero, kin_d76, 2.4)
    assert sat.theta_halide[0] > 1.0 - 1e-6
    assert sat.kappa_tilde[0] < 1e-6 * free.kappa_tilde[0]
    assert sat.lambda_tilde[0] < 1e-6 * free.lambda_tilde[0]


def _x_at(kappa, lam, h, p, t):
    """Developed fraction X(t) of one class with constant kappa and lambda (spec D1 2.4.2,
    2.4.4): dX/dt = kappa (F - X), F = sum_n p_n (1 - exp(-h_n lambda t))."""
    x = 0.0
    for hn, pn in zip(h, p, strict=True):
        b = hn * lam
        if pn == 0.0:
            continue
        if b == 0.0:
            continue
        g = (
            1.0
            - math.exp(-kappa * t)
            - kappa / (kappa - b) * (math.exp(-b * t) - math.exp(-kappa * t))
        )
        x += pn * g
    return x


P_LATENT = np.array([0.1, 0.0, 0.1, 0.2, 0.2, 0.4])  # n = 0 .. 5


def test_bromide_restrains(kin_d76):
    # Rising initial Br-: kappa~, lambda~ and X at a fixed time all fall (spec D1 2.4.13).
    a = np.array([[2e-3], [5e-4]])
    zero = np.zeros(1)
    h = cluster_hazards(kin_d76, 1e-3, 5, 0.7)
    kappas, lams, xs = [], [], []
    for br in (0.0, 1e-4, 1e-3, 1e-2, 1e-1):
        s = surface_state(a, np.array([br]), zero, kin_d76, 2.4)
        kappa, lam = s.kappa_tilde[0], s.lambda_tilde[0]
        kappas.append(kappa)
        lams.append(lam)  # Lambda(t) = lambda~ t at constant chemistry
        xs.append(_x_at(kappa, lam, h, P_LATENT, 300.0))
    assert np.all(np.diff(kappas) < 0.0)
    assert np.all(np.diff(lams) < 0.0)
    assert np.all(np.diff(xs) < 0.0)


def test_lh_fractional_order(kin_d76):
    # Hydroquinone alone, no bromide: d ln kappa~ / d ln C = 1 - theta (spec D1 2.4.3).
    kin = kin_d76
    zero = np.zeros(1)
    for theta_target in (0.4, 0.05, 0.9):
        c = theta_target / (1.0 - theta_target) / kin.K_ads[0]  # K a / (1 + K a) = theta
        h = 1e-4

        def ln_kappa(ln_c):
            a = np.array([[math.exp(ln_c)], [0.0]])
            return math.log(surface_state(a, zero, zero, kin, 2.4).kappa_tilde[0])

        order = (ln_kappa(math.log(c) + h) - ln_kappa(math.log(c) - h)) / (2.0 * h)
        theta = surface_state(np.array([[c], [0.0]]), zero, zero, kin, 2.4).theta[0, 0]
        assert theta == pytest.approx(theta_target, rel=1e-12)
        assert abs(order - (1.0 - theta)) <= 1e-6
    # theta = 0.4 gives James's order 0.6.
    c = 0.4 / 0.6 / kin.K_ads[0]
    up = surface_state(np.array([[c * (1 + 1e-6)], [0.0]]), zero, zero, kin, 2.4).kappa_tilde
    down = surface_state(np.array([[c * (1 - 1e-6)], [0.0]]), zero, zero, kin, 2.4).kappa_tilde
    order = (math.log(up[0]) - math.log(down[0])) / (math.log1p(1e-6) - math.log1p(-1e-6))
    assert order == pytest.approx(0.6, abs=1e-6)


def _kappa_alone_and_together(kin, a_r, a_e):
    zero = np.zeros(1)
    both = surface_state(np.array([[a_r], [a_e]]), zero, zero, kin, 2.4).kappa_tilde[0]
    r = surface_state(np.array([[a_r], [0.0]]), zero, zero, kin, 2.4).kappa_tilde[0]
    e = surface_state(np.array([[0.0], [a_e]]), zero, zero, kin, 2.4).kappa_tilde[0]
    return both, r, e


def test_superadditivity_emerges(kin_d76):
    # Order (hydroquinone, metol). A strongly adsorbing ETA with slow desorption blocks its
    # own sites; a fast regenerator frees them (spec D1 2.4.3).
    kin = _metol_hq(
        kin_d76,
        k_et=[0.005, 0.5],
        K_ads=[10.0, 1e3],
        k_des=[np.nan, 0.01],
        k_reg=[1e4, np.nan],
    )
    both, hq, metol = _kappa_alone_and_together(kin, 1e-3, 1e-3)
    assert both > hq + metol
    assert both / (hq + metol) > 2.0
    # Without regeneration the combination is not superadditive.
    kin0 = dataclasses.replace(kin, k_reg=np.array([0.0, np.nan]))
    both0, hq0, metol0 = _kappa_alone_and_together(kin0, 1e-3, 1e-3)
    assert both0 <= hq0 + metol0
    # k_reg = 0 and Phi_E -> 0: equal to the sum. The agents still compete for sites in Xi,
    # so equality holds only at negligible coverage (K a = 1e-13 here; spec issue).
    kin_lim = dataclasses.replace(
        kin0, K_ads=np.array([1e-10, 1e-10]), k_des=np.array([np.nan, 1e15])
    )
    both, hq, metol = _kappa_alone_and_together(kin_lim, 1e-3, 1e-3)
    assert abs(both - (hq + metol)) <= 1e-12 * (hq + metol)


def test_coverages_bounded(kin_d76):
    rng = np.random.default_rng(20261121)
    n = 4000
    a = 10.0 ** rng.uniform(-12.0, 1.0, size=(2, n))
    a[:, :10] = 0.0
    br = 10.0 ** rng.uniform(-12.0, 1.0, size=n)
    io = 10.0 ** rng.uniform(-12.0, -1.0, size=n)
    kin = _metol_hq(kin_d76, K_ads=[10.0, 3e3], k_reg=[1e3, np.nan], k_des=[np.nan, 0.01])
    s = surface_state(a, br, io, kin, 0.5)
    assert np.all((s.theta >= 0.0) & (s.theta <= 1.0))
    assert np.all((s.theta_halide >= 0.0) & (s.theta_halide <= 1.0))
    assert np.all(s.theta.sum(axis=0) + s.theta_halide < 1.0)
    assert np.all(s.kappa_tilde >= 0.0) and np.all(s.lambda_tilde >= 0.0)
    g = s.regenerated[1]
    assert np.all((g >= 0.0) & (g <= 1.0))
    positive = s.kappa_tilde > 0.0
    np.testing.assert_allclose(s.partition[:, positive].sum(axis=0), 1.0, rtol=1e-14)
    np.testing.assert_array_equal(s.partition[:, ~positive], 0.0)
    assert np.any(~positive)  # the all-zero columns
    for name in (
        "theta",
        "theta_halide",
        "kappa_tilde",
        "lambda_tilde",
        "partition",
        "regenerated",
    ):
        assert not getattr(s, name).flags.writeable


def test_vectorised_equals_scalar(kin_d76):
    rng = np.random.default_rng(7)
    n = 257
    a = 10.0 ** rng.uniform(-8.0, -1.0, size=(2, n))
    br = 10.0 ** rng.uniform(-9.0, -1.0, size=n)
    io = 10.0 ** rng.uniform(-9.0, -3.0, size=n)
    batch = surface_state(a, br, io, kin_d76, 2.4)
    fields = ("theta", "theta_halide", "kappa_tilde", "lambda_tilde", "partition", "regenerated")
    for j in range(n):
        one = surface_state(a[:, j : j + 1], br[j : j + 1], io[j : j + 1], kin_d76, 2.4)
        for name in fields:
            np.testing.assert_array_equal(getattr(one, name)[..., 0], getattr(batch, name)[..., j])
    # Same result for inputs shifted by 8 bytes in memory (spec D1 3.6).
    buf = np.empty(2 * n + 1)
    shifted = buf[1:].reshape(2, n)
    shifted[...] = a
    buf_br = np.empty(n + 1)
    br_shift = buf_br[1:]
    br_shift[...] = br
    moved = surface_state(shifted, br_shift, io.copy(), kin_d76, 2.4)
    for name in fields:
        np.testing.assert_array_equal(getattr(moved, name), getattr(batch, name))


@pytest.mark.parametrize("agent", [0, 1])
def test_arrhenius_scaling(kfile, d76, species, agent):
    # Spec D1 2.4.13 conditions at the surface level: one agent, speciation independent of
    # temperature (the same active concentration at both temperatures), Delta H_Br = 0, no
    # bromide. The time to X = 0.5 scales as exp(Ea/R Delta(1/T)) within 1e-3.
    kf = _with(kfile, "halides", "adsorption_enthalpy", 0.0, "kJ/mol")
    a = np.zeros((2, 1))
    a[agent, 0] = 1e-3
    zero = np.zeros(1)
    times, rates = [], []
    for t_c in (20.0, 24.0):
        kin = resolve_kinetics(kf, dataclasses.replace(d76, temperature_c=t_c), species)
        s = surface_state(a, zero, zero, kin, d76.ionic_strength)
        kappa, lam = s.kappa_tilde[0], s.lambda_tilde[0]
        h = cluster_hazards(kin, 1e-3, 5, 0.7)
        rates.append((kappa, lam))

        def residual(t, kappa=kappa, lam=lam, h=h):
            return _x_at(kappa, lam, h, P_LATENT, t) - 0.5

        times.append(brentq(residual, 1e-6, 1e9, xtol=1e-14, rtol=1e-15))
    name = ("hydroquinone", "metol")[agent]
    ea = kfile.agents[name].activation_energy.value * 1e3
    expected = math.exp(ea / GAS_CONSTANT * (1.0 / 293.15 - 1.0 / 297.15))
    assert times[0] / times[1] == pytest.approx(expected, rel=1e-3)
    assert rates[1][0] / rates[0][0] == pytest.approx(expected, rel=1e-12)
    assert rates[1][1] / rates[0][1] == pytest.approx(expected, rel=1e-12)


def test_bromide_raises_the_apparent_activation_energy(kfile, d76, species):
    # Exothermic halide adsorption (Delta H_Br < 0) frees sites as T rises, so with Br-
    # the apparent Ea exceeds the agent's own (Sugai & Furuichi, spec D1 2.4.5).
    a = np.array([[0.0], [1e-3]])  # metol only

    def apparent_ea(br):
        k = []
        for t_c in (18.0, 24.0):
            kin = resolve_kinetics(kfile, dataclasses.replace(d76, temperature_c=t_c), species)
            k.append(surface_state(a, np.array([br]), np.zeros(1), kin, 2.4).kappa_tilde[0])
        return GAS_CONSTANT * math.log(k[1] / k[0]) / (1.0 / 291.15 - 1.0 / 297.15) / 1e3

    assert apparent_ea(0.0) == pytest.approx(61.8, abs=1e-6)
    assert apparent_ea(1e-2) > 61.8 + 1.0


def test_surface_state_without_agents(kin_d76):
    empty = dataclasses.replace(
        kin_d76,
        agent_ids=(),
        roles=(),
        charges=np.zeros(0),
        k_et=np.zeros(0),
        K_ads=np.zeros(0),
        lambda0=np.zeros(0),
        k_des=np.zeros(0),
        k_reg=np.zeros(0),
        sulfonated=np.zeros(0, dtype=bool),
    )
    s = surface_state(np.zeros((0, 3)), np.full(3, 1e-3), np.zeros(3), empty, 1.0)
    np.testing.assert_array_equal(s.kappa_tilde, 0.0)
    np.testing.assert_array_equal(s.lambda_tilde, 0.0)
    assert s.theta.shape == (0, 3) and s.partition.shape == (0, 3)


def test_surface_state_rejects_bad_input(kin_d76):
    a = np.full((2, 3), 1e-3)
    ok = np.zeros(3)
    with pytest.raises(ValueError, match=r"active must have shape \(2, N\)"):
        surface_state(np.zeros(3), ok, ok, kin_d76, 1.0)
    with pytest.raises(ValueError, match=r"active must have shape \(2, N\)"):
        surface_state(np.zeros((3, 3)), ok, ok, kin_d76, 1.0)
    with pytest.raises(ValueError, match=r"bromide must have shape \(3,\)"):
        surface_state(a, np.zeros(2), ok, kin_d76, 1.0)
    with pytest.raises(ValueError, match=r"iodide must have shape \(3,\)"):
        surface_state(a, ok, np.zeros((3, 1)), kin_d76, 1.0)
    with pytest.raises(ValueError, match="active must be >= 0"):
        surface_state(-a, ok, ok, kin_d76, 1.0)
    with pytest.raises(ValueError, match="bromide must be >= 0"):
        surface_state(a, _col(0.0, -1e-15, 0.0), ok, kin_d76, 1.0)
    with pytest.raises(ValueError, match="iodide must be finite"):
        surface_state(a, ok, _col(0.0, np.nan, 0.0), kin_d76, 1.0)
    for bad in (0.0, -1.0, math.inf, math.nan):
        with pytest.raises(ValueError, match="ionic_strength"):
            surface_state(a, ok, ok, kin_d76, bad)


# ---------------------------------------------------------------------------------------
# cluster_hazards, positive part, C_sat
# ---------------------------------------------------------------------------------------


def test_cluster_hazards_values(kin_d76, kfile):
    u_li = kfile.latent.u_li.value
    zeta = kfile.solvent.latent_attrition.value
    h = cluster_hazards(kin_d76, 2e-3, 5, 0.7)
    assert h.shape == (6,) and not h.flags.writeable
    assert h[0] == 2e-3
    for n in (1, 2, 3, 4):
        assert h[n] == pytest.approx(math.exp(-u_li * (5 - n)) * math.exp(-zeta * 0.7), rel=1e-15)
    assert h[5] == 1.0
    longer = cluster_hazards(kin_d76, 0.0, 8, 0.0)
    np.testing.assert_array_equal(longer[5:], 1.0)
    assert longer[0] == 0.0


@pytest.mark.parametrize("seed", range(5))
def test_hazards_monotone_in_cluster_size(kin_d76, seed):
    rng = np.random.default_rng(seed)
    kin = dataclasses.replace(
        kin_d76, u_li=float(rng.uniform(0.0, 5.0)), n_sat=int(rng.integers(2, 8))
    )
    cap = int(rng.integers(2, 10))
    h = cluster_hazards(kin, float(rng.uniform(0.0, 0.1)), cap, float(rng.uniform(0.0, 1.0)))
    # n >= 1; h_0 = h_fog is a separate fog factor and may exceed h_2.
    assert np.all(np.diff(h[1:]) >= 0.0)
    assert np.all((h[1:] > 0.0) & (h[1:] <= 1.0))


def test_cluster_cap_below_n_sat(kin_d76):
    h = cluster_hazards(kin_d76, 0.0, 3, 0.0)
    assert h.shape == (4,)
    assert h[3] == pytest.approx(math.exp(-kin_d76.u_li * 2.0), rel=1e-15)


def test_sulfite_attrition_slows_only_small_centres(kin_d76):
    low = cluster_hazards(kin_d76, 1e-3, 5, 0.05)
    high = cluster_hazards(kin_d76, 1e-3, 5, 0.8)
    assert np.all(high[1:5] < low[1:5])
    assert high[5] == low[5] == 1.0 and high[0] == low[0]


def test_cluster_hazards_rejects_bad_input(kin_d76):
    for cap in (1, 0, True, 2.0, "5"):
        with pytest.raises(ValueError, match="cluster_cap"):
            cluster_hazards(kin_d76, 0.0, cap, 0.1)
    for fog in (-1e-3, math.nan, math.inf):
        with pytest.raises(ValueError, match="fog_ratio"):
            cluster_hazards(kin_d76, fog, 5, 0.1)
    for so3 in (-1e-6, math.nan):
        with pytest.raises(ValueError, match="free_sulfite_bulk"):
            cluster_hazards(kin_d76, 0.0, 5, so3)


def test_positive_part():
    c = np.array([-1.0, -1e-9, 0.0, 1e-9, 1e-6, 1.0])
    p = positive_part(c)
    assert np.all(p > 0.0)
    assert p[2] == 0.5 * POSITIVE_PART_EPS
    big = c[4:]
    assert np.all(p[4:] - big <= POSITIVE_PART_EPS**2 / (4.0 * big) * (1 + 1e-6))
    assert p[0] < 1e-17  # a negative value maps to nearly 0
    eps = POSITIVE_PART_EPS
    direct = 0.5 * (c[1] + math.sqrt(c[1] ** 2 + eps**2))  # no cancellation at this scale
    assert p[1] == pytest.approx(direct, rel=1e-12)
    assert positive_part(np.array([-1e3]))[0] > 0.0
    np.testing.assert_array_equal(positive_part(5e-10), positive_part(np.array(5e-10)))


def test_silver_saturation(kin_d76, kfile):
    sol = kfile.solvent
    c_ref = sol.c_sat_ref.value
    nu_s, nu_b = sol.sulfite_exponent.value, sol.bromide_exponent.value
    so3 = np.array([0.1, 0.4, 0.4, 0.4])
    br = np.array([1e-3, 1e-3, 0.0, -1e-12])
    c = silver_saturation(so3, br, kin_d76.solvent)

    def expected(s, b):
        # Spec D1 2.4.6 with C+ of each concentration.
        sp, bp = positive_part(np.array([s]))[0], positive_part(np.array([b]))[0]
        return c_ref * (sp / 0.1) ** nu_s * (1e-3 / math.sqrt(bp * bp + 1e-10)) ** nu_b

    for j in range(4):
        assert c[j] == pytest.approx(expected(so3[j], br[j]), rel=1e-13)
    assert c[1] / c[0] == pytest.approx(4.0**nu_s, rel=1e-12)
    # Bromide-free: capped by the 1e-5 M floor, finite.
    assert c[2] == pytest.approx(c_ref * 4.0**nu_s * 100.0**nu_b, rel=1e-6)
    assert np.isfinite(c[3]) and c[3] == pytest.approx(c[2], rel=1e-6)
