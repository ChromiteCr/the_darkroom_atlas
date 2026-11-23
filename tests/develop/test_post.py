"""Stop, fix, wash and dry of the bw-mechanistic model (spec D1 2.6, 2.0.3, 2.0.4)."""

import importlib.util
import itertools
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from darkroom_atlas.data.library import Library
from darkroom_atlas.data.models import Bath, Process, ProcessStep, StepKind
from darkroom_atlas.develop import post
from darkroom_atlas.develop.post import stop_update
from darkroom_atlas.imaging.silver import product_density, silver_density, transmittance_factors
from darkroom_atlas.pipeline import RunContext, StepContext, process
from darkroom_atlas.readout.products import visual_product
from darkroom_atlas.registry import get_step
from darkroom_atlas.spectral import DEFAULT_GRID
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
SHAPE = (11, 13)


@pytest.fixture(scope="module")
def data_lib():
    return Library(syn.DATA)


@pytest.fixture(scope="module")
def ctx(data_lib):
    return RunContext(library=data_lib)


@pytest.fixture(scope="module")
def film():
    return syn.material()


def _developed(film, shape=SHAPE, seed=1, **kw):
    rng = np.random.default_rng(seed)
    f = rng.uniform(0.05, 1.0, (*shape, K))
    x = f * rng.uniform(0.0, 0.9, (*shape, K))
    y = x * rng.uniform(0.3, 1.0, (*shape, K))
    kw.setdefault("growth_rate", rng.uniform(0.0, 0.05, shape))
    kw.setdefault("halide", rng.uniform(0.5, 3.0, shape))
    kw.setdefault("physical", rng.uniform(0.0, 0.1, shape))
    return syn.state(film, shape, 6.0, k=K, initiated=f, developed=x, developed_m2=y, **kw)


def _step(kind, bath, time):
    return ProcessStep(kind=kind, bath=bath, time=time)


def _run(handler_kind, state, step, ctx, data_lib, index=0):
    bath = None if step.bath is None else data_lib.bath(step.bath)
    handler = get_step("bw-mechanistic", handler_kind).fn
    return handler(state, step, StepContext(ctx, index, bath))


# ---- registration and data ---------------------------------------------------------------


@pytest.mark.parametrize("kind", ["stop", "fix", "wash", "dry"])
def test_registered_pointwise(kind):
    entry = get_step("bw-mechanistic", kind)
    assert entry.interaction_range_um == 0.0
    assert entry.fn is getattr(post, f"{kind}_bw")


@pytest.mark.xfail(
    strict=True,
    reason="darkroom_atlas.models does not import darkroom_atlas.develop.post yet; the owner "
    "of models/__init__.py must add the import (then remove this marker)",
)
def test_steps_registered_through_models_import():
    # A real run imports the models package, not post itself: the four steps must be
    # registered by that import alone (spec 2.0.4).
    code = (
        "import sys\n"
        "import darkroom_atlas.models, darkroom_atlas.pipeline, darkroom_atlas.runner\n"
        "from darkroom_atlas.registry import get_step\n"
        "for kind in ('stop', 'fix', 'wash', 'dry'):\n"
        "    assert get_step('bw-mechanistic', kind).interaction_range_um == 0.0, kind\n"
        "assert 'darkroom_atlas.develop.post' in sys.modules\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=False, timeout=120
    )
    assert result.returncode == 0, result.stderr


def test_bath_constants(data_lib):
    # Spec 2.6: acetic stop 2 s, water stop 15 s, fixer 20 s / 60 s, wash 120 s / 300 s.
    expected = {
        "stop-acetic": {"neutralisation_time": 2.0},
        "stop-water": {"neutralisation_time": 15.0},
        "fixer-thiosulfate": {"clearing_time_constant": 20.0, "dye_removal_time": 60.0},
        "water": {"wash_time_constant": 120.0, "wash_dye_time_constant": 300.0},
    }
    for bath_id, constants in expected.items():
        bath = data_lib.bath(bath_id)
        assert isinstance(bath, Bath)
        for name, value in constants.items():
            quantity = getattr(bath.constants, name)
            assert quantity.to("s") == value
            assert quantity.source.value == "assumed"


def test_bath_files_pass_check_data(data_lib):
    problems = [line for line in data_lib.validate() if "baths" in line]
    assert problems == []


# ---- stop ----------------------------------------------------------------------------------


def test_stop_y_closed_form():
    f, x0, y0 = 0.8, 0.3, 0.15
    kappa_end, tau_s = 0.7, 4.0

    def rhs(t, v):
        kappa = kappa_end * math.exp(-t / tau_s)
        return [kappa * (f - v[0]), 2.0 * kappa * (v[0] - v[1])]

    for time_s in (0.5, 2.0, 10.0, 30.0):
        sol = solve_ivp(rhs, (0.0, time_s), [x0, y0], method="DOP853", rtol=1e-13, atol=1e-15)
        x, y = stop_update(
            np.array(x0), np.array(y0), np.array(f), np.array(kappa_end), tau_s, time_s
        )
        assert abs(float(x) - sol.y[0, -1]) <= 1e-10
        assert abs(float(y) - sol.y[1, -1]) <= 1e-10


def test_stop_limits():
    rng = np.random.default_rng(0)
    f = rng.uniform(0.0, 1.0, 50)
    x0 = f * rng.uniform(0.0, 1.0, 50)
    y0 = x0 * rng.uniform(0.2, 1.0, 50)
    kappa = rng.uniform(0.01, 1.0, 50)
    # K = 0: unchanged.
    x, y = stop_update(x0, y0, f, kappa, 2.0, 0.0)
    np.testing.assert_allclose(x, x0, rtol=0, atol=1e-15)
    np.testing.assert_allclose(y, y0, rtol=0, atol=1e-15)
    x, y = stop_update(x0, y0, f, np.zeros(50), 2.0, 30.0)
    np.testing.assert_allclose(x, x0, rtol=0, atol=1e-15)
    # K -> infinity: X and Y tend to F (every initiated grain fully developed).
    x, y = stop_update(x0, y0, f, np.full(50, 1e3), 1.0, 100.0)
    np.testing.assert_allclose(x, f, rtol=0, atol=1e-14)
    np.testing.assert_allclose(y, f, rtol=0, atol=1e-14)
    # t -> infinity saturates at K = kappa tau_s.
    xa, _ = stop_update(x0, y0, f, kappa, 2.0, 1e4)
    k = kappa * 2.0
    np.testing.assert_allclose(xa, f - (f - x0) * np.exp(-k), rtol=1e-14)


def test_stop_update_bad_input():
    one = np.ones(3)
    with pytest.raises(ValueError, match="tau_s"):
        stop_update(one, one, one, one, 0.0, 1.0)
    with pytest.raises(ValueError, match="time_s"):
        stop_update(one, one, one, one, 1.0, -1.0)
    with pytest.raises(ValueError, match="kappa_end"):
        stop_update(one, one, one, -one, 1.0, 1.0)
    with pytest.raises(ValueError, match="finite"):
        stop_update(one * np.nan, one, one, one, 1.0, 1.0)


def test_stop_step_fields(film, ctx, data_lib):
    s = _developed(film)
    step = _step("stop", "stop-acetic", "30s")
    out = _run("stop", s, step, ctx, data_lib, index=1)
    assert not out.has("emulsion", "growth_rate")
    assert out.get("emulsion", "initiated") is s.get("emulsion", "initiated")
    assert out.get("emulsion", "halide") is s.get("emulsion", "halide")
    assert out.get("emulsion", "silver_physical") is s.get("emulsion", "silver_physical")
    const_cap = syn.class_capacity(film, K)
    layer = film.emulsion_layers[0]
    from darkroom_atlas.emulsion.classes import size_classes

    classes = size_classes(layer, K)
    kappa_end = (
        s.get("emulsion", "growth_rate")[..., None] * classes.surface_um2 / classes.volume_um3
    )
    x0 = s.get("emulsion", "silver") / const_cap
    x, y = stop_update(
        x0, s.get("emulsion", "silver_m2"), s.get("emulsion", "initiated"), kappa_end, 2.0, 30.0
    )
    np.testing.assert_allclose(out.get("emulsion", "silver"), x * const_cap, rtol=1e-14)
    np.testing.assert_allclose(out.get("emulsion", "silver_m2"), y, rtol=1e-14)
    assert not out.get("emulsion", "silver").flags.writeable
    # More silver than before wherever grains were initiated but unfinished.
    assert np.all(out.get("emulsion", "silver") >= s.get("emulsion", "silver") * (1 - 1e-15))
    (entry,) = out.meta["post"]
    assert dict(entry) == {
        "index": 1,
        "kind": "stop",
        "bath": "stop-acetic",
        "time_s": 30.0,
        "neutralisation_time_s": 2.0,
    }


def test_stop_extra_development_monotone(film, ctx, data_lib):
    s = _developed(film)
    step = _step("stop", "stop-acetic", "30s")
    totals = []
    for tau in (0.5, 2.0, 15.0, 60.0):
        bath = data_lib.bath("stop-acetic").model_copy(
            update={
                "constants": data_lib.bath("stop-acetic").constants.model_copy(
                    update={
                        "neutralisation_time": data_lib.bath(
                            "stop-acetic"
                        ).constants.neutralisation_time.model_copy(update={"value": tau})
                    }
                )
            }
        )
        out = post.stop_bw(s, step, StepContext(ctx, 0, bath))
        totals.append(out.get("emulsion", "silver"))
    for a, b in itertools.pairwise(totals):
        assert np.all(b >= a) and np.sum(b) > np.sum(a)
    # A water stop (15 s) develops more than an acetic stop (2 s).
    acetic = _run("stop", s, step, ctx, data_lib)
    water = _run("stop", s, _step("stop", "stop-water", "30s"), ctx, data_lib)
    assert np.all(silver_density(water) >= silver_density(acetic))
    assert np.sum(silver_density(water)) > np.sum(silver_density(acetic))


def test_stop_needs_bath_and_growth_rate(film, ctx, data_lib):
    s = _developed(film)
    with pytest.raises(ValueError, match="names no bath"):
        post.stop_bw(s, _step("stop", None, "30s"), StepContext(ctx, 0, None))
    with pytest.raises(ValueError, match="neutralisation_time"):
        _run("stop", s, _step("stop", "water", "30s"), ctx, data_lib)
    with pytest.raises(KeyError, match="growth_rate"):
        _run(
            "stop",
            s.drop("emulsion/growth_rate"),
            _step("stop", "stop-acetic", "30s"),
            ctx,
            data_lib,
        )
    bad = s.with_fields({"emulsion/growth_rate": -np.ones(SHAPE)})
    with pytest.raises(ValueError, match="negative"):
        _run("stop", bad, _step("stop", "stop-acetic", "30s"), ctx, data_lib)


# ---- fix, wash, dry ------------------------------------------------------------------------


def test_fix_closed_form(film, ctx, data_lib):
    s = _developed(film)
    out = _run("fix", s, _step("fix", "fixer-thiosulfate", "90s"), ctx, data_lib, index=2)
    # f_fix = 2 in the fixture film, tau_fix = 20 s, tau_dye = 60 s.
    np.testing.assert_allclose(
        out.get("emulsion", "halide"),
        s.get("emulsion", "halide") * math.exp(-90.0 / 40.0),
        rtol=1e-15,
    )
    dye = out.get("emulsion", "residual_dye")
    assert dye.strides == (0, 0)  # stays a broadcast view
    np.testing.assert_allclose(dye, math.exp(-90.0 / 60.0), rtol=1e-15)
    fixer = out.get("frame", "fixer_residual")
    assert fixer.strides == (0, 0)
    np.testing.assert_array_equal(fixer, 1.0)
    for name in ("silver", "silver_m2", "initiated", "growth_rate"):
        assert out.get("emulsion", name) is s.get("emulsion", name)
    (entry,) = out.meta["post"]
    assert dict(entry) == {
        "index": 2,
        "kind": "fix",
        "bath": "fixer-thiosulfate",
        "time_s": 90.0,
        "clearing_time_constant_s": 20.0,
        "dye_removal_time_s": 60.0,
        "fixing_factor": 2.0,
    }
    # f_fix = 1 clears twice as fast.
    trix_like = syn.material(development={"fixing_factor": 1.0})
    s1 = _developed(trix_like)
    out1 = _run("fix", s1, _step("fix", "fixer-thiosulfate", "90s"), ctx, data_lib)
    np.testing.assert_allclose(
        out1.get("emulsion", "halide"), s1.get("emulsion", "halide") * math.exp(-90.0 / 20.0)
    )


def test_wash_closed_form(film, ctx, data_lib):
    s = _developed(film)
    fixed = _run("fix", s, _step("fix", "fixer-thiosulfate", "5min"), ctx, data_lib)
    washed = _run("wash", fixed, _step("wash", "water", "10min"), ctx, data_lib, index=3)
    np.testing.assert_allclose(washed.get("frame", "fixer_residual"), math.exp(-600.0 / 120.0))
    np.testing.assert_allclose(
        washed.get("emulsion", "residual_dye"),
        math.exp(-300.0 / 60.0) * math.exp(-600.0 / 300.0),
        rtol=1e-15,
    )
    assert washed.get("frame", "fixer_residual").strides == (0, 0)
    assert washed.get("emulsion", "halide") is fixed.get("emulsion", "halide")
    assert [e["kind"] for e in washed.meta["post"]] == ["fix", "wash"]
    assert washed.meta["post"][1]["wash_time_constant_s"] == 120.0
    assert washed.meta["post"][1]["wash_dye_time_constant_s"] == 300.0
    # Without a fix step there is no fixer to wash out: the field stays absent.
    plain = _run("wash", s, _step("wash", "water", "1min"), ctx, data_lib)
    assert not plain.has("frame", "fixer_residual")
    # A non-uniform residual dye is scaled pixel by pixel.
    rng = np.random.default_rng(5)
    dye = rng.uniform(0.0, 1.0, SHAPE)
    s_dye = s.with_fields({"emulsion/residual_dye": dye})
    out = _run("wash", s_dye, _step("wash", "water", "1min"), ctx, data_lib)
    np.testing.assert_allclose(out.get("emulsion", "residual_dye"), dye * math.exp(-60.0 / 300.0))


def test_dry_records_only(film, ctx, data_lib):
    s = _developed(film)
    out = _run("dry", s, ProcessStep(kind="dry", time="10min"), ctx, data_lib, index=4)
    assert set(out.fields) == set(s.fields)
    for key in s.fields:
        assert out.fields[key] is s.fields[key]
    (entry,) = out.meta["post"]
    assert dict(entry) == {"index": 4, "kind": "dry", "bath": None, "time_s": 600.0}


def test_fix_and_wash_need_constants(film, ctx, data_lib):
    s = _developed(film)
    with pytest.raises(ValueError, match="clearing_time_constant"):
        _run("fix", s, _step("fix", "stop-acetic", "1min"), ctx, data_lib)
    with pytest.raises(ValueError, match="names no bath"):
        post.wash_bw(s, _step("wash", None, "1min"), StepContext(ctx, 0, None))
    with pytest.raises(ValueError, match="wash_time_constant"):
        _run("wash", s, _step("wash", "fixer-thiosulfate", "1min"), ctx, data_lib)
    with pytest.raises(KeyError, match="halide"):
        _run(
            "fix",
            s.drop("emulsion/halide"),
            _step("fix", "fixer-thiosulfate", "1min"),
            ctx,
            data_lib,
        )


def test_fix_needs_emulsion_layer(film, ctx, data_lib):
    layers = [
        layer.model_copy(update={"id": "other"}) if layer.id == "emulsion" else layer
        for layer in film.layers
    ]
    odd = film.model_copy(update={"layers": layers})
    s = _developed(film)
    bare = EmulsionState(odd, s.geometry, s.grid, s.seed, {})
    with pytest.raises(ValueError, match="no emulsion layer"):
        _run("fix", bare, _step("fix", "fixer-thiosulfate", "1min"), ctx, data_lib)


def test_skip_fix_leaves_halide_density(film, ctx, data_lib):
    s = _developed(film)
    stopped = _run("stop", s, _step("stop", "stop-acetic", "30s"), ctx, data_lib)
    unfixed = _run("wash", stopped, _step("wash", "water", "10min"), ctx, data_lib)
    fixed = _run(
        "wash",
        _run("fix", stopped, _step("fix", "fixer-thiosulfate", "5min"), ctx, data_lib),
        _step("wash", "water", "10min"),
        ctx,
        data_lib,
    )
    product = visual_product(DEFAULT_GRID)
    window = (slice(0, SHAPE[0]), slice(0, SHAPE[1]))

    def density(state):
        return product_density(
            10.0 ** -silver_density(state), transmittance_factors(state, *window), product
        )

    halide = unfixed.get("emulsion", "halide")
    cleared = unfixed.with_fields({"emulsion/halide": np.zeros(SHAPE)})
    # Skipping the fix adds CP_hal * halide (CP_hal = 0.3 in the fixture).
    np.testing.assert_allclose(density(unfixed) - density(cleared), 0.3 * halide, rtol=1e-12)
    assert np.all(density(unfixed) - density(fixed) > 0.3 * halide * 0.99)


# ---- tiling --------------------------------------------------------------------------------


def _process(steps):
    return Process.model_validate(
        {"schema_version": 1, "id": "test-post", "name": "post", "family": "bw", "steps": steps}
    )


@pytest.mark.parametrize("tile_px", [1, 3, 5, 64])
def test_post_steps_pointwise_tiling(film, data_lib, tile_px):
    rng = np.random.default_rng(9)
    s = _developed(film, residual_dye=rng.uniform(0.5, 1.0, SHAPE))
    proc = _process(
        [
            {"kind": "stop", "bath": "stop-acetic", "time": "30s"},
            {"kind": "fix", "bath": "fixer-thiosulfate", "time": "5min"},
            {"kind": "wash", "bath": "water", "time": "10min"},
            {"kind": "dry", "time": "10min"},
        ]
    )
    whole = process(s, proc, RunContext(library=data_lib))
    tiled = process(s, proc, RunContext(library=data_lib, tile_px=tile_px))
    assert set(whole.fields) == set(tiled.fields)
    assert "emulsion/growth_rate" not in whole.fields
    for key in whole.fields:
        assert np.array_equal(whole.fields[key], tiled.fields[key]), key
    assert whole.meta["post"] == tiled.meta["post"]
    assert [e["kind"] for e in whole.meta["post"]] == ["stop", "fix", "wash", "dry"]
    assert [r.kind for r in whole.history] == ["stop", "fix", "wash", "dry"]
    # The uniform fixer residual stays a zero-byte view through the tiled assembly.
    assert tiled.fields["frame/fixer_residual"].strides == (0, 0)


def test_stop_kind_enum():
    assert StepKind("stop") is StepKind.STOP
