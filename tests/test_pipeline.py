import dataclasses
from pathlib import Path

import numpy as np
import pytest
from scipy.ndimage import gaussian_filter

from darkroom_atlas import registry
from darkroom_atlas.data.library import Library
from darkroom_atlas.data.models import StepKind
from darkroom_atlas.light import Light
from darkroom_atlas.pipeline import (
    Readout,
    RunContext,
    StepContext,
    expose,
    process,
    readout,
    unexposed_state,
)
from darkroom_atlas.registry import RegistryError, register_exposure, register_step
from darkroom_atlas.rng import pixel_uniform
from darkroom_atlas.scene.camera import scene_light
from darkroom_atlas.scene.io import read_image
from darkroom_atlas.spectral import DEFAULT_GRID, SpectralGrid, illuminant
from darkroom_atlas.state import EmulsionState, FrameGeometry
from darkroom_atlas.tiling import iter_tiles

TESTS = Path(__file__).resolve().parent
SAMPLE = TESTS / "assets" / "sample.png"
FIXTURES = TESTS / "fixtures" / "library"
TIME_S = 1 / 125


@pytest.fixture
def scratch():
    """Registrations made in a test are undone afterwards, even if the test fails."""
    with registry.temporary_registrations():
        yield


@pytest.fixture(scope="module")
def lib():
    return Library()


@pytest.fixture(scope="module")
def fx():
    return Library(FIXTURES)


@pytest.fixture(scope="module")
def ideal(lib):
    return lib.material("ideal")


@pytest.fixture(scope="module")
def identity(lib):
    return lib.process("identity")


@pytest.fixture(scope="module")
def light():
    return scene_light(read_image(SAMPLE).linear_rgb, time_s=TIME_S)


@pytest.fixture(scope="module")
def exposed(ideal, light, lib):
    return expose(ideal, light, TIME_S, RunContext(lib, seed=11))


def _as_model(state, model):
    """The same state, but for a copy of its material that runs under ``model``."""
    return dataclasses.replace(state, material=state.material.model_copy(update={"model": model}))


def _assert_states_equal(a, b):
    assert sorted(a.fields) == sorted(b.fields)
    for key in a.fields:
        np.testing.assert_array_equal(a.fields[key], b.fields[key], err_msg=key)
        assert a.fields[key].dtype == b.fields[key].dtype
    assert a.geometry == b.geometry
    assert a.history == b.history


# ----------------------------------------------------------------------------------------
# expose
# ----------------------------------------------------------------------------------------


def test_expose_sets_geometry_seed_meta(ideal, light, lib):
    ctx = RunContext(lib, seed=5)
    state = expose(ideal, light, TIME_S, ctx)
    assert state.geometry.pixel_pitch_um == 36000 / 96
    assert state.geometry.shape_hw == (64, 96)
    assert state.geometry.origin_px == (0, 0)
    assert state.geometry.frame_shape_hw == (64, 96)
    assert state.seed == ctx.seed == 5
    assert state.grid == ctx.grid
    assert state.material is ideal
    assert state.history == ()
    assert state.meta["light"]["iso"] == 100.0
    assert dict(state.meta["light"]) == dict(light.meta)

    medium = expose(ideal, light, TIME_S, RunContext(lib, long_edge_mm=60.0))
    assert medium.geometry.pixel_pitch_um == pytest.approx(60000 / 96, rel=1e-15)


def test_expose_rejects_grid_mismatch(ideal, lib):
    coarse = SpectralGrid(380.0, 780.0, 10.0)
    light = scene_light(read_image(SAMPLE).linear_rgb, grid=coarse, time_s=TIME_S)
    with pytest.raises(ValueError, match="grid"):
        expose(ideal, light, TIME_S, RunContext(lib))
    # A run on the light's own grid works.
    state = expose(ideal, light, TIME_S, RunContext(lib, grid=coarse))
    assert state.grid == coarse


def test_expose_rejects_nonpositive_time(ideal, light, lib):
    ctx = RunContext(lib)
    for bad in (0.0, -1.0, float("nan"), float("inf")):
        with pytest.raises(ValueError, match="time_s"):
            expose(ideal, light, bad, ctx)


def test_expose_unknown_model(exposed, light, lib):
    material = exposed.material.model_copy(update={"model": "nope"})
    with pytest.raises(RegistryError, match="ideal"):
        expose(material, light, TIME_S, RunContext(lib))


def test_expose_checks_model_result(scratch, ideal, light, lib):
    @register_exposure("wronggeom")
    def wrong_geometry(material, light, time_s, ctx):
        geometry = FrameGeometry.from_long_edge(light.shape_hw, 24.0)
        return EmulsionState(material, geometry, ctx.grid, ctx.seed)

    @register_exposure("wrongseed")
    def wrong_seed(material, light, time_s, ctx):
        return dataclasses.replace(unexposed_state(material, light, ctx), seed=ctx.seed + 1)

    @register_exposure("notastate")
    def not_a_state(material, light, time_s, ctx):
        return {}

    ctx = RunContext(lib)
    for model in ("wronggeom", "wrongseed"):
        material = ideal.model_copy(update={"model": model})
        with pytest.raises(ValueError, match=model):
            expose(material, light, TIME_S, ctx)
    with pytest.raises(TypeError, match="notastate"):
        expose(ideal.model_copy(update={"model": "notastate"}), light, TIME_S, ctx)


# ----------------------------------------------------------------------------------------
# process
# ----------------------------------------------------------------------------------------


def test_process_records_history(exposed, identity, lib):
    state = process(exposed, identity, RunContext(lib))
    assert len(state.history) == 1
    record = state.history[0]
    assert record.kind == "identity"
    assert record.index == 0
    assert record.bath is None
    assert record.time_s is None
    for layer in ("x", "y", "z"):
        np.testing.assert_array_equal(state.get(layer, "signal"), state.get(layer, "exposure"))
    # The input state is left alone.
    assert exposed.history == ()
    assert not exposed.has("x", "signal")


def test_process_unknown_step_for_model(exposed, fx):
    with pytest.raises(RegistryError, match="develop") as info:
        process(exposed, fx.process("test-bw"), RunContext(fx))
    assert "ideal" in str(info.value)


def test_process_rejects_grid_mismatch(exposed, identity, lib):
    with pytest.raises(ValueError, match="grid"):
        process(exposed, identity, RunContext(lib, grid=SpectralGrid(380.0, 780.0, 10.0)))


def test_tiled_equals_untiled_pointwise(exposed, identity, lib):
    untiled = process(exposed, identity, RunContext(lib, tile_px=None))
    tiled = process(exposed, identity, RunContext(lib, tile_px=17))
    _assert_states_equal(tiled, untiled)


def test_tiled_equals_untiled_with_range(scratch, exposed, identity, lib):
    pitch = exposed.geometry.pixel_pitch_um

    def blur(state, step, ctx):
        blurred = gaussian_filter(
            state.get("x", "exposure"), sigma=2, truncate=3, mode="nearest"
        )  # radius 6 px
        return state.with_fields({"x/signal": blurred})

    register_step("tiletest", StepKind.identity, interaction_range_um=6 * pitch)(blur)
    register_step("shorttest", StepKind.identity, interaction_range_um=1 * pitch)(blur)

    noise = np.random.default_rng(3).random(exposed.geometry.shape_hw)
    base = exposed.with_fields({"x/exposure": noise})
    state = _as_model(base, "tiletest")
    untiled = process(state, identity, RunContext(lib, tile_px=None))
    tiled = process(state, identity, RunContext(lib, tile_px=16))
    difference = np.max(np.abs(tiled.get("x", "signal") - untiled.get("x", "signal")))
    assert difference < 1e-12

    # The overlap really comes from interaction_range_um: too short a range leaves seams.
    short = _as_model(base, "shorttest")
    seamed = process(short, identity, RunContext(lib, tile_px=16))
    reference = process(short, identity, RunContext(lib, tile_px=None))
    assert np.max(np.abs(seamed.get("x", "signal") - reference.get("x", "signal"))) > 1e-6


def test_tiled_rng_consistent(scratch, exposed, identity, lib):
    calls = []

    @register_step("rngtest", StepKind.identity, interaction_range_um=0.0)
    def noisy(state, step, ctx):
        geometry = state.geometry
        (r0, c0), (h, w) = geometry.origin_px, geometry.shape_hw
        calls.append(geometry.shape_hw)
        noise = pixel_uniform(
            state.seed, "test-noise", (r0, r0 + h), (c0, c0 + w), geometry.frame_shape_hw[1], (2,)
        )
        return state.with_fields({"frame/noise": noise})

    state = _as_model(exposed, "rngtest")
    untiled = process(state, identity, RunContext(lib, tile_px=None))
    assert calls == [(64, 96)]
    calls.clear()
    tiled = process(state, identity, RunContext(lib, tile_px=16))
    assert len(calls) == len(list(iter_tiles((64, 96), 16, 0))) == 24
    assert all(shape == (16, 16) for shape in calls)
    np.testing.assert_array_equal(tiled.get("frame", "noise"), untiled.get("frame", "noise"))
    assert np.ptp(untiled.get("frame", "noise")) > 0.5
    _assert_states_equal(tiled, untiled)


def test_global_step_not_tiled(scratch, exposed, identity, lib):
    received = []

    @register_step("globaltest", StepKind.identity, interaction_range_um=None)
    def whole_tank(state, step, ctx):
        received.append(state.geometry)
        return state.with_fields({"x/signal": state.get("x", "exposure")})

    state = _as_model(exposed, "globaltest")
    result = process(state, identity, RunContext(lib, tile_px=16))
    assert len(received) == 1
    assert received[0] == exposed.geometry
    assert received[0].shape_hw == (64, 96)
    assert received[0].origin_px == (0, 0)
    np.testing.assert_array_equal(result.get("x", "signal"), exposed.get("x", "exposure"))


def test_tiled_step_may_drop_fields(scratch, exposed, identity, lib):
    @register_step("droptest", StepKind.identity, interaction_range_um=0.0)
    def replace_all(state, step, ctx):
        signal = 2.0 * state.get("y", "exposure")
        return dataclasses.replace(state, fields={"y/signal": signal}).with_meta(dropped=True)

    state = _as_model(exposed, "droptest")
    untiled = process(state, identity, RunContext(lib, tile_px=None))
    tiled = process(state, identity, RunContext(lib, tile_px=16))
    assert sorted(untiled.fields) == ["y/signal"]
    _assert_states_equal(tiled, untiled)
    assert tiled.meta["dropped"] is True
    assert "light" in tiled.meta


def test_bath_resolved_into_context(scratch, exposed, fx):
    seen = []

    def record(state, step, ctx):
        assert isinstance(ctx, StepContext)
        seen.append((ctx.step_index, step.kind, None if ctx.bath is None else ctx.bath.id))
        assert ctx.run is run
        return state

    for kind in (StepKind.develop, StepKind.stop, StepKind.fix, StepKind.wash):
        register_step("bathtest", kind)(record)

    run = RunContext(fx, tile_px=None)
    result = process(_as_model(exposed, "bathtest"), fx.process("test-bw"), run)
    assert seen == [
        (0, StepKind.develop, "test-bw-dev"),
        (1, StepKind.stop, None),
        (2, StepKind.fix, "test-fix"),
        (3, StepKind.wash, None),
    ]
    assert [r.kind for r in result.history] == ["develop", "stop", "fix", "wash"]
    assert [r.index for r in result.history] == [0, 1, 2, 3]
    assert [r.bath for r in result.history] == ["test-bw-dev", None, "test-fix", None]
    assert [r.time_s for r in result.history] == [585.0, 30.0, 300.0, 600.0]


def test_run_context_validated(lib):
    with pytest.raises(ValueError):
        RunContext(lib, tile_px=0)
    with pytest.raises(ValueError):
        RunContext(lib, seed=-1)
    with pytest.raises(ValueError):
        RunContext(lib, long_edge_mm=0.0)
    with pytest.raises(TypeError):
        RunContext("data")


# ----------------------------------------------------------------------------------------
# readout
# ----------------------------------------------------------------------------------------


def test_readout_unknown_device(exposed, identity, lib):
    ctx = RunContext(lib)
    state = process(exposed, identity, ctx)
    with pytest.raises(RegistryError, match="ideal_colorimetric"):
        readout(state, "nope", ctx)


def test_readout_as_next_light(ideal, lib):
    ctx = RunContext(lib)
    h, w, b = 4, 5, DEFAULT_GRID.n_bands
    flux = np.broadcast_to(1e6 * illuminant("D50", DEFAULT_GRID), (h, w, b)).copy()
    uniform = Light(flux, DEFAULT_GRID)

    amount = np.linspace(0.0, 1.0, h * w).reshape(h, w)
    negative = Readout(
        kind="spectral_transmittance",
        data=np.repeat(amount[..., None], b, axis=-1),
        grid=DEFAULT_GRID,
    )
    printed = expose(ideal, uniform.filtered(negative.data), 1.0, ctx)
    reference = expose(ideal, uniform, 1.0, ctx)
    for layer in ("x", "y", "z"):
        np.testing.assert_allclose(
            printed.get(layer, "exposure"),
            amount * reference.get(layer, "exposure"),
            rtol=1e-12,
            atol=0.0,
        )


def test_readout_validated():
    with pytest.raises(ValueError):
        Readout(kind="nope", data=np.zeros((2, 2, 3)))
    with pytest.raises(ValueError):
        Readout(kind="spectral_transmittance", data=np.zeros((2, 2, 81)))  # no grid
    with pytest.raises(ValueError):
        Readout(kind="spectral_transmittance", data=np.zeros((2, 2, 80)), grid=DEFAULT_GRID)
    with pytest.raises(ValueError):
        Readout(kind="image_linear_srgb", data=np.zeros((2, 2, 4)))
    image = Readout(kind="image_linear_srgb", data=np.zeros((2, 2, 3)), meta={"a": 1})
    assert not image.data.flags.writeable
    with pytest.raises(TypeError):
        image.meta["a"] = 2
