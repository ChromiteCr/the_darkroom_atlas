import dataclasses
import functools
from pathlib import Path

import numpy as np
import pytest
from scipy.ndimage import gaussian_filter

from darkroom_atlas import pipeline, registry
from darkroom_atlas.data.library import Library
from darkroom_atlas.data.models import Process, ProcessStep, StepKind
from darkroom_atlas.fieldstore import FieldStore
from darkroom_atlas.light import Light, LowRankLight
from darkroom_atlas.pipeline import (
    Readout,
    RunContext,
    StepContext,
    expose,
    process,
    readout,
    segments,
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
    # Seed 0, the RunContext default: process and readout require state.seed == ctx.seed.
    return expose(ideal, light, TIME_S, RunContext(lib))


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


@pytest.mark.parametrize(
    ("first", "rest", "fragment"),
    [
        (np.float32, np.float64, "float32"),  # only the tile at the origin differs
        ((), (2,), "(2,)"),  # trailing dimensions differ
    ],
)
def test_tiled_rejects_inconsistent_tile_fields(
    scratch, exposed, identity, lib, first, rest, fragment
):
    # Each tile must return each field with the same dtype and trailing shape; stitching
    # used to restart the field from zeros, silently losing the tiles pasted before.
    @register_step("mixedtest", StepKind.identity, interaction_range_um=0.0)
    def mixed(state, step, ctx):
        signal = state.get("y", "exposure")
        kind = first if state.geometry.origin_px == (0, 0) else rest
        if isinstance(kind, tuple):
            signal = np.broadcast_to(
                signal[(...,) + (np.newaxis,) * len(kind)], (*signal.shape, *kind)
            )
        else:
            signal = signal.astype(kind)
        return state.with_fields({"y/signal": signal})

    state = _as_model(exposed, "mixedtest")
    with pytest.raises(ValueError, match="y/signal") as caught:
        process(state, identity, RunContext(lib, tile_px=16))
    assert fragment in str(caught.value)
    assert "mixedtest" in str(caught.value)


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


# ----------------------------------------------------------------------------------------
# D1: run context, segments, fused streaming, seeds and low-rank light
# ----------------------------------------------------------------------------------------


def _proc(*kinds):
    """A process of bath-free steps of ``kinds``, each lasting 10 s (identity: no time)."""
    steps = [
        ProcessStep(kind=kind, time=None if kind is StepKind.identity else 10.0) for kind in kinds
    ]
    return Process(schema_version=1, id="test-seg", name="segments", family="test", steps=steps)


def _stencil(field_in, field_out, radius, scale):
    """A step adding ``scale`` times the box sum over ``radius`` px (edge-padded), bit-exact."""

    def step(state, step, ctx):
        a = state.fields[field_in]
        padded = np.pad(a, radius, mode="edge")
        h, w = a.shape
        total = np.zeros_like(a)
        for dr in range(2 * radius + 1):
            for dc in range(2 * radius + 1):
                total = total + padded[dr : dr + h, dc : dc + w]
        return state.with_fields({field_out: a + scale * total})

    return step


def test_run_context_new_fields(lib, tmp_path):
    ctx = RunContext(lib)
    assert ctx.film_format is None
    assert ctx.profile == "image"
    assert ctx.workers == 1
    assert ctx.memory_budget_gb == 6.0
    assert ctx.scratch_dir is None

    ctx = RunContext(
        lib,
        film_format="135",
        profile="calibration",
        workers=3,
        memory_budget_gb=2,
        scratch_dir=str(tmp_path),
    )
    assert ctx.film_format == "135"
    assert ctx.profile == "calibration"
    assert ctx.workers == 3
    assert ctx.memory_budget_gb == 2.0 and isinstance(ctx.memory_budget_gb, float)
    assert ctx.scratch_dir == tmp_path and isinstance(ctx.scratch_dir, Path)
    assert RunContext(lib, profile="reference").profile == "reference"

    bad = [
        {"film_format": ""},
        {"film_format": 135},
        {"profile": "fast"},
        {"workers": 0},
        {"workers": 1.5},
        {"workers": True},
        {"memory_budget_gb": 0.0},
        {"memory_budget_gb": float("inf")},
        {"memory_budget_gb": True},
        {"scratch_dir": 3},
    ]
    for kwargs in bad:
        with pytest.raises((ValueError, TypeError), match=next(iter(kwargs))):
            RunContext(lib, **kwargs)


def test_segments(scratch):
    model = "segtest"
    register_step(model, StepKind.identity, interaction_range_um=5.0)(_stencil("x/a", "x/a", 1, 1))
    register_step(model, StepKind.presoak, interaction_range_um=0.0)(lambda s, st, c: s)
    register_step(model, StepKind.stop, interaction_range_um=None)(lambda s, st, c: s)
    register_step(model, StepKind.wash, interaction_range_um=0.0)(lambda s, st, c: s)
    register_step(model, StepKind.dry, interaction_range_um=2.0)(lambda s, st, c: s)
    I, P, S, W, D = (StepKind.identity, StepKind.presoak, StepKind.stop, StepKind.wash,
                     StepKind.dry)  # fmt: skip
    assert segments(_proc(I, P, S, W, D), model) == [(0, 1), (2,), (3, 4)]
    assert segments(_proc(S), model) == [(0,)]
    assert segments(_proc(S, S), model) == [(0,), (1,)]
    assert segments(_proc(I, P, W, D), model) == [(0, 1, 2, 3)]
    assert segments(_proc(S, I, S), model) == [(0,), (1,), (2,)]
    with pytest.raises(RegistryError, match="stabilize"):
        segments(_proc(I, StepKind.stabilize), model)
    with pytest.raises(TypeError):
        segments("identity", model)


def _fusion_model(scratch_model, calls=None):
    """Steps mixing a finite range, pointwise, global, rng and field-dropping handlers."""

    def noisy(state, step, ctx):
        geometry = state.geometry
        (r0, c0), (h, w) = geometry.origin_px, geometry.shape_hw
        noise = pixel_uniform(
            state.seed, "fusion", (r0, r0 + h), (c0, c0 + w), geometry.frame_shape_hw[1]
        )
        return state.with_fields({"x/a": state.fields["x/a"] * (1.0 + noise)})

    def drop(state, step, ctx):
        return state.drop("x/b").with_fields({"y/c": np.abs(state.fields["x/a"]) ** 0.5})

    def tank(state, step, ctx):
        # Whole-frame coupling: a global mean.
        return state.with_fields({"x/a": state.fields["x/a"] - 0.5 * state.fields["x/a"].mean()})

    pitch = 36000 / 96
    register_step(scratch_model, StepKind.identity, interaction_range_um=2 * pitch)(
        _stencil("x/a", "x/a", 2, 0.01)
    )
    register_step(scratch_model, StepKind.presoak, interaction_range_um=0.0)(noisy)
    register_step(scratch_model, StepKind.stop, interaction_range_um=None)(tank)
    register_step(scratch_model, StepKind.wash, interaction_range_um=0.5 * pitch)(
        _stencil("x/a", "x/b", 1, -0.02)
    )
    register_step(scratch_model, StepKind.dry, interaction_range_um=0.0)(drop)
    return _proc(StepKind.identity, StepKind.presoak, StepKind.stop, StepKind.wash, StepKind.dry)


def _fusion_state(exposed, model):
    rng = np.random.default_rng(7)
    base = dataclasses.replace(exposed, fields={}).with_fields(
        {"x/a": rng.random((64, 96)), "z/keep": rng.random((64, 96, 2))}
    )
    return _as_model(base, model)


def test_segment_fusion_equals_unfused(scratch, exposed, lib):
    proc = _fusion_model("fusiontest")
    assert segments(proc, "fusiontest") == [(0, 1), (2,), (3, 4)]
    state = _fusion_state(exposed, "fusiontest")

    untiled = process(state, proc, RunContext(lib, tile_px=None))
    for tile_px in (16, 23, 200):
        fused = process(state, proc, RunContext(lib, tile_px=tile_px))
        _assert_states_equal(fused, untiled)
        # Step by step: each step its own one-step process, tiled.
        stepwise = state
        for step in proc.steps:
            single = proc.model_copy(update={"steps": [step]})
            stepwise = process(stepwise, single, RunContext(lib, tile_px=tile_px))
        assert sorted(stepwise.fields) == sorted(fused.fields)
        for key in fused.fields:
            np.testing.assert_array_equal(fused.fields[key], stepwise.fields[key], err_msg=key)
            assert fused.fields[key].dtype == stepwise.fields[key].dtype

    assert sorted(untiled.fields) == ["x/a", "y/c", "z/keep"]
    assert [r.kind for r in untiled.history] == ["identity", "presoak", "stop", "wash", "dry"]
    assert [r.index for r in untiled.history] == [0, 1, 2, 3, 4]


def test_segment_halo_is_sum_of_ranges(scratch, exposed, identity, lib):
    pitch = exposed.geometry.pixel_pitch_um
    seen = []

    def first(state, step, ctx):
        seen.append(("first", state.geometry.origin_px, state.geometry.shape_hw))
        return _stencil("x/a", "x/a", 2, 0.1)(state, step, ctx)

    def second(state, step, ctx):
        seen.append(("second", state.geometry.origin_px, state.geometry.shape_hw))
        return _stencil("x/a", "x/a", 3, 0.1)(state, step, ctx)

    register_step("halotest", StepKind.identity, interaction_range_um=1.5 * pitch)(first)
    register_step("halotest", StepKind.wash, interaction_range_um=3 * pitch)(second)
    proc = _proc(StepKind.identity, StepKind.wash)
    state = _fusion_state(exposed, "halotest")
    untiled = process(state, proc, RunContext(lib))
    seen.clear()
    tiled = process(state, proc, RunContext(lib, tile_px=16))
    _assert_states_equal(tiled, untiled)
    # Halo ceil(1.5) + 3 = 5 px; both steps run on the same window, one tile after another.
    tiles = list(iter_tiles((64, 96), 16, 5))
    assert len(seen) == 2 * len(tiles)
    for i, tile in enumerate(tiles):
        rows, cols = tile.outer
        window = ((rows.start, cols.start), (rows.stop - rows.start, cols.stop - cols.start))
        assert seen[2 * i] == ("first", *window)
        assert seen[2 * i + 1] == ("second", *window)
    assert seen[2 * 7][2] == (26, 26)  # an interior tile: 16 + 2 * 5


def test_fused_step_sees_earlier_step_in_history(scratch, exposed, lib):
    histories = []

    @register_step("histtest", StepKind.wash, interaction_range_um=0.0)
    def record(state, step, ctx):
        histories.append(tuple(r.kind for r in state.history))
        return state

    register_step("histtest", StepKind.identity, interaction_range_um=0.0)(lambda s, st, c: s)
    proc = _proc(StepKind.identity, StepKind.wash)
    state = _fusion_state(exposed, "histtest")
    for tile_px in (None, 32):
        histories.clear()
        result = process(state, proc, RunContext(lib, tile_px=tile_px))
        assert histories and set(histories) == {("identity",)}
        assert [r.kind for r in result.history] == ["identity", "wash"]


def test_fused_segment_passes_untouched_fields_through(scratch, exposed, lib):
    @register_step("passtest", StepKind.wash, interaction_range_um=0.0)
    def double(state, step, ctx):
        return state.with_fields({"x/a": 2.0 * state.fields["x/a"]})

    state = _fusion_state(exposed, "passtest")
    result = process(state, _proc(StepKind.wash), RunContext(lib, tile_px=16))
    # The field no step touched is the input array itself, not a copy.
    assert result.fields["z/keep"] is state.fields["z/keep"]
    np.testing.assert_array_equal(result.fields["x/a"], 2.0 * state.fields["x/a"])
    assert not result.fields["x/a"].flags.writeable


def test_fused_segment_field_untouched_only_in_some_tiles(scratch, exposed, lib):
    @register_step("sometest", StepKind.wash, interaction_range_um=0.0)
    def partly(state, step, ctx):
        if state.geometry.origin_px == (0, 0):
            return state  # this tile leaves x/a alone
        return state.with_fields({"x/a": state.fields["x/a"] + 1.0})

    state = _fusion_state(exposed, "sometest")
    result = process(state, _proc(StepKind.wash), RunContext(lib, tile_px=16))
    expected = state.fields["x/a"] + 1.0
    expected[:16, :16] = state.fields["x/a"][:16, :16]
    np.testing.assert_array_equal(result.fields["x/a"], expected)


def test_fused_segment_keeps_uniform_fields_broadcast(scratch, exposed, lib):
    @register_step("unitest", StepKind.wash, interaction_range_um=0.0)
    def uniform(state, step, ctx):
        one = np.array([[True]])
        one.flags.writeable = False  # so the broadcast view is shared, not copied
        flag = np.broadcast_to(one, state.geometry.shape_hw)
        return state.with_fields({"frame/presoaked": flag})

    state = _fusion_state(exposed, "unitest")
    result = process(state, _proc(StepKind.wash), RunContext(lib, tile_px=16))
    presoaked = result.fields["frame/presoaked"]
    assert presoaked.shape == (64, 96) and presoaked.dtype == np.bool_
    assert presoaked.strides == (0, 0)
    assert presoaked.all()


def test_fused_segment_writes_to_scratch_dir(scratch, exposed, lib, tmp_path, monkeypatch):
    # Map every assembled field to a file (threshold 0) to see the scratch directory used.
    monkeypatch.setattr(
        pipeline, "FieldStore", functools.partial(FieldStore, ram_threshold_bytes=0)
    )
    proc = _fusion_model("scratchtest")
    state = _fusion_state(exposed, "scratchtest")
    reference = process(state, proc, RunContext(lib))
    result = process(state, proc, RunContext(lib, tile_px=16, scratch_dir=tmp_path))
    _assert_states_equal(result, reference)
    mapped = result.fields["y/c"]
    assert isinstance(mapped, np.memmap) and mapped.mode == "r"
    assert any(tmp_path.rglob("*.dat"))
    # Without a scratch directory the fields stay in memory.
    in_ram = process(state, proc, RunContext(lib, tile_px=16))
    assert not isinstance(in_ram.fields["y/c"], np.memmap)


def test_tiled_rejects_inconsistent_fields_in_fused_step(scratch, exposed, lib):
    register_step("mixfused", StepKind.identity, interaction_range_um=0.0)(lambda s, st, c: s)

    @register_step("mixfused", StepKind.wash, interaction_range_um=0.0)
    def mixed(state, step, ctx):
        if state.geometry.origin_px == (0, 0):
            return state.with_fields({"y/extra": state.fields["x/a"]})
        return state

    state = _fusion_state(exposed, "mixfused")
    with pytest.raises(ValueError, match="'wash' handler of model 'mixfused'"):
        process(state, _proc(StepKind.identity, StepKind.wash), RunContext(lib, tile_px=16))


def test_process_and_readout_reject_seed_mismatch(exposed, identity, lib):
    assert exposed.seed == 0
    with pytest.raises(ValueError, match="seed"):
        process(exposed, identity, RunContext(lib, seed=1))
    state = process(exposed, identity, RunContext(lib))
    with pytest.raises(ValueError, match="seed"):
        readout(state, "ideal_colorimetric", RunContext(lib, seed=3))
    seeded = expose(exposed.material, scene_light(read_image(SAMPLE).linear_rgb, time_s=TIME_S),
                    TIME_S, RunContext(lib, seed=3))  # fmt: skip
    ctx = RunContext(lib, seed=3)
    assert readout(process(seeded, identity, ctx), "ideal_colorimetric", ctx).data.shape[:2] == (
        64,
        96,
    )


def test_expose_accepts_low_rank_light(scratch, ideal, lib):
    h, w, b = 6, 9, DEFAULT_GRID.n_bands
    rng = np.random.default_rng(5)
    basis = 1e6 * rng.random((3, b))
    low = LowRankLight(rng.random((h, w, 3)), basis, DEFAULT_GRID, {"source": "test"})
    weights = rng.random((b, 1))

    @register_exposure("lowranktest")
    def project_only(material, light, time_s, ctx):
        signal = light.project(weights)[..., 0] * time_s
        return unexposed_state(material, light, ctx).with_fields({"x/exposure": signal})

    material = ideal.model_copy(update={"model": "lowranktest"})
    ctx = RunContext(lib)
    state = expose(material, low, 0.5, ctx)
    assert state.geometry.shape_hw == (h, w)
    assert dict(state.meta["light"]) == {"source": "test"}
    full = expose(material, low.to_light(), 0.5, ctx)
    np.testing.assert_allclose(
        state.get("x", "exposure"), full.get("x", "exposure"), rtol=1e-12, atol=0.0
    )
    with pytest.raises(TypeError, match="light"):
        expose(material, np.zeros((h, w, b)), 0.5, ctx)
