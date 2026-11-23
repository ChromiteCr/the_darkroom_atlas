"""Development and presoak steps of the bw-mechanistic model and the strip path (spec D1 2.0.3,
2.4.10, 3.6-3.8).

Frames are small and development short. The canonical block of the spec (256 cells plus a
halo of about 44 cells) costs about a minute per pass, so the default tests run
:func:`develop_frame` (or the registered step with ``halo_um`` patched) on small blocks;
only ``test_halo_convergence`` (slow) uses the halo of spec D1 3.6.
"""

import dataclasses
import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import pytest

from darkroom_atlas.data.library import Library
from darkroom_atlas.data.models import Agitation, Process, ProcessStep, StepKind
from darkroom_atlas.data.provenance import Quantity
from darkroom_atlas.develop import step as dstep
from darkroom_atlas.fieldstore import FieldStore
from darkroom_atlas.pipeline import RunContext, StepContext, process
from darkroom_atlas.registry import get_step
from darkroom_atlas.spectral import DEFAULT_GRID
from darkroom_atlas.state import EmulsionState, FrameGeometry


def _load_helpers():
    # The helpers live in this directory's conftest.py; load it by path under its own name.
    name = "d1c_develop_conftest"
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name("conftest.py"))
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


_h = _load_helpers()
EXPOSURE_TIME_S = _h.EXPOSURE_TIME_S

SMALL_TANK = Agitation(method="inversion", initial_s=15.0, interval_s=30.0, duration_s=5.0)
ROLL_CM2 = _h.ROLL_AREA_UM2 * 1e-8  # a 135-36 roll in cm2
DEV_TIME = "60s"
N_Z = 2  # image profile
OUT_KEYS = {
    "emulsion/silver",
    "emulsion/initiated",
    "emulsion/silver_m2",
    "emulsion/silver_physical",
    "emulsion/halide",
    "emulsion/growth_rate",
    "emulsion/residual_dye",
}
DEVELOP_META = {
    "bath",
    "dilution",
    "temperature_c",
    "time_s",
    "agitation",
    "tank",
    "profile",
    "chem_cell_um",
    "halo_um",
    "block_cells",
    "n_steps",
    "min_scaled_concentration",
    "prepass",
    "bulk_initial",
    "bulk_final",
    "zero_exposure",
}


@pytest.fixture(scope="module")
def lib():
    return Library()


@pytest.fixture(scope="module")
def film():
    return Library(_h.FIXTURES).material(_h.FILM_ID)


@pytest.fixture(autouse=True)
def test_rates(monkeypatch):
    """The engine tests' kinetics (conftest TEST_RATES): the packaged placeholders hardly
    develop. Patched in the main process; workers receive the resolved kinetics."""
    real = dstep.resolve_kinetics
    monkeypatch.setattr(dstep, "resolve_kinetics", lambda *a: _h.with_test_rates(real(*a)))


def develop_step(time=DEV_TIME, agitation=SMALL_TANK, tank=None, **kw):
    if tank is None:
        tank = {"volume_ml": 300.0, "film_area_cm2": ROLL_CM2}
    return ProcessStep(
        kind=StepKind.DEVELOP,
        bath="d-76",
        dilution=kw.pop("dilution", "stock"),
        time=time,
        temperature=Quantity(value=20.0, unit="C", source="assumed"),
        agitation=agitation,
        tank=tank,
        **kw,
    )


def make_state(material, ev, es, pitch_um, profile="image"):
    meta = {
        "time_s": EXPOSURE_TIME_S,
        "temperature_c": 23.0,
        "n_z": ev.shape[2],
        "profile": profile,
        "film_format": "135",
    }
    one = np.ones((1, 1))
    one.flags.writeable = False
    fields = {
        "emulsion/exposure_volume": ev,
        "emulsion/exposure_surface": es,
        "emulsion/residual_dye": np.broadcast_to(one, ev.shape[:2]),
    }
    geometry = FrameGeometry(ev.shape[:2], pitch_um)
    return EmulsionState(
        material, geometry, DEFAULT_GRID, 0, fields=fields, meta={"exposure": meta}
    )


def scene(shape=(40, 44)):
    """Exposure ramp (log-spaced e_V over x) with an unexposed patch, e_S small."""
    h, w = shape
    x = np.linspace(-1.0, 2.5, w)
    ev = np.repeat(np.repeat((10.0**x)[None, :, None], h, axis=0), N_Z, axis=2)
    ev[h // 4 : h // 2, w // 3 : w // 2] = 0.0
    ev[:, :, 1] *= 0.7
    return ev, 0.01 * ev


def silver_total(state):
    return state.fields["emulsion/silver"].sum(axis=2) + state.fields["emulsion/silver_physical"]


def _ctx(lib, **kw):
    return RunContext(lib, **kw)


def _step_ctx(lib, ctx=None):
    ctx = ctx or _ctx(lib)
    return StepContext(ctx, 0, lib.bath("d-76"))


SMALL_HALO_UM = 3 * 24.0  # 3 cells


def _process(state, ctx, step):
    proc = Process(schema_version=1, id="test/develop", name="t", family="bw", steps=[step])
    return process(state, proc, ctx)


@pytest.fixture(scope="module")
def tiling_reference(lib, film):
    with pytest.MonkeyPatch.context() as mp:
        real = dstep.resolve_kinetics
        mp.setattr(dstep, "resolve_kinetics", lambda *a: _h.with_test_rates(real(*a)))
        mp.setattr(dstep, "BLOCK_CELLS", 8)
        mp.setattr(dstep, "halo_um", lambda *a: SMALL_HALO_UM)
        state = make_state(film, *scene(), 12.0)
        return _process(state, _ctx(lib), develop_step())


def test_registered_steps():
    dev = get_step(dstep.MODEL, StepKind.DEVELOP)
    pre = get_step(dstep.MODEL, StepKind.PRESOAK)
    assert dev.fn is dstep.develop_bw and dev.interaction_range_um is None
    assert pre.fn is dstep.presoak_bw and pre.interaction_range_um == 0.0


def test_develop_fields_and_meta(tiling_reference):
    out = tiling_reference
    assert set(out.fields) == OUT_KEYS
    k = out.fields["emulsion/silver"].shape[2]
    for key in ("emulsion/silver", "emulsion/initiated", "emulsion/silver_m2"):
        assert out.fields[key].shape == (40, 44, k)
    for key in ("emulsion/silver_physical", "emulsion/halide", "emulsion/growth_rate"):
        assert out.fields[key].shape == (40, 44)
    meta = out.meta["develop"]
    assert set(meta) == DEVELOP_META
    assert meta["chem_cell_um"] == 24.0 and meta["block_cells"] == 8
    assert meta["prepass"]["bins"] == 64 and meta["prepass"]["passes"] >= 1
    assert len(meta["prepass"]["delta"]) == meta["prepass"]["passes"]
    m = silver_total(out)
    assert np.all(np.isfinite(m)) and np.all(m >= 0.0)
    # More exposure, more silver (along the ramp, away from the unexposed patch).
    row = m[-1]
    assert row[-1] > row[0] + 0.05
    # The bulk lost developer and gained bromide.
    assert meta["bulk_final"]["bromide"] > meta["bulk_initial"]["bromide"]
    assert len(meta["zero_exposure"]["silver"]) == k
    assert out.history[-1].kind == "develop"


@pytest.mark.parametrize(
    ("tile_px", "workers", "mapped"),
    [(None, 1, False), (7, 1, False), (16, 2, False), (None, 2, True)],
)
def test_develop_independent_of_tiling_and_workers(
    lib, film, tiling_reference, monkeypatch, tile_px, workers, mapped, tmp_path
):
    # Spec D1 3.6: bit-identical whatever tile_px and the number of workers; (None, 1) is a
    # rerun of the reference (determinism).
    monkeypatch.setattr(dstep, "BLOCK_CELLS", 8)
    monkeypatch.setattr(dstep, "halo_um", lambda *a: SMALL_HALO_UM)
    state = make_state(film, *scene(), 12.0)
    if mapped:
        # The run's store maps every field: workers write straight into the outputs and the
        # exposure files are released after development (spec D1 4.4, 4.5).
        store = FieldStore(tmp_path / "store", ram_threshold_bytes=0)
        ev = store.allocate("expose-x:emulsion/exposure_volume", (40, 44, N_Z))
        es = store.allocate("expose-x:emulsion/exposure_surface", (40, 44, N_Z))
        ev[...], es[...] = scene()
        state = make_state(
            film,
            store.freeze(ev_name := "expose-x:emulsion/exposure_volume"),
            store.freeze("expose-x:emulsion/exposure_surface"),
            12.0,
        )
        del ev, es
        ctx = _ctx(lib, tile_px=tile_px, workers=workers, store=store)
    else:
        ctx = _ctx(lib, tile_px=tile_px, workers=workers, scratch_dir=tmp_path)
    out = _process(state, ctx, develop_step())
    if mapped:
        assert ev_name not in store
        assert store.is_mapped(store.name_of(out.fields["emulsion/silver"]))
    for key, value in tiling_reference.fields.items():
        assert np.array_equal(out.fields[key], value), key
    assert out.meta["develop"] == tiling_reference.meta["develop"]


def test_uniform_frame_matches_strip(lib, film):
    # Interior of a uniformly exposed frame vs develop_strip of the same exposure, within the
    # threshold of test_image_profile_uniform_patch (spec D1 3.8): |dM| <= 0.014 g/m2.
    ev = np.full((48, 48, N_Z), 30.0)
    ev[:, :, 1] = 20.0
    es = 0.01 * ev
    step = develop_step()
    out = dstep.develop_frame(
        make_state(film, ev, es, 24.0), step, _step_ctx(lib), block_cells=12, halo_cells=3
    )
    strip = dstep.develop_strip(
        film, "135", (ev[0, :1], es[0, :1]), EXPOSURE_TIME_S, step, _ctx(lib)
    )
    m_strip = float(strip.silver.sum() + strip.silver_physical[0])
    centre = silver_total(out)[16:32, 16:32]
    assert m_strip > 0.1
    assert np.max(np.abs(centre - m_strip)) <= 0.014


def test_strip_image_area_fraction(lib, film):
    # The strip returns only its own steps; the unexposed remainder still depletes the tank.
    ev, es = _h.strip_exposures(3, N_Z)
    full = develop_step(tank={"volume_ml": 300.0, "film_area_cm2": ROLL_CM2})
    one = develop_step(
        tank={
            "volume_ml": 300.0,
            "film_area_cm2": ROLL_CM2,
            "image_area_fraction": {"value": 1.0, "unit": "1", "source": "assumed"},
        }
    )
    a = dstep.develop_strip(film, "135", (ev, es), EXPOSURE_TIME_S, full, _ctx(lib))
    b = dstep.develop_strip(film, "135", (ev, es), EXPOSURE_TIME_S, one, _ctx(lib))
    assert a.silver.shape == b.silver.shape == (3, a.silver.shape[1])
    assert not np.array_equal(a.silver, b.silver)
    c = dstep.develop_strip(film, "135", (ev, es), EXPOSURE_TIME_S, full, _ctx(lib))
    assert np.array_equal(a.silver, c.silver)


def test_presoak_writes_broadcast_flag(lib, film):
    state = make_state(film, *scene((4, 5)), 12.0)
    step = ProcessStep(kind=StepKind.PRESOAK, time="60s")
    out = dstep.presoak_bw(state, step, StepContext(_ctx(lib), 0, None))
    flag = out.fields["frame/presoaked"]
    assert flag.dtype == bool and flag.shape == (4, 5) and flag.strides == (0, 0)
    assert bool(flag.all())
    assert dstep._presoaked(out) and not dstep._presoaked(state)


def test_halo_um_caps_and_grows(lib, film, d76, kin_test):
    f = _h.film_column_model(film, "135", EXPOSURE_TIME_S, kin_test, d76, "image")
    stand = develop_step(time="60min", agitation=Agitation(method="stand", initial_s=30.0))
    cap = 4.6 * f.grid.delta_max_um
    assert dstep.halo_um(f, d76, kin_test, stand) == pytest.approx(cap, rel=1e-12)
    short = dstep.halo_um(f, d76, kin_test, develop_step(time="60s"))
    assert 0.0 < short <= cap
    no_agitation = develop_step(time="60s", agitation=None)
    assert dstep.halo_um(f, d76, kin_test, no_agitation) >= short


def test_develop_bad_input(lib, film):
    state = make_state(film, *scene((6, 6)), 12.0)
    sctx = _step_ctx(lib)
    with pytest.raises(ValueError, match="exposure fields"):
        dstep.develop_frame(state.drop("emulsion/exposure_volume"), develop_step(), sctx)
    with pytest.raises(ValueError, match="exposure"):
        dstep.develop_frame(state.with_meta(exposure={}), develop_step(), sctx)
    with pytest.raises(ValueError, match="film_area_cm2"):
        dstep.develop_frame(state, develop_step(tank={"volume_ml": 300.0}), sctx)
    no_temp = develop_step().model_copy(update={"temperature": None})
    with pytest.raises(ValueError, match="temperature"):
        dstep.develop_frame(state, no_temp, sctx)
    with pytest.raises(ValueError, match="bath"):
        dstep.develop_frame(state, develop_step(), StepContext(_ctx(lib), 0, None))
    with pytest.raises(ValueError, match="block_cells"):
        dstep.develop_frame(state, develop_step(), sctx, block_cells=0)
    with pytest.raises(ValueError, match="halo_cells"):
        dstep.develop_frame(state, develop_step(), sctx, halo_cells=-1)
    window = dataclasses.replace(state, geometry=FrameGeometry((6, 6), 12.0, (1, 1), (8, 8)))
    with pytest.raises(ValueError, match="whole frame"):
        dstep.develop_frame(window, develop_step(), sctx)
    rolls = develop_step(tank={"volume_ml": 300.0, "film_area_cm2": 10.0, "rolls_before": 1})
    with pytest.raises(NotImplementedError):
        dstep.develop_frame(state, rolls, sctx)
    ev, es = _h.strip_exposures(2, N_Z)
    with pytest.raises(ValueError, match="n_z"):
        dstep.develop_strip(
            film, "135", (ev[:, :1], es[:, :1]), EXPOSURE_TIME_S, develop_step(), _ctx(lib)
        )
    with pytest.raises(ValueError, match="shape"):
        dstep.develop_strip(film, "135", (ev, es[:1]), EXPOSURE_TIME_S, develop_step(), _ctx(lib))


# ---------------------------------------------------------------------------------------
# Halo convergence (spec D1 3.6, 3.8; slow)
# ---------------------------------------------------------------------------------------

LARGE_TANK = Agitation(method="lift", initial_s=22.5, interval_s=60.0, duration_s=7.5)
STAND = Agitation(method="stand", initial_s=15.0)


@pytest.mark.slow
@pytest.mark.parametrize(
    "agitation", [SMALL_TANK, LARGE_TANK, STAND], ids=["small", "large", "stand"]
)
@pytest.mark.parametrize("layers_um", [(50.0, 300.0), (100.0, 500.0)], ids=["typical", "upper"])
def test_halo_convergence(lib, film, monkeypatch, agitation, layers_um):
    # Halo x 1.5 changes the developed silver by <= 3e-3 g/m2 (spec D1 3.8). One block holds
    # the frame (an exposure edge next to unexposed film outside the frame); 60 s.
    delta_agit, delta_max = layers_um
    real = dstep.resolve_kinetics

    def kinetics(*a):
        kin = _h.with_test_rates(real(*a))
        return _h.with_transport(kin, delta_agit=(delta_agit, "um"), delta_max=(delta_max, "um"))

    monkeypatch.setattr(dstep, "resolve_kinetics", kinetics)
    ev = np.zeros((8, 8, N_Z))
    ev[:, 4:] = 100.0
    state = make_state(film, ev, 0.01 * ev, 24.0)
    step = develop_step(agitation=agitation)
    base = dstep.develop_frame(state, step, _step_ctx(lib), block_cells=8)
    halo = base.meta["develop"]["halo_um"]
    wide = dstep.develop_frame(
        state, step, _step_ctx(lib), block_cells=8, halo_cells=math.ceil(1.5 * halo / 24.0)
    )
    assert np.max(np.abs(silver_total(wide) - silver_total(base))) <= 3e-3
