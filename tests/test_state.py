import math
import time

import numpy as np
import pytest

from darkroom_atlas.data.library import Library
from darkroom_atlas.fieldstore import FieldStore, open_field
from darkroom_atlas.spectral import DEFAULT_GRID
from darkroom_atlas.state import (
    EmulsionState,
    FieldSpec,
    FrameAssembler,
    FrameGeometry,
    StepRecord,
)
from darkroom_atlas.tiling import iter_tiles

H, W = 50, 70


@pytest.fixture(scope="module")
def material():
    return Library().material("ideal")


@pytest.fixture
def empty(material):
    return EmulsionState(material, FrameGeometry.from_long_edge((H, W)), DEFAULT_GRID, seed=3)


def test_geometry_from_long_edge():
    geom = FrameGeometry.from_long_edge((4000, 6000), 36.0)
    assert geom.pixel_pitch_um == 6.0
    assert geom.um_to_px(12) == 2.0
    assert geom.shape_hw == (4000, 6000)
    assert geom.origin_px == (0, 0)
    assert geom.frame_shape_hw == (4000, 6000)
    assert FrameGeometry.from_long_edge((6000, 4000)).pixel_pitch_um == 6.0


def test_geometry_rejects_bad_values():
    with pytest.raises(ValueError):
        FrameGeometry((0, 10), 6.0)
    with pytest.raises(ValueError):
        FrameGeometry((10, 10), 0.0)
    with pytest.raises(ValueError):
        FrameGeometry((10, 10), 6.0, origin_px=(5, 0), frame_shape_hw=(12, 10))


def test_with_fields_readonly_and_immutable(empty):
    source = np.ones((H, W))
    state = empty.with_fields({"x/exposure": source})
    assert state is not empty
    assert not empty.has("x", "exposure")
    assert state.has("x", "exposure")
    stored = state.get("x", "exposure")
    assert not stored.flags.writeable
    with pytest.raises(ValueError):
        stored[0, 0] = 2.0
    assert source.flags.writeable  # the caller's array is left alone ...
    source[0, 0] = 7.0
    assert stored[0, 0] == 1.0  # ... and was copied
    with pytest.raises(TypeError):
        state.fields["x/exposure"] = source
    later = state.with_fields({"y/exposure": np.zeros((H, W))})
    assert later.get("x", "exposure") is stored  # existing read-only fields are shared
    assert not state.has("y", "exposure")


def test_field_shape_validated(empty):
    with pytest.raises(ValueError):
        empty.with_fields({"x/exposure": np.ones((H, W + 1))})
    with pytest.raises(ValueError):
        empty.with_fields({"x/exposure": np.ones(H)})
    with pytest.raises(ValueError, match="nope"):
        empty.with_fields({"nope/exposure": np.ones((H, W))})
    with pytest.raises(ValueError):
        empty.with_fields({"exposure": np.ones((H, W))})
    assert empty.with_fields({"frame/mask": np.zeros((H, W, 2), bool)}).has("frame", "mask")


def test_get_missing_lists_keys(empty):
    state = empty.with_fields({"x/exposure": np.ones((H, W)), "frame/mask": np.ones((H, W))})
    with pytest.raises(KeyError) as info:
        state.get("y", "exposure")
    message = str(info.value)
    assert "y/exposure" in message
    assert "x/exposure" in message
    assert "frame/mask" in message


def test_crop_tracks_origin(empty):
    data = np.arange(H * W, dtype=np.float64).reshape(H, W)
    state = empty.with_fields({"x/exposure": data})
    once = state.crop((slice(10, 30), slice(5, 25)))
    assert once.geometry.shape_hw == (20, 20)
    assert once.geometry.origin_px == (10, 5)
    assert once.geometry.frame_shape_hw == (H, W)
    assert once.geometry.pixel_pitch_um == state.geometry.pixel_pitch_um
    np.testing.assert_array_equal(once.get("x", "exposure"), data[10:30, 5:25])
    twice = once.crop((slice(2, 12), slice(3, 20)))
    assert twice.geometry.shape_hw == (10, 17)
    assert twice.geometry.origin_px == (12, 8)
    assert twice.geometry.frame_shape_hw == (H, W)
    np.testing.assert_array_equal(twice.get("x", "exposure"), data[12:22, 8:25])
    assert state.geometry.origin_px == (0, 0)
    with pytest.raises(ValueError):
        state.crop((slice(10, 10), slice(0, 5)))
    with pytest.raises(ValueError):
        state.crop((slice(0, 10, 2), slice(0, 5)))


def test_paste_roundtrip(empty):
    rng = np.random.default_rng(0)
    full = empty.with_fields(
        {
            "x/exposure": rng.random((H, W)),
            "y/counts": rng.integers(0, 100, (H, W, 3)),
            "frame/mask": rng.random((H, W)) > 0.5,
        }
    )
    out = empty
    for tile in iter_tiles((H, W), 16, 5):
        out = out.paste(full.crop(tile.outer), tile)
    assert set(out.fields) == set(full.fields)
    for key, value in full.fields.items():
        assert out.fields[key].dtype == value.dtype
        np.testing.assert_array_equal(out.fields[key], value)
        assert not out.fields[key].flags.writeable
    assert out.geometry == empty.geometry
    assert empty.fields == {}


def test_paste_into_existing_field(empty):
    base = empty.with_fields({"x/exposure": np.zeros((H, W))})
    tile = next(iter(iter_tiles((H, W), 16, 4)))
    patch = base.crop(tile.outer).with_fields({"x/exposure": np.ones((20, 20))})
    out = base.paste(patch, tile)
    expected = np.zeros((H, W))
    expected[tile.inner] = 1.0
    np.testing.assert_array_equal(out.get("x", "exposure"), expected)
    np.testing.assert_array_equal(base.get("x", "exposure"), np.zeros((H, W)))


@pytest.mark.parametrize(
    "patch_value",
    [np.ones((20, 20), np.float32), np.ones((20, 20, 2))],
    ids=["dtype", "trailing-shape"],
)
def test_paste_rejects_incompatible_field(empty, patch_value):
    # Restarting the field from zeros would silently discard what was already assembled.
    base = empty.with_fields({"x/exposure": np.full((H, W), 5.0)})
    tile = next(iter(iter_tiles((H, W), 16, 4)))
    patch = base.crop(tile.outer).with_fields({"x/exposure": patch_value})
    with pytest.raises(ValueError, match="x/exposure"):
        base.paste(patch, tile)


def test_paste_rejects_mismatched_tile(empty):
    tiles = list(iter_tiles((H, W), 16, 4))
    wrong = empty.crop(tiles[1].outer)
    with pytest.raises(ValueError):
        empty.paste(wrong, tiles[0])


def test_history_and_meta(empty):
    record = StepRecord(index=0, kind="identity", bath=None, time_s=None)
    state = empty.with_history(record)
    assert state.history == (record,)
    assert empty.history == ()
    assert record.note == ""
    again = state.with_history(StepRecord(1, "wash", None, 300.0, note="running water"))
    assert [r.index for r in again.history] == [0, 1]
    tagged = again.with_meta(light={"iso": 100.0}, run="a")
    assert tagged.meta["light"]["iso"] == 100.0
    assert tagged.meta["run"] == "a"
    assert "run" not in again.meta
    assert tagged.with_meta(run="b").meta == {"light": {"iso": 100.0}, "run": "b"}
    with pytest.raises(TypeError):
        tagged.meta["run"] = "c"
    with pytest.raises(TypeError):
        tagged.meta["light"]["iso"] = 200.0
    assert tagged.history == again.history
    assert tagged.seed == 3


# --------------------------------------------------------------------------------------
# drop and FrameAssembler (spec §4.4)
# --------------------------------------------------------------------------------------


def test_drop_fields(empty):
    a = np.ones((H, W))
    b = np.zeros((H, W, 2))
    state = empty.with_fields({"x/a": a, "x/b": b, "frame/c": a}).with_meta(note="m")
    dropped = state.drop("x/a", "frame/c")
    assert set(dropped.fields) == {"x/b"}
    assert dropped.fields["x/b"] is state.fields["x/b"]
    assert dropped.meta == state.meta and dropped.history == state.history
    assert dropped.geometry == state.geometry
    assert set(state.fields) == {"x/a", "x/b", "frame/c"}  # the original is unchanged
    assert state.drop().fields == state.fields
    with pytest.raises(KeyError, match="available fields"):
        state.drop("x/missing")


def _write_all(assembler, key, full, tile_px):
    for tile in iter_tiles(full.shape[:2], tile_px, 0):
        assembler.write(key, tile.inner, full[tile.inner])


def test_assembler_each_pixel_once(tmp_path):
    geom = FrameGeometry.from_long_edge((H, W))
    full = np.random.default_rng(0).random((H, W, 3))
    tiles = list(iter_tiles((H, W), 16, 0))

    asm = FrameAssembler(geom, [FieldSpec("x/a", (3,))])
    _write_all(asm, "x/a", full, 16)
    out = asm.finish()
    np.testing.assert_array_equal(out["x/a"], full)
    assert not out["x/a"].flags.writeable

    # A missed tile raises in finish.
    asm = FrameAssembler(geom, [FieldSpec("x/a", (3,))])
    for tile in tiles[1:]:
        asm.write("x/a", tile.inner, full[tile.inner])
    with pytest.raises(ValueError, match="exactly once"):
        asm.finish()

    # A rewritten (or overlapping) region raises at once.
    asm = FrameAssembler(geom, [FieldSpec("x/a", (3,))])
    asm.write("x/a", tiles[0].inner, full[tiles[0].inner])
    with pytest.raises(ValueError, match="overlaps"):
        asm.write("x/a", tiles[0].inner, full[tiles[0].inner])
    with pytest.raises(ValueError, match="overlaps"):
        asm.write("x/a", (slice(10, 20), slice(10, 20)), full[10:20, 10:20])

    # Uniform writes are tracked the same way.
    asm = FrameAssembler(geom, [FieldSpec("x/u")])
    asm.write("x/u", tiles[0].inner, np.broadcast_to(1.0, (16, 16)))
    with pytest.raises(ValueError, match="overlaps"):
        asm.write("x/u", tiles[0].inner, np.broadcast_to(1.0, (16, 16)))
    with pytest.raises(ValueError, match="exactly once"):
        asm.finish()

    # A field never written raises in finish.
    asm = FrameAssembler(geom, [FieldSpec("x/a", (3,)), FieldSpec("x/b")])
    _write_all(asm, "x/a", full, 16)
    with pytest.raises(ValueError, match="x/b"):
        asm.finish()

    # mark: regions written straight into the mapped outputs count like write() (§4.4).
    with FieldStore(tmp_path, ram_threshold_bytes=0) as store:
        asm = FrameAssembler(geom, [FieldSpec("x/a", (3,))], store, name_prefix="s4:")
        with pytest.raises(ValueError, match="descriptors"):
            asm.mark(tiles[0].inner)  # nothing handed out yet
        (desc,) = asm.descriptors().values()
        assert asm.store_names == {"x/a": "s4:x/a"}
        writer = open_field(desc, "r+")
        for tile in tiles[1:]:
            writer[tile.inner] = full[tile.inner]
            asm.mark(tile.inner)
        writer.flush()
        with pytest.raises(ValueError, match="overlaps"):
            asm.mark(tiles[1].inner)  # marked twice
        with pytest.raises(ValueError, match="overlaps"):
            asm.write("x/a", tiles[2].inner, full[tiles[2].inner])  # marked, then written
        with pytest.raises(ValueError, match="exactly once"):
            asm.finish()  # tile 0 missing
        asm.write("x/a", tiles[0].inner, full[tiles[0].inner])
        with pytest.raises(ValueError, match="overlaps"):
            asm.mark(tiles[0].inner)  # written, then marked

        # reset: forget the pass, keep the outputs, and overwrite them in a second pass.
        asm.reset()
        with pytest.raises(ValueError, match="exactly once"):
            asm.finish()
        assert asm.descriptors() == {"x/a": desc}  # same allocation
        second = 2.0 * full
        for tile in tiles:
            writer[tile.inner] = second[tile.inner]
            asm.mark(tile.inner)
        writer.flush()
        del writer
        out = asm.finish()["x/a"]
        assert isinstance(out, np.memmap) and out.mode == "r"
        np.testing.assert_array_equal(out, second)
        with pytest.raises(ValueError, match="finished"):
            asm.reset()
        with pytest.raises(ValueError, match="finished"):
            asm.mark(tiles[0].inner)
        with pytest.raises(ValueError, match="finished"):
            asm.descriptors()

    # reset also works for written (not marked) fields and for uniform ones.
    asm = FrameAssembler(geom, [FieldSpec("x/a", (3,)), FieldSpec("x/u")])
    _write_all(asm, "x/a", full, 16)
    for tile in tiles:
        shape = full[tile.inner].shape[:2]
        asm.write("x/u", tile.inner, np.broadcast_to(1.0, shape))
    asm.reset()
    _write_all(asm, "x/a", full[::-1, ::-1], 25)
    for tile in iter_tiles((H, W), 25, 0):
        shape = full[tile.inner].shape[:2]
        asm.write("x/u", tile.inner, np.broadcast_to(2.0, shape))
    out = asm.finish()
    np.testing.assert_array_equal(out["x/a"], full[::-1, ::-1])
    assert out["x/u"].strides == (0, 0) and (out["x/u"] == 2.0).all()

    # Fields in RAM have no descriptor; mark needs one.
    asm = FrameAssembler(geom, [FieldSpec("x/a", (3,))])
    assert asm.descriptors() == {}
    with pytest.raises(ValueError, match="descriptors"):
        asm.mark(tiles[0].inner)


def test_assembler_uniform_fields_stay_broadcast():
    geom = FrameGeometry.from_long_edge((H, W))
    one = np.ones((1, 1))
    one.flags.writeable = False
    asm = FrameAssembler(geom, [FieldSpec("x/u"), FieldSpec("x/v", (2,)), FieldSpec("x/m")])
    for i, tile in enumerate(iter_tiles((H, W), 16, 0)):
        shape = (tile.inner[0].stop - tile.inner[0].start, tile.inner[1].stop - tile.inner[1].start)
        asm.write("x/u", tile.inner, np.broadcast_to(one, shape))
        asm.write("x/v", tile.inner, np.broadcast_to(np.array([[[1.0, 2.0]]]), (*shape, 2)))
        # Different scalars in different tiles: allocated and filled in correctly.
        asm.write("x/m", tile.inner, np.broadcast_to(float(i % 2), shape))
    out = asm.finish()
    assert asm.store_names.keys() == {"x/m"}
    for key in ("x/u", "x/v"):
        assert out[key].strides[:2] == (0, 0) and out[key].base is not None
        assert not out[key].flags.writeable
    assert out["x/u"].shape == (H, W) and (out["x/u"] == 1.0).all()
    assert out["x/v"].shape == (H, W, 2) and (out["x/v"] == [1.0, 2.0]).all()
    expected = np.zeros((H, W))
    for i, tile in enumerate(iter_tiles((H, W), 16, 0)):
        expected[tile.inner] = i % 2
    np.testing.assert_array_equal(out["x/m"], expected)
    # Into an EmulsionState without materialising the broadcast view.
    state = EmulsionState(
        Library().material("ideal"), geom, DEFAULT_GRID, seed=0, fields={"frame/u": out["x/u"]}
    )
    assert state.fields["frame/u"] is out["x/u"]


def test_assembler_uniform_then_real_write_backfills():
    geom = FrameGeometry.from_long_edge((H, W))
    full = np.random.default_rng(2).random((H, W))
    tiles = list(iter_tiles((H, W), 16, 0))
    asm = FrameAssembler(geom, [FieldSpec("x/a")])
    asm.write("x/a", tiles[0].inner, np.broadcast_to(7.0, full[tiles[0].inner].shape))
    for tile in tiles[1:]:
        asm.write("x/a", tile.inner, full[tile.inner])
    out = asm.finish()["x/a"]
    expected = full.copy()
    expected[tiles[0].inner] = 7.0
    np.testing.assert_array_equal(out, expected)


def test_assembler_uses_store_memmap(tmp_path):
    geom = FrameGeometry.from_long_edge((H, W))
    full = np.random.default_rng(3).random((H, W, 2))
    with FieldStore(tmp_path, ram_threshold_bytes=0) as store:
        asm = FrameAssembler(geom, [FieldSpec("x/a", (2,))], store=store)
        _write_all(asm, "x/a", full, 20)
        out = asm.finish()["x/a"]
        assert isinstance(out, np.memmap) and out.mode == "r"
        np.testing.assert_array_equal(out, full)
        assert asm.store_names["x/a"] in store
        # Two assemblers producing the same field key share a store under distinct
        # prefixes (steps use f"s{index}:", spec §4.3); the same name twice raises.
        asm2 = FrameAssembler(geom, [FieldSpec("x/a", (2,))], store=store, name_prefix="s7:")
        _write_all(asm2, "x/a", full, 25)
        assert asm2.store_names == {"x/a": "s7:x/a"}
        np.testing.assert_array_equal(asm2.finish()["x/a"], out)
        asm3 = FrameAssembler(geom, [FieldSpec("x/a", (2,))], store=store, name_prefix="s7:")
        with pytest.raises(ValueError, match="already allocated"):
            _write_all(asm3, "x/a", full, 25)


def test_assembler_write_state_matches_paste(empty):
    rng = np.random.default_rng(4)
    source = empty.with_fields({"x/a": rng.random((H, W, 2)), "x/b": rng.random((H, W))})
    pasted = empty
    asm = FrameAssembler(empty.geometry, [FieldSpec("x/a", (2,)), FieldSpec("x/b")])
    for tile in iter_tiles((H, W), 16, 5):
        tile_state = source.crop(tile.outer)
        pasted = pasted.paste(tile_state, tile)
        asm.write_state(tile_state, tile)
    out = asm.finish()
    for key in ("x/a", "x/b"):
        np.testing.assert_array_equal(out[key], pasted.fields[key])


def test_assembler_rejects_bad_input(empty):
    geom = empty.geometry
    with pytest.raises(ValueError, match="duplicate"):
        FrameAssembler(geom, [FieldSpec("x/a"), FieldSpec("x/a")])
    with pytest.raises(ValueError):
        FrameAssembler(geom, [])
    with pytest.raises(ValueError):
        FieldSpec("x/a", (0,))
    with pytest.raises(ValueError):
        FieldSpec("x/a", dtype="U3")
    asm = FrameAssembler(geom, [FieldSpec("x/a", (2,))])
    with pytest.raises(KeyError):
        asm.write("x/b", (slice(0, 2), slice(0, 2)), np.zeros((2, 2, 2)))
    with pytest.raises(ValueError, match="shape"):
        asm.write("x/a", (slice(0, 2), slice(0, 2)), np.zeros((2, 2)))
    with pytest.raises(ValueError):
        asm.write("x/a", (slice(0, 2), slice(0, W + 1)), np.zeros((2, W + 1, 2)))
    with pytest.raises(ValueError, match="cast"):
        asm.write("x/a", (slice(0, 2), slice(0, 2)), np.zeros((2, 2, 2), dtype=complex))
    tile = next(iter_tiles((H, W), 16, 2))
    with pytest.raises(KeyError, match="lacks"):
        asm.write_state(empty.crop(tile.outer), tile)
    with pytest.raises(ValueError, match="geometry"):
        asm.write_state(empty.with_fields({"x/a": np.zeros((H, W, 2))}), tile)
    asm.write("x/a", (slice(0, H), slice(0, W)), np.zeros((H, W, 2), dtype=np.float32))
    asm.finish()
    with pytest.raises(ValueError, match="finished"):
        asm.finish()
    with pytest.raises(ValueError, match="finished"):
        asm.write("x/a", (slice(0, 2), slice(0, 2)), np.zeros((2, 2, 2)))


class _CountingArray(np.ndarray):
    """An array that counts the elements assigned into it."""

    assigned = 0

    def __setitem__(self, key, value):
        _CountingArray.assigned += np.asarray(self[key]).size
        super().__setitem__(key, value)


class _CountingStore(FieldStore):
    """A RAM store whose fields count every element written into them."""

    def allocate(self, name, shape, dtype="float64"):
        array = super().allocate(name, shape, dtype).view(_CountingArray)
        self._entries[name].array = array
        return array


def _assemble(blocks_per_side, block_px, values, store=None):
    side = blocks_per_side * block_px
    asm = FrameAssembler(
        FrameGeometry.from_long_edge((side, side)), [FieldSpec("x/a", values.shape[2:])], store
    )
    for tile in iter_tiles((side, side), block_px, 0):
        asm.write("x/a", tile.inner, values)
    return asm.finish()["x/a"]


def test_assembler_linear_cost():
    # Spec §10 D1a: assembling costs O(pixels), not O(blocks x pixels) as D0's repeated
    # paste did. Deterministic part: every pixel is copied into the frame exactly once, for
    # 4, 16 and 64 blocks (D0's paste copied the whole frame once per block).
    block_px, trailing = 64, 6
    values = np.ones((block_px, block_px, trailing))
    for blocks_per_side in (2, 4, 8):
        _CountingArray.assigned = 0
        out = _assemble(blocks_per_side, block_px, values, _CountingStore(None))
        assert _CountingArray.assigned == out.size

    # Timing part: 64 and 16 blocks against 4 (area ratios 16 and 4), each within 1.5 x the
    # ratio of a plain block-by-block NumPy copy into a new frame of the same size, which is
    # the O(pixels) floor. Normalising by that floor rather than by the area ratio cancels
    # the cache and memory-bandwidth effects of the larger frames, which a loaded machine
    # amplifies; D0's paste would still come out 16 times above it at 64 blocks. Best of
    # many repetitions, sizes interleaved.
    sizes = (2, 4, 8)
    best = {(n, kind): math.inf for n in sizes for kind in ("assembler", "copy")}
    for _ in range(15):
        for n in sizes:
            start = time.perf_counter()
            _assemble(n, block_px, values)
            middle = time.perf_counter()
            _plain_copy(n, block_px, values)
            end = time.perf_counter()
            best[n, "assembler"] = min(best[n, "assembler"], middle - start)
            best[n, "copy"] = min(best[n, "copy"], end - middle)
    relative = {n: best[n, "assembler"] / best[n, "copy"] for n in sizes}
    growth = {n: relative[n] / relative[2] for n in sizes}
    assert growth[8] <= 1.5, (growth, best)
    assert growth[4] <= 1.5, (growth, best)


def _plain_copy(blocks_per_side, block_px, values):
    side = blocks_per_side * block_px
    out = np.zeros((side, side, *values.shape[2:]))
    for tile in iter_tiles((side, side), block_px, 0):
        out[tile.inner] = values
    return out
