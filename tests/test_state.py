import numpy as np
import pytest

from darkroom_atlas.data.library import Library
from darkroom_atlas.spectral import DEFAULT_GRID
from darkroom_atlas.state import EmulsionState, FrameGeometry, StepRecord
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
