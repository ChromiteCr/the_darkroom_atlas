import numpy as np
import pytest

from darkroom_atlas._frozen import readonly_array
from darkroom_atlas.data.library import Library
from darkroom_atlas.fieldstore import FieldStore
from darkroom_atlas.spectral import DEFAULT_GRID
from darkroom_atlas.state import EmulsionState, FrameGeometry


def test_fieldstore_memmap_readonly_shared(tmp_path):
    with FieldStore(tmp_path, ram_threshold_bytes=1024) as store:
        big = store.allocate("x/exposure", (40, 30, 2))  # 19200 bytes > threshold
        assert isinstance(big, np.memmap) and big.mode == "w+"
        big[...] = np.arange(big.size, dtype=float).reshape(big.shape)
        frozen = store.freeze("x/exposure")
        assert isinstance(frozen, np.memmap) and frozen.mode == "r"
        assert store.freeze("x/exposure") is frozen
        np.testing.assert_array_equal(frozen.ravel(), np.arange(frozen.size, dtype=float))

        # readonly_array shares the read-only map and its views without copying.
        assert readonly_array(frozen) is frozen
        view = frozen[5:10, 3:7]
        assert readonly_array(view) is view
        plain = np.asarray(frozen)
        assert readonly_array(plain) is plain
        with pytest.raises(ValueError):
            frozen[0, 0, 0] = 1.0
        with pytest.raises(ValueError):
            view[0, 0, 0] = 1.0

        # The path the pipeline takes: the map goes into an EmulsionState unchanged.
        material = Library().material("ideal")
        state = EmulsionState(
            material,
            FrameGeometry.from_long_edge((40, 30)),
            DEFAULT_GRID,
            seed=0,
            fields={"frame/exposure": frozen},
        )
        stored = state.fields["frame/exposure"]
        assert stored is frozen
        cropped = state.crop((slice(2, 9), slice(4, 20))).fields["frame/exposure"]
        assert np.shares_memory(cropped, frozen)

        small = store.allocate("x/small", (4, 4))
        assert not isinstance(small, np.memmap) and not store.is_mapped("x/small")
        assert store.is_mapped("x/exposure")


def test_readonly_array_copies_writable_maps(tmp_path):
    path = tmp_path / "m.dat"
    for mode in ("w+", "r+", "c"):
        m = np.memmap(path, dtype=float, mode=mode, shape=(3, 3))
        out = readonly_array(m)
        assert not np.shares_memory(out, m)
        assert not out.flags.writeable
        del m


def test_fieldstore_ram_without_root():
    store = FieldStore(None, ram_threshold_bytes=0)
    a = store.allocate("x/a", (100, 100))
    assert type(a) is np.ndarray and a.flags.writeable and not a.any()
    frozen = store.freeze("x/a")
    assert frozen is a and not frozen.flags.writeable
    assert readonly_array(frozen) is frozen


def test_fieldstore_threshold_is_strict(tmp_path):
    store = FieldStore(tmp_path, ram_threshold_bytes=800)
    store.allocate("x/at", (10, 10))  # exactly 800 bytes: stays in RAM
    store.allocate("x/above", (10, 10, 2))
    assert not store.is_mapped("x/at")
    assert store.is_mapped("x/above")


def test_fieldstore_release_and_close(tmp_path):
    store = FieldStore(tmp_path, ram_threshold_bytes=0)
    store.allocate("x/a", (8, 8))
    store.allocate("x/b", (8, 8))
    files = sorted(tmp_path.iterdir())
    assert len(files) == 2
    assert all("/" not in f.name for f in files)
    store.freeze("x/a")
    store.release("x/a")
    assert "x/a" not in store
    assert len(list(tmp_path.iterdir())) == 1
    store.allocate("x/a", (2, 2))  # the key is free again
    store.close()
    assert list(tmp_path.iterdir()) == []
    assert tmp_path.exists()  # the root belongs to the caller
    with pytest.raises(RuntimeError):
        store.allocate("x/c", (2, 2))


def test_fieldstore_context_manager_removes_files(tmp_path):
    root = tmp_path / "scratch"
    with FieldStore(root, ram_threshold_bytes=0) as store:
        store.allocate("x/a", (4, 4))
        assert len(list(root.iterdir())) == 1
    assert list(root.iterdir()) == []


def test_fieldstore_rejects_bad_use(tmp_path):
    store = FieldStore(tmp_path)
    store.allocate("x/a", (2, 2))
    with pytest.raises(ValueError, match="already allocated"):
        store.allocate("x/a", (2, 2))
    with pytest.raises(KeyError, match="x/a"):
        store.freeze("x/missing")
    with pytest.raises(KeyError):
        store.release("x/missing")
    with pytest.raises(ValueError):
        store.allocate("x/b", (2, -1))
    with pytest.raises(ValueError):
        store.allocate("", (2, 2))
    with pytest.raises(ValueError):
        FieldStore(tmp_path, ram_threshold_bytes=-1)
    with pytest.raises(TypeError):
        FieldStore(tmp_path, ram_threshold_bytes=1.5)
