import pickle
from pathlib import Path

import numpy as np
import pytest

from darkroom_atlas._frozen import readonly_array
from darkroom_atlas.data.library import Library
from darkroom_atlas.fieldstore import FieldDescriptor, FieldStore, descriptor_of, open_field
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
    with pytest.raises(ValueError, match="already allocated"):
        store.allocate("x/a", (2, 2))  # names stay used after release (spec §4.3)
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


def test_fieldstore_unique_backing_files(tmp_path):
    # Spec §4.3, §10 D1a: two steps writing the same field key get different files, and the
    # second allocation does not truncate the first field's read-only map.
    with FieldStore(tmp_path, ram_threshold_bytes=0) as store:
        first = store.allocate("s2:emulsion/silver", (6, 5, 3))
        first[...] = np.arange(first.size, dtype=float).reshape(first.shape)
        frozen = store.freeze("s2:emulsion/silver")
        second = store.allocate("s3:emulsion/silver", (6, 5, 3))
        second[...] = -1.0
        store.freeze("s3:emulsion/silver")
        paths = {store.descriptor(n).path for n in ("s2:emulsion/silver", "s3:emulsion/silver")}
        assert len(paths) == 2
        for path in paths:
            assert path.parent == tmp_path and path.suffix == ".bin"
            assert "/" not in path.name and ":" not in path.name
            stem, _, tag = path.stem.rpartition("-")
            assert stem in ("s2_emulsion_silver", "s3_emulsion_silver")
            assert len(tag) == 12 and all(c in "0123456789abcdef" for c in tag)
        np.testing.assert_array_equal(frozen.ravel(), np.arange(frozen.size, dtype=float))
        with pytest.raises(ValueError, match="already allocated"):
            store.allocate("s2:emulsion/silver", (6, 5, 3))
        # The same name in two stores of one directory still gives two files.
        with FieldStore(tmp_path, ram_threshold_bytes=0) as other:
            other.allocate("s2:emulsion/silver", (6, 5, 3))
            assert other.descriptor("s2:emulsion/silver").path not in paths
        np.testing.assert_array_equal(frozen.ravel(), np.arange(frozen.size, dtype=float))


def test_fieldstore_descriptor_roundtrip(tmp_path):
    with FieldStore(tmp_path, ram_threshold_bytes=0) as store:
        out = store.allocate("s1:x/out", (4, 3, 2), "float32")
        desc = store.descriptor("s1:x/out")
        assert desc == FieldDescriptor(desc.path, "float32", (4, 3, 2), 0)
        assert pickle.loads(pickle.dumps(desc)) == desc
        assert len(pickle.dumps(desc)) < 1024
        # A worker writes through its own "r+" map; the owner sees the data.
        writer = open_field(desc, "r+")
        assert isinstance(writer, np.memmap) and writer.mode == "r+"
        writer[1:3] = 7.0
        writer.flush()
        del writer
        assert out[1:3].tolist() == np.full((2, 3, 2), 7.0).tolist()
        frozen = store.freeze("s1:x/out")
        reader = open_field(store.descriptor("s1:x/out"), "r")
        with pytest.raises(ValueError):
            reader[0, 0, 0] = 1.0
        np.testing.assert_array_equal(reader, frozen)
        # A whole frozen map (e.g. a state field) has a descriptor of its own; views do not.
        assert descriptor_of(frozen) == store.descriptor("s1:x/out")
        with pytest.raises(ValueError, match="whole"):
            descriptor_of(frozen[1:])
        with pytest.raises(ValueError, match="whole"):
            descriptor_of(np.zeros((2, 2)))
        # RAM fields have no descriptor.
        store_ram = FieldStore(None, ram_threshold_bytes=0)
        store_ram.allocate("s1:x/ram", (2, 2))
        with pytest.raises(ValueError, match="RAM"):
            store_ram.descriptor("s1:x/ram")
        with pytest.raises(KeyError):
            store.descriptor("s9:x/missing")
        with pytest.raises(ValueError, match="mode"):
            open_field(desc, "w+")  # type: ignore[arg-type]
        with pytest.raises(TypeError):
            open_field(str(desc.path), "r")  # type: ignore[arg-type]
        with pytest.raises(ValueError):
            FieldDescriptor(Path("x"), "float64", (2,), -1)
        # The store itself never crosses a process boundary.
        with pytest.raises(TypeError, match="pickled"):
            pickle.dumps(store)


def test_fieldstore_rejects_non_native_byte_order(tmp_path):
    # A descriptor carries the dtype name, which has no byte order: a non-native field
    # would be described (and reopened by a worker) as native, so it is refused.
    swapped = np.dtype("float64").newbyteorder()
    with pytest.raises(ValueError, match="byte order"):
        FieldDescriptor(Path("x"), swapped.str, (2,))
    with FieldStore(tmp_path, ram_threshold_bytes=0) as store:
        with pytest.raises(ValueError, match="byte order"):
            store.allocate("s1:x/swapped", (2, 2), swapped.str)
        assert "s1:x/swapped" not in store
        foreign = np.memmap(tmp_path / "foreign.bin", dtype=swapped, mode="w+", shape=(3,))
        with pytest.raises(ValueError, match="byte order"):
            descriptor_of(foreign)
        del foreign
    # The native byte order spelled out explicitly is the same dtype.
    native = np.dtype("float64").newbyteorder("=")
    assert FieldDescriptor(Path("x"), native.str, (2,)).dtype == "float64"


def test_fieldstore_name_of(tmp_path):
    # Spec §4.3: release goes by store name; name_of finds the name of a field handed on
    # (for example a state field a later step deletes).
    with FieldStore(tmp_path, ram_threshold_bytes=64) as store:
        mapped = store.allocate("s1:x/big", (4, 4))
        assert store.name_of(mapped) == "s1:x/big"
        frozen = store.freeze("s1:x/big")
        store.allocate("s1:x/small", (2, 2))  # 32 bytes: in RAM
        ram = store.freeze("s1:x/small")
        assert not store.is_mapped("s1:x/small")
        assert store.name_of(frozen) == "s1:x/big"
        assert store.name_of(frozen[1:3]) == "s1:x/big"  # views of a field
        assert store.name_of(ram) == "s1:x/small"
        assert store.name_of(ram[:1]) == "s1:x/small"
        # A map of the same file opened elsewhere (a worker's, say) names the field too.
        reopened = open_field(store.descriptor("s1:x/big"), "r")
        assert store.name_of(reopened) == "s1:x/big"
        del reopened
        store.release(store.name_of(frozen))
        assert "s1:x/big" not in store
        assert list(tmp_path.iterdir()) == []
        with pytest.raises(KeyError, match="not a field of this store"):
            store.name_of(frozen)  # released
        with pytest.raises(KeyError, match="not a field of this store"):
            store.name_of(np.zeros((2, 2)))
        with pytest.raises(KeyError, match="not a field of this store"):
            store.name_of(ram.copy())
        with pytest.raises(TypeError):
            store.name_of([1.0, 2.0])  # type: ignore[arg-type]
