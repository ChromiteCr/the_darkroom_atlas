import struct
import zlib
from pathlib import Path

import imageio.v3 as iio
import numpy as np
import pytest
import tifffile

from darkroom_atlas.scene.io import ImageInput, read_image, srgb_decode, srgb_encode, write_png16

SAMPLE = Path(__file__).resolve().parents[1] / "assets" / "sample.png"
EIGHT_BIT_WARNING = (
    "8-bit display-referred input: highlight latitude of negative film will not be exercised"
)


def _png16_rgb(path: Path, image: np.ndarray) -> None:
    """Write a 16-bit RGB PNG byte by byte, independently of any imaging library."""

    def chunk(tag: bytes, data: bytes) -> bytes:
        crc = zlib.crc32(tag + data) & 0xFFFFFFFF
        return struct.pack(">I", len(data)) + tag + data + struct.pack(">I", crc)

    h, w, _ = image.shape
    big_endian = image.astype(">u2")
    raw = b"".join(b"\x00" + big_endian[r].tobytes() for r in range(h))
    ihdr = struct.pack(">IIBBBBB", w, h, 16, 2, 0, 0, 0)
    path.write_bytes(
        b"\x89PNG\r\n\x1a\n"
        + chunk(b"IHDR", ihdr)
        + chunk(b"IDAT", zlib.compress(raw))
        + chunk(b"IEND", b"")
    )


def test_srgb_roundtrip():
    x = np.linspace(0.0, 1.0, 256)
    np.testing.assert_allclose(srgb_encode(srgb_decode(x)), x, rtol=0, atol=1e-12)


def test_srgb_known_values():
    assert srgb_decode(0.5) == pytest.approx(0.21404, abs=1e-5)
    assert srgb_decode(0.04045) == pytest.approx(0.0031308, abs=1e-7)
    assert srgb_decode(0.0) == 0.0
    assert srgb_decode(1.0) == pytest.approx(1.0, abs=1e-15)
    assert srgb_encode(1.0) == pytest.approx(1.0, abs=1e-15)


def test_read_sample_png():
    image = read_image(SAMPLE)
    assert isinstance(image, ImageInput)
    rgb = image.linear_rgb
    assert rgb.shape == (64, 96, 3)
    assert rgb.dtype == np.float64
    assert rgb.min() >= 0.0
    assert rgb.max() <= 1.0
    assert image.bit_depth == 8
    assert image.encoding == "srgb"
    assert EIGHT_BIT_WARNING in image.warnings
    assert not rgb.flags.writeable
    encoded = iio.imread(SAMPLE)
    np.testing.assert_allclose(rgb, srgb_decode(encoded / 255.0), rtol=0, atol=1e-15)


def test_read_grayscale(tmp_path):
    grey = np.tile(np.arange(0, 250, 10, dtype=np.uint8), (4, 1))
    path = tmp_path / "grey.png"
    iio.imwrite(path, grey)
    rgb = read_image(path).linear_rgb
    assert rgb.shape == (*grey.shape, 3)
    np.testing.assert_array_equal(rgb[..., 0], rgb[..., 1])
    np.testing.assert_array_equal(rgb[..., 0], rgb[..., 2])
    np.testing.assert_allclose(rgb[..., 0], srgb_decode(grey / 255.0), rtol=0, atol=1e-15)


def test_read_rgba_drops_alpha(tmp_path):
    rng = np.random.default_rng(2)
    rgba = rng.integers(0, 256, (6, 7, 4), dtype=np.uint8)
    rgba[..., 3] = 0  # fully transparent: the colour must still be read as is
    path = tmp_path / "rgba.png"
    iio.imwrite(path, rgba)
    image = read_image(path)
    assert image.linear_rgb.shape == (6, 7, 3)
    np.testing.assert_allclose(
        image.linear_rgb, srgb_decode(rgba[..., :3] / 255.0), rtol=0, atol=1e-15
    )
    assert any("alpha" in w for w in image.warnings)


def test_read_png16(tmp_path):
    values = np.array([[[65535, 0, 1], [12345, 40000, 257]]], dtype=np.uint16)
    path = tmp_path / "rgb16.png"
    _png16_rgb(path, values)
    image = read_image(path)
    assert image.bit_depth == 16
    assert image.encoding == "srgb"
    assert EIGHT_BIT_WARNING not in image.warnings
    assert image.linear_rgb[0, 0, 0] == 1.0
    assert image.linear_rgb[0, 0, 1] == 0.0
    # Full 16-bit precision: no library may silently drop the low byte.
    np.testing.assert_allclose(image.linear_rgb, srgb_decode(values / 65535.0), rtol=0, atol=1e-15)
    assert image.linear_rgb[0, 0, 2] > 0.0


def test_read_float_tiff_linear(tmp_path):
    rng = np.random.default_rng(3)
    values = (rng.random((5, 6, 3)) * 4.0).astype(np.float32)  # scene-linear, above 1 allowed
    path = tmp_path / "linear.tif"
    tifffile.imwrite(path, values, photometric="rgb")
    image = read_image(path)
    assert image.encoding == "linear"
    assert image.bit_depth == 32
    assert image.linear_rgb.dtype == np.float64
    np.testing.assert_array_equal(image.linear_rgb, values.astype(np.float64))
    assert EIGHT_BIT_WARNING not in image.warnings


def test_read_rejects_nonfinite(tmp_path):
    values = np.ones((4, 4, 3), dtype=np.float32)
    values[1, 2, 0] = np.nan
    path = tmp_path / "broken_scene.tif"
    tifffile.imwrite(path, values, photometric="rgb")
    with pytest.raises(ValueError, match="broken_scene.tif"):
        read_image(path)
    values[1, 2, 0] = np.inf
    tifffile.imwrite(path, values, photometric="rgb")
    with pytest.raises(ValueError, match="broken_scene.tif"):
        read_image(path)


def test_write_png16_roundtrip(tmp_path):
    rng = np.random.default_rng(4)
    linear = rng.random((9, 11, 3))
    linear[0, 0] = (0.0, 1.0, 0.5)
    path = tmp_path / "out.png"
    write_png16(path, linear)
    image = read_image(path)
    assert image.bit_depth == 16
    assert image.linear_rgb.shape == linear.shape
    assert np.max(np.abs(image.linear_rgb - linear)) < 1e-4


def test_write_png16_clips_and_is_deterministic(tmp_path):
    linear = np.array([[[-0.5, 2.0, 0.25]]])
    first, second = tmp_path / "a.png", tmp_path / "b.png"
    write_png16(first, linear)
    write_png16(second, linear)
    assert first.read_bytes() == second.read_bytes()
    back = read_image(first).linear_rgb
    assert back[0, 0, 0] == 0.0
    assert back[0, 0, 1] == 1.0
    with pytest.raises(ValueError):
        write_png16(tmp_path / "c.png", np.full((2, 2, 3), np.nan))
