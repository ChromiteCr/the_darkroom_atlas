import struct
import zlib
from pathlib import Path

import imageio.v3 as iio
import numpy as np
import pytest
import tifffile
from PIL import Image, ImageOps

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


# ----------------------------------------------------------------------------------------
# Photometric interpretations, colour modes, orientation and unreadable files
# ----------------------------------------------------------------------------------------

_PALETTE_RGB = np.array([[255, 0, 0], [0, 255, 0], [0, 0, 255], [255, 255, 255]], np.uint8)
_INDICES = np.array([[0, 1], [2, 3]], np.uint8)


def _expected(encoded_uint8: np.ndarray) -> np.ndarray:
    return srgb_decode(encoded_uint8 / 255.0)


@pytest.mark.parametrize("scale", [257, 256, 1], ids=["16-bit", "8-bit-times-256", "raw-8-bit"])
def test_read_palette_tiff(tmp_path, scale):
    # TIFF colour maps are 16-bit; writers fill them as v * 257 (the standard), v * 256
    # (Pillow) or, against the standard, with raw 8-bit values (which libtiff also accepts).
    colormap = np.zeros((3, 256), np.uint16)
    colormap[:, :4] = _PALETTE_RGB.T.astype(np.uint16) * scale
    path = tmp_path / "palette.tif"
    tifffile.imwrite(path, _INDICES, photometric="palette", colormap=colormap)
    image = read_image(path)
    np.testing.assert_allclose(image.linear_rgb, _expected(_PALETTE_RGB[_INDICES]), atol=1e-15)
    assert image.bit_depth == 8
    assert not any("alpha" in w for w in image.warnings)


def test_read_palette_tiff_written_by_pillow(tmp_path):
    palette = Image.new("P", (2, 2))
    palette.putpalette(_PALETTE_RGB.ravel().tolist())
    palette.putdata(_INDICES.ravel().tolist())
    path = tmp_path / "pillow_palette.tif"
    palette.save(path)
    expected = np.asarray(Image.open(path).convert("RGB"))
    np.testing.assert_allclose(read_image(path).linear_rgb, _expected(expected), atol=1e-15)


def test_read_palette_tiff_16bit_colours(tmp_path):
    colormap = np.zeros((3, 256), np.uint16)
    colormap[:, 1] = (1000, 30000, 65535)
    path = tmp_path / "palette16.tif"
    tifffile.imwrite(path, np.array([[0, 1]], np.uint8), photometric="palette", colormap=colormap)
    image = read_image(path)
    assert image.bit_depth == 16
    expected = srgb_decode(np.array([[[0, 0, 0], [1000, 30000, 65535]]]) / 65535.0)
    np.testing.assert_allclose(image.linear_rgb, expected, atol=1e-15)


def test_read_min_is_white_tiff(tmp_path):
    grey = np.array([[0, 255, 100]], np.uint8)
    path = tmp_path / "miniswhite.tif"
    tifffile.imwrite(path, grey, photometric="miniswhite")
    image = read_image(path)
    np.testing.assert_allclose(image.linear_rgb[..., 0], _expected(255 - grey), atol=1e-15)
    np.testing.assert_array_equal(image.linear_rgb[..., 0], image.linear_rgb[..., 2])


def test_read_min_is_white_bilevel_tiff(tmp_path):
    path = tmp_path / "fax.tif"
    tifffile.imwrite(path, np.array([[True, False, True]]), photometric="miniswhite")
    expected = np.asarray(Image.open(path).convert("L"))  # 1 is black in MinIsWhite
    np.testing.assert_allclose(read_image(path).linear_rgb[..., 1], _expected(expected), atol=0)


def test_read_jpeg_compressed_ycbcr_tiff(tmp_path):
    colour = np.full((16, 16, 3), (200, 30, 100), np.uint8)
    path = tmp_path / "ycbcr.tif"
    tifffile.imwrite(path, colour, photometric="rgb", compression="jpeg")
    assert tifffile.TiffFile(path).pages.first.photometric == tifffile.PHOTOMETRIC.YCBCR
    np.testing.assert_allclose(read_image(path).linear_rgb, _expected(colour), atol=0.01)


def _cmyk_tiff(path):
    tifffile.imwrite(path, np.full((2, 2, 4), (50, 0, 0, 200), np.uint8), photometric="separated")


def _ycbcr_tiff(path):
    Image.new("RGB", (4, 4), (200, 30, 100)).convert("YCbCr").save(path)


def _lab_tiff(path):
    Image.new("RGB", (4, 4), (200, 30, 100)).convert("LAB").save(path)


def _cmyk_jpeg(path):
    Image.new("RGB", (8, 8), (200, 30, 100)).convert("CMYK").save(path, format="JPEG")


@pytest.mark.parametrize(
    ("name", "make", "fragment"),
    [
        ("cmyk.tif", _cmyk_tiff, "CMYK"),
        ("ycbcr.tif", _ycbcr_tiff, "YCBCR"),
        ("lab.tif", _lab_tiff, "CIELAB"),
        ("cmyk.jpg", _cmyk_jpeg, "CMYK"),
    ],
)
def test_read_rejects_unsupported_colour_models(tmp_path, name, make, fragment):
    path = tmp_path / name
    make(path)
    with pytest.raises(ValueError, match=fragment) as caught:
        read_image(path)
    message = str(caught.value)
    assert name in message
    assert "alpha" not in message
    assert len(message.splitlines()) == 1


def test_read_jpeg_and_palette_formats(tmp_path):
    colour = np.random.default_rng(5).integers(0, 256, (6, 5, 3), dtype=np.uint8)
    jpeg = tmp_path / "photo.jpg"
    Image.fromarray(colour).save(jpeg, quality=95)
    image = read_image(jpeg)
    np.testing.assert_allclose(
        image.linear_rgb, _expected(np.asarray(Image.open(jpeg))), atol=1e-15
    )
    assert EIGHT_BIT_WARNING in image.warnings

    palette = Image.new("P", (2, 2))
    palette.putpalette(_PALETTE_RGB.ravel().tolist())
    palette.putdata(_INDICES.ravel().tolist())
    gif = tmp_path / "palette.gif"
    palette.save(gif)
    np.testing.assert_allclose(
        read_image(gif).linear_rgb, _expected(_PALETTE_RGB[_INDICES]), atol=1e-15
    )


@pytest.mark.parametrize(
    ("name", "content", "fragment"),
    [
        ("empty.png", b"", "empty"),
        ("text.png", b"this is not an image\n", "not a recognised image"),
        ("tiny.png", b"ab", "not a recognised image"),
        ("fake.exr", b"abc", "not a recognised image"),
        ("broken.png", b"\x89PNG\r\n\x1a\n\x00\x00", "cannot decode"),
        ("broken.tif", b"II*\x00\xff\xff\xff\xff", "cannot decode"),
    ],
)
def test_read_unreadable_files_one_line(tmp_path, name, content, fragment):
    path = tmp_path / name
    path.write_bytes(content)
    with pytest.raises(ValueError, match=fragment) as caught:
        read_image(path)
    message = str(caught.value)
    assert message.startswith(str(path))
    assert len(message.splitlines()) == 1
    assert "pip" not in message


def _asymmetric_rgb() -> np.ndarray:
    """A 6 x 10 image whose eight orientations all differ."""
    rows, cols = np.mgrid[0:6, 0:10]
    return np.stack([rows * 40, cols * 25, (rows * 10 + cols) % 256], axis=-1).astype(np.uint8)


def _displayed(path) -> np.ndarray:
    """What a viewer shows: Pillow's decode with the EXIF orientation applied."""
    return np.asarray(ImageOps.exif_transpose(Image.open(path)).convert("RGB"))


@pytest.mark.parametrize("orientation", range(1, 9))
@pytest.mark.parametrize("suffix", [".jpg", ".tif", ".png"])
def test_read_applies_orientation(tmp_path, orientation, suffix):
    stored = _asymmetric_rgb()
    path = tmp_path / f"oriented{suffix}"
    if suffix == ".tif":
        tifffile.imwrite(path, stored, photometric="rgb", extratags=[(274, "H", 1, orientation)])
    else:
        exif = Image.Exif()
        exif[0x0112] = orientation
        Image.fromarray(stored).save(path, exif=exif.tobytes(), quality=100)
    displayed = _displayed(path)
    if orientation in (5, 6, 7, 8):
        assert displayed.shape[:2] == (10, 6)
    image = read_image(path)
    np.testing.assert_allclose(image.linear_rgb, _expected(displayed), atol=1e-15)
    assert not any("orientation" in w.lower() for w in image.warnings)


def test_read_png_with_broken_exif_warns(tmp_path):
    data = SAMPLE.read_bytes()
    garbage = b"not exif at all"
    chunk = (
        struct.pack(">I", len(garbage))
        + b"eXIf"
        + garbage
        + struct.pack(">I", zlib.crc32(b"eXIf" + garbage) & 0xFFFFFFFF)
    )
    iend = data.rindex(b"IEND") - 4
    path = tmp_path / "broken_exif.png"
    path.write_bytes(data[:iend] + chunk + data[iend:])
    image = read_image(path)
    np.testing.assert_array_equal(image.linear_rgb, read_image(SAMPLE).linear_rgb)
    assert any("EXIF" in w for w in image.warnings)
