"""Reading input images and writing output images, and the sRGB transfer functions.

Integer images (8 or 16 bits per sample) are taken to be sRGB encoded and are decoded to
linear RGB; floating-point images (TIFF) are taken to be scene-linear already. Grey images
are copied to three channels, alpha channels are dropped, and non-finite samples are an
error. Palette images are expanded through their colour map, MinIsWhite TIFFs are inverted,
and the EXIF or TIFF Orientation tag is applied, so the array is the image as a viewer shows
it. Colour models other than grey, palette and RGB (CMYK, YCbCr stored as such, CIELab, raw
sensor data) are rejected rather than misread. Outputs are 16-bit sRGB PNG.

PNG files are decoded and encoded with imagecodecs (libpng), because Pillow reads 16-bit
RGB PNG as 8 bits per sample without warning and cannot write it. TIFF files are read with
tifffile, everything else (JPEG and the other formats Pillow reads) with Pillow.
"""

from collections.abc import Callable, Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from numbers import Integral
from pathlib import Path
from typing import Literal

import imagecodecs
import numpy as np
import tifffile
from PIL import Image, UnidentifiedImageError

from darkroom_atlas._frozen import readonly_array
from darkroom_atlas.backend import xp_of

__all__ = [
    "EIGHT_BIT_WARNING",
    "ImageInput",
    "read_image",
    "srgb_decode",
    "srgb_encode",
    "write_png16",
]

EIGHT_BIT_WARNING = (
    "8-bit display-referred input: highlight latitude of negative film will not be exercised"
)

_PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"
_TIFF_SIGNATURES = (b"II*\x00", b"MM\x00*", b"II+\x00", b"MM\x00+")
_ORIENTATION_TAG = 274  # EXIF and TIFF Orientation
_UNREADABLE_EXIF = "unreadable EXIF data: the Orientation tag, if any, is ignored"

# How to turn the stored array into the displayed one for each Orientation value (EXIF 2.3,
# TIFF 6.0): 2-4 mirror or turn by 180 degrees, 5-8 also swap rows and columns.
_ORIENTATIONS: dict[int, Callable[[np.ndarray], np.ndarray]] = {
    1: lambda a: a,
    2: lambda a: a[:, ::-1],
    3: lambda a: a[::-1, ::-1],
    4: lambda a: a[::-1],
    5: lambda a: np.swapaxes(a, 0, 1),
    6: lambda a: np.rot90(a, -1),
    7: lambda a: np.swapaxes(a[::-1, ::-1], 0, 1),
    8: lambda a: np.rot90(a, 1),
}

# Pillow modes read as they are, and those converted first (palette modes are handled apart).
_PILLOW_MODES = frozenset({"1", "L", "LA", "RGB", "RGBA", "I;16", "I;16L", "I;16B", "F"})
_PILLOW_CONVERSIONS = {"La": "LA", "RGBa": "RGBA", "RGBX": "RGB", "PA": "RGBA"}

# TIFF photometric interpretations that tifffile returns as grey or RGB samples. YCbCr is
# converted to RGB by tifffile only when it decodes JPEG-compressed data.
_PHOTOMETRIC = tifffile.PHOTOMETRIC
_GREY_OR_RGB = frozenset({_PHOTOMETRIC.MINISBLACK, _PHOTOMETRIC.RGB})

# IEC 61966-2-1 piecewise transfer function.
_DECODE_THRESHOLD = 0.04045
_ENCODE_THRESHOLD = 0.0031308
_LINEAR_SLOPE = 12.92
_OFFSET = 0.055
_GAMMA = 2.4


# --------------------------------------------------------------------------------------
# sRGB transfer functions (Array API)
# --------------------------------------------------------------------------------------


def _as_float(x):
    try:
        xp = xp_of(x)
        was_scalar = False
    except TypeError:
        was_scalar = np.ndim(x) == 0
        x = np.asarray(x, dtype=np.float64)
        xp = xp_of(x)
    if not xp.isdtype(x.dtype, "real floating"):
        x = xp.astype(x, xp.float64)
    return xp, x, was_scalar


def srgb_decode(encoded):
    """Decode sRGB-encoded values to linear (IEC 61966-2-1).

    ``c <= 0.04045`` maps to ``c / 12.92``, larger values to ``((c + 0.055) / 1.055) ** 2.4``.
    Negative values are mirrored (``f(-c) == -f(c)``). Scalars give a ``float``.
    """
    xp, c, was_scalar = _as_float(encoded)
    a = xp.abs(c)
    linear = xp.where(
        a <= _DECODE_THRESHOLD, a / _LINEAR_SLOPE, ((a + _OFFSET) / (1.0 + _OFFSET)) ** _GAMMA
    )
    result = xp.where(c < 0.0, -linear, linear)
    return float(result) if was_scalar else result


def srgb_encode(linear):
    """Encode linear values with the sRGB transfer function (IEC 61966-2-1).

    ``l <= 0.0031308`` maps to ``12.92 * l``, larger values to ``1.055 * l ** (1/2.4) - 0.055``.
    Negative values are mirrored. Scalars give a ``float``.
    """
    xp, x, was_scalar = _as_float(linear)
    a = xp.abs(x)
    encoded = xp.where(
        a <= _ENCODE_THRESHOLD,
        a * _LINEAR_SLOPE,
        (1.0 + _OFFSET) * a ** (1.0 / _GAMMA) - _OFFSET,
    )
    result = xp.where(x < 0.0, -encoded, encoded)
    return float(result) if was_scalar else result


# --------------------------------------------------------------------------------------
# Reading
# --------------------------------------------------------------------------------------


@dataclass(frozen=True, eq=False)
class ImageInput:
    """An input image as linear RGB.

    ``linear_rgb`` is a read-only ``float64`` array of shape ``(H, W, 3)``. ``bit_depth`` is
    the number of bits per sample in the file as decoded, ``encoding`` says whether the file
    was sRGB encoded (integer files) or linear (floating-point files), and ``warnings``
    collects everything the reader had to tell the user.
    """

    linear_rgb: np.ndarray
    bit_depth: int
    encoding: Literal["srgb", "linear"]
    warnings: tuple[str, ...]

    def __post_init__(self) -> None:
        rgb = readonly_array(self.linear_rgb, dtype=np.float64)
        if rgb.ndim != 3 or rgb.shape[2] != 3:
            raise ValueError(f"linear_rgb must have shape (H, W, 3), got {rgb.shape}")
        if self.encoding not in ("srgb", "linear"):
            raise ValueError(f"encoding must be 'srgb' or 'linear', got {self.encoding!r}")
        object.__setattr__(self, "linear_rgb", rgb)
        object.__setattr__(self, "bit_depth", int(self.bit_depth))
        object.__setattr__(self, "warnings", tuple(self.warnings))


def read_image(path: str | Path) -> ImageInput:
    """Read an image file as linear RGB.

    * 8- and 16-bit integer images are sRGB decoded after scaling by 255 or 65535; 8-bit
      input adds :data:`EIGHT_BIT_WARNING` to the warnings.
    * Floating-point images (TIFF) are scene-linear and kept as they are, except that
      negative samples are clipped to 0 (with a warning); NaN or inf raises ``ValueError``.
    * Grey images are copied to three channels; an alpha channel is dropped (with a warning
      when it is not fully opaque).
    * Palette images are expanded through their colour map (8-bit when the map holds 8-bit
      colours), MinIsWhite TIFFs are inverted, and the Orientation tag (EXIF, or TIFF tag
      274) is applied.
    * CMYK, YCbCr (other than JPEG-compressed TIFF, which is decoded to RGB), CIELab and the
      other colour models raise ``ValueError``.

    Problems with the file's content, including an empty or unrecognised file, raise a
    one-line ``ValueError`` naming the file; a missing file raises ``FileNotFoundError``.
    """
    path = Path(path)
    with path.open("rb") as handle:
        head = handle.read(8)
    warnings: list[str] = []
    raw = _decode(path, head, warnings)
    colour = _rgb_channels(path, raw, warnings)

    kind = colour.dtype.kind
    if kind == "b":
        colour = colour.astype(np.uint8) * np.uint8(255)
        kind = "u"
    if kind == "u" and colour.dtype.itemsize in (1, 2):
        bit_depth = 8 * colour.dtype.itemsize
        if bit_depth == 8:
            warnings.insert(0, EIGHT_BIT_WARNING)
        full_scale = float(np.iinfo(colour.dtype).max)
        linear = srgb_decode(colour.astype(np.float64) / full_scale)
        encoding: Literal["srgb", "linear"] = "srgb"
    elif kind == "f":
        bit_depth = 8 * colour.dtype.itemsize
        linear = colour.astype(np.float64)
        bad = ~np.isfinite(linear)
        if bad.any():
            raise ValueError(
                f"{path}: {int(bad.sum())} non-finite sample(s) (NaN or inf) in the image; "
                "replace them before running"
            )
        negative = linear < 0.0
        if negative.any():
            warnings.append(
                f"{int(negative.sum())} negative linear sample(s) clipped to 0 "
                "(light cannot be negative)"
            )
            linear = np.maximum(linear, 0.0)
        encoding = "linear"
    else:
        raise ValueError(
            f"{path}: unsupported sample type {colour.dtype}; expected 8- or 16-bit unsigned "
            "integers or floating point"
        )
    return ImageInput(linear, bit_depth, encoding, tuple(warnings))


def _decode(path: Path, head: bytes, warnings: list[str]) -> np.ndarray:
    """Decode ``path`` to grey, grey + alpha, RGB or RGBA samples, as a viewer shows it."""
    if not head:
        raise ValueError(f"{path}: the file is empty")
    if head.startswith(_PNG_SIGNATURE):
        array, orientation = _decode_png(path, warnings)
    elif head[:4] in _TIFF_SIGNATURES:
        array, orientation = _decode_tiff(path)
    else:
        array, orientation = _decode_pillow(path, warnings)
    return _orient(array, orientation, warnings)


@contextmanager
def _decoding(path: Path) -> Iterator[None]:
    """Turn any failure of a decoder into a one-line ``ValueError`` naming ``path``.

    Decoders raise many exception types on malformed files (``struct.error``, ``EOFError``,
    ``SyntaxError``, ...), so everything except a missing file and lack of memory is caught.
    """
    try:
        yield
    except (FileNotFoundError, MemoryError):
        raise
    except UnidentifiedImageError:
        raise ValueError(
            f"{path}: not a recognised image file (expected PNG, TIFF, JPEG or another format "
            "that Pillow reads)"
        ) from None
    except Exception as exc:
        message = " ".join(str(exc).split()) or type(exc).__name__
        raise ValueError(f"{path}: cannot decode image: {message}") from exc


def _decode_png(path: Path, warnings: list[str]) -> tuple[np.ndarray, object]:
    data = path.read_bytes()
    with _decoding(path):
        array = np.asarray(imagecodecs.png_decode(data))
    try:
        orientation = _png_orientation(data)
    except Exception:  # noqa: BLE001 - a broken eXIf chunk must not stop the image
        warnings.append(_UNREADABLE_EXIF)
        orientation = None
    return array, orientation


def _png_orientation(data: bytes) -> object:
    """The Orientation value of the PNG's ``eXIf`` chunk, or ``None`` if it has none."""
    pos = len(_PNG_SIGNATURE)
    while pos + 8 <= len(data):
        length = int.from_bytes(data[pos : pos + 4], "big")
        tag = data[pos + 4 : pos + 8]
        if tag == b"eXIf":
            exif = Image.Exif()
            exif.load(data[pos + 8 : pos + 8 + length])
            return exif.get(_ORIENTATION_TAG)
        if tag == b"IEND":
            break
        pos += 12 + length
    return None


def _decode_tiff(path: Path) -> tuple[np.ndarray, object]:
    with _decoding(path), tifffile.TiffFile(path) as tif:
        page = tif.pages.first
        array = page.asarray()
        photometric = page.photometric
        compression = page.compression
        bits = page.bitspersample
        colormap = page.colormap
        orientation = page.tags.valueof(_ORIENTATION_TAG)
        if array.ndim == 3 and page.planarconfig == tifffile.PLANARCONFIG.SEPARATE:
            array = np.moveaxis(array, 0, -1)

    if photometric in _GREY_OR_RGB:
        return array, orientation
    if photometric == _PHOTOMETRIC.YCBCR and compression == tifffile.COMPRESSION.JPEG:
        return array, orientation
    if photometric == _PHOTOMETRIC.MINISWHITE:
        return _invert_min_is_white(path, array, bits), orientation
    if photometric == _PHOTOMETRIC.PALETTE:
        return _expand_palette(path, array, colormap), orientation
    try:
        name = _PHOTOMETRIC(photometric).name
    except ValueError:
        name = str(int(photometric))
    if name == "SEPARATED":
        name = "SEPARATED (CMYK)"
    raise ValueError(
        f"{path}: unsupported TIFF photometric interpretation {name}; convert the image to "
        "RGB or greyscale first"
    )


def _invert_min_is_white(path: Path, array: np.ndarray, bits: int) -> np.ndarray:
    """MinIsWhite grey (0 is white) as MinIsBlack; extra samples such as alpha are kept."""
    if array.dtype == np.bool_:
        flipped = np.logical_not(array)
    elif array.dtype.kind == "u":
        flipped = ((1 << int(bits)) - 1 - array.astype(np.int64)).astype(array.dtype)
    else:
        raise ValueError(
            f"{path}: MinIsWhite TIFF with {array.dtype} samples is not supported; expected "
            "unsigned integer or bilevel samples"
        )
    if array.ndim == 3:  # grey plus extra samples: invert only the grey one
        flipped = np.concatenate([flipped[..., :1], array[..., 1:]], axis=-1)
    return flipped


def _expand_palette(path: Path, array: np.ndarray, colormap: np.ndarray | None) -> np.ndarray:
    """Palette indices ``(H, W)`` as RGB through the TIFF colour map ``(3, N)``."""
    if colormap is None:
        raise ValueError(f"{path}: palette TIFF without a colour map")
    if array.ndim != 2 or array.dtype.kind != "u":
        raise ValueError(
            f"{path}: unsupported palette TIFF: expected one unsigned-integer index per pixel, "
            f"got samples of shape {array.shape[2:]} and type {array.dtype}"
        )
    colours = np.asarray(colormap, dtype=np.uint16)
    if colours.ndim != 2 or colours.shape[0] != 3 or int(array.max()) >= colours.shape[1]:
        raise ValueError(
            f"{path}: palette index {int(array.max())} lies outside the colour map of shape "
            f"{colours.shape}"
        )
    # Colour maps hold 16-bit values. Writers store 8-bit colours as v * 257 (full scale) or
    # v * 256 (Pillow), and some store the raw 8-bit values, which libtiff also accepts when
    # no entry exceeds 255; those maps are read back as 8-bit colours.
    if colours.max() <= 255:
        colours = colours.astype(np.uint8)
    elif not np.any(colours % 257):
        colours = (colours // 257).astype(np.uint8)
    elif not np.any(colours % 256):
        colours = (colours // 256).astype(np.uint8)
    return np.moveaxis(colours[:, array], 0, -1)


def _decode_pillow(path: Path, warnings: list[str]) -> tuple[np.ndarray, object]:
    with _decoding(path), Image.open(path) as image:
        mode, image_format = image.mode, image.format
        try:
            orientation = image.getexif().get(_ORIENTATION_TAG)
        except Exception:  # noqa: BLE001 - broken EXIF data must not stop the image
            warnings.append(_UNREADABLE_EXIF)
            orientation = None
        if mode == "P":
            has_alpha = "transparency" in image.info
            array = np.asarray(image.convert("RGBA" if has_alpha else "RGB"))
        elif mode in _PILLOW_CONVERSIONS:
            array = np.asarray(image.convert(_PILLOW_CONVERSIONS[mode]))
        elif mode in _PILLOW_MODES:
            array = np.asarray(image)
        else:
            array = None
    if array is None:
        raise ValueError(
            f"{path}: unsupported colour mode {mode} in this {image_format} image; convert the "
            "image to RGB or greyscale first"
        )
    return array, orientation


def _orient(array: np.ndarray, orientation: object, warnings: list[str]) -> np.ndarray:
    """Apply an Orientation tag value (``None``: no tag) to the decoded array."""
    if orientation is None:
        return array
    if isinstance(orientation, Integral) and int(orientation) in _ORIENTATIONS:
        return _ORIENTATIONS[int(orientation)](array)
    warnings.append(f"unknown Orientation tag value {orientation!r} ignored")
    return array


def _rgb_channels(path: Path, raw: np.ndarray, warnings: list[str]) -> np.ndarray:
    """Return ``raw`` as ``(H, W, 3)``: grey is copied to three channels, alpha dropped."""
    if raw.ndim == 2:
        raw = raw[..., np.newaxis]
    if raw.ndim != 3 or raw.shape[0] < 1 or raw.shape[1] < 1:
        raise ValueError(f"{path}: unsupported image shape {raw.shape}; expected (H, W[, C])")
    channels = raw.shape[2]
    if channels not in (1, 2, 3, 4):
        raise ValueError(f"{path}: unsupported number of channels {channels}")
    if channels in (2, 4):
        alpha = raw[..., -1]
        opaque = np.iinfo(raw.dtype).max if raw.dtype.kind in "ui" else 1
        if raw.dtype.kind == "b" or np.any(alpha != opaque):
            warnings.append("alpha channel dropped: transparency is ignored")
        raw = raw[..., :-1]
    if raw.shape[2] == 1:
        raw = np.repeat(raw, 3, axis=2)
    return raw


# --------------------------------------------------------------------------------------
# Writing
# --------------------------------------------------------------------------------------


def write_png16(path: str | Path, linear_rgb) -> None:
    """Write linear RGB ``(H, W, 3)`` as a 16-bit sRGB PNG.

    Values are clipped to ``[0, 1]``, sRGB encoded and rounded to 16 bits. The encoding is
    deterministic: the same input always gives the same bytes. Non-finite values raise
    ``ValueError``.
    """
    rgb = np.asarray(linear_rgb, dtype=np.float64)
    if rgb.ndim != 3 or rgb.shape[2] != 3 or rgb.shape[0] < 1 or rgb.shape[1] < 1:
        raise ValueError(f"linear_rgb must have shape (H, W, 3), got {rgb.shape}")
    if not np.all(np.isfinite(rgb)):
        raise ValueError("linear_rgb must be finite (found NaN or inf)")
    encoded = srgb_encode(np.clip(rgb, 0.0, 1.0))
    quantised = np.round(encoded * 65535.0).astype(np.uint16)
    Path(path).write_bytes(imagecodecs.png_encode(quantised))
