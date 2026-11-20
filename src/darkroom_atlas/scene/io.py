"""Reading input images and writing output images, and the sRGB transfer functions.

Integer images (8 or 16 bits per sample) are taken to be sRGB encoded and are decoded to
linear RGB; floating-point images (TIFF) are taken to be scene-linear already. Grey images
are copied to three channels, alpha channels are dropped, and non-finite samples are an
error. Outputs are 16-bit sRGB PNG.

PNG files are decoded and encoded with imagecodecs (libpng), because Pillow reads 16-bit
RGB PNG as 8 bits per sample without warning and cannot write it. TIFF files are read with
tifffile, everything else through imageio.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import imagecodecs
import imageio.v3 as iio
import numpy as np
import tifffile

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

    Problems with the file's content raise ``ValueError`` naming the file; a missing file
    raises ``FileNotFoundError``.
    """
    path = Path(path)
    with path.open("rb") as handle:
        head = handle.read(8)
    raw = _decode(path, head)
    warnings: list[str] = []
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


def _decode(path: Path, head: bytes) -> np.ndarray:
    try:
        if head.startswith(_PNG_SIGNATURE):
            return np.asarray(imagecodecs.png_decode(path.read_bytes()))
        if head[:4] in _TIFF_SIGNATURES:
            with tifffile.TiffFile(path) as tif:
                page = tif.pages.first
                array = page.asarray()
                if array.ndim == 3 and page.planarconfig == tifffile.PLANARCONFIG.SEPARATE:
                    array = np.moveaxis(array, 0, -1)
                return array
        return np.asarray(iio.imread(path, index=0))
    except (OSError, ValueError, TypeError, RuntimeError, imagecodecs.PngError) as exc:
        if isinstance(exc, FileNotFoundError):
            raise
        raise ValueError(f"{path}: cannot decode image: {exc}") from exc


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
