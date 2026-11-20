"""Generate ``sample.png``, the deterministic 96 x 64 8-bit RGB test image.

Layout (width 96, height 64):

* Left 48 columns: horizontal grey ramp, column 0 = 0 and column 47 = 255.
* Right 48 columns: 2 rows x 4 columns of 12 x 32 px patches, row-major:
  red, green, blue, white / cyan, magenta, yellow, black.

Run with ``uv run python tests/assets/make_sample.py``; the file is written next to this
script regardless of the working directory.
"""

from pathlib import Path

import numpy as np
from PIL import Image

WIDTH = 96
HEIGHT = 64
RAMP_WIDTH = 48
PATCH_ROWS = 2
PATCH_COLS = 4

PATCH_COLOURS = (
    (255, 0, 0),  # red
    (0, 255, 0),  # green
    (0, 0, 255),  # blue
    (255, 255, 255),  # white
    (0, 255, 255),  # cyan
    (255, 0, 255),  # magenta
    (255, 255, 0),  # yellow
    (0, 0, 0),  # black
)

OUTPUT = Path(__file__).parent / "sample.png"


def make_sample() -> np.ndarray:
    """Return the sample image as a ``(64, 96, 3)`` uint8 array."""
    img = np.zeros((HEIGHT, WIDTH, 3), dtype=np.uint8)

    ramp = np.round(np.linspace(0.0, 255.0, RAMP_WIDTH)).astype(np.uint8)
    img[:, :RAMP_WIDTH, :] = ramp[np.newaxis, :, np.newaxis]

    patch_w = (WIDTH - RAMP_WIDTH) // PATCH_COLS
    patch_h = HEIGHT // PATCH_ROWS
    for index, colour in enumerate(PATCH_COLOURS):
        row, col = divmod(index, PATCH_COLS)
        x0 = RAMP_WIDTH + col * patch_w
        y0 = row * patch_h
        img[y0 : y0 + patch_h, x0 : x0 + patch_w, :] = colour

    return img


def main() -> None:
    Image.fromarray(make_sample()).save(OUTPUT, format="PNG")
    print(OUTPUT)


if __name__ == "__main__":
    main()
