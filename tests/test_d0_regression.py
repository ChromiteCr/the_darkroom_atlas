"""Regression guard: the D0 ideal-material sample must stay bit-identical through D1.

The hashes were computed from the output of the D0 code at commit ``b180f08`` (before the
pipeline gained segment fusion and the frame assembler), on ``tests/assets/sample.png``
with material ``ideal`` and process ``identity`` and every other ``run_image`` setting at
its default (spec §10, general completion criteria; ``test_d0_ideal_regression``).
"""

import hashlib
from pathlib import Path

import imagecodecs
import numpy as np
import pytest

from darkroom_atlas.data.library import Library
from darkroom_atlas.pipeline import RunContext, expose, process, readout
from darkroom_atlas.runner import DEFAULT_ISO, run_image
from darkroom_atlas.scene.camera import scene_light
from darkroom_atlas.scene.io import read_image

SAMPLE = Path(__file__).resolve().parent / "assets" / "sample.png"

# sha256 of the decoded uint16 pixels of the written PNG, shape (64, 96, 3).
D0_PIXELS_SHA256 = "3e766f19c13297fee0923618f2bb3d69cc68accccc773e7ed679a54438c6592b"
# sha256 of the float64 ``ideal_colorimetric`` readout data, shape (64, 96, 3).
D0_READOUT_SHA256 = "5f01b10a43f412ebaa8ba8319b34d20c52a302050f7547e86f83eb24ad47e1fb"


def _sha256(array: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


@pytest.mark.parametrize("tile_px", [None, 7, 16])
def test_d0_ideal_regression(tmp_path, tile_px):
    out, _ = run_image(
        SAMPLE, tmp_path, material_id="ideal", process_id="identity", tile_px=tile_px
    )
    pixels = imagecodecs.png_decode(out.read_bytes())
    assert pixels.dtype == np.uint16
    assert pixels.shape == (64, 96, 3)
    assert _sha256(pixels) == D0_PIXELS_SHA256


@pytest.mark.parametrize("tile_px", [None, 7])
def test_d0_ideal_regression_readout(tile_px):
    lib = Library()
    material = lib.material("ideal")
    ctx = RunContext(lib, tile_px=tile_px)
    time_s = 1 / 125
    light = scene_light(
        read_image(SAMPLE).linear_rgb,
        grid=ctx.grid,
        illuminant="D65",
        iso=DEFAULT_ISO,
        ev=0.0,
        time_s=time_s,
    )
    state = process(expose(material, light, time_s, ctx), lib.process("identity"), ctx)
    result = readout(state, material.default_readout, ctx)
    assert result.data.dtype == np.float64
    assert _sha256(result.data) == D0_READOUT_SHA256
