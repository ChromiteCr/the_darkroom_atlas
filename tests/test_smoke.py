def test_version():
    import darkroom_atlas
    assert darkroom_atlas.ATLAS_VERSION == "D0"


def test_xp_of_numpy_and_strict():
    import array_api_strict as xps
    import numpy as np

    from darkroom_atlas.backend import xp_of
    assert xp_of(np.zeros(2)).__name__.endswith("numpy")
    assert xp_of(xps.zeros(2)) is xps or "strict" in xp_of(xps.zeros(2)).__name__


def test_sample_asset():
    from pathlib import Path

    import imageio.v3 as iio
    img = iio.imread(Path(__file__).parent / "assets" / "sample.png")
    assert img.shape == (64, 96, 3) and img.dtype.name == "uint8"
