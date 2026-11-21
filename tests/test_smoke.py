import re

import darkroom_atlas


def test_version():
    assert re.fullmatch(r"D\d+([a-z]\d*)?", darkroom_atlas.ATLAS_VERSION)
