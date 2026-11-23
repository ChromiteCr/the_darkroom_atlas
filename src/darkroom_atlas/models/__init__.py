"""Built-in material models. Importing this package registers them with the registry.

``ideal`` (D0) and ``bw-mechanistic`` (D1): its exposure model (:mod:`.bw`), its process steps
(:mod:`darkroom_atlas.develop.step`, :mod:`darkroom_atlas.develop.post`) and its readout
devices (:mod:`darkroom_atlas.readout.densitometer`, :mod:`darkroom_atlas.readout.scan`).
"""

from darkroom_atlas.develop import post, step
from darkroom_atlas.models import bw, ideal
from darkroom_atlas.readout import densitometer, scan

__all__ = ["bw", "densitometer", "ideal", "post", "scan", "step"]
