"""Virtual sensitometer results (spec D1 section 6.1).

:class:`StripResult` holds the densities read off one step wedge at a series of development
times: ``log_h`` (J,) in log10(lux*s), ``times_s`` (T,), visual and Type 2 densities (T, J),
and the densities of the unexposed patch (T,). All arrays are float64 and read-only.
"""

from dataclasses import dataclass

import numpy as np

from darkroom_atlas._frozen import readonly_array
from darkroom_atlas.validate.conditions import Condition

__all__ = ["StripResult"]


@dataclass(frozen=True, eq=False)
class StripResult:
    """Densities of a sensitometric strip over a development series (spec D1 6.1)."""

    condition: Condition
    log_h: np.ndarray
    times_s: np.ndarray
    visual: np.ndarray
    type2: np.ndarray
    base_fog_visual: np.ndarray
    base_fog_type2: np.ndarray

    def __post_init__(self) -> None:
        for name in ("log_h", "times_s", "visual", "type2", "base_fog_visual", "base_fog_type2"):
            object.__setattr__(self, name, readonly_array(getattr(self, name), np.float64))
        if self.log_h.ndim != 1 or self.log_h.size < 2:
            raise ValueError(
                f"StripResult.log_h must be 1-D with >= 2 steps, got {self.log_h.shape}"
            )
        if not (np.all(np.isfinite(self.log_h)) and np.all(np.diff(self.log_h) > 0.0)):
            raise ValueError("StripResult.log_h must be finite and strictly increasing")
        if self.times_s.ndim != 1 or self.times_s.size < 1:
            raise ValueError(
                f"StripResult.times_s must be 1-D, non-empty, got {self.times_s.shape}"
            )
        n_t, n_j = self.times_s.size, self.log_h.size
        for name in ("visual", "type2"):
            if getattr(self, name).shape != (n_t, n_j):
                raise ValueError(
                    f"StripResult.{name} must have shape {(n_t, n_j)}, "
                    f"got {getattr(self, name).shape}"
                )
        for name in ("base_fog_visual", "base_fog_type2"):
            if getattr(self, name).shape != (n_t,):
                raise ValueError(
                    f"StripResult.{name} must have shape {(n_t,)}, got {getattr(self, name).shape}"
                )
