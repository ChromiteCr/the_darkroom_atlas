"""The ``ideal`` model: a linear recorder of CIE XYZ that checks the spectral pipeline.

It is not a physical material. Each emulsion layer records the energy it receives,
weighted by its spectral sensitivity; nothing is lost, saturated, scattered or random.

* Exposure (``"ideal"``): ``"<layer>/exposure" = sum_b(N_b * E_b * S_b)`` in J per um^2,
  where ``N_b`` is the photon count per um^2 in band ``b`` (flux times exposure time),
  ``E_b`` the photon energy at the band centre and ``S_b`` the layer's spectral sensitivity
  (resolved with :meth:`Library.resolve_curve`).
* Step ``("ideal", identity)``: ``"<layer>/signal" = "<layer>/exposure"``, pointwise.
* Readout ``"ideal_colorimetric"``: needs layers ``x y z`` with a ``signal``;
  ``XYZ = signal / k`` with ``k = state.meta["light"]["energy_scale_j_per_um2"]``, then
  linear sRGB ``= XYZ_TO_LINEAR_SRGB @ XYZ``.

Why this inverts the camera: scene light carries ``N_b * E_b = k * e_b`` (relative energy
``e = reflectance * illuminant``), so with the CIE 1931 colour-matching functions as
sensitivities the signal is ``k * XYZ`` of the scene under its illuminant, whatever the
ISO, EV or exposure time. Dividing by ``k`` gives that XYZ (an ideal white has ``Y == 1``).
"""

import math
from collections.abc import Mapping
from numbers import Real
from typing import Any

import numpy as np

from darkroom_atlas.data.models import Material, ProcessStep, StepKind
from darkroom_atlas.light import Light
from darkroom_atlas.pipeline import Readout, RunContext, StepContext, unexposed_state
from darkroom_atlas.registry import register_exposure, register_readout, register_step
from darkroom_atlas.spectral import XYZ_TO_LINEAR_SRGB
from darkroom_atlas.state import EmulsionState
from darkroom_atlas.units import photon_energy_j

__all__ = ["MODEL", "READOUT", "expose_ideal", "identity_ideal", "read_ideal_colorimetric"]

MODEL = "ideal"
"""Name of the exposure model and of the step handlers' model."""

READOUT = "ideal_colorimetric"
"""Name of the colorimetric readout device."""

_XYZ_LAYERS = ("x", "y", "z")


@register_exposure(MODEL)
def expose_ideal(material: Material, light: Light, time_s: float, ctx: RunContext) -> EmulsionState:
    """Record, for each emulsion layer, the sensitivity-weighted energy received (J/um^2)."""
    grid = light.grid
    photon_energy = photon_energy_j(grid.wavelengths)
    fields = {}
    for layer in material.emulsion_layers:
        if layer.spectral_sensitivity is None:
            raise ValueError(
                f"material {material.id!r}: emulsion layer {layer.id!r} has no "
                "spectral_sensitivity, which the ideal model needs"
            )
        sensitivity = ctx.library.resolve_curve(layer.spectral_sensitivity, grid)
        fields[f"{layer.id}/exposure"] = (
            light.photon_flux @ (photon_energy * sensitivity)
        ) * time_s
    return unexposed_state(material, light, ctx).with_fields(fields)


@register_step(MODEL, StepKind.IDENTITY, interaction_range_um=0.0)
def identity_ideal(state: EmulsionState, step: ProcessStep, ctx: StepContext) -> EmulsionState:
    """Copy each emulsion layer's exposure to its signal (no chemistry)."""
    return state.with_fields(
        {
            f"{layer.id}/signal": state.get(layer.id, "exposure")
            for layer in state.material.emulsion_layers
        }
    )


@register_readout(READOUT)
def read_ideal_colorimetric(state: EmulsionState, ctx: RunContext, **params: Any) -> Readout:
    """Linear sRGB of the recorded XYZ, normalised by the scene's energy scale."""
    if params:
        raise TypeError(f"{READOUT} takes no parameters, got {', '.join(sorted(params))}")
    missing = [f"{layer}/signal" for layer in _XYZ_LAYERS if not state.has(layer, "signal")]
    if missing:
        raise ValueError(
            f"{READOUT} needs the fields x/signal, y/signal and z/signal (an ideal material "
            f"after processing); missing: {', '.join(missing)}"
        )
    light = state.meta.get("light")
    scale = light.get("energy_scale_j_per_um2") if isinstance(light, Mapping) else None
    if (
        isinstance(scale, bool)
        or not isinstance(scale, Real)
        or not math.isfinite(scale)
        or scale <= 0.0
    ):
        raise ValueError(
            f"{READOUT} needs a positive, finite state.meta['light']['energy_scale_j_per_um2'] "
            f"(recorded by scene_light), got {scale!r}"
        )
    xyz = np.stack([state.get(layer, "signal") for layer in _XYZ_LAYERS], axis=-1) / scale
    rgb = xyz @ XYZ_TO_LINEAR_SRGB.T
    return Readout(kind="image_linear_srgb", data=rgb, meta={"device": READOUT})
