"""Light transport in the film during exposure (spec §2.1): depth response, scatter, halation."""

from darkroom_atlas.optics.halation import (
    HalationKernel,
    halation_fraction,
    halation_kernel,
    halation_otf,
)
from darkroom_atlas.optics.kubelka import DepthResponse, depth_response
from darkroom_atlas.optics.psf import scatter_kernel, scatter_otf
from darkroom_atlas.optics.stack import BaseOptics, OpticalStack

__all__ = [
    "BaseOptics",
    "DepthResponse",
    "HalationKernel",
    "OpticalStack",
    "depth_response",
    "halation_fraction",
    "halation_kernel",
    "halation_otf",
    "scatter_kernel",
    "scatter_otf",
]
