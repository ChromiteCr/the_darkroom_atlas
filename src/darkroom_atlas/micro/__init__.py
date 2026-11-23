"""Microscope mode (spec D1 2.7.2): explicit grains, coverage rendering and apertures."""

from darkroom_atlas.micro.grains import GrainSet, place_grains
from darkroom_atlas.micro.measure import aperture_transmittances
from darkroom_atlas.micro.render import render_coverage

__all__ = ["GrainSet", "aperture_transmittances", "place_grains", "render_coverage"]
