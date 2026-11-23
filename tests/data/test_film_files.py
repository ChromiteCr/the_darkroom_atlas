"""The two Kodak film files of D1 (spec D1 1.1, 2.0, 2.1.8, 2.2.7, 2.4.11, 2.5.3, 2.6, 5.4, 7.1)."""

import copy
import math

import numpy as np
import pytest
from pydantic import ValidationError

from darkroom_atlas.data.library import Library
from darkroom_atlas.data.models import (
    BaseLayer,
    EmulsionLayer,
    Material,
    OvercoatLayer,
    bw_mechanistic_problems,
    iter_quantities,
)
from darkroom_atlas.data.provenance import Quantity, Source
from darkroom_atlas.emulsion import size_classes
from darkroom_atlas.light import LowRankLight
from darkroom_atlas.optics.kubelka import depth_response
from darkroom_atlas.optics.stack import optical_stack
from darkroom_atlas.pipeline import RunContext, expose
from darkroom_atlas.spectral import DEFAULT_GRID

FILMS = ("kodak/tri-x-400", "kodak/t-max-100")
FORMATS = ("135", "120")

# Parameters the calibration fits (spec D1 7.1 blocks B and F), as field paths below the
# emulsion layer, with the bounds and prior shape of the spec's parameter tables (2.0.2,
# 2.1.8, 2.2.7, 2.4.11, 2.5.3). "(对数)" / "log 先验" is a lognormal prior, everything else
# uniform. T-Max d_50 is parameterised by (t_g, AR) and its [0.6, 1.5] bounds are the
# one-sided penalty of 7.1.
_COMMON_FITTED = {
    "grain.sigma_log": ((0.2, 0.7), "uniform"),
    "grain.silver_coverage": ((2.0, 6.0), "lognormal"),
    "optics.intrinsic_absorption.alpha0": ((0.1, 10.0), "lognormal"),
    "optics.intrinsic_absorption.edge_energy": ((2.5, 2.9), "uniform"),
    "optics.intrinsic_absorption.urbach_energy": ((0.02, 0.15), "uniform"),
    "optics.dye_bands[0].amplitude": ((0.0, 0.1), "uniform"),
    "optics.dye_bands[0].peak_wavenumber": ((1.4, 2.0), "uniform"),
    "optics.dye_bands[0].width": ((0.03, 0.3), "uniform"),
    "optics.dye_bands[1].amplitude": ((0.0, 0.1), "uniform"),
    "optics.dye_bands[1].peak_wavenumber": ((1.4, 2.0), "uniform"),
    "optics.dye_bands[1].width": ((0.03, 0.3), "uniform"),
    "optics.scattering_550": ((0.05, 0.4), "uniform"),
    "optics.psf_weight": ((0.0, 1.0), "uniform"),
    "latent.efficiency": ((0.05, 1.0), "uniform"),
    "latent.tau0_sigma_ln": ((0.0, 4.0), "uniform"),
    "latent.ionic_time": ((1.0e-7, 1.0e-4), "lognormal"),
    "development.fog_hazard_ratio": ((1.0e-6, 1.0e-1), "lognormal"),
    "development.covering_power": ((0.3, 1.0), "uniform"),
}
FITTED = {
    "kodak/tri-x-400": {
        **_COMMON_FITTED,
        "grain.mean_diameter": ((0.4, 2.0), "lognormal"),
        "optics.psf_length": ((0.5, 5.0), "uniform"),
        "latent.tau0_median": ((0.01, 2.0), "lognormal"),
    },
    "kodak/t-max-100": {
        **_COMMON_FITTED,
        "grain.mean_diameter": ((0.6, 1.5), "lognormal"),
        "grain.thickness": ((0.08, 0.15), "uniform"),
        "optics.psf_length": ((0.3, 1.6), "uniform"),
        "latent.tau0_median": ((0.05, 10.0), "lognormal"),
    },
}

# Prior centres the spec names (2.0.2, 2.2.7, 5.8, 6.6).
NAMED_CENTRES = {
    "kodak/tri-x-400": {
        "grain.mean_diameter": 0.9,
        "grain.sigma_log": 0.45,
        "grain.silver_coverage": 4.0,
        "optics.psf_length": 2.0,
        "latent.tau0_median": 0.22,
        "development.covering_power": 0.6,
    },
    "kodak/t-max-100": {
        "grain.silver_coverage": 4.0,
        "grain.thickness": 0.12,
        "latent.tau0_median": 0.89,
        "development.covering_power": 0.6,  # 6.6 names it for Tri-X; used for both films
    },
}

# Fixed (`assumed` or `literature`) values of the spec's parameter tables: (value, bounds).
FIXED = {
    "kodak/tri-x-400": {
        "gelatin_coverage": (5.0, (3.0, 8.0)),
        "swelling_ratio": (3.0, (2.0, 5.0)),
        "grain.iodide_fraction": (0.04, (0.02, 0.10)),
        "optics.gelatin_index": (1.54, (1.53, 1.54)),
        "optics.scattering_exponent": (1.0, (0.5, 2.0)),
        "latent.efficiency_size_exponent": (0.0, (-1.0, 1.0)),
        "latent.competing_sites": (3.0, (2.0, 10.0)),
        "latent.hole_attack": (0.0, (0.0, 1.0)),
        "development.covering_power_size_exponent": (0.5, (0.0, 1.0)),
        "development.physical_covering_ratio": (0.5, (0.2, 1.0)),
        "development.halide_density": (0.3, (0.1, 0.6)),
        "development.dye_density": (0.0, (0.0, 0.1)),
        "development.fixing_factor": (1.0, None),
    },
}
FIXED["kodak/t-max-100"] = {
    **FIXED["kodak/tri-x-400"],
    "grain.iodide_fraction": (0.03, (0.01, 0.06)),
    "development.fixing_factor": (2.0, None),
}

# Base per format (spec D1 2.1.8): (thickness um, source, grey density, backing density).
BASE = {
    "kodak/tri-x-400": {
        "135": (127.0, Source.ASSUMED, 0.10, 0.0),
        "120": (91.0, Source.ASSUMED, 0.0, 0.5),
    },
    "kodak/t-max-100": {
        "135": (127.0, Source.DATASHEET, 0.10, 0.0),
        "120": (119.0, Source.DATASHEET, 0.0, 0.5),
    },
}


@pytest.fixture(scope="module")
def lib():
    return Library()


@pytest.fixture(scope="module", params=FILMS)
def film(request, lib) -> Material:
    return lib.material(request.param)


def _emulsion(material: Material) -> EmulsionLayer:
    (layer,) = material.emulsion_layers
    return layer


def _get(obj, path: str):
    for part in path.split("."):
        name, _, index = part.partition("[")
        obj = getattr(obj, name)
        if index:
            obj = obj[int(index.rstrip("]"))]
    return obj


def test_both_films_load_as_bw_mechanistic(lib):
    for film_id in FILMS:
        material = lib.material(film_id)
        assert material.id == film_id
        assert material.model == "bw-mechanistic"
        assert material.kind == "bw_negative"
        assert material.default_format == "135"
        assert set(material.formats) == set(FORMATS)


def test_bw_mechanistic_problems_report_nothing(film):
    assert bw_mechanistic_problems(film) == []


def test_check_data_has_no_problem_in_the_film_files(lib):
    problems = [p for p in lib.validate() if "films/kodak/" in p]
    assert problems == []


def test_layer_structure_and_grain_morphology(film):
    # Spec D1 2.0.1: overcoat, emulsion, base; Tri-X conventional, T-Max tabular (1.1).
    assert [(type(layer), layer.id) for layer in film.layers] == [
        (OvercoatLayer, "overcoat"),
        (EmulsionLayer, "emulsion"),
        (BaseLayer, "base"),
    ]
    grain = _emulsion(film).grain
    expected = "conventional" if film.id == "kodak/tri-x-400" else "tabular"
    assert grain.morphology == expected
    assert grain.halide == "AgBrI"
    assert grain.n_classes == 6
    assert (grain.thickness is not None) == (expected == "tabular")
    assert len(_emulsion(film).optics.dye_bands) == 2  # J = 2 (spec D1 2.1.8)
    assert _emulsion(film).latent.cluster_cap == 5


def test_every_fitted_parameter_has_the_spec_bounds_and_prior(film):
    emulsion = _emulsion(film)
    for path, (bounds, prior) in FITTED[film.id].items():
        quantity = _get(emulsion, path)
        assert isinstance(quantity, Quantity), path
        assert quantity.bounds == pytest.approx(bounds, rel=1e-12), path
        assert quantity.prior == prior, path
        # Not fitted yet: assumed at the prior centre until the D1e write-back (7.8).
        assert quantity.source is Source.ASSUMED, path
        assert quantity.fit is None, path
        assert "PLACEHOLDER" in (quantity.note or ""), path


def test_named_prior_centres(film):
    emulsion = _emulsion(film)
    for path, value in NAMED_CENTRES[film.id].items():
        assert _get(emulsion, path).value == pytest.approx(value, rel=1e-12), path


def test_unnamed_prior_centres_follow_the_midpoint_rule(film):
    # Geometric midpoint of the bounds for a lognormal prior, arithmetic otherwise, unless the
    # spec names the centre or the note gives another reason (T-Max dye amplitude: A_550).
    emulsion = _emulsion(film)
    for path, ((lo, hi), prior) in FITTED[film.id].items():
        if path in NAMED_CENTRES[film.id] or "peak_wavenumber" in path:
            continue
        quantity = _get(emulsion, path)
        if film.id == "kodak/t-max-100" and path in {
            "grain.mean_diameter",
            "optics.dye_bands[0].amplitude",
            "optics.dye_bands[1].amplitude",
        }:
            continue
        centre = math.sqrt(lo * hi) if prior == "lognormal" else 0.5 * (lo + hi)
        assert quantity.value == pytest.approx(centre, rel=1e-12), path


def test_tmax_grain_respects_the_aspect_ratio_box():
    # Spec D1 2.0.2, 7.1: AR = d_50/t_g in [8, 19], d_50 in [0.6, 1.5] um; centre AR = sqrt(152).
    grain = _emulsion(Library().material("kodak/t-max-100")).grain
    d50, t_g = grain.mean_diameter.to("um"), grain.thickness.to("um")
    assert 8.0 <= d50 / t_g <= 19.0
    assert 0.6 <= d50 <= 1.5
    assert d50 / t_g == pytest.approx(math.sqrt(8.0 * 19.0), rel=1e-12)


def test_fixed_parameters(film):
    emulsion = _emulsion(film)
    for path, (value, bounds) in FIXED[film.id].items():
        quantity = _get(emulsion, path)
        assert quantity.source is Source.ASSUMED, path
        assert quantity.value == pytest.approx(value, abs=1e-15), path
        if bounds is None:
            assert quantity.bounds is None, path
        else:
            assert quantity.bounds == pytest.approx(bounds), path
    assert film.layers[0].thickness.value == 1.0
    assert film.layers[0].thickness.bounds == (0.5, 2.0)


def test_base_per_format(film):
    base = film.layers[2]
    assert base.refractive_index.value == 1.48
    assert base.refractive_index.source is Source.LITERATURE
    for name, (thickness, source, grey, backing) in BASE[film.id].items():
        fmt = base.formats[name]
        assert fmt.thickness.to("um") == thickness
        assert fmt.thickness.source is source
        assert fmt.grey_density.value == grey
        assert fmt.backing_density.value == backing
        assert fmt.grey_density.source is Source.ASSUMED
        assert fmt.backing_density.source is Source.ASSUMED


def test_datasheet_only_for_what_kodak_states(film):
    # ISO speed for both films, T-Max base thicknesses; everything else is not datasheet.
    allowed = {"iso_speed"}
    if film.id == "kodak/t-max-100":
        allowed |= {"layers[2].formats.135.thickness", "layers[2].formats.120.thickness"}
    datasheet = set()
    for loc, quantity in iter_quantities(film):
        if quantity.source is Source.DATASHEET:
            assert quantity.ref, loc
            datasheet.add(loc)
        if quantity.source is Source.ASSUMED:
            assert quantity.note, f"{loc}: an assumed value needs a note saying where it is from"
    assert datasheet == allowed
    expected_iso = 400 if film.id == "kodak/tri-x-400" else 100
    assert film.iso_speed.value == expected_iso


@pytest.mark.parametrize("film_format", FORMATS)
def test_layer_absorptance_at_550_within_the_constraint(film, film_format):
    # Spec D1 2.1.8, 7.1: A_550 in [0.2, 0.6] at the prior centres (two-stream transport and
    # the single-pass Beer value both).
    classes = size_classes(_emulsion(film), 6)
    stack = optical_stack(film, film_format, classes, DEFAULT_GRID)
    response = depth_response(stack, 8)
    band = int(np.argmin(np.abs(DEFAULT_GRID.wavelengths - 550.0)))
    two_stream = float(response.absorbed[band])
    single_pass = 1.0 - math.exp(-stack.mu_a_per_um[band] * stack.emulsion_thickness_um)
    assert 0.2 <= two_stream <= 0.6
    assert 0.2 <= single_pass <= 0.6


@pytest.mark.parametrize("film_format", FORMATS)
def test_expose_uniform_frame_gives_finite_positive_fields(lib, film, film_format):
    shape = (6, 7)
    light = LowRankLight(
        np.ones((*shape, 3)), np.full((3, DEFAULT_GRID.n_bands), 20.0), DEFAULT_GRID
    )
    ctx = RunContext(lib, long_edge_mm=0.3, film_format=film_format)
    state = expose(film, light, 0.01, ctx)
    assert state.meta["exposure"]["film_format"] == film_format
    for name in ("emulsion/exposure_volume", "emulsion/exposure_surface"):
        field = np.asarray(state.fields[name])
        assert field.shape == (*shape, 2)
        assert np.all(np.isfinite(field))
        assert np.all(field > 0.0)


def test_missing_field_is_reported():
    material = Library().material("kodak/tri-x-400")
    data = material.model_dump(mode="python", exclude_none=True)
    del data["layers"][1]["latent"]
    broken = Material.model_validate(data)
    assert bw_mechanistic_problems(broken) == [
        "layers[1].latent: required by the bw-mechanistic model"
    ]


def test_value_outside_the_bounds_is_rejected():
    material = Library().material("kodak/t-max-100")
    data = copy.deepcopy(material.model_dump(mode="python", exclude_none=True))
    data["layers"][1]["grain"]["thickness"]["value"] = 0.2  # bounds [0.08, 0.15]
    with pytest.raises(ValidationError, match="outside bounds"):
        Material.model_validate(data)


def test_spectral_sensitivity_curve_files_resolve_and_load(lib, film):
    # Spec D1 5.4: the emulsion's spectral_sensitivity points at the film's training curve.
    emulsion = next(layer for layer in film.layers if isinstance(layer, EmulsionLayer))
    curve = emulsion.spectral_sensitivity
    assert curve is not None
    assert curve.source is Source.DATASHEET and curve.ref
    assert curve.file == f"reference/{film.id}/spectral-sensitivity-d0.3.csv"
    assert lib.curve_files(film) == [curve.file]
    reference = lib.reference_curve(curve.file.removeprefix("reference/"))
    assert reference.material == film.id
    assert reference.ref == curve.ref
    values = lib.resolve_curve(curve, DEFAULT_GRID)
    assert values.shape == (DEFAULT_GRID.n_bands,)
    assert np.all(np.isfinite(values))
    assert values.max() > 0.0
