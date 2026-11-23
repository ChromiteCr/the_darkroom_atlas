"""Calibration parameter table, transforms and library overlay (spec D1 7.1, 7.3, 7.5, 5.1)."""

import hashlib
import math
import shutil
from pathlib import Path

import numpy as np
import pytest

from darkroom_atlas.calibrate.parameters import (
    LibraryOverlay,
    ParameterError,
    ParameterSpec,
    current_values,
    from_theta,
    parameter_specs,
    parameter_value,
    theta_bounds,
    to_theta,
)
from darkroom_atlas.data import DataError, Library
from darkroom_atlas.data.models import iter_quantities

TX = "films/kodak/tri-x-400:emulsion"
TM = "films/kodak/t-max-100:emulsion"
KIN = "chemistry/kinetics"

STAGE_A = [
    f"{KIN}:activity.delta_eps.sulfite",
    f"{KIN}:activity.delta_eps.borate",
    f"{KIN}:activity.delta_eps.organic",
    f"{KIN}:activity.delta_eps.amine",
    "developers/d-76:components.boric-oxide.amount",
    "developers/t-max-dev:components.sulfur-dioxide.amount",
    "developers/hc-110-pre2019:components.potassium-hydroxide.amount",
    "calibration:delta_cc",
]


def _spectral(prefix: str) -> list[str]:
    return [
        f"{prefix}.optics.intrinsic_absorption.alpha0",
        f"{prefix}.optics.intrinsic_absorption.edge_energy",
        f"{prefix}.optics.intrinsic_absorption.urbach_energy",
        f"{prefix}.optics.dye_bands[0].amplitude",
        f"{prefix}.optics.dye_bands[0].peak_wavenumber",
        f"{prefix}.optics.dye_bands[0].width",
        f"{prefix}.optics.dye_bands[1].amplitude",
        f"{prefix}.optics.dye_bands[1].peak_wavenumber",
        f"{prefix}.optics.dye_bands[1].width",
    ]


STAGE_B = _spectral(TX) + _spectral(TM)

STAGE_C_TRI_X = [
    f"{TX}.grain.mean_diameter",
    f"{TX}.grain.sigma_log",
    f"{TX}.latent.efficiency",
    f"{TX}.development.covering_power",
    f"{TX}.grain.silver_coverage",
    f"{TX}.development.fog_hazard_ratio",
    f"{TX}.optics.scattering_550",
    "calibration:kodak/tri-x-400/s_abs",
    f"{TX}.optics.psf_weight",
    f"{TX}.optics.psf_length",
    f"{TX}.latent.tau0_median",
    f"{TX}.latent.tau0_sigma_ln",
    f"{TX}.latent.ionic_time",
]
STAGE_C_T_MAX = [
    "calibration:kodak/t-max-100/aspect_ratio",
    f"{TM}.grain.sigma_log",
    f"{TM}.latent.efficiency",
    f"{TM}.development.covering_power",
    f"{TM}.grain.silver_coverage",
    f"{TM}.development.fog_hazard_ratio",
    f"{TM}.optics.scattering_550",
    "calibration:kodak/t-max-100/s_abs",
    f"{TM}.optics.psf_weight",
    f"{TM}.optics.psf_length",
    f"{TM}.latent.tau0_median",
    f"{TM}.latent.tau0_sigma_ln",
    f"{TM}.latent.ionic_time",
    f"{TM}.grain.thickness",
]
STAGE_C = STAGE_C_TRI_X + STAGE_C_T_MAX

STAGE_D = [
    f"{KIN}:agents.metol.k_et",
    f"{KIN}:agents.metol.K_ads",
    f"{KIN}:agents.metol.lambda0",
    f"{KIN}:agents.metol.k_des",
    f"{KIN}:agents.metol.activation_energy",
    f"{KIN}:agents.hydroquinone.k_et",
    f"{KIN}:agents.hydroquinone.K_ads",
    f"{KIN}:agents.hydroquinone.lambda0",
    f"{KIN}:agents.hydroquinone.k_reg",
    f"{KIN}:agents.hydroquinone.activation_energy",
    f"{KIN}:agents.pyrazolidone.k_et",
    f"{KIN}:agents.pyrazolidone.K_ads",
    f"{KIN}:agents.pyrazolidone.lambda0",
    f"{KIN}:agents.pyrazolidone.k_des",
    f"{KIN}:agents.pyrazolidone.pka",
    f"{KIN}:agents.ascorbate.k_et",
    f"{KIN}:agents.ascorbate.K_ads",
    f"{KIN}:agents.ascorbate.lambda0",
    f"{KIN}:agents.ascorbate.k_reg",
    f"{KIN}:halides.K_Br",
    f"{KIN}:barriers.u_gel",
    f"{KIN}:barriers.u_halide",
    f"{KIN}:latent.u_li",
    f"{KIN}:solvent.latent_attrition",
    f"{KIN}:solvent.k_diss",
    f"{KIN}:solvent.c_sat_ref",
    f"{KIN}:solvent.sulfite_exponent",
    f"{KIN}:solvent.k_pd",
    f"{KIN}:transport.delta_max",
    f"{KIN}:transport.delta_agit",
    "developers/t-max-dev:components.hydroquinone.amount",
    "developers/t-max-dev:components.dimezone-s.amount",
    "developers/t-max-dev:components.sodium-bromide.amount",
    "developers/hc-110-pre2019:components.hydroquinone.amount",
    "developers/hc-110-pre2019:components.methylphenidone.amount",
    "developers/hc-110-pre2019:components.potassium-bromide.amount",
]


@pytest.fixture(scope="module")
def lib() -> Library:
    return Library()


def _by_name(specs: list[ParameterSpec]) -> dict[str, ParameterSpec]:
    return {s.name: s for s in specs}


def _data_hashes(root: Path) -> dict[str, str]:
    return {
        p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(root.rglob("*"))
        if p.is_file()
    }


def _copy_library(tmp_path: Path) -> Path:
    src = Library().root
    root = tmp_path / "data"
    for part in ("films", "chemistry", "developers"):
        shutil.copytree(src / part, root / part)
    return root


def _replace_in(path: Path, old: str, new: str) -> None:
    text = path.read_text(encoding="utf-8")
    assert text.count(old) == 1, old
    path.write_text(text.replace(old, new), encoding="utf-8")


# ---- stage lists (spec D1 7.1, 7.3) --------------------------------------------------------


@pytest.mark.parametrize(
    ("stage", "expected"),
    [
        ("A", STAGE_A),
        ("B", STAGE_B),
        ("C", STAGE_C),
        ("D", STAGE_D),
        ("E", STAGE_C + STAGE_D),
    ],
)
def test_stage_parameter_list(lib: Library, stage: str, expected: list[str]) -> None:
    specs = parameter_specs(lib, stage)
    assert [s.name for s in specs] == expected
    expected_stage = {"A": "A", "B": "B", "C": "C", "D": "D"}
    for spec in specs:
        if stage != "E":
            assert spec.stage == expected_stage[stage]


def test_stage_counts_match_spec(lib: Library) -> None:
    # Spec D1 7.1: A 4 + 4; B 9 per film (two dye bands); C Tri-X 13, T-Max 14; D 19 + 11 + 6.
    assert len(STAGE_A) == 8
    assert len(STAGE_B) == 18
    assert (len(STAGE_C_TRI_X), len(STAGE_C_T_MAX)) == (13, 14)
    assert len(STAGE_D) == 36
    films = parameter_specs(lib, "C", films=["kodak/tri-x-400"])
    assert [s.name for s in films] == STAGE_C_TRI_X


def test_stage_record_names_are_accepted(lib: Library) -> None:
    for short, long in [("A", "A_activity"), ("B", "B_spectral"), ("C", "C_film")]:
        assert parameter_specs(lib, long) == parameter_specs(lib, short)
    assert parameter_specs(lib, "D_chemistry") == parameter_specs(lib, "D")
    assert parameter_specs(lib, "E_polish") == parameter_specs(lib, "E")


@pytest.mark.parametrize("stage", ["F", "G", "a", "", None, "C_chemistry"])
def test_unknown_stage_raises(lib: Library, stage: object) -> None:
    with pytest.raises(ValueError, match="unknown stage"):
        parameter_specs(lib, stage)  # type: ignore[arg-type]


def test_blocks(lib: Library) -> None:
    a = _by_name(parameter_specs(lib, "A"))
    assert {s.block for s in a.values()} == {"A"}
    assert {s.block for s in parameter_specs(lib, "B")} == {"B"}
    assert {s.block for s in parameter_specs(lib, "C")} == {"F"}
    assert {s.block for s in parameter_specs(lib, "D")} == {"C"}


# ---- spec fields read from the data files ------------------------------------------------


def test_spec_fields_from_data(lib: Library) -> None:
    specs = _by_name(parameter_specs(lib, "A") + parameter_specs(lib, "C"))
    specs.update(_by_name(parameter_specs(lib, "D")))

    tau0 = specs[f"{TX}.latent.tau0_median"]
    assert (tau0.unit, tau0.transform, tau0.prior) == ("s", "log10", "lognormal")
    assert (tau0.lower, tau0.upper) == (0.01, 2.0)
    assert tau0.prior_mu == pytest.approx(math.log10(0.22), rel=1e-15)
    assert tau0.prior_tau == math.inf
    assert (tau0.block, tau0.stage) == ("F", "C")

    eps = specs[f"{KIN}:activity.delta_eps.sulfite"]
    assert (eps.transform, eps.prior, eps.prior_mu, eps.prior_tau) == (
        "linear",
        "normal",
        0.0,
        0.15,
    )
    assert (eps.lower, eps.upper) == (-0.45, 0.45)

    # Spec D1 2.4.11 gives the metol activation energy a normal prior and no box.
    ea = specs[f"{KIN}:agents.metol.activation_energy"]
    assert (ea.prior, ea.prior_mu, ea.prior_tau) == ("normal", 61.8, 10.0)
    assert (ea.lower, ea.upper) == (-math.inf, math.inf)

    efficiency = specs[f"{TX}.latent.efficiency"]
    assert efficiency.transform == "logit"
    assert specs[f"{TX}.optics.psf_weight"].transform == "logit"
    assert specs[f"{KIN}:solvent.sulfite_exponent"].transform == "linear"
    assert specs[f"{TM}.grain.thickness"].transform == "linear"

    # Bath nuisances: uniform inside the SDS bounds (spec D1 5.5).
    boric = specs["developers/d-76:components.boric-oxide.amount"]
    assert (boric.unit, boric.prior, boric.lower, boric.upper) == ("g/L", "uniform", 0.0, 3.8)

    cc = specs["calibration:delta_cc"]
    assert (cc.transform, cc.lower, cc.upper, cc.prior, cc.prior_mu) == (
        "linear",
        -0.5,
        0.5,
        "uniform",
        0.0,
    )
    s_abs = specs["calibration:kodak/tri-x-400/s_abs"]
    assert (s_abs.transform, s_abs.lower, s_abs.upper, s_abs.prior_mu) == ("linear", -1, 1, 0)
    ar = specs["calibration:kodak/t-max-100/aspect_ratio"]
    assert (ar.transform, ar.lower, ar.upper, ar.block, ar.stage) == ("log10", 8, 19, "F", "C")
    grain = lib.material("kodak/t-max-100").layers[1].grain
    ratio = grain.mean_diameter.value / grain.thickness.value
    assert ar.prior_mu == pytest.approx(math.log10(ratio), rel=1e-15)
    assert ratio == pytest.approx(math.sqrt(8 * 19), rel=1e-12)


def test_every_lognormal_prior_is_log10(lib: Library) -> None:
    for spec in parameter_specs(lib, "A") + parameter_specs(lib, "B") + parameter_specs(lib, "E"):
        if spec.prior == "lognormal":
            assert spec.transform == "log10", spec.name


def test_missing_prior_raises_naming_file_and_field(tmp_path: Path) -> None:
    root = _copy_library(tmp_path)
    film = root / "films" / "kodak" / "tri-x-400.yaml"
    _replace_in(film, "bounds: [0.01, 2.0], prior: lognormal,", "bounds: [0.01, 2.0],")
    library = Library(root)
    with pytest.raises(DataError, match=r"tri-x-400\.yaml: emulsion\.latent\.tau0_median.*prior"):
        parameter_specs(library, "C")
    assert len(parameter_specs(library, "B")) == 18  # other stages are unaffected


def test_missing_bounds_raises_naming_file_and_field(tmp_path: Path) -> None:
    root = _copy_library(tmp_path)
    kinetics = root / "chemistry" / "kinetics.yaml"
    _replace_in(
        kinetics,
        'bounds: [1.0e-6, 1.0e-2], prior: lognormal, note: "PLACEHOLDER (log centre of the bounds); KIN 2.2"}\n  c_sat',
        'prior: lognormal, note: "x"}\n  c_sat',
    )
    with pytest.raises(DataError, match=r"kinetics\.yaml: solvent\.k_diss.*bounds"):
        parameter_specs(Library(root), "D")


def test_missing_nuisance_listing_raises(tmp_path: Path) -> None:
    root = _copy_library(tmp_path)
    _replace_in(root / "developers" / "d-76.yaml", "nuisance: [boric-oxide]\n", "")
    with pytest.raises(DataError, match=r"d-76\.yaml: nuisance: 'boric-oxide'"):
        parameter_specs(Library(root), "A")


# ---- transforms ------------------------------------------------------------------------


def _sample_values(specs: list[ParameterSpec], rng: np.random.Generator) -> dict[str, float]:
    values = {}
    for spec in specs:
        lo, hi = spec.lower, spec.upper
        if not math.isfinite(lo) or not math.isfinite(hi):
            lo, hi = spec.prior_mu - 3 * spec.prior_tau, spec.prior_mu + 3 * spec.prior_tau
        if spec.transform != "linear":
            lo = max(lo, 1e-6)
        if spec.transform == "logit":
            hi = min(hi, 1 - 1e-6)
        if spec.transform == "log10":
            values[spec.name] = 10 ** rng.uniform(math.log10(lo), math.log10(hi))
        else:
            values[spec.name] = rng.uniform(lo, hi)
    return values


def test_parameter_transform_round_trip(lib: Library) -> None:
    specs = parameter_specs(lib, "A") + parameter_specs(lib, "B") + parameter_specs(lib, "E")
    rng = np.random.default_rng(0)
    cases = [current_values(lib, specs)] + [_sample_values(specs, rng) for _ in range(50)]
    for values in cases:
        back = from_theta(specs, to_theta(specs, values))
        assert list(back) == [s.name for s in specs]
        for name, v in values.items():
            assert abs(back[name] - v) <= 1e-14 * abs(v), (name, v, back[name])


def test_transform_values() -> None:
    def spec(transform: str, lo: float, hi: float) -> ParameterSpec:
        return ParameterSpec("x", "1", transform, lo, hi, "uniform", 0.0, math.inf, "F", "C")

    assert to_theta([spec("log10", 1e-3, 1e3)], {"x": 100.0})[0] == pytest.approx(2.0, abs=1e-15)
    assert to_theta([spec("logit", 0.0, 1.0)], {"x": 0.5})[0] == 0.0
    assert to_theta([spec("logit", 0.0, 1.0)], {"x": 0.75})[0] == pytest.approx(math.log(3))
    assert to_theta([spec("linear", -1.0, 1.0)], {"x": -0.25})[0] == -0.25
    small = from_theta([spec("logit", 0.0, 1.0)], np.array([-700.0]))["x"]
    assert small == pytest.approx(math.exp(-700.0), rel=1e-12)


def test_theta_bounds_round_trip(lib: Library) -> None:
    specs = parameter_specs(lib, "C")
    lo, hi = theta_bounds(specs)
    for spec, a, b in zip(specs, lo, hi, strict=True):
        assert a <= b
        if spec.transform == "logit" and spec.upper == 1.0:
            assert b == math.inf
        for t, bound in ((a, spec.lower), (b, spec.upper)):
            if math.isfinite(t):
                value = from_theta([spec], np.array([t]))[spec.name]
                assert spec.lower <= value <= spec.upper
                assert value == pytest.approx(bound, rel=1e-14, abs=1e-300)


def test_to_theta_rejects_bad_values(lib: Library) -> None:
    specs = parameter_specs(lib, "A")
    values = current_values(lib, specs)
    missing = dict(values)
    missing.pop("calibration:delta_cc")
    with pytest.raises(ParameterError, match="missing"):
        to_theta(specs, missing)
    with pytest.raises(ParameterError, match="unknown"):
        to_theta(specs, {**values, "calibration:nope": 0.0})
    with pytest.raises(ParameterError, match="outside"):
        to_theta(specs, {**values, "calibration:delta_cc": 0.6})
    with pytest.raises(ParameterError, match="finite"):
        to_theta(specs, {**values, "calibration:delta_cc": math.nan})
    psf = [s for s in parameter_specs(lib, "C") if s.name == f"{TX}.optics.psf_weight"]
    with pytest.raises(ParameterError, match="not finite"):
        to_theta(psf, {psf[0].name: 0.0})


def test_from_theta_rejects_bad_theta(lib: Library) -> None:
    specs = parameter_specs(lib, "A")
    with pytest.raises(ParameterError, match="shape"):
        from_theta(specs, np.zeros(len(specs) + 1))
    with pytest.raises(ParameterError, match="finite"):
        from_theta(specs, np.full(len(specs), np.inf))


# ---- overlay -----------------------------------------------------------------------------


def _quantities(model: object) -> dict[str, object]:
    return {path: q.value for path, q in iter_quantities(model)}


def test_override_changes_exactly_the_target(lib: Library) -> None:
    before = _data_hashes(lib.root)
    name = f"{TX}.latent.tau0_median"
    overlay = LibraryOverlay(lib, {name: 0.5})
    base_q = _quantities(lib.material("kodak/tri-x-400"))
    over_q = _quantities(overlay.material("kodak/tri-x-400"))
    changed = {p for p in base_q if base_q[p] != over_q[p]}
    assert changed == {"layers[1].latent.tau0_median"}
    assert over_q["layers[1].latent.tau0_median"] == 0.5
    assert parameter_value(overlay, name) == 0.5
    assert parameter_value(lib, name) == 0.22
    # Other files are the base's.
    assert overlay.material("kodak/t-max-100") == lib.material("kodak/t-max-100")
    assert overlay.kinetics() == lib.kinetics()
    assert overlay.bath("d-76") == lib.bath("d-76")
    # The rest of the quantity (bounds, prior, source, note) is unchanged.
    q = overlay.material("kodak/tri-x-400").layers[1].latent.tau0_median
    q0 = lib.material("kodak/tri-x-400").layers[1].latent.tau0_median
    assert q.model_dump(exclude={"value"}) == q0.model_dump(exclude={"value"})
    assert _data_hashes(lib.root) == before


def test_override_kinetics_and_bath(lib: Library) -> None:
    before = _data_hashes(lib.root)
    overrides = {
        f"{KIN}:agents.metol.k_et": 0.2,
        f"{KIN}:activity.delta_eps.borate": -0.1,
        f"{KIN}:agents.metol.activation_energy": 70.0,
        "developers/t-max-dev:components.dimezone-s.amount": 5.0,
    }
    overlay = LibraryOverlay(lib, overrides)
    kin, kin0 = overlay.kinetics(), lib.kinetics()
    changed = {p for p, v in _quantities(kin0).items() if _quantities(kin)[p] != v}
    assert changed == {
        "agents.metol.k_et",
        "activity.delta_eps.borate",
        "agents.metol.activation_energy",
    }
    for name, value in overrides.items():
        assert parameter_value(overlay, name) == value
    bath, bath0 = overlay.bath("t-max-dev"), lib.bath("t-max-dev")
    changed = {p for p, v in _quantities(bath0).items() if _quantities(bath)[p] != v}
    assert changed == {"components[4].amount"}
    assert bath.components[4].species == "dimezone-s"
    assert lib.kinetics().agents["metol"].k_et.value == 0.05
    assert _data_hashes(lib.root) == before


def test_overlay_from_theta(lib: Library) -> None:
    specs = parameter_specs(lib, "E")
    values = current_values(lib, specs)
    theta = to_theta(specs, values) + 0.01
    overlay = LibraryOverlay(lib, from_theta(specs, theta))
    assert current_values(overlay, specs) == pytest.approx(from_theta(specs, theta), rel=1e-15)
    assert current_values(lib, specs) == values


def test_aspect_ratio_sets_mean_diameter(lib: Library) -> None:
    ar = "calibration:kodak/t-max-100/aspect_ratio"
    t_g = f"{TM}.grain.thickness"
    overlay = LibraryOverlay(lib, {ar: 18.0, t_g: 0.15})
    grain = overlay.material("kodak/t-max-100").layers[1].grain
    assert grain.thickness.value == 0.15
    assert grain.mean_diameter.value == pytest.approx(18.0 * 0.15, rel=1e-15)
    # 2.7 um lies outside the [0.6, 1.5] penalty box, which is not a hard bound (spec D1 7.1).
    assert grain.mean_diameter.bounds == (0.6, 1.5)
    assert parameter_value(overlay, ar) == 18.0
    # Base thickness when only AR is given.
    grain = LibraryOverlay(lib, {ar: 10.0}).material("kodak/t-max-100").layers[1].grain
    assert grain.mean_diameter.value == pytest.approx(1.2, rel=1e-15)
    base = lib.material("kodak/t-max-100").layers[1].grain
    assert base.mean_diameter.value == 1.4794593607125541
    assert parameter_value(lib, ar) == pytest.approx(math.sqrt(152), rel=1e-12)


def test_aspect_ratio_conflicts_with_mean_diameter(lib: Library) -> None:
    with pytest.raises(ParameterError, match="both set"):
        LibraryOverlay(
            lib,
            {"calibration:kodak/t-max-100/aspect_ratio": 10.0, f"{TM}.grain.mean_diameter": 1.0},
        )


def test_s_abs_scales_absorption(lib: Library) -> None:
    alpha0 = f"{TX}.optics.intrinsic_absorption.alpha0"
    overlay = LibraryOverlay(lib, {"calibration:kodak/tri-x-400/s_abs": 0.3, alpha0: 2.0})
    optics = overlay.material("kodak/tri-x-400").layers[1].optics
    optics0 = lib.material("kodak/tri-x-400").layers[1].optics
    factor = 10.0**0.3
    assert optics.intrinsic_absorption.alpha0.value == pytest.approx(2.0 * factor, rel=1e-15)
    for band, band0 in zip(optics.dye_bands, optics0.dye_bands, strict=True):
        assert band.amplitude.value == pytest.approx(band0.amplitude.value * factor, rel=1e-15)
        assert band.peak_wavenumber == band0.peak_wavenumber
        assert band.width == band0.width
    assert optics.scattering_550 == optics0.scattering_550
    assert overlay.material("kodak/t-max-100") == lib.material("kodak/t-max-100")
    assert parameter_value(lib, "calibration:kodak/tri-x-400/s_abs") == 0.0


def test_calibration_values(lib: Library) -> None:
    overlay = LibraryOverlay(lib, {"calibration:delta_cc": -0.2})
    assert dict(overlay.calibration_values) == {"calibration:delta_cc": -0.2}
    assert parameter_value(overlay, "calibration:delta_cc") == -0.2
    assert parameter_value(lib, "calibration:delta_cc") == 0.0
    nested = LibraryOverlay(overlay, {f"{TX}.grain.sigma_log": 0.3})
    assert parameter_value(nested, "calibration:delta_cc") == -0.2
    assert nested.material("kodak/tri-x-400").layers[1].grain.sigma_log.value == 0.3


@pytest.mark.parametrize(
    "name",
    [
        "nonsense",
        "films/kodak/tri-x-400",
        "films/kodak/tri-x-400:",
        f"{TX}.latent.no_such_field",
        "films/kodak/tri-x-400:no-such-layer.latent.tau0_median",
        "films/kodak/no-such-film:emulsion.latent.tau0_median",
        "papers/kodak/tri-x-400:emulsion.latent.tau0_median",
        "chemistry/kinetics:agents.no-such-agent.k_et",
        "chemistry/kinetics:agents.metol",
        "chemistry/other:agents.metol.k_et",
        "developers/d-76:components.no-such-species.amount",
        f"{TX}.optics.dye_bands[5].amplitude",
        f"{TX}.optics.dye_bands[x].amplitude",
        f"{TX}.optics.intrinsic_absorption.alpha0.value",
        f"{KIN}:transport.water_viscosity.viscosity",
        "calibration:no_such",
        "calibration:s_abs",
        "calibration:kodak/tri-x-400/delta_cc",
        "calibration:kodak/no-such-film/s_abs",
        "calibration:kodak/tri-x-400/aspect_ratio",
    ],
)
def test_unknown_name_raises(lib: Library, name: str) -> None:
    with pytest.raises(ParameterError):
        LibraryOverlay(lib, {name: 1.0})


@pytest.mark.parametrize(
    ("name", "value"),
    [
        (f"{TX}.latent.tau0_median", 2.5),
        (f"{TX}.latent.tau0_median", 0.0),
        (f"{KIN}:activity.delta_eps.sulfite", 0.5),
        ("developers/d-76:components.boric-oxide.amount", -0.1),
        ("calibration:delta_cc", 0.51),
        ("calibration:kodak/tri-x-400/s_abs", -1.5),
        ("calibration:kodak/t-max-100/aspect_ratio", 7.9),
        (f"{TX}.latent.tau0_median", math.nan),
        (f"{TX}.latent.tau0_median", math.inf),
        (f"{TX}.latent.tau0_median", "0.5"),
        (f"{TX}.latent.tau0_median", True),
    ],
)
def test_out_of_bounds_or_bad_value_raises(lib: Library, name: str, value: object) -> None:
    with pytest.raises(ParameterError):
        LibraryOverlay(lib, {name: value})  # type: ignore[dict-item]


def test_overlay_delegates_every_read_to_base(lib: Library) -> None:
    class Base(Library):
        def reference_tables(self, material_id: str) -> dict:
            return {"from": "base", "id": material_id}

        def developer_facts(self, developer: str) -> dict:
            return {"from": "base", "developer": developer}

        def reference_curve(self, relpath: str) -> tuple[str, str]:
            return ("base", relpath)

    overlay = LibraryOverlay(Base(lib.root), {f"{TX}.grain.sigma_log": 0.3})
    assert overlay.reference_tables("kodak/tri-x-400") == {"from": "base", "id": "kodak/tri-x-400"}
    assert overlay.developer_facts("d-76") == {"from": "base", "developer": "d-76"}
    assert overlay.reference_curve("x.csv") == ("base", "x.csv")
    assert overlay.root == lib.root
    assert overlay.ids("material") == lib.ids("material")
    assert overlay.species() == lib.species()


def test_overlay_needs_a_library() -> None:
    with pytest.raises(TypeError):
        LibraryOverlay("data", {})  # type: ignore[arg-type]
