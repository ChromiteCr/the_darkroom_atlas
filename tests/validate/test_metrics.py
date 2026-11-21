"""Sensitometric metrics (spec D1 6.3-6.5, 6.9 step 3, 6.10, 8.2, 12.9)."""

import math

import numpy as np
import pytest
from scipy.interpolate import PchipInterpolator

from darkroom_atlas.data.models import Agitation, Tank
from darkroom_atlas.validate.conditions import Condition
from darkroom_atlas.validate.metrics import (
    MetricValue,
    SpeedResult,
    average_gradient,
    average_gradient_detail,
    chord_contrast_index,
    chord_contrast_index_detail,
    ci_at_times,
    contrast_index,
    contrast_index_detail,
    d_min_star,
    energy_exposure_erg_cm2,
    iso_bin,
    iso_speed,
    speed_point,
    time_to_ci,
    time_to_ci_detail,
)
from darkroom_atlas.validate.sensitometer import StripResult


def _condition(times_s: tuple[float, ...], material: str = "kodak/tri-x-400") -> Condition:
    return Condition(
        material_id=material,
        film_format="135",
        developer="d-76",
        dilution="stock",
        temperature_c=20.0,
        agitation=Agitation(method="inversion", initial_s=30.0, interval_s=60.0, duration_s=5.0),
        times_s=times_s,
        tank=Tank(volume_ml=3785.0),
    )


def _hinge(x: np.ndarray, d_min: float, gradient: float, x0: float) -> np.ndarray:
    return d_min + gradient * np.maximum(x - x0, 0.0)


def _strip_from_curves(curves, times_s, material="kodak/tri-x-400") -> StripResult:
    """Resample reference curves onto their common log H range (step 0.005)."""
    lo = max(float(x[0]) for x, _ in curves)
    hi = min(float(x[-1]) for x, _ in curves)
    grid = np.append(np.arange(lo, hi, 0.005), hi)
    visual = np.array([PchipInterpolator(x, y)(grid) for x, y in curves])
    fog = np.array([float(y[0]) for _, y in curves])
    return StripResult(
        condition=_condition(tuple(times_s), material),
        log_h=grid,
        times_s=np.asarray(times_s, dtype=float),
        visual=visual,
        type2=visual.copy(),
        base_fog_visual=fog,
        base_fog_type2=fog,
    )


# ---- 6.3 ISO speed ---------------------------------------------------------------------


@pytest.mark.parametrize(
    ("fmt", "expected_log_h_m", "expected_delta_d"),
    [
        ("35mm", -2.73, (0.585, 0.714, 0.852, 1.017)),
        ("120", -2.69, (0.569, 0.690, 0.817, 0.964)),
    ],
)
def test_iso_speed_on_kodak_tri_x(read_reference, fmt, expected_log_h_m, expected_delta_d):
    minutes = (6, 8, 10, 12)
    curves = [read_reference(f"tri-x-400/char-{fmt}-d76-{m:02d}min.csv")[1:] for m in minutes]
    d_min_star = np.array([float(y[0]) for _, y in curves])
    # Spec 12.9: delta D per curve, read on each raw curve with its own left end.
    for (x, y), expected in zip(curves, expected_delta_d, strict=True):
        assert speed_point(x, y, float(y[0]))[1] == pytest.approx(expected, abs=0.01)
    strip = _strip_from_curves(curves, [60.0 * m for m in minutes])
    result = iso_speed(strip, density="visual", base_fog=d_min_star)
    assert result.reason == ""
    assert result.log_h_m == pytest.approx(expected_log_h_m, abs=0.03)
    assert result.delta_d == pytest.approx(0.80)
    assert result.speed == pytest.approx(0.80 / 10.0**result.log_h_m)
    assert result.speed_degrees == pytest.approx(1.0 + 10.0 * math.log10(result.speed))
    assert result.iso_bin == "400"


def test_iso_speed_default_uses_type2_and_zero_exposure_patch():
    x = np.linspace(-4.0, 0.0, 401)
    gradients = (0.5, 0.7)
    t2 = np.array([_hinge(x, 0.05, g, -3.0) for g in gradients])
    strip = StripResult(
        condition=_condition((300.0, 600.0)),
        log_h=x,
        times_s=np.array([300.0, 600.0]),
        visual=t2 + 0.5,  # would give another answer if used
        type2=t2,
        base_fog_visual=np.full(2, 0.55),
        base_fog_type2=np.full(2, 0.05),
    )
    result = iso_speed(strip)
    # Speed point 0.10 above fog: x_m = -3 + 0.1/G; delta D = 1.30 G crosses 0.80 at G = 8/13.
    w = (0.8 / 1.3 - 0.5) / 0.2
    expected = (1 - w) * (-3.0 + 0.1 / 0.5) + w * (-3.0 + 0.1 / 0.7)
    assert result.log_h_m == pytest.approx(expected, abs=1e-6)
    assert result.delta_d == pytest.approx(0.80)


def test_iso_speed_nan_when_not_bracketed():
    x = np.linspace(-4.0, 0.0, 401)
    curves = np.array([_hinge(x, 0.05, g, -3.0) for g in (0.3, 0.4)])
    strip = StripResult(
        condition=_condition((300.0, 600.0)),
        log_h=x,
        times_s=np.array([300.0, 600.0]),
        visual=curves,
        type2=curves,
        base_fog_visual=np.full(2, 0.05),
        base_fog_type2=np.full(2, 0.05),
    )
    result = iso_speed(strip)
    assert isinstance(result, SpeedResult)
    assert math.isnan(result.log_h_m) and math.isnan(result.speed)
    assert result.iso_bin is None
    assert "bracket" in result.reason


def test_iso_speed_rejects_bad_options():
    x = np.linspace(-4.0, 0.0, 41)
    c = _hinge(x, 0.05, 0.6, -3.0)[None, :]
    strip = StripResult(_condition((300.0,)), x, [300.0], c, c, [0.05], [0.05])
    with pytest.raises(ValueError, match="density"):
        iso_speed(strip, density="blue")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="base_fog"):
        iso_speed(strip, base_fog=np.array([0.1, 0.2]))


def test_speed_point_nan_when_curve_starts_above_point_or_ends_early():
    x = np.linspace(-3.0, 0.0, 31)
    y = _hinge(x, 0.3, 0.6, -2.5)
    assert all(math.isnan(v) for v in speed_point(x, y, 0.1))  # already above 0.2
    x_short = np.linspace(-3.0, -1.5, 16)
    x_m, delta_d = speed_point(x_short, _hinge(x_short, 0.1, 0.6, -2.5), 0.1)
    assert math.isfinite(x_m) and math.isnan(delta_d)


def test_iso_bin_matches_table_anchors():
    # Spec 6.3 step 6: ISO 400 is -2.75 ... -2.66, ISO 100 is -2.15 ... -2.06.
    for v in (-2.75, -2.70, -2.66):
        assert iso_bin(v) == "400"
    assert iso_bin(-2.65) == "320"
    assert iso_bin(-2.76) == "500"
    for v in (-2.15, -2.10, -2.06):
        assert iso_bin(v) == "100"
    assert iso_bin(-2.05) == "80"
    assert iso_bin(-2.16) == "125"
    assert iso_bin(math.nan) is None
    assert iso_bin(-20.0) is None


# ---- 6.4 average gradient --------------------------------------------------------------


def test_average_gradient_on_straight_line():
    x = np.linspace(-4.0, 1.0, 501)
    assert average_gradient(x, _hinge(x, 0.1, 0.62, -3.0), 0.1) == pytest.approx(0.62, abs=1e-9)


# ---- 6.5 contrast index ----------------------------------------------------------------


@pytest.mark.parametrize("gradient", [0.45, 0.62, 0.9, 1.4])
def test_ci_straight_line(gradient):
    x = np.round(np.arange(-4.0, 1.0 + 1e-9, 0.01), 10)
    y = _hinge(x, 0.25, gradient, -2.5)
    assert contrast_index(x, y, 0.25) == pytest.approx(gradient, abs=1e-9)


@pytest.mark.parametrize(
    ("minutes", "ci", "g_iso", "log_h_m"),
    [
        ("06", 0.542, 0.547, -2.011),
        ("07.5", 0.693, 0.701, -2.129),
        ("10", 0.895, 0.879, -2.212),
    ],
)
def test_ci_matches_report_on_tmax_d76(read_reference, minutes, ci, g_iso, log_h_m):
    _, x, y = read_reference(f"t-max-100/char-d76-{minutes}min.csv")
    d_min = float(y[0])
    assert contrast_index(x, y, d_min) == pytest.approx(ci, abs=0.01)
    assert average_gradient(x, y, d_min) == pytest.approx(g_iso, abs=0.01)
    assert speed_point(x, y, d_min)[0] == pytest.approx(log_h_m, abs=0.02)
    # The chord cross-check agrees with the collinear construction to a few hundredths.
    assert chord_contrast_index(x, y, d_min) == pytest.approx(ci, abs=0.05)


def test_ci_nan_on_short_curve():
    # Neither the log H span nor the distance from any centre reaches 2.2.
    x = np.linspace(-3.0, -1.5, 151)
    assert math.isnan(contrast_index(x, _hinge(x, 0.2, 0.3, -2.8), 0.2))
    x2 = np.linspace(-3.0, -0.9, 211)
    assert math.isnan(contrast_index(x2, np.full_like(x2, 0.2), 0.2))


def test_ci_nan_on_flat_curve():
    x = np.linspace(-4.0, 1.0, 501)
    assert math.isnan(contrast_index(x, np.full_like(x, 0.2), 0.2))


def test_chord_ci_on_straight_line_and_short_curve():
    x = np.linspace(-4.0, 1.0, 501)
    assert chord_contrast_index(x, _hinge(x, 0.2, 0.7, -3.0), 0.2) == pytest.approx(0.7, abs=1e-9)
    x_short = np.linspace(-4.0, -2.5, 151)
    assert math.isnan(chord_contrast_index(x_short, _hinge(x_short, 0.2, 0.7, -3.5), 0.2))


@pytest.mark.parametrize(
    "fn", [contrast_index, chord_contrast_index, average_gradient], ids=lambda f: f.__name__
)
def test_curve_metrics_reject_bad_input(fn):
    x = np.linspace(-3.0, 0.0, 11)
    y = np.linspace(0.2, 2.0, 11)
    with pytest.raises(ValueError, match="same length"):
        fn(x, y[:-1], 0.2)
    with pytest.raises(ValueError, match="increasing"):
        fn(x[::-1], y, 0.2)
    with pytest.raises(ValueError, match="log_h must be finite"):
        fn(np.where(x > -1, np.nan, x), y, 0.2)
    with pytest.raises(ValueError, match="at least 2"):
        fn(x[:1], y[:1], 0.2)


@pytest.mark.parametrize(
    "fn",
    [contrast_index_detail, chord_contrast_index_detail, average_gradient_detail],
    ids=lambda f: f.__name__,
)
def test_curve_metrics_nan_with_reason_on_undefined_input(fn):
    # Spec 8.2: a failed strip (non-finite densities) or an undefined D_min gives NaN with a
    # reason, never an exception, so the metric is counted as failed rather than skipped.
    x = np.linspace(-4.0, 1.0, 501)
    y = _hinge(x, 0.2, 0.7, -3.0)
    result = fn(x, np.where(x > -1, np.nan, y), 0.2)
    assert isinstance(result, MetricValue)
    assert math.isnan(result.value) and "non-finite" in result.reason
    result = fn(x, y, math.nan)
    assert math.isnan(result.value) and "not finite" in result.reason
    good = fn(x, y, 0.2)
    assert math.isfinite(good.value) and good.reason == ""


def test_ci_detail_reason_on_short_curve():
    x = np.linspace(-3.0, -1.5, 151)
    result = contrast_index_detail(x, _hinge(x, 0.2, 0.3, -2.8), 0.2)
    assert math.isnan(result.value) and "collinear" in result.reason
    result = chord_contrast_index_detail(x, _hinge(x, 0.2, 0.3, -2.8), 0.2)
    assert math.isnan(result.value) and "2.0 log H" in result.reason


# ---- 6.3 step 2 D_min* -----------------------------------------------------------------


def _sim_strip(times_s, x_lo, x_hi, fog_patch):
    """A smooth synthetic strip on [x_lo, x_hi] whose unexposed patch differs from D(x_lo)."""
    x = np.append(np.arange(x_lo, x_hi, 0.01), x_hi)
    curves = []
    for g in (0.48, 0.58, 0.68, 0.78):
        w = 0.25
        curves.append(0.24 + g * w * np.logaddexp(0.0, (x + 2.9) / w))
    curves = np.array(curves)
    fog = np.full(len(times_s), fog_patch)
    return StripResult(
        condition=_condition(tuple(times_s)),
        log_h=x,
        times_s=np.asarray(times_s, dtype=float),
        visual=curves,
        type2=curves.copy(),
        base_fog_visual=fog,
        base_fog_type2=fog,
    )


def test_dmin_star_same_operator(read_reference):
    minutes = (6, 8, 10, 12)
    times = [60.0 * m for m in minutes]
    refs = [read_reference(f"tri-x-400/char-35mm-d76-{m:02d}min.csv")[1:] for m in minutes]
    x_min = min(float(x[0]) for x, _ in refs)
    x_max = max(float(x[-1]) for x, _ in refs)
    # Spec 6.1: the strip spans the reference range widened by 0.30 on each side.
    sim = _sim_strip(times, x_min - 0.30, x_max + 0.30, fog_patch=0.31)

    def two_sided(shift):
        shifted = [(x, y + shift) for x, y in refs]
        stars = [d_min_star(x, y, sim.log_h, sim.visual[i]) for i, (x, y) in enumerate(shifted)]
        ref_star = np.array([s[0] for s in stars])
        sim_star = np.array([s[1] for s in stars])
        ref_speed = iso_speed(
            _strip_from_curves(shifted, times), density="visual", base_fog=ref_star
        )
        sim_speed = iso_speed(sim, density="visual", base_fog=sim_star)
        ci_diff = [
            contrast_index(x, y, ref_star[i])
            - contrast_index(sim.log_h, sim.visual[i], sim_star[i])
            for i, (x, y) in enumerate(shifted)
        ]
        return ref_star, sim_star, ref_speed.log_h_m - sim_speed.log_h_m, np.array(ci_diff)

    ref_star, sim_star, d_log_h_m, d_ci = two_sided(0.0)
    for i, (x, y) in enumerate(refs):
        # Both sides are read at the reference's leftmost log H, not at the zero patch.
        assert ref_star[i] == y[0]
        expected = float(PchipInterpolator(sim.log_h, sim.visual[i])(x[0]))
        assert sim_star[i] == pytest.approx(expected, abs=1e-15)
        assert abs(sim_star[i] - sim.base_fog_visual[i]) > 0.01
    assert math.isfinite(d_log_h_m) and np.all(np.isfinite(d_ci))

    ref_star2, sim_star2, d_log_h_m2, d_ci2 = two_sided(0.02)
    np.testing.assert_allclose(ref_star2, ref_star + 0.02, rtol=0, atol=1e-12)
    np.testing.assert_array_equal(sim_star2, sim_star)
    assert abs(d_log_h_m2 - d_log_h_m) <= 1e-9
    np.testing.assert_allclose(d_ci2, d_ci, rtol=0, atol=1e-9)


def test_dmin_star_nan_outside_simulated_range_and_bad_input():
    x_ref = np.linspace(-3.0, 0.0, 31)
    y_ref = _hinge(x_ref, 0.2, 0.6, -2.5)
    x_sim = np.linspace(-2.5, 0.5, 31)
    ref, sim = d_min_star(x_ref, y_ref, x_sim, _hinge(x_sim, 0.2, 0.6, -2.5))
    assert ref == y_ref[0] and math.isnan(sim)
    x_sim = np.linspace(-3.5, 0.5, 41)
    y_sim = _hinge(x_sim, 0.2, 0.6, -2.5)
    assert math.isnan(d_min_star(x_ref, y_ref, x_sim, np.full_like(x_sim, np.nan))[1])
    with pytest.raises(ValueError, match="density_ref"):
        d_min_star(x_ref, np.where(x_ref > -1, np.nan, y_ref), x_sim, y_sim)
    with pytest.raises(ValueError, match="same length"):
        d_min_star(x_ref, y_ref[:-1], x_sim, y_sim)
    with pytest.raises(ValueError, match="increasing"):
        d_min_star(x_ref, y_ref, x_sim[::-1], y_sim)


def test_iso_speed_nan_with_reason_on_failed_strip():
    x = np.linspace(-4.0, 0.0, 401)
    curves = np.full((2, x.size), np.nan)
    strip = StripResult(
        _condition((300.0, 600.0)), x, [300.0, 600.0], curves, curves, [0.05, 0.05], [0.05, 0.05]
    )
    result = iso_speed(strip)
    assert math.isnan(result.log_h_m) and "non-finite" in result.reason
    good = np.array([_hinge(x, 0.05, g, -3.0) for g in (0.5, 0.7)])
    strip = StripResult(
        _condition((300.0, 600.0)), x, [300.0, 600.0], good, good, [0.05, 0.05], [0.05, 0.05]
    )
    result = iso_speed(strip, base_fog=np.array([math.nan, math.nan]))
    assert math.isnan(result.log_h_m) and "base_fog is not finite" in result.reason


# ---- 6.10 time to CI -------------------------------------------------------------------


def test_time_to_ci_interpolates_monotone_series():
    t = np.array([300.0, 420.0, 540.0, 660.0])
    ci = 0.30 + 0.001 * t
    assert time_to_ci(t, ci, 0.80) == pytest.approx(500.0, abs=1e-9)
    assert time_to_ci(t, ci, 0.60) == pytest.approx(300.0)
    assert time_to_ci(t, ci, 0.96) == pytest.approx(660.0, abs=1e-9)


def test_time_to_ci_nan_outside_range_and_skips_nan():
    t = np.array([300.0, 420.0, 540.0, 660.0])
    ci = np.array([math.nan, 0.5, 0.6, 0.7])
    assert math.isnan(time_to_ci(t, ci, 0.75))
    assert math.isnan(time_to_ci(t, ci, 0.40))
    assert time_to_ci(t, ci, 0.55) == pytest.approx(480.0, abs=1e-6)
    assert math.isnan(time_to_ci(t, np.full(4, math.nan), 0.5))


def test_time_to_ci_rejects_bad_input():
    with pytest.raises(ValueError, match="increasing"):
        time_to_ci(np.array([2.0, 1.0]), np.array([0.5, 0.6]), 0.55)
    with pytest.raises(ValueError, match="same length"):
        time_to_ci(np.array([1.0, 2.0]), np.array([0.5]), 0.55)
    with pytest.raises(ValueError, match="target"):
        time_to_ci(np.array([1.0, 2.0]), np.array([0.5, 0.6]), math.nan)


def test_time_to_ci_does_not_bridge_undefined_values():
    t = np.array([300.0, 420.0, 540.0, 660.0, 780.0])
    ci = np.array([0.40, 0.50, math.nan, 0.70, 0.80])
    result = time_to_ci_detail(t, ci, 0.60)
    assert math.isnan(result.value) and "undefined" in result.reason
    assert time_to_ci(t, ci, 0.75) == pytest.approx(720.0, abs=1e-6)
    assert time_to_ci(t, ci, 0.45) == pytest.approx(360.0, abs=1e-6)
    assert "does not reach" in time_to_ci_detail(t, ci, 0.9).reason
    assert "above" in time_to_ci_detail(t, ci, 0.3).reason
    assert "no finite" in time_to_ci_detail(t, np.full(5, math.nan), 0.5).reason


# ---- 6.10 / 7.2 CI at tabulated times --------------------------------------------------


def test_ci_at_times_interpolates_monotone_series():
    t = np.arange(0.0, 900.0 + 1.0, 15.0)
    ci = 0.9 * (1.0 - np.exp(-t / 400.0))
    at = np.array([[15.0, 397.5], [450.0, 900.0]])
    out = ci_at_times(t, ci, at)
    assert out.shape == at.shape and out.dtype == np.float64
    assert not out.flags.writeable
    assert out[0, 0] == pytest.approx(ci[1], abs=1e-15)
    assert out[1, 1] == pytest.approx(ci[-1], abs=1e-15)
    np.testing.assert_allclose(out, 0.9 * (1.0 - np.exp(-at / 400.0)), rtol=0, atol=1e-5)
    # Monotone (PCHIP keeps the data's monotonicity) on a fine grid.
    fine = ci_at_times(t, ci, np.linspace(0.0, 900.0, 5001))
    assert np.all(np.diff(fine) >= 0.0)
    # Linear data is reproduced exactly.
    lin = ci_at_times(t, 0.2 + 0.001 * t, np.array([7.0, 333.3]))
    np.testing.assert_allclose(lin, [0.207, 0.5333], rtol=0, atol=1e-12)
    assert ci_at_times(t, ci, 400.0).shape == ()


def test_ci_at_times_nan_outside_and_across_undefined_values():
    t = np.array([300.0, 420.0, 540.0, 660.0, 780.0])
    ci = np.array([0.40, 0.50, math.nan, 0.70, 0.80])
    out = ci_at_times(t, ci, np.array([299.0, 360.0, 420.0, 480.0, 540.0, 600.0, 700.0, 781.0]))
    assert math.isnan(out[0]) and math.isnan(out[-1])
    assert out[1] == pytest.approx(0.45) and out[2] == pytest.approx(0.50)
    assert all(math.isnan(v) for v in out[3:6])
    assert out[6] == pytest.approx(0.70 + 0.1 * 40.0 / 120.0)
    single = ci_at_times(t, np.array([math.nan, 0.5, math.nan, 0.7, 0.8]), np.array([420.0, 430.0]))
    assert single[0] == 0.5 and math.isnan(single[1])


def test_ci_at_times_rejects_bad_input():
    t = np.array([300.0, 420.0])
    with pytest.raises(ValueError, match="increasing"):
        ci_at_times(t[::-1], np.array([0.5, 0.6]), np.array([350.0]))
    with pytest.raises(ValueError, match="same length"):
        ci_at_times(t, np.array([0.5]), np.array([350.0]))
    with pytest.raises(ValueError, match="at_s"):
        ci_at_times(t, np.array([0.5, 0.6]), np.array([math.nan]))
    with pytest.raises(ValueError, match="empty"):
        ci_at_times(np.array([]), np.array([]), np.array([1.0]))


# ---- 6.9 step 3 energy conversion ------------------------------------------------------


def test_spectral_energy_conversion():
    # Spec 6.9 step 3 (LAT 3.3): at 550 nm, 1 erg/cm2 = 2769 photons/um2, relative 1e-3.
    photons_per_erg_cm2 = 1.0 / energy_exposure_erg_cm2(1.0, 550.0)
    assert photons_per_erg_cm2 == pytest.approx(2769.0, rel=1e-3)
    assert energy_exposure_erg_cm2(2769.0, 550.0) == pytest.approx(1.0, rel=1e-3)
    out = energy_exposure_erg_cm2(np.array([1.0, 2.0]), np.array([400.0, 400.0]))
    assert out.shape == (2,) and out[1] == pytest.approx(2.0 * out[0])
    assert energy_exposure_erg_cm2(1.0, 400.0) == pytest.approx(550.0 / 400.0 / 2768.8, rel=1e-3)


# ---- StripResult -----------------------------------------------------------------------


def test_strip_result_is_read_only_and_checks_shapes():
    x = np.linspace(-3.0, 0.0, 5)
    d = np.zeros((2, 5))
    strip = StripResult(_condition((60.0, 120.0)), x, [60.0, 120.0], d, d, [0.0, 0.0], [0.0, 0.0])
    assert not strip.visual.flags.writeable
    assert strip.visual.dtype == np.float64
    with pytest.raises(ValueError, match="visual"):
        StripResult(_condition((60.0, 120.0)), x, [60.0, 120.0], d[:1], d, [0, 0], [0, 0])
    with pytest.raises(ValueError, match="base_fog_type2"):
        StripResult(_condition((60.0, 120.0)), x, [60.0, 120.0], d, d, [0, 0], [0])
    with pytest.raises(ValueError, match="log_h"):
        StripResult(_condition((60.0, 120.0)), x[::-1], [60.0, 120.0], d, d, [0, 0], [0, 0])
    with pytest.raises(ValueError, match="times_s"):
        StripResult(_condition((60.0, 120.0)), x, [120.0, 60.0], d, d, [0, 0], [0, 0])
    with pytest.raises(ValueError, match="times_s"):
        StripResult(_condition((60.0, 120.0)), x, [60.0, math.nan], d, d, [0, 0], [0, 0])
