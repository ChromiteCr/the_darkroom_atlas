"""Sensitometric metrics (spec D1 6.3-6.5, 6.10, 12.9)."""

import math

import numpy as np
import pytest
from scipy.interpolate import PchipInterpolator

from darkroom_atlas.data.models import Agitation, Tank
from darkroom_atlas.validate.conditions import Condition
from darkroom_atlas.validate.metrics import (
    SpeedResult,
    average_gradient,
    chord_contrast_index,
    contrast_index,
    iso_bin,
    iso_speed,
    speed_point,
    time_to_ci,
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
    with pytest.raises(ValueError, match="finite"):
        fn(x, np.where(x > -1, np.nan, y), 0.2)
    with pytest.raises(ValueError, match="at least 2"):
        fn(x[:1], y[:1], 0.2)
    with pytest.raises(ValueError, match="finite"):
        fn(x, y, math.nan)


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
