"""Reference observations for calibration and validation (spec D1 7.2, 7.6, 7.7).

An :class:`Observation` is one item's reference values with the sigma and weights of spec 7.2
(``sigma`` already holds the sqrt(N_pts/N_eff), sigma_floor and sqrt(n_class) inflations).
:func:`training_observations` reads only through a
:class:`~darkroom_atlas.validate.splits.TrainingLibrary` and only items of ``TRAINING``;
every id it returns passes :func:`~darkroom_atlas.validate.splits.assert_trainable`.
:func:`within_check_observations` gives the within-pair checks (never in the objective) and
:func:`held_out_observations` the held-out items, only for a frozen calibration record.

Item ids are the split registry's. An observation built from several files (the speed point
of a development series) joins their ids with ``" + "`` (:func:`observation_sources`). A row
of a developer publication (``developers/*.yaml``) is one registry item but gives one
observation per cell: those share ``item_id`` and differ by ``condition.temperature_c``. pH
anchors come from the bath files, outside the reference registry, as
``bath:<bath id>#ph_anchors/<dilution>``; a colorimetric (Carlton & Crabtree) anchor carries
the suffix ``#delta_cc`` (:func:`takes_delta_cc`): its residual is pH_sim + delta_CC - pH_ref
with one offset for the whole series (spec 7.2). The bath files hold one such anchor so far
(``d-76-published`` 8.5); with one member delta_CC absorbs it exactly, so the series gives no
information until the other Carlton & Crabtree variants (8.1, 8.9, 9.2, 8.1, 8.6) are anchors
too. XTOL Part A (10.19, sigma 0.3) is not a dilution of a bath and is not an anchor yet.

Held-out observations carry ``sigma = NaN`` and ``in_objective = False``: they are scored
against the pre-registered thresholds of spec 8.3, never fitted.
"""

import dataclasses
import itertools
import math
from collections import Counter
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import Any, Literal

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

from darkroom_atlas._frozen import readonly_array
from darkroom_atlas.data.library import Library
from darkroom_atlas.data.models import CalibrationRecord, ReferenceCurve
from darkroom_atlas.validate.conditions import (
    Condition,
    ConditionedItem,
    reference_conditions,
)
from darkroom_atlas.validate.metrics import ci_at_times, contrast_index, d_min_star, speed_point
from darkroom_atlas.validate.splits import (
    LeakageError,
    TrainingLibrary,
    assert_trainable,
    classify,
    enumerate_items,
)

__all__ = [
    "CI_N_EFFECTIVE",
    "DENSITY_N_EFFECTIVE",
    "KINDS",
    "PH_PREFIX",
    "Observation",
    "ObservationKind",
    "ci_sigma",
    "density_sigma",
    "held_out_observations",
    "observation_counts",
    "observation_sources",
    "sigma_floor",
    "table_ci_sigma",
    "takes_delta_cc",
    "training_observations",
    "within_check_observations",
]

ObservationKind = Literal[
    "density", "ci", "table_ci", "speed", "ph", "spectral", "mtf", "granularity", "reciprocity"
]
KINDS: tuple[str, ...] = (
    "density",
    "ci",
    "table_ci",
    "speed",
    "ph",
    "spectral",
    "mtf",
    "granularity",
    "reciprocity",
)
PH_PREFIX = "bath:"
DELTA_CC_SUFFIX = "#delta_cc"
_JOIN = " + "

_TX, _TM = "kodak/tri-x-400", "kodak/t-max-100"
_FILMS = (_TX, _TM)

# ---- spec 7.2 constants ------------------------------------------------------------------
# Characteristic-curve density: sigma 0.025 (CAL 2.3) times sqrt(N_pts/N_eff), 15 resampled
# points from D_ref = D_min* + 0.05 to the last point, about 6 independent numbers per curve
# (7.5); points with D < D_min* + 0.1 weigh 0.5.
DENSITY_SIGMA_0 = 0.025
DENSITY_POINTS = 15
DENSITY_N_EFFECTIVE = 6.0
DENSITY_START_RISE = 0.05
DENSITY_LOW_RISE = 0.10
DENSITY_LOW_WEIGHT = 0.5
# CI on a CI(t) curve: sqrt(0.025^2 + sigma_floor^2) * sqrt(5/3), 5 points over the span,
# about 3 independent numbers (CAL 3.2).
CI_SIGMA_0 = 0.025
CI_POINTS = 5
CI_N_EFFECTIVE = 3.0
# sigma_floor (7.2): Tri-X 0 (assumed); T-Max half the CI difference between the 7.5 min
# D-76 characteristic curve and the small-tank D-76 CI(t) curve at 7.5 min, recomputed here
# from the two training CSVs (7.6: the fact item is not read).
SIGMA_FLOOR_TRIX = 0.0
_FLOOR_CHAR = f"{_TM}/char-d76-07.5min.csv"
_FLOOR_CI = f"{_TM}/ci-small-tank-d76.csv"
_FLOOR_TIME_S = 7.5 * 60.0
# Table cell CI (7.2): sqrt(sigma_tab^2 + (s_ref * 0.25/sqrt(12) min)^2) * sqrt(n_class);
# sigma_tab = 0.02 (assumed); XTOL stock non-20 C cells times 2 (KIN 1.6).
SIGMA_TAB = 0.02
ROUNDING_SD_MIN = 0.25 / math.sqrt(12.0)
XTOL_STOCK_FACTOR = 2.0
_REFERENCE_F = 68  # the 20 C column of the Kodak tables (spec 6.1)
# pH anchors (7.2): the anchor's own uncertainty; colorimetric (Carlton & Crabtree) 0.10.
PH_SIGMA_COLORIMETRIC = 0.10
# Spectral shape (7.2, 6.9 step 4): 0.03 log S at 400-700 nm every 20 nm inside the reference
# span, each curve minus its maximum over 450-600 nm.
SPECTRAL_SIGMA = 0.03
SPECTRAL_BANDS_NM = tuple(float(b) for b in range(400, 701, 20))
SPECTRAL_SHAPE_WINDOW_NM = (450.0, 600.0)
SPECTRAL_REPORT_STEP_NM = 5.0  # 6.9 step 4: the report uses every 5 nm in the span
# MTF tail (7.2, 6.7 steps 6, 8): 2 percentage points, Tri-X at 50 c/mm and above 3.
MTF_SIGMA = 2.0
MTF_SIGMA_TRIX_HIGH = 3.0
MTF_TRIX_HIGH_C_MM = 50.0
MTF_FREQUENCIES_C_MM = (2.5, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0)
MTF_FREQUENCIES_TMAX_EXTRA = (80.0, 100.0, 120.0, 150.0)
MTF_FULL_RESPONSE = 100.0
# RMS granularity (7.2, 6.6): G_exp - 8 for T-Max at net density 1.0, sigma 0.6.
GRANULARITY_SIGMA = 0.6
GRANULARITY_NET_DENSITY = 1.0
# Reciprocity (7.2): difference in stops, sigma 0.15.
RECIPROCITY_SIGMA = 0.15
# CI_aim of the ISO speed series (6.3 step 5).
_ISO_DELTA_D = 0.80
_XTOL = 1e-12


@dataclass(frozen=True, eq=False)
class Observation:
    """One reference item's values, sigma and weights (spec 7.2, 10 D1e).

    ``x``: log H (density), times in s (ci, table_ci, speed), nm (spectral), c/mm (mtf),
    net density (granularity), indicated exposure time in s (reciprocity), temperature in C
    (ph). ``y``: density, CI (table_ci: CI_aim), log H_m* (speed), pH, log S shape, MTF %,
    granularity, stops.
    """

    kind: ObservationKind
    item_id: str
    condition: Condition | None
    x: np.ndarray
    y: np.ndarray
    sigma: np.ndarray
    weight: np.ndarray
    n_effective: float
    d_min_ref: float | None = None
    in_objective: bool = True

    def __post_init__(self) -> None:
        if self.kind not in KINDS:
            raise ValueError(f"Observation.kind must be one of {KINDS}, got {self.kind!r}")
        if not self.item_id:
            raise ValueError("Observation.item_id must not be empty")
        for name in ("x", "y", "sigma", "weight"):
            value = np.atleast_1d(np.asarray(getattr(self, name), dtype=np.float64))
            object.__setattr__(self, name, readonly_array(value, np.float64))
        n = self.x.size
        for name in ("x", "y", "sigma", "weight"):
            arr = getattr(self, name)
            if arr.ndim != 1 or arr.size != n or n == 0:
                raise ValueError(
                    f"Observation.{name} must be 1-D with {n} > 0 values, got {arr.shape}"
                )
        if not (np.all(np.isfinite(self.x)) and np.all(np.isfinite(self.y))):
            raise ValueError(f"{self.item_id}: Observation x and y must be finite")
        if not (np.all(np.isfinite(self.weight)) and np.all(self.weight >= 0.0)):
            raise ValueError(f"{self.item_id}: Observation.weight must be finite and >= 0")
        if self.in_objective and not (np.all(np.isfinite(self.sigma)) and np.all(self.sigma > 0.0)):
            raise ValueError(f"{self.item_id}: an objective term needs finite sigma > 0")
        if not (math.isfinite(self.n_effective) and self.n_effective > 0.0):
            raise ValueError(f"{self.item_id}: n_effective must be positive")


def takes_delta_cc(obs: Observation) -> bool:
    """True for a Carlton & Crabtree pH anchor, whose residual adds delta_CC (spec 7.2)."""
    return obs.kind == "ph" and obs.item_id.endswith(DELTA_CC_SUFFIX)


def observation_sources(obs: Observation) -> tuple[str, ...]:
    """Registry ids an observation was built from; ``()`` for pH anchors (bath files)."""
    if obs.item_id.startswith(PH_PREFIX):
        return ()
    return tuple(obs.item_id.split(_JOIN))


def observation_counts(observations: Iterable[Observation]) -> dict[tuple[str, str], int]:
    """``{(material or "bath", kind): count}`` for a summary."""
    counts: Counter[tuple[str, str]] = Counter()
    for obs in observations:
        owner = obs.condition.material_id if obs.condition is not None else "bath"
        counts[(owner, obs.kind)] += 1
    return dict(sorted(counts.items()))


# ---- sigma (spec 7.2) ----------------------------------------------------------------------


def density_sigma() -> float:
    """0.025 * sqrt(15/6) (spec 7.2)."""
    return DENSITY_SIGMA_0 * math.sqrt(DENSITY_POINTS / DENSITY_N_EFFECTIVE)


def ci_sigma(floor: float) -> float:
    """sqrt(0.025^2 + floor^2) * sqrt(5/3) (spec 7.2)."""
    return math.hypot(CI_SIGMA_0, floor) * math.sqrt(CI_POINTS / CI_N_EFFECTIVE)


def table_ci_sigma(slope_per_min: float, n_class: int, *, xtol_stock_off_20c: bool) -> float:
    """sqrt(0.02^2 + (s_ref * 0.072 min)^2) * sqrt(n_class), XTOL stock off 20 C x 2."""
    if n_class < 1:
        raise ValueError(f"n_class must be >= 1, got {n_class}")
    sigma = math.hypot(SIGMA_TAB, slope_per_min * ROUNDING_SD_MIN) * math.sqrt(n_class)
    return sigma * (XTOL_STOCK_FACTOR if xtol_stock_off_20c else 1.0)


def sigma_floor(material_id: str, library: Library) -> float:
    """sigma_floor of spec 7.2: 0 for Tri-X, recomputed from training CSVs for T-Max."""
    if material_id == _TX:
        return SIGMA_FLOOR_TRIX
    if material_id != _TM:
        raise ValueError(f"no sigma_floor for {material_id!r} (spec 7.2)")
    char = library.reference_curve(_FLOOR_CHAR)
    ci_curve = library.reference_curve(_FLOOR_CI)
    ci_char = contrast_index(char.x, char.y, float(char.y[0]))  # D_min* (6.5)
    ci_t = float(ci_at_times(ci_curve.x * 60.0, ci_curve.y, np.array([_FLOOR_TIME_S]))[0])
    floor = abs(ci_char - ci_t) / 2.0
    if not math.isfinite(floor):
        raise ValueError("sigma_floor of T-Max is not finite")
    return floor


# ---- builders ------------------------------------------------------------------------------


def _first_crossing(curve: PchipInterpolator, xs: np.ndarray, level: float) -> float:
    values = curve(xs) - level
    if values[0] >= 0.0:
        return float(xs[0])
    hits = np.flatnonzero((values[:-1] < 0.0) & (values[1:] >= 0.0))
    if hits.size == 0:
        raise ValueError(f"the curve never reaches {level:.3f}")
    i = int(hits[0])
    if values[i + 1] == 0.0:
        return float(xs[i + 1])
    return float(brentq(lambda t: float(curve(t)) - level, xs[i], xs[i + 1], xtol=_XTOL))


def _density(item: ConditionedItem, curve: ReferenceCurve, *, held_out: bool) -> Observation:
    d_min, _ = d_min_star(curve.x, curve.y, curve.x, curve.y)  # spec 6.3 step 2
    pchip = PchipInterpolator(curve.x, curve.y, extrapolate=False)
    start = _first_crossing(pchip, curve.x, d_min + DENSITY_START_RISE)
    x = np.linspace(start, float(curve.x[-1]), DENSITY_POINTS)
    y = pchip(x)
    weight = np.where(y < d_min + DENSITY_LOW_RISE, DENSITY_LOW_WEIGHT, 1.0)
    sigma = math.nan if held_out else density_sigma()
    return Observation(
        "density",
        item.item_id,
        item.condition,
        x,
        y,
        np.full(x.size, sigma),
        weight,
        DENSITY_N_EFFECTIVE,
        d_min_ref=d_min,
        in_objective=not held_out,
    )


def _ci(item: ConditionedItem, curve: ReferenceCurve, floor: float | None) -> Observation:
    x = np.asarray(item.condition.times_s, dtype=np.float64)
    y = ci_at_times(curve.x * 60.0, curve.y, x)
    sigma = math.nan if floor is None else ci_sigma(floor)
    return Observation(
        "ci",
        item.item_id,
        item.condition,
        x,
        y,
        np.full(x.size, sigma),
        np.ones(x.size),
        CI_N_EFFECTIVE,
        in_objective=floor is not None,
    )


def _speed(
    items: Sequence[ConditionedItem], curves: dict[str, ReferenceCurve]
) -> Observation | None:
    """log H_m* of one development series at Delta D = 0.80 (spec 6.3 steps 2-5)."""
    rows = []
    for item in sorted(items, key=lambda i: i.condition.times_s[0]):
        c = curves[item.item_id]
        log_h_m, delta_d = speed_point(c.x, c.y, float(c.y[0]))
        rows.append((item.condition.times_s[0], log_h_m, delta_d))
    for (t0, h0, d0), (t1, h1, d1) in itertools.pairwise(rows):
        if not all(math.isfinite(v) for v in (h0, d0, h1, d1)):
            continue
        if (d0 - _ISO_DELTA_D) * (d1 - _ISO_DELTA_D) <= 0.0 and d0 != d1:
            f = (_ISO_DELTA_D - d0) / (d1 - d0)
            t_star, h_star = t0 + f * (t1 - t0), h0 + f * (h1 - h0)
            ids = _JOIN.join(sorted(i.item_id for i in items))
            first = items[0].condition
            cond = dataclasses.replace(
                first, times_s=tuple(sorted(i.condition.times_s[0] for i in items)), source=ids
            )
            return Observation(
                "speed",
                ids,
                cond,
                [t_star],
                [h_star],
                [math.nan],
                [0.0],
                1.0,
                d_min_ref=None,
                in_objective=False,
            )
    return None


_CI_SLUG = {"d-76": "d76", "xtol": "xtol", "hc-110-pre2019": "hc110-b", "t-max-dev": "tmax-dev"}


def _slope_per_min(
    library: Library, material_id: str, bath: str, dilution: str, aim: float
) -> float:
    """dCI/dt (1/min) of the pair's small-tank 20 C CI(t) curve at CI_aim (spec 7.2)."""
    suffix = "-1-1" if dilution == "1+1" else ""
    rel = f"{material_id}/ci-small-tank-{_CI_SLUG[bath]}{suffix}.csv"
    curve = library.reference_curve(rel)
    pchip = PchipInterpolator(curve.x, curve.y)
    deriv = pchip.derivative()
    lo, hi = float(curve.y.min()), float(curve.y.max())
    if aim <= lo:
        return float(deriv(curve.x[0]))  # outside the span: end-point slope (7.2)
    if aim >= hi:
        return float(deriv(curve.x[-1]))
    t = _first_crossing(pchip, curve.x, aim)
    return float(deriv(t))


def _class_key(item: ConditionedItem) -> tuple[str, str, str, float]:
    c = item.condition
    return (c.developer, c.dilution, str(item.tank), float(item.ci_aim or 0.0))


def _table_ci(
    items: Sequence[ConditionedItem], library: Library, *, mode: Literal["train", "within", "held"]
) -> list[Observation]:
    counts: Counter[tuple[str, str, str, float]] = Counter(
        _class_key(i) for i in items if i.temperature_f != _REFERENCE_F
    )
    slopes: dict[tuple[str, str, float], float] = {}
    out = []
    for item in items:
        c = item.condition
        assert item.ci_aim is not None
        sigma = math.nan
        if mode != "held":
            key = (c.developer, c.dilution, item.ci_aim)
            if key not in slopes:
                slopes[key] = _slope_per_min(library, c.material_id, *key)
            off_20 = item.temperature_f != _REFERENCE_F
            n_class = counts[_class_key(item)] if (off_20 and mode == "train") else 1
            xtol_stock = c.developer == "xtol" and c.dilution == "stock" and off_20
            sigma = table_ci_sigma(slopes[key], n_class, xtol_stock_off_20c=xtol_stock)
            # One class counts about one independent number; the XTOL factor 2 is a quarter.
            n_eff = 1.0 / (n_class * (XTOL_STOCK_FACTOR**2 if xtol_stock else 1.0))
        else:
            n_eff = 1.0
        out.append(
            Observation(
                "table_ci",
                item.item_id,
                c,
                [c.times_s[0]],
                [item.ci_aim],
                [sigma],
                [1.0],
                n_eff,
                in_objective=mode == "train",
            )
        )
    return out


def _spectral(
    item: ConditionedItem, curve: ReferenceCurve, *, absolute: bool
) -> Observation | None:
    lo, hi = float(curve.x[0]), float(curve.x[-1])
    pchip = PchipInterpolator(curve.x, curve.y, extrapolate=False)
    if absolute:
        # 6.9 step 6 (within check): absolute log S on 450-600 nm, every 5 nm in the span.
        a, b = SPECTRAL_SHAPE_WINDOW_NM
        bands = np.arange(a, b + 0.5 * SPECTRAL_REPORT_STEP_NM, SPECTRAL_REPORT_STEP_NM)
        bands = bands[(bands >= lo) & (bands <= hi)]
        if bands.size == 0:
            return None
        return Observation(
            "spectral",
            item.item_id,
            item.condition,
            bands,
            pchip(bands),
            np.full(bands.size, math.nan),
            np.ones(bands.size),
            float(bands.size),
            in_objective=False,
        )
    bands = np.array([b for b in SPECTRAL_BANDS_NM if lo <= b <= hi])
    a, b = SPECTRAL_SHAPE_WINDOW_NM
    window = np.arange(max(a, lo), min(b, hi) + 1e-9, SPECTRAL_REPORT_STEP_NM)
    if bands.size == 0 or window.size == 0:
        return None
    shape = pchip(bands) - float(np.max(pchip(window)))
    return Observation(
        "spectral",
        item.item_id,
        item.condition,
        bands,
        shape,
        np.full(bands.size, SPECTRAL_SIGMA),
        np.ones(bands.size),
        float(bands.size),
    )


def _mtf_split(material_id: str, curve: ReferenceCurve) -> tuple[np.ndarray, np.ndarray]:
    """(tail, peak) frequencies of spec 6.7 step 8 inside the reference span."""
    freqs = MTF_FREQUENCIES_C_MM + (MTF_FREQUENCIES_TMAX_EXTRA if material_id == _TM else ())
    f = np.array([v for v in freqs if curve.x[0] <= v <= curve.x[-1]])
    above = np.flatnonzero(curve.y >= MTF_FULL_RESPONSE)
    if above.size == 0:
        cut = -math.inf
    else:
        i = int(above[-1])
        if i == curve.x.size - 1:
            cut = math.inf
        else:
            x0, x1, y0, y1 = curve.x[i], curve.x[i + 1], curve.y[i], curve.y[i + 1]
            cut = float(x0 + (MTF_FULL_RESPONSE - y0) * (x1 - x0) / (y1 - y0))
    return f[f > cut], f[f <= cut]


def _mtf(item: ConditionedItem, curve: ReferenceCurve, *, peak: bool) -> Observation | None:
    tail, peak_f = _mtf_split(item.condition.material_id, curve)
    f = peak_f if peak else tail
    if f.size == 0:
        return None
    pchip = PchipInterpolator(curve.x, curve.y, extrapolate=False)
    sigma = np.full(f.size, MTF_SIGMA)
    if item.condition.material_id == _TX:
        sigma[f >= MTF_TRIX_HIGH_C_MM] = MTF_SIGMA_TRIX_HIGH
    return Observation(
        "mtf",
        item.item_id,
        item.condition,
        f,
        pchip(f),
        sigma,
        np.ones(f.size),
        float(f.size),
        in_objective=not peak,
    )


def _fact(tables: dict, fact_id: str) -> dict[str, Any] | None:
    for fact in tables.get("facts") or []:
        if fact.get("id") == fact_id:
            return fact
    return None


def _granularity(item: ConditionedItem, library: Library, *, train: bool) -> Observation | None:
    fact = _fact(library.reference_tables(item.condition.material_id), "rms-granularity")
    if fact is None or not isinstance(fact.get("value"), int | float):
        return None
    return Observation(
        "granularity",
        item.item_id,
        item.condition,
        [GRANULARITY_NET_DENSITY],
        [float(fact["value"])],
        [GRANULARITY_SIGMA if train else math.nan],
        [1.0],
        1.0,
        in_objective=train,
    )


def _reciprocity(item: ConditionedItem, library: Library, *, train: bool) -> Observation | None:
    tables = library.reference_tables(item.condition.material_id)
    t = item.condition.exposure_time_s
    for row in (tables.get("reciprocity") or {}).get("rows") or []:
        if float(row["indicated_time_s"]) == t and row.get("aperture_stops") is not None:
            return Observation(
                "reciprocity",
                item.item_id,
                item.condition,
                [t],
                [float(row["aperture_stops"])],
                [RECIPROCITY_SIGMA],
                [1.0],
                1.0,
                in_objective=train,
            )
    return None


def _ph(library: Library) -> list[Observation]:
    """pH anchors of the bath files (spec 7.2); Carlton & Crabtree colorimetric sigma 0.10."""
    out = []
    for bath_id in library.ids("bath"):
        bath = library.bath(bath_id)
        for anchor in getattr(bath, "ph_anchors", None) or []:
            colorimetric = anchor.method == "colorimetric"
            sigma = (
                PH_SIGMA_COLORIMETRIC if colorimetric else float(anchor.ph.uncertainty)  # type: ignore[arg-type]
            )
            suffix = DELTA_CC_SUFFIX if colorimetric else ""
            out.append(
                Observation(
                    "ph",
                    f"{PH_PREFIX}{bath_id}#ph_anchors/{anchor.dilution}{suffix}",
                    None,
                    [float(anchor.temperature.value)],
                    [float(anchor.ph.value)],  # type: ignore[arg-type]
                    [sigma],
                    [1.0],
                    1.0,
                )
            )
    return out


def _film_observations(
    material_id: str,
    library: Library,
    items: Sequence[ConditionedItem],
    *,
    mode: Literal["train", "within", "held"],
) -> list[Observation]:
    out: list[Observation] = []
    curves: dict[str, ReferenceCurve] = {}
    for item in items:
        if item.kind in ("density", "ci", "mtf", "spectral"):
            curves[item.item_id] = library.reference_curve(item.item_id)
    floor = sigma_floor(material_id, library) if mode == "train" else None
    series: dict[tuple, list[ConditionedItem]] = {}
    for item in items:
        if item.kind == "density":
            out.append(_density(item, curves[item.item_id], held_out=mode == "held"))
            c = item.condition
            series.setdefault((c.developer, c.dilution, c.film_format, item.tank), []).append(item)
    if mode in ("train", "held"):
        # A series that does not bracket Delta D = 0.80 has no log H_m* (8.3: not applicable).
        for members in series.values():
            speed = _speed(members, curves)
            if speed is not None:
                out.append(speed)
    for item in items:
        if item.kind == "ci":
            out.append(_ci(item, curves[item.item_id], floor))
    out.extend(_table_ci([i for i in items if i.kind == "table_ci"], library, mode=mode))
    for item in items:
        obs: Observation | None = None
        if item.kind == "spectral":
            obs = _spectral(item, curves[item.item_id], absolute=False)
        elif item.kind == "mtf":
            obs = _mtf(item, curves[item.item_id], peak=False)
        elif item.kind == "granularity":
            obs = _granularity(item, library, train=mode == "train")
        elif item.kind == "reciprocity":
            obs = _reciprocity(item, library, train=mode == "train")
        if obs is not None:
            out.append(obs)
    return out


def training_observations(library: TrainingLibrary) -> list[Observation]:
    """Every spec 7.2 observation buildable from the training items (spec 7.6).

    Raises ``TypeError`` for anything but a :class:`TrainingLibrary` and
    :class:`LeakageError` if a non-training id would be returned.
    """
    if not isinstance(library, TrainingLibrary):
        raise TypeError(
            f"training_observations reads only through a TrainingLibrary (spec 7.7), "
            f"got {type(library).__name__}"
        )
    out: list[Observation] = []
    for material_id in _FILMS:
        items = reference_conditions(material_id, library, split="train")
        out.extend(_film_observations(material_id, library, items, mode="train"))
    out.extend(_ph(library))
    for obs in out:
        for source in observation_sources(obs):
            assert_trainable(source)
    return out


def within_check_observations(library: Library) -> list[Observation]:
    """Within-pair checks of spec 7.6 (``WITHIN_CHECKS``), never in the objective.

    Adds the training-file checks the registry does not list as items: the MTF peak band
    (6.7 step 8) and the T-Max absolute spectral level (6.9 step 6).
    """
    view = library if isinstance(library, TrainingLibrary) else TrainingLibrary(library)
    out: list[Observation] = []
    for material_id in _FILMS:
        items = reference_conditions(material_id, view, split="within")
        out.extend(_film_observations(material_id, view, items, mode="within"))
        for item in reference_conditions(material_id, view, split="train"):
            if item.kind == "mtf":
                obs = _mtf(item, view.reference_curve(item.item_id), peak=True)
            elif item.kind == "spectral" and material_id == _TM:
                obs = _spectral(item, view.reference_curve(item.item_id), absolute=True)
            else:
                continue
            if obs is not None:
                out.append(obs)
    for obs in out:
        for source in observation_sources(obs):
            if classify(source).split == "held_out":
                raise LeakageError(f"{source!r} is held out (7.6)")
    return out


def held_out_observations(library: Library, record: CalibrationRecord) -> list[Observation]:
    """Held-out evaluated items (spec 7.6); only for a record with ``frozen: true`` (7.7).

    Items without a Kodak process file get agitation and tank from the processes of the same
    tank type (:func:`~darkroom_atlas.validate.conditions.reference_conditions`). Includes the
    T-Max rows of ``developers/*.yaml`` and the log H_m* of the development series (8.3).
    Raises ``ValueError`` listing every evaluated held-out item that gives no observation and
    is not a null cell (8.2: a metric with reference data is never skipped).
    """
    if not bool(getattr(record, "frozen", False)):
        raise LeakageError("held-out observations need a frozen calibration record (7.7, 8.4)")
    if isinstance(library, TrainingLibrary):
        raise LeakageError("held-out observations cannot be read through a TrainingLibrary")
    out: list[Observation] = []
    nulls: set[str] = set()
    for material_id in _FILMS:
        items = reference_conditions(
            material_id, library, split="held_out", require_process=False, nulls=nulls
        )
        out.extend(_film_observations(material_id, library, items, mode="held"))
    covered: set[str] = set()
    for obs in out:
        for source in observation_sources(obs):
            c = classify(source)
            if c.split != "held_out" or not c.evaluated:
                raise LeakageError(f"{source!r} is not an evaluated held-out item")
            covered.add(source)
    evaluated = {
        i
        for i in enumerate_items(library.root / "reference")
        if (c := classify(i)).split == "held_out" and c.evaluated
    }
    missing = sorted(evaluated - covered - nulls)
    if missing:
        raise ValueError(
            "evaluated held-out items without an observation (spec 8.2):\n" + "\n".join(missing)
        )
    return out
