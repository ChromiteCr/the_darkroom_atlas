"""Agitation periods and the shared, event-aligned time-step table (spec D1 2.4.8, 3.4).

Every column, member, pre-pass and image block uses the same table. Breakpoints are 0, the
start and stop of every agitation period, every record (snapshot) time and the end. Within
each segment the step starts at ``h0``, grows by ``ratio`` and is capped at ``h_max``; the
last step is cut to the segment end, and folded into the previous one if shorter than
``h0 / 2`` (so a merged step can exceed ``h_max`` by less than ``h0 / 2``). The step restarts
at ``h0`` after every breakpoint, because agitation is a discontinuity of the data.
"""

import itertools
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np

from darkroom_atlas._frozen import FrozenDict, readonly_array
from darkroom_atlas.data.models import Agitation

__all__ = ["AgitationPeriod", "Schedule", "agitation_periods", "step_schedule"]

_ALWAYS = ("rotary", "continuous", "tray")
_PERIODIC = ("inversion", "lift")


@dataclass(frozen=True)
class AgitationPeriod:
    """One agitation period ``[start_s, stop_s)`` in seconds from immersion."""

    start_s: float
    stop_s: float


@dataclass(frozen=True, eq=False)
class Schedule:
    """A time-step table: ``times_s`` (S+1,), ``agitated`` (S,) and record-time indices.

    ``record_steps[t]`` is the index ``i`` with ``times_s[i] == t`` exactly, i.e. the number
    of steps completed at the snapshot (a record at the end maps to ``S``).
    """

    times_s: np.ndarray
    agitated: np.ndarray
    record_steps: Mapping[float, int]

    @property
    def n_steps(self) -> int:
        """Number of steps S."""
        return int(self.agitated.shape[0])

    @property
    def steps_s(self) -> np.ndarray:
        """Step lengths, shape (S,)."""
        return np.diff(self.times_s)


def agitation_periods(
    agitation: Agitation | None, duration_s: float
) -> tuple[AgitationPeriod, ...]:
    """Agitation periods of ``agitation`` within ``[0, duration_s]`` (spec D1 2.4.8 table).

    ``inversion`` and ``lift``: ``[0, initial_s]`` plus ``[k interval_s, k interval_s +
    duration_s]`` for every integer ``k >= 0`` with ``k interval_s < duration_s``.
    ``rotary``, ``continuous``, ``tray``: the whole development. ``stand``: ``[0, initial_s]``
    only. ``None`` means no agitation at all. Periods are clipped to the development, those
    that overlap or touch are merged and empty ones dropped. ``ValueError`` on a
    non-positive duration or a periodic method without ``interval_s`` or ``duration_s``.
    """
    _positive_finite("duration_s", duration_s)
    duration_s = float(duration_s)
    if agitation is None:
        return ()
    raw: list[tuple[float, float]] = []
    method = agitation.method
    if method in _ALWAYS:
        raw.append((0.0, duration_s))
    else:
        raw.append((0.0, float(agitation.initial_s)))
        if method in _PERIODIC:
            if agitation.interval_s is None or agitation.duration_s is None:
                raise ValueError(f"agitation method {method!r} needs interval_s and duration_s")
            interval = float(agitation.interval_s)
            length = float(agitation.duration_s)
            k = 0
            while k * interval < duration_s:
                raw.append((k * interval, k * interval + length))
                k += 1
        elif method != "stand":
            raise ValueError(f"unknown agitation method {method!r}")
    return _union(raw, duration_s)


def step_schedule(
    periods: Sequence[AgitationPeriod],
    duration_s: float,
    *,
    h0_s: float,
    ratio: float,
    h_max_s: float,
    record_s: Sequence[float] = (),
) -> Schedule:
    """Build the event-aligned geometric step table (spec D1 3.4).

    ``periods`` must be sorted, non-overlapping, non-empty and inside ``[0, duration_s]``
    (as :func:`agitation_periods` returns them). Every period boundary and record time is a
    grid time, hit exactly. ``ValueError`` on ``h0_s <= 0``, ``ratio < 1``, ``h_max_s < h0_s``,
    a non-positive duration, a malformed period or a record time outside ``[0, duration_s]``.
    """
    _positive_finite("duration_s", duration_s)
    _positive_finite("h0_s", h0_s)
    _positive_finite("h_max_s", h_max_s)
    _positive_finite("ratio", ratio)
    if ratio < 1.0:
        raise ValueError(f"ratio must be >= 1, got {ratio!r}")
    if h_max_s < h0_s:
        raise ValueError(f"h_max_s ({h_max_s!r}) must be >= h0_s ({h0_s!r})")
    duration_s = float(duration_s)

    previous_stop = 0.0
    spans: list[tuple[float, float]] = []
    for period in periods:
        start, stop = float(period.start_s), float(period.stop_s)
        if not (math.isfinite(start) and math.isfinite(stop)) or not start < stop:
            raise ValueError(f"agitation period must have start < stop, got {period!r}")
        if start < previous_stop or stop > duration_s:
            raise ValueError(
                f"agitation periods must be sorted, disjoint and within [0, {duration_s}], "
                f"got {period!r}"
            )
        spans.append((start, stop))
        previous_stop = stop

    records = [float(t) for t in record_s]
    for t in records:
        if not math.isfinite(t) or t < 0.0 or t > duration_s:
            raise ValueError(f"record time {t!r} is outside [0, {duration_s}]")

    breakpoints = sorted({0.0, duration_s, *records, *(x for span in spans for x in span)})
    times = [0.0]
    for start, end in itertools.pairwise(breakpoints):
        times.extend(_segment(start, end, float(h0_s), float(ratio), float(h_max_s)))
    times_s = np.array(times, dtype=np.float64)

    midpoints = 0.5 * (times_s[:-1] + times_s[1:])
    agitated = np.zeros(midpoints.shape, dtype=bool)
    for start, stop in spans:
        agitated |= (midpoints >= start) & (midpoints < stop)

    index = {t: i for i, t in enumerate(times)}
    record_steps = FrozenDict({t: index[t] for t in records})
    return Schedule(
        times_s=readonly_array(times_s),
        agitated=readonly_array(agitated),
        record_steps=record_steps,
    )


def _segment(start: float, end: float, h0: float, ratio: float, h_max: float) -> list[float]:
    """Step end times in ``(start, end]``; the last one is exactly ``end``."""
    length = end - start
    steps: list[float] = []
    elapsed = 0.0
    h = h0
    while elapsed < length:
        step = min(h, h_max, length - elapsed)
        steps.append(step)
        elapsed += step
        h *= ratio
    if len(steps) > 1 and steps[-1] < 0.5 * h0:
        tail = steps.pop()
        steps[-1] += tail
    out = list(start + np.cumsum(steps))
    out[-1] = end
    return [float(t) for t in out]


def _union(raw: list[tuple[float, float]], duration_s: float) -> tuple[AgitationPeriod, ...]:
    clipped = sorted(
        (max(0.0, a), min(b, duration_s)) for a, b in raw if min(b, duration_s) > max(0.0, a)
    )
    merged: list[list[float]] = []
    for a, b in clipped:
        if merged and a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    return tuple(AgitationPeriod(a, b) for a, b in merged)


def _positive_finite(name: str, value: float) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, float, np.floating, np.integer)):
        raise TypeError(f"{name} must be a real number, got {value!r}")
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be positive and finite, got {value!r}")
