import numpy as np
import pytest

from darkroom_atlas.data.models import Agitation
from darkroom_atlas.develop.schedule import (
    AgitationPeriod,
    agitation_periods,
    step_schedule,
)

IMAGE = {"h0_s": 0.1, "ratio": 1.6, "h_max_s": 5.0}
CALIBRATION = {"h0_s": 0.05, "ratio": 1.4, "h_max_s": 3.0}
COARSE = {"h0_s": 0.2, "ratio": 2.0, "h_max_s": 10.0}


def _periods(pairs):
    return tuple(AgitationPeriod(a, b) for a, b in pairs)


def _boundaries(periods):
    return {x for p in periods for x in (p.start_s, p.stop_s)}


# ---- agitation_periods -------------------------------------------------------------------


def test_inversion_kodak_small_tank():
    agitation = Agitation(method="inversion", initial_s=30.0, interval_s=30.0, duration_s=5.0)
    periods = agitation_periods(agitation, 100.0)
    assert periods == _periods([(0.0, 35.0), (60.0, 65.0), (90.0, 95.0)])


def test_inversion_initial_shorter_than_duration():
    agitation = Agitation(method="inversion", initial_s=2.0, interval_s=30.0, duration_s=5.0)
    assert agitation_periods(agitation, 70.0) == _periods([(0.0, 5.0), (30.0, 35.0), (60.0, 65.0)])


def test_lift_large_tank():
    agitation = Agitation(method="lift", initial_s=22.5, interval_s=60.0, duration_s=7.5)
    periods = agitation_periods(agitation, 200.0)
    assert periods == _periods([(0.0, 22.5), (60.0, 67.5), (120.0, 127.5), (180.0, 187.5)])


def test_periods_are_clipped_to_duration():
    agitation = Agitation(method="inversion", initial_s=5.0, interval_s=30.0, duration_s=5.0)
    assert agitation_periods(agitation, 62.0)[-1] == AgitationPeriod(60.0, 62.0)
    # A period that would start exactly at the end is dropped.
    assert agitation_periods(agitation, 60.0)[-1] == AgitationPeriod(30.0, 35.0)


def test_duration_at_least_interval_is_continuous():
    agitation = Agitation(method="inversion", initial_s=0.0, interval_s=10.0, duration_s=10.0)
    assert agitation_periods(agitation, 45.0) == _periods([(0.0, 45.0)])


@pytest.mark.parametrize("method", ["rotary", "continuous", "tray"])
def test_continuous_methods_cover_everything(method):
    assert agitation_periods(Agitation(method=method), 300.0) == _periods([(0.0, 300.0)])


def test_stand_only_initial():
    agitation = Agitation(method="stand", initial_s=30.0, interval_s=60.0, duration_s=10.0)
    assert agitation_periods(agitation, 600.0) == _periods([(0.0, 30.0)])
    assert agitation_periods(Agitation(method="stand"), 600.0) == ()


def test_no_agitation_is_empty():
    assert agitation_periods(None, 100.0) == ()


def test_periodic_without_interval_raises():
    with pytest.raises(ValueError, match="interval_s"):
        agitation_periods(Agitation(method="inversion", initial_s=5.0), 100.0)
    with pytest.raises(ValueError, match="interval_s"):
        agitation_periods(Agitation(method="lift", initial_s=5.0, interval_s=60.0), 100.0)


@pytest.mark.parametrize("duration", [0.0, -1.0, float("nan"), float("inf")])
def test_agitation_bad_duration_raises(duration):
    with pytest.raises(ValueError):
        agitation_periods(Agitation(method="tray"), duration)


# ---- step_schedule -----------------------------------------------------------------------


def _instant_resets(times, duration):
    """Touching periods starting at 0 and at every time: each one is a reset point."""
    starts = sorted({0.0, *times})
    return _periods(list(zip(starts, [*starts[1:], duration], strict=True)))


@pytest.mark.parametrize(
    ("table", "every", "expected"),
    [
        (IMAGE, 30.0, 178),
        (CALIBRATION, 30.0, 274),
        (COARSE, 30.0, 110),
        (IMAGE, 60.0, 136),
    ],
)
def test_prototype_step_counts(table, every, expected):
    # Spec D1 3.4 / 12.6: 400 s, instantaneous liquid reset every 30 s (or 60 s). The 60 s
    # row is reproduced with resets {0, 30, 60, 120, ..., 360}: a 30 s initial period.
    times = list(np.arange(every, 400.0, every))
    if every == 60.0:
        times = [30.0, *times]
    periods = _instant_resets(times, 400.0)
    assert {p.start_s for p in periods} == {0.0, *times}
    schedule = step_schedule(periods, 400.0, **table)
    assert schedule.n_steps == expected


# Spec D1 3.4 table. The small-tank initial agitation is not given there; initial_s = 15 s
# reproduces the table (5 s and 10 s give one step more, 30 s about ten fewer).
SMALL_TANK = Agitation(method="inversion", initial_s=15.0, interval_s=30.0, duration_s=5.0)
LARGE_TANK = Agitation(method="lift", initial_s=22.5, interval_s=60.0, duration_s=7.5)


@pytest.mark.parametrize(
    ("agitation", "minutes", "image", "calibration"),
    [
        (SMALL_TANK, 7.5, 209, 314),
        (SMALL_TANK, 12.0, 335, 503),
        (LARGE_TANK, 8.0, 160, 248),
        (LARGE_TANK, 12.0, 240, 372),
    ],
)
def test_condition_step_counts(agitation, minutes, image, calibration):
    duration = 60.0 * minutes
    periods = agitation_periods(agitation, duration)
    assert step_schedule(periods, duration, **IMAGE).n_steps == image
    assert step_schedule(periods, duration, **CALIBRATION).n_steps == calibration


@pytest.mark.parametrize("agitation", [SMALL_TANK, LARGE_TANK])
@pytest.mark.parametrize("table", [IMAGE, CALIBRATION])
def test_snapshots_add_no_steps(agitation, table):
    # Spec D1 3.4: snapshots are alignment points, so 15 s and 60 s snapshots cost the same.
    duration = 1080.0
    periods = agitation_periods(agitation, duration)
    every_15 = step_schedule(periods, duration, record_s=np.arange(15.0, 1081.0, 15.0), **table)
    every_60 = step_schedule(periods, duration, record_s=np.arange(60.0, 1081.0, 60.0), **table)
    bare = step_schedule(periods, duration, **table)
    assert every_15.n_steps == every_60.n_steps == bare.n_steps
    for r, i in every_15.record_steps.items():
        assert every_15.times_s[i] == r


@pytest.mark.parametrize("table", [IMAGE, CALIBRATION, COARSE])
def test_kodak_inversion_schedule_properties(table):
    agitation = Agitation(method="inversion", initial_s=30.0, interval_s=30.0, duration_s=5.0)
    periods = agitation_periods(agitation, 487.0)
    records = [45.0, 100.25, 487.0, 0.0]
    schedule = step_schedule(periods, 487.0, record_s=records, **table)
    t = schedule.times_s
    assert t[0] == 0.0
    assert t[-1] == 487.0
    assert np.all(np.diff(t) > 0)
    assert schedule.agitated.shape == (t.size - 1,)
    # Every boundary is a grid time, so no step straddles one.
    for b in _boundaries(periods):
        assert b in set(t.tolist())
    for p in periods:
        inside = (t[:-1] >= p.start_s) & (t[1:] <= p.stop_s)
        straddle = (t[:-1] < p.start_s) & (t[1:] > p.start_s) | (t[:-1] < p.stop_s) & (
            t[1:] > p.stop_s
        )
        assert not straddle.any()
        assert schedule.agitated[inside].all()
    # Steps never exceed h_max by h0/2 or more.
    assert np.max(np.diff(t)) < table["h_max_s"] + 0.5 * table["h0_s"]
    assert schedule.agitated.sum() == sum(
        np.count_nonzero((t[:-1] >= p.start_s) & (t[1:] <= p.stop_s)) for p in periods
    )
    # Record times map to exact step indices.
    for r in records:
        assert t[schedule.record_steps[r]] == r
    assert schedule.record_steps[487.0] == schedule.n_steps
    assert schedule.record_steps[0.0] == 0


def test_steps_restart_at_h0_only_at_period_starts():
    schedule = step_schedule(_periods([(0.0, 10.0), (20.0, 25.0)]), 40.0, **IMAGE)
    t = schedule.times_s
    steps = np.diff(t)
    assert steps[0] == pytest.approx(0.1)
    # The stop at 10 s is only aligned to: the step after it is not h0.
    i = int(np.flatnonzero(t == 10.0)[0])
    assert steps[i] > 1.0
    # The start at 20 s restarts the sequence.
    k = int(np.flatnonzero(t == 20.0)[0])
    assert steps[k] == pytest.approx(0.1)
    assert steps[k + 1] == pytest.approx(0.16)


def test_cut_step_at_alignment_does_not_advance_sequence():
    # Steps 0.1, 0.16, 0.256, 0.4096 end at 0.9256; the next step 0.65536 would pass 1.2, so
    # it is cut to 1.2 (0.2744 >= h0/2) and the step after it is 0.65536 again.
    schedule = step_schedule((), 3.0, record_s=[1.2], **IMAGE)
    steps = schedule.steps_s
    np.testing.assert_allclose(steps[:4], [0.1, 0.16, 0.256, 0.4096])
    assert schedule.times_s[5] == 1.2
    assert steps[4] == pytest.approx(1.2 - 0.9256)
    assert steps[5] == pytest.approx(0.65536)
    assert schedule.record_steps[1.2] == 5


def test_short_cut_at_alignment_stretches_previous_step():
    # 0.9256 + 0.65536 = 1.58096; an alignment at 1.6 leaves 0.01904 (< h0/2), so the
    # previous step is stretched to 1.6 and the sequence continues with 1.048576.
    schedule = step_schedule((), 4.0, record_s=[1.6], **IMAGE)
    steps = schedule.steps_s
    assert schedule.times_s[5] == 1.6
    assert steps[4] == pytest.approx(1.6 - 0.9256)
    assert steps[5] == pytest.approx(1.048576)


def test_short_tail_is_merged_into_previous_step():
    # 0.1 + 0.16 + 0.256 = 0.516; a 0.54 s segment would leave a 0.024 s tail (< h0/2).
    schedule = step_schedule((), 0.54, **IMAGE)
    np.testing.assert_allclose(schedule.steps_s, [0.1, 0.16, 0.28])
    # A tail of h0/2 or more is kept.
    assert step_schedule((), 0.6, **IMAGE).n_steps == 4


def test_segment_shorter_than_h0_is_one_step():
    schedule = step_schedule(_periods([(0.0, 0.01)]), 1.0, **IMAGE)
    assert schedule.times_s[1] == 0.01
    assert schedule.agitated[0]
    assert not schedule.agitated[1:].any()


def test_identical_input_identical_schedule():
    agitation = Agitation(method="lift", initial_s=22.5, interval_s=60.0, duration_s=7.5)
    a = step_schedule(agitation_periods(agitation, 540.0), 540.0, record_s=[120.0], **IMAGE)
    b = step_schedule(agitation_periods(agitation, 540.0), 540.0, record_s=[120.0], **IMAGE)
    assert np.array_equal(a.times_s, b.times_s)
    assert np.array_equal(a.agitated, b.agitated)
    assert dict(a.record_steps) == dict(b.record_steps)


def test_schedule_is_read_only():
    schedule = step_schedule(_periods([(0.0, 5.0)]), 20.0, record_s=[10.0], **IMAGE)
    with pytest.raises(ValueError):
        schedule.times_s[0] = 1.0
    with pytest.raises(ValueError):
        schedule.agitated[0] = False
    with pytest.raises(TypeError):
        schedule.record_steps[3.0] = 1  # type: ignore[index]


@pytest.mark.parametrize(
    ("periods", "duration", "kwargs"),
    [
        ((), 0.0, {}),
        ((), -5.0, {}),
        ((), float("inf"), {}),
        ((), 10.0, {"h0_s": 0.0}),
        ((), 10.0, {"ratio": 0.9}),
        ((), 10.0, {"h_max_s": 0.05}),
        ((), 10.0, {"record_s": [11.0]}),
        ((), 10.0, {"record_s": [-1.0]}),
        ((), 10.0, {"record_s": [float("nan")]}),
        (_periods([(5.0, 5.0)]), 10.0, {}),
        (_periods([(6.0, 5.0)]), 10.0, {}),
        (_periods([(0.0, 11.0)]), 10.0, {}),
        (_periods([(-1.0, 2.0)]), 10.0, {}),
        (_periods([(4.0, 6.0), (1.0, 2.0)]), 10.0, {}),
        (_periods([(1.0, 6.0), (5.0, 7.0)]), 10.0, {}),
    ],
)
def test_bad_input_raises(periods, duration, kwargs):
    with pytest.raises(ValueError):
        step_schedule(periods, duration, **{**IMAGE, **kwargs})
