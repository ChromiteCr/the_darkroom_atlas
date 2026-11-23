# Kodak reference processes

One process file per distinct development condition of the Kodak reference curves and time
tables of `kodak/tri-x-400` and `kodak/t-max-100` (spec D1 5.6, 5.8, 6.1, 7.6). The files are
generated from the reference data (`data/reference/kodak/<film>/`, spec D1 5.9); every file lists the reference items that
use it under `refs`, in the item-id format of the split registry (spec D1 7.6,
`validate/splits.py`).

## Naming rule

```
kodak/<film>/<bath>-<dilution>-<temperature>c-<tank>-<time>min
```

| Part | Rule | Examples |
|------|------|----------|
| `<film>` | material id without `kodak/` | `tri-x-400`, `t-max-100` |
| `<bath>` | developer id in `data/developers/` | `d-76`, `xtol`, `hc-110-pre2019`, `t-max-dev` |
| `<dilution>` | explicit ratio of the step (`stock`, `a+b`) | `stock`, `1+1`, `1+4`, `1+31` |
| `<temperature>` | development temperature in °C to one decimal, `.0` dropped | `20`, `18.3`, `21.1`, `22.2`, `23.9` |
| `<tank>` | `small-tank`, `large-tank`, `rotary-tube`, `tray` | |
| `<time>` | development time in minutes, trailing zeros dropped | `8`, `6.75`, `7.8` |

So the spec's example `kodak/tri-x-400/d-76-stock-20c-large-tank-8min` is the 8-minute D-76
characteristic curve of Tri-X. The film format is not part of a process: Tri-X 35 mm and 120
curves developed alike share one file (the format is a run option, `RunContext.film_format`).
Rolls and sheets of the same tank, temperature and time share one file.

## Fixed choices

- Temperature: Kodak's °F converted exactly, `(F - 32) * 5/9` (spec D1 6.1), tagged
  `datasheet` with the table or caption it comes from. The label in the id is rounded; the
  value is not.
- Dilution: HC-110 dilution B is `1+31` of the concentrate with the pre-2019 formulation
  `hc-110-pre2019` (J-24 Dec 2017 p. 2; spec D1 5.5, 5.9); T-MAX Developer is `1+4`
  (J-86 Dec 2017 p. 1). D-76 and XTOL rows are `stock` or `1+1` as printed.
- Agitation (spec D1 2.4.8) and tank volume (spec D1 6.1):

  | Tank | `agitation` | `volume_ml` |
  |------|-------------|-------------|
  | small tank | `inversion`, `initial_s` 5, `interval_s` 30, `duration_s` 5 | 473 |
  | large tank | `lift`, `initial_s` 22.5, `interval_s` 60, `duration_s` 7.5 | 3785 |
  | rotary tube | `rotary` (whole time) | 250 |
  | tray | `tray` (whole time) | 1000 |

  J-109 Table 1 (Tri-X rows of `xtol_contrast_index_times`) states no agitation; the F-4017
  small-tank procedure is used.
- Tank film area 35 cm² and `image_area_fraction` 1.0: a sensitometric strip (spec D1 5.8);
  the virtual sensitometer sets its own strip area, 1 cm² per step (spec D1 6.1).
- After development: `stop-acetic` 30 s, `fixer-thiosulfate` 5 min, `water` 10 min
  (`assumed`, the spec D1 5.8 example; not part of the Kodak conditions).

## Which conditions have a file

- Characteristic curves (`char-*.csv`): every curve except T-MAX RS, which D1 does not use
  (spec D1 7.6). This includes the four held-out T-Max 100 × T-MAX Developer curves, whose
  condition (developer, time, temperature, tank) is in the file name and conditions header.
- Time-table cells (`tables.yaml#development_times`, Tri-X `#xtol_contrast_index_times`,
  T-Max `#development_time_adjustment_factors/d-76`): cells the split registry classifies as
  `training` or `within_checks`, plus the 68 °F cells that are `unused` only because they
  duplicate a CI(t) curve at 20 °C (they give each training CI(t) curve its Kodak 20 °C time,
  and include the granularity conditions of spec D1 6.6: T-Max D-76 small tank 6½ min, Tri-X
  HC-110 B small tank 3¾ min). The cells are classified before their minutes are read; no
  held-out table row (spec D1 7.6, 7.7) is read or written here.
- Adjustment factors, `[d-76]` block, Normal Contrast column: `t = f × 6½ min`, the small-tank
  D-76 stock time at 68 °F (spec D1 7.6 within checks).
- Not here: CI(t) curves as such (their condition is a time series, which a process with one
  `time` cannot express; the sensitometer reads a `Condition` with `times_s`, spec D1 6.1);
  held-out table rows and held-out CI(t) curves; push tables; rack-and-tank and replenishment
  rows; T-MAX RS.
