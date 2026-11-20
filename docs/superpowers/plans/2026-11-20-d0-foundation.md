# D0 底座与约定 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 搭好 The Darkroom Atlas 的底座：单位与光谱约定、带来源标签的数据模型、`expose / process / readout` 三个操作及其注册机制、分块与确定性随机数、运行记录和命令行，并让一张图经"理想材料"走通全流程。

**Architecture:** 物理实现通过注册表接入：感光材料文件里的 `model` 字段选择曝光模型和各处理步骤的处理器，出口设备按名字注册。D0 只注册一个 `ideal` 模型（线性记录 CIE XYZ），用它验证整条光谱链路的单位、形状、分块和可复现性。真实物理在 D1 起以新的注册项加入，不改底座。

**Tech Stack:** Python 3.13、uv、NumPy、SciPy、pydantic v2、PyYAML、typer、imageio + tifffile + pillow、colour-science（只用其 CIE 标准数据和 Mallett 2019 基函数）、array-api-compat；开发依赖 pytest、ruff、array-api-strict。

**Spec:** `docs/superpowers/specs/2026-11-20-darkroom-atlas-design.md`

## Global Constraints

- Python `>=3.13`，环境与依赖一律用 `uv`（`uv add`、`uv run`），不直接用 pip。
- 数组轴顺序固定为 `(H, W, ..., band)`，波段永远是最后一维。内部一律 `float64`。
- 波长网格默认 380–780 nm、步长 5 nm（81 个波段）。光谱量按"波段内已积分"存放，对波长积分就是对最后一维求和。
- 胶片面长度单位 µm；曝光量单位为每 µm² 每波段的光子数；密度为以 10 为底的光学密度。
- 单位字符串只用 ASCII：`um` `nm` `mm` `s` `min` `C` `K` `mL` `L` `g/L` `g/mol` `g/m2` `1`。
- 数据文件每个物理参数都是带来源标签的 `Quantity`；来源为 `datasheet` 或 `literature` 时必须有 `ref`。D0 不写入任何没有出处的真实胶片或药液数据。
- 所有随机数必须经 `rng.py`，由（种子、流名、像素绝对坐标）唯一决定；禁止 `np.random.seed`、无种子的 `default_rng()`。
- 状态对象不可变：操作返回新状态，数组设为只读。
- 文档字符串和标识符用英文；`README.md` 与 `docs/` 用中文。
- 提交信息不得出现任何 AI 或代理的名字，不加任何 `Co-authored-by` 尾注。只做本地提交，不推送，不建远程。
- 不添加 LICENSE 文件。
- 与本计划不一致的地方（如某个库里的名字和这里写的不同）：以能通过本计划测试的最小改动处理，并在任务汇报里写明偏差，不要悄悄改接口。

## Review Focus

设计文档没有明说、但使用者最容易撞上的五类输入。每条的测试已加到对应任务里。

1. **非标准输入图像**：灰度图、带透明通道的图、16 位 PNG、含 NaN/inf 的浮点 TIFF。应分别复制为三通道、丢弃透明通道、按 65535 归一、报错并指出问题。（Task 5）
2. **全黑画面**：平均测光时亮度为 0。应按 18% 灰处理而不是除以零。（Task 5）
3. **分块边界**：画面小于块、块不能整除画面、重叠大于块。结果必须与不分块完全一致。（Task 4、Task 6）
4. **数据文件写错**：多余字段、`datasheet` 来源缺 `ref`、`id` 与路径不符、引用了不存在的药液或物种。应报出文件路径和具体问题，而不是一段 pydantic 堆栈。（Task 3）
5. **命令行找不到东西**：未知的材料或流程名、输出目录不存在。应列出可用名字并以退出码 2 结束；输出目录自动创建。（Task 7）

## File Structure

```
pyproject.toml  .gitignore  README.md
src/darkroom_atlas/
  __init__.py        ATLAS_VERSION
  backend.py         Array API 命名空间入口
  units.py           物理常数、单位换算、时间与稀释比例解析
  spectral.py        SpectralGrid、CIE 数据、标准光源
  rng.py             按像素绝对坐标确定的随机数
  tiling.py          分块迭代与拼接
  light.py           Light
  state.py           FrameGeometry、StepRecord、EmulsionState
  registry.py        曝光模型 / 步骤处理器 / 出口设备的注册表与协议
  pipeline.py        RunContext、expose / process / readout
  manifest.py        RunManifest
  runner.py          run_image：读图到写图与运行记录的整条流程
  cli.py             atlas run / atlas check-data
  data/
    provenance.py    Source、Quantity 等
    models.py        Species、Bath、Material、Process 及其子模型
    library.py       Library：按 id 加载、解析曲线、整库校验
  scene/
    io.py            读写图像、sRGB 编解码
    upsample.py      RGB → 反射率光谱（Mallett 2019）
    camera.py        虚拟相机测光，产出 Light
  models/
    __init__.py      导入即注册
    ideal.py         ideal 曝光模型、identity 步骤、ideal_colorimetric 出口
data/
  chemistry/species.yaml   films/ideal.yaml   processes/identity.yaml
tests/               与源文件一一对应的 test_*.py；assets/；fixtures/library/
docs/physics.md
```

任务依赖：1 →（2 → 3）与 4 并行 → 5 → 6 → 7 与 8 并行。

---

### Task 1: 项目骨架

**Files:**
- Create: `pyproject.toml`, `.gitignore`, `src/darkroom_atlas/__init__.py`, `src/darkroom_atlas/backend.py`
- Create: `tests/test_smoke.py`, `tests/assets/make_sample.py`, `tests/assets/sample.png`

**Interfaces:**
- Produces: `darkroom_atlas.ATLAS_VERSION: str = "D0"`；`darkroom_atlas.backend.xp_of(*arrays)`，返回 `array_api_compat.array_namespace(*arrays)`；`tests/assets/sample.png`（96×64，8 位 RGB）。

- [ ] **Step 1: 建立 uv 项目。** `uv init --lib --name darkroom-atlas --python 3.13`（如生成了示例文件则删去），src 布局，构建后端 hatchling，包版本 `0.0.0`。运行时依赖：`numpy scipy pydantic pyyaml typer imageio tifffile pillow colour-science array-api-compat`。开发依赖组 `dev`：`pytest ruff array-api-strict`。`[project.scripts] atlas = "darkroom_atlas.cli:app"`。`[tool.pytest.ini_options] testpaths = ["tests"]`。`[tool.ruff] line-length = 100`。提交 `uv.lock`。
- [ ] **Step 2: `.gitignore`**：`.venv/ __pycache__/ .pytest_cache/ .ruff_cache/ *.egg-info/ dist/ out/ .DS_Store`。
- [ ] **Step 3: 写失败的测试** `tests/test_smoke.py`：

```python
def test_version():
    import darkroom_atlas
    assert darkroom_atlas.ATLAS_VERSION == "D0"

def test_xp_of_numpy_and_strict():
    import numpy as np, array_api_strict as xps
    from darkroom_atlas.backend import xp_of
    assert xp_of(np.zeros(2)).__name__.endswith("numpy")
    assert xp_of(xps.zeros(2)) is xps or "strict" in xp_of(xps.zeros(2)).__name__

def test_sample_asset():
    import imageio.v3 as iio
    from pathlib import Path
    img = iio.imread(Path(__file__).parent / "assets" / "sample.png")
    assert img.shape == (64, 96, 3) and img.dtype.name == "uint8"
```

- [ ] **Step 4: 运行确认失败**：`uv run pytest tests/test_smoke.py`，预期 ImportError / 文件不存在。
- [ ] **Step 5: 实现** `__init__.py`、`backend.py`；写 `tests/assets/make_sample.py`，确定性地生成 `sample.png`：左半 48 列为水平灰阶渐变（0→255），右半分成 2 行 × 4 列共 8 个色块，依次为红、绿、蓝、白、青、品红、黄、黑。运行它生成 PNG 并提交。
- [ ] **Step 6: 运行确认通过**：`uv run pytest tests/test_smoke.py`，`uv run ruff check .`。
- [ ] **Step 7: 提交**：`git add pyproject.toml uv.lock .gitignore src tests && git commit -m "chore: project skeleton with uv, pytest and sample asset"`

---

### Task 2: 单位与光谱网格

**Files:**
- Create: `src/darkroom_atlas/units.py`, `src/darkroom_atlas/spectral.py`
- Test: `tests/test_units.py`, `tests/test_spectral.py`

**Interfaces:**
- Consumes: `backend.xp_of`
- Produces（`units.py`）:
  - 常数 `H_PLANCK = 6.62607015e-34`、`C_LIGHT = 299792458.0`、`N_AVOGADRO = 6.02214076e23`、`R_GAS = 8.314462618`、`K_M = 683.0`
  - `photon_energy_j(wavelength_nm)`、`transmittance_to_density(t, *, floor=1e-12)`、`density_to_transmittance(d)`（这三个按 Array API 写）
  - `convert(value, from_unit: str, to_unit: str)`：支持长度 `nm um mm cm m`、时间 `s min h`、温度 `C K`、体积 `mL L`、质量浓度 `g/L mg/L`；两个单位相同则原样返回；量纲不符或单位未知抛 `ValueError`
  - `parse_duration_s(text: str | float | int) -> float`
  - `parse_dilution(text: str) -> float`（原液占最终体积的比例）
  - `pixel_pitch_um(long_edge_mm: float, shape_hw: tuple[int, int]) -> float`
- Produces（`spectral.py`）:
  - `@dataclass(frozen=True) class SpectralGrid(start_nm=380.0, stop_nm=780.0, step_nm=5.0)`，属性 `wavelengths -> np.ndarray`、`n_bands -> int`；方法 `resample(wavelengths_nm, values, *, outside: Literal["zero","hold","error"]="error") -> np.ndarray`（线性插值）、`integrate(values, axis=-1)`（对该轴求和，因为光谱量已按波段积分）
  - `DEFAULT_GRID`
  - `cmf_xyz(grid) -> (B, 3)`：CIE 1931 2° 色匹配函数，取自 colour-science
  - `illuminant(name: str, grid) -> (B,)`：`"D65" "D50" "A" "E"`，返回按波段积分的相对能量，归一化到 `sum(spd * ybar) == 1`
  - `blackbody(temperature_k: float, grid) -> (B,)`：普朗克公式，同样归一化
  - `XYZ_TO_LINEAR_SRGB: (3, 3)` 与 `LINEAR_SRGB_TO_XYZ`：IEC 61966-2-1 矩阵

- [ ] **Step 1: 写失败的测试。**

```python
# tests/test_units.py
def test_photon_energy_555nm():      assert photon_energy_j(555.0) == pytest.approx(3.5792e-19, rel=1e-4)
def test_density_roundtrip():        # d = [0, 0.3, 1, 3]; density(transmittance(d)) == d (abs 1e-12)
def test_density_floor():            assert transmittance_to_density(0.0) == pytest.approx(12.0)
def test_density_array_api_strict(): # 输入 array_api_strict 数组，返回同命名空间数组
def test_convert():                  assert convert(1.0, "mm", "um") == 1000.0; convert(20.0, "C", "K") == 293.15; convert(2.0, "min", "s") == 120.0; convert(5, "g/L", "g/L") == 5
def test_convert_rejects():          # convert(1, "mm", "s") 与 convert(1, "furlong", "m") 抛 ValueError
def test_parse_duration():           # "9:45"→585.0  "1:02:03"→3723.0  "30s"→30.0  "8min"→480.0  90→90.0  "1.5min"→90.0
def test_parse_duration_rejects():   # "", "abc", "-5", "9:75" 抛 ValueError
def test_parse_dilution():           # "stock"→1.0  "1+1"→0.5  "1+50"→1/51  "1+0"→1.0
def test_parse_dilution_rejects():   # "1:1"（有歧义）, "0+1", "x" 抛 ValueError
def test_pixel_pitch():              assert pixel_pitch_um(36.0, (4000, 6000)) == 6.0; pixel_pitch_um(36.0, (6000, 4000)) == 6.0
```

```python
# tests/test_spectral.py
def test_default_grid():             # n_bands == 81, wavelengths[0] == 380, wavelengths[-1] == 780
def test_grid_rejects_non_integer(): # SpectralGrid(380, 780, 7) 抛 ValueError；step<=0、stop<=start 同样
def test_resample_linear_exact():    # 对 y = 2x+1 在 10 nm 采样后重采样到默认网格，误差 < 1e-9
def test_resample_outside_modes():   # 输入只覆盖 400–700："error" 抛 ValueError；"zero" 在 380 处为 0；"hold" 在 380 处等于 400 处的值
def test_cmf_shape_and_peak():       # shape (81, 3)；ybar 在 555 nm 处 ≈ 1.0 (abs 2e-3)
def test_illuminant_normalised():    # 对 D65, D50, A, E：sum(spd * ybar) == approx(1.0, rel 1e-12)
def test_illuminant_chromaticity():  # D65 xy≈(0.3127, 0.3290)，A xy≈(0.4476, 0.4074)，E xy≈(1/3, 1/3)，各 abs 1e-3
def test_blackbody_2856_matches_A(): # xy 与 illuminant("A") 相差 < 2e-3
def test_unknown_illuminant():       # illuminant("X", grid) 抛 ValueError，信息里列出可用名字
def test_srgb_matrices_inverse():    # XYZ_TO_LINEAR_SRGB @ LINEAR_SRGB_TO_XYZ ≈ I (abs 1e-6)
```

- [ ] **Step 2: 运行确认失败**：`uv run pytest tests/test_units.py tests/test_spectral.py`
- [ ] **Step 3: 实现 `units.py`。** `parse_duration_s` 接受 `m:ss`、`h:mm:ss`、带 `s`/`min`/`h` 后缀的数字、裸数字（秒）；分秒位须小于 60。
- [ ] **Step 4: 实现 `spectral.py`。** CIE 数据从 colour-science 取并重采样到网格，不手抄数表。光源的"按波段积分"值 = 光谱功率密度 × `step_nm`，再做归一化。
- [ ] **Step 5: 运行确认通过**，并 `uv run ruff check src tests`。
- [ ] **Step 6: 提交**：`git add src/darkroom_atlas/units.py src/darkroom_atlas/spectral.py tests/test_units.py tests/test_spectral.py && git commit -m "feat: units, constants and spectral grid with CIE data"`

---

### Task 3: 数据模型、来源标签与数据库加载

**Files:**
- Create: `src/darkroom_atlas/data/__init__.py`, `provenance.py`, `models.py`, `library.py`
- Create: `data/chemistry/species.yaml`, `data/films/ideal.yaml`, `data/processes/identity.yaml`
- Create: `tests/fixtures/library/`（见 Step 2）
- Test: `tests/data/test_provenance.py`, `tests/data/test_models.py`, `tests/data/test_library.py`

**Interfaces:**
- Consumes: `units.convert`、`units.parse_duration_s`、`units.parse_dilution`、`spectral.SpectralGrid`、`spectral.cmf_xyz`
- Produces（`provenance.py`，全部 `StrEnum` / pydantic 模型，`frozen=True, extra="forbid"`）:
  - `Source`: `datasheet literature fitted assumed`
  - `FormulaConfidence`: `published patent msds-estimated`
  - `VerificationStatus`: `validated predicted failed`
  - `CurveOrigin`: `datasheet measured`
  - `Quantity(value: float | list[float], unit: str, source: Source, ref: str | None = None, uncertainty: float | None = None, note: str | None = None)`；校验：来源为 `datasheet`/`literature` 时 `ref` 非空；`uncertainty >= 0`。方法 `to(unit: str) -> float | np.ndarray`。
- Produces（`models.py`，每个顶层文件模型带 `schema_version: Literal[1]`）:
  - 枚举 `SpeciesRole`（`developing_agent color_developing_agent preservative alkali buffer restrainer solvent acid fixer bleach coupler other`）、`BathKind`（`developer first_developer color_developer stop fixer bleach bleach_fix reversal wash stabilizer toner`）、`MaterialKind`（`bw_negative color_negative color_reversal bw_paper color_paper ideal`）、`StepKind`（`identity presoak develop stop fix bleach bleach_fix reversal wash stabilize tone dry`）
  - `Species(id, name, formula=None, cas=None, molar_mass: Quantity, roles: list[SpeciesRole], pka: list[Quantity] = [], diffusion_coefficient: Quantity | None = None, parameters: dict[str, Quantity] = {})`；`SpeciesFile(schema_version, species: list[Species])`，id 不得重复
  - `Bath(schema_version, id, name, kind: BathKind, formula_confidence: FormulaConfidence, components: list[BathComponent], ph: Quantity | None = None, standard_dilutions: list[str] = ["stock"], refs: list[str] = [], note=None)`；`BathComponent(species: str, amount: Quantity)`；`standard_dilutions` 每项须能被 `parse_dilution` 解析
  - `SpectralCurve(quantity: str, source: Source, ref: str | None = None, inline: InlineCurve | None = None, file: str | None = None, builtin: str | None = None)`，`inline / file / builtin` 恰好给一个；`InlineCurve(wavelengths_nm: list[float], values: list[float])` 等长且波长严格递增
  - 层（按 `type` 判别的联合，列表顺序为从上到下）：`OvercoatLayer(type="overcoat", id, thickness=None)`、`EmulsionLayer(type="emulsion", id, thickness: Quantity | None = None, sensitized_to: str | None = None, grain: Grain | None = None, spectral_sensitivity: SpectralCurve | None = None, couplers: list[Coupler] = [], parameters: dict[str, Quantity] = {})`、`FilterLayer(type="filter", id, spectral_density: SpectralCurve | None = None, removed_in_processing: bool = True)`、`InterlayerLayer(type="interlayer", id, thickness=None, scavenger: bool = False)`、`AntihalationLayer(type="antihalation", id, kind: Literal["undercoat","dyed_base","remjet","backing"], spectral_density=None, removed_in_processing: bool = True)`、`BaseLayer(type="base", id, material: str, thickness: Quantity | None = None, refractive_index: Quantity | None = None, spectral_density=None)`
  - `Grain(morphology: Literal["conventional","tabular","cubic"], mean_diameter: Quantity, sigma_log: Quantity, silver_coverage: Quantity, halide: str | None = None)`
  - `Coupler(id, dye: Literal["cyan","magenta","yellow"], kind: Literal["image","dir","colored_mask"], dye_spectral_density: SpectralCurve | None = None, parameters: dict[str, Quantity] = {})`
  - `Material(schema_version, id, maker, name, kind: MaterialKind, model: str, iso_speed: Quantity | None = None, default_readout: str | None = None, layers: list[Layer], refs: list[str] = [], note=None)`；校验：层 id 唯一、至少一个 emulsion 层、`kind != ideal` 时恰好一个 base 层；属性 `emulsion_layers`
  - `Agitation(method: Literal["inversion","rotary","stand","continuous","tray"], initial_s: float = 0.0, interval_s: float | None = None, duration_s: float | None = None)`、`Tank(volume_ml: float, film_area_cm2: float | None = None)`
  - `ProcessStep(kind: StepKind, bath: str | None = None, dilution: str = "stock", time: str | float | None = None, temperature: Quantity | None = None, agitation: Agitation | None = None, tank: Tank | None = None, note=None)`，属性 `time_s -> float | None`、`dilution_fraction -> float`；校验：`develop fix bleach bleach_fix tone` 必须有 `bath`；除 `identity` 外必须有 `time` 且大于 0
  - `Process(schema_version, id, name, family: str, steps: list[ProcessStep], refs=[], note=None)`，至少一步
  - `ReferenceCurve(material: str, quantity: str, x_unit: str, y_unit: str, origin: CurveOrigin, ref: str, conditions: str | None, x: np.ndarray, y: np.ndarray)`
- Produces（`library.py`）:
  - `class DataError(Exception)`（信息以文件路径开头）、`class DataNotFound(DataError)`（信息列出可用 id）
  - `default_data_root() -> Path`：环境变量 `DARKROOM_ATLAS_DATA`，否则仓库根的 `data/`（`Path(__file__).resolve().parents[3] / "data"`）
  - `class Library(root: Path | None = None)`，方法 `material(id)`（`films/<id>.yaml`，没有再找 `papers/<id>.yaml`）、`bath(id)`（`developers/<id>.yaml`，没有再找 `baths/<id>.yaml`）、`process(id)`（`processes/<id>.yaml`）、`species() -> dict[str, Species]`（文件缺失返回空字典）、`reference_curve(relpath) -> ReferenceCurve`（`reference/<relpath>`）、`resolve_curve(curve: SpectralCurve, grid) -> np.ndarray`（`builtin` 支持 `cie1931_x cie1931_y cie1931_z`；`inline` 与 `file` 重采样到网格，范围外取 0）、`ids(kind: Literal["material","bath","process"]) -> list[str]`、`path_of(kind, id) -> Path`、`sha256(path) -> str`、`validate() -> list[str]`
  - id 可含斜杠（`kodak/tri-x-400`），文件里的 `id` 必须等于由路径得出的 id。
  - 参考曲线 CSV 格式：文件开头若干行 `# key: value`（`material quantity x_unit y_unit origin ref conditions`），随后是无表头的 `x,y` 数据行。

- [ ] **Step 1: 写数据文件。**

```yaml
# data/films/ideal.yaml
schema_version: 1
id: ideal
maker: none
name: Ideal linear XYZ recorder
kind: ideal
model: ideal
default_readout: ideal_colorimetric
note: Not a real material. Records CIE XYZ linearly; used to verify the spectral pipeline.
layers:
  - {type: emulsion, id: x, spectral_sensitivity: {quantity: relative_energy_sensitivity, source: assumed, builtin: cie1931_x}}
  - {type: emulsion, id: y, spectral_sensitivity: {quantity: relative_energy_sensitivity, source: assumed, builtin: cie1931_y}}
  - {type: emulsion, id: z, spectral_sensitivity: {quantity: relative_energy_sensitivity, source: assumed, builtin: cie1931_z}}
```

```yaml
# data/processes/identity.yaml
schema_version: 1
id: identity
name: Identity (no chemistry)
family: ideal
steps:
  - {kind: identity}
```

`data/chemistry/species.yaml`：`schema_version: 1` 与 `species: []`。

- [ ] **Step 2: 写测试用的假数据库** `tests/fixtures/library/`，名字一律以 `test-` 开头、来源一律 `assumed`，以证明结构能同时容纳黑白与彩色：`chemistry/species.yaml`（3 个物种：一个显影剂、一个保护剂、一个彩色显影剂）、`developers/test-bw-dev.yaml`、`developers/test-color-dev.yaml`、`baths/test-bleach.yaml`、`baths/test-fix.yaml`、`films/test/bw-neg.yaml`（overcoat + 一个带 `grain` 的 emulsion + antihalation + base）、`films/test/color-neg.yaml`（overcoat、感蓝 emulsion、黄 filter、感绿 emulsion、interlayer、感红 emulsion、antihalation、base；每个 emulsion 带一个 image coupler，其中一层再带 dir 和 colored_mask 各一个，光谱敏感度用 inline 曲线）、`processes/test-bw.yaml`（develop 带 agitation 与 tank、stop、fix、wash）、`processes/test-c41.yaml`（develop、bleach、fix、wash、stabilize）、`reference/test/bw-neg/characteristic.csv`。

- [ ] **Step 3: 写失败的测试。**

```python
# tests/data/test_provenance.py
def test_quantity_ok():                 # Quantity(value=1.5, unit="um", source="assumed").to("nm") == 1500.0
def test_quantity_list_to_array():      # value=[1, 2], unit="mm" → to("um") 为 np.ndarray [1000, 2000]
def test_ref_required():                # source="datasheet" 无 ref → ValidationError；source="literature", ref="" 同样
def test_negative_uncertainty():        # uncertainty=-1 → ValidationError
def test_extra_field_forbidden():       # Quantity(..., foo=1) → ValidationError
def test_frozen():                      # 赋值 q.value = 2 抛错

# tests/data/test_models.py
def test_layer_union_dispatch():        # {"type": "base", ...} 解析为 BaseLayer；未知 type → ValidationError
def test_material_layer_ids_unique():   # 重复层 id → ValidationError
def test_material_requires_emulsion():  # 没有 emulsion 层 → ValidationError
def test_real_material_requires_base(): # kind=bw_negative 且无 base → ValidationError；kind=ideal 无 base 通过
def test_spectral_curve_exactly_one():  # 三者全空或给两个 → ValidationError
def test_inline_curve_monotonic():      # 波长不递增或长度不等 → ValidationError
def test_step_requires_bath():          # kind=develop 无 bath → ValidationError
def test_step_requires_time():          # kind=wash 无 time → ValidationError；kind=identity 无 time 通过
def test_step_time_and_dilution():      # time="9:45", dilution="1+1" → time_s == 585.0, dilution_fraction == 0.5
def test_bad_dilution_in_bath():        # standard_dilutions=["1:1"] → ValidationError
def test_species_ids_unique():

# tests/data/test_library.py   （fx = Library(tests/fixtures/library)；repo = Library()）
def test_repo_library_valid():          assert repo.validate() == []
def test_fixture_library_valid():       assert fx.validate() == []
def test_load_ideal():                  # repo.material("ideal").model == "ideal"；3 个 emulsion 层
def test_load_color_negative():         # fx.material("test/color-neg") 有 3 个 emulsion 层、1 个 filter、1 个 base；couplers 齐全
def test_ids_listing():                 # fx.ids("material") == ["test/bw-neg", "test/color-neg"]
def test_not_found_lists_available():   # fx.material("nope") 抛 DataNotFound，信息含 "test/bw-neg"
def test_id_must_match_path(tmp_path):  # 文件 films/a.yaml 里 id: b → DataError，信息含文件路径
def test_unknown_field_reports_path(tmp_path):   # 多一个字段 → DataError，信息含文件路径和字段名
def test_missing_ref_reports_path(tmp_path):     # datasheet 来源缺 ref → DataError，信息含文件路径
def test_validate_dangling_bath(tmp_path):       # 流程引用不存在的药液 → validate() 返回含该药液 id 的一条
def test_validate_dangling_species(tmp_path):    # 药液引用不存在的物种 → validate() 返回一条
def test_validate_missing_curve_file(tmp_path):  # SpectralCurve.file 指向不存在的文件 → validate() 返回一条
def test_resolve_builtin():             # resolve_curve(builtin cie1931_y) 等于 cmf_xyz(grid)[:, 1]
def test_resolve_inline_zero_outside(): # 只覆盖 400–700 的 inline 曲线在 380 处为 0
def test_reference_curve():             # 读出 origin == CurveOrigin.DATASHEET，x、y 等长
def test_sha256_stable():               # 同一文件两次相同，64 个十六进制字符
def test_env_override(monkeypatch):     # DARKROOM_ATLAS_DATA 指向 fixtures 时 default_data_root() 返回它
```

- [ ] **Step 4: 运行确认失败**：`uv run pytest tests/data`
- [ ] **Step 5: 实现 `provenance.py`、`models.py`。**
- [ ] **Step 6: 实现 `library.py`。** pydantic 的 `ValidationError` 一律包成 `DataError`，信息形如 `<path>: <字段路径>: <问题>`，每个问题一行。`validate()` 遍历库内全部文件，收集加载错误和交叉引用错误，不在第一个错误处停下。
- [ ] **Step 7: 运行确认通过**，并 `uv run ruff check src tests`。
- [ ] **Step 8: 提交**：`git add src/darkroom_atlas/data data tests/data tests/fixtures && git commit -m "feat: data models with provenance and library loader"`

---

### Task 4: 确定性随机数与分块

**Files:**
- Create: `src/darkroom_atlas/rng.py`, `src/darkroom_atlas/tiling.py`
- Test: `tests/test_rng.py`, `tests/test_tiling.py`

**Interfaces:**
- Produces（`rng.py`）:
  - `stream_key(seed: int, stream: str) -> int`：`blake2b(f"{seed}:{stream}".encode(), digest_size=8)` 的整数值
  - `pixel_uniform(seed: int, stream: str, rows: tuple[int, int], cols: tuple[int, int], frame_width: int, trailing: tuple[int, ...] = ()) -> np.ndarray`：返回形状 `(rows[1]-rows[0], cols[1]-cols[0], *trailing)`、取值 `[0, 1)` 的 `float64`。`rows`、`cols` 是画面绝对坐标的半开区间。
- Produces（`tiling.py`）:
  - `@dataclass(frozen=True) class Tile(outer: tuple[slice, slice], inner: tuple[slice, slice], inner_in_outer: tuple[slice, slice])`
  - `iter_tiles(shape_hw: tuple[int, int], tile_px: int, overlap_px: int) -> Iterator[Tile]`
  - `map_tiles(fn: Callable[[np.ndarray, Tile], np.ndarray], array: np.ndarray, tile_px: int, overlap_px: int) -> np.ndarray`

- [ ] **Step 1: 写失败的测试。**

```python
# tests/test_rng.py
def test_deterministic():            # 同参数两次调用 array_equal
def test_range_and_mean():           # 全部在 [0, 1)；200×300 的均值在 0.5 ± 0.01
def test_tile_equals_crop_of_full(): # full = pixel_uniform(7, "s", (0, 50), (0, 80), 80, (3,))；
                                     # part = pixel_uniform(7, "s", (10, 30), (20, 60), 80, (3,))；array_equal(part, full[10:30, 20:60])
def test_streams_and_seeds_differ(): # 换 stream 或换 seed 结果不同
def test_rows_independent():         # 第 0 行与第 1 行不同
def test_rejects_bad_ranges():       # cols 超出 frame_width、区间反向、负数 → ValueError

# tests/test_tiling.py
def test_inner_tiles_partition():    # 对 (50, 70) 配 tile=16, overlap=5：用计数数组验证每个像素恰好被 inner 覆盖一次
def test_outer_clipped_and_contains_inner():   # outer 不越界；outer 在内部边上比 inner 多出 overlap
def test_inner_in_outer_consistent():          # array[outer][inner_in_outer] 与 array[inner] 相等
def test_frame_smaller_than_tile():  # (10, 12) 配 tile=64 → 恰好一个块，outer == inner == 全画面
def test_overlap_larger_than_tile(): # tile=8, overlap=20 仍正确分割
def test_rejects_bad_args():         # tile_px < 1 或 overlap_px < 0 → ValueError
def test_map_identity():             # map_tiles(lambda a, t: a, x, 16, 3) 与 x 相等
def test_map_seamless_blur():        # x 为 (60, 90) 随机数组；f = gaussian_filter(sigma=2, truncate=3, mode="nearest")（半径 6）
                                     # map_tiles(f, x, 16, 6) 与 f(x) 的最大绝对差 < 1e-12
def test_map_insufficient_overlap_differs():   # 同上但 overlap=1，最大绝对差 > 1e-6
def test_map_trailing_dims():        # fn 把 (h, w, 3) 变成 (h, w, 5)，输出形状 (H, W, 5)
```

- [ ] **Step 2: 运行确认失败**：`uv run pytest tests/test_rng.py tests/test_tiling.py`
- [ ] **Step 3: 实现 `rng.py`。** 为了在不同 NumPy 版本间保持逐位一致，只使用位生成器的原始输出，不用 `Generator` 的分布方法：对每个绝对行号 `r`，

```python
bitgen = np.random.PCG64(np.random.SeedSequence([seed, stream_key(seed, stream), r]))
raw = bitgen.random_raw(frame_width * prod(trailing))           # uint64
row = ((raw >> np.uint64(11)) * 2.0**-53).reshape(frame_width, *trailing)[cols[0]:cols[1]]
```

`seed` 必须是非负整数，否则 `ValueError`。

- [ ] **Step 4: 实现 `tiling.py`。** inner 块按 `tile_px` 网格划分画面；outer = inner 向四周各扩 `overlap_px` 并裁到画面内。`map_tiles` 用第一个块的输出确定结果的尾部维度和 dtype。
- [ ] **Step 5: 运行确认通过**，并 `uv run ruff check src tests`。
- [ ] **Step 6: 提交**：`git add src/darkroom_atlas/rng.py src/darkroom_atlas/tiling.py tests/test_rng.py tests/test_tiling.py && git commit -m "feat: coordinate-keyed RNG and seamless tiling"`

---

### Task 5: Light、EmulsionState 与场景

**Files:**
- Create: `src/darkroom_atlas/light.py`, `src/darkroom_atlas/state.py`
- Create: `src/darkroom_atlas/scene/__init__.py`, `scene/io.py`, `scene/upsample.py`, `scene/camera.py`
- Test: `tests/test_light.py`, `tests/test_state.py`, `tests/scene/test_io.py`, `tests/scene/test_upsample.py`, `tests/scene/test_camera.py`

**Interfaces:**
- Consumes: `units.*`、`spectral.*`、`tiling.Tile`、`data.models.Material`、`data.library.Library`
- Produces（`light.py`）:
  - `@dataclass(frozen=True) class Light(photon_flux: np.ndarray, grid: SpectralGrid, meta: Mapping[str, Any] = {})`；`photon_flux` 形状 `(H, W, B)`，单位 光子 / µm² / s / 波段，非负有限，否则 `ValueError`
  - 属性 `shape_hw`；方法 `filtered(transmittance) -> Light`（`(B,)` 或 `(H, W, B)`，取值须在 `[0, 1]`）、`luminous_exposure_lux_s(time_s: float) -> np.ndarray`（形状 `(H, W)`；`683 × Σ_b 光子数 × 光子能量 × 1e12 × ybar`）
- Produces（`state.py`）:
  - `@dataclass(frozen=True) class FrameGeometry(shape_hw: tuple[int, int], pixel_pitch_um: float, origin_px: tuple[int, int] = (0, 0), frame_shape_hw: tuple[int, int] | None = None)`；`from_long_edge(shape_hw, long_edge_mm=36.0)`；`um_to_px(um) -> float`。`origin_px` 与 `frame_shape_hw` 记录裁剪后的块在整幅画面中的位置和整幅尺寸（未裁剪时 `frame_shape_hw` 等于 `shape_hw`）。
  - `@dataclass(frozen=True) class StepRecord(index: int, kind: str, bath: str | None, time_s: float | None, note: str = "")`
  - `@dataclass(frozen=True) class EmulsionState(material, geometry, grid, seed: int, fields: Mapping[str, np.ndarray] = {}, history: tuple[StepRecord, ...] = (), meta: Mapping[str, Any] = {})`
    - 字段键为 `"<layer_id>/<name>"`；与层无关的用 `"frame/<name>"`。数组前两维必须等于 `geometry.shape_hw`，存入时设为只读。
    - `get(layer_id, name) -> np.ndarray`（缺失时 `KeyError` 列出已有键）、`has(layer_id, name) -> bool`、`with_fields(updates: Mapping[str, np.ndarray]) -> EmulsionState`、`with_history(record) -> EmulsionState`、`with_meta(**kw) -> EmulsionState`
    - `crop(window: tuple[slice, slice]) -> EmulsionState`：裁所有字段，更新 `geometry.shape_hw` 与 `origin_px`（多次裁剪累加）
    - `paste(tile_state: EmulsionState, tile: Tile) -> EmulsionState`：把 `tile_state` 各字段的 `tile.inner_in_outer` 部分写入自身副本的 `tile.inner` 位置；自身没有的字段先按零数组创建（尾部维度和 dtype 取自 `tile_state`）
- Produces（`scene/io.py`）:
  - `srgb_decode(encoded)`、`srgb_encode(linear)`（IEC 61966-2-1 分段函数）
  - `@dataclass(frozen=True) class ImageInput(linear_rgb: np.ndarray, bit_depth: int, encoding: Literal["srgb", "linear"], warnings: tuple[str, ...])`
  - `read_image(path) -> ImageInput`：8/16 位整数图按 sRGB 解码；浮点 TIFF 视为线性。8 位输入追加警告 `"8-bit display-referred input: highlight latitude of negative film will not be exercised"`。
  - `write_png16(path, linear_rgb) -> None`：裁到 `[0, 1]`，sRGB 编码，量化为 16 位。
- Produces（`scene/upsample.py`）: `mallett2019_basis(grid) -> (B, 3)`（取自 colour-science 的 Mallett & Yuksel 2019 sRGB 基函数）、`rgb_to_reflectance(linear_rgb, grid) -> (H, W, B)`
- Produces（`scene/camera.py`）:
  - `scene_light(linear_rgb, *, grid=DEFAULT_GRID, illuminant="D65", iso=100.0, ev=0.0, time_s=1/125, white_exposure_lux_s: float | None = None) -> Light`
  - 相对能量 `e = reflectance × illuminant(...)`（已按波段积分，理想白的 `Σ e·ybar == 1`）。能量标度 `k`（J/µm²，对应相对能量 1）：给了 `white_exposure_lux_s` 时 `k = white_exposure_lux_s / 683 × 1e-12`；否则按平均测光，使全画面平均光照量等于 `(10 / iso) × 2**ev` lux·s（ISO 2721）。画面平均相对亮度小于 `1e-9` 时按 0.18 计。
  - `photon_flux = k × e / photon_energy / time_s`
  - `meta` 含 `illuminant iso ev time_s metering`（`"average"` 或 `"fixed"`）与 `energy_scale_j_per_um2`（即 `k`）

- [ ] **Step 1: 写失败的测试。**

```python
# tests/test_light.py
def test_luminous_exposure_anchor():   # 1×1 画面，仅 555 nm 波段有 4090.7 光子/µm²/s，time_s=1 → 1.000 lux·s (rel 2e-3)
def test_filtered_scalar_band():       # transmittance 全 0.5 → photon_flux 减半；meta 保留
def test_filtered_rejects_out_of_range():   # 1.2 或 -0.1 → ValueError
def test_rejects_bad_flux():           # 负值、NaN、形状不是 3 维、波段数与网格不符 → ValueError
def test_filtered_array_api_strict():  # 按 Array API 写的核函数接受 array_api_strict 输入

# tests/test_state.py   （material = Library().material("ideal")）
def test_geometry_from_long_edge():    # (4000, 6000), 36 mm → pitch 6.0；um_to_px(12) == 2.0
def test_with_fields_readonly_and_immutable():   # 返回新对象；原对象不变；数组不可写
def test_field_shape_validated():      # 前两维不符 → ValueError
def test_get_missing_lists_keys():     # KeyError 信息含已有键
def test_crop_tracks_origin():         # crop((slice(10, 30), slice(5, 25))) 后 shape (20, 20)，origin (10, 5)；再 crop 累加；frame_shape_hw 不变
def test_paste_roundtrip():            # 把整幅状态按 iter_tiles 逐块 crop 再 paste 回空状态，各字段与原状态相等
def test_history_and_meta():

# tests/scene/test_io.py
def test_srgb_roundtrip():             # encode(decode(x)) ≈ x，x 取 0..1 的 256 个点 (abs 1e-12)
def test_srgb_known_values():          # decode(0.5) ≈ 0.21404 (abs 1e-5)；decode(0.04045) ≈ 0.0031308 (abs 1e-7)
def test_read_sample_png():            # shape (64, 96, 3)，float64，范围 [0, 1]，bit_depth 8，有 8 位警告
def test_read_grayscale(tmp_path):     # 灰度 PNG → 三通道相等
def test_read_rgba_drops_alpha(tmp_path):
def test_read_png16(tmp_path):         # 值 65535 → 1.0；bit_depth 16；无 8 位警告
def test_read_float_tiff_linear(tmp_path):   # float32 TIFF 原值读出，encoding "linear"
def test_read_rejects_nonfinite(tmp_path):   # 含 NaN 的 TIFF → ValueError，信息含文件名
def test_write_png16_roundtrip(tmp_path):    # 写后读回与原线性值相差 < 1e-4

# tests/scene/test_upsample.py
def test_basis_partition_of_unity():   # 三条基函数之和在每个波段 ≈ 1 (abs 1e-3)
def test_white_is_flat():              # rgb=(1,1,1) → 反射率 ≈ 1；rgb=0 → 0
def test_roundtrip_through_cmf():      # 对 8 个色块：反射率 × D65 → XYZ → 线性 sRGB，与输入相差 < 5e-3

# tests/scene/test_camera.py
def test_average_metering_hits_iso():  # 任意图，iso=100 → 平均 luminous_exposure == 0.1 lux·s (rel 1e-9)；iso=400 → 0.025
def test_ev_doubles():                 # ev=+1 → 平均光照量加倍
def test_fixed_white_exposure():       # white_exposure_lux_s=1.0，纯白像素 → 1.0 lux·s (rel 1e-9)；18% 灰 → 0.18
def test_time_does_not_change_exposure():    # 改 time_s，flux 反比变化，光照量不变
def test_black_frame_metered_as_grey():      # 全黑图不抛错，flux 全 0，meta["energy_scale_j_per_um2"] 有限且为正
def test_meta_recorded():
def test_tungsten_shifts_spectrum():   # illuminant="A" 时白色像素的长波光子数多于短波
```

- [ ] **Step 2: 运行确认失败。**
- [ ] **Step 3: 实现 `light.py`、`state.py`。**
- [ ] **Step 4: 实现 `scene/`。** colour-science 里 Mallett 2019 数据集的确切名字以安装的版本为准，先在解释器里确认再写。
- [ ] **Step 5: 运行确认通过**，并 `uv run ruff check src tests`。
- [ ] **Step 6: 提交**：`git add src/darkroom_atlas/light.py src/darkroom_atlas/state.py src/darkroom_atlas/scene tests/test_light.py tests/test_state.py tests/scene && git commit -m "feat: light, emulsion state and scene-to-light conversion"`

---

### Task 6: 注册表、三个操作与理想模型

**Files:**
- Create: `src/darkroom_atlas/registry.py`, `src/darkroom_atlas/pipeline.py`, `src/darkroom_atlas/models/__init__.py`, `src/darkroom_atlas/models/ideal.py`
- Test: `tests/test_registry.py`, `tests/test_pipeline.py`, `tests/models/test_ideal.py`

**Interfaces:**
- Consumes: Task 2–5 的全部产出
- Produces（`registry.py`）:
  - `class RegistryError(LookupError)`（信息列出已注册的名字）
  - 装饰器 `register_exposure(model: str)`、`register_step(model: str, kind: StepKind, *, interaction_range_um: float | None = 0.0)`、`register_readout(name: str)`；重复注册同一键抛 `ValueError`
  - `get_exposure(model)`、`get_step(model, kind) -> RegisteredStep`（带属性 `fn` 与 `interaction_range_um`）、`get_readout(name)`
  - 签名约定：曝光 `fn(material, light, time_s, ctx: RunContext) -> EmulsionState`；步骤 `fn(state, step: ProcessStep, ctx: StepContext) -> EmulsionState`；出口 `fn(state, ctx: RunContext, **params) -> Readout`
  - `interaction_range_um` 的含义：`0.0` 表示逐像素；正数表示横向作用范围（µm），管线据此取分块重叠；`None` 表示该步骤有全局耦合（如有限体积药液槽），管线不替它分块，由处理器自己用 `tiling.map_tiles` 处理。
- Produces（`pipeline.py`）:
  - `@dataclass(frozen=True) class RunContext(library: Library, grid: SpectralGrid = DEFAULT_GRID, seed: int = 0, tile_px: int | None = None, long_edge_mm: float = 36.0)`
  - `@dataclass(frozen=True) class StepContext(run: RunContext, step_index: int, bath: Bath | None)`
  - `@dataclass(frozen=True) class Readout(kind: Literal["image_linear_srgb", "spectral_transmittance", "spectral_reflectance", "density"], data: np.ndarray, grid: SpectralGrid | None = None, meta: Mapping[str, Any] = {})`
  - `expose(material, light, time_s, ctx) -> EmulsionState`：校验 `light.grid == ctx.grid`、`time_s > 0`；按 `material.model` 分派；保证返回状态的 `geometry` 来自 `FrameGeometry.from_long_edge(light.shape_hw, ctx.long_edge_mm)`、`seed == ctx.seed`，并把 `light.meta` 存入 `state.meta["light"]`
  - `process(state, proc, ctx) -> EmulsionState`：逐步执行；按 `(state.material.model, step.kind)` 取处理器；有 `bath` 时从库里取出放入 `StepContext`；`ctx.tile_px` 不为空且处理器的 `interaction_range_um` 不为 `None` 时分块执行（重叠 `ceil(range_um / pixel_pitch_um)`），用 `state.crop` 与 `paste` 拼回；每步追加一条 `StepRecord`
  - `readout(state, device, ctx, **params) -> Readout`
- Produces（`models/ideal.py`，导入 `darkroom_atlas.models` 即完成注册；`pipeline.py` 导入它）:
  - 曝光模型 `"ideal"`：对每个 emulsion 层写入字段 `"<layer>/exposure"` = `Σ_b 光子数_b × 光子能量_b × 敏感度_b`（J/µm²，敏感度由 `library.resolve_curve` 得到）
  - 步骤 `("ideal", identity)`：`"<layer>/signal"` = `"<layer>/exposure"`，`interaction_range_um=0.0`
  - 出口 `"ideal_colorimetric"`：要求层 `x y z` 都有 `signal`；`XYZ = signal / state.meta["light"]["energy_scale_j_per_um2"]`，再乘 `XYZ_TO_LINEAR_SRGB`，返回 `Readout(kind="image_linear_srgb")`

- [ ] **Step 1: 写失败的测试。**

```python
# tests/test_registry.py
def test_ideal_registered():            # 导入 darkroom_atlas.pipeline 后 get_exposure("ideal")、get_step("ideal", StepKind.identity)、get_readout("ideal_colorimetric") 都存在
def test_unknown_lists_registered():    # get_exposure("nope") 抛 RegistryError，信息含 "ideal"
def test_duplicate_rejected():          # 同键注册两次 → ValueError
def test_step_range_stored():           # register_step(..., interaction_range_um=12.0) 后 get_step(...).interaction_range_um == 12.0

# tests/models/test_ideal.py   （lib = Library()；ctx = RunContext(lib)；img = read_image(sample.png).linear_rgb）
def test_exposure_fields():             # expose 后有 "x/exposure" "y/exposure" "z/exposure"，形状 (64, 96)
def test_y_exposure_matches_luminous(): # "y/exposure" × 683 × 1e12 等于 light.luminous_exposure_lux_s(time_s) (rel 1e-9)
def test_end_to_end_roundtrip():        # scene_light → expose → process(identity) → readout：与 img 的最大绝对差 < 5e-3，并打印实际误差
def test_roundtrip_independent_of_iso_ev_time():   # 改 iso、ev、time_s，输出不变 (abs 1e-9)
def test_tungsten_scene_is_warm():      # illuminant="A" 时灰色像素读出的 R > B

# tests/test_pipeline.py
def test_expose_sets_geometry_seed_meta():   # pitch == 36000/96；seed == ctx.seed；meta["light"]["iso"] 存在
def test_expose_rejects_grid_mismatch():     # light 用 10 nm 网格 → ValueError
def test_expose_rejects_nonpositive_time():
def test_process_records_history():          # history 长度 1，kind "identity"，index 0
def test_process_unknown_step_for_model():   # 对 ideal 材料跑 fixtures 里的 test-bw 流程 → RegistryError
def test_tiled_equals_untiled_pointwise():   # tile_px=17 与 None 的结果各字段 array_equal
def test_tiled_equals_untiled_with_range():  # 测试内临时注册 ("tiletest", identity)：对 "x/exposure" 做 gaussian_filter(sigma=2, truncate=3, mode="nearest")，
                                             # interaction_range_um = 6 × pixel_pitch_um；tile_px=16 与 None 的结果最大绝对差 < 1e-12
def test_tiled_rng_consistent():             # 测试内临时注册的步骤用 pixel_uniform（rows/cols 取自 geometry.origin_px，frame_width 取自 geometry.frame_shape_hw）
                                             # 写入一个字段；tile_px=16 与 None 的结果 array_equal
def test_global_step_not_tiled():            # interaction_range_um=None 的步骤即使 tile_px=16 也只被调用一次，收到整幅状态
def test_bath_resolved_into_context():       # 用 fixtures 库和测试内注册的 develop 处理器，ctx.bath.id == "test-bw-dev"
def test_readout_unknown_device():           # RegistryError 列出 "ideal_colorimetric"
def test_readout_as_next_light():            # 用 Light.filtered 把一个 (H, W, B) 透射率作用到均匀光上，再 expose 到 ideal，y/exposure 按透射率成比例
```

测试内的临时注册须在结束后清理，不得污染全局注册表。

- [ ] **Step 2: 运行确认失败。**
- [ ] **Step 3: 实现 `registry.py`、`pipeline.py`。**
- [ ] **Step 4: 实现 `models/ideal.py`。**
- [ ] **Step 5: 运行确认通过。** 若端到端往返误差超过 5e-3，不要放宽阈值，先查出原因（最可能是光源归一化或矩阵白点不一致）并在汇报里写明。
- [ ] **Step 6: 运行全部测试** `uv run pytest` 与 `uv run ruff check .`。
- [ ] **Step 7: 提交**：`git add src/darkroom_atlas/registry.py src/darkroom_atlas/pipeline.py src/darkroom_atlas/models tests/test_registry.py tests/test_pipeline.py tests/models && git commit -m "feat: registry, expose/process/readout pipeline and ideal model"`

---

### Task 7: 运行记录与命令行

**Files:**
- Create: `src/darkroom_atlas/manifest.py`, `src/darkroom_atlas/runner.py`, `src/darkroom_atlas/cli.py`
- Test: `tests/test_manifest.py`, `tests/test_runner.py`, `tests/test_cli.py`

**Interfaces:**
- Consumes: Task 2–6 的全部产出
- Produces（`manifest.py`）:
  - `FileRef(path: str, sha256: str)`
  - `RunManifest`（pydantic）：`schema_version: Literal[1]`、`atlas_version`、`created_utc`、`input: FileRef`、`material: DataRef`、`process: DataRef`、`baths: list[DataRef]`、`species_file: FileRef | None`、`readout: ReadoutRef(device, params)`、`exposure: dict`（即 `Light.meta`）、`geometry: dict`（`shape_hw long_edge_mm pixel_pitch_um`）、`spectral_grid: dict`、`seed`、`tile_px`、`environment: dict`（`python numpy scipy colour platform`）、`git_commit: str | None`、`output: FileRef`、`warnings: list[str]`；`DataRef(id: str, file: FileRef)`
  - `recipe_hash() -> str`：对规范化 JSON 取 sha256，排除 `created_utc`、`environment`、`git_commit`、`tile_px` 以及 `input.path`、`output.path`（保留它们的 sha256）
  - `write(path)`、`RunManifest.read(path)`
- Produces（`runner.py`）: `run_image(input_path: Path, output_dir: Path, *, material_id: str, process_id: str, readout_device: str | None = None, illuminant: str = "D65", iso: float | None = None, ev: float = 0.0, time_s: float = 1/125, seed: int = 0, tile_px: int | None = None, long_edge_mm: float = 36.0, library: Library | None = None) -> tuple[Path, RunManifest]`。`iso` 为空时取材料的 `iso_speed`，再没有则 100。`readout_device` 为空时取材料的 `default_readout`，再没有则抛 `ValueError`。输出文件名 `<输入文件主名>__<material>__<process>.png`（id 中的 `/` 换成 `-`），运行记录为同名加 `.run.json`。输出目录不存在则创建。
- Produces（`cli.py`）: typer 应用 `app`
  - `atlas run INPUT -m/--material ID -p/--process ID -o/--output DIR [-r/--readout NAME] [--illuminant D65] [--iso N] [--ev X] [--time S] [--seed 0] [--tile N] [--long-edge-mm 36] [--data PATH]`：成功时打印输出路径、运行记录路径和警告；`DataNotFound`、`RegistryError`、输入文件不存在时向 stderr 打印一行信息并以退出码 2 结束，不打印堆栈
  - `atlas check-data [--data PATH]`：打印 `Library.validate()` 的每条问题；有问题退出码 1，否则打印 `OK` 加各类文件数量并退出 0

- [ ] **Step 1: 写失败的测试。**

```python
# tests/test_runner.py
def test_run_creates_outputs(tmp_path):      # 输出 sample__ideal__identity.png 与 .run.json 都存在；输出目录原本不存在
def test_output_roundtrip(tmp_path):         # 读回输出，与输入线性值最大绝对差 < 5e-3
def test_rerun_bit_identical(tmp_path):      # 同种子跑到两个目录，PNG 字节相同，recipe_hash 相同
def test_tiled_bit_identical(tmp_path):      # tile_px=32 与不分块的 PNG 字节相同，recipe_hash 相同
def test_seed_changes_recipe_hash(tmp_path): # seed 不同 → recipe_hash 不同
def test_manifest_contents(tmp_path):        # atlas_version == ATLAS_VERSION；material.id == "ideal"；material.file.sha256 等于库里该文件的 sha256；
                                             # exposure["iso"] == 100；geometry["shape_hw"] == [64, 96]；warnings 含 8 位输入的警告；output.sha256 等于实际文件
def test_slash_ids_in_filename(tmp_path):    # 用 fixtures 库里复制出的 test/ideal 材料 → 文件名含 "test-ideal"
def test_no_default_readout_errors(tmp_path):

# tests/test_manifest.py
def test_write_read_roundtrip(tmp_path)
def test_recipe_hash_ignores_volatile():     # 改 created_utc、environment、git_commit、tile_px、路径 → 不变；改 seed 或 material.file.sha256 → 改变

# tests/test_cli.py   （typer.testing.CliRunner）
def test_run_ok(tmp_path):                   # 退出码 0，stdout 含输出文件名
def test_run_unknown_material(tmp_path):     # 退出码 2，输出含 "ideal"（可用名字），不含 "Traceback"
def test_run_unknown_process(tmp_path):      # 退出码 2
def test_run_missing_input(tmp_path):        # 退出码 2
def test_check_data_ok():                    # 退出码 0，含 "OK"
def test_check_data_reports(tmp_path):       # --data 指向含一个坏文件的目录 → 退出码 1，输出含该文件名
def test_entry_point_subprocess(tmp_path):   # subprocess 运行 `uv run atlas run tests/assets/sample.png --material ideal --process identity -o <tmp>`，退出码 0
```

- [ ] **Step 2: 运行确认失败。**
- [ ] **Step 3: 实现 `manifest.py`、`runner.py`、`cli.py`。** `git_commit` 用 `git rev-parse HEAD`，失败时为 `None`。
- [ ] **Step 4: 运行全部测试** `uv run pytest` 与 `uv run ruff check .`。
- [ ] **Step 5: 手工验证。**

Run: `uv run atlas run tests/assets/sample.png --material ideal --process identity -o out/`
Expected: 打印 `out/sample__ideal__identity.png` 和 `.run.json` 的路径，以及 8 位输入的警告。

Run: `uv run atlas check-data`
Expected: `OK`，退出码 0。

- [ ] **Step 6: 提交**：`git add src/darkroom_atlas/manifest.py src/darkroom_atlas/runner.py src/darkroom_atlas/cli.py tests/test_manifest.py tests/test_runner.py tests/test_cli.py && git commit -m "feat: run manifest, runner and atlas CLI"`

---

### Task 8: README、物理文档骨架与版本同步检查

**Files:**
- Create: `README.md`, `docs/physics.md`
- Test: `tests/test_version_sync.py`

**Interfaces:**
- Consumes: `darkroom_atlas.ATLAS_VERSION`

- [ ] **Step 1: 写失败的测试** `tests/test_version_sync.py`：

```python
def test_readme_badge_matches_version():      # README.md 含 f"https://img.shields.io/badge/version-{ATLAS_VERSION}-blue"
def test_badge_is_first_line_after_title():   # 第一行以 "# " 开头；其后第一个非空行以 "![version]" 开头
def test_top_table_row_matches_version():     # "| 版本 | 日期 | 变更内容 | 类型 |" 表头及分隔行之后的第一行，其第一格等于 ATLAS_VERSION
```

- [ ] **Step 2: 运行确认失败。**
- [ ] **Step 3: 写 `README.md`（中文）。** 标题 `# The Darkroom Atlas`；标题下第一行只放一枚徽章 `![version](https://img.shields.io/badge/version-D0-blue)`（仓库没有远程，其余三枚 GitHub 动态徽章不加）。内容依次为：一段话说明项目是什么；"现在能做什么"（D0：底座与理想材料，尚无真实胶片）；快速开始（`uv sync`、`uv run pytest`、`uv run atlas run ...`、`uv run atlas check-data`）；数据目录与来源标签规则的简述；路线（D0–D5，各一行）；指向 `docs/physics.md` 和设计文档的链接；版本记录表：

```markdown
| 版本 | 日期 | 变更内容 | 类型 |
|------|------|----------|------|
| D0 | 2026-11-20 | 完成底座：单位与光谱约定、带来源标签的数据模型、expose/process/readout 管线与注册表、分块与确定性随机数、运行记录、命令行；理想材料走通全流程 | milestone |
```

- [ ] **Step 4: 写 `docs/physics.md`（中文）骨架。** 章节：1 约定（单位、光谱网格与"按波段积分"、数组轴顺序、字段命名 `<layer>/<name>`、随机数与可复现性、分块与 `interaction_range_um` 的三种取值）；2 链路总览（设计文档里的七个环节，各一小段，标明"D0 已实现"或"计划于 Dn"）；3 D0 已实现的物理（测光：ISO 2721 的 `H = 10/S`；光度量到光子数的换算公式；Mallett 2019 光谱还原及其局限；理想材料的定义和它验证了什么）；4 已知的不物理之处；5 参考文献（取自设计文档）。公式用行内代码或 LaTeX 均可，须与代码一致。
- [ ] **Step 5: 运行确认通过**：`uv run pytest tests/test_version_sync.py`
- [ ] **Step 6: 提交**：`git add README.md docs/physics.md tests/test_version_sync.py && git commit -m "docs: README with version record and physics conventions"`

---

## 完成标准

- `uv run pytest` 全绿，`uv run ruff check .` 无报错。
- `uv run atlas run tests/assets/sample.png --material ideal --process identity -o out/` 产出图像与运行记录；同种子重跑逐位相同；分块与不分块逐位相同。
- `uv run atlas check-data` 返回 OK。
- README 的 version 徽章、版本记录表首行与 `ATLAS_VERSION` 一致。
