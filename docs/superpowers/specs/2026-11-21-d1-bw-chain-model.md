# D1 黑白物理链路：模型规格

日期：2026-11-21。适用阶段：D1（D1a–D1f）。上位文档：`docs/superpowers/specs/2026-11-20-darkroom-atlas-design.md`（总体设计）、`docs/physics.md`（D0 约定）。本文件固定 D1 黑白链路的方程、参数、数值方法、数据格式、标定方案与模块接口；D1 的实现以本文件为准，偏离时先改本文件。

本规格写作期间做过的数值核对（脚本不入库）列在第 12 节，正文引用时写"（§12.x）"。

---

## 0 阅读说明

### 0.1 引用缩写

研究报告都在 `docs/research/d1/`，正文按"缩写 §节号"引用：

| 缩写 | 文件 |
|------|------|
| DS | `datasheets.md` |
| DEV | `developer-formulas.md` |
| LAT | `latent-image.md` |
| KIN | `development-kinetics.md` |
| OPT | `emulsion-optics.md` |
| CAL | `calibration-and-validation.md` |
| TMX | `t-max-developer.md` |

暂存的数字化数据在 `tools/digitize/staged/`，正文写作 `staged/kodak/tri-x-400/char-35mm-d76-08min.csv` 这样的相对路径。它们并入 `data/reference/` 后文件名不变（§5.9）。

### 0.2 来源类别

数据文件里每个物理参数的 `source` 取四类之一（D0 规定）。研究报告另有 opened / secondary / recalled / estimate 四个可信度标签，换算规则如下，本规格每个参数都按此标注：

| 报告标签 | 数据文件中的 `source` | 条件 |
|----------|----------------------|------|
| opened（厂商手册） | `datasheet` | 写明出版物、版次、页码 |
| opened（文献、专利、标准） | `literature` | 写明出处 |
| estimate，且全部输入都是 opened | 与输入相同的类别 | `ref` 写"derived: <算法>, from <出处>" |
| estimate 的输入里有 recalled 或 secondary | `assumed` | `note` 写明来历 |
| recalled、secondary | `assumed` | `note` 写"recalled/secondary per <报告 §>" |
| 由 D1 标定得到 | `fitted` | 带 `fit` 记录（§5.1） |

"先验"指标定时对参数施加的约束（区间或正态分布），先验本身也有来源类别，写在参数表"先验/界"一列。

### 0.3 记号与单位

- 胶片面坐标 `x = (x, y)`，单位 µm；像素间距 `p`（µm）。
- 乳剂内的深度用分数深度 `ζ ∈ [0, 1]`，`ζ = 0` 是乳剂上表面（靠保护层），`ζ = 1` 是乳剂与片基的界面。曝光（干膜）与显影（溶胀膜）用同一组分数深度单元，颗粒在分数深度上不动。
- 深度单元 `z = 1 … n_z`，`z = 1` 在最上面，中心 `ζ_z = (z − ½)/n_z`。
- 液体边界层内的高度 `h ≥ 0`，从保护层外表面量起。
- 波段 `b = 1 … B`（B = 81），尺寸档 `c = 1 … K`，团簇大小 `n`，输运组分 `k`。
- 方程写在代码块里，用通用数学记号；字段名、单位串、函数签名一律 ASCII。
- 单位串沿用 D0 的写法（`um` `s` `mol/L` 等），D1 新增的单位见 §5.1。
- 浓度一律指溶胀层（或液体）单位体积中的物质的量，`mol/L`。颗粒数密度用 `um^-3`（体）或 `um^-2`（面）。

---

## 1 范围、非目标与既定决策

### 1.1 范围

D1 实现一条完整、机理化的黑白负片链路：场景 → 光谱曝光（D0 已有）→ 乳剂内光传输 → 潜影 → 显影反应–扩散 → 银密度 → 停显、定影、水洗 → 颗粒 → 读出；以及虚拟感光仪、标定、验证报告和显微镜模式。

| 对象 | 本阶段内容 |
|------|------------|
| 胶片 | `kodak/tri-x-400`（传统多分散颗粒）、`kodak/t-max-100`（扁平颗粒），各含 135 与 120 两种画幅 |
| 显影液 | `d-76`（成品，msds-estimated，B₂O₃ 为有界干扰参数）、`d-76-published`（公开配方 2/100/5/2）、`xtol`（US 5,756,271 例 1）、`hc-110-pre2019`、`hc-110-2019`、`t-max-dev`（msds-estimated，世代 I）、`rodinal`（仅预测） |
| 其他药液 | 醋酸停显液、清水停显、硫代硫酸盐定影液、清水 |
| 出口 | 漫射视觉密度计、ISO Type 2 近似密度计、标准扫描 |

### 1.2 非目标

- 彩色、反转、相纸、印相（D2、D3）。
- 染色显影液（Pyro 类）的染色密度。
- Callier 效应、聚光扫描。标准扫描定义为漫射几何（§2.8）。
- 补充液、长期老化、空气氧化的定量模型（接口留出，D1 取 0）。
- GPU 后端（D5）。
- 推冲表（Kodak 不给推冲的 CI 目标，见 DS §2.2），只作定性对照。

### 1.3 既定决策（复述）

1. 胶片与显影液如 §1.1 表。稀释比一律写成显式比例（`stock`、`1+1`、`1+31`）。
2. 标定分工。胶片参数只用 D-76 特性曲线（加上胶片自己的光谱感光度形状、MTF；T-Max 100 另加 D-76 下测得的 RMS 颗粒度 8）。化学常数为各胶片共用，在训练组合上联合标定：Tri-X×D-76（曲线、CI(t)、时间表），Tri-X×XTOL 与 Tri-X×HC-110（只有 CI(t) 与时间表），Tri-X×T-MAX Developer（特性曲线），T-Max 100×D-76。留出组合从不参与任何拟合或设计选择：T-Max 100×XTOL、T-Max 100×HC-110（反差动力学检验：CI(t)、到达 CI 的时间，含 J-109 表 1 的 T-Max 行与 F-4016 的 HC-110 温度系数），T-Max 100×T-MAX Developer（完整特性曲线检验：趾部、形状、感光度）。Tri-X 的颗粒度 17 是在 HC-110 B 下测的，作为跨显影液检查，不是胶片标定目标。倒易律：潜影常数用 Kodak 表格的一部分标定，按 LAT §5 的建议留出若干行。训练组合内部也留出部分稀释/温度行（例如 1:1 的温度行）作为输运检查。
3. 验证阈值在任何留出评估之前预注册，暂用 CAL §4.2 的建议值，所有者看到第一份报告后确认或修改。每个胶片×显影液组合带状态：`validated`、`predicted` 或 `failed`。标定得到的常数一律标 `fitted`，可信度的主指标是留出组合的预测误差。
4. 物理结构（细化，不替换）：(1) 乳剂内光传输：逐层光谱吸收、散射作为点扩散函数、片基全反射造成的光晕；光学 PSF 不直接拟合手册 MTF；(2) 潜影：按尺寸分档的颗粒群体、泊松光子吸收、成核–生长马尔可夫链给出每档的团簇大小分布，倒易律失效由此产生，并带 LAT 认为必需的陷阱深度分布；(3) 显影：物种库、适用于离子强度约 2 M 的活度处理、带溴化物/碘化物竞争吸附的 Langmuir–Hinshelwood 速率律、Arrhenius、依赖潜影大小的诱导、通过再生的超加和性、亚硫酸盐作为清除剂和卤化银溶剂、pH 作为带缓冲的输运场（明胶为不动缓冲剂）、各物种在凝胶与水中的扩散、带搅动事件的液体边界层、有限体积的充分混合药液槽；(4) 成像：Nutting 公式，各档有遮盖力；(5) 停显、定影、水洗为简单一阶步骤，真实改变状态；(6) 颗粒：全幅统计模式按像素按档抽样，带 OPT 推出的方差修正 κ(D)；小块显微镜模式（布尔模型）作参照；(7) 出口：漫射视觉密度计（ISO Type 2 近似，标 `assumed`），以及规则只定义一次的标准扫描。
5. 软件约束（D0）：数组 `(H, W, ..., band)`、float64、380–780 nm 5 nm 波段积分、µm、每 µm² 每波段光子数、不可变 `EmulsionState`（字段 `"<layer>/<name>"`）、随机数一律经 `rng.py`、有限作用范围的步骤分块与不分块逐位相同、每次运行写运行记录、模型经注册表接入、YAML 数据带来源标签、uv/pytest/ruff。
6. D0 留下的两个工程问题由本规格解决：(a) 内存（全幅 `(H, W, 81)` 光数组、重复 paste 的平方代价）；(b) 药液槽耗尽造成的全幅耦合与分块的矛盾。
7. 计算预算：标定需要成千上万次感光条带的正向计算，因此要有快速感光路径（均匀色块不需要横向维度，每级一根深度柱），与图像路径共用完全相同的物理代码；24 MP 全幅在 14 核、48 GB 笔记本上可以花几分钟，不能花几小时；只用 NumPy/SciPy。
8. 分六个模块交付，各自可测并带一次中版本号：D1a 光传输；D1b 潜影；D1c 化学与显影；D1d 成像、颗粒统计、停显定影与出口；D1e 虚拟感光仪、指标、标定与验证报告（含胶片与显影液数据文件）；D1f 显微镜模式。

对其中几条的澄清和异议在第 11 节。

### 1.4 本规格做出的主要设计选择（索引）

| 选择 | 位置 | 理由（一句话） |
|------|------|----------------|
| 低秩光 `LowRankLight`：系数 `(H, W, 3)` 加基 `(3, B)` | §2.1、§4 | D0 的 RGB→光谱是线性的三基函数组合，`(H, W, 81)` 从不需要存在 |
| 光学先投影后卷积：PSF 与光晕核光谱中性，6 次卷积覆盖全部波段 | §2.1 | 卷积与逐波段加权可交换 |
| 两位点阻塞闭合加 `S_eff` 个竞争位点、对数正态 τ₀ | §2.2 | 单位点不能同时给出 LIRF 的斜率与 HIRF 的幅度（§12.2、§12.3） |
| 可分离的引发风险 λ = λ̃(化学)·h_n(团簇) | §2.4 | Sugai/James 的诱导期形式本身可分离；每像素只需一个"引发时钟" |
| 一阶生长 dx/dt = κ(1 − x)，跟踪 x 与 x² 的群体矩 | §2.4 | 群体方程封闭，每档两个状态；x² 供颗粒统计 |
| SIT 型活度，ε 只用 pH 锚点标定 | §2.3 | 只保留 Debye–Hückel 项已把高离子强度锚点拉到 0.1 pH 以内（§12.4） |
| 化学单元（固定 24 µm）与成员（像素）两层网格；成员按核分解移出时间积分，输出时核双线性插值 | §2.4.4、§3.8 | 化学场在凝胶交换长度 ≈15 µm 以下是平滑的；成本由化学单元数（画幅面积）决定，与像素数几乎无关 |
| 显影步骤内部用固定的规范块，与 `tile_px` 无关 | §3.6 | 带隐式横向扩散的时间推进无法在可变分块下逐位一致 |
| 药液槽：曝光直方图上的柱集合预计算 + 事后不动点检查 | §3.7 | 槽的耗尽只依赖全幅平均源项，与感光条带同一套代码 |
| 固定的、与事件对齐的几何时间步表 | §3.4 | 比自适应步长少一半步数而误差更小（§12.6），且给标定平滑的雅可比 |
| 大字段放在只读内存映射文件里 | §4 | 24 MP 的显影输出约 4 GB，内存只放工作块 |

---

## 2 物理链路

### 2.0 共用部分

#### 2.0.1 层结构

D1 的两种胶片都用同一种层结构（层 id 固定，字段名由它们组成）：

| 层 id | 类型 | 内容 |
|-------|------|------|
| `overcoat` | `overcoat` | 保护层明胶，厚度 `thickness`（干） |
| `emulsion` | `emulsion` | 唯一的感光层：颗粒群体、明胶、增感染料 |
| `base` | `base` | 片基；按画幅（`135`、`120`）给厚度、灰染密度、背涂层密度 |

120 的防光晕背涂层不单列一层，作为 `base.formats["120"].backing_density` 处理（§2.1.5）。T-Max 400 那样的双乳剂层结构不在 D1 范围内。

#### 2.0.2 颗粒尺寸分档

颗粒等效圆直径（ECD）按个数服从对数正态分布：

```
ln d ~ N(μ_d, σ_d²),  μ_d = ln d_50
档界：ln d_c,lo = μ_d + σ_d·(−3 + 6(c−1)/K),  ln d_c,hi = μ_d + σ_d·(−3 + 6c/K)
个数分数：w_c = [Φ(t_c,hi) − Φ(t_c,lo)] / [Φ(3) − Φ(−3)],  t = (ln d − μ_d)/σ_d
档中心：d_c = exp((ln d_c,lo + ln d_c,hi)/2)
```

`Φ` 是标准正态分布函数。分布在 ±3σ 截断后重新归一（截掉的个数分数为 0.27 %，`assumed`）。

单颗几何（`conventional` 按球，`tabular` 按厚度为 `t_g` 的圆片平躺）：

| 量 | conventional | tabular | 单位 |
|----|--------------|---------|------|
| 体积 `V_c` | π d_c³/6 | π d_c² t_g/4 | um³ |
| 表面积 `S_c` | π d_c² | π d_c²/2 + π d_c t_g | um² |
| 投影面积 `A_c` | π d_c²/4 | π d_c²/4 | um² |

银量与个数：

```
M_AgX = (1 − x_I)·187.77 + x_I·234.77                         [g/mol]
m_Ag,c = 1e-12 · ρ_AgX · V_c · 107.87 / M_AgX                 [g per grain]   (ρ 以 g/cm³ 计，1 g/cm³ = 1e-12 g/um³)
N_tot = 1e-12 · m_Ag,coat / Σ_c w_c m_Ag,c                    [um^-2]        (m_Ag,coat 以 g/m² 计)
ν_c = w_c · N_tot                                              [um^-2]
干乳剂厚度：L_e = m_Ag,coat·(M_AgX/107.87)/ρ_AgX + m_gel/ρ_gel   [um]       (g/m² 除以 g/cm³ 正好是 µm)
```

颗粒在分数深度上均匀分布（`assumed`），深度单元 `z` 中 c 档的面密度为 `ν_c/n_z`，溶胀后的体密度为 `n_c,z = ν_c/(n_z·Δ_z)`（`Δ_z` 为溶胀单元厚度，µm）。c 档在单元 z 中的银量浓度（以 Ag 计）：

```
C_c,z = 1e15 · n_c,z · m_Ag,c / 107.87        [mol/L]   (1 L = 1e15 um³)
```

ρ_AgX 取 AgBr 的 6.473 g/cm³（LAT §2.3，recalled → `assumed`），碘化物对密度的影响忽略；AgI 的摩尔质量 234.77 由原子量算出（`assumed`）。L_e 的公式见 LAT §2.3、OPT §1.4（estimate → `assumed`）。

参数：

| 符号 | 字段 | 单位 | 类别 | 先验/界 | 出处 |
|------|------|------|------|---------|------|
| d_50 | `grain.mean_diameter` | um | fitted | Tri-X 0.4–2.0；T-Max 0.6–1.5（t_g 也拟合，所以 T-Max 不直接拟合 d_50：按 (t_g, AR = d_50/t_g ∈ [8, 19]) 参数化，[0.6, 1.5] 为单侧罚项，§7.1） | LAT §2.1、§2.2 |
| σ_d | `grain.sigma_log` | 1 | fitted | 0.2–0.7 | LAT §2.1（教科书范围 0.3–0.6，recalled） |
| t_g | `grain.thickness` | um | fitted（T-Max；先验中心 0.12） | 0.08–0.15，且 d_50/t_g ≥ 8 | LAT §2.2；S/V ≈ 2/t_g + 4/d 是 T-Max 唯一的胶片侧显影速度杠杆（§2.4.4、§9 第 31 条） |
| x_I | `grain.iodide_fraction` | 1 | assumed | Tri-X 0.02–0.10（默认 0.04）；T-Max 0.01–0.06（默认 0.03） | LAT §2.1、§2.2 |
| m_Ag,coat | `grain.silver_coverage` | g/m2 | fitted（每种胶片；先验中心 4.0，对数） | 2–6 | LAT §2.3（无一手数据）、KIN §4.3；可辨识性见 §2.5.3、§7.5 |
| m_gel | `gelatin_coverage` | g/m2 | assumed | 3–8，默认 5.0 | LAT §2.3（Kofron 4.84） |
| K | `grain.n_classes` 与档位覆盖 | 1 | 离散化参数 | 标定与图像 6，参考 12 | §3.1 |

#### 2.0.3 状态字段总表

所有字段都是 float64，前两维为 `(H, W)`。"广播"表示只读的 `np.broadcast_to` 视图，不占内存（§4.3）。

| 字段 | 形状 | 单位 | 写入 | 读取 | 说明 |
|------|------|------|------|------|------|
| `emulsion/exposure_volume` | (H, W, n_z) | um^-3 | 曝光 | 显影 | e_V,z，每单位 AgX 体积吸收的光子数 |
| `emulsion/exposure_surface` | (H, W, n_z) | um^-2 | 曝光 | 显影 | e_S,z，每单位颗粒表面吸收的光子数 |
| `emulsion/residual_dye` | (H, W) 广播 | 1 | 曝光（1.0）、定影、水洗 | 出口 | 残留染料分数 |
| `frame/presoaked` | (H, W) 广播，bool | — | 预湿 | 显影 | 是否预湿 |
| `emulsion/silver` | (H, W, K) | g/m2 | 显影、停显 | 出口 | 化学显影银，按档、深度合并 |
| `emulsion/initiated` | (H, W, K) | 1 | 显影、停显 | 出口 | 已引发颗粒的个数分数 F_c |
| `emulsion/silver_m2` | (H, W, K) | 1 | 显影、停显 | 出口 | 各档颗粒显影分数平方的个数平均 ⟨x²⟩ |
| `emulsion/silver_physical` | (H, W) | g/m2 | 显影 | 出口 | 溶液物理显影沉积的银 |
| `emulsion/halide` | (H, W) | g/m2 | 显影、定影 | 出口 | 剩余卤化银（以 Ag 计） |
| `emulsion/growth_rate` | (H, W) | s^-1 | 显影 | 停显 | 显影结束时的 κ̃，停显后删除 |
| `frame/fixer_residual` | (H, W) 广播 | 1 | 定影、水洗 | — | 残留定影剂分数（无光学作用） |

显影步骤删除两个 `exposure_*` 字段；停显步骤删除 `growth_rate`。`state.meta` 的新键：`"exposure"`（§2.1.8）、`"develop"`（§2.4.10）、`"post"`（§2.6）。

#### 2.0.4 注册项总表

材料文件的 `model` 为 `bw-mechanistic`。

| 种类 | 名字 | `interaction_range_um` | 模块 |
|------|------|------------------------|------|
| 曝光模型 | `bw-mechanistic` | —（内部规范块，§3.6） | `models/bw.py`（D1a） |
| 步骤 | (`bw-mechanistic`, `presoak`) | 0.0 | `develop/step.py`（D1c） |
| 步骤 | (`bw-mechanistic`, `develop`) | None（药液槽全幅耦合；内部规范块） | `develop/step.py`（D1c） |
| 步骤 | (`bw-mechanistic`, `stop`) | 0.0 | `develop/post.py`（D1d） |
| 步骤 | (`bw-mechanistic`, `fix`) | 0.0 | `develop/post.py`（D1d） |
| 步骤 | (`bw-mechanistic`, `wash`) | 0.0 | `develop/post.py`（D1d） |
| 步骤 | (`bw-mechanistic`, `dry`) | 0.0 | `develop/post.py`（D1d） |
| 出口 | `densitometer_visual` | — | `readout/densitometer.py`（D1d） |
| 出口 | `densitometer_type2` | — | `readout/densitometer.py`（D1d） |
| 出口 | `scan_standard` | — | `readout/scan.py`（D1d） |
| 出口 | `silver_density` | — | `readout/densitometer.py`（D1d，诊断用） |

### 2.1 光传输（阶段 1）

#### 2.1.1 状态与输入

输入是胶片面上的光：每像素每波段的光子通量 `Φ_b(x)`（photons um^-2 s^-1），曝光量 `N_b = Φ_b·t_exp`（photons um^-2）。场景光用低秩表示（§4.2）：

```
Φ_b(x) = Σ_r C_r(x) · Q_r,b ,   r = 1…R（RGB 场景 R = 3）
Q_r,b = k · B_b,r · S_b / (E_b · t_exp)      （D0 的 Mallett 基 B、光源 S、能量尺度 k、光子能量 E_b）
```

输出是每个深度单元的两个吸收描述量 `e_V,z(x)`、`e_S,z(x)`（§2.1.4）。光传输只在曝光期间发生，乳剂是干的（厚度 L_e）。

#### 2.1.2 颗粒的吸收截面

```
σ_c,b = V_c · α_b + S_c · s_b                                  [um²]
α_b = α_0                              当 E_ph,b ≥ E_0
α_b = α_0 · exp((E_ph,b − E_0)/E_U)    当 E_ph,b < E_0          [um^-1]，E_ph,b = 1239.84/λ_b [eV]
s_b = Σ_j A_j · exp(−½ ((ν̃_b − ν̃_j)/w_j)²)                     [1]，ν̃_b = 1000/λ_b [um^-1]
层的吸收系数：μ_a,b = Σ_c (ν_c / L_e) · σ_c,b                   [um^-1]
```

`α` 是 AgX 本征吸收（Urbach 尾，体积项；只在 α·d ≪ 1 时成立，粗颗粒在蓝紫段的自遮蔽忽略，`assumed`，LAT §3.1）；`s` 是增感染料每单位颗粒表面的吸收截面（表面项，LAT §3.2–3.3 的参数化）。两者的形状用胶片自己的光谱感光度曲线标定（§7.3 阶段 B），绝对水平由 D-76 特性曲线和感光度确定。

#### 2.1.3 深度方向：两流（Kubelka–Munk）

乳剂内的光在 ℓ* 量级（2–15 µm，OPT §2.3）内变成漫射光，与层厚相当，用两流模型（`assumed`）。从上表面向下的深度 `s ∈ [0, L_e]`：

```
dI⁺/ds = −(K + S) I⁺ + S I⁻
dI⁻/ds =  (K + S) I⁻ − S I⁺
K_b = 2 μ_a,b,  S_b = S_550 · (550/λ_b)^q_S                    [um^-1]
```

边界条件（每波段，`N` 为入射光子数，`J` 为从片基返回、从下表面进入的漫射光子数）：

```
I⁺(0)   = T_top · N + r_i · I⁻(0)
I⁻(L_e) = R_bot · I⁺(L_e) + J
T_top = 1 − R_top,  R_top = ((n_gel − 1)/(n_gel + 1))²
r_i   = 1 − T_top / n_gel²          （明胶–空气界面对内部漫射光的反射，含全反射）
R_bot = 1 − (n_b / n_gel)²          （乳剂–片基界面对向下漫射光的全反射部分）
```

解析解（β 与 R_∞ 用数值稳定的形式）：

```
β   = sqrt(K (K + 2S))
R_∞ = S / (K + S + β)           （= a − b，a = 1 + K/S，b = sqrt(a² − 1)；向上增长的本征解有 I⁺/I⁻ = R_∞）
I⁺(s) = P e^{−βs} + R_∞ Q e^{−β(L_e − s)}
I⁻(s) = R_∞ P e^{−βs} + Q e^{−β(L_e − s)}
```

`(P, Q)` 由两个边界条件给出 2×2 线性方程。指数项都 ≤ 1，不会溢出；S → 0 时 R_∞ → 0、β → K，上式连续地退化为两束独立衰减的光，不需要特殊处理。

光子通量密度（标量通量）`ψ(s) = 2 (I⁺ + I⁻)`（photons um^-2），单位体积吸收 `μ_a ψ = K (I⁺ + I⁻)`。分别取 `(N, J) = (1, 0)` 与 `(0, 1)` 得到两组响应，在深度单元上积分平均：

```
g⁺_b,z = (1/Δs) ∫_cell ψ ds   （N = 1, J = 0）
g⁻_b,z = (1/Δs) ∫_cell ψ ds   （N = 0, J = 1）
τ_b    = (1 − R_bot) I⁺(L_e)  （N = 1, J = 0）：进入片基的漫射光
```

能量守恒（测试用）：`T_top·N + J = ∫K(I⁺ + I⁻)ds + (1 − r_i) I⁻(0) + (1 − R_bot) I⁺(L_e)`。

#### 2.1.4 吸收描述量

c 档颗粒在单元 z 中吸收的平均光子数是两个描述量的线性组合：

```
m_c,z(x) = V_c · e_V,z(x) + S_c · e_S,z(x)                      [photons per grain]
e_V,z(x) = Σ_b α_b · ψ̄_b,z(x)                                   [um^-3]
e_S,z(x) = Σ_b s_b · ψ̄_b,z(x)                                   [um^-2]
ψ̄_b,z(x) = g⁺_b,z · (P ⊛ N_b)(x) + g⁻_b,z · (H ⊛ (τ_b N_b))(x)
```

`P` 是层内散射 PSF（§2.1.6），`H` 是光晕核（§2.1.5）。两个核都按光谱中性处理（`assumed`：PSF 与深度、波长无关；灰染与背涂层中性），因此先在每个像素上做光谱投影，再做卷积：

```
e_X,z = P ⊛ [Σ_b W⁺_X,z,b N_b] + H ⊛ [Σ_b W⁻_X,z,b N_b]
W⁺_X,z,b = w^X_b g⁺_b,z ,   W⁻_X,z,b = w^X_b g⁻_b,z τ_b ,   w^V = α, w^S = s
```

低秩光可以改为先对 R 个系数通道卷积，再组合：`e_X,z = t_exp Σ_r [(P ⊛ C_r)·(Q W⁺_X,z)_r + (H ⊛ C_r)·(Q W⁻_X,z)_r]`（`Q` 为 R×B，`Q W` 为长 R 的向量）。实现取卷积次数少的一种：低秩光 2R = 6 次，全谱光 2·2n_z 次。

#### 2.1.5 光晕核

片基内的光线按斯涅尔不变量 `n_b sin θ ≥ 1` 判断在片基–空气界面是否全反射，与中间层折射率无关（OPT §3）。向下的朗伯漫射光从一点进入片基，返回乳剂的位置在半径 `r = 2 t_b tan θ` 处。单位入射通量产生的返回辐照度：

```
H(r) = R(θ) · T_g(θ) · cos⁴θ / (4π t_b²)                        [um^-2]
θ = arctan(r / (2 t_b))
T_g(θ) = 10^(−2 (D_g + D_back) / cos θ)
R(θ) = 非偏振菲涅耳反射率（片基 n_b → 空气），θ ≥ θ_c = arcsin(1/n_b) 时为 1
返回分数：f_H = ∫₀^{π/2} 2 sin θ cos θ R(θ) T_g(θ) dθ
```

推导：朗伯源每立体角通量 `cos θ / π`，环带面积 `dA = 4 t_b² tan θ sec²θ dθ dφ`，两者相除得 `cos⁴θ` 形式（与 OPT §3 一致）。只算一次返回（`assumed`；二次返回的权重约为 `f_H · R_∞`，灰片基下约为一次的 7 %，§9）。

离散化：像素积分用 8×8 超采样，跨过 `r_c = 2 t_b tan θ_c` 的像素用 32×32；核截断在累计返回量达到 `(1 − 1e-3) f_H` 的半径 `R_H`，上限 3000 µm，截掉的部分记入 `meta["exposure"]["halation_truncated"]`。数值（§12.7）：TAC n_b = 1.475、t_b = 127 µm、D_g = 0.10 时 θ_c = 42.7°，r_c = 234 µm，f_H = 0.222，`R_H`(99.9 %) ≈ 1.83 mm；透明片基 f_H = 0.581，R_H ≈ 10.5 mm（D1 的胶片都有吸收层，不会用到）。

#### 2.1.6 散射点扩散函数

Frieser 指数线扩散函数对应的二维 PSF（OPT §2.1、§2.3）：

```
P(r) = (1 − w_P) δ(r) + w_P · K₀(r/a_P) / (2π a_P²)             [um^-2]
P̂(ν) = 1 − w_P + w_P / (1 + (2π a_P ν)²)                       （ν 以 cycles/um 计）
```

`K₀` 为第二类修正贝塞尔函数，原点处对数奇异但可积。离散化：每像素 16×16 超采样积分，中心像素用等面积圆盘的解析式 `1 − (ρ/a_P) K₁(ρ/a_P)`；截断半径 `R_P = 10.7 a_P`（尾部质量 < 1e-4，§12.7），截断后归一到 1。PSF 的两个参数通过完整正向模型（含显影）对手册 MTF 标定（§6.7、§7.3），不直接拟合 MTF 曲线。

#### 2.1.7 卷积的执行

卷积在固定的规范块上用 FFT 计算（§3.6）：内块 512 像素，晕边 `h = ceil(max(R_H, R_P)/p)` 像素，FFT 尺寸取 `scipy.fft.next_fast_len(512 + 2h)`，对所有块相同；画幅以外补零（片门之外没有光，这是物理条件，不是近似）。块的结果只取内块写回（§4.4）。

#### 2.1.8 参数

| 符号 | 字段 | 单位 | 类别 | 先验/界 | 出处 |
|------|------|------|------|---------|------|
| n_gel | `emulsion.optics.gelatin_index` | 1 | assumed（unverified） | 1.53–1.54，默认 1.54 | OPT §1.1 |
| n_b | `base.refractive_index` | 1 | literature 1.48 | 1.47–1.48 | OPT §1.1（US 6,398,371 B1） |
| t_b | `base.formats.<f>.thickness` | um | T-Max 135：127，120：119（datasheet）；Tri-X 135：127，120：91（assumed，沿用旧 TX） | — | OPT §1.2 |
| D_g | `base.formats.<f>.grey_density` | 1 | assumed（estimate 的解释未核实；不拟合，§7.1） | 135：0.10（0–0.2）；120：0 | OPT §1.3 |
| D_back | `base.formats.<f>.backing_density` | 1 | assumed | 120：0.5（0.2–1.5）；135：0 | OPT §1.3 |
| α_0 | `emulsion.optics.intrinsic_absorption.alpha0` | um^-1 | fitted | 0.1–10（log 先验） | LAT §3.1 |
| E_0 | `…intrinsic_absorption.edge_energy` | eV | fitted | 2.5–2.9 | LAT §3.1、OPT §1.1（边缘约 2.7 eV） |
| E_U | `…intrinsic_absorption.urbach_energy` | eV | fitted | 0.02–0.15 | LAT §3.1 |
| A_j, ν̃_j, w_j | `emulsion.optics.dye_bands[j]` | 1, um^-1, um^-1 | fitted（J = 2，必要时 3） | A ≤ 0.1；ν̃ 1.4–2.0；w 0.03–0.3 | LAT §3.2–3.3 |
| S_550 | `emulsion.optics.scattering_550` | um^-1 | fitted | 0.05–0.4 | OPT §2.3（S ≈ 0.75/ℓ*） |
| q_S | `emulsion.optics.scattering_exponent` | 1 | assumed 1.0 | 0.5–2 | OPT §2.3 Mie 表 |
| w_P | `emulsion.optics.psf_weight` | 1 | fitted | 0–1 | OPT §2.3 |
| a_P | `emulsion.optics.psf_length` | um | fitted | Tri-X 0.5–5；T-Max 0.3–1.6 | OPT §2.4（尾部拟合给上界） |
| — | 层吸收率 A_550 | 1 | 对 α_0、A_j 的 assumed 约束（§7.1 的单侧罚残差，σ 0.02） | 0.2–0.6 | LAT §2.3（约 6 层颗粒）、§3.2 |

A_550 约束的用处：吸收截面的绝对水平与潜影效率 η 只以乘积 η·σ 进入潜影（泊松稀疏化），二者只能通过深度效应区分；没有这条约束 η 不可辨识（§7.5）。分工：阶段 B 只定 α_0、A_j 的相对形状，绝对尺度由胶片块的 s_abs（§7.1）在阶段 C 定，A_550 罚项把它限制在 0.2–0.6 之内。

#### 2.1.9 状态、注册与作用范围

- 曝光模型 `bw-mechanistic`（`models/bw.py`）：从 `unexposed_state` 出发，写 `emulsion/exposure_volume`、`emulsion/exposure_surface`（形状 `(H, W, n_z)`，`n_z` 取当前离散化档位，§3.1）与广播字段 `emulsion/residual_dye = 1.0`。
- `meta["exposure"]` = `{"time_s", "temperature_c": 23.0, "n_z", "profile", "film_format", "halation_fraction", "halation_radius_um", "halation_truncated", "psf_radius_um", "block_px", "halo_px"}`。23 °C 是 ISO 6 §5.3.1 的曝光环境（CAL §1.1）。
- 曝光不经流水线分块，内部规范块与 `RunContext.tile_px` 无关，结果对 `tile_px` 与进程数逐位不变。

#### 2.1.10 极限情况（D1a 测试）

| 测试 | 条件 | 断言 |
|------|------|------|
| `test_km_energy_balance` | 随机 100 组 (K, S, L_e, r_i, R_bot)，N = 1, J ∈ {0, 1} | 能量守恒相对误差 ≤ 1e-12 |
| `test_km_beer_lambert_limit` | S = 1e-12，r_i = R_bot = 0 | 吸收分数（相对进入的 T_top·N）= 1 − e^{−K L_e}，误差 ≤ 1e-10 |
| `test_km_thick_layer_reflectance` | K = 0.01, S = 1, L_e = 1e4 | 漫反射 → R_∞，误差 ≤ 1e-9 |
| `test_scatter_kernel_normalised` | a_P ∈ {0.3p, p, 4p} | 核之和 = 1 ± 1e-12 |
| `test_scatter_kernel_otf` | 同上，ν ≤ 0.25/p，沿 x 轴 | 离散核的传递函数（DTFT）与像素盒滤波并混叠的 OTF `Σ_{k=−8}^{8} P̂(ν + k/p) · sinc(νp + k)` 之差 ≤ 2e-3（盒的 y 向因子在 k' ≠ 0 处为 0；直接与 P̂(ν) 比较时差 0.002–0.018 × w_P，必然失败） |
| `test_halation_fraction_clear_base` | D_g = D_back = 0 | f_H 与数值积分之差 ≤ 1e-4；全反射部分 = 1 − 1/n_b² ± 1e-4 |
| `test_halation_grey_base_numbers` | n_b = 1.475, t_b = 127, D_g = 0.10 | f_H = 0.222 ± 0.002；r_c = 234 µm ± p |
| `test_low_rank_projection_matches_full` | 64×48 帧，随机 RGB | 低秩与展开为全谱的 `e_V`、`e_S` 相对差 ≤ 1e-12 |
| `test_expose_linear_in_light` | 光加倍 | 字段加倍，相对差 ≤ 1e-12 |
| `test_expose_uniform_interior` | 均匀光 | 距边缘超过晕边的像素相对差 ≤ 1e-12 |
| `test_block_convolution_equals_global_fft` | 300×200 帧，晕边 ≥ 核半径 | 与整幅一次 FFT（同样补零）相对差 ≤ 1e-12 |
| `test_expose_independent_of_tiling_and_workers` | `tile_px` ∈ {None, 37, 256}，`workers` ∈ {1, 3} | 逐位相同 |
| `test_size_classes_conserve_silver` | 随机分布参数 | `1e12 Σ ν_c m_Ag,c` = 涂布量，相对差 ≤ 1e-12 |
| `test_layer_thickness_formula` | 1.75 g/m² Ag + 4.84 g/m² 明胶 | L_e = 4.1 ± 0.05 µm（LAT §2.3） |

### 2.2 潜影（阶段 2）

#### 2.2.1 输入与状态

每个像素、每档 c、每个深度单元 z 的平均吸收光子数 `m_c,z = V_c e_V,z + S_c e_S,z`（§2.1.4），曝光时间 `t_exp`（每帧一个值，恒定照度）。单颗粒的光子吸收率 `Φ_g = m/t_exp`，到达潜影位点的电子率

```
ν = η_c · Φ_g ,   η_c = η_0 · (d_c / 1 um)^ξ_η
```

`η_c` 包含复合损失（LAT §1.2 的 E1–E8 合在一个效率里）。泊松吸收由马尔可夫链本身体现：以恒定速率到达的事件数服从泊松分布，所以链从空状态出发即给出泊松统计（LAT §1.3a）。

#### 2.2.2 两位点阻塞链

每颗粒跟踪两个位点：主位点 1（浓缩原理）与一个汇总的竞争位点 2。状态 `(k₁, k₂)`，`k_s ∈ {0, 1, …, n_cap}` 为该位点的银原子数，`n_cap = 5`（≥5 并为一态）。共 `(n_cap + 1)² = 36` 个状态。

```
q₁ = 1 / (1 + ν τ_i)                       到达主位点的概率（离子步骤占用期间改去别处，M/M/1/1 阻塞）
q₂ = (1 − q₁) / (S_eff − 1)                 到达被跟踪的竞争位点
q_lost = (1 − q₁)(S_eff − 2)/(S_eff − 1)    分散到其余位点、形成会衰减的单原子而丢失

转移（对 s = 1, 2）：
  k_s → k_s + 1   速率 ν q_s          （k_s < n_cap）
  k_s = 1 → 0     速率 1/τ₀ + ρ₁ ν    （单原子热衰减与空穴攻击）
P(t_exp) = expm(Q · t_exp) · e_(0,0)
```

`S_eff = 2` 时没有丢失，退化为纯两位点模型。这个闭合在结构上是 LAT §1.3(c) 的方案 (ii)：`q₁ = 1/(1 + ν τ_i)` 是平均场的稀疏因子（随强度变的效率），离子步骤不是链的状态，被跟踪的竞争位点也没有离子阻塞；它借用了方案 (i) 的多位点图像，但按 LAT 的要求标为估计式闭合（`assumed`）。后果：ν τ_i → ∞ 时各转移率趋于常数比例，高照度失效饱和在 m 的 S_eff − 1 倍（S_eff = 3 时恰为 1 档，`test_hirf_asymptote`），短于约 1e-5 s 的曝光不再继续失效。只用两位点而不加丢失项时，强光下电子几乎全部流向位点 2，分散度反而回落，高强度失效最大只有约 0.45 档，达不到 Tri-X 在 1/100 000 s 的 +1 档；`S_eff = 3`、τ_i 取 1e-5 与 3e-5 s 时，同一结构在 (1/100 000 s, 1/10 000 s) 给出 (+1.14, +0.21) 与 (+1.10, +0.66) 档，表中的 (+1, +½) 落在其间（§12.3）。

#### 2.2.3 陷阱深度分布

单原子的寿命按颗粒（位点）固定，跨颗粒服从对数正态（Katz 1949 要求的准连续陷阱深度分布，LAT §1.2、§1.3b）：

```
τ₀ = τ̄₀ · exp(σ_τ x),   x ~ N(0, 1)   （σ_τ = σ_E / (k_B T_exp)）
p(n | m) = Σ_q ω_q · p(n | m; τ₀(x_q))     概率论意义下的 Gauss–Hermite 节点 x_q、权 ω_q（Σω = 1）
```

节点数 `Q_GH`：标定与图像 9，参考 21。

#### 2.2.4 输出：团簇大小分布

曝光结束后，每个颗粒的潜影大小取最大位点 `n = max(k₁, k₂)`。显影前的存放时间（ISO 6 规定 4 h 至 10 天，CAL §1.1）远长于 τ₀，单原子全部衰减，所以 `n = 1` 并入 `n = 0`（`assumed`；稳定团簇的退化不建模，§9）。输出为每档每单元的分布

```
p_c,z(n),  n ∈ {0, 2, 3, 4, 5}，Σ_n p = 1
```

可显影性不在这里判定：显影步骤用整个分布计算引发风险（§2.4.4），所以不同显影液对"最小可显影团簇"的不同要求由显影动力学产生（LAT §5 第 5 条）。

#### 2.2.5 查找表

每帧的 `t_exp` 固定，每档建一张表：`log10 m` 在 [−6, 4] 上等距 201 点，存 `p(n | m)`；表外下方取第一点（p(n ≥ 2) < 1e-12），上方取最后一点。查表在 `log10 m` 上用三次样条（not-a-knot，C²，建表时对每档每个 n 求好系数）插值，负值截到 0 后重新归一。线性插值在 201 点上误差 2.25e-3，且只有 Lipschitz 连续，参数扰动穿过节点时目标函数出现折点；三次样条在 201 点上误差 9e-6（401 点 5e-7，Poisson r 次命中极限、η = 0.5 实测）。成员的 m 在一次显影中不变，所以每个成员只在聚合 p̄ 与合成输出时各查一次（§2.4.4）。建表只依赖胶片参数与 `t_exp`，与化学无关，进程内按 `(material 哈希, t_exp, profile)` 缓存（CAL §3.3"缓存潜影查找表"）。一张表 = K × 201 × Q_GH 次 36×36 矩阵指数，约 1.1 万次，单核约 0.43 s（§12.8）。

#### 2.2.6 倒易律失效如何产生

- 低照度失效：`ν τ₀` 与 1 可比时，第一个原子在第二个电子到达前衰减。单一 τ₀ 时 1–100 s 段每十倍时间约需 1.3 档，比 Tri-X 表（约 1 档/十倍）陡；σ_τ ≈ 1–3.5 把它压平（§12.2，与 Katz 的结论一致）。失效对较高密度的影响更大，反差上升，与 Kodak 长曝光"减少显影"一致。
- 高照度失效：`ν τ_i ≳ 1` 时电子分散到多个位点，形成多个亚可显影团簇；较高密度需要更多补偿，反差下降，与 Kodak 短曝光"增加显影"一致（§12.3）。
- 两个机理都只依赖 `ν`（即强度）与 `t_exp`。低照度一侧没有经验公式；高照度一侧的 q₁ 是方案 (ii) 型的估计式闭合（§2.2.2），不是由位点状态推出的。

#### 2.2.7 参数

| 符号 | 字段 `emulsion.latent.*` | 单位 | 类别 | 先验/界 | 出处 |
|------|--------------------------|------|------|---------|------|
| η_0 | `efficiency` | 1 | fitted | 0.05–1 | LAT §1.2、§1.4（量子感光度 6–27 光子/颗粒作诊断） |
| ξ_η | `efficiency_size_exponent` | 1 | assumed 0；可辨识时 fitted | −1–1 | LAT §1.3(d) |
| τ̄₀ | `tau0_median` | s | fitted | Tri-X 0.01–2，T-Max 0.05–10（对数先验；中心取 LAT 区间的几何中点 0.22 与 0.89 s） | LAT §1.4 推导（0.05–1 与 0.2–4 s，只用了 1 s 训练行）向下放宽 5 倍、向上放宽 2 倍（`assumed`）；不取 §12.2 的拟合（它用了全部行） |
| σ_τ | `tau0_sigma_ln` | 1 | fitted | 0–4 | LAT §1.2（Katz 要求有分布）；上界：Hailstone & De Keyzer 2004 的边缘陷阱深度 0.2–0.6 eV 作 ±2σ，σ_E ≈ 0.1 eV，除以 23 °C 的 k_B T = 0.0255 eV 得 3.9（LAT 表 E2、§1.3） |
| τ_i | `ionic_time` | s | fitted | 1e-7–1e-4（对数） | LAT §1.4（µs 量级，recalled） |
| S_eff | `competing_sites` | 1 | assumed 3 | 2–10 | LAT §1.2 E10；§12.3 |
| ρ₁ | `hole_attack` | 1 | assumed 0 | 0–1 | LAT §1.3(b) |
| n_cap | `cluster_cap` | 1 | assumed 5 | — | 与 n_sat = 5 一致（KIN §1.1 Ohzeki） |
| T_exp | — | C | assumed 23 | — | CAL §1.1（ISO 6 §5.3.1） |

#### 2.2.8 状态与注册

潜影分布不存为字段：`p_c,z(n)` 是曝光字段与 `t_exp` 的确定函数，显影步骤从 `exposure_*` 字段与查找表现算（存下来 24 MP 需要约 46 GB）。`meta["exposure"]["time_s"]` 是唯一需要传下去的量。D1b 不增加注册项。

#### 2.2.9 极限情况（D1b 测试）

| 测试 | 条件 | 断言 |
|------|------|------|
| `test_poisson_r_hit_limit` | τ̄₀ = ∞，τ_i = 0，S_eff = 2，ρ₁ = 0；m ∈ {0.5, 3, 8, 20}；t ∈ {1e-4, 1, 100} s | P(n ≥ r) = 1 − Σ_{j<r} e^{−ηm}(ηm)^j/j!，r = 2…5，误差 ≤ 1e-10（§12.1） |
| `test_probability_conservation` | 随机参数 | \|Σ_n p − 1\| ≤ 1e-12 |
| `test_monotone_in_exposure` | 查找表网格 | P(n ≥ r) 随 m 不减，r = 2…5 |
| `test_no_reciprocity_failure_without_decay_or_blocking` | τ̄₀ = ∞，τ_i = 0 | p(n \| m, t) 在 1e-5…1e3 s 上不变，误差 ≤ 1e-10 |
| `test_single_site_nucleation_formula` | n_cap = 2，S_eff = 2，τ_i = 0，单一 τ₀ | P(k ≥ 2) 与三态链的双指数闭式解之差 ≤ 1e-10 |
| `test_lirf_raises_contrast` | τ̄₀ = 0.3 s，σ_τ = 2，100 s | 判据 P = 0.4 所需补偿档数大于 P = 0.05 的（§12.2） |
| `test_hirf_asymptote` | τ̄₀ = ∞，S_eff = 3，ν τ_i → 10⁴ | P(n ≥ 2) = 0.5 所需 m 相对互易极限的倍数 = S_eff − 1，即 1 档 ± 0.02 |
| `test_hirf_lowers_contrast` | τ_i = 1e-5 s，S_eff = 3，t = 1e-4 s | 判据 P = 0.4 所需档数大于 P = 0.1 的（§12.3：0.37 对 0.21） |
| `test_gauss_hermite_converged` | σ_τ ≤ 4 | Q_GH = 9 与 21 的 \|Δp\| ≤ 1e-5 |
| `test_latent_table_interpolation` | 网格中点，201 点表 | 与直接计算之差 ≤ 1e-4（三次样条实测 9e-6） |
| `test_latent_table_deterministic` | 相同输入两次 | 数组逐位相同 |

### 2.3 药液化学

#### 2.3.1 配料与组分

两类条目（都在 `data/chemistry/species.yaml`，§5.2）：

- 配料（ingredient）：配方里称量的东西，如 `sodium-sulfite`、`borax-decahydrate`、`metol-sulfate`。每种配料给出溶解后对各输运组分的贡献 `s_i,k`（mol/mol）、对质子总量的贡献 `τ_i`（相对于各组分的参考形态）、阳离子数（Na⁺、K⁺）。
- 组分（component）：溶液中输运的总量，如 `sulfite`（SO₃²⁻ + HSO₃⁻ + H₂SO₃）、`hydroquinone`（H₂Q + HQ⁻ + Q²⁻）。每个组分给出酸碱族（各质子化态的电荷与相对参考态的质子数、pKa）、水中扩散系数、凝胶/水扩散比、活度类别、显影剂角色。

另有两个特殊组分：`toth`（质子过量总量 TOTH，KIN §5.2）与 `gelatin`（只在凝胶中、不输运的缓冲剂）。两者是内置伪组分，不写进 `species.yaml`（D0 的 `Species` 要求 `molar_mass` 与至少一个 `roles`，对它们没有意义）：由 `chemistry/species.py` 的 `component_table` 自动加入，明胶的 g_OH、pK、初始 pH 取自 `kinetics.yaml` 的 `gelatin` 块（§2.4.11），其酸碱态见 §2.3.7。Na⁺、K⁺、硫酸根、二甘醇等只进入离子强度，不输运（`assumed`：其浓度梯度不影响反应）。

质子过量的定义：

```
TOTH = Σ_f Σ_j n_f,j [f, j] + [H⁺] − [OH⁻]
```

`n_f,j` 是族 f 第 j 态相对参考态多出的质子数，求和含凝胶单元里的明胶族。TOTH 是每个单元的守恒库存，但它本身不按菲克律输运：它的通量只由可动形态合成（§2.4.8），明胶结合的质子 `T_gel n̄_gel` 留在原单元。可动质子量记为

```
TOTH_mob = Σ_{f 可动} T_f n̄_f(pH) + [H⁺] − [OH⁻] = TOTH − T_gel n̄_gel(pH)
```

参考态与 n 的取值见 §2.3.7 的表。配料对 TOTH 的贡献等于它溶解后各态的 n 之和，例如 Na₂S₂O₅ → 2 HSO₃⁻，τ = +2；硼砂 → 2 B(OH)₃ + 2 B(OH)₄⁻，τ = −2；米吐尔硫酸盐 → 2 MetH₂⁺，τ = +2（硫酸根作旁观离子）。

#### 2.3.2 由配方得到工作液

药液文件给出"基液"（D-76、XTOL 为原液 stock；HC-110、T-MAX、Rodinal 为浓缩液）中每种配料的 g/L。流程步骤的稀释写成相对基液的 `a+b`，稀释分数 `f = a/(a+b)`（D0 `parse_dilution`）。工作液：

```
T_k   = f · Σ_i (w_i / M_i) · s_i,k        [mol/L]
TOTH  = f · Σ_i (w_i / M_i) · τ_i          [mol/L]
m_cat = f · Σ_i (w_i / M_i) · (n_Na,i + n_K,i)
```

`w_i` 为配料 g/L，`M_i` 为摩尔质量。离子强度在形态分布之后算 `I = ½ Σ_s z_s² [s]`（含旁观离子），从 `I⁰ = m_cat` 起做定点迭代，直到 |ΔI| < 1e-6 M（3–5 次）。显影过程中 I 视为常数（`assumed`：显影释放的 Br⁻ 与 H⁺ 在 10 mM 量级，占 D-76 的 I ≈ 2.4 M 不到 1 %）。

干扰参数（如成品 D-76 的 B₂O₃）就是某个配料的 `w_i`，标定时在界内拟合（§5.5）。

#### 2.3.3 活度：SIT 型修正

离子强度约 2 M 时 Davies 式失效（DEV §5.1）。采用 SIT（specific ion interaction）形式，pH 取氢离子活度（玻璃电极惯例，`assumed`）：

```
log10 γ_i = −z_i² · D(I) + ε_i · m_cat
D(I) = A(T) · √I / (1 + 1.5 √I)
```

对酸碱族第 j 步 `HA^z ⇌ A^(z−1) + H⁺`，混合（条件）常数把浓度比与 a_H 联系起来：

```
pK'_f,j(T, I) = pK°_f,j(T) + (2 z_HA − 1) · D(I) + Δε_f · m_cat
pH = pK'_f,j + log10([A]/[HA])
```

`(2z_HA − 1)` 来自 `log γ_A − log γ_HA`：HSO₃⁻/SO₃²⁻（z_HA = −1）为 −3D，B(OH)₃/B(OH)₄⁻ 为 −D，DEAH⁺/DEA 为 +D。水：`[H⁺] = a_H / γ_H`，`[OH⁻] = K_w/(a_H γ_OH)`，两者 ε 取 0。

- `A(T)` 在 0.5050（20 °C）与 0.5091（25 °C）之间对 T 线性插值、外推；1.5 是 SIT 的惯用常数。两者都是 recalled 的标准值（`assumed`，待核对）。
- `Δε_f` 按族归为四类，标定阶段 A 只用 pH 锚点拟合（§7.3）：`sulfite`、`borate`、`organic`（所有酚、烯醇、螯合剂的酸步）、`amine`。
- 只保留 D(I) 项、Δε = 0 时，D-76 公开配方算得 pH 8.44（实测 8.5），XTOL 例 1 为 8.29（实测 8.20）；理想溶液分别为 8.95、8.89（§12.4）。

温度：

```
van 't Hoff：pK°(T) = pK°(T_r) + ΔH°/(R ln 10) · (1/T − 1/T_r)
硼酸（Manov 等 1944）：pK° = 2237.94/T + 0.016883 T − 3.305
胺：pK°(T) = pK°(25 °C) + (dpK/dT)(T − 298.15)
pK_w：过 14.17（20 °C）与 14.00（25 °C）的 van 't Hoff 式，ΔH_w = 56.9 kJ/mol（由两值导出）
```

有机显影剂的电离焓取 0（`assumed`）：它们的温度效应并入显影的表观活化能（§2.4.5），以免与 Sugai 的表观值重复计算。

#### 2.3.4 形态分布

给定单元内各族总量 `T_f` 与 `TOTH`，解一个标量方程求 pH：

```
F(pH) = Σ_f T_f · n̄_f(pH) + [H⁺] − [OH⁻] − TOTH = 0
n̄_f(pH) = Σ_j n_f,j · α_f,j(pH)，α 由 pK' 给出
```

F 对 pH 单调递减，`−∂F/∂pH = β`（缓冲容量，≥ 0）。牛顿迭代带区间保护（pH ∈ [0, 14.5]），容差 |ΔpH| ≤ 1e-12，上一次的 pH 作初值，最多 30 次，对所有单元向量化。导数：

```
∂pH/∂TOTH = −1/β ,   ∂pH/∂T_f = n̄_f / β
```

明胶族只在凝胶单元出现，`T_gel = c_gel · g_OH`（§2.4.7）。

#### 2.3.5 显影剂的活性形态与种类

| 组分 | 角色 | 活性形态（电荷 z） | 氧化产物 | 产物是否被亚硫酸盐磺化 | 出处 |
|------|------|-------------------|----------|------------------------|------|
| `hydroquinone` | 再生剂 | HQ²⁻（2） | 醌 | 是 | KIN §1.2（Eskes）、§2.2 |
| `ascorbate`（异抗坏血酸） | 再生剂 | HA⁻（1） | 脱氢抗坏血酸 | 否 | KIN §1.2、§4.1 |
| `metol` | 电子传递剂（ETA） | 酚氧负离子（1），`assumed` | 醌单亚胺 | 是（`assumed`） | KIN §1.2（recalled） |
| `dimezone-s`、`phenidone`、`methylphenidone` | ETA，共用 `pyrazolidone` 动力学类 | 阴离子（1） | 羟基吡唑 | 否 | KIN §1.2、§4.1 |
| `p-aminophenol` | 单剂（无再生剂） | 酚氧负离子（1） | 醌亚胺 | 是（Rodinal 亚硫酸盐少） | KIN §1.2、§4.1 |

#### 2.3.6 参数

| 符号 | 字段 | 单位 | 类别 | 先验/界 | 出处 |
|------|------|------|------|---------|------|
| A(T) | 常数 | (kg/mol)^½ | assumed 0.5050/0.5091 | — | recalled 标准值 |
| Δε_sulfite, Δε_borate, Δε_organic, Δε_amine | `kinetics.yaml: activity.delta_eps.*` | kg/mol | fitted（amine 不可辨识时 assumed 0） | N(0, 0.15²)（assumed） | DEV §5.1 建议；§7.3 阶段 A |
| pK_w | 常数 | 1 | assumed（recalled） | 14.17/14.00 | DEV §6.2 |

各组分的酸碱与输运常数在 §2.3.7。

#### 2.3.7 组分数据表

pK° 为 25 °C、I → 0（或报告所给条件）；D_w 为 20 °C 水中扩散系数；r_gel = D_gel/D_w。

| 组分 | 各态（参考态加粗）与 n | pK° | 温度关系 | D_w（um²/s） | r_gel | 出处与类别 |
|------|------------------------|-----|----------|--------------|-------|------------|
| `sulfite` | H₂SO₃ (+2)，HSO₃⁻ (+1)，**SO₃²⁻** (0) | 1.86，7.19 | ΔH₂ = −3.65 kJ/mol | 1030 | 0.5 | pK、ΔH：DEV §6.2（Goldberg & Parker 1985，literature）；D：由 recalled 离子扩散系数按 Nernst–Hartley 算出（Na₂SO₃ 是主体电解质，盐值适用，`assumed`）；r_gel：注 † |
| `borate` | **B(OH)₃** (0)，B(OH)₄⁻ (−1) | Manov 式 | Manov 式 | 1030 | 0.5 | pK：DEV §6.2（literature）；D：未找到，取同 sulfite（`assumed`）；多硼酸根忽略（`assumed`，DEV §6.2） |
| `hydroquinone` | **H₂Q** (0)，HQ⁻ (−1)，Q²⁻ (−2) | 9.85，11.45 | ΔH = 0 | 855 | 0.10 | pK：DEV §6.1（IUPAC，literature；pKa₂ 取 11.39–11.56 的中点）；D：DEV §6.3（estimate → `assumed`）；r_gel：KIN §3.2（Iwano，literature） |
| `metol` | MetH₂⁺ (+1)，**MetH** (0)，Met⁻ (−1) | 5.9，10.25 | ΔH = 0 | 740 | 0.13 | pKaH：DEV §6.1（Fieser，literature）；酚 pKa：DEV §6.1 estimate 10.1–10.4（`assumed`）；D：DEV §6.3；r_gel：Iwano |
| `p-aminophenol` | H₃N⁺ (+1)，**中性** (0)，酚氧负离子 (−1) | 5.48，10.46 | ΔH = 0 | 825 | 0.12 | pK：DEV §6.1（literature）；r_gel：`assumed`（与 HQ、米吐尔同量级） |
| `dimezone-s` | **PH** (0)，P⁻ (−1) | fitted（7.5–9.5） | ΔH = 0 | 555 | 0.2（0.1–0.3） | pKa 冲突 9.48/7.57 未解决（DEV §6.1、§7 第 12 条），按 DEV 建议 fitted；D：DEV §6.3；r_gel：KIN §3.2 先验（`assumed`） |
| `phenidone` | 同上 | 共用 `pyrazolidone` 类 pKa | — | 650 | 0.2 | 同上 |
| `methylphenidone` | 同上 | 共用 | — | 605 | 0.2 | 同上 |
| `ascorbate` | **H₂A** (0)，HA⁻ (−1)，A²⁻ (−2) | 4.04，11.34 | ΔH = 0 | 710 | 0.2（0.1–0.3） | pK：DEV §6.1（Khan & Martell，literature；异抗坏血酸按同值，`assumed`）；D：DEV §6.3 |
| `diethanolamine` | DEAH⁺ (+1)，**DEA** (0) | 8.883 | −0.025 /K | 785 | 0.4（0.2–0.8） | pK：DEV §6.2（literature）；r_gel：`assumed`（参照 Iwano 的 HCHO） |
| `ethanolamine` | MEAH⁺ (+1)，**MEA** (0) | 9.498 | −0.030 /K | 1090 | 0.4 | 同上 |
| `methylethanolamine` | 同上 | 9.82 | −0.028 /K（`assumed`） | 1090（取 MEA，`assumed`） | 0.4 | pK：TMX §2.4（"Approximate"，literature） |
| `dtpa` | **H₅L** (0) … L⁵⁻ (−5) | 1.8，2.6，4.3，8.6，10.5 | ΔH = 0 | 555（`assumed`） | 0.2 | recalled（DEV §6.2）→ `assumed` |
| `hedta` | **H₃L** (0) … L³⁻ (−3) | 2.6，5.4，9.8 | ΔH = 0 | 555（`assumed`） | 0.2 | recalled → `assumed` |
| `bromide` | Br⁻ | — | — | 1816 | 0.5 | Br⁻ 在支持电解质中是痕量离子，取示踪扩散系数 2.080e-5 cm²/s（25 °C，DEV §6.3，recalled）按下式换到 20 °C（`assumed`），不用 NaBr 盐值 1420；r_gel：见注 † |
| `iodide` | I⁻ | — | — | 1786 | 0.5 × 0.94 | D：示踪值 2.045e-5 cm²/s（25 °C，recalled → `assumed`）换到 20 °C；r_gel 再乘 Iwano KI/KBr 比 0.94（KIN §3.2，literature） |
| `silver` | 溶解的 Ag(I)–亚硫酸根配合物 | — | — | 1030（同 sulfite，`assumed`） | 0.5 | `assumed`（注 †） |
| `toth` | — | — | — | 不单独给 D：通量由各可动族的质子量与 H⁺、OH⁻ 合成（§2.4.8） | — | KIN §5.2"明胶进入 TOTH 但不输运" |
| H⁺、OH⁻（水的离子，不作组分） | — | — | — | 8130、4600 | 0.5 | 示踪值 9.31e-5、5.27e-5 cm²/s（25 °C，recalled → `assumed`）换到 20 °C；字段 `transport.proton_diffusion`、`transport.hydroxide_diffusion`；r_gel 取盐类值（`assumed`） |
| `gelatin` | G–NH₃⁺ (+1)，**G–NH₂** (0) | 10.5（9.5–11.5） | ΔH = 0 | 0（不动） | — | 基团含量 4.8e-4 mol/g：KIN §3.2（Iwano，literature）；pKa：`assumed` |

注 †：盐类（sulfite、borate、bromide、iodide、silver、H⁺/OH⁻）的 r_gel 取 0.5（`assumed`，界 0.3–0.8）。Iwano 的 KBr 比值在 8 % 明胶中约 0.8、20 % 中约 0.55（KIN §3.2，literature），而 α_sw = 3 的溶胀乳剂含明胶约 30–45 wt %，所以 0.8 不适用，0.5 是向高浓度的外推。r_gel 不随 α_sw 变（`assumed`；凝胶内扩散不限速，§2.4.12）。HQ、米吐尔的比值是 Iwano 在 α = 3.03 的涂层中测的，保持 literature。

扩散系数的温度与黏度换算（Stokes–Einstein，KIN §3.1；DEV §6.3）：

```
D(T, I) = D_w(T_r) · (T / T_r) · (η_w(T_r) / η_w(T)) / (1 + b_η I),   T_r = 293.15 K
ln η_w 对 1/T 线性插值，节点：1.053 cP（18 °C）、1.002（20）、0.911（24）、0.890（25）
b_η = 0.125 L/mol
```

η_w 的四个值是 KIN §1.6、§3.1 引用的 CRC 值（recalled → `assumed`）；b_η 由 DEV §6.3 的估计"0.7–0.8 M 亚硫酸盐中约 1.2–1.4 cP"在 I ≈ 2.4 M 处取中值 1.3 cP 得到（`assumed`）。

#### 2.3.8 极限情况（D1c 测试）

| 测试 | 断言 |
|------|------|
| `test_speciation_ideal_limit` | I → 0、Δε = 0 时，单一弱酸–共轭碱混合物的 pH 与 Henderson–Hasselbalch 解析式之差 ≤ 1e-10 |
| `test_buffer_titration` | TOTH 增加 δ 后的 pH 变化与 ∫ dTOTH/β 的数值积分之差 ≤ 1e-9 |
| `test_toth_round_trip` | 由形态分布重算的 TOTH 与输入的相对差 ≤ 1e-13 |
| `test_d76_published_ideal_ph` | 理想溶液（D = 0、Δε = 0）pH = 8.95 ± 0.02（DEV §5.1 给 8.94，§12.4） |
| `test_xtol_ideal_ph` | 理想溶液：工作液 8.89 ± 0.02，A 部分（850 mL）11.04 ± 0.03（DEV §5.1：8.88、11.03） |
| `test_sit_debye_huckel_only` | Δε = 0、含 D(I)：D-76 公开配方 8.44 ± 0.02，XTOL 工作液 8.29 ± 0.02（§12.4；实现的回归测试，不是验证） |
| `test_ionic_strength_d76` | D-76 公开配方 I = 2.39 ± 0.03 M |
| `test_boric_acid_pka_temperature` | 20 °C 9.278 ± 0.001，25 °C 9.235 ± 0.001 |
| `test_dilution_raises_ph_in_sulfite_buffer` | 含 D(I) 时 XTOL 1+1 的 pH 高于原液（DEV §5.2 的方向） |
| `test_working_solution_conserves_mass` | 各配料 g/L 换算后，各元素（S、B、Br、显影剂）总量与配方相等，相对差 ≤ 1e-13 |

### 2.4 显影（阶段 3）

#### 2.4.1 几何与状态

一根"柱"（cell）自下而上由以下几层组成，每层是若干有限体积单元（level）：

| 层 | 单元 | 厚度 | 说明 |
|----|------|------|------|
| 乳剂 | z = n_z … 1（自下而上） | Δ_z = α_sw L_e / n_z | 溶胀后；颗粒与明胶在此 |
| 保护层 | 1 | α_sw L_oc | 明胶，无颗粒 |
| 液体边界层 | l = 1 … n_l，分两段：下段 n_a = n_l/2 个在 [0, δ_agit]，上段 n_l − n_a 个在 [δ_agit, δ_max] | 每段内几何递增，比 r_l；段 [a, b] 内 n 个单元的首层 h₁ = (b − a)(r_l − 1)/(r_l^n − 1) | 从胶片表面到 δ_max；h = δ_agit 恰好是一个单元界面，所以所有单元厚度都是 δ_agit、δ_max 的光滑函数 |
| 槽液 | — | — | 充分混合，作为顶端狄利克雷值 |

柱的状态：各组分 k 在各单元的浓度 `C_k,ℓ`（mol/L）；每个乳剂单元的引发时钟 `Λ_z`（1）与溶解计数 `D_z`（mol/L，卤化银溶解量）。

一根柱可以带若干"成员"（member）：感光条带与 MTF 条带里每柱一个成员；图像路径里化学单元边长固定为 Δx_c = 24 µm（与像素间距无关），每个像素是一个成员，按与单元的面积重叠分配（§3.8）。成员的量：

| 量 | 含义 | 每成员个数 |
|------|------|-----------|
| `X_c,z` | c 档在单元 z 中已显影银占该档银量的分数 | K·n_z |
| `Y_c,z` | 同一群颗粒显影分数平方的个数平均 ⟨x²⟩ | K·n_z |
| `P_z` | 物理显影沉积的银占单元总银量的分数 | n_z |

成员的输入是潜影分布 `p_c,z(n)`，由它自己的 `m_c,z` 查表得到（§2.2.5）。X、Y、P 对 p 是线性的，所以成员不进入时间积分：柱只积分按团簇大小 n 分解的核 `G、H、Q`（每单元 3·K·n_z·N_n 个，N_n = 5），成员的 X、Y、P 在显影结束与各快照时刻由核合成（§2.4.4）。

#### 2.4.2 引发（诱导期）

一颗团簇大小为 n 的颗粒开始显影的风险率可分离为化学因子与团簇因子（Sugai 1953 的 τ_ind = C·exp((U_charge − U_LI)/kT) 正是这种乘积形式，KIN §1.4、§1.10）：

```
λ_n(t) = h_n · λ̃_z(t)                                                      [s^-1]
λ̃_z = Σ_i λ0_i · A_i(T) · θ_i · exp(−|z_i| · (u_gel √(I_ref/I) + u_X θ_X))
h_0 = h_fog
h_n = exp(−u_LI (n_sat − n)) · exp(−ζ_s [SO₃²⁻]_bulk / (1 mol/L))          （2 ≤ n < n_sat）
h_n = 1                                                                     （n ≥ n_sat）
Λ_z(t) = ∫₀ᵗ λ̃_z dt'
已引发分数：F_c,z(t) = Σ_n p_c,z(n) · (1 − exp(−h_n Λ_z(t)))
```

- 第一个指数是电荷势垒：明胶势垒随离子强度屏蔽（James 1944，KIN §1.1；`√(I_ref/I)` 形式为 `assumed`，I_ref = 1 mol/L），吸附卤化物加高势垒（KIN §1.4、§1.5）。
- `h_n` 的阶梯：≥ 5 个原子的中心显影快、< 5 个的慢（Ohzeki 1996，KIN §1.1），所以 n_sat = 5（`literature`）。
- 亚硫酸盐的溶解作用使小潜影中心更难引发（KIN §2.2），用槽液初始游离亚硫酸根表示（`assumed` 形式），使 h_n 在一次显影中不随时间变。
- 雾翳是 n = 0 的颗粒以 h_fog 引发（KIN §1.8）。

`F_c,z` 只通过 `Λ_z` 依赖时间，且对 n 可分离，所以直接精确求值，不建表：每个乳剂单元每次 f 求值算 N_n = 5 个 `E_n,z = 1 − exp(−h_n Λ_z)` 与 `∂E_n,z/∂Λ = h_n exp(−h_n Λ_z)`，所有档共用；`F_c,z = Σ_n p_c,z(n) E_n,z`。F 对 Λ 光滑（旧版的 201 × 121 双线性 F 表误差 3.4e-3、且只有 Lipschitz 连续，已删除）。

#### 2.4.3 吸附：竞争 Langmuir–Hinshelwood 与 ETA 的表面周转

```
a_i = α_i,act(pH) · C_i                         活性形态浓度 [mol/L]
Ξ = 1 + Σ_E K_E a_E (1 + Φ_E) + Σ_R K_R a_R + K_Br [Br⁻] + K_I [I⁻]
θ_i = K_i a_i / Ξ ,   θ_X = (K_Br [Br⁻] + K_I [I⁻]) / Ξ
Φ_E = k_et,E A_E / (k_des,E A_E + Σ_R k_reg,R A_R a_R)
g_E = Σ_R k_reg,R A_R a_R / (k_des,E A_E + Σ_R k_reg,R A_R a_R)
```

E 取电子传递剂（米吐尔、吡唑烷酮类、对氨基酚），R 取再生剂（对苯二酚、异抗坏血酸）。被氧化的 ETA 留在表面占位，要么由溶液中的再生剂就地还原（概率 g_E，ETA 不消耗），要么脱附而损失（`assumed`：脱附后的自由基在溶液中的再生忽略）。占位使 ETA 单独使用时周转受脱附速率限制；有再生剂时位点很快腾出，所以合用的速率大于两者单用之和。超加和因子是输出，不是参数（KIN §2.1 的再生机理）。

这一形式在中等覆盖率下给出分数级数 `∂ln κ̃/∂ln C = 1 − θ`（θ ≈ 0.4 时为 James 1946 的 0.6，KIN §1.3），并给出"pH 越高溴化物影响越小"（活性形态浓度高时压过 Br⁻，KIN §1.5）。

#### 2.4.4 生长与成员方程

已引发颗粒按一阶动力学显影：`dx/dt = κ_c (1 − x)`。对一档的颗粒群体求和得到封闭方程（新引发的颗粒 x = 0，对 X、Y 无初始贡献）：

```
κ̃_z = Σ_i k_et,i A_i(T) θ_i                                              [s^-1]
κ_c,z = (ℓ_ref S_c / V_c) · κ̃_z ,   ℓ_ref = 1 um                         （显影速率按比表面积缩放，assumed）
dX_c,z/dt = κ_c,z (F_c,z − X_c,z)
dY_c,z/dt = 2 κ_c,z (X_c,z − Y_c,z)
dP_z/dt   = k_pd · a_dev,z · C_Ag(I),z / C_tot,z
a_dev,z   = Σ_c n_c,z S_c X_c,z              （已显影银的表面，um^-1，assumed 与显影分数成正比）
C_tot,z   = Σ_c C_c,z
```

核分解（精确，不是近似）：给定单元量 κ_c,z(t)、Λ_z(t)、C_Ag(I),z(t)，上面三式对 p_c,z(n) 线性，所以

```
X_c,z = Σ_n p_c,z(n) G_c,z,n ,   Y_c,z = Σ_n p_c,z(n) H_c,z,n ,   P_z = Σ_c Σ_n p_c,z(n) Q_c,z,n
dG_c,z,n/dt = κ_c,z (E_n,z − G_c,z,n)
dH_c,z,n/dt = 2 κ_c,z (G_c,z,n − H_c,z,n)
dQ_c,z,n/dt = k_pd · n_c,z S_c G_c,z,n · C_Ag(I),z / C_tot,z
G = H = Q = 0（t = 0）
```

柱的化学只需要成员的面积权平均 `⟨X_c,z⟩ = Σ_n p̄_c,z(n) G_c,z,n`，其中 `p̄_c,z(n) = Σ_m w_m p_m,c,z(n)` 是单元内成员潜影分布的面积权平均（每个单元在显影开始时算一次）。所以一个单元无论有多少成员，ROS2 只积分它自己的核（3·K·n_z·N_n 个状态）。每柱一个成员时 p̄ = p，与直接积分 X、Y、P 在舍入误差内相同（`test_kernel_reduction_exact`）。

生长律在文献里没有定论（面积律自催化与线性丝状生长之争，KIN §1.7）。一阶律在群体层面封闭、每档只要两个状态；引发风险的分布已经给出单颗粒的 S 形时间曲线。这是一个模型结构选择（`assumed`），其影响用参考解里的"零级生长"变体检查（§9）。

#### 2.4.5 温度

```
A_i(T) = exp(−(Ea_i / R) · (1/T − 1/T_r)),   T_r = 293.15 K
```

同一种显影剂的 k_et、λ0、k_des、k_reg 共用 A_i（`assumed`）。米吐尔与对苯二酚的 Ea 用 Sugai & Furuichi 1953 的值作先验（0.64 与 0.70 eV，即 61.8 与 67.5 kJ/mol，±10 kJ/mol，KIN §1.6、§7 第 2 条）。吡唑烷酮类、异抗坏血酸、对氨基酚的 Ea 没有文献值，且 Kodak 的 XTOL 时间–温度行与 D-76 共用一个斜率（KIN §1.6），不能用来拟合，固定为 65 kJ/mol（`assumed`）。电离焓取 0（§2.3.3），所以这里的 Ea 就是显影剂在固定配方下的表观值。卤化物吸附常数按 van 't Hoff 随温度变（ΔH_Br，§2.4.11），所以溴化物多时表观活化能升高（Sugai & Furuichi：加 KBr 后 0.70 对 0.64 eV，KIN §1.6）可以由模型产生；其余吸附常数与势垒 u 不随温度变（`assumed`）。

#### 2.4.6 反应速率与化学计量

每个乳剂单元（mol L⁻¹ s⁻¹；`⟨·⟩` 表示对柱内成员按面积权平均）：

```
R_dev  = Σ_c C_c,z ⟨κ_c,z (F_c,z − X_c,z)⟩                化学显影
a_undev = Σ_c n_c,z S_c ⟨1 − X_c,z⟩                       未显影卤化银表面 [um^-1]
C_sat  = C_sat,ref · ([SO₃²⁻]⁺/0.1 M)^ν_s · (1 mM / √(([Br⁻]⁺)² + (1e-5 M)²))^ν_b
R_diss = k_diss · a_undev · max(C_sat − C_Ag(I), 0)          溶解
R_pd   = k_pd · ⟨a_dev⟩ · C_Ag(I)                           溶液物理显影
R_e    = R_dev + R_pd                                         电子需求
r_i    = f_i R_e ,  f_i = k_et,i A_i θ_i / Σ_j k_et,j A_j θ_j   各显影剂分担
Q̇      = ½ (Σ_E g_E r_E + Σ_R r_R)                            再生剂被氧化（2 个电子一个）
Ė_ox   = ½ (1 − g_E) r_E                                      ETA 损失
R_sulf = σ_R Q̇ + σ_E Ė_ox,  σ = (产物可磺化) · [SO₃²⁻]/([SO₃²⁻] + c_s)
```

`k_diss`、`k_pd` 的单位是 um/s，乘以 um^-1 的表面积密度得到 s^-1，再乘以 mol/L，量纲一致。`c_s = 1e-4 mol/L`（`assumed`）使低亚硫酸盐显影液（Rodinal）平滑地退出磺化。

光滑正部：所有速率式（本节、§2.4.2、§2.4.3、形态分布的输入）中的浓度一律取 `C⁺ = ½(C + √(C² + ε_c²))`，ε_c = 1e-9 mol/L（`assumed`）；输运、库存与守恒账用 C 本身。ROS2 不保正（k₂ 含 −2k₁），干膜入液、无溴显影液开始时会出现舍入级的负值，C⁺ 使速率对它们光滑且有定义（C⁺ − C ≤ ε_c²/(4C)，C ≫ ε_c 时可忽略）。

组分源项：

| 组分 | 源项（每单位体积每秒） | 依据 |
|------|------------------------|------|
| ETA（每种） | −Ė_ox | 每个电子半个 ETA |
| 再生剂 | −Q̇ | KIN §4.1 |
| `sulfite` | −R_sulf − n_lig R_diss + n_lig R_pd | 磺化、络合溶解、物理显影放出配体 |
| `bromide` | +(1 − x_I) R_dev + R_diss | KIN §4.1 |
| `iodide` | +x_I R_dev | KIN §4.1 |
| `silver` | +R_diss − R_pd | KIN §2.2 |
| `toth` | +R_e − R_sulf | 每个电子放出 1 个 H⁺（相对参考态）；每次磺化放出 1 个 OH⁻（KIN §4.1，HQ 加磺化合计 0.5 H⁺/Ag） |
| `D_z`（计数） | +R_diss | 银守恒的记账 |

溶解式中的 max 在 Jacobian 中用 softplus 平滑（宽度 1e-3 C_sat，`assumed`）。脱氢抗坏血酸的水解（KIN §4.1，最终 1.5 H⁺/Ag）不建模（`assumed`），溶解对颗粒剩余卤化银的削减也不反馈到显影（`assumed`，只计入守恒账）。

#### 2.4.7 明胶缓冲与溶胀

```
c_gel = 1e3 · m_gel / (α_sw L_e)        [g/L]（保护层用 ρ_gel L_oc 代替 m_gel，厚度 α_sw L_oc）
T_gel = c_gel · g_OH                    [mol/L]，g_OH = 4.8e-4 mol/g
```

明胶族只进入凝胶单元的 TOTH 平衡，不输运（KIN §3.2、§5.2）：它的结合质子 T_gel n̄_gel(pH) 计入单元的 TOTH 库存，但不进入 TOTH 的通量（§2.4.8），所以它作为局部缓冲池使 pH 前沿以 D/(1 + β_gel/β_mob) 推进（KIN §3.2 Iwano 的 OH⁻ 观察）。干膜入液时凝胶里只有明胶，pH = pH_gel,0；显影液进入后中和明胶并升 pH，这是起始瞬态。溶胀用扩散系数的爬升表示（KIN §3.2 的建议：Fickian 输运加起始渗透延迟）：

```
D_k,gel(t) = r_k · D_k(T, I) · s(t),   s(t) = 1 − exp(−t ln2 / T_½)   （预湿后 s ≡ 1）
```

几何从 t = 0 起就取溶胀后的厚度（`assumed`）。

#### 2.4.8 输运、边界与搅动

竖直方向（每组分、每单元，有限体积，界面用厚度加权的调和平均 D）：

```
dC_k,ℓ/dt = (J_k,ℓ−½ − J_k,ℓ+½) / Δ_ℓ + S_k,ℓ
J_k,ℓ+½ = −D̄_k,ℓ+½ (C_k,ℓ+1 − C_k,ℓ) / (½(Δ_ℓ + Δ_ℓ+1))     [mol L^-1 um s^-1]
D̄_ℓ+½ = (Δ_ℓ + Δ_ℓ+1) / (Δ_ℓ/D_ℓ + Δ_ℓ+1/D_ℓ+1)
```

上式用于除 `toth` 以外的组分（`gelatin` 不动，没有通量）。TOTH 的通量由可动形态合成，明胶结合的质子不跨面（§2.3.1）：

```
P_f,ℓ = T_f,ℓ · n̄_f(pH_ℓ)                                         族 f 的质子量（相对参考态）
J_toth,ℓ+½ = Σ_{f 可动} −D̄_f,ℓ+½ (P_f,ℓ+1 − P_f,ℓ) / h_ℓ+½
           − D̄_H,ℓ+½ ([H⁺]_ℓ+1 − [H⁺]_ℓ) / h_ℓ+½ + D̄_OH,ℓ+½ ([OH⁻]_ℓ+1 − [OH⁻]_ℓ) / h_ℓ+½
h_ℓ+½ = ½(Δ_ℓ + Δ_ℓ+1)
```

`D̄_f,ℓ+½` 与该族总量 T_f 用的是同一个面上的 D（凝胶中含 r_f 与 s(t)），否则结果会依赖参考态的选取。pH_ℓ 来自该单元的形态分布（§2.3.4）；顶端狄利克雷节点用槽液的形态分布。均匀 pH、均匀总量时 J_toth = 0，与明胶含量无关。

- 底面（片基）零通量；凝胶与液体之间分配系数取 1（`assumed`）。
- 液体单元的 D 不乘 r_k 与 s(t)。
- 横向：每个单元层另有 `D_k,ℓ ∇²_xy C_k,ℓ`（凝胶与液体都有），在条带中没有这一项；TOTH 的横向项作用在可动质子量 TOTH_mob 上（§3.5）。
- 顶端：静置期在 h = δ_max 处取槽液值（狄利克雷，自然对流把扩散层限制在这个厚度，KIN §3.4）；搅动期狄利克雷面移到 h = δ_agit（§2.4.1 保证它是单元界面），下段顶层的通量用半单元距离 `J = −D (C_bulk − C_ℓ) / (½Δ_ℓ)`，上段各单元在搅动期内取槽液值 `C_bulk(t)`（不参与求解，库存计入槽而不计入液层），搅动结束时以槽液值为初值重新参与求解。搅动开始时上段混入槽：它与槽液的库存差 `ΔI_k = Σ_{ℓ∈上段} Δ_ℓ (C_k,ℓ − C_k,bulk)`（mol L⁻¹ µm，静置期积累的 Br⁻ 富余与显影剂亏损）作为一次性的混合通量计入槽（§2.4.9），否则每次搅动丢掉约三分之一个周期的溴化物释放（§2.4.12 的剖面）。于是各档位看到的搅动层厚度都恰好是 δ_agit，模型对 δ_agit、δ_max 连续可微（`test_strip_smooth_in_transport`，§2.4.13）。
- 搅动期由流程的 `agitation` 字段给出（§5.6）：

| `method` | 搅动期 |
|----------|--------|
| `inversion`（小罐） | [0, initial_s]，以及每个 interval_s 的整数倍开始、持续 duration_s 的时段（Kodak：5 s，每 30 s 一次 5 s，DS §2.2、KIN §3.4） |
| `lift`（大罐，新增） | [0, initial_s]（Kodak 前 15–30 s 连续，取 22.5 s），之后每 interval_s（60 s）一次、持续 duration_s（提起、倾斜、沥 5–10 s，取 7.5 s；沥干期间离开药液这一点忽略，`assumed`） |
| `rotary`、`continuous`、`tray` | 全程搅动 |
| `stand` | 只有 [0, initial_s] |

搅动的效果因此只有两个共享参数（δ_agit、δ_max），两种 Kodak 罐型只在间隔和时长上不同。这也处理了 Tri-X 的 D-76 曲线（大罐）与 T-Max 100 的 D-76 曲线（小罐）搅动不一致的问题（DS §0、KIN §3.4）。

初始条件：液体单元为槽液成分；乳剂与保护层单元除明胶外为 0（干膜）或清水（预湿）；X = Y = P = Λ = D = 0。

#### 2.4.9 药液槽

```
dC_k,bulk/dt = −1e-15 · (A_film / V_tank) · J̄_k(t)        [mol L^-1 s^-1]
```

`J̄_k` 为进入胶片的面积平均通量（mol L^-1 um s^-1，由当时狄利克雷面（δ_max 或 δ_agit）的有限体积通量精确得到），`A_film` 以 um² 计，`V_tank` 以 L 计，1e-15 把 um³ 换成 L。空气氧化项留出接口、取 0（`assumed`）。每个搅动开始时另有一次跳变 `ΔC_k,bulk = +1e-15 · (A_film / V_tank) · ΔĪ_k`（ΔĪ_k 为 §2.4.8 上段库存差的面积平均；进入胶片为正的通量取负号，所以富余的 Br⁻ 使槽液 Br⁻ 增加）；在 `flux_mean`、`J̄_blk`（§3.7）中它记作该步的冲量 `−ΔĪ_k / h`，所以重新积分槽 ODE 时自动计入。守恒测试（`test_component_conservation`、`test_tank_mass_balance`）把搅动期上段算作槽的一部分。图像中的面积构成见 §3.7。

#### 2.4.10 输出、注册与作用范围

显影结束时按深度合并写字段（§2.0.3）：

```
silver_c          = 1e12 · m_Ag,c · Σ_z (ν_c/n_z) X_c,z                 [g/m2]
initiated_c       = (1/n_z) Σ_z F_c,z                                    [1]
silver_m2_c       = (1/n_z) Σ_z Y_c,z                                    [1]
silver_physical   = 1e12 · Σ_z P_z Σ_c (ν_c/n_z) m_Ag,c                  [g/m2]
halide            = m_Ag,coat − Σ_c silver_c − 0.10787 · Σ_z D_z Δ_z     [g/m2]
growth_rate       = (1/n_z) Σ_z κ̃_z(t_end)                               [s^-1]
```

（`D_z Δ_z` 的单位 mol L^-1 um，乘 107.87 g/mol × 1e-15 L/um³ × 1e12 um²/m² = 0.10787 得 g/m²。）

- 注册：`register_step("bw-mechanistic", StepKind.DEVELOP, interaction_range_um=None)`。流水线不分块；步骤内部按 §3.6 的规范块与 §3.7 的预计算执行，结果与 `tile_px`、进程数无关。
- 读取 `emulsion/exposure_*`、`frame/presoaked`（没有这个字段时按干膜处理）、`meta["exposure"]`；删除 `exposure_*`。
- `meta["develop"]` = `{"bath", "dilution", "temperature_c", "time_s", "agitation", "tank", "profile", "chem_cell_um", "halo_um", "block_cells", "n_steps", "min_scaled_concentration", "prepass": {"bins", "passes", "delta"（每遍一个）, "contraction"（L̂，有重算时）, "depletion_change", "ph_change"}, "bulk_initial", "bulk_final", "zero_exposure": {...每档 silver/initiated/silver_m2 与 silver_physical、halide}}`。
- 预湿步骤 `presoak`（作用范围 0）写广播字段 `frame/presoaked = True`，对化学没有其他影响。

#### 2.4.11 参数

胶片侧（`emulsion.development.*` 与层厚）：

| 符号 | 字段 | 单位 | 类别 | 先验/界 | 出处 |
|------|------|------|------|---------|------|
| h_fog | `fog_hazard_ratio` | 1 | fitted | 1e-6–1e-1（对数） | KIN §1.8 |
| α_sw | `swelling_ratio` | 1 | assumed 3 | 2–5 | KIN §3.3（Iwano"正常坚膜"） |
| L_oc | overcoat `thickness` | um | assumed 1.0 | 0.5–2 | OPT §1.4 |
| x_I | `grain.iodide_fraction` | 1 | assumed | §2.0.2 | LAT §2.1–2.2 |

共用化学（`data/chemistry/kinetics.yaml`，§5.3）。显影剂类 a ∈ {`metol`, `hydroquinone`, `pyrazolidone`, `ascorbate`, `p-aminophenol`}：

| 符号 | 字段 | 单位 | 类别 | 先验/界 | 出处 |
|------|------|------|------|---------|------|
| k_et,a | `agents.<a>.k_et` | s^-1 | fitted；对氨基酚 assumed 取米吐尔值 | 1e-4–10（对数） | KIN §1.10 |
| K_a | `agents.<a>.K_ads` | L/mol | fitted | 1–1e4（对数） | KIN §1.5、§1.10 |
| λ0_a | `agents.<a>.lambda0` | s^-1 | fitted | 1e-4–10（对数） | KIN §1.4、§1.10 |
| Ea_a | `agents.<a>.activation_energy` | kJ/mol | 米吐尔、HQ fitted（先验 61.8、67.5 ± 10，literature）；其余 assumed 65 | — | KIN §1.6 |
| k_des,E | `agents.<E>.k_des` | s^-1 | fitted（米吐尔、吡唑烷酮） | 1e-3–10（对数） | KIN §2.1 |
| k_reg,R | `agents.<R>.k_reg` | L/(mol s) | fitted（HQ、抗坏血酸） | 1e-2–1e4（对数） | KIN §2.1 |
| pK_pyr | `agents.pyrazolidone.pka` | 1 | fitted | 7.5–9.5 | DEV §6.1 |
| K_Br | `halides.K_Br` | L/mol | fitted | 10–1e4（对数）；θ_Br,sat ≤ 0.18 作合理性检查 | KIN §1.5（Herz） |
| ΔH_Br | `halides.adsorption_enthalpy` | kJ/mol | literature −45（固定） | −60 – −30（只作敏感性） | KIN §1.5（Herz 1980：立方/八面体 AgBr 上 −10.7/−11.5 kcal/mol；标准态有出入）。`K_X(T) = K_X,20 · exp(−(ΔH_Br/R)(1/T − 1/T_r))`，K_I 用同一因子（`assumed`）；18.3–23.9 °C 之间 K_Br 变约 1.4 倍 |
| K_I/K_Br | `halides.iodide_ratio` | 1 | assumed 30 | 10–1000 | KIN §4.1（碘化物抑制更强，无数值） |
| u_gel | `barriers.u_gel` | 1（kT） | fitted | 0–10 | KIN §1.1、§1.4 |
| u_X | `barriers.u_halide` | 1 | fitted | 0–20 | KIN §1.4、§1.10 |
| u_LI | `latent.u_li` | 1 | fitted | 0–5 | KIN §1.10（Ohzeki） |
| n_sat | `latent.n_sat` | 1 | literature 5 | — | KIN §1.1 |
| ζ_s | `solvent.latent_attrition` | 1 | fitted | 0–5 | KIN §2.2（定性） |
| k_diss | `solvent.k_diss` | um/s | fitted | 1e-6–1e-2（对数） | KIN §2.2 |
| C_sat,ref | `solvent.c_sat_ref` | mol/L | fitted | 1e-5–1e-2（对数） | KIN §2.2：Carroll & Hubbard 30 °C 时 AgBr 7.5e-5 M @ 0.010 N（= 0.005 M）Na₂SO₃；按 ν_s 外推到 0.1 M 为 4.5e-4（ν_s = 0.6）至 1.5e-3 M（ν_s = 1），作量级参照（作者称反应温度下"可能大几倍"） |
| ν_s | `solvent.sulfite_exponent` | 1 | fitted（kinetics.yaml 初值 1.0，`assumed`） | 0.6–2，线性，均匀先验 | Carroll & Hubbard 0.6–5 mM 三点的斜率 0.62（KIN §2.2，literature）只在低亚硫酸盐段成立；D1 用在 40–800 mM，Ag(SO₃)ₙ 络合的质量作用律给 1（不加溴）到 2（溴过量），所以由训练集的亚硫酸盐水平（D-76、D-76 1+1、XTOL、T-MAX、HC-110 B）拟合，外推性质记入 §9 第 9 条 |
| ν_b | `solvent.bromide_exponent` | 1 | assumed（外推）0.4 | 0.3–1 | 低亚硫酸盐下加 4e-4 N KBr 溶解度降一半（KIN §2.2）；质量作用律在溴过量时给 1 |
| n_lig | `solvent.ligands` | 1 | assumed 2 | 1–3 | — |
| k_pd | `solvent.k_pd` | um/s | fitted | 1e-6–1e-2（对数） | KIN §2.2 |
| δ_max | `transport.delta_max` | um | fitted | 200–500 | KIN §3.4 |
| δ_agit | `transport.delta_agit` | um | fitted | 20–100 | KIN §3.4（US 5,315,338 推得 6–104 µm） |
| T_½ | `transport.swelling_half_time` | s | assumed 20 | 10–40 | KIN §3.3（secondary） |
| b_η | `transport.viscosity_ionic_coefficient` | L/mol | assumed 0.125 | 0–0.2 | DEV §6.3 |
| D_H、D_OH | `transport.proton_diffusion`、`transport.hydroxide_diffusion` | um²/s（20 °C） | assumed 8130、4600 | — | 示踪值 9.31e-5、5.27e-5 cm²/s（25 °C，recalled）；只用于 TOTH 通量（§2.4.8），换算同 §2.3.7 |
| g_OH | `gelatin.groups_per_gram` | mol/g | literature 4.8e-4 | — | KIN §3.2（Iwano） |
| pK_gel | `gelatin.pka` | 1 | assumed 10.5 | 9.5–11.5 | — |
| pH_gel,0 | `gelatin.initial_ph` | 1 | assumed 5.5 | 5–6.5 | — |

药液配方的干扰参数（成品 D-76 的 B₂O₃、T-MAX Developer 的 HQ/Dimezone S/初始溴化物、旧 HC-110 的各组分）在药液文件里，见 §5.5、§7.1。

#### 2.4.12 D-76（20 °C）各项的量级

按 §12.5 的脚本。银 3 g/m²（`assumed`），高光区 300 s 显影一半；溶胀层取任务要求的 10–20 µm（另列 15 µm），边界层取 50–300 µm：

| 过程 | 量级 | 对离散化与降阶的含义 |
|------|------|----------------------|
| 溶胀层中的银浓度与显影速率 | 10 µm：2.78 M、4.6 mM/s；15 µm：1.85 M、3.1 mM/s；20 µm：1.39 M、2.3 mM/s | 面积通量与层厚无关：J = 46.4 mM·µm/s（KIN §3.6 的 5e-5 mol m⁻² s⁻¹） |
| 凝胶内扩散时间 L²/D | 10/15/20 µm：HQ 1.0/2.3/4.0 s，Br⁻ 0.11/0.25/0.44 s（r_gel = 0.5） | 远短于显影时间：凝胶近乎准稳态，2–4 个深度单元足够 |
| HQ 在凝胶中的 Thiele 模数（0.5 HQ/Ag） | φ² = 0.05/0.08/0.10，有效因子 0.98/0.97/0.97；按 1 HQ/Ag 时 φ² 加倍 | 深度方向浓度差 ≤ 5 % |
| 边界层扩散时间 δ²/D | δ = 50/100/200/300 µm：HQ 2.8/11/46/102 s，Br⁻ 1.4/5.5/22/49 s | 与搅动间隔（30/60 s）可比：边界层必须动态求解，不能准稳态 |
| 边界层两侧的浓度差 | δ = 50/100/200/300 µm：Br⁻ 富余 1.3/2.5/5.1/7.6 mM（KBr 1 g/L = 8.4 mM），HQ 亏损 1.3/2.6/5.3/7.9 mM（3/6/12/17 %） | 补偿与邻接效应主要由边界层决定（KIN §3.6） |
| 胶片表面 pH 降（0.5 H⁺/Ag） | 条件 pKa（混合常数，含 D(I)）：0.04/0.07/0.14/0.21；理想 pKa：0.01/0.03/0.05/0.08 | 缓冲容量对活度模型敏感：pH 8.5 时总 β 由 92 mM/pH（理想：亚硫酸盐 81、硼酸盐 6、HQ 4、米吐尔 0.5）降到 33 mM/pH（§2.3.3 只含 D(I) 的混合常数：17、8、7、1），约 2.8 倍；凝胶单元另有明胶约 2–4 mM/pH；活度标定必须先于动力学标定 |
| 封闭单元全部显影的酸量 | 15 µm 时 0.93 M H⁺，对 β ≈ 0.03–0.09 M/pH | 不与外界交换时 pH 会降几个单位：pH 场是一阶物理 |
| 横向长度：凝胶交换长度 √(D τ_v) | ≈ 15 µm（Br⁻ 与 HQ 都是） | 化学单元 ≤ 24 µm 时子单元邻接仍有正确的饱和行为（§3.8） |
| 横向长度：液体内 30/60 s | Br⁻ 234/331 µm，HQ 162/229 µm | 规范块晕边约 1–1.4 mm（§3.6） |
| 显式稳定性上限 | 竖直（Δz = 2.5 µm）1.7 ms；横向（24 µm）79 ms | 必须隐式 |
| 传输算子的刚性比 | 最快 68 s^-1 / 最慢 0.022 s^-1 ≈ 3e3（§12.6 的原型柱） | L 稳定的 Rosenbrock |

#### 2.4.13 极限情况（D1c 测试）

| 测试 | 条件 | 断言 |
|------|------|------|
| `test_diffusion_into_slab` | 无反应、s ≡ 1、顶端恒浓度；乳剂 n_z = 8/16/32 | 与级数解（有限厚度）之差：tD/L² ∈ {0.01, 0.1, 0.3} 上相邻加密的误差比 ≥ 3（观测阶 ≥ 1.6）；tD/L² ≥ 1 时 n_z = 8 的误差 ≤ 2e-3 C₀（8 单元实测：0.01 时 7.3e-2，0.3 时 3.3e-3，约 0.5 时到 2e-3） |
| `test_reaction_diffusion_effectiveness` | 一阶消耗、表面恒浓度 | 有效因子 = tanh φ / φ，误差 ≤ 1e-3 |
| `test_moving_reaction_front` | 瞬时反应、不动反应物 C_s ≫ C₀；专用细网格（200 个乳剂单元，层厚取前沿最大行程的 2 倍） | 由累计消耗量定义的前沿位置 `x_f = (已消耗的不动反应物)/C_s` 与 √(2 D C₀ t/(ν C_s)) 之差 ≤ 3 %（KIN §5.7） |
| `test_silver_conservation` | 任意条带 | 化学银 + 物理银 + 卤化银 + 膜内与流失的溶解银 = 涂布量，相对差 ≤ 1e-12 |
| `test_component_conservation` | Br、显影剂、亚硫酸根、TOTH | 膜 + 液体 + 槽的库存变化 = 积分源项，相对差 ≤ 1e-12 |
| `test_gel_equilibrium_no_toth_flux` | 未曝光、h_fog = 0、k_diss = k_pd = 0（无溶解与物理显影）、s ≡ 1；凝胶单元的可动总量取槽液值，TOTH = TOTH_mob,bath + T_gel n̄_gel(pH_bath)；显影 300 s | 各单元 pH 与槽液 pH 之差 ≤ 1e-10；每个面的 \|J_toth\| ≤ 1e-14 mol L⁻¹ µm s⁻¹ |
| `test_immobile_buffer_retards_ph_front` | 单一可动缓冲族（D）加明胶、无反应、s ≡ 1、半无限凝胶；顶端 TOTH 小扰动（ΔpH ≤ 1e-3） | pH 扰动剖面与 D_app = D/(1 + β_gel/β_mob) 的 erfc 解之差 ≤ 1e-2 × 扰动幅度 |
| `test_surface_ph_drop_matches_buffer` | 小扰动 | 表面 pH 降 = ΔTOTH/β，误差 ≤ 2 % |
| `test_lh_fractional_order` | 只有 HQ、无 Br | ∂ln κ̃/∂ln C = 1 − θ，误差 ≤ 1e-6；θ = 0.4 时为 0.6 |
| `test_bromide_restrains` | 加初始 Br⁻ | 固定时刻 X 单调下降，Λ 增长变慢 |
| `test_superadditivity_emerges` | 米吐尔 + HQ | k_reg > 0 时合用速率 > 单用之和；k_reg = 0 且 Φ_E → 0 时相等（≤ 1e-12） |
| `test_arrhenius_scaling` | 输运极快、单一显影剂；形态分布与温度无关（全部电离焓与 dpK/dT 取 0，A(T) 与 pK_w 固定在 20 °C 值），ΔH_Br = 0，无溴化物 | 达到 X = 0.5 的时间按 exp(Ea/R·Δ(1/T)) 缩放，误差 ≤ 1e-3 |
| `test_agitation_limits` | δ_agit、δ_max → 0；以及无搅动大 δ_max | 前者与"表面直接取槽液"之差 ≤ 1e-4；后者早期表面亏损 ∝ t^½（指数 0.5 ± 0.02） |
| `test_strip_smooth_in_transport` | D-76 小罐条带，标定档与图像档；δ_agit = 39/41/90 µm、δ_max = 300 µm 附近，扰动 1e-4/1e-5/1e-6（相对） | 各级 X 的二阶差分 / 扰动² 在三个尺度上相差 < 10 %（与 `test_forward_model_smooth` 同一判据）；同一 δ_agit 下三个档位的搅动层厚度都等于 δ_agit |
| `test_dry_film_bromide_free_start` | 干膜、D-76（无溴）、图像档 h₀ = 0.1 s 与标定档 | 无 NaN；全程 min C/s_k ≥ −1e-3（§3.3）；记录值写入 `meta["develop"]`；`test_silver_conservation` 的平衡照样成立 |
| `test_tank_infinite_volume` | V = 1e9 L | 槽液成分变化 ≤ 1e-12 |
| `test_tank_mass_balance` | 有限体积 | 槽内变化 = −∫ J̄ A / V，误差 ≤ 1e-12 |
| `test_zero_exposure_fog_only` | m ≡ 0 | X 只来自 h_fog；h_fog = 0 时 X ≡ 0 |
| `test_long_development_limit` | 强曝光、显影剂充足、5000 s | 可显影档 X → 1，误差 ≤ 1e-3 |
| `test_presoak_removes_delay` | 预湿与干膜 | 干膜早期 X 较小；T_½ → 0 时两者相等（≤ 1e-6） |
| `test_rosenbrock_second_order` | 固定步长逐次减半 | 相邻差之比 3.5–4.5 |
| `test_fixed_schedule_vs_reference` | D-76 条带，小罐与大罐搅动，标定档与图像档时间表（§3.4 的事件规则：只在搅动开始处复位）与参考解 | 标定档 max \|ΔX\| ≤ 1e-3，图像档 ≤ 2e-3（§12.6 的 2.3e-4/4.8e-4 来自瞬时复位的原型，不直接适用于有持续时间的搅动期） |
| `test_member_schur_elimination_exact` | 每柱一个成员 | 与完整耦合线性解之差 ≤ 1e-12 |
| `test_kernel_reduction_exact` | 每柱一个成员；以及一个单元多个成员、最近单元耦合 | 核合成的 X、Y、P 与直接积分成员方程之差 ≤ 1e-12；单元化学只依赖 p̄（成员重排、拆分不改变结果，≤ 1e-13） |
| `test_ensemble_vectorisation` | N 根独立柱 | 与逐根计算逐位相同 |
| `test_develop_independent_of_tiling_and_workers` | `tile_px` ∈ {None, 50}，`workers` ∈ {1, 3} | 逐位相同 |
| `test_uniform_image_matches_strip` | 均匀曝光的小图，宽于两倍晕边加 2 个规范块；无限大槽（两边都是恒定槽液）；两边同一档位与时间表 | 离画幅边缘超过两个晕边的内部像素，与单柱结果之差 ≤ 1e-10（画幅外是未曝光片，边缘附近本来就不同） |

### 2.5 成像（阶段 4）

#### 2.5.1 Nutting 与遮盖力

N 个面积为 a 的不透明颗粒随机放置时 `T = exp(−N a)`，`D = 0.4343 N a`（Nutting 1913；布尔模型的覆盖公式，OPT §4.1）。按档求和，并把遮盖力 CP（每 g/m² 银的密度）作为每档的参数：

```
CP_c = CP_ref · (d_c / d_ref)^(−β_CP),   d_ref = 1 um
a_c = ln10 · 1e12 · CP_c · m_Ag,c                      [um²]（一颗完全显影的 c 档颗粒的有效投影面积）
D_Ag = Σ_c CP_c · silver_c + r_CP · CP_ref · silver_physical
```

部分显影的颗粒的有效面积取 `a_c · x`，即遮盖力在显影过程中不变（`assumed`；Solman & Farnell 1973 指出 CP 沿曲线变化，OPT §4.2）。

#### 2.5.2 光谱透射

```
T_b(x) = 10^−[ D_base,b + D_Ag(x) · n_Ag,b + CP_hal · halide(x) + D_dye,b · residual_dye(x) ]
D_base,b = D_g（135 灰片基，中性）+ 0（120 的背涂层在冲洗中去除，`assumed`）
n_Ag,b = 1（银的光谱中性，assumed）
```

D_Ag 是"挡住的光"的比例，最接近漫射密度（OPT §4.3、CAL §1.4），所以所有标定都在漫射视觉密度上进行。

#### 2.5.3 参数

| 符号 | 字段 `emulsion.development.*` | 单位 | 类别 | 先验/界 | 出处 |
|------|------------------------------|------|------|---------|------|
| CP_ref | `covering_power` | 1/(g/m2) | fitted | 0.3–1.0 | OPT §4.2 |
| β_CP | `covering_power_size_exponent` | 1 | assumed 0.5 | 0–1 | OPT §4.2（CP 随粒径下降） |
| r_CP | `physical_covering_ratio` | 1 | assumed 0.5 | 0.2–1 | OPT §4.2 第 3 点 |
| CP_hal | `halide_density` | 1/(g/m2) | assumed 0.3 | 0.1–0.6 | 无数据；只在定影不足时起作用 |
| D_dye,b | `dye_density` | 1 | assumed 0 | 0–0.1 | OPT §1.3（T-Max 粉红残色，定性） |

m_Ag,coat 与 CP_ref 只在长时间显影的最大密度里以乘积出现（Kodak 曲线不到饱和，CAL §3.4）。其他地方两者可以分开：m_Ag,coat 决定银浓度 C_c,z（§2.0.2），因而决定每单位密度放出的 Br⁻、I⁻、H⁺ 与消耗的显影剂（∝ 1/CP），即局部与槽的耗尽、补偿强度与邻接效应（MTF > 100 %）的幅度；它还经 N_tot 决定层吸收与 L_e。所以标定在胶片块中同时拟合 m_Ag,coat 与 CP_ref（先验 2–6 g/m²），由曲线形状与显影时间系列区分（MTF 峰——T-Max 19 c/mm 处 113 %、Tri-X 10 c/mm 处 111 %，OPT §2.4——也依赖它，但峰段只报告，§7.2，可作事后检查）；阶段 F 的剖面似然判为不可辨识时固定在 4.0，并按 §7.5 做敏感性运行。

#### 2.5.4 测试（D1d）

| 测试 | 断言 |
|------|------|
| `test_nutting_linear_in_silver` | 无颗粒噪声时 D_Ag 对银量线性，相对差 ≤ 1e-12 |
| `test_covering_power_area_identity` | 由 a_c 与 ν_c 算的 0.4343 Σ ν_c a_c X_c 与 Σ CP_c silver_c 相等，≤ 1e-12 |
| `test_transmittance_factorisation` | 中性介质下 T_b 与 b 无关；视觉密度 = Type 2 密度，≤ 1e-12 |

### 2.6 停显、定影、水洗、干燥（阶段 5）

都是逐点的一阶步骤，作用范围 0，真实改变状态：

| 步骤 | 方程 | 改变的字段 | 参数（药液或胶片文件） |
|------|------|------------|------------------------|
| `stop` | 显影剂在 τ_s 内被中和，κ(t) = κ_c,end e^{−t/τ_s}；以下式更新每档的 X、Y（换算成 silver、silver_m2），κ_c,end = (ℓ_ref S_c/V_c)·`growth_rate` | `silver`、`silver_m2`（`initiated` 不变：停显中不再引发）；删除 `growth_rate` | `stop.neutralisation_time` τ_s：醋酸停显 2 s，清水 15 s（`assumed`） |
| `fix` | `halide ← halide · exp(−t/(f_fix τ_fix))`；`residual_dye ← residual_dye · exp(−t/τ_dye)`；`fixer_residual ← 1` | `halide`、`residual_dye`、`frame/fixer_residual` | 定影液 `fix.clearing_time_constant` τ_fix = 20 s，`fix.dye_removal_time` τ_dye = 60 s（`assumed`）；胶片 `development.fixing_factor` f_fix：Tri-X 1、T-Max 2（`assumed`） |
| `wash` | `fixer_residual ← fixer_residual · exp(−t/τ_w)`；`residual_dye ← residual_dye · exp(−t/τ_dye,w)` | `frame/fixer_residual`、`residual_dye` | 清水 `wash.time_constant` τ_w = 120 s，τ_dye,w = 300 s（`assumed`） |
| `dry` | 无（记录） | — | — |

停显的闭式解（F 不变，`dX/dK = F − X`，`dY/dK = 2(X − Y)`，以 `K(t) = κ_c,end τ_s (1 − e^{−t/τ_s})` 为自变量）：

```
X = F − (F − X₀) e^{−K}
Y = F − 2(F − X₀) e^{−K} + (F − 2X₀ + Y₀) e^{−2K}
```

K = 0 时还原为 (X₀, Y₀)，K → ∞ 时 X、Y 都趋于 F（已引发的颗粒全部显影完，x² = 1）。这些常数都没有参考数据，D1 不验证它们，只保证"缩短或跳过有相应后果"（总设计 §5）：跳过定影会留下卤化银密度，跳过停显会多显影。

`meta["post"]` 记录每步的常数与时间。测试（D1d）：`test_stop_extra_development_monotone`（τ_s 越长银越多）、`test_skip_fix_leaves_halide_density`（D 增加 CP_hal·halide）、`test_post_steps_pointwise_tiling`（任意 `tile_px` 逐位相同）、`test_stop_y_closed_form`（与数值积分之差 ≤ 1e-10）。

### 2.7 颗粒（阶段 6）

#### 2.7.1 统计模式（全幅）

计数档与均值档：c 档满足 `p ≥ 4 r_c`（`r_c = √(a_c/π)`，完全显影的盘半径，每帧固定，与密度无关）时进入计数积（集合 𝒞）；否则是"均值档"，只取确定的因子 `exp(−μ_c q_c)`，它的方差全部由 ε_pix 提供。均值档保证 `0 ≤ 1 − q_c` 不被违反，不再报错；哪些档是均值档写入读出结果的 `Readout.meta["grain"]["mean_only_classes"]`（下文的 `meta["grain"]` 都指它）。

每像素每档：

```
μ_c = ν_c · A_pix · F_c                         已显影颗粒数的均值（A_pix = p²）
N_c ~ Poisson(μ_c)                              由 rng.py 的均匀数逆变换得到（只对 c ∈ 𝒞）
x̄_c = X_c / F_c ,   ⟨x²⟩_c = Y_c / F_c         已引发颗粒的显影分数一阶、二阶矩
φ = r_CP · CP_ref · silver_physical / Σ_c CP_c silver_c     物理显影银按比例加大已显影的盘（Σ = 0 时 φ = 0，物理银按均值计入）
q_c = a_c (1 + φ) x̄_c / A_pix
T_pix = Π_{c∈𝒞} (1 − q_c)^N_c · Π_{c∉𝒞} exp(−μ_c q_c) · (1 + ε_pix)  截断到 [0, 1]
```

- 计数部分：均匀放在像素内的 N 个圆盘，任意一点未被覆盖的概率是 `(1 − q)^N`，对泊松 N 求期望得 `exp(−μ q) = exp(−ν A F a x̄/A)`，正好是 Nutting 均值（§2.5），所以统计模式的平均密度与 §2.5 和显微镜模式完全一致。φ 使 `0.4343 Σ_c μ_c q_c = D_Ag`（§2.5.1，含物理银）。物理显影银沉积在已显影颗粒上（§2.4.4），所以它放大盘而不是作为无噪声的均匀密度（`assumed`：按银量统一放大，不按各档表面分配）；片基、卤化银、残留染料只乘平均透射（§2.5.2）。
- 泊松而不是固定总数的二项：颗粒位置是泊松分布的（显微镜模式与布尔模型的前提），已显影颗粒是对泊松颗粒的二项稀疏化，闭式仍是泊松（§11 第 1 条）。
- `ε_pix` 补上计数模型漏掉的方差：圆盘之间的重叠（布尔协方差）、显影分数的分散，以及均值档的全部方差。

```
Var_count  = T̄² · (exp(Σ_{c∈𝒞} μ_c q_c²) − 1)        q_c 用实际抽样的 x̄_c（不是 ⟨x²⟩_c）：E[(1−q)^{2N}] = exp(−2μq + μq²)
Var_bool   = (T̄² / A_pix) · ∫ (exp(γ_mix(τ)) − 1) 2πτ dτ
γ_mix(τ)   = Σ_c ν_c F_c · γ(τ; r̃_c),   r̃_c = √(a_c (1 + φ) ⟨x²⟩_c^{½} / π)，γ 为两圆盘交叠面积（含全部档）
Var(ε_pix) = max(Var_bool − Var_count, 0) / (T̄² + Var_count) ,   ε_pix = √Var(ε_pix) · z,  z ~ N(0, 1)
```

分母是计数部分的二阶矩 `E[T_count²] = T̄² + Var_count`：ε 与计数独立，`Var(T_count (1 + ε)) = Var_count + E[T_count²] Var(ε)`，所以截断前的像素方差正好等于 Var_bool。显影分数的分散（⟨x²⟩ > x̄²）只在 Var_bool 里（通过 r̃_c），由 ε 补上。

Var_bool 是"大像素"（白噪声等效）方差：相互独立、各带这一方差的像素，在任何远大于盘的聚合尺度上（48 µm 孔径、降采样的像素）给出布尔模型的正确方差，与 p 无关；盘与像素可比时，单个像素的方差高于真实的像素方差，像素间的相关不表示（`assumed`）。聚合尺度为直径 D_ap 的圆孔时，真实方差比这一极限低约 `b_ap = 1.15 r̃_eff / D_ap`（常半径盘的数值核对：D_ap = 10/20/40 r 时 11.5/5.8/2.9 %），`r̃_eff = √(Σ_c ν_c F_c (π r̃_c²)² / (π Σ_c ν_c F_c π r̃_c²))`。带孔径的读出与 §6.6 把 b_ap 写入 `meta["grain"]`，b_ap > 0.05 时警告。均值档多时 ε 很大，`(1 + ε)` 的截断会抬高平均透射；截断的像素比例写入 `meta["grain"]["truncated_fraction"]`。

γ 的积分在 [0, 2 r̃_max] 上用 48 点 Gauss–Legendre，常半径时复现 OPT §5.4 的 κ(D) 表（1.055、1.18、1.34、1.89、2.86、4.66，§12.7 逐位复现到三位小数）。像素远大于颗粒时 Var_bool ≈ κ(D)·Var_Siedentopf，κ 是 OPT 推出的修正；忽略它会使 D = 1 处的 σ_D 低 √1.89 ≈ 1.37 倍（OPT §0 第 1 条）。

随机数：流名 `grain/emulsion/<c>` 给第 c 档的泊松逆变换（每像素一个均匀数），`grain/emulsion/overlap` 给 z（每像素一个均匀数，`scipy.special.ndtri`）。泊松逆变换：μ ≤ 30 用累积分布直接查找（向量化，上限 μ + 10√μ + 10）；μ > 30 用 Cornish–Fisher 展开 `k = ⌊μ + √μ z + (z² − 1)/6 + (z³ − 7z)/(72√μ) + ½⌋`。结果只依赖（种子、流名、像素绝对坐标），与分块无关。

#### 2.7.2 显微镜模式（D1f，参照）

小块（≤ 2 mm 见方）内显式放置每一颗颗粒：

1. 每档颗粒在三维乳剂中按泊松过程放置（面密度 ν_c，分数深度均匀）。
2. 每颗粒的吸收光子均值取其所在像素、所在深度单元的 m_c,z；团簇大小 n 按 p_c,z(n) 抽样。
3. 显影：在同一小块上跑显影引擎（化学单元 = 显微网格，§3.8 的参考设置），得到每个深度单元的两个单调累计量 Λ_z(t) 与 `K_z(t) = ∫₀ᵗ κ̃_z dt'`（`run_ensemble(record_cells=64)`：在最接近 64 个等间隔时刻的步上记录，不改时间表；时刻之间用单调 PCHIP 插值）。每颗粒的引发时刻 t₀ 由 `h_n Λ_z(t₀) = E`、`E ~ Exp(1)` 求出，显影分数 `x = 1 − exp(−(ℓ_ref S_c/V_c)(K_z(T) − K_z(t₀)))`，再按停显式延续。记录量：2 mm 小块、4 µm 单元、n_z = 8 时 2 × 64 × 2.5e5 × 8 × 8 B ≈ 2 GB（逐步记录要 10–32 GB）。均匀曝光（跨模式颗粒度测试的 D = 0.3、1.0、2.0）只算一根柱，广播到全部单元。
4. 投影：已显影颗粒在胶片面上是圆盘，面积 `a_c x (1 + φ)`（φ 同 §2.7.1，由所在像素的显影字段算出），中心为颗粒的 (x, y)。三维泊松颗粒的投影仍是平面泊松，所以这是 Newson 等（2017）的布尔模型，只是盘的密度和大小来自模拟的颗粒群体（OPT §6.1）。
5. 渲染：在像素宽 `p_r = min(r̃_min/4, 0.05 um)` 的光栅上按"到圆边距离"计算抗锯齿覆盖率；测量孔径（48 µm 圆）内的透射为未覆盖面积分数。

两种模式共用潜影表、显影引擎、遮盖力与随机数流（显微模式用 `micro/<c>`、`micro/position` 等流，按颗粒所在的显微格坐标作键）。

#### 2.7.3 参数

颗粒统计没有新的自由参数：ν_c、F_c、X_c、Y_c、a_c 都来自前面各阶段。T-Max 100 的 RMS 颗粒度 8 是胶片标定目标（约束粒径与遮盖力的组合），Tri-X 的 17 是跨显影液检查（§7.6）。

#### 2.7.4 测试

D1d：`test_poisson_sampler_moments`（μ ∈ {0.1, 3, 30, 300, 3000}，各 1e6 样本，均值与方差在 4σ 内）、`test_poisson_sampler_deterministic`（逐位可重复，与分块无关）、`test_statistical_mean_is_nutting`（大块均匀区的平均透射与 exp(−Σ ν F a x̄ (1 + φ)) 之差在 3 个标准误差内，含 silver_physical > 0 的情形）、`test_overlap_kappa_table`（常半径 κ(0.1…2.0) 与 OPT 表之差 ≤ 0.005）、`test_siedentopf_limit`（D → 0 时 Var_bool / (T̄² Σ_c μ_c a_c² ⟨x²⟩_c / A_pix²) → 1，D = 0.001 时差 ≤ 0.1 %（该比值即 κ(D)，现有 `kappa(0.001)` = 1.00053）；x ≡ 1 时 Var_bool/Var_count 同样 → 1）、`test_statistical_variance_broad_x`（单档、a_c/A_pix = 0.02（p ≈ 12.5 r_c）、X = F/2、Y = F/3（x 在 [0, 1] 上均匀）、D ∈ {0.3, 1}（D = 1 时 σ_ε ≈ 0.18）：1e6 个像素的样本方差与 Var_bool 之差在 4 个标准误差内，截断比例 < 1e-4；Var_bool 的 Siedentopf 极限与随机盘面积的精确值 T̄² Σ μ a² ⟨x²⟩/A² 的相对差 ≤ 1e-4）、`test_mean_only_classes`（4 r_c > p 的档不抽泊松数、均值仍是 Nutting、Var_count 不含该档）、`test_trix_prior_centre_statistical`（Tri-X 先验中心参数（d_50 0.9 µm、σ_d 0.45、CP_ref 0.6、K = 6）在 p = 6 µm 与 p_g = 4 µm 下 D = 1 的均匀块不报错，`meta["grain"]` 含 mean_only_classes、b_ap、truncated_fraction；这组参数的截断比例大，只检查记录，不检查方差）。D1f：§10 的跨模式测试。

### 2.8 出口（阶段 7）

漫射透射密度，按光谱乘积加权（CAL §1.4）：

```
D_Π(x) = −log10( Σ_b Π_b T_b(x) / Σ_b Π_b )
```

| 设备 | Π_b | 类别 | 出处 |
|------|-----|------|------|
| `densitometer_visual` | S_A,b · ȳ_b（CIE 标准光源 A × CIE 1931 ȳ） | literature | CAL §1.4（ISO 5-3 §8.1 的等价关系，与 ISO 表格相差 ≤ 0.002 log） |
| `densitometer_type2` | 10^{P(λ_b)}，P 为对 ISO 5-3:1995 表 2 五个锚点（360 nm 4.280，400 nm 4.916，430 nm 5.000，500 nm 3.993，540 nm 1.0）在 log 空间的单调三次插值，[380, 540] nm 以外为 0 | assumed（近似；锚点 literature） | CAL §1.2 |
| `scan_standard` | S_D50,b · ȳ_b | assumed（固定约定） | — |
| `silver_density` | 直接返回 D_Ag | 诊断 | — |

参数：`aperture_um: float | None = None`（None：逐像素；给值则先对透射做圆孔平均再取对数）、`grain: bool`（密度计默认 False，用均值透射；扫描默认 True，用统计模式）。

几何是理想漫射（收集全部透射光，不模拟乳白玻璃的互反射，CAL §1.4）。D1 的银、片基、残留染料都按光谱中性处理，所以 Type 2 与视觉密度在模型里相等，Type 2 的近似对 ISO 感光度检验没有影响（§11 第 6 条）。

标准扫描的反相规则（固定，所有胶片与显影液共用，不随图像自适应）：

```
D_net = D_scan − D_g(format)                         只减片基，不减灰雾（灰雾是负片的真实密度）
Y_out = 0.18 · 10^((D_net − D_ref) / G_ref)
G_ref = 0.56,   D_ref = 0.10 + G_ref · log10(10 / 0.8) = 0.714
输出 RGB = (Y_out, Y_out, Y_out)，交给 D0 的 write_png16（截断到 [0, 1]）
```

D_ref 的来历：ISO 速度点在基雾以上 0.10，相机平均测光的胶片面曝光量是 10/S lux·s（D0 §3.3），速度点曝光量是 0.8/S，二者相差 log10(12.5) = 1.10；以 Kodak 的目标 CI 0.56 从速度点直线外推，测光平均落在 0.714，扫描后映射到 18 % 灰（`assumed`）。

测试（D1d）：`test_visual_product_matches_iso_anchors`（log10(A·V) 归一到 570 nm 后，450/500/560/600/650 nm 与 CAL §1.4 的 ISO 表值之差 ≤ 0.002）、`test_type2_anchors`（五个锚点逐点复现）、`test_neutral_film_equal_densities`、`test_scan_inversion_fixed_point`（D_net = D_ref 时 Y_out = 0.18 ± 1e-12）、`test_aperture_average_on_transmittance`（孔径平均先对 T 后取对数）。

---

## 3 数值方法

### 3.1 离散化档位

同一套代码按三个档位运行，`RunContext.profile` 选择（§4.7）。档位只改离散化，不改方程或参数。

| 量 | `reference` | `calibration` | `image` |
|----|-------------|---------------|---------|
| 乳剂单元 n_z | 8 | 4 | 2 |
| 液体单元 n_l，几何比 r_l | 16，1.25 | 8，1.4 | 4，2.0 |
| 尺寸档 K | 12 | 6 | 6 |
| Gauss–Hermite 节点 Q_GH | 21 | 9 | 9 |
| 潜影表节点（三次样条，§2.2.5；F 不建表，§2.4.2） | 401 | 201 | 201 |
| 液体单元分段 n_a + (n_l − n_a)（界面在 δ_agit，§2.4.1） | 8 + 8 | 4 + 4 | 2 + 2 |
| 时间推进 | 自适应 ROS2，rtol 1e-4（列测试另用 SciPy BDF rtol 1e-8 对照） | 固定表 (h₀, r, h_max) = (0.05 s, 1.4, 3 s) | 固定表 (0.1 s, 1.6, 5 s) |
| 横向化学网格 | 条带无横向；MTF/边缘条带 Δx ≤ 2 µm 或每周期 ≥ 32 格 | 无横向；MTF 用小信号的两柱对（§6.7） | 化学单元 Δx_c = 24 µm，固定，与像素间距 p 无关；网格锚定在画幅原点 |
| 成员 | 每柱 1 个 | 每柱 1 个 | 每像素 1 个；按面积重叠聚合到单元（p̄），输出时核在像素中心双线性插值（§3.8） |
| 药液槽 | 精确耦合 | 精确耦合（条带 + 其余胶片） | 预计算 + 事后检查（§3.7） |

### 3.2 算子分裂

每个时间步 Δt 用 Strang 分裂（Strang 1968；Hundsdorfer & Verwer 2003，KIN §5.3，recalled）：

```
y(t + Δt) = L(Δt/2) ∘ C(Δt) ∘ L(Δt/2) · y(t)
```

- `L`：横向扩散，逐单元层、逐组分精确求解（§3.5）。条带没有这一步。
- `C`：柱步。竖直输运、全部反应与核方程（§2.4.4）耦合在一起，用 ROS2 积分（§3.3）。反应与竖直扩散不分裂：显影速率依赖表面附近由扩散同时设定的浓度，分裂会在反应前沿产生误差（KIN §5.3）。
- 搅动的固定（h = δ_agit 界面上的狄利克雷值，其上的液体单元取槽液值，§2.4.8）属于 C，是 C 的边界条件。
- 槽液：条带与预计算中，槽液 ODE 用每步的面积平均通量显式更新一次（槽的时间尺度远大于 Δt）；图像块中槽液轨迹是给定输入，在步内线性插值。

分裂误差 O(Δt²)。横向耦合弱（梯度跨 10–100 个像素），`test_splitting_second_order`（§3.9）检查。

### 3.3 柱步：ROS2

Verwer、Spee、Blom、Hundsdorfer（SIAM J. Sci. Comput. 20, 1456, 1999）的 ROS2，γ = 1 + 1/√2（KIN §5.4，recalled）：

```
(I − γ h J) k₁ = f(y_n)
(I − γ h J) k₂ = f(y_n + h k₁) − 2 k₁
y_{n+1} = y_n + (3/2) h k₁ + (1/2) h k₂
误差估计（自适应时）：e = y_{n+1} − (y_n + h k₁) = (h/2)(k₁ + k₂)
```

ROS2 是 W 方法：对任意 J 近似都保持二阶（同一文献），所以 J 可以近似。

- 未知量：每柱各单元的组分浓度、每个乳剂单元的 Λ_z 与 D_z；该柱的核 G、H、Q（§2.4.4，成员不是未知量）。
- 形态分布在每次 f 求值时对每个单元做一次（§2.3.4，上一步的 pH 作初值）。pH 对总量的导数用 §2.3.4 的解析式代入链式法则。
- F 及其对 Λ 的导数按 §2.4.2 精确求值（每个乳剂单元 5 个指数）。max 用 softplus 平滑。
- J 的结构：凝胶单元是 K_c×K_c 稠密块（反应耦合）加 Λ、D；保护层与液体单元是对角块（只有各组分自己的扩散）；相邻单元间是对角的输运耦合。块三对角，块 Thomas 算法，只在凝胶单元做稠密 LU。
- TOTH 的输运项（§2.4.8）在 J 中只保留对 TOTH 自身的导数：`∂J_toth,ℓ+½/∂TOTH_ℓ = (Σ_{f 可动} D̄_f,ℓ+½ β_f,ℓ + D̄_H β_H,ℓ + D̄_OH β_OH,ℓ) / (β_ℓ h_ℓ+½)`（β_f,ℓ = −T_f dn̄_f/dpH，β_H = ln10 [H⁺]，β_OH = ln10 [OH⁻]，β_ℓ 含明胶），即局部缓冲容量加权的 D；TOTH 通量对其他总量 T_g 的导数略去（W 方法允许），所以保护层与液体单元仍是对角块。原 §2.3.7 的"缓冲容量加权 D"只出现在 J 里，不用于 f。
- 核用 Schur 补消去：核对自身状态的 J 是每个 (c, z, n) 的 3×3 下三角（G、H、Q），`k_u = (I − γhJ_uu)⁻¹(f_u + γh J_uc k_c)`，代回柱方程后柱的稠密块多出 `γ²h² J_cu (I − γhJ_uu)⁻¹ J_uc` 与右端 `γh J_cu (I − γhJ_uu)⁻¹ f_u`；柱方程只经 ⟨X⟩ = Σ_n p̄(n) G_n 等 p̄ 加权和依赖核，两项都是对 n 的向量化求和。图像档的成员只属于自己的单元（最近单元耦合，按面积重叠聚合），柱之间相互独立，同一批柱可以连续做完两个阶段，不需要全块的阶段屏障。
- 负浓度：不拒绝、不拆步。固定表的步长与状态无关（逐位一致与光滑目标函数都依赖这一点）；速率式用光滑正部 C⁺（§2.4.6），所以舍入级负值无害。每步记录 `min C_k,ℓ / s_k`（s_k = max(c_k,bulk(0), 1 mM)），全程最小值写入 `meta["develop"]["min_scaled_concentration"]`，低于 −1e-3 时警告（说明时间表对该条件过粗）。自适应 ROS2（参考档）把 `min C/s_k < −1e-6` 当作一次误差检验失败，按正常步长控制退回，整批柱一起退回。不做截断（截断破坏守恒）。

### 3.4 时间步表与搅动事件

所有柱、所有成员、预计算与各图像块共用同一张时间表（槽液轨迹与逐位一致都需要它）。两类时刻：

- **复位点**：0 与每个搅动期的开始。搅动开始时狄利克雷面从 δ_max 跳到 δ_agit，表面梯度突然变陡，是快瞬态（KIN §5.5），步长重新从 h₀ 开始。
- **对齐点**：每个搅动期的结束、快照时刻、终点。搅动结束时浓度剖面连续（上段原本就是槽液值），只是导数间断，之后的瞬态是向更厚液层的慢扩散；快照不改变方程。所以只让一步恰好落在对齐点上，不复位：若下一步 h_j = min(h₀ r^j, h_max) 会越过对齐点，这一步截到对齐点，下一步仍用 h_j 继续几何序列（截断的那一步不推进 j）；截断步短于 h₀/2 时改为把前一步延长到对齐点。

每段（两个复位点之间）内步长从 h₀ 起按比 r 增长，封顶 h_max。步数（`develop/schedule.py` 的规则；由 §12.6 的核对脚本算出）：

| 条件 | 图像档 (0.1, 1.6, 5) | 标定档 (0.05, 1.4, 3) | 旧规则（起止都复位）图像/标定 |
|------|---------------------|----------------------|------------------------------|
| 小罐 30 s/5 s，7.5 min | 209 | 314 | 293 / 454 |
| 小罐，12 min | 335 | 503 | 473 / 733 |
| 大罐 60 s/7.5 s，8 min | 160 | 248 | 216 / 328 |
| 大罐，12 min | 240 | 372 | 324 / 492 |
| CI(t) 运行到 1.5 × 12 min，每 15 s 快照：小罐 / 大罐 | 407 / 312 | 623 / 492 | 小罐 881 / 1333，大罐 816 / 1212 |

快照不复位以后，每 15 s 的快照不再增加步数（每 60 s 快照的步数相同）。

§12.6 的原型柱（6 个凝胶单元、10 个液体单元、400 s、每 30 s 瞬时复位）：

| 方案 | 步数（含拒绝） | 相对 BDF 参考的 max \|ΔX\| |
|------|----------------|---------------------------|
| 自适应 ROS2，tol 1e-2 | 297（+107 拒绝） | 7.5e-4 |
| 自适应 ROS2，tol 1e-3 | 940（+68） | 1.2e-4 |
| 固定 (0.1, 1.6, 5) | 178 | 4.8e-4 |
| 固定 (0.05, 1.4, 3) | 274 | 2.3e-4 |
| 固定 (0.2, 2.0, 10) | 110 | 1.1e-3 |
| 固定 (0.1, 1.6, 5)，每 60 s 复位 | 136 | 1.8e-4 |

固定表用更少的步数得到更小的误差，并且步长与参数无关，标定时有限差分雅可比不会受自适应步长噪声污染（CAL §3.3）。所以标定档与图像档用固定表，参考档用自适应 ROS2 作为对照。这张表的误差来自瞬时复位的原型；有持续时间的搅动期与上面的对齐规则下的误差由 `test_fixed_schedule_vs_reference`（§2.4.13）界定。

显影系列一次算完：要 6/8/10/12 min 四条曲线时只跑一次 12 min，在 6、8、10 min 记录快照，每个快照各自接停显（§2.6）后读出。CI(t) 曲线同样一次得到。

### 3.5 横向扩散

线性、逐层常系数，所以精确求解：

```
Ĉ(k, t + τ) = Ĉ(k, t) · exp(−D_k,ℓ |k|² τ)
```

- 图像块：DCT-II（正交归一，`scipy.fft.dctn`），块的外边界为零通量（诺伊曼）。块外边界只在晕边之外，误差由 §3.6 的晕边控制。
- MTF 与边缘的周期条带：FFT，周期边界，结果对正弦目标精确。
- 小信号 MTF 柱对（标定档，§6.7）：两根柱 a、b 代表 cos(kx) = +1 与 −1 两个相位，L 步只作用于差模：`C_a,b ← C̄ ± ½(C_a − C_b) · exp(−D_k,ℓ k² τ)`，C̄ = ½(C_a + C_b) 不变。调制幅度 ε → 0 时这就是周期条带的切线性（一阶谐波）响应；±ε 对称，差模的误差为 O(ε³)。
- 凝胶层的 D 取步中点时刻的 `r_k D_k s(t)`；液体层的 D 与时间无关。
- TOTH：横向扩散的是可动质子量 `TOTH_mob = TOTH − T_gel n̄_gel(pH)`（§2.3.1，pH 在 L 步开始时求一次），用逐层常数 `D_toth,ℓ = (Σ_{f 可动} D_f,ℓ β̄_f,ℓ + D_H,ℓ β̄_H,ℓ + D_OH,ℓ β̄_OH,ℓ) / β̄_mob,ℓ`（β̄ 为该层在本块内单元上的平均，β_mob 不含明胶），然后 `TOTH ← TOTH + (TOTH_mob' − TOTH_mob)`。明胶结合的质子因此不横向移动，TOTH 守恒，变换个数不变。单一 D 忽略了各族 D 的差别与总量梯度对质子通量的贡献（`assumed`；横向耦合弱，§3.2）。它的误差由可选的参考档对照 `test_lateral_toth_single_d`（D1c，slow）检查：同一边缘条带改为逐族扩散 P_f = T_f n̄_f 与 [H⁺]、[OH⁻]（各用本族 D），\|ΔD\| ≤ 2e-3。
- 核 G、H、Q、Λ_z、D_z 不横向扩散。

### 3.6 规范块与逐位一致

两个带横向作用的计算都在固定的规范块上进行，块网格锚定在画幅原点，与 `RunContext.tile_px`、进程数、调度顺序都无关：

| 计算 | 块（内） | 晕边 | 块外的填充 |
|------|----------|------|------------|
| 光学卷积（曝光模型内） | 512×512 像素 | `ceil(max(R_H, R_P)/p)` 像素 | 0（片门外无光） |
| 显影（显影步骤内） | B_c×B_c 个化学单元，B_c = 256 | `H_c = ceil(halo_um / Δx_c)` 个单元（Δx_c = 24 µm） | 未曝光的胶片（零曝光成员）：135 画幅之外是片边，物理上就是未曝光的同一种胶片 |

```
halo_um = min( max_k 3 · sqrt(2 · D_k · t_q,max / (1 − ρ_k)),  4.6 · δ_max )
ρ_k = exp(−t_a,min · D_k · (π / (2 δ_agit))²)      （没有周期搅动时 ρ_k = 0）
D_k = 液体中组分 k 的 D_w(T)（只做 T/η 换算，不乘 1/(1 + b_η I)，取上界）；k 也遍及 TOTH_mob 的横向 D（§3.5），按槽液初始成分的形态分布计算（D-76 约等于亚硫酸盐的 D；高 pH 显影液中趋向 D_OH = 4600，此时多由第二项封顶）
t_q,max = 最长的静置间隔（连续搅动时取显影时间，此时第二项起作用）
t_a,min = 最短的周期搅动时长
```

第一项是液体中横向扩散的 3σ。搅动只把 δ_agit 以上的液体复位，δ_agit 以下的子层、保护层与凝胶里的横向结构以速率 D_k(π/2δ_agit)² 弛豫，一次搅动后留下 ρ_k，所以各周期的横向方差按几何级数累加，乘 1/(1 − ρ_k)（`assumed`，估计式；由 `test_halo_convergence` 在先验上界处检验）。第二项来自厚度 δ_max、顶端狄利克雷的液层的横向屏蔽长度 2δ_max/π（衰减到 1e-3 需要约 4.4 δ_max，4.6 留出胶层厚度的余量），对所有阶段都是上界。δ_agit = 100 µm、δ_max = 300 µm、20 °C 时：小罐 30 s（5 s 搅动，Br⁻ 1816 um²/s，ρ = 0.11）1047 µm；大罐 60 s（7.5 s，ρ = 0.035）1380 µm（第二项）；静置 1380 µm。

逐位一致的论证：每个块的结果是（填充后的输入、槽液轨迹、参数、离散化）的纯函数，数组形状对所有块相同（边缘块也按完整形状填充），FFT/DCT 的尺寸固定，归约只在块内按固定形状进行，块之间的汇总（槽液通量、直方图）按块的行优先顺序进行。因此不同的 `tile_px`、不同的进程数得到逐位相同的结果（`test_expose_independent_of_tiling_and_workers`、`test_develop_independent_of_tiling_and_workers`）。物理计算里的 BLAS 归约不得依赖线程数：工作进程以单线程启动（环境变量 `OMP_NUM_THREADS = OPENBLAS_NUM_THREADS = MKL_NUM_THREADS = VECLIB_MAXIMUM_THREADS = 1`），热点代码不用 `np.linalg` 的批量求解（§12.8 也表明手写的向量化消元快 3 倍）。BLAS 一级归约的结果还依赖内存对齐（Accelerate 上 `np.dot` 对同一向量错开 8 字节，2000 次中 1267 次结果不同；`np.sum` 与默认 `einsum` 为 0 次），而工作进程收到的数组对齐取决于分配器，所以逐块物理代码里不得用 `np.dot`、`np.inner`、`np.vdot`、一维 `@`、`np.linalg.norm` 与 `einsum(optimize=True)`；改用逐元素乘加 `np.sum` 或 `einsum(optimize=False)`（`test_block_result_independent_of_alignment`：同一块的输入以故意错开 8 字节的副本重算，结果逐位相同）。

逐位一致不等于与"无限晕边"的解相等：有限晕边是降阶模型的一部分，误差由 `test_halo_convergence`（晕边 × 1.5 时显影银的变化 ≤ 3e-3 g/m²，即 2e-3 D，§3.8）界定。

像素的归属（显影块，Δx_c 与 p 不必整除）：块的 `inner`（单元下标 [i₀, i₁) × [j₀, j₁)）拥有中心落在 `[i₀Δx_c, i₁Δx_c) × [j₀Δx_c, j₁Δx_c)` 内的像素，即像素行 `ceil(i₀Δx_c/p − ½) … ceil(i₁Δx_c/p − ½) − 1`（列同理）；该块写这些像素并以这两个像素切片调用 `mark`，所以每个像素恰好属于一个块。块的成员窗口取与填充后单元区域有重叠的全部像素。画幅末端不足一个单元时（120 的 56 mm 为 2333.3 个单元），单元在画幅外的面积按未曝光片计入 p̄。

逐点步骤（停显、定影、水洗、干燥）沿用 D0 的分块，逐位一致由 D0 的机制保证。

### 3.7 药液槽对全幅的耦合

槽的耗尽只依赖全幅的平均源项。平均源项主要由曝光量的直方图决定，邻接效应对平均值只有二阶影响。所以先用直方图上的柱集合算槽液轨迹，再让所有块在这条轨迹下独立计算，最后检查。

槽中胶片的构成（`assumed`）：总面积 `A_film`（流程 `tank.film_area_cm2`）中，比例 `f_img`（`tank.image_area_fraction`，默认 135：(36×24)/(38×35) = 0.65，120：0.80）是与本帧同分布的画面，其余是未曝光的片边与片间。

1. **直方图**（读曝光字段，按规范块的顺序累加）：键 `s(x) = log10 m̄(x)`，`m̄ = (1/n_z) Σ_z Σ_c w_c (V_c e_V,z + S_c e_S,z)`。在本帧 s 的范围上取 64 个等宽箱，另加一个零曝光箱（m̄ < 1e-6）。每箱记像素数（int64，精确）与曝光描述量之和。
2. **预计算**：65 根柱（图像档离散化），权重为 `f_img × 箱内像素数 / 总像素数`，零曝光箱另加 `1 − f_img`，与槽精确耦合，得到槽液轨迹 `c_bulk(t_n)` 与平均通量。
3. **块计算**：所有块用同一条 `c_bulk(t)`，各自累加内单元的通量 `J̄_blk(t_n)`。
4. **检查（按输出判定）**：用 `f_img J̄_blk + (1 − f_img) J̄_zero` 重新积分槽 ODE 得 `c'_bulk`。把第 2 步的 65 根柱在给定轨迹 `c'_bulk` 下重算一次（不耦合槽，约 1 s），与它们在 `c_bulk` 下的结果比较各箱的显影银总量 `M_b = Σ_c silver_c + silver_physical`（g/m²）：

```
δ = max_b |M_b(c'_bulk) − M_b(c_bulk)| / (M_b(c_bulk) + 0.3 g/m²)
通过：δ ≤ ε_tank = 1e-3
```

   ΔD ≈ CP·ΔM，CP ≈ 0.6–0.8 m²/g，所以 ε_tank 对应 D ≤ 2 时 |ΔD| ≲ 2e-3（与 `test_tank_prepass_check` 一致；0.3 g/m² 的下限使近零曝光箱不按相对差判）。不通过时以 `c'_bulk` 为新轨迹重算各块，各块的输出原地覆盖（`FrameAssembler.reset()`，§4.4），最多再算 2 遍；仍不通过则警告。`meta["develop"]["prepass"]` 记录每遍的 δ、估计的收缩率 `L̂ = δ_{i+1}/δ_i`（有重算时）与按耗尽量归一的槽液变化 `max_{k,n} |c'_k − c_k| / s_k`（`s_k = max(|c_k(0)|, max_n |c_k(t_n) − c_k(0)|, 1e-6 M)`；bromide、iodide、silver 的初值为 0，按初值归一没有意义）以及槽液 pH 的最大变化，后两项只作诊断。

误差界：离真正不动点的输出误差约为 δ/(1 − L)。L 是映射 `c_bulk ↦ c'_bulk` 在输出上的收缩率，含组分之间的交叉项与 pH 通道（`∂ln R/∂pH ≈ ln10·|z|(1 − θ)`，一次性小罐的槽液 pH 降约 0.06–0.12），不作先验断言；预期 L ≈ 0.1–0.5，L̂ > 0.5 时警告。第一遍的失配来自直方图分箱与邻接效应对平均通量的二阶影响，预期通常一遍通过；不通过时每多一遍的代价是一次完整的块计算（§3.10）。感光条带不用预计算：条带的所有柱在同一个集合里与槽精确耦合（与图像预计算是同一段代码）。

### 3.8 图像路径的降阶模型

完整模型（参考档）逐像素做全深度化学，24 MP 需要约 2e14 flop 的每步计算、每像素约 150 个状态（3.6e9 个 float64 ≈ 29 GB），不可行。图像档是同一组方程的较粗离散化，降阶只有以下几处，每处都有界定误差的测试：

| 降阶 | 做法 | 由完整模型怎样得到 | 界定误差的测试 |
|------|------|--------------------|----------------|
| 化学单元共享 | 化学单元固定为 Δx_c = 24 µm（与 p 无关，网格锚定在画幅原点）。每个像素（成员）按与单元的面积重叠 w_m,cell 聚合进单元的 p̄（§2.4.4；p > 24 µm 时一个像素分到几个单元，p < 24 µm 时几个像素合进一个单元）；显影期间单元只受自己的 p̄ 驱动（最近单元耦合）。显影结束与各快照时，像素的 X、Y、P、F 由四个相邻单元的核在像素中心双线性插值后与它自己的 p(n) 合成 | 完整模型令 Δx_c → p 且成员只属于自己的单元即还原；凝胶交换长度 ≈ 15 µm（§2.4.12），比它细的细节本来就处在邻域平均的化学环境中，这正是高频邻接效应的饱和极限。插值输出使全幅银量与单元化学消耗的银量相差一个二阶小量（梯度²·Δx_c²），均匀区为 0 | `test_image_profile_edge`、`test_image_profile_mtf`、`test_resolution_independence` |
| 深度粗化 | n_z = 2，n_l = 4 | 同一有限体积格式，单元变粗 | `test_image_profile_uniform_patch` |
| 固定时间表 | (0.1 s, 1.6, 5 s) | 同一 ROS2，步长表变粗 | 同上，及 `test_fixed_schedule_vs_reference` |
| 有限晕边的规范块 | §3.6 | 晕边 → ∞ 即还原 | `test_halo_convergence` |
| 药液槽预计算 | §3.7 | 不动点迭代收敛即还原 | `test_tank_prepass_check` |

阈值（预注册，D1c/D1d 测试）。D1c 还没有 Nutting 与遮盖力（D1d），所以 D1c 的测试对显影银总量 `M = Σ_c silver_c + silver_physical`（g/m²）断言，密度阈值按 `ΔM = ΔD / CP_nom`、`CP_nom = 0.7 m²/g` 换算（`assumed`，只用于换算阈值）；D1d 起的测试直接对密度断言。

| 测试 | 目标 | 阈值 |
|------|------|------|
| `test_image_profile_uniform_patch`（D1c） | 34 级曝光的条带（`develop_strip`，曝光由 `uniform_exposure` 给出），图像档 vs 参考档 | max \|ΔM\| ≤ 0.014 g/m²，RMS ≤ 0.007 g/m²（即 0.01、0.005 D） |
| `test_image_profile_edge` | 阶跃边 Δlog H = 1.5 与 0.6（1024×64 像素，p = 6 µm）vs 参考条带（Δx = 2 µm，再按像素盒平均） | 各位置 \|ΔD\| ≤ 0.02；边缘过冲幅度相差 ≤ 20 % |
| `test_image_profile_mtf` | 正弦目标 2.5、5、10、20、30、50 c/mm | ≤ 30 c/mm 时 \|ΔMTF\| ≤ 3 个百分点，50 c/mm 时 ≤ 5 个 |
| `test_halo_convergence`（D1c） | 晕边 × 1.5；小罐、大罐、静置各一次，δ_agit、δ_max 取典型值（50、300 µm）与先验上界（100、500 µm） | \|ΔM\| ≤ 3e-3 g/m²（即 2e-3 D） |
| `test_tank_prepass_check`（D1c） | 小罐 1+1、整卷负载；128×192 小帧上"所有块带实时槽"的全耦合解 | 槽液轨迹差按耗尽量 s_k 归一（§3.7）≤ 1e-2，槽液 pH 差 ≤ 1e-3；\|ΔM\| ≤ 3e-3 g/m² |
| `test_resolution_independence` | 同一场景 p = 6、9、12、30 µm（化学单元都是 24 µm；9 µm 是非整数比，30 µm 是像素大于单元的预览） | 120 µm 盒平均密度差 ≤ 0.005；48 µm 孔径颗粒度差 ≤ 5 %（p = 30 µm 只比盒平均密度） |

预注册的退路：`test_image_profile_mtf` 在 24 µm 下不通过时，把 Δx_c 改为 16 µm（化学成本 × 2.25），不改方程。

### 3.9 收敛测试

参考档相对更细的离散化（`slow` 标记）：

| 测试 | 变化 | 阈值（D-76 条带，D1c 用夹具 `test/bw-mechanistic` 的参数；D1c 对 M 断言，换算同 §3.8） |
|------|------|-----------------------------------|
| `test_convergence_depth` | n_z 8 → 16 | \|ΔM\| ≤ 1.4e-3 g/m²（1e-3 D） |
| `test_convergence_liquid` | n_l 16 → 32（两段各加倍） | \|ΔM\| ≤ 1.4e-3 g/m² |
| `test_convergence_time` | rtol 1e-4 → 1e-6 | \|ΔM\| ≤ 1.4e-3 g/m² |
| `test_convergence_classes` | K 12 → 24 | \|ΔM\| ≤ 1.4e-3 g/m² |
| `test_convergence_quadrature` | Q_GH 21 → 41 | \|ΔM\| ≤ 1.4e-4 g/m² |
| `test_calibration_profile_vs_reference` | 标定档 vs 参考档 | \|ΔM\| ≤ 7e-3 g/m²（0.005 D） |
| `test_calibration_profile_vs_reference_metrics`（D1e） | 同上，Tri-X 先验参数，经密度计与 §6 的指标 | \|ΔD\| ≤ 0.005，\|ΔCI\| ≤ 0.005，log H_m 差 ≤ 0.005 |
| `test_splitting_second_order` | 周期条带，固定步长减半 | 相邻差之比 3.5–4.5 |

### 3.10 计算代价

本机单进程实测的核函数耗时（§12.8）：

| 操作 | 耗时 |
|------|------|
| 8×8 稠密消元（向量化，N = 1e5） | 179 ns/柱 |
| `np.linalg.solve` 批量 8×8 | 551 ns/柱（不用） |
| 成员指数更新（K·n_z = 12） | 52 ns/成员（旧设计，成员已移出时间积分） |
| 双线性查表 ×12 | 299 ns/成员（旧设计；现每成员每次显影只查两次潜影表） |
| pH 牛顿（3 次迭代，7 个单元，D-76 各族） | 2.66 µs/柱·f 求值（评审脚本 `column.py`） |
| 两个凝胶单元的 10×10 消元与回代 | 0.51 µs/柱（同上） |
| 柱步多进程吞吐（内存带宽受限的成员核） | 10 进程 ×6.2，14 进程 ×4.8（`scaling2.py`） |
| 布尔方差求积 Var_bool（48 点） | 2.9 µs/像素（`grain_bench.py`） |
| DCT 横向半步（56 个场，340²） | 9.0 ns/场·单元 |
| FFT 卷积（6 通道，896² 块，512² 内块） | 151 ns/内像素 |
| `rng.pixel_uniform` | 4.0 ns/样本 |
| `ndtri` | 8.4 ns/样本 |

**每个条带**（标定档：13 个单元、约 8 个组分、每柱 360 个核，35–40 根柱一次调用）：每步约 3 ms（每步约 200 次 NumPy 调用约 1 ms，加上 13 个单元的 pH 求解与凝胶消元：评审实测两者在 35 柱时已是 1.13 ms/步，350 柱时 6.0 ms/步，所以按条件拼批只省约 2 倍，不做跨条件批处理，每个条件一次 `run_ensemble` 调用）。步数见 §3.4：系列运行 248–503 步，CI(t) 运行 492–623 步，平均约 450 步，**一个条件约 1.4 s**（估计，D1c 实测）。

**标定**（本机 14 核 = 10 性能核 + 4 能效核；10 个工作进程，有效并行按 ×8 计，`assumed`：条带调用的工作集只有几 MB，受 Python 开销与计算限制，不像旧成员核那样受内存带宽限制（后者 10 进程只有 ×6.2）；若实测只有 ×6，下表各项乘 1.3，合计约 3 天。D1c 实测后改写）：

| 项 | 规模 | 时间 |
|----|------|------|
| 阶段 D 的一次目标函数 | 约 130 个条件 × 1.4 s / 8 | 约 23 s |
| 阶段 D 的一次雅可比 | 约 36 列；每列只重算读到该参数的条件（§7.4），约 2500 次条件运行，为全部重算的一半 | 约 7 min |
| 阶段 D 的一次 TRF 迭代 | 雅可比 + 1–2 次求值 | 约 8 min |
| 阶段 D 多起点 | 8 个起点各 10 次迭代，χ² 最好的 2 个再各 20 次（§7.4） | 约 16 h |
| 阶段 C 的一次求值（每种胶片） | 2 个系列运行、CI(t) 运行、MTF 小信号柱对（§6.7，1 次调用）、倒易律 6 行（6 张潜影表 × 0.42 s + 定点求曝光）、光谱 1 张表 + 192 柱、颗粒度 1 柱（闭式，§6.6） | 单核约 25 s |
| 阶段 C 拟合（每种胶片） | 约 13 列雅可比并行；30 次迭代；多起点同上 | 一次拟合约 25 min，多起点约 1.5 h |
| 阶段 E | ≤ 5 轮，每轮 C 与 D 从上一轮热启动各约 10 次迭代 | ≤ 约 10 h |
| 剖面似然（§7.5） | 化学块 4 个参数 × 8 个格点 × 约 5 次热启动迭代（约 21 h）；胶片块每种 6 × 8 × 5 次（约 3 h） | 约 27 h |

合计约 2–2.5 天（10 进程），一次性运行。参考档的潜影表（12 × 401 × 21 次 36×36 矩阵指数，约 3.9 s）只在验证中出现，不进入拟合。

**24 MP 全幅**（6000×4000，p = 6 µm，图像档，D-76 小罐 30 s/5 s，显影 7.5 min，209 步，§3.4）：

| 部分 | 规模 | 单核时间 |
|------|------|----------|
| 光学卷积 | 2400 万像素 | 约 5 s |
| 潜影表 | K × 201 × 9 次 36×36 矩阵指数（F 不建表） | 约 0.4 s |
| 化学单元（含晕边 44 格（1047 µm，§3.6），开销 (256 + 88)²/256² = 1.81） | 1.5e6 × 1.81 × 209 步 × 约 12 µs/单元·步（两个阶段：7 个单元的 pH 求解 2 × 2.66 µs、凝胶消元 2 × 0.51 µs，加反应、雅可比、输运与 180 个核，估计） | 约 6800 s |
| 成员 | 开始时聚合 p̄（2.4e7 × 1.81 × 约 0.5 µs），结束时合成（2.4e7 × 约 1.5 µs）；不进入时间循环 | 约 60 s |
| 横向 DCT | 56 个场 × 2.7e6 × 9 ns × 2 × 209 | 约 570 s |
| 预计算与检查（2 × 65 根柱） | — | 约 3 s |
| 颗粒抽样与读出 | Var_bool 2.4e7 × 2.9 µs，加均匀数、泊松逆变换与读出 | 约 90 s |
| 合计 | | 约 7500 s |

块之间完全独立，按进程并行。默认 6 GB 预算下 9 个工作进程（每个约 0.5 GB，§4.6），有效并行按 ×6 计（`assumed`；以 1.6e4 单元为一批的柱步工作区约 100 MB，超出缓存，与评审实测的内存带宽受限核同类，后者 10 进程 ×6.2），约 **21 min**。同一条件显影 12 min（335 步）约 33 min，大罐 8 min（160 步）约 16 min。符合"几分钟而不是几个小时"，但不是旧版写的 5–10 min。

预览不比全幅便宜多少：化学单元固定为 24 µm，单元数只由画幅面积决定（135 为 1500 × 1000 = 1.5e6），1 MP 预览的显影与全幅一样约 20 min（D-76 小罐 7.5 min），只有光学、成员与颗粒随像素数减少。是否另设一个更粗的、未经验证的预览离散化，留给所有者决定。

D1d 的 `test_benchmark_image_path`（`slow`）：上述条件下（a）单进程算 1 个规范块，（b）以默认进程数同时算 `workers` 个块，实测并行效率 `η = t_a · workers / t_b`；外推 24 MP（24 个块）的时间 `t_b · ceil(24 / workers)` 加上其余各项的实测值，超过 30 min 即失败；单元·步耗时、η 与外推值写入测试日志，并据此改写本节。

---

## 4 内存设计

### 4.1 问题

24 MP（6000×4000）下 D0 的做法：`scene_light` 生成 `(H, W, 81)` float64 光子通量，15.6 GB；`EmulsionState` 的所有字段常驻内存；分块步骤每拼一块就复制整个字段（`EmulsionState.paste` 的 `current.copy()`），总代价随块数线性增长、随画幅平方增长。D1 的显影输出每像素 21 个 float64（K = 6），约 4.0 GB，再加曝光字段 0.8 GB。

### 4.2 低秩光

RGB 场景经 Mallett 基还原是线性的：`Φ_b(x) = Σ_r C_r(x) Q_r,b`（§2.1.1）。`LowRankLight` 只存系数 `(H, W, 3)`（即线性 RGB 本身，576 MB，可与输入数组共享）与基 `(3, B)`，从不生成 `(H, W, B)`。曝光模型只需要投影 `light.project(W)`（`(H, W, M)`，M 为权向量个数）与窗口 `light.window(rows, cols)`（按块生成 `(h, w, B)` 的小块，只给全谱输入与理想模型使用）。全谱输入（光谱图像）仍用 D0 的 `Light`，按块读取。

### 4.3 字段存储

新模块 `fieldstore.py`：

- `FieldStore(root, ram_threshold_bytes=256 MiB)`：大小超过阈值、且 `root` 不为 None 的字段分配为 `np.memmap(mode="w+")`，写完后 `freeze` 重新以 `mode="r"` 打开（只读映射），交给 `EmulsionState`。小字段留在内存里。`root` 是运行的临时目录（runner 在输出目录下建 `.atlas-scratch-<随机后缀>/`，结束时删除，`--keep-scratch` 保留）。
- `_frozen.readonly_array` 修改：`mode == "r"` 的 `np.memmap` 及其视图视为只读、直接共享不复制（D0 现在会把内存映射复制进内存，因为可写映射可能被别处改写；只读映射没有这个问题）。
- 均匀字段（`residual_dye`、`fixer_residual`、`presoaked`）用 1×1 只读数组的 `np.broadcast_to` 视图，零字节。
- 存储内的名字不是字段名：步骤用 `name = f"s{step_index}:{key}"` 分配，重复的 name 抛 `ValueError`；映射文件名为 `<name 中的 / 与 : 换成 _>-<uuid4 前 12 位>.bin`，每次分配唯一。显影与停显都写 `emulsion/silver`、`emulsion/silver_m2`，按字段名命名文件会以 `"w+"` 截断输入状态仍在只读映射的文件。`freeze`、`release`、`descriptor` 都按 name。
- 跨进程只传描述符，不传数组：`FieldStore.descriptor(name) -> FieldDescriptor(path, dtype, shape, offset)`；工作进程用 `open_field(desc, mode)` 自己打开（输入 `"r"`，预分配的输出 `"r+"`）。把内存映射的视图直接传给工作进程会被序列化成完整内容（4 MB 视图序列化为 4 000 163 字节，已核对），所以 `run_blocks` 的 `fn` 只接受描述符与小参数。内存里的小字段（低于阈值）照常序列化。
- 一个运行只有一个 `FieldStore`：由 runner 创建，挂在 `RunContext.store` 上（`field(compare=False, hash=False, repr=False)`，不进入运行记录与 `recipe_hash`），步骤经 `StepContext.run.store` 取用；store 负责 `release`，runner 负责 `close`。没有 store（D0 路径、测试）时字段留在内存。

### 4.4 拼接器取代重复 paste

`state.FrameAssembler`：按 `FieldSpec(key, trailing, dtype)` 为每个输出字段一次性分配（经 `FieldStore`），`write(key, inner, values)` 直接写入内块区域；记录每个块是否已写，`finish()` 检查每个像素恰好写过一次，冻结并返回只读字段。总代价 O(像素数)。`pipeline._run_step` 改用它；`EmulsionState.paste` 保留为公开 API（文档注明每次调用 O(N)），流水线不再调用。

并行时输出不经主进程：工作进程按 `FrameAssembler.descriptors()` 给出的描述符以 `"r+"` 打开输出字段，把自己的内块直接写进去，只返回小的摘要（内块区域、块内通量 `J̄_blk(t_n)`、计数）；主进程用 `mark(inner)` 登记已写区域，`finish()` 照样检查每个像素恰好写过一次。`run_blocks` 按块的顺序返回摘要（列表，每项为 KB 量级），从不收集块的输出数组（24 块 × 每块 176 MB 会在主进程堆出约 4.2 GB）。`reset()` 清空登记、保留已分配的输出，供药液槽检查不通过时整幅重算一遍（§3.7）原地覆盖。

一个块的某字段若是广播视图，且所有块给出同一个标量，拼接结果仍是整幅的广播视图，不分配内存；否则按普通字段分配。

### 4.5 分段融合与流式执行

`pipeline.process` 把相邻的、有限作用范围的步骤合成一段（segment）。一段内对每个块：从输入状态裁出窗口（内存映射的视图，不复制），依次执行段内所有步骤，只把段的最终字段写进拼接器。段的晕边是各步作用范围之和（D1 的逐点步骤都是 0）。全幅耦合的步骤（显影）单独成段，拿到整幅状态（字段是内存映射），内部用规范块。

`run_image` 的流式顺序与每步驻留的内容：

| 阶段 | 内存中 | 内存映射文件 | 释放 |
|------|--------|--------------|------|
| 读图 | 线性 RGB float64（576 MB） | — | — |
| 曝光 | RGB（作低秩光的系数，不复制）+ 光学工作块 | `exposure_volume`、`exposure_surface`（768 MB） | 曝光后丢弃光与 RGB |
| 预湿、显影 | 直方图、槽液轨迹、工作进程（每个约 0.5 GB） | 工作进程读曝光字段、直接写显影输出（4.0 GB） | 显影后删除曝光字段文件 |
| 停显…干燥（融合段） | 当前块 | 写更新后的 `silver`、`silver_m2`、`halide` 等 | 停显后删除 `growth_rate` |
| 读出 | 输出 `(H, W, 3)`（576 MB）+ 按行块的临时量 | 读显影输出 | — |
| 写 PNG | uint16（144 MB） | — | 结束时删除临时目录 |

### 4.6 24 MP 的峰值内存

| 项 | 大小 | 位置 |
|----|------|------|
| 线性 RGB / 低秩光 | 0.58 GB | 内存，曝光后释放 |
| 曝光字段（2·n_z = 4 个 float64/像素） | 0.77 GB | 内存映射 |
| 显影输出（3K + 3 = 21 个 float64/像素） | 4.0 GB | 内存映射 |
| 读出结果 | 0.58 GB | 内存 |
| 显影工作进程，每个（256² 单元 + 44 格晕边 → 344² = 1.18e5 单元） | 单元状态（56 个浓度、Λ、D、pH）约 65 MB；核 G、H、Q（180 个/单元）170 MB；柱步以 1.6e4 单元为一批、连续做完两个阶段（最近单元耦合下没有全块的阶段屏障，§3.3）的工作区约 200 MB；成员不常驻：p̄ 的聚合与结束时的合成按行批次从曝光字段现算（每批约 30 MB）；合计约 0.5 GB | 内存 |
| 主进程（直方图、槽液轨迹、调度、块摘要） | < 0.3 GB | 内存 |

旧版每进程 0.70 GB 的估计漏了 ROS2 的阶段存储（CIC 耦合要求全块先做完阶段 1，成员需要 y、k₁、k₂，约 1.09 GB）；成员移出时间积分后不再需要。

`RunContext.memory_budget_gb`（默认 6）决定显影的工作进程数：`workers = min(cpu_count, max(1, floor((budget − 1.5)/0.5)))`，默认得 9 个；用户给定 `workers` 时以较小者为准并在超预算时警告。按此，24 MP 的常驻内存约 5–6 GB，其余由操作系统的页缓存按需换入换出。

### 4.7 接口改动

| 文件 | 改动 |
|------|------|
| `light.py` | 新增 `LowRankLight`；`Light` 新增 `window(rows, cols)`、`project(weights, rows=None, cols=None)`；协议 `LightField` |
| `scene/camera.py` | `scene_light(..., low_rank: bool = False)`；`True` 时返回 `LowRankLight` |
| `spectral.py` | 光源名增加 `"D55"`（colour-science 的 CIE D55，供 ISO 感光测定日光） |
| `_frozen.py` | `readonly_array` 共享 `mode="r"` 的 `np.memmap` |
| `fieldstore.py`（新） | `FieldStore`（唯一文件名、分配句柄、`descriptor`）、`FieldDescriptor`、`open_field` |
| `blocks.py`（新） | 规范块网格与并行执行；`run_blocks` 只回传小摘要 |
| `state.py` | 新增 `FieldSpec`、`FrameAssembler`（含 `descriptors`、`mark`、`reset`）、`EmulsionState.drop(*keys)` |
| `pipeline.py` | `RunContext` 新字段 `film_format`、`profile`、`workers`、`memory_budget_gb`、`scratch_dir`、`store`（不参与比较与哈希）；`process` 改为分段融合加拼接器；`expose` 接受 `LightField` |
| `runner.py` | `run_image(..., film_format=None, profile="image", workers=None, memory_budget_gb=6.0, scratch_dir=None, keep_scratch=False)`；低秩光；临时目录的生命周期；运行记录新字段 |
| `manifest.py` | 新字段 `film_format`、`model_settings`（JsonDict：档位与离散化）、`kinetics_file`（FileRef \| None）、`calibration_runs`（list[str]）、`workers`；`recipe_hash` 排除 `workers` 与临时目录路径 |
| `cli.py` | `atlas run` 新选项 `--format`、`--profile`、`--workers`、`--memory-gb`、`--scratch`、`--keep-scratch`；新命令 `atlas validate`、`atlas calibrate`（D1e）、`atlas micro`（D1f） |

函数签名见 §10。

---

## 5 数据模式

所有新字段都可选，D0 的 schema_version 仍为 1，现有的理想材料与测试夹具无需改动；`bw-mechanistic` 模型在曝光与显影开始前检查它需要的字段，缺少时报出文件路径与字段路径（D0 `DataError` 的格式）。

### 5.1 Quantity 的扩展与新单位

```python
class FitRecord(FrozenModel):
    run_id: Text                                   # data/calibration/<run_id>.yaml
    model_version: Text                            # 例如 "D1e"
    training_set_sha256: Sha256                    # 训练数据文件的联合哈希
    method: Literal["covariance", "profile", "prior"]
    transform: Literal["linear", "log10", "logit"]
    sigma_transformed: NonNegativeFloat            # 变换空间中的一个标准差
    interval: tuple[float, float] | None = None    # 95 % 区间（原单位）；剖面法给出时必填
    at_bound: bool = False                         # least_squares 的 active_mask
    identifiable: bool = True                      # §7.5 的判定

class Quantity(FrozenModel):                       # D0 字段不变，新增：
    bounds: tuple[float, float] | None = None      # 先验盒，原单位；lo ≤ value ≤ hi
    prior: Literal["uniform", "normal", "lognormal"] | None = None
    fit: FitRecord | None = None
    # 校验：source == "fitted" 当且仅当 fit 不为 None；bounds 存在时 value 在界内
```

`units.py` 新增（ASCII，换算因子为精确值）：

| 量纲 | 单位 |
|------|------|
| 浓度 | `mol/L`、`mmol/L` |
| 比含量 | `mol/g` |
| 二级速率常数、吸附常数 | `L/mol`、`L/(mol*s)`、`kg/mol` |
| 扩散系数 | `um2/s`、`cm2/s`（1 cm²/s = 1e8 um²/s） |
| 速度 | `um/s` |
| 一级速率 | `1/s` |
| 摩尔能量 | `kJ/mol`、`J/mol`、`eV`（1 eV ≙ 96.4853 kJ/mol） |
| 倒长度 | `um-1`、`cm-1` |
| 面积、体积 | `um2`、`um3` |
| 遮盖力 | `1/(g/m2)` |
| 黏度 | `cP` |

参考曲线的轴单位（`log10(lux*s)`、`cycles/mm`、`%`、`log10(cm2/erg)`、`min`）只作标签，不换算。

### 5.2 物种与组分（`data/chemistry/species.yaml`）

`Species` 新增字段：

```python
class AcidBase(FrozenModel):
    charge_most_protonated: int
    reference_protons: int                 # 参考态比最质子化态少的质子数（§2.3.1）
    pka: list[Dimensionless]               # 25 °C、I → 0
    temperature_law: Literal["van_t_hoff", "manov_boric", "linear"] = "van_t_hoff"
    ionisation_enthalpy: list[Quantity] = []   # kJ/mol，与 pka 一一对应；缺省 0
    dpka_dt: list[Quantity] = []               # 1/K，temperature_law == "linear" 时用
    activity_class: Literal["sulfite", "borate", "organic", "amine", "none"]

class AgentRole(FrozenModel):
    role: Literal["eta", "regenerator", "single"]
    kinetic_class: NameId                  # kinetics.yaml 中 agents 的键
    active_state: int                      # 活性形态在 AcidBase 各态中的下标
    product_sulfonated: bool

class Species(FrozenModel):                # D0 字段不变，新增：
    kind: Literal["ingredient", "component"] = "ingredient"
    dissolves_to: dict[NameId, PositiveFloat] = {}    # 配料：组分 id → mol/mol
    toth: float = 0.0                                 # 配料：对 TOTH 的贡献，mol/mol
    cations: dict[Literal["Na", "K"], NonNegativeFloat] = {}
    spectator_charge: dict[str, float] = {}           # 配料：旁观离子（如 SO4: -2）的 mol/mol，供离子强度
    acid_base: AcidBase | None = None                 # 组分
    gel_diffusion_ratio: Dimensionless | None = None  # 组分
    developing_agent: AgentRole | None = None         # 组分
    immobile: bool = False                            # 组分；明胶与 toth 是内置伪组分，不写进 species.yaml（§2.3.1）
```

校验：组分必须有 `diffusion_coefficient`（um2/s，20 °C）与 `gel_diffusion_ratio`（`immobile` 除外）；配料的 `dissolves_to` 只能指向组分；`Library.validate()` 增加这两条检查。组分也可以直接写在药液的 `components` 里（如 `hydroquinone`），等价于一个 `dissolves_to: {自身: 1}`、`toth: 0`、无阳离子的配料；其 `molar_mass` 是参考态的摩尔质量。

### 5.3 共用动力学常数（`data/chemistry/kinetics.yaml`）

```python
class AgentKinetics(FrozenModel):
    k_et: Quantity                 # 1/s
    K_ads: Quantity                # L/mol
    lambda0: Quantity              # 1/s
    activation_energy: Quantity    # kJ/mol
    k_des: Quantity | None = None  # 1/s（ETA）
    k_reg: Quantity | None = None  # L/(mol*s)（再生剂）
    pka: Quantity | None = None    # 动力学类共用的酸 pKa（吡唑烷酮）

class KineticsFile(FrozenModel):
    schema_version: Literal[1]
    model_version: Text
    activity: ActivityParams       # debye_huckel_A_20C, debye_huckel_A_25C, sit_b, delta_eps: dict[class, Quantity]
    agents: dict[NameId, AgentKinetics]
    halides: HalideParams          # K_Br, iodide_ratio, adsorption_enthalpy
    barriers: BarrierParams        # u_gel, u_halide, ionic_strength_ref
    latent: LatentInteraction      # u_li, n_sat
    solvent: SolventParams         # latent_attrition, k_diss, c_sat_ref, sulfite_exponent, bromide_exponent, ligands, k_pd, sulfonation_switch
    transport: TransportParams     # delta_max, delta_agit, swelling_half_time, viscosity_ionic_coefficient, water_viscosity (4 点表), proton_diffusion, hydroxide_diffusion
    gelatin: GelatinParams         # groups_per_gram, pka, initial_ph, dry_density
    calibration: CalibrationRef | None = None    # run_id 与 sha256
```

字段名与 §2.3.6、§2.4.11 的参数表一一对应；每个值都是带来源的 `Quantity`。

### 5.4 胶片

```python
class Grain(FrozenModel):          # D0 字段不变，新增：
    thickness: Length | None = None            # tabular 必填
    iodide_fraction: Dimensionless | None = None
    n_classes: int = Field(default=6, ge=2, le=64)

class DyeBand(FrozenModel):
    amplitude: Dimensionless
    peak_wavenumber: Quantity      # um-1
    width: Quantity                # um-1

class EmulsionOptics(FrozenModel):
    gelatin_index: Dimensionless
    intrinsic_absorption: IntrinsicAbsorption    # alpha0 (um-1), edge_energy (eV), urbach_energy (eV)
    dye_bands: list[DyeBand] = Field(min_length=1, max_length=4)
    scattering_550: Quantity       # um-1
    scattering_exponent: Dimensionless
    psf_weight: Dimensionless
    psf_length: Length

class LatentParams(FrozenModel):
    efficiency: Dimensionless
    efficiency_size_exponent: Dimensionless
    tau0_median: Quantity          # s
    tau0_sigma_ln: Dimensionless
    ionic_time: Quantity           # s
    competing_sites: Dimensionless
    hole_attack: Dimensionless
    cluster_cap: int = 5

class EmulsionDevelopment(FrozenModel):
    fog_hazard_ratio: Dimensionless
    covering_power: Quantity                   # 1/(g/m2)
    covering_power_size_exponent: Dimensionless
    physical_covering_ratio: Dimensionless
    halide_density: Quantity                   # 1/(g/m2)
    dye_density: Dimensionless
    fixing_factor: Dimensionless

class EmulsionLayer(FrozenModel):  # D0 字段不变，新增：
    gelatin_coverage: ArealMass | None = None
    swelling_ratio: Dimensionless | None = None
    optics: EmulsionOptics | None = None
    latent: LatentParams | None = None
    development: EmulsionDevelopment | None = None

class BaseFormat(FrozenModel):
    thickness: Length
    grey_density: Dimensionless
    backing_density: Dimensionless

class BaseLayer(FrozenModel):      # D0 字段不变，新增：
    formats: dict[NameId, BaseFormat] = {}

class Material(FrozenModel):       # D0 字段不变，新增：
    default_format: NameId | None = None       # formats 非空时必填，且须为其中之一
```

画幅 id 用 `135`、`120`（`NameId` 允许以数字开头）。`EmulsionLayer.spectral_sensitivity` 沿用 D0 字段，指向参考曲线文件，只在标定中使用。

### 5.5 药液

```python
class PhAnchor(FrozenModel):
    dilution: Dilution
    temperature: Temperature
    ph: Dimensionless              # 带 uncertainty
    method: Literal["glass_electrode", "colorimetric", "typical"]
    note: str | None = None

class BathConstants(FrozenModel):
    neutralisation_time: Quantity | None = None    # s（停显）
    clearing_time_constant: Quantity | None = None # s（定影）
    dye_removal_time: Quantity | None = None       # s（定影）
    wash_time_constant: Quantity | None = None     # s（水洗）
    wash_dye_time_constant: Quantity | None = None # s（水洗）

class Bath(FrozenModel):           # D0 字段不变，新增：
    variant: Text | None = None                    # "packaged"、"published"、"pre2019"、"2019"、"generation-I"
    base_solution: Literal["stock", "concentrate"] = "stock"
    ph_anchors: list[PhAnchor] = []
    nuisance: list[NameId] = []                    # 在 bounds 内拟合的配料 id
    constants: BathConstants | None = None
```

`components[*].amount` 的 `bounds` 记录 SDS 范围（msds-estimated）或称量公差；`nuisance` 中的配料必须有 `bounds`。HC-110、T-MAX Developer 与 Rodinal 的组分写为浓缩液的 g/L（DEV §3.3：B 稀释 = 浓缩液 1+31），流程里写 `dilution: "1+31"`、`"1+4"`、`"1+25"`。

### 5.6 流程

```python
class Agitation(FrozenModel):      # method 增加 "lift"（大罐提起–倾斜–沥干）
    method: Literal["inversion", "lift", "rotary", "stand", "continuous", "tray"]
    ...                            # initial_s、interval_s、duration_s 含义不变（§2.4.8）

class Tank(FrozenModel):           # D0 字段不变，新增：
    image_area_fraction: Dimensionless | None = None   # §3.7；缺省按画幅
    rolls_before: int = Field(default=0, ge=0)         # 此前已在同一药液中显影的卷数（容量检查）
```

`rolls_before > 0` 时槽液初值取"此前各卷按同一直方图显影后"的成分（用预计算的柱集合依次推进），只用于容量一致性检查（§8.1）。

### 5.7 标定记录（`data/calibration/<run_id>.yaml`）

```python
class StageRecord(FrozenModel):
    name: Literal["A_activity", "B_spectral", "C_film", "D_chemistry", "E_polish", "F_identifiability"]
    parameters: list[Text]
    data_files: list[FileRef]                      # 路径相对数据根，带 sha256
    table_rows: list[Text]                         # 表格行的标识（§7.6）
    optimizer: JsonDict                            # 方法、损失、x_scale、起点数、种子
    chi2: NonNegativeFloat
    n_observations: int                            # 残差个数
    n_effective: NonNegativeFloat                  # §7.2 的独立数（s² 用它，§7.8）
    n_parameters: int
    n_penalised: int = 0                           # 取有界罚值的残差项数（§7.2），最终必须为 0
    converged: bool

class PairStatus(FrozenModel):                     # §8.2
    material_id: Text
    developer_id: Text
    tier: Literal["training", "held_out", "predicted"]
    status: Literal["validated", "failed", "predicted"]
    headline: float | None = None                  # §8.1 第 3 项的标题数字（仅留出）
    baseline_headline: float | None = None         # 朴素基线的同一数字
    failed_metrics: list[str] = []                 # 指标 id 与原因（含"算不出"）
    thresholds_sha256: Sha256
    after_unblinding: bool = False
    post_hoc_thresholds: bool = False

class ParameterResult(FrozenModel):
    name: Text                                     # 如 "films/kodak/tri-x-400:emulsion.latent.tau0_median"
    value: float
    unit: str
    transform: Literal["linear", "log10", "logit"]
    sigma_transformed: float
    interval: tuple[float, float] | None
    method: Literal["covariance", "profile", "prior"]
    conditional_on_other_block: bool = True        # 协方差与剖面以另一块的不动点为条件（§7.8）
    at_bound: bool
    identifiable: bool
    prior_dominated: bool = False                  # §7.5：σ_data > 0.5 τ_prior

class CalibrationRecord(FrozenModel):
    schema_version: Literal[1]
    run_id: Text
    created_utc: AwareDatetime
    model_version: Text
    git_commit: str | None
    profile: Literal["calibration"]
    stages: list[StageRecord]
    parameters: list[ParameterResult]
    singular_values: dict[str, list[float]]        # 每个参数块的加权雅可比奇异值
    training_set_sha256: Sha256
    thresholds: JsonDict = {}                      # §8.3 全表与组合–指标表，冻结时写入
    split_registry_sha256: Sha256 | None = None    # §8.4 第 3 步
    metrics_code_sha256: Sha256 | None = None      # §8.4 第 3 步
    frozen: bool = False                           # 冻结后为 True，之后才允许留出评估；为 True 时上面三项必须非空
    held_out_evaluated: bool = False
    pair_statuses: list[PairStatus] = []           # §8.2，图集读取
    notes: list[str] = []
```

写回规则（§7.8）：每个拟合值写回它所在的数据文件（胶片、`kinetics.yaml`、药液的干扰参数），`source: fitted`，`fit.run_id` 指向本记录。

### 5.8 示例文件

以下数值标注 `PLACEHOLDER` 的都是占位，D1e 标定之前不得当作结果引用；带出处的数值取自研究报告。

`data/chemistry/species.yaml`（节选）：

```yaml
schema_version: 1
species:
  - id: sodium-sulfite
    name: sodium sulfite, anhydrous
    formula: Na2SO3
    cas: 7757-83-7
    molar_mass: {value: 126.04, unit: g/mol, source: literature, ref: "PubChem; DEV 6.2"}
    roles: [preservative]
    kind: ingredient
    dissolves_to: {sulfite: 1}
    toth: 0
    cations: {Na: 2}
  - id: borax-decahydrate
    name: borax decahydrate
    formula: Na2B4O7.10H2O
    cas: 1303-96-4
    molar_mass: {value: 381.37, unit: g/mol, source: literature, ref: "PubChem; DEV 6.2"}
    roles: [alkali, buffer]
    kind: ingredient
    dissolves_to: {borate: 4}
    toth: -2
    cations: {Na: 2}
  - id: sulfite
    name: sulfite (S(IV) total)
    molar_mass: {value: 80.06, unit: g/mol, source: literature, ref: "SO3 2- formula mass"}
    roles: [preservative, buffer]
    kind: component
    acid_base:
      charge_most_protonated: 0
      reference_protons: 2
      pka:
        - {value: 1.86, unit: "1", source: literature, ref: "Goldberg & Parker, J. Res. NBS 90, 341 (1985); DEV 6.2"}
        - {value: 7.19, unit: "1", source: literature, ref: "Goldberg & Parker 1985; DEV 6.2"}
      ionisation_enthalpy:
        - {value: 0.0, unit: kJ/mol, source: assumed, note: "irrelevant above pH 4"}
        - {value: -3.65, unit: kJ/mol, source: literature, ref: "Goldberg & Parker 1985"}
      activity_class: sulfite
    diffusion_coefficient: {value: 1030, unit: um2/s, source: assumed, note: "Nernst-Hartley for Na2SO3 from recalled ionic D (DEV 6.3), scaled to 20 C"}
    gel_diffusion_ratio: {value: 0.5, unit: "1", source: assumed, bounds: [0.3, 0.8], note: "Iwano KBr 0.8 at 8 % and 0.55 at 20 % gelatin (KIN 3.2), extrapolated to the 30-45 wt % swollen emulsion"}
  - id: hydroquinone
    name: hydroquinone
    formula: C6H6O2
    cas: 123-31-9
    molar_mass: {value: 110.11, unit: g/mol, source: literature, ref: "PubChem; DEV 6.1"}
    roles: [developing_agent]
    kind: component
    acid_base:
      charge_most_protonated: 0
      reference_protons: 0
      pka:
        - {value: 9.85, unit: "1", source: literature, ref: "IUPAC set, Baxendale & Hardy 1953; DEV 6.1"}
        - {value: 11.45, unit: "1", source: literature, ref: "IUPAC set 11.39-11.56; DEV 6.1", uncertainty: 0.1}
      activity_class: organic
    diffusion_coefficient: {value: 855, unit: um2/s, source: assumed, note: "Wilke-Chang/Hayduk-Laudie mean at 20 C, DEV 6.3"}
    gel_diffusion_ratio: {value: 0.10, unit: "1", source: literature, ref: "Iwano 1970 Table 5; KIN 3.2"}
    developing_agent: {role: regenerator, kinetic_class: hydroquinone, active_state: 2, product_sulfonated: true}
```

`data/chemistry/kinetics.yaml`（节选）：

```yaml
schema_version: 1
model_version: D1c
activity:
  debye_huckel_A_20C: {value: 0.5050, unit: "1", source: assumed, note: "recalled standard value"}
  debye_huckel_A_25C: {value: 0.5091, unit: "1", source: assumed, note: "recalled standard value"}
  sit_b: {value: 1.5, unit: "1", source: assumed, note: "SIT convention"}
  delta_eps:
    sulfite: {value: 0.0, unit: kg/mol, source: assumed, bounds: [-0.45, 0.45], prior: normal, note: "PLACEHOLDER until stage A"}
agents:
  hydroquinone:
    k_et: {value: 0.05, unit: 1/s, source: assumed, bounds: [1.0e-4, 10.0], prior: lognormal, note: PLACEHOLDER}
    K_ads: {value: 100.0, unit: L/mol, source: assumed, bounds: [1.0, 1.0e4], prior: lognormal, note: PLACEHOLDER}
    lambda0: {value: 0.01, unit: 1/s, source: assumed, bounds: [1.0e-4, 10.0], prior: lognormal, note: PLACEHOLDER}
    activation_energy: {value: 67.5, unit: kJ/mol, source: literature, ref: "Sugai & Furuichi 1953 Table II (0.70 eV); KIN 1.6", uncertainty: 10.0}
    k_reg: {value: 100.0, unit: L/(mol*s), source: assumed, bounds: [1.0e-2, 1.0e4], prior: lognormal, note: PLACEHOLDER}
transport:
  delta_max: {value: 300.0, unit: um, source: assumed, bounds: [200.0, 500.0], prior: uniform, note: "PLACEHOLDER; prior KIN 3.4"}
gelatin:
  groups_per_gram: {value: 4.8e-4, unit: mol/g, source: literature, ref: "Iwano 1970; KIN 3.2"}
```

`data/films/kodak/tri-x-400.yaml`（节选）：

```yaml
schema_version: 1
id: kodak/tri-x-400
maker: Kodak
name: KODAK PROFESSIONAL TRI-X 400 (400TX)
kind: bw_negative
model: bw-mechanistic
iso_speed: {value: 400, unit: "1", source: datasheet, ref: "F-4017 Oct 2021 p. 6-7 (implied by push tables); DS 2.2"}
default_readout: scan_standard
default_format: "135"
layers:
  - {type: overcoat, id: overcoat, thickness: {value: 1.0, unit: um, source: assumed, note: "OPT 1.4"}}
  - type: emulsion
    id: emulsion
    sensitized_to: panchromatic
    gelatin_coverage: {value: 5.0, unit: g/m2, source: assumed, bounds: [3.0, 8.0], note: "LAT 2.3"}
    swelling_ratio: {value: 3.0, unit: "1", source: assumed, bounds: [2.0, 5.0], note: "KIN 3.3"}
    grain:
      morphology: conventional
      mean_diameter: {value: 0.9, unit: um, source: assumed, bounds: [0.4, 2.0], prior: lognormal, note: PLACEHOLDER}
      sigma_log: {value: 0.45, unit: "1", source: assumed, bounds: [0.2, 0.7], note: PLACEHOLDER}
      silver_coverage: {value: 4.0, unit: g/m2, source: assumed, bounds: [2.0, 6.0], prior: lognormal, note: "PLACEHOLDER; fitted in block F (D1 2.5.3); no primary value, LAT 2.3"}
      halide: AgBrI
      iodide_fraction: {value: 0.04, unit: "1", source: assumed, bounds: [0.02, 0.10]}
      n_classes: 6
    optics:
      gelatin_index: {value: 1.54, unit: "1", source: assumed, note: "unverified, OPT 1.1"}
      psf_length: {value: 2.0, unit: um, source: assumed, bounds: [0.5, 5.0], note: PLACEHOLDER}
      # 其余光学参数同 §5.4，略
    latent:
      tau0_median: {value: 0.22, unit: s, source: assumed, bounds: [0.01, 2.0], prior: lognormal, note: "PLACEHOLDER; geometric midpoint of LAT 1.4 range 0.05-1 s"}
      competing_sites: {value: 3.0, unit: "1", source: assumed, bounds: [2.0, 10.0], note: "§12.3"}
      # 其余同 §5.4，略
    spectral_sensitivity:
      quantity: log_spectral_sensitivity
      source: datasheet
      ref: "F-4017 Oct 2021 p. 8; shape only (old 5063 emulsion, level inconsistent, LAT 3.3); files re-assigned against the swapped printed legend, frame-edge points dropped (CSV legend_check, plot_floor)"
      file: reference/kodak/tri-x-400/spectral-sensitivity-d0.3.csv
  - type: base
    id: base
    material: cellulose triacetate
    refractive_index: {value: 1.48, unit: "1", source: literature, ref: "US 6,398,371 B1; OPT 1.1"}
    formats:
      "135":
        thickness: {value: 127, unit: um, source: assumed, note: "5 mil for legacy TX (F-9 2003); 400TX not stated; OPT 1.2"}
        grey_density: {value: 0.10, unit: "1", source: assumed, note: "35 mm - 120 base+fog offset, interpretation unverified; OPT 1.3"}
        backing_density: {value: 0.0, unit: "1", source: assumed}
      "120":
        thickness: {value: 91, unit: um, source: assumed, note: "3.6 mil legacy TX; OPT 1.2"}
        grey_density: {value: 0.0, unit: "1", source: assumed}
        backing_density: {value: 0.5, unit: "1", source: assumed, bounds: [0.2, 1.5], note: "OPT 1.3 recommendation"}
```

`data/developers/d-76.yaml`（成品 D-76）：

```yaml
schema_version: 1
id: d-76
name: KODAK Developer D-76 (packaged)
kind: developer
formula_confidence: msds-estimated
variant: packaged
base_solution: stock
components:
  - {species: sodium-sulfite, amount: {value: 100.0, unit: g/L, source: assumed, bounds: [84.0, 119.0], note: "PSI SDS 2024 80-95 %, kit mass assumed; DEV 1.2"}}
  - {species: hydroquinone, amount: {value: 5.0, unit: g/L, source: assumed, bounds: [3.2, 6.3], note: "DEV 1.2"}}
  - {species: metol-sulfate, amount: {value: 2.0, unit: g/L, source: assumed, bounds: [1.1, 3.8], note: "DEV 1.2"}}
  - {species: borax-pentahydrate, amount: {value: 1.53, unit: g/L, source: assumed, bounds: [1.1, 3.8], note: "2 g decahydrate equivalent; DEV 1.2"}}
  - {species: boric-oxide, amount: {value: 2.0, unit: g/L, source: assumed, bounds: [0.0, 3.8], note: "PLACEHOLDER; nuisance, fitted in stage A; DEV 1.2"}}
  - {species: dtpa, amount: {value: 0.5, unit: g/L, source: assumed, bounds: [0.0, 1.2], note: "DEV 1.2"}}
nuisance: [boric-oxide]
ph_anchors:
  - dilution: stock
    temperature: {value: 25, unit: C, source: datasheet, ref: "Kodak CIS-61 (2010); TMX 10.2"}
    ph: {value: 8.5, unit: "1", source: datasheet, ref: "Kodak CIS-61 (2010), typical value; TMX 10.2", uncertainty: 0.05}
    method: typical
standard_dilutions: ["stock", "1+1"]
refs: ["J-78 (Dec 2017)", "PSI SDS 5160304 v2.0 (2024)", "Kodak CIS-61"]
```

`data/processes/kodak/tri-x-400/d-76-stock-20c-large-tank-8min.yaml`：

```yaml
schema_version: 1
id: kodak/tri-x-400/d-76-stock-20c-large-tank-8min
name: Tri-X 400 in D-76 stock, large tank, 20 C, 8 min (F-4017 characteristic curve)
family: bw
steps:
  - kind: develop
    bath: d-76
    dilution: stock
    time: "8min"
    temperature: {value: 20.0, unit: C, source: datasheet, ref: "F-4017 Oct 2021 p. 9"}
    agitation: {method: lift, initial_s: 22.5, interval_s: 60.0, duration_s: 7.5}
    tank: {volume_ml: 3785.0, film_area_cm2: 35.0, image_area_fraction: 1.0}
  - {kind: stop, bath: stop-acetic, time: "30s"}
  - {kind: fix, bath: fixer-thiosulfate, time: "5min"}
  - {kind: wash, bath: water, time: "10min"}
note: "agitation per F-4017 p. 3 large-tank procedure; tank volume and strip area assumed"
```

`data/calibration/<run_id>.yaml`（节选）：

```yaml
schema_version: 1
run_id: d1e-PLACEHOLDER
created_utc: 2026-12-01T00:00:00Z     # PLACEHOLDER
model_version: D1e
profile: calibration
stages:
  - name: C_film
    parameters: ["films/kodak/t-max-100:emulsion.grain.thickness", "calibration:aspect_ratio"]   # T-Max d_50 = AR · t_g（§7.1）
    data_files:
      - {path: reference/kodak/t-max-100/char-d76-06min.csv, sha256: "PLACEHOLDER0000000000000000000000000000000000000000000000000000000"}
    table_rows: ["t-max-100/tables.yaml#facts/rms-granularity"]
    optimizer: {method: trf, loss: linear, x_scale: jac, starts: 8, seed: 0}
    chi2: 0.0          # PLACEHOLDER
    n_observations: 0  # PLACEHOLDER
    n_effective: 0.0   # PLACEHOLDER
    n_parameters: 0    # PLACEHOLDER
    converged: false
parameters: []
singular_values: {}
training_set_sha256: "PLACEHOLDER0000000000000000000000000000000000000000000000000000000"
frozen: false
held_out_evaluated: false
```

### 5.9 参考数据的位置

数字化团队把 `tools/digitize/staged/kodak/<film>/*.csv` 与 `tables.yaml` 并入 `data/reference/kodak/<film>/`（文件名不变），把 `staged/developers/*.yaml` 并入 `data/reference/developers/`。D0 的 `Library.reference_curve` 继续读 CSV；新增 `Library.reference_tables(material_id) -> dict` 与 `Library.developer_facts(developer) -> dict`（只读 YAML，结构按 `staged/developers/README.md`）。暂存的显影液名到药液 id 的对应：`d-76 → d-76`，`xtol → xtol`，`hc-110 → hc-110-pre2019`（2017 年的 J-24 与 2018 年之前的胶片手册早于配方更换，DEV §7 第 1 条），`t-max-dev → t-max-dev`，`t-max-rs → `（D1 不用）。

---

## 6 虚拟感光仪、密度计与指标

### 6.1 虚拟感光仪

| 项 | 定义 | 类别 | 出处 |
|----|------|------|------|
| 光源 | ISO 感光测定日光：`S_b = D55_b · τ_lens(λ_b)`；D55 取 CIE 数据（colour-science）；τ_lens 用 PCHIP 过 CAL §1.5 引用的 ISO 7589:1984 表 1 锚点（350 nm 0.00，380 0.42，400 0.74，440 0.94，460–640 0.97，690 0.94），690 nm 以上保持 0.94 | D55 literature；τ assumed（由锚点近似，非 ISO 全表） | CAL §1.5 |
| 钨丝灯（T-Max 的 MTF） | CIE 标准光源 A | assumed（"Tungsten" 的解读） | DS §3.1 |
| 曝光时间 | 1/50 s | assumed | DS §5.1、CAL §1.1（F-9 的唯一明示值） |
| 级差 | log H 等距 0.15，范围 [min(x_ref) − 0.30, max(x_ref) + 0.30]，无参考曲线时 [−4.0, +1.0]；外加一个零曝光片 | — | CAL §1.6 |
| 每级曝光 | `k_j = 1e-12 · H_j / 683`（光源已归一到 Σ S ȳ = 1，D0 §3.3 固定曝光），`N_b = k_j S_b / E_b` | — | D0 §3.2 |
| 几何 | 每级一根深度柱，没有横向；光晕对均匀场贡献 f_H（各级宽数毫米，r99 ≈ 1.2 mm 基本落在本级内，`assumed`） | — | CAL §1.6 |
| 药液槽 | 条带面积每级 1 cm²；罐体积：小罐 0.473 L，大罐 3.785 L，转筒 0.25 L，盘 1.0 L（耗尽可忽略） | assumed | — |
| 读数 | 每级取均值透射（不加颗粒），视觉与 Type 2 两种密度 | — | §2.8 |
| 显影系列 | 一次运行，在各显影时间记录快照，分别接停显–定影–水洗 | — | §3.4 |

```python
@dataclass(frozen=True)
class Condition:
    material_id: str
    film_format: str
    developer: str
    dilution: str
    temperature_c: float
    agitation: Agitation
    times_s: tuple[float, ...]
    tank: Tank
    exposure_time_s: float = 0.02
    illuminant: Literal["iso_daylight", "A"] = "iso_daylight"
    source: str = ""                    # 参考文件或表格行的标识（§7.6）

@dataclass(frozen=True, eq=False)
class StripResult:
    condition: Condition
    log_h: np.ndarray                   # (J,) log10(lux*s)
    times_s: np.ndarray                 # (T,)
    visual: np.ndarray                  # (T, J)
    type2: np.ndarray                   # (T, J)
    base_fog_visual: np.ndarray         # (T,) 零曝光片
    base_fog_type2: np.ndarray          # (T,)

def sensitometric_strip(condition: Condition, ctx: RunContext, *, log_h: np.ndarray | None = None) -> StripResult
```

温度用 Kodak 原表的 °F 精确换算（65/68/70/72/75/80 °F = 18.33/20.00/21.11/22.22/23.89/26.67 °C），不用四舍五入后的 °C 标签（KIN §1.6 的更正）。

### 6.2 密度计

见 §2.8。条带上逐级读数，均值透射，`aperture_um = None`。

### 6.3 ISO 感光度（ISO 6:1993）

1. 用 Type 2 近似密度。曲线 D(log H) 用单调 PCHIP 插值。
2. D_bf 取零曝光片（ISO 6 §5.6.2），只用于绝对 ISO 感光度的陈述与诊断。**凡与 Kodak 曲线比较**（拟合残差、§8.3 的全部曲线指标），两侧都用同一个算子 `D_min* = D_ref(x_left)`：x_left 为参考曲线最左点的 log H，模拟曲线也在同一 x_left 处取值作自己的 D_min*（Kodak 没有给出任何基雾值，曲线左端只是基雾的上界，暂存 `tables.yaml#facts/base-plus-fog`；左端处的斜率 0.002–0.08 D/log H，用零曝光片作模拟侧的基雾会使两侧算子不同，log H_m 偏 0.03–0.07）。此时 D_bf := D_min*，得 log H_m*。
3. m 点：D 第一次达到 D_bf + 0.10 处，得 log H_m（§5.6.3）。
4. ΔD = D(log H_m + 1.30) − D(log H_m)。
5. 显影系列 {t_i}：在 (t, ΔD) 上找夹住 0.80 的相邻两点线性内插得 t*，在 (t, log H_m) 上同样内插得 log H_m(t*)（附录 A）。夹不住时返回 NaN 并给出原因。
6. `S = 0.80 / H_m`，`S° = 1 + 10 log10(0.80/H_m)`；按 ISO 6 表 1 归档（ISO 400 ↔ −2.75…−2.66，ISO 100 ↔ −2.15…−2.06），报告同时给 log H_m 的数值（CAL §1.1）。

测试（D1e）：`test_iso_speed_on_kodak_tri_x`：对暂存的 Tri-X D-76 曲线（视觉密度，`iso_speed(..., base_fog=D_min*)`，即各曲线最左点，与 DS §5.3 的做法相同），35 mm 得 log H_m ≈ −2.73，120 得 −2.69，误差 ≤ 0.03（DS §5.3 的值是"≈"的估计，暂存 CSV 的 x 误差约 ±0.02）。这是指标实现的检查，不是对绝对感光度假设的检验（§9 第 22 条）。`test_dmin_star_same_operator`：模拟条带与参考曲线的 D_min* 在同一 x_left 取值；参考曲线整体加 0.02 D 时，log H_m* 与 CI 的两侧差不变（≤ 1e-9）。

### 6.4 ISO 平均斜率

`G_ISO = ΔD / 1.30`，在指定显影时间上计算（CAL §1.3.1）。

### 6.5 反差指数（Kodak H-740 的 0.2/2.2 共线作图）

D 与 log H 用同样的单位长度。D_min：从模拟曲线算 CI、去和 Kodak 特性曲线上读出的 CI 比较时，取 §6.3 第 2 步的 D_min*（两侧同一算子，视觉密度；暂存 CSV 的 `ci_cross_check` 也按曲线左端密度作图）；与 Kodak 的 CI(t) 曲线、时间表格比较时（假定它们是 Kodak 按真实基雾作图的结果，H-740，`assumed`），以及仅预测的条件，取零曝光片（视觉密度）。CI 对 D_min 的敏感度：暂存的 D-76 曲线上 D_min 降 0.01 D，CI 降 0.003–0.007（降 0.02 D 降 0.007–0.016，§12.10），所以两类比较必须各自固定算子，且 F-4016 两张图之间的差（§7.2）有一部分可能来自左端代理偏高。

```
C = (x_C, D_min)
P₁ = 曲线上 x > x_C、与 C 的欧氏距离为 0.2 的第一个点
P₂ = 曲线上 x > x_C、与 C 的欧氏距离为 2.2 的第一个点
求 x_C 使 atan2(D₁ − D_min, x₁ − x_C) = atan2(D₂ − D_min, x₂ − x_C)（brentq，区间 [x_lo − 1, x_hi]）
CI = (D₂ − D₁) / (x₂ − x₁)
```

交点：沿 PCHIP 曲线以 Δx = 0.002 采样找出距离跨过 r 的区间，再用 brentq 精确求解。无解（曲线过短或过平）时返回 NaN 并说明。来源：Kodak H-740（2006）第 9 页，CAL §1.3.2 核实（`literature`）。另报 0.10 + 2.0 弦近似作交叉检查（CAL §1.3.2），不作拟合与判定。Z-133 的控制条公式与控制条目标值不得用于模拟曲线（CAL §1.3.3）。

测试（D1e）：

| 测试 | 断言 |
|------|------|
| `test_ci_straight_line` | 铰链折线 D = D_min + G·max(x − x₀, 0)，P₁、P₂ 都在直线段时 CI = G，误差 ≤ 1e-9 |
| `test_ci_matches_report_on_tmax_d76` | 暂存 `char-d76-{06,07.5,10}min.csv`，以曲线最左点作 D_min：CI = 0.542/0.693/0.895（误差 ≤ 0.01），G_ISO = 0.547/0.701/0.879（≤ 0.01），log H_m = −2.011/−2.129/−2.212（≤ 0.02）。参考值是 CAL §1.3.4 用它自己的图框标定从 PDF 矢量算出的；暂存 CSV 是另一次数字化，表头给出 x 误差约 ±0.02 log H、y 误差约 ±0.016，轴的整体平移直接移动 log H_m 而几乎不影响 CI 与 G_ISO，容差按此设定 |
| `test_ci_nan_on_short_curve` | 曲线跨度 < 2.2 时返回 NaN |

### 6.6 RMS 颗粒度

Kodak 的定义："净漫射密度 1.0 处、48 µm 孔径、12 倍放大读数"，数值为 1000 σ_D（OPT §5.2）。

1. 求均匀曝光 H₁ 使净视觉密度 = 1.00（从条带插值得初值，再在单柱上做固定 4 次割线迭代，不按容差提前停止，§7.4），同一流程显影。
2. 统计模式在 3.2 mm × 3.2 mm 的栅格上渲染，栅格间距 p_g = 4 µm，800² 像素；均匀曝光下化学是均匀的，只需一根柱的显影结果铺满全块，再逐像素抽样。孔径方差是 Var_bool 的大像素极限，与 p_g 无关（§2.7.1）；4 r_c > p_g 的档按均值档处理，b_ap 与截断比例随 G_rms 一并报告，b_ap > 0.05 或截断比例 > 1e-3 时警告（此时 G_rms 偏高）。
   标定前的合理性检查：Siedentopf 估计 `G_S = 1000 · sqrt(κ(1.0) · 0.4343 · ⟨a²⟩/⟨a⟩ / A_ap)`（⟨a²⟩/⟨a⟩ = Σ_c ν_c a_c² / Σ_c ν_c a_c，完全显影；A_ap = π · 24² um²；κ(1.0) = 1.89）。先验中心的 Tri-X（d_50 0.9 µm、CP_ref 0.6）给 ⟨a²⟩/⟨a⟩ ≈ 14 um²、G_S ≈ 80，而 OPT §5.4 由颗粒度 17 推出的有效盘为 0.64 um²；T-Max 的 G_S 与 8 相差 2 倍以上的起点从多起点（§7.4）中剔除；Tri-X 的 17 是检查项（§7.6），只在标定记录中标出 G_S，不用于筛选（泄漏规则，§7.7）。
3. 48 µm 直径圆孔，按 48 µm 间距的方格不重叠排列（66² = 4356 个）；孔径透射为孔内像素透射按面积覆盖权的平均；`D_ap = −log10` 视觉光谱乘积下的孔径透射。
4. `G_rms = 1000 · std(D_ap)`，标准误差 `≈ G_rms / sqrt(2(N − 1))`（约 1.1 %），一并报告。
   拟合目标不用抽样值：泊松逆变换给出整数，固定均匀数下抽样的 G_rms 是参数的分段常数函数，`test_forward_model_smooth` 不能通过，且 1.1 % 的抽样误差（G = 8 时 0.09）占了 σ = 0.6 的一部分。阶段 C 的残差用同一统计模式的期望值（delta 方法）：`G_exp = 1000 · 0.4343 · sqrt(Σ_i a_i² Var_pix) / (T̄ Σ_i a_i)`，a_i 为孔内像素的面积覆盖权，Var_pix = Var_count + (T̄² + Var_count) Var(ε_pix)（= Var_bool，截断前）为 §2.7.1 的逐像素透射方差（均匀曝光下各像素相同）。抽样的 G_rms 只用于报告与跨模式检查（`test_granularity_expected_vs_sampled`：二者之差在 3 个标准误差内）。
5. 显影条件：T-Max 100 为 D-76、20 °C、小罐，显影到 CI 0.56 的时间（F-4016 的 6½ min，Kodak 未给时间，`assumed`）；Tri-X 为 HC-110 B、20 °C、小罐，显影到 CI 0.56（F-4017 的 3¾ min，`assumed`）。
6. 测微密度计的光学（数值孔径）不明（OPT §5.2），按理想孔径处理（`assumed`）。

### 6.7 MTF（有效曝光法）

1. 一个空间周期 Λ = 1/ν 的周期条带，N_x = 32 列（Δx = Λ/32），参考档离散化。
2. 曝光 `H(x) = H₀ (1 + m cos(2πx/Λ))`，m = 0.4（`assumed`），H₀ 取使大面积平均净密度为 1.0 的值（`assumed`，Kodak 未给，OPT §0 第 7 条）。光学在傅里叶空间精确处理：每个深度单元的曝光描述量调制为 `m · M_X,z(ν)`，`M_X,z(ν) = (W⁺_X,z·N P̂(ν) + W⁻_X,z·N Ĥ(ν)) / (W⁺_X,z·N + W⁻_X,z·N f_H)`（§2.1.4），Ĥ 为光晕核的汉克尔变换。
3. 全化学显影（横向扩散用 FFT，周期边界），读出每列的视觉密度（均值，无颗粒），理想测微密度计（`assumed`）。
4. 用同一次模拟的大面积特性曲线把每列密度反查为有效曝光 H_eff(x)。
5. 取 H_eff 的一阶傅里叶系数 a₁ 与均值 a₀，`MTF(ν) = 100 · (a₁/a₀)/m`（%）。
6. 频率取数字化曲线上的 2.5、3、5、7、10、15、20、30、40、50、60、70 c/mm（T-Max 另加 80、100、120、150）。
7. 曝光与显影条件：T-Max 100 为钨丝灯（A）、D-76、20 °C、小罐，显影到 CI 0.56；Tri-X 为 ISO 日光（未给，`assumed`）、D-76、大罐、20 °C（未给，`assumed`），显影到 CI 0.56。画幅：两种胶片都按 135 模拟（暂存 CSV 写 "format=not stated"，`assumed`；135 灰片基与 120 防光晕层的光晕不同，这个选择在冻结前固定）。
8. 频段（预注册）：**尾段**为参考 MTF 最后一次 ≥ 100 % 的频率以上的各点（T-Max 约 19 c/mm 以上，Tri-X 约 10 c/mm 以上，OPT §2.4），**峰段**为其余各点。光学传递 ≤ 1、光晕只降低响应，峰段的超出 100 % 全部来自横向化学（δ_max、δ_agit、r_gel、K_Br、u_X），w_P、a_P 管不到它；所以只有尾段进入阶段 C 的拟合，峰段只作横向化学的内部检查（§7.6），不拟合、不计入组合状态。D1 没有其他数据检验邻接效应与光晕。

标定中（阶段 C）不用上面的非线性周期条带：参考档的自适应步长会给有限差分雅可比带来步长噪声（§3.4），每个频率 32 列加 H₀ 搜索也太贵（§3.10）。阶段 C 用小信号（切线性）响应，标定档、固定时间表、一次 `run_ensemble` 调用：

```
每个频率 ν 一对柱 (a, b)：曝光描述量 e_z · (1 ± ε M_X,z(ν))，ε = 1e-3；L 步只作用于差模（§3.5），波数 k = 2πν
另一对 (a₀, b₀)：e_z · (1 ± ε)，无横向耦合（ν → 0）
MTF(ν) = 100 · (D_a − D_b)(ν) / (D_a₀ − D_b₀)
```

调制 m 本来就是 `assumed`；小信号极限下有效曝光法的反查是线性的，上式就是 a₁/(a₀ m) 的 m → 0 极限。H₀ 由条带插值得初值、固定 4 次割线迭代求得（§7.4）。非线性周期条带（m = 0.4，参考档）只用于验证报告。

测试（D1e）：`test_mtf_linear_system`：把显影换成线性响应（D = γ log H、无横向化学）时，测得的 MTF 等于光学传递函数 `M_X,z(ν)` 的深度平均，误差 ≤ 1e-3（周期条带与小信号柱对都测）。`test_mtf_periodic_exact`：一个周期与两个周期的条带结果之差 ≤ 1e-10。`test_mtf_small_signal_limit`：同一档位、同一时间表下，m = 1e-3 的周期条带与小信号柱对之差 ≤ 0.1 个百分点。`test_mtf_small_signal_vs_periodic`：标定档小信号与参考档 m = 0.4 周期条带，≤ 30 c/mm 时差 ≤ 3 个百分点，结果写入验证报告。

### 6.8 倒易律

| 量 | 定义 |
|----|------|
| 参考 | t_ref = 1/100 s（两张表都在 1/1000–1/10 s 内无曝光修正，DS §5.1）下的标准显影（胶片在其参考流程下显影到 CI 0.56 的时间） |
| 判据密度 | 净视觉密度 0.10 与 0.80 两个（Kodak 未说明判据，CAL §1.7）；拟合残差取两者的平均档数 |
| 光圈行（固定时间 t） | `Δstops(t) = log2(H(t)/H(t_ref))`，H(t) 为在时间 t 达到判据所需的曝光 |
| 修正时间行（固定照度） | 照度 I = H(t_ref)/t_ind，求达到判据的时间 t_adj，与 Kodak 的修正时间比较 log10 |
| 显影修正（Tri-X） | 在 t_ind 曝光后，求使 CI 恢复到参考 CI 的显影时间变化（%），与 Kodak 的 +20 … −30 % 比较，只作方向与幅度检查（±10 个百分点） |

### 6.9 光谱感光度

1. 每个 5 nm 波段做窄带曝光（全部光子落在该波段），曝光时间取手册的有效曝光时间：T-Max 1.4 s、Tri-X 0.5 s（LAT §3.3、§5 第 2 条：这两个时间处在低照度失效的起点，必须按原时间模拟）。
2. 显影：T-Max 为 D-76、20 °C；Tri-X 手册未给，取 D-76 大罐 20 °C（`assumed`）；都显影到 CI 0.56。
3. 求使净视觉密度达到 0.3 与 1.0（T-Max："高于 D-min"；Tri-X："高于总灰雾"）的能量曝光 `E_b = N_b · (h c / λ_b) · 1e7 · 1e8`（erg/cm²；N_b 为 photons/um²，h c/λ_b 以 J 计，即 D0 的 `photon_energy_j`，不是 §2.1.2 以 eV 计的 E_ph,b），`log S = −log10 E_b`。测试 `test_spectral_energy_conversion`（D1e）：550 nm 时 1 erg/cm² = 2769 photons/um²（LAT §3.3），相对差 ≤ 1e-3。
4. 拟合只比较形状：每条曲线减去它在 450–600 nm 上的最大值（LAT §3.3 第 5 步：绝对水平不用来定标）。拟合用 400–700 nm 每 20 nm 一点（16 个波段），但只取落在该参考曲线数字化跨度内的波段（Tri-X 两条曲线触到图框底边的末点已被数字化团队删去，D = 0.3 曲线止于 664.0 nm、D = 1.0 曲线止于 655.3 nm，所以 Tri-X 没有 660 nm 以上的点）；报告用跨度内全部 5 nm。
5. Tri-X 文件的归属：F-4017 印刷的图例把两条线标反了，暂存文件已按物理约束（0.3 曲线必在 1.0 曲线之上）重新归属，`spectral-sensitivity-d0.3.csv` 就是物理上的 D = 0.3 曲线（表头 `legend_check`）；按文件名使用，不再交换。Tri-X 的曲线属于旧乳剂 5063（DS §2.1），绝对水平即使按正确归属也比 ISO 400 预期高约 0.5 log（LAT §3.3），只用形状。
6. T-Max 100 的绝对水平（D = 0.3 与 1.0 高于 D-min，1.4 s）作**不拟合的内部检查**：模拟的绝对 log S 与参考直接比较，不减最大值。它检验从 lux·s 到 erg/cm² 的换算链（ISO 日光光谱、τ_lens、光谱形状在明视范围的积分）与 1.4 s 的低照度失效；η_0 与吸收尺度的乘积已由 D-76 曲线（以 lux·s 计）定出，所以它不检验二者本身。这是训练数据中唯一不由 D-76 曲线拟合出来的绝对感光度证据（§9 第 22 条）。报告 450–600 nm 上的平均差；|平均差| > 0.2 log S 时标出（`assumed`，LAT §3.3 的水平是目读的）。

### 6.10 到达 CI 的时间

对一个条件，一次运行到 t_max（表中最长时间的 1.5 倍），每 15 s 记录快照，得 CI(t)；单调插值后求 CI(t*) = 目标值。超出范围时把 t_max 放大到 1.5 倍再算（最多 3 倍），仍不行返回 NaN。

`time_to_ci` 只用于报告（诊断、图）与 §6.6–§6.9 的"显影到 CI 0.56"条件，不进入拟合残差与 §8.3 的判定。时间表格逐格按"表中时刻的 CI"评分（§7.2），不需要任何斜率换算，只要 CI 可算就有限。

---

## 7 标定方案

### 7.1 参数向量

变换：正参数用 log10，(0, 1) 内的分数用 logit，其余线性。先验作为额外残差 `(g(θ) − μ)/τ` 加入（CAL §3.2），盒约束作为硬界。

| 块 | 参数（字段路径见 §2 与 §5） | 个数 | 约束它的数据 | 阶段 |
|----|------------------------------|------|--------------|------|
| A 活度 | Δε_sulfite、Δε_borate、Δε_organic、Δε_amine | 4 | pH 锚点（§7.2） | A |
| A 干扰 | `d-76` 的 `boric-oxide`；`t-max-dev` 的 S(IV) 配料量；`hc-110-pre2019` 的未申报碱（以 KOH 当量的配料表示）；Carlton & Crabtree 比色系列的整体 pH 偏移 δ_CC（只在标定中存在，不写回数据文件；线性，界 ±0.5，均匀） | 4 | CIS-61 的 D-76 8.5、T-MAX 8.37、HC-110 A 稀释 8.9（各一个锚点）；δ_CC 由该系列自身与其余锚点共同定 | A |
| B 光谱（每种胶片） | α_0、E_0、E_U、两条（必要时三条）染料带的 A、ν̃、w | 9–12 | 胶片自己的光谱感光度两条曲线（形状） | B |
| F 胶片（每种胶片） | d_50（T-Max 改为 AR，见下）、σ_d、η_0、CP_ref、m_Ag,coat、h_fog、S_550、s_abs（吸收绝对尺度，见下）、w_P、a_P、τ̄₀、σ_τ、τ_i；T-Max 另有 t_g；可选 ξ_η | Tri-X 13–14，T-Max 14–15 | D-76 特性曲线、D-76 MTF 尾段（§6.7 第 8 步）、倒易律拟合行；T-Max 另加 RMS 颗粒度 8 | C |
| C 显影剂 | 米吐尔：k_et、K、λ0、k_des、Ea；HQ：k_et、K、λ0、k_reg、Ea；吡唑烷酮：k_et、K、λ0、k_des、pK；抗坏血酸：k_et、K、λ0、k_reg | 19 | 训练组合的全部 CI(t)、时间表、Tri-X×T-MAX 曲线、D-76 曲线 | D |
| C 共用 | K_Br、u_gel、u_X、u_LI、ζ_s、k_diss、C_sat,ref、ν_s、k_pd、δ_max、δ_agit | 11 | 同上；δ 主要靠小罐/大罐/转筒的差别；ν_s 靠亚硫酸盐水平不同的训练显影液（D-76、D-76 1+1、XTOL、T-MAX、HC-110 B） | D |
| C 干扰 | `t-max-dev`：HQ、Dimezone S、初始 NaBr；`hc-110-pre2019`：HQ、methylphenidone、KBr | 6 | Tri-X×T-MAX 曲线与时间表；Tri-X×HC-110 CI(t) 与时间表 | D |

每个拟合参数在上表中只出现一次。派生量的约束：

- 吸收绝对尺度 `s_abs`（log10，界 [−1, 1]，均匀；只在标定中存在）：阶段 B 只定光谱形状（α_0 与各 A_j 的比、E_0、E_U、ν̃_j、w_j），阶段 C 用 s_abs 同时乘 α_0 与全部 A_j；写回时把 s_abs 并入 α_0、A_j 的值并复位为 0。
- 层吸收率 A_550 ∈ [0.2, 0.6]（§2.1.8）：单侧罚残差 `r = max(0, A_550 − 0.6, 0.2 − A_550)/0.02`（`assumed`，区间内为零，不拉向任何值），加入阶段 B、C、E 的目标函数。
- 纵横比 AR = d_50/t_g ≥ 8（T-Max）：t_g 拟合（§2.0.2），所以 T-Max 用 (t_g, AR) 参数化：t_g 线性、界 [0.08, 0.15]；AR 取 log10、界 [log10 8, log10 19]；d_50 = AR · t_g 的界 [0.6, 1.5] µm 用与 A_550 同形式的单侧罚残差（尺度 0.02 µm）。报告给出 T-Max 留出预测对 t_g 的敏感性（§9 第 31 条）。
- D_g、D_back、t_b 不拟合（`assumed`，§2.1.8）：D_g 决定光晕返回 f_H 与均匀场感光度，只能靠两种画幅的速度差来定，而 Kodak 曲线的左端只是基雾上界，不足以分开 D_g 与 h_fog。所以两种画幅的差别完全由 assumed 片基参数给出，§7.6 的格式一致性检查因此不循环。

固定不拟合的（`assumed` 或 `literature`）：§2 各参数表中标为 assumed 的全部（含 ν_b、ΔH_Br、D_H、D_OH），加上吡唑烷酮/抗坏血酸/对氨基酚的 Ea（§2.4.5）、对氨基酚的全部动力学常数（取米吐尔值，Rodinal 只作预测）。

### 7.2 观测量、残差与权重

| 观测量 | 残差 | σ | 采样 | 出处 |
|--------|------|---|------|------|
| 特性曲线密度 | D_sim − D_ref（两侧都在参考曲线的 log H 上取值；D_min* 见 §6.3 第 2 步） | 0.025 × √(15/6) = 0.040（D < D_min*+0.1 的点权重 0.5） | 每条曲线在 log H ∈ [D_ref = D_min*+0.05 处, 最后一点] 上等距重采样 15 点（原 CSV 约 137 点高度相关）；一条曲线约 6 个独立数（§7.5），σ 乘 √(N_pts/N_eff) | CAL §2.3、§3.2 |
| 给定时刻的 CI | CI_sim − CI_ref（模拟侧 D_min 取零曝光片，§6.5） | √(0.025² + σ_floor²) × √(5/3)；σ_floor：T-Max 100 取 0.018，Tri-X 取 0 | 每条 CI(t) 曲线在其时间跨度上等距 5 点（约 3 个独立数） | CAL §3.2；σ_floor 见下 |
| 表格格的 CI | CI_sim(t_cell; T_cell, 稀释, 罐型) − CI_aim(格)；CI_aim：F-4016/F-4017/J-24/J-86 的时间表为 0.56（J-86 印明 0.56；F-4016 不印 CI 目标，暂存 `facts/processing-aim`，取 0.56 为 `assumed`），J-109 表 1 为该列的 CI | `√(σ_tab² + (s_ref·0.072 min)²) × √n_class`；σ_tab = 0.02（`assumed`）；0.072 min = ¼ min 取整的均匀分布标准差 0.25/√12；s_ref 为同一训练组合小罐 20 °C 参考 CI(t) 曲线在 CI_aim 处的斜率（训练 CSV 算一次后固定；CI_aim 超出曲线跨度时取端点斜率）；n_class 为同一（显影液、稀释、罐型、CI_aim）类中入拟合的非 20 °C 格数（20 °C 格取 1），使每类的温度斜率只按约一个独立数计（KIN §1.6：Kodak 的行很可能共用一个时间–温度换算）；XTOL 原液非 20 °C 格再乘 2（KIN §1.6 只把 XTOL 行认作 D-76 斜率的副本） | 逐格 | CAL §4.2 的判据改写为 CI 单位；KIN §1.6 |
| pH 锚点 | pH_sim − pH_ref（Carlton & Crabtree 系列为 pH_sim + δ_CC − pH_ref） | 锚点自带不确定度（0.05–0.1）；Carlton & Crabtree 比色值 0.10（只有相对值可靠，DEV §1.3，整体偏移由 δ_CC 吸收）；XTOL A 部分 0.3（配方中 DTPA 形态不确定，§12.4） | 逐个 | DEV §1.3、§2.1、§5.1；TMX §4.1、§10 |
| 光谱感光度形状 | log S 差 | 0.03 | 20 nm 一点，只取参考曲线跨度内（§6.9 第 4 步） | LAT §3.3 |
| MTF 尾段 | 百分点差 | 2（Tri-X 50 c/mm 以上 3） | §6.7 第 8 步的尾段频率；峰段不入拟合 | OPT §2.4 |
| A_550 罚项 | §7.1 的单侧罚残差 | 0.02 | 每种胶片一个 | §2.1.8 |
| RMS 颗粒度（T-Max） | G_exp − 8（期望值，§6.6；不用抽样值） | 0.6 | 一个 | OPT §5.2 |
| 倒易律 | 档数差 | 0.15 | 拟合行 | LAT §4 |

不重复计数：小罐 CI(t) 曲线已经给出 20 °C 下到达 CI 0.56 的时间，所以同一显影液、同一稀释的小罐 20 °C 时间表格不再入拟合（CAL §3.2）；F-4017 的 Tri-X 转筒表与小罐表逐格相同（暂存 `tables.yaml` 注明 "Every cell equals the small-tank table"），转筒行不入拟合。推冲表（无 CI 目标）、大罐–机器冲洗的时间区间（"5 to 7"）、滤镜系数、日光曝光表不入拟合。ISO 速度点 log H_m 不入拟合（它由同一批 D-76 曲线算出，已在密度残差里，CAL §3.2），只作报告指标。F-4016 的 `development_time_adjustment_factors` 不入拟合：`[d-76]` 块是 T-Max×D-76 时间表的取整指导值（与时间表逐格比较 RMS 3.9 %、最大 7.1 %，§12.10），归内部检查；反差列（−20/+20/+40 %）在任何块都不评估（没有预注册的 CI 目标）。

**σ_floor（F-4016 自相矛盾）。** F-4016 第 8 页的特性曲线读出的 CI 比第 9 页同显影液、同罐型、同时间的 CI(t) 曲线高 0.035–0.065（暂存 `t-max-100/tables.yaml#facts/characteristic-vs-contrast-index-curves`，各 CSV 表头的 `ci_cross_check`；训练对唯一可比的一格是 D-76 小罐 7.5 min：0.694 对 0.658）。同一参数不可能同时贴合两张图到 0.04 CI 以内，所以 T-Max 100 的 CI 残差加数据底 σ_floor = 0.036/2 = 0.018（只由训练对得出），使拟合折中于两图之间。§12.10 只用训练曲线核对了折中点同时满足 §8.3 的训练与留出阈值。这个差有一部分可能来自 Kodak 曲线左端的基雾代理（§6.5），不作修正，只报告。Tri-X 没有这种交叉核对，σ_floor = 0（`assumed`）。

**残差恒有限。** 任何残差为 NaN 或无穷时（CI 作图无解、条带不收敛），该项取有界罚值 `0.5/σ`（`assumed`，相当于 0.5 CI 或 0.5 D 的偏差）并计数；目标函数返回前断言全部有限（`test_residuals_always_finite`）。罚值项的个数写入 `StageRecord.n_penalised`（§5.7），最终结果中必须为 0。

pH 锚点：

| 锚点 | 值 | 用于 |
|------|----|------|
| D-76 公开配方及 Carlton & Crabtree 系列：无硼砂 8.1、2 g/L 8.5、6 g/L 8.9、20 g/L 9.2；米吐尔 4 + HQ 5 + 硼砂 2 → 8.1；米吐尔 2 + 硼砂 2、无 HQ → 8.6 | 比色，±0.1；作者说只有相对值可靠（DEV §1.3），整个系列共用一个偏移 δ_CC | Δε_borate、Δε_organic 与 Δε_sulfite 的相对部分（绝对的 Δε_sulfite 靠 XTOL 8.20 与低权重的 A 部分） |
| 成品 D-76 原液 8.5 ± 0.05（CIS-61，25 °C） | 典型值 | B₂O₃ 干扰参数 |
| XTOL 工作液 8.20 ± 0.05（专利例 1；J-109、CIS-61 相同） | 玻璃电极 | Δε_sulfite、Δε_organic |
| XTOL A 部分（850 mL）10.19 ± 0.05 | 玻璃电极，低权重 | 同上 |
| T-MAX Developer 1+4 工作液 8.37 ± 0.05（CIS-61，25 °C） | 典型值 | T-MAX 的 S(IV) 干扰参数（Δε_amine 不可同时确定） |
| 旧 HC-110 A 稀释 8.9 ± 0.10（CIS-61） | 典型值 | HC-110 的未申报碱 |

### 7.3 分阶段拟合

| 阶段 | 内容 | 固定的 | 数据 |
|------|------|--------|------|
| A | 活度 Δε 与 pH 干扰参数 | 其余全部 | pH 锚点 |
| B | 每种胶片的光谱形状参数，用代理 `S_q(λ) ∝ Σ_c w_c σ_c(λ)`（LAT §3.3 第 3–4 步）；权重 `w_c ∝ ∂D/∂ln m_c`（判据密度处，由完整模型算出）；B 与 C 交替两轮后用完整模型核对一次 | 化学取先验 | 光谱曲线形状 |
| C | 胶片块（每种胶片分别）；倒易律拟合行；MTF 用 §6.7 的小信号柱对（标定档、固定时间表），颗粒度用期望值 G_exp（§6.6） | 化学取当前值 | 该胶片的 D-76 曲线、D-76 MTF 尾段、倒易律拟合行、A_550 罚项；T-Max 的颗粒度 8 |
| D | 化学块与药液干扰参数 | 胶片块；阶段 A 的 Δε 与 δ_CC | 所有训练组合（§7.6），另加全部 pH 锚点残差（阶段 D 改动 HQ、Dimezone S、methylphenidone 配料量与 pK_pyr 时，阶段 A 的干扰参数仍须满足各自锚点；锚点残差几乎不增加计算） |
| E | 块坐标迭代：C（胶片块只看 D-76 与本胶片自己的数据）与 D（化学块看全部训练数据）交替，直到参数变化 < 0.1σ 且 χ² 相对变化 < 1e-3（最多 5 轮）；报告第 1 轮与最终胶片参数的偏移，偏移大说明胶片参数在吸收化学失配（CAL §3.3） | — | 同 C、D |
| F | 可辨识性（§7.5）；把不可辨识或先验主导的参数固定在**阶段 E 的当前拟合值**（不是先验均值：均匀盒没有可用的"先验值"；改标 `assumed`，`note` 写明"在标定 <run_id> 中不可辨识，固定于阶段 E 值 <v>"），重做 E | — | — |
| G | 写回（§7.8），训练报告交所有者确认阈值；D1f 的跨模式测试通过后冻结（代码提交、数据哈希、阈值、注册表与指标代码哈希，标定记录 `frozen: true`），然后做唯一一次留出评估（§8.4） | — | — |

阶段 E 是对 CAL §3.3"联合精修"的修改：精修时胶片参数仍只看 D-76 与本胶片数据，以满足既定决策 2（§11 第 3 条）。

### 7.4 优化器

- `scipy.optimize.least_squares(method="trf", x_scale="jac")`，变换空间；探索阶段 `loss="soft_l1", f_scale=1.0`，最终结果 `loss="linear"`，使 χ² 保持意义（CAL §3.3）。
- 雅可比：二点有限差分，`diff_step = 1e-3`（变换单位），由 `calibrate/objective.py` 自己计算后作为 `jac` 回调交给 `least_squares`：每一列只重算读到该参数的条件（依赖表由条件的显影液、稀释与胶片决定，例如米吐尔的参数只影响 D-76 的条件），其余条件复用基点结果，结果与全部重算逐位相同（`test_sparse_jacobian_matches_dense`）。正向模型的数值噪声须远小于这个步长：固定时间表、牛顿容差 1e-12、表格插值二阶连续（§2.2.5）、F 精确求值（§2.4.2）、搅动面落在单元界面上（§2.4.1）、颗粒度用期望值（§6.6）（`test_forward_model_smooth`：参数扰动 1e-6 时目标函数的二阶差分与扰动平方之比稳定）。
- 内层求解不按容差提前停止：求判据密度的曝光（§6.6、§6.7、§6.8、§6.9、`exposure_for_net_density`）从条带的 PCHIP 插值得初值，固定做 4 次割线迭代；CI 作图的交点与共线点用 brentq，`xtol = 1e-14`、`rtol = 4 eps`，到机器精度。按容差（如 ±0.002）停止会在参数空间留下 0.002 D 量级的跳变，与 1e-3 参数步长的效应同量级。
- 多起点：阶段 C、D 各 8 个 Sobol 点（`scipy.stats.qmc.Sobol(scramble=True, seed=0)`），落在先验盒内；每个起点先做 10 次迭代（`soft_l1`），χ² 最好的 2 个再继续到收敛（最多再 20 次，§3.10 的预算）。
- 全局：只当多起点的最优 χ² 相差超过 3.84 时，在已识别的 ≤ 8 个硬方向上用 `scipy.optimize.differential_evolution(seed=0)`（CAL §3.3）。它不在 §3.10 的预算内：默认种群 15 × 维数、每代一次求值，8 维约 120 次求值一代，阶段 D 每代约 45 min，几十代就要 1–2 天；需要时另行安排。

### 7.5 可辨识性

方法：阶段 E 收敛后，对每个块计算加权雅可比 `W^{1/2} J`（变换空间）的 SVD，**只用数据残差**（去掉先验残差与 A_550 罚项行，`residuals(..., include_prior=False)`）；奇异值相对最大值 < 1e-3 的方向上分量绝对值 > 0.3 的参数标记为"松散"。对打算报告为 fitted 的主要参数做剖面似然（Raue 等 2009），按 §3.10 的预算限于化学块最多 4 个、每种胶片的胶片块最多 6 个（取 SVD 最弱方向上分量最大的参数，以及报告中作为结论的参数；其余参数用 §7.8 的协方差区间与下面的先验主导判据）：在先验盒内取 8 个格点，同一块内其余参数重新优化，另一块保持在阶段 E 的不动点（剖面是"以另一块为条件"的，报告如此标注）；剖面的 χ² 也只用数据残差；Δχ² 在两侧都越过 3.84 的给出 95 % 区间；任一侧在先验盒内 Δχ² < 1 的判为实际不可辨识。先验主导：有正态或对数正态先验的参数，若只用数据的标准差 `σ_data`（数据残差的 (JᵀWJ)⁻¹ 对角元的平方根）> 0.5 τ_prior，标为"先验主导"，与不可辨识同样处理（阶段 F），不得写为 fitted。所有 Δχ² 按 §7.8 的 s² 缩放后再与 3.84 比较。判定规则在任何留出评估之前确定，不得事后调整（CAL §3.6 的泄漏途径）。

计数与预期（"有效观测数"按曲线的自由度估计：一条特性曲线约 6 个独立数——基雾、速度、趾部形状、斜率、肩部、显影时间依赖；一条 CI(t) 约 3–4 个）：

| 阶段 | 参数 | 有效独立观测 | 预期 |
|------|------|--------------|------|
| A | 4 个 Δε + 4 个干扰（含 δ_CC） | 约 11 个锚点（Carlton & Crabtree 的 6 个只给相对值） | Δε_borate、Δε_organic 由系列内部差别定；绝对的 Δε_sulfite 只靠 XTOL 工作液 8.20 与低权重的 A 部分，报告中写明；每个干扰参数恰好由一个锚点定死，没有剩余检验；Δε_amine 与 T-MAX 干扰参数无法同时定，Δε_amine 固定为 0 |
| B（每种） | 9–12 | 两条曲线，平滑，约 8–10 | 染料带位置可定；幅度比部分可定；本征吸收边（蓝端）弱 |
| C，Tri-X | 13–14 | 8 条曲线（同一胶片两种画幅，约 30）+ 6 个倒易律行 + MTF 尾段（约 5）+ A_550 罚项 | 速度、反差、灰雾好（片基参数 assumed）；**d_50 弱**（没有颗粒度目标，粒径与 η、CP 互换）；m_Ag,coat 只靠显影时间系列与邻接幅度，可能弱；τ_i 只靠 1 行高照度数据 |
| C，T-Max | 14–15 | 3 条曲线（约 18）+ 6 行 + MTF 尾段（约 7）+ 颗粒度 1 + A_550 罚项 | d_50 由颗粒度约束；t_g 靠 D-76 6/7.5/10 min 系列的速度（与共用 κ̃ 部分简并）；m_Ag,coat 同 Tri-X；τ_i 只靠 1 行 |
| D | 36 | D-76：2 条 CI(t)、4 种罐型的温度斜率（约 1 个数一种，且可能是 Kodak 惯例）、曲线；XTOL：2 条 CI(t)、J-109 三个 CI 水平；HC-110：1 条 CI(t)；T-MAX：8 条曲线（胶片参数已定，对化学约 20 个数）；合计约 80–100，但高度相关 | D-76 块与共用输运（δ_max、δ_agit）可定；T-MAX 干扰参数部分可定；**XTOL 块（k_des,pyr、K_pyr、pK_pyr 三者互换）与 HC-110 块（msds 配方 + 干扰参数）松散**；u_gel、u_X、λ0 部分简并；溶剂参数（k_diss、C_sat,ref、ν_s、k_pd）只能靠 D-76、D-76 1+1、XTOL、T-MAX、HC-110 B 的速度与最大密度差别，弱 |

m_Ag,coat 与 t_g 若在阶段 F 被判为不可辨识而固定（m_Ag,coat = 4.0、t_g = 0.12），报告必须给出敏感性：m_Ag,coat ∈ {2, 6} g/m²（每种胶片分别）与 t_g ∈ {0.08, 0.15} µm（T-Max）各做一次阶段 C–E 的重拟合，留出预测的差别并入 §8 的预测带（每次约 1 次阶段 E 的代价，§7.9）。

松散不等于预测无约束（Gutenkunst 等 2007，CAL §3.4）：T-Max 在 XTOL 中的 CI(t) 依赖的组合，与 Tri-X 在 XTOL 中的 CI(t) 所约束的组合大部分相同。留出误差是这一点的直接检验。

### 7.6 训练与留出清单

路径相对 `tools/digitize/staged/`（并入后为 `data/reference/`，§5.9）。表格行的标识写作 `<文件>#<块>/<id>/<developer>[/<稀释>][@<温度°F>]`。

**训练：胶片块**（只有 D-76 与胶片自己的数据）

| 胶片 | 文件 / 行 |
|------|-----------|
| Tri-X | `kodak/tri-x-400/char-35mm-d76-{06,08,10,12}min.csv`、`char-120-d76-{06,08,10,12}min.csv`；`mtf-d76-large-tank.csv`；`spectral-sensitivity-d0.3.csv`、`spectral-sensitivity-d1.0.csv`（形状；5063 旧乳剂的说明见 DS §2.1）；`tables.yaml#reciprocity` 的 1/10 000、1/1000、1/100、1/10、1、10 s 行 |
| T-Max 100 | `kodak/t-max-100/char-d76-{06,07.5,10}min.csv`；`mtf-d76-small-tank.csv`；`spectral-sensitivity-d0.3.csv`、`spectral-sensitivity-d1.0.csv`；`tables.yaml#facts/rms-granularity`（8）；`tables.yaml#reciprocity` 的 1/10 000、1/1000、1/100、1/10、1、10 s 行 |

**训练：化学块**（上面的 D-76 曲线也在内）

| 组合 | 文件 / 行 |
|------|-----------|
| Tri-X×D-76 | `ci-small-tank-d76.csv`、`ci-small-tank-d76-1-1.csv`；`tables.yaml#development_times/normal-small-tank/d-76`：原液 65/70/72/75 °F 入拟合，原液与 1+1 的 68 °F 两格不入（与 CI 曲线重复），1+1 的其余温度归内部检查；`normal-large-tank/d-76`：原液全部温度，1+1 只取 68 °F |
| Tri-X×XTOL | `ci-small-tank-xtol.csv`、`ci-small-tank-xtol-1-1.csv`；`normal-small-tank/xtol`：原液 65/70/72/75 °F 入拟合，原液与 1+1 的 68 °F 不入（与 CI 曲线重复），1+1 的其余温度归内部检查；`normal-large-tank/xtol`：原液全部温度，1+1 只取 68 °F；`xtol_contrast_index_times`（J-109 表 1 的 Tri-X 行）只取原液 CI 0.72 与 0.82 的 65/70/75 °F 格（CI 0.82 的 65 °F 为"—"）。J-109 的 CI 0.56 行与 F-4017 小罐 XTOL 行逐格相同（暂存 `developers/README.md` 的交叉核对），不重复入拟合；68 °F 的各格都能从 CI(t) 曲线（跨 CI 0.34–0.98）读出，也不入 |
| Tri-X×HC-110 | `ci-small-tank-hc110-b.csv`；`normal-small-tank/hc-110`（去掉 68 °F）；`normal-large-tank/hc-110` |
| Tri-X×T-MAX Dev | `char-35mm-tmax-dev-{06,07,09,11}min.csv`、`char-120-tmax-dev-{06,07,09,11}min.csv`；`ci-small-tank-tmax-dev.csv`；`normal-small-tank/t-max-dev`（去掉 68 °F） |
| T-Max×D-76 | `ci-small-tank-d76.csv`、`ci-small-tank-d76-1-1.csv`、`ci-large-tank-d76.csv`；`normal-small-tank-rolls/d-76`：原液 65/70/72/75 °F 入拟合，原液与 1+1 的 68 °F 不入（与 CI 曲线重复），1+1 的其余温度归内部检查；`normal-large-tank-rolls/d-76`（68 °F 一格不入：与 `ci-large-tank-d76.csv` 重复）、`normal-large-tank-sheets/d-76`、`normal-tray-sheets/d-76`、`normal-rotary-tube-rolls/d-76`、`normal-rotary-tube-sheets/d-76`（存在的行）。`development_time_adjustment_factors` 的 `[d-76]` 块归内部检查（§7.2） |

页片行（训练与留出都有）：按 `film_format: "135"` 的片基参数模拟（`assumed`；T-Max 页片实为 7 mil ESTAR 片基，OPT §1.2）。这些行只按表格格的 CI 评分，在均匀条带上片基主要经光晕返回 f_H 平移 log H（另有返回光在深度上分布不同的小量），CI 对 log H 的平移不变，所以片基的选择对评分的影响可忽略。

**留出组合**（不进入任何拟合、任何设计选择；冻结后只评估一次）

| 组合 | 文件 / 行 |
|------|-----------|
| T-Max×XTOL | `kodak/t-max-100/ci-small-tank-xtol.csv`、`ci-small-tank-xtol-1-1.csv`、`ci-large-tank-xtol.csv`、`ci-rotary-tube-xtol.csv`、`ci-rotary-tube-xtol-1-1.csv`、`ci-tray-xtol.csv`、`ci-tray-xtol-1-1.csv`；`tables.yaml` 中所有 `developer: xtol` 的 `development_times` 行；`xtol_contrast_index_times`（J-109 表 1 的 T-Max 行，含 CI 0.52/0.62）；`developers/xtol.yaml` 中 `material: kodak/t-max-100` 的行 |
| T-Max×HC-110 | `ci-small-tank-hc110-b.csv`、`ci-large-tank-hc110-b.csv`；`tables.yaml` 中所有 `developer: hc-110` 的行；`development_time_adjustment_factors` 的 `[hc-110]` 块；`developers/hc-110.yaml` 中 T-Max 的行 |
| T-Max×T-MAX Dev | `char-tmax-dev-{06,07,10,12}min.csv`、`ci-small-tank-tmax-dev.csv`；`tables.yaml` 中所有 `developer: t-max-dev` 的行；`development_time_adjustment_factors` 的 `[t-max-dev, t-max-rs]` 块；`developers/t-max-dev.yaml` 中 T-Max 的行 |

**留出的注册粒度。** 三个留出组合在暂存文件中的每一行与每个子项都进 `HELD_OUT`，不论评估与否；每项带 `evaluated` 与原因，冻结前固定。混合了训练与留出内容的条目按子项给标识，例如 `kodak/t-max-100/tables.yaml#facts/base-plus-fog/char-tmax-dev-06min`（留出）、`…/char-d76-06min`（训练，只作基雾上界的参考）、`…/char-tmax-rs-08min`（不用）。`tables.yaml#facts/characteristic-vs-contrast-index-curves` 整项归留出（不评估：设计期参考，含 T-MAX 的数字）；训练用的 σ_floor 由代码从两份训练 CSV 重算（§7.2），不读这一项。

| 不评估的留出项（`evaluated: false`） | 原因 |
|------|------|
| `developers/xtol.yaml` 的 `table2-rotary-tube-seasoned`、`table3-large-tank-seasoned`、`table6-rotary-tube-sheets-replenished`、`table7-large-tank-seasoned-sheets` 中 T-Max 的行 | 陈液或补充液（§1.2 非目标）；表 3、7 另为氮气脉冲搅动，D1 无此搅动类 |
| T-Max `tables.yaml#push_processing` 中 xtol、hc-110、t-max-dev 的行 | 推冲，无 CI 目标 |
| T-Max `tables.yaml#rack_and_tank`、`#replenishment` 中 xtol 的项 | 机器冲洗的时间区间；补充（§1.2） |
| `development_time_adjustment_factors` 的 `[hc-110]`、`[t-max-dev, t-max-rs]` 块的反差列（含 '*' 列） | 没有预注册的 CI 目标（§7.2） |

评估方式（预注册）：表格格一律按"表中时刻的 CI"（§7.2 的残差形式）；温度系数块的"Normal Contrast"列换算成时刻 `t = f × t_std`，t_std 为同一出版物同显影液小罐表在该块标准温度 T_std 的时间（`[hc-110]` 20 °C，`[t-max-dev, t-max-rs]` 24 °C），再按表格格评分。T-Max×T-MAX Dev 的 CI 只对 `ci-small-tank-tmax-dev.csv` 与表格格评分；四条特性曲线只按密度、log H_m*、趾部与左端密度评分，不从曲线读 CI（F-4016 两张图互相矛盾，§7.2 的 σ_floor）。

**训练组合内部的检查**（不入拟合；随训练报告给出，作为"开发期检查"而非主指标；用于设计决定时必须写入标定记录的 `notes`）

| 检查 | 内容 |
|------|------|
| 1:1 温度行（输运检查） | Tri-X 的 D-76 1+1 与 XTOL 1+1 在 65/70/72/75 °F 的小罐与大罐行；J-109 表 1 Tri-X 1+1 的 CI 0.72、0.82 在 70/75 °F 的格（其 CI 0.56 格与 F-4017 的 XTOL 1+1 行相同，只算一次）；T-Max 的 D-76 1+1 小罐 65/70/72/75 °F。Kodak 的 1:1 表观活化能 ≈ 41–43 kJ/mol，原液 ≈ 66–68（KIN §1.6）：模型应由边界层输运自然给出较低的斜率；不符时报告为"模型或 Kodak 惯例"，不判为物理失败 |
| 倒易律未拟合行（非盲） | Tri-X 1/100 000 s 与 100 s（含显影修正 +20 %、−30 %）；T-Max 100 s。阈值 ≤ ½ 档（§8.3 内部检查列）。这些行在设计期已被看过（下文"设计期暴露"），不计入任何组合状态；Tri-X 的表对 TRI-X 320 与 400 只印一次（暂存 `applies_to`） |
| 倒易律的显影修正（Tri-X 拟合行） | 1/10 000、1/1000、1、10 s 行的 +15 %、+10 %、−10 %、−20 %：按 §6.8 的方向与幅度规则（±10 个百分点）检查；设计期只用过它们的方向 |
| MTF 峰段（横向化学） | 两种胶片 MTF 的峰段（§6.7 第 8 步）：百分点差的 RMS 与峰值、峰值频率；不拟合，D1 唯一的邻接效应观测 |
| T-Max 光谱绝对水平 | §6.9 第 6 步：D = 0.3 与 1.0 的绝对 log S，450–600 nm 平均差；Tri-X 不做（旧乳剂，水平不一致） |
| F-4016 `[d-76]` 温度系数 | Normal 列按 `t = f × 6½ min`（小罐 20 °C）换算后按表格格的 CI 评分；系数与时间表本身的差（RMS 3.9 %、最大 7.1 %，§12.10）同时列出 |
| 跨显影液颗粒度 | Tri-X 在 HC-110 B 下的 RMS 17（不计入 Tri-X×HC-110 的状态） |
| 容量一致性 | D-76：每 4 卷/加仑（约 1 卷/L）后时间 +15 %；D-76 1:1 在 16 oz 罐中两卷 +10 %；XTOL：每卷至少 100 mL 原液；XTOL 每升 15 卷、两次 +15 %（KIN §4.2、§7 第 8 条；依赖假设的平均负片显影比例，只作一致性检查） |
| 格式一致性 | Tri-X 35 mm 与 120 的 D-76 曲线：同一组胶片参数同时拟合两种画幅，片基参数（t_b、D_g、D_back）都是 assumed（§7.1），残差的格式差异作检查（DS §7 第 3 条） |

**D1 不用**：所有 T-MAX RS 数据（`char-tmax-rs-*`、`ci-*-tmax-rs.csv`、`developers/t-max-rs.yaml`、表格中的 `t-max-rs` 行）；训练组合的推冲表（留出组合的推冲行归 `HELD_OUT`、不评估，见上）；J-78 中 `material: null` 的旧产品名行；滤镜系数与日光曝光表。

**显影液手册中的 Tri-X 行**（`staged/developers/{d-76,xtol,hc-110,t-max-dev,t-max-rs}.yaml`）：J-109 表 1 的 Tri-X 行已经以 `xtol_contrast_index_times` 收在胶片的 `tables.yaml` 里，按上表的规则使用，这里的副本不再用；J-109 表 2–8、J-24 与 J-86 的 400TX 行中**新鲜药液**的组，凡是与 F-4017 同罐型、同温度的格，暂存 README 已核对为逐格相同，一律"不用（F-4017 的重复）"；J-109 的陈液与补充液组（`table2-rotary-tube-seasoned`、`table3-large-tank-seasoned`）README 没有核对，也不是 F-4017 的重复，归"不用（陈液/补充液，§1.2 非目标）"。F-4017 没有的罐型只有 J-109 表 3、7 的氮气脉冲大罐（每 10 s 冲 2 s）：D1 的搅动模型没有脉冲气体搅动这一类，归"不用"，留到将来加入该搅动方式时再定。J-78 的行全部是 `material: null` 的旧产品名，归"不用"。J-86 的机器冲洗时间区间（"6 to 8"）不用。这些文件中的容量、补充、稀释定义等事实只用于内部的容量一致性检查。`test_split_registry_complete`（§10 D1e）按本段与上表逐项分类。

设计期暴露：

- §12.2、§12.3 的可行性计算用了两张倒易律表的全部行（包括上面的未拟合行）来确认链的结构能否达到表值，`S_eff` 闭合是在看到 Tri-X 1/100 000 s 那一行之后选的；LAT §1.4 用 Tri-X 100 s 一行论证需要 τ₀ 的分布（决策 4(2) 的依据）。所以倒易律未拟合行检验的是拟合常数的外推，不检验模型结构，也不计入任何状态；报告中如实注明。τ̄₀、σ_τ 的界只取自 LAT（§2.2.7），不取 §12.2 的拟合。训练用 1/10 000 s 行，偏离了 LAT §5(a)（只拟合 1/1000–10 s），因为它是 τ_i 的唯一高照度约束（§7.5）。
- 写作中为了解暂存文件的结构，`tables.yaml` 的若干留出部分被显示在屏幕上：T-Max 100 的 `xtol_contrast_index_times`（J-109 表 1 的 T-Max 行）与 `development_time_adjustment_factors` 的全部块（含 `[hc-110]` 与 `[t-max-dev, t-max-rs]`）。
- 规格所依据的研究报告与部分暂存表头也含留出数字，作者在设定阈值与权重之前都已读过：DS §3.1（T-Max 在 XTOL、HC-110、T-MAX 中各罐型的 CI(t) 跨度）；DS §3.2（20 °C 的 T-Max 各显影液行、HC-110 系数例与 T-Max×HC-110 的系数–表格比较）；DS §4（J-109 表 1 的 T-Max 行全文、J-24 的 T-Max 稀释 B 行，虽然 §7.7 第 1 条禁止查看）；CAL §4.2 与 §9 第 4 条（留出曲线的斜率）；暂存 `developers/README.md` 中 F-4016 对 J-109 的一个 T-Max XTOL 格；`t-max-100/tables.yaml#facts/base-plus-fog` 的 T-MAX 左端密度；训练文件 `ci-small-tank-d76.csv`、`char-d76-07.5min.csv` 表头 `ci_cross_check` 中对"每一条特性曲线"的概括（含留出的 T-MAX 曲线）。
- 原来按 CAL §4.2 的"局部 dCI/dt 下的等效时间"判定与"T-Max×HC-110 改判表中时刻的 CI"用了研究组从留出曲线算出的斜率；本规格已改为对全部格统一用"表中时刻的 CI"（§7.2、§8.3），不再需要任何斜率。原 F-4016 系数阈值 15 % 是在看过 `[hc-110]` 块之后定的，现改为只由训练数据导出（§8.3、§12.10）。
- 修订时为核对阈值与 F-4016 自相矛盾的相容性，两条留出的 T-MAX 特性曲线（10、12 min）被打开计算过一次；§12.10 最终只用训练曲线与 fact 中已印出的差值，规格中没有任何数字依赖那次计算，阈值没有改变，σ_floor 只由训练对得出。
- 因此留出组合对代码是盲的，对写规格与定阈值的人不是。报告的留出部分必须带一句："阈值与权重设定之前，作者已看过部分留出行（§7.6 设计期暴露）。"并入 `data/reference/` 时，训练文件表头中由留出数据得出的统计（`ci_cross_check` 的概括句）删去，或把该表头字段归入 `HELD_OUT`。

### 7.7 泄漏规则

由代码强制：

- `validate/splits.py` 定义 `TRAINING`、`WITHIN_CHECKS`、`HELD_OUT`、`UNUSED` 四个互不相交的注册表（按文件路径、表格行与子项标识，§7.6），另有 `HELD_OUT_EVALUATED ⊂ HELD_OUT` 与每个不评估项的原因。标定的数据加载器遇到 `HELD_OUT` 中任何一项抛出 `LeakageError`（`test_calibration_loader_refuses_held_out`）。`calibrate` 包只能通过 `TrainingLibrary`（Library 的过滤视图）读参考数据：`reference_curve`、`reference_tables`、`developer_facts` 对 `HELD_OUT` 的任何访问（含整个文件中混有留出子项时返回的部分）都抛 `LeakageError`，返回的表格已去掉留出子项（`test_training_library_filters_held_out`）。
- 留出评估函数要求标定记录 `frozen: true`，评估后写 `held_out_evaluated: true`；再次评估必须显式给 `after_unblinding=True`，输出标为"解盲后"（`test_held_out_requires_frozen_record`、`test_second_held_out_run_labelled`）。

由流程约束（CAL §3.6 的泄漏途径，逐条适用）：

1. 任何地方都不得使用 J-109 表 1 或 J-24 中 T-Max 100 的行，包括开发期的"随便看看"。
2. 不得在看到留出残差后调整共用物理：搅动/边界层模型、槽体积假设、HC-110 的 msds 配方、XTOL 的专利配方、显影剂 pKa 的选择。
3. 留出的 T-Max×XTOL 每一种搅动制度（小罐、大罐、转筒、盘）在 T-Max×D-76 的训练数据中都有同制度的数据，模型里转筒与盘又是同一个连续搅动类（§2.4.8）：留出检验的是在同一胶片、同一搅动制度下标定后的化学迁移，没有任何留出项单独检验搅动。报告如此注明；各罐型仍分开列出，只是为了看误差的结构。
4. 固定哪些参数（§7.5）不得参考留出结果。
5. 不得在曲线内随机留点；留出以整条曲线、整个温度、整个稀释为单位。
6. 温度依赖不得用 F-4016 的温度系数拟合（HC-110 与 T-MAX 块属于留出；D-76 块只作内部检查）。
7. CI 的作图法、密度计的光谱类型在看到留出结果前固定（本规格已固定：H-740 法、视觉密度、ISO 感光度用 Type 2 近似）。
8. 模型的任何后续修改都是新模型版本，其留出数字标为"解盲后"。

### 7.8 结果与不确定度的写回

- 协方差按块计算：`C_b = s_b² (J_bᵀ W J_b)⁻¹`，J_b 为块 b（胶片块或化学块）对其自身目标函数数据残差的雅可比，另一块固定在阶段 E 的不动点；`s_b² = χ²_b,min/(N_eff,b − p_b)`（变换空间，CAL §3.5），N_eff 按 §7.2 的独立数计（每条特性曲线 6、每条 CI(t) 3、每个表格类按其 σ 放大后的等效格数），不按重采样点数。阶段 E 是两个目标的不动点，不是一个联合最优（§11 第 3 条），所以跨块协方差无定义：参数表把这些协方差与剖面区间标为"以另一块为条件"。剖面法给出的参数用剖面区间。
- 每个拟合值写回所在数据文件：`source: fitted`，`value`，`uncertainty`（原单位的一个标准差，仅线性变换时；对数变换的参数只在 `fit.sigma_transformed` 中给出，`uncertainty` 留空），`fit` 记录（§5.1）。不可辨识或先验主导、按阶段 F 固定于阶段 E 值的参数改为 `source: assumed`，`note` 写明。
- 标定记录 `data/calibration/<run_id>.yaml` 保存全过程（§5.7）。运行记录（manifest）的 `calibration_runs` 列出所用数据文件中出现的全部 `fit.run_id`。
- 参数不确定度是模型偏差存在时的下界（Kennedy & O'Hagan 2001，CAL §3.5），报告以"在模型 D1x 下标定"表述；可信度的主指标是留出预测误差。
- 预测带，两种并列报告（种子 0，各 500 组，重算留出预测，给 2.5 % 与 97.5 % 分位）：(a) 只在可辨识方向上从 N(θ̂, C_b) 抽样（各块独立抽）；(b) 在 (a) 之外，阶段 F 固定的参数（不可辨识或先验主导，如 XTOL 块的 k_des,pyr、K_pyr、pK_pyr 与 HC-110 的干扰参数）从各自的先验盒抽样（正态/对数正态先验按先验分布，均匀盒按均匀分布）。(a) 是"条件"带，(b) 才反映固定参数对留出预测的影响；留出误差落在 (a) 之外而在 (b) 之内时，报告为"受不可辨识参数限制"。

### 7.9 计算预算

见 §3.10（10 进程，有效并行按 ×8 估计；实测只有 ×6 时各项约乘 1.3）：阶段 D 的目标函数一次约 23 s，一次 TRF 迭代约 8 min，多起点约 16 h；每种胶片的阶段 C 约 1.5 h；阶段 E ≤ 约 10 h；剖面似然约 27 h；合计约 2–2.5 天，一次性运行。若 m_Ag,coat、t_g 在阶段 F 被固定，§7.5 的敏感性重拟合最多 6 次（m_Ag,coat 两值 × 两种胶片，t_g 两值），每次从阶段 E 的结果起做一轮 C–E（单起点，约 10 h），另加最多约 2.5 天。旧版的"目标函数约 1 s、半天"假设了跨条件拼批与每次运行 274 步，两者都不成立（§3.4、§3.10）。D1c 实测条带耗时后、D1e 开始前按实测改写本节。

---

## 8 验证报告与预注册阈值

### 8.1 报告内容

`atlas validate` 为每种胶片生成一个目录：`report.md`（中文说明，标识符英文）、每张图一个 SVG（matplotlib，D1e 引入的 `report` 依赖组）、每个指标一个 CSV。内容：

1. 抬头：模型版本、代码提交、全部数据文件哈希、标定记录 run_id 与 `frozen`/`held_out_evaluated` 状态、离散化档位。
2. 胶片×显影液组合状态表（14 个组合，§8.2）。
3. **主指标：留出预测误差**。T-Max×XTOL：各罐型 CI(t) 的逐点误差与 RMS、全部评估表格格（含 J-109 表 1 各 CI 水平×温度×稀释）的"表中时刻的 CI"误差；T-Max×HC-110：CI(t) 与表格格（含 F-4016 系数块的 Normal 列，§7.6）的 CI 误差；T-Max×T-MAX Dev：四条曲线的 RMS 与最大误差、左端密度、log H_m*、趾部密度，CI(t) 与表格格的 CI 误差。每项带两种预测带（§7.8）。
   - **每个留出组合一个标题数字**（预注册）：T-Max×XTOL 与 T-Max×HC-110 为该组合全部评估 CI 观测（CI(t) 的点与表格格合并）的 RMS(CI_sim − CI_ref)；T-Max×T-MAX Dev 为四条曲线合并的密度 RMS（§8.3 的窗口）。报告摘要只引用这三个数字与组合状态。
   - **朴素基线**（只用训练数据，冻结前写好代码）：CI_base(t; T-Max×X, 罐型, T, 稀释) = CI_ref,T-Max×D-76(t · r)，r = t_TriX×D-76(20 °C, 同罐型, 同稀释) / t_TriX×X(T, 同罐型, 同稀释)，t 取 Tri-X 训练表格中 CI 0.56 的时间，T-Max×D-76 用同罐型的参考 CI(t)（没有同罐型的 Tri-X 数据时用小罐）；T-MAX 的曲线基线取 T-Max×D-76 三条曲线在 t·r 处按时间线性插值（逐 log H）。三个标题数字都与基线的同一数字并列；模型不优于基线时如实写明"机理化学没有显示出超过时间比换算的预测能力"。
4. 训练拟合质量，逐组合、逐指标对照阈值。
5. 训练组合内部的检查（§7.6）。
6. 像质：MTF（两种胶片，D-76；尾段为训练拟合，峰段为横向化学的内部检查），T-Max 颗粒度（训练），Tri-X 颗粒度（跨显影液检查），统计/显微模式一致性。
7. 光谱感光度形状、T-Max 光谱绝对水平（内部检查，§6.9 第 6 步）、倒易律曲线（两个判据）与 Kodak 表。
8. 一致性检查：ISO 感光度对标称感光度（D-76 的 log H_m 由拟合的同一批曲线决定，这只检查指标实现，不检验绝对感光度假设；D1 中对感光度唯一非循环的证据是 T-Max 光谱绝对水平与留出的 T-Max×T-MAX log H_m*）、容量、格式。
9. 不确定度：参数表（值、σ 或区间、方法、是否可辨识、是否在界上）与各块奇异值谱。
10. 数值：图像档对参考档、标定档对参考档的偏差（§3.8、§3.9 的测试结果）。
11. 已知的不物理之处（§9）。
12. 仅预测的组合：曲线与 CI(t) 的预测，不带参考。

### 8.2 组合状态

| 状态 | 规则 |
|------|------|
| `validated` | 该组合所属层级（训练或留出）的全部适用指标（§8.3 的组合–指标表，逐项预注册）都在阈值内 |
| `failed` | 有参考数据，且任一适用指标超出阈值，**或任一适用指标有参考数据却算不出（NaN 或无定义，例如 ΔD = 0.80 夹不住、CI 作图无解）**；算不出的指标从不跳过，报告列出原因 |
| `predicted` | D1 没有该组合的参考数据：Tri-X 与 T-Max 各自配 `rodinal`、`hc-110-2019`、`d-76-published`（Kodak 的数据用的是成品 D-76） |

训练组合的 `validated` 写作"validated（训练拟合）"，与留出组合的"validated（预测）"分开显示。状态连同层级（`training`/`held_out`）、run_id、阈值集哈希与是否"解盲后"写入标定记录的 `pair_statuses`（§5.7），图集从那里读取；D0 的 `VerificationStatus` 枚举不变，层级由 `PairStatus.tier` 给出。内部检查（§7.6）不计入任何状态。

### 8.3 预注册阈值（暂定；所有者在训练报告之后、冻结之前确认或修改，§8.4）

与 Kodak 曲线比较的全部指标都用 §6.3 第 2 步的 D_min* 算子（两侧同一 x_left）。

| 指标 | 训练组合 | 留出组合 | 依据 |
|------|----------|----------|------|
| 特性曲线密度 RMS（D_ref = D_min*+0.1 至最后一点，两侧在参考曲线的 log H 上取值） | ≤ 0.03 D | T-Max×T-MAX：≤ 0.06 D | CAL §4.2；留出值按"留出阈值取训练阈值的 2 倍"的规则（本规格提出，与 CAL 对速度的 0.05→0.10 一致） |
| 特性曲线最大误差 | ≤ 0.08 D | ≤ 0.15 D | 同上 |
| 左端密度 \|D_sim(x_left) − D_ref(x_left)\| | ≤ 0.03 D | ≤ 0.05 D | CAL §4.2 的基雾阈值；Kodak 不给基雾，改比左端（暂存 `facts/base-plus-fog`） |
| 基雾上界（单侧） | D_sim(零曝光片) ≤ D_ref(x_left) + 0.02 | 同左 | 左端是基雾的上界；0.02 ≈ 数字化 y 误差 0.016 |
| log H_m*（Type 2 近似，D_bf = D_min*） | ≤ 0.05 | ≤ 0.10 | CAL §4.2（ISO 6 §6.3） |
| 趾部：参考曲线的 log H_m* + 0.3、+ 0.6 处的 D（两侧同一 log H，纯密度比较） | — | ≤ 0.06 D | 本规格（留出曲线检验的"趾部"） |
| CI(t) 曲线上给定时刻的 CI | RMS ≤ 0.03，max ≤ 0.05 | 每点 ≤ 0.07 | CAL §4.2 |
| 表格格的 CI（§7.2：CI_sim(t_cell) − CI_aim，含 J-109 表 1 与 F-4016 系数块 Normal 列换算出的格） | 每个训练组合 RMS ≤ 0.03，max ≤ 0.05 | 每格 ≤ 0.07 | CAL §4.2 的判据直接写成 CI 单位，不需要 dCI/dt；取代原"到达 CI 的时间"与"F-4016 HC-110 温度系数 ≤ 15 %"两行（后者是在看过 `[hc-110]` 块后定的，已撤销；训练侧的 `[d-76]` 块与时间表之差最大 7.1 %，约 0.04 CI，§12.10，在 0.07 之内） |
| RMS 颗粒度 | T-Max：\|G − 8\| ≤ 1 | — | 本规格（Kodak 给整数；±1 单位约 12 %） |
| MTF 尾段（§6.7 第 8 步） | RMS ≤ 3 个百分点 | — | 本规格（OPT §2.4 读数误差 ±1.5，Tri-X 尾部 ±3） |
| 光谱感光度形状（420–650 nm，且在参考跨度内） | RMS ≤ 0.05 log S | — | 本规格 |
| 倒易律拟合行 | ≤ ⅓ 档 | — | 本规格（Kodak 以 ⅓ 档为分辨率） |

内部检查的阈值（只报告与标出，不计入任何状态，§7.6）：倒易律未拟合行 ≤ ½ 档；Tri-X×HC-110 B 颗粒度 \|G − 17\| ≤ 2；MTF 峰段 RMS > 5 个百分点时标出；T-Max 光谱绝对水平 \|平均差\| > 0.2 log S 时标出；1:1 温度行、`[d-76]` 系数、Tri-X 倒易律显影修正（±10 个百分点）按 §6.8 与 §7.6。

**T-Max 100 与 F-4016 的自相矛盾。** F-4016 的特性曲线比同书的 CI(t) 曲线高 0.035–0.065 CI（§7.2）。所以 T-Max 的 CI 指标只对 CI(t) 曲线与表格格评分，特性曲线只按密度、左端、log H_m*、趾部评分，不从曲线读 CI 去比。§12.10 只用训练曲线核对了两图之间各承担一半差的折中：差 0.036（训练对）时 CI 误差 0.018、密度 RMS 0.013–0.015；差 0.065（fact 中最大的差）时 CI 误差 0.033、密度 RMS 0.024–0.027。都在上表相应列的阈值内，留给模型误差的余量约 0.015 D 与 0.01–0.04 CI，所以上表不要求两图之间的贴合比数据本身更紧。

**组合–指标表（预注册）。** "适用"指下表列出、且参考数据本身算得出的指标（例如参考显影系列夹不住 ΔD = 0.80 时，log H_m* 不适用；这一点冻结前由参考数据确定并写入 `validate/splits.py` 的 `PAIR_METRICS`）。适用指标在模拟侧算不出时判 `failed`（§8.2）。

| 组合 | 层级 | 适用指标 |
|------|------|----------|
| Tri-X×D-76 | 训练 | 8 条 D-76 曲线的密度 RMS 与最大误差、左端密度、基雾上界；两个画幅显影系列的 log H_m*；CI(t) 2 条；表格格；胶片层面：MTF 尾段、光谱形状 2 条、倒易律拟合行 |
| Tri-X×XTOL | 训练 | CI(t) 2 条；表格格（含 J-109 表 1 的训练格） |
| Tri-X×HC-110 | 训练 | CI(t) 1 条；表格格 |
| Tri-X×T-MAX Dev | 训练 | 8 条曲线的密度 RMS 与最大误差、左端密度、基雾上界；两个显影系列的 log H_m*；CI(t) 1 条；表格格 |
| T-Max×D-76 | 训练 | 3 条曲线的密度 RMS 与最大误差、左端密度、基雾上界；显影系列的 log H_m*；CI(t) 3 条；表格格；胶片层面：MTF 尾段、光谱形状 2 条、RMS 颗粒度 8、倒易律拟合行 |
| T-Max×XTOL | 留出 | CI(t) 7 条；评估的表格格（§7.6） |
| T-Max×HC-110 | 留出 | CI(t) 2 条；评估的表格格（含 `[hc-110]` Normal 列） |
| T-Max×T-MAX Dev | 留出 | 4 条曲线的密度 RMS 与最大误差、左端密度、基雾上界、趾部；显影系列的 log H_m*；CI(t) 1 条；评估的表格格（含 `[t-max-dev, t-max-rs]` Normal 列） |
| 两种胶片 × `rodinal`、`hc-110-2019`、`d-76-published` | 仅预测 | — |

冻结前，用一组已知 θ 生成的合成参考（结构与上表各项相同，不读 Kodak 的留出文件）跑一遍完整的留出评估流程，要求每个适用指标都有限（`test_held_out_pipeline_finite_on_synthetic`）。
| 统计/显微模式颗粒度 | σ_D 比 ∈ [0.94, 1.06]（D = 0.3、1.0、2.0） | — | 本规格（两个估计各约 1.1–1.7 % 标准误差的 2.5 倍） |

### 8.4 一次性留出评估

1. 训练与开发完成，所有 §3 的数值测试通过；**D1f 的跨模式测试也已通过**（D1f 复用潜影表、显影引擎、遮盖力与随机数流，并首次在 4 µm 化学单元上运行显影引擎，§2.7.2；它在共用代码里发现的问题若在留出评估之后才修，标题数字就永久变成"解盲后"）。所以唯一一次留出评估放在 D1f 完成之后（§10）。
2. 训练报告交所有者；所有者确认或修改 §8.3 的阈值与组合–指标表（决策 3 的"第一份报告"即训练报告）。确认之后阈值不再改动。
3. 冻结：代码提交（工作树必须干净）、数据文件哈希、标定记录 `frozen: true`；同时写入 `thresholds`（§8.3 全表与组合–指标表）、`split_registry_sha256`（`TRAINING`/`WITHIN_CHECKS`/`HELD_OUT`/`UNUSED`/`HELD_OUT_EVALUATED` 与原因的规范化 JSON 的哈希）、`metrics_code_sha256`（`validate/metrics.py`、`validate/reference.py`、`validate/splits.py` 的内容哈希）。
4. `atlas validate --held-out` 运行一次：工作树不干净、提交或上述三个哈希与记录不一致时拒绝运行（`test_held_out_refuses_mismatched_freeze`）。写出报告，把 `held_out_evaluated` 置为 true，并把各组合状态写入 `pair_statuses`。
5. 此后的任何改动（模型、代码或阈值）都是新版本，再评估时报告标"解盲后"。冻结阈值下的状态永久保留；若之后修改阈值，新阈值下的状态并列显示，标"事后阈值"。

---

## 9 已知的不物理之处与风险

| # | 项 | 影响 | 如何发现 |
|---|----|------|----------|
| 1 | RGB → 光谱还原（D0） | 同色异谱 | 不可检测；光谱输入可绕开 |
| 2 | 活度：SIT 形式，A(T) 与 1.5 为 recalled 值，Δε 只由约 11 个 pH 锚点拟合 | 缓冲容量在 D-76 的 pH 下对活度敏感（总 β 92 → 33 mM/pH，其中亚硫酸盐 81 → 17，§2.4.12），直接决定显影中的 pH 降 | 锚点残差；Carlton & Crabtree 硼砂系列的趋势；1:1 温度检查（输运与 pH 共同作用） |
| 3 | 成品 D-76 与 CIS-61 pH 的张力：只保留 D(I) 项时，加入 SDS 范围内的 B₂O₃ 后 pH 为 8.17–8.32，低于 CIS-61 的 8.5（§12.4） | B₂O₃ 可能被推到界上，或 Δε 被这一个锚点拉偏 | 阶段 A 后看 B₂O₃ 是否在界上（`at_bound`）、Δε 是否与 Carlton & Crabtree 系列冲突；冲突时记录并在报告中列出 |
| 4 | 配方不确定：旧 HC-110 的 SDS 数字全部是 secondary 且未核实（DEV §3.1）；T-MAX Developer 世代 I 的"16 % SO₂ 原料"读法是假设（TMX §3.3）；XTOL 商品与专利例 1 不同（DEV §2.2）；Rodinal 的 KOH 与亚硫酸盐未知 | 这些显影液的化学常数被配方误差污染；HC-110 留出检验同时检验了干扰参数的迁移 | Tri-X×HC-110、Tri-X×T-MAX 的训练残差；留出结果；报告中 HC-110 与 T-MAX 的结果与配方可信度并列 |
| 5 | Kodak 的时间–温度表可能是惯例（D-76 与 XTOL 斜率相同，KIN §1.6） | 温度依赖实际上每个稀释类只有一个观测 | 1:1 温度检查；XTOL/HC-110 原液非 20 °C 行降权 |
| 6 | 引发风险可分离为化学因子 × 团簇因子 | 若溴化物对不同大小的潜影中心作用不同，趾部随显影时间与显影液的变化会系统偏离 | D-76 与 T-MAX 两组特性曲线趾部残差的时间趋势；留出的 T-Max×T-MAX 趾部 |
| 7 | 一阶生长律 | 单颗粒 S 形由引发分布代替自催化 | 参考档的"零级生长"变体（每档 8 个银量格的群体平衡）拟合 Tri-X D-76 时间系列；χ² 改善 > 10 时在训练数据上改结构（允许，属于开发期决定，记入标定记录） |
| 8 | 遮盖力在显影中不变、按粒径的幂律缩放（β_CP assumed） | 肩部形状、最大密度 | 肩部残差；T-Max 颗粒度与曲线的联合拟合是否一致 |
| 9 | 亚硫酸盐溶剂作用简化：静态的 h_n 衰减、C_sat 幂律、配体数 2；ν_s、ν_b 的唯一测量在 0.6–5 mM（30 °C），D1 用到 40–800 mM，外推一到两个数量级（0.79 M 处 ν_s = 0.6 与 1 相差 7.6 倍），所以 ν_s 拟合、ν_b 假设 | 溶剂型与非溶剂型显影液的感光度差、最大密度差；D-76 对 HC-110 B 的溶解力比（ν_s = 0.6 时 4.7，ν_s = 1 时 13） | D-76 与 T-MAX 的速度差残差；D-76 1:1 的速度（只有 CI 数据，弱） |
| 10 | ETA 的表面周转机理；脱附自由基的溶液再生忽略 | 超加和的强度 | 没有单剂对照数据，D1 基本不可检验；以 D-76 与 XTOL 的拟合质量间接判断 |
| 11 | 每组分一个扩散系数，不用 Nernst–Planck | 液接电位对离子的输运 | 参考档的 Nernst–Planck 变体（D1c 可选测试 `test_nernst_planck_effect`）：D-76 条带上 \|ΔCI\| > 0.01 时升级 |
| 12 | 凝胶–液体分配系数 1（忽略明胶电荷的 Donnan 排斥） | 阴离子形态（HQ²⁻）在凝胶中的有效浓度 | 不可直接检测，被 K、k 吸收 |
| 13 | 溶胀用扩散系数爬升表示，几何从一开始就是溶胀后的 | 短显影时间的行为 | HC-110 B 的短时间（3¾ min）行残差 |
| 14 | 两种边界层厚度的搅动模型、完全复位、大罐沥干期忽略 | 罐型之间的差别 | Tri-X D-76 大罐曲线对小罐 CI(t)；T-Max D-76 各罐型行 |
| 15 | 光晕只算一次返回；135 灰片基 0.10 的解释未核实；120 背涂层密度假设；没有光晕参考数据 | 高光周围的光晕强度未经验证 | 只有 35 mm 与 120 基雾差的格式检查；MTF 在 2.5 c/mm 处的微弱影响 |
| 16 | PSF 与深度、波长无关；两流近似；准直光入射后立即视为漫射；明胶折射率未核实 | MTF 的高频部分、深度分布 | MTF 拟合残差；光谱 D=0.3 与 D=1.0 两条曲线形状差 |
| 17 | 潜影：两位点阻塞加 S_eff 的闭合（LAT §1.3(c) 方案 (ii) 型的平均场闭合，`assumed`；ν τ_i → ∞ 时高照度失效饱和在 S_eff − 1 倍）；显影前不衰退；曝光温度固定 23 °C | 高照度与低照度失效的形状；短于 1e-5 s 的曝光不再继续失效 | 倒易律留出行（结构已在设计期暴露，§7.6） |
| 18 | Tri-X 的光谱、MTF、颗粒度数据来自旧乳剂（5063/TX） | 胶片参数混用了两代乳剂的数据 | 光谱形状推出的吸收与 D-76 速度是否一致；报告中标注 |
| 19 | T-Max 的扁平颗粒厚度没有直接数据（拟合，界 0.08–0.15 µm），没有 ECD 数据 | 粒径与遮盖力的分配；S/V（§9 第 31 条） | 颗粒度 8 与曲线的联合拟合；AR ≥ 8 约束是否起作用；t_g 的剖面似然 |
| 20 | 统计模式：泊松放置（无堆积排斥）、重叠噪声取高斯、每档一个有效半径；4 r_c > p 的档只取均值、方差全由高斯 ε 提供；像素间相关不表示（白噪声等效，孔径方差偏高约 b_ap）；物理银按比例放大盘 | 颗粒度；粗颗粒、细像素时单像素的分布形状与截断 | D1f 的跨模式测试；`meta["grain"]` 的 b_ap 与截断比例 |
| 21 | 图像档降阶（化学单元、粗深度、晕边、槽预计算） | 图像与条带的差别 | §3.8 的测试 |
| 22 | 理想漫射几何；Type 2 由五个锚点近似；ISO 镜头透射由锚点近似；曝光时间 1/50 s 为假设 | 绝对感光度 | D-76 的 ISO 速度对标称值是循环的（log H_m 来自 η_0 与吸收尺度可自由拟合的同一批 D-76 曲线），只检查指标实现，发现不了这些假设。非循环的证据只有两项：T-Max 光谱感光度的绝对水平（不拟合的内部检查，§6.9 第 6 步，检验 lux·s 到光子的换算链与 1.4 s 的失效）与留出的 T-Max×T-MAX log H_m*。中性影像下 Type 2 近似无影响 |
| 23 | 停显、定影、水洗的常数全部假设 | 只影响跳过或缩短这些步骤时的后果 | 无参考数据，不验证 |
| 24 | 空气氧化、补充液不建模 | 容量与老化 | 容量一致性检查 |
| 25 | 显影中离子强度不变 | 稀显影液（Rodinal 1+50）中的局部活度 | 参考档变 I 的对照（可选） |
| 26 | 有机显影剂电离焓取 0，非文献显影剂的 Ea 固定为 65 kJ/mol | 温度依赖 | XTOL/HC-110 温度行（弱）；1:1 检查 |
| 27 | 雾翳只来自 h_fog（无自发雾、二色雾） | 基雾随时间的增长 | Tri-X D-76 基雾在 6–12 min 内增加 ≤ 0.05（KIN §1.8） |
| 28 | 物理显影银的遮盖力比、碘化物含量为假设 | 最大密度、抑制强度 | 最大密度水平 |
| 29 | XTOL 与 HC-110 化学块松散 | 参数本身不可信，预测可能仍受约束 | SVD、剖面似然；留出预测带的宽度 |
| 30 | Kodak 曲线只"代表生产涂布"（CAL §2.3） | 可达到误差的下限 | 不可检测；阈值已考虑（σ_D 0.02–0.03） |
| 31 | 显影速率按比表面积缩放 κ_c ∝ S_c/V_c（§2.4.4，`assumed`）；扁平颗粒的 S/V ≈ 2/t_g + 4/d 主要由 t_g 决定（t_g 0.08–0.15 时 29–17 µm⁻¹，Tri-X 球形 d = 0.9 µm 时 6.7） | T-Max 相对 Tri-X 的显影速度、T-Max 留出组合的跨颗粒类型迁移 | T-Max × D-76 的 CI(t) 与 6/7.5/10 min 系列残差；报告中给出留出预测对 t_g 的敏感性 |
| 32 | 涂布银量 m_Ag,coat 无一手数据（拟合，先验 2–6 g/m²）；它不只与 CP_ref 以乘积进入 D_max，还决定 C_c,z、单位密度释放的 Br⁻/I⁻/H⁺ 与消耗的显影剂、N_tot 与 L_e | 耗尽、补偿、邻接效应幅度；不可辨识时其误差被共用化学常数（δ、K_Br、u_X）吸收并带进 T-Max 留出 | 剖面似然；不可辨识时的 {2, 6} g/m² 敏感性运行（§7.5） |

---

## 10 模块计划（D1a–D1f）

通用完成标准（每个模块都适用）：

- 本模块列出的测试全部通过，`uv run pytest` 全绿；`uv run ruff check` 与 `uv run ruff format --check` 无报错。
- `docs/physics.md` 补上本模块的章节，方程与代码一致，文中提到的测试名都存在（`test_physics_doc.py`）。
- `README.md` 版本记录加一行，version 徽章与 `darkroom_atlas.ATLAS_VERSION` 同步改为 `D1a`…`D1f`（`test_version_sync.py`）。
- 用理想材料重跑 D0 的样例，输出图像与 D0 的逐位相同（回归保护）。
- 不改 `tools/`（数字化团队负责）。提交信息不出现任何 AI 或代理的名字。

依赖：D1a → D1b → D1c → D1d → D1e → D1f（D1f 还需要 D1d 的统计模式与 D1e 的条件与指标）。

### D1a 光传输（含内存基础设施）

文件：新建 `fieldstore.py`、`blocks.py`、`emulsion/__init__.py`、`emulsion/classes.py`、`optics/__init__.py`、`optics/stack.py`、`optics/kubelka.py`、`optics/psf.py`、`optics/halation.py`、`optics/expose.py`、`models/bw.py`；修改 `light.py`、`scene/camera.py`、`spectral.py`、`_frozen.py`、`state.py`、`pipeline.py`、`units.py`、`data/models.py`（§5.4 的胶片字段与 `Quantity.bounds`、`Quantity.prior`）、`data/library.py`（`validate` 的新检查）、`models/__init__.py`。

接口：

```python
# light.py
class LightField(Protocol):
    @property
    def grid(self) -> SpectralGrid: ...
    @property
    def meta(self) -> Mapping[str, Any]: ...
    @property
    def shape_hw(self) -> tuple[int, int]: ...
    def window(self, rows: slice, cols: slice) -> "Light": ...
    def project(self, weights: np.ndarray, rows: slice | None = None,
                cols: slice | None = None) -> np.ndarray: ...          # weights (B, M) -> (h, w, M)
    def luminous_exposure_lux_s(self, time_s: float) -> np.ndarray: ...

@dataclass(frozen=True, eq=False)
class LowRankLight:
    coefficients: np.ndarray        # (H, W, R) float64, read-only, finite, >= 0
    basis: np.ndarray               # (R, B) photons um^-2 s^-1 per unit coefficient, >= 0
    grid: SpectralGrid
    meta: Mapping[str, Any] = field(default_factory=dict)
    @property
    def shape_hw(self) -> tuple[int, int]: ...
    @property
    def rank(self) -> int: ...
    def window(self, rows: slice, cols: slice) -> Light: ...
    def project(self, weights: np.ndarray, rows: slice | None = None,
                cols: slice | None = None) -> np.ndarray: ...
    def luminous_exposure_lux_s(self, time_s: float) -> np.ndarray: ...
    def filtered(self, transmittance: np.ndarray) -> "LowRankLight": ...   # (B,) only
    def to_light(self) -> Light: ...

class Light:                        # D0 类，新增两个方法
    def window(self, rows: slice, cols: slice) -> "Light": ...
    def project(self, weights: np.ndarray, rows: slice | None = None,
                cols: slice | None = None) -> np.ndarray: ...

# scene/camera.py
def scene_light(linear_rgb, *, grid: SpectralGrid = DEFAULT_GRID, illuminant: str = "D65",
                iso: float = 100.0, ev: float = 0.0, time_s: float = 1 / 125,
                white_exposure_lux_s: float | None = None,
                low_rank: bool = False) -> Light | LowRankLight: ...

# fieldstore.py
class FieldStore:
    def __init__(self, root: Path | None, *, ram_threshold_bytes: int = 256 * 2**20) -> None: ...
    def allocate(self, name: str, shape: tuple[int, ...], dtype: str = "float64") -> np.ndarray: ...
        # name is a store name such as "s3:emulsion/silver", unique per run (ValueError if reused);
        # backing file "<name with / and : -> _>-<uuid4[:12]>.bin" (§4.3)
    def freeze(self, name: str) -> np.ndarray: ...       # read-only (RAM or np.memmap mode "r")
    def release(self, name: str) -> None: ...            # deletes the backing file if any
    def descriptor(self, name: str) -> "FieldDescriptor": ...   # mapped fields only; RAM fields raise
    def close(self) -> None: ...
    def __enter__(self) -> "FieldStore": ...
    def __exit__(self, *exc: object) -> None: ...

@dataclass(frozen=True)
class FieldDescriptor:                    # picklable; what crosses process boundaries (§4.3)
    path: Path
    dtype: str
    shape: tuple[int, ...]
    offset: int = 0

def open_field(desc: FieldDescriptor, mode: Literal["r", "r+"]) -> np.ndarray: ...   # np.memmap

# blocks.py
@dataclass(frozen=True)
class Block:
    index: tuple[int, int]
    inner: tuple[slice, slice]            # absolute, in grid units (pixels or chemistry cells)
    padded_origin: tuple[int, int]        # absolute position of padded[0, 0]; may be negative
    padded_shape: tuple[int, int]         # the same for every block of a grid

def block_grid(shape: tuple[int, int], block: int, halo: int) -> tuple[Block, ...]: ...
def padded_window(source: np.ndarray, block: Block, fill: float = 0.0) -> np.ndarray: ...
def run_blocks(fn: Callable[[Block], T], blocks: Sequence[Block], *, workers: int = 1) -> list[T]: ...
    # results in block order; with workers > 1, fn and T must be small and picklable: fn gets inputs as
    # FieldDescriptor and writes its inner block into "r+" outputs itself; T is a summary, never a block array
def worker_environment() -> dict[str, str]: ...          # single-threaded BLAS (§3.6)

# emulsion/classes.py
@dataclass(frozen=True, eq=False)
class SizeClasses:
    morphology: Literal["conventional", "tabular", "cubic"]   # cubic is treated as conventional
    ecd_um: np.ndarray                # (K,)
    thickness_um: np.ndarray          # (K,)
    volume_um3: np.ndarray            # (K,)
    surface_um2: np.ndarray           # (K,)
    projected_um2: np.ndarray         # (K,)
    number_fraction: np.ndarray       # (K,) w_c
    number_per_um2: np.ndarray        # (K,) nu_c
    silver_g_per_grain: np.ndarray    # (K,) m_Ag,c
    iodide_fraction: float
    dry_thickness_um: float           # L_e
    silver_coverage_g_m2: float
    @property
    def count(self) -> int: ...

def size_classes(layer: EmulsionLayer, n_classes: int | None = None) -> SizeClasses: ...

# optics/stack.py
@dataclass(frozen=True)
class BaseOptics:
    thickness_um: float
    refractive_index: float
    grey_density: float
    backing_density: float

@dataclass(frozen=True, eq=False)
class OpticalStack:
    grid: SpectralGrid
    emulsion_thickness_um: float
    alpha_per_um: np.ndarray          # (B,) alpha_b
    dye_cross_section: np.ndarray     # (B,) s_b
    mu_a_per_um: np.ndarray           # (B,)
    km_k_per_um: np.ndarray           # (B,) 2 mu_a
    km_s_per_um: np.ndarray           # (B,)
    t_top: float
    r_internal_top: float
    r_bottom: float
    base: BaseOptics
    psf_weight: float
    psf_length_um: float

def optical_stack(material: Material, film_format: str, classes: SizeClasses,
                  grid: SpectralGrid) -> OpticalStack: ...

# optics/kubelka.py
@dataclass(frozen=True, eq=False)
class DepthResponse:
    plus: np.ndarray                  # (n_z, B) g+
    minus: np.ndarray                 # (n_z, B) g-
    into_base: np.ndarray             # (B,) tau_b
    escaped_top: np.ndarray           # (B,)
    absorbed: np.ndarray              # (B,)

def depth_response(stack: OpticalStack, n_z: int) -> DepthResponse: ...

# optics/psf.py
def scatter_otf(nu_per_um: np.ndarray, weight: float, length_um: float) -> np.ndarray: ...
def scatter_kernel(pitch_um: float, weight: float, length_um: float, *, tail: float = 1e-4,
                   supersample: int = 16) -> np.ndarray: ...       # odd square, sums to 1

# optics/halation.py
@dataclass(frozen=True, eq=False)
class HalationKernel:
    kernel: np.ndarray                # odd square, returned fraction per pixel for unit flux
    fraction: float                   # f_H, untruncated
    truncated: float                  # fraction - kernel.sum()
    radius_um: float
    critical_radius_um: float

def halation_fraction(base: BaseOptics) -> float: ...
def halation_kernel(pitch_um: float, base: BaseOptics, *, energy_tol: float = 1e-3,
                    max_radius_um: float = 3000.0) -> HalationKernel: ...
def halation_otf(nu_per_um: np.ndarray, base: BaseOptics) -> np.ndarray: ...

# optics/expose.py
def exposure_fields(light: LightField, time_s: float, stack: OpticalStack, response: DepthResponse,
                    geometry: FrameGeometry, *, block_px: int = 512, workers: int = 1,
                    store: FieldStore | None = None) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]: ...
def uniform_exposure(photons_per_band: np.ndarray, stack: OpticalStack,
                     response: DepthResponse) -> tuple[np.ndarray, np.ndarray]: ...    # (n_z,), (n_z,)

# models/bw.py
MODEL = "bw-mechanistic"
@register_exposure(MODEL)
def expose_bw(material: Material, light: LightField, time_s: float, ctx: RunContext) -> EmulsionState: ...

# state.py
@dataclass(frozen=True)
class FieldSpec:
    key: str
    trailing: tuple[int, ...] = ()
    dtype: str = "float64"

class FrameAssembler:
    def __init__(self, geometry: FrameGeometry, specs: Sequence[FieldSpec],
                 store: FieldStore | None = None, *, name_prefix: str = "") -> None: ...
        # store names are f"{name_prefix}{spec.key}"; steps pass name_prefix=f"s{step_index}:" (§4.3)
    def write(self, key: str, inner: tuple[slice, slice], values: np.ndarray) -> None: ...
    def write_state(self, tile_state: EmulsionState, tile: Tile) -> None: ...
    def descriptors(self) -> dict[str, FieldDescriptor]: ...   # outputs for workers to open "r+" (mapped fields)
    def mark(self, inner: tuple[slice, slice]) -> None: ...      # register a region a worker wrote directly
    def reset(self) -> None: ...                                 # forget all writes, keep allocations (re-pass, §3.7)
    def finish(self) -> dict[str, np.ndarray]: ...               # each pixel written (or marked) exactly once since reset

class EmulsionState:                # D0 类，新增
    def drop(self, *keys: str) -> "EmulsionState": ...

# pipeline.py
@dataclass(frozen=True)
class RunContext:                   # D0 字段不变，新增：
    film_format: str | None = None
    profile: Literal["image", "calibration", "reference"] = "image"
    workers: int = 1
    memory_budget_gb: float = 6.0
    scratch_dir: Path | None = None
    store: FieldStore | None = field(default=None, compare=False, hash=False, repr=False)
        # the run's one FieldStore (§4.3); created and closed by the runner; not in the manifest or recipe_hash

def segments(proc: Process, model: str) -> list[tuple[int, ...]]: ...      # step indices per fused segment
```

测试（文件：`tests/optics/`、`tests/emulsion/`、`tests/test_light.py`、`tests/test_fieldstore.py`、`tests/test_state.py`、`tests/test_pipeline.py`）：§2.1.10 的全部，另加：

| 测试 | 断言 |
|------|------|
| `test_low_rank_light_window_and_project` | `window` 与 `project` 同展开后的 `Light` 计算之差 ≤ 1e-12（相对） |
| `test_low_rank_light_memory` | 6000×4000 的 `LowRankLight` 的 nbytes ≤ 0.6 GB（不生成 81 个波段） |
| `test_scene_light_low_rank_flag` | `low_rank=False` 时与 D0 逐位相同 |
| `test_fieldstore_memmap_readonly_shared` | 超过阈值的字段是 `mode="r"` 的 memmap；`readonly_array` 不复制；写入抛错 |
| `test_assembler_each_pixel_once` | 漏写或重写都抛 `ValueError`（`write` 与 `mark` 合计）；`reset()` 之后可以再写一遍 |
| `test_fieldstore_unique_backing_files` | 两个步骤以 `s2:emulsion/silver`、`s3:emulsion/silver` 分配时文件不同，后者不截断前者的只读映射；name 重复抛 `ValueError` |
| `test_run_blocks_workers_write_directly` | `workers = 3` 时工作进程经 `FieldDescriptor` 以 `"r+"` 写输出；主进程收到的结果只有摘要（每项序列化后 < 64 KB）；与 `workers = 1` 逐位相同 |
| `test_assembler_linear_cost` | 64 块与 4 块的拼接总时间之比 ≤ 1.5 × 面积比（不再平方增长） |
| `test_segment_fusion_equals_unfused` | 融合段与逐步执行的结果逐位相同 |
| `test_drop_fields` | `drop` 后字段消失，其余不变 |
| `test_d0_ideal_regression` | 理想材料样例输出与 D0 逐位相同 |

完成标准：通用标准；新增测试夹具 `tests/fixtures/library/films/test/bw-mechanistic.yaml`（虚构数值，全部 `assumed`，`model: bw-mechanistic`，带 D1a 已定义的 §5.4 字段：层、颗粒、光学；`FrozenModel` 是 `extra="forbid"`，所以 `latent`（D1b）与 `development`（D1c/D1d）等字段由引入相应模型的模块加进夹具，夹具随模块增长）；24 MP 的 `LowRankLight` 加 `bw-mechanistic` 曝光（用该夹具）在 8 GB 以内完成，峰值内存由 `tracemalloc` 加映射文件大小记录在测试日志中。D0 的夹具 `test/bw-neg`（`model: test-bw`）保持不变。

### D1b 潜影

文件：新建 `latent/__init__.py`、`latent/chain.py`、`latent/table.py`；修改 `data/models.py`（`LatentParams`）。

接口：

```python
# latent/chain.py
@dataclass(frozen=True)
class ChainParams:
    efficiency: float
    efficiency_size_exponent: float
    tau0_median_s: float              # math.inf allowed (no decay)
    tau0_sigma_ln: float
    ionic_time_s: float
    competing_sites: float            # >= 2
    hole_attack: float = 0.0
    cluster_cap: int = 5
    @classmethod
    def from_layer(cls, layer: EmulsionLayer) -> "ChainParams": ...

def generator(nu_per_s: float, tau0_s: float, params: ChainParams) -> np.ndarray: ...     # (S, S), columns sum to 0
def cluster_distribution(mean_photons: float, exposure_time_s: float, params: ChainParams, *,
                         ecd_um: float = 1.0, quadrature: int = 9) -> np.ndarray: ...  # (cap+1,), entry 1 == 0

# latent/table.py
@dataclass(frozen=True, eq=False)
class LatentTable:
    log10_m: np.ndarray               # (N_m,)
    probabilities: np.ndarray         # (K, N_m, cap+1), read-only
    exposure_time_s: float
    def lookup(self, mean_photons: np.ndarray, size_class: int) -> np.ndarray: ...   # (..., cap+1)

def latent_table(classes: SizeClasses, params: ChainParams, exposure_time_s: float, *,
                 nodes: int = 201, log10_range: tuple[float, float] = (-6.0, 4.0),
                 quadrature: int = 9) -> LatentTable: ...
```

测试（`tests/latent/`）：§2.2.9 的全部。

完成标准：通用标准；一张 K = 6、201 点、9 节点的表在 1 s 内建成（单核）。

### D1c 化学与显影

文件：新建 `chemistry/__init__.py`、`chemistry/species.py`、`chemistry/activity.py`、`chemistry/speciation.py`、`chemistry/solution.py`、`chemistry/kinetics.py`、`develop/__init__.py`、`develop/grid.py`、`develop/schedule.py`、`develop/ensemble.py`、`develop/rosenbrock.py`、`develop/lateral.py`、`develop/tank.py`、`develop/prepass.py`、`develop/step.py`；修改 `data/models.py`（§5.2、§5.3、§5.5、§5.6）、`data/library.py`（`kinetics()`、`validate` 的新检查）；数据：`data/chemistry/species.yaml`（D1 全部配料与组分，§2.3.7 的数值）、`data/chemistry/kinetics.yaml`（先验值，全部 `assumed` 或 `literature`，带 `bounds`）。

接口：

```python
# chemistry/species.py
@dataclass(frozen=True)
class AcidBaseFamily:
    component: str
    charges: tuple[int, ...]              # per state, most protonated first
    protons_rel: tuple[int, ...]          # per state, relative to the reference state
    pka25: tuple[float, ...]
    temperature_law: Literal["van_t_hoff", "manov_boric", "linear"]
    enthalpy_kj_mol: tuple[float, ...]
    dpka_dt: tuple[float, ...]
    activity_class: str
    def pka(self, temperature_k: float) -> np.ndarray: ...

@dataclass(frozen=True)
class ComponentInfo:
    id: str
    family: AcidBaseFamily | None
    d_water_um2_s: float                  # 20 C
    gel_ratio: float
    immobile: bool
    agent: AgentRole | None

def component_table(species: Mapping[str, Species]) -> dict[str, ComponentInfo]: ...

# chemistry/activity.py
def debye_huckel_d(ionic_strength: float, temperature_k: float, params: ActivityParams) -> float: ...
def conditional_pka(family: AcidBaseFamily, temperature_k: float, ionic_strength: float,
                    cation_molality: float, params: ActivityParams) -> np.ndarray: ...
def water_constants(temperature_k: float, ionic_strength: float,
                    params: ActivityParams) -> tuple[float, float, float]: ...   # log10 Kw, log10 gamma_H, log10 gamma_OH

# chemistry/speciation.py
@dataclass(frozen=True, eq=False)
class Speciation:
    families: tuple[AcidBaseFamily, ...]
    pka: np.ndarray                       # (F, J_max), conditional, NaN padded
    log10_kw: float
    log10_gamma_h: float
    log10_gamma_oh: float

def make_speciation(families: Sequence[AcidBaseFamily], temperature_k: float, ionic_strength: float,
                    cation_molality: float, params: ActivityParams) -> Speciation: ...
def solve_ph(spec: Speciation, totals: np.ndarray, toth: np.ndarray, *, ph0: np.ndarray | None = None,
             tol: float = 1e-12, max_iter: int = 30) -> np.ndarray: ...          # totals (F, N), toth (N,)
def fractions(spec: Speciation, ph: np.ndarray) -> np.ndarray: ...             # (F, J_max + 1, N)
def buffer_capacity(spec: Speciation, totals: np.ndarray, ph: np.ndarray) -> np.ndarray: ...

# chemistry/solution.py
@dataclass(frozen=True, eq=False)
class WorkingSolution:
    bath_id: str
    dilution: str
    dilution_fraction: float
    temperature_c: float
    component_ids: tuple[str, ...]
    totals: np.ndarray                    # (K_c,) mol/L
    toth: float
    cation_molality: float
    ionic_strength: float
    speciation: Speciation
    ph: float

def working_solution(bath: Bath, species: Mapping[str, Species], kinetics: KineticsFile, *,
                     dilution: str, temperature_c: float,
                     overrides: Mapping[str, float] | None = None) -> WorkingSolution: ...   # overrides: ingredient id -> g/L

# chemistry/kinetics.py
@dataclass(frozen=True, eq=False)
class Kinetics:                           # numeric constants resolved for one solution and temperature
    agent_ids: tuple[str, ...]
    roles: tuple[str, ...]
    charges: np.ndarray
    k_et: np.ndarray
    K_ads: np.ndarray
    lambda0: np.ndarray
    k_des: np.ndarray                     # NaN for regenerators
    k_reg: np.ndarray                     # NaN for ETAs
    sulfonated: np.ndarray                # bool
    K_br: float                           # at the solution temperature (Delta H_Br, §2.4.5)
    K_i: float                            # same van 't Hoff factor
    u_gel: float
    u_halide: float
    u_li: float
    n_sat: int
    solvent: SolventParams
    transport: TransportParams
    gelatin: GelatinParams

def resolve_kinetics(kinetics: KineticsFile, solution: WorkingSolution, species: Mapping[str, Species]) -> Kinetics: ...

@dataclass(frozen=True, eq=False)
class SurfaceState:
    theta: np.ndarray                     # (A, N)
    theta_halide: np.ndarray              # (N,)
    kappa_tilde: np.ndarray               # (N,) 1/s
    lambda_tilde: np.ndarray              # (N,) 1/s
    partition: np.ndarray                 # (A, N)
    regenerated: np.ndarray               # (A, N), NaN for regenerators

def surface_state(active: np.ndarray, bromide: np.ndarray, iodide: np.ndarray, kin: Kinetics,
                  ionic_strength: float) -> SurfaceState: ...                 # active (A, N) mol/L
def cluster_hazards(kin: Kinetics, fog_ratio: float, cluster_cap: int, free_sulfite_bulk: float) -> np.ndarray: ...

# develop/grid.py
@dataclass(frozen=True)
class ColumnGrid:
    n_z: int
    n_l: int
    thickness_um: np.ndarray              # (L,), bottom (emulsion z = n_z) to top
    kinds: tuple[Literal["emulsion", "overcoat", "liquid"], ...]
    height_um: np.ndarray                 # (L,) centres above the base
    delta_max_um: float
    delta_agit_um: float
    n_a: int                              # liquid cells below delta_agit (n_l // 2)
    @property
    def agit_cell(self) -> int: ...       # index of the cell whose top face is at delta_agit

def column_grid(emulsion_dry_um: float, overcoat_dry_um: float, swelling_ratio: float,
                delta_agit_um: float, delta_max_um: float,
                profile: Literal["image", "calibration", "reference"]) -> ColumnGrid: ...

# develop/schedule.py
@dataclass(frozen=True)
class AgitationPeriod:
    start_s: float
    stop_s: float

@dataclass(frozen=True, eq=False)
class Schedule:
    times_s: np.ndarray                   # (S+1,)
    agitated: np.ndarray                  # (S,) bool
    record_steps: Mapping[float, int]     # snapshot time -> step index

def agitation_periods(agitation: Agitation | None, duration_s: float) -> tuple[AgitationPeriod, ...]: ...
def step_schedule(periods: Sequence[AgitationPeriod], duration_s: float, *, h0_s: float, ratio: float,
                  h_max_s: float, record_s: Sequence[float] = ()) -> Schedule: ...

# develop/ensemble.py
@dataclass(frozen=True, eq=False)
class FilmColumnModel:
    classes: SizeClasses
    latent: LatentTable                   # cubic-spline lookup (§2.2.5); no F table (§2.4.2)
    hazards: np.ndarray                   # (cap+1,); entry 1 unused (n = 1 merged into 0)
    silver_per_cell: np.ndarray           # (K,) C_c,z in mol/L (equal for all z)
    gelatin_total: np.ndarray             # (L,) mol/L, 0 in the liquid
    fog_ratio: float
    grid: ColumnGrid

def film_column_model(material: Material, film_format: str, exposure_time_s: float, kin: Kinetics,
                      solution: WorkingSolution, profile: str) -> FilmColumnModel: ...

@dataclass(frozen=True, eq=False)
class Members:
    exposure_volume: np.ndarray           # "columns": (n_cells, n_z); "raster": (h, w, n_z) padded pixel window, um^-3
    exposure_surface: np.ndarray          # same shape, um^-2
    layout: Literal["columns", "raster"]  # columns: member i is cell i (strips, MTF pairs); raster: pixels over the cell grid
    pixel_um: float | None = None         # raster only: pixel pitch p
    origin_um: tuple[float, float] = (0.0, 0.0)   # raster only: corner of pixel [0, 0] relative to corner of cell [0, 0]
    # raster: separable 1-D area overlaps (aggregation into p-bar, §2.4.4) and bilinear weights at pixel
    # centres (output synthesis, §3.8) are computed per row batch; no (M, 4) index arrays are stored

@dataclass(frozen=True)
class Lateral:
    shape_cells: tuple[int, int]
    spacing_um: float
    boundary: Literal["neumann", "periodic", "pair"]
    wavenumber_per_um: np.ndarray | None = None   # "pair" only: (n_pairs,) k = 2 pi nu; columns ordered a0, b0, a1, b1, ...

@dataclass(frozen=True)
class TankSpec:
    volume_l: float
    film_area_um2: float

@dataclass(frozen=True, eq=False)
class BulkTrajectory:
    times_s: np.ndarray                   # (S+1,)
    totals: np.ndarray                    # (S+1, K_c)
    toth: np.ndarray                      # (S+1,)

@dataclass(frozen=True, eq=False)
class EnsembleResult:
    silver: np.ndarray                    # (M, K) g/m2, depth-merged
    initiated: np.ndarray                 # (M, K)
    silver_m2: np.ndarray                 # (M, K)
    silver_physical: np.ndarray           # (M,)
    halide: np.ndarray                    # (M,)
    growth_rate: np.ndarray               # (M,)
    snapshots: Mapping[float, "EnsembleResult"]
    flux_mean: np.ndarray                 # (S, K_c + 1) area-mean flux into the film, incl. TOTH (mobile-species
                                          # flux, §2.4.8) and the agitation-start mixing impulses (§2.4.9)
    bulk: BulkTrajectory
    min_scaled_concentration: float       # min over steps, cells, components of C / s_k (§3.3)
    cell_history: "CellHistories | None"  # only when record_cells > 0

@dataclass(frozen=True, eq=False)
class CellHistories:                      # monotone cumulative cell quantities on R recorded steps (§2.7.2)
    times_s: np.ndarray                   # (R,) a subset of schedule.times_s
    lambda_clock: np.ndarray              # (R, n_cells, n_z) Lambda_z(t)
    kappa_integral: np.ndarray            # (R, n_cells, n_z) K_z(t) = integral of kappa-tilde
    cell_spacing_um: float

def run_ensemble(film: FilmColumnModel, solution: WorkingSolution, kin: Kinetics, schedule: Schedule,
                 members: Members, *, n_cells: int, cell_weights: np.ndarray,
                 lateral: Lateral | None = None, tank: TankSpec | None = None,
                 bulk: BulkTrajectory | None = None, presoaked: bool = False,
                 method: Literal["fixed", "adaptive"] = "fixed", rtol: float = 1e-4,
                 record_cells: int = 0) -> EnsembleResult: ...
    # 恰好给 tank 与 bulk 之一：tank 表示与槽精确耦合（条带、预计算），bulk 表示给定槽液轨迹（图像块）
    # 成员只在开始（聚合 p-bar）与结束、各快照（合成 X、Y、P、F）时用到；时间循环里只有单元与核（§2.4.4）
    # 输出数组的第一维是成员：columns 为 (n_cells, ...)，raster 为内块像素 (h_in * w_in, ...)，按行优先

# develop/prepass.py
@dataclass(frozen=True, eq=False)
class ExposureHistogram:
    edges: np.ndarray                     # (bins + 1,) log10 m-bar
    counts: np.ndarray                    # (bins,) int64
    mean_volume: np.ndarray               # (bins, n_z)
    mean_surface: np.ndarray              # (bins, n_z)
    zero_count: int

def exposure_histogram(exposure_volume: np.ndarray, exposure_surface: np.ndarray, classes: SizeClasses, *,
                       bins: int = 64, block_px: int = 512) -> ExposureHistogram: ...
def tank_prepass(hist: ExposureHistogram, film: FilmColumnModel, solution: WorkingSolution, kin: Kinetics,
                 schedule: Schedule, tank: TankSpec, image_area_fraction: float,
                 presoaked: bool) -> tuple[BulkTrajectory, np.ndarray]: ...       # trajectory, J-bar_zero

@dataclass(frozen=True)
class PrepassCheck:
    delta: float                          # §3.7 acceptance quantity (developed-silver change of the 65 columns)
    depletion_change: float               # max |c' - c| / s_k, diagnostic
    ph_change: float                      # max |pH'_bulk - pH_bulk|, diagnostic
    trajectory: BulkTrajectory            # c'_bulk, the next trajectory if the check fails

def prepass_check(hist: ExposureHistogram, film: FilmColumnModel, solution: WorkingSolution, kin: Kinetics,
                  schedule: Schedule, tank: TankSpec, image_area_fraction: float, presoaked: bool,
                  current: BulkTrajectory, block_flux: np.ndarray, zero_flux: np.ndarray) -> PrepassCheck: ...
    # block_flux (S, K_c + 1): area-mean over all inner cells, summed from the block summaries in block order

# develop/step.py
@register_step(MODEL, StepKind.DEVELOP, interaction_range_um=None)
def develop_bw(state: EmulsionState, step: ProcessStep, ctx: StepContext) -> EmulsionState: ...
@register_step(MODEL, StepKind.PRESOAK, interaction_range_um=0.0)
def presoak_bw(state: EmulsionState, step: ProcessStep, ctx: StepContext) -> EmulsionState: ...
def develop_strip(material: Material, film_format: str, exposures: tuple[np.ndarray, np.ndarray],
                  exposure_time_s: float, step: ProcessStep, ctx: RunContext, *,
                  record_s: Sequence[float] = ()) -> EnsembleResult: ...      # exposures (M, n_z) each
```

测试（`tests/chemistry/`、`tests/develop/`）：§2.3.8、§2.4.13 的全部，§3.9 中 `test_convergence_*` 与 `test_calibration_profile_vs_reference`、`test_splitting_second_order`，§3.8 中 `test_halo_convergence`、`test_tank_prepass_check`、`test_image_profile_uniform_patch`（这些都对显影银 M 断言，§3.8；对密度、CI、log H_m 的断言在 D1d、D1e），§3.6 的 `test_block_result_independent_of_alignment`，§3.5 的 `test_lateral_toth_single_d`（可选，slow）。耗时超过 30 s 的标 `@pytest.mark.slow`，默认运行时跳过（在 `pyproject.toml` 注册标记）。

完成标准：通用标准；D-76 的 35 级标定档条带（小罐 7.5 min，314 步）单核 ≤ 2 s，实测每步耗时写入 §3.10；核分解与直接积分一致（`test_kernel_reduction_exact`）；图像档与参考档的均匀色块差 ≤ 0.014 g/m²（0.01 D 的银量当量）。

### D1d 成像、颗粒统计、停显定影与出口

文件：新建 `imaging/__init__.py`、`imaging/silver.py`、`grain/__init__.py`、`grain/poisson.py`、`grain/overlap.py`、`grain/statistical.py`、`develop/post.py`、`readout/__init__.py`、`readout/products.py`、`readout/densitometer.py`、`readout/scan.py`；修改 `runner.py`、`manifest.py`、`cli.py`、`data/models.py`（`EmulsionDevelopment`、`BathConstants`）；数据：`data/baths/stop-acetic.yaml`、`stop-water.yaml`、`fixer-thiosulfate.yaml`、`water.yaml`（常数 `assumed`，§2.6）。

接口：

```python
# imaging/silver.py
def covering_powers(classes: SizeClasses, dev: EmulsionDevelopment) -> np.ndarray: ...       # (K,) 1/(g/m2)
def developed_areas_um2(classes: SizeClasses, covering_power: np.ndarray) -> np.ndarray: ...  # (K,) a_c
def silver_density(state: EmulsionState, rows: slice | None = None, cols: slice | None = None) -> np.ndarray: ...

@dataclass(frozen=True, eq=False)
class TransmittanceFactors:
    base_density: np.ndarray              # (B,)
    dye_density: np.ndarray               # (B,)
    halide_density: np.ndarray            # (h, w)
    dye_residual: np.ndarray              # (h, w)

def transmittance_factors(state: EmulsionState, rows: slice, cols: slice) -> TransmittanceFactors: ...
def product_density(silver_transmittance: np.ndarray, factors: TransmittanceFactors,
                    product: np.ndarray) -> np.ndarray: ...        # (h, w); never builds (h, w, B) for neutral silver

# grain/poisson.py
def poisson_inverse(u: np.ndarray, mean: np.ndarray, *, exact_below: float = 30.0) -> np.ndarray: ...

# grain/overlap.py
def disk_intersection_area(tau_um: np.ndarray, radius_um: np.ndarray) -> np.ndarray: ...
def boolean_excess_variance(intensity_per_um2: np.ndarray, radius_um: np.ndarray, mean_t: np.ndarray,
                            pixel_area_um2: float, *, nodes: int = 48) -> np.ndarray: ...
def kappa(net_density: float, *, nodes: int = 48) -> float: ...            # constant radius, OPT 5.4

# grain/statistical.py
def counted_classes(developed_areas_um2: np.ndarray, pixel_um: float) -> np.ndarray: ...  # (K,) bool, p >= 4 sqrt(a_c/pi)
def epsilon_variance(var_bool: np.ndarray, mean_t: np.ndarray, mu: np.ndarray, q: np.ndarray,
                     counted: np.ndarray) -> np.ndarray: ...
    # max(Var_bool - Var_count, 0) / (T^2 + Var_count), Var_count = T^2 (exp(sum_{counted} mu q^2) - 1), q from x-bar
def statistical_transmittance(state: EmulsionState, rows: slice, cols: slice) -> np.ndarray: ...  # (h, w), keyed by state.seed
# boolean_excess_variance 不变：调用方传入 r̃_c = sqrt(a_c (1 + phi) <x^2>^0.5 / pi)，强度 nu_c F_c 含全部档

# develop/post.py
def stop_update(x: np.ndarray, y: np.ndarray, f: np.ndarray, kappa_end: np.ndarray, tau_s: float,
                time_s: float) -> tuple[np.ndarray, np.ndarray]: ...
@register_step(MODEL, StepKind.STOP, interaction_range_um=0.0)
def stop_bw(state: EmulsionState, step: ProcessStep, ctx: StepContext) -> EmulsionState: ...
@register_step(MODEL, StepKind.FIX, interaction_range_um=0.0)
def fix_bw(state: EmulsionState, step: ProcessStep, ctx: StepContext) -> EmulsionState: ...
@register_step(MODEL, StepKind.WASH, interaction_range_um=0.0)
def wash_bw(state: EmulsionState, step: ProcessStep, ctx: StepContext) -> EmulsionState: ...
@register_step(MODEL, StepKind.DRY, interaction_range_um=0.0)
def dry_bw(state: EmulsionState, step: ProcessStep, ctx: StepContext) -> EmulsionState: ...

# readout/products.py
def visual_product(grid: SpectralGrid) -> np.ndarray: ...       # S_A * ybar
def type2_product(grid: SpectralGrid) -> np.ndarray: ...        # anchor interpolation, §2.8
def scan_product(grid: SpectralGrid) -> np.ndarray: ...         # S_D50 * ybar

# readout/densitometer.py
@register_readout("densitometer_visual")
def read_visual(state: EmulsionState, ctx: RunContext, *, aperture_um: float | None = None,
                grain: bool = False) -> Readout: ...            # kind "density", (H, W)
@register_readout("densitometer_type2")
def read_type2(state: EmulsionState, ctx: RunContext, *, aperture_um: float | None = None,
               grain: bool = False) -> Readout: ...
@register_readout("silver_density")
def read_silver_density(state: EmulsionState, ctx: RunContext) -> Readout: ...

# readout/scan.py
SCAN_G_REF: float = 0.56
SCAN_Y_REF: float = 0.18
SCAN_D_REF: float      # = 0.10 + 0.56 * log10(12.5)
@register_readout("scan_standard")
def read_scan(state: EmulsionState, ctx: RunContext, *, grain: bool = True) -> Readout: ...   # image_linear_srgb

# runner.py
def run_image(input_path: Path, output_dir: Path, *, material_id: str, process_id: str,
              readout_device: str | None = None, illuminant: str = "D65", iso: float | None = None,
              ev: float = 0.0, time_s: float = 1 / 125, seed: int = 0, tile_px: int | None = None,
              long_edge_mm: float = 36.0, library: Library | None = None,
              film_format: str | None = None, profile: Literal["image", "calibration", "reference"] = "image",
              workers: int | None = None, memory_budget_gb: float = 6.0, scratch_dir: Path | None = None,
              keep_scratch: bool = False) -> tuple[Path, RunManifest]: ...

# manifest.py：RunManifest 新增
#   film_format: str | None = None
#   model_settings: JsonDict = {}
#   kinetics_file: FileRef | None = None
#   calibration_runs: list[str] = []
#   workers: int | None = None              # 不计入 recipe_hash
```

测试：§2.5.4、§2.6、§2.7.4、§2.8 的全部；§3.8 中 `test_image_profile_edge`、`test_image_profile_mtf`、`test_resolution_independence`；另加：

| 测试 | 断言 |
|------|------|
| `test_run_image_bw_end_to_end` | 夹具 `test/bw-mechanistic` + 新夹具显影液 `test-bw-mech-dev`（虚构配方，组分取 D1c 的真实物种：`Library` 只读一个根，所以把所需的真实物种条目复制进 `tests/fixtures/library/chemistry/species.yaml`，并加一份先验值的 `tests/fixtures/library/chemistry/kinetics.yaml`）与流程 `test-bw-mechanistic`，256×384 样图跑通，写出图像与运行记录；同种子重跑逐位相同；`tile_px=64` 逐位相同；`workers` 1 与 3 逐位相同 |
| `test_manifest_records_model_settings` | 运行记录含档位、离散化、`kinetics_file` 哈希、`calibration_runs`；`workers` 不影响 `recipe_hash` |
| `test_scratch_cleanup` | 运行结束后临时目录被删除；`keep_scratch=True` 时保留 |
| `test_benchmark_image_path`（slow） | D-76 小罐 7.5 min，默认 6 GB 预算与默认进程数：1 个规范块单进程、`workers` 个块并行，实测并行效率并按 §3.10 外推，24 MP 外推值 ≤ 30 min；峰值常驻内存 ≤ 6 GB；单元·步耗时与效率写入日志 |

完成标准：通用标准；`atlas run tests/assets/sample.png -m test/bw-mechanistic -p test-bw-mechanistic --data tests/fixtures/library -o out/` 输出一张有颗粒、有邻接效应的黑白负片扫描图（夹具数据），运行记录完整。

### D1e 虚拟感光仪、指标、标定与验证报告（含数据文件）

文件：新建 `validate/__init__.py`、`validate/conditions.py`、`validate/sensitometer.py`、`validate/metrics.py`、`validate/reference.py`、`validate/splits.py`、`validate/report.py`、`calibrate/__init__.py`、`calibrate/parameters.py`、`calibrate/objective.py`、`calibrate/fit.py`、`calibrate/identifiability.py`、`calibrate/writeback.py`；修改 `data/models.py`（`FitRecord`、`Quantity.fit` 与校验、`CalibrationRecord`）、`data/library.py`（`reference_tables`、`developer_facts`、`calibration_record`）、`cli.py`（`validate`、`calibrate`）、`pyproject.toml`（`report` 依赖组：matplotlib）；数据：`data/films/kodak/tri-x-400.yaml`、`data/films/kodak/t-max-100.yaml`、`data/developers/{d-76,d-76-published,xtol,hc-110-pre2019,hc-110-2019,t-max-dev,rodinal}.yaml`、`data/processes/kodak/...`（Kodak 曲线对应的流程）、`data/calibration/<run_id>.yaml`。参考数据由数字化团队并入 `data/reference/`（§5.9），本模块只读。

接口（`Condition` 定义在 `validate/conditions.py`，`StripResult` 定义在 `validate/sensitometer.py`，字段见 §6.1）：

```python
# validate/conditions.py
def parse_conditions_header(text: str) -> dict[str, str]: ...
def conditions_for(material_id: str, library: Library, *, split: Literal["train", "within", "held_out"]) -> list[Condition]: ...

# validate/sensitometer.py
def sensitometric_strip(condition: Condition, ctx: RunContext, *, log_h: np.ndarray | None = None) -> StripResult: ...
def sensitometric_batch(conditions: Sequence[Condition], ctx: RunContext) -> list[StripResult]: ...
def exposure_for_net_density(condition: Condition, net_density: float, ctx: RunContext) -> float: ...   # lux*s

# validate/metrics.py
@dataclass(frozen=True)
class SpeedResult:
    log_h_m: float
    delta_d: float
    speed: float
    speed_degrees: float
    iso_bin: str | None
    reason: str = ""

def iso_speed(strip: StripResult, *, density: Literal["type2", "visual"] = "type2",
              base_fog: float | np.ndarray | None = None) -> SpeedResult: ...   # base_fog 标量或 (T,)：None 用零曝光片；与 Kodak 曲线比较时传 D_min*（§6.3）
def d_min_star(log_h_ref: np.ndarray, density_ref: np.ndarray, log_h_sim: np.ndarray,
               density_sim: np.ndarray) -> tuple[float, float]: ...   # (参考侧, 模拟侧)，都在 x_left = log_h_ref[0] 处取值
def average_gradient(log_h: np.ndarray, density: np.ndarray, base_fog: float) -> float: ...
def contrast_index(log_h: np.ndarray, density: np.ndarray, d_min: float) -> float: ...
def chord_contrast_index(log_h: np.ndarray, density: np.ndarray, d_min: float) -> float: ...
def time_to_ci(times_s: np.ndarray, ci: np.ndarray, target: float) -> float: ...   # 只用于报告与"显影到 CI 0.56"（§6.10）
def ci_at_times(times_s: np.ndarray, ci: np.ndarray, at_s: np.ndarray) -> np.ndarray: ...   # 单调插值；表格格评分用（§7.2）

@dataclass(frozen=True)
class GranularityResult:
    rms: float                            # sampled (report only)
    standard_error: float
    expected_rms: float                   # G_exp, closed form of §6.6 (the calibration residual)
    net_density: float
    n_apertures: int

def rms_granularity(condition: Condition, ctx: RunContext, *, net_density: float = 1.0,
                    aperture_um: float = 48.0, patch_um: float = 3200.0, pitch_um: float = 4.0,
                    sample: bool = True) -> GranularityResult: ...      # sample=False: rms, standard_error = NaN
def mtf_effective_exposure(condition: Condition, frequencies_c_mm: Sequence[float], ctx: RunContext, *,
                           method: Literal["small_signal", "periodic"] = "small_signal",
                           modulation: float = 0.4, mean_net_density: float = 1.0) -> np.ndarray: ...   # percent
    # small_signal: §6.7 column pairs on ctx.profile with eps = 1e-3 (modulation ignored); periodic: 32-column strip
@dataclass(frozen=True, eq=False)
class ReciprocityResult:
    times_s: np.ndarray
    stops: np.ndarray                     # (n_times, n_criteria)
    adjusted_time_ratio: np.ndarray       # log10(t_adj_sim / t_adj_kodak), NaN where not applicable
    development_change_percent: np.ndarray

def reciprocity(condition: Condition, times_s: Sequence[float], ctx: RunContext, *,
                criteria: Sequence[float] = (0.10, 0.80)) -> ReciprocityResult: ...
def spectral_sensitivity(condition: Condition, ctx: RunContext, *, exposure_time_s: float,
                         criteria: Sequence[float] = (0.3, 1.0),
                         bands_nm: Sequence[float] | None = None) -> np.ndarray: ...   # (n_bands, n_criteria) log10(cm2/erg)

# validate/splits.py
TRAINING: frozenset[str]
WITHIN_CHECKS: frozenset[str]
HELD_OUT: frozenset[str]
UNUSED: frozenset[str]                       # 四者互不相交，覆盖全部参考项与子项（§7.6）
HELD_OUT_EVALUATED: frozenset[str]           # ⊂ HELD_OUT
NOT_EVALUATED_REASON: Mapping[str, str]      # HELD_OUT − HELD_OUT_EVALUATED 与 UNUSED 中每项的原因
PAIR_METRICS: Mapping[tuple[str, str], tuple[str, ...]]   # (material, developer) → 适用指标 id（§8.3 组合–指标表）
class LeakageError(RuntimeError): ...
def assert_trainable(item_id: str) -> None: ...
def registry_sha256() -> str: ...            # 规范化 JSON 的哈希（§8.4）
class TrainingLibrary(Library):              # 过滤视图：对 HELD_OUT 的访问抛 LeakageError，表格返回时去掉留出子项
    def __init__(self, base: Library) -> None: ...

# validate/reference.py
@dataclass(frozen=True, eq=False)
class Observation:
    kind: Literal["density", "ci", "table_ci", "speed", "ph", "spectral", "mtf", "granularity", "reciprocity"]
    item_id: str
    condition: Condition | None
    x: np.ndarray                         # table_ci：t_cell（s）；density：参考曲线的 log H
    y: np.ndarray                         # table_ci：CI_aim
    sigma: np.ndarray                     # 已含 §7.2 的 √(N_pts/N_eff)、σ_floor、√n_class 放大
    weight: np.ndarray
    n_effective: float                    # 本项的独立数（§7.2），s² 用（§7.8）
    d_min_ref: float | None = None        # density：参考侧 D_min*（§6.3）
    in_objective: bool = True             # speed 与 MTF 峰段为 False（只报告）

def training_observations(library: TrainingLibrary) -> list[Observation]: ...      # raises LeakageError on held-out ids
def within_check_observations(library: Library) -> list[Observation]: ...
def held_out_observations(library: Library, record: CalibrationRecord) -> list[Observation]: ...   # requires record.frozen

# calibrate/parameters.py
@dataclass(frozen=True)
class ParameterSpec:
    name: str                             # "<file id>:<field path>"；只在标定中存在的参数（δ_CC、s_abs，§7.1）为 "calibration:<name>"，不写回
                                          # T-Max 的 "calibration:aspect_ratio"（§7.1）写回为 grain.mean_diameter = AR · grain.thickness
    unit: str
    transform: Literal["linear", "log10", "logit"]
    lower: float
    upper: float
    prior: Literal["uniform", "normal", "lognormal"]
    prior_mu: float
    prior_tau: float
    block: Literal["A", "B", "F", "C"]
    stage: str

class LibraryOverlay(Library):            # in-memory parameter overrides, no files written
    def __init__(self, base: Library, overrides: Mapping[str, float]) -> None: ...   # 标定中 base 必须是 TrainingLibrary，全部读取委托给它（§7.7）

def parameter_specs(library: Library, stage: str) -> list[ParameterSpec]: ...
def to_theta(specs: Sequence[ParameterSpec], values: Mapping[str, float]) -> np.ndarray: ...
def from_theta(specs: Sequence[ParameterSpec], theta: np.ndarray) -> dict[str, float]: ...

# calibrate/objective.py
def residuals(theta: np.ndarray, specs: Sequence[ParameterSpec], observations: Sequence[Observation],
              library: Library, ctx: RunContext, *, include_prior: bool = True) -> np.ndarray: ...
    # include_prior=True：数据残差 + 先验残差 + A_550 罚项（拟合用）；False：只有数据残差（§7.5 的 SVD、剖面与 σ_data）
    # 返回值恒有限：NaN/inf 项替换为有界罚值 0.5/σ 并计数（§7.2）

# calibrate/fit.py
def run_stage(stage: str, library: TrainingLibrary, ctx: RunContext, *, starts: int = 8, seed: int = 0) -> StageRecord: ...
def calibrate(library: TrainingLibrary, ctx: RunContext, *, run_id: str) -> CalibrationRecord: ...

# calibrate/identifiability.py
def weighted_svd(jacobian: np.ndarray, sigma: np.ndarray) -> tuple[np.ndarray, np.ndarray]: ...   # singular values, right vectors
@dataclass(frozen=True, eq=False)
class ProfileResult:
    grid: np.ndarray
    chi2: np.ndarray
    interval: tuple[float, float] | None
    identifiable: bool                    # 数据残差的 χ²（include_prior=False），按 s² 缩放；另一块固定于不动点

def profile_likelihood(spec: ParameterSpec, specs: Sequence[ParameterSpec], library: Library,
                       ctx: RunContext, *, points: int = 8) -> ProfileResult: ...    # specs 只含 spec 所在块（§7.5：8 个格点）
def prior_dominated(spec: ParameterSpec, sigma_data_transformed: float) -> bool: ...   # σ_data > 0.5 τ_prior（§7.5）

# calibrate/writeback.py
def write_back(record: CalibrationRecord, library_root: Path) -> list[Path]: ...

# validate/report.py
def validation_report(material_id: str, ctx: RunContext, out_dir: Path, *, held_out: bool = False,
                      after_unblinding: bool = False) -> Path: ...
```

测试（`tests/validate/`、`tests/calibrate/`）：§6 中的全部（`test_iso_speed_on_kodak_tri_x`、`test_ci_straight_line`、`test_ci_matches_report_on_tmax_d76`、`test_ci_nan_on_short_curve`、`test_dmin_star_same_operator`、`test_mtf_linear_system`、`test_mtf_periodic_exact`、`test_mtf_small_signal_limit`、`test_mtf_small_signal_vs_periodic`、`test_granularity_expected_vs_sampled`、`test_spectral_energy_conversion`），§3.9 的 `test_calibration_profile_vs_reference_metrics`，§7.4 的 `test_sparse_jacobian_matches_dense`，§7.7 的泄漏测试，另加：

| 测试 | 断言 |
|------|------|
| `test_conditions_header_parsing` | 全部暂存 CSV 的 conditions 行解析为 `Condition`，温度按 °F 精确换算 |
| `test_split_registry_complete` | 每个参考文件、表格行与子项（§7.6 的子项标识）恰好属于 `TRAINING`、`WITHIN_CHECKS`、`HELD_OUT`、`UNUSED` 之一；`HELD_OUT_EVALUATED ⊂ HELD_OUT`；每个不评估项与 `UNUSED` 项都有原因；J-78 的 `material: null` 行全在 `UNUSED`；`PAIR_METRICS` 覆盖 §8.3 组合–指标表的 14 个组合 |
| `test_training_library_filters_held_out` | 经 `TrainingLibrary` 读任一 `HELD_OUT` 项抛 `LeakageError`；读混合文件（如 `t-max-100/tables.yaml`）时返回值不含留出子项 |
| `test_residuals_always_finite` | 在阶段 C、D 的 8 个 Sobol 起点上，`residuals` 全部有限；人为让一个 CI 无解时该项为 0.5/σ 且 `n_penalised` 计数 |
| `test_held_out_pipeline_finite_on_synthetic` | 用已知 θ 生成、结构与留出各项相同的合成参考（不读 Kodak 留出文件），完整的留出评估流程对 `PAIR_METRICS` 的每个适用指标给出有限值，并给出三个标题数字与基线数字 |
| `test_nan_metric_fails_pair` | 适用指标算不出时组合状态为 `failed`，`failed_metrics` 带原因 |
| `test_held_out_refuses_mismatched_freeze` | 工作树不干净、提交或 `split_registry_sha256`/`metrics_code_sha256` 与记录不符时，`validation_report(..., held_out=True)` 拒绝运行 |
| `test_identifiability_ignores_prior` | 一个只受先验约束的参数（数据雅可比列为零）：数据残差的 SVD 给出零奇异值，判为不可辨识或先验主导，不写为 fitted |
| `test_parameter_transform_round_trip` | `from_theta(to_theta(v)) == v`，相对差 ≤ 1e-14 |
| `test_synthetic_recovery` | 用已知参数生成的合成数据（标定档、加 σ 噪声，种子 0），阶段 C 拟合恢复可辨识参数在 2σ 内 |
| `test_forward_model_smooth` | 参数扰动 1e-6 时目标函数的二阶差分 / 扰动² 在三个尺度上相差 < 10 % |
| `test_writeback_marks_fitted` | 写回后 `source: fitted` 且带 `fit` 记录；不可辨识参数为 `assumed` 并有说明 |
| `test_report_generates` | 用夹具数据生成报告目录，含 `report.md`、SVG 与 CSV |

完成标准：通用标准；数据文件全部通过 `atlas check-data`；阶段 A–F 完成并写回；`atlas validate kodak/tri-x-400` 与 `kodak/t-max-100` 生成训练报告（含朴素基线与 §7.6 的内部检查），交所有者确认阈值。D1e 结束时**不冻结、不做留出评估**：冻结与唯一一次留出评估（§8.4）在 D1f 的跨模式测试通过之后进行，作为 D1f 的完成标准。

### D1f 显微镜模式

文件：新建 `micro/__init__.py`、`micro/grains.py`、`micro/develop.py`、`micro/render.py`、`micro/measure.py`；修改 `cli.py`（`atlas micro`）。

接口：

```python
# micro/grains.py
@dataclass(frozen=True, eq=False)
class GrainSet:
    x_um: np.ndarray                      # (G,)
    y_um: np.ndarray
    depth_index: np.ndarray               # (G,) int
    size_class: np.ndarray                # (G,) int
    cluster_size: np.ndarray              # (G,) int, 0 or 2..cap
    init_time_s: np.ndarray               # (G,), inf if never initiated
    developed_fraction: np.ndarray        # (G,)

def place_grains(classes: SizeClasses, patch_um: tuple[float, float], n_z: int, seed: int,
                 origin_um: tuple[float, float] = (0.0, 0.0)) -> GrainSet: ...    # keyed by micro-cell coordinates

# micro/develop.py（CellHistories 定义在 develop/ensemble.py，由 run_ensemble(record_cells=64) 给出，§2.7.2）
def develop_patch(material: Material, film_format: str, exposure: tuple[np.ndarray, np.ndarray],
                  exposure_time_s: float, step: ProcessStep, ctx: RunContext, *,
                  cell_um: float = 4.0, record: int = 64) -> CellHistories: ...
    # uniform exposure: one column, histories broadcast to all cells (read-only views)
def develop_grains(grains: GrainSet, exposure: tuple[np.ndarray, np.ndarray], histories: CellHistories,
                   film: FilmColumnModel, seed: int) -> GrainSet: ...

# micro/render.py
def render_coverage(grains: GrainSet, areas_um2: np.ndarray, window_um: tuple[float, float, float, float],
                    pitch_um: float) -> np.ndarray: ...       # anti-aliased covered fraction per raster pixel

# micro/measure.py
def aperture_transmittances(grains: GrainSet, areas_um2: np.ndarray, centres_um: np.ndarray, *,
                            aperture_um: float = 48.0, pitch_um: float | None = None) -> np.ndarray: ...
def micro_granularity(condition: Condition, ctx: RunContext, *, net_density: float,
                      patch_um: float = 2000.0, seed: int = 0) -> GranularityResult: ...
```

测试（`tests/micro/`）：

| 测试 | 断言 |
|------|------|
| `test_render_single_disk_area` | 单个圆盘、p_r ≤ r/4 时覆盖面积相对误差 ≤ 1e-4 |
| `test_cell_history_thinning` | 一根 D-76 柱：64 个记录步上的单调 PCHIP 插值给出的颗粒显影分数，与逐步记录的结果之差 ≤ 1e-3 |
| `test_micro_mean_density_matches_nutting` | D = 0.3、1.0、2.0 的大块平均密度与 Nutting 之差 ≤ 0.005 |
| `test_micro_constant_radius_kappa` | 全部颗粒等大且完全显影时，孔径方差与 Siedentopf 之比 = κ(D) ± 5 % |
| `test_micro_vs_statistical_granularity`（slow） | T-Max 与 Tri-X 在 D = 0.3、1.0、2.0：σ_D 之比 ∈ [0.94, 1.06] |
| `test_micro_deterministic` | 同种子逐位相同；不同原点的同一显微格得到同一批颗粒 |

完成标准：通用标准；`atlas micro` 输出一张 0.5 mm 见方的显微渲染图（PNG）与它的颗粒度；跨模式测试通过即完成总设计 D1 的"统计模式与显微镜模式的颗粒度一致"；然后按 §8.4 冻结标定记录并做唯一一次留出评估，生成含留出部分的验证报告（D1 的最后一步）。

---

## 11 对既定决策的异议与澄清

1. **"每像素每档二项抽样"（决策 4(6)）。** 我按"对泊松分布的颗粒数做二项稀疏化"实现，闭式就是泊松抽样（§2.7.1）。显微镜模式与布尔模型的颗粒位置是泊松的，若统计模式用固定总数 N 的二项分布，方差会比参照小 (1 − p) 倍，大颗粒档在 D ≈ 1 时 p 接近 1，跨模式检验会按设计失败。若所有者要的确实是固定 N（例如为了表示颗粒的堆积排斥），显微镜模式也要相应改成硬核点过程，那是另一个模型。
2. **"21 级光楔"（决策 7）。** 21 级 × 0.15 只覆盖 3.0 log H，Kodak 的曲线横跨约 4 log H（−3.45 至 +0.5）。虚拟光楔按参考曲线的范围取级（约 28–34 级，CAL §1.6），每级仍是一根深度柱，成本不变。
3. **胶片参数只用 D-76（决策 2）与联合精修。** CAL §3.3 建议最后对全部训练数据做联合精修，那样 Tri-X 的胶片参数会看到 XTOL、HC-110、T-MAX 的数据，与决策 2 相冲突。我改为块坐标迭代：胶片块只看 D-76 与胶片自己的数据，化学块看全部训练数据，交替到收敛（§7.3 阶段 E）。代价是得到的是两个目标的不动点而不是一个联合最优，报告胶片参数在迭代中的偏移作为诊断。
4. **Tri-X×T-MAX Developer 的训练数据。** 决策 2 只写了它的特性曲线。我把它的 CI(t) 曲线与时间表也列入训练（§7.6），因为它们属于训练组合、能约束 T-MAX 的配方干扰参数。若所有者只想用曲线，从清单删去这两项即可，不影响其他部分。
5. **显影步骤的分块（决策 5、6b）。** 显影步骤注册为全幅耦合（`interaction_range_um=None`），内部用固定在画幅坐标上的规范块。结果对 `tile_px` 与进程数逐位不变，但与"无限晕边"的解不逐位相等，差别由晕边测试界定（§3.6）。用户的 `--tile` 对显影没有作用，文档里写明。
6. **ISO Type 2 近似（决策 4(7)）。** D1 的银、片基与残留染料都按光谱中性处理，模型里的 Type 2 密度与视觉密度恒等，所以 Type 2 的近似对 ISO 感光度检验没有影响；它只在将来出现非中性成分时起作用。ISO 5-3 的 Type 2 全表与 ISO 7589 的镜头透射全表在已缓存的预览里可以读到（CAL §0 第 2 条），是否可以写进数据文件是许可问题，留给所有者。
7. **留出曲线的阈值（决策 3）。** CAL §4.2 没有为留出的特性曲线给阈值（它的方案里没有留出曲线）。我按"留出阈值取训练阈值的 2 倍"给了暂定值（§8.3），与 CAL 对速度的 0.05 → 0.10 一致，请所有者确认。
8. **倒易律留出行在设计期暴露。** 为确认链的结构能达到 Kodak 表值，§12.2–§12.3 用了两张表的全部行做可行性计算（玩具群体，不是胶片拟合）。这些行仍不入拟合，但它们检验的是常数的外推，不再能检验结构选择（`S_eff` 闭合正是在看到 Tri-X 1/100 000 s 那一行之后选的）。若所有者需要一个对结构也是盲的倒易律检验，只能采用 LAT §5 的方案 (b)：只用 Tri-X 标定潜影常数、预测 T-Max 的整张表；但 T-Max 的表我也已经看过，同样不再是盲的。所以这些行归内部检查（§7.6），不计入任何状态；τ̄₀、σ_τ 的界与示例值改为只由 LAT 导出（§2.2.7）。
9. **Rodinal。** 对氨基酚的动力学常数没有任何数据（KIN §7 第 11 条），取米吐尔的值并标 `assumed`；醌亚胺水解释放碱的支路不建模。Rodinal 的结果只是"在这些假设下的预测"，报告中这样表述。
10. **阈值确认与留出评估的时点（决策 3）。** "所有者看到第一份报告后确认或修改"按训练报告理解：确认在训练报告之后、冻结之前；冻结之后改阈值只能并列显示为"事后阈值"，冻结阈值下的状态永久保留（§8.4）。唯一一次留出评估从 D1e 移到 D1f 完成之后，因为 D1f 复用并首次在细网格上运行共用的显影代码（§8.4 第 1 步）。若所有者希望 D1e 就给出留出结果，替代办法是在 D1e 冻结前先跑 D1f 的引擎细网格测试。
11. **MTF 峰段不拟合。** 峰段（> 100 %）只由横向化学产生，而胶片块里能动 MTF 的只有 w_P、a_P（§6.7 第 8 步），所以峰段归横向化学的内部检查。若所有者更想把 MTF 峰段放进阶段 D 拟合 δ、r_gel 等，D1 就不再有任何对邻接效应的检查。

---

## 12 设计期数值核对

脚本放在会话的临时目录，用 `uv run` 在仓库环境中运行，不入库。

### 12.1 泊松 r 次命中极限

两位点链（36 态），τ₀ = ∞、τ_i = 0、`S_eff = 2`、ρ₁ = 0：m ∈ {0.5, 3, 8, 20}、t ∈ {1e-4, 1, 100} s、r = 2…5 的 P(n ≥ r) 与 `1 − Σ_{j<r} e^{−ηm}(ηm)^j/j!` 之差全部 < 1e-10。解析上这是速率为 ην 的纯生过程，事件数在 t 时刻服从均值 ηνt = ηm 的泊松分布，与 t 无关，故没有倒易律失效。

### 12.2 低照度失效与陷阱深度分布

玩具群体（5–7 档对数正态，截面 ∝ d^{1.5}，判据为 P(n ≥ 3)），η = 0.5：

| τ̄₀，σ_τ | 判据 | 1 s | 10 s | 100 s | 1000 s（档，相对 1/100 s） |
|---------|------|-----|------|-------|------------------------------|
| 0.3 s，0 | 0.05 | +0.31 | +1.22 | +2.56 | +4.10 |
| 0.3 s，2 | 0.05 | +0.35 | +0.93 | +1.89 | +3.21 |
| 0.3 s，2 | 0.4 | +0.45 | +1.22 | +2.46 | +3.98 |
| 0.3 s，3 | 0.05 | +0.35 | +0.80 | +1.47 | +2.53 |

单一 τ₀ 在 10–1000 s 段约 1.3–1.5 档/十倍，比 Tri-X 表（约 1 档/十倍）陡；陷阱深度分布把它压平，且较高判据需要更多补偿（反差上升），与 Kodak 长曝光减显影一致。拟合到两张表的全部行（两位点、无丢失）：Tri-X 得 τ̄₀ = 0.074 s、σ_τ = 1.12、τ_i = 7.8e-6 s，各行 (1e-5 … 100 s) 为 0.46/0.15/−0.01/0/0.16/0.73/1.81/3.24 档，对表值 1/0.5/0/0/0/1/2/3；T-Max 100 得 τ̄₀ = 1.92 s、σ_τ = 3.47、τ_i = 1.8e-5 s，为 0.32/0.33/0.02/0/0.07/0.20/0.47/1.03，对表值 —/0.33/0/0/0/0.33/0.5/1。

### 12.3 高照度失效的闭合

同一玩具群体，τ̄₀ = 0.074 s、σ_τ = 1.12，判据 0.1（括号内为判据 0.4）：

| S_eff | τ_i | 1/100 000 s | 1/10 000 s | 1/1000 s |
|-------|-----|-------------|------------|----------|
| 2 | 3e-6 | +0.43（+0.45） | +0.04（+0.09） | −0.02 |
| 2 | 1e-5 | +0.43（+0.32） | +0.19（+0.29） | 0 |
| 3 | 1e-5 | +1.14（+1.13） | +0.21（+0.37） | 0 |
| 3 | 3e-5 | +1.10（+1.05） | +0.66（+0.90） | +0.04 |
| 5 | 1e-5 | +1.88（+1.99） | +0.21（+0.40） | 0 |

S_eff = 2 时最大约 0.45 档（强光下电子全流向位点 2，分散度回落）；S_eff = 3 能覆盖 Tri-X 的 (+1, +½)，并在 1/10 000 s 给出较高判据需要更多补偿（反差下降），与 Kodak 短曝光加显影一致。

### 12.4 pH 锚点

电荷平衡（25 °C），理想溶液与只含 SIT 的 Debye–Hückel 项（Δε = 0）：

| 溶液 | I（M） | 理想 | 含 D(I) | 实测 |
|------|--------|------|---------|------|
| D-76 公开配方 | 2.39 | 8.95 | 8.44 | 8.5（Carlton & Crabtree，比色） |
| 同上 + B₂O₃ 1 / 2 / 3 g/L | — | 8.82 / 8.73 / 8.67 | 8.32 / 8.23 / 8.17 | 成品 8.5 ± 0.05（CIS-61） |
| XTOL 例 1 工作液 | 2.14 | 8.89 | 8.29 | 8.20 ± 0.05 |
| XTOL 例 1 A 部分（850 mL） | 0.32 | 11.05 | 10.76 | 10.19 ± 0.05 |

理想值复现了 DEV §5.1 的 8.94、8.88、11.03。D(I) 在 I = 2.39 M 时为 0.237，亚硫酸氢根的混合 pKa 下降 3D = 0.71。高离子强度的两个锚点已在 0.1 以内，低离子强度的 A 部分仍高 0.57，不是活度能解释的（I 只有 0.3 M），更可能是配方读法（例如 DTPA 以酸而非五钠盐加入），所以它在拟合中低权重（§7.2）。成品 D-76 的张力见 §9 第 3 条。

### 12.5 D-76 的量级

见 §2.4.12 的表。I = 2.39 M 时 D(I) = 0.237；按 §2.3.3 的混合常数（pH 为 a_H，HSO₃⁻/SO₃²⁻ 位移 −3D），亚硫酸盐在 pH 8.5 的缓冲容量 17.1 mM/pH（理想 81.3）；加上硼酸盐 8.4（5.9）、HQ 6.9（4.3）、米吐尔 0.8（0.5），总 β 为 33.3 mM/pH（理想 92.0），§2.4.12 的表面 pH 降按总 β 计（早先只计亚硫酸盐，约大一倍）。初稿的量级脚本误把 γ_H 也算进位移（−4D），得 10.0 mM/pH 与 0.70 的表面 pH 降；已按 −3D 重算并更正。这个差别本身说明 pH 标度的约定会改变缓冲容量近一倍，实现必须严格按 §2.3.3。

### 12.6 原型柱与时间步

一维原型（6 个凝胶单元、10 个液体单元、两组分加一个银量状态，r = k C/(1 + K_B Br)(1 − X)，400 s，每 30 s 复位液体）：SciPy BDF（rtol 1e-8）作参考，4566 步；结果见 §3.4 的表。传输算子的特征值最快 68 s⁻¹、最慢 0.022 s⁻¹，比值 3.1e3。粗离散（n_z = 2、n_l = 4，自适应 tol 3e-3）的平均显影分数 0.2296，细离散 0.2307，相对差 0.5 %。

修订核对（数值评审之后）：按 §2.4.8 的搅动期（每 30 s 一次 5 s、大罐每 60 s 一次 7.5 s）与旧规则（起止与快照都复位），图像档 7.5 min 要 293 步、12 min 473 步，不是原型的 178 步；§3.4 的新规则（只在搅动开始复位，止点与快照只对齐）给 209 与 335 步，表见 §3.4（脚本 `sched2.py`，在会话临时目录）。潜影表的插值（§2.2.5）：Poisson r 次命中极限、η = 0.5，201 点线性插值的中点误差 2.25e-3，not-a-knot 三次样条 9.0e-6，401 点 5.3e-7（`lat.py`）。

### 12.7 光学

- 光晕：TAC（n = 1.475）θ_c = 42.7°；t_b = 127 µm 时 r_c = 234 µm。D_g = 0.10：返回分数 0.222（全反射部分 0.199），99 % 与 99.9 % 能量半径 1167 µm 与 1831 µm。透明片基：0.581，其中全反射部分 0.540 = 1 − 1/n_b²，其余为 θ < θ_c 的菲涅耳反射；99 % 与 99.9 % 半径 3324 µm 与 10532 µm。T-Max 120（119 µm）r_c = 220 µm。与 OPT §3 一致。
- 指数 LSF 的 PSF：尾部质量 (R/a)K₁(R/a) < 1e-3、1e-4、1e-6 分别需要 R = 8.24a、10.66a、15.44a。
- 布尔重叠因子（常半径）：κ(0.1, 0.3, 0.5, 1.0, 1.5, 2.0) = 1.055、1.182、1.336、1.891、2.861、4.655；48 点 Gauss–Legendre 给 1.8913（D = 1）、4.6548（D = 2）。与 OPT §5.4 一致。

### 12.8 核函数耗时

本机单进程（§3.10 的表）。另测：`scipy.linalg.expm` 对 36×36 为 40 µs，一张 K = 6、201 点、9 节点的潜影表（10 854 次）约 0.43 s。

### 12.9 指标算法在训练曲线上的核对

按 §6.3、§6.5 的定义（PCHIP、Δx = 0.002 找交点、brentq 求共线点；以曲线最左点作基雾代理）直接读暂存的训练曲线：

| 曲线 | CI | log H_m | G_ISO | 参照 |
|------|----|---------|-------|------|
| T-Max 100 × D-76，6 min | 0.542 | −2.011 | 0.548 | CAL §1.3.4：0.542、−2.011、0.547 |
| 同上，7.5 min | 0.694 | −2.129 | 0.701 | 0.693、−2.129、0.701 |
| 同上，10 min | 0.895 | −2.213 | 0.879 | 0.895、−2.212、0.879 |

Tri-X D-76 的 ΔD（6/8/10/12 min）：35 mm 为 0.585/0.714/0.852/1.017，120 为 0.569/0.690/0.817/0.964；在 8 与 10 min 之间内插到 ΔD = 0.80，得 log H_m = −2.726（35 mm）与 −2.688（120），DS §5.3 为 −2.73 与 −2.69。差别都 ≤ 0.001，所以 `test_ci_matches_report_on_tmax_d76` 与 `test_iso_speed_on_kodak_tri_x` 的容差（0.01/0.02/0.03）留有余量，余量用于数字化之间的差别。

### 12.10 标定修订的核对（D_min 敏感度、F-4016 两图之差、温度系数）

用 `tools/digitize/contrast.py` 的共线作图，直接读暂存 CSV：

- CI 对 D_min 的敏感度：把 D_min 从曲线左端密度降 0.01 / 0.02 / 0.04 D，T-Max D-76 三条曲线的 CI 降 0.003–0.005 / 0.007–0.011 / 0.018–0.026，Tri-X 八条 D-76 曲线降 0.005–0.007 / 0.011–0.016 / 0.026–0.041。
- F-4016 两图之差的折中（只用训练曲线）：把 T-Max D-76 特性曲线绕使密度 RMS 最小的支点在 log H 上伸缩，使其 CI 从曲线值向下移 f × 差；在 §7.2 的 15 点窗口上算对原曲线的密度误差。差取暂存 fact 中已印出的三个量级 0.036（训练对 D-76 7.5 min）、0.043 与 0.065（fact 里 T-MAX 的两个数值，只借用其大小，作用在训练曲线上）。

| 差 | f = 0.5 时对 CI(t) 的 CI 误差 | 三条 D-76 曲线的密度 RMS / 最大（f = 0.5） | f = 1 时密度 RMS / 最大 |
|----|------------------------------|-------------------------------------------|--------------------------|
| 0.036 | 0.018 | 0.013–0.015 / 0.023–0.026 | 0.026–0.030 / 0.046–0.052 |
| 0.043 | 0.021 | 0.016–0.018 / 0.027–0.031 | 0.032–0.036 / 0.055–0.063 |
| 0.065 | 0.033 | 0.024–0.027 / 0.042–0.047 | 0.048–0.055 / 0.084–0.096 |

f = 0.5 时：差 0.036 下训练阈值（CI RMS 0.03、max 0.05；密度 RMS 0.03、max 0.08）都满足；差 0.065 下留出阈值（CI 每点 0.07；密度 RMS 0.06、max 0.15）都满足。只贴合一张图（f = 1）时训练的密度 RMS 已到 0.026–0.030。所以 §7.2 的 σ_floor = 0.018 使拟合折中，§8.3 的阈值不需要放宽。
- F-4016 `[d-76]` 块 Normal 列对 T-Max D-76 各时间表（小罐、大罐卷片与页片、盘、转筒卷片与页片，18/21/22/24 °C 共 24 格）：系数 / (t_T / t_20) − 1 的 RMS 3.9 %、最大 7.1 %（小罐 24 °C）。在 dCI/dt ≈ 0.09 CI/min、t ≈ 6.5 min 处，7.1 % 约合 0.04 CI。



---

## 13 修订记录（草稿之后）

草稿之后按物理、数值、标定三方面的评审做了一轮修订。下表每行一项，列出改动所在的章节；决策本身见各节正文。

**物理**

| 改动 | 章节 |
|------|------|
| TOTH 不再按菲克律整体输运：它仍是每单元的守恒库存（含明胶结合质子），但面通量由各可动族的质子量 P_f 与 H⁺、OH⁻ 合成，明胶不跨面；新增 D_H = 8130、D_OH = 4600 um²/s；J 只保留 ∂J_toth/∂TOTH；横向扩散作用在 TOTH_mob 上，用逐层单一 D；新增 `test_gel_equilibrium_no_toth_flux`、`test_immobile_buffer_retards_ph_front`、`test_lateral_toth_single_d` | §2.3.1、§2.3.7、§2.4.7、§2.4.8、§2.4.11、§2.4.13、§3.3、§3.5、§5.3 |
| 统计模式方差：Var_count 改用实际抽样的 q_c（x̄ 而非 ⟨x²⟩），Var(ε) 的分母改为 T̄² + Var_count，使截断前像素方差恰为 Var_bool；`test_siedentopf_limit` 改述并改在 D = 0.001 检查；新增 `test_statistical_variance_broad_x`、`counted_classes`、`epsilon_variance` | §2.7.1、§2.7.4、§6.6、§9 第 20 条、§10 D1d |
| 适用条件取消：p < 4 r_c 的档改为均值档（不抽泊松数，方差进 ε），不再报错；记录 b_ap 与截断比例；§6.6 增加标定前的 Siedentopf 合理性检查 G_S；删除 `test_validity_condition`，新增 `test_mean_only_classes`、`test_trix_prior_centre_statistical` | §2.7.1、§2.7.4、§6.6、§9 第 20 条、§10 D1d |
| 物理显影银按比例 φ 放大已显影的盘（统计与显微两种模式），使 0.4343 Σ μ q = D_Ag 含物理银 | §2.7.1、§2.7.2、§2.7.4 |
| 涂布银量 m_Ag,coat 由固定改为每种胶片拟合（对数先验中心 4.0，界 2–6）；不可辨识时做 {2, 6} 的敏感性重拟合 | §2.0.2、§2.5.3、§5.8、§7.1、§7.5、§7.9、§9 第 32 条 |
| T-Max 颗粒厚度 t_g 由假设改为拟合（0.08–0.15）；T-Max 用 (t_g, AR) 参数化，d_50 的 [0.6, 1.5] 改为单侧罚项；新增 S/V 风险条目 | §2.0.2、§7.1、§7.5、§9 第 19、31 条、§10 calibrate/parameters.py |
| 亚硫酸盐溶剂指数：ν_s 进入共用 C 块拟合（0.6–2，初值 1.0），ν_b 改为外推假设 0.4（0.3–1）；C_sat,ref 的锚点改写为 Carroll & Hubbard 30 °C 的原始条件；共用块 10 → 11、化学块 35 → 36 | §2.4.11、§7.1、§7.5、§9 第 9 条 |
| 溴化物 D_w 1420 → 1816 um²/s（示踪值），碘化物 1786；晕边公式加入搅动残留因子 1/(1 − ρ_k) 并遍及 TOTH_mob 的横向 D；小罐晕边 875 → 1047 µm，晕边 37 → 44 格；`test_halo_convergence` 覆盖三种搅动与先验上界 | §2.3.7、§3.6、§3.8、§3.10、§4.6 |
| 盐类的凝胶扩散比 r_gel 0.8 → 0.5（assumed，界 0.3–0.8，向高明胶浓度外推）；Br⁻ 的凝胶扩散时间相应改为 0.11/0.25/0.44 s | §2.3.7（注 †）、§2.4.12、§5.8 |
| 缓冲容量按全部族重算（条件值总 β 33.3 mM/pH，理想 92.0），表面 pH 降与封闭单元 β 范围随之更正；测试数值不变 | §2.4.12、§9 第 2 条、§12.5 |
| 卤化物吸附常数按 van 't Hoff 随温度变，新增 ΔH_Br = −45 kJ/mol（literature，固定） | §2.4.5、§2.4.11、§5.3、§10 chemistry/kinetics.py |
| `test_arrhenius_scaling` 的条件补齐（形态分布与温度无关、ΔH_Br = 0、无溴化物） | §2.4.13 |
| 高照度闭合 q₁ 改称 LAT §1.3(c) 方案 (ii) 型的平均场闭合，写明饱和在 S_eff − 1 倍；方程与测试不变 | §2.2.2、§2.2.6、§9 第 17 条 |
| 光谱感光度的能量换算明确用 D0 的 `photon_energy_j`（J），新增 `test_spectral_energy_conversion` | §6.9 第 3 步、§10 D1e |

**数值**

| 改动 | 章节 |
|------|------|
| 液层分两段几何网格，δ_agit 恰为单元界面；搅动期狄利克雷面在 δ_agit，上段取槽液值并在搅动开始时把库存差混入槽；`column_grid` 增加 `delta_agit_um`；新增 `test_strip_smooth_in_transport` | §2.4.1、§2.4.8、§2.4.9、§2.4.13、§3.1、§3.2、§10 develop/grid.py |
| 成员移出时间积分：柱只积分按团簇大小分解的核 G、H、Q（精确），成员在开始聚合 p̄、在结束与快照由核双线性插值合成；`Members` 改为 columns/raster 两种布局；新增 `test_kernel_reduction_exact` | §2.4.1、§2.4.2、§2.4.4、§2.4.6、§3.1、§3.3、§3.5、§3.8、§4.6、§10 develop/ensemble.py、§1.4 |
| 删除二维 F 表，F 与 ∂F/∂Λ 精确求值；删除 `test_f_table_interpolation` 与 `FilmColumnModel` 的表字段 | §2.4.2、§3.1、§3.3、§10 FilmColumnModel |
| 潜影表改用 not-a-knot 三次样条插值；`test_latent_table_interpolation` 收紧到 1e-4 | §2.2.5、§2.2.9、§12.6 |
| 时间表只在 0 与搅动开始复位，止点、快照与终点只对齐；步数表重算；`test_fixed_schedule_vs_reference` 覆盖小罐与大罐、两个档位 | §3.4、§2.4.13、§12.6 |
| 负浓度不再拒步：速率式用光滑正部 C⁺，记录 min C/s_k；`EnsembleResult.n_halvings` 改为 `min_scaled_concentration`；新增 `test_dry_film_bromide_free_start` | §2.4.6、§2.4.10、§3.3、§10 EnsembleResult |
| 药液槽检查改为按 65 根预计算柱的显影银判定（δ ≤ 1e-3），失败时 `FrameAssembler.reset()` 整幅重算，最多再 2 遍；不再断言 L；新增 `PrepassCheck` 与 `prepass_check` | §3.7、§2.4.10、§4.4、§3.8、§10 develop/prepass.py、§10 D1a state.py |
| 化学单元固定为 24 µm，与像素间距无关；像素按面积重叠聚合并按中心归属块；`test_resolution_independence` 增加 9 与 30 µm | §2.4.1、§3.1、§3.6、§3.8、§4.6 |
| 24 MP 与标定预算重估（约 7500 s 单核，9 进程约 21 min；标定约 2–2.5 天）；`test_benchmark_image_path` 改为按实测并行效率外推、30 min 上限；D1c 条带 ≤ 2 s | §3.10、§7.9、§10 D1c、D1d |
| `FieldStore` 按存储名（`s{i}:{key}`）分配、文件名唯一，新增 `FieldDescriptor`、`open_field`；工作进程直接写输出、只回传摘要；`FrameAssembler` 增加 `name_prefix`、`descriptors`、`mark`、`reset`；`RunContext.store`；每进程约 0.5 GB，默认 9 个工作进程 | §4.3–§4.7、§10 D1a |
| D1c 的测试改对显影银 M 断言（ΔM = ΔD/0.7）；新增 D1e 的 `test_calibration_profile_vs_reference_metrics` | §3.8、§3.9、§10 D1c、D1e |
| 夹具分阶段增长（D1a 只含层、颗粒、光学）；toth 与 gelatin 是内置伪组分，不写进 `species.yaml`；D1d 的夹具库复制所需物种并带先验的 `kinetics.yaml` | §2.3.1、§5.2、§10 D1a、D1d |
| 阶段 C 的 MTF 改用小信号柱对（`Lateral.boundary = "pair"`）；颗粒度残差用闭式期望值 G_exp；内层求解固定迭代次数；稀疏雅可比；多起点 8 × 10 再取最好 2 个；剖面限于 4 + 6 个参数、8 个格点 | §6.6、§6.7、§3.5、§7.2、§7.3、§7.4、§7.5、§10 D1c、D1e |
| `run_ensemble(record_cells=…)` 给出 `CellHistories`（移到 develop/ensemble.py，记 K_z = ∫κ̃ dt 的 64 个时刻）；新增 `test_cell_history_thinning` | §2.7.2、§10 D1c、D1f |
| 逐块物理代码禁用依赖内存对齐的 BLAS 归约（`np.dot` 等）；新增 `test_block_result_independent_of_alignment` | §3.6、§10 D1c |
| 若干测试的定义更正：`test_scatter_kernel_otf`（与盒滤波并混叠的 OTF 比）、`test_diffusion_into_slab`（收敛阶）、`test_moving_reaction_front`（专用细网格）、`test_uniform_image_matches_strip`（只比内部像素，≤ 1e-10） | §2.1.10、§2.4.13 |

**标定**

| 改动 | 章节 |
|------|------|
| 与 Kodak 曲线的一切比较两侧都用 D_min* = D_ref(x_left) 与 log H_m*；基雾行改为左端密度加单侧基雾上界；新增 `d_min_star` 与 `test_dmin_star_same_operator` | §6.3、§6.5、§7.2、§8.3、§10 validate/metrics.py |
| "到达 CI 的时间"残差改为表格格的 CI（`table_ci`，σ 含取整误差与 √n_class）；撤销 F-4016 HC-110 系数 15 % 阈值；新增 `ci_at_times`，`time_to_ci` 只作报告 | §6.10、§7.2、§7.6、§8.1、§8.3、§10 |
| 相关性按 σ 放大处理（曲线 √(15/6)、CI(t) √(5/3)、T-Max 的 σ_floor = 0.018）；s² 按 N_eff；非有限残差取有界罚值并计数；ISO 速度点移出目标函数；`Observation`、`StageRecord` 增加相应字段 | §5.7、§7.2、§7.8、§10 validate/reference.py |
| 组合状态按预注册的组合–指标表判定，算不出的适用指标判 failed；新增 `PAIR_METRICS`、`test_held_out_pipeline_finite_on_synthetic`、`test_nan_metric_fails_pair` | §8.2、§8.3、§10 validate/splits.py |
| 注册表分为 `TRAINING`、`WITHIN_CHECKS`、`HELD_OUT`、`UNUSED` 四个，标识细到子项；新增 `TrainingLibrary`、`HELD_OUT_EVALUATED`、`NOT_EVALUATED_REASON`、`registry_sha256`；不评估的留出项列出原因 | §7.6、§7.7、§10 D1e |
| F-4016 温度系数：`[d-76]` 块改为内部检查，反差列一律不评估，留出块的 Normal 列换算为表格格评分 | §7.2、§7.6、§7.7 第 6 条、§8.3、§12.10 |
| 可辨识性只用数据残差（`residuals(..., include_prior=False)`），新增先验主导判据；剖面与协方差按块、以另一块为条件；阶段 F 固定于阶段 E 值；预测带分两种 | §5.7、§7.3、§7.5、§7.8、§10 calibrate/ |
| 冻结与唯一一次留出评估移到 D1f 之后；`CalibrationRecord` 增加阈值与哈希字段及 `pair_statuses`，新增 `PairStatus`；新增 `test_held_out_refuses_mismatched_freeze` | §5.7、§7.3 阶段 G、§8.2、§8.4、§10 D1e、D1f、§11 第 10 条 |
| 胶片块加入 S_550 与只在标定中存在的吸收尺度 s_abs，D_g 不拟合；A_550 改为单侧罚残差；化学块 A 加入 Carlton & Crabtree 系列的偏移 δ_CC | §2.0.2、§2.1.8、§7.1、§7.2、§7.3、§7.5、§10 calibrate/parameters.py |
| Carlton & Crabtree 残差加 δ_CC、σ 统一为 0.10；全部 pH 锚点残差保留在阶段 D、E | §7.2、§7.3、§7.5 |
| MTF 分为尾段（阶段 C 拟合）与峰段（横向化学的内部检查）；两种胶片的 MTF 都按 135 模拟 | §6.7、§7.1、§7.2、§7.6、§8.1、§8.3、§11 第 11 条 |
| 光谱拟合只取参考曲线的数字化跨度；Tri-X 文件按名使用；新增 T-Max 光谱绝对水平的内部检查 | §5.8、§6.9、§7.2、§7.6、§8.1、§9 第 22 条 |
| τ̄₀、σ_τ 的界只由 LAT 导出，Tri-X 示例 τ̄₀ 0.07 → 0.22 s；倒易律未拟合行与 Tri-X 显影修正归内部检查 | §2.2.7、§5.8、§7.6、§8.3、§11 第 8 条 |
| 每个留出组合预注册一个标题数字，并与只用训练数据的朴素基线并列 | §8.1 第 3 项、§5.7 PairStatus |
| 设计期暴露清单补全，泄漏规则第 3 条改写 | §7.6、§7.7 |
| 新增设计期核对：CI 对 D_min 的敏感度、F-4016 两图之差的折中、`[d-76]` 系数与时间表之差 | §12.10 |

**一致性修订（合并三方改动时）**

| 改动 | 章节 |
|------|------|
| d_50 一行删去"t_g 是 assumed 常数、下界 0.96"的说法，与 t_g 拟合及 (t_g, AR) 参数化一致 | §2.0.2 |
| G_exp 中的逐像素方差改为 Var_count + (T̄² + Var_count) Var(ε)，与 §2.7.1 的新分母一致 | §6.6 |
| `profile_likelihood` 的缺省格点 10 → 8 | §10 calibrate/identifiability.py |
| `Kinetics.K_br`、`K_i` 与 `EnsembleResult.flux_mean` 的注释补上温度与 TOTH 通量的含义 | §10 D1c |
| D1c、D1e 的测试清单补上 `test_lateral_toth_single_d`、`test_dmin_star_same_operator`、`test_spectral_energy_conversion` | §10 D1c、D1e |
| 标定记录示例补 `n_effective`，T-Max 的阶段 C 参数改为 t_g 与 `calibration:aspect_ratio` | §5.8 |
| 第 11 节条目恢复 8、9、10、11 的顺序 | §11 |
