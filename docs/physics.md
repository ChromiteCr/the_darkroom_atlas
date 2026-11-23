# 物理模型与约定

这份文档记录模型用到的约定和方程，以及每一处不物理的地方。方程必须与代码一致；改代码时同步改这里。

第 1 节和第 3 节是 D0 的约定与物理；第 4 节记录 D1 已完成的模块（颗粒分档、光传输与曝光、潜影、药液化学、内存设计）；显影、成像、颗粒、读出和标定完成后再补。

## 1 约定

### 1.1 单位

| 量 | 单位 | 说明 |
|----|------|------|
| 波长 | nm | |
| 胶片面上的长度 | µm（代码里写 `um`） | 扩散长度、光晕半径、颗粒大小都用它 |
| 曝光量 | 每 µm² 每波段的光子数 | |
| 光通量 | 每 µm² 每秒每波段的光子数 | `Light.photon_flux` |
| 密度 | 以 10 为底的光学密度 | `D = -log10(T)` |
| 时间 | s | |
| 温度 | °C 或 K（代码里写 `C`、`K`） | |

单位字符串只用 ASCII：`um` `nm` `mm` `s` `min` `C` `K` `mL` `L` `g/L` `g/mol` `g/m2` `1`。

像素间距由画幅和像素数决定：

```
pixel_pitch_um = long_edge_mm * 1000 / max(H, W)
```

所有物理尺度都以 µm 表达，所以同一画面换一个像素数来算，物理结果不变，只是采样密度不同。

用到的常数（`units.py`）：普朗克常数 `h = 6.62607015e-34 J s`，光速 `c = 299792458 m/s`，阿伏伽德罗常数，气体常数，玻尔兹曼常数，以及最大光视效能 `K_m = 683 lm/W`。

### 1.2 光谱网格与"按波段积分"

默认网格是 380–780 nm、步长 5 nm，共 81 个波段，波段中心在 380、385、…、780 nm。

光谱量按"波段内已积分"存放：数组里的每个数是这个波段里的总量，不是每纳米的密度。因此对波长积分就是对最后一维求和，不需要再乘步长。从每纳米的表格数据（例如标准光源）转过来时，先重采样到波段中心，再乘以步长。

无量纲的光谱量（透射率、反射率、光谱敏感度、色匹配函数）在波段中心取值，不乘步长。

标准光源归一化为"理想白的亮度等于 1"：

```
sum_b( S_b * ybar_b ) = 1
```

其中 `S_b` 是光源在波段 `b` 的相对能量，`ybar_b` 是 CIE 1931 色匹配函数的 y 分量。

### 1.3 数组轴顺序

固定为 `(H, W, ..., band)`：前两维是画面的行和列，波段永远是最后一维。内部一律 `float64`。

### 1.4 状态与字段命名

乳剂的状态是 `EmulsionState`，其中的场用 `"<层 id>/<名字>"` 命名，例如 `y/exposure`。不属于任何一层的场用保留前缀 `frame/`。

状态不可变：每个操作返回新的状态，数组设为只读。状态里带着处理历史（每一步一条记录）。

### 1.5 随机数与可复现性

所有随机数都从 `rng.py` 取。一个随机数只由（种子、流名、像素的绝对行列坐标、尾部下标）决定，与画面怎样分块、各块以什么顺序计算无关。

实现：每个绝对行 `r` 有自己的 PCG64 发生器，种子为 `SeedSequence([seed, stream_key(seed, stream), r])`，其中 `stream_key` 是 `blake2b("{seed}:{stream}")` 的 8 字节摘要。只取发生器的原始 64 位输出，换算为 `[0, 1)` 上的均匀数：

```
u = (raw >> 11) * 2**-53
```

不使用 NumPy 的分布函数，所以数值不随 NumPy 版本变化。

横向比较不同胶片或显影液时，各组合用同一个种子，差异只来自物理。

每次运行写一份运行记录（`.run.json`），含输入与输出文件的哈希、所用数据文件的哈希、版本、种子、曝光参数和运行环境。数据文件包括材料、流程、药液、物种表，以及材料引用的每一个光谱曲线文件（`curve_files`，不论该材料的模型是否读取它）。

其中的 `recipe_hash` 由决定结果的内容算出。它不含创建时间 `created_utc`、运行环境 `environment`、git 提交 `git_commit`、分块大小 `tile_px`，也不含输入和输出文件的路径 `input.path`、`output.path`（两者的哈希保留）。数据文件的路径按相对数据库根目录的形式记录，并计入哈希：同一份数据放在两个检出目录里得到相同的 `recipe_hash`，但在库内移动或改名一个数据文件会改变它。两次运行的 `recipe_hash` 相同就应当得到相同的图像。

### 1.6 分块

大图可以分块计算。每个处理步骤的实现要声明它的横向作用范围 `interaction_range_um`，有三种取值：

| 取值 | 含义 | 分块方式 |
|------|------|----------|
| `0.0` | 逐点运算，像素之间互不影响 | 块之间不需要重叠 |
| 正数 | 影响范围有限（例如扩散、散射） | 块向外扩 `ceil(interaction_range_um / pixel_pitch_um)` 个像素，算完只取中间 |
| `None` | 全画面耦合（例如有限体积的药液槽） | 管线不分块，整幅交给该步骤，由它自己决定怎么算 |

要求：分块与不分块的结果逐位相同。

## 2 链路总览

| 环节 | 内容 | 状态 |
|------|------|------|
| 0 场景 → 光谱曝光 | RGB 还原为光谱；虚拟相机按胶片感光度测光 | D0 已实现 |
| 1 乳剂内的光传输 | 逐层光谱吸收、散射（点扩散函数）、片基内反射造成的光晕 | D1 已实现（第 4.2 节） |
| 2 潜影形成 | 颗粒按尺寸分档；光子吸收服从泊松分布；成核-生长模型，倒易律失效由此产生 | D1 已实现分档与潜影链（第 4.1、4.3 节），尚未接入显影 |
| 3 显影 | 由配方得到的反应-扩散：显影剂离子化、溴化物抑制、温度、pH 场、超加和性、搅动与药液耗尽 | 药液化学已实现（第 4.4 节），显影本身计划于 D1 |
| 4 成像 | 黑白：银量 × 遮盖力 → 密度。彩色：成色、层间效应、色罩、反转 | 黑白计划于 D1，彩色计划于 D2 |
| 5 停显、定影、漂白、水洗 | 一阶模型，作为真实步骤改变状态 | 计划于 D1（漂白于 D2） |
| 6 颗粒 | 全幅用二项分布抽样；小块区域逐颗粒模拟（显微镜模式） | 计划于 D1 |
| 7 出口 | 密度计、标准扫描、标准印相 | 密度计与标准扫描计划于 D1，标准印相计划于 D3 |

D0 的管线已经按这条链路的形状搭好：`expose`、`process`、`readout` 三个操作，加上"材料的 `model` 字段选择物理模型"的注册表。D0 只注册了一个模型，即第 3.4 节的理想材料；D1 另注册了曝光模型 `bw-mechanistic`（第 4.2 节）。

底片读出的光谱透射率可以乘到一束光上（`Light.filtered`），再对下一种材料曝光。暗房印相因此复用同一套引擎。

## 3 D0 已实现的物理

### 3.1 RGB 还原为光谱

输入图像先变成线性 sRGB，再用 Mallett & Yuksel 2019 的三条基光谱还原为反射率：

```
R_b = r * B_b,red + g * B_b,green + b * B_b,blue
```

三条基光谱非负，在每个波段上相加等于 1，所以 `(1, 1, 1)` 给出处处为 1 的平直反射率。这样还原出的光谱在 D65 光源和 CIE 1931 观察者下的颜色等于原来的 sRGB 值，只差一处白点。

D0 测得的往返最大误差约为 `2.7e-4`，全部来自白点的差异：sRGB 标准把 D65 白点写成四舍五入的 xy = (0.3127, 0.3290)，读出矩阵按它算出；而网格上的 D65（CIE 表格光谱与色匹配函数按 5 nm 波段求和）白点是 xy = (0.312721, 0.329031)。所以理想白读回为 (1.000043, 1.000015, 0.999727)。改用网格自身的白点算矩阵，往返误差降到约 `3.5e-8`（见测试 `test_roundtrip_error_is_the_srgb_white_point`）。D0 保留标准矩阵。重采样不引入误差：默认网格的每个波段中心都是表格数据（色匹配函数、标准光源、基光谱）本身的取样点，重采样只是取值，不做插值。

实现上有一处调整：表格里的基光谱在各波段相加只精确到约 `1e-5`，代码把每个波段除以三者之和，使白色严格等于 1。

局限见第 5 节。

### 3.2 光度量到光子数

一个光子的能量：

```
E_b = h * c / lambda_b
```

场景的相对能量是反射率乘光源：`e_b = R_b * S_b`。由于光源已经归一化，理想白的 `sum_b(e_b * ybar_b) = 1`。

设比例系数 `k`（单位 J/µm²，对应相对能量为 1）。某个像素在胶片面上的光度曝光量是

```
H = K_m * 1e12 * k * sum_b( e_b * ybar_b )      [lux s]
```

`K_m = 683 lm/W`；`1e12` 是把每 µm² 换成每 m²。

曝光时间 `t` 内到达的光子通量（每 µm² 每秒每波段）：

```
Phi_b = k * e_b / ( E_b * t )
```

曝光量（光子数）是 `N_b = Phi_b * t`，与 `t` 无关：改曝光时间只改通量，不改曝光量。倒易律失效是胶片对通量的响应，从 D1 的潜影模型开始才会出现。

### 3.3 测光

虚拟相机决定 `k`，与原数码照片的拍摄参数无关。

平均测光（默认）：让全画面 `H` 的平均值等于

```
H_mean = (10 / S) * 2**ev      [lux s]
```

`S` 是胶片的感光度，`ev` 是曝光补偿。`H = 10 / S` 是设计文档从 ISO 2721 取的平均胶片面曝光量关系。这个常数决定了模拟出来的感光度是否对得上标称值，D1 标定感光度时要对照标准原文再核一次。

全黑画面（平均相对亮度低于 `1e-9`）无法测光，按 18% 灰处理。

固定曝光：直接指定理想白得到的曝光量 `H_white`，此时 `k = H_white / K_m * 1e-12`，感光度和曝光补偿不起作用。虚拟感光仪用这种方式。

### 3.4 理想材料

理想材料不是物理模型。它有三个感光层 `x`、`y`、`z`，光谱敏感度分别是 CIE 1931 的三条色匹配函数。

曝光：每层记录它收到的、按敏感度加权的能量

```
exposure_layer = sum_b( N_b * E_b * s_layer,b )      [J/um^2]
```

冲洗（`identity` 步骤）：`signal = exposure`，逐点。

读出（`ideal_colorimetric`）：

```
XYZ = signal / k
RGB_linear = M * XYZ
```

`M` 是 XYZ 到线性 sRGB 的矩阵，由 IEC 61966-2-1 的基色和 D65 白点算出（没有用通常印出来的四位小数版本，那一对矩阵互逆只精确到 `4e-5`；白点用的是标准里四舍五入的值，见 3.1 节）。

为什么这样能还原：场景光带着 `N_b * E_b = k * e_b`，所以三层的信号正好是 `k` 乘以场景在该光源下的 XYZ，与感光度、曝光补偿、曝光时间都无关。除以 `k` 就得到 XYZ，理想白的 `Y = 1`。

它检验了什么：光谱还原、光源归一化、色匹配函数、曝光对波段求和、读出矩阵。这些环节只要有一处不一致，图像就不会原样返回。

它不检验什么：绝对光度量。`k` 在相机里乘上、在读出时除掉；光子能量在相机里除掉、在理想模型里乘回。所以 `K_m`、µm² 与 m² 的换算、测光常数 `10 / S`、曝光时间，以及相机和模型共用的任何逐波段因子（例如写错的光子能量），在往返中都会整体抵消，图像照样原样返回（见测试 `test_roundtrip_is_blind_to_absolute_photometry`）。这些量由单元测试单独锚定：

- 光子能量：`test_photon_energy_555nm`，以及逐波段对照 `h * c / lambda` 的 `test_photon_energy_across_the_grid`；
- 光度换算：`test_luminous_exposure_anchor`（555 nm 处每 µm² 4090.7 个光子对应 1 lux s）；
- 测光：`test_average_metering_hits_iso` 和 `test_fixed_white_exposure`；
- 理想材料的 y 层与光度曝光量一致：`test_y_exposure_matches_luminous`。

它也不检验任何胶片物理。

## 4 D1 已实现的物理

这一节只写代码已经完成的 D1 模块：颗粒分档、光传输与曝光、潜影、药液化学，以及内存设计。显影、成像、停显定影、颗粒、读出和标定还没有完成，等它们完成后再补，这里不写。

参数的来源类别沿用 D0 的四类：`datasheet`、`literature`、`fitted`、`assumed`。有两点要先说明，免得读错：

- **标定还没有做。** 表里标为 `fitted` 的参数，目前只是有先验区间、等着被拟合的量：本文写作时，`data/films/kodak/` 里的 `bw-mechanistic` 胶片文件对这些参数只给了 `assumed` 的起始值，没有一个值被标为 `fitted`；测试用的胶片在 `tests/fixtures/library/` 里。所以下面各节的数值结果检验的是实现与方程是否一致，不是对任何一款真实胶片的预测。
- **测试里的 pH 锚点是实现的回归，不是验证。** 它们来自设计期的手算（规格 §12），用来确认代码算的是文档写的方程。

### 4.1 颗粒尺寸分档

代码：`emulsion/classes.py`，入口 `size_classes(layer, n_classes=None)`，返回 `SizeClasses`（所有数组为只读 float64，长度 `K`）。

这一步做什么：把乳剂层的颗粒群体参数变成 `K` 个尺寸档，给出每档的单颗几何、单颗含银量和单位胶片面积上的个数。光传输（4.2）用它算吸收，潜影（4.3）用每档的直径和体积，后面的显影与颗粒统计也读同一个对象。

颗粒的等效圆直径（ECD）按个数服从对数正态分布，在 ±3σ 处截断后重新归一，再按标准化对数直径等分成 `K` 档：

```
ln d ~ N(ln d_50, sigma_d^2)                       截断于 ±3 sigma
t_k = -3 + 6 k / K,   k = 0 … K                    档界（标准化对数直径）
w_c = [Phi(t_c) - Phi(t_(c-1))] / [Phi(3) - Phi(-3)]       个数分数，c = 1 … K
d_c = d_50 * exp( sigma_d * (t_(c-1) + t_c) / 2 )         档中心（对数直径的中点）
```

`Phi` 是标准正态分布函数，`d_50` 是 ECD 的个数中位数（µm）。截去的个数分数是 0.27 %。

单颗几何（`cubic` 按 `conventional` 处理；`tabular` 是直径 `d_c`、厚度 `t_g` 的圆片）：

| 量 | conventional（球） | tabular（圆片） | 单位 |
|----|-------------------|-----------------|------|
| 体积 `V_c` | π d_c³ / 6 | π d_c² t_g / 4 | um³ |
| 表面积 `S_c` | π d_c² | π d_c² / 2 + π d_c t_g | um² |
| 投影面积 `A_c` | π d_c² / 4 | π d_c² / 4 | um² |

含银量、个数密度和干乳剂厚度：

```
M_AgX   = (1 - x_I) * 187.77 + x_I * 234.77                    [g/mol]
m_Ag,c  = 1e-12 * rho_AgX * V_c * 107.87 / M_AgX               [g，每颗粒中的银]
N_tot   = 1e-12 * m_Ag,coat / sum_c( w_c * m_Ag,c )            [um^-2]
nu_c    = w_c * N_tot                                           [um^-2，c 档的面密度]
L_e     = m_Ag,coat * (M_AgX / 107.87) / rho_AgX + m_gel / rho_gel     [um]
```

`x_I` 是碘化物摩尔分数；`rho_AgX` 单位 g/cm³（1 g/cm³ = 1e-12 g/um³）；`m_Ag,coat`、`m_gel` 单位 g/m²（1 g/m² = 1e-12 g/um²）；g/m² 除以 g/cm³ 正好得 µm。

代码里的常数：

| 常数 | 值 | 类别 | 说明 |
|------|----|------|------|
| `M_AG` | 107.87 g/mol | 原子量 | |
| `M_AGBR` | 187.77 g/mol | 由原子量算出 | |
| `M_AGI` | 234.77 g/mol | `assumed` | 由原子量算出 |
| `RHO_AGX` | 6.473 g/cm³ | `assumed` | AgBr 的密度，碘化物对密度的影响忽略 |
| `RHO_GELATIN` | 1.35 g/cm³ | `assumed` | 干明胶密度；规格没有给数值，这是它的 4.1 µm 算例背后的值 |
| `TRUNCATION_SIGMA` | 3 | `assumed` | 截断位置 |

档数 `K`（`grain.n_classes` 可被参数覆盖，范围 2 到 64）由离散化档位给出：`reference` 12，`calibration` 6，`image` 6（`models.bw.PROFILE_CLASSES`）。

| 符号 | 字段 | 单位 | 类别 | 说明 |
|------|------|------|------|------|
| d_50 | `grain.mean_diameter` | um | fitted | 待拟合；T-Max 按 (t_g, 长径比) 参数化 |
| σ_d | `grain.sigma_log` | 1 | fitted | 待拟合 |
| t_g | `grain.thickness` | um | fitted（仅 tabular） | 没有直接数据，待拟合 |
| x_I | `grain.iodide_fraction` | 1 | assumed | 缺省按 0（纯溴化银） |
| m_Ag,coat | `grain.silver_coverage` | g/m2 | fitted | 没有一手数据；它与其他参数的可辨识性问题见第 5 节 |
| m_gel | `gelatin_coverage` | g/m2 | assumed | |

限制情形与检验它们的测试：

| 测试 | 检验什么 |
|------|----------|
| `test_size_classes_conserve_silver` | `1e12 Σ ν_c m_Ag,c` 等于涂布量，相对差 ≤ 1e-12 |
| `test_layer_thickness_formula` | 1.75 g/m² 银加 4.84 g/m² 明胶得 L_e = 4.1 ± 0.05 µm |
| `test_dry_thickness_with_iodide` | 含碘时 L_e 随 `M_AgX` 变化 |
| `test_number_fractions_sum_to_one_and_match_truncated_lognormal` | 个数分数之和为 1，且等于截断对数正态的各档质量 |
| `test_truncated_mass_is_0_27_percent` | 被截去的个数分数 |
| `test_class_centres_are_geometric_means_of_bounds` | 档中心取档界的几何平均 |
| `test_conventional_geometry`、`test_tabular_geometry`、`test_cubic_is_treated_as_conventional` | 单颗几何 |
| `test_fine_classes_approach_the_continuous_mean_volume` | 档数增加时平均体积收敛到连续分布的值 |
| `test_missing_iodide_counts_as_pure_bromide` | 没有碘化物字段时按 0 |

读写：不读写 `EmulsionState`；读 `EmulsionLayer` 的 `grain` 与 `gelatin_coverage`。

不物理之处：颗粒在分数深度上均匀分布；每档只用一个代表直径，档内的体积分布被丢掉；对数正态与 ±3σ 截断是假设；`cubic` 与 `conventional` 不区分；碘化物只改摩尔质量，不改密度，也不改感光性。

### 4.2 光传输与曝光

代码：`optics/stack.py`（`optical_stack`）、`optics/kubelka.py`（`depth_response`）、`optics/halation.py`（`halation_kernel`、`halation_fraction`）、`optics/psf.py`（`scatter_kernel`）、`optics/expose.py`（`exposure_fields`、`uniform_exposure`）、`models/bw.py`（曝光模型 `bw-mechanistic`，函数 `expose_bw`）。

这一步做什么：输入是胶片面上的光子曝光量 `N_b(x) = Φ_b(x) t_exp`（每 µm² 每波段的光子数），输出是每个深度单元 `z` 的两个吸收描述量。c 档颗粒在单元 `z` 吸收的平均光子数是它们的线性组合：

```
m_c,z(x) = V_c * e_V,z(x) + S_c * e_S,z(x)        [每颗粒的光子数]
e_V,z    [um^-3]：每单位卤化银体积吸收的光子数（本征吸收，体积项）
e_S,z    [um^-2]：每单位颗粒表面吸收的光子数（增感染料，表面项）
```

**吸收与散射系数。** 记 `E_ph,b = 1239.84 / lambda_b`（eV，波长以 nm 计），`nu_b = 1000 / lambda_b`（µm⁻¹）：

```
alpha_b = alpha_0                              E_ph,b >= E_0
alpha_b = alpha_0 * exp((E_ph,b - E_0)/E_U)    E_ph,b <  E_0          [um^-1]，Urbach 尾
s_b     = sum_j A_j * exp( -((nu_b - nu_j)/w_j)^2 / 2 )               [1]，染料吸收带
mu_a,b  = (sum_c nu_c V_c / L_e) * alpha_b + (sum_c nu_c S_c / L_e) * s_b     [um^-1]
K_b = 2 mu_a,b ,      S_b = S_550 * (550/lambda_b)^q_S                [um^-1]
```

`alpha`（体积项）只在 `alpha·d ≪ 1` 时成立，粗颗粒在蓝紫段的自遮蔽被忽略；`s`（表面项）是增感染料每单位颗粒表面的吸收截面。

**深度方向：两流（Kubelka–Munk）模型。** 乳剂内的光走几个微米就变成漫射光，与层厚相当，所以用两流近似。`s ∈ [0, L_e]` 是从上表面向下的深度，`I+`、`I-` 是向下、向上的漫射光子通量：

```
dI+/ds = -(K + S) I+ + S I-
dI-/ds =  (K + S) I- - S I+
I+(0)   = T_top * N + r_i * I-(0)          N：入射光子数
I-(L_e) = R_bot * I+(L_e) + J              J：从片基返回、自下表面进入的漫射光子数
R_top = ((n_gel - 1)/(n_gel + 1))^2 ,   T_top = 1 - R_top
r_i   = 1 - T_top / n_gel^2               明胶–空气界面对内部漫射光的反射（含全反射）
R_bot = max(0, 1 - (n_b / n_gel)^2)       乳剂–片基界面对向下漫射光的全反射部分
```

`n_b >= n_gel` 时没有全反射，`R_bot` 取 0（`test_interface_no_total_reflection_into_denser_base`）。

解析解有两种等价形式，按 `beta L_e` 切换：`beta = sqrt(K (K + 2S))`。`beta L_e > 1` 时用指数形式（`R_inf = S / (K + S + beta)`，所有指数项 ≤ 1，不会溢出）；`beta L_e <= 1` 时用 `u = I+ + I-`、`v = I+ - I-` 的双曲形式，它在 `K/S → 0` 和 `K = 0` 处依然良态（指数形式的两个模态在那里合并）。两种形式在切换点互相一致（`test_km_branches_agree_near_switch`）。各深度单元内的积分是解析的。

分别取 `(N, J) = (1, 0)` 与 `(0, 1)`，得到每个深度单元的标量通量 `psi = 2 (I+ + I-)` 的平均 `g+_b,z`、`g-_b,z`，以及进入片基的漫射光 `tau_b = (1 - R_bot) I+(L_e)`（`DepthResponse.plus`、`.minus`、`.into_base`）。单位吸收为 `mu_a * psi = K (I+ + I-)`；入射光的去向满足 `absorbed + escaped_top + into_base = T_top`。

**光晕。** 片基里的光线在片基–空气界面按 `n_b sin θ ≥ 1` 判断是否全反射。向下的朗伯漫射光从一点进入片基，在半径 `r = 2 t_b tan θ` 处返回乳剂，单位入射通量产生的返回辐照度为

```
H(r) = R(theta) * T_g(theta) * cos^4(theta) / (4 pi t_b^2)            [um^-2]
theta = arctan( r / (2 t_b) )
T_g(theta) = 10^( -2 (D_g + D_back) / cos(theta) )            两次斜穿灰染与背涂层
f_H = integral_0^(pi/2) 2 sin(theta) cos(theta) R(theta) T_g(theta) d theta         返回分数
```

`R(theta)` 是非偏振菲涅耳反射率（片基到空气），`theta ≥ arcsin(1/n_b)` 时为 1。只算一次返回（二次返回约为一次的 `f_H · R_inf`，灰片基下估计约 7 %，被忽略）。离散成像素核时，每像素用 8×8 超采样积分，跨过全反射半径 `r_c = 2 t_b tan(theta_c)` 的像素用 32×32；核截断在累计返回量达到 `(1 - 1e-3) f_H` 的半径 `R_H`，上限 3000 µm，不再归一，被截去的部分记入 `halation_truncated`。

**层内散射 PSF。** Frieser 指数线扩散函数对应的二维 PSF，`K0` 为第二类修正贝塞尔函数：

```
P(r)   = (1 - w_P) * delta(r) + w_P * K0(r/a_P) / (2 pi a_P^2)         [um^-2]
P^(nu) = 1 - w_P + w_P / (1 + (2 pi a_P nu)^2)                         nu：cycles/um
```

离散成像素核时，每像素 16×16 超采样，中心像素用等面积圆盘的解析式 `1 - (rho/a_P) K1(rho/a_P)`；K0 项截断在 `(R_P/a_P) K1(R_P/a_P) = 1e-4` 处（`R_P ≈ 10.66 a_P`），截断后归一到 `w_P`。

**组合与卷积。** 两个核都按光谱中性处理（PSF 与深度、波长无关；灰染与背涂层中性），所以先对每个像素做光谱投影，再卷积：

```
e_X,z = P ⊛ [ sum_b W+_X,z,b N_b ] + H ⊛ [ sum_b W-_X,z,b N_b ],      X ∈ {V, S}
W+_X,z,b = w^X_b * g+_b,z ,    W-_X,z,b = w^X_b * g-_b,z * tau_b ,    w^V = alpha,  w^S = s
```

低秩光（4.5 节）有 `R` 个系数通道，且 `R < 2 n_z` 时，改为先对 `R` 个系数通道做卷积再组合（`2R` 次卷积）；否则对 `4 n_z` 个投影通道做卷积。卷积在固定的规范块网格上用 FFT 计算：内块 512 像素（`BLOCK_PX`），晕边取 `ceil(max(R_H, R_P)/p)` 像素（`p` 为像素间距，且不小于任一核的半宽），FFT 尺寸 `next_fast_len(512 + 2 halo)` 对所有块相同，画幅以外补零。块内只用 FFT 与逐元素运算，所以结果与分块、进程数无关。

`uniform_exposure` 给出横向均匀曝光的一列结果（虚拟感光仪用）：`e_X,z = Σ_b N_b (W+ + f_H W-)`，用的是未截断的 `f_H`。均匀画幅内部与它之差，恰好是被截去的晕尾。

**曝光模型。** `bw-mechanistic` 的曝光（`models.bw.expose_bw`）：深度单元数 `n_z` 由离散化档位决定（`reference` 8、`calibration` 4、`image` 2，见 `develop/grid.py` 的 `PROFILE_LEVELS`），尺寸档数见 4.1 节；片基取 `RunContext.film_format`，缺省取材料的 `default_format`。

参数（字段在 `emulsion.optics.*` 与 `base.*`；类别是规格对将来胶片文件的规定，因为还没有这样的文件入库）：

| 符号 | 字段 | 单位 | 类别 | 说明 |
|------|------|------|------|------|
| n_gel | `emulsion.optics.gelatin_index` | 1 | assumed（未核实） | 缺省 1.54 |
| n_b | `base.refractive_index` | 1 | literature | 1.48（三醋酸纤维素，专利 US 6,398,371 B1） |
| t_b | `base.formats.<f>.thickness` | um | datasheet（T-Max）/ assumed（Tri-X 沿用旧乳剂） | |
| D_g | `base.formats.<f>.grey_density` | 1 | assumed | 对"135 与 120 的基雾差"的解释未核实 |
| D_back | `base.formats.<f>.backing_density` | 1 | assumed | 120 的防光晕背涂层 |
| α_0、E_0、E_U | `emulsion.optics.intrinsic_absorption.*` | um^-1、eV、eV | fitted | 待拟合 |
| A_j、ν_j、w_j | `emulsion.optics.dye_bands[j]` | 1、um^-1、um^-1 | fitted | 待拟合 |
| S_550 | `emulsion.optics.scattering_550` | um^-1 | fitted | 待拟合 |
| q_S | `emulsion.optics.scattering_exponent` | 1 | assumed | 缺省 1.0 |
| w_P、a_P | `emulsion.optics.psf_weight`、`psf_length` | 1、um | fitted | 将对手册 MTF 标定 |

另有常数：`PHOTON_ENERGY_EV_NM = 1239.84`；晕核超采样 8/32、截断份额 1e-3、半径上限 3000 µm；PSF 超采样 16、尾部质量 `PSF_TAIL = 1e-4`；块 `BLOCK_PX = 512`。这些是数值参数，不是物理量。

限制情形与测试：

| 测试 | 检验什么 |
|------|----------|
| `test_km_energy_balance`、`test_km_energy_balance_vanishing_absorption` | 能量守恒，相对误差 ≤ 1e-12（含吸收趋于 0） |
| `test_km_beer_lambert_limit` | `S → 0`、无界面反射时，吸收分数 = `1 - exp(-K L_e)` |
| `test_km_thick_layer_reflectance` | 极厚层的漫反射 → `R_inf` |
| `test_depth_response_balance_and_cells` | `absorbed + escaped_top + into_base = T_top`，且各深度单元可加 |
| `test_interface_reflectances_numbers`、`test_urbach_edge_is_continuous` | 界面反射率；Urbach 尾在 `E_0` 处连续 |
| `test_scatter_kernel_normalised` | PSF 核之和为 1 |
| `test_scatter_kernel_otf` | 离散核的传递函数与像素积分后的连续 OTF 之差 ≤ 2e-3 |
| `test_scatter_kernel_centre_pixel_disc`、`test_tail_radius_matches_spec` | 中心像素的圆盘公式；截断半径 |
| `test_halation_fraction_clear_base`、`test_halation_grey_base_numbers` | 透明片基的 `f_H` 与全反射份额 `1 - 1/n_b²`；灰片基 `f_H = 0.222`、`r_c = 234 µm` |
| `test_halation_kernel_truncation`、`test_halation_kernel_radius_cap`、`test_halation_otf_matches_kernel` | 核的截断与传递函数 |
| `test_low_rank_projection_matches_full` | 低秩路径与展开为全谱的路径相对差 ≤ 1e-12 |
| `test_expose_linear_in_light` | 光加倍则场加倍 |
| `test_expose_uniform_interior`、`test_uniform_frame_equals_uniform_exposure` | 均匀光时远离边缘处是常数，且等于一列的结果 |
| `test_photon_conservation` | 画幅内一块亮区，吸收的光子数守恒 |
| `test_halation_ring_position`、`test_edge_spread` | 点光源的光晕环落在 `2 t_b tan θ_c`；刃边响应等于两个核边缘分布的卷积 |
| `test_block_convolution_equals_global_fft`、`test_block_result_independent_of_alignment` | 分块卷积与整幅 FFT 一致 |
| `test_expose_independent_of_tiling_and_workers` | `tile_px` 与进程数不改变结果（逐位） |
| `test_expose_bw_fields_and_meta`、`test_profile_sets_depth_cells`、`test_film_format_selects_base` | 写出的字段与 `meta`；档位决定 `n_z`；画幅决定片基 |

`test_scatter_kernel_otf_vs_continuous`（核与未经像素积分的连续 `P^(nu)` 比较）被标为预期失败：像素孔径与混叠本身就带来约 1e-2 的差，规格字面的 2e-3 达不到，所以 `test_scatter_kernel_otf` 改为与像素积分后的连续 OTF 比较。

`EmulsionState` 读写：

| 读 | 写 |
|----|----|
| 光（`LightField`，见 4.5）、材料的乳剂层与片基 | `emulsion/exposure_volume`，`(H, W, n_z)`，um^-3 |
| | `emulsion/exposure_surface`，`(H, W, n_z)`，um^-2 |
| | `emulsion/residual_dye`，值 1.0 的广播视图 |
| | `meta["exposure"]`：`time_s`、`temperature_c`（23.0）、`n_z`、`profile`、`film_format`、`halation_fraction`、`halation_radius_um`、`halation_truncated`、`psf_radius_um`、`block_px`、`halo_px` |

`pipeline.expose` 另把 `light.meta` 存进 `meta["light"]`。23 °C 是 ISO 6 §5.3.1 的曝光环境，是固定值，不可调。

不物理之处：乳剂内按两流近似处理，准直入射光一进乳剂就视为漫射；PSF 与深度、波长无关；光晕只算一次返回；灰片基与背涂层中性；明胶折射率未核实；Tri-X 的光谱、MTF 数据来自旧乳剂；吸收的绝对水平与下一节的潜影效率 `η_0` 只以乘积进入模型，需要靠层吸收率的约束和深度效应区分（可辨识性问题，标定时处理）。

### 4.3 潜影

代码：`latent/chain.py`（`ChainParams`、`generator`、`cluster_distribution`）、`latent/table.py`（`latent_table`、`LatentTable`）。

这一步做什么：给定每颗粒平均吸收的光子数 `m`（4.2 节的 `m_c,z`）和曝光时间 `t_exp`，算每颗粒潜影团簇大小 `n` 的分布 `p(n | m)`。泊松吸收、复合损失、单原子衰减（低照度失效）和离子步骤阻塞（高照度失效）都在这个连续时间马尔可夫链里。

**链。** 每个颗粒跟踪两个位点：主位点 1 和一个汇总的竞争位点 2。状态 `(k1, k2)`，`k_s ∈ {0, …, n_cap}`（`n_cap = 5`，大于等于 5 并为一态，共 36 个状态）。电子到达潜影位点的速率

```
nu  = eta_c * m / t_exp ,    eta_c = eta_0 * (d_c / 1 um)^xi_eta
q1  = 1 / (1 + nu tau_i)                       到达主位点（离子步骤占用期间改去别处，M/M/1/1 阻塞）
q2  = (1 - q1) / (S_eff - 1)                   到达被跟踪的竞争位点
q_lost = (1 - q1)(S_eff - 2)/(S_eff - 1)       分散到其余位点，不改变状态

k_s -> k_s + 1    速率 nu q_s          k_s < n_cap
k_s = 1 -> 0      速率 1/tau_0 + rho_1 nu      单原子热衰减与空穴攻击
P(t_exp) = expm(Q t_exp) e_(0,0)
```

`S_eff = 2` 时没有丢失，退化为纯两位点模型。两个位点在给定 `nu` 与 `tau_0` 时互相独立，`Q` 是两个单位点生成元的 Kronecker 和，所以 36 状态的矩阵指数精确地等于两个单位点解的外积（`test_matches_36_state_expm`）。每个单位点的 `expm` 用缩放加泰勒级数再平方的方式计算，全部项非负，对很刚性的衰减速率也不会出现相消。

**陷阱深度分布。** 单原子寿命按颗粒固定，跨颗粒服从对数正态（Katz 要求准连续的陷阱深度分布）：

```
tau_0 = tau_0,med * exp(sigma_tau x),   x ~ N(0, 1)
p(n | m) = sum_q omega_q * p(n | m; tau_0(x_q))          概率论意义的 Gauss–Hermite 节点，sum omega = 1
```

`sigma_tau = 0` 或 `tau_0,med = ∞` 时只用一个节点。

**输出。** 曝光结束时每个颗粒的潜影大小取 `n = max(k1, k2)`。显影前的存放时间远长于 `tau_0`，单原子全部衰减，所以 `n = 1` 并入 `n = 0`（稳定团簇的退化不建模）。`cluster_distribution` 返回长 `n_cap + 1` 的向量，第 1 项恒为 0，其余各项之和为 1。可显影性不在这里判定：显影用整个分布算引发风险，所以不同显影液对"最小可显影团簇"的不同要求由显影动力学产生。

**查找表。** 每帧的 `t_exp` 固定，每个尺寸档建一张表：`log10 m` 在 [−6, 4] 上等距 201 个节点，存 `p(n | m)`。查表时在 `log10 m` 上用三次样条（not-a-knot）插值，表外取端点，负值截到 0 后重新归一；`m = 0` 取第一个节点。表只依赖胶片参数与 `t_exp`，与化学无关。表的缓存由调用者负责。

**倒易律失效如何产生。**

- 低照度失效：`nu tau_0` 与 1 可比时，第一个原子在第二个电子到达前衰减。单一 `tau_0` 时每十倍时间约需 1.3 档，比手册（约 1 档/十倍）陡；`sigma_tau` 取 1 至 3.5 把它压平。失效对较高密度影响更大，反差上升。
- 高照度失效：`nu tau_i ≳ 1` 时电子分散到多个位点，形成多个亚可显影团簇；较高密度需要更多补偿，反差下降。`nu tau_i → ∞` 时各转移率趋于常数比例，失效饱和在 `m` 的 `S_eff − 1` 倍（`S_eff = 3` 时恰为 1 档），所以短于约 1e-5 s 的曝光不再继续失效。这个 `q1` 是平均场的估计式闭合，不是由位点状态推出的。

参数（`ChainParams`，字段在 `emulsion.latent.*`）：

| 符号 | 字段 | 单位 | 类别 | 说明 |
|------|------|------|------|------|
| η_0 | `efficiency` | 1 | fitted | 待拟合 |
| ξ_η | `efficiency_size_exponent` | 1 | assumed（缺省 0） | 档的效率随直径的指数 |
| τ_0,med | `tau0_median` | s | fitted | 允许 ∞（无衰减） |
| σ_τ | `tau0_sigma_ln` | 1 | fitted | 对数标准差 |
| τ_i | `ionic_time` | s | fitted | 微秒量级，来自回忆，不是测量 |
| S_eff | `competing_sites` | 1 | assumed（缺省 3） | ≥ 2 |
| ρ_1 | `hole_attack` | 1 | assumed（缺省 0） | |
| n_cap | `cluster_cap` | 1 | assumed（5） | ≥ 2 |
| T_exp | — | C | assumed（23） | 曝光温度，见 `meta["exposure"]` |

数值参数：表节点数 201，Gauss–Hermite 节点数 `quadrature` 缺省 9。

限制情形与测试：

| 测试 | 检验什么 |
|------|----------|
| `test_poisson_r_hit_limit` | `tau_0 = ∞`、`tau_i = 0`、`S_eff = 2`、`rho_1 = 0` 时，`P(n ≥ r)` 等于泊松 r 次命中公式，误差 ≤ 1e-10 |
| `test_probability_conservation` | 概率之和为 1 |
| `test_monotone_in_exposure` | `P(n ≥ r)` 随 `m` 不减 |
| `test_no_reciprocity_failure_without_decay_or_blocking` | 无衰减、无阻塞时分布与 `t_exp` 无关 |
| `test_single_site_nucleation_formula` | 三态链的双指数闭式解 |
| `test_lirf_raises_contrast`、`test_trap_depth_spread_flattens_lirf` | 长曝光使较高密度的补偿增大；`sigma_tau` 压平低照度失效 |
| `test_hirf_asymptote` | `S_eff = 3`、`nu tau_i → 1e4` 时失效为 1 档 ± 0.02 |
| `test_hirf_lowers_contrast` | 短曝光使较高密度的补偿增大 |
| `test_efficiency_size_exponent` | 大颗粒的效率按 `xi_eta` 缩放 |
| `test_gauss_hermite_converged_small_spread` | `sigma_tau = 0.5` 时 9 节点与 21 节点之差 ≤ 1e-5 |
| `test_latent_table_interpolation` | 表中点的插值与直接计算之差 ≤ 1e-4 |
| `test_latent_table_deterministic`、`test_latent_table_build_time` | 逐位可复现；K = 6、201 节点、9 个求积点在一个核上不到 1 s |

没有达到规格要求的一项：`test_gauss_hermite_converged` 要求 `sigma_tau ≤ 4` 时 9 节点与 21 节点之差 ≤ 1e-5，实测达不到（`sigma_tau` 为 1、2、3、4 时分别约 1.4e-4、5.9e-3、2.0e-2、3.5e-2；21 节点在 4 时本身离收敛也有约 9.4e-3），该测试被标为预期失败。原因是衰减的转折在 `x` 上只有约 `1/sigma_tau` 宽。后果：`sigma_tau` 较大时（标定里 `sigma_tau` 的先验上界到 4），9 个求积点给出的 `p(n | m)` 可能有百分之几的误差，这是个有待处理的数值缺陷。

`EmulsionState` 读写：这一步不读也不写状态字段。潜影分布不存为字段，它是曝光字段与 `t_exp` 的确定函数，显影步骤从 `emulsion/exposure_volume`、`emulsion/exposure_surface` 和 `meta["exposure"]["time_s"]` 现算（存下来 24 MP 需要约 46 GB）。

不物理之处：竞争位点的闭合是平均场式的估计；`sigma_tau` 与 `tau_i` 的物理值没有一手测量；显影前团簇不衰退；曝光温度固定；效率 `eta_c` 把电子–空穴复合、陷阱捕获等多个步骤合成一个数。

### 4.4 药液化学

代码：`chemistry/species.py`（`component_table`、`AcidBaseFamily`）、`chemistry/activity.py`（`conditional_pka` 等）、`chemistry/speciation.py`（`make_speciation`、`solve_ph`）、`chemistry/solution.py`（`working_solution`）。

这一步做什么：由药液文件的配方（每种配料的 g/L）、稀释和温度，算出工作液里各组分的总量、质子过量、离子强度和 pH。显影方程需要这些量，但它们本身与显影无关，所以单独成块。

**配料与组分。** `species.yaml` 有两类条目：配料（ingredient，配方里称量的东西，如 `sodium-sulfite`）和组分（component，溶液里输运的总量，如 `sulfite`＝SO₃²⁻＋HSO₃⁻＋H₂SO₃）。每种配料给出溶解后对各组分的贡献 `s_i,k`（mol/mol）、对质子总量的贡献 `tau_i`、Na⁺ 与 K⁺ 的个数和旁观离子电荷。每个组分带酸碱族（各质子化态的电荷与相对参考态的质子数、pK、温度律）、水中扩散系数 `D_w`（20 °C）、凝胶与水的扩散比 `r_gel`、活度类别和显影剂角色。`component_table` 另外追加两个内置伪组分：`toth`（质子过量，没有自己的扩散系数，通量由可动形态合成）和 `gelatin`（凝胶中不动的缓冲剂，G–NH₃⁺/G–NH₂，pK 取 `kinetics.gelatin.pka`，缺省 10.5）。Na⁺、K⁺、硫酸根等只进入离子强度，不输运。

质子过量：

```
TOTH = sum_f sum_j n_f,j [f, j] + [H+] - [OH-]
```

`n_f,j` 是族 `f` 的第 `j` 态相对参考态多出的质子数（硼砂 `B(OH)3` 为 0、`B(OH)4-` 为 −1，`HSO3-` 为 +1，等）。配料对 TOTH 的贡献等于它溶解后各态的 `n` 之和，例如 `Na2S2O5` → 2 `HSO3-`，`tau = +2`；硼砂 → 2 `B(OH)3` + 2 `B(OH)4-`，`tau = −2`；米吐尔硫酸盐 → 2 `MetH2+`，`tau = +2`。

**工作液。** 药液文件的基液（D-76、XTOL 为原液，HC-110 等为浓缩液）每种配料有 `w_i`（g/L）。稀释 `a+b` 的稀释分数 `f = a/(a+b)`，则

```
T_k   = f * sum_i (w_i / M_i) * s_i,k        [mol/L]，组分 k 的总量
TOTH  = f * sum_i (w_i / M_i) * tau_i        [mol/L]
m_cat = f * sum_i (w_i / M_i) * (n_Na,i + n_K,i)
```

`M_i` 是摩尔质量。直接称量的组分（如 `hydroquinone`）算作 `s = {自身: 1}`、`tau = 0`；它的参考态必须不带电，否则报错。`bromide`、`iodide`、`silver`、`sulfite` 在每个工作液里都有（总量可为 0），因为显影的源项与溶剂作用按名字读它们。`overrides` 可以替换某种配料的 g/L（干扰参数，如成品 D-76 的 B₂O₃）。没有酸碱族的组分的电荷，由含它的配料的电中性推出（`component_charges`）。

**活度。** 离子强度约 2 M 时 Davies 式失效，采用 SIT 形式，pH 取氢离子活度（玻璃电极惯例）：

```
log10 gamma_i = -z_i^2 * D(I) + eps_i * m_cat
D(I) = A(T) * sqrt(I) / (1 + b sqrt(I)),    b = 1.5
A(T)：在 0.5050（20 °C）与 0.5091（25 °C）之间对 T 线性插值、外推
```

对酸碱族第 `j` 步 `HA^z <=> A^(z-1) + H+`，混合（条件）常数把浓度比与氢离子活度联系起来：

```
pK'_f,j(T, I) = pK_f,j(T) + (2 z_HA - 1) * D(I) + d_eps_class * m_cat
pH = pK'_f,j + log10( [A] / [HA] )
[H+] = a_H / gamma_H ,    [OH-] = K_w / (a_H gamma_OH),     log10 gamma_H = log10 gamma_OH = -D(I)
```

`(2 z_HA − 1)` 来自 `log γ_A − log γ_HA`：`HSO3-`/`SO3^2-`（`z_HA = −1`）为 −3D，`B(OH)3`/`B(OH)4-` 为 −D，质子化胺为 +D。`d_eps_class` 按活度类别分四类：`sulfite`、`borate`、`organic`、`amine`，单位 kg/mol；`m_cat` 的单位在代码里按 mol/L 计，忽略了摩尔浓度与质量摩尔浓度的差别。

温度律（`AcidBaseFamily.pka`）：

```
van 't Hoff：pK(T) = pK(T_r) + dH / (R ln 10) * (1/T - 1/T_r),   T_r = 298.15 K
硼酸（Manov 等 1944）：pK(T) = 2237.94/T + 0.016883 T - 3.305
胺：pK(T) = pK(25 °C) + (dpK/dT) (T - 298.15)
pK_w：过 14.17（20 °C）与 14.00（25 °C）的 van 't Hoff 直线，dH_w ≈ 56.9 kJ/mol
```

有机显影剂的电离焓取 0：它们的温度效应并入显影的表观活化能，以免重复计算。

**形态分布与 pH。** 给定各族总量 `T_f` 与 `TOTH`，解一个标量方程：

```
F(pH) = sum_f T_f * nbar_f(pH) + [H+] - [OH-] - TOTH = 0
nbar_f(pH) = sum_j n_f,j * alpha_f,j(pH),    alpha 由 pK' 给出
beta = -dF/dpH = ln10 * ( sum_f T_f var_f(pH) + [H+] + [OH-] )     缓冲容量，mol/(L pH)
```

`var_f` 是族 `f` 的质子数在各态上的方差。`F` 对 pH 单调递减。求根用带区间保护的牛顿法（pH ∈ [0, 14.5]；牛顿步出界或收敛不够快时改为区间二分），容差 `|ΔpH| ≤ 1e-12`，最多 30 次，对所有单元向量化，每个单元收敛后即冻结，所以结果不依赖同批的其他单元。根在区间外时报错。

**离子强度。** `I = ½ Σ z_s² [s]`，求和包含阳离子、旁观离子、每个酸碱态、没有酸碱族的带电组分和水的离子。从 `I⁰ = m_cat` 起做定点迭代，直到 `|ΔI| < 1e-6 M`（IONIC_STRENGTH_TOL），最多 100 次，不收敛则报错。显影过程中 `I` 视为常数。

参数：

| 符号 | 来源 | 类别 | 说明 |
|------|------|------|------|
| A(T) | `kinetics.activity.debye_huckel_A_20C`、`_25C` | assumed | 0.5050、0.5091，回忆的标准值，待核对 |
| b | `kinetics.activity.sit_b` | assumed | 1.5，SIT 惯例 |
| Δε（四类） | `kinetics.activity.delta_eps.*` | 先验 N(0, 0.15²)，assumed | **本文写作时是占位值 0**，等标定阶段 A 用约 11 个 pH 锚点拟合 |
| pK_w | `activity.PKW_20C`、`PKW_25C` | assumed | 14.17、14.00，回忆值 |
| 各族 pK | `species.yaml` | 见下 | |
| D_w、r_gel | `species.yaml` | 见下 | 本节只存储，输运由显影步骤使用 |

`species.yaml` 里组分常数的来源（只列类别，具体出处与数值见该文件）：亚硫酸盐的 pK₁、pK₂ 与 ΔH 是 `literature`（Goldberg & Parker 1985）；硼酸 pK 是 `literature`（Manov 等 1944）；氢醌、对氨基苯酚、二乙醇胺、乙醇胺的 pK 是 `literature`；米吐尔的酚 pK、抗坏血酸的 pK、DTPA 与 HEDTA 的 pK 是 `assumed`；Dimezone S、菲尼酮（phenidone）、甲基菲尼酮共用的吡唑烷酮 pK 是 `assumed` 的占位值 8.5（若 `kinetics` 文件给出该类的 pK，则以它覆盖）；所有水中扩散系数 `D_w` 是 `assumed`（由回忆的离子扩散系数换算或估计）；凝胶扩散比 `r_gel` 里氢醌与米吐尔是 `literature`（Iwano），其余是 `assumed`（盐类取 0.5）。

限制情形与测试：

| 测试 | 检验什么 |
|------|----------|
| `test_speciation_ideal_limit`、`test_speciation_ideal_limit_textbook_form` | `I → 0`、Δε = 0 时，单一弱酸–共轭碱混合物的 pH 与 Henderson–Hasselbalch 解析式之差 ≤ 1e-10 |
| `test_buffer_titration` | TOTH 增加后 pH 的变化与 `∫ dTOTH/β` 的数值积分一致 |
| `test_toth_round_trip`、`test_charge_balance_residual`、`test_mass_balance` | 由形态分布重算的 TOTH 与输入一致；电荷与质量守恒 |
| `test_buffer_capacity_matches_numerical_derivative` | `β` 与 `−dF/dpH` 的数值导数一致 |
| `test_boric_acid_pka_temperature` | 20 °C 9.278 ± 0.001，25 °C 9.235 ± 0.001 |
| `test_van_t_hoff_sulfite`、`test_linear_law_amines`、`test_organic_agents_have_no_temperature_dependence` | 各温度律 |
| `test_d76_published_ideal_ph` | 理想溶液（`D = 0`、Δε = 0）D-76 公开配方 pH = 8.95 ± 0.02（手册估计 8.94） |
| `test_xtol_ideal_ph` | 理想溶液 XTOL 工作液 8.89 ± 0.02，A 部分（850 mL）11.04 ± 0.03 |
| `test_sit_debye_huckel_only` | 只保留 `D(I)` 项（Δε = 0）：D-76 公开配方 8.44 ± 0.02（实测约 8.5），XTOL 工作液 8.29 ± 0.02。这是实现的回归测试，不是验证 |
| `test_ionic_strength_d76` | D-76 公开配方 `I = 2.39 ± 0.03` M |
| `test_ionic_strength_is_a_fixed_point` | 返回的 `I` 是定点 |
| `test_dilution_raises_ph_in_sulfite_buffer` | 含 `D(I)` 时 XTOL 1+1 的 pH 高于原液 |
| `test_working_solution_conserves_mass`、`test_toth_and_cations` | 各配料换算后元素总量守恒；TOTH 与阳离子数的例子 |

XTOL 的理想溶液与 `D(I)` 数值是在 Dimezone S 不电离的读法下算的（手册当时没有它的 pK）；`species.yaml` 的占位 pK 8.5 会使 pH 更低（`test_placeholder_dimezone_pka_lowers_xtol_ph`），测试里把它的用量覆盖为 0 来对应那个读法。

`EmulsionState` 读写：这一步不读也不写状态；输入是 `Bath`、`Species` 表和 `KineticsFile`，输出是不可变的 `WorkingSolution`。

不物理之处：SIT 形式的 `A(T)`、`b`、pK_w 是回忆值；Δε 只由约 11 个 pH 锚点标定，缓冲容量在 D-76 的 pH 下对活度很敏感（设计期估计总 β 从 92 降到 33 mM/pH，其中亚硫酸盐从 81 降到 17），所以显影中的 pH 降直接受它支配；成品 D-76 的 pH 与 B₂O₃ 的用量之间有张力（只保留 `D(I)` 项时算得 8.17 至 8.32，低于测量的 8.5）；配方本身有不确定（HC-110 旧配方的数字为二手，T-MAX 的配料读法是假设）；多硼酸根忽略；离子强度在显影中不变；每个组分一个扩散系数，不用 Nernst–Planck，凝胶–液体分配系数取 1；空气氧化和补充液不建模。

### 4.5 内存设计

代码：`light.py`（`LowRankLight`）、`fieldstore.py`（`FieldStore`）、`blocks.py`、`state.py`（`FrameAssembler`、`EmulsionState.drop`）、`pipeline.py`（`segments`、`process`）。

问题：24 MP（6000×4000）下，`(H, W, 81)` 的 float64 光通量就有 15.6 GB，每拼一块就复制整个字段的做法随块数和画幅平方增长。下面几处改动把常驻内存和总代价压下来；它们不改变任何物理量。

**低秩光。** RGB 场景经 Mallett 基还原是线性的：

```
Phi_b(x) = sum_r C_r(x) * Q_r,b ,     r = 1 … R    （RGB 场景 R = 3）
```

`LowRankLight` 只存系数 `coefficients`（`(H, W, R)`，可与输入 RGB 数组共享）和基 `basis`（`(R, B)`，每单位系数每波段的光子通量），从不生成 `(H, W, B)`。`project(weights)` 返回 `(H, W, M)`，计算为 `coefficients @ (basis @ weights)`；`window(rows, cols)` 只在需要时生成一个块的全谱 `Light`；`filtered(transmittance)` 只改基。`Light` 与 `LowRankLight` 都满足协议 `LightField`，`pipeline.expose` 接受任一种；`scene_light(..., low_rank=True)` 返回前者。检验：`test_low_rank_light_window_and_project`、`test_low_rank_light_memory`、`test_low_rank_light_filtered`、`test_expose_accepts_low_rank_light`。

**字段存储。** `FieldStore(root, ram_threshold_bytes=256 MiB)`：大小超过阈值（严格大于）且 `root` 不为 None 的字段分配为 `np.memmap(mode="w+")`；`freeze` 把它刷盘并以 `mode="r"` 重新打开，交给 `EmulsionState`；其余字段留在内存。映射文件名是 `<序号>-<字段名>.dat`，存储内的名字必须唯一，`release` 删除文件。只读的 `np.memmap` 在 `_frozen.readonly_array` 里被直接共享而不复制（可写映射仍被复制，因为别处可能改写它）。均匀字段（`residual_dye`）用 1×1 数组的 `np.broadcast_to` 视图，零字节。检验：`test_fieldstore_memmap_readonly_shared`、`test_readonly_array_copies_writable_maps`、`test_fieldstore_threshold_is_strict`、`test_fieldstore_release_and_close`。

**拼接器取代重复 paste。** `FrameAssembler`（`FieldSpec(key, trailing, dtype)`）为每个输出字段一次性分配（经 `FieldStore`），`write(key, inner, values)` 把块直接写进内块区域，`finish()` 检查每个像素恰好被写过一次（重叠的写入当场报错），冻结并返回只读字段；所有块写入同一个广播值的字段仍是广播视图，不分配内存。总代价 O(像素数)。`EmulsionState.paste` 保留为公开 API，每次调用复制整个字段，流水线不再用它；`EmulsionState.drop(*keys)` 删除字段。检验：`test_assembler_each_pixel_once`、`test_assembler_uniform_fields_stay_broadcast`、`test_assembler_uses_store_memmap`、`test_assembler_linear_cost`（按时间比较，机器负载高时可能偶发失败）、`test_drop_fields`。

**规范块。** `blocks.block_grid` 给出固定在画幅原点的正方形块网格：网格只由画幅形状、块大小和晕边决定，与 `RunContext.tile_px`、进程数、调度无关；每个块（包括边缘块）的补边后形状相同，使 FFT 尺寸和块内归约的形状固定；`padded_window` 在画幅外补零。`run_blocks(fn, blocks, workers)` 按块的顺序返回结果，`workers > 1` 时用新生成（spawn）的进程，且 BLAS 设为单线程，所以跨块的求和在任意进程数下逐位相同。检验：`test_block_grid_tiles_frame_once`、`test_block_grid_anchored_at_origin`、`test_run_blocks_order_and_identity_across_workers`。

**分段融合。** `pipeline.process` 把相邻的、作用范围有限的步骤（`interaction_range_um` 为有限数）合成一段：对每个块，从输入状态裁出窗口，依次执行段内所有步骤，只把段的最终字段写进拼接器；段的晕边是各步作用范围之和。作用范围为 `None`（全幅耦合，如显影）的步骤单独成段，拿到整幅状态。`pipeline.segments(proc, model)` 返回分段结果。检验：`test_segments`、`test_segment_fusion_equals_unfused`、`test_segment_halo_is_sum_of_ranges`、`test_fused_segment_writes_to_scratch_dir`。

**曝光的内存。** 曝光字段 `exposure_volume`、`exposure_surface`（各 `(H, W, n_z)`）经 `FieldStore` 映射到 `RunContext.scratch_dir` 下。多进程时，卷积源经临时 `.npy` 文件传给工作进程，各块结果返回主进程再写入拼接器。检验：`test_expose_workers_with_scratch_dir`、`test_expose_bw_small_fields_stay_in_ram`；慢速测试 `test_expose_24mp_within_memory_budget` 用 24 MP 的低秩光检查曝光阶段的峰值内存（追踪到的内存加映射文件）在 8 GB 以内。

设计估算（尚未实测）：整个黑白链路在 24 MP 下常驻内存约 5 至 6 GB，其中显影输出约 4.0 GB（每像素 21 个 float64）以内存映射存放；显影还没有实现，这一估算不能当作已经达到的指标。

状态：`EmulsionState` 的 `fields` 是只读数组（可为内存映射或广播视图），不改变语义。`RunContext` 新增的字段 `film_format`、`profile`、`workers`、`memory_budget_gb`、`scratch_dir` 不进入配方哈希。运行器对临时目录的创建与清理、命令行选项的接线还没有完成，这里不写。

### 4.6 各阶段读写的状态字段

| 阶段 | 读 | 写 |
|------|----|----|
| 颗粒分档（4.1） | 乳剂层参数 | 无（返回 `SizeClasses`） |
| 曝光（4.2） | `LightField`、材料 | `emulsion/exposure_volume`、`emulsion/exposure_surface`、`emulsion/residual_dye`、`meta["exposure"]` |
| 潜影（4.3） | 无（显影时由曝光字段与 `meta["exposure"]["time_s"]` 现算） | 无（返回 `LatentTable`） |
| 药液化学（4.4） | 药液、物种表、动力学文件 | 无（返回 `WorkingSolution`） |

## 5 已知的不物理之处

- **RGB → 光谱的还原**。一组 RGB 值对应无穷多种光谱（同色异谱），原始场景是哪一种已经无从得知；基光谱法为每种颜色挑了一条平滑的光谱。胶片的光谱敏感度与人眼不同，所以真实场景里两种看上去一样、光谱不同的颜色，在胶片上可能不一样，这一点模拟不出来。直接输入光谱图像可以绕开这一步。
- **8 位输入**。8 位 sRGB/JPEG 已经被裁切和压缩到六七档亮度，负片的宽容度和高光表现体现不出来。
- **超出 1 的场景值**。场景线性数据里高于 1 的值按"反射率大于 1"处理，相当于光源或镜面高光，光谱形状仍由基光谱决定。
- **乳剂配方是商业机密**（D1 起）。胶片参数是对着手册曲线标定出来的。
- **部分专有显影液的配方是推算值**（D1 起）。
- **显影动力学常数大部分是标定值**（D1 起）。模型的结构是物理的，常数的可信度靠"保留组合"的预测误差来衡量：这些组合从不参与标定，只用来检验。

D1 已实现模块里的不物理之处（章节号指第 4 节）：

- **所有 `fitted` 参数现在都没有被拟合**（4.1 至 4.3）。颗粒分布、吸收、散射、PSF、潜影的参数只有先验区间；`data/films/kodak/` 里的 `bw-mechanistic` 文件只带 `assumed` 的起始值。D1 的测试用的是 `tests/fixtures/` 里的值，检验方程与代码一致，不检验与任何胶片的符合程度。
- **涂布银量没有一手数据**（4.1）。它不只与遮盖力以乘积进入最大密度，还决定每档的银浓度、单位密度释放的溴离子与消耗的显影剂；标定时若不可辨识，它的误差会被别的常数吸收。
- **颗粒群体的简化**（4.1）：对数正态加 ±3σ 截断；颗粒在深度上均匀分布；每档一个代表直径；不区分立方与球。
- **乳剂光学**（4.2）：两流近似，准直光入射后立即视为漫射；PSF 与深度、波长无关；光晕只算一次返回，片基灰染和背涂层视为中性，135 的灰染密度 0.10 与 120 的背涂层密度 0.5 是假设（对"135 与 120 的基雾差"的解释未核实）；明胶折射率未核实；Tri-X 的光谱、MTF 与颗粒度数据来自旧乳剂；吸收的绝对水平与潜影效率只以乘积进入模型，要靠层吸收率的约束与深度效应才能区分。核按一个微小的截断尾量截短（晕核 0.1 %，PSF 0.01 %）。
- **潜影链**（4.3）：高照度失效来自 `q1 = 1/(1 + nu tau_i)` 这个平均场式的估计闭合，不是由位点状态推出的，`nu tau_i → ∞` 时失效饱和在 `m` 的 `S_eff − 1` 倍，短于约 1e-5 s 的曝光不再继续失效；`tau_i`、`sigma_tau` 没有一手测量；显影前团簇不衰退；曝光温度固定为 23 °C；效率把多个微观步骤合成一个数；`sigma_tau` 较大时 9 个 Gauss–Hermite 节点不够收敛（`test_gauss_hermite_converged` 预期失败，见 4.3）。
- **药液化学**（4.4）：SIT 的 `A(T)`、`b`、pK_w 是回忆的标准值；四个 Δε 现在是占位值 0，将来只由约 11 个 pH 锚点拟合，缓冲容量对活度很敏感；成品 D-76 的 pH 与 B₂O₃ 用量之间有张力；吡唑烷酮的 pK 是占位值；`D_w` 全部是 `assumed`；每个组分一个扩散系数，不用 Nernst–Planck，不计明胶电荷对阴离子的 Donnan 排斥；离子强度在显影中不变；多硼酸根、空气氧化、补充液忽略。
- **内存设计**（4.5）：24 MP 下整个黑白链路 5 至 6 GB 的预算是设计估算，只有曝光阶段有测试；性能测试 `test_assembler_linear_cost` 按时间判断，机器负载高时可能偶发失败。

以上各项在数据文件的来源标签和验证报告里都会体现。

## 6 参考文献

- T. H. James (ed.), *The Theory of the Photographic Process*, 4th ed., Macmillan, 1977.（潜影、显影动力学）
- J. C. Dainty, R. Shaw, *Image Science*, Academic Press, 1974.（颗粒度、Nutting 公式、MTF）
- G. Haist, *Modern Photographic Processing*, Wiley, 1979.（显影液化学）
- S. Anchell, B. Troop, *The Film Developing Cookbook*.（公开配方）
- A. Newson, J. Delon, B. Galerne, "A Stochastic Film Grain Model for Resolution-Independent Rendering", *Computer Graphics Forum*, 2017.（颗粒渲染）
- Zhang et al., 2023.（颗粒渲染的解析加速；完整条目待 D1 实现显微镜模式时补全）
- I. Mallett, C. Yuksel, "Spectral Primary Decomposition for Rendering with sRGB Reflectance", EGSR 2019.（D0 使用的光谱还原）
- W. Jakob, J. Hanika, "A Low-Dimensional Function Space for Efficient Spectral Upsampling", *Computer Graphics Forum*, 2019.（光谱还原）
- F. C. Williams, F. R. Clapper, "Multiple Internal Reflections in Photographic Color Prints", *JOSA*, 1953.（相纸反射密度）
- ISO 6:1993，黑白胶片感光度；ISO 2721，照相机自动曝光控制。（测光常数；曝光环境 23 °C 取自 ISO 6）
- R. N. Wolfe, E. W. Marchand, J. J. DePalma, *JOSA* **58**(9), 1245–1256 (1968).（用 Kubelka–Munk 两流理论分析细颗粒溴化银层的体积反射与透射；二手来源，只读到摘要）
- H. Frieser，指数线扩散函数；P. Nowak, "Modelling of the modulation transfer function of silver halide light-sensitive layer including DIR and DAR dye couplers", *Optica Applicata* **28**(2), 95 (1998).（PSF 的函数形式，Nowak 一文核实过；Frieser 1960 年原文的引用只来自二手来源）
- US 6,398,371 B1（Fuji Photo Film）。（三醋酸纤维素片基折射率 1.48）
- E. Katz, "On the Photographic Reciprocity Law Failure and Related Effects. I. The Low Intensity Failure", *J. Chem. Phys.* **17**(11), 1132 (1949).（需要准连续的陷阱深度分布；只读到摘要）
- R. K. Hailstone, J. French, R. De Keyzer, "Latent-image formation in tabular AgBr grains: experimental studies", *Imaging Sci. J.* **52**(3), 151 (2004); R. Hailstone, R. De Keyzer, "Latent-image formation in tabular AgBr grains: simulation studies", *Imaging Sci. J.* **52**(3), 164 (2004).（边缘陷阱深度，量子感光度；只读到摘要）
- R. N. Goldberg, V. B. Parker, "Thermodynamics of solution of SO₂(g) in water and of aqueous sulfur dioxide solutions", *J. Res. NBS* **90**(5), 341 (1985).（亚硫酸的 pK 与焓）
- H. Iwano, "The Rate of Penetration of Photographic Chemicals into Gelatin Layers (2): The Diffusion-retarding Factors in Gelatin", *J. Soc. Photogr. Sci. Technol. Japan* **33**(3), 143 (1970).（氢醌、米吐尔在凝胶中的扩散比）
- G. G. Manov, N. J. DeLollis, S. F. Acree, "Ionization constant of boric acid and the pH of certain borax-chloride buffer solutions from 0° to 60° C", *J. Res. NBS* **33**, 287 (1944).（硼酸 pK 的温度式）
- IEC 61966-2-1，sRGB。
- CIE 1931 标准色度观察者；CIE 标准光源 D65、D50、A、E。
- Kodak、Ilford、Fujifilm 各胶片的技术数据手册（出版物编号以实际下载为准）。
