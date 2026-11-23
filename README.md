# The Darkroom Atlas

![version](https://img.shields.io/badge/version-D0c-blue)
![last commit](https://img.shields.io/github/last-commit/ChromiteCr/the_darkroom_atlas)
![commit activity](https://img.shields.io/github/commit-activity/m/ChromiteCr/the_darkroom_atlas)
![stars](https://img.shields.io/github/stars/ChromiteCr/the_darkroom_atlas)

在计算机里用物理方式模拟胶片的曝光、冲洗和成像，并能处理真实照片。胶片是乳剂数据，显影液是化学配方，两者分开建模；一种胶片配一种显影液会得到什么结果，由模拟算出来，而不是套一条预先画好的曲线。

## 现在能做什么

当前版本是 D0c：底座、Kodak 手册参考数据，以及 D1 黑白链路的大部分部件（曝光、潜影、药液化学、显影动力学、银影像与颗粒、停显定影）。**还不能冲出一张真实胶片的图**：把这些部件串起来的显影引擎、读出和命令行正在开发；两种胶片的参数是标定前的起始值，不代表真实的 Tri-X 或 T-Max。

仓库里唯一的材料是"理想材料"：它线性记录 CIE XYZ，没有颗粒、没有化学、没有损失。它不是物理模型，用途是检验光谱链路（RGB 还原为光谱、光源归一化、色匹配函数、曝光对波段求和、冲洗、读出）的约定是否自洽。一张图经过它之后应当原样返回，误差小于 0.005。测光和光子数换算这些绝对光度量在往返中会整体抵消，它们由单元测试单独锚定，见 `docs/physics.md` 第 3.4 节。

D0 已经定下来的东西：

- 单位与光谱约定：380–780 nm、步长 5 nm；胶片面长度用 µm；曝光量用每 µm² 每波段的光子数。
- 三个对象：感光材料（胶片和相纸同一结构）、药液、流程。
- 三个操作：曝光 `expose`、冲洗 `process`、读出 `readout`。新的物理模型通过注册表接入，不改管线。
- 带来源标签的数据格式，以及数据检查命令。
- 大图分块计算，分块与不分块的结果逐位相同。
- 所有随机数由种子和像素坐标决定，同一种子重跑的结果逐位相同。
- 每次运行附带一份运行记录，写明配方、数据文件哈希、版本和种子。

## 快速开始

需要 [uv](https://docs.astral.sh/uv/) 和 Python 3.13。

```bash
uv sync
```

```bash
uv run pytest
```

```bash
uv run atlas run tests/assets/sample.png --material ideal --process identity -o out/
```

```bash
uv run atlas check-data
```

`atlas run` 在输出目录里写两个文件：`sample__ideal__identity.png`（16 位）和同名的 `.run.json`（运行记录）。输入是 8 位图像时会给出提示，因为 8 位图像只保留了场景的六七档亮度，负片的宽容度体现不出来；以后的真实胶片应当喂场景线性数据（RAW、线性 TIFF）。

常用选项：`--iso`、`--ev`、`--time`（曝光时间，秒）、`--seed`、`--tile`（分块边长，像素）、`--long-edge-mm`（画幅长边，默认 36）、`--data`（另一个数据目录）。

## 数据

数据都是纯文本文件，放在 `data/`：

```
data/chemistry/species.yaml      化学物种
data/developers/<名字>.yaml      显影液配方
data/films/<厂商>/<名字>.yaml    胶片：层结构、颗粒、光谱敏感度
data/papers/                     相纸，与胶片同一结构
data/processes/<名字>.yaml       流程：有序的步骤
data/reference/<材料>/*.csv      从厂商手册数字化的参考曲线
```

每个物理参数都带来源标签，四选一：

| 标签 | 含义 |
|------|------|
| `datasheet` | 取自厂商手册，必须写明出处 |
| `literature` | 取自文献，必须写明出处 |
| `fitted` | 标定得到 |
| `assumed` | 假设值 |

显影液配方另带可信度：`published`（公开配方）、`patent`（专利）、`msds-estimated`（由安全数据表推算）。

曲线只从厂商手册自行数字化，不拷贝其他项目的数据，不使用有版权的汇编数据。没有出处的真实胶片或药液数据不进仓库。

提取工具在 `tools/digitize/`，直接读手册 PDF 里的矢量路径，每张图一份配方文件，写明出处、PDF 哈希和误差估计，可以重跑并逐字节核对。Tri-X 400 与 T-Max 100 的曲线和冲洗时间表在 `data/reference/kodak/`，显影液手册的表格在 `data/reference/developers/`。手册 PDF 本身不进仓库。说明见 [tools/digitize/README.md](tools/digitize/README.md)。

## 路线

| 阶段 | 内容 |
|------|------|
| D0 | 底座与约定（已完成） |
| D1 | 黑白物理链路：潜影、显影反应-扩散、银密度、颗粒；Tri-X 400、T-Max 100 与 D-76、XTOL、HC-110、Rodinal |
| D2 | 彩色物理链路：多层乳剂、成色、层间效应、色罩、反转；Portra 400、Ektachrome E100 |
| D3 | 暗房印相：放大机、相纸、相纸显影液、反射密度 |
| D4 | 图鉴扩充：批量添加胶片、显影液和流程 |
| D5 | 性能与界面：GPU 后端、交互界面 |

## 文档

- [docs/physics.md](docs/physics.md)：约定、模型方程、已知的不物理之处、参考文献。
- [设计文档](docs/superpowers/specs/2026-11-20-darkroom-atlas-design.md)：整体设计与验证方法。

## 版本记录

| 版本 | 日期 | 变更内容 | 类型 |
|------|------|----------|------|
| D0c | 2026-11-23 | 显影动力学与常数文件、显影数值内核、银影像与统计颗粒、停显/定影/水洗/干燥；显微镜模式的颗粒放置与渲染；Tri-X 400 与 T-Max 100 的胶片文件（标定前起始值）和 107 个流程文件；参考数据迁入 data/reference；训练/留出清单与参考观测；标定参数表；physics.md 新增 D1 章节 | feat |
| D0b | 2026-11-21 | D1 模型规格定稿；黑白链路前半段：D1 数据结构、低秩光与按块落盘的内存基础设施、层内光传输与光晕的曝光模型、潜影成核-生长链、带活度修正的酸碱平衡、七份显影液配方、感光度与反差指数指标 | feat |
| D0a | 2026-11-21 | 手册曲线提取工具；Tri-X 400、T-Max 100 的特性曲线、反差指数曲线、MTF、光谱敏感度与冲洗时间表，以及 D-76、XTOL、HC-110、T-MAX 显影液手册的表格（暂存，未接入模型）；D1 调研报告 | feat |
| D0 | 2026-11-20 | 完成底座：单位与光谱约定、带来源标签的数据模型、expose/process/readout 管线与注册表、分块与确定性随机数、运行记录、命令行；理想材料走通全流程 | milestone |
