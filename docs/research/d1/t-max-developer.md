# D1 T-MAX Developer 与 T-MAX RS Developer and Replenisher 的组成

研究日期：2026-11-21。范围：

* KODAK PROFESSIONAL T-MAX Developer（单液浓缩液，1+4 稀释）各版安全数据表（SDS/MSDS）列出的组分、CAS、浓缩液中的质量分数范围。
* 换算到工作液（1+4）的 g/L 与 mol/L，工作液 pH，初始溴化物，缓冲体系。
* 由 SDS 推算的配方有多大不确定度，并给出一个带上下界的中心配方，供 `data/developers/t-max.yaml` 使用（配方可信度 `msds-estimated`）。
* 顺带整理 T-MAX RS Developer and Replenisher（两液）的组成。

本报告沿用姊妹报告 `developer-formulas.md` 的格式和物种数据；凡是那份报告已经给出的 pKa、摩尔质量、CAS，这里只引用它的小节号，不重新推导。

**可信度标签**（与姊妹报告相同，标签本身保留英文）：
* **opened**：我自己打开并读了原文（PDF 用 PyMuPDF 抽取文本层后核对）。
* **secondary**：只来自搜索引擎摘要或片段，原文没有打开。
* **recalled**：凭记忆，本次没有打开原文；每条都写明要去哪里核实。
* **estimate**：我根据 opened 数据自己算的，方法写明。
* 送入 `data/` 时的来源标签按设计规格：`datasheet`、`literature`、`fitted`、`assumed`。

**关于托管站点。** 姊妹报告没有抓取经销商和大学托管的 SDS 副本，只把它们的搜索片段记为 secondary。本次任务没有给出站点白名单，而 T-MAX Developer 的 Kodak / Kodak Alaris SDS 在厂商站点上已经找不到（见下文"未取得"），所以我下载并读了这些第三方托管的副本，在表中逐一注明托管站点。它们在内容上是 Kodak、Kodak Alaris 或 Photo Systems Inc.（PSI）出具的文件（抬头、文档号、版本号齐全），但**我无法证明托管副本未被改动**。唯一来自 Kodak 自己站点的是 CIS-61 和 J-86（经 web.archive.org 存档或 kodakprofessional.com）。如果编排方要求严格执行姊妹报告的站点规则，可把下表中第三方托管的 SDS 一律降为 secondary；结论不变，因为 pH 锚点来自 Kodak 自己的 CIS-61。

**已打开并缓存的来源**（`.cache/datasheets/` 下，已被 git 忽略）：

| 缓存文件 | 文件内容 | 托管站点 |
|---|---|---|
| `tmax_sds/ek2007v1.1_1599844.pdf/.txt` | MSDS "KODAK T-MAX Developer"，Eastman Kodak Company，文档号 000000020762，Version 1.1，Revision Date 10/02/2007，Print Date 01/29/2009，product code 1599844 | web.archive.org 存档（2013-07-21 抓取）的 westliberty.edu 副本 |
| `tmax_sds/ka2014v1.8_1402767.pdf/.txt` | MSDS 同名，Kodak Alaris Inc.，000000020762，Version 1.8，Revision 07.02.2014，Print 08.03.2017，product code 1402767，synonym PCD 5337 | beauphoto.com（经销商） |
| `tmax_sds/ka2014v2.1_1402767.pdf/.txt` | SDS 同名，Kodak Alaris，000000020762，Version 2.1，Revision 02/07/2014，Print 05/27/2015，product code 1402767，PCD 5337 | macodirect.de（经销商） |
| `tmax_sds/ka2014v2.1_5050851.pdf/.txt` | 同上 Version 2.1，Print 07/28/2016，product code 5050851，PCD 5337；除 product code 和打印日期外与上一份逐行相同（diff 已查） | parallaxphotographic.coop（经销商） |
| `tmax_sds/ka2019v11_5050851.pdf/.txt` | SDS 同名，Kodak Alaris，SDS number 5050851，PCD 8021，Version 11，Revision 11-14-2019，Issue 08-15-2016 | archpdfs.lps.org（学区；文件名写作 "YR2000"，与内容不符） |
| `tmax_sds/ka2019_1058718.pdf/.txt` | SDS "KODAK T-MAX Developer"，Kodak Alaris，product 1058718，Version 02，Revision 06-20-2019，Issue 05-17-2019 | freestylephoto.com（经销商） |
| `tmax_sds/psi2024_1058718.pdf/.txt` | SDS "KODAK PROFESSIONAL T-Max Developer"，Photo Systems Inc.（Dexter, MI），product no. 1058718，SDS Version 2.0，2/22/2024，上一版 10/11/2023（1.0） | freestylephoto.com |
| `tmax_sds/agar2013_gb.pdf/.txt` | EU 格式 SDS "KODAK T-Max Developer"，由经销商 Agar Scientific Ltd 自行编制（ChemGes 软件），article AGP9167，Revision 06.08.2013，Version 1 | agarscientific.com |
| `tmax_sds/ka2014v2.1_rsAB_5054184.pdf/.txt` | SDS "KODAK T-MAX RS Developer and Replenisher, Part A"（Z17000000414，v2.1，02/07/2014，product 5054184，PCD 5492）与 "Part B"（PCD F1220）合订 | digitaltruth.com |
| `tmax_sds/ka2014v2.1_rsB_5054184.pdf/.txt` | SDS "… Part B"，Z33000000480，v2.1，02/07/2014，Print 12/14/2017 | freestylephoto.com |
| `tmax_kodak/cis61_1999.html` | Kodak CIS-61 "Specific Gravity and pH Measurements of KODAK Processing Chemicals"，October 1999（HTML 版） | web.archive.org 存档的 kodak.com（1999-12-01 抓取） |
| `tmax_kodak/cis61_2006.pdf/.txt` | CIS-61，February 2003 版 | web.archive.org 存档的 kodak.com（2006-08-13 抓取） |
| `tmax_kodak/cis61_2015.pdf/.txt` | CIS-61，May 2010 版，© Eastman Kodak 2010 | web.archive.org 存档的 kodak.com（2015-10-23 抓取） |
| `tmax_kodak/j86_2005.html` | J-86 "KODAK PROFESSIONAL T-MAX Developers"，Revised October 2002（HTML 版） | web.archive.org 存档的 kodak.com |
| `ds_txt/j86_2017.txt`（已有） | J-86，December 2017，13 页 | kodakprofessional.com（姊妹报告 `datasheets.md` 已打开） |
| `patents/US5756271A.txt`（已有） | US 5,756,271，Table I 与 Examples 2–5 | patents.google.com |
| `patents/US4987060A.html/.txt` | US 4,987,060（Marchesano，Minnesota Mining and Manufacturing），单液 HQ 浓缩显影液；**与 Kodak 无关**，只为排除而打开 | patents.google.com |
| `tmax_pubchem/*.json` | PubChem 记录：CAS 63149-47-3、21049-70-7、2218-94-2、7631-90-5、109-83-1 | pubchem.ncbi.nlm.nih.gov |
| `spotcheck_devf/iupac.csv`、`refs.csv`（已有） | IUPAC Digitized pKa Dataset v2.4，用于 2-(甲氨基)乙醇 | github.com/IUPAC |

**未取得或不存在：**

| 目标 | 结果 |
|---|---|
| Kodak Alaris 官方 SDS 门户（3E Online） | 姊妹报告已查：没有任何 Kodak Professional 黑白显影液。本次未重查。 |
| photosys.com、kodakprofessional.com 上的 T-MAX Developer SDS | 首页与 `/sds` 页没有 SDS 链接；kodakprofessional.com/chemistry 返回 "Page not found"。 |
| Sino Promise | 搜索摘要称 Kodak Alaris 已把相纸与化学品业务卖给 Sino Promise，并指向 kodaksino.com/photo-chemistry.html（secondary）。该页实际只有广告追踪脚本和 "Do Not Sell or Share My Personal Information"，像停放域名；**没有找到任何 Sino Promise 出具的 T-MAX Developer SDS**。搜索摘要提到一份 Sino Promise 出具的 EKTACOLOR RA 显影补充液 SDS，与本题无关。 |
| T-MAX RS Part A 的 puc.edu 副本 | Cloudflare 质询页，未绕过；改用 digitaltruth.com 副本（内容是同一份 Z17000000414 v2.1）。 |
| emsdiasum.com 的 T-MAX Developer SDS（article 74245，搜索摘要称 2023 年 1 月打印，含 hydroquinone 与 2,2'-oxybisethanol） | 返回 HTML 拦截页；只有 secondary 片段，与下文 1058718 世代一致。 |
| 125px.com 的 CIS-61 镜像 | 不需要；已从 kodak.com 存档取得同一文件。 |
| 公开 T-MAX Developer 配方的 Kodak 专利 | 未找到，见 §9。 |

---

## 0. 简答

| 项目 | 结论 | 来源 |
|---|---|---|
| 配方世代 | SDS 记录了**两种完全不同的化学体系**：(I) **胺–二氧化硫浓缩液**（二乙醇胺–SO₂ 加合物 + 亚硫酸氢盐 + 对苯二酚 + Dimezone S），见于 2007、2013、2014、2019（v11）各版；(II) 2019 年起 product 1058718 的**亚硫酸钾/硼砂体系**（K₂SO₃、HQ、二甘醇、硼砂、KOH、NaBr、HEDTA、Dimezone S），见于 Kodak Alaris 2019 v02 与 PSI 2024 v2.0。 | opened（§1–2） |
| Kodak 参考曲线用的是哪种 | 几乎可以肯定是 (I)。CIS-61 给出工作液比重 1.045，按 1+4 反推浓缩液约 1.225，与 (I) 的 1.21–1.23 吻合，与 (II) 的 1.07 不符；F-4016/F-4017 的 T-MAX Developer 曲线早于 2019。 | CIS-61 opened；反推为 estimate |
| 稀释 | 1 份浓缩液 + 4 份水，稀释倍数 5 | J-86 2017 p. 1 与 2002 版（opened） |
| **工作液 pH** | **8.37 ± 0.05，25 °C**，工作液比重 1.045 ± 0.003（新配槽液，25 °C）。1999、2003、2010 三版数值相同。 | CIS-61（opened，Kodak 自己的站点） |
| 浓缩液 pH | (I)：8.3（2007）、8.3–8.8（2014）、8.8（2019 v11）。(II)：9.7（2019 v02）、9.6（PSI 2024）。 | SDS §9（opened） |
| 缓冲体系 | (I)：以 **DEAH⁺/DEA**（pKa 8.883，25 °C）为主，HSO₃⁻/SO₃²⁻ 为辅。工作液理想溶液缓冲容量约 0.37–0.41 mol L⁻¹ pH⁻¹，约为 D-76 的 7 倍（estimate）。(II)：亚硫酸盐 + 硼酸盐，与新版 HC Developer 同类。 | §4 |
| 初始溴化物 | (I) 的任何一版都**没有列出**溴化物；(II) 列 NaBr < 1 %，即工作液 0–21 mM。(I) 的初始溴化物只能 `assumed`/`fitted`。 | §5 |
| 中心配方（工作液，世代 I） | DEA 0.78 M（0.68–0.89），S(IV) 总量 0.32 M（0.24–0.41），HQ 7.4 g/L = 67 mM（22–113），Dimezone S 0.77 g/L = 3.8 mM（1.2–12），Br⁻ 0（0–21 mM），pH 8.37（25 °C）。 | §6 |
| 不确定度 | 各组分界宽从 ±15 %（DEA）到 ×5（HQ）、×10（Dimezone S）不等。最大的风险不在范围宽度，而在**S(IV) 的记账方式**：2014 版若按字面读（SO₂ 15–20 %），S(IV) 为 0.59–0.89 M，理想溶液 pH 只有 6.3–7.0，与实测 8.37 矛盾。2007 版的"二乙醇胺–SO₂ 络合物"若按 Kodak 自己专利里的原料"含 16 % SO₂ 的二乙醇胺"读，理想 pH 为 8.33–8.62，与实测吻合。 | §4、§7 |

**对 D1 的三点影响：**
1. **pH 用 CIS-61 的实测值做锚点**（8.37 ± 0.05，25 °C，`datasheet`）。不要从 SDS 范围推 pH。
2. **T-MAX Developer 是一种以胺为缓冲的中等亚硫酸盐显影液**：S(IV) 约 0.3 M（D-76 为 0.79 M），HQ 与 D-76 同量级，Dimezone S 是 XTOL 的几倍。它的缓冲容量很大，所以显影产酸引起的局部 pH 下降应当比 D-76、XTOL 小得多。如果 D1 用 T-Max 100 × T-MAX Developer 做留出检验（`datasheets.md` 的建议），这一点是检验内容的一部分。
3. **建两个 YAML 文件**：`t-max-amine`（世代 I，参考数据用这个）和 `t-max-2019`（世代 II，不用于标定）。做法与姊妹报告处理 HC-110 新旧两版相同。

---

## 1. 产品编号与配方世代

| 世代 | SDS（出具方、版本、日期） | product code / PCD | 浓缩液 pH | 比重 / 相对密度 | 状态 |
|---|---|---|---|---|---|
| I-a | Eastman Kodak MSDS v1.1，Rev. 10/02/2007 | 1599844 | 8.3 | "> 1" | opened |
| I-b | Agar Scientific（经销商自编）EU SDS v1，Rev. 06.08.2013 | AGP9167 | "Not determined" | "Not determined" | opened；数值都是 I-c 各范围的上限，见 §2.2 |
| I-c | Kodak Alaris MSDS v1.8，Rev. 07.02.2014；SDS v2.1，Rev. 02/07/2014（打印于 2015、2016） | 1402767、5050851；PCD 5337 | 8.3–8.8 | 1.22–1.23 | opened |
| I-d | Kodak Alaris SDS v11，Rev. 11-14-2019，Issue 08-15-2016 | 5050851；**PCD 8021** | 8.8 | 1.21 | opened |
| II | Kodak Alaris SDS v02，Rev. 06-20-2019，Issue 05-17-2019 | 1058718 | 9.7 | 1.07（specific gravity） | opened |
| II | PSI SDS v2.0，2/22/2024（上一版 1.0，10/11/2023） | 1058718 | 9.6 | 1.07（relative density） | opened |

* J-86（2002 与 2017 版）列出的订货号是 CAT 140 2767（配 1 加仑）和 159 9844（配 5 加仑）。它们分别对应 SDS 上的 1402767 和 1599844（opened）。5050851 是 Freestyle 销售页上"1 L 配 5 L"规格的编号（搜索结果，secondary）。1058718 是现行编号，CineStill 与 B&H 的销售页用它（搜索结果，secondary）。
* I-c 与 I-d 文件号相同（都是 5050851），但 PCD 从 5337 变为 8021，组分记账方式也变了（§2.3）。PCD 的含义 SDS 没有说明；我记得它是 Kodak 的内部配方代号（recalled，未核实）。改号**可能**意味着配方调整，也可能只是重新编制。我没有找到说明。
* 世代 II 与世代 I 没有一个共同的组分记账方式（II 不含任何胺–SO₂ 加合物），只能当作新配方。它的组分比例和新版 HC Developer（姊妹报告 §3.2）很像，各组分大约是 HC 的 1/4–1/5，Dimezone S 更低。这只是观察，不是来源里的说法。
* **推断：** F-4016（2016/2018）的 T-MAX Developer 曲线、F-4017 2016 版、J-86 的显影时间（2017 年 12 月）和 US 5,756,271 Table I（1996）都早于世代 II 首次出现（2019 年 5 月），所以对应世代 I。F-4017 的 Oct 2021 版晚于 2019；姊妹报告 `developer-formulas.md` §7 对同一版里的 HC-110 数据推断为沿用早期版本、未见重测，我对其中的 T-MAX Developer 曲线做同样的推断，未核实，可比较 2016 与 2021 版的曲线矢量来确认。CIS-61 比重的一致性（§3.1）支持世代 I。

---

## 2. 各版 SDS 组分（浓缩液，% w/w）

所有 SDS 都写明：以范围表示的浓度是为保护商业秘密而给出的范围，实际值不公开。SDS 只列出按当时法规需申报的组分，**范围之和不构成封闭配方**。

### 2.1 世代 I-a：Eastman Kodak 2007（opened）

| 组分 | CAS | % w/w |
|---|---|---|
| Water | 7732-18-5 | 50–60 |
| Diethanolamine–sulphur dioxide complex | 63149-47-3 | 35–45 |
| Sodium metabisulphite | 7681-57-4 | 1–5 |
| Hydroquinone | 123-31-9 | 1–5 |
| 4-hydroxymethyl-4-methyl-1-phenyl-3-pyrazolidinone（Dimezone S） | 13047-13-7 | 0.1–1 |

§9：liquid，colourless，odour amine，specific gravity "> 1"，volatile fraction by weight 50–60 %，**pH 8.3**，does not flash。

PubChem 把 63149-47-3 登记为 "Ethanol, 2,2'-iminobis-, compd. with sulfur dioxide (1:?)"，即**化学计量未定**（opened）。

### 2.2 世代 I-b：Agar Scientific 2013（opened，权重低）

二乙醇胺 35.0 %，sulphur dioxide 20.0 %，hydroquinone 5.0 %，sodium hydrogensulphite 5.0 %，Dimezone S 1.0 %。描述为 "Mixture of substances listed below with nonhazardous additions"。§9 pH 与密度均为 "Not determined"；"Organic solvents 35.0 %"。

这些数值恰好是 I-c 各范围的**上限**（DEA 30–35、SO₂ 15–20、HQ 1–5、NaHSO₃ 1–5、Dimezone S 0.1–<1），是经销商为危害分类取的保守值，不是配方值。不进入中心配方。

### 2.3 世代 I-c：Kodak Alaris 2014（v1.8 与 v2.1 内容相同，opened）

| 组分 | CAS | % w/w |
|---|---|---|
| Water | 7732-18-5 | 35–45 |
| Diethanolamine | 111-42-2 | 30–35 |
| Sulphur dioxide | 7446-09-5 | 15–20 |
| Sodium bisulphite | 7631-90-5 | 1–5 |
| Hydroquinone | 123-31-9 | 1–5 |
| Dimezone S | 13047-13-7 | 0.1 – < 1 |

§9：liquid，clear，odour amine，**specific gravity 1.22–1.23**，vapour pressure 24 mbar，boiling > 100 °C（estimated），**pH 8.3–8.8**。

与 2007 版比：水从 50–60 % 降到 35–45 %；"络合物 35–45 %"拆成了 DEA 30–35 % 加 SO₂ 15–20 %，合计 45–55 %；焦亚硫酸钠换成了"亚硫酸氢钠"。按中点，2014 版的 DEA : SO₂ 摩尔比为 (32.5/105.14) : (17.5/64.07) = 1.13 : 1（estimate），接近 1:1 加合物。两版是否同一配方，SDS 本身无法判断（§4.2）。

### 2.4 世代 I-d：Kodak Alaris 2019 v11，PCD 8021（opened）

| 组分 | CAS | % w/w |
|---|---|---|
| 2-(methylamino)ethanol, compound with sulphur dioxide | 21049-70-7 | 20–30 |
| Diethanolamine | 111-42-2 | 15–20 |
| Hydroquinone | 123-31-9 | 1–5 |
| Sodium hydrogensulfite | 7631-90-5 | 1–5 |
| Dimezone S | 13047-13-7 | 0.1–1 |

§9：liquid，clear，odour amine-like，**pH 8.8**，vapor pressure 18 mm Hg，**relative density 1.21**，boiling > 212 °F。

PubChem 把 21049-70-7 登记为 "Ethanol, 2-(methylamino)-, compd. with sulfur dioxide (1:?)"，化学计量同样未定（opened）。2-(甲氨基)乙醇（N-甲基乙醇胺，MMEA，CAS 109-83-1，M 75.11，PubChem opened）的 pKaH 为 9.77 和 9.88（25 °C），IUPAC 数据集 Perrin 补编条目，评级均为 "Approximate" [Douhéret & Pariaud, *J. Chim. Phys.* 59, 1013 (1962)；Christensen, Izatt, Wrathall & Hansen, *J. Chem. Soc. A* 1969, 1212]（opened，引用按 `perrin_supp` 列解析）。

### 2.5 世代 II：product 1058718（opened）

| 组分 | CAS | Kodak Alaris 2019 v02 | PSI 2024 v2.0 |
|---|---|---|---|
| Potassium sulphite | 10117-38-1 | 5 – < 10 | "Potassium Sulfite Solution 45 %" 10–15（≙ K₂SO₃ 4.5–6.75，estimate） |
| Diethylene glycol（2,2'-oxydiethanol） | 111-46-6 | 1 – < 3 | 1–3 |
| Hydroquinone | 123-31-9 | 1 – < 3 | 1–3 |
| Sodium tetraborate pentahydrate | 12179-04-3 | 0.1–1 | 1–3 |
| Potassium hydroxide | 1310-58-3 | 0.1–1 | "Potassium hydroxide 45 %" < 1 |
| Sodium bromide | 7647-15-6 | 未列 | < 1 |
| Dissolvine H-40（HEDTA 三钠盐溶液） | 139-89-9 | 未列 | < 1 |
| Dimezone S | 13047-13-7 | 未列 | < 0.1 |

§9：2019 v02 为 liquid，tan，**pH 9.7**，specific gravity 1.07；PSI 2024 为 liquid，clear，odour "None, Amine"，**pH 9.6**，relative density 1.07，boiling 100 °C。

PSI 版的气味写"Amine"，但没有列出任何胺；可能是未申报的胺，也可能是模板残留，无法判断。

### 2.6 版本间分歧汇总

| 组分 | I-a 2007 | I-c 2014 | I-d 2019 v11 | II 2019 / 2024 |
|---|---|---|---|---|
| 胺 | DEA–SO₂ 络合物 35–45 | DEA 30–35（另有 SO₂ 15–20） | MMEA–SO₂ 20–30 + DEA 15–20 | 无（仅 PSI 气味栏写 amine） |
| 亚硫酸盐 | 络合物中的 SO₂ + Na₂S₂O₅ 1–5 | SO₂ 15–20 + NaHSO₃ 1–5 | 络合物中的 SO₂ + NaHSO₃ 1–5 | K₂SO₃ 约 4.5–10 |
| HQ | 1–5 | 1–5 | 1–5 | 1–3 |
| Dimezone S | 0.1–1 | 0.1–<1 | 0.1–1 | 未列 / < 0.1 |
| 溴化物 | 未列 | 未列 | 未列 | 未列 / NaBr < 1 |
| 浓缩液 pH | 8.3 | 8.3–8.8 | 8.8 | 9.7 / 9.6 |
| 比重 | > 1 | 1.22–1.23 | 1.21 | 1.07 |

---

## 3. 1+4 稀释与工作液浓度

### 3.1 稀释、密度与一致性检验

* **稀释：**"mixing one part of the concentrate with four parts water"，特性栏写 "Concentrate mix ratio 1:4"（J-86 2017 p. 1、p. 9；2002 版相同；opened）。J-86 的 "1:4" 指 1 份加 4 份，稀释倍数 5。
* **工作液：**pH 8.37 ± 0.05，新配槽液比重 1.045 ± 0.003，均在 25 °C（77 °F）（CIS-61，1999 年 10 月、2003 年 2 月、2010 年 5 月三版，opened）。CIS-61 说明：pH 按 CIS-121 测量，用校准良好的 pH 计；比重用 ASTM E100-05 比重计，分度 0.001。
* **一致性检验（estimate）：**体积可加时，SG_work ≈ (SG_conc + 4)/5，所以 SG_conc ≈ 5 × 1.045 − 4 = **1.225**（±0.015，来自 ±0.003）。这与世代 I-c 的 1.22–1.23、I-d 的 1.21 一致。世代 II 的 1.07 会给出工作液约 1.014，与 CIS-61 不符。因此 CIS-61 的 pH 8.37 属于世代 I。
* **换算方法（estimate）：**工作液 g/L = w × ρ_conc × 1000 / 5。中心值取 ρ = 1.225 kg/L 和范围中点，界取范围端点与 ρ 的端点（I-c 用 1.22–1.23；I-a 没有给密度，用 CIS-61 反推的 1.21–1.24，标 `assumed`）。忽略混合时的体积收缩（`assumed`，大小未查证）。

### 3.2 世代 I-c（2014）按字面换算

| 组分 | 浓缩液 % | 工作液 g/L（下界 / 中心 / 上界） | 工作液 mol/L |
|---|---|---|---|
| DEA（M 105.14） | 30–35 | 73.2 / 79.6 / 86.1 | 0.696 / 0.757 / 0.819 |
| SO₂（M 64.07） | 15–20 | 36.6 / 42.9 / 49.2 | 0.571 / 0.669 / 0.768 |
| NaHSO₃（M 104.06，PubChem opened） | 1–5 | 2.44 / 7.35 / 12.3 | 0.023 / 0.071 / 0.118 |
| HQ（M 110.11） | 1–5 | 2.44 / 7.35 / 12.3 | 0.022 / 0.067 / 0.112 |
| Dimezone S（M 206.24） | 0.1–<1 | 0.24 / 1.35 / 2.46 | 0.0012 / 0.0065 / 0.0119 |

按字面读，S(IV) 总量为 0.59–0.89 M，比 D-76 的 0.79 M 还多。§4.2 说明这个读法与实测 pH 不相容。

### 3.3 世代 I-a（2007）：两种读法

"二乙醇胺–SO₂ 络合物"的化学计量未定（PubChem "1:?"）。我比较了三种读法，络合物在工作液中均为 84.7 / 98.0 / 111.6 g/L（下界 / 中心 / 上界）：

| 读法 | 依据 | DEA（mM） | 来自络合物的 SO₂（mM） |
|---|---|---|---|
| 1:1 加合物（SO₂ 占 37.9 %） | 化学计量假设 | 501 / 579 / 660 | 501 / 579 / 660 |
| 2:1（"二乙醇胺亚硫酸盐"，SO₂ 占 23.4 %） | 化学计量假设 | 618 / 714 / 814 | 309 / 357 / 407 |
| **"含 16 % SO₂ 的二乙醇胺"** | Kodak 在 US 5,756,271 Examples 2–5 中用的原料名称："Diethanolamine (with 16% sulfur dioxide)"（opened） | **677 / 783 / 892** | **212 / 245 / 279** |

焦亚硫酸钠 1–5 % → 2.42 / 7.35 / 12.4 g/L = 12.7 / 38.7 / 65.2 mM，相当于 25 / 77 / 130 mM HSO₃⁻。

把"络合物"读成 16 % SO₂ 的原料是一个**假设**（`assumed`）。依据只有两点：这种原料在 Kodak 同时期的液体显影液专利里出现过；这种读法是唯一能让 SDS 组成与实测 pH 同时成立的读法（§4.2）。没有任何来源直接说明 T-MAX Developer 用的是这种原料。

### 3.4 世代 I-d（2019 v11），ρ = 1.21，按 1:1 加合物读

MMEA–SO₂ 络合物 48.4–72.6 g/L，即 MMEA 与 SO₂ 各 348–522 mM；DEA 36.3–48.4 g/L（345–460 mM）；HQ 2.42–12.1 g/L（22–110 mM）；NaHSO₃ 23–116 mM；Dimezone S 1.2–11.7 mM（estimate）。

### 3.5 世代 II（PSI 2024），ρ = 1.07

K₂SO₃ 9.6–14.5 g/L（61–91 mM）；HQ 2.1–6.4 g/L（19–58 mM）；二甘醇 2.1–6.4 g/L；硼砂五水合物 2.1–6.4 g/L（7.3–22 mM，即总 B 29–88 mM）；KOH 0–0.96 g/L（0–17 mM）；NaBr 0–2.1 g/L（**0–21 mM Br⁻**）；Dimezone S 0–0.21 g/L（0–1.0 mM）（estimate）。

---

## 4. pH 与缓冲体系

### 4.1 实测锚点（Kodak，opened）

| 溶液 | pH（25 °C） | 比重（新配槽液，25 °C） | 来源 |
|---|---|---|---|
| T-MAX Developer，工作液（1:4） | **8.37 ± 0.05** | **1.045 ± 0.003** | CIS-61 1999 / 2003 / 2010 |
| T-MAX RS Developer and Replenisher，工作液 | **8.5 ± 0.05** | **1.041 ± 0.003** | 同上 |
| T-MAX Developer 浓缩液（世代 I） | 8.3（2007）；8.3–8.8（2014）；8.8（2019 v11） | 1.22–1.23（2014）；1.21（2019 v11） | SDS §9 |
| T-MAX RS Part A 浓缩液 | 8.9 | 1.21 | SDS 2014 v2.1 §9 |
| T-MAX Developer 浓缩液（世代 II） | 9.7（2019）；9.6（2024） | 1.07 | SDS §9 |

SDS 没有写浓缩液 pH 的测量温度和方法。

**工作液 pH 不在任何 SDS 里**，只在 CIS-61。J-86（2002、2017）没有 pH，只把 "Buffered solution" 列为特性（opened）。

### 4.2 理想溶液电荷平衡检验（estimate）

方法与姊妹报告 §5.1 相同：25 °C，活度 = 浓度，pKa 取姊妹报告 §6（DEA 8.883；HSO₃⁻/SO₃²⁻ 7.19；HQ 9.85；pKw 14.00）。SO₂ 和 HSO₃⁻ 都按 S(IV) 计入；在 pH > 5 时忽略 H₂SO₃。脚本不入库，结果如下：

| 读法（工作液） | DEA | S(IV) 总量 | 理想 pH | 与实测 8.37 比 |
|---|---|---|---|---|
| I-c 2014 字面，中点 | 0.757 M | 0.740 M | **6.31** | 低 2.1 |
| I-c 2014 字面，最利于高 pH 的端点（DEA 最多、SO₂ 与 NaHSO₃ 最少） | 0.819 M | 0.595 M | 7.01 | 低 1.4 |
| I-a 2007，1:1 加合物 | 0.579 M | 0.656 M | < 5（HSO₃⁻ 两性体系） | 不可能 |
| I-a 2007，2:1 | 0.714 M | 0.434 M | 7.63 | 低 0.7 |
| **I-a 2007，"16 % SO₂"原料，中点** | **0.783 M** | **0.322 M** | **8.50** | **高 0.13** |
| I-a 2007，"16 % SO₂"原料，所有端点组合 | 0.68–0.89 M | 0.24–0.41 M | **8.33–8.62** | 包住了实测值 |
| I-d 2019 v11，MMEA–SO₂ 按 1:1，各端点组合 | MMEA 0.35–0.52 + DEA 0.35–0.46 | 0.37–0.64 M | 7.24–8.34 | 多数偏低 |
| 世代 II（PSI 2024）中点 | — | K₂SO₃ 0.076 M + 硼 0.059 M | 9.17 | 无工作液实测；浓缩液 9.6 |

此外，"16 % SO₂"读法对浓缩液本身给出理想 pH 8.50，在 SDS 的 8.3–8.8 之内。

结论：
1. 2014 版的字面记账（SO₂ 15–20 %）**不可能是工作液中真实的 S(IV) 量**。用 SDS 的 DEA（0.70–0.82 M）、NaHSO₃（23–118 mM）、实测 pH 8.37 ± 0.05 和 DEA 的 pKa 反解电荷平衡，S(IV) 只能是 **0.28–0.39 M**（中心 0.335 M）；反过来，要容纳 0.74 M S(IV)，DEA 要 1.78 M（浓缩液约 76 % w/w），这也不可能。这与姊妹报告在旧 HC-110 上遇到的矛盾（SDS 中点给出 pH 6.4–7.3）性质相同。可能的解释是 SDS 编制软件把"络合物"按 1:1 化学计量拆成了 DEA 和 SO₂ 两项；这只是**猜测**，没有来源。
   这个排除不依赖理想溶液假设：姊妹报告发现理想模型在高离子强度下**高估** pH；高离子强度使 HSO₃⁻ 的表观 pKa 下降，也是让 pH 更低。所以活度校正只会让 6.3 离 8.37 更远。
2. 2007 版加上"16 % SO₂"原料的读法自洽：浓缩液和工作液的 pH 都对得上，给出的 S(IV) 0.24–0.41 M 与第 1 条由数据表反解的 0.28–0.39 M 重合。两条路径是独立的：一条只用 2014 版的 DEA 与 CIS-61 的 pH，另一条只用 2007 版的组成与专利中的原料名称。
3. 姊妹报告发现，亚硫酸盐/硼酸盐显影液在 I ≈ 2 M 下，理想模型会把 pH 高估 0.4–0.8。这里的缓冲对主要是 DEAH⁺/DEA，属于阳离子酸，解离前后电荷数的平方和不变，所以离子强度对其表观 pKa 的影响一阶上很小（Debye–Hückel 的一般结论，**recalled**；可在 Bates, *Determination of pH* (1973) 核实）。工作液离子强度估计约 0.9 M（estimate）。所以 8.50 对 8.37 的 0.13 偏差可以接受。不要把它当成证据去调 Na₂S₂O₅ 的取值。

### 4.3 缓冲体系与缓冲容量（estimate）

* **主缓冲对：**DEAH⁺/DEA，pKa 8.883（25 °C），dpKa/dT ≈ −0.025 K⁻¹（Bower, Robinson & Bates 1962，经姊妹报告 §6.2 opened）。在 pH 8.37 时，DEA 有 23.5 % 为游离碱，76.5 % 为 DEAH⁺。
* **次缓冲对：**HSO₃⁻/SO₃²⁻（7.19）。pH 8.37 时 S(IV) 的 94 % 为 SO₃²⁻。HQ（9.85）贡献很小。
* **缓冲容量 β：**中心配方在 pH 8.37–8.50 时约 **0.37–0.41 mol L⁻¹ pH⁻¹**。同一方法下 D-76 为 0.054，XTOL 为 0.048（姊妹报告 §5.2）。所以 T-MAX Developer 的缓冲能力约是它们的 7–8 倍。这与 J-86 的 "Buffered solution"、"Less affected by differences in water supplies" 定性一致（opened）。
* **温度：**CIS-61 的 pH 是 25 °C 值，J-86 的主推荐温度是 20 °C。以 DEA 的 pKa(20 °C) = 9.005、pKw(20 °C) = 14.17 重算，理想 pH 升高 **+0.11**（estimate；亚硫酸 pKa 的温度系数按姊妹报告 ΔH° = −3.65 kJ/mol 视为可忽略）。所以 20 °C 下的工作液 pH 约 8.45–8.5，标 `estimate`，宜在 D1 中由温度模型自洽算出，不作为第二个锚点。

### 4.4 世代 II 的缓冲

亚硫酸钾 + 硼砂五水合物 + 少量 KOH，与新版 HC Developer 同类（姊妹报告 §3.2、§5）。工作液 pH 没有任何实测值。理想溶液中点 9.17 是上限估计。世代 II 不用于 D1 标定。

---

## 5. 初始溴化物与未列出的组分

* **世代 I 的所有版本都没有列出溴化物**（2007、2013、2014、2019 v11；opened）。三种可能：没有溴化物；含量低于申报阈值；或者当时的法规不要求申报 KBr/NaBr。SDS 本身无法区分。
* 世代 II 列出 NaBr < 1 %，工作液 0–21 mM（opened 加 estimate）。
* J-86 写着 T-MAX RS 不需要 "starter"（opened），即开槽不另加溴化物起始液。这说明不了配方里本来有没有溴化物。
* **建议：**`bromide_initial` 中心值取 0，界 0–21 mM（上界借用世代 II "< 1 % NaBr" 的写法，属**类比**，标 `assumed`），在 D1 中作为 `fitted` 参数，用 T-MAX Developer 的灰雾和趾部数据标定。
* 其他可能未列出的组分：螯合剂（同类 Kodak 产品 HC Developer、XTOL、D-76 都有 DTPA 或 HEDTA，姊妹报告 §1.2、§2.2、§3.2）、苯并三唑类防灰雾剂（DE 1,118,005 Example 1 有，姊妹报告 §3.1）、二甘醇之类的溶剂。以上都是类比推测，**没有来源证明 T-MAX Developer 含有这些组分**。
* 质量闭合检验（2014 中点）：水 40 + DEA 32.5 + SO₂ 17.5 + NaHSO₃ 3 + HQ 3 + Dimezone S 0.55 = 96.55 %，剩余约 3.5 % 未申报（estimate；范围中点之和本身就有 ±15 % 量级的不确定度，这个数只说明有余量，不说明余量是什么）。

---

## 6. 推荐的中心配方（世代 I，`msds-estimated`）

工作液（1 份浓缩液 + 4 份水），每升。基础读法为 2007 版组成，"二乙醇胺–SO₂ 络合物"按"含 16 % SO₂ 的二乙醇胺"处理，浓缩液密度 1.225 kg/L（1.21–1.24）。

| 字段（建议） | 中心值 | 下界 | 上界 | 单位 | 来源标签与依据 |
|---|---|---|---|---|---|
| `dilution_factor` | 5 | 5 | 5 | — | `datasheet`：J-86 p. 1 |
| `diethanolamine_total` | 82.3（0.783 M） | 71.1（0.677 M） | 93.7（0.892 M） | g/L | `datasheet` 范围 + `assumed`（16 % SO₂ 原料）；2014 字面读法给出 0.70–0.82 M，落在本区间内 |
| `sulfur_dioxide_from_amine_complex` | 15.7（0.245 M） | 13.6（0.212 M） | 17.9（0.279 M） | g/L | 同上 |
| `sodium_metabisulfite` | 7.35（38.7 mM，= 77 mM HSO₃⁻） | 2.42（12.7 mM） | 12.4（65.2 mM） | g/L | `datasheet`：2007 版 1–5 %。2014 版记为 NaHSO₃ 1–5 %，即 2.44–12.3 g/L（23–118 mM HSO₃⁻），在 S(IV) 计量上与此重叠 |
| → `sulfite_total_S_IV` | 0.322 M | 0.237 M | 0.409 M | mol/L | 由上两行相加（estimate）。**主要依据**是电荷平衡反解：2014 版 DEA 范围加 CIS-61 pH 8.37 ± 0.05 要求 0.28–0.39 M（§4.2 第 1 条，`datasheet` + `literature` pKa）；16 % 读法是独立印证。**另一读法**（2014 字面）给出 0.59–0.89 M，与实测 pH 不相容（§4.2），列为已排除 |
| `hydroquinone` | 7.35（66.8 mM） | 2.42（22.0 mM） | 12.4（113 mM） | g/L | `datasheet`：1–5 %，取算术中点 |
| `dimezone_s` | 0.77（3.8 mM） | 0.24（1.2 mM） | 2.48（12.0 mM） | g/L | `datasheet`：0.1–1 %。中心值取几何中点 0.32 %（`assumed`），因为范围跨一个数量级；同时它与 DE 1,118,005 Examples 1、8 的 HQ : 吡唑烷酮 = 10 : 1 质量比一致（姊妹报告 §3.1，opened）。算术中点 0.55 % 会给出 1.35 g/L（6.5 mM） |
| `bromide_initial` | 0 | 0 | 21 | mM | `assumed`，§5；在 D1 中按 `fitted` 处理 |
| `pH_25C` | 8.37 | 8.32 | 8.42 | — | `datasheet`：CIS-61 |
| `pH_20C` | 8.48 | — | — | — | `estimate`：+0.11，§4.3；应由模型算出 |
| `specific_gravity_25C` | 1.045 | 1.042 | 1.048 | — | `datasheet`：CIS-61；用于检查配方总量 |
| `buffer` | DEAH⁺/DEA（主），HSO₃⁻/SO₃²⁻（次） | — | — | — | §4.3 |

物种数据复用姊妹报告：二乙醇胺（§6.2，M 105.14，CAS 111-42-2，pKa 8.883）、对苯二酚（§6.1）、Dimezone S（§6.1，M 206.24，CAS 13047-13-7，pKa 未找到）、焦亚硫酸钠（§6.2，M 190.11，CAS 7681-57-4，K(2HSO₃⁻ ⇌ S₂O₅²⁻) = 0.032 kg/mol）、亚硫酸（§6.2，pKa 1.86 / 7.19）、扩散系数（§6.3，DEA 0.81/0.76，Dimezone S 0.57/0.54 ×10⁻⁵ cm²/s，20 °C，estimate）。新增：NaHSO₃ M 104.06（PubChem opened），MMEA（§2.4）。

**实现建议：**D1 若用电荷平衡算 pH，可以直接用上表的 DEA 和 S(IV)。若改用"pH 固定为锚点"的做法，则把 DEA 和 S(IV) 当作总量参数，质子化分数由 pH 8.37 和 pKa 决定。两种做法下，中心配方的理想 pH 相差 0.13，在可接受范围内。

---

## 7. msds-estimated 配方的不确定度

按影响从大到小排列：

1. **配方世代。**有四种记账方式和至少两种化学体系（§1）。选错世代，S(IV) 差 4 倍以上，缓冲体系完全不同。本报告的选择（世代 I）有 CIS-61 比重和时间顺序两方面的支持，但仍是推断。
2. **S(IV) 的记账。**2007 版加"16 % SO₂"读法给出 0.24–0.41 M；由 2014 版 DEA 和实测 pH 反解给出 0.28–0.39 M；2014 版字面读法给出 0.59–0.89 M。前两者一致且与实测 pH 自洽，第三者被排除（§4.2）。剩余风险在于 DEA 的真实含量：反解结果随 DEA 线性变化。亚硫酸盐浓度直接影响溶解物理显影和颗粒，是 D1 中最敏感的未知量之一。
3. **范围宽度。**HQ 1–5 %（中心值 ±67 %），Dimezone S 0.1–1 %（中心值 ×/÷ 3.2），亚硫酸氢盐 1–5 %。这三项决定活度和超加和性。D1 应把 HQ 与 Dimezone S 的工作浓度设为带先验界的 `fitted` 量，先验界取上表。
4. **未申报组分。**溴化物、螯合剂、防灰雾剂可能存在而未列出（§5）。其中溴化物影响最大。
5. **密度与体积收缩。**ρ_conc 由 SDS（1.21–1.23）和 CIS-61 反推（1.225 ± 0.015）确定到约 ±1.5 %；体积收缩忽略。这一项在总不确定度中最小。
6. **SDS 范围的含义。**范围是"为保护商业秘密而给的区间"，没有声明真实值在区间内均匀分布，也没有声明各版本之间区间是否对应同一真实值。所以本报告取中点和取并集的做法都属于 `assumed`。

粗略总评：DEA 约 ±15 %，S(IV) 约 ±25 %（以 §6 的读法为条件；不加条件时是 ×2–3），HQ 约 ×/÷ 2–3，Dimezone S 约 ×/÷ 3，溴化物完全未知。pH 是唯一精确到 ±0.05 的量。

---

## 8. T-MAX RS Developer and Replenisher

| 部分 | 组分（% w/w） | §9 | 来源 |
|---|---|---|---|
| Part A（product 5054184 – Part A，PCD 5492） | DEA（111-42-2）30–35；SO₂（7446-09-5）15–20；HQ（123-31-9）1–5；NaHSO₃（7631-90-5）1–5；Dimezone S（13047-13-7）0.1–<1；"Additional Components Include: Water" | specific gravity **1.21**，**pH 8.9** | SDS Z17000000414 v2.1，Rev. 02/07/2014，Print 09/24/2015（opened，digitaltruth.com 副本） |
| Part B（5054184 – Part B，PCD F1220） | 二甘醇（111-46-6）90–99；乙酸（64-19-7）1–5；1,4-diphenyl-3-(phenylammonio)-1H-1,2,4-triazolium（2218-94-2，PubChem 名 "Nitron"，M 312.4）0.1–<1 | specific gravity **1.12**，pH "No data available"，colour amber，odour vinegar | SDS Z33000000480 v2.1，Rev. 02/07/2014（opened，两份副本内容一致） |
| 搜索摘要中的 Part A 8446163（PCD 5492） | 同 Part A | — | secondary（puc.edu 副本无法打开）。8446163 与 J-86 的 CAT 844 6163（配 1 加仑）一致 |

* Part A 的记账方式与 T-MAX Developer 2014 版完全相同，只是没有给出水的范围。所以 §4.2 的矛盾同样适用：字面读 SO₂ 15–20 %，与 Part A 自身的 pH 8.9 不相容。
* **工作液：**pH 8.5 ± 0.05，比重 1.041 ± 0.003，25 °C（CIS-61，opened）。J-86 把 RS 描述为 "hydroquinone-based, two-part developer"，混合液直接作槽液或补充液，每卷 135-36 或 120（或每张 8×10 英寸）补充 45 mL（J-86 2017 p. 1–2，opened）。
* **稀释：**J-86 只给"配 1 加仑 / 配 10 加仑"的包装规格，**没有给出 Part A、Part B 与水的体积比**（2002、2017 版都查过，opened）。论坛帖子里有 1+4 和 1+9 两种说法，属 secondary 中的论坛级证据，不采用。按比重估算（estimate）：若 Part B 的体积可以忽略，1 + 0.21 f_A = 1.041，得 Part A 体积分数 f_A ≈ 0.20，相当于约 1+4；Part B 每占 1 % 体积，f_A 约减少 0.006。这只是比重一致性估算，不是配方信息。
* Part B 中 Nitron 的作用不明：我记得它是分析化学中的硝酸根沉淀剂（recalled，未核实），在显影液中的作用没有找到来源。二甘醇作为溶剂载体。都只是记录，不作解释。
* D1 不需要 RS 的完整配方；它只与 T-Max 100 的 T-MAX RS 大槽曲线有关（`datasheets.md`）。如果要用，分类为 `msds-estimated`，pH 锚点 8.5（25 °C）。

---

## 9. 专利检索

* **没有找到直接描述 T-MAX Developer 配方的 Kodak 专利。**检索词包括 "T-MAX developer"、"diethanolamine-sulfur dioxide"、CAS 63149-47-3、"4-hydroxymethyl-4-methyl-1-phenyl-3-pyrazolidone" 与 "diethanolamine" 加 "concentrate" 的组合（WebSearch；Google Patents 的 XHR 查询接口被本会话的权限系统拒绝，未再尝试）。
* **最接近的结构先例：**DE 1,118,005 / GB 958,678（Kodak，1959 年优先权；姊妹报告 §3.1，opened），要求保护羟烷基胺 + 其 SO₂ 加成物 + HQ + 1-苯基-3-吡唑烷酮类的浓缩液。T-MAX Developer 的 SDS 组成符合这个类型，但该专利的实施例是无水浓缩液，吡唑烷酮也不是 4-羟甲基-4-甲基取代物；而 T-MAX Developer 含 35–60 % 水。说 T-MAX Developer 源于此专利是**推断**。
* **US 5,756,271**（Kodak，1996 年优先权，opened）：
  * Examples 2–5 用 "Diethanolamine (with 16% sulfur dioxide)" 作原料，与 Dimezone S 搭配。这是 §3.3 读法的唯一依据。
  * Table I 用 "KODAK T-MAX™ Developer" 作对照，写明它 "containing hydroquinone"，但**没有给出它的组成、稀释或 pH**。所以 Table I 只刻画性能，不刻画组成。
* US 4,987,060（3M，1990，opened）是 HQ/亚硫酸盐无胺浓缩液，与 T-MAX Developer 无关，记此以免重复检索。

**US 5,756,271 Table I 中的 T-MAX Developer 行**（opened；20 °C，8 L 槽，氮气脉冲搅拌每 15 s 1 s，21 级碳灰阶，3–15 min 显影系列内插）。D-76 行列出以便对比。原表列名为 "Time to CI 0.56 (Normal)" 与 "Time to CI 0.82 (+3 Stops)"。

| 胶片（1996 年乳剂） | 显影液 | Log E at ANSI 0.62 | EI | t(CI 0.56) min | t(CI 0.82) min | Dmin @ CI 0.56 | Dmin @ CI 0.82 |
|---|---|---|---|---|---|---|---|
| TMAX 100 | T-MAX | −2.26 | 160 | 6.7 | 11.6 | 0.23 | 0.27 |
| TMAX 100 | D-76 | −2.10 | 100 | 6.7 | 10.2 | 0.22 | 0.27 |
| TMAX 400 | T-MAX | −2.79 | 500 | 6.0 | 9.5 | 0.26 | 0.30 |
| TMAX 400 | D-76 | −2.76 | 500 | 7.6 | 11.0 | 0.24 | 0.25 |
| TMAX P3200 | T-MAX | −3.05 | 1000 | 7.4 | 11.3 | 0.34 | 0.40 |
| TMAX P3200 | D-76 | −3.00 | 800 | 9.6 | 14.0 | 0.33 | 0.38 |
| TRI-X PAN | T-MAX | −2.77 | 500 | 6.7 | 12.0 | 0.32 | 0.36 |
| TRI-X PAN | D-76 | −2.65 | 400 | 6.5 | 12.0 | 0.26 | 0.28 |
| PLUS-X PAN | T-MAX | −2.34 | 160 | 6.2 | — | 0.42 | — |
| PLUS-X PAN | D-76 | −2.20 | 125 | 5.6 | 12.0 | 0.30 | 0.32 |

与 D-76 比，T-MAX Developer 一贯给出 0.03–0.16 log E 的感光度优势，Dmin 也更高（Tri-X 高 0.06，Plus-X 高 0.12）。这与"高缓冲、胺活化"的配方性质定性相符；作为 D1 的趋势检验可以用，但胶片是 2007 年以前的乳剂（姊妹报告 §7 第 11 条），不能当参考曲线。

---

## 10. 顺带发现（影响姊妹报告）

1. **旧 HC-110 的工作液 pH 有 Kodak 实测值了。**CIS-61（1999 / 2003 / 2010）：HC-110 Developer dilution A（1:15）pH **8.9 ± 0.10**，比重 1.015 ± 0.003；dilution B（1:31）pH "—"，比重 1.008 ± 0.003；HC-110 Developer Replenisher dilution A pH 9.1 ± 0.10，比重 1.023；均为 25 °C（opened）。这三版都早于 2019 年的配方更换，所以是旧 HC-110 的值。姊妹报告 `developer-formulas.md` §3.1 写的是"未找到"，并用专利类比值 9.6–10.0 作先验；这条应当更新，HC-110 的 pH 可以从 `fitted` 改为 `datasheet` 锚点（dilution A）。
2. **D-76（成品包装）pH：**CIS-61 给出 D-76 stock **8.5 ± 0.05**，比重 1.078 ± 0.003；1:1 工作液比重 1.039，pH "—"；D-76R stock 比重 1.091（25 °C，opened）。CIS-61 称这些数值为 "Typical pH Values"，测量方法按 CIS-121，属于厂商给出的典型值，不是单次测量记录。它与 Carlton & Crabtree 1929 的比色值 8.5 一致。姊妹报告 §1.3 可以补上这一条。
3. **XTOL：**CIS-61 给出 8.2 ± 0.05 和 1.085 ± 0.003，与 J-109 相同（opened）。
4. **Shopify CDN `0339/5113` 是 CineStill 的商店。**同一路径下的 Df96 说明书标题为 "CINESTILL INC, 100 Latona Road"（搜索结果标题，secondary）。姊妹报告的 HC-110、XTOL、D-76 SDS 都来自这个 CDN，它在报告里写过"无法确定这家商店与 Kodak Alaris 或 PSI 的关系"。现在可以确定这是经销商托管。
5. 本会话里 web.archive.org 可以访问（姊妹报告当时 "temporarily offline"）。kodak.com 的 techPubs 存档（CIS-61、J-86 HTML 版）可以用来补全其他 Kodak 出版物的旧版本。

---

## 11. 分歧与未解决问题

1. **世代 I 的 S(IV) 量：**0.24–0.41 M（2007 版加"16 % SO₂"读法，与 pH 自洽）对 0.59–0.89 M（2014 版字面读法，与 pH 矛盾）。没有一手来源能裁决。
2. **2014 → 2019 v11 的 PCD 变化（5337 → 8021）**和 MMEA–SO₂ 的出现：重新编制，还是配方改变？未知。参考数据的年代（2016–2018 版 F-4016，2016 / 2021 版 F-4017）跨过了 2016 年的 Issue date，所以严格说，参考曲线对应的是 I-c 还是 I-d 不确定；两者 pH 与密度相近（8.3–8.8 对 8.8；1.22–1.23 对 1.21）。
3. **世代 II 的 Dimezone S：**Kodak Alaris 2019 版未列，PSI 2024 版 < 0.1 %。与姊妹报告 XTOL Part A 的 PSI "< 0.01 %" 一样，可能是 SDS 低报。
4. **世代 II 的气味 "Amine"** 与组分表中没有胺之间的矛盾。
5. **初始溴化物**在世代 I 完全没有来源。
6. **T-MAX RS 的混合体积比**没有一手来源。
7. **Dimezone S 的酸 pKa** 仍未找到（姊妹报告 §7 第 12 条）。在 pH 8.4 下它的离子化程度未知，影响超加和模型。
8. **温度：**CIS-61 的 pH 是 25 °C 值；20 °C 下约 +0.11 只是理想溶液估计。
