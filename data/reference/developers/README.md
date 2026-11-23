# 显影液手册数据（暂存）

这个目录放的是从柯达显影液技术手册里整理出来的表格和数值，每种显影液一个文件。内容只收与两种胶片有关的部分：KODAK PROFESSIONAL TRI-X 400（400TX，材料 `kodak/tri-x-400`）和 KODAK PROFESSIONAL T-MAX 100（100TMX，材料 `kodak/t-max-100`）。

| 文件 | 显影液 | 出版物 | 版次 |
|------|--------|--------|------|
| `d-76.yaml` | D-76 | J-78 | 2017 年 12 月（Revised 12-17） |
| `xtol.yaml` | XTOL | J-109 | 2018 年 2 月（Revised 02/2018） |
| `hc-110.yaml` | HC-110 | J-24 | 2017 年 12 月（Revised 12-17） |
| `t-max-dev.yaml` | T-MAX Developer | J-86 | 2017 年 12 月（Revised 12-17） |
| `t-max-rs.yaml` | T-MAX RS Developer and Replenisher | J-86 | 同上 |

## 来源

| 出版物 | 网址 | sha256 |
|--------|------|--------|
| J-78 | https://web.archive.org/web/20211017161519id_/https://imaging.kodakalaris.com/sites/default/files/files/resources/j78.pdf | `d1998e9e3807bb8b4c1969ce44e5f9d7750f8043e5a19142a99b83eb546d9915` |
| J-109 | https://kodakprofessional.com/sites/default/files/wysiwyg/pro/chemistry/J-109_Feb_2018.pdf | `10baa7d8875eb6abceda8a7020e07f67a1fa8281f244ddb489feb67785838f87` |
| J-24 | https://kodakprofessional.com/sites/default/files/wysiwyg/pro/chemistry/j24.pdf | `e66818bc246a63db1f7667a703c32d9eeb9ea4cd4c73908f913023f2828171ce` |
| J-86 | https://kodakprofessional.com/sites/default/files/wysiwyg/pro/chemistry/j86.pdf | `a20c0d249d04e8feaeb1240eca457f782f0e57dcd510e26e8a9b1c73bc13e00c` |

J-78 在 kodakprofessional.com 上已经找不到（原路径返回 404）。这里用的是 Wayback Machine 对 imaging.kodakalaris.com 的存档，抓取时间 2021-10-17，封底印 "Revised 12-17"、© 2017 Kodak Alaris。同一网站更早的路径（`/sites/prod/files/...`，2016 年存档）存的是 © 2002 的旧版，没有用。PDF 放在 `.cache/datasheets/pdf/<出版物>-<版次>.pdf`，不提交。

## 文件结构

每个文件的顶层键：

- `sources`：出版物编号、标题、版次、网址、本地文件名和 sha256。
- `notes`：收录范围，以及没收什么、为什么没收。
- `facts`：单个数值和小表格，比如容量、时间补偿、补充液用量、稀释定义、pH、保存期限。每一项都写明 `source`、印刷页码 `page` 和位置 `where`；值写在 `value` 和 `unit` 里，或者按行写在 `rows` 里。出版物里没有的就写 `null`，`statement` 说明"没有写"。
- `agitation`：出版物自己对搅动方式的定义（小罐、大罐、盘、转筒、气体搅动），写出间隔、时长、翻转次数这些原文给出的量。各个显影时间表的 `agitation_ref` 指向这里的 `id`。不同出版物的定义没有互相借用：J-24 只在表头和脚注里写了间隔，就只记间隔。
- `development_times`：显影时间表。原表按列分组的（比如 J-109 表 1 的"原液"和"1:1"两组），拆成两个条目，`column_group` 记原表的组标题。每行记：
  - `film_as_printed`：表中印的胶片名，照抄（J-24 的 "PROFFESSIONAL" 拼错也照抄）；
  - `material`：对应的材料 id；
  - `format`、`ei`、`ci`：画幅、曝光指数、对比度指数。`ei` 是字符串，一般就是印刷原样；
  - `printed` 和 `minutes`：每个温度一格，`printed` 是印刷原样（"6¾"、"4½†"、"NR"、"—"），`minutes` 是十进制分钟。NR、"—" 和读不出的格记为 `null`；
  - `bold`：只在出版物规定了粗体含义时才记，每个温度一格，`true` 或 `false`。J-109 的粗体表示标称感光度，J-86 的粗体表示建议的起始时间。J-109 的 EI 格是粗体时，另记 `nominal_speed_bold: true`，这时 `bold` 数组一定写出，时间格不是粗体就写 `false`，不省略。

温度列表 `temperature_c` 按原表的列序。

**印刷错误的改正（`_as_printed`）。** 只有别的证据明确指出印错时才改，改法是：字段本身写改正后的值，旁边的 `<字段>_as_printed` 保留印刷原样，`note` 写明证据。目前有两处：

- J-78 转筒表（第 5 页）的表头印成 18、21、20、22、24 °C。按这个表头读，18 行里有 13 行在 21 °C 的时间比 20 °C 长，没有一行更短；把第二、三列的表头对调后，没有一行的时间随温度升高而变长。所以 `temperature_c` 记 [18, 20, 21, 22, 24]（`temperature_f` 同样对调），`temperature_c_as_printed` 记印刷顺序。时间格仍按印刷列序，没有挪动。这张表只有旧产品名（`material: null`），也没有别的出版物能直接印证。
- J-109 表 4（第 13 页）400TX 那一行的 EI 印成 "160"。F-4017 第 6 页 400TX 大罐 EI 1600 的 XTOL 时间在 18、20、21、24 °C 与这一行逐格相同，CI 0.72 也是表 1–3 里 EI 1600 的目标，所以 `ei` 记 "1600"，`ei_as_printed` 记 "160"。

## 怎么读的

这些表格不是配方生成的。步骤如下：

1. 用 `tools/digitize/textrows.py` 把表格页的文字层按行打印。这个脚本把三种分数写法都合并成一个记号：J-109 是分子和分母上下错开的两段，J-78、J-24、J-86 是 "61⁄2" 或 "6 ½"。它还给粗体字（字体名含 Semibold 或 Bold）加 `^B` 标记。
2. 按表头的温度位置把每个数字分到列里。每行的格数必须和表头列数相等，不等就报错。这一步和组装 YAML 用的是没有提交的临时脚本；任何一格都可以用 `textrows` 加渲染页复核，每行都写明了页码和表名。
3. 把所有页面渲染成 PNG，和抽出的每一格逐个对照。

文中提到的出版物里的数值都是这样读的。复现命令：

```bash
uv run --group digitize python -I -m tools.digitize.textrows .cache/datasheets/pdf/j109-2018-02.pdf 6
uv run --group digitize python -I -m tools.digitize render .cache/datasheets/pdf/j109-2018-02.pdf 6 out.png --scale 2
```

J-78 的增感表（第 4 页）和 J-86 的机器冲洗表（第 8 页，时间是"6 to 8"这样的区间）版式特殊，是对照文字层手工录入的，同样逐格核对过渲染页。区间记在 `minutes_range` 里。

## 没有曲线

四份显影液手册里都没有特性曲线，也没有对比度指数曲线。检查方法：逐页查找非水平、非竖直的长矢量路径和嵌入位图，再把每一页渲染出来看。J-109 有三张位图，是第 1、16 页的 Kodak Alaris 标志和封底的黑色横条。J-86 第 3、10 页有一幅手腕翻转搅动的示意图，不是数据曲线。所以这次没有写 `j109-...` 这类曲线 CSV，也没有核对图。

## 各出版物要注意的地方

**J-78（D-76）**

- 表里只有旧产品名："T-MAX 100 Professional" 和 "TRI-X Pan"，没有 100TMX 和 400TX。这些行照录，但 `material: null`，`legacy_name: true`。它们的时间和现行胶片手册对不上：小罐原液 20 °C，J-78 的 T-MAX 100 Professional 是 9 分钟，F-4016 的 100TMX 是 6½ 分钟；J-78 的 TRI-X Pan 是 8 分钟，F-4017 的 400TX 是 6¾ 分钟。D-76 用于这两种胶片时，以胶片手册为准。
- 小罐原液表中 TRI-X Pan 21 °C 那一格，文字层和渲染页上都印着 "7v"，读不出来，记为 `null`。
- 转筒表（第 5 页）的 20 °C 和 21 °C 两列表头印反了，`temperature_c` 已按上文改正，印刷顺序留在 `temperature_c_as_printed`。
- "1:1" 在 J-78 里没有用文字说明是哪两种体积之比。手册也没有配方、配制步骤、pH 和对比度指数目标值。
- 容量表里 "16 (4)23" 末尾的 2、3 是上标脚注号，不是数字的一部分。

**J-109（XTOL）**

- 1:1 的定义是 developer:water（第 2 页），每卷至少要用 100 mL 原液来稀释。
- pH 只出现在一句话里：工作液比重 1.085 ± 0.003（25 ± 0.3 °C 下测），"at pH 8.2 ± 0.05"，指原液工作液。
- 表 4（大罐、新液）400TX 那一行的 EI 印成 "160"（CI 0.72），应该是 1600。`ei` 记改正后的 "1600"，`ei_as_printed` 记 "160"，证据见上文。
- 表 1 中 400TX 只列了 135 画幅，表 2–4 列的是 135/120。
- 大罐的搅动：表 3、7 是氮气搅动，每 10 秒冲 2 秒；表 4、8 是每分钟一次手工搅动，但 J-109 没有描述手工动作本身。转筒没有给搅动参数。

**J-24（HC-110）**

- 原液 = 1 份浓缩液 + 3 份水。B 稀释 = 原液 1:7 = 浓缩液 1:31，意思是 1 份对 7 份水（1 L 里 125 mL 原液加 875 mL 水），两种写法浓度相同。
- 2017 年 12 月版只讲一种浓缩液，没有区分新旧配方。
- 表头 24 °C 那一列印成 "24°F"，下面一行是 "(75°F)"，是印刷错误。
- 配制表里有三处明显的印刷错误：原液表 B 稀释 5 L 一行的 "4.38 mL"（应为 L），浓缩液表 A 稀释 5 L 一行的 "313 L"（旁边水是 4687 mL，应为 mL），以及 A 稀释 13.3 L 一行的 "828 L (28 fl oz)"（应为 mL）。这几个数都没有使用。
- 没有 pH，也没有对比度指数目标值，只说按漫射放大机印放调到合适反差，用聚光放大机时时间减约 30%。
- 只收现行产品的行（PROFESSIONAL T-MAX 100、PROFESSIONAL TRI-X 400 Film / 400TX）。旧产品 TMX、TX 的行没有收。

**J-86（T-MAX Developer、T-MAX RS）**

- 目标对比度指数：T-MAX 400 是 0.60，其他胶片 0.56。
- T-MAX Developer 推荐温度是 24 °C（第 12 页）。它的表里，T-MAX 胶片的粗体格都在 24 °C 列，400TX 没有粗体格。
- 小罐搅动有两种说法：第 3、10 页是每 30 秒 5–7 次翻转，时长 5 秒；第 12 页速查表是每 30 秒约 4 次（T-MAX 胶片 5 秒 5 次）。两种都记下了，没有调和。
- 第 12 页的温度容差印成 "?1.7°C [?3°F]"，符号字形丢了，大概是 ±。
- 时间补偿表把胶片规格印成 "136-36"。
- T-MAX RS 没有给混合比例和未补充时的容量。

## 与胶片手册的交叉核对

拿这里的显影时间去对 `data/reference/kodak/tri-x-400/tables.yaml`（F-4017）和 `data/reference/kodak/t-max-100/tables.yaml`（F-4016）。比较的是同一显影液、同一罐型、两边都有的温度：

XTOL 只和 J-109 的新液时间比：表 1（小罐）、表 2 新液组和 1:1 组（转筒卷片）、表 4（大罐卷片新液）、表 5（盘冲页片）、表 6 新液组和 1:1 组（转筒页片）、表 8（大罐页片新液）。胶片手册的时间也是按新液给的；表 2、3、6、7 里陈液、补充液那几组本来就不同，不算不一致。

**两边都有的组合，逐格一致：**

- 400TX 标准冲洗：小罐的 T-MAX、T-MAX RS、HC-110 (B)、XTOL、XTOL 1:1；大罐的 T-MAX RS、HC-110 (B)、XTOL；转筒的 T-MAX、T-MAX RS、HC-110 (B)、XTOL、XTOL 1:1。
- 400TX 增感：EI 1600 和 EI 3200 小罐、转筒的 T-MAX、T-MAX RS、XTOL、XTOL 1:1；大罐 EI 1600 的 T-MAX RS 和 XTOL（J-109 表 4 那一行 EI 印成 "160"），大罐 EI 3200 的 XTOL。
- 100TMX：小罐卷片的 T-MAX (1:4)、T-MAX RS、XTOL 1:1、HC-110 (B)；大罐卷片和页片的 T-MAX RS、XTOL、HC-110 (B)；盘冲页片的 T-MAX RS、XTOL、XTOL 1:1、HC-110 (B)；转筒卷片的 T-MAX (1:4)、T-MAX RS、XTOL、XTOL 1:1、HC-110 (B)；转筒页片的 T-MAX RS、XTOL、XTOL 1:1、HC-110 (B)。小罐卷片 XTOL 原液只有一格不一致：18 °C，F-4016 第 3 页印 8½，J-109 表 1（CI 0.56）印 8¾。两份出版物各照原样记录，没有改。

比较只覆盖两边都印了时间的温度：J-86 的 T-MAX 显影液表和 T-MAX RS 的大罐、盘、转筒表从 20 °C 起，J-109 没有 22 °C 一列，XTOL 1:1 从 20 °C 起。

**没法比的组合（显影液手册里没有对应的表，不算一致）：**

- D-76 和 D-76 1:1：J-78 只有旧产品名。
- 400TX 和 100TMX 大罐的 T-MAX 显影液：J-86 没有 T-MAX 显影液的大罐表（F-4017 在这里印 NR）。
- 400TX 大罐的 XTOL 1:1：J-109 没有 1:1 的大罐表（表 3、4 只有原液）。
- 400TX 增感的 HC-110 (B)：J-24 没有增感表。
- 400TX 大罐 EI 3200 的 T-MAX RS：J-86 的大罐卷片表里 400TX 只有 400/800 和 1600 两行。
- 100TMX 的 T-MAX 1:7、1:9 和 T-MAX RS 1:7、1:9：J-86 只给 T-MAX 1:4 和 T-MAX RS 原液的表。

## Rodinal

Rodinal 没有收录。D1 的抓取规则只允许从柯达自己的网站（及其 Wayback 存档）下载 PDF，Agfa 或 ADOX 的网站不在允许范围内，所以没有去取它们的技术说明。另外，柯达的出版物里没有任何 Rodinal 的数据。
