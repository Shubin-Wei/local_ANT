# 详细教程：从原始数据到三维横波速度模型

本教程按实际操作顺序写，每一步说明做什么、为什么、命令、输出在哪里、怎样判断结果正常、常见问题。参数的含义、单位、依据和取值范围统一在 [parameters.md](parameters.md)；想让 AI 协助完成全流程，见 [ai_workflow.md](ai_workflow.md)。

适用范围：密集台阵（节点仪或宽频带）的垂直分量连续记录，台间距几百米到十几千米，目标是 0.1–2 s 左右（0.5–10 Hz）的 Rayleigh 波相速度和浅部（地表到约 2–3 km）三维 Vs。面状台阵直接三维反演；线状台阵沿测线反演（窄条网格）。更长周期、更大台阵也能用，但要按 parameters.md 调整频带、时延、网格和深度节点。

---

## 1. 安装

### 1.1 Python 环境和依赖

需要 conda（Miniconda 或 Anaconda）。在仓库根目录：

```bash
git clone https://github.com/Shubin-Wei/local_ANT.git
cd local_ANT
conda env create -f environment.yml      # 创建名为 lant 的环境：Python 3.10、ObsPy、GMT、编译器等
conda activate lant
pip install -e .                         # 安装 lant 命令（代码改动立即生效）
```

`conda env create` 很久（半小时以上）停在 “Solving environment” 时，是旧版 conda（23.10 以前）的求解器太慢：用 `conda env create --solver libmamba -f environment.yml`（先 `conda install -n base conda-libmamba-solver`），或用 mamba / micromamba：`micromamba create -f environment.yml -n lant`。

不用 conda 也可以：Python ≥ 3.9 上 `pip install -e .` 装 Python 依赖，另装 GMT 6（≥ 6.1）、gfortran、gcc/g++、make、git、FFTW（单精度，`libfftw3f`）。

### 1.2 第三方程序

仓库里不放第三方程序的源码，由脚本从官方仓库下载固定版本并编译到 `external/`：

```bash
bash tools/install_external.sh            # DSurfTomo + fwd1d、pDSurfTomo、CC-FJpy 三个都装
bash tools/install_external.sh dsurftomo  # 只装 DSurfTomo（借用参考曲线、不用并行版时足够）
```

| 程序 | 用在哪里 | 必需 |
|---|---|---|
| DSurfTomo（Fang et al. 2015，MIT） | 三维反演；`fwd1d` 用它的正演程序算一维起始模型 | 是 |
| pDSurfTomo（Zhu et al. 2026，GPL-3） | 并行反演，结果与 DSurfTomo 相同，数据量大时快 | 否 |
| CC-FJpy（Li et al. 2021，LGPL-3） | 用本台阵的 F-J 谱得到第一遍参考曲线 | 借用参考曲线时不需要 |

### 1.3 检查安装

```bash
lant check
```

逐项列出 Python 包、GMT、DSurfTomo、fwd1d、CC-FJpy 的状态；缺什么会写出安装命令。

### 1.4 先跑通合成算例（推荐）

```bash
python examples/synthetic/make_synthetic.py       # 生成 16 台、3 天的合成数据（5–20 分钟，约 300 MB）
cd examples/synthetic
lant run all -c synthetic.yaml                    # 全流程约 15–30 分钟（其中反演和检验约 15 分钟）
python compare_truth.py                           # 与已知答案比较，写 output/figures/truth_comparison.png
```

合成数据由一个已知的一维模型（`true_model.txt`）的 Rayleigh 波从各方向随机入射叠加而成；数据故意用 miniSEED 和 SAC 混放、分在不同文件夹、SAC 没有台网代码，台站表 `stations.txt` 没有表头且经度在前，用来演示自动识别。测得的相速度应接近 `true_dispersion.txt`，反演结果应接近真模型（在分辨能力以内）。`compare_truth.py` 打印每个频率测得的中位数与真值之差、每个深度节点反演结果与真值之差。本仓库发布前的检验：相速度中位数在 0.9–4.8 Hz 与真值相差 −0.85 % 到 0（低频端略慢，台距只有 2–3 个波长时远场近似的偏差），第二遍参考曲线 0.9–9 Hz 相差在 0.8 % 以内，0–1.2 km 各深度节点的 Vs 与真值相差在 2 % 以内；16 个台的棋盘恢复有限（3 km 格子斜率 0.03–0.40），起始模型整体乘 1.05 的对照与主结果相差 0.14 km/s，说明台站这么少时结果明显依赖起始模型。

跑完看 `output/figures/` 和 `output/records/`，熟悉每一步的输出后再换成自己的数据。

---

## 2. 准备数据

**只需要两样东西：原始数据所在的目录，和一个台站坐标文件。** 不需要整理目录、改文件名、转格式、切天或合并，这些都由程序识别和处理。

### 2.1 原始记录

放在一个目录下即可，可以有任意层子目录、任意文件名，格式可以混用：miniSEED、SAC、GSE2、SEG-2、SEISAN 等 ObsPy 能读的格式都行（逐个文件自动识别；照片、表格、说明文档之类的非数据文件会被跳过并列在报告里）。例如下面两种组织都可以：

```
raw/                               raw/
├── 2024-06-01/                    ├── A01_whatever.sac
│   ├── YY.A01..HHZ.mseed          ├── A02_day1.sac
│   └── YY.A02..HHZ.mseed          ├── station_A03/
└── 2024-06-02/                    │   ├── part1.miniseed
    └── ...                        │   └── part2.miniseed
                                   └── notes.txt
```

程序从每条记录的记录头读取台网、台站、通道、采样率和起止时间，然后自动完成：

| 识别内容 | 做法 |
|---|---|
| 垂直分量 | 通道代码以 Z 结尾，或 SAC 头的 cmpinc = 0；某台只有一个通道时直接用它（报告里注明）；三分量都在时只用垂直分量 |
| 一台多种传感器或采样率 | 只用数据最多的一种，其余列出 |
| 台站名与台站表对应 | 见 2.2 |
| 处理日的起点 | 大多数长记录从哪个 UTC 整点开始（野外按北京时间切文件时是 16 时） |
| 降采样率 | 50、40、25、20、10、5 Hz 中能整除所有原始采样率的最大者 |
| 重复、重叠、采样网格偏移、缺数 | 清洗步骤处理并计数（步骤 clean） |

不需要事先去仪器响应、合并、切天或重采样。不去仪器响应的原因：互相干是谱归一的，同型号仪器的响应完全抵消；混用不同类型仪器时见第 7 节。

### 2.2 台站表

常见的写法都能识别：

- CSV、制表符或空格分隔的文本，有没有表头都可以；Excel（.xlsx 第一个工作表）；StationXML。
- 列名中英文都认：台网 / network / net，台站 / 台站名 / 站名 / 点号 / station / sta / name，纬度 / latitude / lat / stla，经度 / longitude / lon / stlo，高程 / 海拔 / elevation / elev / stel 等。
- 没有表头时按数值判断：不是数字的那一列是台站名，±90 以内的是纬度、另一个是经度（两列都在 ±90 以内时分不出，程序会要求加表头）。
- 台站名可以分台网和台站两列，也可以写在一起：`YY1001`、`YY.1001`、`YY_1001`，或只写台站号 `1001`。
- 坐标用十进制度或度分秒（`35°07'24.2"`、`35:07:24.2`），可以带 N/S/E/W。高程以米计，可以没有（只用于剖面图）。

台站表里的名字与记录头里的名字自动对应：台网和台站都相同；台站号相同（台站表或记录头没写台网，或台网写法不同）；台站表的名字是“台网 + 台站”连写；台站表的名字去掉前面的字母等于记录头的台站号。有歧义（一个名字对上多个台）的不用并在报告里列出。对应结果写成标准台站表 `01_clean/stations.csv`，此后所有输出里台站名都是 `NET.STA`。

**坐标要逐台核对**：坐标错一个台，它所有台对的走时都会错，且不会报错。SAC 头里有坐标时，程序会把与台站表相差 50 m 以上的台列出来，并写出 `01_clean/stations_from_headers.csv`。

### 2.3 磁盘和计算量

- 清洗后的数据档案约与原始数据一样大（`clean.write_archive: false` 可不写，后续直接读原始文件）。
- `02_ccf` 的逐小时互相关：每天约 `台对数 × 0.2 MB`（160 台、24 天约 60 GB）。全部完成后可以删掉（台站日筛选、出“叠加时长”图要用它）。
- 反演的时间主要取决于网格大小和“台站 × 周期”数：例如 59×71×19 的网格、160 台、53 个周期，DSurfTomo 10 线程约 4–6 小时一个工程；棋盘和对照各是一个工程，可以并行（`inversion.njobs`）。

---

## 3. 写配置

```bash
lant init surface -o my_project.yaml
```

只需要填三个路径：

```yaml
project:
  output_dir: /data/me/my_array_output     # 结果目录（不要放在代码目录里）
data:
  raw_dir: /data/me/my_array_raw
  station_file: /data/me/stations.csv
```

其余参数默认是 `auto` 或通用的默认值。`auto` 的参数在 `lant run clean` 时按书面规则从数据推出，写进 `<输出目录>/records/auto_parameters.json`（每个值都附推导依据），之后每一步都用它们：

| 参数 | 规则 |
|---|---|
| 台阵类型 `project.array_type` | 台站分布主轴长宽比 > 4 为线状，否则面状 |
| 处理日起点 `data.day_start_utc_hour` | 多数长记录开始的 UTC 整点 |
| 降采样率 `correlation.fs`、互相关频带上限 | 能整除所有采样率的最大值；上限 min(12.5, 0.4 fs) Hz |
| 频散频带下限 `dispersion.fmin_hz`、导出最长周期 | 90 % 的台距 ≥ 2 个波长（按 2.5 km/s），取到 1/12 倍频程网格 |
| 时延范围 `correlation.maxlag_s` | 最长台距 / 1.2 km/s + 2 个最长周期，取到 10 s |
| 网格间隔 `inversion.grid_spacing_deg` | 面阵：最近台站间距中位数的一半，取 0.0025° 的倍数；线状：0.005° |
| 棋盘格子 `checks.checkerboards`、`checks.vertical` | 台阵孔径的 0.12 / 0.18 / 0.25 倍；垂向变号深度按可反演深度的比例 |
| 剖面 `plots.sections` | 面阵：过台站中心、沿两个主轴方向 |

规则依据见 [parameters.md](parameters.md)。想改哪一个，在配置里直接写值，写了的值优先于规则。其他常改的参数（剔除台站、借用参考曲线、并行设置、地形底图）在 `config/surface.example.yaml` 的“optional”部分，含义见 [parameters.md](parameters.md)。

检查：

```bash
lant validate -c my_project.yaml          # 配置能否读、台站表怎样被识别、原始目录里有多少候选文件
lant run clean -c my_project.yaml
lant inspect auto -c my_project.yaml      # 从数据推出的参数和依据
lant inspect data -c my_project.yaml      # 台站几何、台间距分布
```

---

## 4. 逐步运行

一次跑完：`lant run all -c my_project.yaml`。建议第一次逐步跑，每步看完结果再往下：

```bash
lant run clean -c my_project.yaml
lant run cc -c my_project.yaml
...
lant run all -c my_project.yaml --from stack --to disp2     # 也可以跑一段
lant status -c my_project.yaml                               # 每一步的状态
lant inspect <步骤名> -c my_project.yaml                    # 每一步的关键数字和要看的图
```

每一步只读前面步骤的输出，所以改了某一步的参数后从那一步重跑即可。每次运行的配置、日志和摘要都记在 `<输出目录>/records/`（第 6 节）。

下表中的路径都相对于输出目录。

### 步骤 clean：原始数据清洗（`01_clean`）

**做什么**：扫描所有文件的记录头，建索引；认出每台的垂直分量；读台站表并与记录头的名字对应（写出 `stations.csv`）；从数据推出 `auto` 参数（`records/auto_parameters.json`）；每台定一种通道和采样率；按处理日把记录放到 UTC 采样网格上、合并、去掉完全相同的重复，内容不同的重叠段两份都不用；写出统一的日文件档案 `archive/NET.STA/NET.STA.LOC.CHA.YYYYDDD.mseed`。不滤波、不插补缺数、不改样点值。

**为什么**：后面的步骤只需面对一种格式、一种组织方式；重复、重叠、采样网格偏移这类问题在这里一次处理并记录，而不是散在各步。

**输出**：`report.txt`（文件数、格式、跳过的文件、用了哪个通道、台站表怎样被识别、名字怎样对应、对不上的台站、坐标不一致、几种仪器的台、推出的参数）、`station_days.csv`（每个台站日的小时数、缺口、平移和插值次数、重叠）、`coverage.png`（台站 × 天的数据覆盖）。

**怎样算正常**：`report.txt` 没有 ERROR；台站表里的台都对上了数据（对不上的会列出）；推出的参数合理（`lant inspect auto`）；`coverage.png` 大部分格子接近 24 小时；`interpolated` 很少（有很多说明记录时间不在采样网格上，可能是授时问题，应查原因）；`conflicting overlaps` 为 0 或很少。

**常见问题**：
- “no station of the table could be matched”：台站表和记录头的名字差得太远，看报告里举的例子，用 `data.rename` 指定对应关系（例如 `{"XX.001": "YY.A01"}`）。
- “ambiguous names”：一个名字对上了多个台（例如两个台网有同号台站），给台站表加上台网列。
- 某台“none recognisable as vertical”：通道名看不出哪个是垂直分量，用 `data.channel` 指定（如 `"*3"`）。
- 读不了的文件很多：多半不是地震数据（照片、表格），在 `ignore_patterns` 里排除即可；真是数据却读不了时指定 `data.format`。

### 步骤 cc：逐小时互相干（`02_ccf`）

**做什么**：每个台站日降采样到 50 Hz（零相位抗混叠），去掉 1 s 授时脉冲干扰；60 s 窗、50 % 重叠，含缺数或峰值超过当日标准差 10 倍的窗不用；每窗去线性趋势、5 % 尖灭、逐频点谱归一（互相干）；每小时平均，存 ±20 s 时延。一天一个文件 `day_YYYYDDD.h5`。

**为什么**：谱归一的互相干对高频、强工业干扰的密集台阵最稳（Xu et al. 2021）；逐小时保存是为了后面能按天、按小时筛选和估计误差，而不必重算。

**判断**：屏幕上每天一行 `ok N stations M pairs`。日文件里 `segments_used / segments_gapfree` 是窗的保留比例，正常应在 95 % 以上；很低说明 `max_over_std` 对这批数据太严或有大量瞬态干扰。

**常见问题**：内存不够（进程被杀）→ 减小 `correlation.nproc`。改了互相关参数后重跑会报 “STALE … written with other settings”：删掉 `02_ccf`（或换输出目录），旧文件不会被悄悄沿用。

### 步骤 stack：叠加和筛选（`03_stack`）

**做什么**：先叠加全部小时（`stack_raw.h5`）；面阵按天比较每个台站日与该台总叠加的相关和时移，剔除异常台站日（`exclude_station_days.txt`）；线状按小时比较（`station_hour_qc/`）；得到干净的叠加 `stack_clean.h5`（含 24 h 或 12 h 子叠加，供刀切）。

**为什么**：钟差、被挖出或移动的仪器、耦合变差的时段会给所有台对带来系统误差，平均之前要去掉；规则只看与自身的一致性，不看结果好不好。

**判断**：看 `station_day_qc/station_day_qc.png` 和 `.txt`。左图是“相关系数减当日中位数”，右图是共同时移；红点是被剔的。正常时被剔的是少数孤立点（一般只占片段的百分之一二）。某台连续多天都被剔，或整体时移 ≥ 0.04 s，多半是仪器问题：考虑把它整台放进 `data.exclude_stations`，并查野外记录。

### 步骤 notch：陷波（`04_notch`）

**做什么**：在原始波形和叠加谱上按规则找持续的窄带谱线（工业干扰），生成陷波表并对叠加做陷波，得到最终互相关成品 `stack.h5`。

**判断**：`notch_report.png` 上每个尖峰都被阴影盖住、宽缓的鼓包没有被陷；`raw_lines.txt` 列出原始波形上的谱线和占比。谱线很多（十几条以上）时对照野外环境（水泵、电机、高压线）核实。

### 步骤 reference：第一遍参考曲线（`05_reference`）

**做什么**：`source: fj` 时用全部台对做 F-J 谱（CC-FJpy），取基阶脊线并平滑成第一遍参考曲线；`source: file` 时读入借用的曲线。

**判断**：看 `fj.png`：白线应沿着最强、最连续的能量脊（基阶），白圈是采用的拾取点。图上若有两条脊（高阶模式）而白线跳到了另一条，或低频端没有能量，说明需要换 `fj_mode`、调 `track_start_band_hz`，或借用参考曲线。参考曲线只是向导，第二遍会换成实测中位数；但它定的周数错了，整条曲线会差一个周期，所以必须看图。

### 步骤 disp1：第一遍测量和第二遍参考（`06_dispersion`、`05_reference`）

**做什么**：用第一遍参考测所有台对（不做偏离参考的筛选），取每个频率被接受的相速度的中位数作第二遍参考（`reference.npz`、`reference_pass2.png`）。借用参考曲线时跳过。

**判断**：`reference_pass2.png` 左图：红线（第二遍）应在灰带（第一遍 16–84 %）中间、平滑；与蓝虚线（F-J）相差一般在几个百分点以内。右图是每个频率的台对数，低于门限的频率没有参考、第二遍不测。中位数与 F-J 差 5 % 以上时先查 F-J 拾取。

### 步骤 disp2：正式测量、质控和导出（`06_dispersion`）

**做什么**：用第二遍参考正式测量；相速度质控（信噪比、波长、旁瓣、偏离参考、刀切、平滑、相群一致、连续性）；群速度另一套质控；冗余路径精度检验；按测量频点的周期导出，每条相速度曲线只留最长的一段（`curves/`，反演读这里）；抽样画单台对图（`pair_figures/`）。

**判断**：
1. `pair_figures/`：随便翻二三十张，黑点应落在相速度图和群速度图的能量脊上，不应跳到旁边的脊；
2. `lant inspect disp2`：每个频率的点数，应该在中间频率最多、两端逐渐减少，没有某一个频率突然掉到零（陷波频点附近少一些是正常的）；
3. `precision.txt`：冗余路径的相对差，代表测量精度（一般为百分之一到几）；
4. `figures/3_dispersion_curves.png`（出图后）：曲线集中在平均线附近。

**常见问题**：点很少 → 先看 `pre_c`（`inspect` 的 `candidate_points_before_phase_qc`）：候选点就少，问题在信噪比或台距（数据太短、台阵太小）；候选多而最终少，看哪条质控挡得最多。**不要为了多出点而放宽门限**：质控的目的是只留可靠的点，宁缺毋滥。

### 步骤 inversion：建反演工程（`07_inversion`）

**做什么**：每个周期取台对中位数得平均曲线（`average_curve.txt`）；按“最长波长的一半”定可反演的深度节点和固定底节点；由平均曲线反演一维起始模型（`vs1d.txt`、`vs1d_fit.txt`）；写主工程 `main/`（`surfdata.dat`、`DSurfTomo.in`、`MOD`），以及棋盘 `checkerboards/` 和对照 `controls/`。

**判断**：`figures/4_start_model.png`（或 `vs1d_fit.txt`）：预测曲线应贴着平均曲线（rms 通常 < 0.02 km/s）；起始模型随深度变化平缓。`main/DSurfTomo.in` 的网格应覆盖全部台站、周期数 ≤ 80。

### 步骤 run：运行 DSurfTomo（`07_inversion/*/`）

```bash
lant run run -c my_project.yaml --what main      # 先只跑主工程
lant run run -c my_project.yaml --what checks    # 再跑棋盘和对照
```

每个工程报告 ok / skipped（已算完且输入未变） / FAILED。中断过的工程会报 FAILED 而不被覆盖：确认后删掉该工程目录里的输出（或整个目录后重跑 inversion）再运行。

### 步骤 evaluate：检验（`07_inversion/evaluation`）

**做什么**：用 DSurfTomo 自己的正演对最终模型和每一轮的模型计算全部走时的残差（`residuals_final.csv`、`evaluation.json`）；棋盘按台站 0.8 km 以内的节点统计每个深度恢复的相关系数和斜率；对照工程与主工程的差值。只报数和出图，不设通过标准。

**判断**：
- 残差：从起始模型到最终模型应明显下降，最后几轮变化很小（`figures/6_residuals.png` 右图）；残差直方图中心在 0 附近。最终残差的大小应与数据误差同量级；比数据误差小很多说明可能过拟合（正则化太弱）。
- 棋盘：`figures/6_checkerboard_*` 和 `6_recovery.png`。斜率表示恢复了多少幅度（1 = 完全恢复），相关系数表示图样是否正确。随深度变差是正常的；只有恢复得好的尺度和深度上的结构才能解释。
- 对照：正则化、起始模型变化引起的模型差（rms 和最大值）。差值大的区域，结构依赖于设置，解释要谨慎。

**不要做的**：不要用检验结果反过来挑参数直到“好看”；正则化按 parameters.md 的常规范围选定后，检验只用来说明结果能信到什么程度。

### 步骤 figures：出图（`figures`）

```bash
lant plot -c my_project.yaml                       # 全部
lant plot depth_slices sections -c my_project.yaml # 只画几张
```

| 文件 | 内容 |
|---|---|
| `1_stations` | 台站分布（有地形时叠晕渲地形） |
| `2_ccf_moveout` | 全部台对的互相关随台距排列，1–5 Hz，两条参考速度线 |
| `2_ccf_stack_time` | 一个台对：叠加波形、逐小时互相关、信噪比随叠加时长 |
| `3_dispersion_pair` | 一个台对的波形、相速度图、群速度图和拾取 |
| `3_dispersion_curves` | 全部相速度曲线、平均和标准差 |
| `4_start_model` | 平均曲线、一维模型的拟合和起始模型 |
| `4_kernels` | 起始模型上几个周期的相速度深度敏感度（DSurfTomo 正演链有限差分；数据在 `07_inversion/kernels.txt`） |
| `5_depth_slices` | Vs 深度切片（浅灰地形上半透明，台站凸包内） |
| `5_sections` | Vs 剖面（面阵按 `plots.sections`；线状沿测线） |
| `6_checkerboard_*`、`6_recovery`、`6_residuals` | 棋盘输入与恢复、恢复随深度、残差 |

图中文字为英文，格式和分辨率由 `plots.formats`、`plots.dpi` 控制，`plots.font_scale` 放大字号（幻灯片用 1.5–2）。色标范围 `plots.vs_range` 对所有图相同，便于比较。地图和剖面用 GMT 画：地形用 `plots.topography`（GMT 网格文件，或 `@earth_relief_01s` 之类的 GMT 远程数据，需要联网下载）。图上不叠加任何地质信息；需要叠加断裂、矿区等时，用输出的模型文件 `07_inversion/main/DSurfTomo.inMeasure.dat`（经度 纬度 深度 Vs）自己画。

---

## 5. 线状台阵

一条测线一个项目（`lant init line`）。与面阵的不同在 parameters.md 第 12 节。要点：

1. **参考曲线和起始模型**。测线的 F-J 频带窄、常被高阶模式或斜入射的噪声主导。同一地区有面阵时，借用面阵的结果：
   ```yaml
   reference: {source: file, file: /path/to/surface_output/05_reference/reference.npz}
   inversion: {start_model: {source: file, file: /path/to/surface_output/07_inversion/vs1d.txt}}
   ```
   没有面阵时用本线的 F-J（`fj_mode: track_fundamental`），务必看 `05_reference/fj.png`。借用时，线状结果与面阵相符不能算独立验证。
2. **起始模型与数据整体偏离时**（借来的模型比本线整体快或慢），DSurfTomo 的残差剔除第一轮会单侧剔掉大量数据，看 `evaluate` 的 `rejected_by_rule`。
3. **分辨**：只沿测线方向有分辨；测线两端台站稀，剖面两端不可靠。噪声斜入射时沿线视速度偏高，高频群速度常落在搜索上限。
4. 同一批仪器在不同测线复用时，每条线单独一个项目、单独一个原始目录，避免同名台站混在一起。

---

## 6. 记录与复现

每个输出目录的 `records/` 里：

| 文件 | 内容 |
|---|---|
| `config.yaml` | 最近一次运行用的完整配置（含全部默认值） |
| `config_history/` | 配置每次改变前的版本 |
| `run_log.jsonl` | 每一步的开始、结束、失败、耗时 |
| `summary/<步骤>.json` | 每一步的关键数字（`lant inspect` 读它） |
| `decisions.md` | 参数选择及其依据（人工或 AI 填写） |

复现：`lant run all -c <输出目录>/records/config.yaml`。所有记录都在输出目录里，代码目录不会被写入。

---

## 7. 常见问题

**只想换一个参数重算后面的步骤？** 改配置后从受影响的那一步重跑（例如改质控门限：`lant run disp2 inversion run evaluate figures`）。改了互相关参数要删掉 `02_ccf`。

**测量阶段很慢？** 增加 `dispersion.nproc`。刀切每组都要重算图像，组数越多越慢。

**混用了不同类型的仪器？** 互相干消掉了振幅响应；同型号仪器的相位响应也完全抵消。不同型号混用时，混合台对有相位偏差，在低频（1 Hz 以下）可达百分之一的量级，本程序不改正，解释低频结果时要考虑。

**F-J 装不上？** 用 `reference.source: file` 借一条参考曲线（例如同一地区已有研究的平均频散、或 `true_dispersion.txt` 那样的文本），其余步骤不变。

**DSurfTomo 报 segmentation fault？** 多半是内存或网格：台站贴近网格边缘、`sparsity` 太小、走时条数很大。看 `07_inversion/main/run.log` 末尾；网格过大时增大 `grid_spacing_deg` 或减少 `grid_refine`。

**模型边缘出现很强的异常？** 网格边界和固定底节点不更新，台站覆盖外的节点没有约束。切片只在台站凸包内着色；解释前看棋盘。

**结果可以直接用于地质解释吗？** 先看 evaluate：只有棋盘恢复得好、对照差值小的尺度和深度上的结构才可解释；Vp、密度由 Vs 按 Brocher (2005) 经验关系换算（DSurfTomo 内部写死），真实 Vp/Vs 偏离时 Vs 有百分之几的系统偏差。
