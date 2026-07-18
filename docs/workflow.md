# 背景噪声成像完整流程(01–06)

本项目实现从原始 miniSEED 数据到横波速度模型的完整背景噪声成像流程。公共流程为 01–05 五步;第 05 步之后按台阵类型二选一进入 **06a 线状反演** 或 **06b 面状反演** 分支。

所有命令都在**项目根目录**执行;每步独立运行,可从任意中间步骤继续。

## 流程概览

```
原始 miniSEED 数据
  │
  ▼ [01] scripts/01_clean.py                数据清洗(NoisePy 命名格式)
  ▼ [02] scripts/02_correlation.py          互相关 + 叠加(NoisePy)
  │        前置:准备台站文件 station.csv
  ▼ [03] scripts/03_stack_qc.py             叠加质检可视化(可选但推荐)
  │      scripts/03b_stack_panels.py        出版级多面板图
  ▼ [04] scripts/04_export_egf.py           导出 EGFAnalysisPy 输入(.dat)
  ▼ [05] scripts/05_dispersion.py           频散分析 + DisperPicker 拾取
  │                                         → qc_result/CDisp.*.txt、GDisp.*.txt
  │
  ├─▶ [06a] scripts/06a_line_inversion/     线状台阵:evodcinv 一维反演 → 伪二维剖面
  │         (详见 docs/line_inversion.md)
  └─▶ [06b] scripts/06b_area_inversion/     面状台阵:DSurfTomo 三维反演
            (详见 docs/area_inversion.md)
```

### 如何选择反演分支

| | 06a 线状反演 | 06b 面状反演 |
|---|---|---|
| 台阵形态 | 台站近似排成一条测线 | 台站覆盖一个二维区域 |
| 方法 | evodcinv:逐台站对一维 Vs(z) 反演(CPSO),按中点拼接 | DSurfTomo:频散直接反演三维 Vs |
| 输入 | 05 步 `qc_result/`(直接使用) | 05 步 `qc_result/`(经 06b-1 转换为 surfdata.dat) |
| 结果 | 伪二维 Vs 深度剖面(section.png/.npz) | 三维 Vs 模型 + 深度切片图 |
| 附加依赖 | evodcinv(装入 noisepy 环境,INSTALL §2e) | DSurfTomo 编译(INSTALL §4)、GMT 6 |

两个分支相互独立,互不依赖,也可以对同一批数据都跑。

### 各步骤环境与输入输出

| 步骤 | 用户入口 | 环境 | 主要输入 | 主要输出 |
|------|----------|------|----------|----------|
| 01 | `scripts/01_clean.py` | noisepy | 原始 miniSEED | 清洗后 `.ms` |
| 02 | `scripts/02_correlation.py` | noisepy | 清洗数据 + station.csv | CCF + Stack(HDF5) |
| 03 | `scripts/03_stack_qc.py`、`03b_stack_panels.py` | noisepy | Stack | moveout 图 |
| 04 | `scripts/04_export_egf.py` | noisepy | Stack | CF `.dat` + stations.txt |
| 05 | `scripts/05_dispersion.py` | egfanalysispy | CF `.dat` | 频散曲线 `qc_result/` |
| 06a | `scripts/06a_line_inversion/01_invert.py`、`02_stitch.py` | noisepy + evodcinv | `qc_result/` | 一维模型集 + 伪二维剖面 |
| 06b | `scripts/06b_area_inversion/01_prepare.py`、`02_invert.py`、`03_plot.py` | noisepy + DSurfTomo 可执行文件 + GMT | `qc_result/` | 三维模型 + 切片图 |

### 旧脚本名 → 新入口对照(迁移指南)

| 旧名称 | 新名称 |
|--------|--------|
| `scripts/clean_mseed_for_noisepy.py` | `scripts/01_clean.py` |
| `scripts/local_cc_workflow.py` | `scripts/02_correlation.py` |
| `scripts/stack_visualization.py` | `scripts/03_stack_qc.py` |
| `scripts/plot_stack_panels.py` | `scripts/03b_stack_panels.py` |
| `scripts/export_egfanalysis_input.py` | `scripts/04_export_egf.py` |
| 手改 `EGFAnalysisPy/config.py` + `python run.py` | `scripts/05_dispersion.py --config configs/05_dispersion.yaml` |
| `lant_invert.py` / `lant_stitch.py`(本地项目版本) | `scripts/06a_line_inversion/01_invert.py` / `02_stitch.py` |
| `configs/local_cc_workflow.example.yaml` 等 | `configs/02_correlation.example.yaml` 等(编号对应) |

配置行为不变:旧配置文件内容可直接复制到新的同义文件中使用。

## 软件包引用

本项目使用了以下开源软件包,感谢原作者的贡献:

- **NoisePy** — 背景噪声互相关计算
  Jiang, C., & Denolle, M. A. (2020). NoisePy: A New High-Performance Python Tool for Ambient-Noise Seismology. *Seismological Research Letters*, 91(3), 1853–1866. https://doi.org/10.1785/0220190364
  GitHub: https://github.com/noisepy/NoisePy

- **EGFAnalysisPy** — 面波频散分析(EGFAnalysisTimeFreq 的纯 Python 重写,内置 DisperPicker)
  Yao, H., van der Hilst, R. D., & de Hoop, M. V. (2006). Surface-wave array tomography in SE Tibet from ambient seismic noise and two-station analysis—I. Phase velocity maps. *Geophysical Journal International*, 166(2), 732–744. https://doi.org/10.1111/j.1365-246X.2006.03028.x
  作者:姚华建(中国科学技术大学)

- **DisperPicker** — 基于深度学习的频散曲线自动提取(已内置于 EGFAnalysisPy)
  Yang, S., Zhang, H., Gu, N., Gao, J., Xu, J., Jin, J., Li, J., & Yao, H. (2022). Automatically Extracting Surface Wave Group and Phase Velocity Dispersion Curves from Dispersion Spectrograms Using a Convolutional Neural Network. *Seismological Research Letters*. https://doi.org/10.1785/0220210280
  GitHub: https://github.com/ShaoboYang-USTC/DisperPicker

- **DSurfTomo** — 面波直接反演三维横波速度结构(仅 06b 分支使用)
  Fang, H., Yao, H., Zhang, H., Huang, Y. C., & van der Hilst, R. D. (2015). Direct inversion of surface wave dispersion for three-dimensional shallow crustal structure based on ray tracing: methodology and application. *Geophysical Journal International*, 201(3), 1251–1263. https://doi.org/10.1093/gji/ggv069

- **evodcinv / disba** — 频散曲线一维演化算法反演(仅 06a 分支使用)
  Luu, K. GitHub: https://github.com/keurfonluu/evodcinv

## 环境准备

详细步骤见 [INSTALL.md](INSTALL.md)。摘要:

```bash
conda env create -f environments/noisepy.yml        # 01–04、06a、06b 入口
conda env create -f environments/egfanalysispy.yml  # 05
cd DSurfTomo && python configure && cd ..           # 仅 06b 需要
conda activate noisepy && python tools/check_environment.py
```

预训练模型已包含在 `EGFAnalysisPy/DisperPicker/saver/`,无需下载。

---

## 01 数据清洗

**目的:** 剔除不完整的天数据,将原始 miniSEED 转换为 NoisePy 命名格式。

**入口:** `scripts/01_clean.py`(环境:noisepy)

**输入:** 原始 miniSEED 目录,按台站分子目录:

```
RawData/
├── CF.17982/
│   ├── 230901000000    # 文件名为北京时间戳 YYMMDDHHmmss
│   └── ...
├── CF.58808/
└── CN.20496/
```

**命令:**

```bash
# 用配置文件(推荐;先 cp configs/01_clean.example.yaml configs/01_clean.yaml 并编辑)
python scripts/01_clean.py --config configs/01_clean.yaml

# 或直接命令行(命令行参数覆盖配置文件)
python scripts/01_clean.py \
    --input-root  TestData/RawData \
    --output-root TestData/CleanData \
    --target-channel Z
```

**参数说明:**

| 参数 | 说明 |
|------|------|
| `--config` | YAML 配置文件,键名与下列参数一致 |
| `--input-root` | 原始 miniSEED 根目录 |
| `--output-root` | 清洗后输出目录 |
| `--min-completeness` | 单日数据最低完整度(0~1,默认 0.8) |
| `--target-channel` | 目标通道/分量,如 `Z`、`BNZ`(默认自动取 Z 分量) |
| `--date-range START END` | 仅处理日期严格在两者之间的文件(两端不含),格式 YYMMDD |
| `--overwrite` | 覆盖已有输出(默认不覆盖) |
| `--dry-run` | 仅打印操作,不写文件 |
| `--verbose` | 逐文件详细日志 |

**输出:** `CleanData/CF17982_BNZ00_2023258.ms ...`,命名 `{台网}{台站}_{通道}{位置码}_{年+儒略日}.ms`

**前置准备(02 步需要):台站坐标 CSV**

```csv
network,station,channel,latitude,longitude,elevation
CF,17982,BNZ,42.323349,118.546867,844
CF,58808,BNZ,42.321987,118.54599,822
```

台网/台站需与 miniSEED 头一致;经纬度十进制度;高程单位米。

---

## 02 互相关计算与叠加

**目的:** 用 NoisePy 计算台站对间背景噪声互相关函数(CCF)并叠加。

**入口:** `scripts/02_correlation.py`(环境:noisepy)

```bash
cp configs/02_correlation.example.yaml configs/my_02.yaml   # 编辑路径与参数
```

**关键配置参数:**

```yaml
data-dir: ../TestData/CleanData       # 01 步输出
station-csv: ../TestData/station.csv  # 台站文件
ccf-dir: ../TestData/CCFData          # 互相关输出
stack-dir: ../TestData/StackData      # 叠加输出

freqmin: 0.2          # 最低频率 (Hz)
freqmax: 10           # 最高频率 (Hz)
sampling-rate: 25.0   # 重采样率 (Hz)
maxlag: 100           # 最大滞后 (s)
inc-hours: 24         # 数据分段长度 (h)
cc-len: 600           # 互相关窗长 (s)
step: 300.0           # 滑动步长 (s)
time-norm: rma        # 时间域归一化
freq-norm: rma        # 频率域归一化
stack-method: all     # linear / pws / robust / all
```

**命令:**

```bash
python scripts/02_correlation.py --config configs/my_02.yaml --validate-only  # 只验证配置
python scripts/02_correlation.py --config configs/my_02.yaml --stage cc       # 只互相关
python scripts/02_correlation.py --config configs/my_02.yaml --stage stack    # 只叠加
python scripts/02_correlation.py --config configs/my_02.yaml                  # 全部
```

---

## 03 叠加质检可视化

**目的:** 检查叠加结果质量,选择后续使用的叠加方法。可选但强烈建议。

**入口:** `scripts/03_stack_qc.py`(环境:noisepy)

```bash
python scripts/03_stack_qc.py \
    --stack-dir TestData/StackData \
    --freqmin 1 --freqmax 5 \
    --component ZZ --stack-method linear --maxlag 10 --style both
# 或
python scripts/03_stack_qc.py --config configs/03_stack_qc.yaml
```

**输出:** moveout 图(默认写入 `--outdir`)。

**出版级多面板图(可选):** `scripts/03b_stack_panels.py` 生成 2×N 多面板 moveout 图(N=叠加方法数,上行 image、下行 wiggle/symmetric):

```bash
python scripts/03b_stack_panels.py --stack-dir TestData/StackData --methods linear pws robust --component ZZ
# 或 --config configs/03b_stack_panels.yaml
```

**叠加方法怎么选:** `linear` 最常用;`pws` 信噪比更高;`robust` 抗异常值;建议 02 步用 `all` 全算,再用本步对比选优。

---

## 04 导出 EGFAnalysisPy 所需数据

**目的:** 将 NoisePy 叠加结果(HDF5)转换为 EGFAnalysisPy 的文本格式(CF `.dat`)。

**入口:** `scripts/04_export_egf.py`(环境:noisepy)

```bash
python scripts/04_export_egf.py \
    --stack-dir  TestData/StackData \
    --stack-method all --component ZZ \
    --output-root TestData/EGFAnalysisInputData
# 或 --config configs/04_export_egf.yaml
```

**参数说明:**

| 参数 | 说明 |
|------|------|
| `--config` | YAML 配置文件 |
| `--stack-dir` | NoisePy 叠加结果目录 |
| `--stack-method` | `linear` / `pws` / `robust` / `all`(默认) |
| `--component` | 分量,默认 `ZZ` |
| `--station` | 可选,仅导出含指定台站的台站对 |
| `--output-root` | 输出根目录 |

**输出目录结构:**

```
EGFAnalysisInputData/
├── Allstack_linear/ZZ/
│   ├── CF17982_CF58808.dat    # 波形(三列:时间, gab, gba)
│   ├── stations.txt           # 台站坐标
│   └── EGFfiles.txt           # 波形文件列表
├── Allstack_pws/ZZ/
└── Allstack_robust/ZZ/
```

`.dat` 文件格式:第 1、2 行为源/接收台站的 `经度 纬度 高程`,第 3 行起为 `时间(s) gab gba`。

---

## 05 频散分析与曲线拾取

**目的:** 对 CF/EGF 做时频分析生成频散能量图,用内置 DisperPicker(CNN)自动拾取群速度、相速度频散曲线,经质量控制后输出。整个过程一个命令完成,无需 MATLAB、无需进入 EGFAnalysisPy 目录改配置。

**入口:** `scripts/05_dispersion.py`(环境:**egfanalysispy**)

```bash
cp configs/05_dispersion.example.yaml configs/my_05.yaml   # 编辑
conda activate egfanalysispy
python scripts/05_dispersion.py --config configs/my_05.yaml
```

**关键配置(`configs/05_dispersion.example.yaml`,参数含义与 EGFAnalysisPy/config.py 一致):**

```yaml
paths:
  waveform_path: /path/to/EGFAnalysisInputData/Allstack_robust/ZZ  # 04 步输出
  result_path:   /path/to/DisperPicker_result
run:
  parallel: false     # true 用多进程(run_parallel.py)
  workers: 4
dispersion:
  StartT: 0.4         # 周期范围与步长——06 步会用到,请记下
  EndT: 3.0
  DeltaT: 0.02
  StartV: 0.5
  EndV: 3.5
  DeltaV: 0.005       # DisperPicker 要求固定 0.005 km/s
```

入口脚本会把 YAML 渲染成运行时配置(`<result_path>/runtime_config.py`,含来源与时间戳,可追溯),再调用 `EGFAnalysisPy/run_parallel.py`;仓库内的 `EGFAnalysisPy/config.py` 不会被修改。内部流程:频散能量图 → DisperPicker 拾取 → 按群速度加窗时变滤波重算相速度图 → 再拾取 → 质量控制 → 输出。

**高级用法:** 也可以直接 `cd EGFAnalysisPy && python run_parallel.py --config 自定义config.py -w 8`(或编辑 `config.py` 后 `python run.py`),行为与原版一致。

**输出:**

```
<result_path>/
├── runtime_config.py          # 本次运行的完整参数(自动生成)
├── qc_result/
│   ├── CDisp.FD01_FD16.txt    # 相速度频散曲线
│   ├── GDisp.FD01_FD16.txt    # 群速度频散曲线
│   └── ...
├── qc_plot/                   # 能量图 + 拾取对比图(按 random_plot_ratio 抽样)
├── all_C.png / mean_C.png
└── all_G.png / mean_G.png
```

频散曲线文件格式(前两行为两台站的 `经度 纬度`,数据列为 `周期 速度 信噪比 概率`):

```
117.43333100    31.81069500
117.48089000    31.66345700
0.20   1.250   6.100   0.950
0.30   1.310   5.800   0.880
```

> `qc_result/` 就是两个反演分支的公共输入:06a 直接读取;06b 由 `01_prepare.py` 转换。台站坐标内嵌在文件头,无需额外 station 文件。

---

## 06a 线状台阵反演(evodcinv)

台站近似排成一条测线时使用。两步:逐台站对一维反演 → 拼接伪二维剖面。

```bash
cp configs/06a_line_inversion/line.example.yaml configs/06a_line_inversion/my_line.yaml
conda activate noisepy    # 需已安装 evodcinv(INSTALL §2e)
python scripts/06a_line_inversion/01_invert.py --config configs/06a_line_inversion/my_line.yaml
python scripts/06a_line_inversion/02_stitch.py --config configs/06a_line_inversion/my_line.yaml
```

详细说明(适用条件、配置逐项解释、结果解读)见 [line_inversion.md](line_inversion.md)。

## 06b 面状台阵反演(DSurfTomo)

台站覆盖二维区域时使用。三步:准备工程 → 运行反演 → GMT 绘图。

```bash
cp configs/06b_area_inversion/dsurftomo.example.yaml configs/06b_area_inversion/my_area.yaml
conda activate noisepy    # 需已编译 DSurfTomo(INSTALL §4);绘图需 GMT 6
python scripts/06b_area_inversion/01_prepare.py --config configs/06b_area_inversion/my_area.yaml
python scripts/06b_area_inversion/02_invert.py  --config configs/06b_area_inversion/my_area.yaml
python scripts/06b_area_inversion/03_plot.py    --config configs/06b_area_inversion/my_area.yaml
```

也可以使用 `DSurfTomo/scripts/` 下的交互式脚本手动分步执行(`PrepareDSurfTomoProject.py` 一体式,或 `extractnew.py` → `GenerateIniMOD.py` → `GenerateDSurfTomoInputFile.py` 分步),详见 [area_inversion.md](area_inversion.md)。

---

## 常见问题

常见问题已集中到 [FAQ.md](FAQ.md)。

## 许可证

各软件包遵循其原始许可证。辅助脚本(`scripts/` 和 `DSurfTomo/scripts/`)随项目发布。
