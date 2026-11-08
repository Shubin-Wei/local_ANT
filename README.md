# local_ANT

**密集台阵高频背景噪声面波成像：从任意组织、任意格式的原始记录到三维横波速度模型**
**High-frequency ambient-noise surface-wave tomography for dense (nodal) arrays: from raw records of any layout and format to a 3-D shear-wave velocity model**

[中文](#中文) | [English](#english)

---

## 中文

### 能做什么

对密集台阵（节点仪或宽频带，台间距几百米到十几千米）的垂直分量连续记录，完成 0.5–10 Hz（0.1–2 s）Rayleigh 波相速度测量和浅部（约 0–3 km）三维 Vs 反演。面状台阵直接三维反演，线状台阵沿测线反演。方法和参数来自密集节点台阵（面阵和测线）的实际研究，每个处理步骤都有文献或开源程序依据，没有文献取值的参数在文档里如实标出。

| 步骤 | 内容 |
|---|---|
| clean | 原始数据识别与清洗：任意目录结构和文件名，miniSEED、SAC、GSE2 等 ObsPy 能读的格式；记录头读台站和时间，自动认垂直分量、识别台站表并对应名字、由数据推出参数；统一采样网格、去重复、处理重叠，写出统一的日文件 |
| cc | 逐小时互相干（1 分钟窗、谱归一），去 1 s 授时脉冲干扰、剔除瞬态窗 |
| stack | 叠加；面阵按台站日、线状按台站小时剔除异常时段 |
| notch | 按规则找持续工业谱线并陷波 |
| reference | 第一遍参考曲线（本台阵 F-J 谱，或借用其他台阵） |
| disp1 / disp2 | 两遍两台法频散测量：相速度、群速度，刀切误差，逐频点质控，按实测周期导出、每对只留一段 |
| inversion / run | 平均曲线、深度节点、一维起始模型，DSurfTomo（或并行版 pDSurfTomo）三维反演 |
| evaluate | 正演残差、水平和垂向棋盘、正则化和起始模型对照 |
| figures | 台站图、互相关、频散、起始模型、深度切片、剖面、棋盘、残差（地图和剖面用 GMT） |

另有一套 **AI 协作流程**（[docs/ai_workflow.md](docs/ai_workflow.md)、[AGENTS.md](AGENTS.md)）：在仓库里启动 Claude Code、Codex 等 AI 工具，把数据位置告诉它，它会按规则检查数据、由数据推出参数、逐步运行和检查，必须由你决定的事先问你，全部记录留在你的输出目录。

### 安装

```bash
git clone https://github.com/Shubin-Wei/local_ANT.git && cd local_ANT
conda env create -f environment.yml && conda activate lant
pip install -e .
bash tools/install_external.sh       # 下载并编译 DSurfTomo、pDSurfTomo、CC-FJpy（固定版本）
lant check                           # 检查安装
```

### 快速开始

```bash
# 1. 先跑合成算例（已知答案）
python examples/synthetic/make_synthetic.py
cd examples/synthetic && lant run all -c synthetic.yaml && cd ../..

# 2. 用自己的数据：只填原始数据目录、台站文件、输出目录三个路径
lant init surface -o my_project.yaml
lant validate -c my_project.yaml
lant run all -c my_project.yaml            # 或逐步：lant run clean -c ...，每步后 lant inspect <步骤>
```

**只需要指出原始数据目录和台站文件**，其余由程序识别：原始数据的目录结构、文件名、格式（miniSEED、SAC 等可混用）、垂直分量、台站名与台站表的对应、采样率、野外切文件的时区都从数据里读出；台站表的写法（有无表头、中英文列名、分隔符、Excel、度分秒）自动识别；台阵类型（面状 / 线状）、频带、时延、网格、棋盘等按书面规则由数据推出，每个值和依据记在输出目录的 `records/auto_parameters.json`，写进配置的值优先。结果、图和全部记录都在 `output_dir` 里。

### 文档

| 文档 | 内容 |
|---|---|
| [docs/tutorial.md](docs/tutorial.md) | 详细教程：安装、数据准备、配置、每一步做什么、输出、怎样判断正常、绘图、线状台阵、常见问题 |
| [docs/parameters.md](docs/parameters.md) | 全部参数：含义、单位、默认值、依据、什么时候改、设置不当的后果 |
| [docs/ai_workflow.md](docs/ai_workflow.md) | 让 AI 协助完成全流程：怎样开始、AI 能自己定什么、必须问什么、各阶段判据 |
| [docs/history.md](docs/history.md) | 本仓库的整理记录 |

### 引用

方法与参数见 [docs/parameters.md](docs/parameters.md) 的参考文献。使用的开源程序：DSurfTomo（Fang et al. 2015, GJI 201:1251）、pDSurfTomo（Zhu et al. 2026, arXiv:2604.11920）、CC-FJpy（Li et al. 2021, SRL 92:3179）、ObsPy、GMT。

### 许可

MIT（见 [LICENSE](LICENSE)）。第三方程序由安装脚本从各自的官方仓库下载，保留各自的许可（DSurfTomo MIT、pDSurfTomo GPL-3、CC-FJpy LGPL-3），不包含在本仓库中。

---

## English

### What it does

Rayleigh-wave phase velocities at 0.5–10 Hz (0.1–2 s) and shallow (about 0–3 km) 3-D Vs from continuous vertical-component records of dense arrays (nodal or broadband, station spacing from a few hundred metres to a few kilometres): direct 3-D inversion for areal arrays, narrow-strip inversion along linear arrays. The method and its parameters come from a study of a dense nodal array (an areal array and lines); every step is based on published methods or open-source programs, and parameters without a published value are marked as such in the documentation.

You only point to the raw-data directory and the station file: file layout and names, formats (miniSEED, SAC, ... mixed), vertical component, station-name matching, sampling rate, field time zone, the form of the station table, array type, frequency band, lag window, grid and checkerboard sizes are recognised or derived from the data by written rules (recorded with their reasons in the output directory).

Steps: **clean** (raw data of any layout and any ObsPy-readable format → recognition, name matching, derived parameters, a uniform day archive), **cc** (hourly cross-coherence), **stack** (station-day / station-hour screening), **notch** (spectral lines by rule), **reference** (F-J curve of the array, or a borrowed one), **disp1/disp2** (two-pass two-station phase and group velocities with jackknife errors and per-frequency QC, exported on the measured periods, one segment per pair), **inversion/run** (average curve, depth nodes, 1-D starting model, DSurfTomo or pDSurfTomo), **evaluate** (forward residuals, horizontal and vertical checkerboards, regularisation and starting-model controls), **figures** (GMT maps and sections, matplotlib figures).

An **AI-assisted workflow** ([docs/ai_workflow.md](docs/ai_workflow.md), [AGENTS.md](AGENTS.md)) lets an AI coding assistant (Claude Code, Codex, ...) run the whole chain on your data: it derives parameters from the data by written rules, checks every step, asks you where a decision is yours, and keeps every record in your output directory.

### Install and quick start

```bash
git clone https://github.com/Shubin-Wei/local_ANT.git && cd local_ANT
conda env create -f environment.yml && conda activate lant && pip install -e .
bash tools/install_external.sh            # DSurfTomo, pDSurfTomo, CC-FJpy at fixed versions
lant check
python examples/synthetic/make_synthetic.py && (cd examples/synthetic && lant run all -c synthetic.yaml)
lant init surface -o my_project.yaml      # only output_dir, raw_dir and station_file are needed
lant run all -c my_project.yaml
```

The detailed documentation ([tutorial](docs/tutorial.md), [parameters](docs/parameters.md), [AI workflow](docs/ai_workflow.md)) is in Chinese; code, configuration keys and figure labels are in English.

### License

MIT. Third-party programs are downloaded from their own repositories by `tools/install_external.sh` and keep their licenses (DSurfTomo MIT, pDSurfTomo GPL-3, CC-FJpy LGPL-3); none of them is included here.
