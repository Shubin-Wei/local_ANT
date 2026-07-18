<h1 align="center">local_ANT</h1>

<p align="center">
  Ambient Noise Surface Wave Tomography — Full Pipeline from Local miniSEED Data
</p>

<p align="center">
  <a href="#中文">🇨🇳 中文</a>&nbsp;&nbsp;|&nbsp;&nbsp;<a href="#english">🇬🇧 English</a>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/python-3.10%20%7C%203.6-blue" alt="Python">
  <img src="https://img.shields.io/badge/license-MIT-green" alt="License">
</p>

---

<!-- ==================== 中文 ==================== -->

<a name="中文"></a>

## 中文

<p align="right"><a href="#english">Switch to English ↗</a></p>

### 简介

本项目实现了从原始 miniSEED 地震数据到横波速度模型的 **完整背景噪声成像流程**。

公共流程(01–05)集成 NoisePy(互相关)与 EGFAnalysisPy(频散分析,内置深度学习拾取工具 DisperPicker 及预训练模型,无需 MATLAB);第 05 步之后按台阵类型二选一:

- **06a 线状台阵反演**:evodcinv 逐台站对一维反演 → 拼接伪二维剖面;
- **06b 面状台阵反演**:DSurfTomo 直接反演三维横波速度结构。

### 流程概览

```
原始 miniSEED 数据
  │
  ▼  01 数据清洗                scripts/01_clean.py
  ▼  02 互相关 + 叠加           scripts/02_correlation.py（NoisePy）
  ▼  03 叠加质检可视化          scripts/03_stack_qc.py / 03b_stack_panels.py
  ▼  04 导出 EGF 输入           scripts/04_export_egf.py
  ▼  05 频散分析 + 曲线拾取     scripts/05_dispersion.py（EGFAnalysisPy + DisperPicker）
  │
  ├── 06a 线状台阵反演          scripts/06a_line_inversion/（evodcinv → 伪二维剖面）
  └── 06b 面状台阵反演          scripts/06b_area_inversion/（DSurfTomo → 三维模型）
```

每个步骤独立运行、可从中间步骤继续;所有入口都从项目根目录执行,配置集中在 `configs/`。

### 快速开始

```bash
# 1. 克隆仓库
git clone https://github.com/your-username/local_ANT.git
cd local_ANT

# 2. 创建环境(详见 docs/INSTALL.md)
conda env create -f environments/noisepy.yml        # 步骤 01–04、06a、06b
conda env create -f environments/egfanalysispy.yml  # 步骤 05

# 3. 编译 DSurfTomo(仅面状反演 06b 需要)
cd DSurfTomo && python configure && cd ..

# 4. 检查环境
conda activate noisepy && python tools/check_environment.py

# 5. 复制示例配置并按需修改,然后按 docs/workflow.md 逐步执行
python scripts/01_clean.py --config configs/01_clean.yaml
```

### 各步骤运行环境

| 步骤 | 用户入口 | 环境 | 主要输入 | 主要输出 |
|------|----------|------|----------|----------|
| 01 数据清洗 | `scripts/01_clean.py` | noisepy | 原始 miniSEED | 清洗后 `.ms` |
| 02 互相关+叠加 | `scripts/02_correlation.py` | noisepy | 清洗数据 + station.csv | CCF + Stack(HDF5) |
| 03 叠加质检 | `scripts/03_stack_qc.py`、`03b_stack_panels.py` | noisepy | Stack | moveout 图 |
| 04 导出 EGF 输入 | `scripts/04_export_egf.py` | noisepy | Stack | CF `.dat` 波形 |
| 05 频散分析拾取 | `scripts/05_dispersion.py` | **egfanalysispy** | CF `.dat` | 频散曲线 `qc_result/` |
| 06a 线状反演 | `scripts/06a_line_inversion/01_invert.py`、`02_stitch.py` | noisepy(含 evodcinv) | 频散曲线 | 一维模型 + 伪二维剖面 |
| 06b 面状反演 | `scripts/06b_area_inversion/01_prepare.py`、`02_invert.py`、`03_plot.py` | noisepy + DSurfTomo 可执行文件 + GMT 6 | 频散曲线 | 三维模型 + 深度切片图 |

### 项目结构

```
local_ANT/
├── scripts/                     # 编号用户入口(从根目录运行)
│   ├── 01_clean.py … 05_dispersion.py
│   ├── 06a_line_inversion/      # 线状反演:01_invert.py, 02_stitch.py
│   ├── 06b_area_inversion/      # 面状反演:01_prepare.py, 02_invert.py, 03_plot.py
│   └── plot_station_distribution.py
├── configs/                     # 全部用户配置(00_project + 各阶段 example)
├── environments/                # Conda 环境定义(noisepy.yml, egfanalysispy.yml)
├── tools/                       # check_environment.py 环境自检
├── docs/                        # 安装、流程、配置、反演、FAQ 文档(中英)
├── EGFAnalysisPy/               # 频散分析 + DisperPicker(含预训练模型 saver/)
├── DSurfTomo/                   # 面波层析反演(Fortran 源码 + 准备/绘图脚本)
├── TODO.md
└── LICENSE
```

### 文档

| 文档 | 说明 |
|------|------|
| [安装指南](docs/INSTALL.md) | 从零安装全部环境(含 DSurfTomo 编译) |
| [完整流程](docs/workflow.md) | 01–06 各步骤详细操作(含旧脚本名迁移对照) |
| [配置说明](docs/CONFIG.md) | 配置分层、优先级与运行时配置机制 |
| [线状反演](docs/line_inversion.md) | 06a 分支:适用场景、配置与运行 |
| [面状反演](docs/area_inversion.md) | 06b 分支:DSurfTomo 编译、配置与运行 |
| [常见问题](docs/FAQ.md) | FAQ |

### 引用

本项目集成了以下开源软件包:

> **NoisePy** — Jiang, C., & Denolle, M. A. (2020). NoisePy: A New High-Performance Python Tool for Ambient-Noise Seismology. *Seismological Research Letters*, 91(3), 1853–1866. [DOI](https://doi.org/10.1785/0220190364)

> **EGFAnalysisPy** (基于 EGFAnalysisTimeFreq 的 Python 重写) — Yao, H., van der Hilst, R. D., & de Hoop, M. V. (2006). Surface-wave array tomography in SE Tibet from ambient seismic noise and two-station analysis—I. Phase velocity maps. *Geophysical Journal International*, 166(2), 732–744. [DOI](https://doi.org/10.1111/j.1365-246X.2006.03028.x)

> **DisperPicker** — Yang, S., Zhang, H., Gu, N., et al. (2022). Automatically Extracting Surface Wave Group and Phase Velocity Dispersion Curves from Dispersion Spectrograms Using a Convolutional Neural Network. *Seismological Research Letters*. [DOI](https://doi.org/10.1785/0220210280)

> **DSurfTomo** — Fang, H., Yao, H., Zhang, H., et al. (2015). Direct inversion of surface wave dispersion for three-dimensional shallow crustal structure based on ray tracing: methodology and application. *Geophysical Journal International*, 201(3), 1251–1263. [DOI](https://doi.org/10.1093/gji/ggv069)

> **evodcinv / disba** — Luu, K. evodcinv: Inversion of dispersion curves using evolutionary algorithms. [GitHub](https://github.com/keurfonluu/evodcinv)

---

<!-- ==================== English ==================== -->

<a name="english"></a>

## English

<p align="right"><a href="#中文">切换到中文 ↗</a></p>

### Overview

This project implements a **complete ambient noise imaging pipeline** — from raw miniSEED seismic data to shear-wave velocity models.

The common stages (01–05) integrate NoisePy (cross-correlation) and EGFAnalysisPy (dispersion analysis, with the deep-learning picker DisperPicker and its pretrained model bundled — no MATLAB required). After step 05 the pipeline splits into two parallel branches, chosen by array geometry:

- **06a Linear-array inversion**: per-pair 1-D inversion with evodcinv → stitched pseudo-2D section;
- **06b Areal-array inversion**: DSurfTomo direct inversion for a 3-D shear-wave model.

### Pipeline

```
Raw miniSEED data
  │
  ▼  01 Data cleaning             scripts/01_clean.py
  ▼  02 Cross-correlation + stack scripts/02_correlation.py (NoisePy)
  ▼  03 Stack QC visualization    scripts/03_stack_qc.py / 03b_stack_panels.py
  ▼  04 Export EGF input          scripts/04_export_egf.py
  ▼  05 Dispersion + picking      scripts/05_dispersion.py (EGFAnalysisPy + DisperPicker)
  │
  ├── 06a Linear-array inversion  scripts/06a_line_inversion/ (evodcinv → pseudo-2D section)
  └── 06b Areal-array inversion   scripts/06b_area_inversion/ (DSurfTomo → 3-D model)
```

Every step runs independently and can be resumed midway; all entries are executed from the repo root, with configuration centralized in `configs/`.

### Quick Start

```bash
# 1. Clone the repository
git clone https://github.com/your-username/local_ANT.git
cd local_ANT

# 2. Create the environments (see docs/INSTALL_en.md)
conda env create -f environments/noisepy.yml        # steps 01–04, 06a, 06b
conda env create -f environments/egfanalysispy.yml  # step 05

# 3. Compile DSurfTomo (only needed for areal inversion 06b)
cd DSurfTomo && python configure && cd ..

# 4. Check the environment
conda activate noisepy && python tools/check_environment.py

# 5. Copy the example configs, edit them, then follow docs/workflow_en.md
python scripts/01_clean.py --config configs/01_clean.yaml
```

### Environment per step

| Step | User entry | Environment | Main input | Main output |
|------|------------|-------------|------------|-------------|
| 01 Data cleaning | `scripts/01_clean.py` | noisepy | raw miniSEED | cleaned `.ms` |
| 02 CC + stacking | `scripts/02_correlation.py` | noisepy | cleaned data + station.csv | CCF + Stack (HDF5) |
| 03 Stack QC | `scripts/03_stack_qc.py`, `03b_stack_panels.py` | noisepy | Stack | moveout figures |
| 04 Export EGF input | `scripts/04_export_egf.py` | noisepy | Stack | CF `.dat` waveforms |
| 05 Dispersion + picking | `scripts/05_dispersion.py` | **egfanalysispy** | CF `.dat` | dispersion curves `qc_result/` |
| 06a Linear inversion | `scripts/06a_line_inversion/01_invert.py`, `02_stitch.py` | noisepy (with evodcinv) | dispersion curves | 1-D models + pseudo-2D section |
| 06b Areal inversion | `scripts/06b_area_inversion/01_prepare.py`, `02_invert.py`, `03_plot.py` | noisepy + DSurfTomo binary + GMT 6 | dispersion curves | 3-D model + depth slices |

### Project Structure

```
local_ANT/
├── scripts/                     # Numbered user entries (run from repo root)
│   ├── 01_clean.py … 05_dispersion.py
│   ├── 06a_line_inversion/      # linear branch: 01_invert.py, 02_stitch.py
│   ├── 06b_area_inversion/      # areal branch: 01_prepare.py, 02_invert.py, 03_plot.py
│   └── plot_station_distribution.py
├── configs/                     # All user configs (00_project + per-stage examples)
├── environments/                # Conda env definitions (noisepy.yml, egfanalysispy.yml)
├── tools/                       # check_environment.py
├── docs/                        # Install, workflow, config, inversion, FAQ (CN & EN)
├── EGFAnalysisPy/               # Dispersion analysis + DisperPicker (pretrained model in saver/)
├── DSurfTomo/                   # Surface-wave tomography (Fortran + prep/plot scripts)
├── TODO.md
└── LICENSE
```

### Documentation

| Document | Description |
|----------|-------------|
| [Installation Guide](docs/INSTALL_en.md) | Install everything from scratch (incl. DSurfTomo build) |
| [Full Workflow](docs/workflow_en.md) | Detailed steps 01–06 (with old→new script-name mapping) |
| [Configuration](docs/CONFIG.md) | Config layering, precedence, runtime-config mechanism |
| [Linear inversion](docs/line_inversion.md) | Branch 06a: scope, config, running |
| [Areal inversion](docs/area_inversion.md) | Branch 06b: DSurfTomo build, config, running |
| [FAQ](docs/FAQ.md) | Frequently asked questions |

### Citations

This project integrates the following open-source packages:

> **NoisePy** — Jiang, C., & Denolle, M. A. (2020). NoisePy: A New High-Performance Python Tool for Ambient-Noise Seismology. *Seismological Research Letters*, 91(3), 1853–1866. [DOI](https://doi.org/10.1785/0220190364)

> **EGFAnalysisPy** (Python rewrite of EGFAnalysisTimeFreq) — Yao, H., van der Hilst, R. D., & de Hoop, M. V. (2006). Surface-wave array tomography in SE Tibet from ambient seismic noise and two-station analysis—I. Phase velocity maps. *Geophysical Journal International*, 166(2), 732–744. [DOI](https://doi.org/10.1111/j.1365-246X.2006.03028.x)

> **DisperPicker** — Yang, S., Zhang, H., Gu, N., et al. (2022). Automatically Extracting Surface Wave Group and Phase Velocity Dispersion Curves from Dispersion Spectrograms Using a Convolutional Neural Network. *Seismological Research Letters*. [DOI](https://doi.org/10.1785/0220210280)

> **DSurfTomo** — Fang, H., Yao, H., Zhang, H., et al. (2015). Direct inversion of surface wave dispersion for three-dimensional shallow crustal structure based on ray tracing: methodology and application. *Geophysical Journal International*, 201(3), 1251–1263. [DOI](https://doi.org/10.1093/gji/ggv069)

> **evodcinv / disba** — Luu, K. evodcinv: Inversion of dispersion curves using evolutionary algorithms. [GitHub](https://github.com/keurfonluu/evodcinv)

---

## License

This project is licensed under the [MIT License](LICENSE). Third-party components (DSurfTomo, EGFAnalysisPy, and the DisperPicker bundled within it) retain their original licenses.
