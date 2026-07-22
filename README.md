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

本项目实现了从原始 miniSEED 地震数据到横波速度模型的 **完整背景噪声面波层析成像流程**，支持面状台阵的三维层析反演和测线型台阵的伪二维剖面反演两条路径。

整个流程集成了 3 个第三方软件包（NoisePy、EGFAnalysisPy、DSurfTomo），并提供了一组辅助脚本将它们串联起来，形成端到端的自动化工作流。其中 **EGFAnalysisPy** 是 EGFAnalysisTimeFreq 的纯 Python 重写，并内置了基于深度学习的频散曲线提取工具 DisperPicker（含预训练模型），在一步内完成频散分析与曲线拾取，无需 MATLAB。

### 流程概览

```
原始 miniSEED 数据
  │
  ▼  Step 1 — 数据清洗（clean_mseed_for_noisepy.py）
  ▼  Step 2 — 准备台站文件（station.csv）
  ▼  Step 3 — 互相关计算与叠加（local_cc_workflow.py → NoisePy）
  ▼  Step 4 — 导出 EGF 输入数据（export_egfanalysis_input.py）
  ▼  Step 5 — 频散分析 + 曲线提取（EGFAnalysisPy/run.py 或 run_parallel.py）
  │
  ├──▼  Step 6A — 面波层析反演（DSurfTomo）        → 三维横波速度模型
  │       面状台阵，射线覆盖能形成二维网格
  │
  └──▼  Step 6B — 1D 反演 + 剖面拼接（lant_invert.py / lant_stitch.py）
          测线型台阵、浅层结构              → 伪二维横波速度剖面
```

Step 6 的两条路径按台阵形态二选一：测线型台阵的射线几乎共线，无法约束二维网格，DSurfTomo 不适用，此时走 Step 6B。

### 快速开始

```bash
# 1. 克隆仓库
git clone https://github.com/your-username/local_ANT.git
cd local_ANT

# 2. 安装环境（详见 docs/INSTALL.md）
conda create -n noisepy -y python=3.10 pip
conda activate noisepy
pip install noisepy-seis obspy numpy pandas scipy matplotlib pyyaml datetimerange h5py

# 3. 编译 DSurfTomo
cd DSurfTomo && python configure && cd ..

# 4. 按照 docs/workflow.md 逐步执行
```

### 项目结构

```
local_ANT/
├── scripts/                 # 辅助脚本（连接各软件包）
│   ├── clean_mseed_for_noisepy.py
│   ├── plot_station_distribution_all.py
│   ├── local_cc_workflow.py
│   ├── stack_visualization.py
│   ├── plot_stack_panels.py
│   ├── export_egfanalysis_input.py
│   ├── lant_invert.py       # Step 6B：逐台站对 1D 反演
│   └── lant_stitch.py       # Step 6B：拼接为伪二维剖面
├── configs/                 # 配置文件示例
├── docs/                    # 文档
│   ├── workflow.md          # 完整流程（中文）
│   ├── workflow_en.md       # 完整流程（English）
│   ├── INSTALL.md           # 安装指南（中文）
│   └── INSTALL_en.md        # 安装指南（English）
├── EGFAnalysisPy/           # 频散分析 + CNN 曲线提取（内置 DisperPicker）
│   ├── run.py               # 一步完成频散分析与拾取
│   ├── run_parallel.py      # 同上，多进程并行版
│   └── DisperPicker/        # 内置 DisperPicker 及预训练模型（saver/）
├── DSurfTomo/               # Step 6A：面波层析反演
│   └── scripts/             # 反演数据准备与绘图
└── LICENSE
```

### 文档

| 文档 | 说明 |
|------|------|
| [安装指南](docs/INSTALL.md) | 从零安装全部环境 |
| [完整流程](docs/workflow.md) | 6 个步骤的详细操作说明 |

### 引用

本项目集成了以下开源软件包：

> **NoisePy** — Jiang, C., & Denolle, M. A. (2020). NoisePy: A New High-Performance Python Tool for Ambient-Noise Seismology. *Seismological Research Letters*, 91(3), 1853–1866. [DOI](https://doi.org/10.1785/0220190364)

> **EGFAnalysisPy** (基于 EGFAnalysisTimeFreq 的 Python 重写) — Yao, H., van der Hilst, R. D., & de Hoop, M. V. (2006). Surface-wave array tomography in SE Tibet from ambient seismic noise and two-station analysis—I. Phase velocity maps. *Geophysical Journal International*, 166(2), 732–744. [DOI](https://doi.org/10.1111/j.1365-246X.2006.03028.x)

> **DisperPicker** — Yang, S., Zhang, H., Gu, N., et al. (2022). Automatically Extracting Surface Wave Group and Phase Velocity Dispersion Curves from Dispersion Spectrograms Using a Convolutional Neural Network. *Seismological Research Letters*. [DOI](https://doi.org/10.1785/0220210280)

> **DSurfTomo** — Fang, H., Yao, H., Zhang, H., et al. (2015). Direct inversion of surface wave dispersion for three-dimensional shallow crustal structure based on ray tracing: methodology and application. *Geophysical Journal International*, 201(3), 1251–1263. [DOI](https://doi.org/10.1093/gji/ggv069)

---

<!-- ==================== English ==================== -->

<a name="english"></a>

## English

<p align="right"><a href="#中文">切换到中文 ↗</a></p>

### Overview

This project implements a **complete ambient noise surface wave tomography pipeline** — from raw miniSEED seismic data to a shear-wave velocity model, supporting both 3D tomographic inversion for areal arrays and pseudo-2D section inversion for linear arrays.

The workflow integrates 3 third-party packages (NoisePy, EGFAnalysisPy, DSurfTomo) and provides a set of helper scripts that glue them together into an end-to-end automated workflow. **EGFAnalysisPy** is a pure-Python rewrite of EGFAnalysisTimeFreq with the deep-learning picker DisperPicker (and its pretrained model) bundled in, performing dispersion analysis and curve extraction in a single step — no MATLAB required.

### Pipeline

```
Raw miniSEED data
  │
  ▼  Step 1 — Data cleaning (clean_mseed_for_noisepy.py)
  ▼  Step 2 — Prepare station file (station.csv)
  ▼  Step 3 — Cross-correlation & stacking (local_cc_workflow.py → NoisePy)
  ▼  Step 4 — Export EGF input (export_egfanalysis_input.py)
  ▼  Step 5 — Dispersion analysis + curve extraction (EGFAnalysisPy/run.py or run_parallel.py)
  │
  ├──▼  Step 6A — Surface wave tomography (DSurfTomo)      → 3D shear-wave velocity model
  │       Areal array whose ray coverage constrains a 2D grid
  │
  └──▼  Step 6B — 1D inversion + stitching (lant_invert.py / lant_stitch.py)
          Linear array, shallow structure          → pseudo-2D velocity section
```

Choose between the two Step 6 paths by array geometry: in a linear array the ray paths are nearly collinear and cannot constrain a 2D grid, so DSurfTomo does not apply — use Step 6B instead.

### Quick Start

```bash
# 1. Clone the repository
git clone https://github.com/your-username/local_ANT.git
cd local_ANT

# 2. Set up the environment (see docs/INSTALL_en.md for details)
conda create -n noisepy -y python=3.10 pip
conda activate noisepy
pip install noisepy-seis obspy numpy pandas scipy matplotlib pyyaml datetimerange h5py

# 3. Compile DSurfTomo
cd DSurfTomo && python configure && cd ..

# 4. Follow docs/workflow_en.md step by step
```

### Project Structure

```
local_ANT/
├── scripts/                 # Helper scripts (glue between packages)
│   ├── clean_mseed_for_noisepy.py
│   ├── plot_station_distribution_all.py
│   ├── local_cc_workflow.py
│   ├── stack_visualization.py
│   ├── plot_stack_panels.py
│   ├── export_egfanalysis_input.py
│   ├── lant_invert.py       # Step 6B: per-pair 1D inversion
│   └── lant_stitch.py       # Step 6B: stitch into a pseudo-2D section
├── configs/                 # Example configuration files
├── docs/                    # Documentation
│   ├── workflow.md          # Full workflow (Chinese)
│   ├── workflow_en.md       # Full workflow (English)
│   ├── INSTALL.md           # Installation guide (Chinese)
│   └── INSTALL_en.md        # Installation guide (English)
├── EGFAnalysisPy/           # Dispersion analysis + CNN picker (DisperPicker bundled)
│   ├── run.py               # One-step dispersion analysis & picking
│   ├── run_parallel.py      # Same work, multiprocessing version
│   └── DisperPicker/        # Bundled DisperPicker with pretrained model (saver/)
├── DSurfTomo/               # Step 6A: Surface wave tomography inversion
│   └── scripts/             # Inversion data preparation & plotting
└── LICENSE
```

### Documentation

| Document | Description |
|----------|-------------|
| [Installation Guide](docs/INSTALL_en.md) | Install all dependencies from scratch |
| [Full Workflow](docs/workflow_en.md) | Detailed step-by-step instructions |

### Citations

This project integrates the following open-source packages:

> **NoisePy** — Jiang, C., & Denolle, M. A. (2020). NoisePy: A New High-Performance Python Tool for Ambient-Noise Seismology. *Seismological Research Letters*, 91(3), 1853–1866. [DOI](https://doi.org/10.1785/0220190364)

> **EGFAnalysisPy** (基于 EGFAnalysisTimeFreq 的 Python 重写) — Yao, H., van der Hilst, R. D., & de Hoop, M. V. (2006). Surface-wave array tomography in SE Tibet from ambient seismic noise and two-station analysis—I. Phase velocity maps. *Geophysical Journal International*, 166(2), 732–744. [DOI](https://doi.org/10.1111/j.1365-246X.2006.03028.x)

> **DisperPicker** — Yang, S., Zhang, H., Gu, N., et al. (2022). Automatically Extracting Surface Wave Group and Phase Velocity Dispersion Curves from Dispersion Spectrograms Using a Convolutional Neural Network. *Seismological Research Letters*. [DOI](https://doi.org/10.1785/0220210280)

> **DSurfTomo** — Fang, H., Yao, H., Zhang, H., et al. (2015). Direct inversion of surface wave dispersion for three-dimensional shallow crustal structure based on ray tracing: methodology and application. *Geophysical Journal International*, 201(3), 1251–1263. [DOI](https://doi.org/10.1093/gji/ggv069)

---

## License

This project is licensed under the [MIT License](LICENSE). Third-party components (DSurfTomo, EGFAnalysisPy, and the DisperPicker bundled within it) retain their original licenses.
