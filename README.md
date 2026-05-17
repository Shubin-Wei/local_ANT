<h1 align="center">local_ANT</h1>

<p align="center">
  Ambient Noise Surface Wave Tomography — Full Pipeline from Local miniSEED Data
</p>

<p align="center">
  <a href="#中文">🇨🇳 中文</a>&nbsp;&nbsp;|&nbsp;&nbsp;<a href="#english">🇬🇧 English</a>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/python-3.10-blue" alt="Python">
  <img src="https://img.shields.io/badge/MATLAB-required-orange" alt="MATLAB">
  <img src="https://img.shields.io/badge/license-MIT-green" alt="License">
</p>

---

<!-- ==================== 中文 ==================== -->

<a name="中文"></a>

## 中文

<p align="right"><a href="#english">Switch to English ↗</a></p>

### 简介

本项目实现了从原始 miniSEED 地震数据到三维横波速度模型的 **完整背景噪声面波层析成像流程**。

整个流程集成了 4 个第三方软件包（NoisePy、EGFAnalysisTimeFreq、DisperPicker、DSurfTomo），并提供了一组辅助脚本将它们串联起来，形成端到端的自动化工作流。

### 流程概览

```
原始 miniSEED 数据
  │
  ▼  Step 1 — 数据清洗（clean_mseed_for_noisepy.py）
  ▼  Step 2 — 准备台站文件（station.csv）
  ▼  Step 3 — 互相关计算与叠加（local_cc_workflow.py → NoisePy）
  ▼  Step 4 — 导出 EGF 输入数据（export_egfanalysis_input.py）
  ▼  Step 5 — 频散分析（EGFAnalysisTimeFreq, MATLAB）
  ▼  Step 6 — 格式转换（prepare_disperpicker_input.py）
  ▼  Step 7 — 频散曲线提取（DisperPicker）
  ▼  Step 8 — 面波层析反演（DSurfTomo）
  │
  ▼
三维横波速度模型
```

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
│   ├── local_cc_workflow.py
│   ├── stack_visualization.py
│   ├── export_egfanalysis_input.py
│   └── prepare_disperpicker_input.py
├── configs/                 # 配置文件示例
├── docs/                    # 文档
│   ├── workflow.md          # 完整流程（中文）
│   ├── workflow_en.md       # 完整流程（English）
│   ├── INSTALL.md           # 安装指南（中文）
│   └── INSTALL_en.md        # 安装指南（English）
├── EGFAnalysisTimeFreq/     # MATLAB 频散分析
├── DisperPicker/            # CNN 频散曲线提取
├── DSurfTomo/               # 面波层析反演
│   └── scripts/             # 反演数据准备与绘图
└── LICENSE
```

### 文档

| 文档 | 说明 |
|------|------|
| [安装指南](docs/INSTALL.md) | 从零安装全部环境 |
| [完整流程](docs/workflow.md) | 8 个步骤的详细操作说明 |

### 引用

本项目集成了以下开源软件包：

> **NoisePy** — Jiang, C., & Denolle, M. A. (2020). NoisePy: A New High-Performance Python Tool for Ambient-Noise Seismology. *Seismological Research Letters*, 91(3), 1853–1866. [DOI](https://doi.org/10.1785/0220190364)

> **EGFAnalysisTimeFreq** — Yao, H., van der Hilst, R. D., & de Hoop, M. V. (2006). Surface-wave array tomography in SE Tibet from ambient seismic noise and two-station analysis—I. Phase velocity maps. *Geophysical Journal International*, 166(2), 732–744. [DOI](https://doi.org/10.1111/j.1365-246X.2006.03028.x)

> **DisperPicker** — Yang, S., Zhang, H., Gu, N., et al. (2022). Automatically Extracting Surface Wave Group and Phase Velocity Dispersion Curves from Dispersion Spectrograms Using a Convolutional Neural Network. *Seismological Research Letters*. [DOI](https://doi.org/10.1785/0220210280)

> **DSurfTomo** — Fang, H., Yao, H., Zhang, H., et al. (2015). Direct inversion of surface wave dispersion for three-dimensional shallow crustal structure based on ray tracing: methodology and application. *Geophysical Journal International*, 201(3), 1251–1263. [DOI](https://doi.org/10.1093/gji/ggv069)

---

<!-- ==================== English ==================== -->

<a name="english"></a>

## English

<p align="right"><a href="#中文">切换到中文 ↗</a></p>

### Overview

This project implements a **complete ambient noise surface wave tomography pipeline** — from raw miniSEED seismic data to a 3D shear-wave velocity model.

The workflow integrates 4 third-party packages (NoisePy, EGFAnalysisTimeFreq, DisperPicker, DSurfTomo) and provides a set of helper scripts that glue them together into an end-to-end automated workflow.

### Pipeline

```
Raw miniSEED data
  │
  ▼  Step 1 — Data cleaning (clean_mseed_for_noisepy.py)
  ▼  Step 2 — Prepare station file (station.csv)
  ▼  Step 3 — Cross-correlation & stacking (local_cc_workflow.py → NoisePy)
  ▼  Step 4 — Export EGF input (export_egfanalysis_input.py)
  ▼  Step 5 — Dispersion analysis (EGFAnalysisTimeFreq, MATLAB)
  ▼  Step 6 — Format conversion (prepare_disperpicker_input.py)
  ▼  Step 7 — Dispersion curve extraction (DisperPicker)
  ▼  Step 8 — Surface wave tomography inversion (DSurfTomo)
  │
  ▼
3D shear-wave velocity model
```

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
│   ├── local_cc_workflow.py
│   ├── stack_visualization.py
│   ├── export_egfanalysis_input.py
│   └── prepare_disperpicker_input.py
├── configs/                 # Example configuration files
├── docs/                    # Documentation
│   ├── workflow.md          # Full workflow (Chinese)
│   ├── workflow_en.md       # Full workflow (English)
│   ├── INSTALL.md           # Installation guide (Chinese)
│   └── INSTALL_en.md        # Installation guide (English)
├── EGFAnalysisTimeFreq/     # MATLAB dispersion analysis
├── DisperPicker/            # CNN-based dispersion curve picker
├── DSurfTomo/               # Surface wave tomography inversion
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

> **EGFAnalysisTimeFreq** — Yao, H., van der Hilst, R. D., & de Hoop, M. V. (2006). Surface-wave array tomography in SE Tibet from ambient seismic noise and two-station analysis—I. Phase velocity maps. *Geophysical Journal International*, 166(2), 732–744. [DOI](https://doi.org/10.1111/j.1365-246X.2006.03028.x)

> **DisperPicker** — Yang, S., Zhang, H., Gu, N., et al. (2022). Automatically Extracting Surface Wave Group and Phase Velocity Dispersion Curves from Dispersion Spectrograms Using a Convolutional Neural Network. *Seismological Research Letters*. [DOI](https://doi.org/10.1785/0220210280)

> **DSurfTomo** — Fang, H., Yao, H., Zhang, H., et al. (2015). Direct inversion of surface wave dispersion for three-dimensional shallow crustal structure based on ray tracing: methodology and application. *Geophysical Journal International*, 201(3), 1251–1263. [DOI](https://doi.org/10.1093/gji/ggv069)

---

## License

This project is licensed under the [MIT License](LICENSE). Third-party components (DSurfTomo, DisperPicker, EGFAnalysisTimeFreq) retain their original licenses.
