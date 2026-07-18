# 安装指南

本指南将帮助你从零开始安装本项目所需的全部软件环境。请按照顺序完成以下步骤。

## 目录

- [系统要求](#系统要求)
- [1. 安装 Conda（如果尚未安装）](#1-安装-conda如果尚未安装)
- [2. 安装 NoisePy 及 Python 依赖](#2-安装-noisepy-及-python-依赖)
- [3. 安装 EGFAnalysisPy](#3-安装-egfanalysispy)
- [4. 编译 DSurfTomo](#4-编译-dsurftomo)
- [5. 安装 GMT（可选，用于绘图）](#5-安装-gmt可选用于绘图)
- [验证安装](#验证安装)
- [常见问题](#常见问题)

---

## 系统要求

| 项目 | 要求 |
|------|------|
| 操作系统 | Linux（推荐 Ubuntu 18.04+）或 macOS |
| Python | 3.8 – 3.10（推荐 3.10） |
| Fortran 编译器 | gfortran（DSurfTomo 编译需要） |
| 磁盘空间 | ≥ 5 GB（含 Conda 环境和软件包） |
| 内存 | ≥ 4 GB |

---

## 1. 安装 Conda（如果尚未安装）

Conda 是一个 Python 环境管理工具，可以帮你创建相互隔离的 Python 环境，避免不同项目之间的依赖冲突。

> 官方安装文档：https://docs.conda.io/en/latest/miniconda.html

如果你已经安装了 Conda（Anaconda 或 Miniconda），可以跳过这一步。

### 检查是否已安装

```bash
conda --version
```

如果显示版本号（如 `conda 24.x.x`），说明已安装，跳到第 2 步。

### 安装 Miniconda（轻量版，推荐）

```bash
# 下载 Miniconda 安装脚本（Linux）
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# 运行安装（按提示操作，一路回车/yes 即可）
bash Miniconda3-latest-Linux-x86_64.sh

# 安装完成后，重新打开终端，或执行：
source ~/.bashrc
```

macOS 用户请将下载链接替换为：
```
https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
```

Apple Silicon (M1/M2/M3) 用户：
```
https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
```

---

## 2. 安装 NoisePy 及 Python 依赖

步骤 01–04 及 06b 的入口脚本都运行在 noisepy 环境中。

> NoisePy 官方安装文档：https://github.com/noisepy/NoisePy

### 2a. 创建环境并安装 NoisePy

```bash
# 创建名为 noisepy 的环境，指定 Python 3.10
conda create -n noisepy -y python=3.10 pip

# 激活环境（每次打开新终端都需要执行这一步）
conda activate noisepy

# 安装 NoisePy
pip install noisepy-seis
```

> **什么是 NoisePy？**
> NoisePy 是一个计算背景噪声互相关函数的 Python 软件包。
> 引用：Jiang, C., & Denolle, M. A. (2020). *Seismological Research Letters*, 91(3), 1853–1866.
> GitHub: https://github.com/noisepy/NoisePy

### 2b. 辅助脚本的依赖（自动带全，无需单独安装）

`pip install noisepy-seis` 会通过依赖关系自动安装本项目辅助脚本需要的全部包：
obspy（经 pyasdf）、numpy、pandas、scipy、matplotlib（经 obspy）、pyyaml、
datetimerange、h5py、tqdm——已在实际环境中逐一核验，**不需要再手动安装**。

可用以下命令确认：

```bash
conda activate noisepy
python -c "import obspy, numpy, pandas, scipy, matplotlib, yaml, datetimerange, h5py, tqdm; print('all deps OK')"
```


### 2c. 关于 MPI

本项目的入口脚本**目前不支持 MPI 并行**（已列入 TODO.md）。NoisePy 本体
支持 MPI（`pip install noisepy-seis[mpi]`），如需自行使用请参阅
NoisePy 官方文档：https://github.com/noisepy/NoisePy

### 2d. （可选）安装 Jupyter Notebook

如果你想在 Jupyter Notebook 中使用 NoisePy：

```bash
pip install ipykernel notebook
python -m ipykernel install --user --name noisepy
```

### 2e. 线状台阵反演环境（evodcinv，独立 Conda 环境）

线状台阵反演（步骤 06a，`scripts/06a_line_inversion/`）使用
[evodcinv](https://github.com/keurfonluu/evodcinv) 做一维频散反演，
其正演内核为 [disba](https://github.com/keurfonluu/disba)。
**请为它单独创建一个 Conda 环境，不要装进 noisepy 环境：**

```bash
conda create -n evodcinv -y python=3.10 pip
conda activate evodcinv
pip install evodcinv          # 官方推荐安装方式，会自动安装 disba
pip install pyyaml matplotlib scipy tqdm   # 06a 入口脚本的其余依赖
```

也可以直接用仓库提供的环境定义一步创建：

```bash
conda env create -f environments/evodcinv.yml
```

验证：

```bash
conda activate evodcinv
python -c "import evodcinv, disba; print('evodcinv OK, disba', disba.__version__)"
```

> **版本说明**：原项目机器上确认过 disba 0.6.1 可用；evodcinv 的具体版本
> 未能从原环境考证（待确认），如遇接口不兼容，可尝试
> `pip install "evodcinv==2.*"` 并参阅其官方文档。

---

## 3. 安装 EGFAnalysisPy

EGFAnalysisPy 是 EGFAnalysisTimeFreq 的纯 Python 重写，用于面波频散分析，并**内置了基于深度学习的频散曲线提取工具 DisperPicker**，一步完成频散分析与曲线拾取。它需要单独的 Python 环境，因为它依赖 TensorFlow 1.x，与 NoisePy 使用的包版本不兼容。

> **什么是 EGFAnalysisPy / DisperPicker？**
> EGFAnalysisPy 从经验格林函数 / 互相关函数计算频散能量图（Yao et al., 2006），并用内置的 DisperPicker（CNN，Yang et al., 2022, *SRL*）自动提取频散曲线。
> DisperPicker GitHub: https://github.com/ShaoboYang-USTC/DisperPicker

### 3a. 创建独立环境

```bash
# EGFAnalysisPy 依赖 TensorFlow 1.13，需要 Python 3.6
conda create -n egfanalysispy -y python=3.6 pip
conda activate egfanalysispy
```

### 3b. 安装依赖

```bash
cd EGFAnalysisPy
pip install -r requirements.txt
```

> **注意：** TensorFlow 1.13.1 不支持 Python 3.7+。如果你遇到安装失败，请确认 Python 版本为 3.6。
> TensorFlow 1.13.1 官方页面：https://pypi.org/project/tensorflow/1.13.1/

### 3c. 验证安装

```bash
conda activate egfanalysispy
cd EGFAnalysisPy
python -c "import tensorflow as tf; print('TensorFlow version:', tf.__version__)"
```

应该输出 `TensorFlow version: 1.13.1`。

### 3d. 预训练模型

DisperPicker 的预训练模型已经包含在 `EGFAnalysisPy/DisperPicker/saver/` 目录中，无需额外下载：

```
EGFAnalysisPy/DisperPicker/saver/
├── -10000.data-00000-of-00001
├── -10000.index
└── -10000.meta
```

---

## 4. 编译 DSurfTomo

DSurfTomo 是一个 Fortran 编写的面波层析成像程序，需要编译才能使用。

> **什么是 DSurfTomo？**
> DSurfTomo 直接将面波频散数据反演为三维横波速度模型，无需中间步骤构建相/群速度图。
> 引用：Fang, H., et al. (2015). *Geophysical Journal International*, 201(3), 1251–1263.
> 官方仓库：https://github.com/HongjianFang/DSurfTomo（本仓库内置修改版，直接编译内置版即可）

### 4a. 安装 Fortran 编译器

```bash
# Ubuntu/Debian
sudo apt update
sudo apt install gfortran

# CentOS/RHEL
sudo yum install gcc-gfortran

# macOS (需要先安装 Homebrew: https://brew.sh)
brew install gcc
```

验证安装：

```bash
gfortran --version
```

应该显示版本信息，如 `GNU Fortran (Ubuntu 9.4.0-...) 9.4.0`。

### 4b. 编译

```bash
cd DSurfTomo

# 运行配置脚本（这是一个 Python 脚本，会调用 make）
python configure
```

该脚本会自动：
1. 创建 `bin/` 目录
2. 进入 `src/` 目录执行 `make`
3. 将编译好的可执行文件 `DSurfTomo` 复制到 `bin/`

### 4c. 验证编译

```bash
ls -la DSurfTomo/bin/DSurfTomo
```

如果文件存在，说明编译成功。

### 4d. （可选）添加到 PATH

为了方便使用，可以将 DSurfTomo 添加到系统路径：

```bash
# 将以下行添加到 ~/.bashrc 或 ~/.zshrc
export PATH="/path/to/local_ANT/DSurfTomo/bin:$PATH"

# 使其生效
source ~/.bashrc
```

### 4e. 编译失败排查

| 错误信息 | 解决方法 |
|----------|----------|
| `gfortran: command not found` | 安装 gfortran（见 4a） |
| `fatal error: cannot find -lgomp` | 安装 OpenMP: `sudo apt install libgomp1` |
| `m64: No such file or directory` | 确认系统为 64 位 |

---

## 5. 安装 GMT（可选，用于绘图）

GMT (Generic Mapping Tools) 用于 DSurfTomo 的结果绘图（深度切片、剖面、射线路径、棋盘测试）。如果你不需要绘图，可以跳过此步。

> 官方网站与安装文档：https://www.generic-mapping-tools.org/ 、https://docs.generic-mapping-tools.org/latest/install.html

### 5a. 安装 GMT 6

```bash
# Ubuntu/Debian
sudo apt update
sudo apt install gmt gmt-dcw gmt-gshhg

# CentOS/RHEL
sudo yum install GMT

# macOS
brew install gmt
```

### 5b. 验证安装

```bash
gmt --version
```

应该显示 `6.x.x`。本项目的 GMT 脚本需要 GMT 6.0 以上版本。

---

## 验证安装

运行以下命令验证所有组件是否安装正确：

```bash
echo "=== 1. NoisePy 环境 ==="
conda activate noisepy
python -c "
import noisepy; print('NoisePy OK')
import obspy; print('ObsPy OK')
import h5py; print('h5py OK')
import numpy; print('NumPy OK')
import pandas; print('Pandas OK')
import yaml; print('PyYAML OK')
import scipy; print('SciPy OK')
import matplotlib; print('Matplotlib OK')
import datetimerange; print('datetimerange OK')
"

echo ""
echo "=== 1b. evodcinv 环境（可选，步骤 06a 需要）==="
conda activate evodcinv
python -c "
import evodcinv; print('evodcinv OK')
import disba; print('disba OK')
import tqdm; print('tqdm OK')
" || echo "evodcinv 环境未创建（仅做面状反演可忽略；见 2e 节）"

echo ""
echo "=== 2. EGFAnalysisPy 环境 ==="
conda activate egfanalysispy
python -c "
import tensorflow as tf; print('TensorFlow', tf.__version__, 'OK')
import numpy; print('NumPy OK')
import scipy; print('SciPy OK')
import pandas; print('Pandas OK')
import matplotlib; print('Matplotlib OK')
"

echo ""
echo "=== 3. DSurfTomo ==="
if [ -f DSurfTomo/bin/DSurfTomo ]; then
    echo "DSurfTomo binary OK"
else
    echo "DSurfTomo binary NOT FOUND"
fi

echo ""
echo "=== 4. gfortran ==="
gfortran --version 2>/dev/null | head -1 || echo "gfortran NOT FOUND"

echo ""
echo "=== 5. GMT (可选) ==="
gmt --version 2>/dev/null || echo "GMT not installed (optional)"
```

---

## 常见问题

### Q: `conda activate noisepy` 提示 "conda: command not found"
请确认 Conda 已正确安装。运行 `source ~/.bashrc`（或 `source ~/.zshrc`）后重试。如果仍然找不到，请检查安装路径是否在 `$PATH` 中。

### Q: `pip install noisepy-seis` 超时或下载慢
使用国内镜像源加速：
```bash
pip install noisepy-seis -i https://pypi.tuna.tsinghua.edu.cn/simple
```

### Q: EGFAnalysisPy 安装 TensorFlow 1.13.1 失败
TensorFlow 1.13.1 仅支持 Python 3.6。请确认：
```bash
python --version  # 应该显示 3.6.x
```
如果版本不对，重新创建环境：`conda create -n egfanalysispy -y python=3.6 pip`

### Q: TensorFlow 报 GPU 相关错误
EGFAnalysisPy（内置 DisperPicker）默认使用 CPU 运行，不需要 GPU。如果看到 CUDA/GPU 相关的警告信息，可以忽略，或设置环境变量：
```bash
export CUDA_VISIBLE_DEVICES=""
```

### Q: DSurfTomo 编译报错 `gfortran: command not found`
安装 gfortran：
```bash
sudo apt install gfortran    # Ubuntu/Debian
sudo yum install gcc-gfortran # CentOS/RHEL
brew install gcc               # macOS
```

### Q: 两个 Conda 环境可以合并吗？
不建议。NoisePy 需要 Python 3.10 和较新的 NumPy，而 EGFAnalysisPy 需要 Python 3.6 和 TensorFlow 1.13。两者无法在同一环境中共存。使用时切换环境即可：
```bash
conda activate noisepy        # 运行 NoisePy 相关脚本
conda activate egfanalysispy          # 运行 EGFAnalysisPy（频散分析 + 拾取）
```

### Q: 运行 DSurfTomo 脚本时报 `ModuleNotFoundError: No module named 'numpy'`
DSurfTomo 的 Python 脚本（`DSurfTomo/scripts/` 下的）需要在 `noisepy` 环境中运行：
```bash
conda activate noisepy
cd DSurfTomo/scripts
python PrepareDSurfTomoProject.py
```

---

## 环境切换速查

| 步骤 | 用户入口 | 使用的环境 |
|------|----------|------------|
| 01 数据清洗 | `scripts/01_clean.py` | `conda activate noisepy` |
| 02 互相关 + 叠加 | `scripts/02_correlation.py` | `conda activate noisepy` |
| 03 叠加质检可视化 | `scripts/03_stack_qc.py`、`scripts/03b_stack_panels.py` | `conda activate noisepy` |
| 04 导出 EGF 输入 | `scripts/04_export_egf.py` | `conda activate noisepy` |
| 05 频散分析 + 拾取 | `scripts/05_dispersion.py` | `conda activate egfanalysispy` |
| 06a 线状反演 | `scripts/06a_line_inversion/01_invert.py`、`02_stitch.py` | `conda activate evodcinv`（独立环境，见 2e） |
| 06b 面状反演 | `scripts/06b_area_inversion/01_prepare.py`、`02_invert.py`、`03_plot.py` | `conda activate noisepy` + 编译好的 `DSurfTomo/bin/DSurfTomo`；绘图需 GMT 6 |
