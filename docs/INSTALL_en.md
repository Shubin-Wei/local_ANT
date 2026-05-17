# Installation Guide

This guide walks you through installing all software required by this project, from scratch. Follow the steps in order.

## Table of Contents

- [System Requirements](#system-requirements)
- [1. Install Conda (if not already installed)](#1-install-conda-if-not-already-installed)
- [2. Install NoisePy and Python Dependencies](#2-install-noisepy-and-python-dependencies)
- [3. Install DisperPicker](#3-install-disperpicker)
- [4. Compile DSurfTomo](#4-compile-dsurftomo)
- [5. Install MATLAB (for EGFAnalysisTimeFreq)](#5-install-matlab-for-egfanalysistimefreq)
- [6. Install GMT (Optional, for Plotting)](#6-install-gmt-optional-for-plotting)
- [Verify Installation](#verify-installation)
- [Troubleshooting](#troubleshooting)

---

## System Requirements

| Item | Requirement |
|------|-------------|
| OS | Linux (Ubuntu 18.04+ recommended) or macOS |
| Python | 3.8 – 3.10 (3.10 recommended) |
| Fortran compiler | gfortran (required to compile DSurfTomo) |
| MATLAB | R2018a or later (required for EGFAnalysisTimeFreq) |
| Disk space | ≥ 5 GB (including Conda environments and packages) |
| RAM | ≥ 4 GB |

---

## 1. Install Conda (if not already installed)

Conda is a Python environment manager that lets you create isolated environments for different projects, avoiding dependency conflicts.

If you already have Conda (Anaconda or Miniconda) installed, skip to Step 2.

### Check if Conda is installed

```bash
conda --version
```

If it prints a version number (e.g. `conda 24.x.x`), you're good — skip to Step 2.

### Install Miniconda (lightweight, recommended)

```bash
# Download the installer (Linux)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Run the installer (follow prompts — press Enter/yes when asked)
bash Miniconda3-latest-Linux-x86_64.sh

# After installation, restart your terminal or run:
source ~/.bashrc
```

macOS (Intel):
```
https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
```

macOS (Apple Silicon — M1/M2/M3):
```
https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
```

---

## 2. Install NoisePy and Python Dependencies

Most Python scripts in this project run in a single Conda environment. We install NoisePy and all helper-script dependencies together.

### 2a. Create the environment and install NoisePy

```bash
# Create an environment named "noisepy" with Python 3.10
conda create -n noisepy -y python=3.10 pip

# Activate it (you need to do this every time you open a new terminal)
conda activate noisepy

# Install NoisePy
pip install noisepy-seis
```

> **What is NoisePy?**
> NoisePy is a Python package for computing ambient noise cross-correlation functions.
> Citation: Jiang, C., & Denolle, M. A. (2020). *Seismological Research Letters*, 91(3), 1853–1866.
> GitHub: https://github.com/noisepy/NoisePy

### 2b. Install helper-script dependencies

While still in the `noisepy` environment:

```bash
conda activate noisepy

# Seismic data processing
pip install obspy

# Data processing and visualization
pip install numpy pandas scipy matplotlib

# YAML config parsing
pip install pyyaml

# Time range handling (needed by local_cc_workflow.py)
pip install datetimerange

# HDF5 file I/O (needed by export_egfanalysis_input.py)
pip install h5py
```

Or install everything in one command:

```bash
pip install obspy numpy pandas scipy matplotlib pyyaml datetimerange h5py
```

### 2c. (Optional) Install MPI support

For large datasets, MPI enables parallel processing:

```bash
# First install the MPI library
# Ubuntu/Debian:
sudo apt install libopenmpi-dev

# macOS:
brew install open-mpi

# Then install NoisePy with MPI support
pip install noisepy-seis[mpi]
```

### 2d. (Optional) Install Jupyter Notebook

```bash
pip install ipykernel notebook
python -m ipykernel install --user --name noisepy
```

---

## 3. Install DisperPicker

DisperPicker uses deep learning (TensorFlow) to automatically extract dispersion curves from spectrograms. It requires a **separate** Python environment because it depends on TensorFlow 1.x, which is incompatible with the packages used by NoisePy.

> **What is DisperPicker?**
> A CNN-based tool for automatic dispersion curve extraction.
> Citation: Yang, S., et al. (2022). *Seismological Research Letters*.
> GitHub: https://github.com/ShaoboYang-USTC/DisperPicker

### 3a. Create a separate environment

```bash
# DisperPicker requires TensorFlow 1.13, which needs Python 3.6
conda create -n disperpicker -y python=3.6 pip
conda activate disperpicker
```

### 3b. Install dependencies

```bash
cd DisperPicker

# Option 1: Use the requirements file (recommended)
pip install -r requirements.txt

# Option 2: Install manually
pip install tensorflow==1.13.1 matplotlib==3.3.4 numpy==1.19.0
```

> **Note:** TensorFlow 1.13.1 does not support Python 3.7 or later. If installation fails, verify your Python version is 3.6.

### 3c. Verify installation

```bash
conda activate disperpicker
cd DisperPicker
python -c "import tensorflow as tf; print('TensorFlow version:', tf.__version__)"
```

Expected output: `TensorFlow version: 1.13.1`

### 3d. Pre-trained model

The pre-trained CNN model is already included in the repository under `saver/` — no extra download needed:

```
DisperPicker/saver/
├── -10000.data-00000-of-00001
├── -10000.index
└── -10000.meta
```

---

## 4. Compile DSurfTomo

DSurfTomo is a Fortran program for surface wave tomography. It must be compiled from source.

> **What is DSurfTomo?**
> DSurfTomo directly inverts surface wave dispersion data for 3D shear-wave velocity models, without the intermediate step of constructing phase/group velocity maps.
> Citation: Fang, H., et al. (2015). *Geophysical Journal International*, 201(3), 1251–1263.
> License: MIT

### 4a. Install the Fortran compiler

```bash
# Ubuntu/Debian
sudo apt update
sudo apt install gfortran

# CentOS/RHEL
sudo yum install gcc-gfortran

# macOS (requires Homebrew: https://brew.sh)
brew install gcc
```

Verify:

```bash
gfortran --version
```

You should see version info, e.g. `GNU Fortran (Ubuntu 9.4.0-...) 9.4.0`.

### 4b. Compile

```bash
cd DSurfTomo

# Run the configure script (it's a Python script that invokes make)
python configure
```

This script will:
1. Create a `bin/` directory
2. Enter `src/` and run `make`
3. Copy the compiled `DSurfTomo` executable to `bin/`

### 4c. Verify

```bash
ls -la DSurfTomo/bin/DSurfTomo
```

If the file exists, compilation succeeded.

### 4d. (Optional) Add to PATH

For convenience, add DSurfTomo to your system path:

```bash
# Add this line to ~/.bashrc or ~/.zshrc
export PATH="/path/to/local_ANT/DSurfTomo/bin:$PATH"

# Apply the change
source ~/.bashrc
```

### 4e. Troubleshooting compilation errors

| Error | Solution |
|-------|----------|
| `gfortran: command not found` | Install gfortran (see 4a) |
| `fatal error: cannot find -lgomp` | Install OpenMP: `sudo apt install libgomp1` |
| `m64: No such file or directory` | Make sure you're on a 64-bit system |

---

## 5. Install MATLAB (for EGFAnalysisTimeFreq)

EGFAnalysisTimeFreq is a MATLAB GUI program for surface wave dispersion analysis.

> **What is EGFAnalysisTimeFreq?**
> Software for surface wave EGF/CF phase velocity dispersion analysis.
> Author: Huajian Yao (University of Science and Technology of China)
> Citation: Yao, H., et al. (2006). *Geophysical Journal International*, 166(2), 732–744.

### 5a. Install MATLAB

- Download and install from the [MathWorks website](https://www.mathworks.com/)
- A valid MATLAB license is required (universities typically provide campus-wide licenses)
- Minimum version: R2018a

### 5b. Verify

Open MATLAB and run:

```matlab
cd('/path/to/local_ANT/EGFAnalysisTimeFreq')
EGFAnalysisTimeFreq
```

If the GUI window opens correctly, the installation is successful.

> **Tip:** For detailed usage instructions, refer to `EGFAnalysisTimeFreq/EGFAnalysisTimeFreq_Manual.pdf`

---

## 6. Install GMT (Optional, for Plotting)

GMT (Generic Mapping Tools) is used by the DSurfTomo plotting scripts for depth slices, cross-sections, ray paths, and checkerboard test results. Skip this if you don't need these plots.

### 6a. Install GMT 6

```bash
# Ubuntu/Debian
sudo apt update
sudo apt install gmt gmt-dcw gmt-gshhg

# CentOS/RHEL
sudo yum install GMT

# macOS
brew install gmt
```

### 6b. Verify

```bash
gmt --version
```

Expected: `6.x.x`. The GMT scripts in this project require GMT 6.0 or later.

---

## Verify Installation

Run the following script to check that all components are properly installed:

```bash
echo "=== 1. NoisePy environment ==="
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
echo "=== 2. DisperPicker environment ==="
conda activate disperpicker
python -c "
import tensorflow as tf; print('TensorFlow', tf.__version__, 'OK')
import numpy; print('NumPy OK')
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
echo "=== 5. GMT (optional) ==="
gmt --version 2>/dev/null || echo "GMT not installed (optional)"
```

---

## Troubleshooting

### Q: `conda activate noisepy` says "conda: command not found"
Make sure Conda is properly installed. Run `source ~/.bashrc` (or `source ~/.zshrc`) and try again. If it still doesn't work, check that the Conda installation directory is in your `$PATH`.

### Q: `pip install noisepy-seis` is slow or times out
Use a mirror (Chinese users):
```bash
pip install noisepy-seis -i https://pypi.tuna.tsinghua.edu.cn/simple
```

### Q: DisperPicker — TensorFlow 1.13.1 installation fails
TensorFlow 1.13.1 only supports Python 3.6. Check your version:
```bash
python --version  # Should show 3.6.x
```
If it's wrong, recreate the environment: `conda create -n disperpicker -y python=3.6 pip`

### Q: TensorFlow prints GPU-related errors
DisperPicker runs on CPU by default — no GPU required. You can safely ignore CUDA/GPU warnings, or suppress them:
```bash
export CUDA_VISIBLE_DEVICES=""
```

### Q: DSurfTomo compilation fails with `gfortran: command not found`
Install gfortran:
```bash
sudo apt install gfortran      # Ubuntu/Debian
sudo yum install gcc-gfortran  # CentOS/RHEL
brew install gcc                # macOS
```

### Q: Can I merge the two Conda environments into one?
No. NoisePy requires Python 3.10 and recent NumPy versions, while DisperPicker requires Python 3.6 and TensorFlow 1.13. They cannot coexist in the same environment. Simply switch between them:
```bash
conda activate noisepy        # For NoisePy-related scripts
conda activate disperpicker   # For DisperPicker
```

### Q: Running DSurfTomo scripts gives `ModuleNotFoundError: No module named 'numpy'`
The Python scripts under `DSurfTomo/scripts/` should be run in the `noisepy` environment:
```bash
conda activate noisepy
cd DSurfTomo/scripts
python PrepareDSurfTomoProject.py
```

---

## Environment Quick Reference

| Step | Environment |
|------|-------------|
| Step 1: Data cleaning | `conda activate noisepy` |
| Step 2: Prepare station file | No Python needed |
| Step 3: Cross-correlation | `conda activate noisepy` |
| Step 4: Export EGF data | `conda activate noisepy` |
| Step 5: MATLAB dispersion | MATLAB |
| Step 6: Prepare DisperPicker input | `conda activate noisepy` |
| Step 7: DisperPicker extraction | `conda activate disperpicker` |
| Step 8: DSurfTomo inversion | `conda activate noisepy` (Python scripts) + `DSurfTomo/bin/DSurfTomo` (executable) |
| GMT plotting | System terminal (no Python environment needed) |
