# local_ANT

**密集台阵背景噪声面波成像：从连续记录到三维横波速度模型。**

[![tests](https://github.com/Shubin-Wei/local_ANT/actions/workflows/tests.yml/badge.svg)](https://github.com/Shubin-Wei/local_ANT/actions/workflows/tests.yml)
[![Python](https://img.shields.io/badge/python-3.11%20%7C%203.12%20%7C%203.13-blue)](pyproject.toml)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)

中文 | [English](README.md)

lant（`local_ANT`）是用于密集（节点）台阵（面状或测线）高频背景噪声成像的命令行工具。它只做文献中的基本处理流程：

- **通用**：ObsPy 能读的任何记录格式、任何目录结构、任何台阵；不为某一批数据加功能。
- **易用**：一个配置文件，每个设置都有注释；`lant init` 按你的台阵给出建议值。
- **透明**：每一步写出所用的设置、推出的每个值及其规则、日志和图。中间结果是普通的 HDF5 和 CSV 文件；反演工程是 DSurfTomo 的格式。
- **简洁**：三步，一步一步运行和检查。
- **正确**：每个方法都有出处（[方法](docs/methods.zh.md)）；单元测试用精确解检验公式。

## 步骤

| 命令 | 做什么 | 依据 |
|---|---|---|
| `lant cc` | 读记录，计算每个台对的互相关（重叠分段的互相干），叠加 | Bensen et al. (2007)；Seats et al. (2012)；Xu et al. (2021) |
| `lant disp` | 参考曲线（台阵 F-J 谱、一维模型或曲线文件）；两台法相速度（相位匹配滤波、刀切误差、质量控制） | Wang et al. (2019)；Levshin & Ritzwoller (2001)；Lin et al. (2008)；Luo et al. (2015) |
| `lant inv` | 一维起始模型；基于射线追踪的走时三维直接反演 | Fang et al. (2015)（DSurfTomo 的方法，Python 实现） |
| `lant checker` | 棋盘分辨率检验 | Lévêque et al. (1993)；Rawlinson & Spakman (2016) |

## 安装

Linux，需要 conda（或 mamba）和 git：

```bash
git clone https://github.com/Shubin-Wei/local_ANT.git
cd local_ANT
conda env create -f environment.yml
conda activate lant
lant check
```

地图由 GMT 6 绘制（环境里已包含）；没有 GMT 时其余功能照常。

## 快速开始

```bash
lant init -o my.yaml --raw /data/raw --stations stations.txt --out /data/out --velocity 0.3 3.0
# 先读一遍 my.yaml，然后：
lant cc   -c my.yaml      # 检查 1_cc/figures，然后
lant disp -c my.yaml      # 尤其检查 2_disp/figures/fj_reference.png，然后
lant inv  -c my.yaml
lant checker -c my.yaml --size 2
```

台站表需要表头：`network station longitude latitude elevation`（可以没有 network 列）。`--velocity` 是本区 Rayleigh 波相速度的大致最慢和最快值。

每一步之后看什么：[教程](docs/tutorial.zh.md)。

## 输出

```
OUT/1_cc/       stack.h5, days/*.h5, stations_used.csv, record_index.csv
OUT/2_disp/     dispersion.csv, measurements.h5, reference.csv, median_vs_reference.csv, fj.h5
OUT/3_inv/      vs_model.csv, start_model.csv, residuals.csv, dsurftomo/
OUT/4_checker_*/recovery.csv
每一步还有：params.yaml, log.txt, summary.txt, figures/
```

## 文档

- [教程](docs/tutorial.zh.md)：每一步怎样运行、检查什么
- [方法](docs/methods.zh.md)：原理、公式、出处和局限
- [参数](docs/parameters.zh.md)：每个设置、默认值及依据
- [版本记录](docs/history.zh.md)

## 引用

请引用所用方法的文献（见[方法](docs/methods.zh.md)），尤其是三维反演的 Fang et al. (2015)，以及本软件（[CITATION.cff](CITATION.cff)）。

## 许可

MIT。反演部分（`lant/inversion/dsurf_*.py`）是 DSurfTomo（Fang et al. 2015；MIT 许可）的 Python 移植。
