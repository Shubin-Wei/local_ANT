# local_ANT

**密集台阵背景噪声面波成像：从连续记录到三维横波速度模型。**

[![tests](https://github.com/Shubin-Wei/local_ANT/actions/workflows/tests.yml/badge.svg)](https://github.com/Shubin-Wei/local_ANT/actions/workflows/tests.yml)
[![Python](https://img.shields.io/badge/python-3.11%20%7C%203.12%20%7C%203.13-blue)](pyproject.toml)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)

中文 | [English](README.md)

lant（`local_ANT`）是用于密集（节点）台阵（面状或测线）高频背景噪声成像的命令行工具。输入垂直分量连续记录和台站表，输出台对间的 Rayleigh 波相速度和三维 Vs 模型。各步分开运行；每一步写出所用的设置、推出的值及其规则、日志和图，检查后再做下一步。

## 步骤

| 命令 | 做什么 | 主要依据 | 结果 |
|---|---|---|---|
| `lant clean` | 按台站把记录切成 UTC 日；去中位数，给了响应文件时去仪器响应；低通后降采样 | Bensen et al. (2007) | `1_clean` |
| `lant cc` | 分段（去掉含缺数和振幅过大的段）、每个台对的互相干、线性叠加 | Bensen et al. (2007)；Seats et al. (2012)；Xu et al. (2021) | `2_cc` |
| `lant disp` | 参考曲线（台阵 F-J 谱、一维模型或曲线文件）；两台法相速度、刀切误差和质量控制 | Wang et al. (2019)；Levshin & Ritzwoller (2001)；Lin et al. (2008) | `3_disp` |
| `lant inv` | 一维起始模型；基于射线追踪的走时三维直接反演 | Fang et al. (2015)（DSurfTomo 的 Python 实现） | `4_inv` |
| `lant checker` | 棋盘分辨率检验 | Rawlinson & Spakman (2016) | `5_checker_*` |

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

需要：

- 连续记录，ObsPy 能读的任何格式（miniSEED、SAC 等），任何目录结构；
- 带表头的台站表 `network station longitude latitude elevation`（可以没有 network 列），或 StationXML；台站名须与记录头里的一致；
- 本区 Rayleigh 波相速度大致的最慢和最快值（km/s）。

```bash
lant init -o my.yaml --raw /data/raw --stations stations.txt --out /data/out --velocity 0.3 3.0
# 先读一遍 my.yaml（lant init 写入了建议值及其规则），然后：
lant clean -c my.yaml     # 检查 1_clean/figures
lant cc    -c my.yaml     # 检查 2_cc/figures
lant disp  -c my.yaml     # 尤其检查 3_disp/figures/fj_reference.png
lant inv   -c my.yaml
lant checker -c my.yaml --size 2
```

每一步之后看什么：[教程](docs/tutorial.zh.md)。

## 结果

```
OUT/1_clean/       days/*.h5（清洗后的台站日）, availability.csv, stations_used.csv, record_index.csv
OUT/2_cc/          stack.h5, days/*.h5（子叠加）
OUT/3_disp/        dispersion.csv, measurements.h5, reference.csv, median_vs_reference.csv, fj.h5
OUT/4_inv/         vs_model.csv, start_model.csv, residuals.csv, dsurftomo/
OUT/5_checker_*/   recovery.csv
每一步还有：       params.yaml, log.txt, summary.txt, figures/
```

大数组用 HDF5，其余是 CSV；`lant export-sac` 把叠加后的互相关写成 SAC 文件。`4_inv/dsurftomo/` 是 DSurfTomo 格式的工程。

## 文档

- [教程](docs/tutorial.zh.md)：每一步怎样运行、检查什么
- [方法](docs/methods.zh.md)：原理、公式、出处和局限
- [参数](docs/parameters.zh.md)：每个设置、默认值及依据
- [版本记录](docs/history.zh.md)

## 引用

请引用所用方法的文献（见[方法](docs/methods.zh.md)），尤其是三维反演的 Fang et al. (2015)，以及本软件（[CITATION.cff](CITATION.cff)）。

## 许可

MIT。反演部分（`lant/inversion/dsurf_*.py`）是 DSurfTomo（Fang et al. 2015；MIT 许可）的 Python 移植。
