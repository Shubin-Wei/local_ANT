# 快速上手

中文 | [English](quickstart.md)

本文介绍安装、用自己的数据建立工程，以及流程停下等你判断的地方。每一步的细节见[教程](tutorial.zh.md)，原理见[方法](methods.zh.md)。

## 1. 安装（约 10 分钟）

环境要求：Linux、conda（或 mamba）、git。

```bash
git clone https://github.com/Shubin-Wei/local_ANT.git && cd local_ANT
conda env create -f environment.yml && conda activate lant   # 同时从本仓库安装 lant
lant check                                                   # 各项均应显示正常
```

`lant check` 会对缺失项给出安装命令。

反演就是当前环境中的 DSurfTomo Python 移植版，不需要另装求解器。与原生程序的对比见 [python_backends.zh.md](python_backends.zh.md)。

## 2. 准备自己的数据

需要两项输入：

- **原始记录**：一个目录，子目录结构和文件名不限，ObsPy 可读的任意格式（miniSEED、SAC 等），可以混合。
- **台站表**：CSV、文本或 Excel，至少包含台站名、经度、纬度；中英文列名均可识别。节点仪器移动过位置时，补充序列号和布设起止时间（[教程](tutorial.zh.md)第 2.2 节）。

先写出配置并检查，再开始运行：

```bash
lant init -o /data/result/config.yaml --raw /data/raw --stations /data/stations.csv --out /data/result
lant validate -c /data/result/config.yaml      # 检查配置、台站表和原始数据目录
lant plan -c /data/result/config.yaml          # 每一步将执行的内容及所用数值
```

有两项设置值得注意：

- `data.velocity_range_kms`：研究区面波相速度的大致范围。默认 0.1–5 km/s 覆盖软土到硬岩；时延窗与最低速度成反比，范围过宽会增加计算时间和磁盘占用。仅在对研究区有独立认识时才收窄。
- `data.response`：有仪器响应文件时给出；没有时，所有台站须为同一型号仪器。

`lant run all --raw /data/raw --stations /data/stations.csv --out /data/result` 会写出同样的配置并立即开始运行。

## 3. 流程在哪里停下等你判断

```bash
lant run all -c /data/result/config.yaml
```

流程在需要科学判断处停止，并给出应查看的文件和继续运行的命令。**每一项都须由你查看图件后决定。**

| 步骤 | 判断事项 | 依据 |
|---|---|---|
| `clean` | 未提供响应文件时：所有台站是否为同一型号仪器？ | `work/01_clean/report.txt` |
| `stack`（列出候选时） | 是否剔除所列的台站子叠加？ | `results/figures/2_station_qc.png` |
| `reference` | F-J 拾取是否沿基阶脊？ | `results/figures/3_fj_phase_velocity.png` |
| `inversion` | 一维起始模型及其对平均曲线的拟合是否合理？ | `results/figures/4_start_model.png` |

记录每项决定及其依据后重新运行，流程从停止处继续：

```bash
lant confirm reference -c /data/result/config.yaml --reason "拾取沿 3_fj_phase_velocity.png 中的基阶脊"
lant run all -c /data/result/config.yaml
lant status -c /data/result/config.yaml        # 各步骤状态及下一条命令
```

运行中遇到的问题见[常见问题](troubleshooting.zh.md)。

## 4. 运行之后

先查看输出目录中的 `REPORT.md`：各步骤状态、关键数字、待你判断的事项和主要图件，以及成果文件的链接。然后检查：

| 文件 | 用途 |
|---|---|
| `.lant/steps/<步骤>.md` | 每一步做了什么、使用的参数（数值和来源）、留给你判断的事项 |
| `results/figures/0_data_coverage.png` | 各台站各天的数据量 |
| `results/figures/3_fj_phase_velocity.png` | 参考曲线是否沿基阶 |
| `results/figures/3_dispersion_curves.png` | 全部频散曲线及参考曲线 |
| `results/figures/5_depth_slices.png`、`5_sections.png`、`6_checkerboard_*.png` | 模型及其可分辨的尺度，须结合起来看 |

各输出文件的说明见[输出文件](outputs.zh.md)。
