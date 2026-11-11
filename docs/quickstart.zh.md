# 快速上手

中文 | [English](quickstart.md)

本文介绍安装、在已知答案的合成数据上试运行，以及用自己的数据建立工程，全程约半小时。每一步的细节见[教程](tutorial.zh.md)，原理见[方法](methods.zh.md)。

## 1. 安装（约 10 分钟）

环境要求：Linux、conda（或 mamba）、git。

```bash
git clone https://github.com/Shubin-Wei/local_ANT.git && cd local_ANT
conda env create -f environment.yml && conda activate lant   # 同时从本仓库安装 lant
lant install-external                                        # 下载并编译 DSurfTomo 和 pDSurfTomo
lant check                                                   # 各项均应显示正常
```

`lant check` 会对缺失项给出安装命令。

## 2. 用演示检验安装（约 5 分钟）

```bash
lant demo --out ~/lant_demo
```

演示生成一个小型合成台阵（已知一维模型中的 Rayleigh 平面波），运行完整流程，并给出报告位置 `~/lant_demo/output/REPORT.md`。演示会自行确认判断点；这只在数据由它自己生成时才成立，因此它拒绝在不是由它创建的目录中运行。

## 3. 逐步运行

如需了解流程在哪些地方停下等待使用者判断，可生成同样的数据并手动运行：

```bash
python -m lant.synthetic --quick --out ~/lant_quick          # 10 个台站、6 小时记录，约 1 分钟
lant run all -c ~/lant_quick/synthetic.yaml
```

流程在需要科学判断处停止，并给出应查看的文件和继续运行的命令。合成数据的答案已知；**处理实际数据时，每一项都须由你查看图件后决定**。

| 步骤 | 判断事项 | 依据 | 合成数据上的答案 |
|---|---|---|---|
| `clean` | 未提供响应文件时：所有台站是否为同一型号仪器？ | `work/01_clean/report.txt` | 是 |
| `stack`（列出候选时） | 是否剔除所列的台站子叠加？ | `work/03_stack/station_qc/station_qc.png` | 没有故障台站，保留 |
| `reference` | F-J 拾取是否沿基阶脊？ | `work/04_reference/fj.png` | 是 |
| `inversion` | 一维起始模型及其对平均曲线的拟合是否合理？ | `results/figures/4_start_model.png` | 是 |

记录每项决定后重新运行，流程从停止处继续：

```bash
lant confirm reference -c ~/lant_quick/synthetic.yaml --reason "已知答案的合成数据"
lant run all -c ~/lant_quick/synthetic.yaml
lant status -c ~/lant_quick/synthetic.yaml                     # 各步骤状态及下一条命令
```

流程结束后与真值比较：

```bash
python examples/synthetic/compare_truth.py --dir ~/lant_quick
```

输出各频率上实测与真实相速度之差，以及反演模型与真实模型之差。6 小时记录得到的数字比完整算例粗略（完整算例为 16 个台站、3 天，约 300 MB，处理约 1 小时；去掉 `--quick` 即可生成），见[验证](validation.zh.md)。

## 4. 处理自己的数据

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

然后开始运行，停止点与第 3 节相同：

```bash
lant run all -c /data/result/config.yaml
```

`lant run all --raw /data/raw --stations /data/stations.csv --out /data/result` 会写出同样的配置并立即开始运行。运行中遇到的问题见[常见问题](troubleshooting.zh.md)。

## 5. 运行之后

先查看输出目录中的 `REPORT.md`：各步骤状态、关键数字、待你判断的事项和主要图件，以及成果文件的链接。然后检查：

| 文件 | 用途 |
|---|---|
| `.lant/steps/<步骤>.md` | 每一步做了什么、使用的参数（数值和来源）、留给你判断的事项 |
| `work/01_clean/coverage.png` | 各台站各天的数据量 |
| `work/04_reference/fj.png` | 参考曲线是否沿基阶 |
| `results/figures/3_dispersion_curves.png` | 全部频散曲线及参考曲线 |
| `results/figures/5_depth_slices.png`、`5_sections.png`、`6_checkerboard_*.png` | 模型及其可分辨的尺度，须结合起来看 |

各输出文件的说明见[输出文件](outputs.zh.md)。
