# 快速上手

中文 | [English](quickstart.md)

这份说明用大约半小时带你走完一遍：安装，在一套已知答案的小合成数据上跑完整个流程，再换成你自己的数据。每一步的原理见 [methods.zh.md](methods.zh.md)，详细操作见 [tutorial.zh.md](tutorial.zh.md)。

## 1. 安装（约 10 分钟）

需要 Linux、conda（或 mamba）、git。

```bash
git clone https://github.com/Shubin-Wei/local_ANT.git && cd local_ANT
conda env create -f environment.yml && conda activate lant     # also installs lant from this checkout
lant install-external      # 下载并编译 DSurfTomo、pDSurfTomo
lant check                          # 全部显示 ok 即可；缺什么会给出安装命令
```

## 2. 在合成数据上试一遍（约 10 分钟）

合成数据由一个已知一维模型的 Rayleigh 平面波噪声场生成，因此可以把结果和真值比较。

**一条命令：**

```bash
lant demo --out ~/lant_demo        # 生成数据、跑完整条流程、打印报告位置
```

demo 会自己确认判断点，这只因为数据是它自己生成的；它拒绝在不是它创建的目录里运行。完成后打开 `~/lant_demo/output/REPORT.md`。

**逐步运行**，看流程在哪里停下来等你判断：

```bash
python -m lant.synthetic --quick --out ~/lant_quick     # 10 台、6 小时，约 1 分钟
lant run all -c ~/lant_quick/synthetic.yaml
```

流程会在需要科学判断的地方停下，屏幕上打印要看的文件和下一条命令。合成数据的答案已知，可以直接确认；**在实际数据上，每一处都要你自己看图后再决定**。

| 停在 | 要回答的问题 | 看什么 | 合成数据上 |
|---|---|---|---|
| clean | 没有响应文件，全部台站是不是同一型号？ | `work/01_clean/report.txt` | 是：`lant confirm clean -c … --reason "synthetic data with a known answer"` |
| stack（有候选时） | 列出的台站子叠加要不要剔除？ | `work/03_stack/station_qc/station_qc.png` | 合成数据没有坏台站，保留：`lant confirm stack -c … --reason "synthetic data with a known answer"` |
| reference | F-J 拾取是不是基阶脊线？ | `work/04_reference/fj.png` | 是：`lant confirm reference -c … --reason "synthetic data with a known answer"` |
| inversion | 一维起始模型和它对平均曲线的拟合是否合理？ | `results/figures/4_start_model.png` | 是：`lant confirm inversion -c … --reason "synthetic data with a known answer"` |

每次确认后再运行一次 `lant run all -c ~/lant_quick/synthetic.yaml`，它会从停下的地方接着跑。任何时候都可以用

```bash
lant status -c ~/lant_quick/synthetic.yaml     # 每一步的状态、原因和下一条命令
```

跑完后与真值比较：

```bash
python examples/synthetic/compare_truth.py --dir ~/lant_quick
```

它给出各频率上实测相速度与真实曲线之差、反演模型与真实模型之差。这套快速数据只有 6 小时记录，数字比完整算例粗（完整算例的数字见 [validation.zh.md](validation.zh.md)）；完整算例去掉 `--quick` 即可（16 台、3 天，约 300 MB，处理约一小时）。

## 3. 处理你自己的数据

只需要两样东西：

- **原始数据目录**：任意子目录结构和文件名，ObsPy 能读的格式（miniSEED、SAC 等）都可以，可以混放。
- **台站表**：CSV、文本或 Excel，至少有台站名、经度、纬度三列（列名中英文都认）。节点搬过站时再加序列号和布设起止时间（[tutorial.zh.md](tutorial.zh.md) 2.2 节）。

```bash
lant run all --raw /data/raw --stations /data/stations.csv --out /data/result
```

这条命令写出 `/data/result/config.yaml` 并开始运行。开始前建议在配置里改两处：

- `data.velocity_range_kms`：本区面波相速度的大致范围。默认 0.1–5 km/s 覆盖软土到硬岩，但时延窗与最慢速度成反比，范围很宽时互相关的计算量和磁盘占用都大得多。
- `data.response`：有仪器响应文件就写上；没有时全部台站必须为同一型号。

改完后 `lant run all -c /data/result/config.yaml`。之后的停点与合成算例相同，但每一处都要你看图判断；可参照 [troubleshooting.zh.md](troubleshooting.zh.md) 处理遇到的问题。结果文件的含义见 [outputs.zh.md](outputs.zh.md)。

## 4. 每次运行后建议看的文件

| 文件 | 为什么 |
|---|---|
| `.lant/steps/<步骤>.md` | 每一步做了什么、用了哪些参数（值和依据）、留给你判断的事项 |
| `work/01_clean/coverage.png` | 每个台站有多少数据 |
| `work/04_reference/fj.png` | 参考曲线是否沿基阶 |
| `results/figures/3_dispersion_curves.png` | 全部频散曲线与参考曲线 |
| `results/figures/5_depth_slices.png`、`5_slices_3d.png`（线阵看 `5_sections.png`）、`6_checkerboard_*.png` | 模型及其可分辨的尺度（两者一起看） |
