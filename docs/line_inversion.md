# 06a 线状台阵速度结构反演 / Linear-Array Inversion

## 适用场景

台站近似排成一条测线(如浅层勘探测线、公路沿线布设)。此时二维/三维层析(DSurfTomo)缺乏垂直测线方向的射线交叉,不适用;改用"逐台站对一维反演 + 沿测线拼接"的策略:

1. **01_invert.py**:对每个台站对的频散曲线,用 [evodcinv](https://github.com/keurfonluu/evodcinv)(CPSO 演化算法,正演内核 disba)反演一维 Vs(depth) 模型,结果挂到该台站对的**中点**;
2. **02_stitch.py**:把全部中点用 PCA 投影到测线主方向得到沿线距离,按距离分箱平均、插值、平滑,拼成**伪二维 Vs 深度剖面**。

> 本分支与 DSurfTomo 无任何依赖关系;反演深度大致受最大周期约束(经验上灵敏深度 ≈ 波长/3 ≈ c·T/3)。

## 运行环境

独立的 evodcinv Conda 环境(安装见 [INSTALL.md](INSTALL.md) §2e):

```bash
conda env create -f environments/evodcinv.yml   # 首次创建
conda activate evodcinv
```

## 输入

第 05 步的 `qc_result/` 目录,直接使用、无需转换:

- `CDisp.<台A>_<台B>.txt` — 相速度曲线;
- `GDisp.<台A>_<台B>.txt` — 群速度曲线;
- 文件头两行是两台站经纬度,数据列为 `周期 速度 信噪比 概率`(速度 ≤0 的点视为无效)。

另需台站坐标 CSV(含高程,用于剖面顶部地形),格式同 02 步的 `station.csv`。

## 配置

模板:`configs/06a_line_inversion/line.example.yaml`,一条测线一份;`01_invert.py` 与 `02_stitch.py` 共用同一份。三段:

**`paths:`** — `input_dir`(qc_result 路径)、`station_csv`、`output_dir`。

**`inversion:`** — `01_invert.py` 使用:

| 键 | 含义 |
|----|------|
| `curve_type` | `phase`(用 CDisp)/ `group`(GDisp)/ `both`(联合) |
| `workers` | 并行进程数(每进程单线程,已自动钉死 BLAS 线程) |
| `min_points` | 有效点少于该数的曲线跳过 |
| `plot` | 是否输出每对的反演图 |
| `zmax` / `dz` | 输出深度网格(km);zmax 建议 ≈ c·Tmax/3 |
| `layers` | 分层搜索空间:每层 `{thickness: [min,max], vs: [min,max]}`(km, km/s),厚度上限之和应≈参数化深度 |
| `optimizer` | evodcinv CPSO 设置:`popsize`、`maxiter`、`maxrun` |

**`stitch:`** — `02_stitch.py` 使用:`misfit_max`(剔除失配大的模型)、`nbins`(沿线分箱数)、`smooth`(高斯平滑 σ)、`contour_step`、`cmap`/`vmin`/`vmax`、`topography`(剖面顶随地形)、`show_stations`、`vertical_exaggeration`、`taper_slope`(倒梯形掩膜)、`depth_guides`、`dpi`/`save_pdf` 等,逐项注释见模板。

## 运行

```bash
cp configs/06a_line_inversion/line.example.yaml configs/06a_line_inversion/my_line.yaml
# 编辑 my_line.yaml

python scripts/06a_line_inversion/01_invert.py --config configs/06a_line_inversion/my_line.yaml
#   可选:--max-files N  只跑前 N 对(冒烟测试)

python scripts/06a_line_inversion/02_stitch.py --config configs/06a_line_inversion/my_line.yaml
#   可选:--output-dir DIR  覆盖 paths.output_dir
```

多条测线:复制多份 yaml(如 `line1_phase.yaml`、`line1_group.yaml`),分别运行即可;相/群速度或不同 `curve_type` 也用独立配置区分。

## 输出

```
<output_dir>/
├── models/         # 每台站对的一维 Vs 模型(挂在中点)
├── figs/           # 每对的反演图(plot: true 时)
├── logs/           # 每个并行 worker 的日志
├── summary.txt     # 全部反演的汇总(含 misfit)
├── section.png     # 伪二维 Vs 深度剖面
├── section.pdf     # save_pdf: true 时
└── section.npz     # 剖面数值结果(距离、深度、Vs 网格)
```

反演失败的台站对会在日志与汇总中标记;若模型已保存但事后诊断失败,会保留模型并打印警告(不中断批处理)。

## 常见调参建议

- 剖面出现异常条带:先收紧 `misfit_max`,再增大 `smooth`;
- 深部不可信:减小 `zmax`(受最大周期限制,宁浅勿深);
- 反演太慢:减小 `optimizer.maxiter`/`popsize`,或增大 `workers`;
- 曲线太少:检查 05 步拾取阈值,或把 `curve_type` 从 `both` 放宽为单一类型。
