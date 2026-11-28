# 教程：从记录到三维模型，一步一步来

中文 | [English](tutorial.md)

按顺序一步一步运行，每一步看完它的图和日志再做下一步。重跑某一步只覆盖它自己的文件夹，它后面的步骤也要跟着重跑。

每一步写到输出目录下自己的文件夹：

| 文件 | 内容 |
|---|---|
| `params.yaml` | 这一步用的设置，以及它推出的每个值和推出规则 |
| `log.txt` | 这一步打印的全部内容（带时间） |
| `summary.txt` | 关键数字 |
| `figures/` | 要检查的图 |

## 0. 输入和配置

需要：

- 连续记录的文件夹，ObsPy 能读的任何格式（miniSEED、SAC 等），任何目录结构；
- 带表头的台站表，例如

  ```
  network station longitude latitude elevation
  XX S001 100.12345 30.12345 1520
  ```

  （逗号或空格分隔；可以没有 network 列；也可以是 StationXML）。名字必须与记录头里的代码一致：lant 严格对应，对不上的会列出来；
- 本区 Rayleigh 波相速度大致的最慢和最快值（km/s），例如软沉积 0.2–1.5，硬岩 1.5–3.5。

```bash
lant init -o my.yaml --raw /data/raw --stations stations.txt --out /data/out --velocity 0.3 3.0
```

`lant init` 读记录头，写出 `my.yaml`，给出采样率、互相关频带、时延窗和段长的建议值，每个建议值上方写明规则和代入的本台阵数字（例如 `2 x v_fast 3 km/s / longest pair 9.8 km = ...`）。打开文件读一遍，每个设置都有注释；你更清楚的就改（例如台站仪器不全相同时填 `data.response`）。[`example.yaml`](../example.yaml) 是 `lant init` 写出的一份完整配置示例。

## 1. 数据清洗：`lant clean -c my.yaml`

做什么：建立记录索引，与台站对应，把每个台站的每个 UTC 日清洗后写出：缺数保持为缺数，去中位数，去仪器响应（给了 `data.response` 时），低通后降采样到 `clean.fs_hz`（方法第 2 节）。这一步不按振幅丢数据，按段剔除在 `lant cc` 里做。

检查：

- **`log.txt`**：台站对应情况（没有记录的台站、表外台站的记录、有多个传感器的台站）、不是地震记录的文件、复用的天。
- **`figures/stations.png`**：台站是否齐全、位置是否正确。
- **`figures/availability.png`** 和 `availability.csv`：每个台站每天有多少小时数据。开始晚、结束早、缺数多的台站在这里能看出来。

结果：`days/YYYY-MM-DD.h5`（每天一个数组：台站 × 样点，没有有效样点处为 NaN）、`stations_used.csv`、`availability.csv`、`record_index.csv`。只有原始文件或设置变了的天才会重新清洗。

## 2. 互相关：`lant cc -c my.yaml`

做什么：读清洗后的数据，分段，去掉含缺数或振幅过大的段，计算每个台对的互相干，叠加（方法第 3 节）。`data.exclude_stations` 里的台站和 `data.start_date` / `end_date` 以外的天直接不用，不必重新清洗。

检查：

- **`log.txt`**：用到的台站、台对和天数，复用的天。
- **`figures/ccf_<频段>Hz.png`**：互相关随台距排列，互相关频带分成几个子频段各画一张（上为图像，下为按台距分组平均的波形）。面波表现为两条虚线（速度范围）之间的时距同相轴。噪声来自各个方向时零时延两侧都有；来自一个方向时主要在一侧。
- **`figures/ccf_time_*.png`**：三个台对：叠加波形及其信号窗（灰）和噪声窗（蓝），按时间排列的子叠加，以及前 n 个子叠加的叠加信噪比。到时应稳定，信噪比应随 n 增大。某个子叠加明显不同，说明当天某台有问题（钟差、仪器）或有强的局部源。可用 `data.exclude_stations` 去掉该台，或用 `data.start_date` / `end_date` 去掉这些天，再重跑这一步。

结果：`stack.h5`（叠加）、`days/*.h5`（每天的子叠加）。`lant export-sac -c my.yaml -o 目录` 把叠加写成 SAC 文件，供其他软件使用。

## 3. 相速度和群速度：`lant disp -c my.yaml`

做什么：参考曲线，然后每个台对的相速度、群速度及其刀切误差，质量控制（方法第 4–5 节）。

按顺序检查：

1. **`figures/fj_reference.png`**：全流程最重要的检查。黑线（参考曲线）必须在**基阶**脊上，即每个频率最慢的那条脊。如果某处跳到了高阶：
   - 把拾取限制在基阶清楚的频段：`disp.fj_band_hz: [下限, 上限]`；或
   - 给一维模型（`disp.reference: model`，`disp.reference_file: model.txt`，两列 `depth_km Vs_km/s`），或人工拾取的曲线（`disp.reference: file`，两列 `frequency_Hz c_km/s`），

   然后重跑 `lant disp`。紫色菱形是实测相速度的中位数，应落在参考曲线上或其附近。

   测线台阵的 F-J 谱常常不清楚（台对少、台距分布不均）。这时用附近面阵或文献给出的一维模型或曲线（`disp.reference: model` 或 `file`）。
2. **`median_vs_reference.csv`**：偏离达 20 % 以上，或在相邻频率间跳动，说明整周数有错（参考曲线有错）；稳定的几个百分点偏离是与参考的真实差别。
3. **`figures/phase_curves.png`**、**`figures/group_curves.png`**：所有接受的曲线及其中位数和标准差；**`figures/rejections.png`**：每个频率接受了多少台对，其余的被拒原因（原因名称见方法第 5 节）。
4. **`figures/pair_*.png`**：三个示例台对：互相关及两个速度之间的窗口，相速度图像（波峰对应各种可能的整周数）和群速度图像。接受的相速度应落在同一条波峰上，群速度应在能量极大上。

不要为了多出点而放宽质控：`disp` 里的数字及其依据见 [parameters.zh.md](parameters.zh.md)。

结果：`dispersion.csv`（接受的值，每行一个台对和频率，`c_kms` 和 `u_kms` 两列，被拒的一项留空；是 `lant inv` 的输入）、`measurements.h5`（全部值及每个被拒值的原因）、`reference.csv`、`fj.h5`。

## 4. 三维反演：`lant inv -c my.yaml`

做什么：平均曲线、深度节点和一维起始模型、网格，以及相速度走时、群速度走时或两者的三维直接反演（`inv.data`；方法第 6 节）。推出的数字（深度节点、网格间距、Vs 上下限）会打印出来并记在 `params.yaml` 里；要改哪个，就写进 `inv` 再重跑。

检查：

- **`figures/start_model.png`**：起始模型是否拟合平均曲线。
- **`figures/residuals.png`** 和 `residuals.csv`：最终模型的残差应比起始模型小；按周期的均方根说明哪里拟合较差。
- **`figures/paths.png`**：数据覆盖在哪里。
- **`figures/slices.png`**、**`figures/sections.png`**：模型的深度切片，以及沿切片上所画两条剖面线的剖面。只在分辨率检验能恢复异常的地方解释。

结果：`vs_model.csv`（经度、纬度、深度、Vs、相对起始模型的 dVs（%）、DWS：0 表示该节点没有数据）、`start_model.csv`、`average_curve.csv`、`residuals.csv`，以及 `dsurftomo/`（DSurfTomo 格式的工程，原版 DSurfTomo 也能直接运行）。

## 5. 分辨率检验：`lant checker -c my.yaml --size 公里数`

格子边长 `--size` 公里、幅度 ±10 %（`checker.amplitude`）的棋盘，噪声与数据相当；`--flip-depth 公里数` 使图样在该深度以下反号。做几个尺寸，例如网格间距的 1、2、3 倍：

```bash
lant checker -c my.yaml --size 1.5
lant checker -c my.yaml --size 3
lant checker -c my.yaml --size 3 --flip-depth 0.5
```

看 `figures/checker.png`（各深度的输入与输出）和 `recovery.csv`（按深度的相关系数和斜率；斜率即恢复的比例）。

## 出错时

每一步出错时会停下，并说明缺什么或哪里不一致、该改什么。完整经过在该步的 `log.txt` 里。数据问题（读不了的文件、对不上的台站）都会列出来，不会悄悄跳过。
