# 参数

中文 | [English](parameters.md)

配置文件里的每个设置、默认值和依据。带同样注释的模板是 `lant/config_template.yaml`（`lant init` 会写出一份填好的）。"约定"表示取值是一种选择，不是理论或文献给出的值，有理由时就改。

## 必填

| 键 | 含义 |
|---|---|
| `output_dir` | 全部结果的文件夹，每步一个子文件夹 |
| `data.raw_dir` | 连续记录的文件夹（ObsPy 能读的任何格式、任何目录结构） |
| `data.station_file` | 带表头的台站表（`network station longitude latitude [elevation]`），或 StationXML |
| `data.velocity_range_kms` | 本区 Rayleigh 波相速度的最慢和最快值（km/s） |
| `clean.fs_hz` | 降采样后的采样率；`lant init` 给建议值（规则见下） |
| `cc.fmin_hz`、`cc.fmax_hz`、`cc.maxlag_s`、`cc.segment_s` | 互相关设置；`lant init` 给建议值（规则见下） |

## data

| 键 | 默认 | 含义与依据 |
|---|---|---|
| `channel` | `"*Z"` | 垂直分量的通道代码（通配符） |
| `response` | `null` | 响应文件（StationXML、RESP、dataless SEED；文件或文件夹）：去响应到速度，预滤波在 `cc.fmin_hz`–`cc.fmax_hz` 内平坦（因此这时 `lant clean` 也要用到频带）。仪器不全相同时需要；同一型号时互相干已使其不必要（方法第 2 节） |
| `start_date`、`end_date` | `null` | 使用的第一天和最后一天（UTC）（`lant clean` 和 `lant cc` 都用：去掉某些天不必重新清洗） |
| `exclude_stations` | `[]` | 不用的台站（台站表里的名字；两步都用） |

## clean（第 1 步）

| 键 | 默认 | 含义与依据 |
|---|---|---|
| `fs_hz` | init | 降采样后的采样率，须能由每个原始采样率整数倍抽取。建议值：满足 `cc.fmax_hz` ≤ 0.4 fs 的最低整数采样率（抗混叠低通在 0.8 倍奈奎斯特频率） |

## cc（第 2 步）

| 键 | 默认 | 含义与依据 |
|---|---|---|
| `fmin_hz` | init | 频带下限。建议值：最长台对 r_max 为最快速度下的 2 个波长，f = 2 v_fast / r_max（两台法远场，Bensen et al. 2007；Luo et al. 2015），再按高斯滤波器半增益宽度下调，× (1 − √(ln 2/α)) |
| `fmax_hz` | init | 频带上限。建议值：空间采样定理，台站间距不超过半个波长（如 Rost & Thomas 2002），Δx 取各台到最近台站距离的中位数：f = v_slow / (2Δx)，× (1 + √(ln 2/α))；时间采样：不超过最低原始采样率的 0.4 倍 |
| `maxlag_s` | init | 时延窗。建议值：r_max / (0.6 v_slow) + 4 / fmin。0.6 为约定：频散层中群速度可远低于相速度；4 个周期留作到时之后的噪声窗 |
| `segment_s` | init | 段长，须大于时延窗的两倍（循环相关）。建议值：max(2.5 maxlag, 10 / fmin) |
| `overlap` | 0.5 | 段的重叠（Welch 平均；Seats et al. 2012） |
| `max_over_std` | 10 | 最大振幅超过该台站日标准差这么多倍的段不用（NoisePy，Jiang & Denolle 2020）；0 为不用 |
| `substack_hours` | 24 | 子叠加长度，即刀切误差的组（至少需要 3 组）。记录不足 3 天时 `lant init` 建议更短的 |
| `max_distance_km` | `null` | 只计算这一距离以内的台对（超大台阵省内存和磁盘） |

## disp（第 3 步）

| 键 | 默认 | 含义与依据 |
|---|---|---|
| `fmin_hz`、`fmax_hz` | `null` | 测量频带；`null`：互相关频带向内收，使端点滤波器的半增益点在带内，f_c (1 ∓ √(ln 2/α)) |
| `per_octave` | 12 | 每倍频程的测量频率数（约定） |
| `reference` | `fj` | `fj`（本台阵 F-J 谱）、`model`（一维 Vs 模型）、`file`（一条曲线），见方法第 4 节 |
| `reference_file` | `null` | `model`：`depth_km Vs_km/s`；`file`：`frequency_Hz c_km/s` |
| `fj_band_hz` | `null` | F-J 拾取的频段（`null`：全部频率） |
| `alpha` | 5 | 高斯滤波器宽度，exp(−α((f′−f)/f)²)（FTAN；局部尺度的约定值） |
| `group_window` | 0.3 | 群到时在参考群速度 ±30 % 内搜索（约定） |
| `snr_min` | 5 | 到时的最低信噪比（文献常用 5–15） |
| `min_wavelengths` | 2 | 最短路径（波长数）（Bensen et al. 2007：3；Luo et al. 2015：1–2） |
| `max_jackknife_rel` | 0.02 | 刀切标准差与 c 之比的上限（约定；`null` 为不用） |
| `reference_error` | 0.1 | 参考曲线最大误差的假定值；只在路径不超过 1/(2 × 此值) 个波长处确定整周数（理论：分支间距的一半） |
| `slip_cycles` | 0.3 | 相位跳动超过这么多周时曲线结束（AFTAN） |
| `min_points` | 6 | 接受曲线的最少连续频点数（约定） |

## inv（第 4 步）

| 键 | 默认 | 含义与依据 |
|---|---|---|
| `data` | `phase` | 反演的走时：`phase`、`group` 或 `both`（`3_disp/dispersion.csv` 里的 Rayleigh 波相速度和/或群速度） |
| `grid_km` | `null` | 水平节点间距；`null`：实测波长中位数 / 3（横向分辨率约一个波长，用三个节点采样） |
| `depths_km` | `null` | 深度节点；`null`：从 0 起步长 max(最短波长/4, 深度/4)，到超过最长波长一半为止（Foti et al. 2018 第 4.2.3 节）；最后一个节点是固定半空间 |
| `start_model_file` | `null` | `depth_km Vs_km/s`；`null`：对所用数据的中位数曲线做一维反演 |
| `per_octave` | 4 | 每倍频程使用的周期数（测量周期的子集） |
| `iterations` | 6 | 三维反演迭代次数 |
| `weight` | 0.5 | 平滑权重（走时均方根拟合差（s）与均方根粗糙度（km/s）之比；数据项和模型项取均值，数值不随数据量和节点数变化；不是 DSurfTomo.in 里的 weight）。约定；加高斯噪声的合成反演中 0.35–0.5 恢复模型最好；越大越平滑 |
| `damping` | 0.1 | 向起始模型的阻尼（与 `weight` 一样按参数取均值），使没有射线的节点留在起始模型附近。约定 |
| `outlier_mad` | 3 | 偏离本周期残差中位数超过 k × 1.4826 MAD 的走时剔除，相速度和群速度分开（Leys et al. 2013：3 为"非常保守"） |
| `vs_min_kms`、`vs_max_kms` | `null` | Vs 上下限；`null`：起始模型最慢 Vs 的 0.6 倍和最快的 1.4 倍（约定） |

## checker

| 键 | 默认 | 含义与依据 |
|---|---|---|
| `amplitude` | 0.1 | 起始模型的 ±10 % |
| `noise` | `auto` | 按周期取三维模型剩余残差的稳健标准差（s）（Rawlinson & Spakman 2016：噪声与数据相当）；或给秒数 |

## plots、cores

| 键 | 默认 | 含义 |
|---|---|---|
| `plots.relief` | `auto` | 地图底下的地形：`auto` 用 GMT 的 `@earth_relief_03s`（需要联网；取不到时画不带地形的图）、地形网格文件，或 `null` |
| `plots.vs_range_kms` | `null` | Vs 色标范围（`null`：台站下方模型的 2–98 % 分位） |
| `plots.slice_depths_km` | `null` | 平面切片的深度（`null`：最多 8 个深度节点） |
| `plots.dpi` | 300 | PNG 图的分辨率 |
| `cores` | `auto` | `auto`：空闲 CPU 的 80 %（共享机器上给别人留余量）；或给数字 |
