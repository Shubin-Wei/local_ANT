# 参数说明

配置文件里每个键的含义、单位、默认值、依据、什么时候要改、设得不合适会怎样。默认值来自密集节点台阵（台间距几百米到十几千米、0.1–2 s）的实际研究。“依据”一栏写的是文献或开源程序；没有文献取值可依的，标为**项目取值**。

`array_type: line` 时会换一组默认值，见第 12 节。配置里没写的键都取默认值；写了不存在的键会报错，所以拼错的键不会被悄悄忽略。完整的默认值在 `lant/config.py`。默认为 `auto` 的键在 `lant run clean` 时按书面规则从数据推出（`lant/auto.py`），结果和依据写在 `<输出目录>/records/auto_parameters.json`；在配置里写了值就用写的值。

## 1. project

| 键 | 默认 | 说明 |
|---|---|---|
| `name` | my_array | 项目名，只用于图和记录 |
| `array_type` | auto | `surface` 面状台阵（台站铺成二维）；`line` 线状台阵（一条测线一个项目），决定一组默认值（第 12 节）。auto：台站分布主轴长宽比 > 4 为线状（面阵一般在 1–2，测线常在 10 以上） |
| `output_dir` | 必填 | 全部结果、图和记录都写在这里，不会写进代码目录 |

## 2. data：原始数据

| 键 | 默认 | 说明 |
|---|---|---|
| `raw_dir` | 必填 | 原始数据目录，递归搜索。目录结构和文件名都不要求，台站、通道、采样率、时间都从记录头读 |
| `format` | auto | ObsPy 能读的格式都可以：MSEED、SAC、GSE2、SEG2、SEISAN 等，逐个文件自动识别。识别很慢或误识别时可指定一种 |
| `file_patterns` | `["*"]` | 要尝试读取的文件名模式 |
| `ignore_patterns` | 常见非数据文件 | 不读的文件（`*.txt`、`*.csv`、`*.xml`、图片、压缩包等） |
| `station_file` | 必填 | 台站表，常见写法都能识别（tutorial 2.2）：CSV / 文本（有无表头、各种分隔符、中英文列名）、Excel、StationXML；台站名分列或连写；十进制度或度分秒。名字与记录头自动对应，结果写成 `01_clean/stations.csv` |
| `channel` | auto | auto：每台的垂直分量（通道以 Z 结尾，或 SAC cmpinc = 0；只有一个通道时用它）。通道名看不出垂直分量时给一个模式，如 `"*3"` |
| `network_default` | 无 | 可选：记录头没有台网代码时强制用的代码（一般不需要，名字会自动对应） |
| `rename` | `{}` | 可选：自动对应失败或有歧义时，指定 `{"记录头的 NET.STA": "台站表的 NET.STA"}` |
| `day_start_utc_hour` | auto | 处理日从这个 UTC 小时开始。只影响按“天”做的事：台站日筛选、刀切分组、出图的年积日。auto：60 % 以上的长记录（> 12 小时）从同一个整点开始时取该小时，否则 0。野外按北京时间切文件时推出 16 |
| `min_trace_length_s` | 0（线状 600） | 短于此长度（秒）的记录片段不用。**项目取值**：线状测线野外记录在布设、回收时有零碎片段 |
| `exclude_stations` | `[]` | 整台剔除的台站（判定仪器或耦合有问题），可写列表或一个文件名（每行一个台站）。这是对仪器的判断，不是规则生成的；依据通常是 `03_stack` 的筛选图、功率谱或与近邻台的相关 |
| `field_note_stations` | `[]` | 有野外备注（被挖出、移动过）的台站，台站日筛选对它们多用一条“持续偏低才截断”的规则（第 4 节） |
| `start_date`、`end_date` | 无 | 只处理这段日期（YYYY-MM-DD） |

## 3. clean：原始数据清洗

| 键 | 默认 | 说明 |
|---|---|---|
| `write_archive` | true | 写出清洗后的日文件（`01_clean/archive/NET.STA/…mseed`），后续步骤只读它。体积约与原始数据相当；磁盘不够时设 false，后续直接读原始文件（此时重叠段不会被剔除） |
| `encoding` | FLOAT32 | 日文件的 miniSEED 编码。整数数据用 STEIM2 可省一半以上空间、无损；FLOAT32 对 2^24 以内的整数也无损 |
| `nproc` | 4 | 并行进程数 |
| `min_free_disk_factor` | 1.2 | 可用磁盘小于估计体积的这么多倍时拒绝写出 |

清洗只做组织层面的事，**不改变样点值**：每台只用一种通道和采样率（数据最多的那种）；记录放到 UTC 采样网格上（偏差不超过 0.05 个采样时只平移起始时间，更大时 Lanczos 插值，两者都计数）；完全相同的重复和重叠去掉；**内容不同的重叠段两份都不用**（当作缺数），不替用户选；缺数保持缺数，不插补、不滤波、不去仪器响应。每一项处理的次数写在 `01_clean/station_days.csv`。

## 4. correlation：逐小时互相干（`02_ccf`）

| 键 | 默认 | 依据 / 说明 |
|---|---|---|
| `fs` | auto | 降采样后的采样率，原始采样率必须是它的整数倍。auto：50、40、25、20、10、5 Hz 中能整除所有原始采样率的最大者（200 Hz 节点仪 → 50 Hz）。要做 12 Hz 以上时手工提高 |
| `seg_len_s` / `seg_step_s` | 60 / 30 s | Xu et al. (2021, GJI 227:875) 第 3 节的 1 分钟窗、50 % 重叠。必须不小于 2 × maxlag |
| `maxlag_s` | auto | 时延范围。要大于最长台距 ÷ 最慢速度再留 2 个最长周期作噪声窗。auto：最长台距 / 1.2 km/s + 2 / dispersion.fmin_hz，向上取到 10 s（例：最长台距 14 km、频带下限 0.5 Hz → 20 s）。`seg_len_s` 自动不小于 2 × maxlag。太小：远台对的面波落在窗外或没有噪声窗 |
| `fmin_hz` / `fmax_hz` | 0.2 / auto | 互相干的频带（两端 0.1 Hz 余弦尖灭）。**项目取值**。auto 上限 = min(12.5, 0.4 fs) Hz |
| `smooth_bins` | 1 | 互相干分母的谱平滑点数，1 = 不平滑（逐频点用两台振幅谱归一，Xu et al. 2021 式 1） |
| `max_over_std` | 10 | 片段峰值超过整个台站日标准差的这么多倍就不用该片段（NoisePy `max_over_std`）。去掉地震、爆破、人为冲击。设太小会丢掉大量正常片段（`02_ccf` 日文件里 `segments_used` / `segments_gapfree` 可看比例，正常只丢百分之一二） |
| `remove_1s_comb` | true | 每台每小时把记录按 1 s 折叠求模板再减掉，去掉节点仪 GPS 授时脉冲造成的 1 Hz 梳状谱线（思路 Takagi et al. 2015；现象 Wang et al. 2017）。对没有这种干扰的仪器几乎没有影响；`04_notch/raw_lines.txt` 里若没有 1、2、3… Hz 的谱线可以关掉 |
| `nproc`、`threads_per_process` | 4、2 | 同时算几天、每个进程的线程数。内存：每个进程约 `台站数 × 17 MB + 台对数 × 0.2 MB`（160 台约 5.5 GB） |

叠加是按有效片段数加权的线性叠加。不用相位加权等改变相位的非线性叠加（它们会改变相位，使频散离散度变大）。

## 5. stack：叠加与台站日 / 台站小时筛选（`03_stack`）

| 键 | 默认 | 说明 |
|---|---|---|
| `block_hours` | 24（线状 12） | 子叠加的长度，刀切误差按子叠加分组。必须整除 24 |
| `station_day_qc.enabled` | true（线状 false） | 面阵按天筛选 |
| `station_day_qc.band_hz` | [1, 4] | 比较用的频带 |
| `station_day_qc.dist_km` | [2, 10] | 参与比较的台对台距 |
| `station_day_qc.corr_drop` | 0.25 | 某台某日与它自己总叠加的相关系数（对伙伴取中位）比当日全台阵中位数低这么多就剔除 |
| `station_day_qc.shift_s` | 0.04 s | 或整体时移不小于此值（钟差、仪器被移动） |
| `station_day_qc.min_partners` / `min_segments` | 10 / 200 | 当日伙伴少于 10 个、或台对当日少于 200 个片段（约 1.7 小时）就不评估（不评估 ≠ 剔除） |
| `station_day_qc.field_drop` | 0.05 | 野外备注台站：从某天起此后每个可评估日都比当日中位低这么多，才从那天截到记录结束 |
| `station_hour_qc.*` | 线状启用 | 线状只有几天记录，按小时筛：该小时互相关与总叠加的零延迟相关（1–4 Hz、±10 s、台距 ≤ 6 km）对伙伴取中位，低于 0.4（`threshold`）或伙伴少于 3 个（无法评估，例如测线收走后个别台还在记）都剔除 |
| `exclude_station_days` | `[]` | 另外手工剔除的台站日 `[台站, YYYYDDD]` |

以上规则都是**项目取值**；用互相关检测时间问题的思路见 Stehly et al. (2007)、Sens-Schönfelder (2008)。判断结果看 `03_stack/station_day_qc/station_day_qc.png`：大多数点应在 0 附近，被剔的点（红）应是孤立的离群点。若某天全台阵都偏低（噪声场变化），规则按“比当日中位低”来判，不会整天剔掉。剔除过多（超过 5 % 的台站日）要先查时间和仪器问题，而不是放宽门限。

## 6. notch：陷波（`04_notch`）

| 键 | 默认 | 说明 |
|---|---|---|
| `enabled` | true | 工业谱线持续存在时不陷波，谱线所在频段的相速度测点会明显减少 |
| `band_hz` | [0.5, 12] | 找谱线的频带 |
| `raw_nsta`、`raw_ndays`、`raw_days` | 20、4、[] | 原始波形谱线搜索随机抽的台数、均匀抽的天数（或直接给天） |
| `raw_min_frac` | 0.15 | 1 小时功率谱（0.01 Hz 分辨）高出 1 Hz 滑动中位 6 dB、且占 15 % 以上台站小时的窄峰算持续谱线 |
| `db` | 6 | 未陷波叠加的台对平均谱上高出 1 Hz 滑动中位 6 dB 的窄峰也陷波 |
| `table_file` | 无 | 用给定的陷波表（每行 `中心频率 半宽`，Hz）代替规则 |

半宽 = 谱线漂移范围的一半 + 0.03 Hz，不小于 0.05 Hz；相距 0.06 Hz 以内的合并。都是**项目取值**。陷波只做一次（成品已陷波，测频散时不再陷）。看 `04_notch/notch_report.png`：阴影应正好盖住谱上的尖峰；宽而缓的鼓包不是谱线，不应陷波。

## 7. reference：参考曲线（`05_reference`）

参考曲线用来引导群速度脊线、定相位的周数、做偏离筛选。

| 键 | 默认 | 说明 |
|---|---|---|
| `source` | fj | `fj`：本台阵的 F-J 谱（CC-FJpy `fj_noise`，频率-Bessel 变换，Wang, Wu & Chen 2019；Li et al. 2021）。`file`：借用别的台阵的曲线（`reference.file`：lant 写的 `reference.npz`，或两列 `频率_Hz 相速度_km/s` 的文本），插值到测量频点，范围外不外推、不测 |
| `fj_mode` | max_peak（线状 track_fundamental） | `max_peak`：每个频率取最大峰（假设基阶最强），去掉陷波频点后在 log f 上三次拟合到 `fj_fit_fmax`，拟合范围外取端点值。`track_fundamental`：从 `track_start_band_hz` 里一个无歧义的峰出发逐频点追基阶脊线（每 1/12 倍频程变化不超过 5 %、随频率上升不超过 3 %），三次拟合并强制随频率不增 |
| `fj_rmax_lambda` | 5 | 每个频率只用台距不超过 5 个（按 `fj_c_nominal` = 2.5 km/s 算的）波长的台对：长路径在高频彼此不相干，F-J 图会碎成斑点 |
| `fj_cmin`/`fj_cmax`/`fj_dc` | 1.2/4.0/0.005 | F-J 扫描速度范围和步长（km/s） |
| `fj_fit_fmax` | 6 Hz | 高频 F-J 图常没有稳定的脊线，高于此频率曲线取 6 Hz 的值；看 `05_reference/fj.png` 判断 |
| `pass2` | true | 第二遍参考：用第一遍参考测量（不做任何偏离筛选），取每个频率被接受的相速度的中位数作第二遍参考，再正式测量。中位数是为了让 ±20 % 偏离筛选居中，不代表它比 F-J 更接近真值；第二遍参考对第一遍参考的选择不敏感 |
| `pass2_min_pairs` | 50（线状 20） | 某频率少于这么多台对就不取中位数；缺口不超过 `pass2_max_gap_octave`（0.25 倍频程）时在 log f 上线性插补，否则那里没有参考、不测量 |

**线状台阵**：测线的 F-J 频带窄，噪声常斜入射，本线的 F-J 曲线可能系统偏快，本线的中位曲线可能有台阶。同一地区有面阵（或更好的测线）时，借用它的第二遍参考（`source: file`，`file: <面阵输出>/05_reference/reference.npz`）通常更好；没有时用本线的 `track_fundamental`，并仔细看 `05_reference/fj.png`。借用时线状与面阵结果相符就不是独立的验证。

## 8. dispersion：频散测量与质控（`06_dispersion`）

测量程序是两台法：正负时延平均的对称分量（Bensen et al. 2007），高斯窄带滤波，群速度动态规划追脊，相速度在群到时处读解析信号相位，相位常数 π/4（Lin et al. 2008 式 7），在 1–3.5 个参考波长处用参考曲线定周数后按实测群慢度向高频递推（AFTAN `phtovel.f`），按子叠加分组刀切估计误差。

| 键 | 默认 | 说明 |
|---|---|---|
| `fmin_hz`/`fmax_hz`/`per_octave` | auto/auto/12 | 测量频点：对数等间隔，每倍频程 12 点（例如 0.5–10.08 Hz 共 53 点）。auto 下限：90 % 的台距不小于 2 个波长（按 2.5 km/s）的频率，向下取到过 1 Hz 的 1/12 倍频程网格（例：9.5 km → 0.5 Hz）；auto 上限 = 0.84 × 互相干上限（10.5 Hz）。DSurfTomo 最多 80 个周期 |
| `vmin`/`vmax`/`dv` | 1.2/4.0/0.01 km/s | 速度搜索范围。要包住本区可能的速度；落在边界上的群速度拾取一般通不过相速度质控（线状噪声斜入射时沿线视速度偏高，落在上边界的更多） |
| `alpha` | 5 | 高斯滤波相对带宽 H = exp(−α((f−fc)/fc)²)，与 FTAN 同（Bensen et al. 2007 第 4 节），经验参数。越大越窄：频率分辨好、时间分辨差 |
| `phase_constant` | pi/4 | 远场 2-D 格林函数相位（各向同性噪声场） |
| `min_dist_km` | 0.3 | 更短的台对不测 |
| `lam_min` | 1.0 | 台距不足 1 个参考波长的频点不测（测量阶段的宽门限；最终门限见 qc.lam_c） |
| `anchor_max_lambda` | 3.5 | 周数只在台距为 1.0–3.5 个参考波长处定。**项目取值**：台距为 N 个波长时相邻周数的速度相差约 1/N，3.5 λ 时约 29 %，远大于横向不均匀。代价是一部分合格频点因所在连续段没有锚点而没有相速度 |
| `njack` | 5（线状 auto） | 刀切组数：子叠加轮流分到各组。auto = 数据量达到最满子叠加 25 % 的子叠加个数 |
| `min_jackknife_groups` | 3 | 有效组少于此值时不给误差（该频点不能通过刀切质控） |
| `nproc`、`chunk` | 8、20 | 并行进程数、每块台对数 |
| `side_off_bands` | "" | 在这些频带不用旁瓣判据（`"lo-hi,lo-hi"`）；默认处处都用 |
| `qc.snr_c` | 6 | 读相位处的信噪比门限（信号窗峰值 ÷ 尾部噪声均方根）。文献 4–10 |
| `qc.lam_c` | 2 | 最终相速度要求台距 ≥ 2 个波长（用该点自己的相速度算）。Bensen et al. (2007) 用 3、Luo et al. (2015) 用 1、Ma et al. (2024)、Jing et al. (2025) 用 2。放到 1 个点数明显增多，但近场偏差增大 |
| `qc.jk_c_rel` | 0.02 | 刀切标准差 ≤ 2 % |
| `qc.dev_c` | 0.20（线状 0） | 相速度偏离参考 ≤ 20 %。0 = 不筛（线状借用面阵参考时实测可能整体偏离参考，±20 % 会单侧截断） |
| `qc.dev_pre` | 0.30 | 读相位前群速度偏离参考 ≤ 30 % |
| `qc.smooth_c` | 0.03 | 相邻频点 \|Δln c\| ≤ 3 % |
| `qc.ucons_c` | 0.25 | 由相速度推出的群速度与实测群速度相差 ≤ 25 % |
| `qc.side_c` | 1.0 | 群速度图上另有一个独立的峰（两峰间的谷低于较矮者的 70 %）不低于所选的峰时，该频点不测相速度。**项目取值**，是挡掉点数最多的一条：关掉它点数显著增加（高频尤甚），但曲线明显变乱，被它挡掉的点走时残差离散度明显大于保留的点 |
| `qc.min_run` | 6 | 通过的频点连成不少于 6 个（半个倍频程）才保留 |
| `group_qc.*` | 见说明 | 群速度只作产品，不进反演：snr_g 8、side_g 0.9、jk_g 3 %、dev_g 25 %、lam_g 3 个**相速度**波长（Bensen 2007 第 5 节）；面阵另要求正负两支群速度相差 ≤ 20 %（`branch_diff`）且弱支信噪比 ≥ 4（`snr_branch_min`），线状改用相邻频点跳变 ≤ 3 %（`step_g`），因为线状噪声方向性强，两支本来不对称 |

判断：看 `06_dispersion/pair_figures/`（抽样台对的波形、相速度图、群速度图和拾取点：拾取点应在能量脊上）、`figures/3_dispersion_curves.png`（曲线应集中，离群曲线少）、`06_dispersion/precision.txt`（端点各在 0.5 km 内的冗余路径的相对差，代表测量精度）、`records/summary/disp2.json` 的每个频点点数。

## 9. export：导出（`06_dispersion/curves`）

| 键 | 默认 | 说明 |
|---|---|---|
| `tmin_s`/`tmax_s` | 0.09/auto | 导出的周期范围（按测量频点的周期 1/f，不插值、不外推）。auto 上限 = 1 / dispersion.fmin_hz |
| `one_segment` | true | 每条相速度曲线只留最长的一段连续曲线（一样长时取平均信噪比高的），依据 EGFAnalysisTimeFreq 手册（Yao 2015）p. 9 “select only the most continuous part with higher SNR”。一般只去掉百分之一二的点 |
| `curve_figures` | 200 | 画单台对频散图的台对数（0 不画，−1 全画；全画很慢） |

文件格式：头两行是两台的 `经度 纬度`，此后每行 `周期_s 速度_km/s 信噪比 误差_km/s`，速度 0 表示该周期没有测量。误差 = max(刀切标准差, 该频带冗余路径相对差 × 速度)。DSurfTomo 不读第 4 列。

## 10. inversion：DSurfTomo 反演（`07_inversion`）

方法是 DSurfTomo（Fang et al. 2015）直接反演三维 Vs，只用 Rayleigh 波相速度。

| 键 | 默认 | 说明 |
|---|---|---|
| `median_min_pairs` | 50（线状 20） | 平均曲线：每个周期取台对中位数，少于此数的周期不要。平均曲线决定反演的最长周期、深度节点和起始模型 |
| `depth_sequence_km` | 0, 0.2, … 3.0（线状多 0.1） | 标准深度节点。不深于平均曲线最长波长一半的节点参与反演，下一个节点固定（DSurfTomo 不更新、不输出最后一个节点，它的速度是边界假设，不解释）。“最长波长的一半”是 Foti et al. (2018) 4.2.3 节的经验探测深度。序列不够深时自动按步长逐次 ×1.3 加深（记录在 inversion 的摘要里） |
| `fine_step_km`/`fine_zmax_km` | 0.1/1.5（线状 0） | 在 0–1.5 km 之间加到每 0.1 km 一个节点（起始模型在这些节点上线性插值）。加密不会改善垂向分辨（分辨由数据决定，看垂向棋盘），只让切片可取更多深度 |
| `grid_spacing_deg` | auto | 水平网格间隔（纬度、经度方向），先按它定范围（台站范围外扩 max(1 %, 3 格)），再按 `grid_refine` 细分。可分辨尺度看棋盘，不由网格决定。auto：面阵细分后的间隔约为最近台站间距中位数的一半，取 0.0025° 的倍数（例：0.45 km → 0.0025°，粗网格 0.005° × `grid_refine` 2）；线状 0.005° |
| `line_strip_cells` | 无（线状 1） | 线状：横向只留台站两侧各约 1 格（DSurfTomo 从第 2 个节点起算走时，宽 ny−3 格；台站贴边不报错但走时算成 0，所以不能再窄） |
| `start_model.source` | invert | `invert`：由平均曲线反演一维起始模型（阻尼平滑高斯-牛顿，正演用 DSurfTomo 自己的正演链；初值 1.087 × 最短周期相速度、平滑 2.0、阻尼 0.05、14 次、每步限幅 ±15 %，都是**项目取值**）。`file`：给一个两列 `深度_km Vs_km/s` 的模型（例如线状借用面阵的 `07_inversion/vs1d.txt`）。不用 DSurfTomo 自带的 GenerateIniMOD.py（它令 Vs 等于相速度，偏慢约 9 %） |
| `weight`/`damp` | 1.0/1.0（线状 1.5/1.5） | 平滑和阻尼。DSurfTomo 自带示例 4.0/1.0；常规范围 weight < 10、多数 < 5，damp 0.5 或 1。做法：在 0.5/0.5、1.0/1.0、1.5/1.5 几组里并排看切片、棋盘、残差后选定。**不要用 L 曲线选**：程序每轮惩罚的是模型更新量，权衡曲线不对应固定的目标函数，L 曲线常给出明显过度平滑的值。程序的平滑和阻尼不随走时条数归一，数据量差很多时同一个值的效果不同 |
| `sublayers` | 3 | 算深度核的子层数（2–5） |
| `vmin`/`vmax` | 1.5/4.5 km/s | 速度上下限；不应有节点碰到 |
| `iterations` | 10 | 手册：10 次通常足够。程序没有收敛判据，看 `evaluate` 的逐轮残差 |
| `sparsity` | 0.05（线状 0.30） | 只影响内存分配：每条射线占全部节点的比例上限。窄网格（线状）要大 |
| `threshold` | 3 | 残差剔除：每轮把带符号残差在 [3×第 25 百分位, 3×第 75 百分位] 之外的走时权重置零。它默认残差以零为中心，**起始模型整体偏离数据时第一轮会单侧剔除大量好数据**，所以起始模型必须与数据整体一致 |
| `noise` | 0.02 | 棋盘合成数据加的噪声（2 %，自带示例值；小于实际残差时棋盘偏乐观） |
| `solver` | dsurftomo | `pdsurftomo`：并行版（Zhu et al. 2026）Mode 2，与 DSurfTomo 给出相同的模型，快但内存大（十几万条走时需十几 GB；安装脚本打了降低内存预留的补丁，不改变结果） |
| `threads`/`njobs` | 4/1 | 每个工程的线程数、同时运行的工程数（主工程 + 棋盘 + 对照）。内存 ≈ njobs × 单个工程峰值 |

## 11. checks 与 plots

| 键 | 默认 | 说明 |
|---|---|---|
| `checks.checkerboards` | auto | 水平棋盘：`[{name, cell_km: [南北, 东西]}]`，格子边长（半波长），0 = 该方向不变；线状为 [沿线, 横向]。幅度 `amplitude` ±0.4 km/s（自带工具默认）。auto：孔径的 0.12 / 0.18 / 0.25 倍（例：孔径 14 km → 1.7 / 2.5 / 3.5 km）；线状为测线长度的 0.165 / 0.24 / 0.34 倍，只沿测线变化 |
| `checks.vertical` | auto | 随深度变号的棋盘：常深度的棋盘看不出垂向拖尾（Lévêque et al. 1993；Rawlinson & Spakman 2016）。auto：格子 0.18 × 孔径（线状 0.34 × 长度），在最深可反演节点的 0.21 / 0.46 / 0.71 倍处变号（例：2.4 km → 0.5 / 1.1 / 1.7 km）；线状另加一个只在 0.46 倍处变号的 |
| `checks.controls` | 正则化、起始模型各两组 | `reg_w<weight>_d<damp>`：只改正则化；`start_x<系数>`：起始模型整体乘系数。只报差值，不用来挑参数 |
| `checks.min_rays_per_node` | 5 | 统计和周期图的着色规则 |
| `checks.iteration_residuals` | true | evaluate 时对每一轮的模型也做一次正演，得到残差随迭代的变化；每轮要一次正演，耗时约等于一轮反演；数据量大时可关掉，只算起始和最终模型 |
| `plots.topography` | 无 | GMT 网格文件，或 GMT 远程网格（如 `@earth_relief_01s`，首次使用会下载，GMT 6.1 以上）；不给就不画地形 |
| `plots.vs_range`/`vs_cpt` | 1.9–3.6/seis | 色标范围和色表（红慢蓝快）。多张图用同一个范围才可比 |
| `plots.slice_depths_km` | 0.2…1.5 | 深度切片（必须是模型节点） |
| `plots.sections` | auto | 面阵剖面：`[{name: A, start: [经度, 纬度], end: [经度, 纬度]}]`；auto：过台站中心、沿两个主轴方向各一条；线状自动沿测线 |
| `plots.section_zmax_km` | 最深可更新节点 | 剖面最大深度 |
| `plots.contour_step` | 0.2 km/s | 剖面等值线间隔 |
| `plots.formats`/`dpi`/`font_scale` | png+pdf/300/1 | 导出格式、分辨率、字号缩放（幻灯片可设 1.5–2） |
| `plots.pair_example` | 自动 | 频散示例图用的台对（`NET.STA_NET.STA`） |

## 12. array_type: line 时改变的默认值

| 键 | 面阵 | 线状 | 原因 |
|---|---|---|---|
| `data.min_trace_length_s` | 0 | 600 | 布设、回收时的零碎记录 |
| 筛选 | 台站日 | 台站小时（门限 0.4） | 每条线只有几天记录 |
| `stack.block_hours` | 24 | 12 | 刀切需要足够多的子叠加 |
| `dispersion.njack` | 5 | auto | 同上 |
| `reference.fj_mode` | max_peak | track_fundamental | 测线上最大峰不一定是基阶（高阶模式、斜入射） |
| `reference.pass2_min_pairs`、`inversion.median_min_pairs` | 50 | 20 | 台对少 |
| `dispersion.qc.dev_c` | 0.20 | 0（不筛） | 借用参考时实测整体偏离参考 |
| 群速度 | 正负两支一致 | 相邻频点跳变 ≤ 3 % | 线状噪声方向性强 |
| `inversion.depth_sequence_km` | 0, 0.2, … | 多一个 0.1 km | 测线台站密、高频数据多，浅部加一个节点 |
| `inversion.fine_step_km`、`grid_refine` | 0.1、2 | 0、1 | |
| `inversion.line_strip_cells` | — | 1 | 窄条网格（Xie et al. 2024 规则网格用于线状台阵） |
| `inversion.weight/damp`、`sparsity` | 1.0/1.0、0.05 | 1.5/1.5、0.30 | 线状数据少、约束弱，正则化稍强；窄网格每条射线占的节点比例大 |
| 棋盘 | 正方形 | 只沿测线变化 | |

## 参考文献

- Bensen et al. (2007) GJI 169:1239 —— 预处理、FTAN、信噪比、3 个波长、时间稳定性
- Brocher (2005) BSSA 95:2081 —— Vp、密度与 Vs 的经验关系（DSurfTomo 内部）
- Fang et al. (2015) GJI 201:1251 —— DSurfTomo；手册 v1.4
- Foti et al. (2018) Bull. Earthquake Eng. 16:2367 —— 探测深度约为最长波长的一半
- Jiang & Denolle (2020) SRL 91:1853 —— NoisePy（片段剔除 max_over_std）
- Jing et al. (2025) Ore Geol. Rev. 184:106718 —— 平均频散作起始模型、2 个波长、棋盘
- Lévêque et al. (1993) GJI 115:313；Rawlinson & Spakman (2016) GJI 205:1221 —— 棋盘检验及其局限
- Li et al. (2021) SRL 92:3179 —— CC-FJpy
- Lin et al. (2008) GJI 173:281 —— 相位常数 π/4
- Luo et al. (2015) GJI 201:652；Ma et al. (2024) Tectonophysics —— 波长门限
- Sens-Schönfelder (2008) GJI；Stehly et al. (2007) GJI —— 用互相关检测时间误差
- Takagi et al. (2015) SRL；Wang et al. (2017) JGR —— 周期性仪器噪声
- Wang, Wu & Chen (2019) JGR 124:3708 —— F-J 方法
- Xie et al. (2024) GJI 239:1775 —— DSurfTomo 用于密集线状台阵
- Xu et al. (2021) GJI 227:875 —— 互相干、1 分钟窗、对数等间隔频点
- Yao (2015) EGFAnalysisTimeFreq 手册 —— 每条曲线只留最连续、信噪比高的一段
- Zhu et al. (2026) arXiv:2604.11920 —— pDSurfTomo
