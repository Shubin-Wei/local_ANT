# 输出文件说明

中文 | [English](outputs.md)

全部结果写在配置的 `project.output_dir` 下，原始数据不会被改动。下面按目录列出每个文件是什么、由哪一步写出、什么时候需要看它。标 ★ 的是每次运行都建议看的文件。

## 工程入口

先打开 `REPORT.md`：每一步的状态和下一条命令、关键数字、留给你判断的事项、主要图件，以及保留成果的链接和各自来自哪次运行。`README.md` 说明各目录和操作规则。`notes.md` 是你的笔记和决定（`lant note`）。

## results/：保留成果

`results/dispersion/` 在 disp 完整完成后替换；`results/model/`、`results/validation/`、`results/reproduce/` 在 evaluate 完整完成后一起替换，所以重建或重跑反演时，上一次评估过的模型一直保留到新模型评估完。这些目录里未登记的文件不会被替换。`lant clean-work` 之后它们仍可读；`work/` 不能整个删除。

| 文件 | 内容 |
|---|---|
| ★ `model/vs_model.csv` | 三维模型：`longitude_deg`、`latitude_deg`（WGS84）、`depth_km`、`vs_km_s`，按深度、纬度、经度排序；数值与 DSurfTomo 写出的完全相同。远离射线的节点保持起始模型：要结合棋盘图看 |
| ★ `dispersion/curves.csv` | 每个台对、每个测量频率：测量值、`accepted_qc`、`used_in_inversion` 和交给 DSurfTomo 的值 |
| `dispersion/disp.h5` | 全部测量，含拒绝原因和诊断量 |
| `dispersion/median_vs_reference.txt` | 实测中位曲线与参考曲线的比较 |
| `figures/` | 成果图（见下） |
| `validation/evaluation.json` | 起始模型和最终模型的残差统计、各深度的棋盘恢复、对照反演 |
| `validation/residuals_start.csv`、`residuals_final.csv` | 每个台对和周期的观测、预测走时，以及是否被 DSurfTomo 的规则保留 |
| `reproduce/` | `effective.yaml`（含全部取值的配置）、`stations.csv`、平均曲线和一维起始模型，以及求解器文件 `DSurfTomo.in`、`surfdata.dat`、`MOD`、`DSurfTomo.inMeasure.dat` |

逐台对的诊断图留在 `work/05_dispersion/pair_figures/`；棋盘、对照工程及其恢复模型留在 `work/06_inversion/`。

## .lant/：运行记录

| 文件 | 内容 |
|---|---|
| ★ `steps/<步骤>.md` | 每一步的报告：做了什么、关键数字、用到的每个参数（值、依据类别、来源）、要看的文件、留给你判断的事项、改什么能改变结果 |
| `summary/<步骤>.json` | 同一报告中全部数字，供程序读取 |
| `manifest.json` | 每一步依据哪些输入、参数和上游结果成立；参数或输入变化后，`lant status` 据此判断哪一步过时 |
| `auto_parameters.json` | 由数据和规则推出的每个参数值及其推导（代入的数字） |
| `config.yaml`、`config_history/` | 最近一次运行的完整配置（含全部默认值）和历次配置 |
| `run_log.jsonl` | 每一步的起止时间、耗时、成败 |
| `provenance/` | 代码版本（含未提交改动的指纹）、依赖版本、第三方程序的指纹 |
| `../notes.md` | 你的判断和理由，用 `lant note` 追加（AI 助手写的是你的话，不是它自己的决定） |
| `lock` | 运行锁；同一输出目录同时只能有一个 `lant` 进程写入 |
| `artifacts.json`、`published.json` | lant 写过的文件（供 `lant status`、`lant clean-work` 使用），以及每个保留成果来自哪次运行 |
| `backups/` | `--rebuild` 保留的步骤输出，以及被中断求解工程的输出 |
| `lant.log` | 各步骤的详细输出（每份 10 MiB，保留三份旧文件） |

## work/01_clean/：原始记录整理

| 文件 | 内容 |
|---|---|
| ★ `report.txt` | 文件数、格式、不能读的文件、名字对应、台站与处理日、推出的参数、资源估计、1 s 信号检验结果 |
| ★ `coverage.png` | 每个台站每个处理日的数据小时数 |
| ★ `stations.csv` | 实际使用的台站：`NET.STA`、经纬度、高程（此后所有步骤的台站名以它为准） |
| `processing_days.txt` | 每个处理日的标签和对应的 UTC 时间窗 |
| `station_days.csv` | 每个台站日：数据小时数、缺数、平移和插值的记录数、重叠、丢弃的短片段等计数 |
| `pulse_1s.csv` | 每个台站日 1 s 折叠检验的统计量、门限、有效自由度 |
| `input_files.csv`、`raw_index.json` | 读到的每个文件和记录 |
| `archive/` | 清洗后的日文件（每台每天一个 miniSEED，`*.done.json` 为续跑用的指纹）；`clean.write_archive: false` 时没有 |
| `archive_index.json` | 后面各步读取记录的索引 |
| `header_cache.json` | 文件头扫描缓存（文件未变时不再读） |

## work/02_ccf/：逐小时互相干

每个处理日一个 `day_YYYYDDD.h5`：`hourly [24, 台对, 时延]`、每小时有效段数 `count`、台对名、设置与输入指纹。这是数据量最大的目录；叠加完成后若不再需要按小时剔除，可以删除，但删除后改动 cc 之前的参数需要重算。

## work/03_stack/：叠加与台站诊断

| 文件 | 内容 |
|---|---|
| `stack.h5` | 叠加互相关 `stack`、子叠加 `daily`、各子叠加段数、台距、方位、台站坐标 |
| ★ `station_qc/station_qc.png` | 每个台站子叠加的相关系数（左）和时移（右），候选标红 |
| `station_qc/station_qc.txt` | 门限的数值和由来、全部候选及原因 |
| `candidates.txt` | 候选清单（确认剔除时使用的格式） |
| `removed.txt` | 实际剔除的台站日和台站小时 |

## work/04_reference/：参考曲线

| 文件 | 内容 |
|---|---|
| ★ `fj.png` | 频率—贝塞尔谱（逐频率归一化）、拾取点、参考曲线、测量频带、你给的速度范围 |
| `fj.txt` | 每个频率的台对数、拾取值、是否采用、参考曲线值 |
| `fj.npz` | F-J 图像数组 |
| `reference.npz` | 截到测量频带的参考曲线（disp 步骤读取） |
| `zones.npz`、`zones.png` | 分区参考（只在设置 `reference.zone_km` 时） |

## work/05_dispersion/：频散测量

| 文件 | 内容 |
|---|---|
| `disp.h5` | 每个台对、每个测量频率的相速度 `vc`、群速度 `vg`、信噪比、刀切误差、整周数、是否接受 `ok_c`、被拒原因位掩码 `reason`、频带边缘比值 `band_amp` |
| ★ `median_vs_reference.txt` | 各频率接受点的中位数与参考曲线之差；系统偏离说明参考曲线或整周数可能有问题 |
| `curves.csv` | 测量样点及质控、反演选择标记；独立反演列保留六位小数周期、五位小数速度 |
| `pair_figures/` | 单台对图（`export.curve_figures` 个按台距均匀抽取的台对，另加被拒点最多和偏离中位数最大的台对）：速度—频率图像、拾取、被拒原因 |

`reason` 的位：1 无群速度拾取，2 信噪比不足，4 无数据，8 群速度偏离参考，16 路径太短，32 旁瓣，64 无法定整周数，128 跳周，512 刀切误差，4096 群速度与相速度不一致，8192 连续段太短，16384 波长规则，32768 无参考曲线，65536 频带边缘。

## work/06_inversion/：反演工程

| 文件 | 内容 |
|---|---|
| `average_curve.txt` | 台阵平均曲线（每个周期的中位数和台对数） |
| `vs1d.txt`、`vs1d_fit.txt` | 一维起始模型及其对平均曲线的拟合 |
| `kernels.txt` | 各周期的深度敏感核 |
| `main/` | 主工程：`surfdata.dat`、`DSurfTomo.in`、`MOD`、`project.json`（网格、周期、稀疏预留）、`precheck.json`、`run.log`、最终模型 `DSurfTomo.inMeasure.dat`、各轮模型 `*.iterNNN` |
| `checkerboards/<名字>/` | 棋盘检验工程 |
| `checkerboard_noise.json` | 棋盘合成走时的噪声（`inversion.noise: auto`）：数值和来历 |
| `controls/<名字>/` | 对照工程（起始模型对照、正则化对照） |
| `evaluation/evaluation.json` | 残差、棋盘恢复（各深度的相关系数和斜率）、对照与主模型之差 |
| `evaluation/residuals_final.csv` | 每条走时：台站坐标、周期、距离、观测与预测走时、是否被剔除规则保留 |

`DSurfTomo.inMeasure.dat` 四列：经度、纬度、深度（km）、Vs（km/s）。

## results/figures/：成果图

| 图 | 内容 | 台阵 |
|---|---|---|
| `1_stations` | 台站叠在彩色晕渲地形上（SRTM，`plots.topography`），`data.exclude_stations` 去掉的台站画红叉 | 全部 |
| `2_ccf_moveout_<频带>` | 叠加互相关随台距排列，每个频带一张（测量频带在整数频率处切开）；台距画到该频带有测量的距离为止；两条虚线是该频带实测的最慢群速度和最快相速度 | 全部 |
| `2_ccf_stack_time_<k>` | 三个台对（近、中、远）：叠加互相关、按时间排列的逐小时互相关、对称分量信噪比随叠加小时数的变化 | 全部 |
| `3_dispersion_pair_<k>` | 同样三个台对：波形；相速度及其他整周数分支和参考曲线；群速度能量图和测得的到时 | 全部 |
| `3_dispersion_curves` | 全部接受的相速度曲线、平均值与标准差 | 全部 |
| `4_start_model`、`4_kernels` | 起始模型与拟合、各周期敏感深度 | 全部 |
| `5_depth_slices` | 各切片深度的 Vs，底图为灰色地形，标出剖面位置 | 面阵 |
| `5_slices_3d` | 五个深度切片的透视图，地形悬在上方 | 面阵 |
| `5_sections` | 剖面，上方为地形曲线，两端标字母和方位；按模型长深比作图（必要时注明垂向夸大） | 全部（线阵：沿测线） |
| `6_checkerboard_*` | 棋盘输入与恢复，红（慢）—白—蓝（快）：面阵为最多七个深度的平面图，线阵为沿测线的输入剖面和恢复剖面 | 全部 |
| `6_recovery` | 各深度的恢复斜率和相关系数 | 全部 |
| `6_residuals` | 起始模型和最终模型的走时残差：直方图、各周期均方根 | 全部 |

画不出的图在 figures 报告里列出原因（例如某频带没有测量结果、离线时取不到地形）。
默认只写 PNG；设置 `plots.formats: [png, pdf]` 可同时输出 PDF。
