# 输出文件说明

中文 | [English](outputs.md)

全部结果写在配置的 `project.output_dir` 下，原始数据不会被改动。下面按目录列出每个文件是什么、由哪一步写出、什么时候需要看它。标 ★ 的是每次运行都建议看的文件。

## 总览

先打开 `REPORT.md`：每一步的状态、耗时和下一条命令，关键数字（含平滑权重），留给你判断的事项、主要图件，以及保留成果的链接和各自来自哪次运行。`README.md` 说明各目录和操作规则。`notes.md` 是你的笔记和决定（`lant note`）。

| 文件 | 内容 |
|---|---|
| `REPORT.md`、`README.md` | 报告和目录说明（lant 写） |
| `AGENTS.md`、`CLAUDE.md` 或 `GEMINI.md` | 给 AI 助手的规则：按 `project.ai_assistant` 所选助手只写一份完整规则（codex 写 AGENTS.md，claude 写 CLAUDE.md，gemini 写 GEMINI.md；未选写 AGENTS.md）；已有文件不覆盖 |
| `notes.md` | 你的笔记和决定（`lant note`） |
| `scratch/` | 可选：人和助手的临时文件；`lant status` 不检查其中文件（只报告大小），可以随时删除 |

这里以及 `work/`、`results/`、`.lant/` 中任何不是 lant 写的文件，`lant status` 都会列为未登记；lant 不会删除或覆盖它们。

## results/：保留成果

`results/dispersion/` 在 disp 完整完成后替换；`results/model/`、`results/validation/` 在 evaluate 完整完成后一起替换，所以重建或重跑反演时，上一次评估过的模型一直保留到新模型评估完。这些目录里未登记的文件不会被替换。`lant clean-work` 之后它们仍可读；`work/` 不能整个删除。

每步完整完成后重写 `results/parameters.yaml`：到目前为止全部取值已确定的配置，以及你的确认记录（每一步实际用的取值在该步报告 `.lant/steps/<step>.md` 和 `.lant/manifest.json` 里）。清洗覆盖、台站诊断和 F-J 图在相应步骤完成时移到 `results/figures/`（只在那里保留一份）；`REPORT.md` 链接它们并注明来源轮次。重跑时关闭某项诊断会保留此前已发布图，仍标为此前的轮次。

| 文件 | 内容 |
|---|---|
| ★ `model/vs_model.csv` | 三维模型：`longitude_deg`、`latitude_deg`（WGS84）、`depth_km`、`vs_km_s`、`dws`（最后一轮该节点的数据敏感度之和：没有射线经过为 0，不参与更新的节点为空），按深度、纬度、经度排序；数值与求解器写出的完全相同。远离射线的节点保持起始模型：要结合棋盘图看 |
| ★ `dispersion/curves.csv` | 每个台对、每个测量频率：测量值、`accepted_qc`、`used_in_inversion` 和交给 DSurfTomo 的值 |
| `parameters.yaml` | 全部取值已确定的配置及确认记录；由最近完成的步骤写出（部分 auto 值在后续步骤才推出） |
| `figures/` | 成果图（见下） |
| `validation/evaluation.json` | 起始模型和最终模型的残差统计、各深度的棋盘恢复、对照反演 |
| `validation/checkerboard_recovery.csv` | 每个棋盘、每个深度：异常的恢复比例（斜率）、图案相关系数 r、节点数 |
| `validation/residuals_start.csv`、`residuals_final.csv` | 每个台对和周期的观测、预测走时，以及是否被 DSurfTomo 的规则保留 |

重跑或重画所需的一切都在 `work/`（不必重新计算）：主模型的求解器文件在 `work/06_inversion/main/`（`DSurfTomo.in`、`surfdata.dat`、`MOD`、最终模型 `DSurfTomo.inMeasure.dat`），起始模型和平均曲线在 `work/06_inversion/`，台站表在 `work/01_clean/stations.csv`。全部测量及其拒绝原因和诊断量（`disp.h5`、`median_vs_reference.txt`）、逐台对的诊断图留在 `work/05_dispersion/`（`results/` 里不再复制）；棋盘、对照工程及其恢复模型留在 `work/06_inversion/`。棋盘图只画当前配置的检验工程（`work/06_inversion/checks.json`，与评估所用的相同）；以前配置留下的工程目录保留在 `work/`，但不再出图。

## .lant/：运行记录

| 文件 | 内容 |
|---|---|
| ★ `steps/<步骤>.md` | 每一步的报告：做了什么、关键数字、用到的每个参数（值、依据类别、来源）、要看的文件、留给你判断的事项、改什么能改变结果 |
| `summary/<步骤>.json` | 同一报告中全部数字，供程序读取 |
| `manifest.json` | 每一步依据哪些输入、参数和上游结果成立；参数或输入变化后，`lant status` 据此判断哪一步过时 |
| `auto_parameters.json` | 由数据和规则推出的每个参数值及其推导（代入的数字） |
| `config.yaml`、`config_history/` | 最近一次运行的完整配置（含全部默认值）和历次配置 |
| `run_log.jsonl` | 每一步的起止时间、耗时、成败 |
| `provenance/` | 代码版本（含未提交改动的指纹）、依赖版本、第三方程序的指纹；进程导入时固定的包内源码 SHA-256 |
| `forward_cache/*.npz` | 经完整依赖和校验和核查的正演预测，仅供预检、自动噪声和残差复用；内部可重建数据，不是结果模型 |
| `../notes.md` | 你的判断和理由，用 `lant note` 追加（AI 助手写的是你的话，不是它自己的决定） |
| `lock` | 运行锁；同一输出目录同时只能有一个 `lant` 进程写入 |
| `artifacts.json`、`published.json` | lant 写过的文件（供 `lant status`、`lant clean-work` 使用），以及每个保留成果来自哪次运行 |
| `backups/` | `--rebuild` 保留的步骤输出，以及被中断求解工程的输出 |
| `lant.log` | 各步骤的详细输出（每份 10 MiB，保留三份旧文件） |

`run_log.jsonl` 中的 `phase_start` / `phase` 记录预检、求解器批次/工程、自动噪声和评估的墙钟时间与结束状态；`forward_evaluation` 另记录命中/未命中、身份检查、真正计算及总时间。并行或嵌套阶段不能直接相加。

## work/01_clean/：原始记录整理

| 文件 | 内容 |
|---|---|
| ★ `report.txt` | 文件数、格式、不能读的文件、名字对应、台站与处理日、推出的参数、资源估计、1 s 信号检验结果 |
| ★ `stations.csv` | 实际使用的台站：`NET.STA`、经纬度、高程（此后所有步骤的台站名以它为准） |
| `processing_days.txt` | 每个处理日的标签和对应的 UTC 时间窗 |
| `station_days.csv` | 每个台站日：数据小时数、缺数、平移和插值的记录数、重叠、丢弃的短片段等计数 |
| `pulse_1s.csv` | 每个台站日 1 s 折叠检验的统计量、门限、有效自由度 |
| `input_files.csv`、`raw_index.json` | 读到的每个文件和记录；首尾两条记录都属于不可能在台站表里的台站的 miniSEED 文件，只读这两条（`not_read_because` 列写明） |
| `archive/` | 清洗后的日文件（每台每天一个 miniSEED，`*.done.json` 为续跑用的指纹）；`clean.write_archive: false` 时没有 |
| `archive_index.json` | 后面各步读取记录的索引 |
| `header_cache.json` | 文件头扫描缓存（文件未变时不再读） |

## work/02_ccf/：逐小时互相干

每个处理日一个 `day_YYYYDDD.h5`：`hourly [24, 台对, 时延]`、每小时有效段数 `count`、台对名、设置与输入指纹。这是数据量最大的目录；叠加完成后若不再需要按小时剔除，可以删除，但删除后改动 cc 之前的参数需要重算。

## work/03_stack/：叠加与台站诊断

| 文件 | 内容 |
|---|---|
| `stack.h5` | 叠加互相关 `stack`、子叠加 `daily`、各子叠加段数、台距、方位、台站坐标 |
| `station_qc/station_qc.txt` | 门限的数值和由来、全部候选及原因 |
| `candidates.txt` | 候选清单（确认剔除时使用的格式） |
| `removed.txt` | 实际剔除的台站日和台站小时 |
| `all_hours.json` | `stack.h5` 为全部小时的叠加时写下：其输入的指纹。`lant confirm stack --exclude-candidates` 之后的那次运行据此直接从中剔除候选，不再把每个小时重新叠加 |

## work/04_reference/：参考曲线

| 文件 | 内容 |
|---|---|
| `fj.txt` | 每个频率的台对数、拾取值、是否采用、参考曲线值 |
| `fj.npz` | F-J 图像数组 |
| `reference.npz` | 截到测量频带的参考曲线（disp 步骤读取） |
| `zones.npz` | 分区参考（只在设置 `reference.zone_km` 时） |

clean（数据覆盖）、stack（台站诊断）和 reference（F-J 谱、分区）的图由步骤写在这里，步骤完成时移到 `results/figures/`（`0_data_coverage`、`2_station_qc`、`3_fj_phase_velocity`、`3_reference_zones`）；步骤失败时留在这里。

## work/05_dispersion/：频散测量

| 文件 | 内容 |
|---|---|
| `disp.h5` | 每个台对、每个测量频率的相速度 `vc`、群速度 `vg`、信噪比、刀切误差、整周数、是否接受 `ok_c`、被拒原因位掩码 `reason` |
| ★ `median_vs_reference.txt` | 各频率接受点的中位数与参考曲线之差；系统偏离说明参考曲线或整周数可能有问题 |
| `curves.csv` | 测量样点及质控、反演选择标记；独立反演列保留六位小数周期、五位小数速度 |
| `pair_figures/` | 单台对图（`export.curve_figures` 个按台距均匀抽取的台对，另加被拒点最多和偏离中位数最大的台对）：速度—频率图像、拾取、被拒原因 |

`reason` 的位：1 无群速度拾取，2 信噪比不足，4 无数据，8 群速度偏离参考，16 路径太短，32 旁瓣，64 无法定整周数，128 跳周，512 刀切误差，4096 群速度与相速度不一致，8192 连续段太短，16384 波长规则，32768 无参考曲线（65536 频带边缘已在 2.0.1 去掉）。

## work/06_inversion/：反演工程

| 文件 | 内容 |
|---|---|
| `average_curve.txt` | 台阵平均曲线（每个周期的中位数和台对数） |
| `vs1d.txt`、`vs1d_fit.txt` | 一维起始模型及其对平均曲线的拟合 |
| `kernels.txt` | 各周期的深度敏感核 |
| `main/` | 主工程：`surfdata.dat`、`DSurfTomo.in`、`MOD`、`project.json`（网格及其间距的由来、周期及每倍频程个数、平滑权重）、`precheck.json`、`run.log`、`solver.json`（实现、线程数、各轮各模块秒数、LSMR 停止代码和被剔除的走时数、峰值内存）、`DWS.dat`（最后一轮每个节点的 DWS）、最终模型 `DSurfTomo.inMeasure.dat`、各轮模型 `*.iterNNN` |
| `checkerboards/<名字>/` | 棋盘检验工程：输入、`MOD.true`、真模型和恢复模型（`Vs_model.real`、`DSurfTomo.inSyn.dat`）、日志。评估后删除这里和 `controls/` 中 lant 不读的 DSurfTomo 文件（`*.iterNNN`、`residualFirst.dat`、`residualLast.dat`、`velmap2dRc.dat`），仅当它们由 lant 写出且未改动；`main/` 保留 |
| `checkerboard_noise.json` | 棋盘合成走时的噪声（`inversion.noise: auto`）：每个周期的标准差（秒）和来历；每个棋盘工程里存为 `noise_seconds.txt` |
| `controls/<名字>/` | 对照工程（起始模型对照、正则化对照） |
| `evaluation/evaluation.json` | 残差、棋盘恢复（各深度的相关系数和斜率）、对照与主模型之差 |
| `evaluation/residuals_final.csv` | 每条走时：台站坐标、周期、距离、观测与预测走时、是否被剔除规则保留 |

`DSurfTomo.inMeasure.dat` 四列：经度、纬度、深度（km）、Vs（km/s）。

## results/figures/：成果图

| 图 | 内容 | 台阵 |
|---|---|---|
| `0_data_coverage` | 清洗后的数据覆盖；clean 完成即发布 | 全部 |
| `2_station_qc` | 台站诊断，候选需用户判断；启用时随 stack 发布 | 全部 |
| `3_fj_phase_velocity`、`3_reference_zones` | 阵列 F-J 谱与可选分区参考；实际计算时随 reference 发布，区别于单台对能量图 | 全部 |
| `1_stations` | 台站叠在彩色轻度晕渲地形上（SRTM，`plots.topography`），限制图高并在图内图例区分有记录台站数与反演台站数；`data.exclude_stations` 去掉的台站画红叉 | 全部 |
| `2_ccf_moveout_<频带>` | 叠加互相关随台距排列，每个频带一张（测量频带在整数频率处切开）；台距画到该频带有测量的距离为止；两条虚线是该频带实测的最慢群速度和最快相速度 | 全部 |
| `2_ccf_stack_time_<k>` | 三个台对（近、中、远）：叠加互相关、按时间排列的逐小时互相关、对称分量信噪比随叠加小时数的变化 | 全部 |
| `3_dispersion_pair_<k>` | 同样三个台对：波形；相位匹配滤波拾取及整周数分支与 F-J 参考曲线（不画传统窄带相位背景）；群速度包络能量图和测得的到时 | 全部 |
| `3_dispersion_curves` | 全部接受的相速度曲线、平均值与标准差 | 全部 |
| `4_start_model` | 中位观测相速度与一维初模拟合预测，含 RMS；连接 Vs 参数节点并标固定底值，不表示探测到层界面 | 全部 |
| `4_kernels` | 独立加密延伸模型上的一维联合敏感性示意（Vs ±5%，Vp/密度随动），单位每 km；不是实际反演 Jacobian，虚线为反演底边界 | 全部 |
| `5_depth_slices` | 各切片深度的 Vs，底图为灰色地形，标出剖面位置 | 面阵 |
| `5_slices_3d` | 五个深度切片的透视图，地形悬在上方 | 面阵 |
| `5_sections` | 剖面，上方为参考地形曲线，稀疏等值线与两端字母/方位；线阵并列相对一维参考模型的百分比异常面板，底部注明垂向放大倍数与未作分辨率遮罩 | 全部（线阵：沿测线） |
| `6_checkerboard_*` | 棋盘输入与恢复，红（慢）—白（零）—蓝（快）：面阵为最多七个深度的平面图，线阵为沿测线的输入剖面和恢复剖面 | 全部 |
| `6_recovery` | 近台站评估区各试验随深度的幅度恢复（100% / 1.0 为理想值）与形态相关系数 r（1.0 为理想值）热图及深度曲线 | 全部 |
| `6_residuals` | 全部观测的毫秒残差直方图、分周期相对 RMS；均值、总体标准差和 RMS 分列，并列保留/全部数量；已有逐轮评估才显示实际迭代曲线 | 全部 |

画不出的图在 figures 报告里列出原因（例如某频带没有测量结果、离线时取不到地形）。
工作流只写 PNG（`plots.formats: [png]`）。图件用统一的期刊风格（衬线字体、刻度朝内、红慢蓝快），输入与恢复异常图使用相同色标范围。

恢复热图支持单深度，检查深度严格递增且各试验对齐，保留负相关与超出 100% 的幅度。色块表示节点采样，不代表可分辨深度区间。剖面明确标为相对一维参考的异常，避免将它误称为非均匀初模的实际更新。

父进程在主求解任务完成时立即发布 `results/main_preview/`，无需等待其他工程；evaluate 步骤把评估过的模型发布到 `results/model/` 时删除它（仅当其中每个文件都仍是 lant 写的原样）。保留模型、求解输入、实际参数；绘图环境可用时给出剖面。`status.json` 明确检查可能未完成，记录绘图失败且保留模型。这是预览，不表示科学验收或完整评估已完成。

`lant clean-work --cache forward` 预览已登记的精确正演 `.npz` 缓存；明确执行时仅在评估完成后删除未改动的已登记文件，不将处理步骤标为缓存已清。默认预览包括 archive、cc、forward。未知或改动过的文件、模型、参数与最终图均受保护。占空间的 clean/cc 缓存同样仅按明确清理请求删除。

主模型预览发布时立即打印路径，并写入新项目 README 的目录说明，使完整 run 报告发布前也可找到。
