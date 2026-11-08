# AI 协作流程：让 AI 用本项目处理你的数据

本文件写给两类读者：**使用者**（想让 AI 助手帮自己从原始数据做到速度模型），和 **AI 助手**（Claude Code、Codex、Cursor 等能在终端里运行命令、读写文件的工具）。AI 助手进入仓库后先读根目录的 `AGENTS.md`，它会指到这里。

## 1. 怎样开始（使用者）

1. 按 [tutorial.md](tutorial.md) 第 1 节装好环境（或者让 AI 照着装：它会运行 `lant check` 看缺什么）。
2. 在仓库目录里启动你的 AI 工具（例如在 `local_ANT/` 下运行 `claude` 或 `codex`）。
3. 把下面这段话发给它，按你的情况改掉尖括号里的内容：

```
请按 docs/ai_workflow.md 用本项目处理我的数据。
- 原始数据：<目录>（<格式，例如 SAC 和 miniSEED 混放；不知道就写不知道>）
- 台站坐标：<文件，或“SAC 头里有坐标”>
- 输出目录：<目录>
- 台阵：<面状 / 线状 / 不确定>，<台数>台，<记录天数>，<仪器型号，是否混用>
- 野外按 <北京时间 / UTC / 其他> 切文件
- 已知问题：<例如 A07 被挖出过、某几天有爆破；没有就写无>
- 目标：<例如 0–1.5 km 的三维 Vs；关心的频带或深度>
- 机器：<核数>核，<内存> GB 内存，<可用磁盘> GB
```

4. AI 会先检查数据、给出参数建议和依据；**凡是下面第 3 节列为“必须确认”的事项，它会停下来问你**。你回答后它继续运行、检查每一步的结果，最后给你一份结果摘要。
5. 全部过程（配置、每一步的日志、参数和依据、检查结果）都留在你的输出目录 `records/` 里，下次可以接着做或复现。

最终得到：`<输出目录>/07_inversion/main/DSurfTomo.inMeasure.dat`（三维 Vs 模型：经度 纬度 深度 Vs）、`06_dispersion/curves/`（频散曲线）、`figures/`（全套图）、`07_inversion/evaluation/`（残差、棋盘、对照）、`records/`（记录）。

## 2. AI 的工作方式（AI 助手必须遵守）

- **只调用项目已有的命令**（`lant ...`）。不另写处理代码，不改 `lant/` 里的源码来迁就数据；发现程序缺陷时停下来告诉使用者。
- **不修改原始数据**，不在仓库目录里写结果。所有输出、记录都在 `project.output_dir`。
- **参数从数据和中间结果推出**，每个不同于默认值的参数都要有依据（数据统计、图、文献），写进 `records/decisions.md`（格式见第 5 节）。没有依据时用默认值；默认值本身也不合适而又推不出来时，问使用者。
- **不为了“结果好看”调参数**：不放宽质控来多出点，不按检验结果反复挑正则化，不删掉“不好看”的台站。台站整台剔除只能基于仪器或数据问题的证据，并须使用者确认。
- **每一步跑完都检查**（第 4 节的判据），发现异常时说明现象、可能原因和依据，再决定调整（有依据、属于可自主决定的）或询问使用者。
- 汇报时只报数字和图，不下“稳定”“可靠”之类的结论；说明哪些部分没有验证或受数据限制。

## 3. 哪些事 AI 可以自己定，哪些必须问

| 类别 | 事项 |
|---|---|
| **程序按规则自动定**（`records/auto_parameters.json`，AI 核对） | 文件格式、垂直分量、台站表写法、名字对应、台阵类型、处理日起点、降采样率、频带、时延、网格间隔、棋盘、剖面 |
| **AI 自主决定**（写入 decisions.md） | 并行进程数、线程数（按内存公式和核数）；`data.format`（有数据文件被跳过时）；深度节点序列（不够深时）；出图的深度；是否写清洗档案（按磁盘） |
| **必须向使用者确认** | 自动推出的台阵类型、处理日起点与使用者的描述不一致时；名字对应有歧义或大量对不上时（`data.rename`）；通道认不出垂直分量时；坐标不一致时；整台剔除任何台站（`exclude_stations`）；野外备注台站（`field_note_stations`）；借用参考曲线或起始模型（用哪个台阵的）；F-J 图上基阶不能确定时；正则化不用默认值时（给出几组的切片、棋盘、残差让使用者选）；需要改质控门限时；混用仪器类型的处理；任何删除已有结果的操作 |
| **AI 不做** | 修改源码、修改原始数据、为了多出点放宽质控、用棋盘或残差反选参数、把没有跑过的检验说成已验证 |

## 4. 各阶段任务与判据

每一阶段：运行命令 → `lant inspect <步骤>` 读关键数字 → 看列出的图（AI 能读图片时直接看；不能时依据数字并请使用者看图）→ 按判据决定继续、调整或询问 → 在 decisions.md 记一笔。

### 阶段 0：环境

```bash
lant check
```
缺 GMT 只影响地图和剖面；缺 CC-FJpy 时只能借用参考曲线；缺 DSurfTomo 不能反演。按提示装好再继续。

### 阶段 1：配置与数据检查

```bash
lant init surface -o <项目>.yaml     # 只填 output_dir、raw_dir、station_file 三个路径
lant validate -c <配置>
lant run clean -c <配置>
lant inspect clean -c <配置>
lant inspect auto -c <配置>
lant inspect data -c <配置>
```

`clean` 自己完成识别（格式、垂直分量、台站表写法、名字对应）并按书面规则推出 `auto` 参数（台阵类型、处理日起点、降采样率、频带、时延、网格、棋盘、剖面），值和依据在 `records/auto_parameters.json`。AI 的任务是**核对**，不是重新推一遍：

| 核对 | 判据 | 不符合时 |
|---|---|---|
| 识别的文件和格式 | `report.txt` 里跳过的文件都不是地震数据 | 有数据被跳过：设 `data.format`（自主） |
| 垂直分量 | 每台用的通道是垂直分量；“not recognisable as vertical”的台 | 问使用者通道含义，设 `data.channel` |
| 名字对应 | 台站表的台都对上了；没有 ambiguous；对上的数目与使用者说的台数一致 | 有歧义或大量对不上：把报告里的例子给使用者看，按其回答设 `data.rename` |
| 坐标 | 没有“differ by more than 50 m”；台站分布图（`lant plot stations`）没有明显离群的台 | **问使用者**以哪个为准 |
| 台阵类型 | 与使用者描述一致 | 不一致时**问使用者**，按回答写 `project.array_type` |
| 处理日起点 | 与使用者说的切文件时区一致（北京时间 → 16） | 不一致时**问使用者** |
| 频带、时延、网格、棋盘 | 规则给出的值落在合理范围（频带下限 0.05–2 Hz；网格节点数 < 200 × 200） | 使用者有明确目标（例如关心更深）时说明规则的限制，改值须使用者同意 |
| 覆盖 | `coverage.png` 大部分台站日接近 24 小时 | 整台很少数据：告诉使用者（不自动剔除） |
| 磁盘 | 可用磁盘 > 1.5 × 原始数据量 | 设 `clean.write_archive: false`（自主） |

另据核数和内存设 `correlation.nproc`（每进程 ≈ 台数 × 17 MB + 台对数 × 0.2 MB）、`dispersion.nproc`、`inversion.threads/njobs`（自主）。

### 阶段 2：互相关与筛选

```bash
lant run cc stack -c <配置>
lant inspect stack -c <配置>
```

判据：`station_day_qc.listed` 占台站日 ≤ 5 % 为正常。同一台被剔超过一半的可评估日，或共同时移持续 ≥ 0.04 s：**列出证据问使用者**是否整台剔除。全台阵某天都偏低不是坏数据（规则已按当日中位数处理）。线状：`station_hour_qc.excluded` 主要集中在布设、回收前后的小时是正常的。

### 阶段 3：陷波

```bash
lant run notch -c <配置>; lant inspect notch -c <配置>
```

看 `notch_report.png`：阴影盖住尖峰即可。谱线超过 15 条时告诉使用者（可能是强工业环境），不自行改门限。

### 阶段 4：参考曲线

```bash
lant run reference -c <配置>; lant inspect reference -c <配置>
```

看 `fj.png`：采用的拾取点（白圈）应在一条连续、最强的能量脊上，曲线随频率平滑下降（正常频散）。

- 拾取点跳到另一条脊（高阶模式）、或低频端没有能量：面阵改 `fj_mode: track_fundamental` 并据图设 `track_start_band_hz`（基阶最清楚的频段）后重跑；仍不清楚 → **问使用者**是否有别的参考曲线（同区面阵、已发表的平均频散）。
- `picks_used` < 10：数据太少，F-J 不可靠 → 问使用者借用参考。
- 线状且同区有面阵结果：**问使用者**是否借用面阵的参考曲线和起始模型（见 tutorial 第 5 节）。

### 阶段 5：频散

```bash
lant run disp1 disp2 -c <配置>
lant inspect disp1 -c <配置>; lant inspect disp2 -c <配置>
```

判据：
- `reference_pass2.png`：第二遍参考在第一遍分布中间、平滑；与 F-J 相差 ≤ 5 %。差得更多 → 回到阶段 4 查 F-J。
- 每频点相速度点数（`phase_points_per_frequency`）：中间多、两端少，没有孤立的零。某频段突然为零：多半是那里没有参考（`frequencies_without_reference`）或陷波。
- `pair_figures/` 抽查 20 张以上：拾取点在能量脊上。系统性地跳到别的脊 → 参考曲线有问题，回阶段 4。
- `precision`：冗余路径相对差（有冗余路径时）一般 1–3 %。
- 点数少：先看 `candidate_points_before_phase_qc`。候选就少 → 数据（记录时长、台距）限制，告诉使用者；**不放宽质控**。

### 阶段 6：反演工程

```bash
lant run inversion -c <配置>; lant inspect inversion -c <配置>
```

网格间隔、棋盘、垂向棋盘已由规则定（`lant inspect auto`），深度节点按“最长波长的一半”自动选取，序列不够深时自动加深（`inspect inversion` 的 `sequence_extended_by`）。AI 只需核对，并处理：

| 参数 | 做法 |
|---|---|
| `weight` / `damp` | 默认（面阵 1.0/1.0，线状 1.5/1.5）。使用者要求比较时，建 0.5/0.5、1.0/1.0、1.5/1.5 三组（`checks.controls` 里加 `reg_w…_d…`），把切片、棋盘、残差给使用者选；不用 L 曲线 |
| 深度节点数 | 加密节点（`fine_step_km`）后超过 DSurfTomo 的上限时 inversion 会报错：增大 `fine_step_km` 或减小 `fine_zmax_km`（自主） |

判据：`vs1d_fit.txt` 的 rms < 0.03 km/s；起始模型无锯齿；`main/DSurfTomo.in` 周期数 ≤ 80。

### 阶段 7：运行与检验

```bash
lant run run -c <配置> --what main
lant run run -c <配置> --what checks
lant run evaluate -c <配置>; lant inspect evaluate -c <配置>
```

运行时间长（大数据量每个工程数小时）：先跑主工程，确认没问题再跑检验工程；按内存设 `njobs`。判据：

- 残差（`residuals.iterations`）：起始到最终明显下降，最后两轮变化 < 1 %。没有下降 → 起始模型或数据有问题，查 `rejected_by_rule`（第一轮剔除超过 30 % 说明起始模型整体偏离数据）。
- 最终相对残差与数据误差同量级；远小于数据误差（< 1/3）提示过拟合，告诉使用者。
- 棋盘：报告每个尺度各深度的斜率和相关系数；不评判“通过”。
- 对照：报告差值 rms 和最大值所在深度。

### 阶段 8：出图与汇报

```bash
lant plot -c <配置>
```

据台站分布给 `plots.sections`（面阵：穿过台阵中心、沿主要构造或台站密集方向各一条），切片深度取可更新节点。汇报给使用者：数据量（台站、天数、台对、曲线、走时）、每步剔除了多少及原因、参数及依据（指向 decisions.md）、残差、棋盘和对照的数字、图的位置、没有验证或受限的部分。

## 5. 记录格式

`records/decisions.md` 每个决定一行：

```
| 2026-11-06 | clean | data.day_start_utc_hour | 16 | 88 个长记录都从 16 h UTC 开始（records/summary/clean.json） | AI |
| 2026-11-06 | stack | data.exclude_stations | + YY.A07 | 18 天里 15 天被剔（03_stack/station_day_qc/station_day_qc.txt），野外记录被挖出 | 使用者确认 |
```

`records/` 里的 `config.yaml`、`run_log.jsonl`、`summary/*.json` 由程序自动写，不用手工维护。

## 6. 异常与处理

| 现象 | 可能原因 | 查什么 | 处理 |
|---|---|---|---|
| clean：台站表的名字对不上或有歧义 | 两边写法差得太远、多个台网同号 | report.txt 里举的名字 | 给使用者看例子，按回答设 `data.rename` |
| clean：大量 interpolated | 记录时间不在采样网格上（授时问题） | station_days.csv | 告诉使用者 |
| cc：`segments_used` 远低于 `segments_gapfree` | 强瞬态干扰多、或 max_over_std 不适合 | 日文件属性 | 告诉使用者；改门限须确认 |
| stack：某台连续多天被剔 | 钟差、仪器移动、耦合 | station_day_qc.png | 问使用者是否整台剔除 |
| reference：F-J 拾取跳脊 | 高阶模式、斜入射噪声 | fj.png | 换 fj_mode / 起始频段（自主，有图为据）；仍不清楚则问 |
| disp2：某频段没有点 | 那里没有参考曲线 | reference_pass2.txt | 说明原因；不改门限 |
| disp2：拾取系统性偏到别的脊 | 参考曲线错 | pair_figures | 回参考曲线阶段 |
| inversion：报“depth nodes … DSurfTomo allows” | 加密后深度节点太多 | inspect inversion | 增大 `fine_step_km`（自主） |
| run：FAILED（output present） | 上次中断 | 工程目录 | 问使用者后删除旧输出 |
| evaluate：第一轮剔除 > 30 % | 起始模型整体偏离 | evaluation.json | 告诉使用者；线状借用起始模型时常见 |

## 7. 示例对话（节选）

> **使用者**：请按 docs/ai_workflow.md 处理 /data/me/raw，台站表 /data/me/sta.csv，输出到 /data/me/out。面阵，100 台，SAC 格式，野外按北京时间切文件，B12 被挖出过。
>
> **AI**：（运行 `lant check`、建配置、`lant run clean`、`lant inspect data`）数据：98 台有记录，2 台（YY.B40、YY.B41）没有数据；台间距 0.3–9.8 km，中位 4.1 km；记录从 16 h UTC 开始，与北京时间切文件一致，`day_start_utc_hour` 取 16。按台距，最长周期取 1.6 s（p90 台距 8.1 km、2 个波长）。B12 放进 `field_note_stations`。这些已写入 decisions.md。请确认：YY.B40、YY.B41 没有数据是否正常？
>
> **使用者**：正常，那两台没回收。
>
> **AI**：（继续运行 cc、stack）台站日剔除 21 个（1.3 %），其中 YY.B12 从第 6 天起持续偏低、按野外备注规则截断……
