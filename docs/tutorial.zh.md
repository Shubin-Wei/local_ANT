# 教程：从原始数据到三维横波速度模型

中文 | [English](tutorial.md)

本文逐步说明每一步的目的、方法、依据和检查要点。初次运行见 [quickstart.zh.md](quickstart.zh.md)。原理和公式见 [methods.zh.md](methods.zh.md)；参数的默认值和依据见 [parameters.zh.md](parameters.zh.md)；每个输出文件见 [outputs.zh.md](outputs.zh.md)；遇到问题见 [troubleshooting.zh.md](troubleshooting.zh.md)；与 AI 助手一起运行见 [ai_workflow.zh.md](ai_workflow.zh.md)；合成检验的数字见 [validation.zh.md](validation.zh.md)。

全部用 Python 运行；反演是 DSurfTomo 的移植版，它与原版的对比见 [python_backends.zh.md](python_backends.zh.md)。

给 AI 助手的规则文件按 `project.ai_assistant` 写：`codex` 写 AGENTS.md，`claude`（Claude Code）写 CLAUDE.md，`gemini`（Gemini CLI）写 GEMINI.md，只写一份完整规则；默认写 AGENTS.md。加载方式与已有文件处理见 [AI 工作流程](ai_workflow.zh.md)。

每步完成后看 `OUT/REPORT.md`：可找到已有主图及 `results/parameters.yaml`。覆盖、台站诊断和阵列 F-J 谱在反演前就可查看；完整模型成果在评估后保留。缓存仍在 `work/`；文件职责见[输出说明](outputs.zh.md)，选择清理前先查看 `lant clean-work --help` 与清理预览。

## 1. 安装

### 1.1 环境

```bash
git clone https://github.com/Shubin-Wei/local_ANT.git && cd local_ANT
conda env create -f environment.yml && conda activate lant     # Python、ObsPy、GMT 6 和 lant 本身
```

`conda env create` 长时间停在 “Solving environment” 时，改用 `--solver libmamba`，或用 mamba / micromamba。不用 conda 也可以：Python ≥ 3.11 上 `pip install -e .`，另装 GMT 6（≥ 6.1）。只支持 Linux。不需要编译器：lant 的每一部分（包括反演）都是 Python（NumPy、SciPy、Numba）。

### 1.2 检查

```bash
lant check          # 逐项列出 Python 包和 GMT 的状态及缺失项的安装命令
```

### 1.3 计算资源

核数不用人定：`resources.cores: auto` 取别人没在用的空闲核的 80 %，max(1, 0.8 × (可用 CPU 数 − 1 分钟负载))，本次运行自己造成的负载不算；每一步、每个反演工程开始时重算。并行任务的估计峰值还不超过可用内存的一半（其余排队）。单台对诊断 PNG 也按此计划在独立进程中绘制并显示完成数；测量进程池结束后才启动绘图。每一步的报告在 COMPUTING RESOURCES 下写明用了多少核和峰值内存。`resources.cores` 写整数则照用。

最终主图组也按此资源计划在独立进程中运行，每进程一个数值线程。只请求一组或预算仅够一个进程时串行执行。记录单图耗时，失败组保留上次完整图。经验内存估计包含 600 MB 解释器／绘图预算及工作数组，仍保留相同的可用内存余量。

## 2. 准备数据

所需输入只有两项：原始数据目录和台站表。无需整理目录、重命名文件、转换格式、按天切分或合并记录。

### 2.1 原始记录

任意层子目录、任意文件名，ObsPy 能读的格式都可以，可以混用；非数据文件被跳过并列在报告里。已用测试覆盖的格式：miniSEED、SAC、GSE2、SEG-Y（ObsPy 写出的文件）；SEG-2 走同一条读取路径，但没有测试（ObsPy 不能写 SEG-2，手头也没有长时间的被动记录）。没有台站码的记录（SEG-Y、SEG-2、部分 SAC）按文件名或目录名命名；一个文件里同一时间有多条这样的记录（每道一个通道）时无法这样命名，报告里列为未使用。程序从记录头读台网、台站、通道、采样率、起止时间，然后：

| 内容 | 做法 |
|---|---|
| 垂直分量 | 通道代码以 Z 结尾，或 SAC 头 cmpinc = 0；某台只有一个通道时用它并在报告里注明；名字看不出来时写 `data.channel` |
| 一台多种传感器或采样率 | 用数据最多的一种，其余列出 |
| 处理日起点 | 多数长记录开始的 UTC 整点（`data.day_start_utc_hour`） |
| 重复、重叠、采样网格偏移、缺数 | 处理并计数；缺数不填补 |
| 记录头没有台站名 | 用文件名或上级目录名里与台站表对得上的那一段 |

### 2.2 台站表

CSV、空格或制表符分隔的文本（有无表头均可）、Excel、StationXML。列名中英文都认。必需：台站名、经度、纬度（十进制度或度分秒；不支持投影坐标）。可选列：

| 列 | 用途 |
|---|---|
| 高程 | 剖面图 |
| 台网 | 与记录头对应 |
| 仪器型号（instrument） | 没有响应文件时用来确认全部台站同一型号 |
| 序列号（serial）、开始、结束（start、end） | 节点搬站：同一台仪器先后放在不同点位时，按序列号和时间段把记录分给各点位。时间可以只写日期（如 `260913`、`2026-09-13`），按 `data.day_label` 的时区取整天。没有序列号列时以台站名为仪器代码 |
| 测线（line） | 多条测线在一个表里时，`data.line` 选一条；`lant run all --raw … ` 会给每条测线建一个项目 |

名字的对应规则：台网和台站都相同；台站号相同；“台网 + 台站”连写；去掉前导字母或补零后相同。有歧义的不用，列在报告里。对应结果写成 `work/01_clean/stations.csv`，此后台站名都是 `NET.STA`。

**坐标要逐台核对**：一个台的坐标错了，它所有台对的走时都错，而且不会报错。SAC 头里有坐标时，与台站表相差 50 m 以上的台会列出来。跨日期变更线、纬度超过 80°、两台相距不到 1 m 的台阵会被拒绝或指出。

### 2.3 仪器响应

- 有响应文件（StationXML、dataless SEED、RESP）：写 `data.response`，clean 步骤去响应到速度（Bensen et al. 2007），预滤波频带由互相关频带推出。
- 没有：全部台站必须是同一型号。互相干消掉振幅响应；同一型号的相位响应在台对之间大体相消，但同型号的传感器在出厂容差内仍有差别，留下约为走时 0.1–1 % 量级的随机相位误差。台站表里有型号列且出现多个型号而没有响应文件时，程序报错；没有型号信息时，clean 之后请你确认。

### 2.4 磁盘和计算量

- 清洗后的日文件与原始数据大小相当（`clean.write_archive: false` 可不写，但去响应需要写）。
- `02_ccf` 每天约为 台对数 × 24 × (2 × 时延窗 × 采样率 + 1) × 4 字节；clean 报告给出这一估计和每个 cc 进程的内存。这里的内存是当前可用的内存，包括作业或容器的限制（Linux cgroup），可回收的文件缓存计为可用。它的一半放不下每个核一个进程时，cc 同时少跑几个；一个进程也放不下时，clean 给出警告，cc 停下。可以释放内存或换更大的机器；写窄 `data.velocity_range_kms` 或用 `correlation.max_distance_km`（只算这个距离内的台对）也能减少用量，但只能依据对本区的独立了解来改，不能为了让程序跑得下而改。
- 运行完成并评估后，`lant clean-work -c CONFIG` 列出可重建的缓存（清洗后的日文件、逐小时互相关）及其大小，`--execute --yes` 才删除。成果和报告仍可读；之后要重建 cc、stack 或重画叠加时间图，需要重新生成这些缓存。
- 反演时间取决于网格节点数和“台站 × 周期”数；棋盘和对照各是一个工程；核数和内存允许多少就同时跑多少。

## 3. 配置

最短的方式不用写配置：

```bash
lant run all --raw /data/raw --stations /data/stations.csv --out /data/result
```

它写出 `/data/result/config.yaml`（三个路径和常改项；完整模板用 `lant init --full` 查看）并开始运行。想先看配置再跑：

```bash
lant init -o my.yaml --raw /data/raw --stations /data/stations.csv --out /data/result
lant validate -c my.yaml        # 检查配置和台站表
lant plan -c my.yaml            # 每一步将做什么、用什么值、值从哪里来
```

开始前建议检查以下键：`data.velocity_range_kms`（默认 0.1–5 km/s，覆盖软土到硬岩；对研究区有独立认识时可以收窄，时延窗和 F-J 扫描随之变小，计算更快）、`data.response`、`data.day_label`（默认 UTC）。其余都有默认值或由数据推出。

## 4. 运行

```bash
lant run all -c CONFIG          # 从第一个不成立的步骤开始，停在需要判断的地方
lant status -c CONFIG           # 每一步的状态、原因、下一条命令（--json 供程序读）
lant confirm STEP -c CONFIG     # 看过报告和图之后确认
lant run all -c CONFIG --from disp     # 改了参数后从某一步重跑
```

步骤运行时有一行进度：在整链中的位置、计数（台站日、天、台对批次、求解工程及其迭代、图件）和已用时间；终端里是每秒刷新的进度条，写入日志时每 30 秒一行；`LANT_PROGRESS=off` 关闭。详细输出写进 `.lant/lant.log`（`-v` 同时显示在屏幕上）。整链完成后打印报告和成果的位置。

运行被中断（Ctrl-C、`kill`、终端关闭；在 `nohup` 下关闭终端不会中断）时，求解进程会停止，步骤记为失败，锁被释放；再运行 `lant run all -c CONFIG` 就继续，被中断的求解工程从输入重新运行（半成品保留在 `.lant/backups/`）。即使是 `kill -9` 或内存耗尽，下一条命令也会把该步骤记为失败，而不是拒绝运行。步骤运行期间 `lant status` 显示它正在运行，`lant note` 照常可用。

每一步的报告存到 `.lant/steps/<步骤>.md`（`-v` 时也打印）：做了什么、关键数字、用到的每个参数（值、依据、来源）、要看的文件、留给你判断的事项。`.lant/manifest.json` 记录每一步依据哪些输入和参数成立；输入或参数变了，这一步及其后的步骤自动标为过时，已做的确认失效。没有需要判断的事项时流程不停。

### 4.1 clean：清洗（`01_clean`）

头段扫描按纳秒文件元数据复用未改动文件，旧签名会重扫一次。全部命中时不启动扫描进程池，也不重写缓存。 清洗台站日使用同一文件身份；给了响应时也包含响应来源。修改响应文件会将 clean 标为过时；升级复用键需要重建一次。

**做什么**：读全部文件头；台站表与记录对应；选垂直通道；切成处理日；把样点放到 UTC 采样网格上（偏移小于 0.05 个样点的直接移，其余 Lanczos 插值）；去掉重复和相同的重叠；给了响应文件就去响应；用按 1 s 折叠的 F 检验查每个台站日有没有每秒一次的相干脉冲（Takagi et al. 2015）；由台站分布判断台阵类型。除去响应外不滤波、不归一化。再次运行时，输入文件和设置都没变的台站日直接复用。

**看什么**：`report.txt`（没对上的名字、没有数据的台站、cc 的资源估计）、`0_data_coverage.png`（每台每天的数据量）、`stations.csv`、`pulse_1s.csv`（检出脉冲的台站）。

**留给你判断的**：名字歧义、坐标不一致、分量方向不明、台阵类型不明确（长宽比 3–10）、没有去响应且不知道型号。

### 4.2 cc：逐小时互相干（`02_ccf`）

**做什么**：每个台站日降采样到处理采样率，分成重叠 50 % 的段（Seats et al. 2012），段内峰值超过台站日标准差 `max_over_std` 倍的丢弃（Jiang & Denolle 2020），每段做互相干（两台的谱各自除以自己的振幅，Xu et al. 2021），按小时平均。同时读入相邻日的记录，跨日界的段不丢失。clean 检出 1 s 脉冲的台站先逐小时扣除折叠波形，其余台站不动。

**参数**：采样率、段长、时延窗、频带都由台站分布和速度范围按规则推出（报告里写出代入数字后的推导）。

### 4.3 stack：叠加和台站诊断（`03_stack`）

**做什么**：全部小时按段数加权叠加；同时存下若干小时一组的子叠加，供刀切误差用（`stack.block_hours`）。诊断：每个台对的每个子叠加与该台对其余子叠加的叠加比较（留一法），取相关系数和亚采样时移（Stehly et al. 2007；Sens-Schönfelder 2008）；每个台站取它与最近若干台站的台对的中位数。候选是这两个分布的离群值：相关系数先做 Fisher z 变换，尺度用中位数和 1.4826 × MAD（Leys et al. 2013），门限按 Chauvenet 准则由样本数给出，不是固定数值。时移只在同一台站各台对的时移一致时列出，即中位时移以其在各台对间的标准误差为单位也超过同一门限：钟差使一个台站的全部台对一起移动，而相关系数低的子叠加里时移是分散的。

**诊断只列出，不剔除**。`candidates.txt` 和 `results/figures/2_station_qc.png` 给你看；候选数超过准则允许的偶然个数（1 个）时流程停下，由你决定：`lant confirm stack`（保留）或 `lant confirm stack --exclude-candidates`（剔除），也可以在 `stack.exclude_station_days` 里自己写。

### 4.4 reference：参考曲线（`04_reference`）

**做什么**：计算台阵的频率—贝塞尔谱（Wang et al. 2019；双精度，逐频率只用若干波长以内的台对），拾取能量脊，在 log f 上拟合一条光滑曲线作为参考相速度曲线；再由这条曲线定测量频带下限和速度轴。参考曲线只用来确定相位的整周数和限制群速度拾取，不进入最后的数值。

**看什么**：`3_fj_phase_velocity.png`。**这是全流程最重要的一个判断：拾取的脊是不是基阶。** 高阶模式能量更强时，`max_peak` 会拾到高阶；这时改用 `track_fundamental`（每隔 1/3 倍频程试一个起点，取能沿脊线跟踪最多频点的，拾取限于 `data.velocity_range_kms` 之内；测线上只取波长不短于平均台间距两倍的频率，`3_fj_phase_velocity.png` 中以虚线标出），或缩小 `data.velocity_range_kms`，或用 `reference.source: file` 给一条曲线。报告同时给出实测速度范围与 `data.velocity_range_kms` 的比较。

**结构横向变化强的大台阵**：`reference.zone_km` 把台阵按路径中点分成方形分区，每区各做一条曲线（`3_reference_zones.png`），拾取不清楚的分区用全台阵曲线。

### 4.5 disp：频散测量、质控、导出（`05_dispersion`）

**做什么**：
1. 每个台对取正负时延的对称分量，在对数等间隔的频率上（每倍频程 12 个，Xu et al. 2021）做高斯窄带滤波（同 FTAN，Bensen et al. 2007），沿参考曲线引导拾取群到时。
2. 相位用相位匹配滤波读取（Herrin & Goforth 1977；Levshin & Ritzwoller 2001）：用实测群到时曲线构造相位改正，消除频散造成的包络内相位变化，再在中心频率读相位；远场相位常数 π/4（Lin et al. 2008）。走时 t = −(Φ − π/4)/ω 加整周数，相速度 c = 台距 / t（两台法，Yao et al. 2006），台距按 WGS84 椭球计算（Vincenty 1975）。
3. 整周数在路径为 1–3.5 个参考波长的频段确定（这里参考曲线不可能差一周），再沿频率连续递推。
4. 误差：把子叠加分组做 delete-a-group 刀切，每组重测一遍。误差未做校准，只反映时间上的可重复性。
5. 质控：每个被拒的点记录原因（位掩码），报告按原因统计。判据包括信噪比、刀切误差、路径不短于 2 个自身波长（Luo et al. 2015）、连续段长度、群速度与相速度的一致性等；每个门限的依据见 parameters.zh.md，标“经验值”的没有文献取值。
6. 导出：每个台对只留最长的一段连续曲线，写成 DSurfTomo 的输入。

**看什么**：`results/figures/3_dispersion_curves.png`（全部曲线与参考）、`pair_figures/`（单台对：速度谱、拾取、被拒原因）、`median_vs_reference.txt`。

**不要为了多出点放宽质控。** 门限改了要用 `lant note` 记下理由。

### 4.6 inversion：建反演工程（`06_inversion`）

**做什么**：台阵平均曲线做一维反演得到起始模型（阻尼最小二乘，平滑约束模型本身，Constable et al. 1987）；深度节点步长取 max(最短波长 / 4, 深度 / 4)，到最长波长的一半（经验法则，Foti et al. 2018）；网格间距取测量波长中位数 / 3；建主工程、四个棋盘工程（Lévêque et al. 1993；三个水平，格子为测量波长中位数的 1、2、3 倍，一个随深度变化，见 methods.zh.md 第 9 节）和默认一个对照工程：换起始模型（`start_grad0.05`；Rawlinson & Spakman 2016）；平滑权重减半和加倍（`reg_x0.5`、`reg_x2`）可在 `checks.controls` 里加。正则化平滑的是模型本身，数据项和模型项都取均值（见 methods.zh.md 第 8.4 节）；平滑权重（1）是在面状合成算例上标定的固定默认值（不用阻尼），检验工程用同样的值。三维反演取测量周期中落在每倍频程 `inversion.per_octave`（4）个网格上的那些，即 12 个测量周期每 3 个取 1 个（methods.md 8.3）；最多 60 个周期（DSurfTomo 的 CalSurfG），超出时这一步停下，并给出合适的 `inversion.per_octave`。

**看什么**：`results/figures/4_start_model.png`、`4_kernels.png`（各周期的敏感深度）。**留给你判断的**：起始模型、最深节点、速度上下限是否合理。

### 4.7 run、evaluate、figures

- **run**：运行 Python 版 DSurfTomo（Fang et al. 2015：直接由走时反演三维 Vs，射线随模型更新；每个工程一个 `python -m lant.inversion.solver`）。运行前检查主工程第一轮被剔除规则保留的数据比例（`checks.precheck`）；主工程无效时停下。主工程和对照工程同时运行；棋盘在主工程结束后再运行，因为它们的噪声取主模型的相对拟合残差（`inversion.noise: auto`；数值写在 run 报告和 `work/06_inversion/checkerboard_noise.json` 里）。`inversion.noise` 写成数字时，所有工程同时运行。
- **evaluate**：只给数字和图，不下结论。走时残差（起始模型、最终模型、保留比例）、棋盘恢复（各深度的相关系数和回归斜率，只用剔除后数据雅可比列非零的节点：DWS > 0）、对照模型与主模型的差及各对照工程的残差历程（残差未降到主工程水平的对照尚未收敛）。写在 `work/06_inversion/evaluation/`。
- **figures**：`results/figures/` 下（清单见 outputs.zh.md）：带地形的台站图、分频带的互相关走时图、三个台对的叠加时长图和频散图、频散曲线、起始模型和核函数、深度切片和三维图（面阵）、剖面、棋盘、残差。比例、距离、时延窗和速度范围都随本项目的数据而定；地形用 SRTM，自动从 GMT 数据服务器取（`plots.topography`）。线阵的模型是沿测线的窄条，不画深度切片和三维图，棋盘画成沿测线的输入剖面和恢复剖面。自定义剖面超出模型范围时按模型范围截断并说明截掉多少。

初模以连接的参数节点表示，最后一个值在三维反演中固定，不表示实测界面。核曲线明确标为独立延伸模型的示意。残差图 (a) 以毫秒显示全部观测，表格分别给出均值、标准差和 RMS；(b) 比较不同周期，不是迭代次数。模型解释还须单独查看棋盘试验。

## 5. 线状台阵

一条测线一个项目。与面阵不同的默认值列在 parameters.zh.md 末尾。

1. 测线的 F-J 谱只有一个方向的台对，常被高阶模式或斜入射噪声影响，默认用 `track_fundamental`；务必看 `3_fj_phase_velocity.png`。同一地区有面阵结果时可借用其参考曲线和起始模型（`reference.source: file`、`inversion.start_model.source: file`）；借用后，测线结果与面阵一致不算独立验证。
2. 测线务必写本区的大致速度范围（`data.velocity_range_kms`）：默认的宽范围下测线 F-J 拾取的频带可能很窄，参考曲线以外的高频点会混入偏快的一支（[troubleshooting.zh.md](troubleshooting.zh.md) disp 一节）。
3. 只沿测线方向有分辨；测线两端射线少。噪声斜入射时沿线视速度偏高，本流程不检验这一点。
4. 仪器在不同测线之间搬动时，用台站表的序列号、起止时间和测线列（2.2 节）。

## 6. 假设与局限

- 噪声场近似各向同性。方向性噪声造成的相速度偏差不检验、不改正。
- 只用垂直分量、只测基阶 Rayleigh 波相速度。
- Vp 和密度由 Vs 按 Brocher (2005) 的经验关系换算（DSurfTomo 及其移植版内部固定）；实际 Vp/Vs 偏离时 Vs 有系统偏差。
- DSurfTomo 的坐标是单精度（约 0.4 m）。观测走时和预测走时用同一组坐标；在检验工程中，进入反演的台对（0.5 km 以上）与双精度距离相差不超过 0.2 %。远短于 0.2 km 的台对误差会到百分之一以上，频散步骤的波长规则通常会把它们去掉。
- 刀切误差未校准。
- 模型只在棋盘恢复得出的尺度和深度上可解释；台站覆盖以外的节点没有约束。

## 7. 记录与复现

`输出目录/.lant/`：`config.yaml`（最近一次运行的完整配置）、`config_history/`、`manifest.json`（每步依据什么成立）、`run_log.jsonl`（每步起止、耗时）、`summary/<步骤>.json`（关键数字）、`steps/<步骤>.md`（报告）、`auto_parameters.json`（规则推出的值及推导）、`provenance/`（代码版本、依赖版本）、`../notes.md`（你的决定和理由）。配置文件末尾的 `confirmed:` 段记录每次确认及其依据的结果指纹。复现：用同一份配置在新的输出目录重跑。棋盘检验的合成噪声使用固定随机种子。

## 8. 遇到问题

按步骤排列的现象、原因和处理见 [troubleshooting.zh.md](troubleshooting.zh.md)；每个输出文件的内容见 [outputs.zh.md](outputs.zh.md)。

父进程在主求解任务完成时立即发布 `results/main_preview/`，无需等待其他工程。保留模型、求解输入、实际参数；绘图环境可用时给出剖面。`status.json` 明确检查可能未完成，记录绘图失败且保留模型。这是预览，不表示科学验收或完整评估已完成。

`lant clean-work --cache forward` 预览已登记的精确正演 `.npz` 缓存；明确执行时仅在评估完成后删除未改动的已登记文件，不将处理步骤标为缓存已清。默认预览包括 archive、cc、forward。未知或改动过的文件、模型、参数与最终图均受保护。占空间的 clean/cc 缓存同样仅按明确清理请求删除。
