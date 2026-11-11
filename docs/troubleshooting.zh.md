# 常见问题与处理

中文 | [English](troubleshooting.md)

先运行 `lant status -c CONFIG`：它给出每一步的状态、停下的原因和下一条命令。报错时加 `--debug` 可看到完整的 Python 调用栈。下面按出现问题的步骤排列。

## 安装与检查

| 现象 | 原因与处理 |
|---|---|
| `lant: command not found` | 没有在仓库目录里执行 `pip install -e .`，或没有激活装了 lant 的环境 |
| `lant check` 列出 DSurfTomo 或 fwd1d 缺失 | 运行 `lant install-external dsurftomo`（需要 git、make、gfortran） |
| `conda env create` 长时间停在 Solving environment | 加 `--solver libmamba`，或用 mamba / micromamba |

## 配置

| 现象 | 原因与处理 |
|---|---|
| `configuration key ... unknown` | 键名写错，或是旧版本的键；报错信息会说明新版本里改成了什么 |
| `removed in 1.4` | `dispersion.qc.dev_c`、`smooth_c` 已删除（在含噪声和第二震相的检验中不起作用），从配置中删去 |
| `data.start_date must be YYYY-MM-DD` | 日期要写成 `2026-03-18`；YAML 中不加引号也可以 |
| 修改了参数，不知从哪一步重跑 | `lant status` 指出第一个过时的步骤，`lant run all` 自动从那里继续 |

## clean

| 现象 | 原因与处理 |
|---|---|
| `records without a station in the table` | 记录头里的台站名在台站表中找不到。检查台站表；名字写法不同时用 `data.rename` 指定对应 |
| `files hold several records without station codes at the same times` | 每道一个通道、没有台站码的 SEG-Y 或 SEG-2 文件：文件名分不出哪一道是哪个台站。转换成每个台站一个文件（或带台站码的 miniSEED），例如用 ObsPy |
| `stations of the table without data` | 台站表里的点位在给定时间内没有记录。节点搬站时检查台站表的布设起止时间 |
| 台站表的布设时间只有日期（如 `260913`） | 按 `data.day_label` 的时区取整天，起始日 0 时到结束日 24 时；时间不准确时写完整的日期时间（UTC） |
| `N files gave reader warnings` | 记录文件中有不是数据的块（采集器预分配的空白区等），读取时跳过并计数；数据本身正常读入 |
| 停在 clean 等确认“是否同一型号仪器” | 没有给响应文件时，只有全部台站为同一型号才能直接使用。确认无误后 `lant confirm clean`；否则在 `data.response` 中给出响应文件 |
| 磁盘不够写归档 | 设 `clean.write_archive: false`，后面的步骤直接读原始文件（去响应和节点搬站需要归档）。这时 clean 只读文件头，不做 1 s 周期信号检验，报告标为未检验；已知采集器有这种信号时设 `correlation.remove_1s_comb: true` |
| 重跑 clean 很慢 | 输入文件和设置未变的台站日会直接复用；改了 `data.*` 中影响清洗的键（如采样、响应、日界）才会重做 |

## cc

| 现象 | 原因与处理 |
|---|---|
| `N processes instead of correlation.nproc` | 内存不够同时跑设定的进程数，自动减少，结果不变，只是更慢 |
| `One process alone does not fit` | 一天全部台对的逐小时互相关放不进内存。最有效的是把 `data.velocity_range_kms` 写成本区的实际范围（时延窗与最慢速度成反比）；大台阵另用 `correlation.max_distance_km` 只算一定距离内的台对 |
| `02_ccf` 占用磁盘很大 | 同上：时延窗决定文件大小。clean 报告中有估计 |

## stack

| 现象 | 原因与处理 |
|---|---|
| 候选很多，而且多数是“shift” | 先看 `station_qc.png`：相关系数低的子叠加时移测不准。时移候选要求同一台站各台对的时移一致（超过其标准误差），真正的钟差会表现为同一台站连续若干子叠加的一致时移。默认速度范围很宽时诊断时窗很长、频带很宽，噪声比例高；写窄 `data.velocity_range_kms` 后重跑 |
| 不想剔除任何数据 | `lant confirm stack`（保留全部）；剔除候选用 `--exclude-candidates`；剔除特定台站日写 `stack.exclude_station_days` |

## reference

| 现象 | 原因与处理 |
|---|---|
| `the F-J ridge could be followed over only N frequencies` | 线状台阵的 F-J 谱上找不到连续的基阶脊线。看 `fj.png`：若肉眼能看到基阶，在 `reference.track_start_band_hz` 写出它清楚的频带；若脊线在你给的速度范围之外，修改 `data.velocity_range_kms`；都不行时借用邻近台阵的曲线（`reference.source: file`） |
| `only N clear picks that follow each other`（面状） | 台对太少或频带内面波能量弱。看 `results/figures/2_ccf_moveout_*.png` 是否有随台距移动的波列 |
| 拾到的是高阶模式 | 曲线整体偏快、频率高端上翘时很可能是高阶。面状台阵改用 `reference.fj_mode: track_fundamental`，或缩小速度范围，或借用曲线 |
| 测量频带下限比预期高 | 频带下限由最长台对等于 `dispersion.qc.lam_c` 个参考波长推出，这是两台法的有效范围，不是错误 |

## disp

| 现象 | 原因与处理 |
|---|---|
| 接受的点很少 | `.lant/steps/disp.md` 按原因统计了被拒的点。原因主要是 snr 时，叠加时间不够或频带内无能量；主要是 wavelength / path_short 时，台距相对波长太短；主要是 off_band 时，该频率超出了噪声源激发的频带 |
| `median_vs_reference.txt` 显示系统偏离 | 参考曲线可能拾错（高阶、跳周），回到 reference 检查 `fj.png` |
| 测线：参考曲线只覆盖测量频带的一小段，高频的中位数随速度范围设置而变 | 参考曲线拾取范围以外，整周数靠沿频率递推，群速度拾取不受参考约束，可能混入一支偏快的点。在一条实测测线上，默认速度范围 0.1–5 km/s 时参考只覆盖 1.7–2.8 Hz，3–13 Hz 的中位数比写成本区范围 1–4 km/s 时高 4–9 %，多出的点比两次都接受的点快 8–16 %，与 F-J 脊线不符。测线务必写本区的速度范围，并看 `fj.png` 中拾取是否沿脊线到达混叠线 |
| 想放宽质控门限多得点 | 不建议。确需修改时在 `.lant/decisions.md` 中写明理由，并在报告中说明 |

## inversion 与 run

| 现象 | 原因与处理 |
|---|---|
| `periods ... more than DSurfTomo allows` | DSurfTomo 最多 60 个周期（pDSurfTomo 80 个）。按提示减小 `dispersion.per_octave`，或缩短导出的周期范围 |
| `inversion.sparsity = ... reserves too little` | 手工给的稀疏预留比例太小，DSurfTomo 会以段错误退出。改为 `auto` |
| DSurfTomo `segmentation fault` | 先看 `work/06_inversion/main/run.log` 末尾。常见原因：稀疏预留不够（`inversion.sparsity: auto` 已按网格和路径估计）、内存不够（增大 `grid_spacing_deg` 或减小 `grid_refine`）、栈空间不够 |
| 工程被预检标为 invalid | 起始模型下 DSurfTomo 的剔除规则会去掉一半以上的数据，或残差系统偏离。说明起始模型与数据不匹配；检查参考曲线与一维模型拟合（`results/figures/4_start_model.png`） |
| 主工程运行很慢 | 运行时间随网格节点数、周期数和台站数增长。用 `inversion.solver: pdsurftomo`（并行），或加大网格间距 |

## evaluate

| 现象 | 读法 |
|---|---|
| 棋盘恢复斜率很低 | 该尺度的异常在模型中被压到这个比例；用更大格子的棋盘看能恢复的尺度。台站稀、台对少时三维反演对起始模型改动很小 |
| 起始模型对照与主模型相差大 | 两个模型之差接近所加的差别时，数据不约束那里的模型；同时看对照工程的残差是否已降到主工程的水平（未降到时差别里有未收敛的部分） |

## 计算资源的粗略估计

- 磁盘：原始数据大小（归档）+ 每天 台对数 × 24 × (2 × 时延窗 × 采样率 + 1) × 4 字节（逐小时互相关）。clean 报告给出估计。
- 内存：cc 每个进程约为一天逐小时互相关大小的 2.5 倍加一天全部台站的波形；DSurfTomo 随网格节点数与走时数增长。
- 时间（参考，16 台合成算例）：clean、cc、stack、reference 各数秒至数十秒，disp 约半分钟，一个 DSurfTomo 工程（10 次迭代）约 5–25 分钟。
