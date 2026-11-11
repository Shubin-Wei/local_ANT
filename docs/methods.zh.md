# 方法：原理、公式与依据

中文 | [English](methods.md)

本文逐步说明 `lant` 每一步计算什么、为什么这样算、依据哪篇文献或哪条理论关系，以及每种做法在什么条件下成立。操作步骤见 [tutorial.zh.md](tutorial.zh.md)，每个参数的默认值和依据类别见 [parameters.zh.md](parameters.zh.md)，合成检验的数字见 [validation.zh.md](validation.zh.md)。

文中“经验值”指没有文献给出取值、由本项目在合成数据上选定的数；它们在 `parameters.zh.md` 里都有标注。

---

## 1. 总体思路

背景噪声成像的物理基础是：在近似各向同性的扩散噪声场中，两台记录的互相关与两台之间的格林函数有确定的关系，其相位给出台间面波的走时。对垂直分量，互相关中最强的是基阶 Rayleigh 波。流程因此分为三段：

| 阶段 | 步骤 | 得到什么 |
|---|---|---|
| 从记录到经验格林函数 | clean → cc → stack | 每个台对的叠加互相关，以及若干时间分组的子叠加 |
| 从格林函数到频散 | reference → disp | 每个台对一条基阶 Rayleigh 波相速度曲线及其刀切误差 |
| 从频散到模型 | inversion → run → evaluate | 三维横波速度模型、走时残差、分辨检验 |

整个流程只测基阶 Rayleigh 波相速度，只用垂直分量。群速度只用来定位相位读取的到时，不作为产品。

---

## 2. 台站与距离

### 2.1 台站与记录的对应

台站表与记录头按以下顺序对应：台网和台站都相同；台站号相同；台网与台站连写相同；去掉前导字母或前导零后相同。有歧义的对应不使用，列在报告里。节点仪器在野外搬站时，台站表给出序列号（或台站名即记录头代码）和布设起止时间，记录按时间窗分给各点位。只给日期时，按 `data.day_label` 的时区取整天（起始日 0 时至结束日 24 时）。

### 2.2 距离

台站坐标是大地坐标（WGS84 椭球上的经纬度）。台间距离和方位角按 WGS84 椭球用 Vincenty (1975) 反解计算（`lant.signal.geodesic_km`）。若按半径 6371 km 的球计算，距离误差随路径方位变化：纬度 30° 处南北向偏长 0.30 %、东西向偏短 0.18 %；在赤道附近可达 +0.56 %。这一误差与测量精度同级，且形态与 2θ 方位各向异性相同，因此必须用椭球距离。

相速度测量使用椭球距离：c = Δ / t。DSurfTomo 内部用球面距离把输入的相速度换成走时并做射线追踪；由于输入的是速度而不是走时，均匀介质下 DSurfTomo 的模型速度就等于输入的速度，不带方位偏差；它内部的几何比真实椭球有不超过 0.6 % 的拉伸，对横向位置的影响远小于网格间距。DSurfTomo 的源码不改。

---

## 3. clean：原始记录的整理

每个台站、每个处理日（从 `data.day_start_utc_hour` 起的 24 小时）：

1. **采样网格**：记录起点放到该采样率的 UTC 整采样网格上。偏差不超过 0.05 个样点时直接平移，否则 Lanczos 插值；两者都计数。
2. **重叠与重复**：完全相同的重复和数值相同的重叠合并；数值不同的重叠不择一保留，重叠段作为缺数。
3. **缺数**：缺数保持为缺数，不插值、不补零。后面的分段只使用无缺数的段。
4. **去仪器响应**（只在给了 `data.response` 时）：按 Bensen et al. (2007) 去响应到速度，预滤波带由互相关频带推出。不给响应文件时，记录按原样使用；这只在全部台站为同一型号时成立（互相干消去振幅响应，同型号仪器的相位响应在台对之间相消）。
5. **1 s 周期信号检验**：部分数据采集器每秒产生一次与 UTC 同步的相干噪声（Takagi et al. 2015），它在互相关零时延附近形成虚假的相干信号。检验方法是把记录按 UTC 秒折叠：设一秒有 p 个样点、共 M 秒，每秒先去均值，再逐相位平均得到折叠波形 s(j)。若记录中没有这样的信号，折叠波形的方差只是记录方差的 1/M。检验统计量

   F = [M Σⱼ s(j)² / (p − 1)] / [(Σ 残差²) / ((p − 1)(M − 1))]

   无信号时其期望为 1，与自由度为 (ν, ν(M − 1)) 的 F 分布的 0.1 % 分位数比较。白噪声 ν = p − 1；实际噪声的相邻样点相关，独立值更少，取等效自由度 ν = (p − 1) / Σₖ ρ(k)²，ρ 为一秒内的自相关（相关变量平方和的等效自由度）。不做这一修正时，250 Hz 采样的红噪声和带限噪声误报 7–15 %（标称 0.1 %）；修正后这些噪声 200 次试验均无误报，每秒一个 0.5 σ 的脉冲 200 次均检出（[validation.zh.md](validation.zh.md) 第 8 节）。一个台站超过一半的处理日检出时，cc 步骤对该台逐小时减去折叠波形。

清洗后的日文件写入 `work/01_clean/archive`（`clean.write_archive: false` 时不写，后面的步骤直接读原始文件）。每个台站日记录其输入文件和设置的指纹；再次运行时，输入和设置都没变的台站日直接复用。

---

## 4. cc：逐小时互相干

### 4.1 预处理与分段

每个台站日零相位低通（8 阶 Butterworth，截止于新奈奎斯特频率的 0.8 倍）后降采样到 `correlation.fs`。按段长 L、步长 L/2 分段（50 % 重叠，Seats et al. 2012）。一段被丢弃，若它含缺数，或其峰值绝对值超过该台站日标准差的 `max_over_std` 倍（NoisePy 的做法，Jiang & Denolle 2020）。每段去线性趋势、加 5 % 余弦窗。同时读入相邻日的记录，跨越归档日界的段不会丢失。

### 4.2 互相干

两台的频谱各自除以自身振幅后相乘：

C_ij(ω) = F_i*(ω) F_j(ω) / (|F_i(ω)| |F_j(ω)|)

即互相干（cross-coherence；Xu et al. 2021 比较了多种归一化后推荐的做法）。只保留 `fmin_hz`–`fmax_hz`，两端各加四分之一倍频程的余弦过渡。同一小时内的段在频域平均后变回时域，保留 ±`maxlag_s`。正时延表示能量从台 i 传到台 j，i、j 按台站名排序。

### 4.3 由速度范围推出的参数

用户只需给出本区相速度的大致范围 [v_slow, v_fast]（`data.velocity_range_kms`），以下各量按理论关系推出（数字写在每一步的报告里）：

| 量 | 关系 | 理由 |
|---|---|---|
| 最低频率 f_lo | 最长台对 = 2 个波长（v_slow 下），再放宽到该频率高斯滤波器的半增益点 | 两台法要求路径不短于约 2 个波长（Luo et al. 2015） |
| 时延窗 | 最长台距 / (0.6 v_slow) + 4 个最长周期 | 群速度可低于相速度；到时之后留两个周期判断信噪比 |
| 段长 | max(2 × 时延窗, 10 个最长周期) | 频域互相关是循环相关，段长不小于时延窗两倍才不混叠 |
| 采样率 | 不超过 `fs_max`、各原始采样率都能整数倍降到的最高整数 | 降采样不插值 |
| 最高测量频率 | 0.4 fs / (1 + √(ln 2 / α)) | 最高频率的高斯滤波器半增益点不超过抗混叠滤波带 |

时延窗和段长随 v_slow 反比增长，计算量和磁盘也随之增长。默认范围 0.1–5 km/s 覆盖软土到硬岩，但在孔径十余千米的台阵上时延窗达 600 s；clean 报告给出每进程内存和逐小时文件的磁盘估计，cc 步骤自动把进程数降到内存装得下的数目。知道本区大致速度时写窄这一范围，是缩短计算最有效的办法。

---

## 5. stack：叠加与台站诊断

### 5.1 叠加

全部小时按有效段数加权平均得到叠加互相关；同时按 `stack.block_hours` 小时一组存子叠加，供刀切误差使用。`block_hours` 取 24、12、8、6、4、3、2、1 小时中能让中位台站至少有 8 个子叠加的最长者（经验值）。

### 5.2 台站诊断

钟差、仪器移位、耦合变差会表现为一个台站的所有台对在某段时间内一起变化（Stehly et al. 2007；Sens-Schönfelder 2008）。诊断方法：

1. 台对的每个子叠加与“该台对其余子叠加的叠加”比较（留一法，避免与自身比较），取相关系数 r 和时移 τ（互相关峰值，样点间抛物线插值），时窗 |t| ≤ 距离 / (0.6 v_slow)。
2. 每个台站、每个子叠加，取它与最近 `partners` 个台站的台对的中位数。
3. r 先做 Fisher 变换 z = atanh(r)（使相关系数近似正态），减去同一子叠加的台阵中位数（噪声水平随时间变化）。
4. 尺度用中位数绝对偏差 σ = 1.4826 MAD（Leys et al. 2013）。门限由 Chauvenet 准则按样本数 N 给出：z 只看低侧，偏离超过 Φ⁻¹(1 − 0.5/N) σ；时移看两侧，超过 Φ⁻¹(1 − 0.25/N) σ，且至少 0.25 个样点。按构造，正常数据中每个量平均约有 0.5 个值超过门限。
5. 时移候选另需一个条件：该台站中位时移的不确定度。低相关的子叠加里互相关峰值位置主要由噪声决定，同一台站各台对的时移相互不一致；真正的钟差使该台站全部台对一起移动。因此取该台站各台对时移的中位数标准误差 se = 1.2533 × 1.4826 MAD / √n（正态样本中位数的标准误差为 √(π/2) σ/√n），要求中位时移与台阵中位数之差同时超过 Φ⁻¹(1 − 0.25/N) × se。

诊断只列候选，不剔除。候选数超过按准则期望的误报数（相关系数和时移各 0.5 个，合计 1 个）时，流程停下由你判断；不超过时只在报告里列出。

---

## 6. reference：参考曲线

参考曲线有三个用途：限定群速度拾取的范围；在短路径上确定相位的整周数；给出测量频带下限和速度轴。它不进入最终的相速度数值。

### 6.1 频率—贝塞尔谱

对每个频率 f，用台对互相关的频谱实部 C(r, ω) 计算

I(f, c) = ∫ C(r, ω) J₀(ω r / c) r dr

（Wang et al. 2019, 式 4），积分用按台距排序后的梯形公式，双精度计算。每个频率只用台距不超过 `fj_rmax_lambda` 个波长的台对：先用全部台对做一遍并拾取曲线，第二遍按这条曲线换算距离。谱上的能量脊对应各模式的相速度。

### 6.2 拾取

| 模式 | 适用 | 做法 |
|---|---|---|
| `max_peak` | 面状台阵 | 每个频率取最强峰；第二强峰低于最强峰 0.7 倍且台对数足够时为“清楚”；取最长一串彼此连续的清楚拾取，在 log f 上拟合三次多项式，拾取范围外取端点值 |
| `track_fundamental` | 线状台阵 | 测线的 F-J 谱只有一个方向的台对，高阶模式或斜入射噪声常更强。从一个起点沿频率跟踪脊线：每个频率在允许范围内取最强的峰（脊线是能量极大）；允许范围是对斜率 d ln c / d ln f 的限制，距上一个拾取每隔一个 1/12 倍频程放宽 5 %，随频率上升不超过 3 %；拾取限于 `data.velocity_range_kms` 之内，并且波长不短于平均台间距的两倍（测线按台间距空间采样，波长更短时 F-J 谱出现混叠的假脊线；不规则采样时决定混叠的是平均采样密度，Landau 1967）。起点自动选择：每隔 1/3 倍频程试一个起始频带，取能跟踪最多频点的；一样多时取速度较低者（同一频率上基阶最慢）。跟踪不到 5 个频点时报错，并给出三条出路：指定起始频带、写窄速度范围、借用一条曲线 |

基阶判断是全流程最关键的人工判断：拾到高阶模式时，整周数会整体错一周，全部相速度系统偏高。reference 步骤因此总是停下，由你看 `fj.png` 确认。

### 6.3 借用曲线

`reference.source: file` 使用别处的曲线（lant 写的 npz，或两列文本“频率 Hz、相速度 km/s”）。只在源曲线自身的有效段内按 log f 插值，不跨越无效段、不外推。

---

## 7. disp：两台法相速度

### 7.1 对称分量与窄带滤波

取正、负时延的平均（对称分量，Bensen et al. 2007），在对数等间隔的频率上（每倍频程 12 个）做高斯窄带滤波

G(ω) = exp[−α ((ω − ω_c) / ω_c)²]，α = 5，

得到各频率的包络。包络除以信号窗之后噪声段的均方根，成为信噪比图像。

### 7.2 群到时

由参考相速度推出参考群速度 U = c / (1 − d ln c / d ln f)（Bensen et al. 2007, 式 7），在 0.65–1.35 U 的范围内用动态规划拾取包络脊线：偏离参考超过 25 % 加罚，相邻频率之间的变化平滑约束 3 %。

### 7.3 相位：相位匹配滤波

在中心频率直接读窄带信号在群到时的相位，当群延迟在滤波带内随频率变化（有频散）时有系统偏差，与 Shapiro & Singh (1999) 指出的群速度系统误差同源。本流程用相位匹配滤波（Herrin & Goforth 1977；FTAN 的做法，Levshin & Ritzwoller 2001）：

Ψ(ω) = ∫ t_g(ω) dω，  Φ(ω_c) = arg Σ_ω X(ω) G(ω) exp{i [Ψ(ω) − Ψ(ω_c)]}，

其中 t_g 为实测群延迟。乘 exp(iΨ) 消去了到时的频散，带内各频率的相位得以对齐。相位走时

t_p = −(Φ − φ₀) / ω_c + N T，  c = Δ / t_p，

φ₀ = π/4 是二维格林函数的远场相位（Lin et al. 2008, 式 7），N 为整周数，T 为周期，Δ 为椭球距离。在精确的 J₀ 互相关上，这一做法的误差中位数在 0.1 % 以内（`validation.zh.md` 第 2 节）。

### 7.4 整周数

整周数在路径为 `dispersion.lam_min`–`anchor_max_lambda`（默认 1–3.5）个参考波长的频段确定：在这样的路径上，参考曲线误差即使达到 10 %，走时误差也不到半个周期。然后在 4 倍加密的频率网格上沿频率递推，用实测群慢度预测相邻频率的相位（AFTAN 的 phtovel 做法），偏离预测超过 0.3 周即视为跳周，停止递推。

### 7.5 误差

把子叠加分成 `njack` 组，每次去掉一组重测（delete-a-group 刀切）：

σ² = (g − 1)/g · Σᵢ (cᵢ − c̄)²，g 为组数。

它只反映测量在时间上的可重复性，不包括参考曲线、整周数和噪声方向性的误差，没有经过校准，DSurfTomo 也不使用它（全部走时等权）。

### 7.6 质量控制

每个点只要有一项不满足就不接受，被拒原因以位掩码记录（`work/05_dispersion/disp.h5` 的 `reason`），报告按原因计数。

| 判据 | 内容 | 默认 | 依据 |
|---|---|---|---|
| `snr_g` | 该频率的最大信噪比 | ≥ 4 | 经验值 |
| `snr_c` | 拾取点的信噪比 | ≥ 6 | 文献常用 5–15（Bensen et al. 2007 等），取较低一端 |
| `lam_c` | 路径不短于自身相速度的这么多个波长 | 2 | Bensen et al. (2007) 用 3；Luo et al. (2015) 表明 1–2 个波长即可；本项目 J0 基准：2–3 λ 误差中位数 0.007 % |
| `dev_pre` | 群速度拾取偏离参考 | ≤ 30 % | 经验值 |
| `side_c` | 竞争包络峰的相对高度 | ≤ 1 | 经验值 |
| `jk_c_rel` | 刀切标准差 / 相速度 | ≤ 2 % | 经验值 |
| `ucons_c` | 实测群速度与由相速度曲线推出的群速度之差 | ≤ 25 % | 经验值 |
| `band_edge` | 相干信号振幅谱在该频率与滤波带内均值之比 | ≥ 0.5（−6 dB） | 频带边缘的通常定义 |
| `min_run` | 连续接受的频点数 | ≥ 6 | 经验值 |
| 跳周 | 递推中相位偏离预测 | ≤ 0.3 周 | AFTAN |

**频带边缘判据的定义**：超出噪声源激发的频带后，高斯滤波器通过的只是相邻频率的能量，测出的是别的频率的速度，而时域信噪比仍然很高。判据比较的是相干信号的振幅谱：每个台对的信号窗能谱减去后段时延（只有非相干噪声）按长度换算后的能谱，全部台对合在一起，再开平方。非相干噪声在两段中相同，相减后消去，因此结论不随叠加时长和噪声水平变化。在一天和三天的合成数据（噪声源在 12 Hz 截止）上，12.7 Hz 的比值为 0.64–0.66，13.45 Hz 为 0.33–0.38，14.25 Hz 为 0.19–0.25；只用一半台对时结论相同。

### 7.7 导出

每个台对只保留最长的一段连续接受的曲线（并列时取平均信噪比高者；Yao 的 EGFAnalysisTimeFreq 手册的建议），写成 DSurfTomo 的输入。周期就是测量频率的倒数，不插值。

---

## 8. inversion：一维起始模型与三维反演

### 8.1 平均曲线与深度节点

每个周期取全部台对相速度的中位数，台对数不少于 `inversion.median_min_pairs`（小台阵取台对数的四分之一，至少 5）。深度节点从 0 起，步长取 max(最短波长 / 4, 深度 / 4)：最短波长限制了能分辨的最薄层，面波的垂向分辨随深度变差。最长波长一半以内的节点参与更新（Foti et al. 2018 第 4.2.3 节的经验探测深度），其下再加一个固定的底部节点（DSurfTomo 不更新最后一个节点）。

### 8.2 一维反演

对 ln Vs 做阻尼、平滑的高斯—牛顿迭代，正演用 DSurfTomo 自己的链路（Brocher 2005 由 Vs 换算 Vp 和密度，surfdisp96 计算频散），通过小程序 `fwd1d` 调用，保证起始模型与三维反演用同一套正演。平滑约束作用于模型本身的曲率（“jumping”，Shaw & Orcutt 1985；Constable et al. 1987），而不是每一步的改变量，因此结果不随迭代次数漂移；按实际节点间距归一化。初值取最短周期相速度的 1.087 倍。

### 8.3 三维反演

DSurfTomo（Fang et al. 2015）由走时直接反演三维 Vs，每轮迭代按当前模型重新追踪射线。lant 写出 `surfdata.dat`（台站对、周期、相速度）、`DSurfTomo.in` 和起始模型 `MOD`，不改动其源码。几点接口约束：

- 周期数：DSurfTomo 最多 60 个，pDSurfTomo 最多 80 个；超出时 inversion 步骤停下并给出合适的 `dispersion.per_octave`。
- 稀疏矩阵预留：DSurfTomo 按 sparsity × 走时数 × 节点数预留非零元，预留不够时它先越界写入再检查，直接以段错误退出。lant 按直线路径穿过的网格数估计填充率，`inversion.sparsity: auto` 取估计值的两倍、不低于 0.05；手工给的值低于估计的 1.2 倍时报错。
- 剔除规则：DSurfTomo 每轮剔除残差超出四分位距 `threshold` 倍的数据。运行前用一轮正演检查起始模型下保留的比例（`checks.precheck`）；保留不到一半或残差中位数偏离超过一个四分位距的工程标为无效、不运行。

---

## 9. evaluate：检验

只给数字和图，不设通过标准。

| 检验 | 做法 | 读法 |
|---|---|---|
| 走时残差 | 起始模型和最终模型的残差均方根、中位数、剔除规则保留的比例 | 三维反演把残差降了多少 |
| 棋盘（Lévêque et al. 1993） | 在一维模型上加 ±10 % 的棋盘（测线上图样关于网格中点反对称，使合成走时残差不偏向一侧，否则 DSurfTomo 按带符号四分位数的剔除规则会去掉大部分数据），按实际台对和周期正演走时、加高斯噪声（`inversion.noise`；auto：主工程最终模型中被 DSurfTomo 剔除规则保留的残差的相对均方根，使合成数据与实际数据反演后剩下的误差一样大——噪声与数据误差相当，Rawlinson & Spakman 2016；拟合残差还包含模型拟合不了的结构，是数据误差的上限估计；这时棋盘在主工程之后运行），再反演；台站覆盖范围内的节点上计算恢复图样与输入图样的相关系数和回归斜率 | 斜率就是该尺度的异常在模型中被压低后剩下的比例。它不仅取决于射线路径，也取决于正则化和所加噪声 |
| 垂向棋盘 | 图样符号随深度翻转 | 不同深度之间的串扰 |
| 起始模型对照（Rawlinson & Spakman 2016） | 起始模型按深度加 ±5 % 的梯度（按走时残差居中），其余不变，再反演 | 两个模型之差与所加差别相当的地方，数据没有约束模型；同时给出对照工程自身的残差历程，残差尚未降到主工程水平时，差别中还含有未收敛的部分 |
| 正则化对照（可选） | 平滑和阻尼同乘 0.5 或 2 | 结构对正则化的依赖 |

---

## 10. 适用条件与局限

- 噪声场近似各向同性。方向性噪声使两台法相速度偏离，本流程不检验、不改正。
- 只测基阶 Rayleigh 波；高阶模式强于基阶时，质控能去掉大部分但不是全部错误点（`validation.zh.md` 第 7 节）。
- Vp 和密度由 Vs 按 Brocher (2005) 的经验关系换算（DSurfTomo 内部固定），实际 Vp/Vs 偏离时 Vs 有系统偏差。
- DSurfTomo 内部坐标为单精度：台距在 0.1 km 以下时距离误差可达 1 %。
- 刀切误差未校准。
- 模型只在棋盘检验能恢复的尺度和深度上可解释；台站覆盖以外的节点没有数据约束。

---

## 参考文献

- Bensen, G. D., Ritzwoller, M. H., Barmin, M. P., Levshin, A. L., Lin, F., Moschetti, M. P., Shapiro, N. M., & Yang, Y. (2007). Processing seismic ambient noise data to obtain reliable broad-band surface wave dispersion measurements. *Geophysical Journal International*, 169, 1239–1260. doi:10.1111/j.1365-246X.2007.03374.x
- Brocher, T. M. (2005). Empirical relations between elastic wavespeeds and density in the Earth's crust. *Bulletin of the Seismological Society of America*, 95, 2081–2092. doi:10.1785/0120050077
- Constable, S. C., Parker, R. L., & Constable, C. G. (1987). Occam's inversion. *Geophysics*, 52, 289–300. doi:10.1190/1.1442303
- Fang, H., Yao, H., Zhang, H., Huang, Y.-C., & van der Hilst, R. D. (2015). Direct inversion of surface wave dispersion for three-dimensional shallow crustal structure based on ray tracing. *Geophysical Journal International*, 201, 1251–1263. doi:10.1093/gji/ggv080
- Foti, S., et al. (2018). Guidelines for the good practice of surface wave analysis. *Bulletin of Earthquake Engineering*, 16, 2367–2420. doi:10.1007/s10518-017-0206-7
- Herrin, E., & Goforth, T. (1977). Phase-matched filters: application to the study of Rayleigh waves. *Bulletin of the Seismological Society of America*, 67, 1259–1275.
- Jiang, C., & Denolle, M. A. (2020). NoisePy: a new high-performance Python tool for ambient-noise seismology. *Seismological Research Letters*, 91, 1853–1866. doi:10.1785/0220190364
- Lévêque, J.-J., Rivera, L., & Wittlinger, G. (1993). On the use of the checker-board test to assess the resolution of tomographic inversions. *Geophysical Journal International*, 115, 313–318.
- Levshin, A. L., & Ritzwoller, M. H. (2001). Automated detection, extraction, and measurement of regional surface waves. *Pure and Applied Geophysics*, 158, 1531–1545. doi:10.1007/978-3-0348-8264-4_11
- Leys, C., Ley, C., Klein, O., Bernard, P., & Licata, L. (2013). Detecting outliers: do not use standard deviation around the mean, use absolute deviation around the median. *Journal of Experimental Social Psychology*, 49, 764–766. doi:10.1016/j.jesp.2013.03.013
- Landau, H. J. (1967). Necessary density conditions for sampling and interpolation of certain entire functions. *Acta Mathematica*, 117, 37–52. doi:10.1007/BF02395039
- Lin, F.-C., Moschetti, M. P., & Ritzwoller, M. H. (2008). Surface wave tomography of the western United States from ambient seismic noise: Rayleigh and Love wave phase velocity maps. *Geophysical Journal International*, 173, 281–298. doi:10.1111/j.1365-246X.2008.03720.x
- Luo, Y., Yang, Y., Xu, Y., Xu, H., Zhao, K., & Wang, K. (2015). On the limitations of interstation distances in ambient noise tomography. *Geophysical Journal International*, 201, 652–661. doi:10.1093/gji/ggv043
- Rawlinson, N., & Spakman, W. (2016). On the use of sensitivity tests in seismic tomography. *Geophysical Journal International*, 205, 1221–1243. doi:10.1093/gji/ggw084
- Seats, K. J., Lawrence, J. F., & Prieto, G. A. (2012). Improved ambient noise correlation functions using Welch's method. *Geophysical Journal International*, 188, 513–523. doi:10.1111/j.1365-246X.2011.05263.x
- Sens-Schönfelder, C. (2008). Synchronizing seismic networks with ambient noise. *Geophysical Journal International*, 174, 966–970. doi:10.1111/j.1365-246X.2008.03842.x
- Shapiro, N. M., & Singh, S. K. (1999). A systematic error in estimating surface-wave group-velocity dispersion curves and a procedure for its correction. *Bulletin of the Seismological Society of America*, 89, 1138–1142.
- Shaw, P. R., & Orcutt, J. A. (1985). Waveform inversion of seismic refraction data and applications to young Pacific crust. *Geophysical Journal of the Royal Astronomical Society*, 82, 375–414.
- Stehly, L., Campillo, M., & Shapiro, N. M. (2007). Traveltime measurements from noise correlation: stability and detection of instrumental time-shifts. *Geophysical Journal International*, 171, 223–230. doi:10.1111/j.1365-246X.2007.03492.x
- Takagi, R., Nishida, K., Aoki, Y., et al. (2015). A single bit matters: coherent noise of seismic data loggers. *Seismological Research Letters*, 86(3), 901–907. doi:10.1785/0220150030
- Vincenty, T. (1975). Direct and inverse solutions of geodesics on the ellipsoid with application of nested equations. *Survey Review*, 23(176), 88–93. doi:10.1179/sre.1975.23.176.88
- Wang, J., Wu, G., & Chen, X. (2019). Frequency-Bessel transform method for effective imaging of higher-mode Rayleigh dispersion curves from ambient seismic noise data. *Journal of Geophysical Research: Solid Earth*, 124, 3708–3723. doi:10.1029/2018JB016595
- Xu, Y., Lebedev, S., Meier, T., Bonadio, R., & Bean, C. J. (2021). Optimized workflows for high-frequency seismic interferometry using dense arrays. *Geophysical Journal International*, 227, 875–897. doi:10.1093/gji/ggab260
- Yao, H., van der Hilst, R. D., & de Hoop, M. V. (2006). Surface-wave array tomography in SE Tibet from ambient seismic noise and two-station analysis – I. Phase velocity maps. *Geophysical Journal International*, 166, 732–744. doi:10.1111/j.1365-246X.2006.03028.x
- Zhu, S., Li, J., Chen, G., Fang, H., & Yao, H. (2026). pDSurfTomo: a high-performance parallel computing package for direct surface wave tomography. arXiv:2604.11920.
