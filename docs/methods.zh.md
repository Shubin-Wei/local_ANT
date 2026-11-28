# 方法：原理、公式与依据

中文 | [English](methods.md)

本文逐步说明每一步计算什么、为什么这样算、依据是什么。配置里每个数的默认值和依据见 [parameters.zh.md](parameters.zh.md)；怎样运行见 [tutorial.zh.md](tutorial.zh.md)。

lant 做背景噪声面波成像的基本流程。凡是约定取值而不是理论推出的，文中都写明。

## 1. 为什么噪声互相关能得到面波

噪声场为扩散场时（噪声源均匀分布在台站周围，或散射足够强），两台 A、B 记录的互相关对时间的导数正比于两台之间的格林函数及其时间反转（Snieder 2004；Wapenaar 2004；Wapenaar et al. 2010）：

  dC_AB(τ)/dτ ∝ −G(x_B, x_A, τ) + G(x_B, x_A, −τ)。

位于两台连线延长线上的源（稳相区）在叠加中同相相加，其他方向的源相互抵消。地表的源主要激发面波；面波沿地表扩散（振幅 ∝ 1/√r），在垂直分量的互相关中占主导，即 lant 测量的基阶 Rayleigh 波。

横向均匀介质中相速度为 c(ω) 时，相距 r 的两台的互相关谱为（Aki 1957；Tsai & Moschetti 2010）

  C(r, ω) ∝ J₀(ω r / c(ω))。

lant 两处用到这一关系：台阵 F-J 谱（第 4 节）和两台法相位（第 5 节，远场形式 J₀(x) ≈ √(2/πx) cos(x − π/4)）。

**适用条件**：噪声源分布不均匀会使走时有偏差（Tsai 2009），lant 不做改正，结果在噪声接近各向同性时成立。软层覆盖硬岩时高阶模式可能强于基阶，必须检查参考曲线（第 4 节）。

## 2. 数据清洗（`lant clean`）

- **索引**：`data.raw_dir` 下每个文件用 ObsPy 读一次头段，`1_clean/record_index.csv` 列出全部记录。
- **台站对应**：记录头的台网、台站代码与台站表一致（台站表没有台网列时只比台站代码）才归到该台站。不猜名字，对不上的列在日志里。
- **分量**：`data.channel`（默认 `*Z`），即垂直分量。
- **按天**：每台每个 UTC 日放到精确的采样网格上（记录偏离网格超过 5 % 个采样点时才用 Lanczos 重采样）；缺数保持为缺数，不插值、不补零；去中位数。
- **去仪器响应**（只在给了 `data.response` 时）：用 ObsPy 去响应到速度，预滤波在互相关频带 `cc.fmin_hz`–`cc.fmax_hz` 内平坦（Bensen et al. 2007 第 2.1 节）。不给响应文件时按原样使用，这只在全部台站为同一型号时成立：互相干（下一节）消去振幅响应，同型号仪器的相位响应在台对之间相消。
- **降采样**到 `clean.fs_hz`：零相位 8 阶 Butterworth 低通（截止于新奈奎斯特频率的 0.8 倍），再整数倍抽取。降采样后的样点只有在其前后 max(1 s, 50 个样点) 内原始数据都完整时才算有效。
- **输出**：每天一个 HDF5 文件 `1_clean/days/YYYY-MM-DD.h5`，数组为台站 × 样点，没有有效样点处为 NaN；每台每天的有效小时数写在 `availability.csv`。这一步不按振幅去掉数据：一段记录用不用与段长有关，所以在 `lant cc` 里决定（第 3 节）。

## 3. 互相关与叠加（`lant cc`）

1. **分段**：段长 `cc.segment_s`，重叠 `cc.overlap`（Welch 平均；Seats et al. 2012）。含缺数的段、最大绝对振幅超过该台站日标准差 `cc.max_over_std` 倍的段不用（地震和突跳；NoisePy 的 max_over_std，Jiang & Denolle 2020）。它代替 Bensen et al.（2007）的时域归一化（one-bit 或滑动平均）：受干扰的段直接不用，而不是被压平。每段去线性趋势、加 5 % 余弦窗。
2. **互相干**（Xu et al. 2021 比较了高频密集台阵的各种归一化后推荐）：

   C_ij(f) = F_i*(f) F_j(f) / (|F_i(f)| |F_j(f)|)，

   即两台频谱都做谱白化（Bensen et al. 2007），只保留相位差。保留 `cc.fmin_hz`–`cc.fmax_hz`，两端各四分之一倍频程余弦过渡。
3. **子叠加**：每 `cc.substack_hours` 小时的段在频域平均，变回时延 −`maxlag`…+`maxlag`。正时延表示能量从台对的第 1 台传到第 2 台（台名排序）。FFT 算的互相关是循环相关，所以段长必须大于时延窗的两倍。
4. **叠加**：全部段的平均（线性叠加；Bensen et al. 2007）。

不用 tf-PWS（Schimmel et al. 2011）：相位加权按数据改变波形，而这里要测的正是相位。

## 4. 参考曲线（`lant disp` 前半）

参考曲线有三个用途：限定群到时的搜索范围；确定相位的 2π 模糊（整周数）；它有定义的频率就是测量的频率。它的数值不进入测量结果。来源有三种，都是文献常用的（`disp.reference`）：

| 来源 | 内容 | 文献用法 |
|---|---|---|
| `fj` | 本台阵的 F-J 谱 | 密集台阵研究（Wang et al. 2019；Li et al. 2021） |
| `model` | 一维 Vs 模型（`depth_km Vs_km/s`）的理论相速度 | Bensen et al.（2007）、Lin et al.（2008）用模型预测值 |
| `file` | 一条曲线 `frequency_Hz c_km/s`，例如在 F-J 图上人工拾取 | EGFAnalysisTimeFreq 自动模式（Yao et al. 2006）要求的参考文件 |

**F-J 谱**（Wang, Wu & Chen 2019）：

  I(ω, c) = ∫ C(r, ω) J₀(ωr/c) r dr，

c 等于某个模式的相速度时最大。与开源的 CC-FJpy（Li et al. 2021）处理噪声数据的做法相同，C 取叠加互相关频谱的实部，积分用按台距排序的梯形公式（双精度，使用全部有数据的台对）。k ≤ π/Δr 时谱没有空间假频，Δr 为台对距离的采样间隔（Wang et al. 2019 式 14）。

**拾取**（`fj`）：在 `disp.fj_band_hz`（默认全部频率）内每个频率取 `data.velocity_range_kms` 范围内的最大值，用 ln f 的三次多项式拟合 ln c，只在拾取覆盖的频段内有定义（不外推）。不管用哪种参考，F-J 谱都会算出来并与所用参考画在一起供核对：**参考曲线必须沿基阶**（每个频率最慢的那条脊）。拾到高阶时，所有整周数都会错，相速度全部系统偏高。这一点 lant 不替你判断，请看 `3_disp/figures/fj_reference.png`。

参考群速度：U = c / (1 − d ln c / d ln f)（Bensen et al. 2007 式 7）。

## 5. 两台法相速度（`lant disp` 后半）

对每个台对、每个测量频率 f（每倍频程 `disp.per_octave` 个）：

1. **对称分量**：正时延与反转后的负时延的平均（Bensen et al. 2007）。
2. **窄带滤波**：高斯滤波器 G(f′) = exp[−α ((f′ − f)/f)²]（相对带宽不变，FTAN；Levshin & Ritzwoller 2001），α = `disp.alpha`。
3. **群到时** t_g：包络在 r/((1+w)U_ref) 到 r/((1−w)U_ref) 之间的最大值，w = `disp.group_window`；最大值落在窗边上时不算到时。**信噪比**：t_g 处的包络除以窗后两个周期起到时延窗末尾（至少是时延窗最后 15 %）的滤波信号均方根。
4. **相位匹配滤波**（Herrin & Goforth 1977；FTAN）：频谱乘以 exp(iΨ)，Ψ(ω) = ∫ t_g dω，消去到时的频散，再在 f 附近按高斯权求和得到相位 Φ(f)。直接在 t_g 读窄带信号的相位，在滤波带内群延迟随频率变化时会有偏差。
5. **相位走时**：

   t_p = −(Φ − π/4) / (2πf) + N / f，  c = r / t_p，

   π/4 是 J₀ 的远场相位（Lin et al. 2008），N 为整周数，r 为 WGS84 椭球上的距离（Vincenty 1975）。
6. **整周数**：在每一段可用频率的最低频率处，取使 c 最接近参考的 N。相邻 N 的相对速度差为 λ/r，参考误差 e 小于其一半时选择唯一，即**路径不超过 1/(2e) 个波长**（`disp.reference_error`，默认 0.1，即 5 个波长）；起点处路径更长的段不测（原因 `cycle_ambiguous`）。然后用实测群延迟把 N 递推到更高频率（AFTAN）：

   ω_k t_k = ω_(k−1) t_(k−1) + ½ (t_g,k + t_g,k−1)(ω_k − ω_(k−1))，

   在比测量网格密 4 倍的频率上递推；相位跳动超过 `disp.slip_cycles` 周时该段结束。
7. **误差**：对 `lant cc` 的子叠加做逐个删除的刀切法（每次去掉一个子叠加；Efron & Stein 1981），σ² = (g−1)/g Σ (c_i − c̄)²。它只反映随时间的可重复性，不包括参考曲线、整周数和噪声方向性的误差。
8. **质量控制**（门限都在 `disp` 里）：

| 判据 | 默认 | 依据 |
|---|---|---|
| 信噪比 ≥ `snr_min` | 5 | 文献常用 5–15（Bensen et al. 2007 对长周期数据用 10） |
| 路径 ≥ `min_wavelengths` 个波长（λ = c/f，用实测 c） | 2 | 远场近似：Bensen et al.（2007）用 3；Luo et al.（2015）表明相速度 1–2 个波长即可；精确 J₀ 互相关的单元测试：3 个波长以上误差在 1.5 % 以内，5 个以上在 0.4 % 以内 |
| σ ≤ `max_jackknife_rel` × c | 2 % | 约定取值；设为 `null` 即不用 |
| 只留一段连续、至少 `min_points` 个频点 | 6 | 每个台对一段（EGFAnalysisTimeFreq 手册，Yao 等）；长度为约定取值 |

每个被拒的值都记下原因（`measurements.h5` 的位掩码 `reason`，名称见 `lant.dispersion.REASONS`）。`median_vs_reference.csv` 逐频率比较接受值的中位数与参考：恒定偏离好几个百分点说明参考曲线（因而整周数）可能有错。

## 6. 三维反演（`lant inv`）

**数据**：落在过 1 Hz、每倍频程 `inv.per_octave` 个的网格上的周期的接受值（相邻测量频率的高斯滤波器重叠 90 % 以上，少用一些周期信息不减，而计算时间与周期数成正比）。

**起始模型**：每个周期取各台对相速度的中位数（该周期至少要有 max(5, min(50, 台对数/4)) 个台对）。深度节点从 0 开始，步长 max(最短波长/4, 深度/4)，到超过最长波长一半的第一个节点为止（常用的探测深度估计，Foti et al. 2018 第 4.2.3 节），最后一个节点是固定的半空间。一维模型用对 ln Vs 的阻尼、平滑高斯—牛顿反演，平滑作用于模型本身的曲率（Shaw & Orcutt 1985；Constable et al. 1987），结果不随迭代次数漂移。也可以给模型文件（`inv.start_model_file`）。

**正演**：Vp 和密度由 Vs 按 Brocher（2005）换算；相速度和深度核来自 DSurfTomo 的 surfdisp96，逐句移植为 Python（`lant/inversion/dsurf_disp96.py`；Herrmann 2013，Fang et al. 2015）。由模型给参考曲线时用的也是这一套。

**三维直接反演**（Fang et al. 2015，DSurfTomo 的方法，Python 移植在 `lant/inversion/dsurf_*.py`）：未知量是三维网格上的 Vs（水平间距 `inv.grid_km`，默认为实测波长中位数的三分之一：面波的横向分辨率约为一个波长）。每轮迭代在当前相速度图里用快速行进法追踪每个周期的射线，计算深度核，求 Vs 的更新量：

  min (1/N) ‖W(G Δm − r)‖² + (w²/n) ‖L(m + Δm − m₀)‖² + (d²/n) ‖m + Δm − m₀‖²

r 为走时残差，G 为走时导数，W 为剔除规则给的 0/1 权，N 为保留的走时条数，n 为参数个数，L 为整个网格的 7 点 Laplacian（每个节点用它实有的邻点），m₀ 为起始模型，w = `inv.weight`，d = `inv.damping`。与 DSurfTomo 有意不同的三点：

- 平滑和阻尼作用于模型本身而不是每轮更新量（Occam 的 "jumping" 形式；Constable et al. 1987），结果与迭代次数无关；
- 数据项和模型项取均值，同一权重对任何数据量和节点数含义相同；
- 每轮最小二乘问题用双精度 LSMR 解到收敛（Fong & Saunders 2011，列归一化见 Paige & Saunders 1982）。

**剔除**：残差偏离本周期残差中位数超过 `inv.outlier_mad` × 1.4826 MAD 的走时权重置 0（Rousseeuw & Croux 1993；k = 3 为"非常保守"，Leys et al. 2013）。按周期分别算，因为相位误差换成走时误差与周期成正比。

`4_inv/dsurftomo/` 里的工程是 DSurfTomo 的文件格式，原版程序可以直接运行。

## 7. 分辨率检验（`lant checker`）

棋盘检验（Lévêque et al. 1993）：起始模型乘以 1 + a sin(πx/s) sin(πy/s)；给 `--flip-depth` 时该深度以下符号反转。用同样的台对和周期正演走时，加高斯噪声（`checker.noise` 为 auto 时按周期取实际反演剩余残差的稳健标准差，使检验与数据的噪声水平相当；Rawlinson & Spakman 2016），用与实际数据相同的设置反演。`recovery.csv` 按深度给出有数据的节点上恢复扰动对真实扰动的相关系数和回归斜率；斜率即该尺寸的异常能恢复的比例。一个尺寸的结果不能说明其他尺寸，请做多个尺寸。

## 8. 局限

- 假定噪声各向同性；方向性噪声使相速度有偏差（Tsai 2009），不做改正。
- 只测基阶 Rayleigh 波；高阶更强的地方须检查参考曲线。
- 不考虑地形，距离按椭球面计算。
- Vp、密度按 Brocher（2005）随 Vs 变化；实际 Vp/Vs 不同时 Vs 有系统偏差。
- 刀切误差不是完整的误差。
- 模型只能在棋盘检验能恢复的尺寸和深度上解释。

## 参考文献


- Aki, K. (1957). Space and time spectra of stationary stochastic waves, with special reference to microtremors. *Bull. Earthq. Res. Inst.*, 35, 415–456.
- Aster, R. C., Borchers, B., & Thurber, C. H. (2013). *Parameter Estimation and Inverse Problems* (2nd ed.). Academic Press.
- Bensen, G. D., et al. (2007). Processing seismic ambient noise data to obtain reliable broad-band surface wave dispersion measurements. *GJI*, 169, 1239–1260. doi:10.1111/j.1365-246X.2007.03374.x
- Brocher, T. M. (2005). Empirical relations between elastic wavespeeds and density in the Earth's crust. *BSSA*, 95, 2081–2092. doi:10.1785/0120050077
- Constable, S. C., Parker, R. L., & Constable, C. G. (1987). Occam's inversion. *Geophysics*, 52, 289–300. doi:10.1190/1.1442303
- Efron, B., & Stein, C. (1981). The jackknife estimate of variance. *Ann. Stat.*, 9, 586–596. doi:10.1214/aos/1176345462
- Fang, H., Yao, H., Zhang, H., Huang, Y.-C., & van der Hilst, R. D. (2015). Direct inversion of surface wave dispersion for three-dimensional shallow crustal structure based on ray tracing. *GJI*, 201, 1251–1263. doi:10.1093/gji/ggv080
- Fong, D. C.-L., & Saunders, M. (2011). LSMR: an iterative algorithm for sparse least-squares problems. *SIAM J. Sci. Comput.*, 33, 2950–2971. doi:10.1137/10079687X
- Foti, S., et al. (2018). Guidelines for the good practice of surface wave analysis. *Bull. Earthq. Eng.*, 16, 2367–2420. doi:10.1007/s10518-017-0206-7
- Herrin, E., & Goforth, T. (1977). Phase-matched filters: application to the study of Rayleigh waves. *BSSA*, 67, 1259–1275.
- Herrmann, R. B. (2013). Computer programs in seismology: an evolving tool for instruction and research. *SRL*, 84, 1081–1088. doi:10.1785/0220110096
- Jiang, C., & Denolle, M. A. (2020). NoisePy: a new high-performance Python tool for ambient-noise seismology. *SRL*, 91, 1853–1866. doi:10.1785/0220190364
- Lévêque, J.-J., Rivera, L., & Wittlinger, G. (1993). On the use of the checker-board test to assess the resolution of tomographic inversions. *GJI*, 115, 313–318.
- Levshin, A. L., & Ritzwoller, M. H. (2001). Automated detection, extraction, and measurement of regional surface waves. *PAGEOPH*, 158, 1531–1545.
- Leys, C., et al. (2013). Detecting outliers: do not use standard deviation around the mean, use absolute deviation around the median. *J. Exp. Soc. Psychol.*, 49, 764–766. doi:10.1016/j.jesp.2013.03.013
- Li, Z., Zhou, J., Wu, G., Wang, J., et al. (2021). CC-FJpy: a Python package for extracting overtone surface-wave dispersion from seismic ambient-noise cross correlation. *SRL*, 92, 3179–3186. doi:10.1785/0220210042
- Lin, F.-C., Moschetti, M. P., & Ritzwoller, M. H. (2008). Surface wave tomography of the western United States from ambient seismic noise: Rayleigh and Love wave phase velocity maps. *GJI*, 173, 281–298. doi:10.1111/j.1365-246X.2008.03720.x
- Luo, Y., Yang, Y., Xu, Y., Xu, H., Zhao, K., & Wang, K. (2015). On the limitations of interstation distances in ambient noise tomography. *GJI*, 201, 652–661. doi:10.1093/gji/ggv043
- Paige, C. C., & Saunders, M. A. (1982). LSQR: an algorithm for sparse linear equations and sparse least squares. *ACM TOMS*, 8, 43–71. doi:10.1145/355984.355989
- Rawlinson, N., & Spakman, W. (2016). On the use of sensitivity tests in seismic tomography. *GJI*, 205, 1221–1243. doi:10.1093/gji/ggw084
- Rost, S., & Thomas, C. (2002). Array seismology: methods and applications. *Reviews of Geophysics*, 40(3), 1008. doi:10.1029/2000RG000100
- Rousseeuw, P. J., & Croux, C. (1993). Alternatives to the median absolute deviation. *JASA*, 88, 1273–1283. doi:10.1080/01621459.1993.10476408
- Schimmel, M., Stutzmann, E., & Gallart, J. (2011). Using instantaneous phase coherence for signal extraction from ambient noise data at a local to a global scale. *GJI*, 184, 494–506. doi:10.1111/j.1365-246X.2010.04861.x
- Seats, K. J., Lawrence, J. F., & Prieto, G. A. (2012). Improved ambient noise correlation functions using Welch's method. *GJI*, 188, 513–523. doi:10.1111/j.1365-246X.2011.05263.x
- Shaw, P. R., & Orcutt, J. A. (1985). Waveform inversion of seismic refraction data and applications to young Pacific crust. *Geophys. J. R. astr. Soc.*, 82, 375–414.
- Snieder, R. (2004). Extracting the Green's function from the correlation of coda waves: a derivation based on stationary phase. *Phys. Rev. E*, 69, 046610. doi:10.1103/PhysRevE.69.046610
- Tsai, V. C. (2009). On establishing the accuracy of noise tomography travel-time measurements in a realistic medium. *GJI*, 178, 1555–1564. doi:10.1111/j.1365-246X.2009.04239.x
- Tsai, V. C., & Moschetti, M. P. (2010). An explicit relationship between time-domain noise correlation and spatial autocorrelation (SPAC) results. *GJI*, 182, 454–460. doi:10.1111/j.1365-246X.2010.04633.x
- Vincenty, T. (1975). Direct and inverse solutions of geodesics on the ellipsoid with application of nested equations. *Survey Review*, 23, 88–93. doi:10.1179/sre.1975.23.176.88
- Wang, J., Wu, G., & Chen, X. (2019). Frequency-Bessel transform method for effective imaging of higher-mode Rayleigh dispersion curves from ambient seismic noise data. *JGR Solid Earth*, 124, 3708–3723. doi:10.1029/2018JB016595
- Wapenaar, K. (2004). Retrieving the elastodynamic Green's function of an arbitrary inhomogeneous medium by cross correlation. *Phys. Rev. Lett.*, 93, 254301. doi:10.1103/PhysRevLett.93.254301
- Wapenaar, K., Draganov, D., Snieder, R., Campman, X., & Verdel, A. (2010). Tutorial on seismic interferometry: Part 1 — Basic principles and applications. *Geophysics*, 75, 75A195–75A209. doi:10.1190/1.3457445
- Xu, Y., Lebedev, S., Meier, T., Bonadio, R., & Bean, C. J. (2021). Optimized workflows for high-frequency seismic interferometry using dense arrays. *GJI*, 227, 875–897. doi:10.1093/gji/ggab260
- Yao, H., van der Hilst, R. D., & de Hoop, M. V. (2006). Surface-wave array tomography in SE Tibet from ambient seismic noise and two-station analysis — I. Phase velocity maps. *GJI*, 166, 732–744. doi:10.1111/j.1365-246X.2006.03028.x
