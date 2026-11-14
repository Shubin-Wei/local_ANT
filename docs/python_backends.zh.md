# DSurfTomo 的 Python 移植

中文 | [English](python_backends.md)

lant 的三维反演是 DSurfTomo `a651faf`（Fang et al. 2015）逐句移植的 Python/Numba 版：主程序（`main.f90`）和 CalSurfG 的正演部分，针对基阶 Rayleigh 相速度。它按 DSurfTomo 自己的文件名和格式读写文件，因此一个工程目录仍可当作 DSurfTomo 工程来读。2.0 以前 lant 调用 Fortran 程序；移植版与原生程序的检验见下文。

## lant 怎样运行它

每个反演工程，lant 在工程目录中启动 `python -m lant.inversion.solver --threads N`。它在 `DSurfTomo.in`、`surfdata.dat`、`MOD` 上运行 `dsurf_main`，然后写出 `solver.json`（实现的身份：依赖包版本和这些模块的哈希；线程数；每次迭代各模块的秒数；峰值内存）和 `DWS.dat`（最后一轮每个节点的 DWS，用于棋盘恢复数字）。只有带当前实现的 `solver.json` 的工程才算完成。线程数按 `resources.cores` 定（[parameters.zh.md](parameters.zh.md)）。预检、棋盘噪声和残差在 lant 进程内计算不反演的第一轮（`lant.inversion.solver.evaluate`）。

单独运行时，`python -m lant.inversion.dsurf_main DSurfTomo.in` 像 Fortran 程序一样在当前目录反演；`--threads N`，默认按 lant 的规则（空闲核的 80 %）。

## 正演部分

| 模块 | 对应的 Fortran 子程序 |
|---|---|
| `dsurf_input` | `main.f90` 读取 `DSurfTomo.in`、`surfdata.dat`、`MOD` 的部分，`delsph.f90` |
| `dsurf_bspline` | `gridder`、`bsplrefine` |
| `dsurf_fmm` | `travel`、`fouds2`、`addtree`、`downtree`、`updtree`、`bilinear` |
| `dsurf_source` | `CalSurfG` 中震源网格加密一段、`srtimes` |
| `dsurf_rays` | `rpaths` |
| `dsurf_rows` | `CalSurfG` 中的行组装（Brocher 链式关系、`ftol`、列顺序） |
| `dsurf_disp96` | `surfdisp96.f`（Rayleigh 相速度、球面、基阶：`sphere`、`gtsolh`、`getsol`、`nevill`、`half`、`dltar4`、`var`、`dnka`、`normc`），以及 CalSurfG 中的 Brocher 换算和 `refineGrid2LayerMdl` |
| `dsurf_kernels` | `depthkernel`、`caldespersion`（在 `dsurf_disp96` 上的逐列循环） |
| `dsurf_forward` | 周期/震源循环；节点核来自 `dsurf_kernels` 或由调用者提供 |

CalSurfG 中写死的设置保持不变（`gdx=gdz=8`、`sgdl=8`、`sgs=8`、地球半径 6371 km、`fom=1`，因此未用到的一阶 `fouds1` 没有移植）。`real=np.float32` 逐步重现 Fortran REAL 运算；`np.float64` 用同样的 REAL 输入、以双精度执行同一算法。参考程序的行为照原样重现，不作修正：射线起点距震源不足 `2*dpl` 或位于加密震源单元内时该行为空；连到震源的最后一段不积分；跨 B 样条单元边界时，前一点在旧单元中的基函数值加到新单元的顶点上；震源四个种子点用角距离除以速度；近震源的接收点走时为 `2d/(vs+vr)`；观测走时用单精度 haversine；`surfdata.dat` 中重复出现的周期会覆盖之前的震源。只接受基阶 Rayleigh 相速度。各震源用 `prange` 并行计算，每个震源有自己的走时场和堆。lant 加了一处：某个震源的粗网格推进，在其接收点走时和射线要读的节点都算完（接收点所在格子的角点，再加最慢速度下两个格子对角线的余量）之后就停止；没算的节点都比这些晚，所以读到的每个值都与原生相同，走时部分的时间约减半。

## 反演

| 模块 | 对应 `main.f90` 的内容 |
|---|---|
| `dsurf_main` | 主程序：读入、合成数据、迭代循环、更新、输出（模块说明逐步列出每一步及其行号） |
| `dsurf_system` | 残差、`getpercentile` 与四分位剔除、DWS、平滑行、残差统计 |
| `dsurf_lsmr` | `dp = kind 4` 的 `lsmrModule.f90`、`aprod.f90`、`lsmrblas.f90` |
| `dsurf_format` | gfortran 的表控输出与 `Fw.d` 输出 |

一切按 Fortran 的做法重现，包括其特殊之处：`main.f90` 声明为 REAL 的地方（残差、剔除、矩阵、LSMR、更新）都用 REAL 运算；LSMR 为单精度，`atol = btol = 1e-6`、`conlim = 100`、`itnlim = 400`，并对最近 10 个向量做局部再正交化（`scipy.sparse.linalg.lsmr` 为双精度且没有局部再正交化，因此不用它）；残差 `obst(i) - dsyn(i)` 按下标相减，前者为文件顺序、后者为 CalSurfG 顺序；剔除 `threshold x [q25, q75]` 之外的残差，四分位数带符号，被剔除的行保留为零行；参数盒子表面上的平滑行只有对角元；速度改变量截断在 ±0.5 km/s，模型截断在速度上下限内；第 *k* 次迭代记录的残差属于第 *k* 次更新之前的模型；合成数据与 CalSurfG 的 `synthetic` 一样在 5 × 5 而不是 8 × 8 加密的传播网格上计算；`velmap2dRc.dat` 沿用 Fortran 错开一格的下标。`spfra` 只决定数组大小：非零元超出时与原生程序一样停止。没有重现的：`noiselevel > 0` 时 gfortran 的随机数序列（改用替代的 Marsaglia 抽样，零噪声消除了随机序列差异，但频散核仍可能使合成观测不同）、屏幕上的射线贴边警告、群速度和 Love 波数据。

节点核照 DSurfTomo 自己的频散代码、按它自己的精度计算：分层模型、地球扁平化和存下的相速度用 REAL，频散函数和求根用双精度，REAL 常数在双精度表达式里保留单精度值，整数次幂按 gfortran 的展开顺序相乘。相速度和 Vs、Vp、密度偏导与原生 depthkernel 逐位相同（以一份存下的原生输出作测试）。`dsurf_kernels` 对相同的列只算一次，只在 CalSurfG 读取的地方计算偏导（内部各列、最后一个深度以上的节点）。未重现：原生 depthkernel 的各 OpenMP 线程共用 getsol 的 SAVE 变量；这里每次调用各用自己的（即单线程运行时的值）。

## 与原版一致（2.0.0）

在四个早先原生运行的输入上——仓库合成数据生成器的面状（16 台）和线状（24 台）工程、pDSurfTomo `example_syn` 的五个周期（合成标志 1）、一条未公开的 55 台测线——移植版十轮迭代的输出与原生 DSurfTomo `a651faf` 逐字节相同：每一轮的 `.iterNNN` 模型、最终模型、`residualFirst.dat`、`residualLast.dat`，合成算例还有 `DSurfTomo.inSyn.dat`、`velmap2dRc.dat` 和 `Vs_model.real`。日志只有第三行不同，即注明本移植的横幅。2.0.0 之前节点核来自 disba（与 surfdisp96 接近，但不逐位相同）；下面的对比都是那个版本的。

## 本轮共同观测对照

本轮对照使用同一份观测、初模、网格及控制文件进行六轮棋盘格对照，双方均以实测模式读取，核对观测字节及 float32 走时相同。最大／RMS 最终模型差为 0.440／0.0506 m/s；双方 9 线程，原版／Python 墙钟时间 524.07／154.57 s。恢复幅度仍弱，不能以实现一致证明分辨率。实测十轮对照中原版已完成，Python 正在计算，尚不报告实测差异。

![共同观测棋盘格](assets/solver_common_checkerboard.png)

以下为上一阶段各自合成观测的历史对照，不能替代上述共同观测比较。

## 与原生程序的对比

上一阶段对比使用与 README 首图相同类型的正弦棋盘格：Vs 扰动 ±6 %、格子 3 km、16 个台站、32 个周期共 2884 个走时，33 × 35 × 11 网格，六次迭代。两个程序读取同一份 `DSurfTomo.in`、`MOD`、`MOD.true`、`surfdata.dat` 的副本。合成观测分别由各自的正演计算，均不加噪声。对比固定平滑权重 = 阻尼 = 1、8 线程；这些是基准设置，不替代 lant 自动选择资源和正则化的默认行为。

![棋盘格：真模型、原生、Python 和差值](assets/solver_checkerboard.png)

图示 0.60 km 深度，输入、原生与 Python 共用 ±6 % 色标，差值单独使用 m/s 色标。在这个格子大小与固定正则化下，恢复幅度偏弱；实现一致不等于分辨率充分。地图显示台站范围，下面的差异数字则包括所有输出节点和深度。

| 最终模型差异，Python − 原生 | km/s |
|---|---:|
| 最大绝对差 | 4.90e-4 |
| 均方根差 | 4.5632e-5 |

（2.0.0 之前。）Disba 与原生频散核接近，但不逐位相同；最终差异也包含合成走时、剔除集合和非线性迭代路径的影响。此前单独的横向结构算例中，Python 驱动使用原生核时，每轮模型数值以及首末轮残差文件与原生一致；这把移植与频散核差异分开检验。那些早期算例上的 disba 单线程与八线程结果也一致。这些结论只对应所测输入，不是普适误差界。

![各模块每次迭代耗时](assets/solver_timing.png)

| 模块 | 原生（秒） | Python（秒） | 原生 / Python |
|---|---:|---:|---:|
| 深度核（SenK） | 76.23 | 25.84 | 2.95 |
| 走时 | 21.81 | 2.17 | 10.04 |
| LSMR | 0.707 | 0.162 | 4.38 |

耗时取六次迭代去掉最快、最慢后的平均。两者在同一服务器、同为八线程预算；原生 SenK 使用 OpenMP，走时与 LSMR 串行。原生程序未修改，通过 gdb 的入口／返回断点计时（含断点开销），Python 使用内部计时。这些模块加速比不是整条流程的加速比。依赖版本、源码／输入哈希及未取整数字见[对比记录](assets/solver_checkerboard_summary.json)。本次图件没有重跑真实记录。

`tools/plot_solver_comparison.py --case VALIDATION_DIR --out FIGURE_DIR` 从已完成的求解器输出重画图，不运行求解器。验证目录需要 `inputs/`、`native/`、`python/` 和 `case.json`；本地路径与运行命令写在交接记录中。

## 来源与许可

- [DSurfTomo](https://github.com/HongjianFang/DSurfTomo)，固定 `a651faf92758c1b4c5eb494a4865ef797b9a4b1c`，MIT；Fang et al. (2015)，doi:10.1093/gji/ggv080。被移植的程序。
- Herrmann (2013)，*Seismol. Res. Lett.* 84(6)，1081–1088，doi:10.1785/0220110096：Computer Programs in Seismology，DSurfTomo 的 surfdisp96.f 由其中的 surf96 修改而来。
- LSMR：Fong & Saunders (2011)，doi:10.1137/10079687X，与 DSurfTomo 的 `lsmrModule.f90` 相同。
- 对比指标沿用 pDSurfTomo（Zhu et al. 2026，arXiv:2604.11920）。
- [PyGMT](https://www.pygmt.org/)，BSD-3-Clause，调用 GMT（LGPL），用于 `plots.backend: pygmt`。
