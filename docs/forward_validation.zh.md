# 一维正演（surfdisp96）

中文 | [English](forward_validation.md)

lant 里层状模型的相速度都由 DSurfTomo 自己的频散代码计算：起始模型及其拟合、核示意图、合成算例，以及三维反演逐列的深度核。`lant.inversion.dsurf_disp96` 把 DSurfTomo `a651faf` 的 `surfdisp96.f`（Rayleigh 相速度、球面地球、基阶）以及 `CalSurfG.f90` 中的 Brocher 换算和 `refineGrid2LayerMdl` 逐句移植，并保持 Fortran 原有的精度。相速度与原生 `depthkernel` 逐位相同（一份存下的原生输出是测试的一部分）。从 2.0.1 起，深度核改为沿 Brocher 关系求导（每个节点 2 次频散计算，原为 6 次），与原生 Vs、Vp、密度三项偏导的合成相差不超过 0.03 %（[python_backends.zh.md](python_backends.zh.md)）。

2.0 之前，一维计算先是 fwd1d（包装 surfdisp96 的小 Fortran 程序），后来是 disba 0.7.0（双精度的 surf96）；disba 的相速度与 surfdisp96 最多差 5e-6 km/s，足以把个别数据推过 DSurfTomo 的剔除界。

## 接口与模型约定

```python
from lant.inversion.start import forward

c = forward([0, 0.1, 0.5, 2], [1.8, 2.0, 2.6, 3.4], [1, 0.1, 0.5])
```

深度单位 km，Vs 和结果单位 km/s，周期单位 s。只移植了 Rayleigh 相速度（`igr=0`，`iwave` 为 0 或 2）。要求至少两个从 0 开始递增的深度节点、正的 Vs、1–80 个正周期（surfdisp96 的 NP）。周期按从小到大求解，与 DSurfTomo 一样（surfdisp96 每个周期的求根从上一个周期的根出发），结果按调用者给的顺序返回。找不到根的周期返回 0，与 surfdisp96 相同；调用者遇到时停止。

与 DSurfTomo 一样，节点深度和 Vs 用 REAL（单精度）。Vp 和密度在节点上按 Brocher (2005) 换算，再插值到各子层中点。厚度为 `h` 的每个区间分成 `int((h + 1e-4) / (h / minthk)) + 1` 个子层；连同半空间最多 200 层。球面扁平化用半径 6370 km、对数厚度、中点速度缩放和密度幂 −2.275（Rayleigh），半空间按临时 1 km 厚度缩放。求根照 surfdisp96：起始速度为最慢层半空间 Rayleigh 速度的 0.95 × 0.90，步长 `dc = 0.005`（REAL 常数：0.004999999888 km/s），Dunkin 复合矩阵频散函数用双精度，用 Neville 插值和二分细化到相对 1e-6；结果存为 REAL。

## 限制

没有移植高阶模式、群速度、Love 波、流体层和各向异性；lant 只用基阶 Rayleigh 相速度。求根失败时，原生 DSurfTomo 打印 "improper initial value in disper" 并以 0 继续；lant 给出这条信息后停止。

## 来源

- [DSurfTomo a651faf](https://github.com/HongjianFang/DSurfTomo/tree/a651faf92758c1b4c5eb494a4865ef797b9a4b1c)，MIT：`surfdisp96.f` 与 `CalSurfG.f90`（`depthkernel`、`caldespersion`、`refineGrid2LayerMdl`）。
- Herrmann (2013)，*Seismol. Res. Lett.* 84(6)，1081–1088，doi:10.1785/0220110096：Computer Programs in Seismology，`surfdisp96.f` 由其中的 surf96 修改而来。
- Brocher (2005)，*BSSA*，[10.1785/0120050077](https://doi.org/10.1785/0120050077)：DSurfTomo 使用的 Vp/Vs/密度经验关系。
- Schwab & Knopoff (1972)，*Methods in Computational Physics* 第 11 卷，以及 Biswas (1972)，*PAGEOPH* 96，61–74：球面扁平化与密度换算，见 surfdisp96 中 `sphere` 的引用。
