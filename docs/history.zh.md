# 版本记录

中文 | [English](history.md)

更早的版本：见 git 标签（`v1.5.0`、`v2.0.0`–`v2.0.2`）及其文档。

## 3.1.0（2026-11-27）

**新增 `lant clean` 一步**，在 `lant cc` 之前：每个台站的每个 UTC 日清洗后（保留缺数、去中位数、给了响应时去响应、低通降采样）写到 `1_clean/days/YYYY-MM-DD.h5`，并写出 `availability.csv` 和 `figures/availability.png`（每台每天的数据小时数）。`lant cc` 读这些天的数据；改设置、台站或天数后重跑 `lant cc` 不必再读原始记录。

**群速度**：`lant disp` 测每个台对的 Rayleigh 波群速度（瞬时频率处的包络极大），有各自的刀切误差和质量控制；`dispersion.csv` 有 `c_kms`、`sigma_kms`、`u_kms`、`u_sigma_kms` 四列。`lant inv` 可反演相速度走时、群速度走时或两者（`inv.data`，默认 `phase`）。

- `cc.fs_hz` 改为 `clean.fs_hz`；`record_index.csv`、`stations_used.csv` 和台站图在 `1_clean`。
- 结果文件夹：`1_clean`、`2_cc`、`3_disp`、`4_inv`、`5_checker_*`。
- 建议的采样率和频带上限（`lant init`）：`clean.fs_hz` 取每个原始采样率都能整数倍抽取到的最高整数采样率，`cc.fmax_hz` = 0.4 × `clean.fs_hz`（奈奎斯特采样定理）。3.0.0 的规则（台对距离的平均间隔）在台对很多时会给出几百赫兹。
- `lant cc` 对清洗后的数据使用 `data.exclude_stations` 和 `data.start_date` / `end_date`；清洗时的采样率不同，或去响应的预滤波窄于互相关频带时，会停下并说明。
- `lant disp` 新增图 `group_curves.png`、`group_rejections.png`。`lant inv` 的 `average_curve.csv`、`start_fit.csv`、`residuals.csv` 增加 `velocity` 列（phase 或 group）。
- `example.yaml`：`lant init` 写出的一份完整配置。
- 某列模型中有 Vs 随深度减小时，频散求根的步长从 0.005 km/s 改为 0.0005 km/s。原步长下可能跳到高阶模式或找不到根（三维反演因此停下）。
- 台对的数据全部来自一个子叠加时，`lant disp` 不再出错停下。

## 3.0.0（2026-11-27）

**怎样运行**：各步依次运行，每步检查后再做下一步：`lant cc`、`lant disp`、`lant inv`，另有分辨率检验 `lant checker`。`lant init` 写出带建议值和规则的配置。每一步写 `params.yaml`（所用设置和推出的值及规则）、`log.txt`、`summary.txt` 和图。不读 2.x 的工程目录。

**方法不变的部分**：读记录与降采样、去响应、分段振幅剔除与互相干、线性叠加；F-J 谱（Wang et al. 2019）；FTAN 群到时、相位匹配滤波、π/4 相位、用群延迟递推整周数；三维反演（DSurfTomo 的正演、射线和深度核，配 lant 的正则化、剔除和求解）。

**改动**

- 恢复向起始模型的阻尼（`inv.damping`，默认 0.1），平滑权重 0.5，两者都按参数个数取均值。
- 参考曲线：来自 F-J 谱（`disp.fj_band_hz` 内每个频率取最大值，三次拟合，不外推）、来自一维模型（`disp.reference: model`，新增）或来自曲线文件。
- 整周数：在路径不超过 1/(2 × `disp.reference_error`) 个波长的最低可用频率处确定（理论：分支间距的一半），代替经验的 1–3.5 个波长锚定频段。
- 刀切：每个子叠加（`cc.substack_hours`，默认每天一个）为一组。
- 质控精简为信噪比、波长、刀切、每个台对一段；去掉了群速度拾取偏离、旁瓣、群速度/相速度一致性和第一道信噪比门限。
- 建议的频带上限按 Wang et al.（2019 式 14）的 F-J 假频条件。
- 直接读原始记录，不再写清洗归档；处理日为 UTC 日。
- 地图用 GMT（能取到时加 GMT 远程地形），其他图用 matplotlib。
- 叠加结果可写成 SAC（`lant export-sac`）。

**删除**：1 秒采集器脉冲的检测与扣除；叠加后的台站诊断（改为子叠加图）；测线专用的 F-J 脊线跟踪和分区 F-J；按序列号和时间分配记录、模糊台站名匹配；工程目录里的确认关卡、清单、溯源文件和步骤报告；起始模型和正则化对照（改设置重跑即可）；合成算例。
