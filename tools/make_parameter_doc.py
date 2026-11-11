#!/usr/bin/env python
"""Write docs/parameters.md and docs/parameters.zh.md from lant/params.py (English texts) and ZH below.

    python tools/make_parameter_doc.py

Run it after changing a default (lant/config.py) or a description (lant/params.py); a key without a Chinese
text stops it.
"""

import os
import sys

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from lant import config, params  # noqa: E402

KIND_ZH = {
    "required": "必填",
    "literature": "文献",
    "theory": "理论",
    "rule": "规则",
    "empirical": "经验值",
    "choice": "选择",
    "technical": "技术",
}

ZH = {
    "project.name": "项目名（只用于图和记录）",
    "project.array_type": "auto：按台站分布判断（台站分布的长宽比，小于 3 为面状，大于 10 为线状，介于两者之间时在 clean 交给你确认）；或直接写 surface / line",
    "project.output_dir": "全部结果、图和记录都写在这里",
    "data.raw_dir": "原始数据目录，递归搜索；目录结构和文件名不限",
    "data.format": "auto：由 ObsPy 逐个文件识别格式（MSEED、SAC、GSE2、SEG2 等）；或写一个 ObsPy 格式名",
    "data.file_patterns": "尝试读取的文件名模式",
    "data.ignore_patterns": "不读的文件名模式",
    "data.station_file": "台站表，常见形式都可以（CSV、文本、Excel、StationXML）",
    "data.component": "Z：垂直分量互相关（Rayleigh 波）；目前只支持这一种",
    "data.channel": 'auto：按通道名或 SAC 头段找每个台的垂直通道；名字看不出来时写模式，如 "*3"',
    "data.network_default": "记录头段没有台网代码时用的台网代码",
    "data.rename": "自动对不上时手工指定：{头段里的名字: 台站表里的名字}",
    "data.day_start_utc_hour": "处理日从这个 UTC 整点开始；auto：长于 12 小时的记录里有 60 % 以上从同一整点开始就取它，否则取 0",
    "data.day_label": '处理日的标签：utc_start，或固定时差如 "UTC+8"',
    "data.velocity_range_kms": "本区 Rayleigh 波相速度的大致范围 [慢, 快]（km/s），决定时延窗、频带和 F-J 扫描范围。默认值覆盖软土到硬岩；范围写窄一些可以省计算时间（reference 步骤会报告实测速度）",
    "data.response": "仪器响应文件或目录（StationXML、dataless SEED、RESP）：给了就在 clean 去响应（Bensen et al. 2007）。不给时全部台站必须是同一型号",
    "data.line": "台站表有测线列时：本项目处理哪条测线（一条测线一个项目）",
    "data.min_trace_length_s": "短于这个长度的记录片段不用",
    "data.exclude_stations": "整个不用的台站（你已知仪器或耦合有问题）",
    "data.start_date": "使用的第一个处理日（YYYY-MM-DD）",
    "data.end_date": "使用的最后一个处理日（YYYY-MM-DD）",
    "clean.write_archive": "写出清洗后的日文件（大小与原始数据相当）；false：后面的步骤直接读原始文件",
    "clean.nproc": "进程数",
    "clean.encoding": "归档 miniSEED 的编码",
    "clean.min_free_disk_factor": "剩余磁盘小于归档估计大小的这个倍数时拒绝写",
    "correlation.fs_max": "最高处理采样率（Hz）；频带上限为 0.4 fs / (1 + sqrt(ln 2 / alpha))：小孔径、高频台阵可调高",
    "correlation.fs": "降采样后的采样率；auto：不超过 fs_max、且每种原始采样率都能整数倍降到的最高整数采样率",
    "correlation.seg_len_s": "分段长度；auto：max(2 × maxlag（循环相关）, 10 个最长周期)",
    "correlation.seg_step_s": "分段步长；auto：段长的一半（50 % 重叠，Seats et al. 2012）",
    "correlation.maxlag_s": "时延窗；auto：最长台距 / (0.6 × 最慢速度) + 4 个最长周期（最慢的到时之后留一段噪声算信噪比）",
    "correlation.fmin_hz": "互相关频带下限；auto：最长台对在最慢速度下为 2 个波长的频率，再放宽到该频率高斯滤波器的半增益点",
    "correlation.fmax_hz": "互相关频带上限；auto：dispersion.fmax_hz 放宽到其高斯滤波器的半增益点，不超过 0.4 fs（抗混叠滤波）",
    "correlation.max_distance_km": "只对不超过这个距离的台对做互相关（不写：全部台对）；用于台对放不进内存的超大台阵",
    "correlation.smooth_bins": "相干分母在频率上的平滑点数；1 = 不平滑（Xu et al. 2021 的互相干）",
    "correlation.max_over_std": "分段峰值超过该台站日标准差的这个倍数就丢弃（NoisePy 的 max_over_std，Jiang & Denolle 2020）",
    "correlation.remove_1s_comb": "逐小时扣除按 1 s 折叠出的波形（记录器每秒一次的相干噪声，Takagi et al. 2015）；auto：只对 clean 检出该信号的台站做（折叠记录的 F 检验）；true：全部台站；false：不做",
    "correlation.nproc": "同时处理的天数（内存：每个进程装一天的全部台站）",
    "correlation.threads_per_process": "每个进程的 FFT 线程数",
    "stack.block_hours": "刀切用的子叠加长度；auto：24、12、8、6、4、3、2、1 小时中能让中位台站至少有 8 个子叠加的最长者",
    "stack.station_qc.enabled": "台站诊断：列出相关系数或时移与其余不同的台站子叠加（不经你同意不会剔除）",
    "stack.station_qc.partners": "每个台站通过它与最近的这么多个台站组成的台对来比较",
    "stack.station_qc.min_partners": "判断一个台站子叠加至少需要的有数据台对数（不足：不判断，并在报告里写明）",
    "stack.exclude_station_days": '你要剔除的台站日：[台站, YYYYDDD]、[台站, "YYYYDDD-"] 等',
    "stack.exclude_station_hours": "你要剔除的台站小时：[台站, YYYYDDD, 小时]",
    "stack.exclude_candidates": "同时剔除诊断列出的候选（由 `lant confirm stack --exclude-candidates` 设置）",
    "reference.source": "fj：本台阵的频率—贝塞尔谱（Wang et al. 2019）；file：你给的曲线（AFTAN 和 Yao 的 EGFAnalysisTimeFreq 也接受参考曲线）",
    "reference.file": "借用的曲线：lant 写的 npz，或文本 “频率_Hz 相速度_km/s”",
    "reference.zone_km": "不写：全台阵一条参考曲线；写边长（km）：按路径中点落在哪个方形分区，每个分区各做一条 F-J 曲线，用于台阵内结构变化强的情况；分区曲线不清楚时用全台阵曲线",
    "reference.fj_mode": "max_peak：取最强的脊（面状台阵）；track_fundamental：从一个无歧义的起点跟踪脊（线状台阵）",
    "reference.normal_dispersion": "track_fundamental：曲线不得随频率升高；有速度倒转、脊确实上升时设 false",
    "reference.fj_rmax_lambda": "每个频率的 F-J 谱只用不超过这么多个波长的台对",
    "reference.fj_c_nominal": "把 fj_rmax_lambda 换算成距离用的速度；auto：先用全部台对做一遍并拾取曲线，第二遍用这条曲线（共两遍）",
    "reference.fj_cmin": "扫描的最慢速度；auto：0.5 × data.velocity_range_kms 的下限",
    "reference.fj_cmax": "扫描的最快速度；auto：1.5 × data.velocity_range_kms 的上限",
    "reference.fj_dc": "扫描的速度步长；auto：不超过 max(fj_cmin / 200, 扫描宽度 / 3000) 的 1-2-2.5-5 系列步长",
    "reference.fj_nfft": "谱的 FFT 长度；auto：不小于 4 倍时延窗的 2 的幂，至少 8192",
    "reference.fj_fit_fmax": "参与拟合的拾取的最高频率；auto：最长一串清楚拾取的上端",
    "reference.track_start_band_hz": "脊跟踪的起点频带；auto：拾取最清楚的那个倍频程",
    "dispersion.fmin_hz": "最低测量频率；auto：最长台对等于 qc.lam_c 个参考曲线波长的频率",
    "dispersion.fmax_hz": "最高测量频率；auto：0.4 fs / (1 + sqrt(ln 2 / alpha))，即高斯滤波器仍在抗混叠频带内的最高频率；数据实际能到多高由质控决定",
    "dispersion.per_octave": "每倍频程的测量频点数，对数等间隔（Xu et al. 2021）",
    "dispersion.vmin": "速度轴下限；auto：0.6 × 频带内参考曲线的最慢速度",
    "dispersion.vmax": "速度轴上限；auto：1.3 × 最快速度",
    "dispersion.dv": "速度轴步长；auto：不超过 vmin / 100 的 1-2-2.5-5 系列步长",
    "dispersion.alpha": "高斯滤波器 exp(-alpha ((f - fc) / fc)^2) 的相对宽度，同 FTAN（Bensen et al. 2007）；越大越窄",
    "dispersion.phase_constant": "二维格林函数远场相位 π/4（Lin et al. 2008 式 7）",
    "dispersion.min_dist_km": "短于此的台对不测（0：只按 lam_min）",
    "dispersion.lam_min": "路径短于这么多个参考波长就不测（Luo et al. 2015：1 个波长可以够）",
    "dispersion.anchor_max_lambda": "在路径长为 lam_min 到这么多个波长的地方定周数（周数少，参考曲线不至于差一周）",
    "dispersion.njack": "子叠加的分组刀切组数；auto：数据量不少于最满子叠加 25 % 的子叠加个数",
    "dispersion.min_jackknife_groups": "有数据的组少于此数：没有误差估计，该点被拒",
    "dispersion.nproc": "进程数",
    "dispersion.chunk": "每个任务的台对数",
    "dispersion.side_off_bands": "不应用旁瓣规则的频带 'lo-hi,...'",
    "dispersion.qc.snr_c": "拾取到时处的信噪比（峰值 / 尾部噪声均方根）；文献用 4–10（Bensen et al. 2007 对全球数据用 10）",
    "dispersion.qc.jk_c_rel": "允许的刀切标准误差，相对于速度",
    "dispersion.qc.ucons_c": "实测群速度与相速度曲线推出的群速度之间允许的不一致（0 = 关）",
    "dispersion.qc.side_c": "竞争包络峰相对于拾取峰的允许高度",
    "dispersion.qc.dev_pre": "群速度拾取相对于参考的允许偏离（0 = 关）",
    "dispersion.qc.lam_c": "接受的点要求路径不短于这么多个自身波长（Bensen et al. 2007 用 3；Luo et al. 2015 表明 1 可以够）",
    "dispersion.qc.min_run": "接受的点必须连成至少这么多个频点的一段",
    "dispersion.qc.band_edge": "相干信号振幅谱（全部台对合算；信号窗减去后段时延，因此与噪声水平无关）在该频率低于其在滤波带内均值的这个比例时不测（0.5 = −6 dB，频带边缘的通常定义）：超出噪声源激发的频带后，滤波器通过的只是别的频率的能量（0 = 不用）",
    "dispersion.qc.snr_g": "第一道门：一个频率的最大信噪比要达到此值",
    "export.tmin_s": "导出的最短周期；auto：1 / dispersion.fmax_hz",
    "export.tmax_s": "导出的最长周期；auto：1 / dispersion.fmin_hz",
    "export.one_segment": "每个台对只留一段连续曲线，取最长的（同 Yao 的 EGFAnalysisTimeFreq 的拾取方式）",
    "export.curve_figures": "画多少张单台对图：按台距均匀选这么多张，外加异常的（0 不画，-1 全画）",
    "inversion.wave": "rayleigh_phase（只支持这一种）",
    "inversion.median_min_pairs": "一个周期进入一维模型平均曲线所需的台对数",
    "inversion.depth_sequence_km": "深度节点；auto：步长 max(最短波长 / 4, 深度 / 4)，到最长波长的一半为止（经验法则，Foti et al. 2018）",
    "inversion.fine_step_km": "在 fine_zmax_km 以浅每隔这么多 km 加一个深度节点（0 = 不加）",
    "inversion.fine_zmax_km": "见 fine_step_km",
    "inversion.grid_spacing_deg": "网格间距 [纬, 经]；auto：加密后约为台站最近邻距离中位数的一半",
    "inversion.grid_refine": "范围定下后每个网格间隔再细分的倍数",
    "inversion.line_strip_cells": "线状台阵：台站两侧沿垂线方向各保留的网格数",
    "inversion.start_model.source": "invert：由本台阵平均曲线反演一维模型；file：深度、Vs 两列的文件",
    "inversion.start_model.file": "你自己的起始模型（两列：深度 km、Vs km/s）",
    "inversion.start_model.smooth": "一维反演：模型曲率的权重",
    "inversion.start_model.damp": "一维反演：阻尼",
    "inversion.start_model.iterations": "一维反演：迭代次数",
    "inversion.start_model.minthk": "一维正演：每个节点间隔的子层数",
    "inversion.weight": "DSurfTomo 平滑权重（Fang et al. 2015）；在常规范围（小于 10，多数小于 5）内人工取值，未经优化；线阵 1.5",
    "inversion.damp": "DSurfTomo 阻尼（Fang et al. 2015）；在常规范围 0.5–1.5 内人工取值，未经优化；线阵 1.5",
    "inversion.sublayers": "DSurfTomo：算核函数时深度节点之间的子层数（手册：2–5）",
    "inversion.vmin": "Vs 下限；auto：0.6 × 起始模型最慢的 Vs",
    "inversion.vmax": "Vs 上限；auto：1.4 × 最快的 Vs",
    "inversion.iterations": "DSurfTomo 迭代次数",
    "inversion.sparsity": "DSurfTomo：为矩阵非零元预留的比例；auto：按路径和网格估计的填充率的两倍，不低于 0.05（预留太少 DSurfTomo 会以段错误退出；多留只多占内存）",
    "inversion.threshold": "DSurfTomo：残差超出这么多个四分位距的数据权重为零",
    "inversion.noise": "棋盘合成走时所加的高斯噪声，取各走时的比例；auto：主工程最终模型中被 DSurfTomo 剔除规则保留的残差的相对均方根，使合成数据与实际数据拟合后剩下的误差一样大（Rawlinson & Spakman 2016）；这时棋盘在主工程之后运行",
    "inversion.solver": "dsurftomo（不超过 60 个周期）或 pdsurftomo（并行，不超过 80 个周期）",
    "inversion.threads": "每个工程的线程数",
    "inversion.njobs": "同时运行的工程数",
    "checks.checkerboards": "水平棋盘 [{name, cell_km}]；auto：格子为孔径的 0.12、0.18、0.25 倍（线阵：沿测线，测线长度的 0.165、0.24、0.34 倍）",
    "checks.vertical": "随深度变化的棋盘；auto：格子为孔径的 0.18 倍，在最深更新节点的 0.21、0.46、0.71 处变号",
    "checks.amplitude": "棋盘扰动幅度，各深度一维 Vs 的比例",
    "checks.controls": "对照反演：reg_w<w>_d<d> 换正则化，start_grad<f> 换起始模型",
    "checks.precheck.min_kept": "DSurfTomo 的剔除规则在第一次迭代保留的数据少于这个比例时，该工程判为无效",
    "checks.precheck.max_median_iqr": "……或残差中位数超出这么多个四分位距时",
    "checks.iteration_residuals": "同时评估每次迭代的模型（每次一遍正演）",
    "plots.formats": "图的格式",
    "plots.dpi": "png 图的分辨率",
    "plots.topography": "地图的地形底图；auto：从 GMT 数据服务器取 SRTM 1 角秒瓦片（Farr et al. 2007），存在 ~/.cache/lant/srtm1（只请求瓦片名；地图宽度不超过 3 度）；也可写网格文件或 GMT 远程网格；null：不画地形",
    "plots.vs_range": "色标范围；auto：台站下方模型的第 2–98 百分位，取到 0.1 km/s",
    "plots.vs_cpt": "GMT 色标",
    "plots.slice_depths_km": "切片深度；auto：从更新过的深度节点里均匀取至多 8 个",
    "plots.sections": "剖面 [{name, start, end}]；auto：过台阵中心沿两个主轴各一条（面状），或沿测线（线状）",
    "plots.section_zmax_km": "剖面深度（默认：最深的更新节点）",
    "plots.contour_step": "等值线间隔；auto：约为色标范围的十分之一",
    "plots.font_scale": "字号系数",
    "plots.pair_example": "单台对图（叠加时长、频散）画哪些台对：一个台对名或列表；默认在曲线较长的台对中按距离的 20 % / 50 % / 80 % 分位数取三对",
}

SECTIONS = [
    ("project", "project", "项目"),
    ("data", "data: raw data", "data：原始数据"),
    ("clean", "clean: cleaning (`01_clean`)", "clean：清洗（`01_clean`）"),
    ("correlation", "correlation: hourly cross-coherence (`02_ccf`)", "correlation：逐小时互相干（`02_ccf`）"),
    ("stack", "stack: stack and station diagnosis (`03_stack`)", "stack：叠加与台站诊断（`03_stack`）"),
    ("reference", "reference: reference curve (`04_reference`)", "reference：参考曲线（`04_reference`）"),
    (
        "dispersion",
        "dispersion: measurement and quality control (`05_dispersion`)",
        "dispersion：频散测量与质量控制（`05_dispersion`）",
    ),
    (
        "export",
        "export: curves passed to the inversion (`05_dispersion/curves`)",
        "export：交给反演的曲线（`05_dispersion/curves`）",
    ),
    ("inversion", "inversion: DSurfTomo (`06_inversion`)", "inversion：DSurfTomo 反演（`06_inversion`）"),
    ("checks", "checks: checkerboards, controls, precheck", "checks：棋盘、对照、预检"),
    ("plots", "plots: figures", "plots：出图"),
]

HEAD = {
    "en": """# Parameters

[中文](parameters.zh.md) | English

Every key of the configuration file, with its default, what it means and what the default rests on. This file is written by `tools/make_parameter_doc.py` from `lant/config.py` (defaults) and `lant/params.py` (texts); `lant init` writes the same information into the configuration file, and the step reports (`OUT/.lant/steps/<step>.md`) show the value each step used and, for derived values, the derivation with your numbers.

## 1. How to read the tables

Only three keys are required: `project.output_dir`, `data.raw_dir`, `data.station_file`. A value written in the configuration always wins. `auto` means: derived from your data by the rule in the description.

The column "basis" says what a default rests on:

| basis | meaning |
|---|---|
| required | only you know it |
| literature | a value or method taken from the cited source |
| theory | follows from a relation of the method (sampling, filter width, far-field phase) |
| rule | derived from your data by a written rule; the rule's own factors are project choices unless a source is cited |
| empirical | experience of this project; no literature value was found. Change it when your data call for it, and say so |
| choice | yours; there is no right value |
| technical | affects speed, memory or files only, not the result |

Where the chain stops for you: after `clean` (ambiguous names, stations without data, coordinate mismatches, unknown component orientation, an array type between areal and line, responses not removed and instrument models unknown), after `stack` (only when the diagnosis lists candidates), after `reference` (is the picked ridge the fundamental mode?), after `inversion` (starting model, deepest node, velocity bounds). A step with nothing to judge does not stop.
""",
    "zh": """# 参数说明

中文 | [English](parameters.md)

配置文件的每一个键：默认值、含义、默认值的依据。本文件由 `tools/make_parameter_doc.py` 根据 `lant/config.py`（默认值）和 `lant/params.py`（说明）生成；`lant init` 把同样的信息写进配置文件，步骤报告（`输出目录/.lant/steps/<步骤>.md`）给出每一步实际用的值，由规则推出的值还写出代入你的数据后的推导。

## 1. 表格怎么读

必填的只有三个键：`project.output_dir`、`data.raw_dir`、`data.station_file`。配置里写了的值总是优先。`auto` 表示按说明里的规则由数据推出。

“依据”一栏说明默认值的来历：

| 依据 | 含义 |
|---|---|
| 必填 | 只有你知道 |
| 文献 | 数值或方法取自所引文献 |
| 理论 | 由方法本身的关系得出（采样、滤波器宽度、远场相位） |
| 规则 | 按写明的规则由你的数据推出；规则里的系数除注明文献者外是本项目的选择 |
| 经验值 | 本项目的经验，没有查到文献取值。数据需要时可以改，改了要说明 |
| 选择 | 由你定，没有对错 |
| 技术 | 只影响速度、内存或文件，不影响结果 |

流程在哪些地方停下来等你：`clean` 之后（名字有歧义、台站没有数据、坐标对不上、分量方向不明、台阵类型介于面状和线状之间、没有去响应且不知道仪器型号）、`stack` 之后（只在诊断列出候选时）、`reference` 之后（拾取的脊是不是基阶）、`inversion` 之后（起始模型、最深节点、速度上下限）。没有需要判断的事项时不停。
""",
}

TAIL = {
    "en": """
## %d. Defaults that differ for line arrays

A line array (`project.array_type: line`, or `auto` on an elongated array) uses:

%s

## References

%s
""",
    "zh": """
## %d. 线状台阵不同的默认值

线状台阵（`project.array_type: line`，或 `auto` 判断为线状）使用：

%s

## 参考文献

%s
""",
}

REFS = [
    "Bensen, G. D., et al. (2007). Processing seismic ambient noise data to obtain reliable broad-band surface wave dispersion measurements. *GJI*, 169, 1239–1260. doi:10.1111/j.1365-246X.2007.03374.x",
    "Brocher, T. M. (2005). Empirical relations between elastic wavespeeds and density in the Earth's crust. *BSSA*, 95, 2081–2092. doi:10.1785/0120050077",
    "Constable, S. C., Parker, R. L., & Constable, C. G. (1987). Occam's inversion. *Geophysics*, 52, 289–300. doi:10.1190/1.1442303",
    "Fang, H., Yao, H., Zhang, H., Huang, Y.-C., & van der Hilst, R. D. (2015). Direct inversion of surface wave dispersion for three-dimensional shallow crustal structure based on ray tracing. *GJI*, 201, 1251–1263. doi:10.1093/gji/ggv080",
    "Farr, T. G., et al. (2007). The Shuttle Radar Topography Mission. *Rev. Geophys.*, 45, RG2004. doi:10.1029/2005RG000183",
    "Foti, S., et al. (2018). Guidelines for the good practice of surface wave analysis. *Bull. Earthquake Eng.*, 16, 2367–2420. doi:10.1007/s10518-017-0206-7",
    "Herrin, E., & Goforth, T. (1977). Phase-matched filters: application to the study of Rayleigh waves. *BSSA*, 67, 1259–1275.",
    "Jiang, C., & Denolle, M. A. (2020). NoisePy: a new high-performance Python tool for ambient-noise seismology. *SRL*, 91, 1853–1866. doi:10.1785/0220190364",
    "Jing, J., et al. (2025). Ambient noise seismic tomography of Tonglushan skarn-type Cu-Fe-Au deposit in Eastern China. *Ore Geology Reviews*, 184, 106718. doi:10.1016/j.oregeorev.2025.106718",
    "Lévêque, J.-J., Rivera, L., & Wittlinger, G. (1993). On the use of the checker-board test to assess the resolution of tomographic inversions. *GJI*, 115, 313–318.",
    "Levshin, A. L., & Ritzwoller, M. H. (2001). Automated detection, extraction, and measurement of regional surface waves. *Pure Appl. Geophys.*, 158, 1531–1545. doi:10.1007/978-3-0348-8264-4_11",
    "Leys, C., Ley, C., Klein, O., Bernard, P., & Licata, L. (2013). Detecting outliers: do not use standard deviation around the mean, use absolute deviation around the median. *J. Exp. Soc. Psychol.*, 49, 764–766. doi:10.1016/j.jesp.2013.03.013",
    "Lin, F.-C., Moschetti, M. P., & Ritzwoller, M. H. (2008). Surface wave tomography of the western United States from ambient seismic noise: Rayleigh and Love wave phase velocity maps. *GJI*, 173, 281–298. doi:10.1111/j.1365-246X.2008.03720.x",
    "Luo, Y., Yang, Y., Xu, Y., Xu, H., Zhao, K., & Wang, K. (2015). On the limitations of interstation distances in ambient noise tomography. *GJI*, 201, 652–661. doi:10.1093/gji/ggv043",
    "Rawlinson, N., & Spakman, W. (2016). On the use of sensitivity tests in seismic tomography. *GJI*, 205, 1221–1243. doi:10.1093/gji/ggw084",
    "Seats, K. J., Lawrence, J. F., & Prieto, G. A. (2012). Improved ambient noise correlation functions using Welch's method. *GJI*, 188, 513–523. doi:10.1111/j.1365-246X.2011.05263.x",
    "Sens-Schönfelder, C. (2008). Synchronizing seismic networks with ambient noise. *GJI*, 174, 966–970. doi:10.1111/j.1365-246X.2008.03842.x",
    "Shapiro, N. M., & Singh, S. K. (1999). A systematic error in estimating surface-wave group-velocity dispersion curves and a procedure for its correction. *BSSA*, 89, 1138–1142.",
    "Shaw, P. R., & Orcutt, J. A. (1985). Waveform inversion of seismic refraction data and applications to young Pacific crust. *Geophys. J. R. astr. Soc.*, 82, 375–414.",
    "Stehly, L., Campillo, M., & Shapiro, N. M. (2007). Traveltime measurements from noise correlation: stability and detection of instrumental time-shifts. *GJI*, 171, 223–230. doi:10.1111/j.1365-246X.2007.03492.x",
    "Takagi, R., Nishida, K., Aoki, Y., et al. (2015). A single bit matters: coherent noise of seismic data loggers. *SRL*, 86(3), 901–907. doi:10.1785/0220150030",
    "Taylor, J. R. (1997). *An Introduction to Error Analysis* (2nd ed.), chapter 6 (Chauvenet's criterion). University Science Books.",
    "Wang, J., Wu, G., & Chen, X. (2019). Frequency-Bessel transform method for effective imaging of higher-mode Rayleigh dispersion curves from ambient seismic noise data. *JGR Solid Earth*, 124, 3708–3723. doi:10.1029/2018JB016595",
    "Xu, Y., Lebedev, S., Meier, T., Bonadio, R., & Bean, C. J. (2021). Optimized workflows for high-frequency seismic interferometry using dense arrays. *GJI*, 227, 875–897. doi:10.1093/gji/ggab260",
    "Yao, H., van der Hilst, R. D., & de Hoop, M. V. (2006). Surface-wave array tomography in SE Tibet from ambient seismic noise and two-station analysis – I. Phase velocity maps. *GJI*, 166, 732–744. doi:10.1111/j.1365-246X.2006.03028.x",
    "Zhu, S., Li, J., Chen, G., Fang, H., & Yao, H. (2026). pDSurfTomo: a high-performance parallel computing package for direct surface wave tomography. arXiv:2604.11920.",
]


def dump(v):
    t = yaml.safe_dump(v, default_flow_style=True, width=10**6).strip().replace("\n...", "")
    return t if len(t) <= 60 else t[:57] + "..."


def main():
    cfg = config.defaults("surface")
    cfg["project"]["array_type"] = "auto"
    keys = params.flat(cfg)
    missing = [k for k, _ in keys if k not in ZH]
    extra = [k for k in ZH if k not in dict(keys)]
    if missing or extra:
        sys.exit("Chinese texts: missing %s; not a key any more %s" % (missing, extra))
    req = {".".join(k) for k in config.REQUIRED}
    for lang, name in (("en", "parameters.md"), ("zh", "parameters.zh.md")):
        L = [HEAD[lang]]
        for n, (sec, en, zh) in enumerate(SECTIONS, start=2):
            L.append("## %d. %s\n" % (n, en if lang == "en" else zh))
            L.append("| key | default | basis | meaning |" if lang == "en" else "| 键 | 默认值 | 依据 | 含义 |")
            L.append("|---|---|---|---|")
            for k, v in keys:
                if k.split(".")[0] != sec:
                    continue
                kind, text = params.DOC[k]
                kind = "required" if k in req else kind
                L.append(
                    "| `%s` | %s | %s | %s |"
                    % (
                        k.split(".", 1)[1],
                        "—" if k in req else "`%s`" % dump(v).replace("|", "\\|"),
                        kind if lang == "en" else KIND_ZH[kind],
                        (text if lang == "en" else ZH[k]).replace("|", "\\|"),
                    )
                )
            L.append("")
        line = "\n".join("- `%s`: `%s`" % (k, dump(v)) for k, v in params.flat(config.LINE_OVERRIDES))
        refs = "\n".join("- " + r for r in REFS)
        L.append(TAIL[lang] % (len(SECTIONS) + 2, line, refs))
        with open(os.path.join(os.path.dirname(HERE), "docs", name), "w") as fh:
            fh.write("\n".join(L).rstrip() + "\n")
        print("wrote docs/%s" % name)


if __name__ == "__main__":
    main()
