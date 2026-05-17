# 本地 miniSEED 数据背景噪声面波层析成像完整流程

本项目实现了从原始 miniSEED 地震数据到三维横波速度模型的完整背景噪声面波层析成像流程。整个流程分为 8 个步骤，依次使用 5 个软件包和若干辅助脚本完成。

## 流程概览

```
原始 miniSEED 数据
  │
  ▼ [Step 1] clean_mseed_for_noisepy.py
清洗后的 miniSEED 数据（NoisePy 格式）
  │
  ▼ [Step 2] 准备台站经纬度文件 station.csv
  │
  ▼ [Step 3] local_cc_workflow.py（NoisePy）
互相关函数 (CCF) + 叠加结果 (Stack)
  │    └─ [可选] stack_visualization.py 可视化叠加结果
  │
  ▼ [Step 4] export_egfanalysis_input.py
EGFAnalysisTimeFreq 输入格式（.dat 波形 + stations.txt + EGFfiles.txt）
  │
  ▼ [Step 5] EGFAnalysisTimeFreq（MATLAB）
频散谱图（C.*.dat / G.*.dat）
  │
  ▼ [Step 6] prepare_disperpicker_input.py
DisperPicker 输入格式（重命名 + 配置更新）
  │
  ▼ [Step 7] DisperPicker（pick_v.py）
自动提取的频散曲线
  │
  ▼ [Step 8] DSurfTomo 反演
    ├── extractnew.py          → surfdata.dat
    ├── GenerateIniMOD.py      → MOD 初始模型
    ├── GenerateDSurfTomoInputFile.py → DSurfTomo.in
    ├── DSurfTomo 程序执行反演
    └── GMT 脚本绘图（深度切片、剖面、射线路径、棋盘测试）
```

## 软件包引用

本项目使用了以下开源软件包，感谢原作者的贡献：

- **NoisePy** — 背景噪声互相关计算
  Jiang, C., & Denolle, M. A. (2020). NoisePy: A New High-Performance Python Tool for Ambient-Noise Seismology. *Seismological Research Letters*, 91(3), 1853–1866. https://doi.org/10.1785/0220190364
  GitHub: https://github.com/noisepy/NoisePy

- **EGFAnalysisTimeFreq** — 面波频散分析（MATLAB）
  Yao, H., van der Hilst, R. D., & de Hoop, M. V. (2006). Surface-wave array tomography in SE Tibet from ambient seismic noise and two-station analysis—I. Phase velocity maps. *Geophysical Journal International*, 166(2), 732–744. https://doi.org/10.1111/j.1365-246X.2006.03028.x
  作者：姚华建（中国科学技术大学）

- **DisperPicker** — 基于深度学习的频散曲线自动提取
  Yang, S., Zhang, H., Gu, N., Gao, J., Xu, J., Jin, J., Li, J., & Yao, H. (2022). Automatically Extracting Surface Wave Group and Phase Velocity Dispersion Curves from Dispersion Spectrograms Using a Convolutional Neural Network. *Seismological Research Letters*. https://doi.org/10.1785/0220210280
  GitHub: https://github.com/ShaoboYang-USTC/DisperPicker

- **DSurfTomo** — 面波直接反演三维横波速度结构
  Fang, H., Yao, H., Zhang, H., Huang, Y. C., & van der Hilst, R. D. (2015). Direct inversion of surface wave dispersion for three-dimensional shallow crustal structure based on ray tracing: methodology and application. *Geophysical Journal International*, 201(3), 1251–1263. https://doi.org/10.1093/gji/ggv069

## 环境准备

### Python 环境

```bash
conda create -n noisepy -y python=3.10 pip
conda activate noisepy
pip install noisepy-seis
pip install obspy numpy pandas pyyaml datetimerange h5py scipy matplotlib
```

### MATLAB

EGFAnalysisTimeFreq 需要 MATLAB 运行环境。

### DisperPicker 环境

```bash
cd DisperPicker
pip install -r requirements.txt
```

### DSurfTomo 编译

```bash
cd DSurfTomo
./configure
make
```

### GMT（可选，用于绘图）

DSurfTomo 的绘图脚本需要 GMT 6+。

---

## Step 1：数据清洗

**目的：** 剔除不满一天的数据，将原始 miniSEED 文件转换为 NoisePy 所需的命名格式。

**脚本：** `scripts/clean_mseed_for_noisepy.py`

**输入：**
- 原始 miniSEED 文件目录，按台站分子目录存放
  ```
  RawData/
  ├── CF.17982/
  │   ├── 230901000000    # 文件名为时间戳 YYMMDDHHmmss
  │   ├── 230902000000
  │   └── ...
  ├── CF.58808/
  └── CN.20496/
  ```

**命令：**

```bash
python scripts/clean_mseed_for_noisepy.py \
    --input-root  TestData/RawData \
    --output-root TestData/CleanData \
    --time-zone beijing \
    --target-channel Z \
    --tolerance-seconds 60
```

**参数说明：**
| 参数 | 说明 |
|------|------|
| `--input-root` | 原始 miniSEED 文件根目录 |
| `--output-root` | 清洗后输出目录 |
| `--time-zone` | 文件名时间戳的时区：`beijing`（默认）或 `utc` |
| `--target-channel` | 目标通道/分量，如 `Z`、`BNZ`、`BHZ` |
| `--tolerance-seconds` | 判定"满一天"的容差（秒），默认 60 |
| `--dry-run` | 仅打印操作不实际写入 |

**输出：**
```
CleanData/
├── CF17982_BNZ00_2023258.ms
├── CF17982_BNZ00_2023259.ms
└── ...
```
文件名格式：`{台网}{台站}_{通道}{位置码}_{年+儒略日}.ms`

---

## Step 2：准备台站信息文件

**目的：** 创建 NoisePy 需要的台站坐标文件。

**格式：** CSV 文件，包含以下列：

```csv
network,station,channel,latitude,longitude,elevation
CF,17982,BNZ,42.323349,118.546867,844
CF,58808,BNZ,42.321987,118.54599,822
CN,20496,BNZ,42.318123,118.543602,786
```

**注意事项：**
- 台站名、台网代码需与 miniSEED 文件头信息一致
- 经纬度使用十进制度数
- 高程单位为米
- 保存为 `TestData/station.csv`

---

## Step 3：互相关计算与叠加

**目的：** 使用 NoisePy 计算台站对之间的背景噪声互相关函数并叠加。

**脚本：** `scripts/local_cc_workflow.py`

### 3a. 编写配置文件

复制并修改示例配置：

```bash
cp configs/local_cc_workflow.example.yaml configs/my_config.yaml
```

**关键配置参数（`configs/local_cc_workflow.example.yaml`）：**

```yaml
# 数据路径
data-dir: ../TestData/CleanData      # Step 1 的输出目录
station-csv: ../TestData/station.csv  # Step 2 的台站文件
ccf-dir: ../TestData/CCFData          # 互相关结果输出目录
stack-dir: ../TestData/StackData      # 叠加结果输出目录

# 处理参数
freqmin: 0.5          # 最低频率 (Hz)
freqmax: 10           # 最高频率 (Hz)
sampling-rate: 50.0   # 采样率 (Hz)
maxlag: 500           # 最大滞后时间 (秒)
inc-hours: 24         # 每段数据长度 (小时)
cc-len: 3600          # 互相关窗口长度 (秒)
step: 1800.0          # 窗口滑动步长 (秒)

# 归一化
time-norm: rma         # 时间域归一化方法
freq-norm: rma         # 频率域归一化方法

# 叠加
stack-method: all      # 叠加方法: linear / pws / robust / all
```

### 3b. 运行互相关和叠加

```bash
# 验证配置（不实际运行）
python scripts/local_cc_workflow.py --config configs/my_config.yaml --validate-only

# 仅运行互相关
python scripts/local_cc_workflow.py --config configs/my_config.yaml --stage cc

# 仅运行叠加
python scripts/local_cc_workflow.py --config configs/my_config.yaml --stage stack

# 一次运行全部
python scripts/local_cc_workflow.py --config configs/my_config.yaml
```

### 3c. 可视化叠加结果（可选）

**脚本：** `scripts/stack_visualization.py`

```bash
python scripts/stack_visualization.py \
    --stack-dir TestData/StackData \
    --freqmin 1 --freqmax 5 \
    --component ZZ \
    --stack-method linear \
    --maxlag 10 \
    --style both
```

也可以使用配置文件（见 `configs/stack_visualization.example.yaml`）：

```bash
python scripts/stack_visualization.py --config configs/my_stack_viz_config.yaml
```

**输出：** 移出（moveout）图，保存在 `--outdir` 指定的目录（默认 `TestData/figures/`）。

---

## Step 4：导出 EGFAnalysisTimeFreq 所需数据

**目的：** 将 NoisePy 叠加结果（HDF5 格式）转换为 EGFAnalysisTimeFreq 需要的文本格式。

**脚本：** `scripts/export_egfanalysis_input.py`

```bash
python scripts/export_egfanalysis_input.py \
    --stack-dir  TestData/StackData \
    --stack-method all \
    --component ZZ \
    --output-root TestData/EGFAnalysisInputData
```

**参数说明：**
| 参数 | 说明 |
|------|------|
| `--stack-dir` | NoisePy 叠加结果目录 |
| `--stack-method` | 叠加方法：`linear` / `pws` / `robust` / `all`（默认） |
| `--component` | 分量，默认 `ZZ` |
| `--station` | 可选，仅导出包含指定台站的台站对 |
| `--output-root` | 输出根目录 |

**输出目录结构：**
```
EGFAnalysisInputData/
├── Allstack_linear/ZZ/
│   ├── CF17982_CF58808.dat    # 波形数据（三列：时间, gab, gba）
│   ├── CF17982_CN20496.dat
│   ├── stations.txt           # 台站坐标
│   └── EGFfiles.txt           # 波形文件列表
├── Allstack_pws/ZZ/
└── Allstack_robust/ZZ/
```

每个 `.dat` 文件的格式：
```
  经度_源台站      纬度_源台站      高程_源台站      (第1行)
  经度_接收台站    纬度_接收台站    高程_接收台站    (第2行)
  时间(s)          gab              gba              (第3行起)
```

---

## Step 5：使用 MATLAB 生成频散谱图

**目的：** 使用 EGFAnalysisTimeFreq 对经验格林函数进行时频分析，生成频散谱图。

**工具：** `EGFAnalysisTimeFreq/EGFAnalysisTimeFreq.m`（MATLAB GUI 程序）

### 操作步骤

1. 打开 MATLAB，运行 `EGFAnalysisTimeFreq.m`
2. 在 GUI 中设置参数（参考 `EGFAnalysisTimeFreq/Feidong_disp/fd_config.png`）
3. **重要参数（需要记住，后续步骤会用到）：**
   - **周期范围**：如 0.2 – 5.0 s
   - **速度范围**：如 0.5 – 4.0 km/s
   - **速度步长**：固定为 0.005 km/s（DisperPicker 要求）
   - **周期步长**：建议 0.1 s
4. 输入数据目录：选择 Step 4 的输出目录（如 `TestData/EGFAnalysisInputData/Allstack_robust/ZZ/`）
5. 选择 **Automatic** 模式，点击 **Start Processing**
6. 程序自动处理所有台站对，输出频散谱图

**输出：**
```
EGFAnalysisOutputData/RobustDisperImage/
├── C.CF17982_CF58808.dat    # 相速度频散谱图
├── C.CF17982_CN20496.dat
├── G.CF17982_CF58808.dat    # 群速度频散谱图
├── G.CF17982_CN20496.dat
└── ...
```

每个 `.dat` 文件是一个矩阵（如 701×49），行方向为速度，列方向为周期。

> **注意：** 详细使用说明请参考 `EGFAnalysisTimeFreq/EGFAnalysisTimeFreq_Manual.pdf`

---

## Step 6：准备 DisperPicker 输入

**目的：** 将 EGFAnalysisTimeFreq 的频散谱图转换为 DisperPicker 所需的格式，并更新 DisperPicker 的配置。

**脚本：** `scripts/prepare_disperpicker_input.py`

```bash
python scripts/prepare_disperpicker_input.py \
    --input-dir  TestData/EGFAnalysisOutputData/RobustDisperImage \
    --picker-dir DisperPicker \
    --start-period 0.2 \
    --end-period 5.0 \
    --min-velocity 0.5 \
    --max-velocity 4.0 \
    --station-file DisperPicker/config/station.txt \
    --clean
```

**参数说明：**
| 参数 | 说明 |
|------|------|
| `--input-dir` | EGFAnalysisTimeFreq 输出的频散谱图目录 |
| `--picker-dir` | DisperPicker 根目录 |
| `--start-period` | 起始周期 (s)，需与 Step 5 一致 |
| `--end-period` | 结束周期 (s)，需与 Step 5 一致 |
| `--min-velocity` | 最小速度 (km/s)，需与 Step 5 一致 |
| `--max-velocity` | 最大速度 (km/s)，需与 Step 5 一致 |
| `--station-file` | 台站坐标文件（可选） |
| `--clean` | 清除 DisperPicker 中已有的数据 |

**该脚本会自动：**
1. 复制 `C.*.dat` → `DisperPicker/data/TestData/phase_image/STA1.STA2.dat`
2. 复制 `G.*.dat` → `DisperPicker/data/TestData/group_image/STA1.STA2.dat`
3. 自动检测矩阵尺寸
4. 更新 `DisperPicker/config/config.py` 中的 `range_T`、`range_V`、`root`

---

## Step 7：DisperPicker 自动提取频散曲线

**目的：** 使用预训练的 CNN 模型自动从频散谱图中提取群速度和相速度频散曲线。

### 7a. 检查配置

确认 `DisperPicker/config/config.py` 已被 Step 6 正确更新：

```python
self.root = '/your/path/to/DisperPicker'
self.range_T = [0.2, 5.0, 49]    # [起始周期, 结束周期, 周期点数]
self.range_V = [0.5, 4.0, 701]   # [起始速度, 结束速度, 速度点数]
```

确认 `DisperPicker/config/station.txt` 包含正确的台站坐标：

```
STA1   经度        纬度         高程
STA2   经度        纬度         高程
```

### 7b. 运行提取

```bash
cd DisperPicker
python pick_v.py
```

### 7c. 质量控制（QC）

```bash
cd DisperPicker/result/process
python process_C.py    # 相速度 QC
python process_G.py    # 群速度 QC
python plot_res.py     # 绘制每条频散曲线对比图
```

**输出：**
```
DisperPicker/result/
├── pick_result/           # 原始提取结果
├── process/
│   └── new/               # QC 后的频散曲线
│       ├── phase_velocity/
│       │   ├── STA1.STA2.txt
│       │   └── ...
│       ├── group_velocity/
│       │   ├── STA1.STA2.txt
│       │   └── ...
│       ├── DPN_C.jpg      # 所有相速度频散曲线汇总图
│       └── DPN_G.jpg      # 所有群速度频散曲线汇总图
```

每条频散曲线文件（如 `STA1.STA2.txt`）的格式：
```
周期(s)  速度(km/s)
0.3      2.15
0.4      2.28
...
```

---

## Step 8：DSurfTomo 面波层析成像反演

**目的：** 将频散曲线数据反演为三维横波速度模型。

**脚本位置：** `DSurfTomo/scripts/`

### 方式一：一体化交互式流程（推荐）

使用 `PrepareDSurfTomoProject.py` 一步完成数据准备、初始模型和配置文件生成：

```bash
cd DSurfTomo/scripts
python PrepareDSurfTomoProject.py
```

该脚本会交互式引导你完成以下步骤：
1. 选择输入类型（DisperPicker 结果或 EGFAnalysisTimeFreq 手动拾取）
2. 设置频散曲线的周期和速度参数
3. 生成 `surfdata.dat`（频散数据格式转换）
4. 生成 `MOD`（初始速度模型）
5. 生成 `DSurfTomo.in`（反演配置文件）
6. 可选：生成 `MOD.true`（棋盘测试真实模型）

### 方式二：分步执行

#### 8a. 转换频散数据 → surfdata.dat

**脚本：** `extractnew.py`

```bash
python extractnew.py \
    --input-dir /path/to/DisperPicker/result/process/new \
    --input-type disperpicker \
    --station-file /path/to/station.txt \
    --output-dir /path/to/project \
    --min-period 0.2 --max-period 5.0 --period-step 0.1 \
    --curve-type both
```

也支持 EGFAnalysisTimeFreq 手动拾取结果：

```bash
python extractnew.py \
    --input-dir /path/to/manual_picks \
    --input-type manual \
    --output-dir /path/to/project \
    --min-period 0.2 --max-period 5.0 --period-step 0.1
```

**输出：** `surfdata.dat` + `surfdata.json`（元数据）

#### 8b. 生成初始模型 → MOD

**脚本：** `GenerateIniMOD.py`

```bash
python GenerateIniMOD.py
```

交互式输入参数：
- 深度网格节点（如 `0 0.2 0.5 1 2 3 5`）
- 水平网格间距或网格点数
- 每层初始速度（默认根据频散数据自动计算）

**输出：** `MOD`（初始模型文件）

#### 8c. 生成配置文件 → DSurfTomo.in

**脚本：** `GenerateDSurfTomoInputFile.py`

```bash
python GenerateDSurfTomoInputFile.py
```

交互式输入参数：
- 反演参数（迭代次数、阻尼因子等）
- 速度约束范围
- 平滑参数

**输出：** `DSurfTomo.in`

#### 8d. 可选：生成棋盘测试模型

**脚本：** `GenerateTrueMOD.py`

```bash
python GenerateTrueMOD.py
```

**输出：** `MOD.true`

#### 8e. 运行 DSurfTomo 反演

```bash
/path/to/DSurfTomo/bin/DSurfTomo
```

程序读取当前目录下的 `DSurfTomo.in`、`surfdata.dat`、`MOD`，执行反演。

**输出文件（部分）：**
- `MOD_iter_*.dat` — 每次迭代的速度模型
- `residual.dat` — 残差
- `raypath_*.dat` — 射线路径

### 8f. 结果可视化（GMT 脚本）

```bash
# 深度切片图
bash plotdepth_subplot.gmt --project-dir /path/to/project

# 剖面图
bash plotcross_subplot.gmt --project-dir /path/to/project \
    --profiles "lon1/lat1/lon2/lat2/label"

# 射线路径图
bash plotray_subplot.gmt --project-dir /path/to/project

# 棋盘测试结果图
bash plotchecker_subplot.gmt --project-dir /path/to/project
```

### 8g. 频散数据统计（可选）

```bash
python plot_dispersion_period_counts.py --surfdata /path/to/surfdata.dat
```

---

## 目录结构说明

```
local_ANT/
├── scripts/                            # 辅助脚本（连接各软件包的胶水代码）
│   ├── clean_mseed_for_noisepy.py      # Step 1: 数据清洗
│   ├── local_cc_workflow.py            # Step 3: NoisePy 互相关工作流
│   ├── stack_visualization.py          # Step 3c: 叠加结果可视化
│   ├── export_egfanalysis_input.py     # Step 4: 导出 EGF 输入
│   └── prepare_disperpicker_input.py   # Step 6: 准备 DisperPicker 输入
├── configs/                            # 配置文件示例
│   ├── local_cc_workflow.example.yaml
│   └── stack_visualization.example.yaml
├── EGFAnalysisTimeFreq/                # MATLAB 频散分析程序
├── DisperPicker/                       # 深度学习频散曲线提取
│   ├── config/                         # 配置文件
│   ├── data/TestData/                  # 输入数据
│   └── result/                         # 输出结果
├── DSurfTomo/                          # 面波层析成像反演
│   ├── bin/                            # 可执行文件
│   ├── scripts/                        # 数据准备和绘图脚本
│   └── src/                            # 源代码
└── TestData/                           # 示例数据
    ├── RawData/                        # 原始 miniSEED
    ├── CleanData/                      # 清洗后数据
    ├── CCFData/                        # 互相关结果
    ├── StackData/                      # 叠加结果
    ├── EGFAnalysisInputData/           # EGF 输入
    ├── EGFAnalysisOutputData/          # 频散谱图
    ├── figures/                        # 可视化图片
    └── station.csv                     # 台站信息
```

## 常见问题

### 周期和速度参数需要一致吗？
是的。Step 5（EGFAnalysisTimeFreq）设置的周期范围、速度范围必须在 Step 6（prepare_disperpicker_input.py）中使用相同的值。DisperPicker 要求速度步长固定为 0.005 km/s。

### 如何选择叠加方法？
- `linear`：线性叠加，最常用
- `pws`：相位加权叠加，信噪比较高
- `robust`：鲁棒叠加，对异常值不敏感
- `all`：同时计算所有方法

建议先用 `all` 全部计算，然后通过 `stack_visualization.py` 对比选择最优结果。

### DisperPicker 的 station.txt 格式是什么？
```
台站名  经度         纬度          高程
FD01   117.4333135  31.81069491  12.939
```

### extractnew.py 支持哪些输入类型？
- `disperpicker`：DisperPicker QC 后的结果（`phase_velocity/` 和 `group_velocity/` 目录下的 `.txt` 文件）
- `manual`：EGFAnalysisTimeFreq 手动拾取的 `CDisp.T.*.dat` 和 `GDisp.*.dat` 文件

---

## 许可证

各软件包遵循其原始许可证。辅助脚本（`scripts/` 和 `DSurfTomo/scripts/`）随项目发布。
