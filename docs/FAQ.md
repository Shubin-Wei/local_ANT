# 常见问题(FAQ)

后续将在此处持续补充项目安装、配置、运行和结果分析过程中遇到的常见问题。

安装类问题(Conda、TensorFlow 1.13、gfortran 编译等)另见 [INSTALL.md](INSTALL.md) 末尾的"常见问题"一节。

---

## 流程与配置

### 周期和速度参数需要各步骤一致吗?

基本一致即可。06b 的 `prepare.period_step` 不必与 05 步的 `DeltaT` 完全相同,只要求 period_step 生成的周期网格点能落在 05 步输出曲线的周期采样中(通常取 `DeltaT` 或其整数倍最稳妥)。DisperPicker 要求速度步长 `DeltaV` 固定为 0.005 km/s。

### 如何选择叠加方法?

- `linear`:线性叠加,最常用;
- `pws`:相位加权叠加,信噪比较高;
- `robust`:鲁棒叠加,对异常值不敏感;
- `all`:同时计算所有方法。

建议 02 步用 `all` 全部计算,再用 03 步(`03_stack_qc.py` / `03b_stack_panels.py`)对比选优。

### 第 05 步之后选 06a 还是 06b?

看台阵形态:台站近似排成一条测线 → **06a**(evodcinv 一维反演拼伪二维剖面);台站覆盖二维区域 → **06b**(DSurfTomo 三维反演)。两分支相互独立,详见 [workflow.md](workflow.md) 的对比表。

### 06b 的数据转换脚本吃什么输入?

`extractnew.py`(以及封装它的 `01_prepare.py`)仅支持 `egfpy` 输入:05 步 `qc_result/` 下的 `CDisp.*.txt`(相速度)与 `GDisp.*.txt`(群速度)。文件头两行为台站经纬度,数据列为 `周期 速度 信噪比 概率`;第 4 列是 DisperPicker 概率,不作为保留标志。台站坐标已内嵌文件头,无需额外 station 文件。

### 输出目录已存在,程序拒绝运行?

这是防静默覆盖的保护:01 步对已存在的输出文件、05 步对非空的 `qc_result/`、06b-1 对非空的工程目录都会报错停止。确认要覆盖时加 `--overwrite`。

## 环境

### 运行 05 步提示 TensorFlow 不可导入?

05 步必须在 `egfanalysispy` 环境运行:`conda activate egfanalysispy`。入口脚本会在调用前检查并给出该提示。

### 运行 06a 提示 `No module named 'evodcinv'`?

evodcinv 环境未创建或未激活。它使用独立的 Conda 环境:`conda env create -f environments/evodcinv.yml && conda activate evodcinv`,见 [INSTALL.md](INSTALL.md) §2e。

### 怎么快速自查环境?

```bash
python tools/check_environment.py
```

在哪个环境里运行就检查哪个环境;与当前环境无关的项报 FAIL 属正常(如在 noisepy 环境里 TensorFlow 报 FAIL)。

### 06b 提示 DSurfTomo executable not found?

DSurfTomo 尚未编译,或 `configs/00_project.yaml` 的 `dsurftomo_bin` 路径不对。编译方法见 [INSTALL.md](INSTALL.md) §4;编译产物默认在 `DSurfTomo/bin/DSurfTomo`。

