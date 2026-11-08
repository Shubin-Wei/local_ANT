# AGENTS.md — 给 AI 助手的说明

本仓库是密集台阵高频背景噪声面波成像的完整流程（原始数据 → 清洗 → 逐小时互相干 → 筛选 → 陷波 → 参考曲线 → 两台法相速度与质控 → DSurfTomo 三维 Vs → 检验和出图），命令行入口是 `lant`。

## 使用者让你处理数据时

严格按 [docs/ai_workflow.md](docs/ai_workflow.md) 做。要点：

1. 先 `lant check`；缺的程序按提示安装（`bash tools/install_external.sh`）。
2. 用 `lant init` 建配置，所有结果写到使用者指定的输出目录；**不在仓库里写结果，不改原始数据，不改 `lant/` 源码**。
3. 逐步运行（`lant run <步骤>`），每步后 `lant inspect <步骤>`、看列出的图，按 ai_workflow.md 第 4 节的判据检查。
4. 参数按 ai_workflow.md 的规则由数据推出，依据写进 `<输出目录>/records/decisions.md`；第 3 节“必须确认”的事项先问使用者。
5. 不为了多出点或结果好看而放宽质控、反选参数；汇报只给数字和图，注明没有验证或受限的部分。

参数的含义和依据：[docs/parameters.md](docs/parameters.md)；操作细节：[docs/tutorial.md](docs/tutorial.md)。

## 修改本项目代码时

- 结构：`lant/`（`clean` `cc` `stack` `notch` `reference` `disp` `qc` `export` 按步骤一个模块，`pipeline.py` 串起来，`inversion/`、`plotting/`），`config/`（配置模板），`docs/`，`tools/install_external.sh`（第三方程序），`examples/synthetic/`（合成算例），`tests/`。
- 默认参数只在 `lant/config.py`；改了默认值或方法要同步 `docs/parameters.md`、`docs/tutorial.md`。
- 代码、配置、图中文字用英文；文档用中文。
- 改动后运行 `python -m pytest tests`，并用合成算例（`examples/synthetic`）跑通受影响的步骤。
