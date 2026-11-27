# AGENTS.md — 给 AI 助手的说明

中文 | [English](AGENTS.md)

本仓库是 `lant`，密集台阵背景噪声面波成像的命令行工具：`lant cc`（互相关）、`lant disp`（参考曲线与相速度）、`lant inv`（三维反演，Fang et al. 2015）、`lant checker`（分辨率检验）。原则：通用、易用、透明、简洁、正确。

## 使用者让你处理数据时

1. 先 `lant check`，缺什么装什么。
2. `lant init -o 配置 --raw 目录 --stations 台站表 --out 目录 --velocity 最慢 最快`；把建议值及其规则（写在配置里）给使用者看。
3. **一步一步**运行（先 `lant cc`，再 `lant disp`，再 `lant inv`，再 `lant checker`）。每步之后读它的 `log.txt` 和 `summary.txt`，把 [docs/tutorial.zh.md](docs/tutorial.zh.md) 列出的图给使用者看，再往下走。最重要的是由使用者判断 `2_disp/figures/fj_reference.png`（参考曲线是不是基阶）。
4. 结果只写到使用者的输出文件夹。不改原始数据，不改 `lant/`。
5. 不为了多出数据或图好看而放宽质控、调参数。汇报给数字和图，并说明哪些没有检验。使用者没要求时不设 `cores`。

参数及依据：[docs/parameters.zh.md](docs/parameters.zh.md)；方法：[docs/methods.zh.md](docs/methods.zh.md)。

## 修改代码时

- 结构：`lant/config.py` + `config_template.yaml`（默认值的唯一位置，每个键都有注释）、`step.py`（每步都写的文件）、`records.py`（索引、对应、台站日）、`correlate.py` + `cc_kernel.py`（第 1 步）、`fj.py`、`reference.py`、`dispersion.py`、`disp_run.py`（第 2 步）、`inversion/`（`run.py` 第 3 步、`checker.py`、`start_model.py`、`project.py`；`dsurf_*.py` 为 DSurfTomo 移植、`model_regularization.py`、`least_squares.py`）、`plotting/`、`init.py`、`export.py`、`cli.py`；`tests/` 单元测试。
- 每个方法或默认值都要有出处（文献或理论），否则在 docs/parameters 里标为"约定"。不为某一批数据加功能。
- 改了设置或方法：同步更新 `config_template.yaml`、`docs/parameters`、`docs/methods`、`docs/tutorial` 和 `docs/history` 的中英两版。
- `dsurf_*` 模块逐句按 DSurfTomo，文档字符串里标 "lant" 的地方除外。
- 提交前：`ruff check lant tests`、`ruff format --check lant tests`、`python -m pytest tests`。
