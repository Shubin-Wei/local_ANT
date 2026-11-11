# AGENTS.md — 给 AI 助手的说明

中文 | [English](AGENTS.md)

本仓库是密集台阵高频背景噪声面波成像的完整流程（原始记录 → 清洗 → 逐小时互相干 → 叠加与台站诊断 → 参考曲线 → 两台法相速度与质量控制 → DSurfTomo 三维 Vs → 检验和出图），命令行入口是 `lant`。

## 使用者让你处理数据时

按 [docs/ai_workflow.zh.md](docs/ai_workflow.zh.md) 做。要点：

1. 先 `lant check`；缺的东西按它列出的命令安装（`lant install-external`）。
2. 用 `lant run all --raw 目录 --stations 台站表 --out 目录` 开始（或 `lant init -o 配置` 后 `lant run all -c 配置`）。所有结果写到使用者的输出目录。**不在仓库里写结果，不改原始数据，不改 `lant/`。**
3. 每次停下后运行 `lant status -c 配置 --json`：它给出每一步的状态、原因和下一条命令。先读 `输出目录/.lant/steps/<步骤>.md` 和它列出的文件，再往下走。
4. 流程停下来等判断时，把报告里 “FOR YOU TO JUDGE” 下的事项连同图交给使用者，**由使用者决定：不得自行运行 `lant confirm`。** 把使用者的决定和理由记到 `输出目录/notes.md`。
5. 不为了多出点或图好看而放宽质量控制、调参数。汇报只给数字和图，并说明哪些没有检验、哪些受限。
6. 只以两种方式建文件：笔记（`lant note`）和 `OUT/scratch/` 下的临时文件；`lant status --json` 会列出其他所有不是 lant 写的文件。

每个参数的含义和依据：[docs/parameters.zh.md](docs/parameters.zh.md)。每一步的细节：[docs/tutorial.zh.md](docs/tutorial.zh.md)；原理和依据：[docs/methods.zh.md](docs/methods.zh.md)；输出文件：[docs/outputs.zh.md](docs/outputs.zh.md)；已知问题与处理：[docs/troubleshooting.zh.md](docs/troubleshooting.zh.md)。

## 修改本项目代码时

- 结构：`lant/`（每步一个模块：`clean` `cc` `stack` `reference` `disp` `qc` `export`；`pipeline.py` 串起来；`manifest.py` 记录哪些结果成立；`auto.py` 是规则；`confirm.py` 是判断事项；`report.py` 是步骤报告；`params.py` 是每个键的含义和依据；另有 `response.py`、`stations.py`（台站表）、`station_match.py`（记录对应台站）、`pulse.py`（1 s 授时脉冲）、`fj.py`（F-J 谱和脊线）、`inversion/`、`plotting/`；工程基础设施：`paths.py`（目录布局）、`records.py`（步骤上下文、日志）、`artifacts.py`（文件登记、写锁上下文、缓存清理）、`delivery.py`（成果、REPORT、工程规则）、`progress.py`、`processes.py`（求解进程组、进程池看门狗）、`resources.py`（内存和磁盘）、`synthetic.py` 与 `demo.py`（已知答案的数据），`_external/`（第三方程序的编译脚本、补丁、fwd1d 源码）），`config/example.yaml`（由 `lant init` 写出），`docs/`，`examples/synthetic/`（与已知答案的比较），`tests/`（每个模块一个文件，公用函数在 `tests/util.py`）。
- 默认值只放在 `lant/config.py`；每个键都要在 `lant/params.py` 里有说明（有测试检查）。改了默认值或方法，要同步改两种语言的 `docs/parameters`（`python tools/make_parameter_doc.py`）、`docs/tutorial` 和 `docs/methods`，并在 `docs/history` 里记下改动。
- 代码、配置、图中文字用英文。每份文档都有英文版（`x.md`）和中文版（`x.zh.md`），两版逐节对应。
- 每个模块一个明确的职责，名字说明它做什么。约 500 行时检查是否混入了几个职责（混入了就拆；默认值表不算）；200 行以内不为凑小文件而拆。
- 运行 `python -m pytest tests`，并用合成算例跑通改动涉及的步骤。
- 不改 DSurfTomo 和 pDSurfTomo（只有安装脚本里 pDSurfTomo 的内存预留补丁）。
- 方法和默认值要有出处：文献、理论，或明确标为“经验值”。不要自己造门限。

## 操作与开发边界

处理数据时，先读 `REPORT.md`、`lant status --json` 及所指证据；程序产物只通过 lant 生成，笔记用 `lant note "内容" -c CONFIG` 追加。用户授权的配置修改需记录理由。科学确认必须来自用户决定，并通过 `lant confirm STEP -c CONFIG --reason "用户决定及依据"` 留痕；终端类型不证明身份。未登记文件只提示，不删除。清理先用 `lant clean-work -c CONFIG` 预览，实际清理需要 `--execute --yes`。处理数据规则不禁止用户明确授权的仓库开发。

开发遵循[十二条规范](docs/engineering.zh.md)。
