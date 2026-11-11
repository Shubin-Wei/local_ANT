# 参与贡献

中文 | [English](CONTRIBUTING.md)

## 报告问题

用 bug 模板新建 issue。模板要的是不需要你的数据也能复现问题的信息：

- `lant --version` 和 `lant check` 的输出；
- 运行的命令和它打印的内容；
- `lant status -c CONFIG --json`，以及失败步骤的报告（`OUT/.lant/steps/<步骤>.md`）；
- `OUT/.lant/lant.log` 的末尾。

不要附上不能公开的原始记录或台站坐标；能复现问题的合成算例（`lant demo`、`python -m lant.synthetic`）最有用。

## 修改代码

遵循[开发规范](docs/engineering.zh.md)和 [AGENTS.zh.md](AGENTS.zh.md) 中关于代码的部分：一次提交一件事；默认值只放在 `lant/config.py`，依据写在 `lant/params.py`；每份文档都有中英文；不修改 DSurfTomo 和 pDSurfTomo；方法和门限都要有来源。

提交 pull request 之前：

```bash
pip install -e ".[dev]"
ruff check lant tests tools examples
ruff format --check lant tests tools examples
python -m pytest tests
python tools/run_synthetic.py --dir DIR      # 改了某一步时，先 python -m lant.synthetic --quick --out DIR
```

改动计算时，要对确定性产物（叠加、测量、求解器输入）与上一版本做同输入比较，并在 `docs/history.zh.md` 中记录。
