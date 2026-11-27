# 参与贡献

中文 | [English](CONTRIBUTING.md)

## 报告问题

用 bug 模板新建 issue：`lant --version`、`lant check`、运行的命令，以及出错那一步的 `log.txt`。不要附上不能公开的数据；能复现问题的小算例最有用。

## 修改代码

- 一次提交一件事；格式化单独提交。
- 一个模块一个职责；约 500 行时检查是否混了几个职责。
- 每个方法或默认值都要有出处（文献或理论），否则在 docs/parameters 里标为"约定"。不为某一批数据加功能。
- 不静默失败：报错说明哪里不对、该改什么。
- 改了设置或方法：同步更新 `lant/config_template.yaml` 和文档（中英两版）以及 `docs/history`。

提交 pull request 之前：

```bash
pip install -e ".[dev]"
ruff check lant tests
ruff format --check lant tests
python -m pytest tests
```
