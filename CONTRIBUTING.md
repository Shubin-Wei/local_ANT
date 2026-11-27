# Contributing

[中文](CONTRIBUTING.zh.md) | English

## Reporting a problem

Open an issue with the bug template: `lant --version`, `lant check`, the command, and the `log.txt` of the step that failed. Do not attach data you may not share; a small example that shows the problem is the most useful.

## Changing the code

- One purpose per commit; formatting in its own commit.
- One responsibility per module; about 500 lines is the point to check whether a module holds several.
- Every method or default needs a source (literature or theory), or is labelled "convention" in docs/parameters.md. Nothing for one particular data set.
- No silent failures: errors say what is wrong and what to change.
- A changed setting or method: update `lant/config_template.yaml` and the documents (English and Chinese), and `docs/history.md`.

Before a pull request:

```bash
pip install -e ".[dev]"
ruff check lant tests
ruff format --check lant tests
python -m pytest tests
```
