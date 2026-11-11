# Contributing

[中文](CONTRIBUTING.zh.md) | English

## Reporting a problem

Open an issue with the bug template. It asks for what makes a problem reproducible without your data:

- `lant --version` and the output of `lant check`;
- the command, and what it printed;
- `lant status -c CONFIG --json`, and the report of the failed step (`OUT/.lant/steps/<step>.md`);
- the end of `OUT/.lant/lant.log`.

Do not attach raw records or station coordinates you may not share; a synthetic example (`lant demo`, `python -m lant.synthetic`) that shows the problem is the most useful.

## Changing the code

Follow the [development rules](docs/engineering.md) and the code section of [AGENTS.md](AGENTS.md): one purpose per commit, defaults only in `lant/config.py` with their basis in `lant/params.py`, every document in English and Chinese, DSurfTomo and pDSurfTomo unmodified, and no method or threshold without a source.

Before a pull request:

```bash
pip install -e ".[dev]"
ruff check lant tests tools examples
ruff format --check lant tests tools examples
python -m pytest tests
python tools/run_synthetic.py --dir DIR      # after python -m lant.synthetic --quick --out DIR, for a changed step
```

A change of a computation needs the same-input comparison of the deterministic products (stack, measurements, solver inputs) against the previous version, and a note in `docs/history.md`.
