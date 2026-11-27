# AGENTS.md — instructions for AI assistants

[中文](AGENTS.zh.md) | English

This repository is `lant`, a command-line tool for ambient-noise surface-wave tomography with dense arrays: `lant cc` (cross-correlation), `lant disp` (reference curve and phase velocities), `lant inv` (3-D inversion, Fang et al. 2015), `lant checker` (resolution test). Principles: general, easy to use, transparent, simple, correct.

## When the user asks you to process data

1. `lant check`; install what is missing.
2. `lant init -o CONFIG --raw DIR --stations FILE --out DIR --velocity SLOW FAST`; show the user the suggested values and their rules (written in CONFIG).
3. Run the steps **one at a time** (`lant cc`, then `lant disp`, then `lant inv`, then `lant checker`). After each, read its `log.txt` and `summary.txt` and show the user the figures listed in [docs/tutorial.md](docs/tutorial.md) before going on. Above all, the user must judge `2_disp/figures/fj_reference.png` (is the reference the fundamental mode?).
4. Write results only to the user's output folder. Do not change the raw data or `lant/`.
5. Do not loosen quality control or tune parameters to get more data or a nicer picture. Report numbers and figures and say what was not checked. Do not set `cores` unless the user asks.

Parameters and their basis: [docs/parameters.md](docs/parameters.md); methods: [docs/methods.md](docs/methods.md).

## When you change the code

- Layout: `lant/config.py` + `config_template.yaml` (the only place of defaults; every key commented), `step.py` (what every step writes), `records.py` (index, matching, station-day), `correlate.py` + `cc_kernel.py` (step 1), `fj.py`, `reference.py`, `dispersion.py`, `disp_run.py` (step 2), `inversion/` (`run.py` step 3, `checker.py`, `start_model.py`, `project.py`; `dsurf_*.py` the DSurfTomo port, `model_regularization.py`, `least_squares.py`), `plotting/`, `init.py`, `export.py`, `cli.py`; `tests/` unit tests.
- Every method or default needs a source (literature or theory) or is labelled "convention" in docs/parameters.md. No function for one particular data set.
- A changed setting or method: update `config_template.yaml`, `docs/parameters.md`, `docs/methods.md`, `docs/tutorial.md` and `docs/history.md` in English and Chinese (`x.zh.md`).
- The `dsurf_*` modules follow DSurfTomo statement by statement except where marked "lant" in their docstrings.
- Before committing: `ruff check lant tests`, `ruff format --check lant tests`, `python -m pytest tests`.
