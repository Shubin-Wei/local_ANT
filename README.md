# local_ANT

[![tests](https://github.com/Shubin-Wei/local_ANT/actions/workflows/tests.yml/badge.svg)](https://github.com/Shubin-Wei/local_ANT/actions/workflows/tests.yml)
[![Python](https://img.shields.io/badge/Python-3.8%2B-blue)](pyproject.toml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green)](LICENSE)
[![Development](https://img.shields.io/badge/status-development-orange)](docs/history.md)

[中文](README.zh.md) | English

High-frequency ambient-noise surface-wave tomography for dense (nodal) arrays, from raw records to a 3-D shear-velocity model. One command line, `lant`; one configuration file; every step writes a report that says what it did, with which values, and what is left for you to judge.

![Synthetic dispersion and shear-velocity model](docs/assets/synthetic-result.png)

*Known-truth synthetic example: dispersion on the left, horizontally averaged Vs on the right. This illustrates a laterally uniform medium; it is not a demonstration of general 3-D resolution.*

Start with `OUT/REPORT.md`. Retained scientific products are under `results/`; processing files are under `work/`; state and rotating logs are under `.lant/`. Manual runs show progress on stderr, and `lant run ... -v` shows details. [Output guide](docs/outputs.md).

```
raw records ─ clean ─ cc ─ stack ─ reference ─ disp ─ inversion ─ run ─ evaluate ─ figures
```

| Step | What it does | Method |
|---|---|---|
| `clean` | reads any directory layout and any format ObsPy reads, matches the station table to the records, cuts processing days, removes instrument responses when response files are given, tests for a 1-s timing pulse | Bensen et al. (2007); Takagi et al. (2015) |
| `cc` | hourly cross-coherence of the vertical components | Xu et al. (2021) |
| `stack` | stack and sub-stacks; station diagnosis (leave-one-out correlation and time shift) | Stehly et al. (2007); Sens-Schönfelder (2008) |
| `reference` | reference phase-velocity curve of the array from its frequency-Bessel spectrum | Wang et al. (2019) |
| `disp` | two-station phase velocities with a phase-matched filter, jackknife errors, quality control | Yao et al. (2006); Levshin & Ritzwoller (2001); Bensen et al. (2007) |
| `inversion`, `run` | 1-D starting model, direct inversion of the travel times for 3-D Vs with DSurfTomo (or pDSurfTomo) | Fang et al. (2015) |
| `evaluate`, `figures` | residuals, checkerboard recovery, control inversions; maps, sections and the other figures | Lévêque et al. (1993); Rawlinson & Spakman (2016) |

## What you give, what the program decides

You give three paths: the raw data, the station table, the output directory. Everything else has a default or is derived from your data by a written rule, and the step reports show each derivation with your numbers. The program stops only where a scientific judgement is yours (for example: is the picked ridge the fundamental mode?) and tells you what to look at and which command continues.

## Scope and limits

- Vertical components only (ZZ, fundamental-mode Rayleigh waves). Areal arrays and line arrays (one project per line).
- The noise field is assumed to be close to isotropic; the bias from directional noise is not tested.
- Linux. Station coordinates in longitude and latitude (no projected coordinates). Arrays across the date line or beyond 80° latitude are refused.
- DSurfTomo accepts at most 60 periods (pDSurfTomo 80); station pairs closer than about 200 m carry distance errors from its single-precision coordinates.
- Without response files every station must be the same instrument model. Sensors of one model still differ within their tolerance, which leaves random phase errors of the order of 0.1–1 % of the travel time.
- Numbers from the synthetic tests are in [docs/validation.md](docs/validation.md).

## Install

```bash
git clone https://github.com/Shubin-Wei/local_ANT.git && cd local_ANT
conda env create -f environment.yml && conda activate lant     # Python, ObsPy, GMT 6, gfortran
pip install -e .
lant install-external                                  # DSurfTomo and pDSurfTomo, fixed versions
lant check                                                      # says what is missing and how to install it
```

## Quick start

```bash
lant run all --raw /data/raw --stations /data/stations.csv --out /data/result
```

writes `/data/result/config.yaml` and runs the chain. When it stops it prints what to look at and the next command, for example

```bash
lant confirm reference -c /data/result/config.yaml --reason "DECISION_AND_EVIDENCE"
lant run all -c /data/result/config.yaml          # continues from the first step that does not stand
```

`lant status -c CONFIG` shows the state of every step and the next command; `lant plan -c CONFIG` shows what every step will do before it runs. To try the chain on data with a known answer first: [docs/quickstart.md](docs/quickstart.md) (about half an hour, installation included).

## Documents

| | |
|---|---|
| [docs/quickstart.md](docs/quickstart.md) | quick start: installation, a complete run on synthetic data, then your own data |
| [docs/tutorial.md](docs/tutorial.md) | tutorial: data preparation, every step: purpose, procedure, parameters, what to look at |
| [docs/methods.md](docs/methods.md) | methods: principles, formulas, sources and conditions of validity of every step |
| [docs/outputs.md](docs/outputs.md) | what every output file holds and when to look at it |
| [docs/troubleshooting.md](docs/troubleshooting.md) | common problems and remedies, step by step |
| [docs/parameters.md](docs/parameters.md) | every configuration key: default, meaning, what the default rests on |
| [docs/ai_workflow.md](docs/ai_workflow.md) | running the chain with an AI assistant |
| [docs/validation.md](docs/validation.md) | synthetic tests and their numbers |
| [docs/history.md](docs/history.md) | changes between versions |
| [AGENTS.md](AGENTS.md) | instructions for AI assistants working in this repository |

## Citation

If you use this workflow, please cite the methods it rests on (table above; full references in [docs/parameters.md](docs/parameters.md)) and the third-party programs DSurfTomo (Fang et al. 2015, GJI 201:1251) and, when used, pDSurfTomo (Zhu et al. 2026, arXiv:2604.11920), as well as ObsPy and GMT.

## License

MIT (see `LICENSE`). The third-party programs are downloaded by `tools/install_external.sh` under their own licences (DSurfTomo MIT, pDSurfTomo GPL-3) and are not part of this repository.
