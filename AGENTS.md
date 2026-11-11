# AGENTS.md — instructions for AI assistants

[中文](AGENTS.zh.md) | English

This repository is a complete workflow for high-frequency ambient-noise surface-wave tomography with dense arrays (raw records → cleaning → hourly cross-coherence → stack and station diagnosis → reference curve → two-station phase velocities with quality control → 3-D Vs with DSurfTomo → checks and figures). The command line is `lant`.

## When the user asks you to process data

Follow [docs/ai_workflow.md](docs/ai_workflow.md). In short:

1. `lant check` first; install what it lists as missing (`bash tools/install_external.sh`).
2. Start with `lant run all --raw DIR --stations FILE --out DIR` (or `lant init -o CONFIG` and `lant run all -c CONFIG`). All results go to the user's output directory. **Do not write results into the repository, do not change the raw data, do not edit `lant/`.**
3. After every stop run `lant status -c CONFIG --json`: it gives the state of each step, the reason, and the next command. Read `OUT/.lant/steps/<step>.md` and the files it lists before going on.
4. Where the chain stops for a judgement, show the user the items under "FOR YOU TO JUDGE" with the figures, and **let the user decide: never run `lant confirm` on your own.** Record the user's decision and reason in `OUT/.lant/decisions.md`.
5. Do not loosen quality control or tune parameters to get more points or a nicer picture. Report numbers and figures, and say what was not checked or is limited.

Meaning and basis of every parameter: [docs/parameters.md](docs/parameters.md). Details of every step: [docs/tutorial.md](docs/tutorial.md); principles and sources: [docs/methods.md](docs/methods.md); output files: [docs/outputs.md](docs/outputs.md); known problems and remedies: [docs/troubleshooting.md](docs/troubleshooting.md).

## When you change the code of this project

- Layout: `lant/` (one module per step: `clean` `cc` `stack` `reference` `disp` `qc` `export`; `pipeline.py` runs them; `manifest.py` keeps the record of what stands; `auto.py` the rules; `confirm.py` the judgements; `report.py` the step reports; `params.py` the meaning and basis of every key; `response.py`, `stations.py` (station tables), `station_match.py` (records to stations), `pulse.py` (1-s timing pulse), `fj.py` (F-J spectrum and ridge), `inversion/`, `plotting/`), `config/example.yaml` (written by `lant init`), `docs/`, `tools/install_external.sh` (third-party programs), `examples/synthetic/` (data with a known answer), `tests/` (one file per module, helpers in `tests/util.py`).
- Defaults live only in `lant/config.py`; every key needs an entry in `lant/params.py` (a test checks it). After changing a default or a method, update `docs/parameters.md` (`python tools/make_parameter_doc.py`), `docs/tutorial.md` and `docs/methods.md` in both languages, and add the change to `docs/history.md`.
- Code, configuration and text in figures are English. Every document exists in English (`x.md`) and Chinese (`x.zh.md`); keep the two in step, section by section.
- One clear responsibility per module, with a name that says it. About 500 lines is the point to check whether a module holds several (split it when it does; tables of defaults are exempt); below 200 lines do not split just to have small files.
- Run `python -m pytest tests`, and the synthetic example for the steps you touched.
- DSurfTomo and pDSurfTomo are not modified (only the memory-reserve patch for pDSurfTomo in the install script).
- A method or a default needs a source: literature, theory, or an explicit label "empirical". Do not invent thresholds.
