# AGENTS.md — instructions for AI assistants

[中文](AGENTS.zh.md) | English

This repository is a complete workflow for high-frequency ambient-noise surface-wave tomography with dense arrays (raw records → cleaning → hourly cross-coherence → stack and station diagnosis → reference curve → two-station phase velocities with quality control → 3-D Vs by direct inversion of the travel times (Fang et al. 2015, DSurfTomo's method; lant's regularization) → checks and figures). The command line is `lant`.

## When the user asks you to process data

Follow [docs/ai_workflow.md](docs/ai_workflow.md). In short:

1. `lant check` first; install what it lists as missing.
2. Start with `lant run all --raw DIR --stations FILE --out DIR` (or `lant init -o CONFIG` and `lant run all -c CONFIG`). All results go to the user's output directory. **Do not write results into the repository, do not change the raw data, do not edit `lant/`.**
3. After every stop run `lant status -c CONFIG --json`: it gives the state of each step, the reason, and the next command. Read `OUT/.lant/steps/<step>.md` and the files it lists before going on.
4. Where the chain stops for a judgement, show the user the items under "FOR YOU TO JUDGE" with the figures, and **let the user decide: never run `lant confirm` on your own.** Record the user's decision and reason in `OUT/notes.md`.
5. Do not loosen quality control or tune parameters to get more points or a nicer picture. Report numbers and figures, and say what was not checked or is limited.
6. Create files only as notes (`lant note`) and as temporary files under `OUT/scratch/`; `lant status --json` lists any other file lant did not write.

Meaning and basis of every parameter: [docs/parameters.md](docs/parameters.md). Details of every step: [docs/tutorial.md](docs/tutorial.md); principles and sources: [docs/methods.md](docs/methods.md); output files: [docs/outputs.md](docs/outputs.md); known problems and remedies: [docs/troubleshooting.md](docs/troubleshooting.md).

## When you change the code of this project

- Layout: `lant/` (one module per step: `clean` (station-day archive worker: `clean_archive.py`) `cc` `stack` `reference` `disp` `qc` `export`; `pipeline.py` runs them; `manifest.py` keeps the record of what stands; `lineage.py` freezes execution parameters and indexes result origins; `auto.py` the rules; `confirm.py` the judgements; `report.py` assembles reports; `step_descriptions.py` describes recorded step results; `params.py` the meaning and basis of every key; `response.py`, `stations.py` (station tables), `station_match.py` (records to stations), `station_qc.py` (station diagnosis of the stack step), `pulse.py` (1-s timing pulse), `cc_kernel.py` (compiled pair sums of cc), `fj.py` (F-J spectrum and ridge), `inversion/` (`run.py` schedules solvers; `model_evaluation.py` reads and evaluates models), `plotting/` (including `figure_batch.py` for isolated main-figure workers); project infrastructure: `record_index.py` (raw discovery and cached headers), `paths.py` (layout), `records.py` (step context, logs), `artifacts.py` (file register, write lock context, cache cleanup), `delivery.py` (results, project README), `run_report.py` (REPORT.md), `assistant_rules.py` (the rules file for AI assistants in a project), `progress.py`, `processes.py` (solver groups, pool watchdog), `resources.py` (memory and disk)), `config/example.yaml` (the full template, as written by `lant init --full`), `docs/`, `tests/` (one file per module, helpers in `tests/util.py`).
- Defaults live only in `lant/config.py`; every key needs an entry in `lant/params.py` (a test checks it). The same generator updates `config/example.yaml`; `python tools/make_parameter_doc.py --check` verifies all three generated references without writing. After changing a default or a method, update `docs/parameters.md` (`python tools/make_parameter_doc.py`), `docs/tutorial.md` and `docs/methods.md` in both languages, and add the change to `docs/history.md`.
- Code, configuration and text in figures are English. Every document exists in English (`x.md`) and Chinese (`x.zh.md`); keep the two in step, section by section.
- One clear responsibility per module, with a name that says it. About 500 lines is the point to check whether a module holds several (split it when it does; tables of defaults are exempt); below 200 lines do not split just to have small files.
- Run `python -m pytest tests` and the smoke test of the chain (`tests/smoke_records.py`, as in the CI); a change of a method also needs the inversion-level synthetic comparison of docs/validation.md.
- The `dsurf_*` modules follow DSurfTomo statement by statement, except the changes lant made on purpose since 2.0.1 (regularization, rejection, least-squares solve, depth-kernel derivative, synthetic noise; listed in docs/python_backends.md "2.0.1 changes" and marked "lant" in the module docstrings). Change their numbers only together with a new comparison against the original, and add any further departure to that list with its basis.
- The program chooses the cores (`resources.cores: auto`); do not set them unless the user asks, and report how many each step used.
- A method or a default needs a source: literature, theory, or an explicit label "empirical". Do not invent thresholds.

## Data processing and repository development

The data-processing rules do not prevent development of the repository when the user explicitly asks for it. A recorded reason, an interactive terminal or `--yes` do not show that the user decided; only the user's decision makes a confirmation valid. Development follows the [development rules](docs/engineering.md).
