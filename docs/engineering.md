# Development rules

[中文](engineering.zh.md) | English

1. Make the smallest sufficient change; leave unrelated code alone.
2. Keep one purpose per commit; separate features, fixes and refactoring.
3. Give each module a clear responsibility and each function one explicit task.
4. Reduce deep nesting and duplication. At about 500 lines, review responsibilities; line count is a reminder, not a limit.
5. Use descriptive names, the existing style and automated formatting.
6. Explain reasons, constraints and unusual handling in comments; keep them current.
7. Preserve error context and tracebacks; never silently swallow failures or log secrets.
8. Validate external input, parameterize queries where applicable and keep credentials out of source control.
9. Run relevant checks; add meaningful regression tests for defects.
10. Review the diff before committing; remove debugging, temporary artifacts and unrelated changes; explain the purpose.
11. Update documentation when interfaces, configuration or startup changes.
12. Solve current needs with proportionate design; avoid speculative abstractions and dependencies.

## Local checks

Install `pip install -e ".[dev]"`, then run `ruff check lant tests tools examples`, `ruff format --check lant tests tools examples` and `python -m pytest tests --cov=lant`; the lint rules and the reasons for the ignored ones are in `pyproject.toml`. Formatting stays in its own commits. Build both distributions with `python -m build`, then verify a wheel outside the checkout. CI runs lint and format, the tests on Python 3.11, 3.12 and 3.13, the quick synthetic example to the inversion project, and an installed wheel.

## Responsibility review, 2026-11-13

The 500-line rule triggered these responsibility splits, kept separate from numerical/performance changes:

| Original module | Current responsibilities and lines |
|---|---|
| `data.py`, 514 | `data.py`, 400: waveform/day access; `record_index.py`: raw discovery and header reuse |
| `report.py`, 823 | `report.py`, 388: parameter provenance and assembly; `step_descriptions.py`, 454: recorded step narratives |
| `clean.py`, 651 | `clean.py`, 471: step orchestration; `clean_archive.py`, 195: one station-day archive and its reuse key |
| `inversion/run.py`, 558 | `run.py`, 286: solver scheduling; `model_evaluation.py`, 292: model I/O, exact predictions and recovery statistics |

Four modules remain above 500 lines after review. `config.py` (716) contains the central default tables and configuration lifecycle; `params.py` (543) is the parameter registry (table exemptions). `cli.py` (628) remains the command adapter/dispatch entry point; `stations.py` reads and normalises station tables, and the matching of record names to the table moved to `station_match.py` (records to stations); `disp.py` (575) keeps one measurement step and its shared worker context. The station diagnosis of the stack step is its own module, `station_qc.py`, so `stack.py` keeps stacking and removal only. Retaining these existing boundaries is deliberate; the line count is a review trigger, not evidence that every file has to fit a hard cap. Reconsider a boundary when adding a distinct responsibility, rather than moving a few lines solely to cross 500.
