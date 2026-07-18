# Configuration Guide

## Layers and precedence

Configuration has three layers, highest priority first:

```
CLI args  >  stage config (configs/0X_*.yaml)  >  project config (configs/00_project.yaml)  >  code defaults
```

- **Stage configs**: one file per step, containing only that step's parameters;
- **Project config** `00_project.yaml`: only paths referenced by several steps (EGFAnalysisPy directory, DSurfTomo executable, etc.). Currently only the step-05 and 06b entries read it; steps 01–04 are driven entirely by their own stage configs;
- Every `*.example.yaml` is a template: copy it without the `.example` suffix (or any name you like) and edit the copy — never the template itself.

```bash
cp configs/02_correlation.example.yaml configs/02_correlation.yaml
python scripts/02_correlation.py --config configs/02_correlation.yaml
```

## File inventory

| Config file | Serves | Notes |
|-------------|--------|-------|
| `00_project.example.yaml` | 05, 06b | project-wide shared paths |
| `01_clean.example.yaml` | `scripts/01_clean.py` | keys mirror the CLI options |
| `02_correlation.example.yaml` | `scripts/02_correlation.py` | NoisePy correlation & stacking parameters |
| `03_stack_qc.example.yaml` | `scripts/03_stack_qc.py` | moveout visualization parameters |
| `03b_stack_panels.example.yaml` | `scripts/03b_stack_panels.py` | panel figure parameters |
| `04_export_egf.example.yaml` | `scripts/04_export_egf.py` | keys mirror the CLI options |
| `05_dispersion.example.yaml` | `scripts/05_dispersion.py` | all dispersion parameters (replaces editing config.py) |
| `06a_line_inversion/line.example.yaml` | both 06a entries | one config per line (shared by invert & stitch) |
| `06b_area_inversion/dsurftomo.example.yaml` | all three 06b entries | one config per inversion project |

## Key naming

For steps 01–04, YAML keys are **spelled the same** as the long CLI options (`input-root`, `stack-dir`, …; underscores also accepted). Explicit CLI flags always override the config file; unknown keys produce a warning and are ignored. Steps 05, 06a and 06b use sectioned structures (`paths:` / `run:` / `prepare:` …) documented inline in each example file.

## Runtime-config mechanism

Where an underlying program only accepts its own fixed format, the entry script does the translation, so users only maintain the root-level YAML:

- **Step 05**: `05_dispersion.py` renders the YAML into `<result_path>/runtime_config.py` (subclassing `EGFAnalysisPy/config.py`, overriding parameters, recomputing derived quantities; the header records the source config and creation time), then hands it to `EGFAnalysisPy/run_parallel.py --config`. The bundled template is never modified; the runtime config in the result directory is a complete parameter archive of that run;
- **06b-1**: `01_prepare.py` translates the `prepare:` section into the command line of `DSurfTomo/scripts/PrepareDSurfTomoProject.py` (generated files such as `DSurfTomo.in` all go into the project directory);
- **06b-2/3**: the executable path comes from `paths.dsurftomo_bin` (stage config), falling back to `dsurftomo_bin` in `00_project.yaml` — a direct instance of the precedence rule.

## Overwrite protection

Entries never overwrite existing results by default: step 01 requires `--overwrite` to rewrite outputs; step 05 aborts on a non-empty `qc_result/`, and 06b-1 aborts on a non-empty project directory unless `--overwrite` is given.

## Path conventions

- `/path/to/...` values in the examples are placeholders — replace them with your real paths;
- Relative paths are resolved against the **directory you run the command from** (always run from the repo root);
- Do not commit real data paths back to the repository — keep personal copies (without `.example`) out of Git or add them to `.gitignore`.
