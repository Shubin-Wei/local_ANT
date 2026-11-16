# Quick start

[中文](quickstart.zh.md) | English

This guide covers installation, the setup of a project with your own data and the points where the chain stops for your judgement. Details of every step are in the [tutorial](tutorial.md), the principles in [methods](methods.md).

## 1. Installation (about 10 minutes)

Requirements: Linux, conda (or mamba) and git.

```bash
git clone https://github.com/Shubin-Wei/local_ANT.git && cd local_ANT
conda env create -f environment.yml && conda activate lant   # also installs lant from this checkout
lant check                                                   # every line should read ok
```

`lant check` names the install command for anything missing.

The inversion is the Python port of DSurfTomo in this same environment; no separate solver installation is needed. Its comparison with the original is in [python_backends.md](python_backends.md).

If an older Conda stays in dependency solving, a machine with Python ≥ 3.11 and GMT 6 can also use an isolated Python environment from the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
lant check
```

GMT and Ghostscript must be installed separately and on PATH; `lant check` checks dependencies. This leaves other Conda environments intact.

## 2. Your data

Two inputs are needed:

- **Raw records**: one directory with any folder structure and file names, in any format ObsPy reads (miniSEED, SAC, …), formats may be mixed.
- **Station table**: CSV, text or Excel with at least station name, longitude and latitude; English and Chinese column names are recognised. For nodes that were moved, add serial numbers and deployment times ([tutorial](tutorial.md), section 2.2).

Write a configuration and review it before the first run:

```bash
lant init -o /data/result/config.yaml --raw /data/raw --stations /data/stations.csv --out /data/result
lant validate -c /data/result/config.yaml      # reads the configuration, the station table and the raw directory
lant plan -c /data/result/config.yaml          # what every step will do, with which values
```

Two settings deserve attention:

- `data.velocity_range_kms`: the approximate range of surface-wave phase velocities in the area. The default 0.1–5 km/s covers soft soil to hard rock; the lag window is inversely proportional to the slowest velocity, so a range wider than necessary costs computing time and disk. Narrow it only on independent knowledge of the area.
- `data.response`: instrument response files, if available; without them all stations must be of one instrument model.

`lant run all --raw /data/raw --stations /data/stations.csv --out /data/result` writes the same configuration and starts at once.

## 3. Where the chain stops for you

```bash
lant run all -c /data/result/config.yaml
```

The chain stops where a scientific judgement is needed and prints the files to look at and the command that continues. **Each decision is yours, made after looking at the figures.**

| Step | Question | Evidence |
|---|---|---|
| `clean` | Without response files: are all stations of the same instrument model? | `work/01_clean/report.txt` |
| `stack` (when candidates are listed) | Remove the listed station sub-stacks? | `results/figures/2_station_qc.png` |
| `reference` | Do the F-J picks follow the fundamental-mode ridge? | `results/figures/3_fj_phase_velocity.png` |
| `inversion` | Are the 1-D starting model and its fit to the average curve reasonable? | `results/figures/4_start_model.png` |

Record each decision with its evidence and run the chain again; it continues where it stopped:

```bash
lant confirm reference -c /data/result/config.yaml --reason "picks follow the fundamental ridge in 3_fj_phase_velocity.png"
lant run all -c /data/result/config.yaml
lant status -c /data/result/config.yaml        # state of every step and the next command
```

Problems on the way are covered in [troubleshooting](troubleshooting.md).

## 4. After a run

Start with `REPORT.md` in the output directory: the state of every step, key numbers, the items left for you and the main figures, with links to the results. Then check:

| File | Purpose |
|---|---|
| `.lant/steps/<step>.md` | What each step did, the parameters it used (value and source), what is left for your judgement |
| `results/figures/0_data_coverage.png` | Data available per station and day |
| `results/figures/3_fj_phase_velocity.png` | Whether the reference curve follows the fundamental mode |
| `results/figures/3_dispersion_curves.png` | All dispersion curves and the reference |
| `results/figures/5_depth_slices.png`, `5_sections.png`, `6_checkerboard_*.png` | The model and the scales it resolves; read them together |

Every output file is described in [output files](outputs.md).

## 5. Edit a step, then redraw its results

The common calculation and display settings are expanded in `config.yaml`; normally only a few need editing. See [parameters](parameters.md) for advanced settings, or generate a commented reference with `lant init --full -o reference.yaml`. Do not maintain two active configurations.

These commands work from any directory after activating the environment; use an absolute configuration path:

```bash
lant plan cc -c /data/result/config.yaml       # this step and changes since its recorded run
lant run cc -c /data/result/config.yaml        # recompute just this step
lant params cc -c /data/result/config.yaml     # recorded versus next settings
```

Each command freezes its configuration at the start; edits during execution belong to subsequent commands. Editing settings does not change existing results. Use `lant status -c CONFIG` to see which downstream steps need recomputing.

Correlation display controls live under `plots`, independently of calculation settings under `correlation`. For example, edit the existing `plots:` section (do not add a duplicate section):

```yaml
plots:
  ccf_bands_hz: [[1.0, 2.0], [2.0, 5.0]]
  ccf_stack_band_hz: [1.0, 5.0]
  ccf_maxlag_s: 10.0
```

These numbers illustrate syntax, not recommended scientific settings. Bands and lags must fit the saved correlation product. `[]` and `null` retain default display choices.

```bash
lant plot ccf_moveout -c /data/result/config.yaml
lant params results/figures/2_ccf_moveout_1-2Hz.png -c /data/result/config.yaml
```

A redraw updates selected figures without recalculating correlations. `lant plot --help` explains the command; `lant plot -c CONFIG` redraws all available figures. See [outputs](outputs.md) for names. After editing `plots.vs_range` or `slice_depths_km`, use `lant plot depth_slices -c CONFIG` for slices alone.

`REPORT.md` links parameter sources. `results/parameters.yaml` indexes actual execution records, upstream sources and display settings for the data and each image; it is not the configuration for the next command. Partial redraws leave other images' records intact. Failed redraws retain the old image and its old record. History preserves parameters, not every overwritten large data file or historical image; missing provenance from older projects is never invented.

The terminal shows the stage, completed work, elapsed time and a percentage when the total is known; unknown work uses a heartbeat. Detailed output goes to `.lant/lant.log`; add `-v` to display it. Start with `REPORT.md` and `results/figures/`, without searching internal logs.
