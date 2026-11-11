# Quick start

[中文](quickstart.zh.md) | English

This guide covers installation, a first run on synthetic data with a known answer, and the setup of a project with your own data; allow about half an hour. Details of every step are in the [tutorial](tutorial.md), the principles in [methods](methods.md).

## 1. Installation (about 10 minutes)

Requirements: Linux, conda (or mamba) and git.

```bash
git clone https://github.com/Shubin-Wei/local_ANT.git && cd local_ANT
conda env create -f environment.yml && conda activate lant   # also installs lant from this checkout
lant install-external                                        # downloads and builds DSurfTomo and pDSurfTomo
lant check                                                   # every line should read ok
```

`lant check` names the install command for anything missing.

## 2. Verify the installation with the demo (about 5 minutes)

```bash
lant demo --out ~/lant_demo
```

The demo generates a small synthetic array (plane Rayleigh waves in a known 1-D model), runs the whole chain and prints the location of the report, `~/lant_demo/output/REPORT.md`. It confirms the judgement points itself, which is legitimate only because it generated the data; it refuses to run in a directory it did not create.

## 3. A run step by step

To see where the chain stops for the user, generate the same data and run them yourself:

```bash
python -m lant.synthetic --quick --out ~/lant_quick          # 10 stations, 6 hours, about 1 minute
lant run all -c ~/lant_quick/synthetic.yaml
```

The chain stops where a scientific judgement is needed and prints the files to look at and the command that continues. With synthetic data the answer is known; **with real data, each decision is yours, made after looking at the figures**.

| Step | Question | Evidence | Answer for the synthetic data |
|---|---|---|---|
| `clean` | Without response files: are all stations of the same instrument model? | `work/01_clean/report.txt` | Yes |
| `stack` (when candidates are listed) | Remove the listed station sub-stacks? | `work/03_stack/station_qc/station_qc.png` | No faulty stations: keep them |
| `reference` | Do the F-J picks follow the fundamental-mode ridge? | `work/04_reference/fj.png` | Yes |
| `inversion` | Are the 1-D starting model and its fit to the average curve reasonable? | `results/figures/4_start_model.png` | Yes |

Record each decision and run the chain again; it continues where it stopped:

```bash
lant confirm reference -c ~/lant_quick/synthetic.yaml --reason "synthetic data with a known answer"
lant run all -c ~/lant_quick/synthetic.yaml
lant status -c ~/lant_quick/synthetic.yaml                     # state of every step and the next command
```

When the chain has finished, compare the result with the truth:

```bash
python examples/synthetic/compare_truth.py --dir ~/lant_quick
```

It reports the difference between measured and true phase velocity per frequency, and between the inverted and the true model. Six hours of records give rougher numbers than the full example (16 stations, 3 days, about 300 MB and an hour of processing; leave out `--quick`); see [validation](validation.md).

## 4. Your own data

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

Then start the chain; the stops are those of section 3:

```bash
lant run all -c /data/result/config.yaml
```

`lant run all --raw /data/raw --stations /data/stations.csv --out /data/result` writes the same configuration and starts at once. Problems on the way are covered in [troubleshooting](troubleshooting.md).

## 5. After a run

Start with `REPORT.md` in the output directory: the state of every step, key numbers, the items left for you and the main figures, with links to the results. Then check:

| File | Purpose |
|---|---|
| `.lant/steps/<step>.md` | What each step did, the parameters it used (value and source), what is left for your judgement |
| `work/01_clean/coverage.png` | Data available per station and day |
| `work/04_reference/fj.png` | Whether the reference curve follows the fundamental mode |
| `results/figures/3_dispersion_curves.png` | All dispersion curves and the reference |
| `results/figures/5_depth_slices.png`, `5_sections.png`, `6_checkerboard_*.png` | The model and the scales it resolves; read them together |

Every output file is described in [output files](outputs.md).
