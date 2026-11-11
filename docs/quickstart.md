# Quick start

[中文](quickstart.zh.md) | English

This guide takes about half an hour: installation, a complete run on a small synthetic data set with a known answer, then the switch to your own data. For the principles of each step see [methods.md](methods.md); for the details of operation see [tutorial.md](tutorial.md).

## 1. Installation (about 10 minutes)

You need Linux, conda (or mamba) and git.

```bash
git clone https://github.com/Shubin-Wei/local_ANT.git && cd local_ANT
conda env create -f environment.yml && conda activate lant
pip install -e .
lant install-external      # downloads and builds DSurfTomo and pDSurfTomo
lant check                          # everything should read ok; missing pieces come with their install command
```

## 2. A first run on synthetic data (about 10 minutes)

The synthetic data are a noise field of plane Rayleigh waves in a known 1-D model, so the results can be compared with the truth.

**In one command:**

```bash
lant demo --out ~/lant_demo        # writes the data, runs the whole chain, prints where the report is
```

The demo confirms the judgement points itself, which is right only because it made the data; it refuses to work in a directory it did not create. Open `~/lant_demo/output/REPORT.md` when it is done.

**Step by step**, to see where the chain stops for you:

```bash
python -m lant.synthetic --quick --out ~/lant_quick     # 10 stations, 6 hours, about 1 minute
lant run all -c ~/lant_quick/synthetic.yaml
```

The chain stops where a scientific judgement is needed and prints the files to look at and the next command. The answer of synthetic data is known, so you can confirm directly; **on real data, decide each point yourself after looking at the figures**.

| Stops at | Question to answer | Look at | On the synthetic data |
|---|---|---|---|
| clean | without response files: are all stations of the same instrument model? | `work/01_clean/report.txt` | yes: `lant confirm clean -c … --reason "synthetic data with a known answer"` |
| stack (when there are candidates) | remove the station sub-stacks listed? | `work/03_stack/station_qc/station_qc.png` | there are no faulty stations, keep them: `lant confirm stack -c … --reason "synthetic data with a known answer"` |
| reference | do the F-J picks follow the fundamental-mode ridge? | `work/04_reference/fj.png` | yes: `lant confirm reference -c … --reason "synthetic data with a known answer"` |
| inversion | are the 1-D starting model and its fit to the average curve reasonable? | `results/figures/4_start_model.png` | yes: `lant confirm inversion -c … --reason "synthetic data with a known answer"` |

After each confirmation run `lant run all -c ~/lant_quick/synthetic.yaml` again; it continues where it stopped. At any time

```bash
lant status -c ~/lant_quick/synthetic.yaml     # state of each step, reason and next command
```

When it is done, compare with the truth:

```bash
python examples/synthetic/compare_truth.py --dir ~/lant_quick
```

It gives the difference between measured and true phase velocity at each frequency, and between the inverted and the true model. These quick data hold only 6 hours of records, so the numbers are rougher than those of the full example (see [validation.md](validation.md)); for the full example leave out `--quick` (16 stations, 3 days, about 300 MB, about an hour of processing).

## 3. Your own data

Only two things are needed:

- **A raw-data directory**: any folder structure and file names, any format ObsPy reads (miniSEED, SAC, …), mixed if you like.
- **A station table**: CSV, text or Excel with at least station name, longitude and latitude (English and Chinese column names are recognised). For nodes that were moved, add serial numbers and deployment start and end ([tutorial.md](tutorial.md), section 2.2).

```bash
lant run all --raw /data/raw --stations /data/stations.csv --out /data/result
```

This writes `/data/result/config.yaml` and starts the run. Two settings are worth changing before you start:

- `data.velocity_range_kms`: the approximate range of surface-wave phase velocities in the area. The default 0.1–5 km/s covers soft soil to hard rock, but the lag window is inversely proportional to the slowest velocity, so a wide range makes the correlations much more expensive in time and disk.
- `data.response`: give the instrument response files if you have them; without them all stations must be of the same model.

Then `lant run all -c /data/result/config.yaml`. The stops are the same as in the synthetic example, but at each one the judgement is yours; [troubleshooting.md](troubleshooting.md) helps with problems on the way, and [outputs.md](outputs.md) explains the result files.

## 4. Files worth a look after every run

| File | Why |
|---|---|
| `.lant/steps/<step>.md` | what each step did, which parameters it used (value and source), what is left for your judgement |
| `work/01_clean/coverage.png` | how much data each station has |
| `work/04_reference/fj.png` | whether the reference curve follows the fundamental mode |
| `results/figures/3_dispersion_curves.png` | all dispersion curves and the reference |
| `results/figures/5_depth_slices.png`, `5_slices_3d.png` (lines: `5_sections.png`), `6_checkerboard_*.png` | the model and the scales it can resolve (look at them together) |
