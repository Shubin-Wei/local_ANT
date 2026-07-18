# Complete Ambient Noise Imaging Workflow (01–06)

This project implements a complete ambient noise imaging pipeline from raw miniSEED data to shear-wave velocity models. Stages 01–05 are common; after step 05 the pipeline splits into two parallel branches — **06a linear-array inversion** or **06b areal-array inversion** — chosen by array geometry.

All commands are executed from the **repository root**; every step runs independently and can be resumed midway.

## Pipeline overview

```
Raw miniSEED data
  │
  ▼ [01] scripts/01_clean.py                data cleaning (NoisePy naming format)
  ▼ [02] scripts/02_correlation.py          cross-correlation + stacking (NoisePy)
  │        prerequisite: station.csv
  ▼ [03] scripts/03_stack_qc.py             stack QC visualization (optional, recommended)
  │      scripts/03b_stack_panels.py        publication-grade panel figure
  ▼ [04] scripts/04_export_egf.py           export EGFAnalysisPy input (.dat)
  ▼ [05] scripts/05_dispersion.py           dispersion analysis + DisperPicker picking
  │                                         → qc_result/CDisp.*.txt, GDisp.*.txt
  │
  ├─▶ [06a] scripts/06a_line_inversion/     linear array: evodcinv 1-D inversion → pseudo-2D section
  │         (see docs/line_inversion.md)
  └─▶ [06b] scripts/06b_area_inversion/     areal array: DSurfTomo 3-D inversion
            (see docs/area_inversion.md)
```

### Choosing the inversion branch

| | 06a linear inversion | 06b areal inversion |
|---|---|---|
| Array geometry | stations roughly along one line | stations covering a 2-D region |
| Method | evodcinv: per-pair 1-D Vs(z) inversion (CPSO), stitched at pair midpoints | DSurfTomo: direct inversion for 3-D Vs |
| Input | step-05 `qc_result/` (used directly) | step-05 `qc_result/` (converted to surfdata.dat by 06b-1) |
| Result | pseudo-2D Vs depth section (section.png/.npz) | 3-D Vs model + depth-slice figures |
| Extra deps | evodcinv (into the noisepy env, INSTALL §2e) | compiled DSurfTomo (INSTALL §4), GMT 6 |

The branches are fully independent; you may also run both on the same dataset.

### Environment and I/O per step

| Step | User entry | Env | Main input | Main output |
|------|------------|-----|------------|-------------|
| 01 | `scripts/01_clean.py` | noisepy | raw miniSEED | cleaned `.ms` |
| 02 | `scripts/02_correlation.py` | noisepy | cleaned data + station.csv | CCF + Stack (HDF5) |
| 03 | `scripts/03_stack_qc.py`, `03b_stack_panels.py` | noisepy | Stack | moveout figures |
| 04 | `scripts/04_export_egf.py` | noisepy | Stack | CF `.dat` + stations.txt |
| 05 | `scripts/05_dispersion.py` | egfanalysispy | CF `.dat` | dispersion curves `qc_result/` |
| 06a | `scripts/06a_line_inversion/01_invert.py`, `02_stitch.py` | noisepy + evodcinv | `qc_result/` | 1-D model set + pseudo-2D section |
| 06b | `scripts/06b_area_inversion/01_prepare.py`, `02_invert.py`, `03_plot.py` | noisepy + DSurfTomo binary + GMT | `qc_result/` | 3-D model + slice figures |

### Old script names → new entries (migration guide)

| Old name | New name |
|----------|----------|
| `scripts/clean_mseed_for_noisepy.py` | `scripts/01_clean.py` |
| `scripts/local_cc_workflow.py` | `scripts/02_correlation.py` |
| `scripts/stack_visualization.py` | `scripts/03_stack_qc.py` |
| `scripts/plot_stack_panels.py` | `scripts/03b_stack_panels.py` |
| `scripts/export_egfanalysis_input.py` | `scripts/04_export_egf.py` |
| editing `EGFAnalysisPy/config.py` + `python run.py` | `scripts/05_dispersion.py --config configs/05_dispersion.yaml` |
| `lant_invert.py` / `lant_stitch.py` (local project copies) | `scripts/06a_line_inversion/01_invert.py` / `02_stitch.py` |
| `configs/local_cc_workflow.example.yaml` etc. | `configs/02_correlation.example.yaml` etc. (matching numbers) |

Configuration behaviour is unchanged: contents of old config files can be copied verbatim into their renamed counterparts.

## Package citations

This project uses the following open-source packages — many thanks to their authors:

- **NoisePy** — ambient noise cross-correlation
  Jiang, C., & Denolle, M. A. (2020). NoisePy: A New High-Performance Python Tool for Ambient-Noise Seismology. *Seismological Research Letters*, 91(3), 1853–1866. https://doi.org/10.1785/0220190364
  GitHub: https://github.com/noisepy/NoisePy

- **EGFAnalysisPy** — surface-wave dispersion analysis (pure-Python rewrite of EGFAnalysisTimeFreq, DisperPicker bundled)
  Yao, H., van der Hilst, R. D., & de Hoop, M. V. (2006). Surface-wave array tomography in SE Tibet from ambient seismic noise and two-station analysis—I. Phase velocity maps. *Geophysical Journal International*, 166(2), 732–744. https://doi.org/10.1111/j.1365-246X.2006.03028.x
  Author: Huajian Yao (USTC)

- **DisperPicker** — deep-learning dispersion-curve extraction (bundled inside EGFAnalysisPy)
  Yang, S., Zhang, H., Gu, N., Gao, J., Xu, J., Jin, J., Li, J., & Yao, H. (2022). Automatically Extracting Surface Wave Group and Phase Velocity Dispersion Curves from Dispersion Spectrograms Using a Convolutional Neural Network. *Seismological Research Letters*. https://doi.org/10.1785/0220210280
  GitHub: https://github.com/ShaoboYang-USTC/DisperPicker

- **DSurfTomo** — direct inversion for 3-D shear-wave structure (branch 06b only)
  Fang, H., Yao, H., Zhang, H., Huang, Y. C., & van der Hilst, R. D. (2015). Direct inversion of surface wave dispersion for three-dimensional shallow crustal structure based on ray tracing: methodology and application. *Geophysical Journal International*, 201(3), 1251–1263. https://doi.org/10.1093/gji/ggv069

- **evodcinv / disba** — 1-D dispersion inversion with evolutionary algorithms (branch 06a only)
  Luu, K. GitHub: https://github.com/keurfonluu/evodcinv

## Environment setup

See [INSTALL_en.md](INSTALL_en.md) for details. Summary:

```bash
conda env create -f environments/noisepy.yml        # entries for 01–04, 06a, 06b
conda env create -f environments/egfanalysispy.yml  # step 05
cd DSurfTomo && python configure && cd ..           # only needed for 06b
conda activate noisepy && python tools/check_environment.py
```

The pretrained model ships in `EGFAnalysisPy/DisperPicker/saver/` — no download needed.

---

## 01 Data cleaning

**Purpose:** drop incomplete days and convert raw miniSEED into NoisePy naming format.

**Entry:** `scripts/01_clean.py` (env: noisepy)

**Input:** raw miniSEED root with one subdirectory per station:

```
RawData/
├── CF.17982/
│   ├── 230901000000    # filename is a Beijing timestamp YYMMDDHHmmss
│   └── ...
├── CF.58808/
└── CN.20496/
```

**Commands:**

```bash
# with a config file (recommended; first cp configs/01_clean.example.yaml configs/01_clean.yaml)
python scripts/01_clean.py --config configs/01_clean.yaml

# or pure CLI (CLI flags override the config file)
python scripts/01_clean.py \
    --input-root  TestData/RawData \
    --output-root TestData/CleanData \
    --target-channel Z
```

**Options:**

| Option | Description |
|--------|-------------|
| `--config` | YAML file whose keys mirror the options below |
| `--input-root` | raw miniSEED root directory |
| `--output-root` | cleaned output directory |
| `--min-completeness` | minimum per-day completeness (0–1, default 0.8) |
| `--target-channel` | channel/component, e.g. `Z`, `BNZ` (default: auto Z) |
| `--date-range START END` | only dates strictly between the two (exclusive), format YYMMDD |
| `--overwrite` | overwrite existing outputs (off by default) |
| `--dry-run` | print actions without writing |
| `--verbose` | per-file logging |

**Output:** `CleanData/CF17982_BNZ00_2023258.ms ...`, named `{net}{sta}_{chan}{loc}_{year+julday}.ms`

**Prerequisite for step 02: station coordinate CSV**

```csv
network,station,channel,latitude,longitude,elevation
CF,17982,BNZ,42.323349,118.546867,844
CF,58808,BNZ,42.321987,118.54599,822
```

Network/station must match the miniSEED headers; decimal degrees; elevation in metres.

---

## 02 Cross-correlation and stacking

**Purpose:** compute and stack ambient-noise cross-correlation functions with NoisePy.

**Entry:** `scripts/02_correlation.py` (env: noisepy)

```bash
cp configs/02_correlation.example.yaml configs/my_02.yaml   # edit paths & params
```

**Key parameters:**

```yaml
data-dir: ../TestData/CleanData       # step-01 output
station-csv: ../TestData/station.csv
ccf-dir: ../TestData/CCFData
stack-dir: ../TestData/StackData

freqmin: 0.2          # Hz
freqmax: 10           # Hz
sampling-rate: 25.0   # Hz
maxlag: 100           # s
inc-hours: 24
cc-len: 600           # s
step: 300.0           # s
time-norm: rma
freq-norm: rma
stack-method: all     # linear / pws / robust / all
```

**Commands:**

```bash
python scripts/02_correlation.py --config configs/my_02.yaml --validate-only
python scripts/02_correlation.py --config configs/my_02.yaml --stage cc
python scripts/02_correlation.py --config configs/my_02.yaml --stage stack
python scripts/02_correlation.py --config configs/my_02.yaml
```

---

## 03 Stack QC visualization

**Purpose:** inspect stack quality and choose the stacking method for later steps. Optional but strongly recommended.

**Entry:** `scripts/03_stack_qc.py` (env: noisepy)

```bash
python scripts/03_stack_qc.py \
    --stack-dir TestData/StackData \
    --freqmin 1 --freqmax 5 \
    --component ZZ --stack-method linear --maxlag 10 --style both
# or
python scripts/03_stack_qc.py --config configs/03_stack_qc.yaml
```

**Output:** moveout figures (written to `--outdir`).

**Publication-grade panels (optional):** `scripts/03b_stack_panels.py` builds a 2×N panel moveout figure (N = number of stack methods; image row + wiggle/symmetric row):

```bash
python scripts/03b_stack_panels.py --stack-dir TestData/StackData --methods linear pws robust --component ZZ
# or --config configs/03b_stack_panels.yaml
```

**Choosing a stack method:** `linear` is the workhorse; `pws` boosts SNR; `robust` resists outliers; run `all` in step 02 and compare here.

---

## 04 Export EGFAnalysisPy input

**Purpose:** convert NoisePy stacks (HDF5) into EGFAnalysisPy's text format (CF `.dat`).

**Entry:** `scripts/04_export_egf.py` (env: noisepy)

```bash
python scripts/04_export_egf.py \
    --stack-dir  TestData/StackData \
    --stack-method all --component ZZ \
    --output-root TestData/EGFAnalysisInputData
# or --config configs/04_export_egf.yaml
```

**Options:**

| Option | Description |
|--------|-------------|
| `--config` | YAML config file |
| `--stack-dir` | NoisePy stack directory |
| `--stack-method` | `linear` / `pws` / `robust` / `all` (default) |
| `--component` | component, default `ZZ` |
| `--station` | optional; only pairs involving this station |
| `--output-root` | output root directory |

**Output layout:**

```
EGFAnalysisInputData/
├── Allstack_linear/ZZ/
│   ├── CF17982_CF58808.dat    # waveform (3 columns: time, gab, gba)
│   ├── stations.txt
│   └── EGFfiles.txt
├── Allstack_pws/ZZ/
└── Allstack_robust/ZZ/
```

`.dat` format: lines 1–2 hold `lon lat elev` of the source/receiver stations; data lines are `time(s) gab gba`.

---

## 05 Dispersion analysis and curve picking

**Purpose:** time–frequency analysis of CF/EGF into dispersion energy diagrams, automatic group/phase-velocity picking with the bundled DisperPicker (CNN), quality control, and output — one command, no MATLAB, no editing files inside EGFAnalysisPy.

**Entry:** `scripts/05_dispersion.py` (env: **egfanalysispy**)

```bash
cp configs/05_dispersion.example.yaml configs/my_05.yaml   # edit
conda activate egfanalysispy
python scripts/05_dispersion.py --config configs/my_05.yaml
```

**Key parameters (`configs/05_dispersion.example.yaml`; meanings follow EGFAnalysisPy/config.py):**

```yaml
paths:
  waveform_path: /path/to/EGFAnalysisInputData/Allstack_robust/ZZ  # step-04 output
  result_path:   /path/to/DisperPicker_result
run:
  parallel: false     # true = multiprocessing (run_parallel.py)
  workers: 4
dispersion:
  StartT: 0.4         # period range & step — step 06 will need these
  EndT: 3.0
  DeltaT: 0.02
  StartV: 0.5
  EndV: 3.5
  DeltaV: 0.005       # DisperPicker requires exactly 0.005 km/s
```

The entry renders the YAML into a runtime config (`<result_path>/runtime_config.py`, with source and timestamp for provenance) and invokes `EGFAnalysisPy/run_parallel.py`; the bundled `EGFAnalysisPy/config.py` is never modified. Internal flow: energy diagram → DisperPicker picking → group-velocity-windowed time-variable filtering and re-computed phase diagram → second picking → QC → output.

**Advanced:** you may still run `cd EGFAnalysisPy && python run_parallel.py --config your_config.py -w 8` (or edit `config.py` and `python run.py`) — behaviour is identical to upstream.

**Output:**

```
<result_path>/
├── runtime_config.py          # full parameters of this run (auto-generated)
├── qc_result/
│   ├── CDisp.FD01_FD16.txt    # phase-velocity curves
│   ├── GDisp.FD01_FD16.txt    # group-velocity curves
│   └── ...
├── qc_plot/                   # energy diagrams + picking comparison (sampled)
├── all_C.png / mean_C.png
└── all_G.png / mean_G.png
```

Curve file format (lines 1–2: `lon lat` of the two stations; data: `period velocity SNR probability`):

```
117.43333100    31.81069500
117.48089000    31.66345700
0.20   1.250   6.100   0.950
0.30   1.310   5.800   0.880
```

> `qc_result/` is the common input of both inversion branches: 06a reads it directly; 06b converts it via `01_prepare.py`. Station coordinates are embedded in the file headers — no extra station file needed.

---

## 06a Linear-array inversion (evodcinv)

Use when stations lie roughly along one line. Two steps: per-pair 1-D inversion → stitch a pseudo-2D section.

```bash
cp configs/06a_line_inversion/line.example.yaml configs/06a_line_inversion/my_line.yaml
conda activate noisepy    # evodcinv required (INSTALL §2e)
python scripts/06a_line_inversion/01_invert.py --config configs/06a_line_inversion/my_line.yaml
python scripts/06a_line_inversion/02_stitch.py --config configs/06a_line_inversion/my_line.yaml
```

Full details (applicability, config walkthrough, result interpretation): [line_inversion.md](line_inversion.md).

## 06b Areal-array inversion (DSurfTomo)

Use when stations cover a 2-D region. Three steps: prepare project → run inversion → GMT plotting.

```bash
cp configs/06b_area_inversion/dsurftomo.example.yaml configs/06b_area_inversion/my_area.yaml
conda activate noisepy    # compiled DSurfTomo required (INSTALL §4); plotting needs GMT 6
python scripts/06b_area_inversion/01_prepare.py --config configs/06b_area_inversion/my_area.yaml
python scripts/06b_area_inversion/02_invert.py  --config configs/06b_area_inversion/my_area.yaml
python scripts/06b_area_inversion/03_plot.py    --config configs/06b_area_inversion/my_area.yaml
```

The interactive scripts under `DSurfTomo/scripts/` remain available for manual step-by-step runs (`PrepareDSurfTomoProject.py` all-in-one, or `extractnew.py` → `GenerateIniMOD.py` → `GenerateDSurfTomoInputFile.py`); see [area_inversion.md](area_inversion.md).

---

## FAQ

Frequently asked questions are collected in [FAQ.md](FAQ.md).

## License

Each package keeps its original license. Helper scripts (`scripts/` and `DSurfTomo/scripts/`) are released with this project.
