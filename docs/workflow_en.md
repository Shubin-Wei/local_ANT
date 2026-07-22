# Complete Workflow: Ambient Noise Surface Wave Tomography from Local miniSEED Data

This project implements a complete ambient noise surface wave tomography pipeline — from raw miniSEED seismic data to a shear-wave velocity model. The workflow consists of 6 steps, using 3 software packages and several helper scripts. **EGFAnalysisPy** is a pure-Python rewrite of EGFAnalysisTimeFreq with the deep-learning picker DisperPicker (and its pretrained model) bundled in, performing dispersion analysis and curve extraction in a single step — no MATLAB required.

Steps 1–5 are shared by every use case: data cleaning → cross-correlation and stacking → dispersion extraction. Step 6 offers **two parallel inversion paths**; pick one based on your array geometry:

| Path | Use when | Product |
|------|----------|---------|
| **Step 6A** DSurfTomo | Areal array whose ray coverage can constrain a 2D grid | 3D shear-wave velocity model |
| **Step 6B** 1D inversion + section stitching | Linear (profile) array, shallow structure | Pseudo-2D shear-wave velocity section |

In a linear array the ray paths are nearly collinear and cannot constrain a 2D grid, so DSurfTomo's tomographic inversion does not apply. Use Step 6B instead: invert each station pair independently in 1D, then stitch the results along the profile.

## Workflow Overview

```
Raw miniSEED data
  │
  ▼ [Step 1] clean_mseed_for_noisepy.py
Cleaned miniSEED data (NoisePy format)
  │
  ▼ [Step 2] Prepare station coordinate file (station.csv)
  │
  ▼ [Step 3] local_cc_workflow.py (NoisePy)
Cross-correlation functions (CCF) + stacked results
  │    └─ [Optional] stack_visualization.py for QC
  │
  ▼ [Step 4] export_egfanalysis_input.py
EGFAnalysisPy input format (.dat waveforms + stations.txt + EGFfiles.txt)
  │
  ▼ [Step 5] EGFAnalysisPy (run.py or run_parallel.py)
Dispersion analysis + DisperPicker auto-picking + QC → dispersion curves (qc_result/CDisp.*.txt, GDisp.*.txt)
  │
  ├─▼ [Step 6A] DSurfTomo inversion (areal array → 3D model)
  │   ├── extractnew.py               → surfdata.dat
  │   ├── GenerateIniMOD.py           → MOD (initial model)
  │   ├── GenerateDSurfTomoInputFile.py → DSurfTomo.in
  │   ├── Run DSurfTomo inversion
  │   └── GMT scripts for plotting (depth slices, cross-sections, ray paths, checkerboard tests)
  │
  └─▼ [Step 6B] 1D inversion + stitching (linear array → pseudo-2D section)
      ├── lant_invert.py              → per-pair 1D inversion, writes models/
      └── lant_stitch.py              → stitch along the profile → section.png
```

## Citations

This project integrates the following open-source software packages. All credit goes to their original authors:

- **NoisePy** — Ambient noise cross-correlation
  Jiang, C., & Denolle, M. A. (2020). NoisePy: A New High-Performance Python Tool for Ambient-Noise Seismology. *Seismological Research Letters*, 91(3), 1853–1866. https://doi.org/10.1785/0220190364
  GitHub: https://github.com/noisepy/NoisePy

- **EGFAnalysisPy** — Surface wave dispersion analysis (pure-Python rewrite of EGFAnalysisTimeFreq, DisperPicker bundled)
  Yao, H., van der Hilst, R. D., & de Hoop, M. V. (2006). Surface-wave array tomography in SE Tibet from ambient seismic noise and two-station analysis—I. Phase velocity maps. *Geophysical Journal International*, 166(2), 732–744. https://doi.org/10.1111/j.1365-246X.2006.03028.x
  Author: Huajian Yao (University of Science and Technology of China)

- **DisperPicker** — CNN-based automatic dispersion curve extraction (bundled within EGFAnalysisPy)
  Yang, S., Zhang, H., Gu, N., Gao, J., Xu, J., Jin, J., Li, J., & Yao, H. (2022). Automatically Extracting Surface Wave Group and Phase Velocity Dispersion Curves from Dispersion Spectrograms Using a Convolutional Neural Network. *Seismological Research Letters*. https://doi.org/10.1785/0220210280
  GitHub: https://github.com/ShaoboYang-USTC/DisperPicker

- **DSurfTomo** — Direct surface wave tomography for 3D shear-wave velocity
  Fang, H., Yao, H., Zhang, H., Huang, Y. C., & van der Hilst, R. D. (2015). Direct inversion of surface wave dispersion for three-dimensional shallow crustal structure based on ray tracing: methodology and application. *Geophysical Journal International*, 201(3), 1251–1263. https://doi.org/10.1093/gji/ggv069

## Prerequisites

### Python Environment

```bash
conda create -n noisepy -y python=3.10 pip
conda activate noisepy
pip install noisepy-seis
pip install obspy numpy pandas pyyaml datetimerange h5py scipy matplotlib
```

### EGFAnalysisPy Environment

EGFAnalysisPy bundles DisperPicker and depends on TensorFlow 1.13, so it needs a separate Python 3.6 environment:

```bash
conda create -n egfanalysispy -y python=3.6 pip
conda activate egfanalysispy
cd EGFAnalysisPy
pip install -r requirements.txt
```

The pretrained model is already included in `EGFAnalysisPy/DisperPicker/saver/` — no extra download needed.

### Compile DSurfTomo

```bash
cd DSurfTomo
python configure          # configure is a Python script; it runs make in src/ and copies the binary to bin/
```

### GMT (Optional, for Plotting)

The DSurfTomo plotting scripts require GMT 6+.

---

## Step 1: Data Cleaning

**Goal:** Remove incomplete (less than one full day) records and convert raw miniSEED files to the naming format required by NoisePy.

**Script:** `scripts/clean_mseed_for_noisepy.py`

**Input:**
- Raw miniSEED files organized by station subdirectory:
  ```
  RawData/
  ├── CF.17982/
  │   ├── 230901000000    # Filename is a timestamp: YYMMDDHHmmss
  │   ├── 230902000000
  │   └── ...
  ├── CF.58808/
  └── CN.20496/
  ```

**Command:**

```bash
python scripts/clean_mseed_for_noisepy.py \
    --input-root  TestData/RawData \
    --output-root TestData/CleanData \
    --time-zone beijing \
    --target-channel Z \
    --tolerance-seconds 60
```

**Parameters:**
| Parameter | Description |
|-----------|-------------|
| `--input-root` | Root directory of raw miniSEED files |
| `--output-root` | Output directory for cleaned files |
| `--time-zone` | Timezone for filename timestamps: `beijing` (default) or `utc` |
| `--target-channel` | Target channel/component, e.g. `Z`, `BNZ`, `BHZ` |
| `--tolerance-seconds` | Tolerance (in seconds) for full-day completeness check (default: 60) |
| `--date-range START END` | Only process files whose date is strictly between START and END (both excluded), format YYMMDD, e.g. `--date-range 260316 260429` |
| `--dry-run` | Print actions without writing any files |

**Output:**
```
CleanData/
├── CF17982_BNZ00_2023258.ms
├── CF17982_BNZ00_2023259.ms
└── ...
```
Filename format: `{network}{station}_{channel}{location}_{year}{day-of-year}.ms`

---

## Step 2: Prepare Station Coordinate File

**Goal:** Create the station catalog file required by NoisePy.

**Format:** CSV file with the following columns:

```csv
network,station,channel,latitude,longitude,elevation
CF,17982,BNZ,42.323349,118.546867,844
CF,58808,BNZ,42.321987,118.54599,822
CN,20496,BNZ,42.318123,118.543602,786
```

**Notes:**
- Station names and network codes must match the miniSEED file headers
- Coordinates use decimal degrees
- Elevation is in meters
- Save as `TestData/station.csv`

---

## Step 3: Cross-Correlation and Stacking

**Goal:** Compute ambient noise cross-correlation functions between station pairs and stack them using NoisePy.

**Script:** `scripts/local_cc_workflow.py`

### 3a. Create Configuration File

Copy and edit the example config:

```bash
cp configs/local_cc_workflow.example.yaml configs/my_config.yaml
```

**Key configuration parameters (`configs/local_cc_workflow.example.yaml`):**

```yaml
# Data paths
data-dir: ../TestData/CleanData      # Output from Step 1
station-csv: ../TestData/station.csv  # Station file from Step 2
ccf-dir: ../TestData/CCFData          # Cross-correlation output
stack-dir: ../TestData/StackData      # Stacking output

# Processing parameters
freqmin: 0.2          # Minimum frequency (Hz)
freqmax: 10           # Maximum frequency (Hz)
sampling-rate: 25.0   # Target sampling rate (Hz)
maxlag: 100           # Maximum lag time (seconds)
inc-hours: 24         # Processing chunk length (hours)
cc-len: 600           # Correlation window length (seconds)
step: 300.0           # Window sliding step (seconds)

# Normalization
time-norm: rma         # Time-domain normalization
freq-norm: rma         # Frequency-domain normalization

# Stacking
stack-method: all      # Stacking method: linear / pws / robust / all
```

### 3b. Run Cross-Correlation and Stacking

```bash
# Validate configuration (dry run)
python scripts/local_cc_workflow.py --config configs/my_config.yaml --validate-only

# Run cross-correlation only
python scripts/local_cc_workflow.py --config configs/my_config.yaml --stage cc

# Run stacking only
python scripts/local_cc_workflow.py --config configs/my_config.yaml --stage stack

# Run both stages
python scripts/local_cc_workflow.py --config configs/my_config.yaml
```

### 3c. Visualize Stacked Results (Optional)

**Script:** `scripts/stack_visualization.py`

```bash
python scripts/stack_visualization.py \
    --stack-dir TestData/StackData \
    --freqmin 1 --freqmax 5 \
    --component ZZ \
    --stack-method linear \
    --maxlag 10 \
    --style both
```

Or use a YAML config file (see `configs/stack_visualization.example.yaml`):

```bash
python scripts/stack_visualization.py --config configs/my_stack_viz_config.yaml
```

**Output:** Moveout plots saved to `--outdir` (default: `TestData/figures/`).

**Publication-quality multi-panel figure (optional):** `scripts/plot_stack_panels.py` produces a 2×N multi-panel moveout figure (N = number of stack methods; top row image, bottom row wiggle/symmetric), suitable for journal figures:

```bash
python scripts/plot_stack_panels.py \
    --stack-dir TestData/StackData \
    --methods linear pws robust \
    --component ZZ
```

A YAML config file can also be used (see `configs/plot_stack_panels.example.yaml`).

---

## Step 4: Export Data for EGFAnalysisPy

**Goal:** Convert NoisePy stacking results (HDF5 format) to the text format required by EGFAnalysisPy (CF `.dat` waveforms).

**Script:** `scripts/export_egfanalysis_input.py`

```bash
python scripts/export_egfanalysis_input.py \
    --stack-dir  TestData/StackData \
    --stack-method all \
    --component ZZ \
    --output-root TestData/EGFAnalysisInputData
```

**Parameters:**
| Parameter | Description |
|-----------|-------------|
| `--stack-dir` | NoisePy stacking result directory |
| `--stack-method` | Stacking method: `linear` / `pws` / `robust` / `all` (default) |
| `--component` | Component, default `ZZ` |
| `--station` | Optional: only export pairs involving this station |
| `--output-root` | Output root directory |

**Output directory structure:**
```
EGFAnalysisInputData/
├── Allstack_linear/ZZ/
│   ├── CF17982_CF58808.dat    # Waveform (3 columns: time, gab, gba)
│   ├── CF17982_CN20496.dat
│   ├── stations.txt           # Station coordinates
│   └── EGFfiles.txt           # Waveform file list
├── Allstack_pws/ZZ/
└── Allstack_robust/ZZ/
```

Each `.dat` waveform file format:
```
  source_longitude      source_latitude      source_elevation      (line 1)
  receiver_longitude    receiver_latitude    receiver_elevation    (line 2)
  time(s)               gab                  gba                   (line 3+)
```

---

## Step 5: Dispersion Analysis & Curve Extraction with EGFAnalysisPy

**Goal:** Compute dispersion energy images from the empirical Green's / cross-correlation functions and automatically pick group and phase velocity dispersion curves with the bundled DisperPicker (CNN), followed by quality control. **This all happens in a single Python process — no MATLAB required.**

**Tool:** `EGFAnalysisPy/run.py`

Internally, `run.py` does: dispersion energy image computation → DisperPicker picking → time-variable filtering (windowed with the picked group velocity) and phase-velocity recomputation → re-picking → quality control → dispersion curve output. The input is the CF `.dat` waveforms exported in Step 4 (identical format to the original EGFAnalysisTimeFreq input).

### 5a. Configure

```bash
cd EGFAnalysisPy
cp config_test.py config.py        # start from the example, or edit config.py directly
```

Edit the key parameters in `config.py`:

```python
self.waveform_path = self.root + '/FeidongCFs'        # ← point to the Step 4 .dat directory
self.result_path   = self.root + '/DisperPicker/result'  # result output directory
self.isEGF  = False        # False for CF, True for EGF
self.StartT = 0.2          # start period (s)
self.EndT   = 5.0          # end period (s)
self.DeltaT = 0.1          # period step (s)
self.StartV = 0.5          # start velocity (km/s)
self.EndV   = 4.0          # end velocity (km/s)
self.DeltaV = 0.005        # velocity step (fixed at 0.005 km/s)
```

> Step 6A (DSurfTomo) reuses the period range and step set here, so note them down. `DeltaV` must be 0.005 km/s (required by DisperPicker).

### 5b. Run

```bash
conda activate egfanalysispy
python run.py
```

The program iterates over every `.dat` file under `waveform_path`, processing each pair. The pretrained model ships with the package (`EGFAnalysisPy/DisperPicker/saver/`) — no extra download needed.

**Parallel version (recommended for large station-pair counts):** `run_parallel.py` performs exactly the same per-file work as `run.py`, but distributes files across worker processes:

```bash
python run_parallel.py                       # min(16, CPU count) workers by default
NWORKERS=32 python run_parallel.py           # set the worker count
python run_parallel.py --config Line2Config.py -w 8   # different config, explicit workers
python run_parallel.py -n 20                 # first 20 files only, for parameter tuning
```

| Option | Env var | Description |
|--------|---------|-------------|
| `-w, --workers` | `NWORKERS` | Worker processes; default `min(16, CPU count)` |
| `-c, --config` | `CONFIG` | Config file path; defaults to `config.py` alongside the script |
| `-n, --max-files` | `MAXFILES` | Process only the first N files; 0 means all |

`--config` accepts any file with the same structure as `config.py`, so multiple profiles can each keep their own config (e.g. `Line1Config.py`, `Line2Config.py`) and be run independently.

> Because TensorFlow 1.x sessions are not picklable, each worker builds its own `Pick()` instance in the pool initializer and reuses it; workers start with `spawn` to avoid the TF fork deadlock; and native thread pools are pinned to 1 per process so parallelism comes from processes rather than oversubscribed BLAS threads. The script handles all of this — normal use requires no attention to it.

### 5c. Output

```
EGFAnalysisPy/DisperPicker/result/
├── qc_result/
│   ├── CDisp.FD01_FD16.txt    # phase velocity dispersion curve
│   ├── GDisp.FD01_FD16.txt    # group velocity dispersion curve
│   └── ...
├── qc_plot/                   # per-pair energy image + picked curve figures
├── all_C.png / mean_C.png     # all / mean phase velocity summary
└── all_G.png / mean_G.png     # all / mean group velocity summary
```

Each dispersion curve file (first two rows are station coordinates; data columns are `period velocity snr probability`):
```
117.43333100    31.81069500       # station A  longitude  latitude
117.48089000    31.66345700       # station B  longitude  latitude
0.20   1.250   6.100   0.950      # period(s)  velocity(km/s)  snr  DisperPicker probability
0.30   1.310   5.800   0.880
...
```

> This output feeds directly into the Step 6A `extractnew.py` `egfpy` input; station coordinates are embedded in the file header, so no separate station file is needed.

---

## Step 6A: DSurfTomo Surface Wave Tomography Inversion

> **One of two parallel inversion paths**, for **areal arrays**. If your stations lie along a line, skip to [Step 6B](#step-6b-1d-inversion-and-pseudo-2d-section-linear-arrays).

**Goal:** Invert dispersion curve data for a 3D shear-wave velocity model.

**Scripts location:** `DSurfTomo/scripts/`

### Option A: All-in-One Interactive Workflow (Recommended)

Use `PrepareDSurfTomoProject.py` to complete data preparation, initial model, and configuration in one step:

```bash
cd DSurfTomo/scripts
python PrepareDSurfTomoProject.py
```

The script interactively guides you through:
1. Entering the Step 5 dispersion curve directory (`EGFAnalysisPy/DisperPicker/result/qc_result`), auto-detected as `egfpy` input
2. Setting period and velocity parameters
3. Generating `surfdata.dat` (dispersion data format conversion)
4. Generating `MOD` (initial velocity model)
5. Generating `DSurfTomo.in` (inversion configuration file)
6. Optional: generating `MOD.true` (checkerboard test true model)

### Option B: Step-by-Step Execution

#### 6a. Convert Dispersion Curves → surfdata.dat

**Script:** `extractnew.py` (interactive — answer the prompts)

```bash
python extractnew.py
# input directory:  /path/to/EGFAnalysisPy/DisperPicker/result/qc_result
# (auto-detected: Detected input source: egfpy)
# curve type (phase/group/both): both
# strict period step (s): 0.1
# wave type (1=Love, 2=Rayleigh): 2
# output surfdata path: surfdata.dat
```

`extractnew.py` reads `CDisp.*.txt` and `GDisp.*.txt` under `qc_result/`; station coordinates are taken from the file headers, so no separate station file is needed. The period step must match the `DeltaT` used in Step 5.

**Output:** `surfdata.dat` + `surfdata.json` (metadata)

#### 6b. Generate Initial Model → MOD

**Script:** `GenerateIniMOD.py`

```bash
python GenerateIniMOD.py
```

Interactive prompts for:
- Depth grid nodes (e.g. `0 0.2 0.5 1 2 3 5`)
- Horizontal grid spacing or grid point count
- Initial velocity per layer (defaults computed from dispersion data)

**Output:** `MOD` (initial model file)

#### 6c. Generate Configuration → DSurfTomo.in

**Script:** `GenerateDSurfTomoInputFile.py`

```bash
python GenerateDSurfTomoInputFile.py
```

Interactive prompts for:
- Inversion parameters (iterations, damping, etc.)
- Velocity constraint range
- Smoothing parameters

**Output:** `DSurfTomo.in`

#### 6d. Optional: Generate Checkerboard Test Model

**Script:** `GenerateTrueMOD.py`

```bash
python GenerateTrueMOD.py
```

**Output:** `MOD.true`

#### 6e. Run DSurfTomo Inversion

```bash
/path/to/DSurfTomo/bin/DSurfTomo
```

The program reads `DSurfTomo.in`, `surfdata.dat`, and `MOD` from the current directory and runs the inversion.

**Output files (partial list):**
- `MOD_iter_*.dat` — velocity model at each iteration
- `residual.dat` — residuals
- `raypath_*.dat` — ray paths

### 6f. Visualize Results (GMT Scripts)

```bash
# Depth slice maps
bash plotdepth_subplot.gmt --project-dir /path/to/project

# Cross-section profiles
bash plotcross_subplot.gmt --project-dir /path/to/project \
    --profiles "lon1/lat1/lon2/lat2/label"

# Ray path maps
bash plotray_subplot.gmt --project-dir /path/to/project

# Checkerboard test results
bash plotchecker_subplot.gmt --project-dir /path/to/project
```

### 6g. Dispersion Data Statistics (Optional)

```bash
python plot_dispersion_period_counts.py --surfdata /path/to/surfdata.dat
```

---

## Step 6B: 1D Inversion and Pseudo-2D Section (Linear Arrays)

> **One of two parallel inversion paths**, for the **shallow structure beneath a linear (profile) array**. If your stations are distributed areally, use [Step 6A](#step-6a-dsurftomo-surface-wave-tomography-inversion) instead.

**Goal:** In a linear array the ray paths are nearly collinear and cannot constrain a 2D velocity grid, so DSurfTomo's tomographic inversion does not apply. This path instead **inverts each station pair's dispersion curve independently in 1D**, yielding a set of 1D shear-wave velocity models, which are then stitched into a pseudo-2D section by projecting each pair's midpoint onto the profile.

**Tools:** `scripts/lant_invert.py` (1D inversion) and `scripts/lant_stitch.py` (stitching)

**Dependency:** the inversion uses [evodcinv](https://github.com/keurfonluu/evodcinv) (CPSO evolutionary algorithm)

```bash
pip install evodcinv
```

### 6B-a. Write the Config File

Both scripts share one YAML config; every parameter lives there and nothing is hard-coded in the scripts:

```bash
cp configs/lant.example.yaml configs/my_lant.yaml
```

Key parameters:

```yaml
paths:
  input_dir:   ../TestData/DisperPicker/qc_result   # ← qc_result from Step 5
  station_csv: ../TestData/station.csv              # ← station file from Step 2
  output_dir:  ../TestData/Inversion                # output directory

inversion:
  curve_type: phase    # phase (CDisp) | group (GDisp) | both (joint)
  workers: 16          # parallel worker processes
  min_points: 8        # skip curves with fewer valid points
  zmax: 0.6            # maximum output depth (km)
  dz: 0.01             # depth step (km)
  layers:              # layered search space; thickness/vs given as [min, max]
    - {thickness: [0.01, 0.04], vs: [0.3, 2.5]}
    ...
```

> **How to choose `zmax`:** surface-wave depth sensitivity is roughly wavelength/3, i.e. `c*Tmax/3`. For example, a maximum period of 0.5 s with a phase velocity around 2.2 km/s constrains only about 0.4 km, so `zmax: 0.6` is appropriate. Parameterizing much deeper buys no resolution — only unconstrained values. The thickness maxima in `layers` should sum to roughly `zmax`.

### 6B-b. Per-Pair 1D Inversion

```bash
python scripts/lant_invert.py --config configs/my_lant.yaml

# Common overrides (CLI takes precedence over the config file)
python scripts/lant_invert.py --config configs/my_lant.yaml \
    --curve-type both --workers 8 --zmax 0.8 --no-plot
python scripts/lant_invert.py --config configs/my_lant.yaml -n 10   # try 10 first
```

**Output:** `output_dir/models/` (one 1D model per station pair), `output_dir/figs/` (fit comparison figures), and `output_dir/summary.txt` (per-pair misfit summary).

Start with `-n 10`, inspect the fits in `figs/` and the misfit distribution in `summary.txt`, tune the `layers` search ranges and `optimizer` settings, then run the full set.

### 6B-c. Stitch into a Pseudo-2D Section

```bash
python scripts/lant_stitch.py --config configs/my_lant.yaml

# Load models from a different output directory without editing the config
python scripts/lant_stitch.py --config configs/my_lant.yaml --output-dir /path/to/other
```

Stitching parameters live under the `stitch:` section of the config:

```yaml
stitch:
  topography: true       # follow local terrain at the top (uses station_csv elevations)
  show_stations: true    # station triangles along the surface
  misfit_max: 0.2        # drop models with misfit above this (null = keep all)
  nbins: 120             # number of along-profile distance bins
  smooth: 1.2            # gaussian smoothing sigma in cells (0 = off)
  vertical_exaggeration: 2.0
  cmap: jet_r            # low-Vs warm / high-Vs cool
  vmax: 4
  taper_slope: 0.1       # mask into an inverted trapezoid (wide top, narrow bottom):
                         # trim this many km per side per km of depth, since edges
                         # get less reliable with depth. 0 = keep the full rectangle
```

**Output:** `output_dir/section.png` (plus a PDF when `save_pdf: true`).

> The horizontal axis is each pair's midpoint projected onto the profile, so this is a **pseudo-2D, not a true 2D, section**: every column comes from an independent 1D inversion with no lateral coupling between columns. Read it as a visual arrangement of 1D results, and do not over-interpret short-wavelength lateral variation.

### 6B-d. Station Distribution Check (Optional)

```bash
python scripts/plot_station_distribution_all.py
```

Plots the map-view distribution of all stations — useful for confirming the array really is linear and for deciding between 6A and 6B.

---

## Directory Structure

```
local_ANT/
├── scripts/                            # Helper scripts (glue code between packages)
│   ├── clean_mseed_for_noisepy.py      # Step 1: Data cleaning
│   ├── plot_station_distribution_all.py # Step 2: Station distribution check
│   ├── local_cc_workflow.py            # Step 3: NoisePy cross-correlation workflow
│   ├── stack_visualization.py          # Step 3c: Stack visualization
│   ├── plot_stack_panels.py            # Step 3c: Publication-quality multi-panel moveout
│   ├── export_egfanalysis_input.py     # Step 4: Export EGF input
│   ├── lant_invert.py                  # Step 6B: Per-pair 1D inversion
│   └── lant_stitch.py                  # Step 6B: Stitch into a pseudo-2D section
├── configs/                            # Example configuration files
│   ├── local_cc_workflow.example.yaml
│   ├── plot_stack_panels.example.yaml
│   ├── stack_visualization.example.yaml
│   └── lant.example.yaml               # Step 6B: 1D inversion + stitching
├── EGFAnalysisPy/                      # Step 5: Dispersion analysis + CNN picker
│   ├── run.py                          # One-step dispersion analysis & picking
│   ├── run_parallel.py                 # Same work, multiprocessing version
│   ├── config.py / config_test.py      # Parameter configuration
│   └── DisperPicker/                   # Bundled DisperPicker with pretrained model (saver/)
├── DSurfTomo/                          # Step 6A: Surface wave tomography inversion
│   ├── bin/                            # Executables
│   ├── scripts/                        # Data preparation and plotting scripts
│   └── src/                            # Source code
└── TestData/                           # Example data
    ├── RawData/                        # Raw miniSEED
    ├── CleanData/                      # Cleaned data
    ├── CCFData/                        # Cross-correlation results
    ├── StackData/                      # Stacking results
    ├── EGFAnalysisInputData/           # EGF input (CF .dat waveforms, generated by Step 4)
    ├── figures/                        # Visualizations
    └── station.csv                     # Station information
```

## FAQ

### Should I use Step 6A (DSurfTomo) or Step 6B (1D inversion)?
It depends on the array geometry. Use **6A** when stations are distributed areally and the station-pair ray paths cross the study region from several azimuths — this yields a 3D model. Use **6B** when stations lie along a single line: all rays are then nearly collinear, the across-profile direction of a 2D grid has no ray constraint at all, and DSurfTomo would return unreliable results. Running `scripts/plot_station_distribution_all.py` first is an easy way to decide.

Note that 6B produces a **pseudo-2D** section: each column comes from an independent 1D inversion with no lateral coupling, so it must not be interpreted as a true 2D inversion result.

### Must the period and velocity parameters be consistent across steps?
Yes. The period range and step set in Step 5 (EGFAnalysisPy `config.py`) must match the period step used in Step 6A (`extractnew.py`). DisperPicker requires a fixed velocity step of 0.005 km/s (`DeltaV`).

### How do I choose a stacking method?
- `linear`: Linear stacking — most commonly used
- `pws`: Phase-weighted stacking — higher SNR
- `robust`: Robust stacking — less sensitive to outliers
- `all`: Compute all methods simultaneously

We recommend computing `all` first, then comparing results with `stack_visualization.py` to select the best.

### What input does extractnew.py take?
Only the `egfpy` input type: the `qc_result/` output of EGFAnalysisPy (`run.py`) — `CDisp.*.txt` (phase) and `GDisp.*.txt` (group). The first two rows hold the station coordinates (lon lat); data columns are `period velocity snr probability`, where the 4th column is the DisperPicker probability and is not used as a keep flag. The input directory may point to the `qc_result` directory or its parent. Station coordinates are embedded in the file headers, so no separate station file is needed.

---

## License

Each software package is distributed under its own original license. The helper scripts (`scripts/` and `DSurfTomo/scripts/`) are distributed with this project.
