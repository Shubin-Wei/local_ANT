# Complete Workflow: Ambient Noise Surface Wave Tomography from Local miniSEED Data

This project implements a complete ambient noise surface wave tomography pipeline — from raw miniSEED seismic data to a 3D shear-wave velocity model. The workflow consists of 8 steps, using 5 software packages and several helper scripts.

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
EGFAnalysisTimeFreq input format (.dat waveforms + stations.txt + EGFfiles.txt)
  │
  ▼ [Step 5] EGFAnalysisTimeFreq (MATLAB)
Dispersion spectrograms (C.*.dat / G.*.dat)
  │
  ▼ [Step 6] prepare_disperpicker_input.py
DisperPicker input format (renamed files + updated config)
  │
  ▼ [Step 7] DisperPicker (pick_v.py)
Automatically picked dispersion curves
  │
  ▼ [Step 8] DSurfTomo inversion
    ├── extractnew.py               → surfdata.dat
    ├── GenerateIniMOD.py           → MOD (initial model)
    ├── GenerateDSurfTomoInputFile.py → DSurfTomo.in
    ├── Run DSurfTomo inversion
    └── GMT scripts for plotting (depth slices, cross-sections, ray paths, checkerboard tests)
```

## Citations

This project integrates the following open-source software packages. All credit goes to their original authors:

- **NoisePy** — Ambient noise cross-correlation
  Jiang, C., & Denolle, M. A. (2020). NoisePy: A New High-Performance Python Tool for Ambient-Noise Seismology. *Seismological Research Letters*, 91(3), 1853–1866. https://doi.org/10.1785/0220190364
  GitHub: https://github.com/noisepy/NoisePy

- **EGFAnalysisTimeFreq** — Surface wave dispersion analysis (MATLAB)
  Yao, H., van der Hilst, R. D., & de Hoop, M. V. (2006). Surface-wave array tomography in SE Tibet from ambient seismic noise and two-station analysis—I. Phase velocity maps. *Geophysical Journal International*, 166(2), 732–744. https://doi.org/10.1111/j.1365-246X.2006.03028.x
  Author: Huajian Yao (University of Science and Technology of China)

- **DisperPicker** — CNN-based automatic dispersion curve extraction
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

### MATLAB

EGFAnalysisTimeFreq requires a MATLAB installation.

### DisperPicker Environment

```bash
cd DisperPicker
pip install -r requirements.txt
```

### Compile DSurfTomo

```bash
cd DSurfTomo
./configure
make
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
freqmin: 0.5          # Minimum frequency (Hz)
freqmax: 10           # Maximum frequency (Hz)
sampling-rate: 50.0   # Target sampling rate (Hz)
maxlag: 500           # Maximum lag time (seconds)
inc-hours: 24         # Processing chunk length (hours)
cc-len: 3600          # Correlation window length (seconds)
step: 1800.0          # Window sliding step (seconds)

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

---

## Step 4: Export Data for EGFAnalysisTimeFreq

**Goal:** Convert NoisePy stacking results (HDF5 format) to the text format required by EGFAnalysisTimeFreq.

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

## Step 5: Generate Dispersion Spectrograms with MATLAB

**Goal:** Perform time-frequency analysis on empirical Green's functions using EGFAnalysisTimeFreq to generate dispersion spectrograms.

**Tool:** `EGFAnalysisTimeFreq/EGFAnalysisTimeFreq.m` (MATLAB GUI)

### Procedure

1. Open MATLAB and run `EGFAnalysisTimeFreq.m`
2. Set parameters in the GUI (reference: `EGFAnalysisTimeFreq/Feidong_disp/fd_config.png`)
3. **Important parameters (you will need these in later steps):**
   - **Period range:** e.g. 0.2 – 5.0 s
   - **Velocity range:** e.g. 0.5 – 4.0 km/s
   - **Velocity step:** must be 0.005 km/s (required by DisperPicker)
   - **Period step:** recommended 0.1 s
4. Set input data directory: select the output from Step 4 (e.g. `TestData/EGFAnalysisInputData/Allstack_robust/ZZ/`)
5. Select **Automatic** mode, click **Start Processing**
6. The program processes all station pairs automatically and outputs dispersion spectrograms

**Output:**
```
EGFAnalysisOutputData/RobustDisperImage/
├── C.CF17982_CF58808.dat    # Phase velocity spectrogram
├── C.CF17982_CN20496.dat
├── G.CF17982_CF58808.dat    # Group velocity spectrogram
├── G.CF17982_CN20496.dat
└── ...
```

Each `.dat` file is a matrix (e.g. 701×49): rows = velocity, columns = period.

> **Note:** For detailed instructions, refer to `EGFAnalysisTimeFreq/EGFAnalysisTimeFreq_Manual.pdf`

---

## Step 6: Prepare DisperPicker Input

**Goal:** Convert EGFAnalysisTimeFreq dispersion spectrograms to DisperPicker's expected format and update its configuration.

**Script:** `scripts/prepare_disperpicker_input.py`

```bash
python scripts/prepare_disperpicker_input.py \
    --input-dir  TestData/EGFAnalysisOutputData/RobustDisperImage \
    --picker-dir DisperPicker \
    --start-period 0.2 \
    --end-period 5.0 \
    --min-velocity 0.5 \
    --max-velocity 4.0 \
    --station-file DisperPicker/config/station.txt \
    --clean
```

**Parameters:**
| Parameter | Description |
|-----------|-------------|
| `--input-dir` | Directory of dispersion spectrograms from EGFAnalysisTimeFreq |
| `--picker-dir` | DisperPicker root directory |
| `--start-period` | Start period (s), must match Step 5 |
| `--end-period` | End period (s), must match Step 5 |
| `--min-velocity` | Minimum velocity (km/s), must match Step 5 |
| `--max-velocity` | Maximum velocity (km/s), must match Step 5 |
| `--station-file` | Station coordinate file (optional) |
| `--clean` | Remove existing data in DisperPicker before copying |

**This script automatically:**
1. Copies `C.*.dat` → `DisperPicker/data/TestData/phase_image/STA1.STA2.dat`
2. Copies `G.*.dat` → `DisperPicker/data/TestData/group_image/STA1.STA2.dat`
3. Detects matrix dimensions
4. Updates `DisperPicker/config/config.py` — sets `range_T`, `range_V`, and `root`

---

## Step 7: Automatic Dispersion Curve Extraction with DisperPicker

**Goal:** Use a pre-trained CNN model to automatically extract group and phase velocity dispersion curves from spectrograms.

### 7a. Verify Configuration

Confirm that `DisperPicker/config/config.py` was correctly updated by Step 6:

```python
self.root = '/your/path/to/DisperPicker'
self.range_T = [0.2, 5.0, 49]    # [start_period, end_period, n_periods]
self.range_V = [0.5, 4.0, 701]   # [start_velocity, end_velocity, n_velocities]
```

Confirm that `DisperPicker/config/station.txt` has correct station coordinates:

```
station_name  longitude     latitude       elevation
FD01          117.4333135   31.81069491    12.939
```

### 7b. Run Extraction

```bash
cd DisperPicker
python pick_v.py
```

### 7c. Quality Control (QC)

```bash
cd DisperPicker/result/process
python process_C.py    # Phase velocity QC
python process_G.py    # Group velocity QC
python plot_res.py     # Plot each curve overlaid on its spectrogram
```

**Output:**
```
DisperPicker/result/
├── pick_result/           # Raw picking results
├── process/
│   └── new/               # QC'd dispersion curves
│       ├── phase_velocity/
│       │   ├── STA1.STA2.txt
│       │   └── ...
│       ├── group_velocity/
│       │   ├── STA1.STA2.txt
│       │   └── ...
│       ├── DPN_C.jpg      # All phase velocity curves summary
│       └── DPN_G.jpg      # All group velocity curves summary
```

Each dispersion curve file (e.g. `STA1.STA2.txt`) format:
```
period(s)  velocity(km/s)
0.3        2.15
0.4        2.28
...
```

---

## Step 8: DSurfTomo Surface Wave Tomography Inversion

**Goal:** Invert dispersion curve data for a 3D shear-wave velocity model.

**Scripts location:** `DSurfTomo/scripts/`

### Option A: All-in-One Interactive Workflow (Recommended)

Use `PrepareDSurfTomoProject.py` to complete data preparation, initial model, and configuration in one step:

```bash
cd DSurfTomo/scripts
python PrepareDSurfTomoProject.py
```

The script interactively guides you through:
1. Selecting input type (DisperPicker results or EGFAnalysisTimeFreq manual picks)
2. Setting period and velocity parameters
3. Generating `surfdata.dat` (dispersion data format conversion)
4. Generating `MOD` (initial velocity model)
5. Generating `DSurfTomo.in` (inversion configuration file)
6. Optional: generating `MOD.true` (checkerboard test true model)

### Option B: Step-by-Step Execution

#### 8a. Convert Dispersion Curves → surfdata.dat

**Script:** `extractnew.py`

```bash
python extractnew.py \
    --input-dir /path/to/DisperPicker/result/process/new \
    --input-type disperpicker \
    --station-file /path/to/station.txt \
    --output-dir /path/to/project \
    --min-period 0.2 --max-period 5.0 --period-step 0.1 \
    --curve-type both
```

Also supports EGFAnalysisTimeFreq manual picks:

```bash
python extractnew.py \
    --input-dir /path/to/manual_picks \
    --input-type manual \
    --output-dir /path/to/project \
    --min-period 0.2 --max-period 5.0 --period-step 0.1
```

**Output:** `surfdata.dat` + `surfdata.json` (metadata)

#### 8b. Generate Initial Model → MOD

**Script:** `GenerateIniMOD.py`

```bash
python GenerateIniMOD.py
```

Interactive prompts for:
- Depth grid nodes (e.g. `0 0.2 0.5 1 2 3 5`)
- Horizontal grid spacing or grid point count
- Initial velocity per layer (defaults computed from dispersion data)

**Output:** `MOD` (initial model file)

#### 8c. Generate Configuration → DSurfTomo.in

**Script:** `GenerateDSurfTomoInputFile.py`

```bash
python GenerateDSurfTomoInputFile.py
```

Interactive prompts for:
- Inversion parameters (iterations, damping, etc.)
- Velocity constraint range
- Smoothing parameters

**Output:** `DSurfTomo.in`

#### 8d. Optional: Generate Checkerboard Test Model

**Script:** `GenerateTrueMOD.py`

```bash
python GenerateTrueMOD.py
```

**Output:** `MOD.true`

#### 8e. Run DSurfTomo Inversion

```bash
/path/to/DSurfTomo/bin/DSurfTomo
```

The program reads `DSurfTomo.in`, `surfdata.dat`, and `MOD` from the current directory and runs the inversion.

**Output files (partial list):**
- `MOD_iter_*.dat` — velocity model at each iteration
- `residual.dat` — residuals
- `raypath_*.dat` — ray paths

### 8f. Visualize Results (GMT Scripts)

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

### 8g. Dispersion Data Statistics (Optional)

```bash
python plot_dispersion_period_counts.py --surfdata /path/to/surfdata.dat
```

---

## Directory Structure

```
local_ANT/
├── scripts/                            # Helper scripts (glue code between packages)
│   ├── clean_mseed_for_noisepy.py      # Step 1: Data cleaning
│   ├── local_cc_workflow.py            # Step 3: NoisePy cross-correlation workflow
│   ├── stack_visualization.py          # Step 3c: Stack visualization
│   ├── export_egfanalysis_input.py     # Step 4: Export EGF input
│   └── prepare_disperpicker_input.py   # Step 6: Prepare DisperPicker input
├── configs/                            # Example configuration files
│   ├── local_cc_workflow.example.yaml
│   └── stack_visualization.example.yaml
├── EGFAnalysisTimeFreq/                # MATLAB dispersion analysis program
├── DisperPicker/                       # Deep learning dispersion curve picker
│   ├── config/                         # Configuration files
│   ├── data/TestData/                  # Input data
│   └── result/                         # Output results
├── DSurfTomo/                          # Surface wave tomography inversion
│   ├── bin/                            # Executables
│   ├── scripts/                        # Data preparation and plotting scripts
│   └── src/                            # Source code
└── TestData/                           # Example data
    ├── RawData/                        # Raw miniSEED
    ├── CleanData/                      # Cleaned data
    ├── CCFData/                        # Cross-correlation results
    ├── StackData/                      # Stacking results
    ├── EGFAnalysisInputData/           # EGF input
    ├── EGFAnalysisOutputData/          # Dispersion spectrograms
    ├── figures/                        # Visualizations
    └── station.csv                     # Station information
```

## FAQ

### Must the period and velocity parameters be consistent across steps?
Yes. The period range and velocity range set in Step 5 (EGFAnalysisTimeFreq) must be passed with the same values in Step 6 (prepare_disperpicker_input.py). DisperPicker requires a fixed velocity step of 0.005 km/s.

### How do I choose a stacking method?
- `linear`: Linear stacking — most commonly used
- `pws`: Phase-weighted stacking — higher SNR
- `robust`: Robust stacking — less sensitive to outliers
- `all`: Compute all methods simultaneously

We recommend computing `all` first, then comparing results with `stack_visualization.py` to select the best.

### What is the station.txt format for DisperPicker?
```
station_name  longitude     latitude       elevation
FD01          117.4333135   31.81069491    12.939
```

### What input types does extractnew.py support?
- `disperpicker`: QC'd results from DisperPicker (`.txt` files under `phase_velocity/` and `group_velocity/`)
- `manual`: Manual picks from EGFAnalysisTimeFreq (`CDisp.T.*.dat` and `GDisp.*.dat` files)

---

## License

Each software package is distributed under its own original license. The helper scripts (`scripts/` and `DSurfTomo/scripts/`) are distributed with this project.
