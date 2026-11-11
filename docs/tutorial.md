# Tutorial: from raw records to a 3-D shear-velocity model

[中文](tutorial.zh.md) | English

Purpose, method, basis and what to look at, step by step. A first run in half an hour: [quickstart.md](quickstart.md). Principles and formulas: [methods.md](methods.md); defaults and what they rest on: [parameters.md](parameters.md); every output file: [outputs.md](outputs.md); problems: [troubleshooting.md](troubleshooting.md); running with an AI assistant: [ai_workflow.md](ai_workflow.md); numbers of the synthetic tests: [validation.md](validation.md).

## 1. Installation

### 1.1 Python environment

```bash
git clone https://github.com/Shubin-Wei/local_ANT.git && cd local_ANT
conda env create -f environment.yml && conda activate lant     # Python, ObsPy, GMT 6, gfortran
pip install -e .
```

If `conda env create` stays at "Solving environment" for a long time, use `--solver libmamba`, or mamba / micromamba. Without conda: `pip install -e .` on Python ≥ 3.9, plus GMT 6 (≥ 6.1), gfortran, gcc, make, git. Linux only.

### 1.2 Third-party programs

The repository holds no third-party source; a script downloads fixed versions from the official repositories and builds them in `external/`:

```bash
bash tools/install_external.sh            # DSurfTomo + fwd1d, pDSurfTomo
bash tools/install_external.sh dsurftomo  # DSurfTomo only
```

| Program | Used for | Required |
|---|---|---|
| DSurfTomo (Fang et al. 2015, MIT) | 3-D inversion; `fwd1d` is its 1-D forward calculation, used for the starting model | yes |
| pDSurfTomo (Zhu et al. 2026, GPL-3) | parallel version for large data sets; at most 80 periods | no |

Neither program is modified, except for one patch that lowers the memory reserve of pDSurfTomo. The frequency-Bessel spectrum is computed by `lant` itself; no further program is needed.

### 1.3 Check

```bash
lant check          # state of the Python packages, GMT, DSurfTomo, fwd1d, and how to install what is missing
```

### 1.4 Run the synthetic example first

```bash
python examples/synthetic/make_synthetic.py              # 16 stations, 3 days; 5–20 minutes, about 300 MB
                                                         # (--quick: 10 stations, 6 hours, about 1 minute; quickstart.md)
lant run all -c examples/synthetic/synthetic.yaml        # stops where a judgement is yours and prints the next command
python examples/synthetic/compare_truth.py               # comparison with the known answer
```

The data are a sum of plane Rayleigh waves of a known 1-D model arriving from random azimuths; miniSEED and SAC files are deliberately mixed in nested folders and the station table has no header. `--layout line` makes a line array, `--model sediment` a slow sedimentary model. These data are laterally uniform, isotropic and fundamental-mode only; a laterally varying medium, a higher mode and added noise are made directly as correlations with `make_stack.py` (see validation.md).

## 2. Preparing the data

You need two things: the raw-data directory and a station table. No sorting into folders, renaming, format conversion, cutting into days or merging.

### 2.1 Raw records

Any depth of sub-folders, any file names, any format ObsPy reads, also mixed; files that are not data are skipped and listed in the report. Formats covered by tests: miniSEED, SAC, GSE2, SEG-Y (files written by ObsPy); SEG-2 goes through the same reading path but has not been tested (ObsPy cannot write it, and no long passive recording was at hand). Records without a station code (SEG-Y, SEG-2, some SAC files) are named from the file or folder name; a file that holds several such records at the same time (one trace per channel) cannot be named that way and is reported as not used. From the headers the program takes network, station, channel, sampling rate and times, and then:

| What | How |
|---|---|
| Vertical component | channel code ending in Z, or SAC header cmpinc = 0; a station with a single channel uses it, and the report says so; write `data.channel` when the names do not tell |
| Several sensors or rates at one station | the one with most data is used, the others are listed |
| Start of the processing day | the full UTC hour at which most long records start (`data.day_start_utc_hour`) |
| Duplicates, overlaps, sample-grid offsets, gaps | handled and counted; gaps are not filled |
| No station name in the header | the part of the file name or of the folder names that matches the station table |

### 2.2 Station table

CSV, text separated by blanks or tabs (with or without header), Excel, StationXML. Column names are recognised in English and Chinese. Required: station name, longitude, latitude (decimal degrees or degrees-minutes-seconds; no projected coordinates). Optional columns:

| Column | Used for |
|---|---|
| elevation | sections |
| network | matching the headers |
| instrument | without response files: confirms that all stations are one model |
| serial, start, end | moved nodes: when one instrument stood at several sites, its records are given to the sites by serial number and time span. Times may be dates only (e.g. `260913` or `2026-09-13`): a whole day in the time zone of `data.day_label`. Without a serial column the station name is the instrument code |
| line | several lines in one table: `data.line` picks one; `lant run all --raw …` makes one project per line |

Names are matched by: network and station equal; station equal; "network + station" written together; equal after dropping leading letters or zero padding. Ambiguous names are not used and are listed. The result is `work/01_clean/stations.csv`; from then on stations are named `NET.STA`.

**Check the coordinates station by station**: one wrong coordinate makes every travel time of that station wrong, without any error message. Where SAC headers hold coordinates, stations more than 50 m from the table are listed. Arrays across the date line, beyond 80° latitude, or with two stations closer than 1 m are refused or pointed out.

### 2.3 Instrument response

- With response files (StationXML, dataless SEED, RESP): set `data.response`; the clean step removes the response to velocity (Bensen et al. 2007), with a pre-filter derived from the correlation band.
- Without: all stations must be one instrument model. Cross-coherence removes the amplitude response; the phase responses of one model largely cancel between the two stations, but sensors of one model still differ within their tolerance, which leaves random phase errors of the order of 0.1–1 % of the travel time. If the table names several models and there are no response files the program stops with an error; if there is no model information you are asked to confirm after clean.

### 2.4 Disk and computing

- The cleaned day files are about as large as the raw data (`clean.write_archive: false` skips them; response removal needs them).
- `02_ccf` takes about pairs × 24 × (2 × lag window × rate + 1) × 4 bytes per day; the clean report gives this estimate and the memory of one cc process. When memory does not hold `correlation.nproc` processes (70 % of the physical memory), the cc step runs fewer; when not even one fits it stops and suggests a narrower `data.velocity_range_kms` or `correlation.max_distance_km` (only pairs within that distance).
- The inversion time depends on the number of grid nodes and of "stations × periods"; every checkerboard and control is one project, and they can run side by side (`inversion.njobs`).

## 3. Configuration

The shortest way needs no configuration file:

```bash
lant run all --raw /data/raw --stations /data/stations.csv --out /data/result
```

It writes `/data/result/config.yaml` (the three paths and common settings; use `lant init --full` for all defaults and their basis) and starts. To look first:

```bash
lant init -o my.yaml --raw /data/raw --stations /data/stations.csv --out /data/result
lant validate -c my.yaml        # checks the configuration and the station table
lant plan -c my.yaml            # what every step will do, with which values, and where they come from
```

Keys worth a look at the start: `data.velocity_range_kms` (default 0.1–5 km/s, from soft soil to hard rock; a narrower range you know makes the lag window and the F-J scan smaller and the run faster), `data.response`, `data.day_label` (default UTC). Everything else has a default or is derived from the data.

## 4. Running

```bash
lant run all -c CONFIG          # from the first step that does not stand; stops where a judgement is yours
lant status -c CONFIG           # state of every step, why, and the next command (--json for programs)
lant confirm STEP -c CONFIG     # after you looked at the report and the figures
lant run all -c CONFIG --from disp     # run again from a step after changing a parameter
```

Every step prints a report and stores it in `.lant/steps/<step>.md`: what it did, key numbers, every parameter it used (value, basis, source), files to look at, items left for you. `.lant/manifest.json` records on which inputs and parameters every step stands; when one of them changes, the step and those after it become stale and confirmations lapse. A step with nothing to judge does not stop.

### 4.1 clean: cleaning (`01_clean`)

**What**: reads all headers; matches the station table to the records; chooses the vertical channel; cuts processing days; puts the samples on the UTC sample grid (offsets below 0.05 sample are shifted, others interpolated with a Lanczos kernel); removes duplicates and identical overlaps; removes the response when response files are given; tests every station-day for a coherent pulse once per second with an F test on the record folded at 1 s (Takagi et al. 2015); judges the array type from the station distribution. Apart from response removal: no filtering, no normalisation. A second run reuses every station-day whose input files and settings did not change.

**Look at**: `report.txt` (unmatched names, stations without data, resource estimate for cc), `coverage.png` (data per station and day), `stations.csv`, `pulse_1s.csv` (stations with the pulse).

**Left to you**: ambiguous names, coordinate mismatches, unknown component orientation, an unclear array type (elongation 3–10), responses not removed and models unknown.

### 4.2 cc: hourly cross-coherence (`02_ccf`)

**What**: every station-day is decimated to the processing rate and cut into segments with 50 % overlap (Seats et al. 2012); segments whose peak exceeds `max_over_std` times the standard deviation of the station-day are dropped (Jiang & Denolle 2020); each segment gives the cross-coherence (each spectrum divided by its own amplitude, Xu et al. 2021), averaged per hour. Records of the neighbouring days are read too, so segments across a day boundary are kept. For stations where clean detected the 1-s pulse, the folded waveform is subtracted hour by hour first; other stations are left alone.

**Parameters**: sampling rate, segment length, lag window and band follow from the geometry and the velocity range by written rules (the report shows each derivation with your numbers).

### 4.3 stack: stack and station diagnosis (`03_stack`)

**What**: all hours stacked, weighted by their segments; sub-stacks of several hours are kept for the jackknife (`stack.block_hours`). Diagnosis: every sub-stack of a pair is compared with the stack of the pair's other sub-stacks (leave one out) by correlation coefficient and sub-sample time shift (Stehly et al. 2007; Sens-Schönfelder 2008); per station the median over its pairs with the nearest stations. Candidates are the outliers of the two distributions: the correlation enters as Fisher z, the scale is the median and 1.4826 × MAD (Leys et al. 2013), and the limit follows from the sample size by Chauvenet's criterion, not from a fixed number. A time shift is listed only when the shifts of the pairs of the station agree, i.e. the median shift is beyond the same limit in units of its standard error over the pairs: a clock error moves every pair of a station alike, while in sub-stacks of low correlation the shifts scatter.

**The diagnosis lists, it does not remove.** `candidates.txt` and `station_qc/station_qc.png` are for you; when there are more candidates than the criterion lets through by chance (1), the chain stops and you decide: `lant confirm stack` (keep) or `lant confirm stack --exclude-candidates` (remove), or write your own list in `stack.exclude_station_days`.

### 4.4 reference: reference curve (`04_reference`)

**What**: the frequency-Bessel spectrum of the array (Wang et al. 2019; double precision, per frequency only the pairs within a few wavelengths), the ridge picked and fitted by a smooth curve in log f: the reference phase-velocity curve. From it follow the lower end of the measurement band and the velocity axis. The reference only fixes the cycle count of the phase and guides the group-arrival pick; it does not enter the final values.

**Look at**: `fj.png`. **This is the most important judgement of the chain: is the picked ridge the fundamental mode?** Where a higher mode is stronger, `max_peak` picks it; then use `track_fundamental` (it follows the ridge from the start, tried every 1/3 octave, that can be followed over the most frequencies, within `data.velocity_range_kms` and, on a line, only at wavelengths of at least twice the mean station spacing, drawn as a dashed line in `fj.png`), narrow `data.velocity_range_kms`, or give a curve with `reference.source: file`. The report also compares the measured velocities with `data.velocity_range_kms`.

**Large arrays with strong lateral variation**: `reference.zone_km` divides the array into square zones by path midpoint, one curve per zone (`zones.png`); zones without a clear ridge use the curve of the whole array.

### 4.5 disp: dispersion measurement, quality control, export (`05_dispersion`)

**What**:
1. The symmetric component of every pair is filtered with Gaussian narrow-band filters (as FTAN, Bensen et al. 2007) at log-spaced frequencies (12 per octave, Xu et al. 2021); the group arrival is picked along a ridge guided by the reference.
2. The phase is read with a phase-matched filter (Herrin & Goforth 1977; Levshin & Ritzwoller 2001): the measured group-arrival curve gives a phase correction that removes the phase variation inside the envelope caused by dispersion, and the phase is then read at the centre frequency; far-field phase constant π/4 (Lin et al. 2008). Travel time t = −(Φ − π/4)/ω plus whole cycles, phase velocity c = distance / t (two-station method, Yao et al. 2006), with distances on the WGS84 ellipsoid (Vincenty 1975).
3. The cycle count is fixed where the path is 1–3.5 reference wavelengths (there the reference cannot be a cycle off) and then followed continuously in frequency.
4. Errors: delete-a-group jackknife over groups of sub-stacks, each measured again. The errors are not calibrated; they show repeatability in time only.
5. Quality control: every rejected sample carries its reasons (bit mask), counted in the report. Criteria include signal-to-noise ratio, jackknife error, a path of at least 2 of the sample's own wavelengths (Luo et al. 2015), run length, agreement of group and phase velocity, and the band edge (frequencies where the amplitude spectrum of the coherent signal, i.e. of the signal window minus the incoherent noise of late lags, over all pairs, is below half its mean inside the filter band are not measured); the basis of every threshold is in parameters.md, and those labelled "empirical" have no literature value.
6. Export: per pair only the longest continuous segment, written as DSurfTomo input.

**Look at**: `results/figures/3_dispersion_curves.png` (all curves and the reference), `pair_figures/` (single pairs: velocity spectrum, picks, reasons), `median_vs_reference.txt`.

**Do not loosen the quality control to get more points.** A changed threshold needs its reason in `.lant/decisions.md`.

### 4.6 inversion: building the projects (`06_inversion`)

**What**: 1-D inversion of the array's average curve for the starting model (damped least squares, smoothness applied to the model itself, Constable et al. 1987); depth nodes in steps of max(shortest wavelength / 4, depth / 4) down to half the longest wavelength (rule of thumb, Foti et al. 2018); grid spacing about half the median nearest-station distance; the main project, checkerboard projects (Lévêque et al. 1993) and by default one control project with another starting model (`start_grad0.05`; regularisation controls can be added in `checks.controls`; Rawlinson & Spakman 2016). The sparse-matrix reserve of DSurfTomo is estimated from the paths (`inversion.sparsity: auto`). DSurfTomo takes at most 60 periods (pDSurfTomo 80); with more, the step stops and names the `dispersion.per_octave` that fits.

**Look at**: `results/figures/4_start_model.png`, `4_kernels.png` (depth sensitivity of the periods). **Left to you**: whether starting model, deepest node and velocity bounds are reasonable.

### 4.7 run, evaluate, figures

- **run**: DSurfTomo (Fang et al. 2015: direct inversion of the travel times for 3-D Vs, rays updated with the model). Before it runs, the share of data its rejection rule keeps in the first iteration is checked (`checks.precheck`); an invalid main project stops the chain. The main project and the controls run side by side; the checkerboards follow when the main project has finished, since their noise is the relative misfit of the main model (`inversion.noise: auto`; the value is in the run report and in `work/06_inversion/checkerboard_noise.json`). A number in `inversion.noise` lets all projects run together.
- **evaluate**: numbers and figures only, no verdicts. Travel-time residuals (starting and final model, share kept), checkerboard recovery (correlation and regression slope by depth), differences between control models and the main model together with the residual history of each control (a control whose residual has not come down to that of the main project is not converged). Written to `work/06_inversion/evaluation/`.
- **figures**: in `results/figures/` (list in outputs.md): station map over relief, correlation moveouts per band, three single pairs (stack time and dispersion), dispersion curves, starting model and kernels, depth slices and a 3-D view (areal arrays), sections, checkerboards, residuals. Scales, distances, lag windows and velocity ranges follow the data of the project; the relief is SRTM, loaded automatically from the GMT data server (`plots.topography`). The model of a line is a narrow strip along it: no depth slices or 3-D view, and the checkerboards are drawn as input and recovered sections along the line. A section of your own that leaves the model is cut to the model, and the note says by how much.

## 5. Line arrays

One project per line. The defaults that differ are listed at the end of parameters.md.

1. The F-J spectrum of a line has pairs in one direction only and is often affected by higher modes or obliquely arriving noise; the default is `track_fundamental`; do look at `fj.png`. Where an areal array of the same area exists, its reference curve and starting model can be borrowed (`reference.source: file`, `inversion.start_model.source: file`); agreement between the line and the areal result is then no independent confirmation.
2. On a line, give the approximate velocity range of the area (`data.velocity_range_kms`): with the wide default the F-J picks of a line can cover a narrow band, and high-frequency points beyond the reference can include a fast population ([troubleshooting.md](troubleshooting.md), section disp).
3. Resolution exists along the line only; the ends have few rays. Oblique noise raises the apparent velocity along the line; this workflow does not test it.
4. For instruments moved between lines, use the serial, start/end and line columns of the station table (section 2.2).

## 6. Assumptions and limits

- The noise field is close to isotropic. The phase-velocity bias from directional noise is neither tested nor corrected.
- Vertical components only; fundamental-mode Rayleigh phase velocities only.
- Vp and density follow Vs by the empirical relations of Brocher (2005), fixed inside DSurfTomo; where the real Vp/Vs differs, Vs is biased.
- DSurfTomo holds coordinates in single precision; pairs closer than about 200 m carry distance errors.
- The jackknife errors are not calibrated.
- The model can be interpreted only at the scales and depths the checkerboards recover; nodes outside the station coverage are unconstrained.

## 7. Records and reproduction

`OUT/.lant/`: `config.yaml` (full configuration of the last run), `config_history/`, `manifest.json` (what every step stands on), `run_log.jsonl` (start, end and duration of every step), `summary/<step>.json` (key numbers), `steps/<step>.md` (reports), `auto_parameters.json` (derived values and derivations), `provenance/` (code version, package versions, fingerprints of the third-party programs), `decisions.md` (your decisions and reasons). The `confirmed:` section at the end of the configuration records each confirmation with the fingerprint of the result it rests on. To reproduce: run the same configuration into a new output directory. The synthetic noise of the checkerboard tests has no fixed random seed; repeated checkerboards differ within the noise level.

## 8. Problems

Symptoms, causes and remedies, step by step: [troubleshooting.md](troubleshooting.md). What each output file holds: [outputs.md](outputs.md).
