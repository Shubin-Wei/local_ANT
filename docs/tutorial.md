# Tutorial: from raw records to a 3-D shear-velocity model

[中文](tutorial.zh.md) | English

This tutorial describes the purpose, method, basis and checks of every step. For a first run see [quickstart.md](quickstart.md). Principles and formulas: [methods.md](methods.md); defaults and what they rest on: [parameters.md](parameters.md); every output file: [outputs.md](outputs.md); problems: [troubleshooting.md](troubleshooting.md); running with an AI assistant: [ai_workflow.md](ai_workflow.md); numbers of the synthetic tests: [validation.md](validation.md).

Everything runs in Python; the inversion is a port of DSurfTomo, compared with the original in [python_backends.md](python_backends.md).

The rules file for AI assistants follows `project.ai_assistant`: `codex` AGENTS.md, `claude` (Claude Code) CLAUDE.md, `gemini` (Gemini CLI) GEMINI.md, one file with the full rules; the default writes AGENTS.md. See [AI workflow](ai_workflow.md) for context loading and existing-file handling.

Open `OUT/REPORT.md` after each step: it links the available main figures and `results/parameters.yaml`. The coverage, station diagnostics and array F-J spectrum are available before inversion; full model products are retained after evaluation. Caches remain under `work/`; see [outputs](outputs.md) and preview `lant clean-work --help` before choosing cache removal.

## 1. Installation

### 1.1 Environment

```bash
git clone https://github.com/Shubin-Wei/local_ANT.git && cd local_ANT
conda env create -f environment.yml && conda activate lant     # Python, ObsPy, GMT 6 and lant itself
```

If `conda env create` stays at "Solving environment" for a long time, use `--solver libmamba`, or mamba / micromamba. Without conda: `pip install -e .` on Python ≥ 3.11, plus GMT 6 (≥ 6.1). Linux only. No compiler is needed: every part of lant, including the inversion, is Python (NumPy, SciPy, Numba).

### 1.2 Check

```bash
lant check          # state of the Python packages and GMT, and how to install what is missing
```

### 1.3 Computing resources

Nobody sets the number of cores: `resources.cores: auto` takes 80 % of the cores other users leave idle, max(1, 0.8 × (usable CPUs − 1-minute load)), where the load this run caused itself does not count; it is taken again at the start of every step and of every solver project. Parallel work is also kept within half of the available memory (the rest waits). Single-pair diagnostic PNGs also render in separate processes under this plan and show their completed count; the measurement pool has finished before they start. Each step report lists the cores and the peak memory it used under COMPUTING RESOURCES. An integer in `resources.cores` is used as given.

Main figure groups also run in separate processes under this resource plan, each with one numerical thread. A single requested group or a one-process budget uses the serial path. Per-figure elapsed times are recorded; failed groups retain their previous complete images. The empirical memory estimate includes a 600 MB interpreter/rendering budget and work arrays, with the same available-memory headroom.

## 2. Preparing the data

Two inputs are needed: the raw-data directory and a station table. The records need no sorting into folders, renaming, format conversion, cutting into days or merging.

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
- `02_ccf` takes about pairs × 24 × (2 × lag window × rate + 1) × 4 bytes per day; the clean report gives this estimate and the memory of one cc process. The memory counted is what is available now, within the limits of a job or container (Linux cgroups), with reclaimable file cache counted as available. When half of it does not hold one process per core, the cc step runs fewer at a time; when not even one fits, the clean step warns and the cc step stops. Free memory or use a larger machine; a narrower `data.velocity_range_kms` or `correlation.max_distance_km` (only pairs within that distance) also reduces it, but change them only on independent knowledge of the area, not to make the run fit.
- After a run is finished and evaluated, `lant clean-work -c CONFIG` lists the rebuildable caches (cleaned day files, hourly correlations) and their size; `--execute --yes` removes them. Results and the report stay readable; rebuilding cc or stack, or redrawing the stack-time figures, then needs the caches again.
- The inversion time depends on the number of grid nodes and of "stations × periods"; every checkerboard and control is one project; as many run side by side as the cores and the memory allow.

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

Keys to review before the first run: `data.velocity_range_kms` (default 0.1–5 km/s, from soft soil to hard rock; a narrower range based on independent knowledge of the area makes the lag window and the F-J scan smaller and the run faster), `data.response`, `data.day_label` (default UTC). Everything else has a default or is derived from the data.

## 4. Running

```bash
lant run all -c CONFIG          # from the first step that does not stand; stops where a judgement is yours
lant status -c CONFIG           # state of every step, why, and the next command (--json for programs)
lant confirm STEP -c CONFIG     # after you looked at the report and the figures
lant run all -c CONFIG --from disp     # run again from a step after changing a parameter
```

While a step runs, one progress line shows its place in the chain, what it counts (station-days, days, pair batches, solver projects with their iterations, figures) and the elapsed time; in a terminal it is a bar updated every second, in a log one line every 30 seconds; `LANT_PROGRESS=off` turns it off. Detailed output goes to `.lant/lant.log` (`-v` shows it as well). A finished chain prints where the report and results are.

An interrupted run (Ctrl-C, `kill`, a closed terminal; under `nohup` a closed terminal does not stop it) stops the solver processes, records the step as failed and releases the lock; `lant run all -c CONFIG` continues, and interrupted solver projects restart from their inputs (their partial output is kept in `.lant/backups/`). Even after `kill -9` or running out of memory the next command records the step as failed instead of refusing to run. While a step runs, `lant status` shows it as running and `lant note` still works.

Every step stores a report in `.lant/steps/<step>.md` (`-v` prints it): what it did, key numbers, every parameter it used (value, basis, source), files to look at, items left for you. `.lant/manifest.json` records on which inputs and parameters every step stands; when one of them changes, the step and those after it become stale and confirmations lapse. A step with nothing to judge does not stop.

### 4.1 clean: cleaning (`01_clean`)

Header scans reuse unchanged files using nanosecond file metadata; old cache signatures are refreshed once. When all headers are cached, no scanning pool starts and the cache is not rewritten. Cleaned station-days use the same file identity, including response sources when supplied. Changing a response file marks clean out of date; upgrading the reuse key requires one rebuild.

**What**: reads all headers; matches the station table to the records; chooses the vertical channel; cuts processing days; puts the samples on the UTC sample grid (offsets below 0.05 sample are shifted, others interpolated with a Lanczos kernel); removes duplicates and identical overlaps; removes the response when response files are given; tests every station-day for a coherent pulse once per second with an F test on the record folded at 1 s (Takagi et al. 2015); judges the array type from the station distribution. Apart from response removal: no filtering, no normalisation. A second run reuses every station-day whose input files and settings did not change.

**Look at**: `report.txt` (unmatched names, stations without data, resource estimate for cc), `0_data_coverage.png` (data per station and day), `stations.csv`, `pulse_1s.csv` (stations with the pulse).

**Left to you**: ambiguous names, coordinate mismatches, unknown component orientation, an unclear array type (elongation 3–10), responses not removed and models unknown.

### 4.2 cc: hourly cross-coherence (`02_ccf`)

**What**: every station-day is decimated to the processing rate and cut into segments with 50 % overlap (Seats et al. 2012); segments whose peak exceeds `max_over_std` times the standard deviation of the station-day are dropped (Jiang & Denolle 2020); each segment gives the cross-coherence (each spectrum divided by its own amplitude, Xu et al. 2021), averaged per hour. Records of the neighbouring days are read too, so segments across a day boundary are kept. For stations where clean detected the 1-s pulse, the folded waveform is subtracted hour by hour first; other stations are left alone.

**Parameters**: sampling rate, segment length, lag window and band follow from the geometry and the velocity range by written rules (the report shows each derivation with your numbers).

### 4.3 stack: stack and station diagnosis (`03_stack`)

**What**: all hours stacked, weighted by their segments; sub-stacks of several hours are kept for the jackknife (`stack.block_hours`). Diagnosis: every sub-stack of a pair is compared with the stack of the pair's other sub-stacks (leave one out) by correlation coefficient and sub-sample time shift (Stehly et al. 2007; Sens-Schönfelder 2008); per station the median over its pairs with the nearest stations. Candidates are the outliers of the two distributions: the correlation enters as Fisher z, the scale is the median and 1.4826 × MAD (Leys et al. 2013), and the limit follows from the sample size by Chauvenet's criterion, not from a fixed number. A time shift is listed only when the shifts of the pairs of the station agree, i.e. the median shift is beyond the same limit in units of its standard error over the pairs: a clock error moves every pair of a station alike, while in sub-stacks of low correlation the shifts scatter.

**The diagnosis lists, it does not remove.** `candidates.txt` and `results/figures/2_station_qc.png` are for you; when there are more candidates than the criterion lets through by chance (1), the chain stops and you decide: `lant confirm stack` (keep) or `lant confirm stack --exclude-candidates` (remove), or write your own list in `stack.exclude_station_days`.

### 4.4 reference: reference curve (`04_reference`)

**What**: the frequency-Bessel spectrum of the array (Wang et al. 2019; double precision, per frequency only the pairs within a few wavelengths), the ridge picked and fitted by a smooth curve in log f: the reference phase-velocity curve. From it follow the lower end of the measurement band and the velocity axis. The reference only fixes the cycle count of the phase and guides the group-arrival pick; it does not enter the final values.

**Look at**: `3_fj_phase_velocity.png`. **This is the most important judgement of the chain: is the picked ridge the fundamental mode?** Where a higher mode is stronger, `max_peak` picks it; then use `track_fundamental` (it follows the ridge from the start, tried every 1/3 octave, that can be followed over the most frequencies, within `data.velocity_range_kms` and, on a line, only at wavelengths of at least twice the mean station spacing, drawn as a dashed line in `3_fj_phase_velocity.png`), narrow `data.velocity_range_kms`, or give a curve with `reference.source: file`. The report also compares the measured velocities with `data.velocity_range_kms`.

**Large arrays with strong lateral variation**: `reference.zone_km` divides the array into square zones by path midpoint, one curve per zone (`3_reference_zones.png`); zones without a clear ridge use the curve of the whole array.

### 4.5 disp: dispersion measurement, quality control, export (`05_dispersion`)

**What**:
1. The symmetric component of every pair is filtered with Gaussian narrow-band filters (as FTAN, Bensen et al. 2007) at log-spaced frequencies (12 per octave, Xu et al. 2021); the group arrival is picked along a ridge guided by the reference.
2. The phase is read with a phase-matched filter (Herrin & Goforth 1977; Levshin & Ritzwoller 2001): the measured group-arrival curve gives a phase correction that removes the phase variation inside the envelope caused by dispersion, and the phase is then read at the centre frequency; far-field phase constant π/4 (Lin et al. 2008). Travel time t = −(Φ − π/4)/ω plus whole cycles, phase velocity c = distance / t (two-station method, Yao et al. 2006), with distances on the WGS84 ellipsoid (Vincenty 1975).
3. The cycle count is fixed where the path is 1–3.5 reference wavelengths (there the reference cannot be a cycle off) and then followed continuously in frequency.
4. Errors: delete-a-group jackknife over groups of sub-stacks, each measured again. The errors are not calibrated; they show repeatability in time only.
5. Quality control: every rejected sample carries its reasons (bit mask), counted in the report. Criteria include signal-to-noise ratio, jackknife error, a path of at least 2 of the sample's own wavelengths (Luo et al. 2015), run length, and agreement of group and phase velocity; the basis of every threshold is in parameters.md, and those labelled "empirical" have no literature value.
6. Export: per pair only the longest continuous segment, written as DSurfTomo input.

**Look at**: `results/figures/3_dispersion_curves.png` (all curves and the reference), `pair_figures/` (single pairs: velocity spectrum, picks, reasons), `median_vs_reference.txt`.

**Do not loosen the quality control to get more points.** A changed threshold needs its reason, recorded with `lant note`.

### 4.6 inversion: building the projects (`06_inversion`)

**What**: 1-D inversion of the array's average curve for the starting model (damped least squares, smoothness applied to the model itself, Constable et al. 1987); depth nodes in steps of max(shortest wavelength / 4, depth / 4) down to half the longest wavelength (rule of thumb, Foti et al. 2018); grid spacing the median measured wavelength / 3; the main project, four checkerboard projects (Lévêque et al. 1993; three horizontal with cells of 1, 2 and 3 median measured wavelengths, and one depth-varying, methods.md section 9) and by default one control project, another starting model (`start_grad0.05`; Rawlinson & Spakman 2016); the smoothing weight halved and doubled (`reg_x0.5`, `reg_x2`) can be added in `checks.controls`. The regularization smooths the model itself, the data and model terms as means (methods.md section 8.4); its weight (1) is a fixed default calibrated on the areal synthetic example (no damping), and the check projects use the same values. The 3-D inversion takes the measured periods on a grid of `inversion.per_octave` (4) per octave, every third of the 12 measured (methods.md 8.3); at most 60 periods (DSurfTomo's CalSurfG); with more, the step stops and names the `inversion.per_octave` that fits.

**Look at**: `results/figures/4_start_model.png`, `4_kernels.png` (depth sensitivity of the periods). **Left to you**: whether starting model, deepest node and velocity bounds are reasonable.

### 4.7 run, evaluate, figures

- **run**: DSurfTomo in Python (Fang et al. 2015: direct inversion of the travel times for 3-D Vs, rays updated with the model; `python -m lant.inversion.solver` per project). Before it runs, the share of data its rejection rule keeps in the first iteration is checked (`checks.precheck`); an invalid main project stops the chain. The main project and the controls run side by side; the checkerboards follow when the main project has finished, since their noise is the relative misfit of the main model (`inversion.noise: auto`; the value is in the run report and in `work/06_inversion/checkerboard_noise.json`). A number in `inversion.noise` lets all projects run together.
- **evaluate**: numbers and figures only, no verdicts. Travel-time residuals (starting and final model, share kept), checkerboard recovery (correlation and regression slope by depth, on nodes with nonzero post-rejection data-Jacobian column support: DWS > 0), differences between control models and the main model together with the residual history of each control (a control whose residual has not come down to that of the main project is not converged). Written to `work/06_inversion/evaluation/`.
- **figures**: in `results/figures/` (list in outputs.md): station map over relief, correlation moveouts per band, three single pairs (stack time and dispersion), dispersion curves, starting model and kernels, depth slices and a 3-D view (areal arrays), sections, checkerboards, residuals. Scales, distances, lag windows and velocity ranges follow the data of the project; the relief is SRTM, loaded automatically from the GMT data server (`plots.topography`). The model of a line is a narrow strip along it: no depth slices or 3-D view, and the checkerboards are drawn as input and recovered sections along the line. A section of your own that leaves the model is cut to the model, and the note says by how much.

The starting model is drawn as connected parameter nodes: the last value is held fixed in the 3-D inversion, not a measured interface. Kernel curves are labelled as a separate extended-model illustration. In the residual figure, panel (a) uses milliseconds for every observation; the table separates mean, standard deviation and RMS. Panel (b) compares periods, not iterations. Read the checkerboards separately before interpreting the model.

## 5. Line arrays

One project per line. The defaults that differ are listed at the end of parameters.md.

1. The F-J spectrum of a line has pairs in one direction only and is often affected by higher modes or obliquely arriving noise; the default is `track_fundamental`; do look at `3_fj_phase_velocity.png`. Where an areal array of the same area exists, its reference curve and starting model can be borrowed (`reference.source: file`, `inversion.start_model.source: file`); agreement between the line and the areal result is then no independent confirmation.
2. On a line, give the approximate velocity range of the area (`data.velocity_range_kms`): with the wide default the F-J picks of a line can cover a narrow band, and high-frequency points beyond the reference can include a fast population ([troubleshooting.md](troubleshooting.md), section disp).
3. Resolution exists along the line only; the ends have few rays. Oblique noise raises the apparent velocity along the line; this workflow does not test it.
4. For instruments moved between lines, use the serial, start/end and line columns of the station table (section 2.2).

## 6. Assumptions and limits

- The noise field is close to isotropic. The phase-velocity bias from directional noise is neither tested nor corrected.
- Vertical components only; fundamental-mode Rayleigh phase velocities only.
- Vp and density follow Vs by the empirical relations of Brocher (2005), fixed inside DSurfTomo (and its port); where the real Vp/Vs differs, Vs is biased.
- DSurfTomo holds coordinates in single precision (about 0.4 m). Observed and predicted travel times use the same coordinates; on the test projects the pairs that reach the inversion (0.5 km and longer) differ from double-precision distances by at most 0.2 %. Pairs much shorter than 0.2 km would carry errors of a percent or more; the wavelength rule of the dispersion step normally removes them.
- The jackknife errors are not calibrated.
- The model can be interpreted only at the scales and depths the checkerboards recover; nodes outside the station coverage are unconstrained.

## 7. Records and reproduction

`OUT/.lant/`: `config.yaml` (full configuration of the last run), `config_history/`, `manifest.json` (what every step stands on), `run_log.jsonl` (start, end and duration of every step), `summary/<step>.json` (key numbers), `steps/<step>.md` (reports), `auto_parameters.json` (derived values and derivations), `provenance/` (code version, package versions), `../notes.md` (your decisions and reasons). The `confirmed:` section at the end of the configuration records each confirmation with the fingerprint of the result it rests on. To reproduce: run the same configuration into a new output directory. The synthetic noise of the checkerboard tests has a fixed random seed.

## 8. Problems

Symptoms, causes and remedies, step by step: [troubleshooting.md](troubleshooting.md). What each output file holds: [outputs.md](outputs.md).

The parent publishes `results/main_preview/` immediately when the main solver future completes, before other projects finish. It retains the model, solver inputs, effective parameters and a section when plotting is available; `status.json` explicitly says checks may be pending and records plotting failure without discarding the model. This is a preview, not scientific acceptance or completed evaluation.

`lant clean-work --cache forward` previews registered exact-forward `.npz` caches; explicit execution removes only unchanged registered files after evaluation, without invalidating processing stages. The default preview includes archive, cc and forward caches. Unknown/modified files, models, parameters and final figures remain protected. Large clean/cc caches are also removed only by an explicit cleanup request.
