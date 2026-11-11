# Output files

[中文](outputs.zh.md) | English

All results are written under `project.output_dir` of the configuration; the raw data are never changed. The tables below list, directory by directory, what each file holds, which step writes it and when you need it. Files marked ★ are worth a look after every run.

## Project entry points

Open `REPORT.md` first. `README.md` explains ownership and operating rules. `notes.md` holds appended notes and decisions. `results/dispersion/` is replaced after a complete disp step. `results/model/` (solver inputs and model), `results/validation/` (finished checks and the evaluation) and `results/reproduce/` (effective configuration, station geometry, starting-model files) are replaced together after a complete evaluate step, so rebuilding or rerunning the inversion keeps the previous evaluated model until the new one is evaluated. `REPORT.md` names the run each of them comes from and marks it when the step has run again since. Scientific products remain available independently of working caches. `work/` is not disposable as a whole.


## .lant/: the run record

| File | Content |
|---|---|
| ★ `steps/<step>.md` | report of each step: what was done, the key numbers, every parameter used (value, source category, origin), files to look at, items left for your judgement, what to change to change the result |
| `summary/<step>.json` | all numbers of the same report, for programs |
| `manifest.json` | the inputs, parameters and upstream results each step stands on; after a change of parameters or inputs `lant status` uses it to tell which steps are stale |
| `auto_parameters.json` | every parameter derived from the data by a rule, with its derivation (the numbers substituted) |
| `config.yaml`, `config_history/` | the complete configuration of the latest run (all defaults included) and the earlier ones |
| `run_log.jsonl` | start, end, duration and outcome of every step |
| `provenance/` | code version (with a fingerprint of uncommitted changes), versions of dependencies, fingerprints of third-party programs |
| `../notes.md` | your judgements and reasons (written by hand, or by an AI assistant from your words) |
| `lock` | run lock; only one `lant` process can write to an output directory at a time |

## work/01_clean/: preparing the raw records

| File | Content |
|---|---|
| ★ `report.txt` | numbers of files, formats, unreadable files, name matching, stations and processing days, derived parameters, resource estimate, result of the 1-s signal test |
| ★ `coverage.png` | hours of data per station and processing day |
| ★ `stations.csv` | the stations used: `NET.STA`, longitude, latitude, elevation (all later steps use these names) |
| `processing_days.txt` | label of each processing day and its UTC window |
| `station_days.csv` | per station-day: hours of data, gaps, shifted and interpolated records, overlaps, dropped short pieces |
| `pulse_1s.csv` | statistic, limit and effective degrees of freedom of the 1-s folding test per station-day |
| `input_files.csv`, `raw_index.json` | every file and record read |
| `archive/` | cleaned day files (one miniSEED per station and day; `*.done.json` are fingerprints for resuming); absent when `clean.write_archive: false` |
| `archive_index.json` | index of the records read by the later steps |
| `header_cache.json` | cache of the header scan (files that did not change are not read again) |

## work/02_ccf/: hourly cross-coherence

One `day_YYYYDDD.h5` per processing day: `hourly [24, pairs, lags]`, the number of valid segments per hour `count`, pair names, settings and input fingerprints. This is by far the largest directory. Once the stack is done and no hours need to be removed any more it can be deleted, but then changing any parameter before cc requires the correlations to be computed again.

## work/03_stack/: stack and station diagnosis

| File | Content |
|---|---|
| `stack.h5` | stacked correlations `stack`, sub-stacks `daily`, segments per sub-stack, distances, azimuths, station coordinates |
| ★ `station_qc/station_qc.png` | correlation (left) and time shift (right) of the sub-stacks of each station, candidates in red |
| `station_qc/station_qc.txt` | values and origin of the limits, all candidates with their reasons |
| `candidates.txt` | list of candidates (in the format used when confirming their removal) |
| `removed.txt` | station-days and station-hours actually removed |

## work/04_reference/: reference curve

| File | Content |
|---|---|
| ★ `fj.png` | frequency–Bessel spectrum (normalised per frequency), picks, reference curve, measured band, your velocity range |
| `fj.txt` | per frequency: pairs used, picked value, whether it was used, reference value |
| `fj.npz` | arrays of the F-J image |
| `reference.npz` | the reference curve cut to the measured band (read by the disp step) |
| `zones.npz`, `zones.png` | zone references (only with `reference.zone_km`) |

## work/05_dispersion/: dispersion measurement

| File | Content |
|---|---|
| `disp.h5` | per pair and measured frequency: phase velocity `vc`, group velocity `vg`, signal-to-noise ratio, jackknife error, number of cycles, accepted `ok_c`, bit mask of rejection reasons `reason`; band-edge ratios `band_amp` |
| ★ `median_vs_reference.txt` | per frequency, median of the accepted points minus the reference; a systematic offset points to a problem with the reference curve or the number of cycles |
| `curves.csv` | Measured samples with QC and inversion-selection flags; separate inversion columns preserve six-decimal periods and five-decimal velocities |
| `pair_figures/` | single-pair figures (`export.curve_figures` pairs spread evenly over distance, plus the pairs with most rejected points and with the largest departure from the median): velocity–frequency image, picks, rejection reasons |

Bits of `reason`: 1 no group pick, 2 signal-to-noise too low, 4 no data, 8 group velocity off the reference, 16 path too short, 32 side lobe, 64 number of cycles undetermined, 128 cycle skip, 512 jackknife error, 4096 group and phase velocities inconsistent, 8192 run too short, 16384 wavelength rule, 32768 no reference curve, 65536 band edge.

## work/06_inversion/: inversion projects

| File | Content |
|---|---|
| `average_curve.txt` | average curve of the array (median and number of pairs per period) |
| `vs1d.txt`, `vs1d_fit.txt` | 1-D starting model and its fit to the average curve |
| `kernels.txt` | depth sensitivity kernels per period |
| `main/` | main project: `surfdata.dat`, `DSurfTomo.in`, `MOD`, `project.json` (grid, periods, sparse reserve), `precheck.json`, `run.log`, final model `DSurfTomo.inMeasure.dat`, models of each iteration `*.iterNNN` |
| `checkerboards/<name>/` | checkerboard projects |
| `checkerboard_noise.json` | the noise of the checkerboard travel times (`inversion.noise: auto`): value and how it was obtained |
| `controls/<name>/` | control projects (starting model, regularisation) |
| `evaluation/evaluation.json` | residuals, checkerboard recovery (correlation and slope per depth), differences between controls and the main model |
| `evaluation/residuals_final.csv` | every travel time: station coordinates, period, distance, observed and predicted time, kept by the rejection rule or not |

`DSurfTomo.inMeasure.dat` has four columns: longitude, latitude, depth (km), Vs (km/s).

## results/figures/: result figures

| Figure | Content | Arrays |
|---|---|---|
| `1_stations` | stations over coloured, hill-shaded relief (SRTM, `plots.topography`), stations left out by `data.exclude_stations` as red crosses | all |
| `2_ccf_moveout_<band>` | stacked correlations sorted by distance, one figure per band (the measurement band cut at round frequencies); distances up to where pairs were measured in the band; the two dashed moveouts are the slowest group and the fastest phase velocity measured in the band | all |
| `2_ccf_stack_time_<k>` | three pairs (short, middle, long distance): stacked correlation, hourly correlations in time order, signal-to-noise ratio of the symmetric component against stacked hours | all |
| `3_dispersion_pair_<k>` | the same three pairs: waveform; phase velocities with the other numbers of cycles and the reference; group-velocity image with the measured arrivals | all |
| `3_dispersion_curves` | all accepted phase-velocity curves, mean and standard deviation | all |
| `4_start_model`, `4_kernels` | starting model and fit, sensitivity depth per period | all |
| `5_depth_slices` | Vs at the slice depths over grey relief, with the section lines | areal |
| `5_slices_3d` | five slices in a perspective view, the relief above them | areal |
| `5_sections` | sections with the elevation profile above, ends marked with letter and direction; at the length-to-depth ratio of the model (a vertical exaggeration is stated when needed) | all (lines: along the line) |
| `6_checkerboard_*` | checkerboard input and recovery, red (slow) – white – blue (fast): maps at up to seven depths (areal), input and recovered section (line) | all |
| `6_recovery` | recovery slope and correlation per depth | all |
| `6_residuals` | travel-time residuals of the starting and the final model: histogram and rms per period | all |

A figure that cannot be drawn is listed in the figures report with the reason (e.g. a band without measurements, or no relief available offline).
PNG is the default; set `plots.formats: [png, pdf]` for vector PDF output as well.
