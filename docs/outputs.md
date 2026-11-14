# Output files

[中文](outputs.zh.md) | English

All results are written under `project.output_dir` of the configuration; the raw data are never changed. The tables below list, directory by directory, what each file holds, which step writes it and when you need it. Files marked ★ are worth a look after every run.

## Overview

Open `REPORT.md` first: the state and the time of every step and the next command, key numbers (including the smoothing weight and damping), the items left for you to judge, the main figures, and links to the retained results with the run each comes from. `README.md` explains the directories and the operating rules. `notes.md` holds your notes and decisions (`lant note`).

| File | Content |
|---|---|
| `REPORT.md`, `README.md` | the report and the guide to the directory (written by lant) |
| `AGENTS.md`, `CLAUDE.md` or `GEMINI.md` | the rules for AI assistants: one file with the full text, for the assistant chosen in `project.ai_assistant` (codex AGENTS.md, claude CLAUDE.md, gemini GEMINI.md; none chosen: AGENTS.md); an existing file is never overwritten |
| `notes.md` | your notes and decisions (`lant note`) |
| `scratch/` | optional: temporary files of people and assistants; not checked by `lant status` (which reports its size), safe to delete |

Any other file here, or in `work/`, `results/` or `.lant/`, that lant did not write is listed by `lant status` as unregistered; lant never deletes or overwrites it.

## results/: retained results

`results/dispersion/` is replaced after a complete disp step; `results/model/` and `results/validation/` are replaced together after a complete evaluate step, so rebuilding or rerunning the inversion keeps the previous evaluated model until the new one is evaluated. Unregistered files in these directories are never replaced. They stay readable after `lant clean-work`; `work/` is not disposable as a whole.

Each completed step rewrites `results/parameters.yaml`: the configuration with every value resolved so far and the record of your confirmations (the values each step used are in its report, `.lant/steps/<step>.md`, and in `.lant/manifest.json`). The cleaning, station-diagnosis and F-J figures are moved into `results/figures/` when their step completes (they are kept only there); `REPORT.md` links them and their originating runs. A rerun that disables a diagnostic keeps its earlier published figure with its earlier run label.

| File | Content |
|---|---|
| ★ `model/vs_model.csv` | the 3-D model: `longitude_deg`, `latitude_deg` (WGS84), `depth_km`, `vs_km_s`, sorted by depth, latitude, longitude; the values exactly as DSurfTomo wrote them. Nodes far from the ray paths keep the starting model: read it with the checkerboards |
| ★ `dispersion/curves.csv` | every pair and measurement frequency: measured values, `accepted_qc`, `used_in_inversion` and the values given to DSurfTomo |
| `parameters.yaml` | every value resolved, with the confirmation record; written by the latest completed step (some auto values resolve only in later steps) |
| `figures/` | result figures (below) |
| `validation/evaluation.json` | residual statistics of the starting and final model, checkerboard recovery per depth, control inversions |
| `validation/residuals_start.csv`, `residuals_final.csv` | travel time observed and predicted per pair and period, kept or rejected by DSurfTomo's rule |

To rerun or redraw, everything is in `work/` (nothing has to be computed again): the solver files of the main model are in `work/06_inversion/main/` (`DSurfTomo.in`, `surfdata.dat`, `MOD`, final model `DSurfTomo.inMeasure.dat`), the starting model and average curve in `work/06_inversion/`, the station table in `work/01_clean/stations.csv`. All measurements with their rejection reasons and diagnostics (`disp.h5`, `median_vs_reference.txt`) and the per-pair diagnostic figures stay in `work/05_dispersion/` (no copy in `results/`); the checkerboard and control projects with their recovered models stay in `work/06_inversion/`. Checkerboard figures are drawn for the check projects of the present configuration (`work/06_inversion/checks.json`), the same ones the evaluation uses; a project directory left from an earlier configuration is kept in `work/` but not drawn.

## .lant/: the run record

| File | Content |
|---|---|
| ★ `steps/<step>.md` | report of each step: what was done, the key numbers, every parameter used (value, source category, origin), files to look at, items left for your judgement, what to change to change the result |
| `summary/<step>.json` | all numbers of the same report, for programs |
| `manifest.json` | the inputs, parameters and upstream results each step stands on; after a change of parameters or inputs `lant status` uses it to tell which steps are stale |
| `auto_parameters.json` | every parameter derived from the data by a rule, with its derivation (the numbers substituted) |
| `config.yaml`, `config_history/` | the complete configuration of the latest run (all defaults included) and the earlier ones |
| `run_log.jsonl` | start, end, duration and outcome of every step |
| `provenance/` | code version (with a fingerprint of uncommitted changes), versions of dependencies, fingerprints of third-party programs; package source SHA-256 frozen at process import |
| `forward_cache/*.npz` | dependency- and checksum-verified forward predictions shared by prechecks, automatic noise and residuals; internal rebuildable data, not a result model |
| `../notes.md` | your judgements and reasons, appended with `lant note` (an AI assistant writes your words, not its own decisions) |
| `lock` | run lock; only one `lant` process can write to an output directory at a time |
| `artifacts.json`, `published.json` | the files lant wrote (for `lant status` and `lant clean-work`) and the run each retained result comes from |
| `backups/` | step outputs kept by `--rebuild`, and the output of interrupted solver projects |
| `lant.log` | detailed output of the steps (10 MiB per file, three older files kept) |

`phase_start` / `phase` entries in `run_log.jsonl` record prechecks, solver batches/projects, automatic noise and evaluation wall times and outcome; `forward_evaluation` records hit/miss, identity checks, actual computation and total time separately. Concurrent or nested phases are not directly additive.

## work/01_clean/: preparing the raw records

| File | Content |
|---|---|
| ★ `report.txt` | numbers of files, formats, unreadable files, name matching, stations and processing days, derived parameters, resource estimate, result of the 1-s signal test |
| ★ `stations.csv` | the stations used: `NET.STA`, longitude, latitude, elevation (all later steps use these names) |
| `processing_days.txt` | label of each processing day and its UTC window |
| `station_days.csv` | per station-day: hours of data, gaps, shifted and interpolated records, overlaps, dropped short pieces |
| `pulse_1s.csv` | statistic, limit and effective degrees of freedom of the 1-s folding test per station-day |
| `input_files.csv`, `raw_index.json` | every file and record read; a miniSEED file whose first and last records are of a station that cannot be in the station table is read only at its ends (`not_read_because` says so) |
| `archive/` | cleaned day files (one miniSEED per station and day; `*.done.json` are fingerprints for resuming); absent when `clean.write_archive: false` |
| `archive_index.json` | index of the records read by the later steps |
| `header_cache.json` | cache of the header scan (files that did not change are not read again) |

## work/02_ccf/: hourly cross-coherence

One `day_YYYYDDD.h5` per processing day: `hourly [24, pairs, lags]`, the number of valid segments per hour `count`, pair names, settings and input fingerprints. This is by far the largest directory. Once the stack is done and no hours need to be removed any more it can be deleted, but then changing any parameter before cc requires the correlations to be computed again.

## work/03_stack/: stack and station diagnosis

| File | Content |
|---|---|
| `stack.h5` | stacked correlations `stack`, sub-stacks `daily`, segments per sub-stack, distances, azimuths, station coordinates |
| `station_qc/station_qc.txt` | values and origin of the limits, all candidates with their reasons |
| `candidates.txt` | list of candidates (in the format used when confirming their removal) |
| `removed.txt` | station-days and station-hours actually removed |
| `all_hours.json` | written when `stack.h5` is the stack of every hour: the fingerprint of its inputs, so that the run after `lant confirm stack --exclude-candidates` removes the candidates from it instead of stacking every hour again |

## work/04_reference/: reference curve

| File | Content |
|---|---|
| `fj.txt` | per frequency: pairs used, picked value, whether it was used, reference value |
| `fj.npz` | arrays of the F-J image |
| `reference.npz` | the reference curve cut to the measured band (read by the disp step) |
| `zones.npz` | zone references (only with `reference.zone_km`) |

The figures of clean (data coverage), stack (station diagnosis) and reference (F-J spectrum, zones) are written here by the step and moved to `results/figures/` (`0_data_coverage`, `2_station_qc`, `3_fj_phase_velocity`, `3_reference_zones`) when it completes; a failed step leaves them here.

## work/05_dispersion/: dispersion measurement

| File | Content |
|---|---|
| `disp.h5` | per pair and measured frequency: phase velocity `vc`, group velocity `vg`, signal-to-noise ratio, jackknife error, number of cycles, accepted `ok_c`, bit mask of rejection reasons `reason`; band-edge ratios `band_amp` |
| ★ `median_vs_reference.txt` | per frequency, median of the accepted points minus the reference; a systematic offset points to a problem with the reference curve or the number of cycles |
| `curves.csv` | every pair and measurement frequency: measured values, QC and inversion flags, and the values given to DSurfTomo (periods to six decimals, velocities to five, as before) |
| `pair_figures/` | single-pair figures (`export.curve_figures` pairs spread evenly over distance, plus the pairs with most rejected points and with the largest departure from the median): velocity–frequency image, picks, rejection reasons |

Bits of `reason`: 1 no group pick, 2 signal-to-noise too low, 4 no data, 8 group velocity off the reference, 16 path too short, 32 side lobe, 64 number of cycles undetermined, 128 cycle skip, 512 jackknife error, 4096 group and phase velocities inconsistent, 8192 run too short, 16384 wavelength rule, 32768 no reference curve, 65536 band edge.

## work/06_inversion/: inversion projects

| File | Content |
|---|---|
| `average_curve.txt` | average curve of the array (median and number of pairs per period) |
| `vs1d.txt`, `vs1d_fit.txt` | 1-D starting model and its fit to the average curve |
| `kernels.txt` | depth sensitivity kernels per period |
| `main/` | main project: `surfdata.dat`, `DSurfTomo.in`, `MOD`, `project.json` (grid, periods, sparse reserve, weight and damping), `precheck.json`, `run.log`, `solver.json` (implementation, threads, seconds per module and iteration, peak memory), `DWS.dat` (DWS of every node, last iteration), final model `DSurfTomo.inMeasure.dat`, models of each iteration `*.iterNNN` |
| `checkerboards/<name>/` | checkerboard projects: inputs, `MOD.true`, true and recovered model (`Vs_model.real`, `DSurfTomo.inSyn.dat`), logs. After the evaluation the DSurfTomo files lant does not read (`*.iterNNN`, `residualFirst.dat`, `residualLast.dat`, `velmap2dRc.dat`) are deleted here and in `controls/`, only when lant wrote them and they are unchanged; `main/` keeps them |
| `checkerboard_noise.json` | the noise of the checkerboard travel times (`inversion.noise: auto`): value and how it was obtained |
| `controls/<name>/` | control projects (starting model, regularisation) |
| `evaluation/evaluation.json` | residuals, checkerboard recovery (correlation and slope per depth), differences between controls and the main model |
| `evaluation/residuals_final.csv` | every travel time: station coordinates, period, distance, observed and predicted time, kept by the rejection rule or not |

`DSurfTomo.inMeasure.dat` has four columns: longitude, latitude, depth (km), Vs (km/s).

## results/figures/: result figures

| Figure | Content | Arrays |
|---|---|---|
| `0_data_coverage` | data coverage after cleaning; published at clean | all |
| `2_station_qc` | station diagnostics; candidates require your decision; published at stack when enabled | all |
| `3_fj_phase_velocity`, `3_reference_zones` | array F-J spectrum and optional zoned references; published at reference when computed, distinct from single-pair energy | all |
| `1_stations` | recorded stations over coloured, lightly shaded relief (SRTM, `plots.topography`), bounded map height with internal legend separating recorded and inversion station counts; stations left out by `data.exclude_stations` as red crosses | all |
| `2_ccf_moveout_<band>` | stacked correlations sorted by distance, one figure per band (the measurement band cut at round frequencies); distances up to where pairs were measured in the band; the two dashed moveouts are the slowest group and the fastest phase velocity measured in the band | all |
| `2_ccf_stack_time_<k>` | three pairs (short, middle, long distance): stacked correlation, hourly correlations in time order, signal-to-noise ratio of the symmetric component against stacked hours | all |
| `3_dispersion_pair_<k>` | the same three pairs: waveform; phase-matched filter picks with integer-cycle branches and the F-J reference (no classical narrow-band phase background); group-velocity envelope image with measured arrivals | all |
| `3_dispersion_curves` | all accepted phase-velocity curves, mean and standard deviation | all |
| `4_start_model` | median observed phase velocities and the prediction used to initialize the 1-D model, with RMS; connected Vs parameter nodes and fixed bottom value, not detected layer interfaces | all |
| `4_kernels` | illustrative 1-D joint sensitivity on a separate finer, extended model (Vs +/-5%, Vp/rho follow), per km; not the actual inversion Jacobian; the dashed line is the inversion bottom boundary | all |
| `5_depth_slices` | Vs at the slice depths over grey relief, with the section lines | areal |
| `5_slices_3d` | five slices in a perspective view, the relief above them | areal |
| `5_sections` | sections with reference elevation profile above, sparse isolines and ends marked with letter and direction; line arrays include the percentage update relative to the 1-D starting model; vertical exaggeration and unmasked resolution noted at the bottom | all (lines: along the line) |
| `6_checkerboard_*` | checkerboard input and recovery, red (slow) – white (zero) – blue (fast): maps at up to seven depths (areal), input and recovered section (line) | all |
| `6_recovery` | test-by-depth heatmaps and depth profiles of amplitude recovery factor (100% / 1.0 ideal) and pattern correlation r (1.0 ideal) within the near-station evaluation region | all |
| `6_residuals` | all-observation residual histogram in ms, relative RMS by period; mean, population standard deviation and RMS are listed separately, with kept/all counts; actual iteration curves only when already evaluated | all |

A figure that cannot be drawn is listed in the figures report with the reason (e.g. a band without measurements, or no relief available offline).
Workflow figures are PNG only (`plots.formats: [png]`). The figure groups use the README overview style; maps of input and recovered anomalies share colour limits.

Recovery heatmaps support a single depth, check aligned increasing depths, and retain negative correlations and amplitude overshoots. Cells display node samples rather than resolved depth intervals. Section anomalies are explicitly relative to the 1-D reference, which need not be a heterogeneous initial model.

The parent publishes `results/main_preview/` immediately when the main solver future completes, before other projects finish; the evaluate step removes it when it publishes the evaluated model in `results/model/` (only when every file in it is still as lant wrote it). It retains the model, solver inputs, effective parameters and a section when plotting is available; `status.json` explicitly says checks may be pending and records plotting failure without discarding the model. This is a preview, not scientific acceptance or completed evaluation.

`lant clean-work --cache forward` previews registered exact-forward `.npz` caches; explicit execution removes only unchanged registered files after evaluation, without invalidating processing stages. The default preview includes archive, cc and forward caches. Unknown/modified files, models, parameters and final figures remain protected. Large clean/cc caches are also removed only by an explicit cleanup request.

Main-preview publication also prints its path immediately and the project README documents it, so it can be found before the complete run report is published.
