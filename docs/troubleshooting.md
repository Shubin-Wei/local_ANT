# Troubleshooting

[中文](troubleshooting.zh.md) | English

Start with `lant status -c CONFIG`: it gives the state of each step, why the chain stopped and the next command. Add `--debug` to an erroring command to see the full Python traceback. The entries below are ordered by the step where the problem appears.

## Installation and check

| Symptom | Cause and remedy |
|---|---|
| `lant: command not found` | the environment with lant is not active (`conda activate lant`), or, without conda, `pip install -e .` was not run in the repository |
| `lant check` lists DSurfTomo or fwd1d as missing | run `lant install-external dsurftomo` (needs git, make, gfortran) |
| `conda env create` hangs at "Solving environment" | add `--solver libmamba`, or use mamba / micromamba |

## Configuration

| Symptom | Cause and remedy |
|---|---|
| `configuration key ... unknown` | misspelt key, or a key of an older version; the message says what it became |
| `removed in 1.4` | `dispersion.qc.dev_c` and `smooth_c` were removed (they had no effect in the tests with noise and a second arrival); delete them from the configuration |
| `data.start_date must be YYYY-MM-DD` | write the date as `2026-03-18`; quotes are not needed in YAML |
| parameters changed, not sure where to restart | `lant status` names the first stale step; `lant run all` continues from there |

## clean

| Symptom | Cause and remedy |
|---|---|
| `records without a station in the table` | station names in the headers are not found in the station table. Check the table; when names are written differently, map them with `data.rename` |
| `files hold several records without station codes at the same times` | SEG-Y or SEG-2 files with one trace per channel and no station codes: the file name cannot tell which trace is which station. Convert them to one file per station (or to miniSEED with station codes), e.g. with ObsPy |
| `stations of the table without data` | sites in the table have no records in the time span. For moved nodes check the deployment times in the table |
| deployment times in the table are dates only (e.g. `260913`) | a whole day is taken in the time zone of `data.day_label`, from 0 h of the first day to 24 h of the last; give full date and time (UTC) when this is not accurate |
| `N files gave reader warnings` | the record files contain blocks that are not data (space pre-allocated by the logger, etc.); they are skipped and counted, the data are read normally |
| clean stops for confirmation that all instruments are of one model | without response files the records can be used directly only when all stations are of the same model. After checking, `lant confirm clean`; otherwise give response files in `data.response` |
| not enough disk for the archive | set `clean.write_archive: false`; the later steps then read the raw files (response removal and moved nodes need the archive). clean then reads only the headers and does not test for the 1-s periodic signal; the report marks it as not checked. When the loggers are known to produce it, set `correlation.remove_1s_comb: true` |
| running clean again is slow | station-days with unchanged input files and settings are reused; they are redone only after a change of the `data.*` keys that affect cleaning (sampling, response, day boundaries) |

## cc

| Symptom | Cause and remedy |
|---|---|
| `N processes instead of correlation.nproc` | memory does not hold the requested number of processes; the number is reduced automatically, the results are the same, only slower |
| `One process alone does not fit` | the hourly correlations of all pairs of one day do not fit in memory. Most effective: write the actual range of the area in `data.velocity_range_kms` (the lag window is inversely proportional to the slowest velocity); on large arrays also use `correlation.max_distance_km` to correlate only pairs within a distance |
| `02_ccf` takes a lot of disk | as above: the lag window sets the file size. The clean report gives an estimate |

## stack

| Symptom | Cause and remedy |
|---|---|
| many candidates, mostly "shift" | look at `station_qc.png` first: time shifts of sub-stacks with low correlation are poorly measured. A shift candidate requires the shifts of the pairs of a station to agree (beyond their standard error); a real clock error shows as a consistent shift of one station over several consecutive sub-stacks. With a very wide default velocity range the diagnosis window is long and the band wide, so the fraction of noise is high; narrow `data.velocity_range_kms` and run again |
| you do not want to remove any data | `lant confirm stack` (keep everything); remove the candidates with `--exclude-candidates`; remove particular station-days with `stack.exclude_station_days` |

## reference

| Symptom | Cause and remedy |
|---|---|
| `the F-J ridge could be followed over only N frequencies` | no continuous fundamental ridge was found in the F-J spectrum of the line. Look at `fj.png`: if you can see the fundamental mode, name the band where it is clear in `reference.track_start_band_hz`; if the ridge lies outside your velocity range, change `data.velocity_range_kms`; if neither works, borrow the curve of a nearby array (`reference.source: file`) |
| `only N clear picks that follow each other` (areal array) | too few pairs, or weak surface-wave energy in the band. Check `results/figures/2_ccf_moveout_*.png` for wave trains that move out with distance |
| a higher mode was picked | a curve that is fast overall and rises at high frequencies is probably a higher mode. On areal arrays use `reference.fj_mode: track_fundamental`, narrow the velocity range, or borrow a curve |
| the lower end of the measured band is higher than expected | it follows from the longest pair being `dispersion.qc.lam_c` reference wavelengths; this is the valid range of the two-station method, not an error |

## disp

| Symptom | Cause and remedy |
|---|---|
| few accepted points | `.lant/steps/disp.md` counts the rejected points by reason. Mostly snr: stacking time too short, or no energy in the band; mostly wavelength / path_short: distances too short relative to the wavelength; mostly off_band: the frequency is outside the band excited by the noise sources |
| `median_vs_reference.txt` shows a systematic offset | the reference curve may be wrong (higher mode, cycle skip); go back to reference and check `fj.png` |
| line: the reference covers only a small part of the measured band, and the medians at high frequencies change with the velocity range | beyond the picks of the reference, the cycle count relies on continuation along frequency and the group pick is not guided by the reference, so a fast population can enter. On one real line, with the default range 0.1–5 km/s the reference covered only 1.7–2.8 Hz, and the medians at 3–13 Hz were 4–9 % higher than with the range of the area, 1–4 km/s; the extra points were 8–16 % faster than the points accepted in both runs and did not agree with the F-J ridge. On a line, always give the velocity range of the area, and check in `fj.png` that the picks follow the ridge up to the aliasing line |
| you want to loosen quality control for more points | not recommended. If it is really needed, record the reason with `lant note` and state it in your report |

## inversion and run

| Symptom | Cause and remedy |
|---|---|
| `periods ... more than DSurfTomo allows` | DSurfTomo takes at most 60 periods (pDSurfTomo 80). Reduce `dispersion.per_octave` as suggested, or shorten the exported period range |
| `inversion.sparsity = ... reserves too little` | the sparse reserve given by hand is too small and DSurfTomo would end with a segmentation fault. Use `auto` |
| DSurfTomo `segmentation fault` | look at the end of `work/06_inversion/main/run.log`. Common causes: sparse reserve too small (`inversion.sparsity: auto` estimates it from the grid and the paths), not enough memory (increase `grid_spacing_deg` or reduce `grid_refine`), stack size too small |
| a project is marked invalid by the precheck | under the starting model the rejection rule of DSurfTomo would remove more than half of the data, or the residuals are systematically offset. The starting model does not match the data; check the reference curve and the 1-D fit (`results/figures/4_start_model.png`) |
| the main project runs very long | the time grows with grid nodes, periods and stations. Use `inversion.solver: pdsurftomo` (parallel), or a larger grid spacing |

## evaluate

| Observation | How to read it |
|---|---|
| low checkerboard recovery slope | anomalies of that size are reduced to this fraction in the model; look at larger checkerboard cells for the scales that can be recovered. With few stations and pairs the 3-D inversion changes the starting model only little |
| large difference between the starting-model control and the main model | where the difference is close to the change applied, the data do not constrain the model there; also check whether the residual of the control project has come down to that of the main project (if not, part of the difference is unconverged) |

## Rough estimate of resources

- Disk: size of the raw data (archive) + per day pairs × 24 × (2 × lag window × sampling rate + 1) × 4 bytes (hourly correlations). The clean report gives an estimate.
- Memory: each cc process needs about 2.5 times the size of one day of hourly correlations plus the waveforms of all stations for one day; DSurfTomo grows with grid nodes and travel times.
- Time (for reference, 16-station synthetic example): clean, cc, stack and reference a few seconds to tens of seconds each, disp about half a minute, one DSurfTomo project (10 iterations) about 5–25 minutes.
