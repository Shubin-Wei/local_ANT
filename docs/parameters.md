# Parameters

[中文](parameters.zh.md) | English

Every setting of the configuration file, its default and its basis. The template with the same comments is `lant/config_template.yaml` (`lant init` writes a filled-in copy). "Convention" marks a value that is a choice, not a result of theory or a literature value; change it when you have a reason.

## Required

| Key | Meaning |
|---|---|
| `output_dir` | folder for all results, one sub-folder per step |
| `data.raw_dir` | folder of continuous records (any format ObsPy reads, any layout) |
| `data.station_file` | station table with a header line (`network station longitude latitude [elevation]`), or StationXML |
| `data.velocity_range_kms` | slowest and fastest Rayleigh phase velocity expected (km/s) |
| `clean.fs_hz` | sampling rate after decimation; `lant init` suggests it (rule below) |
| `cc.fmin_hz`, `cc.fmax_hz`, `cc.maxlag_s`, `cc.segment_s` | correlation settings; `lant init` suggests them (rules below) |

## data

| Key | Default | Meaning and basis |
|---|---|---|
| `channel` | `"*Z"` | channel codes of the vertical component (shell pattern) |
| `response` | `null` | response files (StationXML, RESP, dataless SEED; file or folder): responses are removed to velocity, pre-filter flat over `cc.fmin_hz`–`cc.fmax_hz` (so `lant clean` needs the band then). Needed when the instruments differ; with one instrument model the cross-coherence makes it unnecessary (methods §2) |
| `start_date`, `end_date` | `null` | first and last UTC day used (by `lant clean` and `lant cc`: days can be left out of the correlation without cleaning again) |
| `exclude_stations` | `[]` | stations left out (names as in the table; also by both steps) |

## clean (step 1)

| Key | Default | Meaning and basis |
|---|---|---|
| `fs_hz` | init | sampling rate after decimation; a whole fraction of every raw rate. Suggested: the lowest whole-number rate with `cc.fmax_hz` ≤ 0.4 fs (the anti-alias low-pass is at 0.8 × Nyquist) |

## cc (step 2)

| Key | Default | Meaning and basis |
|---|---|---|
| `fmin_hz` | init | lower band edge. Suggested: the longest pair r_max is 2 wavelengths at the fastest velocity, f = 2 v_fast / r_max (two-station far field, Bensen et al. 2007; Luo et al. 2015), lowered by the half-gain width of the Gaussian filter, × (1 − √(ln 2/α)) |
| `fmax_hz` | init | upper band edge. Suggested: spatial sampling theorem, the station spacing at most half the wavelength (e.g. Rost & Thomas 2002), with Δx the median distance to the nearest station: f = v_slow / (2Δx), × (1 + √(ln 2/α)); time sampling: at most 0.4 × the lowest raw rate |
| `maxlag_s` | init | lag window. Suggested: r_max / (0.6 v_slow) + 4 / fmin. 0.6: convention, group velocities in dispersive layers lie well below the phase velocity; 4 periods give the noise window after the arrival |
| `segment_s` | init | segment length. Must exceed 2 × maxlag (circular correlation). Suggested: max(2.5 maxlag, 10 / fmin) |
| `overlap` | 0.5 | overlap of segments (Welch averaging; Seats et al. 2012) |
| `max_over_std` | 10 | a segment whose peak exceeds this × the std of its station-day is not used (NoisePy, Jiang & Denolle 2020); 0: off |
| `substack_hours` | 24 | length of the sub-stacks, the groups of the jackknife error (at least 3 needed). `lant init` suggests shorter ones for records shorter than 3 days |
| `max_distance_km` | `null` | correlate only pairs up to this distance (saves memory and disk for very large arrays) |

## disp (step 3)

| Key | Default | Meaning and basis |
|---|---|---|
| `fmin_hz`, `fmax_hz` | `null` | measurement band; `null`: the correlation band narrowed so the edge filters keep half gain inside it, f_c (1 ∓ √(ln 2/α)) |
| `per_octave` | 12 | measurement frequencies per octave (convention) |
| `reference` | `fj` | `fj` (F-J spectrum of the array), `model` (1-D Vs model), `file` (a curve) — methods §4 |
| `reference_file` | `null` | `model`: `depth_km Vs_km/s`; `file`: `frequency_Hz c_km/s` |
| `fj_band_hz` | `null` | band of the F-J picks (`null`: all frequencies) |
| `alpha` | 5 | width of the Gaussian filters, exp(−α((f′−f)/f)²) (FTAN; convention for local scales) |
| `group_window` | 0.3 | group arrival searched within ±30 % of the reference group velocity (convention) |
| `snr_min` | 5 | lowest SNR of the arrival (literature 5–15) |
| `min_wavelengths` | 2 | shortest path in wavelengths (Bensen et al. 2007: 3; Luo et al. 2015: 1–2) |
| `max_jackknife_rel` | 0.02 | largest jackknife standard deviation / c (convention; `null`: not used) |
| `reference_error` | 0.1 | assumed largest error of the reference; the cycle count is set only where the path is at most 1/(2 × this) wavelengths (theory: half a branch spacing) |
| `slip_cycles` | 0.3 | a phase jump larger than this fraction of a cycle ends the curve (AFTAN) |
| `min_points` | 6 | shortest accepted curve (consecutive frequencies; convention) |

## inv (step 4)

| Key | Default | Meaning and basis |
|---|---|---|
| `data` | `phase` | travel times inverted: `phase`, `group` or `both` (Rayleigh phase and/or group velocities of `3_disp/dispersion.csv`) |
| `grid_km` | `null` | horizontal node spacing; `null`: median measured wavelength / 3 (lateral resolution ≈ one wavelength, three nodes sample it) |
| `depths_km` | `null` | depth nodes; `null`: from 0 in steps of max(shortest λ/4, depth/4) to past half the longest λ (Foti et al. 2018 §4.2.3); the last node is a fixed half-space |
| `start_model_file` | `null` | `depth_km Vs_km/s`; `null`: 1-D inversion of the median curves of the data used |
| `per_octave` | 4 | periods per octave used (a subset of the measured ones) |
| `iterations` | 6 | iterations of the 3-D inversion |
| `weight` | 0.5 | smoothing weight (rms travel-time misfit in s per rms roughness in km/s; data and model terms are means, so the value does not depend on the number of data or nodes; not the weight of a DSurfTomo.in). Convention; in synthetic inversions with Gaussian noise 0.35–0.5 recovered the model best; larger is smoother |
| `damping` | 0.1 | damping towards the starting model (a mean over the parameters, like `weight`); keeps nodes without rays near the starting model. Convention |
| `outlier_mad` | 3 | rejection beyond k × 1.4826 MAD of the residuals of a period, phase and group apart (Leys et al. 2013: 3 "very conservative") |
| `vs_min_kms`, `vs_max_kms` | `null` | Vs bounds; `null`: 0.6 × slowest and 1.4 × fastest Vs of the starting model (convention) |

## checker

| Key | Default | Meaning and basis |
|---|---|---|
| `amplitude` | 0.1 | ±10 % of the starting model |
| `noise` | `auto` | per period, the robust standard deviation of the residuals of the 3-D model (s) (Rawlinson & Spakman 2016: noise like the data); or a number in seconds |

## plots, cores

| Key | Default | Meaning |
|---|---|---|
| `plots.relief` | `auto` | relief under the maps: `auto` GMT's `@earth_relief_03s` (needs internet; maps without relief when not available), a grid file, or `null` |
| `plots.vs_range_kms` | `null` | colour range of Vs (`null`: 2–98 % of the model under the stations) |
| `plots.slice_depths_km` | `null` | depths of the map slices (`null`: up to 8 depth nodes) |
| `plots.dpi` | 300 | resolution of the PNG figures |
| `cores` | `auto` | `auto`: 80 % of the idle CPUs (others can still work on a shared machine); or a number |
