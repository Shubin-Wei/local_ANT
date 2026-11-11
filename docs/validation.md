# Validation: synthetic tests and their numbers

[中文](validation.zh.md) | English

Numbers only; what they mean for your data is yours to judge. Sections 1–11 use synthetic data with a known answer, each with the version it was run with; section 12 reports trial runs on real data; section 13 lists what has not been tested.

## 1. How the test data are made

| Generator | What it makes | What it can show |
|---|---|---|
| `examples/synthetic/make_synthetic.py` | raw noise records: plane Rayleigh waves of a 1-D model from random azimuths (16 stations, 3 days, aperture 9 km; `--layout line`: 24 stations, 2 days) | the whole chain from raw files; laterally uniform, isotropic, fundamental mode only |
| `examples/synthetic/make_stack.py` | the correlations themselves, J0(2πf t(f)), as a stack file: laterally varying medium (chequer pattern), a second faster arrival, incoherent noise | measurement and inversion with lateral variation, a higher mode, noise |

Comparison scripts: `compare_truth.py` (uniform medium), `compare_hetero.py` (laterally varying medium). The tests in `tests/test_numerics.py` repeat the benchmark of section 2 in a small form.

## 2. Two-station phase velocity on exact J0 correlations

Correlations J0(2πf d / c(f)) of a uniform medium, no noise; relative error of the measured phase velocity by path length in wavelengths:

| path length | median error | 95th percentile of the absolute error |
|---|---|---|
| 2–3 λ | 0.007 % | 2.9 % |
| 3–5 λ | 0.064 % | 0.83 % |
| 5–10 λ | 0.096 % | 0.20 % |
| ≥ 10 λ | 0.085 % | 0.11 % |

The default `dispersion.qc.lam_c` accepts paths from 2 wavelengths.

## 3. Phase reading: phase-matched filter against the earlier instantaneous-frequency correction

Same stack of the uniform example (`make_synthetic.py`), same quality control:

| | earlier method | phase-matched filter |
|---|---|---|
| accepted samples | 5085 | 4955 |
| median difference to the true phase velocity | −0.25 % | +0.07 % |
| at 1.26 Hz | −0.75 % | +0.07 % |
| at 0.89 Hz | −0.83 % | −0.23 % |
| samples above the 12-Hz upper end of the source band | 351 | 234 |

(Version 1.3.0.)

**Band-edge criterion (redefined in version 1.4.0).** The criterion of 1.3.0 used the amplitude spectrum of each pair and did not act on data with a high signal-to-noise ratio: in one day of synthetic data (noise sources cut off at 12 Hz) all 120 pairs were accepted at 12.7, 13.45 and 14.25 Hz. Version 1.4.0 uses the amplitude spectrum of the coherent signal ([methods.md](methods.md), section 7.6):

| frequency (Hz) | amplitude ratio (one day / three days of data) | pairs accepted (one day of data, 1.3.0 → 1.4.0) |
|---|---|---|
| 12.7 | 0.64 / 0.66 | 120 → 120 |
| 13.45 | 0.37 / 0.33 | 120 → 0 |
| 14.25 | 0.25 / 0.19 | 120 → 0 |

With half of the pairs the result is the same.

## 4. Reference curve: frequency-Bessel spectrum computed in lant

Uniform example: the picks differ from those of CC-FJpy (used up to v1.2) by at most 0.39 %; median difference to the true curve +0.30 % (CC-FJpy +0.24 %). With the default velocity range 0.1–5 km/s (scan 0.05–8 km/s, two passes) 43 picks were used between 0.27 and 3.0 Hz; with the range 1.5–3.5 km/s, 40 picks between 0.33 and 3.2 Hz (version 1.3.0).

New examples of version 1.4.0: areal (range 1.5–3.5 km/s) 41 picks at 0.28–3.17 Hz, median difference to the true curve +0.28 %; line (`track_fundamental`) 36 picks at 0.35–2.83 Hz, median +0.15 %, at most 3.5 % (the lowest frequencies). The line picks stop at 2.83 Hz because of the limit of wavelengths to at least twice the mean station spacing (0.775 km for this line); the measurement still reaches 12.7 Hz (section 5), where above the picks the number of cycles is carried along frequency.

## 5. The whole chain on a uniform medium

`make_synthetic.py`, then `lant run all`, then `compare_truth.py`. "Difference" is the median over the pairs against the true curve per frequency; Vs is the horizontal mean of the model against the true model at the depth nodes.

| | areal, range 1.5–3.5 km/s (1.4.0) | line, 24 stations (1.4.0) | areal, default range 0.1–5 km/s (1.3.0) |
|---|---|---|---|
| pairs / curves exported | 120 / 118 | 276 / 215 | 120 / 116 |
| samples accepted by the measurement / points exported | 5908 / 3139 | 9505 / 8039 | 5053 / — |
| phase velocity: median difference over all frequencies | +0.10 % | +0.02 % | +0.04 % |
| phase velocity: largest difference at one frequency | 0.94 % | 0.76 % | 0.79 % |
| Vs at the depth nodes: range of the differences | −0.4 … +0.3 % | −0.6 … +0.2 % | −0.1 … +0.2 % |
| travel-time residual, final (RMS; relative) | 5.6 ms; 0.38 % | 5.0 ms; 0.40 % | — |
| share of data kept by DSurfTomo's rejection rule | 0.86 | 0.86 | — |
| checkerboard (cells of 3 km): noise derived from the main fit; recovery slope over the depths; pattern correlation | 0.11 %; 0.01–0.16; 0.36–0.54 (fixed 2 %: 0.01–0.22; 0.42–0.56) | 0.07 %; 0.00–0.28; 0.21–0.70 (fixed 2 %: 0.00–0.27; 0.28–0.69) | — |
| starting-model control against the main model (RMS / largest) | 0.099 / 0.246 km/s | 0.102 / 0.265 km/s | — |
| station diagnosis: candidates | 0 | 1 | 0 |
| seconds: clean / cc / stack / reference / disp | 6 / 24 / 2 / 5 / 51 | 4 / 29 / 3 / 6 / 64 | 6 / 60 / 73 / 7 / 119 |
| seconds: inversion / run / evaluate / figures | 10 / 1930 + 1786 / 597 / 71 | 12 / 14008 + 8640 / 986 / 41 | — |

The times are from a fully loaded shared server (40 cores, load average 30–100) and give only the order of magnitude. "run" is given in two parts: main project, control and checkerboards with a fixed noise of 2 %, then the checkerboards again with the noise derived from the main fit (`inversion.noise: auto`, the default since then); evaluate and figures are of the second run. Making the data (`make_synthetic.py`) took 33–34 minutes under that load.

Starting-model control of the areal example: the model difference is close to the ±5 % gradient applied; the weighted residual of the control fell from 38.5 to 10.0 ms, that of the main project from 2.7 to 2.0 ms. The control had not converged to the level of the main project, so part of the difference is the starting model itself; 16 stations and 3 days constrain the model only weakly.

In these examples the starting model (1-D inversion of the average curve) already agrees with the true model within the numbers above, and the 3-D inversion changes it by less than 0.002 km/s; they test the measurement and the 1-D step, not lateral resolution.

The quick example (`make_synthetic.py --quick`: 10 stations, 6 hours) took about 36 s to make on the same server and about 6 minutes for the whole chain; one run gave a median phase-velocity difference of +0.33 % and at most 1.35 %.

## 6. Laterally varying medium

(Version 1.3.0.) `make_stack.py --hetero 0.06 --cell 3 --noise 0.02` on the geometry of the areal example (16 stations): Vs pattern of ±6 % with 3-km cells, the same at all depths; then the chain from the reference step, and `compare_hetero.py`.

Phase velocities against the path average of the truth (accepted samples):

| band (Hz) | samples | median difference | 95th percentile of the absolute difference |
|---|---|---|---|
| 0.63–1.37 | 523 | −0.02 % | 0.43 % |
| 1.37–3.00 | 1406 | +0.03 % | 0.63 % |
| 3.00–6.54 | 1508 | +0.07 % | 0.35 % |
| 6.54–14.3 | 1652 | +0.02 % | 0.09 % |
| all | 5089 | +0.03 % | 0.30 % |

Path anomaly relative to the uniform background: standard deviation 1.64 % (truth) and 1.66 % (measured), correlation 0.983, regression slope 0.994.

Model against the pattern, at the 552 grid nodes within one cell of a station, and the checkerboard test of the same run (cells of 3 km, ±10 %) for comparison:

| depth (km) | model: correlation | model: slope | checkerboard: correlation | checkerboard: slope |
|---|---|---|---|---|
| 0.00 | 0.49 | 0.08 | 0.54 | 0.09 |
| 0.12 | 0.50 | 0.17 | 0.57 | 0.20 |
| 0.24 | 0.50 | 0.21 | 0.57 | 0.24 |
| 0.36 | 0.50 | 0.20 | 0.57 | 0.23 |
| 0.48 | 0.49 | 0.16 | 0.57 | 0.18 |
| 0.60 | 0.48 | 0.12 | 0.56 | 0.14 |
| 0.75 | 0.46 | 0.09 | 0.55 | 0.10 |
| 0.94 | 0.45 | 0.06 | 0.54 | 0.07 |
| 1.20 | 0.43 | 0.03 | 0.52 | 0.04 |
| 1.50 | 0.40 | 0.01 | 0.47 | 0.01 |

Travel-time residual: 34.3 ms (1.69 %) for the starting model, 18.9 ms (0.98 %) for the final model; share kept by the rejection rule 0.89. The straight-ray travel times of the test data and the bent rays of DSurfTomo differ; this is part of the numbers.

DSurfTomo's rejection threshold (`inversion.threshold`, in interquartile ranges; default 3) on the same data, model against the pattern:

| threshold | correlation at 0.24 km | slope at 0.24 km | slope at 0.60 km | standard deviation of the weighted residual after the last iteration (ms) |
|---|---|---|---|---|
| 3 | 0.50 | 0.21 | 0.12 | 12.4 |
| 4 | 0.52 | 0.25 | 0.14 | 13.0 |
| 5 | 0.53 | 0.27 | 0.16 | 13.4 |

## 7. Quality-control screens with noise and a second arrival

(Version 1.3.0; this section is the basis for removing `dev_c` and `smooth_c` in 1.4.0.) `make_stack.py` on the same geometry, uniform medium, with incoherent noise (`--noise`, relative to the largest amplitude of each correlation) and a second arrival at 1.25 times the phase velocity above 4 Hz (`--higher-mode`, amplitude relative to the fundamental). Errors are against the true phase velocity; "> 2 %" and "> 5 %" count accepted samples.

| case | screens | accepted | median abs. error | > 2 % | > 5 % | largest |
|---|---|---|---|---|---|---|
| noise 0.05 | all | 5016 | 0.07 % | 38 | 0 | 4.9 % |
| | without `dev_c`, `smooth_c`, `ucons_c` | 5016 | 0.07 % | 38 | 0 | 4.9 % |
| | also without `dev_pre`, `side_c` | 5096 | 0.08 % | 55 | 0 | 4.9 % |
| noise 0.15 | all | 3949 | 0.20 % | 84 | 3 | 6.1 % |
| | without the five | 3949 | 0.20 % | 84 | 3 | 6.1 % |
| second arrival 0.8, noise 0.05 | all | 5069 | 0.09 % | 138 | 5 | 5.9 % |
| | without `dev_c`, `smooth_c`, `ucons_c` | 5080 | 0.09 % | 149 | 6 | 5.9 % |
| | also without `dev_pre`, `side_c` | 5156 | 0.09 % | 179 | 8 | 6.9 % |
| second arrival 1.5, noise 0.05 | all | 2576 | 0.22 % | 99 | 29 | 32.5 % |
| | without `dev_c`, `smooth_c`, `ucons_c` | 2579 | 0.22 % | 101 | 29 | 32.5 % |
| | also without `dev_pre`, `side_c` | 3491 | 0.17 % | 229 | 118 | 32.5 % |

In these four cases `dev_c` and `smooth_c` alone changed nothing; `ucons_c` removed 3–11 samples. With a second arrival stronger than the fundamental, 29 accepted samples are more than 5 % off with all screens on.

## 8. Test for the 1-s periodic signal: false and true detections

The fold F test of `clean` ([methods.md](methods.md), section 3) at a significance of 0.1 %, on two sets of synthetic noise without any pulse:

| noise | sampling rate (Hz) | set A (20 records of 6 h each): 1.3.0 → 1.4.0 | set B (200 records each): 1.3.0 → 1.4.0 |
|---|---|---|---|
| white | 100 / 250 | 0 → 0 / 1 → 1 | 0 → 0 / 0 → 0 |
| red | 100 / 250 | 0 → 0 / 2 → 0 | 3 → 0 / 1 → 0 |
| microseism-like | 100 / 250 | 0 → 0 / 0 → 0 | 0 → 0 / 0 → 0 |
| band-limited | 100 / 250 | 0 → 0 / 3 → 0 | 2 → 0 / 14 → 0 |

For white noise the F test is exact and both versions use the same limit; the one false detection in set A at 250 Hz is the same record in both versions. A pulse of 0.5 times the noise standard deviation once per second: detected in 200 of 200 records of set B (both versions).

## 9. Inter-station distances

Version 1.4.0 computes distances on the WGS84 ellipsoid (Vincenty 1975). The ratio to the distance on a sphere of radius 6371 km depends on azimuth: at latitude 30° +0.30 % north–south and −0.18 % east–west; near the equator up to +0.56 % north–south. On 200 random pairs the relative difference to ObsPy `gps2dist_azimuth` is below 10⁻⁸ (`tests/test_numerics.py`).

## 10. Station diagnosis: injected clock errors

Version 1.4.0. In the stack file of the three-day areal example, every pair of one station in one sub-stack is shifted by some samples (sign by the order of the pair, the full stack updated accordingly), and the station diagnosis is run. Six trials per shift (different stations and sub-stacks):

| shift | detected | other candidates |
|---|---|---|
| 0.5 samples (0.01 s) | 6/6 | 0 in every trial |
| 1 sample | 6/6 | 0 |
| 2 samples | 6/6 | 0 |
| 5 samples | 6/6 | 0 |

Without an injected shift: 0 candidates.

## 11. Instrument responses, SEG-Y input, directional noise

Version 1.4.0, quick example (`make_synthetic.py --quick`: 10 stations, 6 hours), whole chain.

**Two sensor types** (`--responses`): alternate stations recorded through geophones of 2 Hz and 10 Hz natural frequency (damping 0.7), responses given as StationXML.

| | responses removed (`data.response`) | responses not removed |
|---|---|---|
| curves exported: pairs of different / of the same sensor type | 23 of 25 / 20 of 20 | 15 of 25 / 20 of 20 |
| phase velocity: median difference to the true curve; largest at one frequency | +0.39 %; 2.16 % | +0.21 %; 2.23 % |

Without removal the velocities are hardly biased: the phase difference of two sensors shifts the causal and the acausal branch of their correlation the same way in time, and the symmetric component (causal plus time-reversed acausal) cancels it to first order. The pairs with different sensors lose curves instead (15 of 25 against 23 of 25). With noise from one side only this cancellation fails; the rule that stations without response files must be of one model stays.

**SEG-Y without station codes** (unit test `test_segy_input`): three stations, one file each, two hours cut into traces of 300 s (SEG-Y holds at most 32 767 samples per trace), plus one file with three traces at the same time. The 72 traces are named from their file names and joined to 6 station-hours, with and without the archive; the file with three traces at the same time is reported and not used. Before this test such records were taken as one unnamed station and none was cleaned (history.md, v1.4).

**Directional noise** (`--directional F`: a fraction F of the plane waves from azimuth 45° ± 15°, the rest from random azimuths; one run each, 6 hours of record):

| | isotropic (F = 0) | F = 0.4 | F = 0.7 |
|---|---|---|---|
| F-J reference: picks; band | (section 5) | 35; 0.30–2.12 Hz | 11; 0.35–0.63 Hz |
| curves exported | 42 | 36 | 0 |
| phase velocity: median difference to the true curve; largest at one frequency | +0.33 %; 1.62 % | +0.57 %; 2.19 % | — |
| median difference of the pairs at 0–30° / 30–60° / 60–90° from the noise direction | +0.25 / +0.33 / +0.27 % | +0.27 / +0.78 / +0.26 % | — |

With F = 0.7 the frequency-Bessel spectrum gives a ridge only at the lowest frequencies, no number of cycles can be fixed in the measurement band, and nothing is exported: the chain gives no curve rather than doubtful ones. The effect of directional noise on the velocities that are exported depends on the source distribution and was tested for this one geometry only.

## 12. Trial runs on real data

Version 1.4.0, run from the raw records on two unpublished nodal arrays (results kept outside the repository; data and station names are not given here).

**A line**: 55 stations, 15.4 km long, 6 days, 200 Hz vertical component, miniSEED.

| | velocity range left at the default 0.1–5 km/s | velocity range given as 1–4 km/s |
|---|---|---|
| clean | 457 s | (the same cleaning) |
| cc | 943 s (lag window ±600 s, 48 GB of hourly files) | 1015 s |
| stack; station-diagnosis candidates (1.3.0 rule → 1.4.0 rule) | 402 s; 78 → 3 | 83 s; 78 → 4 |
| reference: band of the F-J picks | 1.68–2.83 Hz | 1.06–3.17 Hz |
| disp: curves / points exported | 334 / 5620 | 346 / 5665 |
| inversion: run / evaluate | stopped part-way (see below) | 7995 s (5 projects side by side, 10 iterations; checkerboards at a fixed 2 %), then 7294 s for the 3 checkerboards with the derived noise / 1003 s |
| travel-time residual, starting → final model (RMS; relative); share kept | — | 130.3 → 102.2 ms; 12.3 → 8.6 %; 0.88 |
| checkerboards 2.5 / 3.7 / 5.3 km: noise derived from the main fit; largest recovery slope; pattern correlation | — | 4.81 % (relative RMS of the 4943 of 5623 residuals kept); 0.30 / 0.33 / 0.21; 0.28–0.63 (fixed 2 %: 0.31 / 0.34 / 0.21; 0.31–0.62) |
| starting-model control against the main model (RMS / largest) | — | 0.085 / 0.240 km/s |

- Station diagnosis: under the 1.3.0 rule nearly all 78 candidates were poorly measured time shifts of sub-stacks with low correlation (60 % of the sub-stacks with a correlation below 0.2 were beyond the limit). With the standard-error condition of 1.4.0 what remains are shifts of 0.29–0.61 s that agree between the pairs of a station; two of them appear under both settings.
- Reference curve: with 1–4 km/s, the tracking rule of 1.3.0 broke where the dispersion is steep (a fall of 7–8 % per 1/12 octave) and continued on an aliased ridge about 15 % too slow; the rule of 1.4.0 (slope limit widening with the frequency step, strongest peak, wavelength at least twice the mean station spacing) follows the continuous ridge of the image up to 3.17 Hz.
- Dispersion: of the 2753 points accepted under both settings, 97–99 % have the same number of cycles, with a median residual of 0.02–0.04 cycles. With the default range the reference covers only 1.68–2.83 Hz; the extra points accepted there are 7 / 14 / 16 % faster than the common points at 3.2 / 5 / 8 Hz and do not agree with the F-J ridge at 3.2 Hz, so the medians are 4–9 % higher than with the range of the area. A line needs the velocity range of the area ([troubleshooting.md](troubleshooting.md)).
- Checkerboards of the line: the pattern is recovered only between about 4.5 and 10.5 km along the line, where the stations are dense; at the sparse ends (stations 0.7–1.2 km apart, against 0.08–0.3 km in the middle) nothing is recovered. The slope and correlation are taken over all nodes within 0.3 km of a station, the ends included, which is why they are lower than inside the dense part. The figures are sections along the line (`6_checkerboard_*`).
- On the fully loaded server each pre-check of the default-range project took about 19 minutes, and its inversion was stopped part-way; the inversion was completed only for the project with the velocity range.

**An areal array**: 166 stations, 8 days, 200 Hz, without the archive (raw files read directly). clean 245 s, cc 3954 s (13695 pairs, ±40 s), stack 210 s, reference 203 s, disp 677 s; 4647 curves exported. The correlation candidates of the station diagnosis concentrate on 4 stations, the same ones found abnormal earlier and independently from their noise power spectra. The F-J spectrum shows one continuous fundamental-mode ridge. Against an existing, independent processing of this array (from a longer record), point by point on the common pairs, the median difference of the phase velocities is about +0.1 %; the detailed numbers are not given because the data are unpublished.

## 13. Not tested

- Response removal on real data (both trial data sets are of one instrument model, without response files; section 11 uses synthetic records); moved nodes are covered by unit tests on constructed records only.
- The 1-s periodic signal is not tested when the archive is not written (the report marks it as not checked); the areal trial array is known to carry such a signal, which was not removed in this run.
- SEG-2 input (ObsPy cannot write SEG-2, and no long passive recording was at hand; it goes through the same reading path as SEG-Y); arrays of thousands of stations (the pair-distance limit and the memory estimate were used on small and medium arrays only).
- Directional noise on real data, and a real higher mode with its own dispersion (the second arrival of section 7 is the fundamental curve scaled by 1.25).
- Zone references (`reference.zone_km`) on a medium that needs them.
- The run-through with an AI assistant described in ai_workflow.md.
