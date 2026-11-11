# Methods: principles, formulas and sources

[中文](methods.zh.md) | English

This document explains, step by step, what `lant` computes, why it computes it that way, which publication or theoretical relation it rests on, and under which conditions each choice holds. For how to run the steps see [tutorial.md](tutorial.md); for the default and the source category of every parameter see [parameters.md](parameters.md); for the numbers of the synthetic tests see [validation.md](validation.md).

"Empirical" in this text marks a number for which no publication gives a value and which this project chose on synthetic data; every such number is labelled in `parameters.md`.

---

## 1. Overview

Ambient-noise imaging rests on the relation between the cross-correlation of two records in a nearly isotropic, diffuse noise field and the Green's function between the two stations; its phase gives the surface-wave travel time between them. On the vertical component the strongest arrival in the correlation is the fundamental-mode Rayleigh wave. The workflow therefore has three parts:

| Part | Steps | Result |
|---|---|---|
| From records to empirical Green's functions | clean → cc → stack | the stacked correlation of each pair, and sub-stacks over groups of hours |
| From Green's functions to dispersion | reference → disp | one fundamental-mode Rayleigh phase-velocity curve per pair, with jackknife errors |
| From dispersion to a model | inversion → run → evaluate | a 3-D shear-velocity model, travel-time residuals, resolution tests |

The workflow measures only fundamental-mode Rayleigh phase velocities and uses only the vertical component. Group velocities serve only to locate the time at which the phase is read; they are not a product.

---

## 2. Stations and distances

### 2.1 Matching stations and records

The station table is matched to the record headers in this order: same network and station; same station code; network and station written together; equal after removing leading letters or leading zeros. Ambiguous matches are not used and are listed in the report. When nodal instruments were moved in the field, the table gives the serial number (or the station name is the code in the headers) and the deployment start and end; records are assigned to sites by these time windows. When only dates are given, a whole day is taken in the time zone of `data.day_label` (from 0 h of the first day to 24 h of the last day).

### 2.2 Distances

Station coordinates are geodetic (longitude and latitude on the WGS84 ellipsoid). Distances and azimuths between stations are computed on the WGS84 ellipsoid with the inverse solution of Vincenty (1975) (`lant.signal.geodesic_km`). On a sphere of radius 6371 km the distance error depends on the azimuth of the path: at latitude 30° north–south paths come out 0.30 % too long and east–west paths 0.18 % too short; near the equator the error reaches +0.56 %. This error is of the same size as the measurement precision and has the form of a 2θ azimuthal anisotropy, so ellipsoidal distances are required.

Phase velocities are measured with ellipsoidal distances: c = Δ / t. DSurfTomo turns the input phase velocities back into travel times with spherical distances and traces rays on its sphere. Because the input is a velocity and not a travel time, in a uniform medium the model velocity of DSurfTomo equals the input velocity, without an azimuthal bias; its internal geometry is stretched relative to the true ellipsoid by at most 0.6 %, which moves features by far less than a grid spacing. The source code of DSurfTomo is not changed.

---

## 3. clean: preparing the raw records

For each station and each processing day (24 hours from `data.day_start_utc_hour`):

1. **Sampling grid**: the record start is placed on the UTC grid of whole samples of its rate. Offsets up to 0.05 samples are shifted, larger ones are Lanczos-interpolated; both are counted.
2. **Overlaps and duplicates**: exact duplicates and overlaps with identical values are merged; overlaps with different values are not resolved by choosing one, the overlap becomes a gap.
3. **Gaps**: gaps stay gaps; they are neither interpolated nor zero-filled. The later segmentation uses only segments without gaps.
4. **Instrument response** (only when `data.response` is given): removed to velocity as in Bensen et al. (2007), with a pre-filter band derived from the correlation band. Without response files the records are used as they are; this is valid only when all stations are of the same model (cross-coherence removes the amplitude response, and the phase responses of instruments of one model cancel between the two stations of a pair).
5. **Test for a 1-s periodic signal**: some data loggers produce coherent noise synchronised with UTC once per second (Takagi et al. 2015); it creates a spurious coherent signal near zero lag in the correlations. The test folds the record by UTC seconds: with p samples per second and M seconds, each second has its mean removed, and the samples are averaged phase by phase into the folded waveform s(j). Without such a signal the variance of the folded waveform is 1/M of the record variance. The test statistic

   F = [M Σⱼ s(j)² / (p − 1)] / [(Σ residual²) / ((p − 1)(M − 1))]

   has expectation 1 without a signal and is compared with the 0.1 % quantile of the F distribution with (ν, ν(M − 1)) degrees of freedom. For white noise ν = p − 1; in real noise neighbouring samples are correlated and there are fewer independent values, so the effective number of degrees of freedom ν = (p − 1) / Σₖ ρ(k)² is used, with ρ the autocorrelation within one second (the effective degrees of freedom of a sum of squares of correlated variables). Without this correction red and band-limited noise sampled at 250 Hz give false detections in 7–15 % of the records (nominal 0.1 %); with it, 200 trials of each noise type gave no false detection, and a pulse of 0.5 σ once per second was detected in all 200 trials ([validation.md](validation.md), section 8). When the signal is detected on more than half of the processing days of a station, the cc step subtracts the folded waveform hour by hour from that station.

The cleaned day files are written to `work/01_clean/archive` (not when `clean.write_archive: false`; the later steps then read the raw files). Each station-day records a fingerprint of its input files and settings; on a later run, station-days whose inputs and settings are unchanged are reused.

---

## 4. cc: hourly cross-coherence

### 4.1 Preprocessing and segments

Each station-day is low-pass filtered with zero phase (8th-order Butterworth, corner at 0.8 of the new Nyquist frequency) and decimated to `correlation.fs`. It is cut into segments of length L with step L/2 (50 % overlap, Seats et al. 2012). A segment is dropped if it contains a gap or if its absolute peak exceeds `max_over_std` times the standard deviation of the station-day (the NoisePy practice, Jiang & Denolle 2020). Each segment is detrended and tapered with a 5 % cosine window. Records of the neighbouring days are read as well, so that segments crossing a day boundary of the archive are not lost.

### 4.2 Cross-coherence

The spectra of the two stations are each divided by their own amplitude and then multiplied:

C_ij(ω) = F_i*(ω) F_j(ω) / (|F_i(ω)| |F_j(ω)|)

This is the cross-coherence (recommended by Xu et al. 2021 after comparing several normalisations). Only `fmin_hz`–`fmax_hz` is kept, with a cosine transition of a quarter octave at each end. The segments of one hour are averaged in the frequency domain and transformed back; lags within ±`maxlag_s` are kept. Positive lag means energy travelling from station i to station j, with i and j sorted by station name.

### 4.3 Parameters derived from the velocity range

The user gives only the approximate range of phase velocities in the area, [v_slow, v_fast] (`data.velocity_range_kms`). The following quantities follow from theoretical relations (the numbers are written in each step report):

| Quantity | Relation | Reason |
|---|---|---|
| Lowest frequency f_lo | longest pair = 2 wavelengths (at v_slow), widened to the half-gain point of the Gaussian filter at that frequency | the two-station method needs paths of at least about 2 wavelengths (Luo et al. 2015) |
| Lag window | longest distance / (0.6 v_slow) + 4 longest periods | group velocities can be lower than phase velocities; two periods after the arrival are kept for the signal-to-noise ratio |
| Segment length | max(2 × lag window, 10 longest periods) | frequency-domain correlation is circular; segments at least twice the lag window do not wrap around |
| Sampling rate | the highest whole number not above `fs_max` to which every raw rate can be decimated by an integer factor | decimation without interpolation |
| Highest measured frequency | 0.4 fs / (1 + √(ln 2 / α)) | the half-gain point of the Gaussian filter at the highest frequency stays inside the anti-alias pass band |

The lag window and the segment length grow inversely with v_slow, and so do the computation and the disk space. The default range 0.1–5 km/s covers soft soil to hard rock, but on an array with an aperture of more than ten kilometres it gives a lag window of 600 s. The clean report gives an estimate of the memory per process and of the disk space of the hourly files, and the cc step reduces the number of processes to what fits in memory. Writing a narrower range when the approximate velocities of the area are known is the most effective way to shorten the computation.

---

## 5. stack: stacking and station diagnosis

### 5.1 Stacking

All hours are averaged with weights equal to their numbers of valid segments, giving the stacked correlation; sub-stacks over groups of `stack.block_hours` hours are stored for the jackknife errors. `block_hours` is the longest of 24, 12, 8, 6, 4, 3, 2 and 1 hours for which the median station has at least 8 sub-stacks (empirical).

### 5.2 Station diagnosis

Clock errors, moved instruments or degraded coupling show up as a change of all pairs of one station over some time (Stehly et al. 2007; Sens-Schönfelder 2008). The diagnosis:

1. Each sub-stack of a pair is compared with the stack of the other sub-stacks of that pair (leave one out, so that no sub-stack is compared with itself): correlation coefficient r and time shift τ (peak of the cross-correlation, parabolic interpolation between samples), within |t| ≤ distance / (0.6 v_slow).
2. For each station and sub-stack, the median over its pairs with the nearest `partners` stations is taken.
3. r is Fisher-transformed, z = atanh(r) (approximately normal), and the median of the array for the same sub-stack is subtracted (the noise level changes with time).
4. The scale is the median absolute deviation, σ = 1.4826 MAD (Leys et al. 2013). The limit follows from Chauvenet's criterion for N values: z is tested on the low side only, with a limit of Φ⁻¹(1 − 0.5/N) σ; the time shift is tested on both sides, with Φ⁻¹(1 − 0.25/N) σ and at least 0.25 samples. By construction about 0.5 normal values per quantity exceed the limit on average.
5. A time-shift candidate needs one more condition: the uncertainty of the median shift of the station. In sub-stacks with low correlation the peak position is set mostly by noise, and the shifts of the pairs of one station disagree; a real clock error moves all pairs of the station together. The standard error of the median over the pairs of the station, se = 1.2533 × 1.4826 MAD / √n (the standard error of the median of a normal sample is √(π/2) σ/√n), is therefore computed, and the difference between the median shift and the array median must also exceed Φ⁻¹(1 − 0.25/N) × se.

The diagnosis lists candidates and does not remove anything. When there are more candidates than the number expected by chance under the criterion (0.5 each for correlation and time shift, 1 in total), the chain stops for your judgement; otherwise they are only listed in the report.

---

## 6. reference: the reference curve

The reference curve has three uses: it bounds the group-time picking; it fixes the integer number of cycles of the phase on short paths; it gives the lower end of the measured band and the velocity axis. It does not enter the final phase-velocity values.

### 6.1 Frequency–Bessel spectrum

For each frequency f the real part of the correlation spectra of the pairs, C(r, ω), gives

I(f, c) = ∫ C(r, ω) J₀(ω r / c) r dr

(Wang et al. 2019, eq. 4), integrated with the trapezoidal rule over pairs sorted by distance, in double precision. At each frequency only pairs up to `fj_rmax_lambda` wavelengths are used: a first pass uses all pairs and picks a curve, the second converts the limit into distances with that curve. Ridges of energy in the spectrum mark the phase velocities of the modes.

### 6.2 Picking

| Mode | Used for | Procedure |
|---|---|---|
| `max_peak` | areal arrays | the strongest peak at each frequency; a pick is "clear" when the second peak is below 0.7 of the strongest and enough pairs contribute; the longest run of clear picks that follow each other is fitted with a cubic polynomial in log f, and end values are used outside the picks |
| `track_fundamental` | lines | the F-J spectrum of a line has pairs of one direction only, and higher modes or obliquely incident noise are often stronger. The ridge is followed in frequency from a start: at each frequency the strongest peak within the allowed window is taken (the ridge is the energy maximum); the window limits the slope d ln c / d ln f, widening by 5 % for every 1/12-octave step since the last pick, and the velocity may rise with frequency by at most 3 %; picks are restricted to `data.velocity_range_kms` and to wavelengths of at least twice the mean station spacing (a line samples the wavefield at its station spacing, and shorter wavelengths produce aliased ridges in the F-J spectrum; for irregular sampling the average sampling density is what counts, Landau 1967). The start is chosen automatically: a start band every 1/3 octave is tried, and the one that can be followed over the most frequencies is kept; on a tie the slower one (at one frequency the fundamental mode is the slowest). With fewer than 5 frequencies the step stops with an error and gives three ways forward: name the start band, narrow the velocity range, or borrow a curve |

Whether the fundamental mode was picked is the most important manual judgement of the whole workflow: when a higher mode is picked, the number of cycles is off by one everywhere and all phase velocities are biased high. The reference step therefore always stops for you to check `fj.png`.

### 6.3 Borrowed curves

`reference.source: file` takes a curve from elsewhere (an npz written by lant, or a two-column text file "frequency in Hz, phase velocity in km/s"). It is interpolated in log f only within the valid part of the source curve, not across invalid parts and not beyond the ends.

---

## 7. disp: two-station phase velocities

### 7.1 Symmetric component and narrow-band filters

The positive and negative lags are averaged (symmetric component, Bensen et al. 2007) and filtered with Gaussian narrow-band filters at logarithmically spaced frequencies (12 per octave)

G(ω) = exp[−α ((ω − ω_c) / ω_c)²], α = 5,

giving the envelope at each frequency. The envelope divided by the root mean square of the noise after the signal window gives the signal-to-noise image.

### 7.2 Group times

The reference group velocity follows from the reference phase velocity, U = c / (1 − d ln c / d ln f) (Bensen et al. 2007, eq. 7). The envelope ridge is picked by dynamic programming within 0.65–1.35 U: departures of more than 25 % from the reference are penalised, and changes between neighbouring frequencies are smoothed with a 3 % constraint.

### 7.3 Phase: phase-matched filtering

Reading the phase of the narrow-band signal at the group time at the centre frequency is biased when the group delay changes within the filter band (when there is dispersion); the bias has the same origin as the systematic group-velocity error shown by Shapiro & Singh (1999). This workflow uses phase-matched filtering (Herrin & Goforth 1977; as in FTAN, Levshin & Ritzwoller 2001):

Ψ(ω) = ∫ t_g(ω) dω,  Φ(ω_c) = arg Σ_ω X(ω) G(ω) exp{i [Ψ(ω) − Ψ(ω_c)]},

with t_g the measured group delay. Multiplying by exp(iΨ) removes the dispersion of the arrival, so that the phases of the frequencies within the band line up. The phase travel time is

t_p = −(Φ − φ₀) / ω_c + N T,  c = Δ / t_p,

with φ₀ = π/4 the far-field phase of the 2-D Green's function (Lin et al. 2008, eq. 7), N the integer number of cycles, T the period and Δ the ellipsoidal distance. On exact J₀ correlations the median error of this procedure is within 0.1 % ([validation.md](validation.md), section 2).

### 7.4 Number of cycles

The number of cycles is fixed in the band where the path is `dispersion.lam_min`–`anchor_max_lambda` (default 1–3.5) reference wavelengths long: on such paths even a 10 % error of the reference curve gives less than half a period of travel-time error. It is then carried along frequency on a grid four times finer, predicting the phase at the next frequency from the measured group slowness (the phtovel procedure of AFTAN); a departure of more than 0.3 cycles from the prediction counts as a cycle skip and stops the continuation.

### 7.5 Errors

The sub-stacks are split into `njack` groups and the measurement is repeated with one group left out each time (delete-a-group jackknife):

σ² = (g − 1)/g · Σᵢ (cᵢ − c̄)², g the number of groups.

It reflects only the repeatability of the measurement in time; it does not include errors of the reference curve, of the number of cycles or from directional noise, it is not calibrated, and DSurfTomo does not use it (all travel times have equal weight).

### 7.6 Quality control

A point is rejected as soon as one criterion fails; the reasons are recorded as a bit mask (`reason` in `work/05_dispersion/disp.h5`) and counted in the report.

| Criterion | Meaning | Default | Source |
|---|---|---|---|
| `snr_g` | highest signal-to-noise ratio at that frequency | ≥ 4 | empirical |
| `snr_c` | signal-to-noise ratio at the pick | ≥ 6 | 5–15 common in the literature (Bensen et al. 2007 and others); the low end |
| `lam_c` | path at least this many wavelengths of its own phase velocity | 2 | Bensen et al. (2007) use 3; Luo et al. (2015) show that 1–2 wavelengths suffice; benchmark of this project on J0: median error 0.007 % at 2–3 λ |
| `dev_pre` | departure of the group pick from the reference | ≤ 30 % | empirical |
| `side_c` | relative height of a competing envelope peak | ≤ 1 | empirical |
| `jk_c_rel` | jackknife standard deviation / phase velocity | ≤ 2 % | empirical |
| `ucons_c` | difference between the measured group velocity and the one implied by the phase-velocity curve | ≤ 25 % | empirical |
| `band_edge` | amplitude spectrum of the coherent signal at that frequency relative to its mean within the filter band | ≥ 0.5 (−6 dB) | usual definition of a band edge |
| `min_run` | consecutive accepted frequencies | ≥ 6 | empirical |
| cycle skip | departure of the phase from the prediction during continuation | ≤ 0.3 cycles | AFTAN |

**The band-edge criterion**: beyond the band excited by the noise sources, a Gaussian filter passes only energy from neighbouring frequencies; what is measured is the velocity of another frequency, while the signal-to-noise ratio in the time domain is still high. The criterion compares amplitude spectra of the coherent signal: for each pair the energy spectrum of the signal window minus that of late lags (incoherent noise only), scaled by length; all pairs are summed and the square root is taken. Incoherent noise is the same in both windows and cancels, so the result does not depend on the stacking time or the noise level. On synthetic data of one and three days (noise sources cut off at 12 Hz) the ratio is 0.64–0.66 at 12.7 Hz, 0.33–0.38 at 13.45 Hz and 0.19–0.25 at 14.25 Hz; with half of the pairs the result is the same.

### 7.7 Export

For each pair only the longest run of consecutive accepted frequencies is kept (on a tie the one with the higher mean signal-to-noise ratio; as recommended in the manual of Yao's EGFAnalysisTimeFreq) and written as input for DSurfTomo. The periods are the inverses of the measured frequencies, without interpolation.

---

## 8. inversion: 1-D starting model and 3-D inversion

### 8.1 Average curve and depth nodes

At each period the median of the phase velocities of all pairs is taken, if there are at least `inversion.median_min_pairs` pairs (on small arrays a quarter of the pairs, at least 5). Depth nodes start at 0 with a step of max(shortest wavelength / 4, depth / 4): the shortest wavelength limits the thinnest layer that can be resolved, and the vertical resolution of surface waves degrades with depth. Nodes down to half of the longest wavelength are updated (the empirical depth of investigation in section 4.2.3 of Foti et al. 2018), and one fixed bottom node is added below (DSurfTomo does not update the last node).

### 8.2 1-D inversion

A damped, smoothed Gauss–Newton iteration in ln Vs. The forward problem uses DSurfTomo's own chain (Vp and density from Vs after Brocher 2005, dispersion from surfdisp96) through the small program `fwd1d`, so that the starting model and the 3-D inversion share the same forward problem. The smoothing constrains the curvature of the model itself ("jumping", Shaw & Orcutt 1985; Constable et al. 1987) and not of each update, so the result does not drift with the number of iterations; it is normalised by the actual node spacing. The initial value is 1.087 times the phase velocity at the shortest period.

### 8.3 3-D inversion

DSurfTomo (Fang et al. 2015) inverts travel times directly for 3-D Vs, tracing rays again in the current model at each iteration. lant writes `surfdata.dat` (station pairs, periods, phase velocities), `DSurfTomo.in` and the starting model `MOD`, and does not change the source code. Constraints of the interface:

- Number of periods: at most 60 for DSurfTomo and 80 for pDSurfTomo; beyond that the inversion step stops and suggests a suitable `dispersion.per_octave`.
- Reserve of the sparse matrix: DSurfTomo reserves sparsity × number of travel times × number of nodes non-zero elements; when the reserve is too small it writes beyond it before checking and ends with a segmentation fault. lant estimates the fill from the number of grid cells crossed by straight paths; `inversion.sparsity: auto` takes twice the estimate, at least 0.05; a value given by hand below 1.2 times the estimate is refused with an error.
- Rejection rule: at each iteration DSurfTomo removes data whose residual lies beyond `threshold` interquartile ranges. Before running, one forward pass checks the fraction kept under the starting model (`checks.precheck`); a project that would keep less than half, or whose median residual is off by more than one interquartile range, is marked invalid and not run.

---

## 9. evaluate: checks

Numbers and figures only, without pass criteria.

| Check | Procedure | How to read it |
|---|---|---|
| Travel-time residuals | root mean square, median and fraction kept by the rejection rule, for the starting and final models | how much the 3-D inversion reduced the residuals |
| Checkerboard (Lévêque et al. 1993) | a ±10 % checkerboard added to the 1-D model (on a line the pattern is odd about the middle of the grid, so that the synthetic residuals do not lean to one side; otherwise DSurfTomo's rejection rule on the signed quartiles drops most of the data), travel times computed for the actual pairs and periods, Gaussian noise added (`inversion.noise`; auto: the relative rms of the residuals of the final main model kept by DSurfTomo's rejection rule, so that the synthetic data are as noisy as the real data remain after the inversion — noise of the size of the data errors, Rawlinson & Spakman 2016; the misfit is an upper estimate of those errors, since it also contains structure the model does not fit; the checkerboards then run after the main project), then inverted; correlation and regression slope between recovered and input pattern at the nodes covered by stations | the slope is the fraction of an anomaly of that size that remains in the model. It depends on the regularisation and on the added noise as well as on the paths |
| Vertical checkerboard | the sign of the pattern flips with depth | leakage between depths |
| Starting-model control (Rawlinson & Spakman 2016) | the starting model gets a gradient of ±5 % with depth (centred on the travel-time residuals), everything else unchanged, then inverted | where the difference between the two models is about as large as the change applied, the data do not constrain the model; the residual history of the control project is given as well: as long as its residual has not come down to that of the main project, part of the difference is unconverged |
| Regularisation controls (optional) | smoothing and damping multiplied by 0.5 or 2 | dependence of the structure on the regularisation |

---

## 10. Conditions of validity and limitations

- The noise field is approximately isotropic. Directional noise biases two-station phase velocities; the workflow neither tests nor corrects this.
- Only the fundamental Rayleigh mode is measured; when a higher mode is stronger than the fundamental, quality control removes most but not all wrong points ([validation.md](validation.md), section 7).
- Vp and density follow from Vs through the empirical relations of Brocher (2005) (fixed inside DSurfTomo); a different Vp/Vs biases Vs.
- DSurfTomo stores coordinates in single precision: for pairs closer than 0.1 km the distance error can reach 1 %.
- Jackknife errors are not calibrated.
- The model can be interpreted only at the scales and depths that the checkerboard tests recover; nodes outside the station coverage are not constrained by data.

---

## References

- Bensen, G. D., Ritzwoller, M. H., Barmin, M. P., Levshin, A. L., Lin, F., Moschetti, M. P., Shapiro, N. M., & Yang, Y. (2007). Processing seismic ambient noise data to obtain reliable broad-band surface wave dispersion measurements. *Geophysical Journal International*, 169, 1239–1260. doi:10.1111/j.1365-246X.2007.03374.x
- Brocher, T. M. (2005). Empirical relations between elastic wavespeeds and density in the Earth's crust. *Bulletin of the Seismological Society of America*, 95, 2081–2092. doi:10.1785/0120050077
- Constable, S. C., Parker, R. L., & Constable, C. G. (1987). Occam's inversion. *Geophysics*, 52, 289–300. doi:10.1190/1.1442303
- Fang, H., Yao, H., Zhang, H., Huang, Y.-C., & van der Hilst, R. D. (2015). Direct inversion of surface wave dispersion for three-dimensional shallow crustal structure based on ray tracing. *Geophysical Journal International*, 201, 1251–1263. doi:10.1093/gji/ggv080
- Foti, S., et al. (2018). Guidelines for the good practice of surface wave analysis. *Bulletin of Earthquake Engineering*, 16, 2367–2420. doi:10.1007/s10518-017-0206-7
- Herrin, E., & Goforth, T. (1977). Phase-matched filters: application to the study of Rayleigh waves. *Bulletin of the Seismological Society of America*, 67, 1259–1275.
- Jiang, C., & Denolle, M. A. (2020). NoisePy: a new high-performance Python tool for ambient-noise seismology. *Seismological Research Letters*, 91, 1853–1866. doi:10.1785/0220190364
- Lévêque, J.-J., Rivera, L., & Wittlinger, G. (1993). On the use of the checker-board test to assess the resolution of tomographic inversions. *Geophysical Journal International*, 115, 313–318.
- Levshin, A. L., & Ritzwoller, M. H. (2001). Automated detection, extraction, and measurement of regional surface waves. *Pure and Applied Geophysics*, 158, 1531–1545. doi:10.1007/978-3-0348-8264-4_11
- Leys, C., Ley, C., Klein, O., Bernard, P., & Licata, L. (2013). Detecting outliers: do not use standard deviation around the mean, use absolute deviation around the median. *Journal of Experimental Social Psychology*, 49, 764–766. doi:10.1016/j.jesp.2013.03.013
- Landau, H. J. (1967). Necessary density conditions for sampling and interpolation of certain entire functions. *Acta Mathematica*, 117, 37–52. doi:10.1007/BF02395039
- Lin, F.-C., Moschetti, M. P., & Ritzwoller, M. H. (2008). Surface wave tomography of the western United States from ambient seismic noise: Rayleigh and Love wave phase velocity maps. *Geophysical Journal International*, 173, 281–298. doi:10.1111/j.1365-246X.2008.03720.x
- Luo, Y., Yang, Y., Xu, Y., Xu, H., Zhao, K., & Wang, K. (2015). On the limitations of interstation distances in ambient noise tomography. *Geophysical Journal International*, 201, 652–661. doi:10.1093/gji/ggv043
- Rawlinson, N., & Spakman, W. (2016). On the use of sensitivity tests in seismic tomography. *Geophysical Journal International*, 205, 1221–1243. doi:10.1093/gji/ggw084
- Seats, K. J., Lawrence, J. F., & Prieto, G. A. (2012). Improved ambient noise correlation functions using Welch's method. *Geophysical Journal International*, 188, 513–523. doi:10.1111/j.1365-246X.2011.05263.x
- Sens-Schönfelder, C. (2008). Synchronizing seismic networks with ambient noise. *Geophysical Journal International*, 174, 966–970. doi:10.1111/j.1365-246X.2008.03842.x
- Shapiro, N. M., & Singh, S. K. (1999). A systematic error in estimating surface-wave group-velocity dispersion curves and a procedure for its correction. *Bulletin of the Seismological Society of America*, 89, 1138–1142.
- Shaw, P. R., & Orcutt, J. A. (1985). Waveform inversion of seismic refraction data and applications to young Pacific crust. *Geophysical Journal of the Royal Astronomical Society*, 82, 375–414.
- Stehly, L., Campillo, M., & Shapiro, N. M. (2007). Traveltime measurements from noise correlation: stability and detection of instrumental time-shifts. *Geophysical Journal International*, 171, 223–230. doi:10.1111/j.1365-246X.2007.03492.x
- Takagi, R., Nishida, K., Aoki, Y., et al. (2015). A single bit matters: coherent noise of seismic data loggers. *Seismological Research Letters*, 86(3), 901–907. doi:10.1785/0220150030
- Vincenty, T. (1975). Direct and inverse solutions of geodesics on the ellipsoid with application of nested equations. *Survey Review*, 23(176), 88–93. doi:10.1179/sre.1975.23.176.88
- Wang, J., Wu, G., & Chen, X. (2019). Frequency-Bessel transform method for effective imaging of higher-mode Rayleigh dispersion curves from ambient seismic noise data. *Journal of Geophysical Research: Solid Earth*, 124, 3708–3723. doi:10.1029/2018JB016595
- Xu, Y., Lebedev, S., Meier, T., Bonadio, R., & Bean, C. J. (2021). Optimized workflows for high-frequency seismic interferometry using dense arrays. *Geophysical Journal International*, 227, 875–897. doi:10.1093/gji/ggab260
- Yao, H., van der Hilst, R. D., & de Hoop, M. V. (2006). Surface-wave array tomography in SE Tibet from ambient seismic noise and two-station analysis – I. Phase velocity maps. *Geophysical Journal International*, 166, 732–744. doi:10.1111/j.1365-246X.2006.03028.x
- Zhu, S., Li, J., Chen, G., Fang, H., & Yao, H. (2026). pDSurfTomo: a high-performance parallel computing package for direct surface wave tomography. arXiv:2604.11920.
