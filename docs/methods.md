# Methods: principles, formulas and sources

[中文](methods.zh.md) | English

This page says what every step computes, why, and on what it rests. Each number in the configuration is listed with its basis in [parameters.md](parameters.md); how to run the steps is in [tutorial.md](tutorial.md).

A value that is a choice rather than a result of theory is marked "convention".

## 1. Why correlations of noise give surface waves

If the noise field is diffuse (sources distributed evenly around the stations, or enough scattering), the time derivative of the cross-correlation of the records at two stations A and B is proportional to the Green's function between them and its time reverse (Snieder 2004; Wapenaar 2004; Wapenaar et al. 2010):

  dC_AB(τ)/dτ ∝ −G(x_B, x_A, τ) + G(x_B, x_A, −τ).

Sources in line with the two stations (the stationary-phase zones) add up coherently; other sources cancel out in the stack. Sources at the surface excite mainly surface waves, which spread over a surface (amplitude ∝ 1/√r) and dominate the correlations of vertical records: these are fundamental-mode Rayleigh waves, which lant measures.

In a laterally uniform medium with phase velocity c(ω) the spectrum of the correlation of two stations r apart is (Aki 1957; Tsai & Moschetti 2010)

  C(r, ω) ∝ J₀(ω r / c(ω)).

lant uses this relation twice: for the array F-J spectrum (section 4) and for the phase of two-station measurements, whose far-field form is J₀(x) ≈ √(2/πx) cos(x − π/4) (section 5).

**Limits.** Uneven noise sources bias the travel times (Tsai 2009); lant does not correct this, and the result holds where the noise is close to isotropic. Higher modes can be stronger than the fundamental mode where soft layers lie on hard rock; the reference curve has to be checked for that (section 4).

## 2. Cleaning the records (`lant clean`)

- **Index.** Every file under `data.raw_dir` is read once with ObsPy (headers only); `1_clean/record_index.csv` lists every record.
- **Stations.** A record belongs to a station of the table when its network and station codes match the table (or its station code alone, for a table without networks). Names are never guessed; what does not match is listed in the log.
- **Channel.** `data.channel` (default `*Z`): the vertical component.
- **Day.** Each UTC day of each station is put on its exact sample grid (Lanczos resampling only if a record is off the grid by more than 5 % of a sample). Gaps stay gaps (no interpolation, no zero filling). The median is removed.
- **Instrument response** (only when `data.response` is given): removed to ground velocity with ObsPy, with a pre-filter that is flat inside the correlation band `cc.fmin_hz`–`cc.fmax_hz` (Bensen et al. 2007, sect. 2.1). Without response files the records are used as they are. This is exact only when every station has the same instrument: the cross-coherence (below) removes the amplitude response, and the phase responses of identical instruments cancel in every pair.
- **Decimation** to `clean.fs_hz`: zero-phase 8th-order Butterworth low-pass at 0.8 × the new Nyquist frequency, then every q-th sample (q a whole number). A decimated sample counts as valid only if the raw data within ±max(1 s, 50 samples) were all present.
- **Output**: one HDF5 file per day, `1_clean/days/YYYY-MM-DD.h5`, an array station × sample with NaN where no sample is valid, and the hours of valid data per station (`availability.csv`). Nothing is removed for its amplitude in this step: whether a stretch of record is used depends on the segment length, so it is decided in `lant cc` (section 3).

## 3. Cross-correlation and stack (`lant cc`)

1. **Segments** of `cc.segment_s`, overlapping by `cc.overlap` (Welch averaging; Seats et al. 2012). A segment is not used if it contains a gap, or if its largest |amplitude| exceeds `cc.max_over_std` × the standard deviation of the station-day (earthquakes and glitches; `max_over_std` of NoisePy, Jiang & Denolle 2020). This replaces the temporal normalisation (one-bit or running-mean) of Bensen et al. (2007): a disturbed segment is left out instead of being flattened. Each segment: linear trend removed, 5 % cosine taper.
2. **Cross-coherence** (Xu et al. 2021 compared normalisations for high-frequency dense arrays and recommend it):

   C_ij(f) = F_i*(f) F_j(f) / (|F_i(f)| |F_j(f)|),

   that is, both spectra are whitened (Bensen et al. 2007) and only the phase difference is kept. The band `cc.fmin_hz`–`cc.fmax_hz` is kept, with quarter-octave cosine edges.
3. **Sub-stacks**: the segments of `cc.substack_hours` hours are averaged in the frequency domain and transformed to lags −`maxlag`…+`maxlag`. A positive lag means energy travelling from station 1 to station 2 of the pair (names in sorted order). The segment must be longer than twice the lag window, because the correlation computed with the FFT is circular.
4. **Stack**: the average over all segments (linear stack; Bensen et al. 2007).

Tf-PWS (Schimmel et al. 2011) is not used: phase weighting changes the waveform in a way that depends on the data, and the measured quantity here is the phase itself.

## 4. Reference curve (`lant disp`, part 1)

The reference curve does three things: it sets where the group arrival is looked for, it fixes the 2π ambiguity (cycle count) of the phase, and the frequencies where it is defined are the frequencies measured. Its values do not enter the measurements. Three sources, as in the literature (`disp.reference`):

| Source | What it is | Used by |
|---|---|---|
| `fj` | the F-J spectrum of this array | dense-array studies (Wang et al. 2019; Li et al. 2021) |
| `model` | phase velocities of a 1-D Vs model (`depth_km Vs_km/s`) | Bensen et al. (2007), Lin et al. (2008): predictions of a model |
| `file` | a curve `frequency_Hz c_km/s`, e.g. picked by hand from the F-J figure | the automatic mode of EGFAnalysisTimeFreq (Yao et al. 2006) |

**F-J spectrum** (Wang, Wu & Chen 2019):

  I(ω, c) = ∫ C(r, ω) J₀(ωr/c) r dr,

largest where c is the phase velocity of a mode. As in the open-source CC-FJpy (Li et al. 2021) for noise data, C is the real part of the spectrum of the stack, and the integral is the trapezoidal rule over the sorted pair distances (double precision; every pair with data). The spectrum is free of spatial aliasing for k ≤ π/Δr, where Δr is the sampling interval of the pair distances (Wang et al. 2019, eq. 14).

**Picking** (`fj`): at each frequency of `disp.fj_band_hz` (default: all), the largest value inside `data.velocity_range_kms`; a cubic polynomial of ln c in ln f through these picks; defined over the band of the picks only (no extrapolation). The F-J spectrum is always computed and drawn, together with whichever reference is used, so that it can be checked: **the reference must follow the fundamental mode** (the slowest ridge). If a higher mode is picked, every cycle count and so every phase velocity is wrong. lant does not decide this by itself: look at `3_disp/figures/fj_reference.png`.

Reference group velocity: U = c / (1 − d ln c / d ln f) (Bensen et al. 2007, eq. 7).

## 5. Two-station phase and group velocities (`lant disp`, part 2)

For each pair and each measurement frequency f (`disp.per_octave` per octave):

1. **Symmetric component**: average of the positive lags and the time-reversed negative lags (Bensen et al. 2007).
2. **Narrow-band filters**: Gaussian, G(f′) = exp[−α ((f′ − f)/f)²] (constant relative width, FTAN; Levshin & Ritzwoller 2001), α = `disp.alpha`.
3. **Group arrival** t_g: the maximum of the envelope between r/((1+w)U_ref) and r/((1−w)U_ref), w = `disp.group_window`. A maximum on the edge of this window is not an arrival. **SNR**: the envelope at t_g divided by the rms of the filtered signal from two periods after the window to the end of the lag window (at least its last 15 %).
   **Group velocity** U = r / t_g. The arrival time belongs to the instantaneous frequency of the filtered signal at t_g (the phase derivative there), not to the centre of the filter (Levshin et al. 1972; Bensen et al. 2007); within each run of consecutive arrivals the group times are interpolated (in ln f) from the instantaneous frequencies back onto the measurement frequencies, never across a gap.
4. **Phase-matched filter** (Herrin & Goforth 1977; FTAN): the spectrum is multiplied by exp(iΨ), Ψ(ω) = ∫ t_g dω, which removes the dispersion of the arrival; the Gaussian-weighted sum around f gives its phase Φ(f). Reading the phase of the narrow-band signal at t_g instead is biased when the group delay changes inside the filter band.
5. **Phase travel time**:

   t_p = −(Φ − π/4) / (2πf) + N / f,  c = r / t_p,

   where π/4 is the far-field phase of J₀ (Lin et al. 2008) and N is the cycle count. r is the distance on the WGS84 ellipsoid (Vincenty 1975).
6. **Cycle count.** At the lowest frequency of each run of usable frequencies, N is the value that brings c closest to the reference. Neighbouring values of N differ by λ/r in relative velocity, so the choice is unambiguous when the reference error e is smaller than half of that, i.e. **when the path is at most 1/(2e) wavelengths long** (`disp.reference_error`, default 0.1, i.e. 5 wavelengths). Runs that start on longer paths are not measured (reason `cycle_ambiguous`). N is then followed to higher frequencies with the measured group delay (AFTAN):

   ω_k t_k = ω_(k−1) t_(k−1) + ½ (t_g,k + t_g,k−1)(ω_k − ω_(k−1)),

   on a grid four times denser than the measurement grid; a jump of more than `disp.slip_cycles` of a cycle ends the run.
7. **Error**: delete-one-group jackknife over the sub-stacks of `lant cc` (each sub-stack left out in turn; Efron & Stein 1981): σ² = (g−1)/g Σ (c_i − c̄)², and the same for U. It measures the repeatability over time only: not the errors of the reference, of the cycle count or of uneven noise.
8. **Quality control** (all thresholds are in `disp`), for c and U separately (the cycle count and the phase tracking concern c only):

| Criterion | Default | Basis |
|---|---|---|
| SNR ≥ `snr_min` | 5 | common choice 5–15 (Bensen et al. 2007 use 10 for their long-period data) |
| path ≥ `min_wavelengths` × λ (λ = c/f, measured c) | 2 | far-field approximation: Bensen et al. (2007) use 3; Luo et al. (2015) show that 1–2 wavelengths suffice for phase velocities; the unit test on exact J₀ correlations: within 1.5 % from 3 wavelengths, 0.4 % beyond 5 |
| σ ≤ `max_jackknife_rel` × c (or × U) | 2 % | convention; `null` switches it off |
| one run of at least `min_points` frequencies | 6 | one segment per pair (EGFAnalysisTimeFreq manual, Yao et al.); length a convention |

Every rejected value keeps its reasons (bit masks `reason` for c and `reason_u` for U in `measurements.h5`; names in `lant.dispersion.REASONS`). `median_vs_reference.csv` compares the median of the accepted values with the reference at each frequency. A wrong cycle count changes c by a branch spacing, at least 2 × `disp.reference_error` (20 % by default) on the paths where N is set; a median far from the reference by that order, or jumping between frequencies, points to wrong cycle counts and a wrong reference. An offset of a few percent that does not change with frequency is a difference between the area and the reference, which the measurement allows.

**Accuracy on exact correlations** (J₀ spectra of a dispersive medium, `tests/test_measurement.py`): c within 0.4 % (95 %) beyond 5 wavelengths and 1.5 % at 3–5 wavelengths; U within 1 % (95 %), with a median bias of about +0.5 % at α = 5 (the envelope maximum of a wide filter in a dispersive wave; a larger α reduces it).

## 6. 3-D inversion (`lant inv`)

**Data.** The accepted phase velocities, group velocities or both (`inv.data`: `phase`, `group`, `both`) at the periods on the 1/`inv.per_octave`-octave grid through 1 Hz (neighbouring measurement frequencies overlap by more than 90 % in the Gaussian filters, so fewer periods carry the same information and the computing time scales with their number).

**Starting model.** The median velocity of each period over the pairs, for each kind used (a period needs max(5, min(50, pairs/4)) pairs). Wavelengths (velocity × period) come from the phase velocities, or from the group velocities when only these are used. Depth nodes from 0, in steps of max(shortest wavelength/4, depth/4), to the first node past half the longest wavelength (the usual estimate of the investigation depth, Foti et al. 2018, sect. 4.2.3); the last node is a fixed half-space. The 1-D model is a damped, smoothed Gauss–Newton inversion of ln Vs that fits the average curves used; the smoothing acts on the curvature of the model itself (Shaw & Orcutt 1985; Constable et al. 1987), so the result does not depend on the number of iterations. Or a model file (`inv.start_model_file`).

**Forward calculation.** Vp and density from Vs by Brocher (2005); phase and group velocities and their depth kernels from DSurfTomo's surfdisp96 (Herrmann 2013; Fang et al. 2015), with the same numbers as the original Fortran (`lant/inversion/dsurf_disp96.py`). The root search steps 0.005 km/s; where Vs decreases with depth somewhere in a column, two modes can come closer than that and the search jumps to a higher mode or finds no root, so such columns are solved with a step of 0.0005 km/s. The same calculation gives the reference curve from a model.

**Direct 3-D inversion** (Fang et al. 2015, the method of DSurfTomo; `lant/inversion/dsurf_*.py`): the unknowns are Vs on a 3-D grid (horizontal spacing `inv.grid_km`, default the median measured wavelength/3: lateral resolution of surface waves is about a wavelength). Each iteration traces the rays of every period in the current phase- or group-velocity maps with the fast-marching method, computes the depth kernels, and solves for the update of Vs:

  min (1/N) ‖W(G Δm − r)‖² + (w²/n) ‖L(m + Δm − m₀)‖² + (d²/n) ‖m + Δm − m₀‖²

r: travel-time residuals; G: travel-time derivatives; W: 0/1 weights of the rejection; N kept travel times; n parameters; L: 7-point Laplacian of the whole grid (every node with its neighbours); m₀: starting model; w = `inv.weight`, d = `inv.damping`. Differences from DSurfTomo:

- smoothing and damping act on the model, not on each update (Occam's "jumping" form; Constable et al. 1987), so that the result does not depend on the number of iterations;
- data and model terms are means, so the same weight means the same for any number of data and nodes;
- the least-squares problem is solved to convergence in double precision (LSMR, Fong & Saunders 2011, with column scaling, Paige & Saunders 1982);
- for group travel times the ray and the derivative −ds/U² are taken in the group-velocity map, like the time itself, so that G is the derivative of the computed times; DSurfTomo takes them in the phase-velocity map.

Because of the first two points, the values of `inv.weight` and `inv.damping` are not those of a DSurfTomo.in.

**Rejection**: a travel time further than `inv.outlier_mad` × 1.4826 MAD from the median residual of its period gets weight 0 (Rousseeuw & Croux 1993; k = 3 is "very conservative", Leys et al. 2013). Per period (and for phase and group apart), because a phase error becomes a travel-time error proportional to the period.

The project in `4_inv/dsurftomo/` is in DSurfTomo's file format and runs with the original program (with DSurfTomo's own regularization).

## 7. Resolution test (`lant checker`)

Checkerboard test (Lévêque et al. 1993): the starting model times 1 + a sin(πx/s) sin(πy/s); with `--flip-depth` the sign changes below that depth. Travel times of this model are computed for the same pairs and periods, Gaussian noise is added (`checker.noise` auto: per period, the robust standard deviation of the residuals the real inversion leaves, so the test is about as noisy as the data; Rawlinson & Spakman 2016), and the data are inverted with the same settings. `recovery.csv` gives, per depth, the correlation and the regression slope of recovered on true perturbation at the nodes with data; the slope is the fraction of an anomaly of that size the inversion recovers. One size says nothing about others: run several sizes.

## 8. Limits

- Isotropic noise is assumed; directional noise biases the phase velocities (Tsai 2009) and is not corrected.
- Only the fundamental Rayleigh mode is measured; where a higher mode is stronger, check the reference curve.
- Topography is ignored: distances are on the ellipsoid.
- Vp and density follow Vs by Brocher (2005); a different Vp/Vs biases Vs.
- The jackknife error is not a full error.
- Group velocities carry the bias of the filter width (about +0.5 % at α = 5 on exact data); with `inv.data: both` the phase and group times are then not quite consistent.
- The model can be interpreted only at the sizes and depths the checkerboards recover.

## References

- Aki, K. (1957). Space and time spectra of stationary stochastic waves, with special reference to microtremors. *Bull. Earthq. Res. Inst.*, 35, 415–456.
- Aster, R. C., Borchers, B., & Thurber, C. H. (2013). *Parameter Estimation and Inverse Problems* (2nd ed.). Academic Press.
- Bensen, G. D., et al. (2007). Processing seismic ambient noise data to obtain reliable broad-band surface wave dispersion measurements. *GJI*, 169, 1239–1260. doi:10.1111/j.1365-246X.2007.03374.x
- Brocher, T. M. (2005). Empirical relations between elastic wavespeeds and density in the Earth's crust. *BSSA*, 95, 2081–2092. doi:10.1785/0120050077
- Constable, S. C., Parker, R. L., & Constable, C. G. (1987). Occam's inversion. *Geophysics*, 52, 289–300. doi:10.1190/1.1442303
- Efron, B., & Stein, C. (1981). The jackknife estimate of variance. *Ann. Stat.*, 9, 586–596. doi:10.1214/aos/1176345462
- Fang, H., Yao, H., Zhang, H., Huang, Y.-C., & van der Hilst, R. D. (2015). Direct inversion of surface wave dispersion for three-dimensional shallow crustal structure based on ray tracing. *GJI*, 201, 1251–1263. doi:10.1093/gji/ggv080
- Fong, D. C.-L., & Saunders, M. (2011). LSMR: an iterative algorithm for sparse least-squares problems. *SIAM J. Sci. Comput.*, 33, 2950–2971. doi:10.1137/10079687X
- Foti, S., et al. (2018). Guidelines for the good practice of surface wave analysis. *Bull. Earthq. Eng.*, 16, 2367–2420. doi:10.1007/s10518-017-0206-7
- Herrin, E., & Goforth, T. (1977). Phase-matched filters: application to the study of Rayleigh waves. *BSSA*, 67, 1259–1275.
- Herrmann, R. B. (2013). Computer programs in seismology: an evolving tool for instruction and research. *SRL*, 84, 1081–1088. doi:10.1785/0220110096
- Jiang, C., & Denolle, M. A. (2020). NoisePy: a new high-performance Python tool for ambient-noise seismology. *SRL*, 91, 1853–1866. doi:10.1785/0220190364
- Lévêque, J.-J., Rivera, L., & Wittlinger, G. (1993). On the use of the checker-board test to assess the resolution of tomographic inversions. *GJI*, 115, 313–318.
- Levshin, A. L., Pisarenko, V. F., & Pogrebinsky, G. A. (1972). On a frequency-time analysis of oscillations. *Annales de Géophysique*, 28, 211–218.
- Levshin, A. L., & Ritzwoller, M. H. (2001). Automated detection, extraction, and measurement of regional surface waves. *PAGEOPH*, 158, 1531–1545.
- Leys, C., et al. (2013). Detecting outliers: do not use standard deviation around the mean, use absolute deviation around the median. *J. Exp. Soc. Psychol.*, 49, 764–766. doi:10.1016/j.jesp.2013.03.013
- Li, Z., Zhou, J., Wu, G., Wang, J., et al. (2021). CC-FJpy: a Python package for extracting overtone surface-wave dispersion from seismic ambient-noise cross correlation. *SRL*, 92, 3179–3186. doi:10.1785/0220210042
- Lin, F.-C., Moschetti, M. P., & Ritzwoller, M. H. (2008). Surface wave tomography of the western United States from ambient seismic noise: Rayleigh and Love wave phase velocity maps. *GJI*, 173, 281–298. doi:10.1111/j.1365-246X.2008.03720.x
- Luo, Y., Yang, Y., Xu, Y., Xu, H., Zhao, K., & Wang, K. (2015). On the limitations of interstation distances in ambient noise tomography. *GJI*, 201, 652–661. doi:10.1093/gji/ggv043
- Paige, C. C., & Saunders, M. A. (1982). LSQR: an algorithm for sparse linear equations and sparse least squares. *ACM TOMS*, 8, 43–71. doi:10.1145/355984.355989
- Rawlinson, N., & Spakman, W. (2016). On the use of sensitivity tests in seismic tomography. *GJI*, 205, 1221–1243. doi:10.1093/gji/ggw084
- Rost, S., & Thomas, C. (2002). Array seismology: methods and applications. *Reviews of Geophysics*, 40(3), 1008. doi:10.1029/2000RG000100
- Rousseeuw, P. J., & Croux, C. (1993). Alternatives to the median absolute deviation. *JASA*, 88, 1273–1283. doi:10.1080/01621459.1993.10476408
- Schimmel, M., Stutzmann, E., & Gallart, J. (2011). Using instantaneous phase coherence for signal extraction from ambient noise data at a local to a global scale. *GJI*, 184, 494–506. doi:10.1111/j.1365-246X.2010.04861.x
- Seats, K. J., Lawrence, J. F., & Prieto, G. A. (2012). Improved ambient noise correlation functions using Welch's method. *GJI*, 188, 513–523. doi:10.1111/j.1365-246X.2011.05263.x
- Shaw, P. R., & Orcutt, J. A. (1985). Waveform inversion of seismic refraction data and applications to young Pacific crust. *Geophys. J. R. astr. Soc.*, 82, 375–414.
- Snieder, R. (2004). Extracting the Green's function from the correlation of coda waves: a derivation based on stationary phase. *Phys. Rev. E*, 69, 046610. doi:10.1103/PhysRevE.69.046610
- Tsai, V. C. (2009). On establishing the accuracy of noise tomography travel-time measurements in a realistic medium. *GJI*, 178, 1555–1564. doi:10.1111/j.1365-246X.2009.04239.x
- Tsai, V. C., & Moschetti, M. P. (2010). An explicit relationship between time-domain noise correlation and spatial autocorrelation (SPAC) results. *GJI*, 182, 454–460. doi:10.1111/j.1365-246X.2010.04633.x
- Vincenty, T. (1975). Direct and inverse solutions of geodesics on the ellipsoid with application of nested equations. *Survey Review*, 23, 88–93. doi:10.1179/sre.1975.23.176.88
- Wang, J., Wu, G., & Chen, X. (2019). Frequency-Bessel transform method for effective imaging of higher-mode Rayleigh dispersion curves from ambient seismic noise data. *JGR Solid Earth*, 124, 3708–3723. doi:10.1029/2018JB016595
- Wapenaar, K. (2004). Retrieving the elastodynamic Green's function of an arbitrary inhomogeneous medium by cross correlation. *Phys. Rev. Lett.*, 93, 254301. doi:10.1103/PhysRevLett.93.254301
- Wapenaar, K., Draganov, D., Snieder, R., Campman, X., & Verdel, A. (2010). Tutorial on seismic interferometry: Part 1 — Basic principles and applications. *Geophysics*, 75, 75A195–75A209. doi:10.1190/1.3457445
- Xu, Y., Lebedev, S., Meier, T., Bonadio, R., & Bean, C. J. (2021). Optimized workflows for high-frequency seismic interferometry using dense arrays. *GJI*, 227, 875–897. doi:10.1093/gji/ggab260
- Yao, H., van der Hilst, R. D., & de Hoop, M. V. (2006). Surface-wave array tomography in SE Tibet from ambient seismic noise and two-station analysis — I. Phase velocity maps. *GJI*, 166, 732–744. doi:10.1111/j.1365-246X.2006.03028.x
