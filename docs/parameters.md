# Parameters

[中文](parameters.zh.md) | English

Every key of the configuration file, with its default, what it means and what the default rests on. This file is written by `tools/make_parameter_doc.py` from `lant/config.py` (defaults) and `lant/params.py` (texts); `lant init` writes the same information into the configuration file, and the step reports (`OUT/.lant/steps/<step>.md`) show the value each step used and, for derived values, the derivation with your numbers.

## 1. How to read the tables

Only three keys are required: `project.output_dir`, `data.raw_dir`, `data.station_file`. A value written in the configuration always wins. `auto` means: derived from your data by the rule in the description.

The column "basis" says what a default rests on:

| basis | meaning |
|---|---|
| required | only you know it |
| literature | a value or method taken from the cited source |
| theory | follows from a relation of the method (sampling, filter width, far-field phase) |
| rule | derived from your data by a written rule; the rule's own factors are project choices unless a source is cited |
| empirical | experience of this project; no literature value was found. Change it when your data call for it, and say so |
| choice | yours; there is no right value |
| technical | affects speed, memory or files only, not the result |

Where the chain stops for you: after `clean` (ambiguous names, stations without data, coordinate mismatches, unknown component orientation, an array type between areal and line, responses not removed and instrument models unknown), after `stack` (only when the diagnosis lists candidates), after `reference` (is the picked ridge the fundamental mode?), after `inversion` (starting model, deepest node, velocity bounds). A step with nothing to judge does not stop.

## 2. project

| key | default | basis | meaning |
|---|---|---|---|
| `name` | `my_array` | choice | project name (figures and records only) |
| `array_type` | `auto` | rule | auto: from the station geometry (how many times the station cloud is longer than wide; below 3 areal, above 10 line, in between shown for confirmation); or surface / line |
| `output_dir` | — | required | every result, figure and record is written here |

## 3. data: raw data

| key | default | basis | meaning |
|---|---|---|---|
| `raw_dir` | — | required | raw-data directory, searched recursively; any layout and file names |
| `format` | `auto` | technical | auto: ObsPy detects the format of each file (MSEED, SAC, GSE2, SEG2, ...); or one ObsPy format name |
| `file_patterns` | `['*']` | technical | file-name patterns that are tried |
| `ignore_patterns` | `['*.txt', '*.csv', '*.xml', '*.json', '*.yaml', '*.yml', ...` | technical | file-name patterns that are not read |
| `station_file` | — | required | station table in any common form (CSV, text, Excel, StationXML) |
| `component` | `Z` | choice | Z: vertical-component correlations (Rayleigh waves); the only one supported |
| `channel` | `auto` | rule | auto: the vertical channel of each station from the channel name or the SAC header; or a pattern such as "*3" when the names do not tell |
| `network_default` | `null` | choice | network code for record headers without one |
| `rename` | `{}` | choice | {header name: table name} where the automatic matching fails |
| `day_start_utc_hour` | `auto` | rule | processing days start at this UTC hour; auto: the hour at which >= 60 % of the records longer than 12 h start, else 0 |
| `day_label` | `utc_start` | choice | how a processing day is labelled: utc_start, or a fixed offset such as "UTC+8" |
| `velocity_range_kms` | `[0.1, 5.0]` | choice | approximate Rayleigh-wave phase velocities of the area [slow, fast] in km/s; sets the lag window, the band and the F-J scan.  The default covers soft soil to hard rock; a narrower range saves computing time (the reference step reports the measured velocities) |
| `response` | `null` | choice | instrument response file or directory (StationXML, dataless SEED, RESP): the responses are removed in the clean step (Bensen et al. 2007).  Without it every station must carry the same instrument model |
| `line` | `null` | choice | when the station table has a line column: the line this project processes (one project per line) |
| `min_trace_length_s` | `0.0` | empirical | pieces of record shorter than this are not used |
| `exclude_stations` | `[]` | choice | stations left out entirely (an instrument or coupling problem you know of) |
| `start_date` | `null` | choice | first processing day used (YYYY-MM-DD) |
| `end_date` | `null` | choice | last processing day used (YYYY-MM-DD) |

## 4. clean: cleaning (`01_clean`)

| key | default | basis | meaning |
|---|---|---|---|
| `write_archive` | `true` | technical | write the cleaned day files (about the size of the raw data); false: later steps read the raw files, and the 1-s periodic signal is not tested (response removal and moved nodes need the archive) |
| `nproc` | `4` | technical | processes |
| `encoding` | `FLOAT32` | technical | miniSEED encoding of the archive |
| `min_free_disk_factor` | `1.2` | technical | refuse to write the archive when the free disk is below this factor x its estimated size |

## 5. correlation: hourly cross-coherence (`02_ccf`)

| key | default | basis | meaning |
|---|---|---|---|
| `fs_max` | `50.0` | choice | highest processing rate (Hz); the top of the band is 0.4 fs / (1 + sqrt(ln 2 / alpha)): raise it for small, high-frequency arrays |
| `fs` | `auto` | rule | rate after decimation; auto: the highest whole-number rate <= fs_max that every raw rate reaches by an integer factor |
| `seg_len_s` | `auto` | rule | segment length; auto: max(2 x maxlag (circular correlation), 10 longest periods) |
| `seg_step_s` | `auto` | literature | segment step; auto: half the segment (50 % overlap, Seats et al. 2012) |
| `maxlag_s` | `auto` | rule | lag window; auto: longest distance / (0.6 x slowest velocity) + 4 longest periods (the slowest arrival, then a stretch of noise for the SNR) |
| `fmin_hz` | `auto` | rule | bottom of the correlation band; auto: where the longest pair is 2 wavelengths at the slowest velocity, widened to the half-gain point of the Gaussian filter there |
| `fmax_hz` | `auto` | theory | top of the correlation band; auto: dispersion.fmax_hz widened to the half-gain point of its Gaussian filter, not above 0.4 fs (anti-alias filter) |
| `max_distance_km` | `null` | choice | only pairs up to this distance are correlated (none: all pairs); for very large arrays, where all pairs do not fit into memory |
| `smooth_bins` | `1` | literature | smoothing of the coherence denominator in frequency bins; 1 = none (cross-coherence as in Xu et al. 2021) |
| `max_over_std` | `10.0` | literature | a segment is dropped when its peak exceeds this many standard deviations of the station-day (NoisePy max_over_std, Jiang & Denolle 2020) |
| `remove_1s_comb` | `auto` | literature | subtract the waveform folded modulo 1 s, hour by hour (coherent logger noise at 1-s intervals, Takagi et al. 2015); auto: for the stations where the clean step detects it (F test of the folded records); true: all stations; false: none |
| `nproc` | `4` | technical | days processed at the same time (memory: one day of every station per process) |
| `threads_per_process` | `2` | technical | FFT threads per process |

## 6. stack: stack and station diagnosis (`03_stack`)

| key | default | basis | meaning |
|---|---|---|---|
| `block_hours` | `auto` | rule | sub-stack length of the jackknife; auto: the longest of 24, 12, 8, 6, 4, 3, 2, 1 h that gives the median station at least 8 sub-stacks |
| `station_qc.enabled` | `true` | choice | station diagnosis: sub-stacks of a station that differ from the rest in correlation or time shift are listed (never removed without you) |
| `station_qc.partners` | `20` | empirical | each station is compared through its pairs with this many nearest stations |
| `station_qc.min_partners` | `5` | empirical | pairs with data needed to judge a station sub-stack (fewer: not judged, and reported so) |
| `exclude_station_days` | `[]` | choice | station-days you remove: [station, YYYYDDD], [station, "YYYYDDD-"], ... |
| `exclude_station_hours` | `[]` | choice | station-hours you remove: [station, YYYYDDD, hour] |
| `exclude_candidates` | `false` | choice | also remove the candidates of the diagnosis (set by `lant confirm stack --exclude-candidates`) |

## 7. reference: reference curve (`04_reference`)

| key | default | basis | meaning |
|---|---|---|---|
| `source` | `fj` | literature | fj: frequency-Bessel spectrum of this array (Wang et al. 2019); file: a curve you give (as AFTAN and Yao's EGFAnalysisTimeFreq take a reference curve) |
| `file` | `null` | choice | the borrowed curve: npz written by lant, or text "freq_Hz c_km/s" |
| `zone_km` | `null` | choice | none: one reference curve for the whole array; a size in km: one F-J curve per square zone of that size from the pairs whose path midpoint lies in it, for arrays across which the structure changes strongly; a zone without a clear curve uses that of the whole array |
| `fj_mode` | `max_peak` | empirical | max_peak: the strongest ridge (areal arrays); track_fundamental: the ridge followed from an unambiguous start (lines) |
| `normal_dispersion` | `true` | choice | track_fundamental: the curve may not rise with frequency; false where a velocity reversal makes the ridge really rise |
| `fj_rmax_lambda` | `5.0` | empirical | pairs up to this many wavelengths enter the F-J spectrum of a frequency |
| `fj_c_nominal` | `auto` | rule | velocity that turns fj_rmax_lambda into a distance; auto: the curve picked in a first pass with every pair (two passes) |
| `fj_cmin` | `auto` | rule | slowest scanned velocity; auto: 0.5 x the slowest of data.velocity_range_kms |
| `fj_cmax` | `auto` | rule | fastest scanned velocity; auto: 1.5 x the fastest of data.velocity_range_kms |
| `fj_dc` | `auto` | rule | velocity step of the scan; auto: the 1-2-2.5-5 step not above fj_cmin / 200 |
| `fj_nfft` | `auto` | technical | FFT length of the spectra; auto: power of 2 >= 4 x the lag window, >= 8192 |
| `fj_fit_fmax` | `auto` | rule | highest frequency of the fitted picks; auto: the top of the longest run of clear picks |
| `track_start_band_hz` | `auto` | rule | where the ridge tracking starts; auto: of starts every 1/3 octave, the one followed over the most frequencies (on a tie the slower) |

## 8. dispersion: measurement and quality control (`05_dispersion`)

| key | default | basis | meaning |
|---|---|---|---|
| `fmin_hz` | `auto` | rule | lowest measured frequency; auto: where the longest pair is qc.lam_c wavelengths of the reference curve |
| `fmax_hz` | `auto` | theory | highest measured frequency; auto: 0.4 fs / (1 + sqrt(ln 2 / alpha)), the highest whose Gaussian filter stays inside the anti-alias band; the QC decides how far the data reach |
| `per_octave` | `12` | literature | measurement frequencies per octave, log-spaced (Xu et al. 2021) |
| `vmin` | `auto` | rule | bottom of the velocity axis; auto: 0.6 x the slowest reference velocity in the band |
| `vmax` | `auto` | rule | top of the velocity axis; auto: 1.3 x the fastest |
| `dv` | `auto` | rule | step of the velocity axis; auto: the 1-2-2.5-5 step not above vmin / 100 |
| `alpha` | `5.0` | literature | relative width of the Gaussian filter exp(-alpha ((f - fc) / fc)^2), as in FTAN (Bensen et al. 2007); larger = narrower |
| `phase_constant` | `pi/4` | theory | phase of the far-field 2-D Green's function, pi/4 (Lin et al. 2008 eq. 7) |
| `min_dist_km` | `0.0` | choice | pairs shorter than this are not measured (0: lam_min alone rules) |
| `lam_min` | `1.0` | literature | paths shorter than this many reference wavelengths are not measured (Luo et al. 2015: one wavelength can be enough) |
| `anchor_max_lambda` | `3.5` | empirical | the cycle count is fixed where the path is lam_min to this many wavelengths long (few cycles: the reference cannot be a cycle off) |
| `njack` | `5` | empirical | delete-a-group jackknife groups of the sub-stacks; auto: the sub-stacks holding >= 25 % of the fullest one |
| `min_jackknife_groups` | `3` | empirical | fewer groups with data: no error estimate, the sample is rejected |
| `nproc` | `8` | technical | processes |
| `chunk` | `20` | technical | pairs per task |
| `side_off_bands` | `''` | choice | bands 'lo-hi,...' where the side-lobe rule is not applied |
| `qc.snr_c` | `6.0` | literature | SNR at the picked arrival (peak / rms of the trailing noise); the literature uses 4-10 (Bensen et al. 2007: 10 for global data) |
| `qc.jk_c_rel` | `0.02` | empirical | jackknife standard error allowed, relative to the velocity |
| `qc.ucons_c` | `0.25` | empirical | disagreement allowed between the measured group velocity and the one implied by the phase curve (0 = off) |
| `qc.side_c` | `1.0` | empirical | height of a competing envelope peak allowed, relative to the pick |
| `qc.dev_pre` | `0.3` | empirical | deviation of the group pick allowed from the reference (0 = off) |
| `qc.lam_c` | `2.0` | literature | accepted samples need a path of this many of their own wavelengths (Bensen et al. 2007 use 3; Luo et al. 2015 show 1 can be enough) |
| `qc.min_run` | `6` | empirical | accepted samples must form runs of at least this many frequencies |
| `qc.band_edge` | `0.5` | choice | a frequency is not measured where the amplitude spectrum of the coherent signal (all pairs together; signal window minus the later lags, so that the noise level does not matter) is below this fraction of its mean over the filter band (0.5 = -6 dB, the usual definition of a band edge): beyond the band the noise excites, the filter passes only energy of other frequencies (0 = not applied) |
| `qc.snr_g` | `4.0` | empirical | first gate: the largest SNR of a frequency must reach this |

## 9. export: curves passed to the inversion (`05_dispersion/curves`)

| key | default | basis | meaning |
|---|---|---|---|
| `tmin_s` | `auto` | rule | shortest period passed on; auto: 1 / dispersion.fmax_hz |
| `tmax_s` | `auto` | rule | longest period passed on; auto: 1 / dispersion.fmin_hz |
| `one_segment` | `true` | literature | one continuous segment per pair, the longest (as curves are picked in Yao's EGFAnalysisTimeFreq) |
| `curve_figures` | `30` | choice | single-pair figures drawn: this many spread over distance plus the conspicuous ones (0 none, -1 all) |

## 10. inversion: DSurfTomo (`06_inversion`)

| key | default | basis | meaning |
|---|---|---|---|
| `wave` | `rayleigh_phase` | choice | rayleigh_phase (the only one supported) |
| `median_min_pairs` | `50` | empirical | pairs a period needs to enter the average curve of the 1-D model |
| `depth_sequence_km` | `auto` | rule | depth nodes; auto: steps of max(shortest wavelength / 4, depth / 4) down to half the longest wavelength (rule of thumb, Foti et al. 2018) |
| `fine_step_km` | `0.0` | choice | extra depth nodes every this many km down to fine_zmax_km (0 = none) |
| `fine_zmax_km` | `0.0` | choice | see fine_step_km |
| `grid_spacing_deg` | `auto` | rule | grid spacing [lat, lon]; auto: about half the median nearest-station distance after refinement |
| `grid_refine` | `2` | empirical | every grid interval is divided by this once the extent is fixed |
| `line_strip_cells` | `null` | empirical | lines: cells kept on each side of the stations across the line |
| `start_model.source` | `invert` | choice | invert: 1-D model from this array's average curve; file: depth Vs |
| `start_model.file` | `null` | choice | the starting model of your own (two columns: depth km, Vs km/s) |
| `start_model.smooth` | `0.01` | empirical | 1-D inversion: weight of the model curvature |
| `start_model.damp` | `0.05` | empirical | 1-D inversion: damping |
| `start_model.iterations` | `14` | empirical | 1-D inversion: iterations |
| `start_model.minthk` | `3` | technical | 1-D forward: sublayers per node interval |
| `weight` | `1.0` | empirical | DSurfTomo smoothing weight (Fang et al. 2015); set by hand in the usual range (below 10, mostly below 5), not optimised; lines 1.5 |
| `damp` | `1.0` | empirical | DSurfTomo damping (Fang et al. 2015); set by hand in the usual range 0.5-1.5, not optimised; lines 1.5 |
| `sublayers` | `3` | literature | DSurfTomo: sublayers between depth nodes for the kernels (its manual: 2-5) |
| `vmin` | `auto` | rule | lower Vs bound; auto: 0.6 x the slowest Vs of the starting model |
| `vmax` | `auto` | rule | upper Vs bound; auto: 1.4 x the fastest |
| `iterations` | `10` | empirical | DSurfTomo iterations |
| `sparsity` | `auto` | rule | DSurfTomo: fraction of its matrix reserved for nonzero entries; auto: twice the filling estimated from the paths and the grid, at least 0.05 (too little makes DSurfTomo stop with a segmentation fault; more only costs memory) |
| `threshold` | `3.0` | empirical | DSurfTomo: residuals beyond this many quartiles get zero weight |
| `noise` | `auto` | rule | Gaussian noise added to the synthetic travel times of the checkerboards, a fraction of each time; auto: the relative rms of the residuals of the final main model that DSurfTomo's rejection rule keeps, so that the synthetic data are as noisy as the real data fit (Rawlinson & Spakman 2016); the checkerboards then run after the main project |
| `solver` | `dsurftomo` | choice | dsurftomo (<= 60 periods) or pdsurftomo (parallel, <= 80 periods) |
| `threads` | `4` | technical | threads per project |
| `njobs` | `1` | technical | projects run at the same time |

## 11. checks: checkerboards, controls, precheck

| key | default | basis | meaning |
|---|---|---|---|
| `checkerboards` | `auto` | rule | horizontal checkerboards [{name, cell_km}]; auto: cells of 0.12, 0.18 and 0.25 x the aperture (lines: 0.165, 0.24, 0.34 x the line length, along the line) |
| `vertical` | `auto` | rule | depth-varying checkerboards; auto: cell 0.18 x aperture, sign reversals at 0.21, 0.46, 0.71 x the deepest updated node |
| `amplitude` | `0.1` | empirical | checkerboard perturbation, a fraction of the 1-D Vs of each depth |
| `controls` | `[start_grad0.05]` | empirical | control inversions: reg_w<w>_d<d> other regularisation, start_grad<f> other starting model |
| `precheck.min_kept` | `0.5` | empirical | a project is invalid when DSurfTomo's rejection rule would keep less than this fraction of the data in the first iteration |
| `precheck.max_median_iqr` | `1.0` | empirical | ... or when the median residual is beyond this many interquartile ranges |
| `iteration_residuals` | `false` | choice | also evaluate the model of every iteration (one forward run each) |

## 12. plots: figures

| key | default | basis | meaning |
|---|---|---|---|
| `formats` | `[png]` | choice | figure formats |
| `dpi` | `300` | choice | resolution of the png figures |
| `topography` | `auto` | choice | relief under the maps; auto: SRTM 1-arc-second tiles (Farr et al. 2007) from the GMT data server, kept in ~/.cache/lant/srtm1 (only the tile names are requested; maps up to 3 degrees wide); or a grid file / GMT remote grid; null: no relief |
| `vs_range` | `auto` | rule | colour range; auto: 2nd-98th percentile of the model under the stations, to 0.1 km/s |
| `vs_cpt` | `seis` | choice | GMT colour table |
| `slice_depths_km` | `auto` | rule | depths of the slices; auto: up to 8 updated depth nodes, evenly chosen |
| `sections` | `auto` | rule | sections [{name, start, end}]; auto: two through the array centre along its principal axes (areal) or along the line |
| `section_zmax_km` | `null` | choice | depth of the sections (default: the deepest updated node) |
| `contour_step` | `auto` | rule | contour interval; auto: about a tenth of the colour range |
| `font_scale` | `1.0` | choice | font size factor |
| `pair_example` | `null` | choice | pair(s) of the single-pair figures (stack time, dispersion): a name or a list; default: three pairs at the 20 / 50 / 80 % quantiles of the distances of the pairs with long curves |


## 13. Defaults that differ for line arrays

A line array (`project.array_type: line`, or `auto` on an elongated array) uses:

- `data.min_trace_length_s`: `600.0`
- `reference.fj_mode`: `track_fundamental`
- `dispersion.njack`: `auto`
- `inversion.median_min_pairs`: `20`
- `inversion.grid_refine`: `1`
- `inversion.line_strip_cells`: `1.0`
- `inversion.weight`: `1.5`
- `inversion.damp`: `1.5`

## References

- Bensen, G. D., et al. (2007). Processing seismic ambient noise data to obtain reliable broad-band surface wave dispersion measurements. *GJI*, 169, 1239–1260. doi:10.1111/j.1365-246X.2007.03374.x
- Brocher, T. M. (2005). Empirical relations between elastic wavespeeds and density in the Earth's crust. *BSSA*, 95, 2081–2092. doi:10.1785/0120050077
- Constable, S. C., Parker, R. L., & Constable, C. G. (1987). Occam's inversion. *Geophysics*, 52, 289–300. doi:10.1190/1.1442303
- Fang, H., Yao, H., Zhang, H., Huang, Y.-C., & van der Hilst, R. D. (2015). Direct inversion of surface wave dispersion for three-dimensional shallow crustal structure based on ray tracing. *GJI*, 201, 1251–1263. doi:10.1093/gji/ggv080
- Farr, T. G., et al. (2007). The Shuttle Radar Topography Mission. *Rev. Geophys.*, 45, RG2004. doi:10.1029/2005RG000183
- Foti, S., et al. (2018). Guidelines for the good practice of surface wave analysis. *Bull. Earthquake Eng.*, 16, 2367–2420. doi:10.1007/s10518-017-0206-7
- Herrin, E., & Goforth, T. (1977). Phase-matched filters: application to the study of Rayleigh waves. *BSSA*, 67, 1259–1275.
- Jiang, C., & Denolle, M. A. (2020). NoisePy: a new high-performance Python tool for ambient-noise seismology. *SRL*, 91, 1853–1866. doi:10.1785/0220190364
- Jing, J., et al. (2025). Ambient noise seismic tomography of Tonglushan skarn-type Cu-Fe-Au deposit in Eastern China. *Ore Geology Reviews*, 184, 106718. doi:10.1016/j.oregeorev.2025.106718
- Lévêque, J.-J., Rivera, L., & Wittlinger, G. (1993). On the use of the checker-board test to assess the resolution of tomographic inversions. *GJI*, 115, 313–318.
- Levshin, A. L., & Ritzwoller, M. H. (2001). Automated detection, extraction, and measurement of regional surface waves. *Pure Appl. Geophys.*, 158, 1531–1545. doi:10.1007/978-3-0348-8264-4_11
- Leys, C., Ley, C., Klein, O., Bernard, P., & Licata, L. (2013). Detecting outliers: do not use standard deviation around the mean, use absolute deviation around the median. *J. Exp. Soc. Psychol.*, 49, 764–766. doi:10.1016/j.jesp.2013.03.013
- Lin, F.-C., Moschetti, M. P., & Ritzwoller, M. H. (2008). Surface wave tomography of the western United States from ambient seismic noise: Rayleigh and Love wave phase velocity maps. *GJI*, 173, 281–298. doi:10.1111/j.1365-246X.2008.03720.x
- Luo, Y., Yang, Y., Xu, Y., Xu, H., Zhao, K., & Wang, K. (2015). On the limitations of interstation distances in ambient noise tomography. *GJI*, 201, 652–661. doi:10.1093/gji/ggv043
- Rawlinson, N., & Spakman, W. (2016). On the use of sensitivity tests in seismic tomography. *GJI*, 205, 1221–1243. doi:10.1093/gji/ggw084
- Seats, K. J., Lawrence, J. F., & Prieto, G. A. (2012). Improved ambient noise correlation functions using Welch's method. *GJI*, 188, 513–523. doi:10.1111/j.1365-246X.2011.05263.x
- Sens-Schönfelder, C. (2008). Synchronizing seismic networks with ambient noise. *GJI*, 174, 966–970. doi:10.1111/j.1365-246X.2008.03842.x
- Shapiro, N. M., & Singh, S. K. (1999). A systematic error in estimating surface-wave group-velocity dispersion curves and a procedure for its correction. *BSSA*, 89, 1138–1142.
- Shaw, P. R., & Orcutt, J. A. (1985). Waveform inversion of seismic refraction data and applications to young Pacific crust. *Geophys. J. R. astr. Soc.*, 82, 375–414.
- Stehly, L., Campillo, M., & Shapiro, N. M. (2007). Traveltime measurements from noise correlation: stability and detection of instrumental time-shifts. *GJI*, 171, 223–230. doi:10.1111/j.1365-246X.2007.03492.x
- Takagi, R., Nishida, K., Aoki, Y., et al. (2015). A single bit matters: coherent noise of seismic data loggers. *SRL*, 86(3), 901–907. doi:10.1785/0220150030
- Taylor, J. R. (1997). *An Introduction to Error Analysis* (2nd ed.), chapter 6 (Chauvenet's criterion). University Science Books.
- Wang, J., Wu, G., & Chen, X. (2019). Frequency-Bessel transform method for effective imaging of higher-mode Rayleigh dispersion curves from ambient seismic noise data. *JGR Solid Earth*, 124, 3708–3723. doi:10.1029/2018JB016595
- Xu, Y., Lebedev, S., Meier, T., Bonadio, R., & Bean, C. J. (2021). Optimized workflows for high-frequency seismic interferometry using dense arrays. *GJI*, 227, 875–897. doi:10.1093/gji/ggab260
- Yao, H., van der Hilst, R. D., & de Hoop, M. V. (2006). Surface-wave array tomography in SE Tibet from ambient seismic noise and two-station analysis – I. Phase velocity maps. *GJI*, 166, 732–744. doi:10.1111/j.1365-246X.2006.03028.x
- Zhu, S., Li, J., Chen, G., Fang, H., & Yao, H. (2026). pDSurfTomo: a high-performance parallel computing package for direct surface wave tomography. arXiv:2604.11920.
