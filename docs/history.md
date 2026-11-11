# History

[中文](history.zh.md) | English

## v1.4 (2026-11): corrections after trial runs on real data, complete documentation

After trial runs from raw records on one line and one areal array of real data, the problems found were fixed; every change comes with a synthetic test ([validation.md](validation.md)).

**Correctness of the computation**
- **Records without a station code (SEG-Y, SEG-2, some SAC files)**: they were named from their file names only after the vertical channel had been chosen station by station, so all of them were taken as one station without a name and nothing was cleaned. The names are now taken first, and the records are read back with the code of the file. A file with several such records at the same time (one trace per channel) cannot be named from its file name and is reported as not used.
- **Checkerboard noise from the data** (`inversion.noise: auto`): the relative rms of the residuals of the final main model kept by DSurfTomo's rejection rule, in place of a fixed 2 %; the checkerboards then run after the main project. On the real line the kept residuals have a relative rms of 4.8 % (all residuals 8.6 %), so a fixed 2 % made the checkerboards look better than the data allow.
- **Inter-station distances on the WGS84 ellipsoid** (inverse solution of Vincenty 1975). On the former sphere of 6371 km the distance error depends on azimuth (−0.18 % to +0.30 % at latitude 30°) and has the form of a 2θ azimuthal anisotropy. DSurfTomo receives velocities, so no extra coordinates are needed and its source is unchanged.
- **Band-edge criterion redone**: it used the amplitude spectrum of each pair, and at high signal-to-noise ratio every point beyond the band of the noise sources was accepted. It now uses the amplitude spectrum of the coherent signal (energy spectrum of the signal window minus that of late lags with incoherent noise only, all pairs together, computed once); the result does not depend on the stacking time or the noise level.
- **1-s pulse test**: the F test uses effective degrees of freedom. Red and band-limited noise sampled at 250 Hz gave 7–15 % false detections before (nominal 0.1 %); none in 200 trials now.
- cc reads the records of neighbouring days, so segments across a day boundary of the archive are no longer lost; `remove_1s_comb` with a non-integer sampling rate is an error.
- Periods are written with 6 decimals (exported curves, DSurfTomo.in) instead of being rounded together.
- **Sparse reserve of DSurfTomo**: with too small a reserve DSurfTomo writes out of bounds and crashes. `inversion.sparsity` defaults to `auto`, estimated from the grid cells crossed by straight paths; a value given by hand that is too small is an error.
- Unquoted dates in YAML no longer break the step reports.
- Without the archive (`clean.write_archive: false`) the 1-s periodic signal is not tested, which the report did not say; it is now marked as not checked.

**Simpler methods**
- `dispersion.qc.dev_c` and `smooth_c` removed: leaving them out changed nothing in four tests with noise and a second arrival.
- The checkerboard pattern of a line is odd about the middle of the grid: on a real line of 15 km the former pattern, starting at the grid edge, gave synthetic residuals leaning to one side for two of the three cell sizes, and the precheck marked them invalid.
- One default control project only, the starting-model control `start_grad0.05` (areal and line alike).
- Station diagnosis: when there are no more candidates than expected by chance under Chauvenet's criterion (1), the chain does not stop; they are listed in the report. A time-shift candidate must also exceed the standard error of the median shift over the pairs of the station, so poorly measured shifts of low-correlation sub-stacks are no longer listed in bulk.
- On a line with `data.velocity_range_kms` left at the default, the reference step lists this among the items to judge (on the real line the default gave a reference over a narrow band and medians 4–9 % higher).
- F-J ridge tracking for lines: a start every 1/3 octave is tried and the one followed over the most frequencies is kept (the slower on a tie); at each frequency the strongest peak within the allowed window is taken, the window widening with the frequency distance from the last pick (a limit on the slope of the dispersion); picks are restricted to `data.velocity_range_kms` and to wavelengths of at least twice the mean station spacing (sampling theorem: shorter wavelengths give aliased ridges in the F-J spectrum of a line); fewer than 5 frequencies is an error with ways forward, instead of a reference curve from one or two picks. The automatic start band was rounded and could exclude the frequency where the start had been found; fixed.

**Usability**
- Sections of areal arrays: the elevation profile above a section is interpolated between all stations (inside their convex hull), not drawn only between the few stations projected onto the section.
- Synthetic examples: `--responses` records alternate stations with two geophone types and writes their responses; `--directional F` sends a fraction F of the plane waves from one direction.
- Lower bounds of the Python packages set to versions that were run: from Python 3.8 with NumPy 1.19, SciPy 1.5, ObsPy 1.2, h5py 2.10, Matplotlib 3.5, PyYAML 6.0, openpyxl 3.1 to Python 3.10 with NumPy 2.2, SciPy 1.15, ObsPy 1.5, h5py 3.16, Matplotlib 3.10.
- clean can resume: station-days with unchanged input files and settings are reused; the header scan is cached; reader warnings about non-data blocks are only counted.
- Changing only computing resources (processes, threads, pairs per task, disk margin) no longer makes finished steps out of date or voids their confirmation.
- The clean report estimates the memory per cc process and the disk space of the hourly correlations; cc reduces the number of processes when memory is short.
- Station table: deployment times may be dates only (a whole day in the time zone of `data.day_label`); more Chinese column names are recognised; with start and end times but no serial number, the station name serves as the serial.
- Sections are drawn at the true length-to-depth ratio of the model (within limits; a vertical exaggeration is stated otherwise).
- **Figures redone in the layouts usual in ambient-noise papers**: relief (SRTM 1", downloaded once and cached) under the station map, the slices and the checkerboards; correlation moveouts per frequency band, with distances, lag window and reference velocities taken from the measurements of the band; three single pairs instead of one; the phase panel of the dispersion figure shows the measured phase velocity with its other cycle branches and the reference (the classical narrow-band phase image disagrees with the phase-matched measurement where the dispersion is strong); checkerboards red – white – blue with boxed labels; a 3-D view of the slices for areal arrays; section ends marked with letter and direction, the elevation profile taken from the relief; the residual figure compares the starting and the final model.
- Figures of lines: the checkerboards are drawn as input and recovered sections along the line, and depth slices are left out (the model is a narrow strip; drawn as maps they were some fifteen times taller than wide); the station map is widened to at least 0.4 of its length.
- The evaluate report lists starting and final residuals, the recovery slope and correlation range of each checkerboard, and the model difference of the control project with its residual history (numbers only).
- Quick synthetic example `make_synthetic.py --quick` (10 stations, 6 hours, about 6 minutes for the whole chain); continuous integration runs the tests and the quick example.

**Code**
- Modules split by responsibility (about 500 lines as the point to check, AGENTS.md): station matching and the geometry checks moved from `clean.py` to `station_match.py`; the 1-s pulse test and its removal from `clean.py` and `cc.py` to `pulse.py`; the F-J spectrum and ridge picking from `reference.py` to `fj.py`; the step reports (`report.py`) and the figure set (`plotting/figures.py`) are written by one function per step or group of figures. The tests are split by module (`tests/test_<module>.py`), with new tests of the figure choices (bands, windows, example pairs, relief tiles) and of the checkerboard noise.

**Documentation**
- New [methods.md](methods.md) (principles, formulas, sources), [outputs.md](outputs.md) (output files), [troubleshooting.md](troubleshooting.md), [quickstart.md](quickstart.md); the tutorial was rewritten for this version.

## v1.3 (2026-11): three paths to start, a literature basis for every method

- **Only three required paths** (raw data, station table, output directory). `lant run all --raw … --stations … --out …` writes the configuration and runs; `lant init` no longer has separate templates for areal and line arrays. The array type follows from the station distribution (elongation below 3 areal, above 10 line, in between you confirm); day labels are UTC by default; the velocity range defaults to 0.1–5 km/s and the reference step reports the measured velocities.
- **Instrument response**: with response files (`data.response`) the clean step removes the response (Bensen et al. 2007); without, one instrument model is required and the remaining error is stated in the documents and reports. The formerly required `data.instrument` is gone.
- **The phase is read with a phase-matched filter** (Herrin & Goforth 1977; Levshin & Ritzwoller 2001) instead of the instantaneous-frequency correction. New band-edge criterion `dispersion.qc.band_edge` (reason `off_band`).
- **The F-J spectrum is computed by lant itself in double precision** (Wang et al. 2019); CC-FJpy is no longer needed. Without a given velocity there are two passes: the first takes every pair and picks a curve, the second limits the pairs by the wavelengths of that curve. New zone references `reference.zone_km`.
- **One station diagnosis instead of two**: leave-one-out correlation and sub-sample time shift (Stehly et al. 2007; Sens-Schönfelder 2008), limits from the sample size by Chauvenet's criterion and the MAD (Leys et al. 2013); candidates are listed, removal is yours to confirm. The former station-day and station-hour screens are gone.
- **1-s pulse**: the clean step tests every station-day with an F test on the record folded at 1 s (Takagi et al. 2015); the cc step subtracts it only at stations where it was found. **The notch step is gone**, and so is the second reference pass (pass2).
- **Stops only for scientific judgements**: a step with nothing to judge does not stop; `lant status --json` gives state, reason and next command; `lant plan` lists what every step will do before it runs.
- **Moved nodes and several lines**: serial, start/end and line columns in the station table; `data.line`; one project per line. Station names are taken from file and folder names when the headers have none.
- **Large arrays**: `correlation.max_distance_km` limits the pair distance; memory is estimated before cc; the stack is written day by day.
- Sections of your own that leave the model are cut to it, with a note. The synthetic noise of the checkerboard tests no longer has a fixed random seed (DSurfTomo is no longer patched for it).
- Every configuration key carries the kind of its basis (literature, theory, rule, empirical, choice, technical); `docs/parameters` is generated by a program; citations that could not be verified were removed.
- All documents exist in Chinese and English; new document `docs/validation` with the synthetic tests (phase-matched filter, F-J, J0 benchmark, laterally varying medium, line array).

## v1.2 (2026-11): a shorter chain, confirmed by a person at the key points

v1.1 left every parameter to the program and added default steps for particular data problems; the chain grew complicated and nobody looked when an automatic value was off. This version was settled item by item after two reviews (science, simplicity):

- **The probe step was removed.** Parameters before the correlation follow from a rough velocity range `data.velocity_range_kms` by theoretical relations; the lower band limit and the velocity axis follow from the measured F-J curve. The probe's "group repeatability + extrapolated signal-to-noise ratio" had no precedent in the literature.
- **Required keys**: array type, velocity range, time zone, instrument model (reduced to three paths in v1.3).
- **Four confirmation points**: after clean, stack, reference, inversion; `lant confirm <step>` writes the accepted values into the `confirmed:` section of the configuration; a confirmation lapses when its result changes. Priority: user's value > confirmed value > rule.
- **Diagnose only, a person decides**: station-day / station-hour screens, notch and 1-s pulse removal are off by default and listed as candidates by the stack step.
- **One measurement pass**: disp1 and disp2 merged into disp; the second reference pass became optional.
- **Group-velocity products and redundant-path precision removed**: group velocity only places the phase reading; only phase velocities are exported, with jackknife errors (not calibrated). Every rejected sample records its reason (bit mask).
- **Checks**: regularisation controls at 0.5 and 2 times the main values; the starting-model control is a linear gradient (+5 % at the surface to −5 % at the deepest node), centred by the travel-time residuals. The former overall × 1.05 control lost all its data to DSurfTomo's rejection rule in the first iteration; a precheck (`checks.precheck`) now stops an invalid main project and marks invalid checks. DSurfTomo is not changed.
- **Engineering**: `lant status` compares resolved parameters and shows pending confirmations; the write lock uses flock.
- The configuration template lists the dozen keys changed most often.

## v1.1 (2026-11): every parameter derived from the data

Many parameters of v1.0 were still the values of the author's one array. Changes:

- New step **probe**: correlation and F-J on a sample of the data with a wide velocity range, from which band, sampling rate, lag window, segment length, velocity search range, F-J parameters and screening settings were derived (removed again in v1.2).
- Sub-stack length from the record length; grid spacing from the station spacing.
- Depth nodes from the average curve (steps of the larger of a quarter of the shortest wavelength and a quarter of the depth, down to half the longest wavelength); velocity bounds from the starting model; colour scale, slice depths and contour interval from the final model.
- Absolute quantities became relative ones (taper of the coherence band, error bands per octave, checkerboard statistics by station spacing, and others).
- Fixes: a fixed lower band limit of 0.2 Hz and a fixed shortest export period could contradict the derived band; only six candidate sampling rates; several GB of memory for one station-day at high rates; a nearest-distance calculation growing with the cube of the station number. The 1-D forward calculation passed periods to surfdisp96 in the given order although it needs increasing periods: with decreasing periods it followed other roots in strongly layered models; periods are now sorted first.
- The smoothness criterion `qc.smooth_c` became relative to the reference curve (the old form limited the dispersion slope and rejected true curves in strongly dispersive areas).
- Fixes after a code review: coordinates in scientific notation and seconds of arc; long station names cut in miniSEED; `data.rename`; single horizontal channels; ambiguous table entries; non-integer decimation ratios; NaN samples; station-day screening with sub-day blocks; circular shift search; a silently lowered minimum of jackknife groups; rejected phase velocities used in the group-velocity wavelength rule; group SNR taken from the whole image; borrowed curves interpolated across invalid spans; a possible endless loop with manual depth nodes; the half-space missing from the layer limit; an unreadable error above 80 periods; completion judged from log text only; check projects with stale input; control models compared without checking the grids; the deepest updated layer missing in figures; incomplete configuration checks; non-standard NaN in records; configuration history overwritten; `validate` exit code 0 on problems; thread numbers overridden by environment variables.
- Checkerboard amplitude ±10 % of the background Vs (was ±0.4 km/s); the starting-model control scales updated nodes only; `data.component` (Z only); sensors pointing down are flipped.
- Traceability: code version and environment per step (`records/provenance/`), list of input files, one run per output directory (`records/lock`).
- The upper band limit is the highest frequency the processing rate `correlation.fs_max` allows; the quality control decides how far the data reach. The correlation band is widened to the half-gain points of the Gaussian filters at its ends.
- Step-by-step running: every step prints a report and stores it as `records/steps/<step>.md`; `lant status` compares the parameters of the last run with the configuration.
- The reference curve constrains the measurement only where picks support it; beyond, as in FTAN, the group velocity is the envelope maximum over the whole velocity range and the phase velocity continues the cycle count.
- The lower band limit follows from theory: the frequency where the longest pair is 2 wavelengths.
- Day labels can be configured with `data.day_label`; `start_date` / `end_date` select by label; `01_clean/processing_days.txt` lists the UTC window of every label.
- The smoothing of the 1-D starting model constrains the model itself (jumping, Constable et al. 1987), with the curvature computed from the real depth intervals; before, the result drifted with the number of iterations.
- Normal dispersion in `track_fundamental` became the switch `reference.normal_dispersion`.
- Fixed random seed for the synthetic noise of the checkerboards (removed again in v1.3). More residual statistics in `evaluate`.
- The stack-time figure uses the real signal and noise windows.
- The parameter document says for every item whether it is a literature value or a project value.
- The synthetic generator got a soft-sediment model (`--model sediment`).

## v1.0 (2026-11): rewritten as a general version of the high-frequency workflow

**Why**: the original local_ANT (NoisePy correlations + neural-network picking with EGFAnalysisPy/DisperPicker + DSurfTomo) had been replaced in the author's dense-array work by a workflow with a basis for each step. This version turned that workflow into a general project: no dependence on the author's directories, file names or machine; parameters in the configuration or derived from the data; recognition and cleaning of raw data in any format; checks, figures, documents and an AI procedure.

**Versions brought together**:

| Version | Content | Treatment |
|---|---|---|
| old `main` (6f70b3e) | original workflow + PR #1 | tag `legacy-main`; method retired, not merged |
| old `refactor/project-structure` (55db2f8) | numbered restructuring of the original workflow | tag `legacy-refactor`, branch deleted; not merged |
| the author's unpushed local changes | plotting improvements and project configurations of the original workflow | not merged; sections and slices are implemented anew in `lant/plotting` |
| the author's current high-frequency workflow | source of the methods of the new version | ported |

**Check of the port**: on a subset of real dense-array data (areal array: 30 stations, 3 days; line: 30 stations, 4 days) the author's scripts and the new version were run with the same parameters and compared step by step:

| Stage | Result |
|---|---|
| hourly cross-coherence | first 23 hours of each day agree (relative difference < 3e-8); in the 24th hour the new version uses one more 60-s segment, because it reads all files of a time window and the old scripts read one file per day |
| station-day screen, notch table | identical |
| F-J reference, first and second pass | bit-identical (same input) |
| both measurement passes and quality control (areal, line) | all values and flags bit-identical (same input) |
| export | curve files byte-identical |
| grid, `surfdata.dat`, `DSurfTomo.in`, 1-D starting model | same as the original tools |
| parameters derived from the data (`lant/auto.py`) | for the full data of the original study: the same day start, sampling rate, band, lag, grid and checkerboards as the study used |

**Deliberate changes**:

- Line station-hour screen: stations rejected as a whole are removed before scoring.
- Line export as for areal arrays: measured periods, one segment per pair.
- A line can use its own reference curve and starting model, or borrowed ones.
- Checkerboard cells are given in kilometres.
- Kernel curves by finite differences with fwd1d (DSurfTomo's own forward chain, Vp and density following Vs after Brocher); CPS330 is not distributed with this project.
- Third-party programs are not in the repository; `tools/install_external.sh` downloads fixed versions and builds them; the lower memory reserve of pDSurfTomo is a patch.
- The memory queue for a shared server became the two parameters `inversion.njobs` and `threads`.
- The user gives only the raw-data directory and the station file; format, layout, vertical component, table form and name matching are recognised by the program, and the remaining parameters follow from the data by written rules (`records/auto_parameters.json`).

## Productization (unreleased)

New projects use `work/`, `results/` and `.lant/`, with README and total report entry points. Solver inputs, models and checks are retained separately from working caches. Internal archive paths and configuration snapshots are portable. Progress uses stderr and bounded detail logs; inspect/status no longer write inspection snapshots. No measurement or solver method changes.

Compact `lant init` is now the default; `--full` retains all parameter descriptions. PNG is the default presentation format; add PDF explicitly with `plots.formats`. The dispersion CSV preserves legacy solver rounding and selection semantics.
