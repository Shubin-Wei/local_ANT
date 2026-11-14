# History

[中文](history.zh.md) | English

## 2.0.0 (2026-11-14): Python only, automatic resources, calibrated regularization

- Inversion: the fast marching of each source stops once the nodes its receivers read are done (outputs unchanged bit for bit; travel times about twice as fast).
- Clean: a miniSEED file whose first and last records belong to a station that no matching rule can give a station of the table is not read in full (its header codes stay in the index for the matching); a 55-station line in a 181-station raw directory: 478 → 144 s, the archive byte-identical.
- Cross-correlation: the linear detrend of each segment is compiled (the same LAPACK least squares as SciPy, bit for bit) and the segments of an hour run in threads; one day of a large areal array 168 → 103 s with four threads, outputs identical.
- Stack: the hourly sums run compiled and in parallel over the pairs (the same numbers bit for bit; 8 days of a large areal array 293 → 54 s), and the run after `lant confirm stack --exclude-candidates` reuses the stack of every hour of the last run when its inputs are unchanged (`all_hours.json`) instead of stacking every hour again.

- Inversion: the depth kernels are DSurfTomo's own surfdisp96 and depthkernel, ported statement by statement in their own precision (`dsurf_disp96`), instead of disba. The port's outputs are now byte-identical to native DSurfTomo on four test cases (python_backends.md); before, the small kernel differences moved a few data across the rejection bound and the models differed by up to a few m/s. The kernels are also faster (a 55-station line, 20 threads: 34.6 → 26.6 s per iteration). The 1-D forward calculation (starting model, synthetic examples, kernel figure) uses the same port; disba is no longer a dependency, so the starting models and the synthetic examples change in the last digits.

- Cores: `resources.cores: auto` takes 80 % of the idle cores (was half), max(1, 0.8 × (usable CPUs − 1-minute load of others)).

- Regularization: the smoothing weight and the damping are fixed numbers, weight 1 and damping 0.1 for areal arrays and lines (`auto` is refused), calibrated on lant's areal and line synthetic examples against the true model (methods.md 8.4). The former damping 1 (lines 1.5) held every update to a few hundredths of a km/s, so that the checkerboards kept 3–12 % of their amplitude; the automatic weight (largest first-iteration change 0.3 km/s) reached its target on neither example and is removed with `regularization.json`. The combined rule w = d = max(f_B, f_M) with its discrepancy target from the uncalibrated jackknife errors is removed, with `surfdata_dc.txt`. Controls `reg_x<f>` (weight times f; e.g. `reg_x0.5`, `reg_x2`, not run by default) show how much the model depends on the weight.

- Grid and checkerboards: the grid spacing is the median measured wavelength c T / 3 in both directions, set by the inversion step (was half the median nearest-station distance, set by the clean step); `checks.checkerboards: auto` gives three horizontal tests with cells of 1, 2 and 3 median wavelengths, and the depth-varying test has the 2-wavelength cell (was one test of max(wavelength, 2 station spacings)).

- Figures: checkerboards and anomaly sections are red (slow) – white (zero) – blue (fast) (GMT polar, reversed); Vs keeps `plots.vs_cpt`. The note under a section mentions topography only when it is drawn, and no longer overlaps a single colour bar.

- REPORT.md: the steps table gives the time of each step; the key numbers give the smoothing weight and damping and count the figures actually in `results/figures`; the figure folder is no longer listed file by file. Step reports show project paths relative to the output directory, and the resource lines name the projects of each solver batch. Wording only; no processing changes.

- Move the matching of record names to the station table (`names_from_files`, `assign_by_serial`, `match_names`) from `lant.stations` into `lant.station_match`; `lant.stations` keeps reading and normalising the tables. The moved definitions are AST-identical; no processing changes.

- Move the station diagnosis of the stack step from `lant.stack` into `lant.station_qc` (stacking and removal stay in `lant.stack`). Every moved definition is AST-identical; no processing changes.

- The forward calculation of the solver no longer holds the derivative field of every travel time at once (rows x nx x ny single precision, several GB for large arrays and fine grids): sources are taken in blocks of at most 512 MB of fields (`dsurf_forward.FIELD_BLOCK_BYTES`) and each block is turned into its matrix rows before the next, as DSurfTomo turns each receiver's field into its row. Models, residuals and DWS of the checkerboard and the measured-line comparison are byte-identical, also with 2 MB blocks.

- A 3-D project again carries at most 60 periods, the NP of CalSurfG in DSurfTomo a651faf (as in 1.5.0); 2.0.0 had raised it to the 80 of surfdisp96, which the native solver does not accept, so such a project could not be compared with DSurfTomo. The 1-D forward calculation of the starting model keeps 80.

- Move raw-file discovery and cached header indexing into `lant.record_index`; retain the `lant.data` entry points. This responsibility split does not change header parsing or numerical processing.

- Separate step-specific report narratives into `lant.step_descriptions`, leaving parameter provenance and report assembly in `lant.report`; report wording and scientific decisions remain unchanged.

- Move the station-day archive worker and its existing reuse key into `lant.clean_archive`; `lant.clean` retains whole-step orchestration. Processing statements and archive version are unchanged in this refactoring.

- Separate model I/O, forward evaluation, residuals and recovery into `lant.inversion.model_evaluation`; retain public imports from `lant.inversion.run`, which now focuses on solver scheduling. Numerical statements are preserved; a redundant local JSON import is removed and cache tests patch the new owning module.

- Expose the selected automatic-regularization trial's existing LSMR stop code in its JSON and step report. The native implementation reuses code 3 for damped least-squares tolerance and a condition limit; the report preserves that ambiguity and labels the machine-precision condition and iteration limits separately; the solver controls, search, selected weights and numerical operations are unchanged. Clarify that the manual's 0.2–0.4 guidance refers to the first few iterations.

- Reuse unchanged header indexes without starting a worker pool or rewriting the cache. Cache signatures now include nanosecond modification/change times and file identity; old signatures trigger one rescan. Same-size edits within one second no longer reuse stale headers. This is an engineering cache rule, not a raw-data content checksum.

- Apply the same metadata identity to cleaned station-day reuse and clean-input provenance. Instrument-response inventories and cleaned archives also depend on response-source identity (including SHA256 contents), so replacing a response at the same path invalidates them and marks clean stale. Unchanged files keep their existing processing; upgraded cache keys cause one rebuild.

- Render independent main figure groups in bounded spawned processes under the automatic core/memory plan, with one numerical thread each. Preserve private Matplotlib/GMT state, staged file publication and old images on group failure. Record per-figure elapsed times. The existing 17-image example produced byte-identical PNGs; this changes scheduling, not figure content or scientific calculations.

Output directories of earlier versions are not continued: run the chain again into a new output directory (lant
says so when it meets one). A configuration with a removed key stops with what replaced it.

- **One solver, in Python.** The inversion is the Python port of DSurfTomo (`lant.inversion.dsurf_main`, run per
  project by `python -m lant.inversion.solver`). DSurfTomo and pDSurfTomo are no longer built or called (`lant
  install-external` and the key `inversion.solver` are gone), nor fwd1d and fwdtimes: the starting model and all
  phase velocities come from disba 0.7.0. The Fortran sources, the PyKonal candidate and every tool that compared
  with or linked native code are removed from the repository; the comparison of the port with native DSurfTomo is
  documented in [python_backends.md](python_backends.md) with figures of a synthetic case. Only Rayleigh phase
  velocities are supported.
- **Wavelet L1 removed.** A sparsity (wavelet L1, IRLS) regularization was tried as an alternative to DSurfTomo's
  smoothing: in four of six synthetic cases its model was further from the truth than the starting model. It is
  removed; the inversion keeps DSurfTomo's L2 smoothing and damping.
- **Weight and damping chosen from the data** (`inversion.weight` and `damp`: auto). On the first iteration only:
  w = d = f with f = max(f_B, f_M), f_B where the largest velocity change is 0.3 km/s (DSurfTomo manual: 0.2-0.4),
  f_M where the residual equals the travel-time errors of the jackknife uncertainties (discrepancy principle,
  Morozov 1966); the check projects inherit it. A number in the configuration is used as given.
- **Cores and memory chosen by the program.** `resources.cores` (auto) replaces `clean.nproc`,
  `correlation.nproc` and `threads_per_process`, `dispersion.nproc`, `inversion.threads` and `njobs`: a share of
  the cores other users leave idle (80 % from 2.0, see above), taken again at every step and solver project, and parallel work within half of the
  available memory (was 70 %). Every step report gives the cores and the peak memory it used.
- Checkerboard recovery numbers use nodes with nonzero post-rejection data-Jacobian column support (DWS > 0; the solver writes
  `DWS.dat`); the maps stay whole.
- The built-in synthetic example (`lant demo`, `python -m lant.synthetic`) uses a ±6 % sinusoidal checkerboard on a 1-D background (`true_checkerboard.json`); `--uniform` disables it. This follows the user’s choice to keep checkerboard examples and the README overview style; the short-lived slow/fast-body example is replaced. Straight-ray synthesis is approximate and independent of the inversion forward model.
- One conda environment (`environment.yml`), Python 3.11 or newer, NumPy 2; no compilers, make, Cython or PyKonal.
  Linux only.
- A failing parallel batch no longer cancels other batches (each has its own cancellation).
- **Faster hourly correlations.** The pair products of the segment spectra are summed per hour by a
  compiled kernel (`lant/cc_kernel.py`; the same single-precision products added in double precision in segment
  order, spread over the threads of the day's process) instead of NumPy copies of the spectra of every pair; the
  hourly array is written to the day file hour by hour instead of kept in memory; the stations of a day are read
  and decimated by up to four threads straight into one array. The product and summation precisions and segment order are preserved; NumPy and Numba complex products can
  differ by floating-point roundoff, so bitwise identity is not promised. A failed write removes its temporary
  file and keeps the previously completed day. The memory estimate follows the streamed arrays.

- **Faster dispersion measurement.** Full and delete-group stacks share one batched phase-matched FFT and filter
  weights. Analytic filter banks and daily sub-stacks are processed one pair at a time; the worker estimate includes
  these arrays. Picking, jackknife groups and quality-control rules are unchanged.

- **Figure review.** Keep the README overview style (STIX serif, inward ticks, slow-red/fast-blue anomaly palette). Short panel titles, accepted and rejected phase candidates and the reference wavelength gate make QC visible. Recovery captions name nonzero data-Jacobian support (DWS > 0). Workflow output is PNG only; `plots.formats` rejects other formats.

- **Documentation corrections.** State the signed-quartile rejection rule, the distinct correlation/EGF phase conventions, jackknife use in automatic regularization, and the engineering basis of depth spacing and combined regularization. Correct the F-J equation reference and DWS support description. Numerical solver statements are unchanged.

- **Parallel pair diagnostics.** The unchanged set of single-pair PNGs is rendered in spawned processes under the automatic core/memory plan, with progress and worker-memory reporting. Each process owns its Matplotlib state; figure content and QC decisions are unchanged.

- **Cold-start worker fix.** Retuning a worker keeps Numba's import-time thread-pool ceiling and changes its active mask. This prevents a first compilation from rejecting a changed `NUMBA_NUM_THREADS` after the runtime has started; warmed caches had hidden the error.

## 2.0.0 development step: DSurfTomo-compatible inversion driver

`python -m lant.inversion.dsurf_main DSurfTomo.in` ports DSurfTomo's main program: the iterations around the compatible forward, the quartile rejection, the smoothing rows, single-precision LSMR with local reorthogonalization (`dsurf_lsmr`), the velocity update and the native output files in gfortran's formats (`dsurf_format`), including the synthetic branch (`ifsyn = 1`, on CalSurfG's 5 × 5 synthetic grid). With native depth kernels every output file is byte-identical to native DSurfTomo's on the test project. The depth kernels (`dsurf_kernels`) give the same numbers as before with less work. `tools/dsurf_native_kernels.f90`, `tools/dsurf_native_timing.py`, `tools/compare_dsurf_inversion.py` and `tools/plot_dsurf_inversion.py` compare with native DSurfTomo. Not connected to `inversion.solver`; no default or configuration key changes.

## 2.0.0 development step: DSurfTomo-compatible forward port

`lant.inversion.dsurf_forward` and its kernel modules port the forward part of DSurfTomo's CalSurfG (B-spline grids, refined-source fast marching, receiver times, ray paths, row assembly) to Python/Numba, reproducing the reference's REAL arithmetic and its known behaviour near the source. `tools/dsurf_native_dump.f90`, `tools/compare_dsurf_port.py` and `tools/plot_dsurf_port.py` compare it with the unmodified native routines. `sensitivity.node_partials` exposes the separate Vp, Vs and density partials. No default, configuration key or inversion driver changes.

## 2.0.0 development step: results hold what the user reads

- One `results/parameters.yaml` (every value resolved, with the confirmation record), rewritten by each completed step, in place of one yaml per step; the values each step used remain in its report and the manifest.
- `results/reproduce/` is no longer written: every file in it was a copy of `work/06_inversion/main/` or `work/`, where the inputs to rerun stay. An unchanged `results/reproduce/` or per-step `results/parameters/` of an earlier version is removed at the next publication.
- The data-coverage, station-diagnosis and F-J figures are moved, not copied, to `results/figures/` when their step completes.
- `median_vs_reference.txt` stays in `work/05_dispersion/`; `results/dispersion/` holds `curves.csv` only.
- After the evaluation, the DSurfTomo files lant never reads (`*.iterNNN`, `residualFirst.dat`, `residualLast.dat`, `velmap2dRc.dat`) are deleted from the checkerboard and control projects, only when lant registered them and they are unchanged; the main project keeps them.
- The project rules file for assistants is written by its own module, `lant/assistant_rules.py`.

## 2.0.0 development step: checkerboard cells of at least two station spacings

The automatic checkerboard cell is max(median measured wavelength, 2 x median nearest-station distance); the factor 2 is empirical (twice the sampling limit, so that the pattern is sampled by at least two stations per cell). With one station spacing, the synthetic areal example of 16 stations recovered its 1.5 km cells with a regression slope of at most 0.01.

## 2.0.0 development step: tidier results

`results/dispersion/` no longer holds a copy of `disp.h5` (the measurements stay in `work/05_dispersion/`); `results/reproduce/` holds the inputs of the main project only (lant's `checks.json` and the kernel figure data `kernels.txt` are no longer copied); `results/main_preview/` is removed once the evaluated model is published, when nothing in it was changed; checkerboard figures are drawn only for the check projects of the present configuration, the ones the evaluation uses. No processing data are removed.

## 2.0.0 development step: one rules file for the chosen assistant

A new project gets only the rules file of the assistant in `project.ai_assistant`, with the full text from one template: codex `AGENTS.md`, claude `CLAUDE.md`, gemini `GEMINI.md`; null keeps writing `AGENTS.md`. Before, `AGENTS.md` was always written and Claude/Gemini got a second one-line file importing it. Existing files are never overwritten.

## 2.0.0 development step: two checkerboards sized from the data

`checks.checkerboards: auto` and `checks.vertical: auto` now give exactly one horizontal and one depth-varying checkerboard, each still a full nonlinear inversion, in place of three horizontal sizes (fractions of the aperture) and one or two depth-varying tests. Both are set by the inversion step from the main project: cell = max(median measured wavelength c T, median nearest-station distance) (Yao et al. 2006; sampling theorem), one sign reversal at the median wavelength / 3 (rule of thumb, Fang et al. 2015). On the real line of validation.md the former defaults ran 5 checkerboards and 1 control besides the main project. Further tests can still be listed in the configuration (Rawlinson & Spakman 2016 recommend several scales). The recovery figure accepts tests whose depths differ (blank cells); depth-varying maps show depths spread over the model.

## 2.0.0 development step: failed wavelet attempts keep earlier records

A non-converged or interrupted experimental wavelet run writes `wavelet_irls.failed.json` and `DSurfTomo.in.wavelet.failed.log`; a retained earlier model keeps its own `wavelet_irls.json` and log instead of being paired with the failed attempt's record.

## 2.0.0 development step: truthful IRLS stopping

Experimental IRLS rejects invalid iteration budgets and tolerances, records inner convergence, and stops with a failure reason on LSQR/LSMR condition or iteration limits. Repeated limited updates can no longer be reported as convergence. Stop codes follow the [SciPy LSMR documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.lsmr.html).

## 2.0.0 development step: complete recovery ranges

Recovery heatmaps support a single depth, check aligned increasing depths, and retain negative correlations and amplitude overshoots. Cells display node samples rather than resolved depth intervals. Section anomalies are explicitly relative to the 1-D reference, which need not be a heterogeneous initial model.

## 2.0.0 development step: station map, section, recovery and spectrum layout

Station maps bound Mercator height for narrow line arrays, soften hill-shading under station triangles, compensate GMT colour-bar font scaling and place an internal legend separating recorded and inversion station counts. Line sections add the percentage update relative to the 1-D starting model alongside absolute Vs, thin contour labelling, and move vertical exaggeration and unmasked-resolution notes out of the title into the caption. Recovery figures combine test-by-depth heatmaps with depth profiles, readable labels and explicit ideal values (100% amplitude, 1.0 correlation). Array F-J spectra use the unified energy palette with a colour bar and explicit array-spectrum note; single-pair figures note integer-cycle branches alongside the published array spectrum.

## 2.0.0 development step: diagnostic figure meaning

Initial Vs nodes are connected rather than drawn as fictitious interfaces, and the fixed bottom is explained as a boundary condition. Extended fine-grid joint kernels are explicitly illustrative, with method metadata, rather than labelled as the inversion Jacobian. Residuals use a full-tail millisecond histogram and separate mean/std/RMS and kept/all counts; period and actual iteration panels are distinguished. Pair-figure documentation now describes cycle branches rather than a nonexistent phase-energy image. Numerical models and measurement operators are unchanged.

## 2.0.0 development step: experimental Fang 2015 D4 wavelet L1/IRLS candidate

`lant.inversion.wavelet_irls` and `tools/compare_wavelet_irls.py` implement the paper's 3-D orthogonal Daubechies D4 wavelet transform (`W W^T = W^T W = I` across dyadic, even non-dyadic and odd grid axes), augmented-residual IRLS (`R_ii = 1 / max(|r_i|, epsilon)`), inner LSQR/LSMR updates and documented L-curve/discrepancy parameter rules on top of the existing forward/Jacobian core. Controlled synthetic comparisons against the L2/LSMR baseline and pinned DSurfTomo cover blocky/sparse, smooth, noisy/outlier and uniform- vs non-uniform-depth cases with held-out traveltimes. Default solvers and configuration keys are unchanged.

## 2.0.0 development step: figures and parameters at each step

Cleaning coverage, station diagnostics and existing F-J reference figures now appear in retained results as soon as their step completes, alongside per-step effective configurations and report links with run origins. File replacement is atomic and preserves unregistered files. Solver/model publication remains tied to evaluation; work caches require explicit cleanup.

## 2.0.0 development step: assistant entry files

New output projects always receive shared `AGENTS.md` rules. `project.ai_assistant` selects an optional Gemini CLI or Claude Code import file; the default no longer creates `CLAUDE.md` unconditionally. Existing instruction files are preserved.

## 2.0.0 development step: candidate per-period propagation mesh reuse

Each Python-candidate forward model stably groups independent source fields by period, interpolates one full propagation mesh per period and retains only the current period's immutable mesh. Predictions and CSR rows remain aligned with the original observations; each source owns its writable velocity copy, traveltimes, FMM heap and refined patch. A changed model or parameters prepares new meshes; there is no approximate reuse across models. Interpolation, source seeds, rays, sensitivity cutoffs and default solver are unchanged. Source propagation remains serial and candidate scientific acceptance is still pending.

## 2.0.0 development step: inversion phase timing

`.lant/run_log.jsonl` records starts/completions for prechecks, solver batches, individual projects, automatic noise and evaluation, including failures, interruption and finished-project skips. Timings use a monotonic wall clock; parent-process CPU time is not reported as solver CPU time, and overlapping/nested phase times must not be summed as total elapsed time. Forward requests separately record cache hits and actual computation time so avoided work can be distinguished from scheduling overhead.

## 2.0.0 development step: exact forward prediction reuse

Starting-model prechecks, automatic main-model noise and final/iteration residuals share a content-addressed prediction cache. Keys cover the full model and depths, complete grid/period/physical inputs, original observation row order, solver and native runtime library hashes, frozen process source, threads and numerical environment; mtime never establishes equivalence. Only successful results that pass row-count, finite-value and row-order checks are published atomically and verified by checksum when read. Failure/interruption never publishes a result, and stochastic synthetic prechecks are not reused. Automatic noise no longer has an independent final-model mtime shortcut. Physics, rejection rules, default solver and third-party sources are unchanged. Each request records hit/miss, total, identity and actual computation timings; a hit means one complete evaluation was avoided, not a predicted whole-project speedup from a tiny case.

## 2.0.0 development step: source identity for long tasks

The first lant import freezes SHA-256 hashes of all package Python source files. Each step records that identity, process start and PID; a changed source tree fails at step/solver boundaries, and Python solver subprocesses verify the launching process fingerprint. A new on-disk Git head does not mean imported functions were updated. Use a fixed wheel or dedicated checkout for long production tasks and a separate development worktree; this check cannot make a mutable environment immutable.

## 2.0.0 development step: experimental Python implementations

Opt-in disba starting-model fitting and a shallow fundamental-Rayleigh Python tomography core now use disba, PyKonal and SciPy. PyGMT C-API workers cover maps/sections with Rasterio relief mosaics. Numerical identities, frozen synthetic observations, bounded CPU threads and normal process-group recovery are recorded. New tools compare native depth kernels, raw Jacobians, updates, common-observation models and final residuals. Defaults and QC remain unchanged; scientific replacement review is pending. See [Python implementations](python_backends.md).

## Earlier development: experimental 1-D component validation

An explicitly called disba adapter adds spherical model mapping and strict period/root checks. A comparison tool records per-period differences, failures, versions and timings against a fixed `fwd1d`, with optional pysurf96 comparison. See [forward validation](forward_validation.md). This is exploratory development: the production forward problem, fitting, solver, dependencies and scientific defaults are unchanged.

## v1.5 (2026-11): readable results, visible progress, safe interruption, installation and maintenance

No measurement, quality-control or inversion method and no scientific default changed. On the same synthetic input as v1.4.0 the stack and measurement datasets (34), the solver dispersion curves (35) and the solver inputs are identical bit for bit, and the main model is identical byte for byte.

**Changed behaviour**
- New project layout: `results/` (what you use), `work/` (processing files), `.lant/` (state and logs), with `REPORT.md` and `README.md` at the top. Projects of earlier versions are not read.
- `results/model/vs_model.csv` (longitude, latitude, depth, Vs; DSurfTomo's values with a header) replaces the solver project in the results; the per-pair dispersion text files are replaced by `curves.csv` (measured values, QC and inversion flags, the values given to DSurfTomo with their former rounding).
- `lant init` writes a compact configuration (`--full`: every key with its basis); figures are PNG by default (`plots.formats` adds PDF).
- `lant confirm` needs `--reason` with the user's decision; third-party programs are built with `lant install-external`.

**Results and report**
- The model, checks and reproduction files are replaced only after a complete evaluation, so a rebuilt or interrupted inversion keeps the previous evaluated model; `REPORT.md` opens with the state, key numbers, the items to judge and the main figures, and names the run each result comes from.

**Running**
- One progress line with the place in the chain, counts, the observed iteration of every running solver project (DSurfTomo logs unbuffered) and the elapsed time; detailed output in `.lant/lant.log` (`-v` shows it); a closing summary.
- Ctrl-C, `kill` and a closed terminal stop the solver processes and record the step as failed; `lant run all` continues and restarts interrupted solver projects from their inputs. After `kill -9` or running out of memory the next command recovers instead of refusing; a pool worker killed for lack of memory stops the step instead of leaving it waiting.
- Memory checks count job and container limits; `lant clean-work` removes rebuildable caches of a finished project on request; `lant note` works while a step runs.

**Getting started**
- `lant demo` runs a synthetic array with a known answer through the whole chain in a few minutes; `lant --version`; complete help; suggestions for misspelt configuration keys; messages that name the next command.
- `conda env create -f environment.yml` installs lant itself; installed wheels carry the solver build resources.

**AI assistants**
- Every output directory holds `AGENTS.md` with the rules (and a `CLAUDE.md` importing it); notes go through `lant note`, temporary files under `scratch/`; `lant status` lists every file lant did not write, at the top of the directory too; files lant did not write are never overwritten or deleted.

**Maintenance**
- The repository is formatted and linted with ruff; CI runs lint, the tests on Python 3.8, 3.10 and 3.12, the quick synthetic example and an installed wheel; contributing guides and an issue template.

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

Research comparisons now require a fresh output directory, hash all Python/tool sources and fixture inputs, verify integrity after computation, and retain per-case numerical arrays and inner solver diagnostics. The paired depth-grid control samples one linear-in-depth physical model on the same 0–3 km domain, uses common dense-grid observations and physical-depth quadrature. This measures combined parameterization/forward-grid effects, not an isolated wavelet effect.

The parent publishes `results/main_preview/` immediately when the main solver future completes, before other projects finish. It retains the model, solver inputs, effective parameters and a section when plotting is available; `status.json` explicitly says checks may be pending and records plotting failure without discarding the model. This is a preview, not scientific acceptance or completed evaluation.

`lant clean-work --cache forward` previews registered exact-forward `.npz` caches; explicit execution removes only unchanged registered files after evaluation, without invalidating processing stages. The default preview includes archive, cc and forward caches. Unknown/modified files, models, parameters and final figures remain protected. Large clean/cc caches are also removed only by an explicit cleanup request.

`tools/check_jacobian.py` repeats the original six-iteration weak-regularization native/Python comparison, then measures raw Jv against centered finite differences and one-sided linearization remainders for smooth, local and edge-adjacent directions at both the initial and heterogeneous final model. Perturbation scales are declared in an input protocol; raw row arrays, geometry flags and matrices are retained. It diagnoses approximation errors and does not infer an acceptance tolerance from its results.

Saved wavelet comparisons can be rendered again with `python tools/redraw_wavelet_comparison.py --experiment SAVED_EXPERIMENT --out NEW_DIRECTORY`. The numerical package must match the recorded experiment. Shared anomaly ranges include overshoots; horizontal slices label latitude and longitude coefficient indices correctly. The output records both original numerical inputs and current plotting source.

Measured fundamental-Rayleigh DSurfTomo evaluations can use the installed `fwdtimes` adapter. It links unmodified pinned routines and preserves native float32 input rounding, 8-fold propagation/source refinement, observation order and rejection weights while omitting unused depth partials, rays, Jacobian and LSMR. Builds are tied to both solver and adapter/source hashes; linked libraries enter the prediction cache identity. Older installations, other solvers and stochastic checks retain the existing evaluator. Rebuild external tools to install the adapter.

Experimental IRLS preassembles a small transformed system only when three dense buffers fit its declared 64 MiB budget; larger systems remain matrix-free. Positive regularization permits warm-starting the same unique weighted subproblem. Convergence requires both relative coefficient change and the normalized gradient of the epsilon-floor Huber envelope to meet the declared tolerance. Inner limits and outer budget exhaustion remain explicit failures to converge. Parameter-selection records now retain this quality; a small update alone is insufficient.

The Python candidate treats receivers within one minimum coarse propagation-cell distance using a straight local segment and midpoint slowness for both time and path derivative. This removes singular-gradient ray failures and near-source detours without dropping observations. Uniform-map scaling gives the independent check dT/dc=-T/c. The diagnostic `tools/verify_source_neighborhood.py` calls unmodified native rpaths and exposes its zero derivative for near receivers and omitted terminal segment; the reference itself is not altered.

Native prediction source groups run in separate processes within the existing `inversion.threads` budget, one OpenMP thread each. Workers own all Fortran module state, retain original row IDs, and merge before a single native rejection/formatting pass. Process groups are reaped on failure or interruption. A conservative allocation estimate bounds workers to 70% of effective available memory (the same engineering reserve used elsewhere); it also bounds padded receiver and propagation-grid allocations. This changes scheduling, not the numerical operator.

- Experimental wavelet projects now fail without publishing a final model when IRLS convergence is not verified; stopping diagnostics are retained.

- The controlled wavelet comparison retains holdout/discrepancy stopping histories and enforces its declared convergence requirement for every candidate and fitted model.
