# Version history

[中文](history.zh.md) | English

Earlier versions: see the git tags (`v1.5.0`, `v2.0.0` – `v2.0.2`) and their documents.

## 3.1.0 (2026-11-27)

**New step `lant clean`** before `lant cc`: every UTC day of every station is cleaned (gaps kept, median removed, response removed if given, low-pass and decimation) and written to `1_clean/days/YYYY-MM-DD.h5`, with `availability.csv` and `figures/availability.png` (hours of data per station and day). `lant cc` reads these days; it can be run again with other settings, stations or days without reading the raw records.

**Group velocities.** `lant disp` measures the Rayleigh group velocity of every pair (envelope maximum at the instantaneous frequency) with its own jackknife error and quality control; `dispersion.csv` has the columns `c_kms`, `sigma_kms`, `u_kms`, `u_sigma_kms`. `lant inv` inverts the phase travel times, the group travel times or both (`inv.data`, default `phase`).

**Synthetic example** (`examples/synthetic`): noise records over two 1-D models (measured against theoretical dispersion) and travel times through a 3-D model with a basin, a fault zone and a buried body (true against recovered model), with a spike test and a checkerboard.

- `cc.fs_hz` is now `clean.fs_hz`; `record_index.csv`, `stations_used.csv` and the station map are in `1_clean`.
- Result folders: `1_clean`, `2_cc`, `3_disp`, `4_inv`, `5_checker_*`.
- Suggested sampling rate and upper frequency (`lant init`): `clean.fs_hz` is the highest whole-number rate reachable from every raw rate, `cc.fmax_hz` = 0.4 × `clean.fs_hz` (Nyquist theorem). The rule of 3.0.0 (mean interval of the pair distances) gave hundreds of hertz for arrays with many pairs.
- `lant cc` applies `data.exclude_stations` and `data.start_date` / `end_date` to the cleaned days; it stops if the days were cleaned at another rate or with a response pre-filter narrower than the correlation band.
- Figures redrawn in the style of journal figures: correlations per sub-band with a lag axis fitted to the array (`ccf_<band>Hz.png`); per pair the stack, sub-stacks and signal-to-noise ratio against stacking (`ccf_time_*.png`); dispersion curves with median and standard deviation (`phase_curves.png`, `group_curves.png`); one `rejections.png` for phase and group velocities; example pairs with the correlation, phase- and group-velocity images; paths, slices, sections (`sections.png`, replacing `profiles.png`) and checkerboards with GMT. `average_curve.csv`, `start_fit.csv` and `residuals.csv` of `lant inv` have a column `velocity` (phase or group).
- `example.yaml`: a complete configuration written by `lant init`.
- Dispersion of a model column in which Vs decreases with depth: root search in steps of 0.0005 km/s instead of 0.005 km/s. With the coarser step the search could jump to a higher mode or find no root (the 3-D inversion then stopped).
- `lant disp` no longer stops on a pair whose data all come from one sub-stack.

## 3.0.0 (2026-11-27)

**How to run.** Steps run one after the other, each checked before the next: `lant cc`, `lant disp`, `lant inv`, and `lant checker` for the resolution test. `lant init` writes the configuration with suggested values and their rules. Every step writes `params.yaml` (settings used and values derived, with rules), `log.txt`, `summary.txt` and its figures. Project folders of 2.x are not read.

**Unchanged methods**: reading and decimating the records, response removal, segments with amplitude rejection and cross-coherence, linear stack; the F-J spectrum (Wang et al. 2019); FTAN group arrival, phase-matched filter, π/4 phase, cycle count followed with the group delay; the 3-D inversion (DSurfTomo's forward, rays and kernels with lant's regularization, rejection and solver).

**Changed**

- Damping towards the starting model is back (`inv.damping`, default 0.1), with smoothing weight 0.5, both as means over the parameters.
- Reference curve: from the F-J spectrum (largest value per frequency in `disp.fj_band_hz`, cubic fit, no extrapolation), from a 1-D model (`disp.reference: model`, new), or from a curve file.
- Cycle count: set at the lowest usable frequency where the path is at most 1/(2 × `disp.reference_error`) wavelengths (theory: half a branch spacing), replacing the empirical 1–3.5 wavelength anchor band.
- Jackknife: each sub-stack (`cc.substack_hours`, default one per day) is one group.
- Quality control reduced to SNR, wavelength, jackknife, one run per pair; the group-pick deviation, side-lobe, group/phase consistency and first SNR gate were removed.
- Suggested upper frequency from the F-J aliasing rule of Wang et al. (2019, eq. 14).
- Records are read directly; there is no cleaned archive. Processing days are UTC days.
- Maps with GMT (relief from GMT's remote data when available); other figures with matplotlib.
- Results of the stack can be written as SAC (`lant export-sac`).

**Removed**: detection and removal of 1-s logger pulses; station diagnosis after the stack (replaced by the sub-stack figure); line-specific F-J ridge tracking and zoned F-J; assignment of records by serial number and time, fuzzy station-name matching; confirmation gates, manifests, provenance files and step reports in project folders; starting-model and regularization controls (rerun with other settings instead); synthetic examples.
