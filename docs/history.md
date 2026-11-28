# Version history

[中文](history.zh.md) | English

Earlier versions: see the git tags (`v1.5.0`, `v2.0.0` – `v2.0.2`) and their documents.

## 3.1.0 (2026-11-27)

**Data cleaning is a step of its own**, `lant clean`, before `lant cc`. It writes every UTC day of every station cleaned (gaps kept, median removed, response removed if given, low-passed and decimated) to `1_clean/days/YYYY-MM-DD.h5`, with `availability.csv` and `figures/availability.png` (hours of data per station and day). The cleaned data can be checked before correlating, and `lant cc` can be run again with other settings, stations or days without reading the raw records again. The processing itself is unchanged.

- `cc.fs_hz` is now `clean.fs_hz`; `record_index.csv`, `stations_used.csv` and the station map moved to `1_clean`.
- Result folders renumbered: `1_clean`, `2_cc`, `3_disp`, `4_inv`, `5_checker_*`.
- Suggested upper frequency (`lant init`): the spatial sampling theorem, station spacing at most half the wavelength, f = v_slow / (2 × median station spacing). The rule of 3.0.0 (mean interval of the pair distances) gave hundreds of hertz for arrays with many pairs.
- `lant cc` uses `data.exclude_stations` and `data.start_date` / `end_date` on the cleaned days; it stops if they were cleaned at another rate or with a response pre-filter narrower than the correlation band.

## 3.0.0 (2026-11-27)

lant was rebuilt around the basic processing chain of the literature, so that every step can be followed and checked.

**How to run.** Three steps run one after the other, each checked before the next: `lant cc`, `lant disp`, `lant inv`, plus `lant checker` for the resolution test. `lant init` writes the configuration with suggested values and their rules. There is no `run all`, no confirmation gates, no hidden state: every step writes `params.yaml` (settings used and values derived, with rules), `log.txt`, `summary.txt` and its figures. Old project folders are not read.

**Kept** (methods unchanged): reading and decimating the records, response removal, segments with amplitude rejection and cross-coherence, linear stack; the F-J spectrum (Wang et al. 2019); FTAN group arrival, phase-matched filter, π/4 phase, cycle count followed with the group delay; the Python port of DSurfTomo (forward, rays, kernels) with lant's regularization, rejection and solver.

**Changed**

- Damping towards the starting model is back (`inv.damping`, default 0.1), with smoothing weight 0.5, both as means over the parameters.
- Reference curve: from the F-J spectrum (largest value per frequency in `disp.fj_band_hz`, cubic fit, no extrapolation), from a 1-D model (`disp.reference: model`, new), or from a curve file.
- Cycle count: set at the lowest usable frequency where the path is at most 1/(2 × `disp.reference_error`) wavelengths (theory: half a branch spacing), replacing the empirical 1–3.5 wavelength anchor band.
- Jackknife: each sub-stack (`cc.substack_hours`, default one per day) is one group.
- Quality control reduced to SNR, wavelength, jackknife, one run per pair; the group-pick deviation, side-lobe, group/phase consistency and first SNR gate were removed (no literature basis).
- Suggested upper frequency from the F-J aliasing rule of Wang et al. (2019, eq. 14).
- Records are read directly; there is no cleaned archive. Processing days are UTC days.
- Maps with GMT (relief from GMT's remote data when available); other figures with matplotlib.
- Results of the stack can be written as SAC (`lant export-sac`).

**Removed**: detection and removal of 1-s logger pulses; station diagnosis after the stack (replaced by the sub-stack figure); line-specific F-J ridge tracking and zoned F-J; assignment of records by serial number and time, fuzzy station-name matching; confirmation gates, manifests, provenance, step reports and AI rule files written into projects; starting-model and regularization controls (rerun with other settings instead); synthetic examples, end-to-end smoke tests and the synthetic validation documents (unit tests of the formulas remain).
