# Tutorial: from records to a 3-D model, one step at a time

[中文](tutorial.zh.md) | English

Run the steps one after the other and look at each step's figures and log before the next one. A step that is run again overwrites its own folder only; the steps after it then have to be run again as well.

Every step writes into its own folder of the output directory:

| File | Content |
|---|---|
| `params.yaml` | the settings the step used, and every value it derived with the rule that gave it |
| `log.txt` | everything the step printed (with times) |
| `summary.txt` | the key numbers |
| `figures/` | the figures to check |

## 0. Input and configuration

You need:

- a folder of continuous records in any format ObsPy reads (miniSEED, SAC, ...), in any layout;
- a station table with a header line, for example

  ```
  network station longitude latitude elevation
  XX S001 100.12345 30.12345 1520
  ```

  (CSV or whitespace separated; the network column may be left out; or a StationXML file). The names must be the codes in the record headers: lant matches them exactly and lists what does not match;
- the slowest and fastest Rayleigh phase velocity you expect in the area (km/s), e.g. 0.2–1.5 on soft sediments, 1.5–3.5 on hard rock.

```bash
lant init -o my.yaml --raw /data/raw --stations stations.txt --out /data/out --velocity 0.3 3.0
```

`lant init` reads the record headers and writes `my.yaml` with suggested numbers for the sampling rate, the correlation band, the lag window and the segment length. Above each suggestion it writes the rule and the numbers of your array that gave it (for example `2 x v_fast 3 km/s / longest pair 9.8 km = ...`). Open the file and read it: every setting has a comment. Change what you know better (e.g. `data.response` if your stations do not all have the same instrument). [`example.yaml`](../example.yaml) shows a complete configuration written by `lant init`.

## 1. Cleaning the records: `lant clean -c my.yaml`

What it does: indexes the records, matches them to the stations, and writes every UTC day of every station cleaned: gaps kept as gaps, median removed, instrument response removed (if `data.response` is given), low-passed and decimated to `clean.fs_hz` (methods, section 2). Nothing is thrown away for its amplitude here; that is done segment by segment in `lant cc`.

Check:

- **`log.txt`**: the station matching (stations without records, records of stations not in the table, stations with several sensors), files that are not seismic records, days that were reused.
- **`figures/stations.png`**: are all stations there, in the right place?
- **`figures/availability.png`** and `availability.csv`: hours of data per station and day. Stations that started late, stopped early or have many gaps show up here.

Results: `days/YYYY-MM-DD.h5` (one array per day: station × sample, NaN where there is no valid sample), `stations_used.csv`, `availability.csv`, `record_index.csv`. A day is cleaned again only when its raw files or the settings change.

## 2. Cross-correlation: `lant cc -c my.yaml`

What it does: reads the cleaned days, cuts them into segments, leaves out segments with gaps or large amplitudes, correlates every pair (cross-coherence), and stacks (methods, section 3). The stations in `data.exclude_stations` and the days outside `data.start_date` / `end_date` are left out without cleaning again.

Check:

- **`log.txt`**: stations, pairs and days used, days that were reused.
- **`figures/ccf_<band>Hz.png`**: correlations against distance, one figure per sub-band of the correlation band (image above, wiggles averaged in distance bins below). Surface waves are a moveout between the two dashed lines (the velocity range). Seen on both sides of zero lag when the noise comes from all directions; mostly on one side when it comes from one direction.
- **`figures/ccf_time_*.png`**: three pairs: the stack with its signal window (grey) and noise window (blue), the sub-stacks in time order, and the signal-to-noise ratio of the stack of the first n sub-stacks. The arrival should be stable and the ratio should grow with n. A sub-stack that differs points at a station problem on that day (clock, instrument) or a strong local source. Leave such a station out with `data.exclude_stations`, or such days with `data.start_date` / `end_date`, and run the step again.

Results: `stack.h5` (the stack), `days/*.h5` (sub-stacks of each day). `lant export-sac -c my.yaml -o DIR` writes the stack as SAC files for other programs.

## 3. Phase and group velocities: `lant disp -c my.yaml`

What it does: the reference curve, then the phase and group velocities of every pair with their jackknife errors, and quality control (methods, sections 4–5).

Check, in this order:

1. **`figures/fj_reference.png`**: the most important check of the whole chain. The black line (the reference) must lie on the **fundamental-mode** ridge, the slowest ridge at each frequency. If it jumps to a higher mode somewhere:
   - limit the picking to the frequencies where the fundamental mode is clear: `disp.fj_band_hz: [lo, hi]`, or
   - give a 1-D model (`disp.reference: model`, `disp.reference_file: model.txt` with columns `depth_km Vs_km/s`) or a curve picked by hand (`disp.reference: file`, columns `frequency_Hz c_km/s`),

   then run `lant disp` again. The purple diamonds are the median of the measured phase velocities: they should lie on or near the reference. The [synthetic example](../examples/synthetic/README.md) shows a first higher mode that is stronger than the fundamental over a band of frequencies, and a curve picked along the fundamental.

   On a line array the F-J spectrum is often unclear (few pairs, unevenly spread distances). Then use a 1-D model or a curve from a nearby areal array or from the literature (`disp.reference: model` or `file`).
2. **`median_vs_reference.csv`**: a deviation of 20 % or more, or one that jumps between frequencies, means wrong cycle counts (the reference is wrong); a steady offset of a few percent is a real difference from the reference.
3. **`figures/phase_curves.png`**, **`figures/group_curves.png`**: every accepted curve with the median and standard deviation; **`figures/rejections.png`**: how many pairs are accepted at each frequency, and why the others were rejected (reason names in methods, section 5).
4. **`figures/pair_*.png`**: three example pairs: the correlation with the window between the two velocities, the phase-velocity image (its crests are the possible cycle counts) and the group-velocity image. The accepted phase velocities should lie on one crest and the group velocities on the energy maximum.

Do not loosen the quality control to get more points: the numbers in `disp` have the bases listed in [parameters.md](parameters.md).

Results: `dispersion.csv` (the accepted values, one row per pair and frequency, columns `c_kms` and `u_kms`, a rejected one left empty; the input of `lant inv`), `measurements.h5` (all values and the reasons of every rejection), `reference.csv`, `fj.h5`.

## 4. 3-D inversion: `lant inv -c my.yaml`

What it does: the average curves, the depth nodes and the 1-D starting model, the grid, and the direct 3-D inversion of the phase travel times, the group travel times or both (`inv.data`; methods, section 6). The derived numbers (depth nodes, grid spacing, Vs bounds) are printed and kept in `params.yaml`; to change one, write it into `inv` and run again.

Check:

- **`figures/start_model.png`**: does the starting model fit the average curve?
- **`figures/residuals.png`** and `residuals.csv`: the residuals of the final model should be smaller than those of the starting model; the rms per period tells where the model fits worse.
- **`figures/paths.png`**: where the data are.
- **`figures/slices.png`**, **`figures/sections.png`**: the model in depth slices and along the two sections drawn on the slices. Interpret it only where the resolution test recovers the anomalies.

Results: `vs_model.csv` (lon, lat, depth, Vs, dVs in % of the starting model, DWS: 0 means no data at that node), `start_model.csv`, `average_curve.csv`, `residuals.csv`, and `dsurftomo/` (the project in DSurfTomo's format, which the original DSurfTomo also runs).

## 5. Resolution test: `lant checker -c my.yaml --size KM`

Checkerboard of cells `--size` km, ±10 % (`checker.amplitude`), with noise as in the data; `--flip-depth KM` reverses the pattern below that depth. Run several sizes, e.g. 1, 2 and 3 times the grid spacing:

```bash
lant checker -c my.yaml --size 1.5
lant checker -c my.yaml --size 3
lant checker -c my.yaml --size 3 --flip-depth 0.5
```

Check `figures/checker.png` (input and output at each depth) and `recovery.csv` (correlation and slope per depth; the slope is the fraction of the anomaly recovered).

## When something goes wrong

Each step stops with a message that says what is missing or inconsistent and what to change. The log of the step (`log.txt`) has the full history. Problems with the data (unreadable files, stations that do not match) are listed, never silently skipped.
