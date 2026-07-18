# 06a Linear-Array Inversion

## When to use

Use this branch when stations lie roughly along one line (e.g. a shallow exploration profile). 2-D/3-D tomography (DSurfTomo) lacks crossing rays perpendicular to the line in that geometry, so instead we use "per-pair 1-D inversion + along-line stitching":

1. **01_invert.py**: for each station pair's dispersion curve, invert a 1-D Vs(depth) model with [evodcinv](https://github.com/keurfonluu/evodcinv) (CPSO evolutionary algorithm; forward solver disba) and assign it to the pair **midpoint**;
2. **02_stitch.py**: project all midpoints onto the array's principal direction (PCA) to get along-line distance, then bin, average, interpolate and smooth the 1-D columns into a **pseudo-2D Vs depth section**.

> This branch has no dependency on DSurfTomo. Inversion depth is roughly limited by the maximum period (rule of thumb: sensitivity depth ≈ wavelength/3 ≈ c·T/3).

## Environment

A dedicated evodcinv Conda environment (installation: [INSTALL_en.md](INSTALL_en.md) §2e):

```bash
conda env create -f environments/evodcinv.yml   # first time only
conda activate evodcinv
```

## Input

The step-05 `qc_result/` directory, used directly with no conversion:

- `CDisp.<staA>_<staB>.txt` — phase-velocity curves;
- `GDisp.<staA>_<staB>.txt` — group-velocity curves;
- lines 1–2 are the two stations' lon/lat; data columns are `period velocity SNR probability` (points with velocity ≤ 0 are invalid).

Also required: a station CSV with elevations (shapes the section top topography), same format as the step-02 `station.csv`.

## Configuration

Template: `configs/06a_line_inversion/line.example.yaml`, one per line; `01_invert.py` and `02_stitch.py` share the same file. Three sections:

**`paths:`** — `input_dir` (qc_result path), `station_csv`, `output_dir`.

**`inversion:`** — used by `01_invert.py`:

| Key | Meaning |
|-----|---------|
| `curve_type` | `phase` (CDisp) / `group` (GDisp) / `both` (joint) |
| `workers` | parallel processes (BLAS threads pinned to 1 per worker) |
| `min_points` | skip curves with fewer valid points |
| `plot` | save the per-pair inversion figure |
| `zmax` / `dz` | output depth grid (km); zmax ≈ c·Tmax/3 recommended |
| `layers` | layered search space: per layer `{thickness: [min,max], vs: [min,max]}` (km, km/s); the sum of thickness maxima should ≈ the parameterized depth |
| `optimizer` | evodcinv CPSO settings: `popsize`, `maxiter`, `maxrun` |

**`stitch:`** — used by `02_stitch.py`: `misfit_max` (drop poorly fitting models), `nbins` (distance bins), `smooth` (Gaussian sigma), `contour_step`, `cmap`/`vmin`/`vmax`, `topography` (section top follows terrain), `show_stations`, `vertical_exaggeration`, `taper_slope` (inverted-trapezoid mask), `depth_guides`, `dpi`/`save_pdf`, etc. — all annotated inline in the template.

## Running

```bash
cp configs/06a_line_inversion/line.example.yaml configs/06a_line_inversion/my_line.yaml
# edit my_line.yaml

python scripts/06a_line_inversion/01_invert.py --config configs/06a_line_inversion/my_line.yaml
#   optional: --max-files N   (only the first N pairs; smoke test)

python scripts/06a_line_inversion/02_stitch.py --config configs/06a_line_inversion/my_line.yaml
#   optional: --output-dir DIR   (override paths.output_dir)
```

Multiple lines: keep multiple YAML copies (e.g. `line1_phase.yaml`, `line1_group.yaml`) and run them separately; use separate configs for phase/group or different `curve_type` runs as well.

## Output

```
<output_dir>/
├── models/         # per-pair 1-D Vs models (assigned to midpoints)
├── figs/           # per-pair inversion figures (when plot: true)
├── logs/           # one log per parallel worker
├── summary.txt     # summary of all inversions (incl. misfit)
├── section.png     # pseudo-2D Vs depth section
├── section.pdf     # when save_pdf: true
└── section.npz     # numeric section (distance, depth, Vs grid)
```

Failed pairs are flagged in the logs and the summary; if a model was already saved but a post-inversion diagnostic fails, the model is kept and a warning is printed (the batch is not interrupted).

## Tuning tips

- Odd stripes in the section: tighten `misfit_max` first, then increase `smooth`;
- Untrustworthy deep part: reduce `zmax` (limited by the maximum period — shallower is safer);
- Too slow: lower `optimizer.maxiter`/`popsize` or raise `workers`;
- Too few curves: revisit step-05 picking thresholds, or relax `curve_type` from `both` to a single type.
