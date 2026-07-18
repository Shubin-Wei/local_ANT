# 06b Areal-Array Inversion (DSurfTomo)

## When to use

Use this branch when stations cover a 2-D region (areal array) so that ray paths cross inside it. **DSurfTomo** (Fang et al., 2015) inverts surface-wave dispersion data directly for a 3-D shear-wave velocity model, without constructing intermediate phase/group velocity maps.

> DSurfTomo belongs to this branch only; for linear arrays use the [06a linear inversion](line_inversion_en.md).

## Prerequisites

1. **Compile DSurfTomo** (see [INSTALL_en.md](INSTALL_en.md) §4):

   ```bash
   cd DSurfTomo && python configure && cd ..
   # configure enters src/, runs make, and copies the executable to bin/DSurfTomo
   ```

2. **GMT 6+** (plotting, INSTALL §5);
3. The step-05 `qc_result/` dispersion results.

Set `dsurftomo_bin: DSurfTomo/bin/DSurfTomo` in `configs/00_project.yaml` (or `paths.dsurftomo_bin` in the stage config, which takes precedence).

## Numbered entries (daily use)

One config per inversion project:

```bash
cp configs/06b_area_inversion/dsurftomo.example.yaml configs/06b_area_inversion/my_area.yaml
# edit my_area.yaml: input_dir, project_dir, period_step, grid, depths, ...
conda activate noisepy

# 1. Prepare the project: qc_result -> surfdata.dat + MOD + DSurfTomo.in (into project_dir)
python scripts/06b_area_inversion/01_prepare.py --config configs/06b_area_inversion/my_area.yaml

# 2. Run the inversion (executes inside project_dir; console output also saved to dsurftomo_run.log)
python scripts/06b_area_inversion/02_invert.py --config configs/06b_area_inversion/my_area.yaml

# 3. GMT plotting (depth slices by default; figures go to project_dir/figs)
python scripts/06b_area_inversion/03_plot.py --config configs/06b_area_inversion/my_area.yaml
```

Keys in the `prepare:` section mirror the CLI of `DSurfTomo/scripts/PrepareDSurfTomoProject.py` (translated automatically). Highlights:

| Key | Meaning |
|-----|---------|
| `curve_type` | `both` / `phase` / `group` |
| `period_step` | period grid step (s); it need not equal the step-05 DeltaT exactly — the generated period points only have to exist in the step-05 output curves (DeltaT or an integer multiple is the safe choice) |
| `grid_spacing` or `grid_count` | horizontal grid, exactly one of: `[dvxd,dvzd]` (degrees) or `[nx,ny]` |
| `depths` | depth nodes (km, comma-separated string) |
| `minvel` / `velgrad` / `t1` / `t2` / `model_vel` | initial-model estimation (null = auto) |
| `weight` / `damp` / `sublayers` / `nsrcmax` | smoothing, damping, etc. |

## Internal scripts (for debugging)

The original scripts under `DSurfTomo/scripts/` are the underlying implementation called by the numbered 06b entries; all of them are preserved. Use the numbered entries day-to-day, and run these manually only for step-by-step debugging:

```bash
cd DSurfTomo/scripts

# All-in-one (interactive prompts; recommended for manual runs)
python PrepareDSurfTomoProject.py

# Or step by step:
python extractnew.py                    # qc_result -> surfdata.dat (+ surfdata.json)
python GenerateIniMOD.py                # -> MOD initial model
python GenerateDSurfTomoInputFile.py    # -> DSurfTomo.in
python GenerateTrueMOD.py               # optional -> MOD.true (checkerboard test)

# Run the inversion (inside the project directory)
cd /path/to/project && /path/to/local_ANT/DSurfTomo/bin/DSurfTomo
```

`extractnew.py` reads `CDisp.*.txt` / `GDisp.*.txt` from `qc_result/`; station coordinates come from the file headers, no extra station file. Use the step-05 `DeltaT` (or an integer multiple) as the period step.

## Output

Produced in the project directory (verified against a real run):

- `DSurfTomo.inMeasure.dat` — final model (four columns: lon lat depth velocity; the plotting script reads this by default);
- `DSurfTomo.inMeasure.dat.iterNNN` — model after each iteration;
- `residualFirst.dat` / `residualLast.dat` — first/last iteration residuals;
- `DSurfTomo.in.log` — program log; `dsurftomo_run.log` — full console output saved by the numbered entry.

## Plotting

All plotting is driven by `03_plot.py`: list the scripts to run under `plot.scripts` in the config (executed in order; `--project-dir`/`--output-dir` are passed automatically; scripts needing extra arguments use the `name:` + `args:` form):

```yaml
plot:
  scripts:
    - plotdepth_subplot.gmt              # depth slices (gets --depths automatically)
    - name: plotcross_subplot.gmt        # cross-sections
      args: ["--profiles", "lon1/lat1/lon2/lat2/label"]
    - plotray_subplot.gmt                # ray paths
    - plotchecker_subplot.gmt            # checkerboard test
```

All four GMT scripts share the same `--project-dir`/`--output-dir` interface; see each script's `--help` for further options. Dispersion-data statistics can be plotted separately with `python DSurfTomo/scripts/plot_dispersion_period_counts.py --surfdata <project_dir>/surfdata.dat`.

## Notes

- Keep the project directory (`project_dir`) outside the repository; never commit inversion results to Git;
- `01_prepare.py` refuses to run into a non-empty project directory unless `--overwrite` is given;
- Checkerboard tests: generate `MOD.true` with the internal script `GenerateTrueMOD.py` and set the forward-modeling switches per the DSurfTomo documentation.
