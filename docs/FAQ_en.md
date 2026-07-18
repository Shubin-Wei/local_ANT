# Frequently Asked Questions (FAQ)

Common questions about installation, configuration, running, and result analysis will continue to be collected here.

For installation issues (Conda, TensorFlow 1.13, gfortran builds, etc.) also see the "FAQ" section at the end of [INSTALL_en.md](INSTALL_en.md).

---

## Workflow and configuration

### Do period/velocity parameters have to be consistent across steps?

Roughly, yes. `prepare.period_step` in 06b does not have to equal the step-05 `DeltaT` exactly — the period grid it generates only needs to fall on periods that exist in the step-05 output curves (using `DeltaT` or an integer multiple of it is the safe choice). DisperPicker requires the velocity step `DeltaV` to be exactly 0.005 km/s.

### How do I choose the stacking method?

- `linear`: the workhorse, most common;
- `pws`: phase-weighted stacking, higher SNR;
- `robust`: robust stacking, insensitive to outliers;
- `all`: compute all of them.

Recommended: run `all` in step 02, then compare with step 03 (`03_stack_qc.py` / `03b_stack_panels.py`).

### After step 05, should I pick 06a or 06b?

By array geometry: stations roughly along a line → **06a** (evodcinv 1-D inversion stitched into a pseudo-2D section); stations covering a 2-D region → **06b** (DSurfTomo 3-D inversion). The branches are independent; see the comparison table in [workflow_en.md](workflow_en.md).

### What input does the 06b data converter accept?

`extractnew.py` (and the `01_prepare.py` wrapper around it) accepts only the `egfpy` input type: `CDisp.*.txt` (phase) and `GDisp.*.txt` (group) from the step-05 `qc_result/`. Lines 1–2 are station lon/lat, data columns are `period velocity SNR probability`; column 4 is the DisperPicker probability and is not a keep flag. Station coordinates are embedded in the headers — no extra station file.

### The output directory exists and the program refuses to run?

That is the silent-overwrite protection: step 01 for existing output files, step 05 for a non-empty `qc_result/`, and 06b-1 for a non-empty project directory all stop with an error. Add `--overwrite` when you really mean to overwrite.

## Environment

### Step 05 says TensorFlow is not importable?

Step 05 must run inside the `egfanalysispy` environment: `conda activate egfanalysispy`. The entry script checks this before invoking and prints this hint.

### Step 06a says `No module named 'evodcinv'`?

The evodcinv environment has not been created or activated. It lives in its own Conda environment: `conda env create -f environments/evodcinv.yml && conda activate evodcinv` — see [INSTALL_en.md](INSTALL_en.md) §2e.

### How do I quickly self-check the environment?

```bash
python tools/check_environment.py
```

It checks whatever environment it runs in; FAIL for items irrelevant to the current environment is normal (e.g. TensorFlow failing inside the noisepy env).

### 06b says "DSurfTomo executable not found"?

DSurfTomo has not been compiled, or `dsurftomo_bin` in `configs/00_project.yaml` is wrong. Build instructions: [INSTALL_en.md](INSTALL_en.md) §4; the executable lands in `DSurfTomo/bin/DSurfTomo` by default.

