# Frequently Asked Questions (FAQ)

Common questions about installation, configuration, running, and result analysis will continue to be collected here.

For installation issues (Conda, TensorFlow 1.13, gfortran builds, etc.) also see the "FAQ" section at the end of [INSTALL_en.md](INSTALL_en.md).

---

## Workflow and configuration

### Do period/velocity parameters have to be consistent across steps?

Yes. The period range and step set in step 05 (`dispersion:` section of `configs/05_dispersion.yaml`) must match `prepare.period_step` in 06b (06a reads the curve files directly and has no such requirement). DisperPicker requires the velocity step `DeltaV` to be exactly 0.005 km/s.

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

evodcinv is not installed. `conda activate noisepy && pip install evodcinv tqdm` — see [INSTALL_en.md](INSTALL_en.md) §2e.

### How do I quickly self-check the environment?

```bash
python tools/check_environment.py
```

It checks whatever environment it runs in; FAIL for items irrelevant to the current environment is normal (e.g. TensorFlow failing inside the noisepy env).

### 06b says "DSurfTomo executable not found"?

DSurfTomo has not been compiled, or `dsurftomo_bin` in `configs/00_project.yaml` is wrong. Build instructions: [INSTALL_en.md](INSTALL_en.md) §4; the executable lands in `DSurfTomo/bin/DSurfTomo` by default.

## Migration

### I used the old script names (local_cc_workflow.py etc.)?

Functionality and configuration are fully compatible — only names and locations changed. See the "Old script names → new entries" table in [workflow_en.md](workflow_en.md); old config file contents can be copied verbatim into the renamed counterparts.

### I used to edit EGFAnalysisPy/config.py or LineXConfig.py by hand?

Now you only maintain `configs/05_dispersion.yaml` (one copy per line/batch); the entry generates the runtime config automatically and you never enter the EGFAnalysisPy directory. The manual way still works (`run_parallel.py --config your_config.py`).
