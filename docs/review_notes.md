# Productization review notes

[中文](review_notes.zh.md) | English

Review the changes after `23b2361` on `main`. Commits separate infrastructure, exports, installation, tests and documentation. These are local changes prepared for review; no tag, release or package publication is included. The package version remains 1.4.0.

## Implemented scope

- New projects use `work/`, retained `results/` and internal `.lant/`, with README and REPORT entry points. There is no legacy-project reader or migration layer.
- Managed logs, stderr progress, observed solver iterations and completed figure counts replace lengthy default terminal reports. Unknown totals use heartbeats.
- Registered-file audits, project locks, reasoned confirmations, notes and preview-first cache cleanup make operations inspectable. Unregistered files in a step's output directory block that step from overwriting them.
- `curves.csv` replaces per-pair exports while preserving solver rounding. Compact configuration, PNG by default, bundled compiler resources and installed-package external installation improve delivery.
- Bilingual documentation, development rules, citation metadata, scoped lint/format checks and distribution checks support maintenance.

## Validation performed

| Check | Result |
|---|---|
| Full unit suite, original Python 3.8 environment | 68 passed; 76 upstream/deprecation warnings |
| Full unit suite, isolated newer dependencies on Python 3.8 | 68 passed; one warning; unit-test coverage 53% |
| Same-input comparison with `23b2361` | 34 stack/dispersion HDF5 datasets and 35 solver-ready curves exactly equal |
| Main solver inputs | `DSurfTomo.in`, `surfdata.dat` and `MOD` byte-identical |
| Quick areal synthetic case | All nine stages completed, including the main DSurfTomo solve, evaluation and figures; checkerboard, vertical and control projects disabled |
| Cache cleanup in that isolated case | 75 retained result files unchanged by SHA-256; all stages still complete; cleaned-cache flags and file audit correct |
| Synthetic line case | Through inversion-project construction; station, dispersion and start-model figures rendered |
| Installation | Wheel and source distribution built; wheel installed outside the checkout; packaged resources used to compile DSurfTomo/fwd1d; required installation checks passed |
| Static checks | Selected Ruff lint rules, formatting of new infrastructure and `git diff --check` passed |

The separate workspace review directory holds configuration files, logs and a read-only comparison script. Tests cover lock contention, unknown files, symlinks, relocation, cleanup, failed publication, failed figure rendering, CLI forwarding and solver-child interruption/progress.

## Limits and review priorities

This is engineering validation, not a new scientific validation campaign. The full line solver, random checkerboards, pDSurfTomo and production-scale resource stress were not rerun. Remote CI has not run for these local commits. Unit coverage does not include the separate synthetic runs.

Review result publication and failure paths (`delivery.py`, `plotting/figures.py`), file ownership and cleanup (`artifacts.py`), locks and state (`records.py`, `manifest.py`, `provenance.py`), CSV precision (`export.py`), process cancellation (`processes.py`, `inversion/run.py`) and installed-package resources first.

File auditing uses size and nanosecond modification/change times, not content hashes or authorship. Same-filesystem relocation was tested; copying files may change metadata and trigger audit warnings. A terminal flag and a recorded reason are not authentication. Resource estimates cannot guarantee future load or unknown storage quotas. Figure staging preserves previous files when rendering fails; a whole figure group is not a multi-file atomic transaction.

Useful review commands from the checkout:

```sh
git log --reverse --oneline 23b2361..HEAD
git diff --stat 23b2361..HEAD
git diff 23b2361..HEAD -- lant tests
python -m pytest tests -q
```
