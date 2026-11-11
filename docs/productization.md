# Productization contract

[中文](productization.zh.md) | English

Review baseline: `23b2361` (v1.4.0). Changes are separate, reviewable commits on the existing branch; no publication is implied.

## Output ownership

| Product | Consumer | Retention |
|---|---|---|
| Clean archive, hourly correlations | cc, stack | Rebuildable cache only when explicitly selected for cleanup; raw input must remain available |
| Stack, reference, measurement HDF5 | measurement, QC, diagnostics | Retained by default |
| Exported inversion curves | starting model, inversion | Keep measured-period ordering, missing values and legacy numeric rounding |
| Solver inputs, final models, evaluation, checks | evaluation, figures, scientific review | Publish a complete retained copy; never classify the entire inversion directory as disposable |
| Manifest, configuration history, provenance, decisions | status, reproducibility | Retain; use one state authority |

New projects use `work/`, `results/`, `.lant/`, `README.md`, `REPORT.md`, a configuration and optional `notes.md`. Layout identity is versioned. This is a new version with no legacy-project compatibility or migration layer. Internal project paths are portable where feasible; unavailable external inputs must be reported rather than silently accepted for recomputation.

## Scientific and operational invariants

No measurement or solver method, scientific default, QC criterion or third-party random implementation changes in this work. PNG-only is a presentation default for new configurations; explicit formats remain supported. Compact configuration templates do not change effective defaults.

Progress goes to stderr; machine output remains valid JSON. A terminal is a display capability, not proof of human authorization. Scientific confirmations bind a result fingerprint and record reasons and invocation metadata. Assistants must obtain the user's scientific decision and never invent one.

Cleanup previews by default, holds the project lock when executing, removes only explicitly registered cache files, rejects unsafe paths and preserves deliverables and unknown files. Published results remain readable after cache cleanup. No legacy layout reader or migration command is introduced.

## Acceptance

Use the same raw synthetic inputs for old/new comparison. Compare measurement grids, accepted masks, rejection codes and solver input semantics. Random solver checks use scientific diagnostics, not bitwise comparisons or seed patches. Test relocation, missing outputs, locks, cleanup, interrupt recovery, JSON isolation and installation outside the source tree. Record actual validation and limitations with the final review handoff.

## Progress and logs

`lant run ... -v` shows detailed output. Normally, step stdout goes to the managed `lant.log` (10 MiB per file, three backups); stderr carries progress and the final summary. `LANT_PROGRESS=off` disables progress. Long stages without a reliable total emit a heartbeat; ETA is approximate and shown only after counted work has started.

## Notes, confirmation and cleanup

`lant note "text" -c CONFIG` appends a timestamped note. `lant confirm STEP -c CONFIG --reason "decision and evidence"` records the result basis, reason and invocation mode under the same project lock. Status reports unregistered, missing and modified managed files without attributing an author.

`lant clean-work -c CONFIG [--cache archive|cc]` previews registered unchanged files. Only `--execute --yes` deletes them, after current evaluation has finished and raw data remain available. Final products and unknown files stay intact. Rebuilding cc or stack after cleanup requires recreating upstream caches; status exposes `cache_cleaned`. `--from` and `--rebuild` display their impact and require explicit acknowledgement; rebuild backups are unique and never replace earlier backups.

## Resource preflight

Clean estimates memory and storage before writing the archive; cc rechecks immediately before allocating workers. Effective available memory includes readable Linux cgroup limits. The existing 70% memory reserve is retained; only worker count changes. Storage estimates are conservative full-output bounds including replacement space; quotas and future load cannot be guaranteed. A broad line-array velocity prior is flagged without adding a confirmation gate or changing the prior.

## Dispersion table

`work/05_dispersion/curves.csv` is the solver-facing table, retained under `results/dispersion/`. It stores measured values, QC acceptance and inversion selection separately. The inversion columns deliberately preserve the former text export's rounding (period: six decimals, velocity: five, SNR: three, uncertainty: five, coordinates: eight). `disp.h5` retains rejection codes and full diagnostics. One-segment selection and missing-value semantics are unchanged.
