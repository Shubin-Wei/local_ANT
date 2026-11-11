# Productization contract

[中文](productization.zh.md) | English

Review baseline: `23b2361` (v1.4.0). Changes are separate, reviewable commits on the existing branch; no publication is implied.

Implementation and validation are recorded in the [review notes](review_notes.md), including checks not performed.

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

Cleanup previews by default, holds the project lock when executing, removes only registered cache files of the names the steps write (`*.mseed`, `*.done.json` in the archive, `day_*.h5`), rejects unsafe paths and preserves deliverables and unknown files. Published results remain readable after cache cleanup. No legacy layout reader or migration command is introduced.

## Acceptance

Use the same raw synthetic inputs for old/new comparison. Compare measurement grids, accepted masks, rejection codes and solver input semantics. Random solver checks use scientific diagnostics, not bitwise comparisons or seed patches. Test relocation, missing outputs, locks, cleanup, interrupt recovery, JSON isolation and installation outside the source tree. Record actual validation and limitations with the final review handoff.

## Progress and logs

`lant run ... -v` shows detailed output. Normally, step stdout goes to the managed `lant.log` (10 MiB per file, three backups); stderr carries progress and the final summary. `LANT_PROGRESS=off` disables progress. Long stages without a reliable total emit a heartbeat; ETA is approximate, shown only after counted work has started, and not shown for solver projects, whose durations differ widely.

## Notes, confirmation and cleanup

`lant note "text" -c CONFIG` appends a timestamped note. `lant confirm STEP -c CONFIG --reason "decision and evidence"` records the result basis, reason and invocation mode under the same project lock. Status reports unregistered, missing and modified managed files without attributing an author; the audit compares size and modification time (kept by `cp -a` and `rsync -a`), not content. While a writer holds the lock, status shows its step as `running` and defers the file audit. Every write command (steps, note, confirm, rebuild, cleanup) takes the lock through one context that registers its outputs; a command that ended without cleanup (SIGKILL, out of memory) leaves `.lant/pending.json`, and the next write command registers its outputs and records its step as failed instead of refusing to run. SIGTERM and SIGHUP are handled like Ctrl-C (solver groups stopped, step recorded as failed); a signal ignored by `nohup` stays ignored.

`lant clean-work -c CONFIG [--cache archive|cc]` previews registered unchanged files. Only `--execute --yes` deletes them, after current evaluation has finished and raw data remain available. Final products and unknown files stay intact. Rebuilding cc or stack after cleanup requires recreating upstream caches; status exposes `cache_cleaned`. `--from` and `--rebuild` display their impact and require explicit acknowledgement; rebuild backups are unique and never replace earlier backups.

## Resource preflight

Clean estimates the cc memory and storage and warns (clean itself needs neither, and its report is the evidence for the velocity range); cc checks again before allocating workers and refuses what does not fit. Effective available memory includes readable Linux cgroup limits, with the inactive page cache counted as available as MemAvailable does. The existing 70% memory reserve is retained; only worker count changes. Storage estimates are conservative full-output bounds including replacement space; quotas and future load cannot be guaranteed. A broad line-array velocity prior is flagged without adding a confirmation gate or changing the prior.

## Dispersion table

`work/05_dispersion/curves.csv` is the solver-facing table, retained under `results/dispersion/`. It stores measured values, QC acceptance and inversion selection separately. The inversion columns deliberately preserve the former text export's rounding (period: six decimals, velocity: five, SNR: three, uncertainty: five, coordinates: eight). `disp.h5` retains rejection codes and full diagnostics. One-segment selection and missing-value semantics are unchanged.

## Installation resources

`lant install-external [dsurftomo|pdsurftomo]` uses resources bundled in the wheel and source distribution. Installed packages default to `$XDG_CACHE_HOME/lant/external` (or `~/.cache/lant/external`); a source checkout uses its `external/` directory. `LANT_EXTERNAL` overrides both. The source wrapper remains a convenience entry point. Compiler inputs and the existing memory-only patch moved into `lant/_external/` without scientific changes.

Solver progress includes the last observed iteration from the bounded DSurfTomo log tail when available. Buffered logs fall back to elapsed-time heartbeats; observed iterations are not a guarantee of percentage or remaining time. Figure tasks report their completed count.
