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

New projects use `work/`, `results/`, `.lant/`, `README.md`, `REPORT.md`, a configuration and optional `notes.md`. Layout identity is versioned. Existing projects keep their original paths and lock. Internal project paths are portable where feasible; unavailable external inputs must be reported rather than silently accepted for recomputation.

## Scientific and operational invariants

No measurement or solver method, scientific default, QC criterion or third-party random implementation changes in this work. PNG-only is a presentation default for new configurations; explicit formats remain supported. Compact configuration templates do not change effective defaults.

Progress goes to stderr; machine output remains valid JSON. A terminal is a display capability, not proof of human authorization. Scientific confirmations bind a result fingerprint and record reasons and invocation metadata. Assistants must obtain the user's scientific decision and never invent one.

Cleanup previews by default, holds the project lock when executing, removes only explicitly registered cache files, rejects unsafe paths and preserves deliverables and unknown files. Published results remain readable after cache cleanup. Old layouts are not implicitly migrated or rerun.

## Acceptance

Use the same raw synthetic inputs for old/new comparison. Compare measurement grids, accepted masks, rejection codes and solver input semantics. Random solver checks use scientific diagnostics, not bitwise comparisons or seed patches. Test legacy access, relocation, missing outputs, locks, cleanup, interrupt recovery, JSON isolation and installation outside the source tree. Record actual validation and limitations with the final review handoff.
