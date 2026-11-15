# Method-basis update, 2026-11-13

[中文](method_basis_update.zh.md) | English

This supplements the historical 170-item local audit; it neither changes its counts nor implies that every old threshold was verified. It records the basis of this optimization phase and derives three existing calculations directly from their implementation. Numerical methods and defaults are unchanged.

## Changes and their basis

| Choice | Implementation | Basis and scope |
|---|---|---|
| Separate indexing, archive workers, report narratives and model evaluation | `record_index`, `clean_archive`, `step_descriptions`, `inversion/model_evaluation` | Engineering responsibility boundaries. Existing processing functions retain their statements; evaluation drops one redundant JSON import. Public entry points remain available. |
| Skip fully cached header scans and cache rewrites | `record_index.build_index` | Engineering optimization. An unchanged-file regression refuses any worker-pool or raw-reader call. |
| Raw-header/archive/input identity | `record_index.file_signature`, `clean_archive._task_key`, `manifest.raw_inputs` | Engineering metadata key: path, size, nanosecond mtime/ctime, device/inode. Filesystem metadata resolution limits detection; this is not a content checksum of large waveforms. Same-second changed-mtime edits are covered by regressions. |
| Response inventory/archive invalidation | `response.identity/load`, `clean.instruments` | Engineering metadata plus SHA256 contents for the smaller response inputs. A regression preserves size and mtime while changing contents; both inventory and clean-input provenance must change. |
| Main-figure process isolation and concurrency | `plotting.figure_batch`, `plotting.figures` | Engineering scheduling under the existing automatic resource plan. The user requested speed with server headroom. Private renderer state, source-identity checks, staged publication and failed-group preservation are retained; 17 PNGs match byte for byte. |
| Figure memory estimate | `plotting.figure_batch.memory_estimate` | Explicitly empirical: 600 MB interpreter/rendering allowance, eight float64 buffers for stack/hourly arrays, existing pair-image work estimate, and model-file parsing allowance. A scheduling estimate, not a guaranteed RSS bound. |
| Common-observation solver comparison | external comparison record; `tools/plot_solver_comparison.py` | Same measured-mode files and matching parsed observed times isolate inversion differences from observation generation. The checkerboard observations are serialized once from the existing native zero-noise calculation. Grid, weights and iteration counts are fixed benchmark choices; their values do not become production defaults. |
| Checkerboard and figure appearance | synthetic example and comparison figures | User decision, 2026-11-13: retain checkerboards and one journal figure style. This display preference is not a resolution criterion. |

## Weighted linear stacking

For a particular pair, hour h holds n_h accepted segment correlations and stores their average x_h. The implementation uses

$$X=\frac{\sum_h n_h x_h}{\sum_h n_h}
=\frac{\sum_h\sum_{s=1}^{n_h}x_{hs}}{\sum_h n_h}.$$

Thus count weighting represents the mean over accepted segments rather than giving a sparse hour the same weight as a full hour. This is the algebraic basis of `stack.build_stack`, not a claim that a publication prescribes every grouping interval. Float32 hourly storage and subsequent accumulation introduce the documented rounding choices. Zero-count hours contribute nothing.

## Starting-model least squares and curvature

Let m = ln(Vs / (1 km/s)), r = observed phase velocity − g(m), and J approximate the derivative of g with respect to m. `start.invert1d` solves

$$\min_{\delta}\|J\delta-r\|^2+s^2\|L(m+\delta)\|^2+d^2\|\delta\|^2,$$

using A = [J; sL; dI] and b = [r; −sLm; 0]. This follows directly by expanding the three residual blocks. It penalizes curvature of the updated model and size of the step. The relative perturbation quotient divides by ε = 0.01, a first-order approximation to a log-velocity derivative; its denominator is not ln(1+ε). The finite difference, weights, 14 iterations and clipping δ to ±0.15 are engineering choices. Exponentiating ±0.15 gives approximately +16.18%/−13.93%, rather than exactly ±15%.

For adjacent depth spacings h₁, h₂ and total depth range Z, `start.curvature_operator` uses the nonuniform second-derivative coefficients

$$\left[\frac{2}{h_1(h_1+h_2)},-\frac{2}{h_1h_2},\frac{2}{h_2(h_1+h_2)}\right].$$

Multiplication by Z² converts to derivatives in ξ = (z−z₀)/Z. Multiplication by √((h₁+h₂)/(2Z)) gives each interior node its representative integration width, so the squared norm approximates integrated squared curvature over the represented interior intervals. It is not a proof of layer resolution or exact boundary quadrature. Existing tests check zero curvature for a linear function on unequal spacings and approximate consistency for a cubic on two grids. These algebraic derivations support the implementation without claiming that unopened Shaw/Constable papers prescribe its exact discretization.

## Travel-time errors and regularization boundary

The first-order conversion σ_t ≈ d σ_c/c² (t = d/c at fixed distance) and the combined rule w = d, f = max(f_B, f_M) that used it were replaced in 2.0: the smoothing weight and the damping are fixed values calibrated on the synthetic examples (methods.md 8.4).

## What this establishes

The new comparisons address implementation agreement; the known-truth experiment addresses recovery under its stated geometry/settings. Real-data fit alone does not establish a correct Earth model. The historical audit and the explicitly unverified literature are kept with the development records outside the repository; the present document does not relabel them as source-verified.
