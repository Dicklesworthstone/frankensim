# CONTRACT: fs-cmaes-viz-wasm

Status: **v1 — live surface.** Browser boundary for CMA-ES optimization
internals visualization. Own nested workspace; wasm-bindgen confined to
wasm32; native builds stay dependency-clean.

## Purpose

Emit the FULL internal state stream of a deterministic, seeded CMA-ES run
(Hansen 2016 couplings, mirroring the explainer site's TypeScript fallback
math) so a browser can visualize covariance ellipsoids, evolution paths, and
population flow per generation — including an honest PCA marginal for
dim > 3.

v0.2.0: the active update is the canonical Hansen 2016 form — negative
log-weights on the worst-ranked offspring, scaled by the alpha bounds and
Mahalanobis-rescaled per candidate (n/‖C^{-1/2}y‖²) — followed by a spectral
floor-and-rebuild repair, so every snapshot's `eigvals` describe a genuine
positive-definite covariance. (v0.1.0's simplified mirrored-weight heuristic
could drive C indefinite and reported the raw spectrum.)

v0.2.1: the cumulative-path `h_sigma` normalizer uses the canonical
`sqrt(1 - (1 - c_s)^(2g))` recurrence. The previous linear generation
multiplier left the square-root domain after a few generations and silently
disabled rank-one covariance-path updates for the remainder of a run. The
public `sx`, `sz`, and `sf` population streams are now reordered together so
rank, elite flag, decision vector, sample vector, and fitness stay aligned.

v0.3.0: the step-size damping is the Hansen 2016 default
`1 + 2 max(0, sqrt((mu_eff - 1)/(n + 1)) - 1) + c_s`; v0.2.1 omitted the
inner `- 1` and therefore over-damped every UI-reachable configuration. Each
generation now emits the post-update covariance eigensystem alongside its
post-update mean, sigma, and paths, rather than mixing two optimizer states.
The seeded stream uses the site's exact 32-bit LCG transition and paired
Box-Muller consumption. Reflect-repaired phenotypes are ranked and displayed,
while the latent Gaussian preimages drive mean and covariance adaptation.

No-claims: teaching/viz surface. NOT fs-dfo's production `cmaes` (no BIPOP
restarts, no identity ledgers, no adversarial-refusal hardening). For
production optimization use `fs_dfo::cmaes`.

## Public surface (wasm32 exports)

| Export | Signature | Returns |
|---|---|---|
| `cmaes_viz_run` | 18 scalars (dim, x0_0..x0_5, sigma0, lambda, active, seed, generations, landscape, noise, bounds_enabled, bound_min, bound_max, f_target) -> String | JSON envelope |
| `cmaes_viz_kernel_version` | () -> String | `"fs-cmaes-viz-wasm 0.3.0"` |

## Envelope contract (frozen v1)

- success: `{"ok":{"kernel","dim","landscape","stop_reason","best_f","best_x","total_evals","generations":[...],"pca_basis","pca_center","pca_pool_eigvals"}}`
- refusal: `{"refusal":{"code","message","ranked_repairs"}}`

Each `generations[i]`:
`{"g","mean","sigma","eigvals"(asc),"eigvecs"(row-major, col j ↔ eigvals[j]),
"cond","best_f","evals","proj_mean"[3],"proj_eigvals"[3](asc),
"proj_eigvecs"[9],"sx"(λ·n, rank order),"sz"(λ·n white-noise, rank order),
"sf"(noisy),"se"(0/1 elite),"p_sigma","p_c"}`

## Refusal codes

`dim-out-of-range` (2..=6) · `x0-non-finite` · `sigma0-non-positive` ·
`lambda-out-of-range` (4..=48) · `generations-out-of-range` (1..=200) ·
`landscape-unknown` (0..=4) · `noise-invalid` · `bounds-inverted` ·
`f-target-invalid` · `eigen-decomposition-failed` · `non-finite-objective`

Nothing is silently clamped; no traps cross the boundary.

## Determinism

Single 32-bit LCG stream (1664525/1013904223 — transition and scaling shared
with the site's TS engine) feeding paired Box–Muller draws. Same inputs ⇒
byte-identical envelope, native and wasm (`bitwise_replay_same_seed` test).
`f64::total_cmp` for all orderings; no wall-clock, no entropy.

## Landscapes (minimization)

0 sphere · 1 rosenbrock · 2 cigar · 3 rastrigin · 4 elli

## Phase-space honesty

dim ≤ 3: projection is the identity frame (direct coordinates).
dim ≥ 4: basis = top-3 eigenvectors of the FINAL covariance (pooled frame);
per-generation `proj_*` are the projected 3×3 marginal covariances'
eigendecompositions — a faithful marginal, never a fake 3D ellipsoid.
`pca_pool_eigvals` carries the full spectrum for variance-explained display.

## CI gates

`cargo check --locked --target wasm32-unknown-unknown` then
`wasm-pack build --target web --release`; nested `Cargo.lock` must not be
mutated by the pack step. Native tests: `cargo test` (exact-value and invariant
tests, including canonical damping, RNG consumption, coherent snapshots, and
latent reflection adaptation).
