# CONTRACT: fs-lenia-wasm

Status: **v1 — live surface.** Browser boundary for the CMA-ES explainer
site's continuous Lenia field. Own nested workspace; wasm-bindgen confined to
wasm32; native builds stay dependency-clean.

## Purpose

Replace the site's O(N²·R²) direct-convolution TypeScript Lenia (capped at
96²) with the SAME model computed through a radix-2 FFT — O(N² log N) exact
toroidal convolution — so a 256² or 512² display field with a proportionally
larger ring kernel steps faster than the 96² fallback, phones included. Also
hosts the reduced-resolution snapshot-seeded fitness rollouts the site's
CMA-ES search evaluates.

v0.2.0 performance structure (identical math, verified by the same tests):
precomputed twiddle/bit-reversal tables per size, cache-blocked in-place
transposes over row-padded planes (power-of-two row strides alias cache sets
badly enough under wasm to triple the 512² step cost), a transpose-sparing
convolution path that multiplies in the transposed spectral layout (valid
because the distance-only kernel's spectrum is transpose-symmetric), and a
growth LUT cached across steps while (μ, σ) are unchanged. Measured in wasm
(bun/JSC, M-series): ~1.8 ms per 256² step, ~8.6 ms per 512² step.

No-claims: teaching/viz surface, not a general spectral PDE solver. The TS
fallback remains the behavioral reference at 96²; this kernel is the same
mathematics at higher resolution, verified FFT-vs-direct to < 1e-10.

## Model (frozen v1; mirrors the site's TS fallback)

- Kernel: ring K(r) = exp(−((r/R − 0.5)/0.18)²/2) for r ≤ R, unit-normalized,
  toroidal. R = rel_radius · size (fractional radii honored exactly).
- Growth: G(u) = 2·exp(−((u − μ)/σ)²/2) − 1, applied through a 2048-entry
  per-step LUT with linear interpolation.
- Update: a ← clamp(a + Δt·G(K∗a), 0, 1).
- Metrics: interface = fraction of cells strictly in (0.08, 0.92);
  mass = mean activation.
- Fitness: mean over steps of (interface − 2·|mass − 0.25|), rolled out from
  a box-averaged snapshot at eval resolution.

## Public surface (wasm32 exports)

| Export | Signature | Returns |
|---|---|---|
| `lenia_init` | (size: u32, eval_size: u32, rel_radius: f64) -> String | `{"ok":{"kernel","size","kernelRadius"}}` |
| `lenia_clear` | () | — |
| `lenia_seed_ring` | (cx, cy, radius, ring_frac, width, intensity) | — (additive hollow gaussian ring, wrapped) |
| `lenia_step` | (mu, sigma, dt, steps: u32) -> String | `{"ok":{"interface","mass"}}` (last step) |
| `lenia_render` | () | — (colormaps the field into the RGBA buffer) |
| `lenia_rgba_ptr` / `lenia_rgba_len` | () -> u32 | RGBA buffer location in wasm memory; allocated once per init, never reallocated — the page rewraps per frame in case memory grows |
| `lenia_snapshot_eval` | () | — (freeze the display field, box-averaged, as the eval seed) |
| `lenia_eval` | (mu, sigma, dt, steps: u32) -> String | `{"ok":{"score"}}` from the frozen seed |
| `lenia_version` | () -> String | `"fs-lenia-wasm 0.2.0"` |

## Refusal codes

`size-out-of-range` (power of two, 64..=512) · `eval-size-invalid` (power of
two ≥ 32, divides size) · `rel-radius-out-of-range` (0.01..=0.25) ·
`input-non-finite` · `sim-not-initialized`

Nothing is silently clamped; no traps cross the boundary.

## Determinism

Pure functions of the field state and scalar inputs: no wall-clock, no
entropy. The same call sequence replays bitwise within a build
(`deterministic_replay`, `eval_scores_are_finite_and_replayable`).

## CI gates

`cargo check --locked --target wasm32-unknown-unknown` then
`wasm-pack build --target web --release`. Native tests: `cargo test`
(8 tests, including FFT-vs-direct convolution equivalence at < 1e-10).
