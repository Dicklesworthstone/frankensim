# CONTRACT: fs-heatmap-wasm

Status: **v1 — live surface.** Browser boundary for the CMA-ES explainer
site's objective-landscape heatmaps. Own nested workspace; wasm-bindgen
confined to wasm32; native builds stay dependency-clean.

## Purpose and layer

Layer L6 browser rasterizer for the explainer's objective-landscape heatmaps.

## Public types and semantics

The wasm32 exports `heatmap_render`, `heatmap_rgba_ptr`, `heatmap_rgba_len`,
and `heatmap_version` have the buffer-lifetime and envelope semantics stated
below.

## Invariants

Rasterization preserves the shared pixel mapping, normalization modes, and
truncating color-ramp arithmetic used by the site's reference template.

## Error model

Unknown fields or modes, invalid scale, size, domain, or non-finite inputs
return the listed refusal codes; the boundary does not silently clamp or trap.

## Determinism class

Pure deterministic rasterization: identical scalar inputs produce identical
buffer bytes.

## Cancellation behavior

None: rendering is synchronous and this crate publishes no cancellation API.

## Unsafe boundary

No unsafe boundary is claimed; the manifest forbids unsafe code.

## Feature flags

None are declared in the current manifest.

## Conformance tests

The native, wasm32, and site parity checks listed below cover registry,
orientation, refusal, determinism, and reference-template agreement.

## No-claim boundaries

This is a teaching/visualization surface, not a general plotting library; the
site's JavaScript template remains its behavioral reference and fallback.

## Purpose

Seven site components rasterize a 2D objective into a background canvas with
one shared template: pixel → (x, y) → v = field(x, y) → normalization →
linear color ramp → RGBA. The JS loops burn 500k–1M transcendental
evaluations on the main thread at mount and again on every landscape or
strategy switch. This kernel hosts the field registry and rasterizes the
identical template into an in-memory RGBA buffer the page blits with zero
per-pixel JS work. Measured against the site's original loops: pixel-exact
(max channel diff 0 over 19 full-resolution configs) at roughly 2–9× the
speed for polynomial fields; trig-dense fields (rastrigin, ackley) run near
parity with JS but still land off the paint path via the site's deferred
builder.

No-claims: teaching/viz surface, not a plotting library. The site's JS
template (app/lib/frankensimHeatmap.ts renderHeatmapJs) is the behavioral
reference and fallback; formulas here must stay identical to it.

## Template (frozen v1)

- x = xmin + (px / W)·(xmax − xmin); y = ymax − (py / H)·(ymax − ymin)
  (row 0 is the TOP of the domain).
- Normalization modes: `log10p1` clamp01(log10(1+v)/k) · `tanh` tanh(v/k)
  (unclamped, mirrors the JS) · `linear` clamp01(v/k) · `sqrt`
  clamp01(√v/k) · `log10eps` clamp01(log10(max(1e-4, v+1e-4))/k).
- Ramp: r = trunc(r0 + rk·n), g = trunc(g0 + gk·(1−n)),
  b = trunc(b0 + bk·(1−n)), a = 255 — truncation toward zero, matching the
  JS `|0`.

## Field registry (ids shared with app/lib/frankensimHeatmap.ts)

`rosenbrock100` · `rastrigin` · `ackley` · `cigar-y1000` · `himmelblau` ·
`step-ridge` · `rosenbrock10` · `rot-cigar80` · `cigar-x100` · `sphere` ·
`banana-canyon` · `bowl-ripple` · `box-quad-clamp` · `box-quad-reflect` ·
`box-quad-logit`

## Public surface (wasm32 exports)

| Export | Signature | Returns |
|---|---|---|
| `heatmap_render` | (field, width, height, xmin, xmax, ymin, ymax, norm_mode, norm_k, r0, rk, g0, gk, b0, bk) -> String | `{"ok":{"kernel","width","height"}}` |
| `heatmap_rgba_ptr` / `heatmap_rgba_len` | () -> u32 | RGBA buffer location in wasm memory; valid until the next render call — the page copies immediately |
| `heatmap_version` | () -> String | `"fs-heatmap-wasm 0.1.0"` |

## Refusal codes

`field-unknown` · `norm-unknown` · `norm-scale-zero` · `size-out-of-range`
(8..=4096 per axis, ≤ 4,194,304 pixels) · `domain-degenerate` ·
`input-non-finite`

Nothing is silently clamped; no traps cross the boundary.

## Determinism

Pure function of the scalar inputs: no wall-clock, no entropy. Same call ⇒
byte-identical buffer (`deterministic_replay`).

## Parity

The site repo's bun script (scratchpad `test-heatmap-parity.mjs`) diffs this
kernel against verbatim copies of all seven components' original loops at
production resolution: 19 configs, max per-channel diff 0, zero mismatched
pixels. Rerun after touching either side.

## CI gates

`cargo check --locked --target wasm32-unknown-unknown` then
`wasm-pack build --target web --release`. Native tests: `cargo test`
(6 tests: registry coverage, hand-computed template pixel, orientation,
refusals, determinism, constraint-repair image flatness).
