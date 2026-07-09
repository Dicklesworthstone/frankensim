# CONTRACT: fs-neuroshape-e2e

NeuroShapeCert — a PROVEN neural implicit shape. Layer L5 (LUMEN).

## Purpose and layer

Composes `fs-rep-neural` (Lipschitz + IBP), `fs-viz` (isocontour + Morse),
`fs-evidence` (Verified). Deps point downward.

## Public types and semantics

- `blob_sdf_net() -> MlpSdf` — the `tanh`-MLP `f = Σ tanh(3(±coord−0.7)) + 3`.
- `run_campaign(&MlpSdf, ring_r, inner) -> NeuroShapeReport` — certifies the
  Lipschitz bound, a no-tunnel sphere-trace radius, an interval topology
  certificate (inside box + a closed boundary frame), a Morse single-minimum
  cross-check, and localizes the zero set.

## Invariants

- `safe_radius = |f|/L` under-estimates the distance to the NEAREST surface point
  (sound sphere tracing — no tunneling).
- TOPOLOGY: a certified-inside central box (`hi < 0`) enclosed by FOUR edge
  strips (`lo > 0`) that tile the box boundary into a CLOSED frame proves the
  interior is non-empty and BOUNDED → `Verified`; a too-small box (frame overlaps
  the surface) yields `Estimated`. Boundedness rests on a closed barrier, not on
  discrete spot checks that could leave gaps.
- A single Morse minimum (`classify_hessian → Minimum`) is EVIDENCE (not proof)
  that the bounded region is a single component.
- Deterministic (fixed net + grid; no RNG).

## Error model

Total on the demo net; `eval_interval`/`classify_hessian` are total.

## Determinism class

Fully deterministic (G5).

## Cancellation behavior

None (a synchronous batch).

## Unsafe boundary

None; `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/neuroshape.rs` (3): topology certified (Lipschitz, sound safe radius,
inside box + all 4 boundary strips certified, Morse minimum, surface inside the
frame); a too-small box yields no certificate; determinism.

## No-claim boundaries

2-D demo net; the Lipschitz bound is the (loose) product-of-spectral-norms; the
topology certificate rigorously proves non-emptiness + boundedness (a closed
interval barrier) but NOT single-connectedness — that is only Morse-evidenced,
not proven, and the full homeomorphism type is out of scope. The Hessian is a
finite-difference estimate.
