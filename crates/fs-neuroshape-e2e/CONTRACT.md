# CONTRACT: fs-neuroshape-e2e

NeuroShapeCert — certified facts about a neural implicit shape. Layer L5
(LUMEN).

## Purpose and layer

Composes `fs-rep-neural` (Lipschitz + IBP), `fs-viz` (isocontour + Hessian
classification),
`fs-evidence` (Verified). Deps point downward.

## Public types and semantics

- `blob_sdf_net() -> MlpSdf` — the spectral-normalized `tanh`-MLP whose
  effective field is approximately `2.12·Σ tanh(3(±coord−0.7)) + 6.5`.
- `run_campaign(&MlpSdf, ring_r, inner) -> NeuroShapeReport` — certifies the
  Lipschitz bound, a no-tunnel sphere-trace radius, an interval topology
  certificate (inside box + a closed boundary frame), an origin-Hessian
  curvature cross-check, and localizes the sampled zero set.
- `CertifiedEnclosedComponentExists` — a private-field witness constructed only
  when the central box is interval-negative, all four strips of the closed frame
  are interval-positive, the intervals are finite and ordered, and the central
  box lies strictly inside the frame. It proves one negative component exists
  and is enclosed by the frame.
- `ComponentCountEvidence` — non-exhaustive typed state: `Unknown` has lower
  bound zero; `LowerBound(CertifiedEnclosedComponentExists)` has lower bound one.
  `exact_count()` is always `None` in this tranche.

## Invariants

- `safe_radius = |f|/L` under-estimates the distance to the NEAREST surface point
  (sound sphere tracing — no tunneling).
- TOPOLOGY: a certified-inside central box (`hi < 0`) enclosed by FOUR edge
  strips (`lo > 0`) that tile the box boundary into a CLOSED frame proves that
  the connected component meeting the central box exists and cannot cross the
  frame. `MlpSdf` is continuous (affine maps composed with `tanh`), so the
  connected negative central square lies in one negative component and every
  path from it to the exterior crosses the positive frame. Therefore the global
  component count is at least one. This does not bound the whole negative set
  and does not exclude disconnected components either inside or outside the
  frame.
- `component_enclosure_color` is `Verified` only when the typed witness exists;
  the color describes the enclosure candidate and never an exact component
  count. A too-small/open/invalid frame yields `Estimated` and typed `Unknown`.
- `boundary_frame_certified` says exactly that all four frame strips are
  certified positive. It replaces the ambiguous former field `bounded`.
- A positive-definite finite-difference Hessian at the origin is curvature
  evidence only. Without a certified zero gradient it does not establish a
  critical point or minimum, much less uniqueness or a component count.
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

`tests/neuroshape.rs` (3): G0 typed lower-bound state and private witness payload
for the certified frame, including explicit refusal to return an exact count;
Lipschitz/safe-radius/enclosure checks; an open frame yields typed `Unknown`; G5
determinism includes the typed topology evidence.

## No-claim boundaries

2-D demo net; the Lipschitz bound is the (loose) product-of-spectral-norms. The
interval-frame certificate proves that at least one enclosed negative component
exists. It does NOT prove the full negative set is bounded, exclude exterior or
additional interior components, establish a finite upper component-count bound,
or certify any exact component count. The sampled contour is localization only;
the finite-difference Hessian is not a critical-point or global Morse/Conley
certificate. There is no
complete admitted domain cover, exterior sign certificate, unresolved-cell
accounting, cubical homology witness, refinement-stability witness, sheaf-glued
coverage proof, cancellation protocol, durable replay receipt, or source-bound
exact-topology identity in this tranche. The constructor-sealed witness itself
has no source/field identity, units, budget, schema, or authenticated issuer and
is therefore campaign-local candidate data, not a portable authority receipt.
Those are required before an
`ExactComponentCount` state may exist.
