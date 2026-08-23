# CONTRACT: fs-neuroshape-e2e

NeuroShapeCert — certified facts about a neural implicit shape. Layer L5
(LUMEN).

## Purpose and layer

Composes `fs-rep-neural` (Lipschitz + IBP + interval sign-margin safe step) and
`fs-viz` (isocontour + Hessian classification). Deps point downward. The campaign
no longer depends on `fs-evidence`: the only enclosed-component state it
publishes is the typed `ComponentCountEvidence`, so there is no color left to
mint here.

## Public types and semantics

- `blob_sdf_net() -> MlpSdf` — the spectral-normalized `tanh`-MLP whose
  effective field is approximately `2.12·Σ tanh(3(±coord−0.7)) + 6.5`.
- `run_campaign(&MlpSdf, ring_r, inner) -> NeuroShapeReport` — certifies the
  Lipschitz bound, a no-tunnel sphere-trace radius, an interval topology
  certificate (inside box + a closed boundary frame), an origin-Hessian
  curvature cross-check, and localizes the sampled zero set. Panics on
  inadmissible input.
- `try_run_campaign(...) -> Result<NeuroShapeReport, CampaignError>` — the same
  campaign with a structured, non-trapping refusal for untrusted boundaries:
  a non-2-D network, or a non-finite / non-positive `ring_r` / negative `inner`.
- `NeuroShapeReport::safe_step: SafeStepDerivation` — the replayable no-tunnel
  step derivation from `fs_rep_neural::derive_safe_step`, carrying the origin
  enclosure, the certified `magnitude_lower_bound` on `|f(0)|`, the Lipschitz
  upper bound, the downward-rounded `radius`, and a `SafeStepStatus`.
  `origin_value` remains alongside it as a nominal visualization value with no
  certificate authority.
- `NeuroShapeReport::field_identity` plus the activation/safe-step semantics
  versions, names, and ULP budget — the arithmetic a replay must reproduce.
- `CertifiedEnclosedComponentExists` — a private-field witness constructed only
  when the central box is interval-negative, all four strips of the closed frame
  are interval-positive, the intervals are finite and ordered, and the central
  box lies strictly inside the frame. It proves one negative component exists
  and is enclosed by the frame.
- `ComponentCountEvidence` — non-exhaustive typed state: `Unknown` has lower
  bound zero; `LowerBound(CertifiedEnclosedComponentExists)` has lower bound one.
  `exact_count()` is always `None` in this tranche.
- `NEUROSHAPE_LOCALIZATION_SCHEMA_VERSION = 1` — freezes every wire code of
  the typed zero-set localization vocabulary below
  (`SurfaceLocalizationStatus` codes `1..=8`, `LocalizationStage` codes
  `1..=2`, `LocalizationDiagnostic` codes `1..=18`). Version-aware consumers
  must reject codes they do not implement; no display-string is ever parsed.
- `NeuroShapeReport::surface_localization: SurfaceLocalization` — the
  AUTHORITATIVE outcome of sampled zero-set localization:
  `Localized { crossings, max_radius, nearest_radius }`, `ValidEmpty`, or a
  refusal (`InvalidInput`, `Unrepresentable`, `ResourceRefused`,
  `Cancelled`, `AllocationRefused`, `InternalFault`) carrying its exact
  producing stage plus bounded structured detail (offender indices/packed
  edge endpoints, exact scalar bits, required-vs-limit, cancellation kind,
  stable checkpoint phase). Every `fs_viz::Grid2Error` and
  `fs_viz::IsoContourError` variant maps to exactly one documented outcome;
  nothing collapses to an undifferentiated failure.
- `NeuroShapeReport::surface_crossings` / `max_crossing_radius` /
  `nearest_surface_radius` — DERIVED, non-authoritative compatibility views
  of `surface_localization`. Zero crossings never distinguishes valid-empty
  from a refusal; `NaN` sentinels alone never carry status.
- `iso_contour_resource_code(IsoContourResource) -> u32` — the stable
  resource ordinal used as auxiliary detail for plan-overflow refusals.

## Invariants

- SOUND SPHERE TRACING: `safe_step.radius()` under-estimates the distance to the
  NEAREST surface point (no tunneling). Its authority is the INTERVAL sign margin
  at the origin — `magnitude_lower_bound` is an endpoint of the degenerate IBP
  enclosure `eval_interval([0,0], [0,0])`, hence a certified lower bound on
  `|f(0)|` — divided by the certified Lipschitz upper bound `L` and rounded DOWN.
  The nominal `origin_value` is NOT the certificate: `|origin_value|/L` is an
  ordinary round-to-nearest forward pass whose own evaluation error is
  unaccounted for, and it can exceed the true `|f(0)|/L`. An enclosure that does
  not exclude zero yields `radius = 0` with a `SafeStepStatus` that says so; no
  step is ever published without a certified sign margin.
- TOPOLOGY: a certified-inside central box (`hi < 0`) enclosed by FOUR edge
  strips (`lo > 0`) that tile the box boundary into a CLOSED frame proves that
  the connected component meeting the central box exists and cannot cross the
  frame. `MlpSdf` is continuous (affine maps composed with `tanh`), so the
  connected negative central square lies in one negative component and every
  path from it to the exterior crosses the positive frame. Therefore the global
  component count is at least one. This does not bound the whole negative set
  and does not exclude disconnected components either inside or outside the
  frame.
- `component_count_evidence` is `LowerBound(witness)` only when the typed
  interval-frame witness exists, and it never carries an exact count. A
  too-small/open/invalid frame yields typed `Unknown` with lower bound zero.
- `boundary_frame_certified` says exactly that all four frame strips are
  certified positive. It replaces the ambiguous former field `bounded`.
- A positive-definite finite-difference Hessian at the origin is curvature
  evidence only. Without a certified zero gradient it does not establish a
  critical point or minimum, much less uniqueness or a component count.
- Deterministic (fixed net + grid; no RNG).

Total on the demo net; `eval_interval`/`classify_hessian` are total.
`try_run_campaign` returns a typed `CampaignError` for a wrong input dimension or
a non-finite/out-of-range geometric parameter; `run_campaign` panics on the same
inputs. Untrusted boundaries (the WASM export) must call the fallible form.
Grid-construction and isocontour failures are NEVER erased: they surface inside
the report as the typed `SurfaceLocalization` outcome described above, so a
malformed grid and an empty contour are distinguishable at every boundary.

## Determinism class

Fully deterministic (G5).

## Cancellation behavior

The campaign itself is a synchronous batch with no ambient execution context;
its compat isocontour path therefore cannot observe a mid-run cancellation.
The localization vocabulary still carries the full typed `Cancelled` state
(stable kind + `'static` checkpoint phase) so producers that run the same
extraction under a caller-owned `Cx` publish the identical record shape.

## Unsafe boundary

None; `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/neuroshape.rs` (12): G0 pins component-evidence schema version 1, typed
lower-bound state, the localization schema version and every stable status /
stage / diagnostic code, and the private witness payload for the certified
frame, including explicit refusal to return an exact count; Lipschitz /
interval sign-margin safe-step / enclosure checks; an open frame yields typed
`Unknown`; admission refuses a non-2-D net and non-finite/out-of-range
geometry; G5 determinism includes the field identity and safe-step bits. The
typed localization battery (G0/G3/G4) maps EVERY `Grid2Error` and
`IsoContourError` variant — including plan overflow with its frozen resource
ordinals and all seven fs-exec cancellation kinds — to its documented outcome,
proves a localized campaign agrees bit-for-bit with its derived legacy views,
and proves reachable live outcomes: an identically-zero field reports
`Unrepresentable`/coincident-edge instead of a silent `None`, a strictly
positive field is `ValidEmpty` (never a refusal), and NaN samples name their
first offending node.

## No-claim boundaries

2-D demo net; the Lipschitz bound is the (loose) product-of-spectral-norms. The
interval-frame certificate proves that at least one enclosed negative component
exists. It does NOT prove the full negative set is bounded, exclude exterior or
additional interior components, establish a finite upper component-count bound,
or certify any exact component count. The sampled contour is diagnostic
localization only;
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
The typed `SurfaceLocalization` record is diagnostic evidence about the sampled
visualization only: its status codes classify what the grid/extraction kernels
refused, and they carry no topology, distance, or certificate authority.
