# CONTRACT: fs-asbuilt

As-built ingestion — reality is just another chart (plan addendum,
Proposal 11): register scan data to the design and emit a validated-color δ.

## Purpose and layer

Layer L2 (representation/geometry). Depends only on `fs-evidence` (the `Color`
+ `ValidityDomain`). Pure, deterministic; a closed-form 2-D rigid fit (no SVD).

## Public types and semantics

- `Point2 { x, y }`; `Fiducial { design, measured }` — a design datum and its
  scanned location.
- `register(&[Fiducial]) -> Result<Registration, RegError>` — the rigid
  rotation+translation best mapping design → measured (2-D Umeyama/Procrustes
  closed form). Requires `>= MIN_FIDUCIALS` (3) non-collinear fiducials; carries
  `residual_rms` forward (the registration uncertainty). `Registration::apply`
  maps a design point into measured coordinates.
- `well_posed(&Registration, certified_deviation) -> bool` — the R8 gate: true
  iff the registration residual is BELOW the deviation being certified.
- `as_built_diff(&Registration, design, scanned, design_tolerance,
  measurement_noise, calibration_cert) -> Result<AsBuiltDiff, RegError>` — the
  per-point δ after registration; `within_tolerance`, `above_noise_floor`, and
  a `Color::Validated` anchored (dataset) to the calibration certificate, its
  regime tagged with the residual + measurement noise.
- `RegError` — `TooFewFiducials` / `CollinearFiducials` / `LengthMismatch` /
  `Empty`.

## Invariants

- Registration RECOVERS a ground-truth rigid transform (residual → 0 on clean
  fiducials) and CARRIES its error forward on noisy ones.
- Well-posedness needs `>= 3` non-collinear fiducials (rank-2 design scatter);
  collinear/too-few is refused.
- R8: `well_posed` is false when the residual meets/exceeds the certified
  deviation (signal below the noise floor) or the deviation is non-positive.
- The as-built δ is always `Validated` and anchored to the calibration cert.

## Error model

Structured `RegError` values; no panics.

## Determinism class

Fully deterministic: the fit, gate, and δ are pure functions of the inputs.

## Cancellation behavior

None (synchronous pure functions).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/asbuilt.rs` (Proposal 11, 8 cases): exact rigid-transform recovery +
apply; noisy residual carried forward; too-few + collinear fiducials refused;
the R8 signal-vs-noise gate; the validated δ anchored to the cert; a
below-noise-floor deviation flagged; malformed-input rejection; determinism.

## No-claim boundaries

- v1 is 2-D rigid registration (rotation + translation) with KNOWN
  correspondences; 3-D (Kabsch/SVD), scale, and correspondence-free ICP are
  follow-ons.
- Registration is treated as an optimization whose error is carried forward;
  writing it (and the as-built δ) to the design ledger is fs-ledger's
  integration, and the fiducial/datum PRIMITIVES at design time are fs-geom's
  (this crate consumes the correspondences).
- The scan is modeled as sampled points; admitting a full CT voxel grid /
  point cloud as a representation type with restriction maps to interface trace
  spaces extends fs-rep-voxel + fs-geom's chart zoo.
- The δ reuses the deviation metric directly; the full sheaf δ / watertightness
  machinery is the geometry layer's.
