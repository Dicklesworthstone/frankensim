# CONTRACT: fs-asbuilt

As-built ingestion — reality is just another chart (plan addendum,
Proposal 11): register scan data to the design and emit an honestly colored
as-built candidate.

## Purpose and layer

Layer L2 (representation/geometry). Depends on `fs-evidence` (`Color` and
`ValidityDomain`), `fs-exec` (explicit `Cx`, execution mode, and budgets),
`fs-ivl` (outward-rounded observability enclosures), and the native `fs-blake3`
content-identity primitive. The scientific calculation is deterministic and
uses a closed-form 2-D rigid fit (no SVD).

## Public types and semantics

- `Point2::new(x, y) -> Result<Point2, RegError>` constructs finite points and
  canonicalizes signed zero; coordinates are private and available through
  `x()` / `y()`.
- `Fiducial::new(design, measured)` pairs already-valid typed points; fields
  are private with read-only accessors.
- `register(&[Fiducial], &fs_exec::Cx<'_>) -> Result<Registration, RegError>` —
  the rigid rotation+translation best mapping design → measured (2-D
  Umeyama/Procrustes closed form). Requires `>= MIN_FIDUCIALS` (3)
  non-collinear fiducials and a numerically observable rotation objective;
  retains `residual_rms` as a global fit diagnostic, not registration
  uncertainty. Registration fields
  are private; `Registration::new` is fallible, accessors are read-only, and
  `Registration::apply` refuses non-finite arithmetic overflow.
- `well_posed(&Registration, certified_deviation) -> bool` — the R8 gate: true
  iff the supplied deviation is finite and positive and the registration
  residual is BELOW the deviation being certified.
- `as_built_diff(&Registration, design, scanned, design_tolerance,
  measurement_noise, calibration_candidate, &fs_exec::Cx<'_>) ->
  Result<AsBuiltDiff, RegError>` — the
  per-point δ after registration; `within_tolerance`, `above_noise_floor`, and
  a `Color::Estimated` whose domain-separated BLAKE3 identity binds every
  scientific function input plus the documented execution subset below.
  `proposed_regime` carries residual/noise/tolerance bounds for later
  authenticated review. The calibration string is a bounded, structurally
  valid candidate identity, never authority. Result fields are private and
  exposed through read-only accessors, so callers cannot forge or mutate an
  authenticated-looking `AsBuiltDiff` value.
- `RegError` covers cancellation with exact phase/progress,
  work-plan overflow, insufficient/collinear/unobservable/interval-unresolved/oversized points, length mismatch,
  empty data, non-finite or negative numeric inputs, malformed calibration
  identity, arithmetic overflow, and bounded-allocation failure.

## Invariants

- Registration RECOVERS a ground-truth rigid transform (residual → 0 on clean
  fiducials) and retains the noisy-fit RMS as an advisory diagnostic.
- Well-posedness needs `>= 3` non-collinear fiducials (rank-2 design scatter).
  Centered design and measured coordinates are normalized by their finite
  point-set extents before the relative squared-scatter rank gate and rotation
  objective, so their product expressions do not underflow/overflow merely
  because a representable configuration is uniformly rescaled;
  collinear/too-few is refused. Registration and diff inputs are capped at
  `MAX_AS_BUILT_POINTS`.
- Registration separately requires measured spread and an outward-rounded
  proof that at least one component of the centered cross-covariance vector is
  nonzero. Centroids and both cross sums carry `fs-ivl::Interval` enclosures;
  both enclosures must have finite endpoints. A collapsed measured set returns
  `UnobservableRotation`; if spread exists but both components can contain
  zero, `RotationCertificationUnresolved` reports the distinct numerical
  ambiguity and refuses fail-closed. This rejects collapsed scans and
  reflection/cancellation configurations instead
  of publishing `atan2(0, 0)`'s arbitrary zero-angle convention, without an
  epsilon heuristic.
- Public point, registration, and as-built result fields cannot be forged; all numeric inputs
  are finite, residual/tolerance/noise are non-negative, and unrecoverable
  non-finite arithmetic or a non-finite final result is refused.
- Rotation-plus-translation components preserve the ordinary binary64
  evaluation whenever it is finite and use scaled three-term evaluation when
  a rotation sum overflows before a cancelling finite translation. Recovery is
  fail-closed unless an outward-rounded interval proves the original real
  affine sum remains inside the finite binary64 range. Residual RMS uses scaled
  sum-of-squares accumulation, so a finite RMS is not rejected merely because
  squaring an individual finite distance would overflow.
- R8: `well_posed` requires a finite positive supplied deviation and is false
  when the residual meets or exceeds it (signal below the noise floor).
- The default as-built δ is always `Estimated`. Its bounded identity uses
  length-framed canonical fields followed by a domain-separated native BLAKE3
  digest, preventing delimiter and prefix collisions. Numeric identity fields
  canonicalize `-0.0` to `+0.0`, matching their mathematical equality.
- A well-formed string such as `forged-calibration-claim` cannot promote a
  result: this crate has no validated-promotion API.
- Constant-time preflight declares exactly `6n` point visits for registration
  (extrema/running-mean and anchored-normalized passes for each centroid,
  followed by scatter and residual) and `3n` point visits for a diff
  (deviation, maximum, and identity). Work is checked in `u128` before the
  initial checkpoint. Point scans poll every 256 points, with additional phase
  and final-publication checkpoints; cancellation never publishes a partial
  registration or diff.
- The `asbuilt-diff-v4` identity binds execution mode, every field of the
  ambient `fs_exec::Budget`, work-plan version and exact `3n` shape, poll-policy
  version and 256-point stride, plus all scientific and provenance inputs.
  `StreamKey` is intentionally not part of this identity. Registration has no
  retained execution identity in this crate.

## Error model

Structured `RegError` values; hostile numeric/identity inputs return errors.
`WorkPlanOverflow` refuses an unrepresentable plan, and `Cancelled` retains the
stable phase plus exact completed/planned point visits.
Deviation allocation uses `try_reserve_exact`; no public path intentionally
panics.

## Determinism class

The fit, gate, and δ are deterministic functions of their semantic inputs.
G5 tests lock that mode, budget, work-plan, poll-policy version, and stride move
the retained diff identity without changing the numerical result.

## Cancellation behavior

Synchronous and cancellation-aware. Both public long-running entry points take
an explicit `Cx`; preflight precedes the initial poll, long scans poll at the
fixed 256-point stride, and a final checkpoint gates publication. Cancellation
returns `RegError::Cancelled` with exact progress and no partial output.

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/asbuilt.rs`: exact/noisy registration, fiducial well-posedness, R8,
estimated diff semantics, proposed regime, empty/length errors, NaN/infinity/
negative rejection, invalid registration, arithmetic overflow, malformed and
forged calibration identities, delimiter-collision resistance, bounded
identity, signed-zero canonicalization, scale-invariant rank admission, and
deterministic replay; G4 pre-cancel, exact stride-boundary, mid-phase, and
publication cancellation; and G5 execution/work/poll identity separation.

## No-claim boundaries

- v1 is 2-D rigid registration (rotation + translation) with KNOWN
  correspondences; 3-D (Kabsch/SVD), scale, and correspondence-free ICP are
  follow-ons.
- Registration is treated as an optimization whose global fit RMS diagnostic
  is propagated into advisory screens and the proposed regime. That residual
  is not transform covariance or a pointwise spatial uncertainty bound.
  Writing it (and the as-built δ) to the design ledger is fs-ledger's
  integration, and the fiducial/datum PRIMITIVES at design time are fs-geom's
  (this crate consumes the correspondences).
- The scan is modeled as sampled points; admitting a full CT voxel grid /
  point cloud as a representation type with restriction maps to interface trace
  spaces extends fs-rep-voxel + fs-geom's chart zoo.
- The δ reuses the deviation metric directly; the full sheaf δ / watertightness
  machinery is the geometry layer's.
- `well_posed`, `within_tolerance`, and `above_noise_floor` are advisory
  residual/dispersion screens, not pointwise uncertainty bounds, statistical
  significance tests, or tolerance certificates.
- Calibration authenticity is an explicit no-claim. A future Validated
  promotion must inject a typed verifier, verify retained calibration artifact
  bytes/content hash under a declared policy, and bind that verification
  receipt. No such API exists in the current crate.
- Point-visit work is a deterministic logical accounting unit, not an
  instruction count or a guarantee about wall-clock latency, memory pressure,
  deadline enforcement, drain behavior, or a 200-microsecond cancellation
  bound. Registration also makes no retained provenance claim about the `Cx`
  under which it ran.
- The retained diff identity is a replay/integrity binding, not authenticated
  provenance. In addition to `StreamKey`, it excludes arena identity,
  cancel-gate state, scheduler state, and other internal `Cx` state.
