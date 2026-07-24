# CONTRACT: fs-asbuilt

As-built ingestion — reality is just another chart (plan addendum,
Proposal 11): register scan data to the design and emit an honestly colored
as-built candidate.

## Purpose and layer

Layer L2 (representation/geometry). Depends on `fs-evidence` (`Color` and
`ValidityDomain`), `fs-exec` (explicit `Cx`, execution mode, and budgets),
`fs-ivl` (outward-rounded observability enclosures), `fs-la` (the
deterministic one-sided Jacobi SVD and symmetric eigendecomposition used by
the 3-D modules), and the native `fs-blake3` content-identity primitive. The
legacy scientific calculation is deterministic and uses a closed-form 2-D
rigid fit (no SVD). The additive `uncertainty` module refits that transform
under an explicit calibrated covariance model and keeps its stronger decision
semantics separate from the residual-RMS screen. The additive `rigid3` module
provides closed-form 3-D Kabsch rigid and Umeyama similarity registration
plus a calibrated 6-dof pose covariance; the additive `datum` module provides
datum-priority (3-2-1) registration with per-datum residuals and a
datum-versus-global diagnostic. The additive `propagate` module carries the
calibrated pose covariance into fs-evidence's eight-term engineering budget
as correlated cross-QoI geometry terms.

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
- `registration_invocation_resources(point_count)` returns the checked typed
  work, poll, cost, evaluation, memory, and retained-output grant for one
  registration. `register_budgeted(fiducials, &mut ChildBudget)` consumes that
  grant from an affine parent-issued child rather than reconstructing an
  allowance from the ambient `Cx`.
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
  authenticated-looking `AsBuiltDiff` value. `max_deviation_index()` retains
  the last input-order index attaining the maximum, including deterministic
  ties, so a composed workflow need not rescan the deviation payload.
- `as_built_diff_invocation_resources(...)` returns the checked typed grant and
  conservative retained-payload shape for one diff.
  `as_built_diff_budgeted(..., &Cx, &mut ChildBudget)` reserves that live-memory
  envelope before allocation, spends through the same affine child, and
  publishes retained output only on success.
- `RegError` covers cancellation with exact phase/progress,
  work-plan overflow, insufficient/collinear/unobservable/interval-unresolved/oversized points, length mismatch,
  empty data, non-finite or negative numeric inputs, malformed calibration
  identity, arithmetic overflow, bounded-allocation failure, and typed
  invocation-budget refusal.
- `uncertainty::{Covariance2, CrossFiducialModel, HuberPolicy, BiasBound,
  MetrologyModel}` declares strictly positive-definite, heteroscedastic
  per-fiducial covariance; within-pair x/y covariance; either independent or
  symmetric-principal-factor equicorrelated standardized fiducials; a finite
  radial bound on the total registered-inspection systematic vector error over
  the complete query domain (or explicit unbounded state); and a bounded
  deterministic robust policy. A raw fiducial/scanner bias is not accepted as
  that already-propagated bound. Unknown cross-fiducial dependence refuses.
- `uncertainty::estimate_calibrated_registration(fiducials, model, &Cx)`
  globally solves the fixed-weight constrained system in
  `(tx, ty, cos(theta), sin(theta))`, refitting after every Huber weight update.
  It returns `CalibratedRegistration` with the full bit-symmetric 3x3
  first-order covariance, exact `2n-3` degrees of freedom, final robust weights,
  explicit outlier dispositions/standardized residuals, full-model leverage
  diagnostics, and a domain-separated model identity. Ambiguous global
  rotations refuse. Huber covariance is a frozen-weight sandwich and is marked
  conditional; it cannot issue a finite-sample tolerance decision.
- `uncertainty::assess_calibrated_as_built(...)` propagates pose uncertainty as
  `G Cov(tx,ty,theta) G^T`, adds each independent inspection covariance exactly
  once, and applies a familywise Chebyshev-plus-union radial bound. It returns
  `DecisionState::{WithinTolerance, ExceedsTolerance, Indeterminate}` with
  lower/upper maximum-deviation bounds, confidence, family size, and a stable
  reason. Unknown registration/inspection overlap, unbounded bias, or adaptive
  weights produces an explicit bound-unavailable `Indeterminate` result.
- `uncertainty::{EvidenceReceipt, EvidenceVerifier,
  AuthenticatedAsBuiltEvidence}` separates a full content identity from
  authority. The opaque wrapper is constructible only after an injected
  verifier accepts the exact candidate/receipt under the receipt-bound policy;
  the default verifier denies everything. Authentication proves lineage, not
  physical validation or the calibration assumptions.
- `rigid3::{Point3, Fiducial3}` mirror the 2-D primitives: finite-checked,
  signed-zero-canonical, private fields with read-only accessors.
- `rigid3::register3(&[Fiducial3], &fs_exec::Cx<'_>)` solves the closed-form
  3-D Kabsch rigid fit (design → measured) through the deterministic Jacobi
  SVD of the extent-normalized weighted cross-covariance. Right-handed
  canonicalization of both singular frames makes the optimum exactly
  `V * U^T` in every admitted case, including coplanar rank-2
  cross-covariances. The result retains the rotation matrix, translation,
  advisory residual RMS, and a `RegistrationCondition` payload with the
  design/measured spectra, cross singular values, coplanarity flags, and a
  reflection-preference diagnostic for mirrored data.
- `rigid3::register3_similarity(&[Fiducial3], scale_tolerance, &Cx)` adds the
  Umeyama scale. The scale is reported, never silent: `ScaleAssessment`
  carries the estimate, a first-order standard error under an isotropic
  homoscedastic residual model, the caller-declared tolerance (no default),
  and a `UnitSuspicion` naming the nearest common unit-conversion ratio when
  the estimate leaves the declared band.
- `rigid3::{Covariance3, CrossFiducialModel3, MetrologyModel3}` declare
  strictly positive-definite per-fiducial 3x3 covariance, independence or
  explicit unknown dependence (which refuses), a bounded deterministic Huber
  policy, and a structurally valid calibration identity.
- `rigid3::estimate_calibrated_rigid3(fiducials, model, &Cx)` publishes the
  scalar-weighted closed-form Kabsch estimate (base weights
  `3 / trace(Sigma_i)`, deterministic Huber multipliers refreshed against
  declared-covariance standardized residual norms, re-solving after every
  refresh including the last) together with the full first-order sandwich
  covariance of `(tx, ty, tz, rx, ry, rz)`, exact `3n - 6` degrees of
  freedom, hat-block leverage traces summing to the parameter dimension,
  outlier dispositions, and a domain-separated model identity. The rotation
  block is a left rotation-vector perturbation about the weighted design
  centroid image. For isotropic models the sandwich reduces exactly to the
  generalized-least-squares covariance.
- `datum::DatumSystem` declares the drawing-style datum hierarchy over
  fiducial indices: A (plane, at least three targets), B (direction, at least
  two), C (one point), pairwise disjoint.
- `datum::register3_datum(&[Fiducial3], &DatumSystem, &Cx)` aligns A then B
  then C, each constraint consuming only the degrees of freedom its priority
  allows: B's out-of-plane information and C's non-axial components are
  discarded by construction. It reports signed per-datum residuals
  (out-of-plane for A, off-line for B, along-line for C), every fiducial's
  residual norm, and a `DatumGlobalComparison` carrying the rotation/
  translation delta against the embedded global Kabsch fit plus per-fiducial
  residual-norm deltas — the difference between global and datum results is
  itself a published diagnostic.
- `propagate::{QoiSensitivity, QoiEvaluator, CoveragePolicy}` declare
  first-order QoI pose sensitivities (in the calibrated covariance's exact
  parameterization and centroid pivot), an optional true-map evaluator for
  the linearization spot-check, and the caller's coverage policy (no default
  factor: the half-width multiplier is declared engineering policy, not a
  distributional claim).
- `propagate::propagate_pose_covariance(registration, sensitivities,
  coverage, evaluator, tolerance, unchecked_reason, &Cx)` emits one
  content-addressed `GeometryPropagation` record: the full cross-QoI
  covariance `G Sigma G^T`, per-QoI standard deviations, correlation
  accessors, the registration model identity it consumed, and a
  `PropagationMethod` stating exactly how the numbers were earned. With an
  evaluator, twelve deterministic one-sigma sigma points probe the true map;
  exceeding the declared tolerance downgrades the method to
  `LinearizationRejected`. Without one, the caller's non-empty reason is
  recorded as `LinearizedUnchecked`.
- `propagate::GeometryPropagation::geometry_term(ordinal)` builds the
  `Geometry` budget term for one QoI: a zero-mean `Distribution` with the
  propagated deviation and declared coverage when the linearization stands,
  and an `Unknown` with a named reason when it was rejected — a confidently
  wrong half-width is unrepresentable at this API. Every term cites the
  shared record identity as provenance and replay, which is how the cross-QoI
  correlation structure travels with per-QoI budgets.
- `field::SurfaceSample` is one SUPPLIED correspondence: a design-frame
  nominal point, the outward unit normal declared there, and the measured
  point asserted to correspond to it. A non-unit normal is refused at
  declaration, because a normal of length `l` would silently scale every
  deviation computed against it by `l`.
- `field::ProbeModel` is the instrument's declared first-order error model:
  a measured-frame unit probe direction, the smallest admitted incidence
  cosine, and one-sigma lateral and along-probe position uncertainties.
- `field::DeviationField::extract` returns the signed normal deviation at
  every admitted station with three separately reported one-sigma
  contributions — range noise `range_sigma * c`, registration
  `sqrt(g' C g)`, and correspondence ambiguity
  `lateral_sigma * sqrt(1 - c^2) / c` for `c = |n . v|` — combined in
  quadrature and scaled by the declared `CoveragePolicy`. Stations at or below
  the grazing floor are REFUSED and retained by supplied index in
  `grazing()`; `sample_count()` accounts for every supplied sample.
- `field::DeviationField::measurement_term(role)` projects the field into the
  `Measurement` slot of the eight-term engineering budget as an
  `IntervalBound` bracketing the per-point half-widths. It is deliberately NOT
  a `Geometry` term: the pose contribution reaches a budget through
  `propagate`, and emitting it twice would inflate the budget while looking
  more rigorous.
- `field::ThicknessField::extract` pairs opposing-face stations into local
  thickness `nominal + d_a + d_b`. Its registration term is the quadratic form
  of the SUMMED pose sensitivity, not a quadrature over the two faces: both
  faces ride one pose, so their pose errors are perfectly correlated and
  largely cancel on opposing normals. A station is refused if EITHER face is
  grazing.
- `field::fit_form` fits a total-degree polynomial (capped at
  `MAX_FORM_ORDER`) to the deviations over supplied in-plane stations,
  returning coefficients in an internally normalized frame plus the RMS and
  maximum residual and the fitted surface's peak-to-valley span — the warpage
  statement. A design matrix that cannot determine the declared order refuses
  with `RankDeficientFit` rather than returning an arbitrary member of the
  solution family.
- `field::profile_statistics` returns `Ra`, `Rq`, `Rt`, and segment-averaged
  `Rz` after removing a DECLARED mean-line form. These are UNFILTERED; see the
  no-claim boundaries.

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
- Each typed planner maps logical work to a distinct `WorkUnits` field and an
  equal-valued, representable `CostUnits` field, declares one scientific
  `EvaluationUnits`, and computes its poll and payload byte shapes before a
  child can spend. Registration declares no live-memory payload and only its
  fixed retained result. A diff declares the same conservative byte envelope
  for live memory and retained output. These dimensions are never converted
  into one another by the lower-layer API.
- Budgeted entry points accept only a borrowed, non-cloneable `ChildBudget`.
  Work, cost, evaluation, every poll, memory, and output are charged to that
  child; scientific refusals are latched into the parent receipt. Unused
  capacity returns only when the child is consumed by `finish`, so sibling
  stages cannot mint fresh authority through this crate.
- The `asbuilt-diff-v4` identity binds execution mode, every field of the
  ambient `fs_exec::Budget`, work-plan v2 and exact `3n` shape, poll-policy v2
  and its 256-point/256-byte strides, plus all scientific and provenance inputs.
  `StreamKey` is intentionally not part of this identity. Registration has no
  retained execution identity in this crate.
- Spatial covariance uses the rigid Jacobian ordered as `(tx, ty, theta)` and
  retains every translation/rotation cross term. Fiducial covariance factors
  are symmetric principal square roots, so the declared standardized
  equicorrelation is not an axis-order-dependent Cholesky convention. The
  equicorrelation domain is strict `-1/(n-1) < rho < 1`; boundaries are never
  clamped. Robust weighting is supported only for independent fiducials. Each
  fixed-weight transform is the global unit-circle trust-region minimum after
  eliminating translation; hard cases with multiple minima refuse. The local
  sensitivity includes the global solver's trust multiplier, and covariance is
  symmetrized once and revalidated as positive definite before publication.
- The calibrated model never converts `Registration::residual_rms` into pose or
  pointwise uncertainty. Absolute calibrated fiducial covariance determines
  parameter covariance; scaling it again by residual scatter and then adding
  residual RMS pointwise would double count the same fit error.
- For a disjoint inspection family of size `M`, confidence `1-alpha`, total
  point covariance `S_j`, and finite radial bias `b`, the simultaneous radius
  is `b + sqrt(trace(S_j) * M / alpha)`. The maximum lower bound is
  `max_j max(0, observed_j-radius_j)` and the upper bound is
  `max_j(observed_j+radius_j)`. The union bound assumes no independence among
  inspection points, but it does require calibrated covariance upper models
  and disjointness from the registration measurements. Rotation sine/cosine,
  affine mapping, pose trace, inspection trace, observed norm, radius, and
  final lower/upper arithmetic use `fs-ivl` outward enclosures so
  round-to-nearest equality cannot false-accept.
- Registration-model identity v1 binds the factor/correlation/robust/bias
  model, calibration identity, every ordered fiducial and covariance, final
  transform/covariance, standardized residual, weight, outlier disposition,
  leverage diagnostic, global-solver semantics, and degrees-of-freedom
  semantics. Spatial-evidence identity v1 additionally binds every inspection
  pair/covariance, relation, tolerance, confidence, point bound, and tri-state
  output. Both canonicalize signed zero and are tamper-evident addresses only.
- 3-D well-posedness is spectral and relative: design/measured scatter and
  the cross-covariance are classified at a stated `1e-12` relative rank gate
  after extent normalization. Coincident and collinear configurations refuse
  with a geometric diagnosis; coplanar configurations are admitted and
  flagged. A reflection-preferring cross-covariance with coincident trailing
  singular values refuses as ambiguous instead of publishing one of several
  optimal rotations; a reflection preference with a clear trailing gap is
  admitted and surfaced as a mirrored-data diagnostic on the condition
  payload.
- The scalar-weighted Kabsch fit is the exact global minimizer of its
  declared objective. The 6-dof covariance is the sandwich
  `H^{-1} (sum w^2 J^T Sigma J) H^{-1}` for that estimator under the declared
  per-fiducial covariances, symmetrized once and revalidated positive
  definite before publication; hat-block leverage traces sum to the fitted
  parameter dimension. The `rigid3` model identity binds the schema version,
  every ordered fiducial and covariance, the cross/robust model, calibration
  identity, transform, condition payload, covariance, weights, leverage, and
  outlier diagnostics with length-framed, signed-zero-canonical fields under
  a domain-separated native BLAKE3 digest.
- Datum registration never lets a lower-priority datum contradict a higher
  one: the A constraint consumes the plane orientation and normal offset at
  the fitted feature, B consumes only its in-plane projection, and C consumes
  only the along-line translation. Per-datum residuals are signed components
  in the published orthonormal constraint frame.
- The pose sensitivity of a signed normal deviation is `[n, y x n]` with
  `n = R * normal` and `y = R (point - pivot)`, in the SAME parameterization
  and about the SAME pivot as the published covariance. The pivot — the
  weighted design centroid the covariance's rotation block pivots about — is
  published by `CalibratedRigid3Registration::rotation_pivot` precisely so a
  caller never reconstructs it from the base weights and Huber multipliers and
  gets it subtly wrong; `tests/field.rs` pins the gradient against a finite
  difference of an explicitly reconstructed perturbed pose, so a wrong pivot
  fails as a first-order error.
- Every supplied field sample is accounted for exactly once: admitted stations
  and grazing refusals partition the input, refusals retain the SUPPLIED index
  (so they can be traced back to the scan), and the refusal set is bound into
  the record identity — two scans agreeing on every admitted point but
  disagreeing on what they refused are different records.
- A thickness station's pose uncertainty is the quadratic form of the summed
  face sensitivities, so a rigid motion — which cannot change the distance
  between two faces of one part — contributes no first-order thickness
  uncertainty.
- A form fit either determines its declared order or refuses. Rank is checked
  on the QR diagonal before the solve, so a rank-deficient station set never
  reaches `solve_ls` to divide by a zero pivot and return non-finite
  coefficients.

## Error model

Structured `RegError` values; hostile numeric/identity inputs return errors.
`WorkPlanOverflow` refuses an unrepresentable plan, and `Cancelled` retains the
stable phase plus exact completed/planned point visits.
`InvocationBudget` preserves the underlying typed deadline, cancellation, or
resource refusal from `fs-exec`; a scientific preflight or domain refusal is
also latched fail-closed into that invocation.
`BudgetRefused` (bead sj31i.6) retains the ambient accountant's typed refusal
verbatim: the plain `register`/`as_built_diff` entry points admit `cx.budget()`
plus the preflighted work plan through `fs_exec::AdmittedBudget` before any
work (expired deadlines - `Budget::ZERO` included - deadlines without an
ambient time source, and over-quota cost plans refuse at admission), enforce
cancellation/deadline/poll quota at every checkpoint, and charge completed
work as cost at checkpoint boundaries. Real cancellation keeps the structured
`Cancelled` shape.
Deviation allocation uses `try_reserve_exact`; no public path intentionally
panics.
`uncertainty::SpatialUncertaintyError` separately names malformed covariance,
correlation, confidence, geometry, dependence, arithmetic, allocation, and
cancellation failures. Unknown scientific dependence is never silently
converted to independence.
`field::FieldError` separately names empty/oversized sample sets, index-pairing
length mismatches, malformed scalars, a refusing registration, an all-refused
scan (`NoAdmittedSamples`, which carries the grazing count so the caller learns
WHY nothing survived), underdetermined and rank-deficient form fits, an
over-cap form order, a too-short profile, arithmetic overflow, allocation
failure, cancellation, and a refused uncertainty term. Refusals never degrade
into a silently empty field: a scan whose every station was edge-on is an
error, not an empty success.

## Determinism class

The fit, gate, δ, and calibrated spatial model are deterministic functions of
their semantic inputs.
G5 tests lock that mode, budget, work-plan, poll-policy version, and stride move
the retained diff identity without changing the numerical result.
The calibrated module uses fixed iteration counts, ordered scans, symmetric
covariance factors, canonical binary64 identity fields, and no scheduling-
dependent reduction.
The 3-D modules inherit `fs-la`'s deterministic Jacobi sweep orders and
tie-breaks and add no scheduling-dependent reduction of their own, so replay
on one platform is bitwise. No cross-ISA bit-identity is claimed for the 3-D
paths: the Jacobi kernels evaluate plain binary64 expressions whose
reconstruction accuracy (~1e-13), not bit pattern, is the portable contract.
`field` extraction is a fixed-order loop over supplied samples with no
scheduling-dependent reduction, canonical signed-zero identity fields, and
monomial basis entries built by repeated multiplication rather than `powi`, so
the fitted coefficients cannot move with the build mode. Its least-squares
paths inherit `fs-la`'s deterministic Householder QR, so the same platform
replays bitwise under the same cross-ISA boundary stated above.

## Cancellation behavior

Synchronous and cancellation-aware. Both public long-running entry points take
an explicit `Cx`; preflight precedes the initial poll, long scans poll at the
fixed 256-point stride, and a final checkpoint gates publication. Cancellation
returns `RegError::Cancelled` with exact progress and no partial output.
The budgeted forms poll the child authority, which checks its absolute clock
and originating cancellation gate before spending each poll. Typed output is
not published after a deadline, cancellation, resource, or scientific refusal.
The calibrated registration and spatial assessment also take an explicit `Cx`,
poll at bounded 256-point scan boundaries plus finalization, and publish no
partial result after cancellation. They do not yet have affine `ChildBudget`
entry points; this absence is a no-claim, not declared resource enforcement.
The `rigid3`, `datum`, and `propagate` entry points follow the same pattern:
explicit `Cx`, 256-point scan strides, a final publication checkpoint,
structured `Cancelled { phase }` refusals with no partial output, and no
affine `ChildBudget` forms yet (the same no-claim). The 3x3
SVD/eigendecomposition and 6x6 factorization calls between checkpoints are
constant-bounded work.
`field::DeviationField::extract` and `field::ThicknessField::extract` follow
that same pattern: explicit `Cx`, 256-point strides, a final publication
checkpoint, `FieldError::Cancelled { phase }` with no partial field, and no
affine `ChildBudget` form. `fit_form` and `profile_statistics` are pure
synchronous functions over already-extracted data and take no `Cx`; their cost
is bounded by the caller's station count and the capped basis width.

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
deterministic replay; typed resource planning, affine budgeted registration and
diff execution, retained last-maximum index ties, and receipt integrity; G4
pre-cancel, exact stride-boundary, mid-phase, and publication cancellation; and
G5 execution/work/poll identity separation.

`tests/field.rs`: G0 deviation recovery of injected analytic offsets under
identity and general poses (including a tangential slide that must NOT
register as deviation); the composed half-width equalling the declared
coverage-scaled quadrature exactly; opposite monotonic movement of the range
and ambiguity terms with incidence; the INDEPENDENT pivot oracle — the pose
sensitivity checked against a finite difference of the deviation under an
explicitly reconstructed perturbed pose, at a tolerance that admits only
second-order truncation and therefore rejects any wrong pivot; thickness
recovery of an injected bond line and first-order cancellation of rigid pose
error on opposing faces; adversarial near-edge-on refusal, counted by supplied
index, plus all-grazing refusal; form-fit recovery of an injected quadratic
bow, residual exposure at too low an order, and rank-deficient/underdetermined/
over-order refusals; closed-form square-wave `Ra`/`Rq`/`Rt`/`Rz`, declared
form removal cancelling a ramp, the deliberate blindness of `Ra` to a smooth
residual inside the ripple amplitude with `Rt` discriminating instead,
unfiltered statistics demonstrably absorbing waviness, and trailing-remainder
segment folding; and G5 replay equality, identity movement under a
one-nanometre content change, identity binding of the REFUSAL set, and
cancellation publishing no partial field.

`tests/spatial_uncertainty.rs`: G0 analytic independent/equicorrelated
cardinal-geometry covariance and leverage, covariance/correlation/rank refusal,
ambiguous-rotation refusal, direct-construction Huber validation,
robust-outlier disposition/downweighting with conditional no-claim,
pose-plus-inspection propagation without residual double counting, far-point
rotational leverage, outward tolerance equality at zero and nonzero rotation,
total-bias application, family-size widening, all three decision states,
overlap/bias no-claims, G3
heteroscedastic off-diagonal unit/order metamorphisms, G5 semantic identity
movement/replay, receipt mutation/policy refusal, and pre-cancel publication
refusal.

`tests/rigid3.rs`: G0 exact recovery of a general rotation and a coplanar
configuration, degenerate coincident/collinear refusals on both sides,
mirrored-data reflection preference and the symmetric ambiguous refusal, G3
rigid-conjugation invariance, similarity scale recovery with unit-suspicion
firing on a seeded 25.4 unit error and staying silent near unity, the
analytic axis-configuration 6-dof covariance, Monte-Carlo covariance
agreement on synthetic noise, Huber outlier downweighting and fit
improvement, typed model refusals, oversized-input refusal, G4 structured
cancellation, and G5 bitwise replay with input-sensitive model identity.

`tests/datum.rs`: the hand-worked block fixture recovering a seeded pose
exactly, structural hierarchy invariances (B out-of-plane and C transverse
perturbations provably cannot move the datum pose while the global fit
moves), seeded-deviation exposure with the datum-versus-global delta, typed
system/degeneracy/orientation refusals, the noisy scan-like e2e lane logging
poses, covariances, per-datum residuals, and the delta diagnostic, G5 bitwise
replay, and G4 structured cancellation.

`tests/propagate.rs`: the analytic translation marginal on the axis
configuration, exact +/-1 correlations for parallel/antiparallel gradients
with bilinear covariance scaling, spot-check acceptance on an exactly linear
map and rejection-with-downgrade on a strongly curved one, synthetic-truth
sampling recovery of the propagated variance, typed declaration/evaluator
refusals, G5 record-identity determinism and input sensitivity, G4
structured cancellation, and the register → calibrate → propagate →
eight-term-budget e2e lane logging the correlation structure.

## No-claim boundaries

- Registration requires KNOWN correspondences in both the 2-D and 3-D
  modules. Correspondence-free ICP remains an explicit [F] follow-on and is
  not smuggled in behind the Kabsch path. CT VOLUME registration (volumetric
  intensity alignment) is staged scope for the voxel layer, not this crate.
- FIELD EXTRACTION INHERITS THAT RULE. `field` consumes supplied
  correspondences and never discovers them: it does not project onto a chart,
  search for a closest point, or verify that the declared nominal point is the
  one actually measured. A wrong correspondence produces a confident wrong
  deviation, and nothing in this module can detect it. Chart-driven
  correspondence search is the follow-on that would need `fs-geom`/`fs-query`;
  this module deliberately takes no geometry-query dependency.
- The field uncertainty composition ASSUMES the station measurement errors are
  independent of the fiducial errors that fixed the pose. That is false when
  one instrument measured both, which is the normal case. The three reported
  terms are therefore a declared decomposition, not a joint covariance, and
  the composed half-width is not a confidence interval.
- The grazing floor is an ADMISSION SCREEN on measurement geometry, not a
  correctness certificate. A station above the floor can still be wrong for
  reasons this module does not model: multipath, edge effects, penetration
  into translucent material, or a correspondence error. Passing the screen
  means the declared first-order model was not obviously inapplicable.
- The ambiguity term models lateral footprint elongation at oblique incidence
  as `1/c`. It is a declared engineering model, not a derivation from an
  instrument's point-spread function, and it does not model SURFACE CURVATURE
  interacting with a lateral offset (a second-order `kappa d^2 / 2` effect
  that this module ignores because the curvature is not supplied).
- `field` handles PLATE-LIKE and simple mating geometries: the deviation is a
  scalar along a declared normal and the form fit is a low-order polynomial
  over two supplied in-plane coordinates. FREE-FORM surfaces, where no single
  in-plane parameterization exists, are explicit no-claim — the caller would
  have to supply a chart parameterization that this module does not define.
- WARPAGE is the peak-to-valley span of the FITTED low-order surface, not of
  the part. What the declared order did not capture is reported as the fit
  residual and is not folded into the warpage number; a warpage value quoted
  without its residual is a statement about the fit, not about the part.
- ROUGHNESS IS UNFILTERED AND IS THEREFORE NOT ISO 4287. `Ra`/`Rq`/`Rz`/`Rt`
  are defined there on a profile separated into roughness and waviness by a
  phase-correct Gaussian filter at a declared cutoff (ISO 16610-21), over a
  stated number of sampling lengths. NO SUCH FILTER IS IMPLEMENTED. These are
  the same arithmetic on a merely form-removed profile, so waviness inside the
  supplied trace is counted as roughness and the values generally read HIGH
  against a filtering instrument (`tests/field.rs` demonstrates that inflation
  rather than merely asserting it). They are a self-consistent relative
  statistic and a legitimate input declaration; they are not a conformance
  value and must not be reported as `Ra` without the unfiltered qualifier.
  `Rz` here averages equal INDEX segments of the supplied trace, which is not
  the ISO sampling length unless the caller made it so.
- No MATERIAL-STATE field (density, porosity) is extracted. Those require
  CT-class volumetric data, which no admitted representation in this crate
  carries; the CT lane stays staged scope, and this module reserves no schema
  for it.
- The extracted field is NOT bound into physics. A spatially varying interface
  resistance in `fs-conduction` remains one face-constant `R''` per named
  surface; consuming a thickness or roughness MAP there is a separate L3
  change that this L2 module enables but does not perform.
- Every field record identity is a domain-separated integrity address over the
  typed content and the refusal set. It is not authentication: it proves
  neither that the scan is genuine nor that the cited registration is the one
  that actually produced these coordinates.
- The 3-D calibrated estimator is the scalar-weighted Kabsch fit, not the
  generalized-least-squares optimum under anisotropic per-fiducial
  covariances; the sandwich covariance is correct for the estimator actually
  published, and the efficiency gap under anisotropy is accepted, not hidden.
  The standardized equicorrelation shortcut of the 2-D module is not offered
  in 3-D; unknown dependence refuses.
- The similarity scale standard error is a first-order diagnostic under an
  isotropic homoscedastic residual model estimated from the fit itself; it is
  not a calibrated bound, and no calibrated covariance is offered for the
  7-parameter similarity pose.
- The datum path publishes no retained execution identity and no pose
  covariance; a calibrated datum-pose uncertainty (fitted features are not
  exact constraints) is future work, and the 2-D module's simultaneous
  decision machinery is not extended over the 3-D or datum poses here.
- Geometry propagation is FIRST-ORDER only: when the sigma-point spot-check
  rejects the linearization, this crate downgrades to an explicit unknown
  instead of offering a sampling-based propagation — that machinery belongs
  to the product UQ lane. The propagation record is a tamper-evident
  identity, not authenticated authority; the scenario-side as-built
  placement binding and the nominal-versus-as-built report rendering are the
  product layers' integrations and are NOT provided here.
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
- The calibrated module provides evidence-bearing tri-state bounds, but the
  legacy boolean API remains for compatibility until downstream consumers
  migrate. Those booleans are not projections of the calibrated bounds and
  must not be promoted.
- Spatial evidence remains first-order and conditional on the supplied
  calibrated covariance/correlation and a bound on total systematic error over
  the queried domain. Raw sensor/fiducial bias is not automatically a spatial
  registration-bias bound. Huber sandwich covariance does not cover
  data-dependent weight selection, so its decision is deliberately
  unavailable. No Gaussian, exact nonlinear confidence,
  unknown-dependence, or high-leverage asymptotic claim is made.
- `EvidenceVerifier` authenticates retained lineage/policy binding only. It
  does not independently prove calibration artifact contents, the declared
  noise law, physical validation, or coverage. A lying injected verifier is an
  explicit composition-root trust failure; `NoEvidenceVerifier` admits
  nothing.
- Registration/inspection sample reuse needs retained cross-covariance and
  influence terms that v1 does not accept. Unknown or overlapping input is
  `Indeterminate` with no numeric bound rather than a zero-correlation guess.
- Point-visit work is a deterministic logical accounting unit, not an
  instruction count or a guarantee about wall-clock latency, memory pressure,
  deadline enforcement, drain behavior, or a 200-microsecond cancellation
  bound. Registration also makes no retained provenance claim about the `Cx`
  under which it ran.
- Typed planner byte counts are conservative semantic payload envelopes, not
  allocator-overhead or process-RSS measurements. `CostUnits` is abstract and
  is not a wall-time, currency, or energy certificate. A planner describes a
  grant but does not itself admit an invocation; the parent `fs-exec` issuer
  owns admission, the absolute deadline, and the terminal receipt.
- The retained diff identity is a replay/integrity binding, not authenticated
  provenance. In addition to `StreamKey`, it excludes arena identity,
  cancel-gate state, scheduler state, and other internal `Cx` state.
