# fs-motion CONTRACT

Certified rigid motion for the MORPH layer: a motor path over a time
domain becomes a checkable object (a *tube*) instead of a per-evaluation
pose, and moving geometry exposes time-aware queries instead of
pretending to satisfy the timeless `Chart` contract.

## Purpose and layer

Layer **L2 (MORPH)**. `fs-ga` owns instantaneous SE(3) motor algebra and
`fs-scenario` owns frame trees; nothing previously bound a motor *path*
to a chart. `fs-motion` provides that binding. Dependencies: `fs-ga`,
`fs-geom`, `fs-ivl`, `fs-exec`, `fs-evidence`, `fs-math`, `fs-query`. This crate
must NEVER depend on `fs-scenario` or any higher layer; higher layers
lower their motions into tubes through [`LowerToMotorTube`].

## Public types and semantics

- `CertifiedMotorTube` — a piecewise enclosure of a motor path
  `t ↦ M(t)` over a closed time interval. Each segment stores univariate
  Taylor models (`fs_ivl::TaylorModel1`) for the sixteen PGA multivector
  components of the motor (only the eight even slots are ever nonzero),
  together with a rigorously computed **versor-defect bound**: an upper
  bound of `‖M(t) M̃(t) − 1‖∞` over the segment domain, derived from the
  same Taylor models by rigorous polynomial arithmetic. This is the
  "component models with unit-motor/Study-constraint residual bounds"
  representation: the unit-norm and Study constraints are not assumed,
  they are *measured* and reported.
- `MotorTubeSegment` — one segment: component models, domain, defect.
- `MotorPath` — the point-evaluation view of a tube: `motor_at(t)`
  returns an `fs_ga::Motor` (midpoint of the component enclosures at
  `t`) plus honesty data (`PathSample`: max enclosure width, defect).
- `PointActionEnclosure` / `BoxActionEnclosure` — interval enclosures of
  `M(t)·x` (respectively `M(t)·B` for an AABB `B`) over a time
  subinterval, tagged with an [`EnclosureClass`] and the defect bound.
  The homogeneous weight is divided out as an interval; a weight
  enclosure containing zero refuses. Because the sandwich divides by
  the homogeneous weight, uniform versor scaling cancels: the enclosure
  covers the action of the *constructed component path* exactly.
- `EnclosureClass` — `Certified` (Taylor-model enclosure with rigorous
  remainder) or `FalsifiedOnly` (checked only by sampling). Every
  action/evaluation output states its class.
- `SpacetimeChart<C>` — a base `Chart` plus a body-to-world tube.
  `snapshot(t)` returns an immutable `MotionSnapshot` that implements
  `fs_geom::Chart` with time and path provenance frozen.
  `eval_over(x, span, cx)` returns a certified interval enclosing the
  base field value along the pulled-back trajectory of `x` for base
  charts that claim `TraceStepClaim::ExactDistance` (global
  1-Lipschitz theorem for exact signed distance); all other base
  claims refuse with a typed error.
- `SweptChart<C>` — the sign-correct swept-union implicit
  `ψ(x) = inf_t φ(M(t)⁻¹x)`. `evaluate` uses a deterministic complete
  binary cover of the time domain. Every active cell supplies a rigorous
  lower bound through `SpacetimeChart::eval_over`; every evaluated feasible
  time supplies a rigorous upper bound. Its `SweepReceipt` therefore encloses
  the global infimum even when the accuracy decision is `Unknown`. The
  implicit has the union's inside/outside sign but is **not generally a signed
  distance**, so the `Chart` façade publishes no numerical-distance,
  gradient, Lipschitz, topology, or ray-step claim.
- `EnvelopeChart<C, O>` — a `SweptChart` plus a chart-bound
  `EnvelopeOracle<C>`. The oracle must rigorously enclose `F` and `dF/dt`,
  certify characteristic-root existence (interval Newton, root count, or
  validated implicit-manifold continuation), certify rank margin, trimming,
  and visibility. A zero-containing interval alone never admits a branch.
  Interior regular branches, parameter endpoints, rank singularities,
  trimming, occlusion, and `Unknown` are distinct outcomes. The trace receipt
  carries the containing swept-field band; it is not an exact envelope SDF.
- `WankelApexPoint` and `WankelSealCircle` keep three constructions distinct:
  the ideal apex-point epitrochoid, a declared finite seal's center/contact
  loci, and the actual bore. The first two are derived directly from
  `base · T(e cos α, e sin α) · R(α/3 + β₀)`; the actual bore exists only as
  a visibility/trim-validated envelope of the seal-circle family. Rotor-flank
  conjugacy is a separate problem.
- `separation_over` — deterministic branch-and-bound in time over a
  chart-bound `ClearanceOracle`. Every retained cell has a complete-cover lower
  bound; a feasible fixed-time configuration supplies the optional upper
  bound. `ClearanceRange` says `TwoSided` only when both exist and otherwise
  remains explicitly `LowerOnly`. `ClearanceErrors` separately records chart
  conversion, spatial discretization, motion-tube/model, and optimization
  uncertainty in metres; both range endpoints are inflated outward.
- `SpherePairClearanceOracle` — an analytic production oracle for declared
  spherical body-frame proxies. Certified motor-tube point-action enclosures
  bound center distance over time cells; nonzero proxy errors bind those
  spheres back to the source charts. It can deliberately disable witnesses to
  exercise the one-sided lane.
- `OverlapInradiusWitness` — a deepest-common-interior point and certified ball
  radius from two exact-SDF bands. It is deliberately not named penetration
  depth, minimum translation distance, EPA depth, or SE(3) separating motion.
- `ChamberDefinition`, `ChamberChartFamily`, and `ChamberVolumeFunction` — a
  named, closure-gated `V(theta)` path. The family must construct the complete
  exact-distance chamber chart bounded by every named surface. Only
  `ProofState::Proven` closure reaches `fs-query::geometric_moments`; the result
  records the raw certified spatial quadrature plus chart-conversion,
  motion-tube, boundary-closure, and model-form volume errors.
- `IdealWankelVolumeOracle` — a declared ideal nominal G1 comparison formula
  with explicit generating radius, eccentricity, parallel transfers, depth,
  minimum volume, and phase. It is not a chart and cannot grant volume
  authority to an actual bore/flank/seal model.
- Analytic constructors: `screw_tube` (constant-twist screw about an
  axis line through a center, with translation along the axis) and
  `wankel_tube` (Wankel rotor **pose**: eccentric-center orbit at crank
  rate plus rotor phasing at 1/3 rate; the epitrochoid housing curve is
  the derived apex locus, deliberately NOT constructed here).
- `LowerToMotorTube` — the builder trait higher layers implement to
  lower their motion descriptions (frame trees, MBD trajectories) into
  tubes without this crate importing their types.

Trigonometric component models are built with the identity
`cos u = 1 − 2·sin²(u/2)` so no irrational constant enters a polynomial
coefficient; every transcendental evaluation goes through
`TaylorModel1::sin` with its rigorous remainder.

The double cover is fixed deterministically: at construction the
component models are sign-canonicalized so the scalar component's
midpoint at the domain anchor is positive (falling back to the first
even component exceeding a fixed tolerance in fixed blade order;
refusing as ambiguous when all are tiny). Piecewise tubes validate
chart transitions at every interior boundary — component enclosures of
adjacent segments must intersect and their representative vectors must
have positive dot product — BEFORE any consumer takes a logarithm.

## Invariants

- Geometric-product structure constants are extracted at runtime from
  `fs-ga` basis products (`OnceLock`), never transcribed by hand; the
  conformance suite falsifies the extracted table against
  `Motor::transform_point` on dense batteries.
- All sixteen component models of a segment share one domain and one
  order.
- Sign canonicalization and transition validation are deterministic
  functions of the inputs.
- `PointActionEnclosure` contains every real point `M(t)·x` for the
  constructed component path, for all `t` in the queried subinterval;
  AABB enclosures contain the image of every point of the box (rigid
  action is affine in `x`, so corner enclosures hull the box image).
- The versor-defect bound is an upper bound over the whole segment
  domain, not a sample statistic.
- At every point, the minimum lower endpoint over the current swept time-cell
  cover is a lower bound on the true infimum. The smallest upper endpoint at
  any evaluated feasible time is an upper bound. Splitting replaces one
  parent by two closed children whose union is the parent; no sample is used
  to establish completeness.
- Swept and envelope work selection uses total f64 ordering and fixed endpoint
  tie breaks. It is independent of worker count and scheduler order.
- An envelope branch is admitted as regular only after root existence,
  positive rank margin, in-trim status, and visibility are all `Proven` over
  its retained time enclosure. Endpoint surfaces are classified separately
  and never smuggled through the interior `dF/dt = 0` condition.
- The minimum clearance lower endpoint is the minimum of the lower bounds on a
  complete closed time-cell cover. Its upper endpoint, when present, comes
  from a feasible time and configuration. Splitting replaces exactly one
  parent by two closed children, and total-f64 ordering plus time endpoints
  fixes every work/tie decision.
- A chamber receipt exists only after named boundary closure is `Proven`, the
  family supplies an exact-distance chart, the integration domain contains its
  support, and `fs-query` completes a certified whole-region moments pass.
  Additional volume errors inflate the spatial band once and are logged
  individually.

## Error model

All fallible operations return `Result<_, MotionError>`:
`NonFiniteInput`, `EmptyTimeDomain`, `InvalidSegments`, `MixedModelShape`,
`Taylor(TaylorModelError)` (propagated fs-ivl refusals),
`DegenerateWeight` (homogeneous weight enclosure contains zero),
`DoubleCoverAmbiguous`, `ChartTransition` (adjacent segments fail the
overlap or sign test), `OutOfDomain`, `UnsupportedBaseClaim`
(`eval_over` on a base chart without `ExactDistance`),
`UnboundedSupport` (a swept construction cannot certify a finite support),
`InvalidConfiguration` (non-finite/negative tolerance), `InvalidEvidence`
(malformed, contradictory, or unavailable certificate evidence),
`InconsistentEnclosure` (independent lower and feasible-upper evidence
contradicts), `InvalidGeometry` / `PointActionFailed` (declared machine
geometry or its finite PGA action refuses), `Query(QueryError)` (propagated
certified moments/query refusal), and
`Cancelled` (cooperative cancellation observed). Panics are reserved
for programmer errors (violated internal invariants).

## Determinism class

Deterministic: identical inputs produce bit-identical component
coefficients, remainders, enclosures, and defect bounds on the same
ISA. All transcendentals route through `fs-ivl` Taylor models or
`fs_math::det`; no platform libm, no scheduler dependence (the G5
bit-replay conformance test pins this). Structure-constant extraction
is a fixed traversal of fixed basis products.

## Cancellation behavior

Loops over segments, box corners, clearance time cells, spatial quadrature
cells, and dense falsification samples poll
`cx.checkpoint()` at bounded strides and return
`MotionError::Cancelled` promptly. Single-segment scalar evaluations
are bounded-time and do not poll internally. Clearance minimization checks
immediately before and after every lower-bound or witness provider call, so a
pre-cancelled operation never enters the oracle and a cancellation observed
during provider work cannot publish or trigger the next provider call.
`SpacetimeChart::eval_over` also checks immediately after its injected base
chart evaluation, before constructing a certified field enclosure.

## Unsafe boundary

None. `#![forbid(unsafe_code)]`.

## Feature flags

None. Everything here is `[S]`-ambition machinery. Validated events remain a
separate bead.

## Validated events (bead 6b8h)

The `events` module implements the PRESCRIBED-PATH rung of doctrine
D3: no silent event misses. A [`GuardModel`] pairs a Taylor model of a
guard with a Taylor model whose band rigorously encloses the guard's
TRUE time derivative (the caller's class invariant; discharged by
construction for `plane_crossing_guard`, which weight-clears
`⟨M(t)·x, n⟩ − offset` after certifying the homogeneous weight band
strictly positive and builds the derivative by the sandwich product
rule from `screw_tube_with_derivative`'s rate tube — same constants,
differentiated basis, so both models describe one real function).
`scan_events` classifies a deterministic left-first bisection cover:
every leaf ends root-free (band excludes zero, or certified-monotone
with equal rigorous endpoint signs), a certified unique crossing
(certified-monotone, opposite endpoint signs, bisection-refined
without weakening the certificate), or an explicit `PossibleEvent`
(grazing, resolution floor, or budget) — nothing is dropped, and the
`RootCountCertificate` verdict is set-valued whenever anything stayed
possible. `confirmed` is then a certified lower bound, while
`possible_windows` counts unresolved regions rather than bounding their root
multiplicity. Each guard-model domain must exactly match its corresponding
tube-segment domain; both scan lanes refuse a mismatch before interval
evaluation. Zeno and subdivision budgets surface as typed verdicts, and the
Zeno cap never admits an extra certified event: untouched current and future
segment windows remain explicit possible events.
`estimated_scan` is the classical dense lane: LABELED Estimated, used
as the independent falsifier. `enumerate_simultaneous` groups
overlap-connected certified windows across guards and enumerates only orders
that preserve the proven precedence of every disjoint pair (set-valued; groups
above `MAX_ENUMERATED_GROUP` stay explicitly unordered). Piecewise double-cover
continuity is chained at construction: the first segment uses the anchor rule
and every later segment matches the previous segment's representative at the
shared junction, so rotations through π cannot tear the cover.

## Conformance tests

`tests/conformance.rs`:

- `mt_001` structure-table falsification: constant-motor sandwich
  enclosures contain `Motor::transform_point` results across rotor /
  translator / composed batteries.
- `mt_002` screw-tube containment: dense time sampling of pointwise
  fs-ga motors falsifies point-action enclosures (with a stated
  cross-implementation rounding tolerance).
- `mt_003` double-cover determinism: sign-flipped base poses produce
  bit-identical canonical components; a deliberate interior sign flip
  refuses with `ChartTransition`.
- `mt_004` residual honesty: exact-axis screws report tiny defects; a
  deliberately non-unit axis reports a large defect (the residual
  machinery must DETECT the broken construction, not mask it).
- `mt_005` `eval_over` containment against dense sampling of a moving
  exact-distance sphere (sampling falsifies, never proves).
- `mt_006` G5 bit replay across reconstruction.
- `mt_007` Wankel pose falsification against pointwise composition.
- `mt_008` box-action containment under sampled interior points.
- `mt_009` snapshot chart agreement with pulled-back base evaluation
  and support transport containment.

`tests/swept_envelope.rs`:

- cancellation requested from inside the final base-chart evaluation cannot
  publish a `FieldEnclosure`;
- `po_5` checks the deterministic exhaustive-subdivision accounting and
  falsifies every returned band against the independently derived translation
  capsule distance. The complete interval partition is the proof; closed-form
  point comparisons are only falsifiers.
- budget exhaustion returns `Unknown` with a still-sound interval, while the
  `Chart` façade retains `NoClaim` distance/ray semantics.
- exact parameter endpoints are classified separately from interior
  `dF/dt = 0` branches, and zero-containing field/derivative intervals with no
  root-existence proof remain `Unknown`.
- a derived rolling-rack family
  `F = x cos u + y sin u - r_b u` independently yields the rotated involute
  from `F = F_u = 0`; its validated fixture exercises regular branch
  admission and detailed rejection counters.
- rank-degenerate and invisible characteristic fixtures return `Unknown` or a
  typed rejection rather than a regular envelope branch.
- the derived Wankel apex formula is falsified against the independently
  composed pointwise motor and certified tube action. Finite seal center and
  contact loci remain distinct and no test equates either with the bore.

`tests/clearance_volume.rs`:

- a rotating eccentric spherical cam and stationary follower produce a
  two-sided minimum-clearance enclosure; the feasible witness is independently
  evaluated and all four length-error sources are logged;
- disabling witness production yields an explicit lower-only `Unknown`
  receipt, never a fabricated enclosure;
- a pre-cancelled clearance request refuses before entering either oracle
  method;
- exact-SDF bands yield a common-interior inradius witness without relabeling
  it as global penetration or pose displacement;
- certified quadrature encloses a slider-crank cylinder's independently
  derived G1 volume at several crank angles, with every additional volume
  error logged;
- the ideal nominal Wankel sinusoid is checked through an explicitly
  equivalent-volume manufactured quadrature fixture. A separately named actual
  Wankel chamber with unproven bore/flank/seal closure refuses before chart
  construction. The manufactured fixture makes no Wankel geometry claim.

`tests/events.rs` (bead 6b8h, printed measurements): known analytic
root counts found and certified with directions (ev-001); top-level
no-event certificates (ev-002); tangency yields `PossibleEvent`
grazing windows, never a certified crossing — the Sev-0 gate (ev-003);
a 100×-resolution dense-scan falsifier finds no sign change outside
certified/possible windows (ev-004); duplicated-guard simultaneous
groups enumerate both admissible orders while disjoint crossings stay
separate (ev-005); subdivision-budget exhaustion drains into visible
possible windows with a typed verdict (ev-006); bitwise-deterministic
replay (ev-007); certified windows refine to sub-1e-3 localization
around analytic roots (ev-008); guard/tube domain mismatch refuses in both scan
lanes (ev-009); receipt width includes certified leaves (ev-010); the Zeno cap
is strict and preserves the untouched tail as possible (ev-011); chain-overlap
order enumeration retains precedence between disjoint members (ev-012).

## No-claim boundaries

- **Rigid paths only.** No deformable sweeps, no scaling, no shear.
- Event completeness is claimed ONLY for the Taylor-pair guard class
  (guard + certified derivative band describing one real function). A
  zero-containing band alone never admits a certificate; black-box
  guards get `estimated_scan`'s labeled Estimated outputs. Certified
  windows locate SIGN CROSSINGS of the weight-cleared guard; roots at
  segment junctions or bisection cut points may surface as
  possible-event windows instead of certified crossings. Even-order
  touch points (grazing) are always `PossibleEvent`, never counted.
  Root-count completeness is relative to the scan span. With no possible
  windows, `confirmed` is the exact root count; otherwise it is only a
  certified lower bound, and the unresolved windows carry no finite
  root-multiplicity upper bound in v1. No interval-Newton tightening is
  performed in v1; the mode-ledger and true-flow (ValidatedStep) rungs live in
  fs-time bead ow2o, not here.
- Analytic constructors enclose their **constructed component path**;
  the deviation of that path from the ideal real-number screw or
  Wankel motion is bounded by the reported versor defect, not
  separately certified. A caller needing `defect = 0` semantics must
  check the reported bound against its own tolerance.
- Sampling-based tests falsify; they never certify. Only Taylor-model
  enclosures carry `EnclosureClass::Certified`.
- `MotionSnapshot` deliberately claims nothing for ray stepping
  (`TraceStepClaim::NoClaim`), reports `NumericalCertificate::no_claim`
  for sample error, and drops gradient/Lipschitz data: transporting
  those claims through an approximate motor is future work. Its
  `topology_hint` passes through the base chart's bounds (an invertible
  near-rigid pull-back is a homeomorphism of the zero set).
- `eval_over` supports only `ExactDistance` base charts in this
  version; Lipschitz-implicit transport refuses rather than guessing.
- A `SweptChart` interval certifies the infimum of the transported base
  implicit and hence the swept-union sign. It makes no distance-accuracy claim
  and registers no Rep Router conversion edge in v1. Future swept-to-mesh or
  swept-to-voxel edges must preserve the enclosure class and account for their
  own discretization error.
- `EnvelopeOracle` is the trusted certifier boundary for characteristic-root
  existence, derivatives, rank, trim, and visibility. Opaque charts are never
  finite-differenced and promoted. `Unknown` rank/visibility/existence remains
  `Unknown`; mutually contradictory oracle fields also demote to `Unknown`.
  No global watertightness, uniqueness, self-intersection freedom, or exact
  boundary-distance claim follows from a trace receipt.
- The Wankel apex expression is an ideal point locus. `WankelSealCircle`
  exposes a declared center and a contact point for a supplied unit normal;
  supplying that normal does not prove envelope visibility. The actual bore
  remains an `EnvelopeChart` problem, and the rotor flank remains a separate
  conjugate-envelope design problem. The rotor is not claimed to be an exact
  constant-width triangle, and no rolling-seal assumption is made.
- A spherical clearance proxy certifies the source bodies only to the extent of
  its explicit conversion/spatial/motion error budget. Generic nonconvex
  minimum translation, EPA/GJK penetration, and pose-space separating motion
  remain `fs-query` upgrades; local contact response must use local gaps and
  normals rather than `OverlapInradiusWitness` as a global displacement.
- Chamber volume is authority for the exact named chart family and seal
  convention only. The ideal Wankel formula is a G1 oracle, not evidence that a
  finite-seal bore and conjugate rotor flank close. Ports, clearances,
  deformation, recesses not included in `minimum_volume_m3`, and as-built
  geometry require their own named charts and certified quadrature errors.
- No claim of tightness: enclosures are sound, not minimal; segment
  count and Taylor order are the caller's accuracy budget.
- Double-cover canonicalization is bit-deterministic for identical
  inputs, but a tube built from `−M` canonicalizes to enclosures that
  agree with the `M`-built tube only at roundoff scale (≤ 1e−13
  absolute per endpoint at unit component scale): the Taylor-model
  negation's outward remainder rounding is not perfectly
  sign-symmetric. Exact bit-equality ACROSS the double cover is not
  claimed.
