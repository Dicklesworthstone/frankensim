# fs-motion CONTRACT

Certified rigid motion for the MORPH layer: a motor path over a time
domain becomes a checkable object (a *tube*) instead of a per-evaluation
pose, and moving geometry exposes time-aware queries instead of
pretending to satisfy the timeless `Chart` contract.

## Purpose and layer

Layer **L2 (MORPH)**. `fs-ga` owns instantaneous SE(3) motor algebra and
`fs-scenario` owns frame trees; nothing previously bound a motor *path*
to a chart. `fs-motion` provides that binding. Dependencies: `fs-ga`,
`fs-geom`, `fs-ivl`, `fs-exec`, `fs-evidence`, `fs-math`. This crate
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

## Error model

All fallible operations return `Result<_, MotionError>`:
`NonFiniteInput`, `EmptyTimeDomain`, `InvalidOrder`,
`Taylor(TaylorModelError)` (propagated fs-ivl refusals),
`DegenerateWeight` (homogeneous weight enclosure contains zero),
`DoubleCoverAmbiguous`, `ChartTransition` (adjacent segments fail the
overlap or sign test), `OutOfDomain`, `UnsupportedBaseClaim`
(`eval_over` on a base chart without `ExactDistance`),
`UnboundedSupport` (snapshot support transport of an infinite box), and
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

Loops over segments, box corners, and dense falsification samples poll
`cx.checkpoint()` at bounded strides and return
`MotionError::Cancelled` promptly. Single-segment scalar evaluations
are bounded-time and do not poll internally.

## Unsafe boundary

None. `#![forbid(unsafe_code)]`.

## Feature flags

None. Everything here is `[S]`-ambition machinery; swept envelopes,
clearance volumes, and validated events are separate beads/crates.

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

## No-claim boundaries

- **Rigid paths only.** No deformable sweeps, no scaling, no shear.
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
- The Wankel constructor produces the rotor **pose** path; the
  epitrochoid apex locus and swept/envelope geometry belong to the
  envelope bead (`fs-motion` swept-volume follow-up) and are not
  claimed here.
- No claim of tightness: enclosures are sound, not minimal; segment
  count and Taylor order are the caller's accuracy budget.
