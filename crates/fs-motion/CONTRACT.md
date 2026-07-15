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

## Error model

All fallible operations return `Result<_, MotionError>`:
`NonFiniteInput`, `EmptyTimeDomain`, `InvalidSegments`, `MixedModelShape`,
`Taylor(TaylorModelError)` (propagated fs-ivl refusals),
`DegenerateWeight` (homogeneous weight enclosure contains zero),
`DoubleCoverAmbiguous`, `ChartTransition` (adjacent segments fail the
overlap or sign test), `OutOfDomain`, `UnsupportedBaseClaim`
(`eval_over` on a base chart without `ExactDistance`),
`UnboundedSupport` (a swept construction cannot certify a finite support),
`InvalidConfiguration` (non-finite/negative tolerance),
`InconsistentEnclosure` (independent lower and feasible-upper evidence
contradicts), `InvalidGeometry` / `PointActionFailed` (declared machine
geometry or its finite PGA action refuses), and
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

None. Everything here is `[S]`-ambition machinery. Clearance volumes and
validated events remain separate beads/crates.

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
- No claim of tightness: enclosures are sound, not minimal; segment
  count and Taylor order are the caller's accuracy budget.
- Double-cover canonicalization is bit-deterministic for identical
  inputs, but a tube built from `−M` canonicalizes to enclosures that
  agree with the `M`-built tube only at roundoff scale (≤ 1e−13
  absolute per endpoint at unit component scale): the Taylor-model
  negation's outward remainder rounding is not perfectly
  sign-symmetric. Exact bit-equality ACROSS the double cover is not
  claimed.
