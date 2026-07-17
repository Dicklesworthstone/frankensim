# fs-geocon — CONTRACT

## Purpose and layer

L4 (ASCENT-facing geometry). Differentiable/certifiable geometric
constraint primitives (plan §7.6): the manufacturability and
design-intent constraints ASCENT enforces, each FIRST-CLASS — a value,
a derivative story, a certificate story, and its fs-constraint kind
mapping. Geometry supplies values and derivatives; ASCENT owns the
constraint SEMANTICS.

## Public types and semantics

- `min_thickness_soft` / `min_thickness_soft_clipped` — the anti-paperclip
  constraint: mean p-norm aggregation of fs-query thickness samples.
  The default API requires a finite chart support; the clipped API is an
  explicitly LOCAL report over the caller's finite AABB and never promotes
  that result into a global minimum-thickness claim. `soft_min` is the C¹
  optimizer value (an over-approximation converging DOWN to the
  minimum as p grows; exact on uniform samples, which keeps lever
  derivatives clean); `hard_min` and the LOCALIZED violation list are sampled
  estimates, not ledger certificates. `ThicknessReport::authority` preserves
  the generic oracle's explicit `Estimate` class. The two values are reported
  side by side, never conflated. Only geometric local misses are counted as
  skips; malformed samples/arithmetic and cancellation propagate as refusals,
  and an empty/all-skipped aggregate returns `NoThicknessSamples`.
- `draft_violations` — normals versus the pull direction for the mold
  half being assessed: smooth hinge² penalty (C¹) plus EXACT violating
  regions; normals opposing the pull within the half's own reach are
  UNDERCUTS, flagged separately; faces beyond `n·d < −0.5` belong to
  the other mold half (the v1 parting model is the plane perpendicular
  to the pull).
- `QuotientChart` + `SymmetryGroup` (`ReflectX`, `Cyclic{n}`,
  `Periodic{period}`) — symmetry ENFORCED BY CONSTRUCTION: evaluation
  points fold into the fundamental domain, so the shape is invariant
  for ARBITRARY inner designs and lever values; gradients chain back
  through the fold. Reflection folds bitwise; cyclic/periodic at fp
  scale. Periodic repetition publishes an honest x-unbounded extended AABB
  while preserving the inner chart's transverse bounds. The seam makes the
  chart C0 (declared), clears sample Lipschitz authority, and demotes valid
  finite inner evidence to `Estimate` while preserving its full numerical band
  (`NoClaim` and malformed or nominal-excluding certificates absorb): symmetry
  alone is not an abstract-distance theorem. Cyclic-orbit support radii are rounded outward to
  preserve support containment. `SymmetryGroup::cyclic`/`periodic` validate new
  values; directly constructed invalid public enum values and malformed inner
  supports fail closed before fold/support arithmetic.
- `envelope_violation` — containment (`φ_allowed ≤ 0` on design
  boundary samples) and keep-outs (`flip = true`): the sampled worst
  plus a SUM-FORM log-sum-exp aggregate that is `≥ worst` within
  `ln(n)/β` — CONSERVATIVE, so driving the smooth value to zero drives
  the true worst to zero (never stops short).
- `volume_certified` / `volume_smooth` — over an EXPLICIT integration
  domain (fixed independently of levers, so derivatives see the shape
  change, not grid realignment): the box is admitted through
  `fs_geom::SamplingDomain`, `h` is a finite positive maximum cell width,
  and normalized cell-center placement encloses exact rational centers with
  directed rounding. A RIGOROUS enclosure (sure-inside cells vs an
  outward-rounded L1 radius covering width and center-placement error)
  requires `TraceStepClaim::ExactDistance` plus a finite
  rigorous `Exact`/`Enclosure` certificate at every cell center. A local
  Lipschitz sample, `Estimate`, `NoClaim`, or malformed rigorous certificate
  refuses rather than authorizing `[lo, hi]`. This sits beside the deliberately
  non-certifying smoothed-Heaviside estimate whose lever derivative matches
  Hadamard on fixtures. `VolumeError` reports all preflight, authority,
  evaluation, and cancellation failures; `VOLUME_MAX_CELLS` is the shared
  deterministic work cap.
- `GeoPrimitive::descriptor()` — the declared table: differentiability
  class (fs-opt `Class`), certificate story (`Enclosure` /
  `SmoothEstimate` / `ExactByConstruction`), and the fs-constraint
  `ConstraintKind` mapping; `proof_escalation()` names the interval
  obligations where they exist.

## Invariants

1. Thickness: soft over-approximates converging down with p;
   violations localize to exactly the thin samples; the lever FD
   derivative matches the analytic 2 on the neck; a toy descent DRIVES
   the neck to feasibility (gcp-001). Unresolved extended support propagates
   as `UnboundedSupport` rather than becoming a skipped sample, while an
   explicit finite clip yields a deliberately local report (gcp-001b).
2. Draft: a 10° wall passes 5° and fails 15° with all samples
   localized; vertical walls violate any positive draft; a mushroom
   shoulder is flagged as an UNDERCUT, not low draft; the smooth
   penalty's FD derivative matches the analytic hinge slope (gcp-002).
3. Symmetry: quotient shapes are invariant under their groups for
   ARBITRARY asymmetric inner designs (bitwise for reflection,
   fp-scale for cyclic/periodic) across random levers × groups ×
   orbits; folded gradients match finite differences off-seam; periodic
   support is infinite along x with finite transverse bounds — violation is
   structurally impossible (gcp-003).
4. Envelopes: containment/keep-out match analytic penetrations; the
   LSE aggregate is conservative within `ln(n)/β`; the FD derivative
   is right; the descent returns an escaping design fully inside
   (gcp-004).
5. Volume: certified enclosures bracket the analytic sphere volume at
   two resolutions and tighten with h; the smoothed volume's lever
   derivative matches Hadamard `4πr²` within 2%; a descent meets a
   volume cap (gcp-005). Malformed/unbounded domains, invalid spacing,
   checked-count overflow, and grids beyond `VOLUME_MAX_CELLS` refuse before
   chart evaluation, while a normal finite grid remains usable (gcp-005b).
6. The descriptor table is total, symmetry is `ExactByConstruction`,
   volume is `Enclosure`, thickness maps to Fabrication, and proof
   escalations are declared exactly where they exist (gcp-006).

## Error model

fs-query's teaching errors carry through thickness queries (`SamplingDomain`,
`NoGradient`, `NotOnBoundary`, `Cancelled`, …); envelope assessment is total.
Local thickness failures are counted, but domain-admission failures propagate
and cannot be hidden as skipped samples. Volume APIs return `VolumeError`:
domain admission and finite-positive `h`/`epsilon` validation precede count
math; per-axis counts and their product are checked; the deterministic work cap
is enforced before evaluation; weak chart theorems, weak or malformed
per-sample certificates, non-representable cell measures, and non-finite chart
samples are structured failures rather than false enclosures or NaN results.

## Determinism class

Fully deterministic: fixed sample sets, canonical iteration, no
randomness. Identical inputs give identical reports bitwise, CROSS-ISA:
every transcendental routes through `fs_math::det` (bead frankensim-lyms;
platform libm is not correctly rounded and differs across ISAs), and the
crate is registered in the `check-libm` doctrine lint that keeps raw
libm from reappearing. `sqrt` stays primitive (IEEE-754 correct
rounding).

## Cancellation behavior

Volume integrations poll `cx.checkpoint()` at most every 256 completed cells
and once before publication, returning `VolumeError::Cancelled` with the exact
completed-cell count; thickness aggregation polls through the fs-query oracle
and returns its carried `QueryError::Cancelled` teaching error.

## Unsafe boundary

None. `#![forbid(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None.

## Conformance tests

`tests/conformance.rs`, aggregate cases gcp-001, gcp-001b, gcp-002,
gcp-003, gcp-004, gcp-005, gcp-005b, and gcp-006 emit canonical
`fs_obs::EventKind::ConformanceCase` records after their direct checks complete.
The records use `Info`/`Error` severity, pass failure-record linting and wire
validation before printing, and print before the aggregate assertion. Seeded
quotient campaign gcp-003 carries its literal LCG input seed
`0x1001_2026_0707_0023`; all fixed aggregate cases use zero. Cases executed
through `with_cx` record the fixed Cx execution seed `0x6C0` separately in
detail rather than presenting it as input randomness. The existing
volume/Hadamard and descriptor-table `Custom` companion events remain
object-shaped, wire-validated, and printed; the volume event records the fixed
execution seed separately. Assertions and expectations before aggregate
emission remain ordinary Rust test diagnostics and can terminate a case
without publishing a verdict. In particular, fail-closed matrices gcp-001c
and gcp-005c are direct-assertion tests and intentionally emit no aggregate
record. Any reimplementation must pass the suite unchanged.

## No-claim boundaries

- Derivatives here are FD demonstrations against analytic truths;
  exact adjoints join the gradient-stack bead.
- The draft parting model is the plane perpendicular to the pull;
  parting-line OPTIMIZATION and multi-pull molds are follow-ups.
- Symmetry groups are the plan's named trio about fixed axes/planes;
  arbitrary group generators and lattice groups follow with fs-ga.
- Clipped thickness evidence covers only the recorded finite AABB. It makes no
  claim about thinner features elsewhere on an unbounded chart.
- Certified volume requires the global `ExactDistance` theorem plus a finite
  rigorous trace enclosure at every cell center. Partition spans, widths,
  center-placement error, conservative L1 cell radii, cell measure, and final
  count products use directed outward rounding; local Lipschitz samples or
  generic implicit fields cannot authorize the enclosure. Tighter interval
  cell geometry remains an optimization, not a missing soundness premise.
- Mass is volume × constant density; density fields join fs-material.
