# fs-query — CONTRACT

## Purpose and layer

L2 (MORPH). Geometry queries (plan §7.4): the interrogation layer
every consumer calls constantly (FLUX embedding, ASCENT constraints,
LUMEN), UNIFORM across chart types — every query speaks `&dyn Chart`,
so the same call runs against analytic fixtures, F-rep CSG, dense SDF
grids, and mesh charts, and the conformance battery holds their
answers to the MULTI-CHART AGREEMENT discipline (same abstract region
⇒ same answers within composed certificates).

## Public types and semantics

- `closest_point` / `closest_point_clipped`: damped Newton projection along the chart gradient;
  the post-projection RESIDUAL is measured and reported, never
  assumed. Charts that honestly decline gradients (mesh charts near
  edges) fall back to a central FD on the signed distance. Analytic gradients
  need no sampling-domain admission; an unbounded chart's FD fallback refuses
  unless the caller supplies a finite clip from which to derive the stencil
  scale. Every source point, nominal value, gradient, FD subtraction, robustly
  scaled Newton update, and final residual must remain finite; producer-side
  cancellation is checked immediately after each evaluation. The residual
  still carries the honesty.
- `raycast`: conservative sphere tracing from each sample's rigorous
  trace-value enclosure and certified local Lipschitz bound. The endpoint
  actually produced by floating-point ray evaluation must remain inside the
  certified safe ball. Clean misses are `None` only after the caller's `tmax`
  endpoint is classified; grazing or rounding stalls return an explicit
  `UnresolvedTrace` error. A Lipschitz-implicit chart's `|f|/L` certifies only
  the step radius; without separate proximity evidence, only a rigorous exact
  zero can authorize `RayHit` — incomplete, never unsafe.
- `geometric_moments` → `GeometricMoments` (`MomentEnclosure` volume,
  first moments, `SecondMoments` about the origin): certified
  unit-density integrals over a chart's region by exact-distance cell
  classification with outward rounding — sure-inside cells add exact
  closed-form box integrals, straddling cells add conservative
  brackets, and everything degraded refuses (capability, domain,
  spacing, per-sample evidence class, `MAX_MOMENT_CELLS` work,
  cancellation). `com_enclosure` divides only through a proven-positive
  volume; `translated` applies the exact translation-covariance /
  parallel-axis laws with outward rounding. The domain must contain the
  chart's support box: moments are whole-region claims.
- `ConvexSupportMap` (`ConvexSphere`, `ConvexBox`) and
  `convex_separation` → `ConvexSeparation`: certified `[lo, hi]`
  distance enclosures between compact convex sets via deterministic
  Frank-Wolfe on the Minkowski difference. The upper bound is the
  outward-rounded norm of a realized support-point difference plus the
  shapes' declared support slack; the lower bound is the
  outward-rounded support-plane separation minus slack, clamped at
  zero. `separation_proven` is exactly `lo > 0`; every early stop
  (iteration cap, nonsmooth `1/k` residual) widens the bracket, never
  falsifies it. Constructors validate geometry; degraded arithmetic
  refuses typed.
- `FeatureComplex` (`Feature::{Vertex, Edge, Face}`) and
  `ccd_candidates`: the typed vertex/edge/face decomposition of a
  triangle boundary with outward-rounded per-feature boxes, and
  conservative CCD candidate enumeration through a deterministic
  median-split BVH. Each side's boxes are inflated by its declared
  motion bound with upward rounding, so every feature pair that could
  approach within the combined window is INCLUDED; the output is a
  lexicographically sorted candidate SUPERSET (a pure function of the
  inputs), and exceeding the caller's pair cap refuses rather than
  truncates. Degenerate triangles, bad indices, non-finite positions
  or inflations, and the `MAX_COMPLEX_FEATURES` ceiling refuse typed.
- `ImplicitGapOracle` → `GapSample`: pointwise SDF-pair contact
  queries over two exact-distance charts (weaker claims refuse at
  construction, reusing `SeparationRequiresExactDistance`). The
  outward-rounded `[sum_lo, sum_hi]` encloses `φ_A(p) + φ_B(p)`;
  `separation_upper` (triangle inequality) exists exactly when the
  point is certified outside both bodies; `overlap_inradius` is a
  certified common-ball radius exactly when `max(φ_A, φ_B)` is
  certified negative — a pointwise witness, never a penetration depth;
  `normal` is an uncertified Estimate-class contact axis that is
  absent whenever a gradient is honestly declined. Malformed evidence
  refuses per sample; cancellation is checked after each chart call.
- `OffsetChart` / `minkowski_ball`: dilation/erosion as a chart
  wrapper (`φ − r`); the ball case of the Minkowski sum IS the offset
  (bitwise), which is the fillet/clearance workhorse. Construction is fallible
  and rejects non-finite radii. The wrapper validates the inner certificate after every
  sample, outward-translates its full numerical band, and caps the transformed
  evidence at `Estimate`; malformed, nominal-excluding, overflowing, or
  `NoClaim` evidence remains `NoClaim`. The shifted nominal is included in the
  transformed hull and is never paired with stale or collapsed inner bounds.
- `ClearanceField` (`c(p) = φ_A⁺ + φ_B⁺`) + `separation` /
  `separation_clipped`: `ClearanceField::value` is a nominal convenience only.
  A rigorous bracket requires `TraceStepClaim::ExactDistance` from both charts
  and validates each sample's finite `Exact`/`Enclosure` trace certificate.
  Grid lower endpoints minus the exact-distance field's 2-Lipschitz
  nearest-node slack give `lower_bound`; the smallest certificate upper
  endpoint gives `observed`. Each raw support is validated before union,
  finite-domain admission and checked `(n + 1)^3` arithmetic happen before
  evaluation, and the public `SEPARATION_MAX_CHART_SAMPLES` cap refuses
  representable but excessive work before evaluation. `SeparationScope`
  distinguishes complete finite-support authority from a caller-clipped LOCAL
  bracket.
- `thickness_at` / `min_thickness` and their `_clipped` counterparts: the THICKNESS ESTIMATOR —
  inward-normal march + bisection to the opposite wall; per-sample
  failures are SKIPPED AND COUNTED only after the chart support resolves to a
  finite domain. Domain failures propagate rather than becoming skips; an
  explicit finite clip enables deliberately local queries on unbounded charts.
  March scale and termination use the inward ray's AABB-exit distance, so a
  large transverse clip extent cannot jump over an otherwise nearby opposite
  wall.
  Generic implicit-field magnitudes have no no-skipped-zero theorem, so both
  `Thickness` and `ThicknessMinimum` explicitly carry
  `NumericalKind::Estimate`; no caller may describe them as certified.
  Non-finite values, gradients, query points, and march arithmetic refuse.
  Positive parameters must also move to a distinct representable world point
  during marching and bisection; an unresolvable translated/thin wall refuses
  rather than returning the original boundary as its own opposite.
  The estimates respond smoothly to design
  levers where the walls are smooth (the battery FD-differentiates
  min-thickness through an F-rep neck radius and reads the analytic
  subgradient 2).
- `medial_poles`: interior circumcenters of the Delaunay of a
  boundary sample set, λ-filtered by local sample spacing — the
  medial-axis approximation that cross-checks the oracle (2·pole
  radius ≈ local thickness). Boundary/tetrahedron loops poll cancellation;
  chart-requested cancellation is rechecked immediately after evaluation and
  again before publication, and non-finite pole samples refuse.
- `curvature`: mean/Gaussian/principal from central stencils on the
  reported signed field (shape operator = tangent-restricted Hessian for a
  unit-gradient distance field), with
  a caller-supplied positive finite `h` with a representable squared stencil
  scale that also drives the gradient FD fallback, and
  a PER-CHART ACCURACY CLASS (`CurvatureClass`): `SecondOrder`
  (analytic/F-rep — O(h²), measured), `GridLimited` (C¹ grids — error
  floors at the grid's own interpolation error), `Estimate`
  (mesh fields are non-smooth across facets and may themselves carry only
  estimate/no-claim distance authority). Stencil points and samples,
  differences, Hessian projections, invariants, and final principals are all
  checked finite; this prevents an `Ok` result containing NaN/Inf but does not
  upgrade an accuracy class into a numerical certificate.

## Invariants

1. Closest points agree with analytic truth across all four chart
   families within each chart's OWN certificate (exact/F-rep at fp,
   tiled at its declared bound, mesh at faceting scale), residuals
   are honest, and answers are translation-equivariant. Malformed producer
   samples, overflowing FD points/Newton arithmetic, and producer-requested
   cancellation refuse before publication (gq-001..gq-001a).
2. Raycasts match analytic hits across chart types; tangent rays
   never tunnel (grazes land on the surface or report unresolved); the CSG
   tracer never claims a hit past a dense oracle, and every sample including
   `tmax` revalidates its local Lipschitz and rigorous trace evidence, and
   cancellation requested inside either chart producer wins before hit/miss
   authority, and a loose valid implicit-field `L` cannot promote a small
   normalized residual into a geometric hit (gq-002..gq-002e).
3. Offsets of spheres are exactly spheres of the summed radius;
   erosion shrinks exactly; `minkowski_ball` is BITWISE the offset;
   offset charts retain closest-point and other differential queries; generic
   raycast remains `NoClaim` until a reach/proximity theorem is supplied
   (gq-003).
4. Exact-distance separation brackets hold across shrinking gaps (truth in
   `[lower_bound, observed]`); local Lipschitz/enclosure fields cannot upgrade
   `NoClaim`, malformed per-sample evidence refuses, and cancellation requested
   by either producer wins before bracket authority (gq-004..gq-004b).
5. The thickness estimator reads the graded slab analytically (1% rel),
   finds the dumbbell neck (2× neck radius, zero skips), agrees with
   the medial-pole cross-check, and differentiates through a design
   lever with the analytic subgradient. Every result remains explicitly
   `Estimate`; malformed samples, empty aggregates, and producer cancellation
   fail closed (gq-005..gq-005a).
6. Curvature converges at measured order ≈2 on SecondOrder charts,
   torus principals hit 1/r and 1/(R+r), classes are documented per
   family, grid-limited charts land within their own scale, and
   curvature scalars are rotation-invariant. Malformed samples, overflowing
   stencil arithmetic, and producer cancellation fail closed rather than
   publishing NaN/Inf scalars (gq-006..gq-006a).
7. Support-derived samplers refuse malformed, unresolved-unbounded,
   degenerate, or overflowing domains before source evaluation or span/count
   loops; analytic closest-point queries bypass admission, explicit curvature
   `h` drives its FD gradient, clipped separation is marked local, excessive
   work refuses under a public cap, ratio-first extreme finite coordinates
   remain representable, and thickness aggregates propagate domain failures
   (gq-007).

## Error model

`QueryError` teaching errors: `NoGradient` (with the location),
`SamplingDomain` (the structured `fs-geom` admission refusal),
`InvalidOffsetRadius`, `InvalidFiniteDifferenceStep`, `InvalidPointSample`,
`InvalidPointArithmetic`, `InvalidBoundaryIndex`, `SamplingGridTooLarge`,
`SamplingWorkLimitExceeded`, `InvalidSeparationArithmetic`,
`NoLipschitz`, `NoTraceClaim`, `SeparationRequiresExactDistance`, `InvalidRay`,
`InvalidTraceSample` (with the location), `InvalidThicknessSample`,
`InvalidThicknessArithmetic`, `NoThicknessSamples`,
`UnresolvedTrace` (with the location and sample count), `NotOnBoundary` (with
the sd found and the advice to project first), `NoOppositeWall`,
`MomentsUncertifiedChart` (moments demand the exact-distance capability),
`MomentsInvalidDomain`, `MomentsInvalidSpacing`, `MomentsExcessiveWork`
(the deterministic `MAX_MOMENT_CELLS` ceiling), `MomentsInvalidSample`
(Estimate/NoClaim-class evidence cannot feed a certified integral),
`MomentsVolumeUnproven` (COM needs a positive certified volume lower
bound), `ConvexInvalidShape` (non-finite/degenerate convex geometry or
a zero iteration budget), `ConvexInvalidSupport` (non-finite support
evaluation or bound arithmetic), `FeatureComplexTooLarge`,
`FeatureInvalidInflation`, `FeatureTooManyPairs` (cap refusal, never
silent truncation), `Cancelled`,
`Mesh` (fs-mesh refusals carried through). Honest gaps refuse; nothing guesses.

## Determinism class

Fully deterministic: fixed iteration counts, canonical grid orders,
no randomness. Identical inputs give identical answers bitwise.

## Cancellation behavior

`raycast` polls before each sample and again after `eval` and
`trace_value_enclosure`; `separation` polls immediately before and after every
producer call and once before publishing; thickness polls around every chart
evaluation and before publishing local or aggregate estimates. Producer-
requested cancellation therefore wins over evidence returned by the same call.
Closest-point and curvature queries likewise poll around every analytic or FD
sample and once before publishing. All return `Cancelled` teaching errors.
Other point queries are O(iterations) and non-blocking.

## Unsafe boundary

None. `#![forbid(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None.

## Conformance tests

`tests/conformance.rs`, cases gq-001..gq-007 (+ typed trace refusal checks)
— JSON-line verdicts, seeded LCG randomness, fs-obs events for the
thickness estimator and curvature convergence tables. Any
reimplementation must pass the suite unchanged.
`tests/moments.rs`, cases gm-001..gm-005 — certified geometric moments:
box/sphere closed-form containment with bounded widths, COM enclosures,
the translation-covariance metamorphic (outward-rounded law vs direct
recomputation must overlap), capability/input/work refusals, and
cancellation. `geometric_moments` is Enclosure-class: every returned
bracket contains the true unit-density integral of the chart's region.
`tests/convex.rs`, cases gc-001..gc-005 — certified convex separation:
analytic sphere/box distance containment, nonsmooth box-box honesty,
touching/overlap never claiming separation, bit-identical replay, and
constructor/budget/cancellation refusals.
`tests/features.rs`, cases gf-001..gf-005 — feature complexes and CCD
candidates: canonical construction, BVH-equals-brute-force oracle
agreement, guaranteed inclusion of the closest pair inside the motion
window, translation invariance and window monotonicity metamorphics,
identical replay, and refusal drills.
`tests/gap.rs`, cases gg-001..gg-004 — implicit gap oracle: analytic
containment and true-distance upper bounding on disjoint spheres,
certified common-ball witnesses on overlap, claim/evidence/point/
cancellation refusals, and no normal claim without gradients.

## No-claim boundaries

- General Minkowski sums (non-ball structuring elements, max-plus /
  FFT-assisted convolution) are deferred; the exact ball case is the
  v1 surface.
- The medial-axis approximation is pole-based (filtered Delaunay
  circumcenters of boundary samples); full filtered-Voronoi medial
  complexes with angle criteria and stability guarantees are the
  follow-up.
- Thickness values and subgradients are Estimate-authority demonstrations; a
  certified inward marcher and exact adjoints join later proof lanes.
- Separation deliberately accepts only exact Euclidean signed-distance charts.
  A safe broader theorem for Lipschitz-implicit fields and interval cell bounds
  could widen coverage and tighten the current global 2-Lipschitz slack.
- `separation_clipped` and clipped thickness queries cover only their recorded
  finite domain; a clip never upgrades a local answer into a global claim over
  an unbounded chart.
- Curvature on mesh charts is an ESTIMATE class by design; discrete
  curvature operators (cotan/normal-cycle) on the half-edge mesh are
  a separate surface.
- `geometric_moments` is Enclosure-class only under
  `TraceStepClaim::ExactDistance`; weaker claims refuse rather than
  downgrade. Band brackets are conservative, not convergence-optimal:
  no width-vs-`h` rate is claimed. Moments are unit-density geometric
  integrals; densities, inertia tensors, and material identity live
  downstream (fs-matdb consumers), never here. Rotation covariance and
  spatially-varying weighting are deferred surfaces —
  [`GeometricMoments::translated`] covers translation only. The bead
  rjnd remainder (codimensional thickness, deformation hooks) is not
  yet claimed by this crate.
- `ImplicitGapOracle` is pointwise: `separation_upper` bounds the
  distance between the two bodies from above at ONE probe, and
  `overlap_inradius` witnesses a common ball at ONE probe. Neither
  aggregates into a global gap field, a penetration depth, or a
  time-of-impact claim, and the `normal` estimate carries no
  certificate whatsoever.
- `ccd_candidates` is a broad phase only: candidates are a
  conservative SUPERSET under the caller's declared motion bounds, and
  no narrow-phase, time-of-impact, or contact claim is made. The
  motion bound is a caller assertion (a certified radius over the CCD
  window); the routine cannot detect an understated bound. Boxes are
  axis-aligned in the input frame — rotation-aware swept bounds are a
  deferred surface.
- `convex_separation` proves separation only (`lo > 0`); an enclosure
  containing zero claims NOTHING about overlap, and penetration-depth
  certificates (EPA-style lower bounds) are a deferred rjnd surface.
  No convergence rate is claimed: nonsmooth pairs may return wide,
  honest brackets at the iteration cap. Support maps must satisfy the
  declared slack contract; the routine cannot detect a lying
  `support_slack()`.
- Chart-native fast paths (mesh BVH closest-point dispatch instead of
  generic Newton) are perf-lane work; answers here are correct first.
