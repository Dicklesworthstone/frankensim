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

- `ContactInflation`: a private-radius, proof-bearing absolute-error carrier
  for contact-adjacent queries. `from_conversion` accepts only an unforgeable
  `Certified<T>` receipt, `from_motion` does the same for absolute motion-model
  error, and `from_route` independently revalidates a router `ChainOutcome` as
  finite, ordered, rigorous, nonnegative absolute-error evidence containing
  its QoI. Composition is outward-rounded and refuses overflow. `exact_zero`
  is reserved for native geometry with no unrecorded conversion or motion
  uncertainty and is a bitwise identity.
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
  `geometric_moments_with_inflation` adds the retained absolute error to every
  cell-classification radius, so uncertainty can move sure cells into the
  conservative boundary band but can never tighten an enclosure.
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
  refuses typed. `convex_separation_with_inflation` subtracts the composed
  retained radius from `lo`, adds it to `hi`, and re-derives the verdict.
  `contained_ball_radius` is an optional strict-overlap capability (default:
  no claim). `convex_overlap_witness` seals the smaller positive ball admitted
  by both maps. `convex_penetration_depth` revalidates that witness against the
  concrete pair, seeds an inner Minkowski octahedron, and performs bounded,
  deterministic EPA-style support expansion. The inner hull's closest-face
  distance raises a monotone penetration lower bound; any retained support
  plane lowers a monotone upper bound. A cap returns the honest open bracket,
  while touching or unbound/mismatched witnesses refuse.
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
  `ccd_candidates_with_inflation` adds each side's retained error to its own
  motion window before BVH construction; larger evidence radii can only add
  candidates.
- `ImplicitGapOracle` → `GapSample`: pointwise SDF-pair contact
  queries over two exact-distance charts (weaker claims refuse at
  construction, reusing `SeparationRequiresExactDistance`). The
  outward-rounded `[sum_lo, sum_hi]` encloses `φ_A(p) + φ_B(p)`;
  `separation_upper` (triangle inequality) exists exactly when the
  point is certified outside both bodies; `overlap_inradius` is a
  certified common-ball radius exactly when `max(φ_A, φ_B)` is
  certified negative — a pointwise witness, never by itself a penetration
  depth. `overlap_witness()` carries the probe and radius in a sealed token;
  EPA must still revalidate the ball against its concrete support maps;
  `normal` is an uncertified Estimate-class contact axis that is
  absent whenever a gradient is honestly declined. `new_with_inflation`
  widens the sum by both radii, tests outside/inside against each side's own
  radius, and shrinks overlap witnesses rather than laundering them.
  Malformed evidence refuses per sample; cancellation is checked after each
  chart call.
- `CodimThickness` / `codim_gap` / `codim_gap_from_separation` →
  `CodimGap` with `CodimVerdict::{ProvenClear, ProvenContact,
  Unresolved}`: shell/rod contact as an outward-rounded effective-gap
  bracket `dist(midsurfaces) − (t_a + t_b)` over a caller-certified
  midsurface distance enclosure. Soundness is directional: clear
  proofs need only a true lower bound (convex-hull separations
  qualify), contact proofs need the upper bound realized between
  actual midsurface points — the separation composition DOWNGRADES
  contact to `Unresolved` when the caller cannot vouch for the
  witnesses. `codim_gap_with_inflation` widens the midsurface-distance
  enclosure before thickness subtraction. A straddling bracket claims
  nothing.
- `DeformationMap` / `DeformedChart`: the deformable-adapter hook —
  a certified pull-back into a reference configuration presents a
  reference exact-distance chart as a CURRENT-configuration
  `LipschitzImplicit` chart (composed bound `1 × L_pullback`, rounded
  up). Sign and zero set transfer; distance does not, and the wrapper
  never claims it. Gradients are declined (no certified Jacobian yet),
  reference enclosures pass through rigorously or collapse to
  no-claim, and a broken pull-back yields a no-claim sample every
  typed consumer refuses. Construction refuses non-exact-distance
  references and unusable Lipschitz bounds.
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

- `extract_bore` -> `BoreExtraction`: the mesh->A(x) tube extractor
  (music ingest lane). Pipeline: `medial_poles` cloud -> deterministic
  thinning -> SPINAL filter (a pole survives only when its radius is
  near-maximal within its own-radius neighborhood, which deletes the
  disc-shaped medial sheets at flat tube ends) -> Euclidean MST ->
  diameter-path spine with a reach-based BRANCH refusal -> radius^2-weighted
  pole centroids binned by segment-projection arc length along the
  smoothed, end-trimmed spine -> per-station polar plane-sections
  (bisected boundary rays; area is the exact polygonal polar quadrature)
  -> cross-checks against `pi*(pole radius)^2` and `thickness_at` -> axial
  `integral A ds` vs the certified `geometric_moments` volume enclosure.
  Open tubes are extended to their physical end faces along averaged
  interior tangents, with end stations INSET a tenth of the local radius
  (a station on the face clips its own section plane) and face-slab
  volumes added back; closed loops (a full torus) are detected by
  spine-endpoint proximity and integrated around the wrap. `BoreStation`
  carries per-station diagnostics (s, area, equivalent radius, pole
  radius, deviations, min/max ray); `BoreExtraction::debug_lines` renders
  them as JSON lines; `VolumeClosure` is `Certified` only for
  `TraceStepClaim::ExactDistance` charts and otherwise RECORDS
  `Unavailable { claim }` - never a downgraded number. Authority of every
  non-closure output is `Estimate` by construction (gb-001..gb-007).

- `extract_crook_delta` -> `ValveChartDelta` + `certify_slide_range`
  -> `SlideRange` (music bead 3ez8g.3.2): per-valve inserted-length
  receipts minted by the bore extractor and CERTIFIED against the
  CAD's declared centerline arc length — outside the authored band the
  extraction REFUSES by name (`DeltaLengthOutOfBand`: a crook that
  measures wrong is a data error, not a detuned horn). Port areas are
  read at near-end INTERIOR stations (the extractor's cut-face end
  stations carry the disclosed tangent-tilt term) and checked against
  the DECLARED junction step: an undeclared mismatch refuses
  (`JunctionAreaMismatch`) while the same geometry with the step
  declared passes — the acoustic-feature-vs-extraction-error
  distinction. Slide ranges certify from two endpoint extractions;
  runtime switching lifts belong to the performance images
  (gk-001..gk-003).

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


The bore extractor's authority invariant: `BoreExtraction.authority` is
always `Estimate`; the only rigorous field is the `Certified` closure
enclosure, inherited from `geometric_moments` under its own exact-distance
gate; a deliberately biased area sweep fails the closure comparison (the
executed gb-002 falsifier).

## Error model

`ContactInflationError` refuses non-rigorous, negative, non-finite, inverted,
QoI-inconsistent, or arithmetically overflowing conversion/router evidence.
`QueryError::InvalidContactInflation` carries that refusal at query boundaries.
Other `QueryError` teaching errors: `NoGradient` (with the location),
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
`ConvexOverlapUnproven` (strict positive overlap was absent, touching-only,
or did not revalidate against the supplied pair),
`FeatureInvalidInflation`, `FeatureTooManyPairs` (cap refusal, never
silent truncation), `CodimInvalidThickness`, `CodimInvalidDistance`,
`DeformationRequiresExactDistance`, `DeformationInvalidMap`,
`Cancelled`,
`Mesh` (fs-mesh refusals carried through). Honest gaps refuse; nothing guesses.
The bore extractor has its own typed `BoreError`: `Query` (any underlying
`QueryError`), `InvalidConfig`, `TooFewPoles`, `BranchedLumen` (with the
attachment point and side-arm length - segment-wise extraction between
authored cut planes is the supported workflow), `DegeneratePoleChain`,
`CenterOutsideLumen`, and `RayNoExit` (with the station center, tangent,
and exhausted cap for log-only diagnosis).

## Determinism class

Fully deterministic: fixed iteration counts, canonical grid orders,
no randomness. Identical inputs give identical answers bitwise. Exact-zero
inflation delegates with the original floating-point inputs unchanged;
positive composition and endpoint changes use fixed outward rounding.
Penetration hull faces, horizon edges, and equal-distance ties use canonical
index order; budgeted runs are deterministic prefixes whose `lo` cannot fall
and whose `hi` cannot rise.

## Cancellation behavior

`raycast` polls before each sample and again after `eval` and
`trace_value_enclosure`; `separation` polls immediately before and after every
producer call and once before publishing; thickness polls around every chart
evaluation and before publishing local or aggregate estimates. Producer-
requested cancellation therefore wins over evidence returned by the same call.
Closest-point and curvature queries likewise poll around every analytic or FD
sample and once before publishing. All return `Cancelled` teaching errors.
Other point queries are O(iterations) and non-blocking.
EPA polls before and after every support evaluation and at a fixed iteration
stride; no penetration bracket is published after cancellation.

## Unsafe boundary

None. `#![forbid(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None.

- `felt` module (music bead 3ez8g.3.4): felt thickness-field chart
  ingestion — the missing field -> chart + receipt side of the piano
  hammer's contact geometry (the sampler side already existed in
  fs-contact). A validated station profile (crowned t(theta) or a
  manufacturer strip t(s)) with an EXPLICIT units declaration
  (undocumented units refuse, never guess), an ingest receipt (FNV-1a
  digest over exact bits, thickness statistics min/max/crown, largest
  uncertainty half-width, authority = Estimate ALWAYS), and
  `FeltThicknessChart`: the hammer solid (core cylinder + felt +
  width-rounding sag, slab-bounded) as a `Chart` with sign-exact
  nominal, `Estimate` certificates whose band is the station
  half-width, and `TraceStepClaim::NoClaim` — the AxisymmetricChart
  authority precedent, so gap sampling runs under AllowEstimate and
  the receipt authority is never promoted. LINEAGE (corpus law): the
  committed manufacturer fixture
  `data/felt/steinway-us5125310-strip-taper.tsv` is TRANSCRIBED (not
  digitized) from the public text of US patent US5125310 (Steinway,
  1992) and records exactly the stated endpoint values; the linear
  interpolation between stations is the module's declared rule.

## Conformance tests

`tests/conformance.rs`, cases gq-001..gq-007 (+ typed trace refusal checks)
— every aggregate verdict is a canonical fs-obs `ConformanceCase` event.
The randomized gq-001, gq-002, gq-003, and gq-004 campaigns carry their
literal LCG input seeds (`0x1001_2026_0707_0011` through
`0x1001_2026_0707_0014`); all fixed-input cases use zero. The suite's fixed
`Cx` stream seed `0x9e4` is execution provenance, not a substitute for an
input seed; mixed LCG/`Cx` cases key the aggregate verdict to the LCG input
and record the fixed stream separately in detail. The typed fs-obs companion
events for the thickness estimator and curvature convergence tables remain
independently wire-validated. Assertions and expectations reached before an
aggregate verdict remain ordinary Rust test diagnostics and emit no verdict.
Any reimplementation must pass the suite unchanged.
`tests/moments.rs`, cases gm-001..gm-007 — certified geometric moments:
box/sphere/torus/hollow-shell closed-form containment with bounded widths,
COM enclosures, the translation-covariance metamorphic (outward-rounded law
vs direct recomputation must overlap), capability/input/work refusals,
cancellation, and the open-mesh capability-routing refusal. The canonical
fs-obs suite/session is `fs-query/moments`. Scopes gm-001, gm-002, gm-003,
gm-005, gm-006, and gm-007 each emit one `ConformanceCase` at sequence zero.
Scope gm-004 retains ten ordered refusal rows at sequences 0..9:
`weaker-claim-refuses`, `estimate-sample-refuses`,
`zero-spacing-refuses`, `nan-spacing-refuses`,
`inverted-domain-refuses`, `non-finite-domain-refuses`,
`support-excluding-domain-refuses`, `per-axis-work-refuses`,
`total-work-refuses`, and `unproven-volume-refuses-com`. All fixtures are
fixed-input, so every aggregate `seed` is zero. That input seed is deliberately
distinct from `Cx` execution provenance: gm-001..gm-004, gm-006, and gm-007
use stream `(seed=0x60e5, kernel=9, tile=0, iteration=0)`, while gm-005 uses
`(seed=0x60e5, kernel=10, tile=0, iteration=0)`. Every event detail records
the applicable stream. The suite requires no feature flag. Assertions and
expectations reached before an aggregate remain ordinary Rust diagnostics;
each gm-004 row is linted, wire-validated, and printed before its terminal
assert, so a refusal failure retains the exact emitted prefix. No Custom
metric/forensic rows are produced. `geometric_moments` is Enclosure-class:
every returned bracket contains the true unit-density integral of the chart's
region.
`tests/convex.rs`, cases gc-001..gc-005 — certified convex separation:
analytic sphere/box distance containment, nonsmooth box-box honesty,
touching/overlap never claiming separation, bit-identical replay, and
constructor/budget/cancellation refusals. Every aggregate outcome uses a
canonical fs-obs `ConformanceCase` with the shared `0xC0F` execution seed;
direct assertion failures before those events remain ordinary Rust test
diagnostics.
`tests/penetration.rs`, cases gp-001..gp-006 — analytic sphere/box depth
containment, sealed gap-witness revalidation and cross-pair refusal, touching
and zero-budget refusal, monotone budget prefixes, bit-identical replay, and
cancellation, plus malformed support/slack refusal.
`tests/features.rs`, cases gf-001..gf-005 — feature complexes and CCD
candidates: canonical construction, BVH-equals-brute-force oracle
agreement, guaranteed inclusion of the closest pair inside the motion
window, translation invariance and window monotonicity metamorphics,
identical replay, and refusal drills. Aggregate outcomes use canonical fs-obs
`ConformanceCase` events; evaluated gf-002..gf-005 carry the shared `0xFEA7`
execution seed and constructor-only gf-001 uses zero. Direct assertion
failures before those events remain ordinary Rust test diagnostics.
`tests/gap.rs`, cases gg-001..gg-004 — implicit gap oracle: analytic
containment and true-distance upper bounding on disjoint spheres,
certified common-ball witnesses on overlap, claim/evidence/point/
cancellation refusals, and no normal claim without gradients. Every
aggregate outcome uses a canonical fs-obs `ConformanceCase` with the shared
`0x6A9` execution seed. Direct assertion failures before those events remain
ordinary Rust test diagnostics.
`tests/codim.rs`, cases gd-001..gd-004 — codimensional gating:
bracket containment with all three verdicts on their exact sides,
shelled-sphere composition against analytic offset-body geometry with
hull-only contact downgrades, thickening monotonicity, and refusal
drills. Aggregate outcomes use canonical fs-obs `ConformanceCase` events;
gd-002 carries its `0xC0D1` execution seed and the three deterministic cases
use zero. Direct assertion failures before those aggregate events remain
ordinary Rust test diagnostics.
`tests/deform.rs`, cases gh-001..gh-004 — deformation hook: sign and
safe-step transfer under scale+shift against the analytic deformed
sphere, rigorous enclosure pass-through with Estimate collapse,
construction/sample refusals, and no distance-claim laundering into
exact-distance consumers. Aggregate outcomes use canonical fs-obs
`ConformanceCase` events; gh-001..gh-003 carry the shared `0xDEF0` execution
seed and constructor-only gh-004 uses zero. Direct assertion failures before
those events remain ordinary Rust test diagnostics.
`tests/inflation.rs` — conversion/router receipt refusal, exact-zero bit
neutrality, outward widening and witness shrinking, and monotonicity across
convex, implicit-gap, codimensional, CCD, and moments consumers.

`tests/bore.rs`, cases gb-001..gb-007 - the bore-extraction battery
(fs-obs suite `fs-query/bore`, fixed-input seed 0, execution seed 0xB0E5,
kernels 30/31): straight F-rep cylinder oracle (area vs pi*r^2 at 7.4e-4,
axial volume vs analytic at 3.7e-4, honest `Unavailable` closure for the
`LipschitzImplicit` composite), full-torus closed-loop with the CERTIFIED
volume closure and the executed biased-sweep falsifier, 90-degree
bent-tube oracle with a face-distance two-tier area gate (tight at two
radii or more from a cut face; disclosed 6% tilt envelope beside it),
cone-frustum taper tracking on an `AxisymmetricChart`, refusal battery
(Y-junction branch, too-few-poles, config), pre-cancelled context, and
the OBJ-cylinder end-to-end through the fs-io quarantine (24-gon section
truth, `NoClaim` closure recorded, cap-less variant refused at
promotion).

`tests/crook.rs`, cases gk-001..gk-003 — the three-valve fixture
(bent-torus lumens with analytic arc lengths; worst CAD deviation
0.0148 inside the 0.10 band; delta ordering certified; slide range +
inverted-endpoint refusal), the band refusals by name (lying CAD
length; undeclared vs declared port step), and THE COMPOSED ORACLE:
inserting the extracted delta into a cylinder duct moves the TMM first
impedance peak -936.1 cents vs the -931.6 quarter-wave-plus-end-
correction prediction (4.5 cents inside the authored 12-cent band) —
the emergent-pitch precondition, executed.

`src/felt.rs` in-module battery, cases gf-001..gf-004 — refusals by
name (non-positive/non-finite thickness, unsorted stations, missing
lineage, units mismatch, bad chart geometry), deterministic digest +
stats (a one-bit thickness change moves the digest), chart signs with
the shoulder clamp and Estimate certificates containing the nominal,
and the committed Steinway strip ingesting with the patent's exact
stated values (1 in -> 1/8 in outer, 1/4 in -> 3/32 in under). The
COMPOSED check lives at L3 (fs-couple felt_hammer_island): ingested
field -> AllowEstimate gap sampling -> convergence-checked response
curve -> FiniteGapPoint, with the constant-thickness control matching
the exact cylinder+width-sag closed form at 1e-8 and grid refinement
at rel 0.016 (band 0.05).

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
  [`GeometricMoments::translated`] covers translation only. The fixed-fixture
  battery makes no randomized-coverage claim, and gm-007 checks an honest
  capability refusal for one open mesh rather than certifying a general
  watertightness detector.
- The overlap side of an implicit-gap query is strictly a pointwise common-ball
  witness: `max(φ_A, φ_B)` certified negative proves that ball and nothing
  else. It authorizes convex penetration only after the sealed token is
  revalidated through both concrete support maps' positive contained-ball
  capability. `ClearanceField`'s positive-part sum remains authority-free and
  can never authorize penetration.
- `DeformedChart` transfers sign, zero set, and a Lipschitz theorem —
  never distances, gradients, or support tightness. Pull-back
  bijectivity (no self-interpenetration) and the declared Lipschitz
  bound are the adapter's certificates; a lying map produces sound
  refusals only where non-finiteness betrays it. Certified Jacobian
  chains and swept/rotating deformations are deferred surfaces.
- `codim_gap` verdicts inherit the caller's midsurface-distance
  certificate: a lying enclosure or an overstated thickness radius
  produces confident nonsense the routine cannot detect. Thickness is
  a single isotropic radius — anisotropic shells, tapered rods, and
  swept/rotating codimensional bodies are deferred surfaces.
- `ImplicitGapOracle` is pointwise: `separation_upper` bounds the
  distance between the two bodies from above at ONE probe, and
  `overlap_inradius` witnesses a common ball at ONE probe. It does not
  aggregate into a global gap field or time-of-impact claim; only the sealed
  pointwise token, revalidated against matching convex maps, can seed the
  separately bounded EPA route. The `normal` estimate carries no certificate.
- `ccd_candidates` is a broad phase only: candidates are a
  conservative SUPERSET under the caller's declared motion bounds, and
  no narrow-phase, time-of-impact, or contact claim is made. The
  motion bound is a caller assertion (a certified radius over the CCD
  window); the routine cannot detect an understated bound. Boxes are
  axis-aligned in the input frame — rotation-aware swept bounds are a
  deferred surface.
- `convex_separation` proves separation only (`lo > 0`); an enclosure
  containing zero claims NOTHING about overlap. `convex_penetration_depth`
  is a separate strictly-admitted route; it does not upgrade a zero-containing
  separation bracket.
  No convergence rate is claimed: nonsmooth pairs may return wide,
  honest brackets at the iteration cap. Support maps must satisfy the
  declared slack contract; the routine cannot detect a lying
  `support_slack()` or `contained_ball_radius()`. Penetration is minimum pure
  translation depth for compact convex sets only: no rigid-pose separating
  motion, nonconvex decomposition, time-of-impact, or response impulse is
  claimed.
- Chart-native fast paths (mesh BVH closest-point dispatch instead of
  generic Newton) are perf-lane work; answers here are correct first.
- The plain query entry points assert that no conversion/motion uncertainty
  exists outside the chart's own retained sample evidence. A caller that
  obtained a chart through a router/converter, or that has a separate motion
  error, must use the corresponding `*_with_inflation` entry point. The
  carrier validates receipts and arithmetic; it cannot detect a caller that
  deliberately pairs a valid receipt with the wrong chart.
- The bore extractor is NOT a medial-axis transform and NOT a certificate:
  centerline, A(s), and both cross-checks are `Estimate` (the underlying
  marcher has no no-tunneling theorem). v1 topology is single unbranched
  tubes (open or closed-loop); branched lumens refuse rather than guess,
  and the segment-wise-between-cut-planes workflow belongs to the caller.
  Station sections within about two local radii of a cut face carry a
  disclosed tangent-tilt inflation (measured 3.3% on the bend oracle);
  the thickness cross-check can skip stations whose diametral chords end
  on a tight support box (`thickness_at` correctly refuses there), which
  is counted, not hidden. The certified volume closure exists only for
  exact-distance charts; imported meshes (`MeshChart`, `NoClaim`) record
  `VolumeClosure::Unavailable` - there is no certified-volume route for a
  mesh anywhere in the workspace, and this module does not invent one.
- The crook certificate is CONSISTENCY against CAD intent, not a
  measurement of a physical horn; it inherits the bore extractor's
  Estimate authority and tube-only topology. Continuous-slide
  sub-sample delay interpolation is the wind epic's fractional-delay
  bead (a dependency of claim, not of code).
- The felt chart is contact GEOMETRY ingest only: the felt
  CONSTITUTIVE card is bead 87zbd's (the licensed-coupon absence row
  records that hunt); felt-as-absorber (JCA) is a different bead and
  claim; a 3-D felt continuum would be an offline authority lab. The
  width-rounding sag is an authored geometric model (real hammers are
  crowned across the width); the strip fixture's linear interpolation
  between the patent's stated endpoints is the ingest rule, not a
  patent statement.
