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

- `closest_point`: damped Newton projection along the chart gradient;
  the post-projection RESIDUAL is measured and reported, never
  assumed. Charts that honestly decline gradients (mesh charts near
  edges) fall back to a central FD on the signed distance — usable,
  with the residual still carrying the honesty.
- `raycast`: conservative sphere tracing on the chart's certified
  Lipschitz bound (`|φ| ≤ dist` makes stepping by `φ/L` safe); clean
  misses are `None`; grazing rays may exhaust the step budget while
  approaching — incomplete, never unsafe.
- `OffsetChart` / `minkowski_ball`: dilation/erosion as a chart
  wrapper (`φ − r`); the ball case of the Minkowski sum IS the offset
  (bitwise), which is the fillet/clearance workhorse.
- `ClearanceField` (`c(p) = φ_A⁺ + φ_B⁺`) + `separation`: grid
  minimization with descent polish, then a RIGOROUS lower bound from
  the field's Lipschitz constant — the true separation lies in
  `[lower_bound, observed]`. Collision margins as certified fields.
- `thickness_at` / `min_thickness`: the THICKNESS ORACLE —
  inward-normal march + bisection to the opposite wall; per-sample
  failures are SKIPPED AND COUNTED. Values respond smoothly to design
  levers where the walls are smooth (the battery FD-differentiates
  min-thickness through an F-rep neck radius and reads the analytic
  subgradient 2).
- `medial_poles`: interior circumcenters of the Delaunay of a
  boundary sample set, λ-filtered by local sample spacing — the
  medial-axis approximation that cross-checks the oracle (2·pole
  radius ≈ local thickness).
- `curvature`: mean/Gaussian/principal from central stencils on the
  signed distance (shape operator = tangent-restricted Hessian), with
  a PER-CHART ACCURACY CLASS (`CurvatureClass`): `SecondOrder`
  (analytic/F-rep — O(h²), measured), `GridLimited` (C¹ grids — error
  floors at the grid's own interpolation error), `Estimate`
  (exact-distance mesh charts — non-smooth across facets).

## Invariants

1. Closest points agree with analytic truth across all four chart
   families within each chart's OWN certificate (exact/F-rep at fp,
   tiled at its declared bound, mesh at faceting scale), residuals
   are honest, and answers are translation-equivariant (gq-001).
2. Raycasts match analytic hits across chart types; tangent rays
   never tunnel (grazes land on the surface or approach); the CSG
   tracer never claims a hit past a dense oracle (gq-002).
3. Offsets of spheres are exactly spheres of the summed radius;
   erosion shrinks exactly; `minkowski_ball` is BITWISE the offset;
   offset charts remain fully queryable (gq-003).
4. Separation brackets hold across shrinking gaps (truth in
   `[lower_bound, observed]`) and the clearance field dominates the
   separation everywhere (gq-004).
5. The thickness oracle reads the graded slab analytically (1% rel),
   finds the dumbbell neck (2× neck radius, zero skips), agrees with
   the medial-pole cross-check, and differentiates through a design
   lever with the analytic subgradient (gq-005).
6. Curvature converges at measured order ≈2 on SecondOrder charts,
   torus principals hit 1/r and 1/(R+r), classes are documented per
   family, grid-limited charts land within their own scale, and
   curvature scalars are rotation-invariant (gq-006).

## Error model

`QueryError` teaching errors: `NoGradient` (with the location),
`NoLipschitz`, `NotOnBoundary` (with the sd found and the advice to
project first), `NoOppositeWall`, `Cancelled`, `Mesh` (fs-mesh
refusals carried through). Honest gaps refuse; nothing guesses.

## Determinism class

Fully deterministic: fixed iteration counts, canonical grid orders,
no randomness. Identical inputs give identical answers bitwise.

## Cancellation behavior

`separation` polls per grid slab; `min_thickness` polls every 64
samples; both return `Cancelled` teaching errors. Point queries are
O(iterations) and non-blocking.

## Unsafe boundary

None. `#![forbid(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None.

## Conformance tests

`tests/conformance.rs`, cases gq-001..gq-006 (+ refusal spot checks)
— JSON-line verdicts, seeded LCG randomness, fs-obs events for the
thickness oracle and curvature convergence tables. Any
reimplementation must pass the suite unchanged.

## No-claim boundaries

- General Minkowski sums (non-ball structuring elements, max-plus /
  FFT-assisted convolution) are deferred; the exact ball case is the
  v1 surface.
- The medial-axis approximation is pole-based (filtered Delaunay
  circumcenters of boundary samples); full filtered-Voronoi medial
  complexes with angle criteria and stability guarantees are the
  follow-up.
- Thickness subgradients are FD demonstrations; exact adjoints
  through the oracle join the gradient-stack bead.
- Separation bounds use a global Lipschitz constant; local bounds
  (interval arithmetic over cells, fs-ivl) would tighten the slack.
- Curvature on mesh charts is an ESTIMATE class by design; discrete
  curvature operators (cotan/normal-cycle) on the half-edge mesh are
  a separate surface.
- Chart-native fast paths (mesh BVH closest-point dispatch instead of
  generic Newton) are perf-lane work; answers here are correct first.
