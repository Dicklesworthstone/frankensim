# fs-rep-nurbs — CONTRACT

Rational B-spline charts (plan §7.2): EXACT spline algebra, trimmed
patches with certified classification, certified closest-point brackets,
and the HONEST Boolean position.

Ambition tags: exact algebra / trims / closest-point [F per bead label;
what is shipped is tested to [S] discipline]; direct B-rep Booleans
deliberately NOT shipped (see the Boolean position).

## Purpose and layer

Layer **L2** (MORPH). Runtime deps: `std`, fs-ivl, fs-math. Consumers:
fs-iga (geometry basis = analysis basis), fs-render NURBS tracing
(shares the clipping/Newton machinery), the NURBS↔SDF converter beads
(wqd.11/12).

## Public types and semantics

- `Rat` — exact i128 rational scalar (gcd-reduced, overflow-CHECKED:
  leaving the exactness domain is a named panic, never wraparound).
- `Scalar` — the field abstraction; the SAME generic basis/curve/surface
  code runs at `f64` (fast) and `Rat` (exact). The conformance suite
  also instantiates it at a test dual, so derivative checks flow through
  the identical code path.
- `KnotVector` (clamped, validated), Cox–de Boor basis (Piegl–Tiller
  A2.2 shape).
- `NurbsCurve<S, DIM>` — homogeneous de Boor evaluation; f64 derivatives
  to arbitrary order (homogeneous differencing + rational Leibniz);
  EXACT Boehm `insert_knot`; EXACT `remove_knot` (reconstruction checked
  with scalar EQUALITY — in `Rat` a proof); Bézier decomposition; EXACT
  `elevate_degree` (per-segment binomial elevation reassembled on a
  full-multiplicity knot vector — valid and evaluation-identical;
  minimal-multiplicity elevation is a follow-up).
- `NurbsSurface<S>` — tensor-product evaluation, EXACT directional knot
  insertion, first partials via isocurve nets, per-span control boxes.
- `TrimLoop`/`TrimmedPatch` — trim curves in EXACT RATIONAL form
  (closure validated by rational equality). `classify` is CERTIFIED:
  outside every Bézier span hull ⇒ the exactly-computed control-polygon
  winding equals the curve winding (hull-homotopy argument); ambiguous
  points subdivide EXACTLY (rational midpoints) and become `Boundary`
  after the budget — never a guessed in/out.
- `closest_point_curve`/`closest_point_surface` — best-first
  branch-and-bound over rational-Bézier hull boxes (positive weights
  survive de Casteljau, so every sub-segment hull bounds it), one local
  split per iteration, split junctions as free upper-bound samples,
  Newton polish on curves. Returns `CertifiedDistance { lower, upper }`
  with hulls epsilon-inflated (1e−9·scale) against f64 rounding.
- `boolean(op, policy)` — THE BOOLEAN POSITION: always a structured
  `BooleanRefusal` in v0. Default policy routes through SDF (convert →
  implicit CSG → re-fit); `DirectCertificateGated` refuses pending the
  sheaf watertightness certificate (wqd.13). An attempt, never a promise.

- `sdf` module (plan §7.3 edge 3, bead wqd.11; [F], behind the
  `nurbs-sdf` feature until its Gauntlet tier is green): the CERTIFIED
  NURBS → SDF converter. `ShellSdf::distance` brackets the unsigned
  distance by convex-hull branch-and-bound (hull distance is a rigorous
  lower bound because the patch lies inside its control hull; any
  evaluated point is a rigorous upper bound) with a damped Gauss–Newton
  polish that can only TIGHTEN the upper bound — certification never
  depends on Newton converging. Trim classification of the winning
  parameter DOWNGRADES the certificate (infinite upper bound) when the
  closest point is trimmed away or in the boundary band — the edge
  label (certified vs measured-only) is decided per input.
  `ShellSdfChart` presents the field as a Chart: sign from DECLARED
  orientation (du × dv outward), `nurbs-sdf/unsigned` without one;
  1-Lipschitz; enclosure certificates. `generate_tile` is budget-aware
  (P4): refinement fires within two cell diagonals of the surface,
  tight claims apply to genuinely adjacent cells, and achieved widths
  plus branch-and-bound splits are ledgered per tile.

- `refit` module (plan §7.3 edge 4, bead wqd.12; [F], behind the
  `nurbs-refit` feature — which enables `nurbs-sdf`, since the CSG
  acceptance closes the loop through the forward converter and the
  sheaf certificate): the SDF → NURBS RE-FIT, the edge that makes the
  honest Boolean policy work (§7.2: Booleans route through F-rep, then
  re-fit when a spline chart is required). v1 pipeline for star-shaped
  domains: radial-bisection sampling ON THE FIELD ITSELF, tensor-product
  B-spline least squares with discrete thin-plate (control-lattice
  Laplacian) regularization, exact G⁰ seam closure by control-column
  tying (G¹ measured). ERROR HONESTY: the spline→SDF direction is
  PROMOTED to a certificate — `sup |sdf(S)| ≤ max sampled + (L_u+L_v)·h/2`
  with per-direction hodograph Lipschitz bounds `L ≤ max‖ΔC‖·(n−p)` and
  the SDF's 1-Lipschitz assumption — while coverage (SDF→spline) stays a
  measured estimate; the report records which is which. Thin features
  below patch resolution produce STRUCTURED WARNINGS with parameter and
  world locations, never silent smoothing. The patch-density knobs are
  the ErrBudget trade, ledgered.

## Invariants

1. **Refinement exactness (the definitive test)**: knot insertion,
   degree elevation, and surface insertion leave evaluation EXACTLY
   equal (rational `==`) at common parameters; insert+remove recovers
   the control net exactly.
2. **G0 laws in exact arithmetic**: partition of unity sums to exactly
   1; clamped endpoints interpolate exactly.
3. **Derivatives match duals**: the generic evaluator instantiated at a
   dual scalar equals the analytic derivative pipeline.
4. **Trim classification never lies**: adversarial fixtures (nested
   square/diamond/island, near-tangent vertices, sliver holes) classify
   correctly; on-boundary points report `Boundary`.
5. **Closest-point brackets contain the truth**: dense-sampling oracles
   (1e5 curve samples, 300² surface grid) always fall inside
   `[lower, upper]`.

## Error model

`NurbsError`: `Structure`, `Domain`, `Exactness`. `Rat` overflow panics
with a named message (exactness-domain exit — a documented boundary, not
a data path).

## Determinism class

**D0**: exact arithmetic on the `Rat` path; fixed iteration orders
everywhere; best-first ties broken by queue order.

## Cancellation behavior

All loops carry explicit budgets (`max_subdivision`, `max_splits`);
P7 by boundedness.

## Unsafe boundary

Zero `unsafe`.

## Feature flags

All OFF by default per the Ambition-Tag rule:
- `nurbs-sdf` [F] — the NURBS→SDF converter; disabled until its
  Gauntlet tier is green. Gates the `sdf` integration target.
- `nurbs-refit` [F] — the SDF→NURBS refit (implies `nurbs-sdf`; the CSG
  acceptance closes the loop through the forward converter + sheaf
  certificates). Gates the `refit` integration target.

## Conformance tests

`tests/conformance.rs` (JSON verdicts, suite `fs-rep-nurbs/conformance`):
nb-001 exact G0 laws + derivative-vs-dual through the generic evaluator;
nb-002 six random rational curves: insertion/elevation evaluation-exact +
insert/remove lossless round trip; nb-003 adversarial trim battery
(nesting, tangency, slivers, boundary honesty); nb-004 certified curve +
surface brackets vs dense oracles (iteration stats logged); nb-005 the
Boolean policy refusals with teaching routes; nb-006 exact surface
refinement + partials vs central differences.

## No-claim boundaries

- **No direct B-rep Booleans** — by design, not omission. The SDF route
  is the product; the certificate-gated path unlocks with wqd.13.
- **Degree elevation emits full-multiplicity knot vectors** (valid,
  evaluation-identical); minimal-multiplicity reassembly is follow-up.
- **Closest-point rigor is hull-based with an epsilon-inflated f64
  margin**, not interval-arithmetic end-to-end; the fs-ivl interval
  polish tier (Bézier-clipping-seeded interval Newton) is the certified
  upgrade path this bead names for LUMEN's tracing.
- **Greville/Gauss quadrature tables for IGA** land with the fs-iga
  consumer (tfz.9), which owns the quadrature accuracy claims.
- **`Rat` is i128-bounded**: deep repeated refinement can exceed the
  exactness domain; the failure is loud and named.
- Ray intersection shares this machinery but ships with the LUMEN
  chart-backend bead (qfx.2).

## No-claim boundaries (sdf converter)

- Lower bounds come from control-hull boxes; near MEDIAL AXES (many
  equidistant patches) and pole-degenerate parameterizations, bracket
  widths converge slowly with splits — the documented budget trade
  (~1e-3 at 2000 splits/cell, 2.6e-4 at 8000 on the unit-sphere tile).
  Interval-Newton (Krawczyk) contraction of the projection equations is
  the upgrade path when fs-ivl grows 2-D machinery.
- Sign is TRUSTED from the declared B-rep orientation; the winding-style
  fallback for imperfect shells lands with the quarantine/census beads
  (fs-io owns mesh-side honesty).
- Trim downgrades widen the certificate; distance-to-kept-region
  (excluding trimmed area from the B&B itself) is future work — the
  current lower bound remains rigorous for the UNTRIMMED surface.

## No-claim boundaries (refit)

- v1 parameterization is RADIAL (star-shaped domains around the given
  center; the bracket failure is a structured teaching error).
  General-topology segmentation over the dual-contoured mesh (wqd.10)
  is the upgrade path.
- The promotion assumes the input field is 1-Lipschitz with an exact
  zero set within its own certificate (true SDFs; min/max CSG of true
  SDFs is 1-Lipschitz but only a distance LOWER bound away from the
  surface — the promoted bound remains valid one-sidedly).
- Coverage (surface → spline) is sampled at the projection grid: a
  feature no ray hits is invisible — density is the caller's knob, and
  the warning channel reports what the samples DID see.
- G¹ cross-seam continuity is measured, not enforced; pole rows are
  least-squares-collapsed, not pinned.
