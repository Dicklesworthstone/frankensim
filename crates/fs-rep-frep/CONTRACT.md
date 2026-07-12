# fs-rep-frep — CONTRACT

## Purpose and layer

L2 (MORPH). Function-representation charts (plan §7.2): CSG **DAGs**
over implicit primitives with R-function differentiable Booleans, plus
evaluators AUTO-DERIVED from the same DAG — value+gradient (exact chain
rule), Lipschitz bound (per-node composition), and interval range over
a box (per-node inclusion). The function IS the shape: the abstract
region is the sublevel set `{ p : f(p) < 0 }`, and every claim below is
about THAT region — no silent promotion to "exact distance".

## Public types and semantics

- `FrepBuilder` → `Frep`: arena-style DAG construction. Node ids are a
  topological order; SHARING a subexpression is reusing its id.
- Primitives (L = 1): exact `sphere`, `half_space`, `box_prim`, and
  `cylinder` (infinite, +z); `torus` (+z axis) is exact for ring geometry
  (`major > minor`) and otherwise retains only its exact sign/zero set.
  Unbounded supports are
  reported as ±`UNBOUNDED_HALF` boxes; intersections shrink them back.
- Transforms: `translate`, `rotate` (axis-angle Rodrigues; GA motors
  join with fs-ga), `scale` (uniform, SDF-preserving `s·f(p/s)`),
  `offset` (exact sign/zero set; magnitude remains conservative until a reach
  certificate proves exact-distance preservation).
- `boolean(op, style, a, b)` with `BoolOp::{Union, Intersect,
  Difference}` × `BoolStyle::{Hard, Blend{radius}}`. Every op routes
  through ONE smooth/hard min via sign flips (difference is
  `−min(−a, b)`). `Hard` is exact min/max — and its crease derivative
  discontinuity POISONS shape optimization, which is why it exists only
  LABELED, next to the C¹ alternative. `Blend` is the quadratic
  R-function smooth min `min(a,b) − r·h²/4`, `h = max(r−|a−b|, 0)/r`:
  C¹ everywhere, convex weights, radius exposed as a design lever (the
  `(fillet :r 3mm)` IR shape).
- Auto-derived evaluators:
  - `value(p)` / `value_grad(p)` — exact chain rule; `None` gradient
    propagates honestly from medial points; hard creases return the
    selected branch's subgradient (ties to the left operand).
  - `lipschitz()` — primitives 1; rigid/scale/offset preserve; Booleans
    take `max(La, Lb)` (blend weights are convex). Valid EVERYWHERE.
  - `interval(box)` — exact distance ranges for sphere/half-space,
    conservative interval arithmetic elsewhere; rotated boxes covered
    by corner AABBs; Booleans use monotonicity of min/smin.
- Design levers: `params()` enumerates every numeric in the DAG as
  `(ParamId, name, value)`; `set_param` validates like the builder;
  `d_value_d_param` is the Jacobian action (symmetric FD v1 — see
  no-claims).
- `Chart` impl: composed Lipschitz bound in every sample; certificate
  honesty — pure sphere/half-space/cylinder/box/valid-ring-torus chains under
  rigid/uniform-scale transforms stamp `Exact`;
  spindle tori, offsets, and anything with a Boolean stamp `Estimate` (the value is a conservative
  bound: exact SIGN, `|f(p)| ≤ dist(p, ∂Ω)` — exactly the
  sphere-tracing safety contract). The matching typed
  `Chart::trace_step_claim()` is `ExactDistance` for the exact chains and
  `LipschitzImplicit` for composites; the latter certifies safe steps and the
  zero set, not a geometric-distance upper bound. `differentiability()`
  reports C1 only for kink-free DAGs (no hard Booleans, no box edges).

## Invariants

1. G0 containment: `interval(box)` contains `f(p)` for every sampled
   `p` in the box, on random DAGs mixing all node kinds (frep-001).
2. The composed Lipschitz bound is never violated under adversarial
   near/far pair sampling — and it is TIGHT (observed ratios ≈ 1), so
   the certificate is not vacuous (frep-002).
3. R-function blends are C¹ at seams: analytic gradients match
   crease-straddling central differences for union/intersect/difference,
   while the SAME probe exhibits an O(1) discontinuity on the hard
   union; blend weights are a convex pair (frep-003).
4. A DAG with shared subexpressions evaluates BITWISE identically
   (value and gradient) to its expanded-tree rewrite (frep-004).
5. Sphere tracing with the composed field and Lipschitz bound NEVER
   tunnels: zero safety violations against a dense-scan + bisection
   oracle over random DAGs and rays; certificate kinds and C-class
   report as declared (frep-005).
6. Metamorphic algebra: hard idempotence and commutativity BITWISE;
   blend self-union equals dilation by exactly r/4 BITWISE; rotations
   round-trip and dyadic translations are equivariant to 1e-12; the
   radius/offset levers differentiate exactly (−1), and the
   blend-radius lever is exactly 0 outside its zone and −1/4 on the
   seam (frep-006).

## Error model

Construction and lever edits return `FrepError` teaching errors:
`NonPositive` (radii/scales/blend radii), `ZeroVector` (degenerate
directions), `BadNode` (unknown id), `BadParam` (unknown slot) — each
says what to fix. Evaluation itself is total: honest gaps surface as
`None` gradients and `Estimate` certificates, never as panics.

## Determinism class

Fully deterministic: plain `f64` expression evaluation, no
parallelism, no iteration over unordered containers, no time or
randomness. Identical inputs give bitwise-identical values, gradients,
intervals, and supports on a given target (frep-004's bitwise law is
the regression trip-wire).

## Cancellation behavior

Per-query evaluation is O(DAG) and non-blocking; like the fs-geom
fixture charts, `Chart::eval` does not poll `cx.checkpoint()` — tiled
consumers (fs-rep-sdf builds, contouring) poll at their own strides.

## Unsafe boundary

None. `#![forbid(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None.

## Conformance tests

`tests/conformance.rs`, cases frep-001..frep-006 — JSON-line verdicts,
seeded LCG randomness, per-case fs-obs Custom events (interval
tightness, seam gradient stats, ray-safety counts). Any
reimplementation must pass the suite unchanged.

## No-claim boundaries

- The composite field's MAGNITUDE is a one-sided conservative bound,
  not the exact distance; samples say `Estimate` on purpose. Per-query
  rigorous sd enclosures join with fs-ivl. Accordingly a renderer's
  `|f|/L` termination is a normalized-residual hit, not a certified Euclidean
  distance-to-boundary enclosure.
- The local interval kit rounds to nearest, not outward (frep-001
  carries 1e-9 fp slack); unification with fs-ivl's outward-rounded
  types is the tightening path.
- `d_value_d_param` is symmetric finite difference; exact parameter
  adjoints (chain rule through the DAG) join with fs-xform.
- Revolved/extruded fs-cheb profiles ("revolve THIS function") join
  once fs-cheb's profile evaluators land; the node set here is the
  closed-form primitive zoo.
- Shared subexpressions are re-evaluated per reference (correctness
  proven bitwise; memoized evaluation plans are fs-plan's).
- No self-intersection/validity certificates; `topology_hint` is
  `unknown()` (Betti certificates are the sheaf bead's).
