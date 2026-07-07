# CONTRACT: fs-xform

> Status: ACTIVE (parameterizations v1). Owns the θ → warp contract and
> the FFD / RBF / level-set-velocity / density levers with exact Jacobian
> actions and composition. Manifold harmonics, neural implicits, and
> F-rep procedural programs are follow-up levers (atlas §7.6); Helmholtz
> filtering is topo-simp's.

## Purpose and layer

Design parameterizations — the levers ASCENT pulls (plan §7.6): a
differentiable map θ → space warp shipped with its Jacobian action
(`δT(x) = ∂T/∂θ·δθ`, the boundary velocity field shape-gradient assembly
consumes) and spatial Jacobian (fold-over monitoring, composition chain
rule). Layer: L2 (MORPH). Runtime deps: `std`, fs-geom.

## Public types and semantics

- `Parameterization` trait — `dof`, `apply`, `jacobian_action`,
  `spatial_jacobian`; the contract law is
  `T(θ+ε·δθ)(x) − T(θ)(x) = ε·jacobian_action + o(ε)` (exact for the
  linear levers here; conformance-tested).
- `FfdLattice` — trivariate Bernstein FFD (Sederberg–Parry): linear in
  control displacements, exact basis-contraction Jacobians, analytic
  spatial gradients; points outside the box pass through (boundary
  continuity is the caller's control-point responsibility — documented).
- `RbfMorph` — Wendland-C2 compactly supported handle morphs: radial
  (frame-equivariant, the G3 law), locality by support radius, smooth at
  handle centers (φ′(0) = 0).
- `VelocityBand` — narrow-band scalar velocity DOFs (trilinear, band-
  masked), Appendix C's `xform.level-set-velocity`; plus `advect_sdf`,
  a first-order Godunov upwind step proving the lever drives a level set.
- `DensityField` — raw SIMP densities with clamp diagnostics naming the
  offending component; perturbation = δθ (identity Jacobian).
- `Composed` — `then ∘ first` with chain-rule Jacobian action through the
  spatial Jacobian and matrix-product spatial Jacobians.
- `detect_foldover` — det(∂T/∂x) ≤ 0 at any probe → structured
  `FoldOver { at, det }` refusal with a suggested fix.

## Invariants

1. Every lever's `jacobian_action` is linear in δθ and consistent with
   `apply` (exact secants for linear levers; second-order FD convergence
   for compositions — both conformance-tested).
2. Compact supports are exact: outside an FFD box / RBF support / band
   mask, displacement and Jacobian action are exactly zero.
3. Composition satisfies the chain rule:
   `δT = δT_then(y) + J_then(y)·δT_first(x)`, `J = J_then·J_first`.
4. Refusals name the defect (DOF counts, component index, fold location).

## Error model

`XformError` (`DofMismatch`, `OutOfBounds`, `FoldOver`) — structured,
teaching, never panics across the boundary.

## Determinism class

Pure straight-line f64 arithmetic on the inputs — bit-deterministic
across runs, thread counts, and ISAs. No RNG, no global state.

## Cancellation behavior

All calls are O(dof) or O(grid) loops with no I/O; `advect_sdf` is one
bounded sweep. `Cx` tiling arrives when these run inside fs-exec kernels
(the levers are pure functions, trivially tileable).

## Unsafe boundary

None. Safe Rust only.

## Feature flags

None. All v1 levers are `[S]` default-path.

## Conformance tests

`tests/conformance.rs`: exact-secant + linearity battery for FFD/RBF
(xf-001); dual-number JVP gate on the FFD warp — the solvers'
gradient-gate discipline (xf-002); RBF rotation equivariance + fold-over
refusal on a violent collapse probed through the compression side
(xf-003); composition chain rule with measured second-order FD
convergence (xf-004); level-set velocity driving a Godunov advection step
that grows a sphere at unit speed (xf-005); density/DOF diagnostics
content (xf-006). Module unit tests cover Bernstein partition-of-unity
and derivative-sum laws, Wendland shape/support, band masking, identity
lattices, and validity refusals.

## No-claim boundaries

- Manifold harmonics, neural-implicit, spline-control-point, and F-rep
  procedural levers: follow-up beads (the trait is their contract).
- FFD boundary continuity is not enforced (documented caller
  responsibility); no automatic C¹ blending.
- `advect_sdf` is first-order upwind on a full grid with a frozen
  boundary layer — WENO, narrow-band storage, and redistancing belong to
  topo-levelset.
- Fold-over detection is sample-based monitoring, not a certificate
  (certified invertibility would need fs-ivl interval Jacobians —
  a natural follow-up).
- SIMP chain rule to physics (filtering, projection) is downstream.
