# CONTRACT: fs-xform

> Status: ACTIVE (parameterizations v1). Owns the θ → warp contract and
> the FFD / RBF / level-set-velocity / density levers with exact Jacobian
> actions and composition. Manifold harmonics and elastic metrics are
> feature-gated; neural implicits and F-rep procedural programs remain
> follow-up levers (atlas §7.6). Helmholtz filtering is topo-simp's.

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

- `harmonics` module (plan §7.6, bead wqd.20; [F], behind
  `manifold-harmonics`): the Laplace–Beltrami shape spectrum.
  `cotan_laplacian` retains signed finite-element cotangents, including
  negative weights on obtuse triangles: `uᵀLu` is the assembled
  piecewise-linear Dirichlet energy. `ManifoldBasis::compute` uses
  lumped-mass symmetrization `B = M^{-1/2} L M^{-1/2}`, scales B by a
  Gershgorin bound, and applies fs-la LOBPCG to the dimensionless
  shifted operator with the known component nullspaces removed. One
  constant per connected component is re-admitted analytically, ordered
  by lowest vertex index; uniform inflation is a design direction.
  Up to 64 vertices, the existing dense Jacobi kernel permits the complete
  spectrum. Oversampling obeys the numerical kernel's block and work caps;
  when the requested count itself exceeds them, dense Jacobi is a bounded
  fallback up to 128 vertices. Larger requests must pass LOBPCG admission.
  `k` includes constants and zero requests return an empty basis.
  Eigenvalues and recomputed whitened residuals retain inverse-square
  geometry units; no dimensional shift or eigenvalue cutoff removes
  legitimate low modes after a unit change. Sparse whitened weights are
  prepared once, outside iterative matvecs. Deterministic
  sign/ordering (P2: coordinate j means the same thing every run),
  `displace` (normal displacement Σ θⱼψⱼn̂), `project`/`transfer`
  (M-weighted coefficient transfer across basis refreshes),
  `needs_refresh` (bbox-relative drift criterion),
  `dirichlet_energy` (== λⱼ for M-orthonormal modes).
  Geometry must have finite, representable positive triangle areas,
  cotangents, vertex masses, component masses, and operator scale;
  this trusted constructor panics on inadmissible geometry or counts,
  and refuses a rank-deficient projected search block instead of returning
  fewer than `k` modes. A returned mode count is not a convergence or
  spectral-completeness certificate.

- `elastic` module (plan Bet 6 / §9, bead wqd.27; [F], behind
  `elastic-shapes`, implies `manifold-harmonics`): elastic shape-space
  metrics. Curves: SRV transform (`q = c′/√|c′|`), closed-form L²
  geodesics (`srv_geodesic` — the midpoint bisects the distance),
  reparameterization-invariant `elastic_distance` via slope-constrained
  DP matching, `karcher_mean` with per-iteration shift evidence.
  Surfaces (same connectivity): `deformation_energy`
  (membrane + bending via the cotan machinery), `crossfade_path` /
  `straighten_from` / `straighten_path` (monotone-energy discrete
  geodesics; for affinely related shapes the crossfade IS the fixed
  point). `pullback_metric` — the elastic form through the harmonics
  basis is DIAGONAL with weights λ + λ² — and `within_trust_region`
  (step size = bounded deformation, not bounded coefficients).

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

All OFF by default per the Ambition-Tag rule (the v1 FFD/RBF/velocity/
density levers remain unflagged `[S]`):
- `manifold-harmonics` [F] — manifold-harmonic deformation basis
  (`dep:fs-la`); gates the `harmonics` integration target.
- `elastic-shapes` [F] — elastic shape parameterization (implies
  `manifold-harmonics`); gates the `elastic` integration target.

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

- Neural-implicit, spline-control-point, and F-rep
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

## No-claim boundaries (harmonics)

- The signed cotangent form preserves finite-element Dirichlet energy,
  including on obtuse triangles. It does not impose an intrinsic Delaunay
  triangulation or claim a discrete maximum principle. An open surface
  uses natural Neumann boundary conditions; essential boundary constraints
  require a separate constrained operator. Chordal meshes still approximate
  a smooth surface's intrinsic geometry.
- Meshes only in this tier: the SDF narrow-band (implicit-surface)
  Laplacian variant the plan sketches is follow-up scope — the input is
  deliberately a plain (positions, triangles) surface so no VOLUMETRIC
  body-fitted mesh is ever required, which is the doctrine the bead
  encodes.
- Eigenvalues are numerical Jacobi/LOBPCG estimates with recomputed
  `||M^-1/2 L ψ - λ M^1/2 ψ||₂` residuals, not certified enclosures or a
  convergence guarantee at a fixed iteration count. Repeated eigenspaces
  do not give uniquely identifiable individual eigenvectors.
- `transfer` assumes refresh-in-place (shared vertex set); remeshing
  transfers need a resampling front end.

`tests/harmonics.rs` adds focused G0/G3 checks: exact affine-field energy
on an obtuse triangle (the former clamp produced 3.0 instead of 0.25),
the complete analytic single-triangle spectrum, spectrum and residual
unit scaling by `2^-300` and `2^300`, independent original-pencil residual
recomputation, and distinct constant offsets on disconnected components.

## No-claim boundaries (elastic)

- Curve SRV distances are computed over translations + monotone
  reparameterizations; rotation/scale alignment and CLOSED-curve seam
  search are follow-up scope.
- Surface paths require shared connectivity (refresh-in-place); SRNF
  with remeshing correspondence is the [F] growth path the plan names.
- The DP matching is slope-constrained (1/3..3) on a fixed lattice —
  invariance holds to the documented ~5% grid tolerance, not exactly.
- Trust-region superiority is demonstrated on the seeded fixture
  (median-calibrated, ledgered), not claimed universally.
