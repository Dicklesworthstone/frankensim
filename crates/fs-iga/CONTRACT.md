# CONTRACT: fs-iga

Isogeometric analysis (1D core): Galerkin directly on B-spline spaces.

## Purpose and layer

Layer L3 (FLUX). No dependencies — pure Rust (B-splines, Gauss quadrature,
dense solve).

## Public types and semantics

- `BsplineSpace::clamped_uniform(degree, elements)` — a clamped uniform space on
  `[0, 1]`; `num_basis` = `degree + elements`; `basis(i, x)` / `basis_deriv(i,
  x)` (Cox–de Boor). Panics on `degree == 0` or `elements == 0`.
- `solve_poisson(&BsplineSpace, g) -> Result<Solution, IgaError>` — solves
  `−u'' = g` on `[0, 1]` with `u(0) = u(1) = 0` (assemble `Kᵢⱼ = ∫ Nᵢ'Nⱼ'`,
  `fᵢ = ∫ Nᵢ g` by 4-point Gauss quadrature per knot span; clamp the two
  boundary DOFs; solve the interior block).
- `Solution` — `coeffs`, `eval(x)`, `l2_error(exact)`.
- `IgaError` — `TooFewDofs` / `Singular`.

## Invariants

- The B-spline basis is a PARTITION OF UNITY (`Σ Nᵢ(x) = 1`).
- Galerkin CONSISTENCY: a polynomial solution representable in the degree-`p`
  space is reproduced EXACTLY (to roundoff).
- k-REFINEMENT: raising the degree sharply reduces the error on a smooth
  non-polynomial solution (fewer DOFs, higher accuracy — the IGA superpower).
- The solution satisfies the homogeneous Dirichlet boundary conditions.

## Error model

Structured `IgaError`; the only panic is a nonsensical space request.

## Determinism class

Fully deterministic: assembly and solve are pure functions of the space + `g`.

## Cancellation behavior

None (synchronous pure functions).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/iga.rs` (7 cases): partition of unity; DOF count = degree + elements; a
polynomial solution reproduced exactly; k-refinement convergence on a smooth
solution; too-small-space rejection; the Dirichlet BCs are satisfied;
determinism.

## No-claim boundaries

- v0 is 1D scalar Poisson on a CLAMPED UNIFORM B-spline space with 4-point Gauss
  quadrature and homogeneous Dirichlet BCs. NURBS weights, Bézier extraction for
  fs-feec element-loop compatibility, multi-patch MORTAR coupling (with sheaf
  interface certificates), KIRCHHOFF–LOVE shells on spline surfaces, reduced
  quadrature, and 2D/3D are staged with this basis/assembly core.
- The dense assemble-and-solve is for the 1D fixture; the production path is
  fs-feec's batched small-dense assembly + the sparse solver stack.
- Inhomogeneous / typed BCs via fs-scenario on patch boundaries are a downstream
  integration.
