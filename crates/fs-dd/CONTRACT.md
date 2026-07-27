# fs-dd — CONTRACT

Domain decomposition (plan §8.9, bead tfz.11): BDDC substructuring with
corners-primal coarse spaces, nodal rho-scaled interface averaging,
sheaf-derived edge enrichment, and CCD-aligned partitioning metrics.

## Purpose and layer
Layer L3. The extreme-core-count solver path; single-machine p-MG covers
the committed scale (the plan's staging note).

## Public types and semantics
- `Decomposition`: the structured 2-D model problem — `s × s` square
  subdomains of `m × m` cells, variable per-cell coefficients
  (`uniform`, `checkerboard` fixtures), 5-point stencil with averaged
  edge coefficients, Dirichlet rim. `apply_global` is the whole-system
  oracle.
- `Bddc::new(decomp, with_edge_modes)`: substructures, factors interior
  blocks, local free-boundary Schur blocks, and the coarse Galerkin
  matrix ONCE; `with_edge_modes` adds one average constraint per open
  interface (the sheaf-derived enrichment) and deflates the local
  corrections against them so the two spaces do not fight.
- `schur_apply` (never assembled), `precondition` (nodal rho-scaled local
  corrections plus an unweighted continuous coarse correction),
  `solve_cg -> Result<CgReport, CgError>` (preconditioned CG with explicit
  `Converged` / `IterationLimit` / `Breakdown` termination, last finite
  relative residual, and an optional `LanczosConditionEstimate`),
  `ccd_locality` (the topological locality metric for island-aligned
  partitions), `gamma_len`/`coarse_dim`. At each interface node, the local
  rho-scaling numerator is the largest coefficient among that subdomain's
  incident cells; a max-normalized denominator forms a finite partition of
  unity over all touching subdomains.
- `LanczosConditionEstimate { kappa, ritz_min, ritz_max, krylov_dimension }`
  describes the extreme Ritz values of the finite RHS-dependent Krylov
  projection. One accepted step is transparently a one-dimensional `κ=1`
  projection; zero steps have no spectral evidence.

## Invariants
- The preconditioner is SPD and symmetric on uniform and checkerboard-1e6
  probes, with and without the sheaf-edge coarse modes (dd-001).
- After explicit convergence, condition estimates track the BDDC signature:
  `κ/(1+log(H/h))²` stays in a bounded band across H/h ∈ {4, 8, 16} (dd-002).
- On the deterministic seed-5, tolerance-1e-8 checkerboard sweep, the
  RHS-dependent Ritz estimate stays in 2.63..3.71 across contrasts
  1..1e6 (a 1.42× band, guarded by a <4× limit). The 1e6 case takes
  12 iterations versus 14 for uniform, preserving the existing ≤3×
  gate. This is evidence for subdomain-aligned jumps, not arbitrary
  coefficient geometry (dd-003).
- The MEASURED seed-11, tolerance-1e-8 coarse-space comparison (dd-004)
  is: uniform corners 19 iterations / estimate 15.58 / dimension 9,
  sheaf 14 / 3.74 / 33; checkerboard-1e6 corners 8 / 7.64 / 9,
  sheaf 11 / 2.63 / 33. The sheaf enrichment strictly improves the
  RHS-dependent estimate on the adversarial jump fixture, but the jump
  iteration count rises from 8 to 11 and the coarse dimension rises from
  9 to 33. The test records that trade without claiming iteration or
  wall-clock dominance.
- Nodal rho weights form a partition of unity. Uniform coefficients reduce
  exactly to multiplicity scaling; a spatial falsifier proves a nonincident
  high-coefficient cell cannot influence another interface node's weight.
- Deflation lesson (recorded): NAIVE edge enrichment DEGRADES BDDC —
  local corrections must be projected away from the coarse edge space.

## Error model
Construction panics on non-SPD local blocks (a mesh/coefficient bug, not a
runtime condition). `CgError` refuses a dimension mismatch, the first
non-finite RHS component, and non-finite/negative tolerance before any solve
shortcut. A zero iteration budget is valid: a zero RHS or tolerance-admitted
initial guess is `Converged`, while any other RHS is `IterationLimit` with
relative residual 1 and no spectral estimate. A zero tolerance requires the
verified true residual norm itself to be exactly zero; division rounding cannot
manufacture exact convergence.

Finite nonzero RHS values are max-normalized before dot products. CG checks
preconditioned residual products, Schur curvature, step/recurrence scalars, and
solution/residual/direction components before use. Non-positive SPD quantities
or non-finite arithmetic produce a `Breakdown` report retaining the last finite
relative residual and only a condition estimate from a completely valid
Lanczos prefix. Convergence and iteration-limit reports recompute the true
normalized residual `b - Sx`; a false recursive convergence restarts from that
verified residual and resets the Lanczos chain.

Lanczos reconstruction refuses empty/mismatched, non-finite, non-positive, or
overflowing coefficients and non-positive/non-finite Ritz bounds. It never
clamps malformed evidence into `κ=1`.

## Determinism class
Fully deterministic: fixed assembly order, dense factorizations, fixed-seed
test probes, max-normalized rho sums in row-major subdomain order, normalized
CG recurrences, true-residual checks, and scaled Sturm bisection. Subdomain
iteration order is row-major.

## Cancellation behavior
Each subdomain factorization/solve is an independent unit (the scoped
resumability contract); the fs-exec tile integration lands when this
crate leaves its feature gate.

## Unsafe boundary
No `unsafe` anywhere in this crate.

## Feature flags
`bddc` ([M], default OFF) gates everything; `sheaf-coarse` (implies
`bddc`, pulls fs-geom) gates the Bet-11 cross-checks.

## Conformance tests
tests/bddc.rs — CG G0 admission/zero-budget/one-step-exhaustion diagnostics;
G3 RHS sign and power-of-two scaling equivalence; dd-001 G0 preconditioner
properties with uniform/jump symmetry and convergence-gated solve claims;
dd-002 log²(H/h) scaling; dd-003 the seven-point subdomain-aligned contrast
sweep; dd-004 the measured sheaf-vs-corners comparison; dd-005 the sheaf
cross-check (the 2×2 subdomain adjacency is a 4-cycle whose 1-D harmonic space
is exactly the mode the corner constraint pins — the sheaf explains WHY
corners are primal) + CCD-locality metrics. Private unit tests verify the
rho-scaling partition of unity and spatial locality, reject malformed Lanczos
histories, and recover a known 2×2 Ritz ratio.

## No-claim boundaries
- The model problem is the structured 2-D 5-point Laplacian; unstructured
  meshes, elasticity, and FETI-DP variants are follow-up scope (the
  bead names FETI-DP as infrastructure-sharing, not shipped here).
- The demonstrated coefficient-contrast independence is limited to the
  subdomain-aligned checkerboard family. Nodal rho-max scaling is defined for
  arbitrary positive per-cell fields, but no robustness claim is made for
  within-subdomain channels or other adversarial coefficient geometry.
  Face-block deluxe scaling and adaptive constraints are follow-up work.
- The sheaf-harmonic coarse space is the EDGE-AVERAGE enrichment framed
  and cross-checked by Bet 11's machinery; adaptive eigenmode
  enrichment (adaptive-BDDC class) is the [M] growth path.
- `ccd_locality` is a topological metric, NOT a wall-clock claim — no
  performance numbers without benchmarks (AGENTS perf rules).
- Condition diagnostics are extreme-Ritz ratios of an RHS-dependent Krylov
  projection, not certified enclosures and not necessarily the full operator's
  condition number. G0/G2 condition claims require `Converged`; an
  `IterationLimit`, `Breakdown`, or one-dimensional `κ=1` projection is not
  affirmative conditioning evidence.
