# Adaptive optimization and enrichment preconditioning

The `adaptive_elastic_sdf3` example prepares numerical factors once per evaluated
optimization density, and separately once per enriched goal-load family. Each
preparation is shared by that density's independent loads. Geometry, density
pullbacks, OC acceptance and the DWR estimator retain their original mathematics.
Bare library backends still select identity, without additional setup work.

```sh
# Arguments: rounds, updates, outer iteration budget, enrichment mode,
# then optimization mode (jacobi by default; two-level also supported).
cargo run --locked -p fs-topopt --features cutfem-marquee --release \
  --example adaptive_elastic_sdf3 -- 2 3 250000 two-level jacobi
cargo run --locked -p fs-topopt --features cutfem-marquee --release \
  --example adaptive_elastic_sdf3 -- 2 3 250000 two-level two-level
# Compare enrichment modes with the same optimization preparation:
cargo run --locked -p fs-topopt --features cutfem-marquee --release \
  --example adaptive_elastic_sdf3 -- 2 3 250000 identity jacobi
```

Use the repository's RCH execution lane where required.

## Repeated optimization solves

Construct `fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3`
with `jacobi(fine, max_diagonal_contributions)` or
`two_level(fine, &coarse, AdaptiveSolveOptions3, checkpoint)`, then pass the
result to the existing `CutDensityStudy3::new`. The same evaluator and optimizer
accept this backend; it is not a separate optimization implementation.

The backend owns the fine operator and its geometric interpolation, so a matrix
cannot be substituted just because its dimension happens to match. Only scales
can change. Interpolation is constructed once per grid; the constrained diagonal
and coarse factor are rebuilt at **every evaluated density**, including rejected
trials. They are shared across all loads of that density, never carried forward
as stale factors. The current generic solver clones/transposes the retained CSR
interpolation at preparation; this is not a zero-allocation cache.

The example's optional optimization correction uses a fixed one-cell slab space.
It is a correction space only, not a replacement coarse physics solve or a
claim of sufficient resolution. The coarsest space is not suitable for every
implicit geometry. Failed support, transfer or factor admission refuses rather
than constructing missing geometry or silently falling back. The existing
Jacobi option has no Galerkin setup applications and may be cheaper on small
load families. Preconditioning can change roundoff and subsequent optimizer
trajectories; only the original identity path promises its original arithmetic.

## Enriched goal solves

Programmatic callers select `GoalRefinementOptions3::preconditioner`:
`Identity` retains the old default, `Jacobi` uses the exact condensed diagonal,
and `TwoLevel` adds the geometric Galerkin correction. Both bare and prepared
adaptive studies pass their identical underlying physical operator and accepted
fields to this estimator. Supply explicit diagonal and `TwoLevelBudget` limits.
The default budget admits 384 coarse coordinates; the generic coarse solve is
hard-capped at 512. This is not a recursive multigrid hierarchy.

The coarse matrix is `P^T A_fine P`, not the separately integrated coarse
elasticity operator. Hanging-node cross terms and density-dependent ghost terms
are included in the diagonal. Every preconditioner is fixed during CG; its
borrow prevents a mutable stiffness update while its numerical data are live.

`SolveWork::linear_iterations` counts outer Krylov iterations.
`preconditioner_operator_applications` separately retains all Galerkin setup
applications, including cancelled/failed and rejected-trial work. Setup limits
apply per preparation, separately from the cumulative Krylov cap. Count both
before comparing work; neither measures wall-clock time or total arithmetic.
A refused preparation, interrupted load family or rejected trial cannot replace
accepted scales or fields. The true-residual gate still checks the actual PDE.

No continuum certificate, asymptotic scalability, mesh-independent iteration
count, or allocation-free/real-time claim is made. Native Rust execution is
pending; numerical development references do not substitute for Rust tests.
