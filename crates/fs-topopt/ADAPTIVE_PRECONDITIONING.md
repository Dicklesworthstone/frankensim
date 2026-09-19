# Adaptive enrichment preconditioning

The `adaptive_elastic_sdf3` example now prepares a bounded geometric two-level
preconditioner once per enriched load family. The optimization, geometry,
density pullbacks and DWR estimator are unchanged. Only enriched linear solves
use the selected preconditioner; ordinary optimization solves still use their
existing identity-preconditioned path.

```sh
cargo run --locked -p fs-topopt --features cutfem-marquee --release \
  --example adaptive_elastic_sdf3 -- 2 3 250000 two-level
# Compare with the same workload and stopping controls:
cargo run --locked -p fs-topopt --features cutfem-marquee --release \
  --example adaptive_elastic_sdf3 -- 2 3 250000 jacobi
cargo run --locked -p fs-topopt --features cutfem-marquee --release \
  --example adaptive_elastic_sdf3 -- 2 3 250000 identity
```

Use the repository's RCH execution lane where required.

Programmatic callers select `GoalRefinementOptions3::preconditioner`:
`Identity` retains the old default, `Jacobi` uses the exact condensed diagonal,
and `TwoLevel` adds the geometric Galerkin correction. Supply explicit diagonal
contribution and `TwoLevelBudget` limits. Coarse-size exhaustion refuses rather
than silently switching methods. The example admits 384 coarse coordinates;
the generic coarse solve is hard-capped at 512. This bounded coarse factor is
not a recursive large-scale multigrid hierarchy.

The coarse matrix is `P^T A_fine P`, NOT the separately integrated coarse
elasticity operator. Hanging-node cross terms and density-dependent ghost
terms are included in the diagonal. All preconditioners remain fixed during
CG; the direct coarse solve has no residual-dependent inner iteration count.
Preparation borrows the fine density state, preventing stale reuse after a
normal mutable material update. Prepare again for a changed design or grid.

`SolveWork::linear_iterations` continues to count outer Krylov iterations.
`preconditioner_operator_applications` separately retains all Galerkin setup
applications, including cancelled/failed work. Its limits come from the setup
policy, not the Krylov cap. Count both before comparing work; neither count is
wall-clock time or total arithmetic. One setup is reused across all loads.

No continuum certificate, asymptotic scalability, mesh-independent iteration
count, or allocation-free/real-time claim is made. Native Rust execution is
pending; numerical development references do not substitute for Rust tests.
