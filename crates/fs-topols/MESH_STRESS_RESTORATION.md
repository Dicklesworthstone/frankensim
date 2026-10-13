# Multi-load stress restoration with independent mesh checks

The P2 raw-SDF structural-design lane now combines explicit sampled-stress
restoration with `MultiLoadProjectedOptimizer::advance_one_resolution_polling`.
It uses the existing fixed-material-area search, CutFEM equilibrium, per-load
stress probes, checkpoint and lifetime case-solve budget. There is no substitute
solver, new allowable or hidden relaxation of the requested resolution policy.

## Native study declaration

Start from `examples/marquee/bracket-multi-load-mesh-2d.fsim` and explicitly add
`:stress-restoration-reduction 0.01` inside its `:load-family (independent ...)`
block. Keep `:mesh-check (refinement ...)` and supply the sampled-stress limit
appropriate to the declared screening problem. The existing example's very high
limit is deliberately nonbinding; it does not become a calibrated engineering
allowable merely because restoration is enabled.

Run with the existing command surface:

```sh
cargo run --release -p fs-cli -- study study.fsim study.db
cargo run --release -p fs-cli -- study --resume study-<hash> study.db
cargo run --release -p fs-cli -- report study-<hash> study.db
```

The canonical study retains both immutable policies. Resume uses the same
measured admission rules as live execution, including the full independent load
family and every declared grid. Existing strict studies retain their old rules.

## Admission and interpretation

While the accepted optimization-grid design is overstressed, each infeasible
grid must either reach the original admitted stress bound or strictly reduce its
worst **unweighted** sampled-stress excess by more than the declared fraction.
A zero-objective-weight load still governs stress. A grid already within its
bound cannot lose feasibility. Each load must continue to satisfy the original
observed compliance-resolution and material-area tolerances; stress restoration
is not permission to ignore an underresolved solution.

Compliance may increase during restoration. It is not reported as compliance
improvement, and an accepted partial repair remains explicitly stress-infeasible.
A candidate that reaches the bound on the optimization grid must also reach it
on every finer grid before leaving restoration. Subsequent updates use ordinary
all-grid stress-preserving compliance descent. An initially coarse-feasible but
fine-overstressed baseline remains `mesh-unresolved`; this lane does not invent a
fine-grid stress derivative to repair that different situation.

Every started independent solve, including interrupted and refused fine-grid
work, spends the existing allowance. Publication remains atomic with respect to
accepted geometry, phase and update ordinal. Resume re-solves the retained
endpoints under the separate recovery budget; it does not replay update prefixes.

This is a compliance-direction **heuristic** for sampled-stress restoration. It
can stop with `no-feasible-descent`, `mesh-unresolved` or `budget-exhausted` while
still infeasible. Neither an accepted update nor an exhausted iteration budget
establishes convergence, a continuous stress bound, physical validation or an
optimal design.

## Focused regressions

The topology suite covers real two-grid CutFEM restoration with a governing
zero-weight load, independent re-evaluation, exact checkpoint replay, cancellation
and actual solve accounting. Native-study tests cover explicit canonical input,
phase-aware retained-report admission, fine-grid regressions, premature coarse
feasibility, cancellation/resume and unchanged feasible compliance descent.

```sh
cargo test --release -p fs-topols --lib projected
cargo test --release -p fs-cli --lib study::elasticity
```
