# Repeated combined-load refinement and minimum-material refitting

```sh
cargo run --release -p fs-topopt --features cutfem-marquee \
  --example worst_load_sdf3 -- --load-space --adapt-material 1500 200 250000 0.02 2
```

The arguments are an authored physical compliance cap, accepted optimizer
updates per grid (default 200), cumulative Krylov iterations (250000), the
normalized full-Gram change tolerance (0.02), and the maximum successful grid
refits (2). The cap is an illustrative dimensionless requirement, not a
baseline-derived rating or a promise of feasibility. Existing fixed-grid,
binary-material and one-refit commands are unchanged.

This command repeats the existing complete-load-set assessment and constrained
refit. It never substitutes separate basis-load checks for the worst combination
`f = F q, ||q||_2 <= 1`, nor reuses a frozen adverse direction. Diagnostic grids
inherit physical eroded stiffness; new optimization grids inherit raw controls,
rebuild their filters and loads, and restore the same requirement. The cap,
load laws, material model, smoothing and filter radius stay fixed.

Every active source cell is refined. Basis-wise DWR masses do not localize the
whole load-space matrix error, so this route does not promote those indicators
into selective full-set error estimates. Tree depth (4) and background leaf
count (4096) remain explicit limits. A large requested refit count does not
remove them. All geometry builds share one quadrature allowance, including
failures. All restoration, gradient, optimization, recheck and assessment solves
share the original physical work context.

## Resumable library owner

`load_space::design::limits::campaign::LoadSpaceCapCampaign3` consumes a completed
`LoadSpaceCapRefit3`, its background, the original reference laws and work
context, a caller-owned geometry builder and `LoadSpaceCapCampaignOptions3`.
The builder must preserve physical identity, poll its supplied checkpoint and
retain a cumulative geometry budget. The campaign owns its current design and
immutable policies across run segments; `run(n)` permits at most n promotions.
`run(0)` does no work. A positive segment assesses the final installed grid even
at the lifetime refit cap, and completed probes are not repeated on resume.

Successful candidates remain pending until the publication checkpoint. A stop
there preserves the candidate, so resuming does not rebuild or reoptimize it.
A stop after publication retains the entire new tree/operator/field family.
Failed candidates remain separately inspectable, never installed, and consume
a lifetime candidate attempt. Physical work is never refunded. Old complete
operators are released on promotion; compact grid/probe summaries retain the
history. `current()`, `current_probe()`, `pending()` and `rejected()` distinguish
these states, and `into_current()` returns the installed tree and its owner.
This is an in-memory session, not a disk checkpoint or an AL-iteration resume.

Only `StableAndFeasible` passes: the current material must satisfy the physical
and smooth combined-load requirements on both grids and the norm of the entire
normalized Gram-matrix difference must meet the authored tolerance. Refinement
limits, attempt limits, candidate stops and pauses are not accuracy passes.
The command exports installed fields before propagating later failure. Its
initial and intermediate outputs remain grid-labeled; no rejected candidate is
presented as the current design.

The public `SolveControl::checkpoint` permits command clients to use the same
cancellation callback between geometry and solve phases. It does not charge
iterations, reset accounting, or independently prove a budget is available.

## Validation boundary

This adds orchestration, not another elasticity solver, eigensolver, error
estimator or optimizer. It provides discrete two-grid evidence, not continuum
error bounds, dilated-volume error estimates, stress/buckling qualification,
zero-void certification, manufactured feature sizes or cross-grid descent.
The new Rust tests and command have not been compiled or run in the editing
environment, which lacks Cargo/rustc.

```sh
cargo test --release -p fs-topopt --features cutfem-marquee --lib \
  sdf3::robust::envelope::load_space::design::limits::campaign
cargo test --release -p fs-topopt --features cutfem-marquee --example worst_load_sdf3 \
  load_space::refinement::adaptive::tests
cargo test --release -p fs-topopt --lib control::tests::phase_checkpoints
cargo test --release -p fs-topopt --doc
```
