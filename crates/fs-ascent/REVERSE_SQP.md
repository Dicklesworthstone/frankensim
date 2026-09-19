# Constrained reverse-mode studies

`ReverseSqpStudy` solves smooth, small-dense Euclidean problems directly from
`fs_opt::ReverseProblem`. Unlike `ReverseStudy`, it retains and enforces both
`EqZero` and `LeZero` constraints. Users do not supply callbacks or coordinate
finite-difference steps.

```rust,ignore
use fs_ascent::ReverseSqpStudy;
use fs_opt::ReverseProblem;
use fs_opt::reverse::ReverseLimits;

let oracle = ReverseProblem::new(&problem, ReverseLimits {
    max_nodes: 10_000,
    max_scalar_slots: 100_000,
})?;
let mut study = ReverseSqpStudy::new(&oracle, &initial_point, 128, Some(&cx))?;
let report = study.run(1e-7, 80, Some(&cx))?;
let checkpoint = study.clone();
let point = study.optimizer().point();
```

The `128` above is an explicit upper bound on the sum of decision coordinates
and declared constraints, checked before the initial evaluation or dense
Jacobian allocation. It limits this dense path, not total process memory.
The reverse compiler's limits separately bound its tape. Dense SQP retains an
n-by-n BFGS model and explicit constraint Jacobians; this is not a sparse,
matrix-free constrained optimizer.

## Values, derivatives and multipliers

One sample attempt evaluates all objective and constraint roots on one shared
primal tape. The weighted objective gradient uses one reverse sweep. One further
reverse sweep per declared constraint forms the Jacobian rows required by the
existing dense SQP kernels. There are no separate primal evaluations to compute
those rows, update curvature, or report the KKT residual.

Objective weights and minimize/maximize senses come from the problem. Negative
satisfied-inequality residuals are retained, not clipped. Equality and inequality
rows are grouped internally for SQP. `report.constraint_multipliers` maps them
back to original declaration order and can be passed directly to the oracle's
`lagrangian_gradient`; `report.solution.lambda` and `nu` remain grouped.
The KKT residual is a local numerical diagnostic, not an interval certificate
or a proof of global optimality or nonlinear infeasibility.

## Budget, continuation and errors

The immutable problem's evaluation limit counts complete sample *attempts*,
including initialization, rejected trials and callback failures. Each attempt
uses at most one primal plus `1 + constraint_count` reverse sweeps. Evaluation
budgets do not count QP pivots or individual reverse sweeps. The existing QP
has its own deterministic pivot bound. A cheaper residual-only line-search
path and a separate reverse-work budget are not implemented here.

`run_with_budget(tolerance, additional_steps, cumulative_cap, cx)` permits a
stricter external ceiling but cannot raise the problem limit. It checks before
every trial, including backtracking. The final permitted trial can be accepted.
An exhausted budget takes stop-attribution priority over convergence. Reports
use the retained point and duals without spending another graph evaluation.

Clone retains the accepted sample, BFGS model, penalty, multipliers, history and
work counts. The shared borrowed oracle fixes problem meaning: resume cannot
silently swap variables, objectives, weights or constraints. Work counts are
per checkpoint, not a shared budget pool across forks.

Non-finite initial arithmetic is a typed failure. Non-finite trial arithmetic
is an explicit unavailable sample, so the merit search can shorten its step;
`domain_rejections()` and `last_rejection()` expose those failures. Other errors
propagate. No placeholder derivative is ever accepted. Cancellation preserves
the accepted checkpoint without refunding attempted work. Partial searches are
not serialized; resume may repeat their probes. Splits at complete accepted
steps retain the numerical trajectory and accounting. Cx polls bracket dense
QP/BFGS phases and run inside existing reverse sweeps; dense kernels and Packing
are still whole phases, with no bounded wall-clock cancellation guarantee.

Non-Euclidean variables are refused. Kinks, unavailable physics/UQ execution and
unsupported problem tags retain `ReverseProblem`'s refusals. The legacy `sqp`,
`Study`, and `ReverseStudy` paths are unchanged.

## Runnable sizing study and focused checks

```sh
cargo run -p fs-ascent --example reverse_sqp_beam
cargo run -p fs-ascent --example reverse_sqp_beam -- 10000000
cargo test -p fs-ascent --test sqp_state --test reverse_sqp
cargo test -p fs-ascent --example reverse_sqp_beam
cargo test -p fs-ascent --test sqp_globalization --test constrained_battery
```

The example minimizes the mass of a rectangular cantilever with a fixed 2:1
height/width aspect ratio, stress and deflection limits, and lower size bounds.
It declares load, length, elastic modulus, density and reference dimensions in
SI, then solves explicitly dimensionless ratios. Changing allowable stress
switches the controlling constraint. The two example tests compare the designs
with closed-form solutions of this same algebraic model. This is illustrative
linear beam algebra, not a mesh-based physics solve or a validated engineering
design. Rust build/tests must run in the project toolchain; mathematical
reference translations are not native execution evidence.
