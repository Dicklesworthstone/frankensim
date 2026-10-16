# Optimize the live IR with reverse-mode L-BFGS

`ReverseStudy` connects `fs_opt::ReverseProblem` to the existing limited-memory
L-BFGS engine. No objective closure, finite-difference step size, or dense
Jacobian is required. Each objective callback evaluates one shared primal tape
and computes its signed, weighted objective gradient in one reverse sweep.

```rust,ignore
use fs_ascent::{ReverseStudy, StopRule};
use fs_opt::ReverseProblem;
use fs_opt::reverse::ReverseLimits;

let oracle = ReverseProblem::new(&problem, ReverseLimits {
    max_nodes: 100_000,
    max_scalar_slots: 1_000_000,
})?;
let mut study = ReverseStudy::new(&oracle, &packed_initial_point, 17, None)?;
let report = study.run(&StopRule::GradNorm(1e-7), 200, None)?;
let solution = &study.optimizer().x;
let checkpoint = study.clone();
// Continue without re-evaluating the initial point or discarding curvature.
let next = study.run(&StopRule::GradNorm(1e-9), 100, None)?;
```

Use `Some(&cx)` instead of `None` to observe an explicit `fs_exec::Cx` cancellation
context. The compiled oracle is borrowed and shared across checkpoints. Resume
does not accept a replacement problem, so cached values cannot accidentally be
used with different objective weights, expressions, or variable meanings. Work
counts belong to each study/checkpoint; cloning does not install a shared
cross-branch budget pool.

## Budgets and returned state

The sealed problem's `EvalLimit` is a cumulative **objective-callback** ceiling,
including the initial evaluation and all rejected or failed trial attempts.
Every callback performs at most one primal and one reverse sweep. Gradient
cost therefore does not add two objective evaluations per decision coordinate.
Every `StopRule::Budget` leaf is also a hard ceiling, even inside `All`.

The limit is checked before each line-search probe, including zoom. A trial
that uses the final evaluation may still be accepted; the report then returns
`Budget` with that valid accepted point. A budget exhausted on a rejected trial
returns the previous accepted point. Resuming under an exhausted problem budget
performs no further evaluations. A zero-iteration run only inspects cached state.
The explicit positive/unlimited problem limit always permits the one constructor
evaluation; a tighter run-time rule cannot retroactively undo that evaluation.

The report exposes the existing stop reason, objective, gradient norm and
cumulative work. `optimizer()` is read-only and also exposes the point, gradient,
accepted history and iteration count. Objective weights and minimization versus
maximization senses are taken directly from `ReverseProblem`.

## Domain failures and cancellation

An invalid initial value or derivative is a typed error. During line search,
non-finite primal, adjoint, or scalarization results instead make that trial
unavailable. The search uses the existing positive-infinity barrier convention
to shorten the step. Its zero placeholder gradient is never accepted as a
solution gradient. `rejected_trials()` and `last_rejection()` retain the count
and latest original numerical error. Other failures, including cancellation,
propagate without being recast as domain rejections.

Cancellation is observed at iteration boundaries and inside reverse sweeps.
The accepted point, gradient and curvature remain a valid checkpoint; spent
trial evaluations are not refunded. Cancellation inside a search discards that
unfinished search, so continuation may repeat trials and is not a claim of
identical total work. Splits between complete accepted iterations retain both
trajectory and accounting. Packing, two-loop recursion, and existing binding
validation remain whole non-interruptible phases. There is no process-restart
serialization or bounded wall-clock cancellation-latency claim.

## Scope and runnable example

This driver admits **unconstrained Euclidean variables only**. Constraints are
refused rather than dropped, and Sphere/SO(3)/Stiefel variables are refused rather
than updated as flat coordinates. Those need constrained or Riemannian engine
integration. The lower-level reverse adapter already provides their derivative
products, but this driver does not claim to solve them. Kinks, external physics
and uncertainty nodes, and unexecuted problem tags retain the adapter's refusals.

Existing `Study` and legacy `LbfgsState::new/run` trajectories are unchanged.
The new `LbfgsState::try_new/try_run` APIs are also usable with independent
fallible callbacks. The underlying Wolfe loop is shared, not a second numerical
implementation. Graph limits govern the compiled tape; this is not a certified
bound on total L-BFGS memory, time, roundoff, or global optimization error.

```sh
cargo run -p fs-ascent --example reverse_study
cargo test -p fs-ascent --test lbfgs_fallible --test reverse_study
cargo test -p fs-ascent --lib wolfe::tests
cargo test -p fs-opt --test reverse_problem
cargo test -p fs-ascent --example reverse_ir_optimize
```

The runnable example minimizes a 256-coordinate quadratic with three budgeted
reverse evaluations. The regression suite also covers nonlinear Rosenbrock,
weighted multi-variable objectives, logarithmic-domain recovery, checkpoint
splits, strict budgets, typed errors, and cancellation/resumption with real Cx
values. Rust tests must be executed in the project toolchain before their
success can be claimed; development mathematical references are not substitutes.
