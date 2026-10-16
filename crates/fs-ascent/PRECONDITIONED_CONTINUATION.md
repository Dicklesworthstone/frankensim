# Preconditioned constrained Newton and continuation

The existing inverse-action kernel is connected to `AugmentedNewtonStudy`
and `NewtonContinuation`. The original unpreconditioned entrypoints retain their
numerical paths. It introduces no linear solver, derivative approximation or
external dependency.

For a single constrained model:

```rust,ignore
let result = study.run_preconditioned(
    &mut inverse, tolerance, additional_iterations,
    optimizer_budget, absolute_inverse_limit, Some(&cx),
)?;
```

The inverse receives `(accepted_point, residual, frozen_merit, constraints)`.
`frozen_merit.penalty()` and `.multipliers()` are the actual current PHR model,
not the candidate multipliers in a KKT report. Constraint values belong to the
same accepted point and use declaration order. This permits rebuilding penalty-
and active-set-aware inverse maps without repeating a primal physics solve.
The inverse must be fixed, linear and positive definite during each Krylov
solve. It may change between outer iterations. Positive-product checks are not
an SPD proof, and callback-owned caches are not cloned with the study.

For multiple immutable physical models:

```rust,ignore
let result = continuation.run_preconditioned(
    &mut stage_inverse, additional_iterations,
    optimizer_budget, absolute_inverse_limit, Some(&cx),
)?;
```

The callback adds a leading zero-based stage index. It receives the same frozen
merit and accepted constraints. A cache keyed only to coordinates is insufficient:
two stages can share coordinates but have different physics and curvature.

`preconditioner_evals` in the study report is cumulative over its multiplier
phases. `ContinuationWork::preconditioner_evals` includes retired stages and the
currently active one. The allowance is absolute across all runs and stages;
initializing another stage does not grant a fresh inverse allowance. Exhaustion
returns `PreconditionerBudget` (wrapped in `StageStopped` for continuation).
Failures retain attempted work, including a caught inverse panic. An increased
allowance or repaired callback permits retry from the accepted checkpoint.
Switching temporarily to the unpreconditioned entrypoint does not refund work.

The original objective/Hessian/multiplier ceilings, provider-owned native solve
budgets, trust metric and four KKT tests remain independent and unchanged.
Cancellation uses the existing Cx-aware kernel. Callback side effects are not
undone; whole-process memory and wall-time limits are not introduced.

Native verification commands:

```sh
cargo test --release -p fs-ascent --test augmented_preconditioner \
  --test augmented_newton --test newton_continuation
cargo run --release -p fs-ascent --example preconditioned_continuation
```

The example is explicitly synthetic diagonal algebra, not a measured PDE speedup.
Eight new regressions cover exact inverse acceleration, current multiplier and
constraint metadata, absolute ceilings, faults, cancellation, multi-stage
continuation and split replay. They have NOT been compiled or run in the
local authoring environment, which has no DSR, RCH, Cargo or rustc.
