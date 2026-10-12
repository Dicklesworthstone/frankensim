# Matrix-free constrained Newton over the live IR

`reverse_newton::AugmentedNewtonStudy` runs Powell–Hestenes–Rockafellar
multiplier updates around the existing Newton–Krylov kernel. It accepts
Euclidean variables, equality constraints and inequalities expressed as
`c <= 0`. Both objective and constraints must be dimensionless: normalize
physical outputs explicitly in the IR before choosing a penalty.

```rust,ignore
use fs_ascent::reverse_newton::{AugmentedNewtonStudy, AugmentedNewtonOptions,
    AugmentedNewtonBudget};
let mut study = AugmentedNewtonStudy::new(
    &oracle, &initial_point, AugmentedNewtonOptions::default(), Some(&cx))?;
let report = study.run(1e-6, 100, AugmentedNewtonBudget {
    objective_evals: 500,
    hessian_products: 2000,
    multiplier_updates: 50,
}, Some(&cx))?;
```

The complete augmented Hessian is the Lagrangian Hessian **plus**
`Jᵀ mu_active J`, computed by one directional sweep and one reverse pullback.
It never constructs a dense constraint Jacobian or parameter Hessian.
Bound physics uses its own explicit Hessian capability. There is no hidden
finite-difference fallback. The PHR penalty is continuously differentiable;
at an inequality's shifted zero, the inner kernel uses the inactive-side
**generalized** Hessian, not a claimed classical second derivative.

## Stop semantics

The report contains the original objective, original constraints, and all
four numerical KKT residuals. Its candidate multipliers are the projected
`lambda + mu*c`. The retained augmented gradient is the Lagrangian gradient
for that candidate pair. `frozen_multipliers()` exposes the different values
currently defining the inner merit. Only all four KKT residuals strictly
below tolerance produce `Converged`.

Inner steps may be infeasible. A small penalty-gradient norm by itself does
not imply feasible convergence. Objective/Hessian/multiplier budgets, inner
stalls, the call's Newton-iteration cap and the maximum penalty each have
separate stopping reasons. Reaching the penalty limit is not a proof that
the constraints are infeasible. KKT residuals are numerical diagnostics,
not global-optimality or continuum/derivative-error certificates.

## Work, cancellation and continuation

The initial point counts as one objective evaluation. All run ceilings are
absolute cumulative limits. Objective and augmented-Hessian attempts include
failed or discarded searches; each Hessian action includes the extra Jv and
transpose work. Each attempted multiplier retarget is charged separately
before it allocates/reweights the tape. Provider-owned physical solve budgets
remain additional independent limits.

Retargeting reuses the accepted primal tape: it does not run physics again.
A new merit resets radius, history and numerical stop state, because old
model ratios are no longer applicable. Accepted point, cumulative work and
`KrylovOptions` survive. A failed/cancelled retarget keeps the old merit but
not a refund. Cloning retains all state for in-process continuation; provider
caches and budgets remain shared according to their binding contract.

`run(..., additional_iterations = 0, ...)` evaluates only retained KKT state.
Splitting calls at completed Newton boundaries preserves the numerical path
with unchanged tolerances and adequate cumulative allowances. No serialized
process checkpoint, full memory admission or wall-clock bound is claimed.

For an externally controlled multiplier strategy,
`ReverseNewtonStudy::new_augmented` exposes a single **frozen** augmented
subproblem, including the existing preconditioned entry points. Its stop
rules and `snapshot().f` refer to the augmented merit, never KKT convergence.
The original `ReverseNewtonStudy::new` still refuses constrained problems.
