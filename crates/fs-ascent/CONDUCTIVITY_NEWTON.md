# Native thermal Newton design

Run the complete physical example:

```sh
cargo run --release -p fs-ascent --features conduction-design --example conductivity_newton
cargo test --release -p fs-ascent --features conduction-design --test conductivity_newton
```

`ConductivityNewtonModel` binds the existing `ConductivityDesign` to the live
`fs-opt` physics IR. Its scalar is a dimensionless-weighted sum of free-node
temperatures, in kelvin. Declare one `Rn(element_count)` variable with
**dimensionless log-conductivity coordinates** and a `PdeResidual` with the
same study name and temperature dimensions. Use `PhysicsBinding` with
`ReverseProblem::new_with_physics`, then pass a real `Cx` to
`ReverseNewtonStudy::new` and `run` (or `run_preconditioned`).

The example optimizes a normalized temperature-rise objective with a quadratic
log-parameter penalty on a synthetic anisotropic heated cube. It is not a
material-card selector or a continuum-optimal heatsink claim. Objective
normalization is explicit IR algebra, so its complete chain-rule Hessian is
retained. The optional feature uses the already-declared conduction dependency;
the default optimizer dependency graph does not change.

## What is computed

The model maps `rho = exp(x)` and uses the actual assembled P1 thermal operator.
It retains a residual-checked primal/adjoint pair and sparse operator at one
exact-bit-keyed point. A nonzero Hessian action requires two sensitivity solves,
not new primal assembly, a dense design Hessian, or finite differences. Both the
implicit thermal-state derivative and the log-map curvature are included:

```text
g_x    = rho * g_rho
H_x v  = rho * (H_rho (rho * v)) + rho * g_rho * v
```

Fixed nonzero Dirichlet temperatures and fixed Robin conditions participate in
the derivative. Temperature-dependent conductivity, moving meshes, contact
parameters and nonlinear boundary feedback are outside this model's contract.
Physical domain failures remain typed refusals, not convergence.

## Work and continuation

A cache miss reserves two native linear-solve slots; a nonzero Hessian action
reserves two more. These are conservative **reservations**, not reported counts
of successful solves or Krylov iterations. Failed/cancelled admitted work is not
refunded. Read `reserved_solves()` and increase `set_solve_limit()` to retry a
provider-budget refusal. The IR's Hessian/objective budgets remain independent;
all native solves also obey the design's iteration/residual policy.

Repeated products reuse the current point. A rejected trial can evict it, but
an older accepted tape re-prepares its own exact point before applying a Hessian.
Failed cache replacement preserves the preceding linearization. Optimizer
clones share this provider's work budget rather than resetting it. This is an
in-process retained-state contract, not serialized process-restart support.
Cancellation is checked between native Krylov iterations and element tiles;
whole sparse matvec, assembly and allocation phases retain their existing
bounds. There is no claimed global byte limit, wall-clock guarantee, certified
Hessian error, or continuum temperature/error certificate.
