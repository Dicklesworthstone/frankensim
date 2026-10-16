# Minimum material under response limits

`minimum_material_sdf3` minimizes numerical projected material volume subject
to explicit upper/lower embedded-reaction and displacement-integral limits:

```sh
cargo run -p fs-topopt --features cutfem-marquee --release \
  --example minimum_material_sdf3 -- 80 250000 0.0012 0.003
```

Use the repository's DSR/RCH execution lane where available. Arguments are
accepted updates, cumulative Krylov iterations, minimum tensile reaction and
maximum tensile reaction. The example independently solves extension and
compression, with corresponding signed reaction windows and displacement
ceilings. It prints every original value/bound/residual, the smooth constraint,
accepted optimization history and retained raw/projected densities. All units
are dimensionless reference units; the bounds are design requirements, not
calibrated specimen measurements. A fixed reference mass density makes volume
minimization equivalent to mass minimization; heterogeneous costs are not modeled.

## Inequalities, not soft targets

Use `sdf3::response::limits::{ResponseLimitCase3, ResponseLimit3,
BoundedResponse3, ResponseLimitSense3}`. A displacement is `q^T u`, with `q`
assembled on the current grid by the existing volume/surface integrators. A
reaction uses the existing embedded numerical traction mode. Force and moment
modes follow the original reaction convention: traction exerted on the solid,
not actuator work or reactions on additional strongly clamped box nodes.

Each `AtMost` residual is `(value-bound)/scale`; `AtLeast` reverses that sign.
Zero and negative bounds are valid, but the declared physical scale must be
finite and positive. Supply two limits for a window. Limits are not averaged
or prioritized with weights. All independent experiments must pass.

The optimizer uses the unnormalized shifted log-sum-exp upper approximation:

`G = max(c) + log(sum(exp(beta*(c_i-max(c)))))/beta`.

In exact arithmetic, `max(c) <= G <= max(c)+log(m)/beta`. Subtracting `log(m)`
would allow safe observations to average away a violated limit and is not done.
The finite-arithmetic result is numerical, not an interval certificate; the
original maximum is checked as well as `G` in `feasible(tolerance)`. Tolerance
is dimensionless and corresponds to each limit's supplied scale. Smoothing is
fixed during a session. It can conservatively exclude a tightly constrained
window even when the original inequalities have a feasible point.

## Reused physics and optimizer

`CutDensityStudy3::evaluate_response_limits` shares one current-density
preparation across one primal/aggregate-adjoint pair per experiment. It retains
reaction direct material derivatives, the prescribed-motion load derivative,
and the original full stiffness, protected-region, projection and filter
pullbacks. No iteration differentiation or per-limit dense design Jacobian is
introduced. Reaction traces are reintegrated during the adjoint pass to avoid
retaining a state/material Jacobian for every scalar limit. Callbacks must be
pure, density-independent reference laws. Geometry and nodal observations stay
fixed during optimization; rebuilding them on another grid is explicit.

`limits::design::MinimumMaterialStudy3` connects that evaluator to the existing
linear-storage projected augmented-Lagrangian state. Bound constraints are
projected directly; there is no dense KKT matrix. The objective is volume
fraction divided by a fixed positive `volume_scale`; reported volume is not
rescaled. Split `run` calls retain multiplier, penalty, spectral step and work.

AL accepted steps can violate response limits. The driver keeps the latest
accepted physical state and the least-volume accepted feasible endpoint
separately. `into_best_feasible` consumes the optimizer, installs the latter's
matching material and returns its fields; it errors when no such endpoint
exists. It does not pretend an earlier incumbent has the latest iterate's KKT
residual. Rejected/cancelled probes do not change accepted material and all
spent work remains counted. Outer/linear/setup budgets retain their existing
owners. The example emits retained evidence before propagating a late error.

This is fixed-reference linear elasticity with numerical cut quadrature. It
does not provide continuum guarantees, follower loads, shape derivatives,
global optimality or response-limit-specific DWR/adaptive refitting. The existing
squared-response estimator is not silently reused for this different objective.

## Focused verification

The `response_limits_sdf3` target includes six physical-evaluator regressions
and four minimum-material optimizer regressions; two aggregation unit tests
cover conservative bias and numerical refusal. They cover coordinate finite
differences, reaction/lifting terms, independent cases, protected regions,
accepted-field alignment, split-run replay, cancellation and infeasible stops.
Native Rust execution is pending in the authoring environment: DSR, RCH and
Cargo are absent. Independent NumPy/SciPy derivative/optimization references
do not execute these Rust implementations and are not performance claims.
