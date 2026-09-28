# Transient cooling with temperature-dependent conductivity

The existing `cooling-network` command can use the declared scalar `k(T)`
material tables during a transient, rather than freezing conductivity at the
previous temperature. Enable the nonlinear endpoint solver explicitly:

```bash
frankensim --json cooling-network \
  examples/cooling-network/nonlinear-contact-pulse.json
```

The example combines two temperature-dependent solids, different heat
capacities, a finite-resistance matching contact, localized component heating,
a bypass/mixing air network, and a two-speed fan schedule. Its data are
illustrative declarations, not measured hardware or validated materials.

## Explicit endpoint policy

Keep the existing `solid.materials[*].conductivity_curve` declaration and add
this object inside `transient`:

```json
"nonlinear": {
  "max_iterations": 32,
  "residual_rtol": 1e-10,
  "residual_atol_j": 1e-10,
  "armijo_c": 1e-4,
  "shrink": 0.5,
  "max_backtracks": 24
}
```

Every field is required. A k(T) transient without this policy refuses before
hydraulic/thermal execution. Ordinary scalar/tensor constant-material requests
without the policy retain the existing linear timestep path. An explicit
policy also admits constant materials for comparison. A sampled table with
constant ordinates still carries a temperature-validity interval; it is not
silently converted into an unbounded scalar material.

Each endpoint solves, on the free degrees of freedom,

```text
F(T) = C (T - T_old) + dt [A(T) T - b] = 0
J(T) = C + dt J_steady(T)
```

`C` is the existing consistent P1 capacity matrix. `J_steady` is the existing
production Jacobian including K'(T) and the contact operator. Newton directions
use FGMRES with the SPD Picard block as preconditioner. Trial residuals always
re-evaluate the actual material law. Leaving a table's temperature range causes
backtracking, not extrapolation; an inadmissible initial state refuses.

The stopping threshold is `residual_atol_j + residual_rtol * norm(F(T_old))`.
It is fixed during that endpoint solve and denominated in joules. The large
absolute-temperature load `C*T_old` is NOT used as a permissive residual scale.
Neither small updates, exhausted budgets nor a global energy balance alone
can replace this free-residual condition.

`max_iterations` limits Newton updates per solid callback. The existing
`budgets.linear_iterations` caps TOTAL inner Krylov iterations across all
Newton corrections of that callback; it does not reset per correction.
`max_backtracks` limits rejected trials per update. Assembly, Newton updates,
backtracking and bounded Krylov cycles check the existing cancellation context.
No failed endpoint is returned as an accepted temperature field.

## Coupling, adaptive steps and sizing

Every Robin-reference trial retains the same immutable `T_old`. Only after
solid residual, solid energy, interface heat and coupled energy checks pass
does a temperature field become physical history. Constant prescribed
Dirichlet temperatures include their consistent-capacity reactions in the
energy balance; instantaneous prescribed-temperature jumps remain unsupported.

Adaptive coarse/half-step trials, repeated duty cycles and fan/workload-sizing
candidates use the same nonlinear policy. Discarded trials count as solver
work but never as accumulated heat or accepted history. `transient.nonlinear`
in the output records the policy, solid calls, Newton updates, Krylov iterations,
backtracks and worst accepted residual/threshold ratio for that reported cycle.
Outer repeated-cycle/design summaries retain their existing total-work fields.

The independent dense FEM calculation for the example's 2-second steps gives
an illustrative sampled peak of approximately 306.16510336 K at 30 seconds,
and a final peak of 301.41523251 K at 150 seconds. The constant-20/2 W/(m K)
control gives approximately 306.34315850 K and 301.46723363 K respectively.
These are independent numerical reference values, NOT recorded executions of
the Rust implementation. FrankenSim itself (measured 2026-09-25, release build, after 05db922bf made the transient capacitance row-sum lumped; the reference above uses the consistent P1 mass, so the two differ on this 12-tet mesh) gives 304.0634748 K at 30 seconds and a final peak
of 301.6885607 K, and 304.1166334 K and 301.7305465 K for the constant-conductivity
control. The actual-binary regressions now pin these measured values.

## Focused checks and remaining boundaries

```bash
cargo test -p fs-conduction --test nonlinear_backward_euler
cargo test -p fs-cli --test cooling_nonlinear_transient
```

The kernel tests include a hand-derived nonuniform tetrahedral endpoint,
constant-material comparison, heterogeneous capacities, Dirichlet reactions,
immutable history, validity/budget/cancellation refusals, and an intentionally
loose nonlinear tolerance that must still fail the energy gate. The CLI tests
exercise real cooling commands, contacts, replay and adaptive trial accounting.
Their presence does not assert successful execution on a particular checkout.

Heat capacity and contact resistance remain temperature independent. There is
no enthalpy/phase-change model, fluid storage, temperature-dependent air
properties, transient adjoint, continuous-time peak bound or experimental
validation. The existing UQ command still accepts steady base requests only.
The linear theta-method APIs retain their existing constant-conductivity scope;
this addition is the nonlinear backward-Euler endpoint and its cooling consumer.
