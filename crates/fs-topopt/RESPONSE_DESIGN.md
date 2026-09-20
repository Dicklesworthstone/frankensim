# Prescribed-motion response design

`motion_response_sdf3` fits two independently imposed boundary-motion experiments
with explicit material-volume and raw-density constraints. It uses the existing
small-dense SQP engine, not the fixed-force compliance OC update.

```sh
cargo run -p fs-topopt --features cutfem-marquee --release \
  --example motion_response_sdf3 -- 40 250000
# Supply four target integral values instead of generating synthetic targets:
cargo run -p fs-topopt --features cutfem-marquee --release \
  --example motion_response_sdf3 -- 40 250000 T00 T01 T10 T11
```

Use DSR/RCH where available. The first arguments limit accepted SQP steps and
cumulative Krylov iterations. The example declares dimensionless geometry and
loads; its observation integrals and response scales are explicit in the source.
Absent targets, a disclosed synthetic forward design generates them under the
same solve budget. This is not experimental calibration or unique density
identification. Replace the `T..` placeholders with finite numeric target values.

## Programmatic problem

The public types are in `fs_topopt::sdf3::response`. Each `ResponseCase3` carries a
fixed external nodal force, an optional pure prescribed-displacement law on the
retained embedded support, and `ResponseTarget3` observations. An observation is
`q^T u`, with a target, positive normalization scale and nonnegative weight. Use
the actual operator's reference/body integrators to obtain volume or surface
observation vectors. Rebuild these vectors explicitly after a geometry change.

`evaluate_responses` evaluates the sum of weighted, normalized squared response
errors, optionally adding `ResponseOptions3::volume_weight * V`. The derivative
uses an aggregate adjoint per case and **both** terms
`z^T db_g/ds - z^T dK/ds u`, then the existing SIMP/projection/filter pullback.
Nitsche and density-dependent ghost terms are included. Primal and adjoint solves
share one preparation per evaluated density and the same work/cancellation
control. Neither augmented-RHS work nor external-force work is the fitting
objective, an actuator-energy calculation, or a reaction-force measurement.

The evaluator restores incoming material scales after success or returned error;
its complete result retains the evaluated raw/projected design and fields.
`ResponseDesignStudy3` installs only a fully accepted SQP evaluation. It borrows
the study, experiment family and solve control for its lifetime, preventing an
ordinary mesh, target or material substitution beneath cached optimizer state.
Calls to `run` resume the same BFGS/multiplier checkpoint. Failed trials preserve
accepted evidence and spent work, while a repeated line search may spend more.

## Constraints and stopping

The SQP inequalities are `V-volume_cap <= 0`, `density_floor-rho_i <= 0`, and
`rho_i-1 <= 0`, with exact discrete Jacobians. Constraint-violating but physically
admissible starts and intermediate merit-accepted iterates are allowed. Inspect
`constraint_violation()` and the returned KKT residuals before treating an output
as feasible or stationary. An iteration limit or stalled search is not convergence.
The example exits nonzero on a physical/numerical stop or an infeasible endpoint,
after printing any accepted fields/design data it has.

Default dense KKT dimension 256 admits at most 85 cell variables with these
constraints. A hard dimension ceiling of 1024 bounds this small-problem solver
path; it is not a scalable replacement for sparse topology optimization. Setup,
filter, primal and adjoint work remain in `SolveControl`; SQP's callback ceiling
also counts unavailable or rejected trials. Dense QP phases are not preemptible.

Nonzero-motion DWR, shape/follower-load derivatives, geometric changes during an
optimization call, reaction/actuator work, uniqueness and continuum certification
are not implemented by this path. Native Rust compilation/tests remain pending
in the authoring environment; independent numerical references are not Rust runs.
