# Tetrahedral thermal bounds

This is the linear, averaged-functional increment of q61wp.11 / sj31i.25.
It does not close their nonlinear, contact, point-maximum or product-report scope.
The existing 1D speculation APIs and their evidence authorities are unchanged.

## Real consumer

```sh
cargo run -p fs-conduction --features equilibrated-mean --example conduction_mean -- 4 2.0
cargo test -p fs-conduction --features equilibrated-mean --lib verified_mean::
cargo test -p fs-verify --features certified-speculation --lib tet::
```

`conduction::solve_with_mean_bound` consumes the existing `ConductionProblem`
and explicit `MeanSolveConfig`. It returns the original FEM solution/report,
an actual unit-source homogeneous-boundary dual solution/report, and a bound
on whole-domain volume-mean temperature. No second thermal solver is implemented.
The example uses a declared unit slab with unit conductivity, 300 K end faces,
adiabatic sides, and adjustable constant heating; its analytic comparison is
300 + source/12 K. The example emits the actual produced endpoints, not a
precomputed interval. Native test success must be established by execution.

## Mathematical scope

`tet::energy_bound` requires a conforming non-overlapping tetrahedral domain,
element-constant positive scalar conductivity and forcing, affine Dirichlet data,
constant outward Neumann flux, and positive face-constant Robin h with affine
reference data. The candidate must match its Dirichlet trace. Binary64 inputs
are interpreted as their exact real values. The graph solve proposes fluxes;
reverse forest elimination defines an exactly conservative RT0 flux enclosed
by outward arithmetic. Shared face fluxes have opposite signs. Returned flux
boxes are correlated enclosures, not independent conservative degrees of freedom.
The majorant includes both the volume flux defect and Robin trace defect.

For a cell-weighted integral J, `tet::goal_bound` constructs the same-operator
dual with source equal to the exact declared weights. It encloses
`J(v) + R_v(z_h)` with radius `primal_energy_bound * dual_energy_bound`.
The residual correction includes source, Neumann, Robin and volume terms;
it cannot be discarded merely because a solve reported convergence. Both
algebraic and discretization errors remain covered. `tet::mean_bound` uses a
unit-source dual and divides by outward geometric volume only at the end.

## Boundaries

Global mesh non-overlap and fidelity to a CAD surface are caller obligations;
local incidence/orientation checks are not a global geometry certificate.
The conduction adapter refuses unsupported anisotropy, nonlinear k(T), bounded
material-temperature validity without a range proof, nonconstant element source,
nonconstant Neumann normal trace and nonconstant face h. It never averages those
inputs into a different PDE. Matching/mortar contact, nonlinear radiation,
material uncertainty and model discrepancy are outside this increment. A mean
bound does not imply a maximum-temperature bound or a safety decision. No
capability level, persisted verifier authority or evidence colour is promoted.
