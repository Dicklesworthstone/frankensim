# Tetrahedral thermal bounds

This is the linear, averaged-functional increment of q61wp.11 / sj31i.25.
It does not close their nonlinear, contact, point-maximum or product-report scope.
The existing 1D speculation APIs and their evidence authorities are unchanged.

## Real consumer

```sh
cargo run -p fs-verify --features thermal-conduction --example conduction_mean -- 4 2.0
cargo run -p fs-verify --features thermal-conduction --example conduction_mean -- 4 8.0 4.0 2.0 1.0
cargo test -p fs-verify --features thermal-conduction --lib conduction::
cargo test -p fs-verify --features certified-speculation --lib tet::
```

`fs_conduction::verification::solve_with_mean_bound`, enabled by
`fs-conduction/thermal-verification`, consumes the existing `ConductionProblem`
and explicit `MeanSolveConfig`. It returns the original FEM solution/report,
an actual unit-source homogeneous-boundary dual solution/report, and a bound
on whole-domain volume-mean temperature. No second thermal solver is implemented.
The exact elementwise conductivity tensors, material assignment and boundary
operator are shared by the primal, dual and verifier; no isotropic substitution
is used. The example supports optional directional `kx ky kz` conductivities.
It uses a unit slab, 300 K end faces, adiabatic sides and constant heating;
the analytic comparison is `300 + source/(12*kx)` K. It emits produced interval
endpoints, not a precomputed interval. Native success requires execution.

`fs_conduction::verification::bound_temperature_mean` bounds an existing P1 temperature field
without repeating its primal solve. Pass the original `ConductionProblem`,
the full nodal temperatures, explicit dual `SolveConfig` and `FluxBudget`.
The field need not be converged: the bound includes its remaining algebraic
error. Invalid lengths, non-finite values and a mismatched Dirichlet trace
refuse before the dual solve. The returned `MeanFieldBound` contains the actual
dual report and mean enclosure; it does not invent a primal solver report.

The native adapter lives with the conduction solver; `fs-verify::conduction`
retains the actual fields/reports while verifying their mean. The production
dependency direction is `fs-conduction -> fs-adjoint -> fs-verify` (with the
optional direct verification edge). Only the example and cross-crate tests
consume the native solver through a dev dependency, avoiding a Cargo cycle.

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

`TensorTetProblem` and `tensor_energy_bound`, `tensor_goal_bound`, and
`tensor_mean_bound` extend this same construction to elementwise symmetric
positive-definite 3x3 conductivity tensors in mesh coordinates. Normalized
principal minors establish definiteness with outward arithmetic, and the full
inverse tensor weights the quadratic flux defect. Every off-diagonal term is
retained in the residual and exact simplex moments. Inconclusive definiteness
or unbounded inverse arithmetic refuses rather than clipping, symmetrizing or
replacing the operator. `validate_conductivity` offers the same coefficient
preflight before expensive primal work; it is not a geometry/field certificate.

For a cell-weighted integral J, `tet::goal_bound` constructs the same-operator
dual with source equal to the exact declared weights. It encloses
`J(v) + R_v(z_h)` with radius `primal_energy_bound * dual_energy_bound`.
The residual correction includes source, Neumann, Robin and volume terms;
it cannot be discarded merely because a solve reported convergence. Both
algebraic and discretization errors remain covered. `tet::mean_bound` uses a
unit-source dual and divides by outward geometric volume only at the end.
The tensor counterparts retain these same functional and residual semantics.

## Affine heat sources without changing the PDE

The native APIs above also accept `ScalarField::Nodal`: the original four
source values on each tetrahedron go to `AffineSourceTetProblem`, not a cell
average. Constant rows retain the original verifier arithmetic. The existing
example accepts an explicit source ramp without changing its original invocation:

```sh
cargo run -p fs-verify --features thermal-conduction --example conduction_mean -- \
  --affine 4 0.0 24.0 4.0 2.0 1.0
cargo test -p fs-verify --features certified-speculation --test affine_source
```

Here `f(x)=24*x` W/m^3 and `kx=4` W/(m K), so the exact nominal slab solution
is `T(x)=300+x-x^3` K, whose volume mean is `300.25` K. The example emits both
source endpoints and the actual solver/verifier results; it does not precompute
an enclosure. For endpoints f0,f1, its analytical reference is
`300+(f0+f1)/(24*kx)`. Numerical equality of that mean with an averaged-source
problem does NOT make their temperature fields or error bounds interchangeable.

The extension uses the existing conservative face-flux construction for the
exact cell mean and adds the local polynomial

```text
b = (1/4) sum_i f_i * lambda_i * (x - vertex_i).
div(b) = f - mean(f),    b dot n = 0 on every face.
```

Thus the full flux equilibrates the actual affine source. Its normal continuity
and every Neumann/Robin trace are unchanged. The energy defect is degree four
and is integrated by outward simplex moments, including full tensor-inverse
cross terms. The goal correction uses the exact P1-source/P1-dual product;
using separate averages would shift the interval center. No data-oscillation
term is silently discarded and no auxiliary PDE solve is introduced.

Regression controls include a zero-mean nonzero source whose averaged substitute
has zero forcing, tensor cross terms, manufactured cubic mean/energy bounds,
source-order changes, native primal/dual parity, retained incomplete fields,
constant-source compatibility and original boundary refusals. The authoring
environment had no Rust toolchain: these native tests remain unverified there.
Independent polynomial/quadrature checks are not native execution or an
interval-rounding proof. The existing thermal workflow runs the numerical
sources and integration regressions, then the actual FEM consumer when its
original dependency bootstrap succeeds.

## Boundaries

Global mesh non-overlap and fidelity to a CAD surface are caller obligations;
local incidence/orientation checks are not a global geometry certificate.
The conduction adapter refuses nonlinear k(T), bounded material-temperature
validity without a range proof, nonconstant Neumann normal trace and
nonconstant face h. It never averages those inputs into a
different PDE. Exact tensor symmetry is required; an approximately symmetric
material tensor is not silently repaired. Matching/mortar contact, nonlinear
radiation, material uncertainty and model discrepancy are outside this increment.
A mean bound does not imply a maximum-temperature bound or a safety decision.
No capability level, persisted verifier authority or evidence colour is promoted.
