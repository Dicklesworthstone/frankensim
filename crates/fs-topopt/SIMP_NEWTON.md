# Filtered density design through constrained Newton

Enable `fs-topopt/newton-design`. `simp_newton::SimpNewtonModel` supplies
compliance and projected physical-volume views of the **same** raw-density
point. Bind either view as an objective or constraint in the existing live IR,
then use `fs_ascent::reverse_newton::AugmentedNewtonStudy`. The native elasticity,
Helmholtz filter, projection, SIMP map and Newton--Krylov/augmented-Lagrangian
engines are reused; no second solver or dense design Hessian is introduced.

```rust,ignore
let physics = SimpNewtonModel::new(
    &pipeline, &elasticity, &loads, &physical_cell_volumes,
    solve_options, lifetime_budget,
)?;
let raw_c = physics.compliance("compliance")?;
let raw_v = physics.volume("physical-volume")?;
let chart = LogisticBox::new(&lower, &upper, Dims::NONE)?;
let compliance = chart.bind(&raw_c)?;
let volume = chart.bind(&raw_v)?;
let initial_logits = chart.encode(&initial_raw_density, Some(&cx))?;
```

`LogisticBox` is the reusable coordinate adapter in
`fs_opt::reverse::physics::box_map`. Each interval's limits carry the original
provider's physical dimensions; the new coordinates are dimensionless. Study
name and output dimensions are preserved. Bind the resulting views to PDE nodes
with matching names, coordinate count and residual dimensions. Compliance is
force times displacement (joules for SI model data); normalize it explicitly in
the IR before combining it with dimensionless volume/constraints.

## Derivatives and shared work

The raw compliance Hessian includes both native stiffness curvature and the
second derivative of the projection/SIMP composition. Volume is weighted by the
provided **physical cell volumes**, not an unweighted mean or SIMP stiffness.
The coordinate adapter adds the complete second chain:

```text
y = lower + (upper-lower) sigmoid(x)
a = dy/dx; b = d2y/dx2
g_x = a * g_y
H_x v = a * H_y(a * v) + b * g_y * v
```

Both scalar views borrow one exact-point cache, immutable physical inputs and
one owner of lifetime work. A fresh point conservatively reserves `loads + 3`
solves: the forward filter, each independent load, and two transpose filters.
Nonzero compliance products reserve `positive-weight loads + 2`; nonlinear
volume products reserve two filter solves. Actual started solves and cumulative
Krylov iterations are reported separately from reservations. Cached samples and
cached zero products do not spend additional work; zero-sharpness raw-volume
curvature needs no solves. A nonlinear coordinate map can still give that
linear raw-volume function nonzero coordinate curvature.

Reservations and measured work survive failures. A failed preparation preserves
the previous point; a successful trial may evict it, and an older accepted tape
re-prepares its own exact point. Cloning an optimizer or a view cannot refund
provider work. Raise the absolute physics ceilings with `set_budget` and resume
from the retained optimizer state. All native solves share the passed `Cx`
cancellation scope; `None` explicitly disables cooperative cancellation.

## Runnable path and regression commands

```bash
cargo run --release -p fs-topopt --features newton-design --example constrained_simp
cargo test --release -p fs-opt --lib reverse::physics::box_map::tests
cargo test --release -p fs-topopt --features newton-design \
  --test compliance_hessian --test compliance_newton --test simp_hessian --test simp_newton
```

The example has 48 independently varying cells on an unequal-volume tetrahedral
mesh, two separately applied loads, and a projected-volume allowance. It prints
the actual termination, numerical KKT residuals, physical work and raw densities,
and checks the candidate with a fresh first-order physical solve. The two
constrained regression problems additionally use a uniform-cell reduction to
check a known inverse-SIMP optimum, in both objective/constraint orientations;
this reduction is test-only, not imposed by the library or example.

Native compilation and these new regressions were **not run** in the authoring
runtime, which had no DSR, RCH, Cargo or rustc. The repository's DSR policy remains
authoritative; these Cargo commands describe a native fallback, not proof of a
passing gate.

## Boundaries

This is fixed-mesh, fixed-load, homogeneous-support linear-elastic SIMP design.
The logistic chart covers an **open** box, not a closed-box active-set algorithm.
Endpoints and numerically saturated coordinates are refused before calling the
physical provider. Chart-coordinate KKT is not a boundary-KKT or global-optimum
certificate. The native Hessian also refuses projected densities outside its
smooth `(1e-12,1)` interval: a box does not waive that restriction, especially
under sharp projection. Keep the initial point and continuation settings in the
resolved smooth domain. Continuation requires rebuilding immutable physics
models at each new parameter setting; this API does not silently mutate a bound
model or transfer its lifetime budget to a new owner.

No shape/stress derivatives, binary-topology guarantee, material calibration,
continuum certificate, global memory admission or process-restart serialization
is implied. Existing first-order-only physics remains first-order-only through
the generic chart; there is no numerical-difference fallback.
