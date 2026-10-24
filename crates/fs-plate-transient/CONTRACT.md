# fs-plate-transient CONTRACT

Spatial DKT plate trajectories and physical-parameter adjoints through the
shared generalized-alpha integrator in `fs-time`.

## Why a separate crate
Until 2026-09-27 this was `fs-plate`'s optional `transient` feature. Its
`fs-time`/`fs-solver` dependency closed the package cycle
`fs-feec -> fs-couple -> fs-plate -> fs-solver -> fs-feec`: Cargo resolves the
lockfile with every workspace feature, so an optional edge still counts. Every
workspace cargo command failed. Here the dependency points one way only
(`fs-plate-transient -> fs-plate`), and nothing under `fs-couple` depends back.

## Surface
`PlateDynamics`, its budgets, parameters and load traits (`src/transient.rs`,
moved here from `fs-plate` on 2026-09-28). The mathematical contract (M/K/C parameterization, VJP ordering,
cancellation and refusal behaviour) is stated in `fs-plate/CONTRACT.md` under
the transient bullet and is unchanged by the move. The file reaches `fs-plate`
only through the public `PlateModel`.

## Tests
`tests/plate_transient.rs`, `tests/plate_preconditioning.rs` and the
`moving_plate` / `preconditioned_plate` examples live in this crate.
`fs-ascent`'s `plate_calibration` example includes `examples/moving_plate.rs`.

The `preconditioned_plate` example accepts mesh size 8 (147 free DOFs) or 12
(363 free DOFs), comparing physical Jacobi scaling against the default identity
hook under identical short-restart controls. These are measured fixed-mesh,
fixed-step cases, not a mesh-independent convergence or timing guarantee.
`tests/plate_preconditioning.rs` checks the effective M/C/K
diagonal, both bounded larger-mesh solves, a 147-DOF dense-LU endpoint, and a
12-step sampled trajectory with checkpoint replay and all five total physical
parameter derivatives against independent dense differences. Forward and
adjoint solves retain explicit restart/cycle limits; initialization and load
chain rules are included in the gradient.

## Unsafe boundary
None. Workspace `unsafe_code = "deny"`.

## Moving finite-area loads

`moving_load::MovingPlateLoad` implements the original `PlateLoad` trait for
a translating finite footprint and a piecewise-linear physical load history.
Explicit positive area quadrature is located through `fs-plate`'s shared
`PlatePatchStencil`; the reference mesh, supports and material matrices stay
owned by the original spatial model. Three physical load parameters per knot
(Fz, Mx, My, in that order) participate in the existing trajectory adjoint,
including the consistent-initial-acceleration chain rule. Knot times, path
coordinates and quadrature measures are fixed for these derivatives.

Admission bounds the knot/point/mesh scans and checks every knot footprint.
Each later query also locates the actual footprint, so a path leaving a
nonconvex mesh refuses without clipping or partially updating the force/VJP
buffer. The explicitly declared time horizon is closed; extrapolation refuses.
Loads and observations are stateless under retry, replay and out-of-order time
queries. Observed power is physical load times material velocity; a translating
probe over a motionless displaced surface contributes no mechanical power.
This layer does not integrate energy, solve contact, infer complete dynamic
support reactions or certify the continuously swept footprint.

The `moving_patch` example shares the original forward/adjoint experiment with
`moving_plate`. `tests/moving_patch.rs` covers point-quadrature equivalence,
all knot-load VJPs, physical work, failed-query rollback and an independently
assembled dense trajectory/gradient reference. These authored regressions
require native execution for qualification; focused source type checks alone
do not establish whole-trajectory behavior.
