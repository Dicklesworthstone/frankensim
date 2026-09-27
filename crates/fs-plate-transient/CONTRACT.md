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
`PlateDynamics`, its budgets, parameters and load traits, re-exported from the
source file `crates/fs-plate/src/transient.rs`, which is compiled here through
`#[path]`. The mathematical contract (M/K/C parameterization, VJP ordering,
cancellation and refusal behaviour) is stated in `fs-plate/CONTRACT.md` under
the transient bullet and is unchanged by the move. The file reaches `fs-plate`
only through the public `PlateModel`.

## Tests
`crates/fs-plate/tests/plate_transient.rs` and the `moving_plate` example stay
in `fs-plate`, which takes this crate as a dev-dependency (a legal cycle).
`fs-ascent`'s `plate_calibration` example and test consume it the same way.

The `preconditioned_plate` example accepts mesh size 8 (147 free DOFs) or 12
(363 free DOFs), comparing physical Jacobi scaling against the default identity
hook under identical short-restart controls. These are measured fixed-mesh,
fixed-step cases, not a mesh-independent convergence or timing guarantee.
`crates/fs-plate/tests/plate_preconditioning.rs` checks the effective M/C/K
diagonal, both bounded larger-mesh solves, a 147-DOF dense-LU endpoint, and a
12-step sampled trajectory with checkpoint replay and all five total physical
parameter derivatives against independent dense differences. Forward and
adjoint solves retain explicit restart/cycle limits; initialization and load
chain rules are included in the gradient.

## Unsafe boundary
None. Workspace `unsafe_code = "deny"`.
