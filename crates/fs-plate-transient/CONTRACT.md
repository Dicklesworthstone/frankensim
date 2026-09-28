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
