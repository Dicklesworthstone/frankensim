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

## Reciprocal finite-mass impact

`PlateImpact` appends one physical body translation to the complete reduced
spatial plate state. A fixed `PlatePatchStencil` supplies a work-conjugate row
`b`: closure is `z - b*q_plate`, and the same row applies the opposite
contact force to the plate. Constructor admission checks the patch's support
indexing against the bound model, including equally sized but differently
indexed models. Mesh geometry and the fixed footprint remain caller-owned.

The adapter reuses the admitted scalar `fs-dcontact::Obstacle` and its
`ContactStorage` potential, Hunt–Crossley loss, non-attractive unloading cap,
and analytic tangents. Original law gaps, quadrature weights, stiffness,
exponent, internal loss and provenance remain unchanged. Footprint area
distributes a **total** contact law; it does not multiply that law again.
The contact-only inner storage is zero, so plate elastic energy is counted once.

Displacement- and velocity-dependent resistance participates in the shared
implicit generalized-alpha solve, consistent initial acceleration and discrete
adjoints. There is no lagged force injection, modal truncation or second time
integrator. Plate and applied-load parameters retain their original order;
body mass, contact-law controls, footprint and geometry are fixed in these
derivatives. Contact switching and unloading-cap corners use the original
owner's selected one-sided derivatives.

The probe reports plate/body kinetic energy, plate/contact potential energy,
closure and relative velocity, penetration, normal force and separate physical
contact/plate dissipation. These are instantaneous physical quantities.
Generalized-alpha may also dissipate numerically and is not an exact discrete
contact-energy ledger. The footprint is an averaged normal port: this does not
claim local surface conformity, tangential friction, changing contact topology
or rigid-body rotations.

Run `cargo run -p fs-plate-transient --example plate_impact` for a clamped
DKT plate and a moving body. The example prints CSV of both motions, contact
force, penetration and energy/dissipation, using explicitly illustrative SI
inputs. `tests/plate_impact.rs` covers reciprocal power, unchanged law measure,
unloading/separation branches, tangent/transpose finite differences, momentum
transfer and rebound, coupled step derivatives, and binding/budget refusals.

## Source-bound material impact

`PlateImpact::from_source` consumes the existing
`fs-contact::normal_patch::source::SourcePointContact` and retains that exact
owner, including both bulk states, the ordered interface/model card, geometry,
applicability envelope and original receipts. The adapter lowers its original
Hertz coefficient and Hunt–Crossley coefficient into the existing distributed
contact owner with unit point weight. It neither accepts replacement contact
coefficients nor adds a second constitutive formula. `source()` exposes the
retained owner for original pressure/footprint observations and identity checks.

The shared structural time owner calls the new model admission hooks with the
actual timestep, starting state, implicit stage and proposed endpoint. Invalid
stages never reach the contact constitutive callbacks. A failed solve or
inadmissible endpoint leaves displacement, velocity, acceleration, time, step
count and history unchanged. Consistent initialization and adjoints enforce
the same state domain. Direct contact actions, transposed actions, parameter
pullbacks and probes also check the original source's indentation/rate limits.
Instantaneous actions use the source's maximum duration only for a state-only
box check; actual integration duration is checked independently by the driver.

An arbitrary assembled `PlateModel` does not retain its source material
identity. Callers of the general constructor must therefore assemble the plate
from the matching source surface and keep the depth/layer geometry consistent.
The `material_plate_impact` example demonstrates this complete path: it derives
plate Young modulus, Poisson ratio and density from the exact retained
surface-B state, then uses the same source for contact. Its two explicitly
synthetic cases change both plate response and contact response. A material
change requires re-resolution and reassembly; the existing dynamics scale
parameters are not derivatives of the material-card resolver.

Run `cargo run -p fs-plate-transient --example material_plate_impact`.
The example uses the existing material/interface evaluators and retains
`SyntheticFixture` authority. It does not claim experimental calibration,
thermal evolution, three-dimensional contact conformity, or independent
material identity verification for arbitrary preassembled matrices.
