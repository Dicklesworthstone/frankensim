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

## Consistent physical inertia

`PlateDynamics::new` admits the complete symmetric positive-definite sparse
mass, including the consistent panel/beam and displacement/slope cross terms.
The existing diagonal mass takes its original linear-time path. Nondiagonal
mass uses the existing sparse LDLT owner once for positive-definiteness
admission; it is neither diagonalized nor replaced by the admission factor.
Symmetry is required within 64 eps times the geometric mean of the two mass
diagonals, admitting assembly rounding but refusing materially nonsymmetric
inertia. Every original stored entry continues through the primal action,
explicit transpose, Rayleigh damping, energy and parameter derivatives.

`new_with_mass_options` accepts the existing `LdltOptions` for explicit pivot
and numeric factor-byte ceilings. Original operator/DOF ceilings apply before
scanning, and cancellation is checked around the synchronous sparse analysis
and factorization. Singular, indefinite and nonfinite mass refuse before a
model is returned. No factorization or temporary matrix is added to time steps.
The allocation-free mass Jacobi preconditioner is exact only for lumped mass;
consistent-mass initialization and its adjoint require adequate explicitly
supplied Krylov controls. This is admission and integration of the chosen
discretization, not an independent continuum accuracy certificate.

`tests/consistent_mass.rs` supplies a manufactured displacement/slope pencil
with cross inertia. Its G1/G3 trajectory and all five physical/load parameter
derivatives are compared with dense LU initialization and stepping (state
tolerance 2e-10 relative with a 1e-3 component scale floor, gradient tolerance
2e-6 relative). A midpoint free-vibration test checks assembled energy at
2e-9 relative and dissipative decay; G4 checks cancellation and mass/factor
refusals. Full crate and measured shell/plate qualification remain separate
from these bounded mathematical regressions.

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
integrator. By default, plate and applied-load parameters retain their original
order and body/contact controls are fixed. Authored models can opt into the
contact-calibration suffix described below. Contact switching and unloading-cap
corners use the original owner's selected one-sided derivatives.

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

## Physical contact calibration

`PlateImpact::with_contact_calibration` enables physical-parameter derivatives
for an explicitly authored impact law. It appends
`[indenter_mass_kg, contact_stiffness, internal_loss, gap_0, ...]` after the
unchanged plate/load prefix. Every gap retains its original point index and SI
units. The time-independent contact controls contribute to the structural
residual; the applied-load VJP writes exact zeros to this suffix.

The mass pullback includes the indenter's actual acceleration. Contact
coefficient and gap pullbacks come from the original `fs-dcontact::Obstacle`
owner's derivative of the complete positive resisting force, including the
unloading cap and original point measures. The shared initialization and
trajectory adjoints therefore include these controls without a second solver
or a finite-difference production gradient. Forward behavior is unchanged by
enabling the derivative coordinates.

The caller rebuilds a model at each new admissible parameter point. Exponent,
collocation, footprint and physical quadrature weights remain fixed. Contact
switching uses the existing one-sided convention; derivatives do not account
for changes to contact topology. Source-bound material models refuse this
opt-in because free contact coefficients would bypass their material/interface
resolution. This API does not differentiate the material-card resolver or
establish that all fitted controls are independently identifiable.

Run `cargo run -p fs-plate-transient --example calibrate_plate_impact` for a
bounded fit of authored contact stiffness to synthetic velocity observations
on a full DKT plate. It uses the shared sampled trajectory adjoint, chains the
consistent initial acceleration and parameterized load, and keeps all other
controls fixed. `tests/impact_calibration.rs` compares the step and sampled
trajectory gradients to finite differences that rebuild the physical model
and its initial acceleration. It also checks the unchanged default parameter
layout, source-bound refusal and reduction of the actual observation loss.

## Simultaneous finite-body contacts

`PlateImpactSet` composes prepared `PlateImpact` bodies on one full spatial
plate. Its state order is `[plate DOFs, body_0, body_1, ...]`. Build every body
from a copy of the same immutable `PlateDynamics`: admission checks the actual
model and load identities, physical plate parameters and cached load parameter
count. Equal-sized independent plates or different load owners cannot silently
be treated as one specimen. The set separately caps body count, combined DOFs,
total original contact quadrature points and all parameter coordinates.

Every original plate mass, stiffness and damping action is evaluated once.
Each body's original contact owner sees that same trial plate state and its
own finite-mass coordinate; only its reciprocal contact contribution is
scattered into the shared plate prefix and corresponding body entry. Contact
forces at different footprints therefore interact through the actual spatial
plate during one existing generalized-alpha/Newton solve. No sequential contact
solver, repeated plate inertia, or independently advanced plate copies are used.

Each body retains its own footprint, original quadrature measures, material
source, loss law, mass and calibration choice. All source step/state limits
participate in joint admission, with the refusing body index retained by the
structural callback. An inadmissible trial or endpoint cannot publish a partial
set of accepted body states. The original shared time owner provides rollback,
cancellation and trajectory replay. Body-to-body contact detection, changing
footprints and tangential/rotational motion remain outside this normal model.

The common plate/load parameter prefix appears once. Enabled authored body
suffixes follow in constructor order; `contact_parameter_range` exposes their
exact ranges, including empty ranges for fixed/source-bound bodies. Full state
transposes preserve the same gather/scatter maps, and contact parameter
pullbacks reuse the original owner. The exact block mass inverse supports the
shared consistent initialization and its pullback. Source coefficients retain
their original immutable ownership and are never silently made fit variables.

`with_body_forcing` attaches independently prescribed normal forces to the
ordered body coordinates through the existing pure time-dependent `PlateLoad`
protocol. Its output has one force in newtons per body, positive in the body's
positive coordinate direction. The shared plate load remains in the plate
prefix. Force-history parameters append after every contact suffix, and
`body_forcing_parameter_range` exposes their range. These controls contribute
only through the forcing VJP; their fixed-state residual partials are zero.
Consistent initialization and checkpointed trajectory pullbacks therefore
include the original force history without duplicating parameter contributions.
Source-owned contact coefficients remain fixed even when excitation is fitted.

The original aggregate parameter cap also bounds these controls. Rebinding a
history replaces its parameter tail without moving contact ranges. Changing
the provider's parameter count after admission, nonfinite forces or pullbacks,
and callbacks that leave output entries unwritten are refused. Each physical
parameter point uses an immutable provider; state-dependent forces belong in
the coupled residual. Applied body work is external work and must be included
when interpreting driven-system energy changes.

`PlateImpactSetProbe` reports plate storage and Rayleigh loss once, sums the
separate body/contact storage and contact dissipation, and preserves ordered
per-body observations. These are instantaneous physical observations rather
than an assertion of exact generalized-alpha step energy conservation.

Run `cargo run -p fs-plate-transient --example material_plate_multi_impact` for
two material-sourced finite bodies on the same DKT plate. The comparison changes
only the first body's initial velocity and reports the resulting second-body
response through shared plate compliance. A third scenario applies separate
smooth force pulses and reports their amplitudes and signed sampled drive work.
`tests/plate_multi_impact.rs` covers singleton equivalence, shared momentum and
coupling, body permutation, local and parameter derivatives, force-history
trajectory pullbacks, original source refusal and aggregate admission.
